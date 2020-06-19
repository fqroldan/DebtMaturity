using QuantEcon, Interpolations, Optim, PlotlyJS, ColorSchemes, ForwardDiff, LinearAlgebra, Printf

using ORCA

include("reporting_routines.jl")

mutable struct DebtMat{K}
	pars::Dict{Symbol, Number}

	gr::Dict{Symbol, Vector{Float64}}
	grlong::Dict{Symbol, Vector{Int64}}

	agg::Dict{Symbol,Array{Float64,K}}

	P::Matrix{Float64}

	ϕ::Dict{Symbol,Array{Float64,K}}

	vf::Array{Float64,K}
end

function DebtMat(;
	β = 0.96,
	γ = 0.9,
	ψ = 1/3.3,
	Nb = 6,
	Nd = 6,
	Nθ = 7,
	ρθ = 0.7,
	σθ = 0.01,
	τ = 0.3,
	)

	# ψ > 1-τ || throw(error("ψ too low, should be at least (1-τ) = $(1-τ)")) 

	r_star = 1/β - 1
	ρ = 0.2 # Target average maturity of 7 years: ~0.05 at quarterly freq
	κ = ρ + r_star
	pars = Dict(:β=>β, :γ=>γ, :ψ=>ψ, :Nb=>Nb, :Nd=>Nd, :Nθ=>Nθ, :κ=>κ, :ρ=>ρ, :τ=>τ)

	bgrid = range(-0.5,2,length=Nb)
	dgrid = range(-0.5,2,length=Nd)

	mc = tauchen(Nθ, ρθ, σθ, 0, 1)
	θgrid = mc.state_values .+ 0.1
	P = mc.p

	gr = Dict(:b => bgrid, :d => dgrid, :θ => θgrid)
	Jgrid = gridmake(1:Nb, 1:Nd, 1:Nθ)
	grlong = Dict(
		:b=>[Jgrid[js,1] for js in 1:size(Jgrid,1)],
		:d=>[Jgrid[js,2] for js in 1:size(Jgrid,1)],
		:θ=>[Jgrid[js,3] for js in 1:size(Jgrid,1)],
	)

	agg = Dict([sym => ones(Nb, Nd, Nθ) for sym in [:b′, :d′, :Uc, :qb, :qd]])

	ϕ = Dict([sym => ones(Nb, Nd, Nθ)*0.2 for sym in [:c, :n, :b, :d, :τ, :g]])
	vf = zeros(Nb, Nd, Nθ)

	K = length(gr)
	return DebtMat{K}(pars, gr, grlong, agg, P, ϕ, vf)
end

w(dd::DebtMat, gc) = w(gc,dd.pars[:γ])
function w(gc,γ)
	u = 0.0
	gmin = 1e-8
	if gc < gmin
		return w(gmin,γ) + (gc-gmin) * (gmin)^-γ
	end
	if γ != 1
		# u = 1*gc
		u = (gc.^(1-γ))./(1-γ)
	else
		u = log.(gc)
	end
	return u
end

utility(dd::DebtMat,cc,nn) = utility(cc,nn,dd.pars[:γ],dd.pars[:ψ])
function utility(cc, nn, γ, ψ)
	u = 0.0
	cmin = 1e-8
	Lmin = 1e-8

	L = 1 - nn

	if cc < cmin
		return utility(cmin, nn, γ, ψ) + (cc-cmin) * (cmin)^-γ
	end
	if L < Lmin
		return utility(cc, 1-Lmin, γ, ψ) + ψ * (L-Lmin) * (Lmin)^-γ
	end
	if γ != 1
		u = (cc.^(1-γ))./(1-γ) + ψ * (L.^(1-γ))./(1-γ)
	else
		u = log.(cc) .+ ψ * log.(L)
	end
	return u
end

∇u(dd::DebtMat,cc,nn) = ∇u(cc,nn,dd.pars[:γ],dd.pars[:ψ])
∇u(cc,nn,γ,ψ) = ForwardDiff.gradient(x->utility(x[1],x[2],γ,ψ), [cc,nn])

Uc(dd::DebtMat,cc,nn) = ∇u(cc,nn,dd.pars[:γ],dd.pars[:ψ])[1]
Uc(cc,nn,γ,ψ) = cc.^(-γ)

Un(dd::DebtMat,cc,nn) = ∇u(cc,nn,dd.pars[:γ],dd.pars[:ψ])[2]
Un(cc,nn,γ,ψ) = ∇u(cc,nn,γ,ψ)[2]

""" Computes τ from household intratemporal FOC """
labor_wedge(dd::DebtMat,cc,nn) = labor_wedge(cc,nn,dd.pars[:γ],dd.pars[:ψ])
function labor_wedge(cc,nn,γ,ψ)
	Uc, Un = ∇u(cc,nn,γ,ψ)
	return 1 + Un / Uc
end

n(dd::DebtMat, k) = n(k, dd.pars[:γ], dd.pars[:ψ])
n(k, γ, ψ) = (-k/ψ).^(1/(1-γ))

P(dd::DebtMat,B,D) = P(B,D,dd.pars[:κ])
P(B,D,κ) = B + κ*D

R(dd::DebtMat, B′,D′,D,qb,qd) = R(B′,D′,D,qb,qd,dd.pars[:ρ])
R(B′,D′,D,qb,qd,ρ) = qb*B′ + qd*(D′-ρ*D)


function iterate_q(dd::DebtMat, itp_U, itp_qd)
	β = dd.pars[:β]
	Jgrid = dd.grlong
	
	qd = similar(dd.agg[:qd])
	qb = similar(dd.agg[:qb])

	for js in 1:size(dd.grlong[:b],1)
		jθ = first([dd.grlong[sym][js] for sym in [:θ]])

		pθ = dd.P[jθ,:]

		uc1 = dd.agg[:Uc][js]

		bpv = dd.agg[:b′][js]
		dpv = dd.agg[:d′][js]

		Eqd = 0.0
		Eqb = 0.0
		for (jθp, θpv) in enumerate(dd.gr[:θ])
			prob = pθ[jθp]
			uc2 = itp_U(bpv, dpv, θpv)
			repayment = dd.pars[:κ] + (1-dd.pars[:ρ]) * itp_qd(bpv, dpv, θpv)

			Eqd += prob * uc2 * repayment
			Eqb += prob * uc2
		end

		qd[js] = β * Eqd / uc1
		qb[js] = β * Eqb / uc1
	end

	return qd, qb
end

function update_q!(dd::DebtMat; maxiter::Int64=2500, tol::Float64=1e-8, verbose::Bool=false)
	knots = (dd.gr[:b], dd.gr[:d], dd.gr[:θ])
	itp_U = interpolate(knots, dd.agg[:Uc], Gridded(Linear()))

	iter, dist = 0, 1+tol
	t0 = time()
	while iter < maxiter && dist > tol
		iter += 1

		old_q = copy(dd.agg[:qd])
		itp_qd = interpolate(knots, old_q, Gridded(Linear()))

		dd.agg[:qd], dd.agg[:qb] = iterate_q(dd, itp_U, itp_qd)

		dist = sum( (old_q - dd.agg[:qd]).^2 ) / sum( old_q.^2 )
	end

	# dd.agg[:Qd] = dd.agg[:Uc] .* dd.agg[:qd]
	# dd.agg[:Qb] = dd.agg[:Uc] .* dd.agg[:qb]

	message = "Done in $iter iterations ($(time_print(time()-t0)))"
	if dist < tol
		message *= " ✓"
	else
		message *= ". Dist = $(@sprintf("%0.3g", dist))"
	end
	if verbose
		println(message)
	end
	return (dist < tol)
end

resource_constraint_g(dd::DebtMat, cc, nc) = nc - cc

function budget_constraint(dd::DebtMat, cv, bpv, dpv, τv, state, pθ, itp_U, itp_qd)
	γ, ψ, κ, ρ = [dd.pars[key] for key in [:γ, :ψ, :κ, :ρ]]
	bv, dv = [state[key] for key in [:b, :d]]

	nv = 1 - cv * ( ψ / (1-τv) )^(1/γ)
	# nv = max(1e-6, nv)
	gv = resource_constraint_g(dd, cv, nv)
	# gv = max(0, gv)

	Ucv = max(1e-10,cv)^(-γ)
	qbv, qdv = debt_price(dd, bpv, dpv, pθ, itp_U, itp_qd, Ucv)

	LHS = gv + bv + κ * dv

	RHS = τv * nv + qbv*bpv + qdv * (dpv - (1-ρ)*dv)

	return LHS, RHS, cv, nv, gv, qbv, qdv
end

function budget_surplus(dd::DebtMat, cv, bpv, dpv, τv, state, pθ, itp_U, itp_qd)
	LHS, RHS, cv, nv, gv, qbv, qdv = budget_constraint(dd, cv, bpv, dpv, τv, state, pθ, itp_U, itp_qd)

	return RHS - LHS
end

function eval_vf(dd::DebtMat, τv, bpv, dpv, state, pθ, itp_U, itp_qd, itp_v)
	β, ψ, τ, γ = [dd.pars[sym] for sym in [:β, :ψ, :τ, :γ]]
	θv = state[:θ]

	max_c = 1.0

	res = Optim.optimize(cv -> budget_surplus(dd, cv, bpv, dpv, τv, state, pθ, itp_U, itp_qd)^2, 0, max_c)

	cc = res.minimizer

	_, _, cv, nv, gv, qbv, qdv = budget_constraint(dd, cc, bpv, dpv, τv, state, pθ, itp_U, itp_qd)

	ut = utility(dd, cc, nv) + θv * w(dd, gv)
	
	Eu = 0.0
	for (jθp, θpv) in enumerate(dd.gr[:θ])
		Eu += pθ[jθp] * itp_v(bpv, dpv, θpv)
	end
	return ut + β * Eu
end

n_FOC(dd::DebtMat, τv, cc) = (dd.pars[:ψ] / (1-τv))^(1/dd.pars[:γ]) * cc

# function budget_constraint(dd::DebtMat, cc, bpv, dpv, state, qbv, qdv)
# 	bv, dv, θv = [state[sym] for sym in [:b, :d, :θ]]

# 	γ, ρ, κ, ψ = [dd.pars[sym] for sym in [:γ, :ρ, :κ, :ψ]]

# 	RHS = 1 + (qbv*bpv + qdv * (dpv - (1-ρ)*dv) - bv - κ*dv) / cc

# 	RHS = max(1e-8, RHS)

# 	nv = ψ * RHS^(1/(1-γ)) * cc
# 	# nv = max(nv, cc+1e-2)

	
# 	τv = 1 - ψ * (cc / nv)^γ
# 	# τv = max(0, min(1, τv))
# 	# nv = max(0.1, min(nv, 0.9))

# 	return τv, nv
# end

function debt_price(dd::DebtMat, bp, dp, pθ, itp_U, itp_qd, Ucv)
	κ, ρ, β = [dd.pars[sym] for sym in [:κ, :ρ, :β]]

	qb, qd = zeros(2)
	for (jθp, θpv) in enumerate(dd.gr[:θ])
		prob = pθ[jθp]
		
		qb += prob * itp_U(bp, dp, θpv)
		qd += prob * itp_U(bp, dp, θpv) * (κ + (1-ρ) * itp_qd(bp, dp, θpv))
	end

	qb *= β / Ucv
	qd *= β / Ucv

	return qb, qd
end

function solve_opt_value(dd::DebtMat, guess::Dict{Symbol, Float64}, state::Dict{Symbol, Float64}, pθ, itp_U, itp_qd, itp_v)
	bv, dv, θv = [state[sym] for sym in [:b, :d, :θ]]

	bmin, bmax = minimum(dd.gr[:b])+1e-4, maximum(dd.gr[:b])-1e-4
	dmin, dmax = minimum(dd.gr[:d])+1e-4, maximum(dd.gr[:d])-1e-4
	τmin, τmax = 0.01, min(1,1-dd.pars[:ψ])
	
	wrap_vf(x) = -eval_vf(dd, x[1], x[2], x[3], state, pθ, itp_U, itp_qd, itp_v)

	res = Optim.optimize(wrap_vf, [τmin, bmin, dmin], [τmax, bmax, dmax], [guess[sym] for sym in [:τ, :b, :d]], Fminbox(NelderMead()))

	if Optim.converged(res)
	else
		println("COULDN'T FIND OPTIMUM AT b=$bv, d=$dv, θ=$θv")
	end
	τv, bp, dp = res.minimizer
	vf = first(-wrap_vf(res.minimizer))


	res_inner = Optim.optimize(cv -> budget_surplus(dd, cv, bp, dp, τv, state, pθ, itp_U, itp_qd)^2, 0, 1.0)
	cc = res_inner.minimizer
	_, _, cv, nv, gv, qbv, qdv = budget_constraint(dd, cc, bp, dp, τv, state, pθ, itp_U, itp_qd)

	ϕv = Dict(:c=>cv, :n=>nv, :b=>bp, :d=>dp, :τ=>τv, :g=>gv)
	return ϕv, vf
end

function optim_step(dd::DebtMat, itp_U, itp_qd, itp_v)
	Jgrid = gridmake(1:dd.pars[:Nb], 1:dd.pars[:Nd], 1:dd.pars[:Nθ])

	new_v = similar(dd.vf)
	new_ϕ = Dict(key => similar(val) for (key,val) in dd.ϕ)

	for js in 1:size(Jgrid,1)
		state = Dict(sym => dd.gr[sym][dd.grlong[sym][js]] for sym in [:b, :d, :θ])
		
		jθ = dd.grlong[:θ][js]
		pθ = dd.P[jθ,:]

		guess = Dict(key=>dd.ϕ[key][js] for key in keys(dd.ϕ))
		ϕv, vv = solve_opt_value(dd, guess, state, pθ, itp_U, itp_qd, itp_v)

		for sym in keys(new_ϕ)
			new_ϕ[sym][js] = ϕv[sym]
		end
		new_v[js] = vv
	end
	return new_v, new_ϕ
end

function vfi_iter(dd::DebtMat, itp_qd)
	knots = (dd.gr[:b], dd.gr[:d], dd.gr[:θ])
	itp_v = interpolate(knots, dd.vf, Gridded(Linear()));
	itp_U = interpolate(knots, dd.agg[:Uc], Gridded(Linear()));

	new_v, new_ϕ = optim_step(dd, itp_U, itp_qd, itp_v)

	return new_v, new_ϕ
end

function update_agg!(dd::DebtMat)
	dd.agg[:Uc] = dd.ϕ[:c].^(-dd.pars[:γ])
	dd.agg[:b′] = dd.ϕ[:b]
	dd.agg[:d′] = dd.ϕ[:d]
	nothing
end

function vfi!(dd::DebtMat; maxiter::Int64=250, tol=Float64=1e-6)
	iter, dist = 0, 1+tol
	
	knots = (dd.gr[:b], dd.gr[:d], dd.gr[:θ]);
	itp_qd = interpolate(knots, dd.agg[:qd], Gridded(Linear()));

	upd_η = 0.2

	t0 = time()
	while iter < maxiter && dist > tol
		iter += 1
		println("Iteration $iter")
		old_v = copy(dd.vf)

		dd.vf, dd.ϕ = vfi_iter(dd, itp_qd);

		update_agg!(dd)

		dist = sum( (old_v - dd.vf).^2 ) / sum( old_v.^2 )
		println("dist_v = $dist at ||v|| = $(norm(old_v))")
	end
	return dist
end

function equil!(dd::DebtMat; maxiter::Int64=250, tol::Float64=1e-4)
	iter, dist = 0, 1+tol

	t0 = time()
	while iter < maxiter && dist > tol
		iter += 1
		println("Outer iteration $iter")
		old_qs = Dict(sym => copy(dd.agg[sym]) for sym in [:qd, :qb])

		update_q!(dd)

		dist_v = vfi!(dd, maxiter = 2)

		dist_q = maximum([sum((old_qs[key] - dd.agg[key]).^2) / sum( old_qs[key].^2 ) for key in keys(old_qs)])
		println("dist_q = $dist at ||q|| = $(norm(dd.agg[:qd]))")

		dist = max(dist_q, 10dist_v)

		println()
	end
end

