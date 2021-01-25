using QuantEcon, Interpolations, Optim, PlotlyJS, ColorSchemes, ForwardDiff, LinearAlgebra, Printf, Random, JLD, Dates, ProgressBars

# using ORCA

include("type_def.jl")
include("reporting_routines.jl")
include("simul.jl")

w(dd::DebtMat, gc) = w(gc,dd.pars[:γ])
function w(gc,γ)
	u = 0.0
	gmin = 5e-2
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
	cmin = 5e-2
	Lmin = 1e-3

	L = 1 - nn

	if cc < cmin
		return utility(cmin, nn, γ, ψ) + (cc-cmin) * (cmin)^-γ
	end
	if L < Lmin
		return utility(cc, Lmin, γ, ψ) + ψ * (L-Lmin) * (Lmin)^-γ
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


function iterate_q!(qd, qb, dd::DebtMat, itp_U, itp_qd)
	β = dd.pars[:β]
	Jgrid = dd.grlong

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
	# println(extrema(qd))
end

function update_q!(dd::DebtMat; maxiter::Int64=2500, tol::Float64=1e-8, verbose::Bool=false)
	knots = (dd.gr[:b], dd.gr[:d], dd.gr[:θ])
	itp_U = interpolate(knots, dd.agg[:Uc], Gridded(Linear()))

	qd = similar(dd.agg[:qd])
	qb = similar(dd.agg[:qb])

	upd_η = 0.75

	iter, dist = 0, 1+tol
	t0 = time()
	while iter < maxiter && dist > tol
		iter += 1

		itp_qd = interpolate(knots, dd.agg[:qd], Gridded(Linear()))

		iterate_q!(qd, qb, dd, itp_U, itp_qd)

		dist = norm( qd - dd.agg[:qd] ) / (1+norm(dd.agg[:qd]))

		dd.agg[:qd] = dd.agg[:qd] * (1-upd_η) + upd_η * qd
		dd.agg[:qb] = dd.agg[:qb] * (1-upd_η) + upd_η * qb
	end

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

function debt_price(dd::DebtMat, bp, dp, pθ, itp_U, itp_qd)
	κ, ρ, β = (dd.pars[sym] for sym in [:κ, :ρ, :β])

	qb, qd = zeros(2)
	for (jθp, θpv) in enumerate(dd.gr[:θ])
		prob = pθ[jθp]

		qb += prob * itp_U(bp, dp, θpv)
		qd += prob * itp_U(bp, dp, θpv) * (κ + (1-ρ) * itp_qd(bp, dp, θpv))
	end

	qb *= β
	qd *= β

	return qb, qd
end

function debt_price(dd::DebtMat, bp, dp, pθ, itp_U, itp_qd, Ucv)
	Eub, Eud = debt_price(dd, bp, dp, pθ, itp_U, itp_qd)

	return Eub / Ucv, Eud / Ucv
end	

function post_tax_wage(dd, cv, bpv, dpv, state, pθ, qbv::Number, qdv::Number)
	ρ, κ = (dd.pars[sym] for sym in [:ρ, :κ])
	bv, dv = (state[key] for key in [:b, :d])

	# wage_priv = n * (1-τ)
	wage_priv = qbv * bpv + qdv * (dpv - (1-ρ) * dv) - bv - κ * dv + cv
end

function IC_given_cbd(dd::DebtMat, cv, bpv, dpv, state, pθ, itp_U, itp_qd)
	Ucv = Uc(dd, cv, 0.0)
	qbv, qdv = debt_price(dd, bpv, dpv, pθ, itp_U, itp_qd, Ucv)

	IC_given_cbd(dd, cv, bpv, dpv, state, pθ, qbv, qdv)
end

function IC_given_cbd(dd::DebtMat, cv, bpv, dpv, state, pθ, qbv::Number, qdv::Number)
	ψ, γ = (dd.pars[sym] for sym in [:ψ, :γ])
	
	wage_priv = post_tax_wage(dd, cv, bpv, dpv, state, pθ, qbv, qdv)
	τv = 1 - wage_priv + ψ * cv
	# τv = max(0, min(1-1e-4, τv))
	nv = 1 - ψ/(1-τv) * cv

	# nv = max(1e-4, nv)
	gv = nv - cv

	return τv, nv, gv
end

function eval_vf(dd::DebtMat, cv, bpv, dpv, state, pθ, Eucb::Number, Euqd::Number, itp_v)
	β, ψ, γ, ρ, κ = (dd.pars[sym] for sym in [:β, :ψ, :γ, :ρ, :κ])
	θv = state[:θ]

	Ucv = Uc(dd, cv, 0.0)
	qbv = Eucb / Ucv
	qdv = Euqd / Ucv

	τv, nv, gv = IC_given_cbd(dd, cv, bpv, dpv, state, pθ, qbv, qdv)

	ut = utility(dd, cv, nv) + θv * w(dd, gv)

	Eu = 0.0
	for (jθp, θpv) in enumerate(dd.gr[:θ])
		Eu += pθ[jθp] * itp_v(bpv, dpv, θpv)
	end
	return ut + β * Eu
end

function find_C_bounds(dd, bpv, dpv, state, pθ, itp_U, itp_qd)

	bv = state[:b]
	dv = state[:d]

	κ, ρ, β = (dd.pars[key] for key in [:κ, :ρ, :β])

	Eucb, Euqd = debt_price(dd, bpv, dpv, pθ, itp_U, itp_qd)

	numer = bv + κ * dv
	denom = 1 + β * (Eucb + Euqd)

	return 1e-4, numer / denom, Eucb, Euqd
end

function eval_vf_c(dd, bpv, dpv, state, pθ, itp_U, itp_qd, itp_v, save_r)

	cmin, cmax, Eucb, Euqd = find_C_bounds(dd, bpv, dpv, state, pθ, itp_U, itp_qd)

	if cmax > cmin
		res = Optim.optimize(cv -> -eval_vf(dd, cv, bpv, dpv, state, pθ, Eucb, Euqd, itp_v), cmin, cmax, GoldenSection())

		cv = res.minimizer
		save_r[1] = cv
		save_r[2] = bpv
		save_r[3] = dpv
		vf = res.minimum
	else
		vf = 1e3 * (cmax - cmin) - 1e4
	end
	return vf
end

function solve_opt_value(dd::DebtMat, guess::Dict{Symbol, Float64}, state::Dict{Symbol, Float64}, pθ, itp_U, itp_qd, itp_v)

	bmin, bmax = minimum(dd.gr[:b])+1e-4, maximum(dd.gr[:b])-1e-4
	dmin, dmax = minimum(dd.gr[:d])+1e-4, maximum(dd.gr[:d])-1e-4

	save_r = [guess[key] for key in [:c, :b, :d]]

	wrap_vf(x) = -eval_vf_c(dd, x[1], x[2], state, pθ, itp_U, itp_qd, itp_v, save_r)

	res = Optim.optimize(wrap_vf, [bmin, dmin], [bmax, dmax], [guess[sym] for sym in [:b, :d]], Fminbox(BFGS()))

	cv, bp, dp = save_r

	if Optim.converged(res)
	else
		bv, dv, θv = (state[sym] for sym in [:b, :d, :θ])
		println("COULDN'T FIND OPTIMUM AT b=$bv, d=$dv, θ=$θv")
	end
	# bp, dp = res.minimizer
	vf = first(-wrap_vf(res.minimizer))


	# cmin, cmax = 5e-4, 1 - 5e-4
	# res_c = Optim.optimize(cv -> -eval_vf(dd, cv, bp, dp, state, pθ, itp_U, itp_qd, itp_v), cmin, cmax, GoldenSection())
	# cv = res_c.minimizer

	τv, nv, gv = IC_given_cbd(dd, cv, bp, dp, state, pθ, itp_U, itp_qd)

	ϕv = Dict(:c=>cv, :n=>nv, :b=>bp, :d=>dp, :τ=>τv, :g=>gv)
	return ϕv, vf
end

function optim_step!(new_v, new_ϕ, dd::DebtMat, itp_U, itp_qd, itp_v)
	Jgrid = gridmake(1:dd.pars[:Nb], 1:dd.pars[:Nd], 1:dd.pars[:Nθ]);

	# Threads.@threads for js in 1:size(Jgrid,1)
	for js in ProgressBar(1:size(Jgrid,1))
	# for js in 1:size(Jgrid, 1)
		# println(js)
		state = Dict(sym => dd.gr[sym][dd.grlong[sym][js]] for sym in [:b, :d, :θ])

		jθ = dd.grlong[:θ][js]
		pθ = dd.P[jθ,:]

		# guess = Dict(key=>dd.ϕ[key][js] for key in keys(dd.ϕ))
		guess = Dict(:c => 0.5, :b => maximum(dd.gr[:b])*.9 + mean(dd.gr[:b])*.1, :d=>maximum(dd.gr[:d])*.9 + mean(dd.gr[:d])*.1)
		ϕv, vv = solve_opt_value(dd, guess, state, pθ, itp_U, itp_qd, itp_v)

		for sym in keys(new_ϕ)
			new_ϕ[sym][js] = ϕv[sym]
		end
		new_v[js] = vv
	end
end

function vfi_iter(dd::DebtMat, itp_qd)
	knots = (dd.gr[:b], dd.gr[:d], dd.gr[:θ])
	itp_v = interpolate(knots, dd.vf, Gridded(Linear()));
	itp_U = interpolate(knots, dd.agg[:Uc], Gridded(Linear()));

	new_v = similar(dd.vf);
	new_ϕ = Dict(key => similar(val) for (key,val) in dd.ϕ);

	optim_step!(new_v, new_ϕ, dd, itp_U, itp_qd, itp_v)

	return new_v, new_ϕ
end

function update_agg!(dd::DebtMat; upd_η = 0.5)
	dd.agg[:Uc] = [Uc(dd, cv, 0.0) for cv in dd.ϕ[:c]] * upd_η + (1-upd_η) * dd.agg[:Uc]
	dd.agg[:b′] = dd.ϕ[:b] * upd_η + (1-upd_η) * dd.agg[:b′]
	dd.agg[:d′] = dd.ϕ[:d] * upd_η + (1-upd_η) * dd.agg[:d′]
	nothing

	dd.agg[:Uc] = max.(dd.agg[:Uc], 100)
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

		dist = norm( old_v - dd.vf ) / max(1, norm( old_v ))
		println("dist_v = $(@sprintf("%.3g", dist)) at ‖v‖ = $(@sprintf("%.3g", norm(old_v)))")
	end
	return dist
end

function equil!(dd::DebtMat; maxiter::Int64=250, tol::Float64=1e-4)
	iter, dist = 0, 1+tol

	t0 = time()
	tol_vfi = 1e-2
	while iter < maxiter && dist > tol
		iter += 1
		print("\nOuter iteration $iter at $(Dates.format(now(), "HH:MM"))\n")

		old_qs = Dict(sym => copy(dd.agg[sym]) for sym in [:qd, :qb])
		
		update_q!(dd)

		dist_v = vfi!(dd, tol = tol_vfi)
		update_agg!(dd)

		update_q!(dd)

		if maximum(isnan.(dd.agg[:qb])) == 1 || maximum(isnan.(dd.agg[:qd])) == 1 #|| max(maximum(dd.agg[:qb]), maximum(dd.agg[:qd])) > 1e10
			dd.agg[:qb] = ones(size(dd.agg[:qb]))
			dd.agg[:qd] = ones(size(dd.agg[:qd]))

			dist_q = 1.0
		else
			dist_q = maximum([norm(old_qs[key] - dd.agg[key]) / norm( old_qs[key] ) for key in keys(old_qs)])
		end

		println("dist_q = $(@sprintf("%.3g", dist_q)) at ‖q‖ = $(@sprintf("%.3g", norm(dd.agg[:qd])))")

		dist = max(dist_q, dist_v)
		tol_vfi = max(1e-6, tol_vfi * 0.95)
	end
	s = "Done in $(time_print(time()-t0))"
	dist <= tol ? s *= " ✓" : s *= ". Failed to converge."
	return dist <= tol
end
