function iter_simul!(pp::Path{T}, dd::DebtMat, t, itp_ϕb, itp_ϕd, itp_ϕτ, itp_ϕc, itp_U, itp_qd, itp_pθ) where T

	bt = getfrompath(pp, t, :b)
	dt = getfrompath(pp, t, :d)
	jθ = Int(getfrompath(pp, t, :jθ))

	θt = dd.gr[:θ][jθ]

	pθ = itp_pθ(θt, dd.gr[:θ])
	pθ = pθ / sum(pθ)

	state = Dict(:b=>bt, :d=>dt, :θ=>θt)

	cc = itp_ϕc(bt, dt, θt)
	bp = itp_ϕb(bt, dt, θt)
	dp = itp_ϕd(bt, dt, θt)


	# _, _, cv, nv, gv, qbv, qdv = budget_constraint(dd, cc, bp, dp, τp, state, pθ, itp_U, itp_qd)
	τv, nv, gv = IC_given_cbd(dd, cc, bp, dp, state, pθ, itp_U, itp_qd)

	Ucv = Uc(dd, cc, nv)
	qb, qd = debt_price(dd, bp, dp, pθ, itp_U, itp_qd, Ucv)

	# Fill values of equilibrium at t
	fill_path!(pp, t, Dict(:c => cc, :n => nv, :g => gv, :qb => qb, :qd => qd))

	# Draw shocks for t+1
	probθ = cumsum(pθ)
	jθp = findfirst(probθ .> rand())
	θp = dd.gr[:θ][jθp]

	if t < T
		fill_path!(pp, t+1, Dict(:jθ=>jθp, :θ=>θp, :b=>bp, :d=>dp))
	end
	nothing
end

function simul(dd::DebtMat, simul_length=10000, burn_in=1000)
	Random.seed!(1)
	T = simul_length + burn_in

	pp = Path(T=T)

	knots = (dd.gr[:b], dd.gr[:d], dd.gr[:θ])
	itp_ϕb = interpolate(knots, dd.ϕ[:b], Gridded(Linear()))
	itp_ϕd = interpolate(knots, dd.ϕ[:d], Gridded(Linear()))
	itp_ϕτ = interpolate(knots, dd.ϕ[:τ], Gridded(Linear()))
	itp_ϕc = interpolate(knots, dd.ϕ[:c], Gridded(Linear()))
	itp_U  = interpolate(knots, dd.agg[:Uc], Gridded(Linear()))
	itp_qd = interpolate(knots, dd.agg[:qd], Gridded(Linear()))
	itp_pθ = interpolate((dd.gr[:θ], dd.gr[:θ]), dd.P, Gridded(Linear()))

	fill_path!(pp, 1, Dict(:jθ=>1, :θ=>dd.gr[:θ][1], :b =>  2*mean(dd.gr[:b]), :d=> 2*mean(dd.gr[:d])))

	for jt in 1:T
		iter_simul!(pp, dd, jt, itp_ϕb, itp_ϕd, itp_ϕτ, itp_ϕc, itp_U, itp_qd, itp_pθ)
	end
	return trim_path(pp, burn_in)
	# pp
end