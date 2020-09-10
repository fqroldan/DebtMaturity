mutable struct DebtMat{K}
	pars::Dict{Symbol, Number}

	gr::Dict{Symbol, Vector{Float64}}
	grlong::Dict{Symbol, Vector{Int64}}

	agg::Dict{Symbol,Array{Float64,K}}

	P::Matrix{Float64}

	ϕ::Dict{Symbol,Array{Float64,K}}

	vf::Array{Float64,K}
end

abstract type AbstractPath
end
mutable struct Path{T} <: AbstractPath
	data::Dict{Symbol, Vector{Float64}}
end
function Path(; T::Int64 = 1)
	data = Dict( key => Vector{Float64}(undef, T) for key in [:b, :d, :θ, :jθ, :c, :n, :g, :qb, :qd])
	return Path{T}(data)
end

periods(pv::Vector{T}) where T <: AbstractPath = sum([periods(pp) for pp in pv])
periods(p::Path{T}) where T = T

function check_periods(p::Path, t::Int64)
	0 < t <= periods(p) || throw("t out of bounds")
	nothing
end

getfrompath(p::Path, t::Int64, sym::Symbol) = p.data[sym][t]
getfrompath(p::Path, t::AbstractArray, sym::Symbol) = [p.data[sym][tv] for tv in t]
getfrompath(p::Path, t::Int) = Dict(key => p.data[key][t] for key in keys(p.data))
getfrompath(p::Path, sym::Symbol) = p.data[sym]
series(p::Path, sym::Symbol) = getfrompath(p,sym)
getmean(p::Path, sym::Symbol) = getmean([p], sym)
getmean(pv::Vector{T}, sym::Symbol) where T <: AbstractPath = mean(vcat([series(pp, sym) for pp in pv]...))

function fill_path!(p::Path, t::Int64, d::Dict=Dict())
	check_periods(p,t)
	missing_keys = 0
	for (key, val) in d
		if haskey(p.data, key)
			p.data[key][t] = val
		else
			missing_keys += 1
		end
	end

	if missing_keys > 0
		print("WARNING: $missing_keys missing keys")
	end
	nothing
end

function trim_path(p::Path{T}, t0::Int64) where T
	check_periods(p,t0)

	return Path{T-t0}(Dict(key => val[t0+1:end] for (key, val) in p.data))
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
	# ρ = 0.2 # Target average maturity of 7 years: ~0.05 at quarterly freq
	ρ = 0.0 # For a consol
	κ = ρ + r_star
	pars = Dict(:β=>β, :γ=>γ, :ψ=>ψ, :Nb=>Nb, :Nd=>Nd, :Nθ=>Nθ, :κ=>κ, :ρ=>ρ, :τ=>τ)

	bgrid = range(-0.2,0.4,length=Nb)
	dgrid = range(-0.2,0.4,length=Nd)

	mc = tauchen(Nθ, ρθ, σθ, 0, 1)
	θgrid = mc.state_values .+ 0.1*1.3
	P = mc.p

	gr = Dict(:b => bgrid, :d => dgrid, :θ => θgrid)
	Jgrid = gridmake(1:Nb, 1:Nd, 1:Nθ)
	grlong = Dict(
		:b=>[Jgrid[js,1] for js in 1:size(Jgrid,1)],
		:d=>[Jgrid[js,2] for js in 1:size(Jgrid,1)],
		:θ=>[Jgrid[js,3] for js in 1:size(Jgrid,1)],
	)

	agg = Dict([sym => ones(Nb, Nd, Nθ)*mean(bgrid) for sym in [:b′, :d′, :Uc, :qb, :qd]])

	ϕ = Dict([sym => ones(Nb, Nd, Nθ)*0.2 for sym in [:c, :n, :b, :d, :τ, :g]])
	vf = zeros(Nb, Nd, Nθ)

	K = length(gr)
	return DebtMat{K}(pars, gr, grlong, agg, P, ϕ, vf)
end
