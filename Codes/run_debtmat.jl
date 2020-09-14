include("debtmat.jl")

function run_debtmat()

	βvec = range(0.94, 0.99, length=11)

	for (jβ, βv) in enumerate(βvec)

		dd = DebtMat(β = βv)

		equil!(dd)

		save("../Output/debtmat$(jj).jld", "dd", dd)
	end
end


run_debtmat()