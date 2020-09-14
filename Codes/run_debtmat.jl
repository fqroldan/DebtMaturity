include("debtmat.jl")

function run_debtmat()

	βvec = range(0.94, 0.99, length=11)

	for (jβ, βv) in enumerate(βvec)

		dd = DebtMat(β = βv)
		write("../output.txt", "Created DebtMat with β = $βv")
		equil!(dd)
		write("../output.txt", "Solved DebtMat with β = $βv")
		save("../Output/debtmat$(jj).jld", "dd", dd)
		write("../output.txt", "Saved DebtMat with β = $βv")
	end
end


run_debtmat()