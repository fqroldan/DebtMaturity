include("debtmat.jl")

function run_debtmat()

	βvec = range(0.94, 0.99, length=11)
	
	dd = DebtMat()

	for (jβ, βv) in enumerate(βvec)

		dd.pars[:β] = βv
		write("../output.txt", "Created DebtMat with β = $βv")

		flag = equil!(dd)
		write("../output.txt", "Solved DebtMat with β = $βv $(ifelse(flag, "✓","NO"))")

		if flag
			save("../Output/debtmat$(jβ).jld", "dd", dd)
			write("../output.txt", "Saved DebtMat with β = $βv")
		end
	end
end


run_debtmat()