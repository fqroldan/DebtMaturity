""" Define styles """

def_style = let
	axis = attr(showgrid = true, gridcolor="#e2e2e2", gridwidth=0.5, zeroline=false)
	layout = Layout(xaxis = axis, yaxis=axis)
	Style(layout=layout)
end

slides_def = let
	layout = Layout(plot_bgcolor="#fafafa", paper_bgcolor="#fafafa",
		width=1920*0.45, height=1080*0.45, font_size=16, font_family="Lato",
		legend = attr(orientation = "h", x=0.05))
	Style(def_style, layout=layout)
end

dark_bg = let
	axis = attr(gridcolor="#1b1b1b")
	layout = Layout(plot_bgcolor="#020202", paper_bgcolor="#020202", font_color="white", xaxis=axis,yaxis=axis)
	Style(layout=layout)
end
slides_dark = Style(slides_def, dark_bg)

paper = let
	layout = Layout(width = 1920 * 0.5, height = 1080 * 0.35, font_size=16, font_family = "Linux Libertine",
		legend = attr(orientation = "h", x=0.05))
	Style(def_style, layout=layout)
end

function time_print(tfloat::Float64)
	t = floor(Int, tfloat)
	if t < 1
		t_print = "no time"
	elseif t < 60
		t_print = "$(t) second"
		t == 1 ? nothing : t_print = t_print * "s"
	elseif t < 3600
		t_print = "$(floor(Int,t/60)) minute"
		floor(Int,t/60) == 1 ? nothing : t_print = t_print * "s"
		if t % 60 != 0
			t_print = t_print * " and $(t%60) second"
			t % 60 == 1 ? nothing : t_print = t_print * "s"
		end
	else
		t_print = "$(floor(Int,t/3600)) hour"
		floor(Int, t/3600) == 1 ? nothing : t_print = t_print * "s"
		t = t % 3600
		t_print = t_print * " and $(floor(Int,t/60)) minute"
		floor(Int,t/60) == 1 ? nothing : t_print = t_print * "s"
	end
	return t_print
end

function makeplot_simul(dd::DebtMat, pp::Path; style::Style=slides_def)
	T = periods(pp)
	data = [
	scatter(x=1:T, y=series(pp,:b), name="<i>b", yaxis="y1")
	scatter(x=1:T, y=series(pp,:d), name="<i>d", yaxis="y1")
	scatter(x=1:T, y=series(pp,:c), name="<i>c", yaxis="y2")
	# scatter(x=1:T, y=series(pp,:n), name="<i>n", yaxis="y2")
	scatter(x=1:T, y=series(pp,:g), name="<i>g", yaxis="y3")
	scatter(x=1:T, y=series(pp,:θ), name="<i>θ", yaxis="y4")
	]

	layout = Layout(
	yaxis1 = attr(domain = [0.76, 1]),
	yaxis2 = attr(domain = [0.51, 0.74]),
	yaxis3 = attr(domain = [0.26, 0.49]),
	yaxis4 = attr(domain = [0, 0.24]),
	xaxis = attr(anchor = "y4"),
	)

	plot(data, layout, style=style)
end
