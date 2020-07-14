using MAT, PlotlyJS, ORCA, CSV, DataFrames, ColorSchemes

""" DEFINE STYLES """
def_style = let
	axis = attr(gridcolor="#e2e2e2", gridwidth=0.5, zeroline=false)
	layout = Layout(xaxis = axis, yaxis=axis)
	Style(layout=layout)
end

slides_def = let
	layout = Layout(plot_bgcolor="#fafafa", paper_bgcolor="#fafafa",
		width=1920*0.45, height=1080*0.45, font_size=18, font_family="Lato",
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


function makeplot(psi)
	colpal = ColorSchemes.lajolla
	file = matopen("../../hpc/tenedorLoop_psi$(floor(Int,psi)).mat");

	M = read(file)["BOPTFINAL"];

	b01, b02, delta = M[1,:], M[2,:], M[3,:];

	data = [scatter(x=delta, y=b01, line_width=3, marker_color=get(colpal, 0.25), name="b<sub>0</sub><sup>1")
       scatter(x=delta, y=b02, line_width=3, marker_color=get(colpal, 0.75), name="b<sub>0</sub><sup>2")]

    layout = Layout(title="Debt at time <i>0</i> for <i>ψ = $(psi)", xaxis_title="<i>δ", yaxis_range=[-0.025, 0.475], height=1080*.35)

	plot(data, layout, style=slides_dark)
end

function makeCAs()

	colpal = ColorSchemes.romaO

	df = CSV.read("../../Data/CAs.csv")

	df.:group_1 = [findfirst(unique(df.:group).==df.group[jj]) for jj in 1:length(df.:group)]
	rename!(df, :group_1=>:gr)

	fiscal = [sum(df[(df.gr.==jg) .& (df.year.==year),:].:fiscal_gr) for jg in 1:length(unique(df.gr)), (jt, year) in enumerate(unique(df.year))]

	countrynames = ["USA", "DEU", "Euro Area Debtors", "CHN", "JPN", "Oil Exporters", "AEs", "EMs"]

	data = [
		[scatter(yaxis = "y1", xaxis = "x1", x=df[df.gr.==jg,:].year, y=fiscal[jg,:], showlegend=true, mode="lines", marker_color=get(colpal, (jg-1)/(length(countrynames))), name=countrynames[jg]) for jg in 1:length(unique(df.gr))]
		# [scatter(yaxis = "y1", x=df[(df.gr.==jg).&(df.surp.==1),:].year, y=df[(df.gr.==jg).&(df.surp.==1),:].fiscal_gr, showlegend=false, marker_color=get(colpal, (jg-1)/(length(countrynames))), name=countrynames[jg]) for jg in 1:length(unique(df.gr))]
		[scatter(yaxis = "y2", xaxis="x2",x=df[(df.gr.==jg).&(df.surp.==0),:].year, y=df[(df.gr.==jg).&(df.surp.==0),:].bca_gdp_gr*100, showlegend=false, mode="lines", marker_color=get(colpal, (jg-1)/(length(countrynames))), name=countrynames[jg]) for jg in 1:length(unique(df.gr))]
		[scatter(yaxis = "y2", xaxis="x2",x=df[(df.gr.==jg).&(df.surp.==1),:].year, y=df[(df.gr.==jg).&(df.surp.==1),:].bca_gdp_gr*100, showlegend=false, marker_color=get(colpal, (jg-1)/(length(countrynames))), name=countrynames[jg]) for jg in 1:length(unique(df.gr))]
	]

	annotations = [
		attr(text="Cyclically-adjusted fiscal balance", x = 0.225, xanchor="center", y=1, xref = "paper", yref="paper", showarrow=false, font_size=20)
		attr(text="Current account", x = 0.775, xanchor="center", y=1, xref = "paper", yref="paper", showarrow=false, font_size=20)
	]

	layout = Layout(
		xaxis1 = attr(domain=[0, 0.45]),
		xaxis2 = attr(domain=[0.55, 1]),
		yaxis1 = attr(zeroline=true, anchor="x1", title="% of group GDP"),
		yaxis2 = attr(zeroline=true, anchor="x2", title="% of group GDP"),
		width=1920*0.6, height=1080*0.6,
		annotations = annotations,
		title="Fiscal balance and current accounts (WEO Live)"
	)


	plot(data, layout, style=slides_dark)
end