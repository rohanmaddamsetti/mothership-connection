### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 9c4cf388-bacd-11ed-3a30-0f0349e5f60c
begin
	using Plots
	using Plots.PlotMeasures
	using PlutoUI
	using DifferentialEquations
	using LinearAlgebra
	using Images
	using StatsPlots
	using DataFrames
	using ColorSchemes
	using CSV
end

# ╔═╡ 067221ad-e5d7-48fa-a6b1-6a89ffc33b87
md"""
# **_Toy model of ARG duplication and transposition dynamics_**

**linear ODE model version 6: julia 1.8+**.

Rohan Maddamsetti, Teng Wang, Yi Yao, and Lingchong You.
"""


# ╔═╡ 8706c369-4696-4edb-bc80-199d3fd01af3
load("../results/diagrams/toy-model-v3.jpg")

# ╔═╡ 9e545059-b2c5-4134-b977-94984b01dc66
md"""
## **Model description**

I built a toy model to illustrate why multiple identical copies of a protein sequence in a genome may reveal recent positive selection. A diagram of the model is shown above.

There are three subpopulations of bacteria. Each cell contains a chromosome and a multi-copy plasmid. Each chromosome and plasmid may contain an antibiotic resistance gene (ARG). An ARG on a chromosome is shown as a large black bar, and an ARG on the plasmid is shown as a small black bar.

We are interested in the dynamics of the three subpopulations due to growth and mutation (duplication, loss, and transposition dynamics of the ARG). I roughly follow the modeling framework used by Lopatkin et al. (2017) "Persistence and reversal of plasmid-mediated antibiotic resistance" and by Yao et al. (2022) "Intra-and interpopulation transposition of mobile genetic elements driven by antibiotic selection" and described in those papers' Supplementary Information.

**Model Assumptions**

*Growth dynamics*. We assume that there is a steady inflow of nutrients and antibiotic, and a steady outflow of depleted media and cells, reflected by a constant dilution rate, $D$. This assumption allows the population to grow continuously at a steady-state population size. We normalize the number of cells by the carrying capacity, such that each state variable represents the percentage of carrying capacity that is taken up by the subpopulation-- note that this is *not* the relative frequency of cells in the population, because the total population may be at a steady state that is less than carrying capacity. The growth rate of each subpopulation is modeled by  growth functions $f_i > 0$, that we describe in greater detail below.

*Mutation dynamics.* We define a mutation as a transition from one state to another due to gene duplication by transposition. Each transition occurs at a constant rate $\eta$. We assume that transposon excision rates rates are negligible, such that duplication events leave the original copy unchanged in the chromosome.

These assumptions lead to a system of differential equations of the form:

$\frac{dx_i}{dt} = f_ix_i(1 - \Sigma x_i) - Dx_i + Q_i$ where the first term reflects logistic growth at a rate $f_i$ when carrying capacity has not been reached, the second term reflects constant dilution due to a fixed outflow rate, and the third term wraps up all the state transitions (mutation dynamics). 


**Model Equations**

$\frac{dx_1}{dt} = f_1x_1 (1 - \Sigma x_i) - Dx_1 -  2\eta x_1$

$\frac{dx_2}{dt} = f_2x_2 (1 - \Sigma x_i) - Dx_2 + \eta x_1$

$\frac{dx_3}{dt} = f_3x_3 (1 - \Sigma x_i) - Dx_3 + \eta x_1$


**Growth functions**

$f_i = (1-c)^x \frac{K_j^n}{K_j^n + A^n}$ where $A$ is antibiotic concentration and $K_i$ is the concentration of antibiotic that reduces growth by 50%, $n$ is a Hill coefficient, $c$ is the cost of expressing the ARG, and $x$ is the physical number of ARGs in the cell.

We assume that the plasmid has a copy number of y, with values ranging from 0 to 4. We assume a Hill cofficient $n = 3$. We also assume that $0 < c < 1$, and that $A > 0$. $K$ varies depending on the configuration of genes on chromosome or plasmid:

$f_1 = (1-c) \frac{1^3}{1^3 + A^3}$

$f_2 = (1-c)^2 \frac{2^3}{2^3 + A^3}$

$f_3 = (1-c)^{(1+y)} \frac{(1+y)^3}{(1+y)^3 + A^3}$
"""

# ╔═╡ 1cfeb2d1-e9f9-43cb-a700-cc758944cdc1
function calc_f(Aₜ, k, n, c)
		(1 - c)^k * k^n/(k^n + Aₜ^n)
end

# ╔═╡ 5943c034-659a-4ea7-8781-ac22d13d0baf
md""" My code follows this tutorial: 
[https://diffeq.sciml.ai/stable/tutorials/ode_example](https://diffeq.sciml.ai/stable/tutorials/ode_example)
"""

# ╔═╡ 757af6dd-ed81-48d4-b703-bad3f5d77e71
function dynamics!(du, u, p, t)	
		x1, x2, x3 = u
		xtotal = sum(u)

		η, antibiotic_conc_func, c, D, plasmid_copy_number = p
		
		k1, k2, k3 = [1 2 (1+plasmid_copy_number)]
		n = 3 # Hill coefficient
	
		f1 = calc_f(antibiotic_conc_func(t), k1, n, c)
		f2 = calc_f(antibiotic_conc_func(t), k2, n, c)
		f3 = calc_f(antibiotic_conc_func(t), k3, n, c)
		
    	du[1] = dx1 = f1*x1*(1 - xtotal) - D*x1 - 2η*x1
		du[2] = dx2 = f2*x2*(1 - xtotal) - D*x2 + η*x1
		du[3] = dx3 = f3*x3*(1 - xtotal) - D*x3 + η*x1
end

# ╔═╡ 31940f75-496f-4f5d-a027-f6bc82f9bea3
md""" 
For Figure 1, I use the following parameter settings:  

Antibiotic Concentration = 2.0

Duplication Cost = 0.1

Transposition Rate = 0.0002  

Dilution Rate = 0.1.

Plasmid copy number = 2.
"""

# ╔═╡ d787a932-5cbe-459d-a193-a14261c46eb5
md""" Antibiotic Concentration Slider"""

# ╔═╡ 6ee60bf2-c926-47ba-a3b7-8f3b7fd9c3e6
@bind AntibioticConcentration Slider(0:0.01:5, default=2.0, show_value=true)

# ╔═╡ e53719fe-5949-458e-94b8-ecca22ad05b9
md""" Duplication Cost Slider"""

# ╔═╡ c7b81a78-9045-4004-ac11-aa6420df7dea
@bind DuplicationCost Slider(0:0.01:0.5, default=0.1, show_value=true)

# ╔═╡ 26e1c2e3-063c-4116-a0ac-b07fd6efb896
md""" Transposition Rate Slider"""

# ╔═╡ 1589175b-fa55-4592-b8a8-36e7cd239101
@bind TranspositionRate Slider(0:0.00001:0.0002, default=0.0002, show_value=true)

# ╔═╡ 73848f0b-6f2e-4aba-92e9-97ebaa97f5c6
md""" Dilution Rate Slider"""

# ╔═╡ 0a59dd38-1478-4e3b-bde8-514aafc112c6
@bind DilutionRate Slider(0:0.01:0.1, default=0.1, show_value=true)

# ╔═╡ d32f0e11-b46b-4474-b617-2f0c33128b0d
md""" Plasmid copy number Slider"""

# ╔═╡ af18e81a-c1c8-4965-ba85-8c973608980e
@bind PlasmidCopyNumber Slider(0:1:4, default=2, show_value=true)

# ╔═╡ 3e58b2fc-6052-4602-8991-a4935e97e4a9
begin	
	## initial conditions
	x1, x2, x3 = 1, 0, 0
	u₀ = [x1, x2, x3]
	## time interval
	tspan = (0.0,200.0)
end

# ╔═╡ fabbfa13-8b8e-4268-a48a-7ecac577bbb4
begin 

	## transposition rate.
	η = TranspositionRate
	
	## Dilution rate.
	D = DilutionRate

	## Plasmid copy number
	y = PlasmidCopyNumber
	
	## Antibiotic concentration as a function of time.
	##A₀ = t->AntibioticConcentration ## constant function
	A₀ = t->AntibioticConcentration
	Apulse = t->ifelse(t<5000, AntibioticConcentration,0) ## step function
	
	## fitness cost of duplication c: can be anywhere between 0 and 1.
	c = DuplicationCost
	
	match_Dynetica = false
	if (match_Dynetica)
		η = 3e-5
		D = 0.05
		c = 0.1
		A₀ = t->4
		Apulse = t->ifelse(t<200, 4,0)
	end
	
	## bundle parameters into a vector.
	antibiotic_treatment = [η, A₀, c, D, y]
	pulse_antibiotic_treatment = [η, Apulse, c, D, y]
end

# ╔═╡ 2eec564a-522d-4c0a-b038-cfd5e92125d0
md""" ### plots of the fitness functions. 

fitnesses of x1, x2, x3. 

"""

# ╔═╡ 1f6bd6e1-0750-4624-8b39-e6f5616f8e0e
let ## local scope block
	
	k1, k2, k3 = [1 2 (1+PlasmidCopyNumber)]
	n = 3 # Hill coefficient
	
	f1 = calc_f(AntibioticConcentration, k1, n, c)
	f2 = calc_f(AntibioticConcentration, k2, n, c)
	f3 = calc_f(AntibioticConcentration, k3, n, c)
	fitnesses = [f1, f2, f3]
	
	bar(fitnesses,legend=:bottomright)
	
end

# ╔═╡ ac2b4292-c087-449b-a414-d193b007b3cd
function ConcAndCostToXTotal(antibiotic_conc::Float64, my_cost::Float64, fixed_parameters)
	
	η, D, plasmid_copy_number = fixed_parameters
	antibiotic_conc_func = t->antibiotic_conc

	my_parameters = [η, antibiotic_conc_func, my_cost, D, plasmid_copy_number]
	
	my_prob = ODEProblem(dynamics!, u₀, tspan, my_parameters)
	my_sol = solve(my_prob)
	return sum(my_sol[end])
end

# ╔═╡ 823bad9a-5c39-4492-9164-3e8e74f1ef22
let
	fixed_parameters = [η, D, y]
	
	p = plot()
	for cost in 0:0.01:0.08
	antibiotic_concs = [x for x in 0:0.02:4]
	xtotals = [ConcAndCostToXTotal(x, cost, fixed_parameters) for x in antibiotic_concs]
	my_label = "cost = $cost"
	plot!(antibiotic_concs, xtotals, label = my_label,
			legend = false, palette = :YlGnBu_9,
			ylabel = "Total population size",
			xlabel="Antibiotic Concentration", grid = false,
			fontfamily="Helvetica")
	end
	p
end

# ╔═╡ 9d277a7a-8677-45e8-a8c2-a31561ade308
antibiotic_prob = ODEProblem(dynamics!, u₀, tspan, antibiotic_treatment);	

# ╔═╡ 5cc61696-e79b-4426-baee-99b5df4c0ce7
antibiotic_sol = solve(antibiotic_prob);

# ╔═╡ 46f89e44-2a2a-4b18-81d0-2f4c93db2fc2
let
	Fig1B = plot(antibiotic_sol,linewidth=2,xaxis="Time", size=(3.5*72,3*72),
					legend = false, fontfamily = "Helvetica", grid = false,
					yaxis = "Biomass")
	savefig(Fig1B, "../results/linear-ODE-model-figures/Fig1B-pop-dynamics.pdf")
	## save Source Data for Fig1B.
	CSV.write("../results/Source-Data/Fig1B-Source-Data.csv", antibiotic_sol, header=["Time", "TypeI", "TypeII","TypeIII"])
	Fig1B
end

# ╔═╡ 3ada8460-59d2-4650-8170-e1330be1182c
antibiotic_sol_array = transpose(hcat(antibiotic_sol.u...))

# ╔═╡ 832c8291-e35d-4b4f-a97c-8302193a2975
let
	## IMPORTANT NOTE: THE TIME UNITS ARE MESSED UP! DON'T USE THIS FIGURE!
	## This is *probably* because the timestep is not a fixed interval of time,
	## but I haven't checked this to confirm.
oldFig1B = groupedbar(antibiotic_sol_array, bar_position = :stack, size=(3.5*72,3*72),
					bar_width = 1, legend = false, fontfamily = "Helvetica", lw = 0,
					grid = false,
					ylabel = "Total biomass",
					xlabel = "Time")
savefig(oldFig1B, "../results/linear-ODE-model-figures/old-Fig1B-pop-dynamics.pdf")
oldFig1B
end

# ╔═╡ a5035e2d-b3d8-4f50-942f-ab6ec2aa91e6
pulse_antibiotic_prob = ODEProblem(dynamics!, u₀, tspan, pulse_antibiotic_treatment);

# ╔═╡ 1fbc5140-008d-42d3-a7ed-5980bfdf77a0
pulse_antibiotic_sol = solve(pulse_antibiotic_prob);

# ╔═╡ d2c5f306-d784-47eb-bf15-6c9848979e06
begin
p1 = plot([Apulse(t) for t in 1:10000], linewidth=1,label="Antibiotic")
p2 = plot(pulse_antibiotic_sol,linewidth=2,xaxis="t")
p3 = plot(p1, p2, layout = (2, 1))
p3
end

# ╔═╡ 682f67a5-2bb8-46a9-be45-96c108b2c2d7
savefig(p3, "../results/linear-ODE-model-figures/toy-model-dynamics-v0.5.pdf")

# ╔═╡ 427ae3f1-3555-42dc-9481-4e40f044d064
md""" __Duplication Index calculation__

DI  (duplication index) = $(x_2 + x_3)/ x_{total}$, that is, the fraction of the gene that will be duplicated (for a certain ARG cost).

"""

# ╔═╡ 5c3d5841-d7bb-4d29-a717-509937e046df
function FinalDuplicationIndex(sol)
	## fraction of population containing duplicates
	return (sol[end][2] + sol[end][3])/sum(sol[end])
end

# ╔═╡ 9babc592-5fc3-4840-b231-c93a26e7c54d
function ConcAndCostToFinalDuplicationIndex(antibiotic_conc::Float64, my_cost::Float64, fixed_parameters)
	
	η, D, y = fixed_parameters
	antibiotic_conc_func = t->antibiotic_conc

	my_parameters = [η, antibiotic_conc_func, my_cost, D, y]
	
	my_prob = ODEProblem(dynamics!, u₀, tspan, my_parameters)
	my_sol = solve(my_prob)
	return FinalDuplicationIndex(my_sol)
end

# ╔═╡ d4f12532-dd2b-4fb8-997b-d8b2740ac710
let
	fixed_parameters = [η, D, y]

	p = plot(size=(3.5*72,3*72))

	## Create arrays to store the Source Data for this figure.
	antibiotic_concs_data = Float64[]
	dup_indices_data = Float64[]
	cost_values_data = Float64[]
	
	for cost in 0.05:0.05:0.25
		antibiotic_concs = [x for x in 0.25:0.001:1.2]
		dup_indices = [ConcAndCostToFinalDuplicationIndex(x, cost, fixed_parameters) for x in antibiotic_concs]

		## Append data to the arrays for Source Data.
    	append!(antibiotic_concs_data, antibiotic_concs)
    	append!(dup_indices_data, dup_indices)
    	append!(cost_values_data, fill(cost, length(antibiotic_concs)))
		
		my_label = "cost = $cost"
		plot!(antibiotic_concs, dup_indices, label=my_label,
			legend = false,
			ylabel="Duplication Index",
			palette = :Hokusai3,
			xlabel="Antibiotic Concentration",
			fontfamily="Helvetica",
			grid = false)
	end

	## Create a DataFrame to store the Source Data
	Fig1C_source_data_df = DataFrame(
    	Antibiotic_Concentration = antibiotic_concs_data,
    	Duplication_Index = dup_indices_data,
    	Cost = cost_values_data
	)
	## Save Fig1C Source Data to file.
	CSV.write("../results/Source-Data/Fig1C-Source-Data.csv", Fig1C_source_data_df)
	
	savefig(p, "../results/linear-ODE-model-figures/Fig1C-DI-versus-selection.pdf")
	p
end

# ╔═╡ 0fb01ad1-3e36-42ee-93f2-ee8714cd4909
md"""Fig. 1D: Keep fitness cost and [A] constant. Show the time trajectories of the Duplication index for different rates of duplication."""

# ╔═╡ 817f5ea0-16e2-417b-8bbf-b163974b9c6c
function DuplicationIndexOverTime(sol)
	## fraction of population containing duplicates over time
	duplication_index_vec = [(v[2] + v[3])/sum(v) for v in sol.u]
	return duplication_index_vec
end

# ╔═╡ 4f859c18-0a3a-48c1-8fc6-e3d70324ddb0
struct DuplicationIndexTimeSeries
    t::Vector{Float64}
    v::Vector{Float64}
end

# ╔═╡ 625b835c-7bb1-47ef-93ab-f954741a55b3
function TimeSeriesDuplicationIndex(antibiotic_conc::Float64, my_cost::Float64, η, D, y)
	antibiotic_conc_func = t->antibiotic_conc
	my_parameters = [η, antibiotic_conc_func, my_cost, D, y]
	my_prob = ODEProblem(dynamics!, u₀, tspan, my_parameters)
	my_sol = solve(my_prob) 
	DI_TimeSeries = DuplicationIndexTimeSeries(my_sol.t, DuplicationIndexOverTime(my_sol))
	return DI_TimeSeries
end

# ╔═╡ 022c37b3-6761-4958-8a43-5434dcccdedc
let
	cost = 0.1
	antibiotic_conc = 2.0
	fig1D = plot(size=(3.5*72,3*72))

	## Create arrays to store the Source Data for this figure.
	time_data = Float64[]
	duplication_index_data = Float64[]
	duplication_rate_data = Float64[]
	
	for duplication_rate in [0, 2e-7, 2e-6, 2e-5, 2e-4]
		DI_timeseries = TimeSeriesDuplicationIndex(antibiotic_conc, cost, duplication_rate, D, y)

		## Append data to the arrays for Source Data.
    	append!(time_data, DI_timeseries.t)
    	append!(duplication_index_data, DI_timeseries.v)
    	append!(duplication_rate_data, fill(duplication_rate,length(DI_timeseries.t)))
		
		my_label = "cost = $cost"
		plot!(DI_timeseries.t, DI_timeseries.v, label=my_label,
			legend = false,
			ylabel="Duplication Index",
			palette = :Hokusai3, #:tol_bright,
			xlabel="Time",
			fontfamily="Helvetica",
			grid = false)
	end

	## Create a DataFrame to store the Source Data
	Fig1D_source_data_df = DataFrame(
    	Time = time_data,
   		Duplication_Index = duplication_index_data,
    	Duplication_Rate = duplication_rate_data
	)
	## Save Fig1D Source Data to file.
	CSV.write("../results/Source-Data/Fig1D-Source-Data.csv", Fig1D_source_data_df)
	
	savefig(fig1D, "../results/linear-ODE-model-figures/Fig1D-DI-versus-TranspositionRate.pdf")
	fig1D
end

# ╔═╡ c67c9b24-3fad-44d5-a0a4-69716a5b8aad
md"""Figure 1E: Keep fitness cost and [A] constant. Vary both [A] and transposition rate and show the final Duplication index."""

# ╔═╡ c65feae3-8b17-435d-9103-21c8ce054937
function make_DI_matrix_entry(antibiotic_conc, transposition_rate)
	## helper function for making the heatmap.
	DI_timeseries = TimeSeriesDuplicationIndex(antibiotic_conc, c, transposition_rate, D, y)
	final_DI = DI_timeseries.v[end]
	return final_DI
end

# ╔═╡ d120286a-a043-4080-9030-9c599ad11684
let

	fig1E = plot(size=(3.5*72,3*72))
	
	antibiotic_concs_vec = [x for x in 0.0:0.1:1.2]
	transposition_rates_vec = [1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4]

	## fill the DI_matrix.
	## IMPORTANT: Julia is column-major, so antibiotics concentrations are columns (y-axis) while transposition rates are rows (x-axis).
	DI_matrix = [make_DI_matrix_entry(i, j) for i in antibiotic_concs_vec, j in transposition_rates_vec]

	## Create arrays to store the Source Data for Figure 1E.
	antibiotic_conc_column = Float64[]
	transposition_rate_column = Float64[]
	DI_column = Float64[]

	for i in antibiotic_concs_vec
		for j in transposition_rates_vec
			final_DI = make_DI_matrix_entry(i, j)
			## Append data to the arrays for Source Data.
    		append!(antibiotic_conc_column, i)
    		append!(transposition_rate_column, j)
    		append!(DI_column, final_DI)
		end
	end

	## Create a DataFrame to store the Source Data
	Fig1E_source_data_df = DataFrame(
    	AntibioticConcentration = antibiotic_conc_column,
		TranspositionRate = transposition_rate_column,
   		Duplication_Index = DI_column
	)
	## Save Fig1E Source Data to file.
	CSV.write("../results/Source-Data/Fig1E-Source-Data.csv", Fig1E_source_data_df)

	## Make Figure 1E.
	fig1E = heatmap(
		DI_matrix, 
		ylabel="Antibiotic Concentration", 
		xlabel="log(Transposition Rate)", 
		title="Duplication Index", 
		fontfamily="Helvetica", 
		c=:viridis,
		bottom_margin = 3mm,
		colorbar=true
	)
	
	# Set tick mark labels on the X-axis
	xticks!(fig1E, 1:4:9, ["-12", "-8", "-4"])

	# Set tick mark labels on the Y-axis
	yticks!(fig1E, 1:6:13, ["0","0.6","1.2"])
	
	# Customize font size
	plot!(fig1E, titlefontsize=16, xlabelfontsize=14, ylabelfontsize=14, tickfontsize=14)
	
	savefig(fig1E, "../results/linear-ODE-model-figures/Fig1E-DI-heatmap.pdf")
	fig1E
	
end

# ╔═╡ 36934771-b7a6-482e-b316-85ec429c5574
findcolorscheme("cvd")

# ╔═╡ Cell order:
# ╠═9c4cf388-bacd-11ed-3a30-0f0349e5f60c
# ╠═067221ad-e5d7-48fa-a6b1-6a89ffc33b87
# ╠═8706c369-4696-4edb-bc80-199d3fd01af3
# ╠═9e545059-b2c5-4134-b977-94984b01dc66
# ╠═1cfeb2d1-e9f9-43cb-a700-cc758944cdc1
# ╟─5943c034-659a-4ea7-8781-ac22d13d0baf
# ╠═757af6dd-ed81-48d4-b703-bad3f5d77e71
# ╟─31940f75-496f-4f5d-a027-f6bc82f9bea3
# ╠═d787a932-5cbe-459d-a193-a14261c46eb5
# ╠═6ee60bf2-c926-47ba-a3b7-8f3b7fd9c3e6
# ╠═e53719fe-5949-458e-94b8-ecca22ad05b9
# ╠═c7b81a78-9045-4004-ac11-aa6420df7dea
# ╠═26e1c2e3-063c-4116-a0ac-b07fd6efb896
# ╠═1589175b-fa55-4592-b8a8-36e7cd239101
# ╠═73848f0b-6f2e-4aba-92e9-97ebaa97f5c6
# ╠═0a59dd38-1478-4e3b-bde8-514aafc112c6
# ╠═d32f0e11-b46b-4474-b617-2f0c33128b0d
# ╠═af18e81a-c1c8-4965-ba85-8c973608980e
# ╠═3e58b2fc-6052-4602-8991-a4935e97e4a9
# ╠═fabbfa13-8b8e-4268-a48a-7ecac577bbb4
# ╠═2eec564a-522d-4c0a-b038-cfd5e92125d0
# ╠═1f6bd6e1-0750-4624-8b39-e6f5616f8e0e
# ╠═ac2b4292-c087-449b-a414-d193b007b3cd
# ╠═823bad9a-5c39-4492-9164-3e8e74f1ef22
# ╠═9d277a7a-8677-45e8-a8c2-a31561ade308
# ╠═5cc61696-e79b-4426-baee-99b5df4c0ce7
# ╠═46f89e44-2a2a-4b18-81d0-2f4c93db2fc2
# ╠═3ada8460-59d2-4650-8170-e1330be1182c
# ╠═832c8291-e35d-4b4f-a97c-8302193a2975
# ╠═a5035e2d-b3d8-4f50-942f-ab6ec2aa91e6
# ╠═1fbc5140-008d-42d3-a7ed-5980bfdf77a0
# ╠═d2c5f306-d784-47eb-bf15-6c9848979e06
# ╠═682f67a5-2bb8-46a9-be45-96c108b2c2d7
# ╠═427ae3f1-3555-42dc-9481-4e40f044d064
# ╠═5c3d5841-d7bb-4d29-a717-509937e046df
# ╠═9babc592-5fc3-4840-b231-c93a26e7c54d
# ╠═d4f12532-dd2b-4fb8-997b-d8b2740ac710
# ╠═0fb01ad1-3e36-42ee-93f2-ee8714cd4909
# ╠═817f5ea0-16e2-417b-8bbf-b163974b9c6c
# ╠═4f859c18-0a3a-48c1-8fc6-e3d70324ddb0
# ╠═625b835c-7bb1-47ef-93ab-f954741a55b3
# ╠═022c37b3-6761-4958-8a43-5434dcccdedc
# ╠═c67c9b24-3fad-44d5-a0a4-69716a5b8aad
# ╠═c65feae3-8b17-435d-9103-21c8ce054937
# ╠═d120286a-a043-4080-9030-9c599ad11684
# ╠═36934771-b7a6-482e-b316-85ec429c5574
