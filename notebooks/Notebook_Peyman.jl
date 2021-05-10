### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 3d583f2a-9fce-4409-8f0d-83716b908bad
using Plots, PlutoUI, Distributions

# ╔═╡ 65a16bbf-a6e9-4199-9b4b-7bf3c99b5f38
md""" 

Problem setup: We want to estimate the mean $\mu$ of a gaussian distribution 		with known variance $\sigma^2$, from one observed data point $x_0$

Observed data x₀ $(@bind Obs Slider(LinRange(-5,5,101), show_value=true, default=0)) 

Variance σ² $(@bind σ² Slider(LinRange(1,5,51), show_value=true, default=1)) 
"""

# ╔═╡ edcec37a-a2ea-4104-baf4-c76a0cb26cd1
begin
		function likelihood(θ, x)
			return exp(-(θ-x)^2/(2*σ²))
		end
		
		function xminmax(α)
			s = sqrt(-2*σ²*log(α))
			return [Obs-s, Obs+s]	
		end
	x = LinRange(-12,12,200)
end;

# ╔═╡ b8d27d20-e226-4ed5-8721-d84ae962acc6
md""" $\alpha$ and $\beta$ are related to one another and control the risk in our prediction 

α $(@bind αp Slider(LinRange(0.01,0.99,51), show_value=true, default=0.3))

Choose β $(@bind bul CheckBox())

β $(@bind βp Slider(LinRange(0.01,0.97,51), show_value=true, default=0.3))

""" 

# ╔═╡ 8c704d70-a318-4763-be57-43117f9b3b19
begin
	get_β(alf) = 1-cdf(Chisq(1), -2*log(alf))
	get_α(bet) = exp(quantile(Chisq(1), 1-bet)/(-2))
	α = ifelse(bul, get_α(βp), αp)
	β = ifelse(!bul, get_β(αp), βp)
end;

# ╔═╡ 61958cf2-c345-4344-ab01-46d9b88d9017
begin
	function risk(β)
	return σ²*quantile(Chisq(1), 1-β)
	#return -8*σ²*log(get_α(β))  	
	end
	p1 = plot(LinRange(0, 1, 101), get_β.(LinRange(0, 1, 101)) , xlabel = "α", ylabel = "β")
	scatter!([α], [β])
	p2 =plot(x, likelihood.(x, Obs))
	v = xminmax(α)
	plot!(v, [α, α], lw =1)
	vline!(v)
	#p3 = plot(LinRange(0.01, 0.99, 201),risk_vals, xlabel = "β", ylabel = "risk")
	p3 = plot(LinRange(0.01, 0.99, 201),risk.(LinRange(0.01, 0.99, 201)),xlabel = "β", ylabel = "risk")

	scatter!([β], [risk(β)])
	plot(p1, p2, p3, layout = (3,1), size = (400, 600), legend = false)
end

# ╔═╡ Cell order:
# ╠═3d583f2a-9fce-4409-8f0d-83716b908bad
# ╟─edcec37a-a2ea-4104-baf4-c76a0cb26cd1
# ╟─65a16bbf-a6e9-4199-9b4b-7bf3c99b5f38
# ╟─b8d27d20-e226-4ed5-8721-d84ae962acc6
# ╟─61958cf2-c345-4344-ab01-46d9b88d9017
# ╟─8c704d70-a318-4763-be57-43117f9b3b19
