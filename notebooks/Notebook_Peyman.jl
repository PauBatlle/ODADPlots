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
using Plots, PlutoUI, Distributions, BoundingSphere, Contour

# ╔═╡ 65a16bbf-a6e9-4199-9b4b-7bf3c99b5f38
md""" 

Problem 1 setup: We want to estimate the mean $\mu$ of a gaussian distribution 		with known variance $\sigma^2$, from one observed data point $x_0$

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

# ╔═╡ 61958cf2-c345-4344-ab01-46d9b88d9017
begin
	get_β(alf) = 1-cdf(Chisq(1), -2*log(alf))
	get_α(bet) = exp(quantile(Chisq(1), 1-bet)/(-2))
	α = ifelse(bul, get_α(βp), αp)
	β = ifelse(!bul, get_β(αp), βp)
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

# ╔═╡ 0ef04a53-5d9f-4e39-9203-4a2b6947df70
md""" 

Problem 2 setup: We want to estimate the probabilities $\theta_1$ and $\theta_2$ of two different coins to land heads after observing $h_1$ heads and $t_1$ tails from coin 1 and $h_2$ heads and $t_2$ tails from coin 2
"""

# ╔═╡ 26bbc0d9-a602-4ceb-9e2f-ceefeab674dc
md"""
h₁: number of heads coin 1 $(@bind h₁ Slider(1:10, show_value=true, default=1))

t₁: number of tails coin 1 $(@bind t₁ Slider(1:10, show_value=true, default=1))

h₂: number of heads coin 2 $(@bind h₂ Slider(1:10, show_value=true, default=1))

t₂: number of tails coin 2 $(@bind t₂ Slider(1:10, show_value=true, default=1))
"""

# ╔═╡ 5c1071c4-12d5-426a-bd5e-4f339a2a6554
md"""
α $(@bind α₂ Slider(LinRange(0.01,0.99,101), show_value=true, default=0.2))
"""

# ╔═╡ 780edd65-9acf-44cc-a281-2c8a4b3ac973
begin
	function likelihoodcoins(θ₁,θ₂)
		num = θ₁^h₁ * (1-θ₁)^t₁ * θ₂^h₂ * (1-θ₂)^t₂
		den = (h₁/(h₁+t₁))^h₁ * (t₁/(h₁+t₁))^t₁ *  (h₂/(h₂+t₂))^h₂ * (t₂/(h₂+t₂))^t₂
		return num/den
	end
	x₂ = LinRange(0,1,100)
	z = [likelihoodcoins(xi,yi) for xi in x₂, yi in x₂];
	MLE_X = h₁/(h₁+t₁) 
	MLE_Y = h₂/(h₂+t₂)
	cc = contours(x₂,x₂,z,[α₂]);
	line=lines(levels(cc)[1])[1]
	xs = coordinates(line)[1]
	ys = coordinates(line)[2]
	points_boundary = [[xs[i], ys[i]] for i in 1:length(xs)]
	cent, rad = boundingsphere(points_boundary)
	rad_ap = round(rad,digits=3)
end;

# ╔═╡ 8c704d70-a318-4763-be57-43117f9b3b19
begin
	function circle(h, k, r)
	t = LinRange(0, 2*π, 500)
	h .+ r*sin.(t), k.+ r*cos.(t)
	end
	Plots.contour(x₂,x₂,likelihoodcoins, levels=[α₂], aspectratio= 1, xlims = (-0.2,1.2), ylims = (-0.2,1.2), color = "black", xlabel = "θ₁", ylabel = "θ₂", legend = :outerleft)
	plot!(circle(cent[1],cent[2],rad), label = "Bounding circle") 
	plot!([cent[1], cent[1]], [cent[2]-rad, cent[2]], label = "Radius $rad_ap")
	scatter!([cent[1]], [cent[2]], label = "Circumcenter")
	scatter!([MLE_X], [MLE_Y], label = "MLE")
end

# ╔═╡ Cell order:
# ╠═3d583f2a-9fce-4409-8f0d-83716b908bad
# ╟─edcec37a-a2ea-4104-baf4-c76a0cb26cd1
# ╟─65a16bbf-a6e9-4199-9b4b-7bf3c99b5f38
# ╟─b8d27d20-e226-4ed5-8721-d84ae962acc6
# ╟─61958cf2-c345-4344-ab01-46d9b88d9017
# ╟─0ef04a53-5d9f-4e39-9203-4a2b6947df70
# ╟─780edd65-9acf-44cc-a281-2c8a4b3ac973
# ╟─26bbc0d9-a602-4ceb-9e2f-ceefeab674dc
# ╟─5c1071c4-12d5-426a-bd5e-4f339a2a6554
# ╟─8c704d70-a318-4763-be57-43117f9b3b19
