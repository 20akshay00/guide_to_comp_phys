### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 7141c6f9-25a9-4253-bae6-d194e48e1e7d
using PlutoUI, Plots

# ╔═╡ b19ec43d-b84e-4f6d-8aa7-5595dcb835a7
md"""
### Implementing the forward Euler method

---

"""

# ╔═╡ 9d0a01f0-e9f9-491e-8bfe-a34c6ac91113
md"""
Let us begin with a first-order differential equation involving a function of only one variable and then proceed to see how we can extend the same code using broadcasting in Julia, to numerically integrate a system of first-order differential equations.
"""

# ╔═╡ 8d43f02c-93af-4a01-b0c5-dcee10335071
md"""
Below is the code to solve an equation of the type,

$\dot x = f(x, t)$
"""

# ╔═╡ 2fb572d1-2fd4-4dc1-a100-d7d1d0dccd45
function forwardEuler(f, xi, span, dt = 1.0e-3)
    ti, tf = span 
    steps = floor((tf-ti)/dt) |> Int #obtains the number of points in the discretized domain 
	
	# Allocate empty arrays to store x and t values
    x = Vector{Float64}(undef, steps) 
	t = Vector{Float64}(undef, steps)
	
	# Set the initial values of the arrays
    x[1], t[1] = xi, 0.0
	
	# Euler method iteration
    for i in 1:(steps-1)
        x[i+1] = x[i] + dt * f(x[i],t[i])
		t[i+1] = t[i] + dt
    end
	
    return x, t
end

# ╔═╡ 63e917df-6ff8-412d-a60a-16e774504e18
md"""!!! info "Input parameters"
	`f` : which is the function giving the value of the derivative at x
	
	`x0` : this is the initial condition for the problem
	
	`span` : this is the span of the domain $[t_i, t_f]$ over which we wish to integrate the equation.
	
	`dt` : this is the step size, more often than not, $10^{-3}$ works well enough. But lower the step size, higher the accuracy!
"""

# ╔═╡ c01147de-1f99-4674-b856-f8b8e1865256
md"""
### Examples

---

"""

# ╔═╡ 78280809-36f2-41ca-ad5c-8c360396471c
md"""
##### 1) Exponential growth

We would like to study the solutions of this differential equation:

$\frac{dx}{dt} = kx$


"""

# ╔═╡ 05fde7ba-43c1-4e76-9412-e0f7046a07d5
md"""
Let us solve this for some random initial conditions:

1) Let `k = 10`,  `x(0.1) = 1` and we would like to solve for the solution in `t in [0.1, 1.0)`.
"""

# ╔═╡ 8518b0a6-7541-49cf-963e-73812f42d24e
let
	f(x, t) = 10*x
	x0 = 1
	span = (0.1, 1.0)
	solution, t = forwardEuler(f, x0, span, 1e-2)
	
	plot(t, solution, label = "Interpolated solution", legend = :topleft, lw = 2, color = :blue)
	scatter!(t, solution, label = "Euler method points", legend = :topleft, lw = 2, color = :red, markersize = 3)
end

# ╔═╡ 7ebef684-108b-4faf-9b39-19c5576e528e
md"""
As expected, we obtain an exponential growth. Try plotting the known analytical solution and plot the error in our numerical estimate.

---

"""

# ╔═╡ 321fe129-19e2-406d-b4ee-915d94b73fee
md"""
2) Let `k = -10`,  `x(0.1) = 1` and we would like to solve for the solution in `t in [0.1, 1.0)`.
"""

# ╔═╡ abe33b06-9378-40b7-b3f0-f630034b61f5
let
	f(x, t) = -10*x
	x0 = 1
	span = (0.1, 1.0)
	solution, t = forwardEuler(f, x0, span, 1e-2)
	
	plot(t, solution, label = "Interpolated solution", legend = :topleft, lw = 2, color = :blue)
	scatter!(t, solution, label = "Euler method points", legend = :topleft, lw = 2, color = :red, markersize = 3)
end

# ╔═╡ a498094e-16c9-4dc2-97af-5fee1bc44a18
md"""
By changing the sign of the exponent, we now obtain an exponentially decaying solution. This serves as a decent sanity-check for whether our integration scheme is working or not.
"""

# ╔═╡ f7fb6f63-5a16-43e0-b155-921262b23406
md"""

---

"""

# ╔═╡ 56cc30ca-c1cc-45cd-a5cd-2b3103e1ea02
md"""
##### 2. Logistic growth
This is a simple toy model that can be used to simulate population growth when there is a factor that limits their growth (a population cap). The model is given by the differential equation,

$\frac{dx}{dt} = \gamma \cdot x \cdot (1 - \frac{x}{K})$

-$x$ : number of individuals

-$\gamma$ : growth rate

-$K$ : population cap
"""

# ╔═╡ 121a761e-c5bc-4e6f-82b1-8a8499f57ee0
md"""
1) Let the model parameters be `gamma = 0.95`,  `K = 100`. Let `x(0.0) = 1` and we would like the solution for `t in [0.0, 11.0)`.

"""

# ╔═╡ 8da33672-5d94-41d5-a370-3414d8cb9ed7
let
	logistic(x,t) = 0.95*x*(1 - (x/100))
	x0 = 1
	span = (0, 11)
	solution, t = forwardEuler(logistic, x0, span, 1e-1)
	plot(t, solution, label = "Interpolated solution", legend = :top, lw = 2, color = :blue)
	scatter!(t, solution, label = "Euler method points", legend = :top, lw = 2, color = :red, markersize = 3)
end

# ╔═╡ bd331fdf-9606-445b-9533-83a48b81cbb3
md"""
We can now play around with the parameter values and initial conditions to study the behaviour of this differential equation.
"""

# ╔═╡ 74ab82e0-bdea-11eb-2998-793ccc59b7a5
html"<style>
main {
    max-width:75%;
    padding: unset;
    margin-right: unset !important;
    margin-bottom: 20px;
    align-self: center !important;
}
pluto-helpbox {
    visibility: hidden;
}
</style>"

# ╔═╡ Cell order:
# ╠═7141c6f9-25a9-4253-bae6-d194e48e1e7d
# ╟─b19ec43d-b84e-4f6d-8aa7-5595dcb835a7
# ╟─9d0a01f0-e9f9-491e-8bfe-a34c6ac91113
# ╟─8d43f02c-93af-4a01-b0c5-dcee10335071
# ╠═2fb572d1-2fd4-4dc1-a100-d7d1d0dccd45
# ╟─63e917df-6ff8-412d-a60a-16e774504e18
# ╟─c01147de-1f99-4674-b856-f8b8e1865256
# ╟─78280809-36f2-41ca-ad5c-8c360396471c
# ╟─05fde7ba-43c1-4e76-9412-e0f7046a07d5
# ╠═8518b0a6-7541-49cf-963e-73812f42d24e
# ╟─7ebef684-108b-4faf-9b39-19c5576e528e
# ╟─321fe129-19e2-406d-b4ee-915d94b73fee
# ╠═abe33b06-9378-40b7-b3f0-f630034b61f5
# ╟─a498094e-16c9-4dc2-97af-5fee1bc44a18
# ╟─f7fb6f63-5a16-43e0-b155-921262b23406
# ╟─56cc30ca-c1cc-45cd-a5cd-2b3103e1ea02
# ╟─121a761e-c5bc-4e6f-82b1-8a8499f57ee0
# ╠═8da33672-5d94-41d5-a370-3414d8cb9ed7
# ╟─bd331fdf-9606-445b-9533-83a48b81cbb3
# ╟─74ab82e0-bdea-11eb-2998-793ccc59b7a5
