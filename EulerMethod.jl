### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ d82d7324-62e2-4b72-ac67-cabd6b9a5f49
using Plots, PlutoUI

# ╔═╡ ec8b53fd-ec1a-4493-b517-444e9e0f8a8e
md"""
# Introduction & Motivation
"""

# ╔═╡ 648015f7-ad40-40ef-8d7e-11dd1ab67379
md"""
Consider the following initial value problem (IVP): (perhaps from a perspective where $x(t)$ is a dynamical variable, and $t$ is the time elapsed.)

$\frac{dx}{dt} = f(x,t)\ \ ; \hspace{1cm} x(t_0) = x_0$

Our aim is to 'solve' this problem. Analytically, this means that we would like to find a closed form expression for the function $x(t) \ \ \forall t \in [0, \infty)$. However, we find that most differential equations may not have simple analytical expressions for their solutions. This motivates us to look for a numerical solution instead where we trade the knowledge of an analytical expression in exchange for the values of the function $x(t)$ at a finite number of time-points $t$.

To begin with, this involves breaking our domain into a discrete set of values, $T_{discrete} = \{t_1,t_2,...,t_N\} \subset [0, \infty)$. In general, the spacing between $t_i$ need not be constant, but for starters, let us say they are evenly spaced. This produces a set of discrete values constructed like so (an arithmetic progression): 

$\boxed{t_k = t_1 + (k-1)\Delta t}$

Hence, solving this problem is equivalent to obtaining the values $x(t) = \{x_1,x_2,...,x_N\} \equiv \{x(t_1), x(t_2),...,x(t_N)\}$
$\forall t_i \in T_{discrete}.$ 

Euler's method is our first (and simplest) attempt to numerically approximate the solution of this differential equation. After obtaining these values, we can apply specialized interpolation schemes to construct our version of the true solution.
"""

# ╔═╡ 011b25ed-7873-48f0-8dcd-86828e3de60f
md"""
## The Euler method
"""

# ╔═╡ e186e1f6-0311-4fe7-94d3-9b11f27bd7f5
md"""
From the first principle of derivatives, we have 

$$\frac{dx}{dt} \bigg |_{t=t_n} = \lim_{\Delta t \to 0} \frac{x(t_n + \Delta t) - x(t_n)}{\Delta t}$$

We can discretize this, keeping in mind that $\Delta t \to 0$; 

$$\frac{dx}{dt} \bigg |_{t=t_n} = \frac{x_{n+1} - x_n}{t_{n+1} - t_n} = f(x_n,\ t_n)$$

Denoting $t_{n+1} - t_n = h$, the equation can be rearranged to give the value of $x_{n+1}$ in terms of $x_n$ like so;

$x_{n+1} = x_n + h \cdot f(x_n,\ t_n) \ \ ; \hspace{0.5cm} \forall n \in {1, 2, ..., N}$

The above equation gives an iterative scheme, called the forward Euler method. 
"""

# ╔═╡ 16edb287-443a-48bd-b657-308825adfb06
md"""!!! note "A recurring theme"
	This idea of approximating derivatives by replacing the limiting value $dt$ with a finite value $\Delta t$ in addition to being the basis of other integrations schemes for olving **ODEs**, also forms the central idea of the [method of finite differences](https://en.wikipedia.org/wiki/Finite_difference) used to solve partial differential equations, like those governing diffusion, wave propagation, etc.
"""

# ╔═╡ 95927608-2f5d-4eb5-8adf-98202077531b
md"""
## A graphical perspective
"""

# ╔═╡ b9f4be19-1d9a-4d92-ba0f-7cbb3471c4ff
md"""
With a plot, this process can be visualised in a more intuitive manner. In the graph shown below, there is a smooth blue curve which is the solution to some IVP $(\dot{x} = f(x, \ t))$. The dots on the curve represent the points $t_i \in T_{discrete}$, at which we would like to numerically approximate the solution. Using the iteration scheme shown above, the approximate solution at the specified points are shown in red."""

# ╔═╡ f67b2ea1-28af-4dee-b75c-3755649967a4
begin
	h = 0.5
	
	t = 1:0.01:7
	ticks = 1.2:h:6.2
	
	y_approx, y_exact = exp.(t) ./ 100, exp.(ticks)./100
	
	t0, y0 = ticks[1], y_approx[1]
	
	lines = [[t0, y0]]
	
	for i in 1:length(ticks)
		w, e = lines[i]
		foo_t = w + h
		foo_y = e + h * e
		push!(lines, [foo_t, foo_y])
	end	
end

# ╔═╡ ffea36f2-2059-4573-8e02-e95a343a2697
begin
	plot(t, y_approx, xticks = ticks, ylim = (-2,5), xlim = (-0.5,8), label = false)
	scatter!(ticks, y_exact, label = "Analytical solution", color = :blue, markersize = 2)
	scatter!(getindex.(lines,1), getindex.(lines,2), label = "Numerical solution", markercolor = :red, legend = :top, markersize = 2.5, xlabel = "t", ylabel = "x(t)")
end

# ╔═╡ 566f148c-0ec8-4fa6-b1bc-788ce145d599
md"""
Notice that the numerical solution diverges from the anaytical solution for larger $t$. This is due to a poorly chosen spacing, $\Delta t$ between the discrete points, which results in an accumulation of error (ideally, we would require $\Delta t \to 0$). 
"""

# ╔═╡ 7c6b7839-a4bb-491c-8813-dd6049dbe796
md"""
## Numerical errors
"""

# ╔═╡ 5a3908a7-573a-4eac-ae44-63b6f5539456
md"""
While performing any numerical operation, it is necessary to obtain some information about the numerical errors that we introduce by using these approximate iterative schemes to gauge the utility of our approximate solution. For the euler method with fixed discretization, the error analysis involved is fairly straightforward.

We begin with the taylor expansion of our function $x(t)$ about $t = t_n$, truncated upto first order;

$$x(t_n + \Delta t) = x(t_n) + \left(\frac{dx}{dt}\bigg |_{t = t_n}\right)\Delta t + \mathcal{O}(\Delta t^2)$$

Locally, we see that the error drops off quadratically with $\Delta t$. However, we perform this iterative scheme $N$ times, where $N = \frac{t_f - t_0}{\Delta t}$. 

So, globally the error accumulated is $\mathcal{O}(\Delta t)$.
"""

# ╔═╡ ace793b9-2dfd-4540-9980-0ff1a3a451b8
md"""!!! tip "Taylor expansions"
	Similar methods will often be useful to derive and obtain the error bars for various higher-order integration schemes.

"""

# ╔═╡ e3a0f46f-fba7-496b-bf08-a847958427d6
md"""

---
## Implementation
Check out these links to find language-specific implementations of the euler method. 
- [Julia](EulerMethod_Julia.html) version
- [Python](EulerMethod_Python.html) version
"""

# ╔═╡ 32186aa5-2128-4713-ba0d-a4d6e238cfba
TableOfContents()

# ╔═╡ b032b264-0b4b-4ac3-9c53-631f5c0c02fd
html"<style>
main {
    max-width:70%;
    margin-right: 375px !important;
}

</style>"

# ╔═╡ Cell order:
# ╟─ec8b53fd-ec1a-4493-b517-444e9e0f8a8e
# ╟─648015f7-ad40-40ef-8d7e-11dd1ab67379
# ╟─011b25ed-7873-48f0-8dcd-86828e3de60f
# ╟─e186e1f6-0311-4fe7-94d3-9b11f27bd7f5
# ╟─16edb287-443a-48bd-b657-308825adfb06
# ╟─95927608-2f5d-4eb5-8adf-98202077531b
# ╟─b9f4be19-1d9a-4d92-ba0f-7cbb3471c4ff
# ╟─f67b2ea1-28af-4dee-b75c-3755649967a4
# ╟─ffea36f2-2059-4573-8e02-e95a343a2697
# ╟─566f148c-0ec8-4fa6-b1bc-788ce145d599
# ╟─7c6b7839-a4bb-491c-8813-dd6049dbe796
# ╟─5a3908a7-573a-4eac-ae44-63b6f5539456
# ╟─ace793b9-2dfd-4540-9980-0ff1a3a451b8
# ╟─e3a0f46f-fba7-496b-bf08-a847958427d6
# ╟─d82d7324-62e2-4b72-ac67-cabd6b9a5f49
# ╟─32186aa5-2128-4713-ba0d-a4d6e238cfba
# ╟─b032b264-0b4b-4ac3-9c53-631f5c0c02fd
