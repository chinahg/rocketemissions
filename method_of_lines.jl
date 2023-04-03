# # import Pkg
# # Pkg.add("MethodOfLines")
# # Pkg.add("ModelingToolkit")
# # Pkg.add("DomainSets")
# # Pkg.add("OrdinaryDiffEq")
# # Pkg.add("Plots")

# using Pkg, MethodOfLines, ModelingToolkit, DomainSets, OrdinaryDiffEq, Plots
# Pkg.status()
# @parameters t, x
# @variables u(..)
# Dt = Differential(t)
# Dx = Differential(x)
# Dxx = Dx^2

# α = 1.1

# eq = Dt(u(t,x)) ~ α*Dxx(u(t,x))

# domain = [x ∈ Interval(0.0,10.0), t ∈ Interval(0.0,1.0)]

# ic_bc = [u(0.0,x) ~ exp(-(x-4.0)^2) + exp(-(x-6.0)^2),
#         u(t,0.0) ~ 0.0,
#         u(t,10.0) ~ 0.0]

# @named sys = PDESystem(eq,ic_bc,domain,[t,x],[u(t,x)])

# dx = 0.1
# dy = 1.0

# discretization = MOLFiniteDifference([x => dx], t, approx_order=2)

# prob = discretize(sys, discretization)
# sol = solve(prob, Tsit5(), saveat = 0.05)
# grid = get_discrete(sys,discretization)

# anim = @animate for (i, t_disc) in enumerate(sol[t])
    
#     plot(grid[x], map(d -> sol[d][i], grid[u(t, x)]), ylim = (0.,1.), label = "u", title = "t = $t_disc")

# end

# gif(anim, "heatrod.gif", fps = 10)

using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomainSets

# Parameters, variables, and derivatives
@parameters t x
@variables u(..)
Dt = Differential(t)
Dxx = Differential(x)^2

# 1D PDE and boundary conditions
eq  = Dt(u(t, x)) ~ Dxx(u(t, x))
bcs = [u(0, x) ~ cos(x),
        u(t, 0) ~ exp(-t),
        u(t, 1) ~ exp(-t) * cos(1)]

# Space and time domains
domains = [t ∈ Interval(0.0, 1.0),
           x ∈ Interval(0.0, 1.0)]

# PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, x], [u(t, x)])

# Method of lines discretization
dx = 0.1
order = 2
discretization = MOLFiniteDifference([x => dx], t)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization)

# Solve ODE problem
using OrdinaryDiffEq
sol = solve(prob, Tsit5(), saveat=0.2)

# Plot results and compare with exact solution
discrete_x = sol[x]
discrete_t = sol[t]
solu = sol[u(t, x)]

using Plots
plt = plot()

anim = @animate for i in 1:length(discrete_t)
    plot!(discrete_x, solu[i, :], label="Numerical, t=$(discrete_t[i])")
end
gif(anim,"test_anim.gif", fps=10)