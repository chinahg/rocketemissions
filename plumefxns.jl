import Pkg

Pkg.add("Interpolations")
Pkg.add("Plots")
Pkg.add("PyCall")
Pkg.add("OrdinaryDiffEq")
Pkg.add("DataFrames")


using Interpolations, Plots, PyCall, OrdinaryDiffEq, DataFrames
ct = pyimport("cantera")

"Composite type to hold ambient conditions, parameters and initial conditions"
struct AmbientConditions
    "Ambient velocity [m/s]"
    u_a::AbstractFloat
    "Ambient temperature [K]"
    T_a::AbstractFloat
    "Ambient passive tracer mixing ratio"
    χ_a::AbstractFloat
    "Prandtl number"
    Pr::AbstractFloat
    "Lewis number"
    Le::AbstractFloat
    "Specific gas constant of air [J/kg]"
    R::AbstractFloat
    "Ambient pressure [Pa]"
    p::AbstractFloat
    "Initial temperature [K]"
    T0::AbstractFloat
    "Initial velocity"
    u0::AbstractFloat
    "Initial passive tracer mixing ratio"
    χ0::AbstractFloat
end

"""
    second_order_central(x, i, Δ)
Computes the second order central difference approximation to the second
derivative of the variable x w.r.t. its discretization dimension (e.g. x is
the state, and its entries x[i] are the values for different radial dimensions)
using the step size Δ, at location i.
# Inputs
- `x` State vector
- `i` Index to compute second order difference for 
- 'Δ' The spacing between the vector elements
"""
function second_order_central(x::AbstractVector, i::Integer, Δ::AbstractFloat)
    return ((x[i+1] - 2 * x[i] + x[i-1])/Δ^2)
end

"""
    construct_rhs(u, T, y, Δψ, Δϕ, ambient; χ=nothing)
Constructs the right-hand-side vector for the system `Ax=b` resulting
from applying the Crank-Nicholson scheme to the governing equations.
# Inputs 
- `u` axial velocities
- `T` temperature
- `y` y-locations of grid points 
- `Δψ` Spacing between grid points in ϕ-ψ space
- `Δϕ` Step size in ϕ-ψ space
- `ambient` Ambient conditions
- `χ` Mixing ratio of passive tracer
"""
function construct_rhs(u::AbstractVector, T::AbstractVector, y::AbstractVector,
     Δψ::AbstractVector, Δϕ::AbstractFloat, ambient::AmbientConditions;
    χ::AbstractVector=nothing)

    if χ !== nothing
        b = zeros(size(u)[1] + size(T)[1] + size(χ)[1])
        include_χ = true 
    else
        b = zeros(size(u)[1] + size(T)[1])
        include_χ = false 
    end

    ## Boundary conditions

    # Neumann conditions at y=0
    b[1] = 0.
    b[size(u)[1] + 1] = 0.0

    # Dirichlet conditions at y -> infinity 
    b[size(u)[1]] = ambient.u_a 
    b[size(u)[1] + size(T)[1]] = ambient.T_a

    # Definition for convenience 
    λ = ambient.u0 / ((ambient.R/ambient.p)^2 * ambient.T0^2) 
    Le = ambient.Le 

    # Loop through b 
    for i=2:size(u)[1]-1
        b[i] = u[i]/Δϕ +  0.5 * λ * y[i]^2 * second_order_central(u, i, Δψ[i])
        b[i+size(u)[1]] = ( T[i]/Δϕ + 0.5 * λ * y[i]^2 
                    * second_order_central(T, i, Δψ[i]) / ambient.Pr )
    end

    if include_χ
        for i=2:size(u)[1]-1
            b[i+2*size(u)[1]] = (χ[i]/Δϕ +  Le * 0.5 * λ * y[i]^2
             * second_order_central(χ, i, Δψ[i]) / ambient.Pr )
        end
        b[size(u)[1] + size(T)[1] + 1] = 0.0
        b[end] = ambient.χ_a
    end

    return b
end

"""
    construct_tridiagonal_matrix(n, Δψ, Δϕ, y, ambient; include_χ=false)
Constructs the matrix `A` for the system `Ax=b` resulting
from applying the Crank-Nicholson scheme to the governing equations.
# Inputs 
- `n` Number of grid points in radial direction
- `Δψ` Spacing between grid points in ϕ-ψ space
- `Δϕ` Step size in ϕ-ψ space
- `y` y-locations of grid points 
- `ambient` Ambient conditions
- `include_χ` indicates if passive tracer is simulated as well
"""
function construct_tridiagonal_matrix(n::Integer, Δψ::AbstractVector,
    Δϕ::AbstractFloat, y::AbstractVector, ambient::AmbientConditions;
    include_χ::Bool=false) 

    if include_χ
        A = zeros((3*n, 3*n))
    else
        A = zeros((2*n, 2*n))
    end
    ## Boundary conditions 

    # Neumann conditions at y=0
    A[1, 1] = (-1/(Δψ[1]))
    A[1, 2] =  1/(Δψ[1])
    A[n+1, n+1] = -1/Δψ[1]
    A[n+1, n+2] =  1/Δψ[1] 

    # Dirichlet conditions at y -> infinity 
    A[n,n] = 1.0 
    A[2*n,2*n] = 1.0

    # Definition for convenience 
    λ = ambient.u0 / ((ambient.R/ambient.p)^2 * ambient.T0^2) 
    
    Pr = ambient.Pr
    Le = ambient.Le 
    # Loop through matrix 
    for i=2:n-1 

        # For u 
        A[i,i] = 1/Δϕ  + λ * y[i]^2 /(Δψ[i]^2)
        A[i, i-1] = - λ * y[i]^2 / (2 * Δψ[i]^2)
        A[i, i+1] = - λ * y[i]^2 / (2 * Δψ[i]^2)

        # For T 
        j = i + n
        A[j, j] = 1/Δϕ  + λ * y[i]^2 / (Pr * Δψ[i]^2)
        A[j, j-1] = - λ * y[i]^2 / (2 * Pr * Δψ[i]^2)
        A[j, j+1] = - λ * y[i]^2 / (2 * Pr * Δψ[i]^2)
    end
    
    if include_χ
        A[2*n+1:end, 2*n+1:end] = Le .* A[n+1:2*n, n+1:2*n]
        A[2*n+1, 2*n+1] = -1/Δψ[1]
        A[2*n+1, 2*n+2] =  1/Δψ[1]
        A[3*n, 3*n] = 1.0
    end

    return A
end

"""
    compute_y(u, T, Δψ, R, p)
Computes y locations corresponding to particular states, based on the 
grid spacing in ϕ-ψ space. 
# Inputs 
- `u` axial velocities
- `T` temperature
- `Δψ` Spacing in ϕ-ψ space
- `R` Specific gas constant for air in J/kg
- `p` Ambient pressure in Pa
"""
function compute_y(u::AbstractVector, T::AbstractVector, Δψ::AbstractVector,
    R::AbstractFloat, p::AbstractFloat)
    return (2*cumsum((R/p)*T ./ u) .* Δψ).^0.5
end

"""
   compute_x(ϵ, Δϕ)
Computes x locations of the grid based on turbulent mixing coefficient ϵ
and step sizes Δϕ
# Inputs 
- `ϵ` Values of the turbulent mixing coefficient
- `Δϕ` Step sizes
"""
function compute_x(ϵ::AbstractVector, Δϕ::AbstractVector)
    x = zeros(size(ϵ))
    x[2:end] = cumsum(Δϕ[1:end-1] ./ ϵ[1:end-1])
    return x 
end

"""
    compute_ψ(u, T, y, p, R)
Computes the transformed coordinate ψ based on the state and y locations.
# Inputs
- `u` axial velocities
- `T` temperature
- `y` y-locations of grid points 
- `R` Specific gas constant for air in J/kg
- `p` Ambient pressure in Pa
"""
function compute_ψ(u::AbstractVector, T::AbstractVector, y::AbstractVector,
    p::AbstractFloat, R::AbstractFloat)
    return cumsum(p./(R*T) .* (y[2:end] .- y[1:end-1]) .* u)
end

"""
    get_ϵ(κ, u, y)
Computes the turbulent mixng coefficient based on current velocity profile.
# Inputs
- `κ` Empirical parameter
- `u` axial velocities
- `y` y-locations of grid points 
"""
function get_ϵ(κ::AbstractFloat, u::AbstractVector, y::AbstractVector)
    # First find the half-width 
    fractional_change = (u .- u[end])./(u[1] - u[end])
    idx = argmin(abs.(fractional_change .- 0.5))
    half_width = y[idx] 
    return κ * (u[1] - u[end]) * half_width
end

"""
    solve_exhaust_flow(u_init, T_init, ambient, n, Δϕ, Δψ; χ_init=nothing)
Solves for the velocity, temperature and optionally passive tracer profiles
within an aircraft engine exhaust jet.
# Inputs
- `u_init` Initial velocity profile
- `T_init` Initial temperature profile
- `ambient` AmbientConditions struct
- `n` Number of grid points in the radial direction
- `Δϕ` Grid spacing in ϕ-ψ space.
- `Δψ` Grid spacing in ϕ-ψ space.
- `χ_init` Initial passive tracer profile.
"""
function solve_exhaust_flow(u_init::AbstractVector, T_init::AbstractVector,
    ambient::AmbientConditions, n::Integer, Δϕ::AbstractVector, 
    Δψ::AbstractVector; χ_init::AbstractVector=nothing)
    
    u_mem = zeros((size(u_init)[1], n))
    T_mem = zeros((size(u_init)[1], n))
    y_mem = zeros((size(u_init)[1], n))
    u_mem[:,1] = u_init
    T_mem[:,1] = T_init
    ϵ = zeros(n)

    include_χ = false 
    if χ_init !== nothing
        include_χ = true
        χ_mem = zeros((size(χ_init)[1], n))
        χ_mem[:,1] = χ_init
    end

    @inbounds for i=1:n-1
        y_mem[:,i] = compute_y(@view(u_mem[:,i]), @view(T_mem[:,i]), Δψ, ambient.R, ambient.p)
        
        if include_χ
            A = construct_tridiagonal_matrix(size(u_init)[1], Δψ, Δϕ[i],
                    @view(y_mem[:,i]), ambient, include_χ=include_χ)
            b = construct_rhs(@view(u_mem[:,i]), @view(T_mem[:,i]),
                    @view(y_mem[:,i]), Δψ, Δϕ[i], ambient, χ=@view(χ_mem[:,i]))
            sol = A\b
            u_mem[:, i+1] = sol[1:size(u_init)[1]]
            T_mem[:, i+1] = sol[size(u_init)[1]+1:size(u_init)[1] + size(T_init)[1]]
            χ_mem[:, i+1] = sol[size(u_init)[1]+size(T_init)[1] + 1:end]
        else
            A = construct_tridiagonal_matrix(size(u_init)[1], Δψ, Δϕ[i], @view(y_mem[:,i]), ambient)
            b = construct_rhs(@view(u_mem[:,i]), @view(T_mem[:,i]), @view(y_mem[:,i]), Δψ, Δϕ[i], ambient)
            sol = A\b 
            u_mem[:, i+1] = sol[1:size(u_init)[1]]
            T_mem[:, i+1] = sol[size(u_init)[1]+1:end]
        end
        ϵ[i] = get_ϵ(0.02, @view(u_mem[:,i]), @view(y_mem[:,i]))
    end

    x = compute_x(ϵ, Δϕ)
    if include_χ
        return x, y_mem, u_mem, T_mem, χ_mem, ϵ
    else
        return x, y_mem, u_mem, T_mem 
    end
end

"""
    regrid_solution(x, y, u, T, y_spacing)
Maps back solution from a grid in ϕ-ψ space to a grid in x-y space.
# Inputs
- `x` x coordinates corresponding to the grid
- `y` y coordinates corresponding to the grid
- `u` Velocity values of solution
- `T` Temperature values of solution
- `y_spacing` Desired spacing in `y` direction for output grid
"""
function regrid_solution(x::AbstractVector, y::AbstractMatrix, u::AbstractMatrix,
    T::AbstractMatrix, χ::AbstractMatrix, y_spacing::AbstractFloat)

    yy = 0:y_spacing:maximum(y)

    u_gridded = zeros((size(yy)[1], size(x)[1]))
    T_gridded = zeros((size(yy)[1], size(x)[1]))
    χ_gridded = zeros((size(yy)[1]), size(x)[1])

    for j=1:size(x)[1]
        u_interp_extrap = LinearInterpolation(y[:,j], u[:,j], extrapolation_bc=Line())
        u_gridded[:,j] = u_interp_extrap(yy)
        T_interp_extrap = LinearInterpolation(y[:,j], T[:,j], extrapolation_bc=Line())
        T_gridded[:,j] = T_interp_extrap(yy)
        χ_interp_extrap = LinearInterpolation(y[:,j], χ[:,j], extrapolation_bc=Line())
        χ_gridded[:,j] = χ_interp_extrap(yy)
    end

    return x, yy, u_gridded, T_gridded, χ_gridded
end

function show_profiles(x, y, var, cols)
    labels = ["x = $xx m" for xx in x[cols]]
    plot(var[:,cols], y[:,cols], dpi=300, label=permutedims(labels))
end

"""
    plot_heatmap(x, y, var, xlabel, ylabel, clabel, colormap)
Plots gridded solution. 
# Inputs
- `x` x location of grid points
- `y` y location of grid points
- `var` Solution variable values at grid points
- `xlabel` Label for horizontal axis
- `ylabel` Label for vertical axis
- `clabel` Label for colormap
- `colormap` Colormap to use (e.g. :viridis)
"""
function plot_heatmap(x, y, var, xlabel, ylabel, clabel, colormap, h_string, species)
    # This is kind of hard-coded, to avoid any issues with the first entry of
    # `x` which typically is zero. 
    # If you use more horizontal grid points, the value `22` doesn't make much
    # sense either

    xlabel!(xlabel)
    ylabel!(ylabel)
    xticks!([-2, -1, 0, 1, 2, 3], ["0.01", "0.1", "1", "10", "100", "1000"])
    plot = heatmap(log10.(x[2:22]), 10. .* range(0, stop=1, length=size(y)[1]), clim = (0,Inf),
    var[:,2:22], colorbar_title=clabel, size = (700, 500), dpi=300, c=colormap)

    return(plot)
    
end


function plot_aircraft_geometry(x, y, var, x_engine, engine_spacing, xmax, ymax, amb_value)


    # First do only half of the image 
    xx = range(0, stop=xmax, length=100)
    yy = range(0, stop=ymax, length=100)
    engine_idx = argmin(abs.(yy .- engine_spacing))
    values = amb_value .* ones((size(yy)[1], size(xx)[1]))
    

    for (i, xi) in enumerate(xx)

        if xi > x_engine
            x_idx = argmin(abs.(x.-(xi -x_engine)))
            interp = Interpolations.LinearInterpolation(y[:,x_idx], var[:,x_idx],
             extrapolation_bc=Line())
            values[engine_idx:end, i] .= interp.(yy[engine_idx:end] .- yy[engine_idx])
            values[1:engine_idx, i] .= interp.(yy[engine_idx] .- yy[1:engine_idx])
        end
    end

    # Now mirror 
    values = vcat(reverse(values; dims=1), values)
    return xx, vcat(-1 .* reverse(yy; dims=1), yy), values
end


struct AmbientConditionsχ
    "Ambient velocity [m/s]"
    u_a::AbstractFloat
    "Ambient temperature [K]"
    T_a::AbstractFloat
    "Ambient passive tracer mixing ratio"
    χ_a::Array
    "Prandtl number"
    Pr::AbstractFloat
    "Lewis number"
    Le::AbstractFloat
    "Specific gas constant of air [J/kg]"
    R::AbstractFloat
    "Ambient pressure [Pa]"
    p::AbstractFloat
    "Initial temperature [K]"
    T0::AbstractFloat
    "Initial velocity"
    u0::AbstractFloat
    "Initial passive tracer mixing ratio"
    χ0::Array
end


function solve_reaction(χ_h0, T, Δϕ, ϵ, u, gas, j) #there is O2 and CO where it shouldnt be... 
    #temp drops but should have ambient amounts instead of initial amounts
    P = 101325 #Pa
    ω = zeros(size(χ_h0))
    reactor = ct.IdealGasConstPressureReactor(gas)  
    
    #FOR NO REACTIONS
    gas.set_multiplier(0)
    
    for l = 1:50
        for k = 1:53
            if χ_h0[l,k] < 0
                println("negative species input at y=", l, " X=", k)
            end
        end
    end
    
    for i=1:length(χ_h0[:,1]) #index through all "y"s
        
        try
            gas.TPX = T[i], P, χ_h0[i,:] #GO THROUGH TO SEE IF SETTING AND UPDATING GAS CORRECTLY
        catch
           println(i)
           print(gas.report())
        end

        reactor.syncState()
        reactorNet = ct.ReactorNet([reactor])
        
        t_final = Δϕ/(u[i]*abs(ϵ))
        t = 0
        
        reactorNet.advance(t_final, apply_limit=false)

        χ_1[i,:] = 10^6 .*reactor.thermo.X #mole fraction to ppm state.X[len,:] 10^9 #kmol/m^3s, assume 1 m^3, to ppm #rates for specific y (i) and all species UNITS
        
    end
    return χ_1 #for all y and all species [50,n_species]

end

function construct_tridiagonal_matrix_χ(n::Integer, Δψ::AbstractVector,
    Δϕ::AbstractFloat, y::AbstractVector, ambient::AmbientConditionsχ) 

    A = zeros((n, n))

    ## Boundary conditions 
    # Neumann conditions at y=0 #change
    A[1, 1] = -1/(Δψ[1])
    A[1, 2] =  1/(Δψ[1])

    # Dirichlet conditions at y -> infinity #change
    A[n,n] = 1.0 

    # Definition for convenience
    Pr = ambient.Pr
    Le = ambient.Le 
    λ = (Le/Pr)*((ambient.u0*ambient.p^2)/((ambient.R*ambient.T0)^2))
    
    # Loop through matrix 
    for j=2:n-1 
        # For χ
        A[j, j] = 1/Δϕ  + λ * y[j]^2 * 2 / (Δψ[j]^2)
        A[j, j-1] = - λ * y[j]^2 / (Δψ[j]^2)
        A[j, j+1] = - λ * y[j]^2 / (Δψ[j]^2)
    end
    return A
end

function construct_rhs_χ(u, T, y, Δψ, Δϕ, ambient::AmbientConditionsχ, χ, j)

    b = zeros(size(χ)[1])

    ## Boundary conditions
    # Neumann conditions at y=0 #change
    b[1] = 0.

    # Dirichlet conditions at y -> infinity #change
    b[size(χ)[1]] = ambient.χ_a[j]

    # Definition for convenience 
    Pr = ambient.Pr
    Le = ambient.Le 
    λ = (Le/Pr)*((ambient.u0*ambient.p^2)/((ambient.R*ambient.T0)^2))

    # Loop through b 
    #print(size(χ)[1])
    for i=2:size(χ)[1]-1
        b[i] = χ[i]/Δϕ + λ * y[i]^2 * second_order_central(χ, i, Δψ[i])
        #took out 0.5 factor on lambda (not 100% sure why was there) ^^^
    end

    return b
end

function solve_exhaust_flow_χ(u_mem, T_mem, ambient::AmbientConditionsχ, n::Integer, 
        Δϕ,  Δψ::AbstractVector, χ_init, i, j)
   
    y_mem = zeros(size(u_mem))
    
    χ_mem = χ_init
    
    y_mem = compute_y(u_mem, T_mem, Δψ, ambient.R, ambient.p)

    A = construct_tridiagonal_matrix_χ(size(u_mem)[1], Δψ, Δϕ, y_mem, ambient)
    #display(A)
    b = construct_rhs_χ(u_mem, T_mem, y_mem, Δψ, Δϕ, ambient, χ_mem, j)

    sol = A\b
    
    χ_mem = sol
    
    return χ_mem
end