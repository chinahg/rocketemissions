# import Pkg

# Pkg.add("Interpolations")
# Pkg.add("PyPlot")
# Pkg.add("PyCall")
# Pkg.add("OrdinaryDiffEq")
# Pkg.add("YAML")
# Pkg.add("DelimitedFiles")
# Pkg.add("CSV")
# Pkg.add("HDF5")
# Pkg.add("StructArrays")
# Pkg.add("Random")

using Interpolations, PyCall, OrdinaryDiffEq,
YAML, DelimitedFiles, CSV, HDF5, StructArrays, Random, 
NBInclude, PyPlot

ct = pyimport("cantera")

ENV["PYTHON"] = "/opt/conda/envs/lae2020/bin/python"
# Pkg.build("PyCall")

include("plumefunctions.jl")
include("plume_plotting.jl")

struct gas_type
    gas
end

println("started")
n_species = 53
upper = 20000 #[m]
lower = 16000 #[m]
space = convert(Int, (upper - lower) / 2000 + 1) #250
h = [16000] #Int.(LinRange(16000, 20000, space))
g = 1

n = parse(Int, ARGS[1]) #n steps in x dir
s = parse(Int, ARGS[2]) #s steps in y direction
ψ_init= parse(Float64, ARGS[3])
ϕ_init= parse(Float64, ARGS[4])
ψ_mult= parse(Float64, ARGS[5])
ϕ_mult= parse(Float64, ARGS[6])
job_id = ARGS[7]
set_T = ARGS[8] #imported or ambient
set_u = ARGS[9] #imported or ambient

println("Reading arguments from bash: s = ", s, " n = ", n)
println("Reading arguments from bash: ψ_init = ", ψ_init, " ϕ_init = ", ϕ_init)
println("Reading arguments from bash: ψ_mult = ", ψ_mult, " ϕ_mult = ", ϕ_mult)

T_save = zeros(length(h), s, n)
u_save = zeros(length(h), s, n)
χ_save = zeros(length(h), s, n, n_species)
y_save = zeros(length(h),s,n)
x_save = zeros(length(h),n)

gas_g = StructArray{gas_type}(undef,s,n,length(h))
gas_g .= [gas_type(0)]
h_string = string(h[1])

println("started populating phi and psi")
Δϕ = ϕ_init * ones(n) #step size in phi
for i = 2:n
    Δϕ[i] = ϕ_mult * Δϕ[i-1] #enlarge with each step by mult
end

Δψ = ψ_init * ones(s) #s vertical grid points
for i = 2:s
    Δψ[i] = ψ_mult * Δψ[i-1] #enlarge with each step by mult
end

println("started altitude for loop")
for m = 1:lastindex(h)
    h_string = string(h[m])

    ### IMPORT SHOCK EXIT CONDITIONS ###

    p_all = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/plot_data.h5", h_string * "m/P")
    p = convert(AbstractFloat, p_all[2])

    u_a = 1.11849E-19 * big(h[m])^5 - 1.14814E-14 * big(h[m])^4 + 4.22542E-10 * big(h[m])^3 - 6.92322E-06 * big(h[m])^2 + 6.58761E-02 * big(h[m]) + 5.37920E+01

    #curve fit #a = ambient vel [m/s] (speed of rocket) 
    T_a = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/plot_data.h5", h_string * "m/T_a")
    
    if set_T == "ambient"
        T0 = T_a #for testing
    else
        T0 = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/plot_data.h5", h_string * "m/T") #initial plume temperature
        T0 = convert(AbstractFloat, T0[2])
    end

    if set_u == "ambient"
        u0 = u_a #for testing
    else
        u0 = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/plot_data.h5", h_string * "m/u")#initial plume velocity
        u0 = convert(AbstractFloat, u0[2])
    end
 
    ### CALCULATE VELOCITY AND TEMPERATURE FIELDS (NO CHEMISTRY) ###

    Pr = 1.0 #prandtl number
    Le = 1.0 #lewis number
    u_init = u0 .* ones(s)
    T_init = T0 .* ones(s)

    R = 287.0
    y_init = compute_y(u_init, T_init, Δψ, R, p)

    #Geometry of plume
    radius = 1.147
    u_init[y_init.>radius] .= convert(AbstractFloat, u_a)
    T_init[y_init.>radius] .= convert(AbstractFloat, T_a)
   
    ambient = AmbientConditions(u_a, T_a, 1.0, 1.0, R, p, T0, u0)

    #Solve for T and u at all steps: NO CHEMISTRY
    x, y, u, T, ϵ = solve_exhaust_flow(u_init, T_init, ambient, n, Δϕ, Δψ)

    T_save[m, :, :] = T
    u_save[m, :, :] = u
    x_save[m, :] = x
    y_save[m, :, :] = y

    ### CHEMISTRY STARTS HERE

    #SET INITIAL PLUME CONDITIONS
    #import concentration data from upstream combustion
    χ0_full = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/plot_data.h5", h_string * "m/X")
    χ0 = zeros(n_species)
        #for d = 1: length(χ0_full[:,2])-1
        #    χ0[d] = χ0_full[d,2]

    #gas = ct.Solution("gri30.yaml")
    #gas.TPX = T0, p, "Ar:1"
    #println(gas.density*u0*1.186)
    #end
    #same ambient and plume conditions
    #χ0[48] = 780790 #ppm N2
    #χ0[4] = 209445 #ppm O2
    χ0[49] = 100000 #ppm Ar #CHANGE BACK AFTER DEBUG TO 9339
    #χ0[16] = 9765 #ppm CO2 426

    #ppm must sum to a million

    #Redefine species initial conditions for multiple species
    χ_a = zeros(n_species)

    χ_a[48] = 780790 #ppm N2
    χ_a[4] = 209445 #ppm O2
    χ_a[49] = 0 #ppm Ar #CHANGE BACK AFTER DEBUG TO 9339
    χ_a[16] = 9765 #ppm CO2 426

    χ_init = zeros(s, n_species)

    ### COMPUTE MOLE FRACTIONS IN PLUME ###
    println("at radius allocation")
    i = y_init .>= radius
    for j = 1:lastindex(i) 
        if i[j] == 0 #if y position is less than the initial plume radius, assign plume conditions
            for k = 1:n_species
                χ_init[j, k] = χ0[k] #initial plume conditions
            end
        else #if y position is greater than the initial plume radius, assign ambient conditions
            χ_init[j, 48] = χ_a[48]
            χ_init[j, 4] = χ_a[4]
            χ_init[j, 49] = χ_a[49]
            χ_init[j, 16] = χ_a[16]
        end
    end

    # Save ambient conditions to a struct to call later easily
    ambient = AmbientConditionsχ(u_a, T_a, χ_a, 1.0, 1.0, R, p, T0, u0, χ0)

    # Create 3D array to save concentrations at each y,x location
    χ = zeros(length(u[:, 1]), length(u[1, :]), n_species) #[y,x,species]

    χ[:, 1, :] = χ_init # at x=1 assign concentrations for all species at all y

    # Initialize arrays for Strang-Splitting
    χ_h0 = zeros(size(χ_init))
    χ_1 = zeros(size(χ_init))
    
    # Create gas object to store Reactor output gas object state
    gas = ct.Solution("gri30.yaml")

    # Create a dummy reactor to establish a global variable
    dummy_reactor = ct.IdealGasReactor(gas)
    
    println("starting splitting")
    for i = 1:n-1 #x
        
        for j = 1:n_species #species
            #calculate f0 at half step 0.5*Δϕ (x)
            χ_h0[:, j] = solve_exhaust_flow_χ(u[:, i], T[:, i], ambient, n, 0.5 * Δϕ[i], Δψ, χ[:, i, j], i, j)
            # concentration of species j at x = i and y = all
        end

        save_tuple = solve_reaction(χ_h0, T[:, i], Δϕ[i], ϵ[i], u[:, i], gas, i, χ_1, s, n_species, @view gas_g[:,i,m])
        gas_g.gas[:,i+1,m] .= save_tuple[2]
        χ_1 = save_tuple[1]

        for j = 1:n_species-1 #species
            #calculate f0 at full step Δϕ
            χ[:, i+1, j] = solve_exhaust_flow_χ(u[:, i], T[:, i], ambient, n, 0.5 * Δϕ[i], Δψ, χ_1[:, j], i, j)
        end

        i = i + 1
    end
    χ_save[m, :, :, :] = χ

    println(m, " of ", length(h), " altitudes done!\n")
end
println("done computing! 💕")

fid = h5open("plume_" * h_string * "m_" *job_id* ".h5", "w")
create_group(fid, job_id)

fid["n"] = n
fid["s"] = s
fid["x"] = x_save
fid["y"] = y_save
fid["T"] = T_save
fid["X"] = χ_save
fid["u"] = u_save
fid["delPhi"] = Δϕ
fid["delPsi"] = Δψ
fid["job_id"] = job_id
fid["altitude"] = h_string
close(fid)

plot_flag = true

if plot_flag
    plotting(job_id)
end