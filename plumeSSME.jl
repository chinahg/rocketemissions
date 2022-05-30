import Pkg

Pkg.add("Interpolations")
Pkg.add("Plots")
Pkg.add("PyCall")
Pkg.add("OrdinaryDiffEq")
Pkg.add("yaml")


using Interpolations, Plots, PyCall, OrdinaryDiffEq
ct = pyimport("cantera")

loader=Loader

### IMPORT SHOCK EXIT CONDITIONS ###
#Load YAML file to append new data
stream = open(str(h)+"_altitude.yaml", 'r')
dictionary = yaml.safe_load(stream)

h = dictionary["Altitude [m]"]
u0 = dictionary["Shocks Exit Velocity [m/s]"] #initial plume velocity
T0 = dictionary["Shocks Exit Temperature [K]"] #initial plume temperature
χ0 = 4.e4 #tracer species
p = dictionary["Pressure [Pa]"]
u_a = 1E-19*h^5 - 1E-14*h^4 + 4E-10*h^3 - 7E-06*h^2 + 0.0659*h + 53.792 #curve fit #a = ambient vel [m/s] (speed of rocket) UPDATE 
T_a = dictionary["Temperature [K]"] 
χ_a = 70. #passive tracer mixing ratio [ppm]

### CALCULATE VELOCITY AND TEMPERATURE FIELDS (VIN) ###

n = 100 #n steps in x dir
Δϕ = 0.001 * ones(n) #step size in phi
for i=2:n
    Δϕ[i] = 2*Δϕ[i-1] #double step size each step forward
end
Δψ = 0.5*ones(50) #50 vertical grid points
for i=2:50
    Δψ[i] = 1.1*Δψ[i-1] #enlarge with each step by 1.1
end

Pr = 1. #prandtl number
Le = 1.
u_init = u0 .* ones(50)
T_init = T0 .* ones(50)
χ_init = χ0 .* ones(50)

R = 287.
y_init = compute_y(u_init, T_init, Δψ, R, p)

#Geometry of plume
u_init[y_init .> 0.2225] .= u_a
T_init[y_init .> 0.445] .= T_a
χ_init[y_init .> 0.445] .= χ_a

ambient = AmbientConditions(u_a, T_a, χ_a, 1., 1., R, p, T0, u0, χ0);

#Solve for T and u at all steps 
x, y, u, T, χ_vin, ϵ = solve_exhaust_flow(u_init, T_init, ambient, n, Δϕ, Δψ, χ_init=χ_i)

### COMPUTE MOLE FRACTIONS IN PLUME ###
#O2 and CO
#ppm must sum to a million

#Redefine species initial conditions for multiple species
n_species = 53 #for gri30
χ_a = zeros(n_species)

#add iimport of species here !!!!!!!!!!!!!!!!!!
χ_a[48] = 780790 #ppm N2
χ_a[4] = 209445 #ppm O2
χ_a[49] = 9339 #ppm Ar
χ_a[16] = 426 #ppm CO2

χ0 = zeros(n_species)
#From nozzle exit (update with nozzle exit)
χ0[1] = 0.0035593*10^6 #initial ppm H2
χ0[2] = 0.10857*10^6 #initial ppm H
χ0[3] = 0.74269*10^6 #initial ppm O
χ0[4] = 0.097942*10^6 #initial ppm O2
χ0[5] = 0.045682*10^6 #initial ppm OH
χ0[6] = 0.001399*10^6 #initial ppm H2O
χ0[7] = 0.00015422*10^6 #initial ppm HO2
χ0[8] = 0.00000080817*10^6 #initial ppm H2O2

χ_init = zeros(50,n_species)

i=y_init .>= 0.445
for j=1:length(i)
    if i[j] == 0
        χ_init[j,1] = χ0[1] #ppm H2
        χ_init[j,6] = χ0[6] #ppm H2O
        χ_init[j,4] = χ0[4] #ppm O2
    else
        χ_init[j,48] = χ_a[48]
        χ_init[j,4] = χ_a[4] 
        χ_init[j,49] = χ_a[49]
        χ_init[j,16] = χ_a[16]
    end
end

ambient = AmbientConditionsχ(u_a, T_a, χ_a, 1., 1., R, p, T0, u0, χ0)

χ = zeros(length(χ_vin[:,1]), length(χ_vin[1,:]), n_species) #[y,x,species]

χ[:,1,:] =χ_init #first x value (all y and x)
χ_h0 = zeros(size(χ_init))
χ_1 = zeros(size(χ_init))
ω = zeros(size(χ_init))
gas = ct.Solution("gri30.yaml")

println("started solver :)")

for i=1:n-1 #x 
    for j=1:n_species #species
        #calculate f0 at half step 0.5*Δϕ (x)
        χ_h0[:,j] = solve_exhaust_flow_χ(u[:,i], T[:,i], ambient, n, 0.5*Δϕ[i], Δψ, χ[:,i,j], i, j)
        
        for n = 1:length(χ_h0[:,j])
            if χ_h0[n,j] < 0
                χ_h0[n,j] = 0
            end
        end
    end
    
    χ_1 = solve_reaction(χ_h0, T[:,i], Δϕ[i], ϵ[i], u[:,i], gas, i)

    for j=1:n_species-1 #species
        #calculate f0 at full step Δϕ
        χ[:,i+1,j] = solve_exhaust_flow_χ(u[:,i], T[:,i], ambient, n, 0.5*Δϕ[i], Δψ, χ_1[:,j], i, j)
        
        k = i+1
        for n = 1:length(χ[:,i+1,j])
            if χ[n,k,j] < 0
                χ[n,k,j] = 0
            end
        end
    end
    
    i = i+1  
end
println("Splitting complete")
