import Pkg

Pkg.add("Interpolations")
Pkg.add("Plots")
Pkg.add("PyCall")
Pkg.add("OrdinaryDiffEq")
Pkg.add("YAML")
Pkg.add("DelimitedFiles")
Pkg.add("CSV")
using CSV


using Interpolations, Plots, PyCall, OrdinaryDiffEq, YAML, DelimitedFiles
ct = pyimport("cantera")

include("plumefxns.jl")

### IMPORT SPECIES NAMES FROM MECHANISM ###
#gri30_species = YAML.load_file("gri30_julia.yaml")
#gri30_species = collect(keys(gri30_species))
#n_species = length(gri30_species) #for gri30
#print(n_species)
#Assign each species to an integer corresponding to its mechanism index
#d = 1
#while d <= n_species
#    string2var(gri30_species[d],d)
#    d = d+1
#end
#print(O2)
n_species =53
### CHOOSE ALTITUDE [m] ###
h = 20000
h_string = string(h)

### IMPORT SHOCK EXIT CONDITIONS ###
#Load YAML file to save shock
dictionary = YAML.load_file("rockettests/"*h_string*"m/"*h_string*"_altitude.yaml")


u0 = dictionary[27]["Shocks Exit Velocity [m/s]"][1] #initial plume velocity
T0 = dictionary[25]["Shocks Exit Temperature [K]"][1] #initial plume temperature
χ0 = 4.e4 #tracer species
p = dictionary[1]["Pressure [Pa]"][1]
u_a = 1E-19*h^5 - 1E-14*h^4 + 4E-10*h^3 - 7E-06*h^2 + 0.0659*h + 53.792 #curve fit #a = ambient vel [m/s] (speed of rocket) UPDATE 
T_a = dictionary[1]["Temperature [K]"][1]
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
x, y, u, T, χ_vin, ϵ = solve_exhaust_flow(u_init, T_init, ambient, n, Δϕ, Δψ, χ_init=χ_init)

### COMPUTE MOLE FRACTIONS IN PLUME ###
#O2 and CO
#ppm must sum to a million

#Redefine species initial conditions for multiple species
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

display(u_g)
#REGRID SOLUTION
p = 1
xx, yy, u_g, T_g, χ_g =  regrid_solution(x, y, u, T, χ[:,:,p], 0.01)
display(u_g)
### PLOTTING RESULTS ###
##Reshape data into 2D arrays to be stored in CSV
big_length = length(xx)*length(yy)#*n_species

u_g = reshape(u_g,(length(xx)*length(yy),1))
T_g = reshape(T_g,(length(xx)*length(yy),1))
χ_g = reshape(χ_g,(big_length,1))
x_new = zeros(big_length,1)
y_new = zeros(big_length,1)
u_new = zeros(big_length,1)
T_new = zeros(big_length,1)

x_new[1:length(xx),1] = xx
y_new[1:length(yy),1] = yy
u_new[1:length(xx)*length(yy),1] = u_g
T_new[1:length(xx)*length(yy),1] = T_g

#Save data in dataframe
df = DataFrame(altitude = h, 
               x = vec(x_new),
               y = vec(y_new),
               u = vec(u_new),
               T = vec(T_new),
               χ = vec(χ_g)
               )

##Save results to CSV (plot in Jupyter)
CSV.write("test.csv",  df)