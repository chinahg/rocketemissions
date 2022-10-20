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

println("Reading arguments from bash: s = ", s, " n = ", n)
println("Reading arguments from bash: ψ_init = ", ψ_init, " ϕ_init = ", ϕ_init)
println("Reading arguments from bash: ψ_mult = ", ψ_mult, " ϕ_mult = ", ϕ_mult)

T_save = zeros(length(h), s, n)
u_save = zeros(length(h), s, n)
χ_save = zeros(length(h), s, n, n_species)
x_save = zeros(n)
y_save = zeros(s)
gas_g = StructArray{gas_type}(undef,s,n,length(h))
gas_g .= [gas_type(0)]

println("started populating phi and psi")
Δϕ = ϕ_init * ones(n) #step size in phi
for i = 2:n
    Δϕ[i] = ϕ_mult * Δϕ[i-1] #enlarge with each step by 1.1
end

Δψ = ψ_init * ones(s) #s vertical grid points
for i = 2:s
    Δψ[i] = ψ_mult * Δψ[i-1] #enlarge with each step by 1.1
end

println("started altitude for loop")
for m = 1:lastindex(h)
    h_string = string(h[m])

    ### IMPORT SHOCK EXIT CONDITIONS ###

    u0 = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/plot_data.h5", h_string * "m/u")#initial plume velocity
    u0 = convert(AbstractFloat, u0[2])

    T0 = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/plot_data.h5", h_string * "m/T") #initial plume temperature
    T0 = convert(AbstractFloat, T0[2])

    p_all = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/plot_data.h5", h_string * "m/P")
    p = convert(AbstractFloat, p_all[2])

    u_a = 1.11849E-19 * big(h[m])^5 - 1.14814E-14 * big(h[m])^4 + 4.22542E-10 * big(h[m])^3 - 6.92322E-06 * big(h[m])^2 + 6.58761E-02 * big(h[m]) + 5.37920E+01

    #curve fit #a = ambient vel [m/s] (speed of rocket) 
    T_a = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/plot_data.h5", h_string * "m/T_a")
    #T0 = T_a #for testing see vincent messages
    #u0 = u_a #for testing see vincent messages
 
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
    x_save = x
    y_save = y
    T_save[m, :, :] = T
    u_save[m, :, :] = u

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
    #gas = ct.Solution("gri30.yaml")

    # Create a dummy reactor to ????
    #dummy_reactor = ct.IdealGasReactor(gas)
    
    println("starting splitting")
    Threads.@threads for i = 1:n-1 #x
        
        for j = 1:n_species #species
            #calculate f0 at half step 0.5*Δϕ (x)
            χ_h0[:, j] = solve_exhaust_flow_χ(u[:, i], T[:, i], ambient, n, 0.5 * Δϕ[i], Δψ, χ[:, i, j], i, j)
            # concentration of species j at x = i and y = all
        end

        #save_tuple = solve_reaction(χ_h0, T[:, i], Δϕ[i], ϵ[i], u[:, i], gas, i, χ_1, s, n_species, @view gas_g[:,i,m])
        #gas_g.gas[:,i+1,m] .= save_tuple[2]
        #χ_1 = save_tuple[1]

        for j = 1:n_species-1 #species
            #calculate f0 at full step Δϕ
            χ[:, i+1, j] = solve_exhaust_flow_χ(u[:, i], T[:, i], ambient, n, 0.5 * Δϕ[i], Δψ, χ_h0[:, j], i, j)
        end

        i = i + 1
    end
    χ_save[m, :, :, :] = χ

    println(m, " of ", length(h), " altitudes done!\n")
end
println("done computing! 💕")

plot_flag = true

if plot_flag
    close("all")
    # Define altitude range to plot over
    #h = [16000] #Int.(LinRange(16000, 20000, space)) #make an array and put in loop for all alts
    h_string_arr = string(h)

    for m = 1:lastindex(h) #For all specified altitudes

        # Initialize arrays to save results
        T = T_save[m, :, :] #temperature
        u = u_save[m, :, :] #velocity
        χ = χ_save[m, :, :, :] #concentrations
        x = x_save
        y = y_save
        h_string = string(h[m])
        P_atm = 101325 #placeholder, need altitude dependent
        R = 287 #[J/kgK] placeholder

        #REGRID SOLUTION
        #xx, yy, u_g, T_g, χ_gO2 =  regrid_solution(x, y, u, T, χ[:,:,4], 0.01)
        #xx, yy, u_g, T_g, χ_gN2 =  regrid_solution(x, y, u, T, χ[:,:,48], 0.01)
        #xx, yy, u_g, T_g, χ_gNO2 = regrid_solution(x, y, u, T, χ[:, :, 37], 0.01)
        #xx, yy, u_g, T_g, χ_gN2O = regrid_solution(x, y, u, T, χ[:, :, 38], 0.01)
        #xx, yy, u_g, T_g, χ_gNO = regrid_solution(x, y, u, T, χ[:, :, 36], 0.01)
        #xx, yy, u_g, T_g, χ_gAr = regrid_solution(x, y, u, T, χ[:, :, 49], 0.01)

        #PLOT 2D MAPS AND SAVE
        x_max = 20 #[m]

        #CALCULATE NO EI FOR ALTITUDE
        # Sum of all NOx species as a function of x
        sumNO = zeros(n)
        sumNO2 = zeros(n)
        sumN2O = zeros(n)
        Xarea = zeros(s,n)
        EI_Ar = zeros(s,n)
        i = 1
        j = 1

        #x_lim = length(xx)
        #y_lim = length(yy)

        MW_Ar = 39.95 #[kg/mol]
        mdot_fuel = 67.35 #[kg/s]
        rho_tot = zeros(s,n)
        MFar = zeros(s,n)
        MFn2 = zeros(s,n)
        MFo2 = zeros(s,n)

        mdot_Ar_sum = zeros(n)
        mdot_N2_sum = zeros(n)
        mdot_O2_sum = zeros(n)
        #mdot_Ar_sumTrunc = zeros(n)
        #mdot_N2_sumTrunc = zeros(n)
        #mdot_O2_sumTrunc = zeros(n)
        mdot_Ar = zeros(s,n)
        mdot_N2 = zeros(s,n)
        mdot_O2 = zeros(s,n)
        ppm_ar = zeros(n) 
        ppm_n2 = zeros(n)
        ppm_o2 = zeros(n)
        mdot_tot = zeros(n)
        ppm_tot = zeros(n)

        ϕ = cumsum(Δϕ)
        ψ = cumsum(Δψ)

        for i = 1:n, j = 1:s-1
            rho_tot[j,i] = P_atm/(R*T[j,i]) #[kg/m^3] #need to convert to xx yy coordinates
            MFar[j,i] = (χ[j,i,49]*(1/28.97)*(39.9))/1000000 #[kg/kg tot]
            MFn2[j,i] = (χ[j,i,48]*(1/28.97)*(28))/1000000 #[kg/kg tot]
            MFo2[j,i] = (χ[j,i,4]*(1/28.97)*(32))/1000000 #[kg/kg tot]
            Xarea[j,i] = pi*(ψ[j+1]^2 - ψ[j]^2)

            mdot_Ar[j,i] = Xarea[j,i]*MFar[j,i]*rho_tot[j,i]*u[j,i]
            mdot_N2[j,i] = Xarea[j,i]*MFn2[j,i]*rho_tot[j,i]*u[j,i]
            mdot_O2[j,i] = Xarea[j,i]*MFo2[j,i]*rho_tot[j,i]*u[j,i]

            ppm_ar[i] += χ[j,i,49]
            ppm_n2[i] += χ[j,i,48]
            ppm_o2[i] += χ[j,i,4]

            #TODO: populate the last entry of sum arrays (stops at s-1)
            mdot_O2_sum[i] += mdot_O2[j,i]
            mdot_N2_sum[i] += mdot_N2[j,i]
            mdot_Ar_sum[i] += mdot_Ar[j,i]
            
            #if j < 150
            #    mdot_O2_sumTrunc[i] += mdot_O2[i,j]
            #    mdot_N2_sumTrunc[i] += mdot_N2[i,j]
            #    mdot_Ar_sumTrunc[i] += mdot_Ar[i,j]
            #end
            mdot_tot[i] += mdot_Ar[j,i] + mdot_N2[j,i] + mdot_O2[j,i]
        end

        ppm_tot = ppm_ar + ppm_n2 + ppm_o2
        mdot_Ar[s,:] = mdot_Ar[s-1,:]
        mdot_N2[s,:] = mdot_N2[s-1,:]
        mdot_O2[s,:] = mdot_O2[s-1,:]

        #fig,axc = plt.subplots()
        #axc.plot(ϕ, u[s,:], label="5000/5000")
        #axc.plot(ϕ, u[s-1,:], label="4999/5000")
        #axc.plot(ϕ, u[s-50,:], label="4950/5000")
        #axc.plot(ϕ, u[s-100,:], label="4900/5000")
        #axc.plot(ϕ, u[s-150,:], label="4850/5000")
        #axc.plot(ϕ, u[s-200,:], label="4800/5000")
        #axc.plot(ϕ, u[s-250,:], label="4750/5000")
        #axc.plot(ϕ, u[s-1000,:], label="4000/5000")
    
        #axc.set_xlabel("ϕ")
        #axc.set_ylabel("m/s at ψ = 0.1")
        #axc.set_xscale("log")
        #legend = axc.legend(loc = "upper right")
        #fig.suptitle("Centerline Velocity")
        #savefig("/home/chinahg/GCresearch/rocketemissions/rockettests/" * h_string * "m/" * job_id * "_u_cent.png")

        fig,axR = plt.subplots(2,2)
        axR[1,1].plot(ϕ, mdot_Ar[1,:])
        axR[1,1].set_xlabel("ϕ")
        axR[1,1].set_ylabel("kg/s Ar at ψ = 0.1")
        axR[1,1].set_xscale("log")

        axR[2,1].plot(ϕ, χ[1,:,49])
        axR[2,1].set_xlabel("ϕ")
        axR[2,1].set_ylabel("ppm Ar at ψ = 0.1")
        axR[2,1].set_xscale("log")

        axR[2,2].plot(ϕ, T[1,:,1])
        axR[2,2].set_xlabel("ϕ")
        axR[2,2].set_ylabel("T at ψ = 0.1")
        axR[2,2].set_xscale("log")
        fig.suptitle("Single Ring, 0.1< ψ <= 0.105")
        fig.tight_layout()
        savefig("/home/chinahg/GCresearch/rocketemissions/rockettests/" * h_string * "m/" * job_id * "_ring1.png")

        fig,axT = plt.subplots()
        im = axT.imshow(T[:,:,m], cmap="summer")
        axT.set_ylabel("ψ")
        axT.set_xlabel("ϕ")
        axT.set_title("Temperature")
        cbar = fig.colorbar(im,ax=axT)
        savefig("/home/chinahg/GCresearch/rocketemissions/rockettests/" * h_string * "m/" * job_id * "_T.png")

        fig,axu = plt.subplots()
        im = axu.imshow(u[:,:,m], cmap="viridis")
        axu.set_ylabel("ψ")
        axu.set_xlabel("ϕ")
        cbar = fig.colorbar(im,ax=axu)
        axu.set_title("Velocity")
        savefig("/home/chinahg/GCresearch/rocketemissions/rockettests/" * h_string * "m/" * job_id * "_u.png")


        #fig,axs = plt.subplots(2,2)
        #axs[1,1].plot(ϕ, mdot_tot)
        #axs[1,1].set_xlabel("ϕ")
        #axs[1,1].set_ylabel("Total mass flow rate [kg/s]")
        #axs[1,1].set_title("Total mdot")
        #axs[1,1].set_xscale("log")
    #
        #axs[1,2].plot(ϕ, mdot_N2_sumTrunc)
        #axs[1,2].set_xlabel("ϕ")
        #axs[1,2].set_ylabel("N2 mass flow rate [kg/s]")
        #axs[1,2].set_title("N2")
        #axs[1,2].set_xscale("log")
    #
        #axs[2,1].plot(ϕ, mdot_O2_sumTrunc)
        #axs[2,1].set_xlabel("ϕ")
        #axs[2,1].set_ylabel("O2 mass flow rate [kg/s]")
        #axs[2,1].set_xscale("log")
    #
        #axs[2,2].plot(ϕ, mdot_Ar_sumTrunc)
        #axs[2,2].set_xlabel("ϕ")
        #axs[2,2].set_ylabel("Ar mass flow rate [kg/s]")
        #axs[2,2].set_title("Ar")
        #axs[2,2].set_xscale("log")
        #
        #fig.suptitle("Truncated Mass Flow Rates, s = 250/300")
        #fig.tight_layout()
        #savefig("/home/chinahg/GCresearch/rocketemissions/rockettests/" * h_string * "m/" * job_id * "_mdot_trunc.png")

        fig,axX = plt.subplots(2,2)
        axX[1,1].plot(ϕ, mdot_tot)
        axX[1,1].set_xlabel("ϕ")
        axX[1,1].set_ylabel("Total mass flow rate [kg/s]")
        axX[1,1].set_title("Total mdot")
        axX[1,1].set_xscale("log")

        axX[1,2].plot(ϕ, mdot_N2_sum)
        axX[1,2].set_xlabel("ϕ")
        axX[1,2].set_ylabel("N2 mass flow rate [kg/s]")
        axX[1,2].set_title("N2")
        axX[1,2].set_xscale("log")

        axX[2,1].plot(ϕ, mdot_O2_sum)
        axX[2,1].set_xlabel("ϕ")
        axX[2,1].set_ylabel("O2 mass flow rate [kg/s]")
        axX[2,1].set_xscale("log")

        axX[2,2].plot(ϕ, mdot_Ar_sum)
        axX[2,2].set_xlabel("ϕ")
        axX[2,2].set_ylabel("Ar mass flow rate [kg/s]")
        axX[2,2].set_title("Ar")
        axX[2,2].set_xscale("log")

        fig.suptitle("Total Mass Flow (all rings)")
        fig.tight_layout()
        savefig("/home/chinahg/GCresearch/rocketemissions/rockettests/" * h_string * "m/" * job_id * "_mdot_tot.png")


        fig,axs2 = plt.subplots(2,2)
        axs2[1,1].plot(ϕ, ppm_tot)
        axs2[1,1].set_xlabel("ϕ")
        axs2[1,1].set_ylabel("Total ppm")
        axs2[1,1].set_title("Total ppm")
        axs2[1,1].set_xscale("log")

        axs2[1,2].plot(ϕ, ppm_n2)
        axs2[1,2].set_xlabel("ϕ")
        axs2[1,2].set_ylabel("N2 ppm")
        axs2[1,2].set_title("N2")
        axs2[1,2].set_xscale("log")

        axs2[2,1].plot(ϕ, ppm_o2)
        axs2[2,1].set_xlabel("ϕ")
        axs2[2,1].set_ylabel("O2 ppm")
        axs2[2,1].set_xscale("log")

        axs2[2,2].plot(ϕ, ppm_ar)
        axs2[2,2].set_xlabel("ϕ")
        axs2[2,2].set_ylabel("Ar ppm")
        axs2[2,2].set_title("Ar")
        axs2[2,2].set_xscale("log")
        fig.suptitle("Species ppm (all rings)")
        fig.tight_layout()
        savefig("/home/chinahg/GCresearch/rocketemissions/rockettests/" * h_string * "m/" * job_id * "_ppm_tot.png")

        fig,axXg = plt.subplots(2,2)
        custom_xlim = (0,250)
        plt.setp(axXg, xlim=custom_xlim)
        axXg[1,1].plot(xx[1:250], mdot_tot_g[1:250])
        axXg[1,1].set_xlabel("x")
        axXg[1,1].set_ylabel("Total mass flow rate [kg/s]")
        axXg[1,1].set_title("Total mdot")


        axXg[1,2].plot(xx[1:250], mdot_N2_sum_g[1:250])
        axXg[1,2].set_xlabel("x")
        axXg[1,2].set_ylabel("N2 mass flow rate [kg/s]")
        axXg[1,2].set_title("N2")


        axXg[2,1].plot(xx[1:250], mdot_O2_sum_g[1:250])
        axXg[2,1].set_xlabel("x")
        axXg[2,1].set_ylabel("O2 mass flow rate [kg/s]")


        axXg[2,2].plot(xx[1:250], mdot_Ar_sum_g[1:250])
        axXg[2,2].set_xlabel("x")
        axXg[2,2].set_ylabel("Ar mass flow rate [kg/s]")
        axXg[2,2].set_title("Ar")

        fig.suptitle("Total Mass Flow (all rings)")
        fig.tight_layout()
        savefig("/home/chinahg/GCresearch/rocketemissions/rockettests/" * h_string * "m/" * job_id * "_mdot_tot_g.png")

        println("Δ mdot Ar total is ", mdot_Ar_sum[n]-mdot_Ar_sum[1])
        println("Δ mdot Ar/Δ mdot total is ", (mdot_Ar_sum[n]-mdot_Ar_sum[1])/(mdot_tot[n]-mdot_tot[1]))
        println(mdot_Ar_sum[n], mdot_tot[n])

        println("plotted ", h_string[m], " altitude!\n")
        #println(mdot_Ar_sum[1])
        #println(mdot_Ar[:,1])
        #println(yy[473])
    end

    println("done plotting!")
end