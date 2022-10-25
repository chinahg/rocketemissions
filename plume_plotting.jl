using Interpolations, PyCall, OrdinaryDiffEq,
YAML, DelimitedFiles, CSV, HDF5, StructArrays, Random, 
NBInclude, PyPlot

function plotting(job_id::String)

    include("plumefunctions.jl")

    h_string = "16000"

    close("all")
    h5open("/home/chinahg/GCresearch/rocketemissions/slurm/plume_" *h_string* "m_" *job_id* ".h5")

    x = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/slurm/plume_" *h_string* "m_" *job_id* ".h5", "x")
    y = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/slurm/plume_" *h_string* "m_" *job_id* ".h5", "y")
    T = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/slurm/plume_" *h_string* "m_" *job_id* ".h5", "T")
    u = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/slurm/plume_" *h_string* "m_" *job_id* ".h5", "u")
    χ = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/slurm/plume_" *h_string* "m_" *job_id* ".h5", "X")
    s = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/slurm/plume_" *h_string* "m_" *job_id* ".h5", "s")
    n = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/slurm/plume_" *h_string* "m_" *job_id* ".h5", "n")
    Δϕ = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/slurm/plume_" *h_string* "m_" *job_id* ".h5", "delPhi")
    Δψ = HDF5.h5read("/home/chinahg/GCresearch/rocketemissions/slurm/plume_" *h_string* "m_" *job_id* ".h5", "delPsi")

    m = 1
    # Initialize arrays to save results
    T = T[m, :, :] #temperature
    u = u[m, :, :] #velocity
    χ = χ[m, :, :, :] #concentrations
    ϕ = cumsum(Δϕ)
    ψ = cumsum(Δψ)

    P_atm = 101325 #placeholder, need altitude dependent
    R = 287 #[J/kgK] placeholder

    #REGRID SOLUTION
    xx, yy, u_g, T_g, χ_gO2 =  regrid_solution(x[m,:], y[m,:,:], u, T, χ[:,:,4])
    xx, yy, u_g, T_g, χ_gN2 =  regrid_solution(x[m,:], y[m,:,:], u, T, χ[:,:,48])
    #xx, yy, u_g, T_g, χ_gNO2 = regrid_solution(x, y, u, T, χ[:, :, 37])
    #xx, yy, u_g, T_g, χ_gN2O = regrid_solution(x, y, u, T, χ[:, :, 38])
    #xx, yy, u_g, T_g, χ_gNO = regrid_solution(x, y, u, T, χ[:, :, 36])
    xx, yy, u_g, T_g, χ_gAr = regrid_solution(x[m,:], y[m,:,:], u, T, χ[:, :, 49])

    fig,axC = plt.subplots(2,1)
    axC[1,1].plot(yy[:,1], ψ)
    axC[1,1].set_xlabel("y")
    axC[1,1].set_ylabel("ψ")
    axC[2,1].plot(xx, ϕ)
    axC[2,1].set_xlabel("x")
    axC[2,1].set_ylabel("ϕ")
    fig.suptitle("Transformation Conversion")
    fig.tight_layout()
    savefig("/home/chinahg/GCresearch/rocketemissions/rockettests/" * h_string * "m/" * job_id * "_conversion.png")

    Xarea = zeros(s,n)
    i = 1
    j = 1
    rho_tot = zeros(s,n)
    MFar = zeros(s,n)
    MFn2 = zeros(s,n)
    MFo2 = zeros(s,n)
    mdot_Ar_sum = zeros(n)
    mdot_N2_sum = zeros(n)
    mdot_O2_sum = zeros(n)
    mdot_Ar = zeros(s,n)
    mdot_N2 = zeros(s,n)
    mdot_O2 = zeros(s,n)
    ppm_ar = zeros(n) 
    ppm_n2 = zeros(n)
    ppm_o2 = zeros(n)
    mdot_tot = zeros(n)
    ppm_tot = zeros(n)

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

    #IN XY
    r = length(yy)
    k = length(xx)

    mdot_Ar_sum_g = zeros(k)
    mdot_N2_sum_g = zeros(k)
    mdot_O2_sum_g = zeros(k)
    mdot_Ar_g = zeros(r,k)
    mdot_N2_g = zeros(r,k)
    mdot_O2_g = zeros(r,k)
    Xarea_g = zeros(r,k)
    rho_tot_g = zeros(r,k)
    MFar_g = zeros(r,k)
    MFn2_g = zeros(r,k)
    MFo2_g = zeros(r,k)
    mdot_tot_g = zeros(k)
    println("started xy plot")

    # IN XY COORDINATES
    for i = 1:k, j = 1:(r-1)
        rho_tot_g[j,i] = P_atm/(R*T_g[j,i]) #[kg/m^3] #need to convert to xx yy coordinates
        MFar_g[j,i] = (χ_gAr[j,i]*(1/28.97)*(39.9))/1000000 #[kg/kg tot]
        MFn2_g[j,i] = (χ_gN2[j,i]*(1/28.97)*(28))/1000000 #[kg/kg tot]
        MFo2_g[j,i] = (χ_gO2[j,i]*(1/28.97)*(32))/1000000 #[kg/kg tot]
        Xarea_g[j,i] = pi*(yy[j+1]^2 - yy[j]^2)
        mdot_Ar_g[j,i] = Xarea_g[j,i]*MFar_g[j,i]*rho_tot_g[j,i]*u_g[j,i]
        mdot_N2_g[j,i] = Xarea_g[j,i]*MFn2_g[j,i]*rho_tot_g[j,i]*u_g[j,i]
        mdot_O2_g[j,i] = Xarea_g[j,i]*MFo2_g[j,i]*rho_tot_g[j,i]*u_g[j,i]

        #TODO: populate the last entry of sum arrays (stops at s-1)
        mdot_O2_sum_g[i] += mdot_O2_g[j,i]
        mdot_N2_sum_g[i] += mdot_N2_g[j,i]
        mdot_Ar_sum_g[i] += mdot_Ar_g[j,i]
        mdot_tot_g[i] += mdot_Ar_g[j,i] + mdot_N2_g[j,i] + mdot_O2_g[j,i]
    end

    println("populating last gridded index")
    mdot_Ar_g[length(yy),:] = mdot_Ar_g[length(yy)-1,:]
    mdot_N2_g[length(yy),:] = mdot_N2_g[length(yy)-1,:]
    mdot_O2_g[length(yy),:] = mdot_O2_g[length(yy)-1,:]
    println("populated last gridded index")
    
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
    
    println("plotted ", h_string[m], " altitude!\n")

    println("done plotting!")
end