using PyCall
ct = pyimport("cantera")

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