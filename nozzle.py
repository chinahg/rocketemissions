import sys
import os

sys.path.insert(0,"/home/chinahg/GCresearch/cantera/build/python")

import cantera as ct
ct.add_directory('/user/chinahg')
ct.__file__

import numpy as np
import time
import math as math
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize
import scipy.integrate as integrate

def Area_Ratio_Mach(x, A_throat, A_exit, gamma):
    #returns mach number given the turning angle
    return((-A_exit/A_throat) + ((0.5*(gamma+1))**(-(gamma+1)/(2*gamma-2))) * ((1+((gamma-1)*0.5*x**2))**((gamma+1)/(2*gamma-2)))/x)

class ReactorOde:
    def __init__(self, gas):
        # Parameters of the ODE system and auxiliary data are stored in the
        # ReactorOde object.
        self.gas = gas
        
    def __call__(self, t, y):
        """the ODE function, y' = f(t,y) """
        nsp = 53 #number of species in mechanism
        
        drdx = 0.0900547
        dAdx = 3.1415*2*drdx

        mdot = 67.35 + 404.79
        A_in = 0.0599
        
        # State vector is [T, Y_1, Y_2, ... Y_K]
        self.gas.TDY = y[1], y[0], y[2:nsp+2]
        
        rho = self.gas.density_mass
        T = self.gas.T
        Y = self.gas.Y
        
        #converging
        #create new function to find dAdx and A etc
        A = A_in+dAdx*t
        
        MW_mix = self.gas.mean_molecular_weight
        Ru = ct.gas_constant
        R = Ru/MW_mix
        nsp = 53 #nSpecies(gas)
        vx = mdot/(rho*A)
        P = rho*R*T
        
        MW = self.gas.molecular_weights
        h_k = self.gas.partial_molar_enthalpies/self.gas.molecular_weights #J/kg
        h = self.gas.enthalpy_mass #average enthalpy of mixture [J/kg]
        w = self.gas.net_production_rates # 1/s
        Cp = self.gas.cp_mass #J/kg

        #self.gas.set_multiplier(0)

        #--------------------------------------------------------------------------
        #---F(1), F(2) and F(3:end) are the differential equations modelling the---
        #---density, temperature and mass fractions variations along a plug flow---
        #-------------------------reactor------------------------------------------
        #--------------------------------------------------------------------------
        
        dDdx = ((1-R/Cp)*((rho*vx)**2)*(1/A)*(dAdx) + rho*R/(vx*Cp)*sum(MW*w*(h_k-(MW_mix*Cp*T/MW))))/ (P*(1+(vx**2)/(Cp*T)) - rho*vx**2)
       
        dTdx = (vx*vx/(rho*Cp))*dDdx + vx*vx*(1/A)*(dAdx)/Cp - (1/(vx*rho*Cp))*sum(h_k*w*MW)
    
        dYkdx = w[0:nsp]*MW[0:nsp]/(rho*vx) #check if kg/s
            
        return np.hstack((dDdx,dTdx,dYkdx))

def nozzle(T_Noz1, P_Noz1, comp_Noz1, A_throat, A_exit, L_Noz, mdot_ox, mdot_f):

    ### ISENTROPIC CONVERGING SECTION ###

    #Initial conditions from CC
    T = T_Noz1
    P = P_Noz1
    gamma = 1.1

    #Calculate incoming mach number
    M_CC = sp.optimize.newton(Area_Ratio_Mach, 0.7, args=(A_throat, A_exit, gamma))

    #Calculate stagnation properties
    T_t = T*(1+(gamma-1)*0.5*M_CC**2)
    P_t = P*(1+(gamma-1)*0.5*M_CC**2)**(gamma/(gamma-1))

    #Calculate properties at throat
    P_throat = P_t*(2*gamma-1)**(-gamma/(gamma-1))
    T_throat = T_t*(1/(2*gamma-1))

    #print("\nMACH NUMBER AT NOZZLE INLET: ", M_CC)
    #print("\nTHROAT TEMPERATURE: ", T_throat)
    #print("\nTHROAT PRESSURE: ", P_throat)

    #Call diverging section, return final gas state (end of nozzle)
    states_final = nozzle_div(T_throat, P_throat, comp_Noz1, A_throat, A_exit, L_Noz, mdot_ox,mdot_f)
    return(states_final) 

def nozzle_div(T_Noz1, P_Noz1, comp_Noz1, A_throat, A_exit, L_Noz, mdot_ox, mdot_f):

    ### SET INITIAL CONDITIONS ###
    # Temperature of gas, in K
    T0 = T_Noz1 #CC temp

    # Pressure of gas, in Pa
    P0 = P_Noz1

    # Import the gas phase, read out key species indices:
    gas = ct.Solution('gri30.yaml')

    nsp = 53 #nSpecies(gas)

    # Set the initial state and then equilibrate for a given enthalpy and pressure:
    gas.TPX = T0,P0,comp_Noz1
    gas.equilibrate('HP')

    # Initial condition (spark!)
    y0 = np.hstack((gas.density, gas.T, gas.Y))

    ## REACTOR PROPERTIES ###
    # Length of the reactor, in m
    L = L_Noz

    ### SOLVE REACTOR ###

    # Integrate the equations, keeping T(t) and Y(k,t)
    states = ct.SolutionArray(gas, 1, extra={'x': [0.0]})

    ode = ReactorOde(gas)
    solver = integrate.ode(ode)
    solver.set_integrator('vode', method='bdf', with_jacobian=True, atol=0.0000000000001)

    y0 = np.hstack((gas.density, gas.T, gas.Y))
    solver.set_initial_value(y0, solver.t)

    i = 0
    dx = 0.001

    ### SOLVE ODES FOR INDIVIDUAL REACTOR AND SAVE STATE ###
    while solver.t < L:
    
        solver.integrate(solver.t + dx)
        gas.TDY = solver.y[1], solver.y[0], solver.y[2:nsp+2]

        states.append(gas.state, x=solver.t)

        i = i+1
        #print(i, solver.t, gas.report())

    return(states)