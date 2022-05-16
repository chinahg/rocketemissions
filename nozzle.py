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

class ReactorOde:
    def __init__(self, gas):
        # Parameters of the ODE system and auxiliary data are stored in the
        # ReactorOde object.
        self.gas = gas
        

    def __call__(self, t, y):
        """the ODE function, y' = f(t,y) """
        nsp = 53 #number of species in mechanism
        
        dAdx = 1.325
        mdot = 67.35 + 404.79
        A_in = 0.0599
        
        # State vector is [T, Y_1, Y_2, ... Y_K]
        self.gas.TDY = y[1], y[0], y[2:nsp+2]
        
        rho = self.gas.density
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
        w = self.gas.net_production_rates
        Cp = self.gas.cp_mass

        #--------------------------------------------------------------------------
        #---F(1), F(2) and F(3:end) are the differential equations modelling the---
        #---density, temperature and mass fractions variations along a plug flow---
        #-------------------------reactor------------------------------------------
        #--------------------------------------------------------------------------
        
        dDdx = ((1-R/Cp)*((rho*vx)**2)*(1/A)*(dAdx) + 0*rho*R*sum(MW*w*(h_k-MW_mix*Cp*T/MW))/(vx*Cp) )/ (P*(1+vx**2/(Cp*T)) - rho*vx**2)
       
        dTdx = (vx*vx/(rho*Cp))*dDdx + vx*vx*(1/A)*(dAdx)/Cp - (1/(vx*rho*Cp))*sum(h_k*w)
    
        dYkdx = w[0:nsp]*MW[0:nsp]/(rho*vx)
        
        i = 0
        while i < nsp:
            if dYkdx[i] < 0:
                dYkdx[i] = 0
            i = i+1
            
        return np.hstack((dDdx,dTdx,dYkdx))

def nozzle_react(T_Noz1, P_Noz1, comp_Noz1, A_throat, A_exit, L_Noz, mdot_ox, mdot_f):

    ### SET INITIAL CONDITIONS ###
    # Temperature of gas, in K
    T0 = T_Noz1 #CC temp

    # Pressure of gas, in Pa
    P0 = P_Noz1

    # Import the gas phase, read out key species indices:
    gas = ct.Solution('gri30.yaml')
    ih2 = 0 #gas.speciesIndex('CH4')
    io2  = 4 #gas.speciesIndex(gas,'O2')

    nsp = 53 #nSpecies(gas)
    x = np.zeros(nsp) #array of length number of species

    # Set initial composition
    x[ih2] = 1
    x[io2] = 6

    # Set the initial state and then equilibrate for a given enthalpy and pressure:
    gas.TPY = T0,P0,comp_Noz1
    gas.equilibrate('HP')

    # Initial condition (spark!)
    y0 = np.hstack((gas.density, gas.T, gas.Y))

    ## REACTOR PROPERTIES ###
    # The Dimensions and conditions of the reactor are given below

    # Inlet Area, in m^2
    A_in = A_throat
    # Exit Area, in m^2
    A_out = A_exit
    # Length of the reactor, in m
    L = L_Noz
    # The whole reactor is divided into n small reactors
    n = 20
    # Mass flow rate into the reactor, in kg/s
    mdot_calc = mdot_ox + mdot_f #ox+fuel kg/s

    nsp = 53 #gas.nSpecies()

    ### SOLVE REACTOR ###

    # Integrate the equations, keeping T(t) and Y(k,t)
    states = ct.SolutionArray(gas, 1, extra={'x': [0.0]})

    ode = ReactorOde(gas)
    solver = integrate.ode(ode)
    solver.set_integrator('vode', method='bdf', with_jacobian=True, atol=0.0001)

    y0 = np.hstack((gas.density, gas.T, gas.Y))
    solver.set_initial_value(y0, solver.t)

    i = 0
    dx = 0.01

    ### SOLVE ODES FOR INDIVIDUAL REACTOR AND SAVE STATE ###
    while solver.t < L:
    
        solver.integrate(solver.t + dx)
        gas.TDY = solver.y[1], solver.y[0], solver.y[2:nsp+2]
        states.append(gas.state, x=solver.t)

        i = i+1
        #print(i, solver.t, states.T[i])

    return(states)

""" def nozzle_geometry(z, A_CC, A_throat, A_exit, CC_to_throat, throat_to_exit):
    #cone nozzle geometry
    if 0 <= z <= CC_to_throat:
        nozzle_radius = math.sqrt(A_CC/math.pi) - ((math.sqrt(A_CC/math.pi)-math.sqrt(A_throat/math.pi))/CC_to_throat)*z
        nozzle_area = math.pi*nozzle_radius**2
        return(nozzle_area) """