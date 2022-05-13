import sys
import os

sys.path.insert(0,"/home/chinahg/GCresearch/cantera/build/python")
print(sys.path)

import cantera as ct
ct.add_directory('/user/chinahg')
ct.__file__

import numpy as np
import time
import math as math
#import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize


class NozzleReactor(ct.DelegatedIdealGasReactor):

    def __init__(self, *args, **kwargs): #takes arguments the user specifies when creating Reactor
        super().__init__(*args, **kwargs)
        self.pressure = 20000000 #PLACEHOLDER for initial pressure
    
    def after_initialize(self,t0):
        self.n_vars = self.n_vars + 2
    
    def after_update_state(self, y):
    #save solution to state vector y
        self.pressure = 0

gas_Noz = ct.Solution('gri30.yaml')
gas_Noz.TPX = 2000, 101325, 'H2:1'
#Create Reactor
Noz_reactor = NozzleReactor(gas_Noz)
Noz_reactorNet = ct.ReactorNet([Noz_reactor])
#SKIP REACTOR IMPLEMENTATION, JUST ASSUME WORKS FOR NOW

def Isentropic_Mach(x, A, A_throat, gamma):
    #returns mach number given the area ratio and gamma
    return((1/(x)*((1+((gamma-1)/2)*x**2)/(1+((gamma-1)/2)))**((gamma+1)/(2*gamma-2)))*(A/A_throat))

# Isentropic Pressure Profile #
def Isentropic_Pressure(Noz_reactorNet, A, i):  
    
    #solve for mach number at new area ratio
    gamma = Noz_reactorNet.thermo.cp/Noz_reactorNet.cv
    M = sp.optimize.newton(Isentropic_Mach, 2, args=(A, A_throat, gamma))
    a = math.sqrt(gamma*8.314*Noz_reactorNet.T) #use current reactor temperature to calculate a
    u = M*a
    
    #solve for new pressure at new area ratio
    P_static = Noz_reactorNet.P #pressure after react
    P_dyn = 0.5*Noz_reactorNet.density*u**2 #use velocity from new mach CHECK
    P_tot = P_static + P_dyn
    
    #new pressure with updated area
    P = P_tot*(1 + (gamma-1)/2 *M**2)^(-gamma/(gamma-1))
    
    dpdz = (P-Noz_reactorNet.P)/(delta_z)
    return dpdz

def nozzle_geometry(z, A_CC, A_throat, A_exit, CC_to_throat, throat_to_exit):
    #cone nozzle geometry
    if 0 <= z <= CC_to_throat:
        nozzle_radius = math.sqrt(A_CC/math.pi) - ((math.sqrt(A_CC/math.pi)-math.sqrt(A_throat/math.pi))/CC_to_throat)*z
        nozzle_area = math.pi*nozzle_radius**2
        return(nozzle_area)
                                                   
    elif CC_to_throat < z <= CC_to_throat+throat_to_exit:
        print('past throat!')
        nozzle_radius = math.sqrt(A_throat/math.pi) + ((math.sqrt(A_exit/math.pi)-math.sqrt(A_throat/math.pi))/throat_to_exit)*z
        nozzle_area = math.pi*nozzle_radius**2
        return(nozzle_area)
    
# Advance Reactor #
def nozzle_react(T_Noz1, P_Noz1, comp_Noz1, mdot_Noz, A_CC, A_throat, A_exit, CC_to_throat, throat_to_exit, u_CC_exit, mdot_ox, mdot_f):
    gas_Noz = ct.Solution('gri30.yaml') #MOVE OUTSIDE OF FXN?
    gas_Noz.TPX = T_Noz1, P_Noz1, comp_Noz1

    #Create Reactor
    Noz_reactor = NozzleReactor(gas_Noz)
    Noz_reactorNet = ct.ReactorNet([Noz_reactor])
    
    t = 0
    dt = 0.01
    n1 = 0
    z = 0
    u_Noz = u_CC_exit

    while A_Noz > A_throat:
        #need area at n1 for advance for custom function
        Noz_reactorNet.advance(t+dt) #advance to the next time until area function returns throat area
        
        z = z + u_Noz * dt #Move forward to new position z, equivalent to movement during delta t

        #Save the state
        Noz_states.append(Noz_reactor.thermo.state, z = z, u_Noz = u_Noz)

        #Calculate new pressure due to change in A
        #calculate new area based on delta t and velocity
        A_Noz = nozzle_geometry(z, A_CC, A_throat, A_exit, CC_to_throat, throat_to_exit) #returns new cross sectional area
        #calculate new dpdz
        Noz_reactorNet.dpdz = Isentropic_Pressure(Noz_reactorNet, A_Noz, i)
        #recalculate speed
        u_Noz = (mdot_ox + mdot_f)/(Noz_reactorNet.thermo.density*A_Noz)

        n1 = n1+1

    #save final state
    Noz_states = ct.SolutionArray(gas_Noz, extra=['z'])
    
    return(Noz_states)