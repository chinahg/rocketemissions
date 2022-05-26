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
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize

def combustion_chamber(P_CC0, V_CC, mdot_ox, mdot_f):

    #WSR

    #Create Gas Mixture
    gasExhaust = ct.Solution('gri30.yaml')
    gasExhaust.TPY = 300, P_CC0, 'O2: 6, H2:1'
    gasExhaust.equilibrate('HP')
    #print(gasExhaust.report())
    
    gasFuel = ct.Solution('gri30.yaml') #ct.Hydrogen()
    gasFuel.TPX = 20.15, 227527, 'H2:1' #https://science.ksc.nasa.gov/shuttle/technology/sts-newsref/et.html#:~:text=LIQUID%20HYDROGEN%20TANK,-The%20liquid%20hydrogen&text=Its%20operating%20pressure%20range%20is,to%20the%20left%20aft%20umbilical.
    gasOx = ct.Solution('gri30.yaml') #ct.Oxygen()
    gasOx.TPX = 90.15, 140000, 'O2:1'
    
    #Fuel and ox through low pressure turbopump and heat exchanger
    gasFuel.TP = 533.15, P_CC0 #1000000
    gasOx.TP = 145.9, 2900000 #weighted average of inlet ox temps
    
    #Create Reactor Infrastructure
    tankFuel = ct.Reservoir(gasFuel)
    tankOx = ct.Reservoir(gasOx)
    CC_exhaust1 = ct.Reservoir(gasExhaust)
    
    CC_reactor = ct.IdealGasReactor(gasExhaust) #create CC with random gas initially filling it
    
    #Connect tanks to CC
    flow_controller_fuel = ct.MassFlowController(upstream=tankFuel,downstream=CC_reactor,mdot=mdot_f)
    flow_controller_ox = ct.MassFlowController(upstream=tankOx,downstream=CC_reactor,mdot=mdot_ox)
    
    #flow_controller_fuel2 = ct.MassFlowController(upstream=CC_reactor,downstream=CC_exhaust1,mdot=mdot_f+mdot_ox)
    #flow_controller_ox2 = ct.MassFlowController(upstream=CC_reactor,downstream=CC_exhaust1,mdot=mdot_ox)

    #Connect CC to exhaust reservoir
    press_controller_ox = ct.PressureController(upstream=CC_reactor, downstream=CC_exhaust1, master=flow_controller_ox, K=0.001)
    press_controller_fuel = ct.PressureController(CC_reactor, CC_exhaust1, master=flow_controller_fuel, K=0.001)

    CC_reactorNet = ct.ReactorNet([CC_reactor])

    #React!
    t = 0
    state = ct.SolutionArray(gasExhaust, extra=['t','volume','enthalpy','eox','efuel'])

    while t < 0.2:
        #t = t + del_t
        CC_reactorNet.step()
        t = CC_reactorNet.time
        #CC_reactorNet.advance(t)
        
        # Extract the state of the reactor
        state.append(CC_reactor.thermo.state, volume=CC_reactor.volume, t=t, enthalpy=gasExhaust.enthalpy_mole, efuel=gasFuel.enthalpy_mole, eox=gasOx.enthalpy_mole)
    
    A = 0.159 #CC cross sectional area
    n = len(state.T)
    u = (mdot_f+mdot_ox)/(gasExhaust.density*A)
    print("\nFinal CC Composition: ", gasExhaust.report(),"\n", u)
    
    return(state)
