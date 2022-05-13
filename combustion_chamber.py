import sys
import os

sys.path.insert(0,"/home/chinahg/GCresearch/cantera/build/python")

import cantera as ct
ct.add_directory('/user/chinahg')
ct.__file__

import numpy as np
import time
import math as math
#import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize

def combustion_chamber(P_CC0, V_CC, mdot_ox, mdot_f):

    #WSR

    #Create Gas Mixture
    gasExhaust = ct.Solution('gri30.yaml')
    gasExhaust.TPX = 1000, P_CC0, 'O2: 6, H2:1'
    gasExhaust.equilibrate('HP')
    
    gasFuel = ct.Hydrogen()
    gasFuel.TP = 20.15, 227527 #https://science.ksc.nasa.gov/shuttle/technology/sts-newsref/et.html#:~:text=LIQUID%20HYDROGEN%20TANK,-The%20liquid%20hydrogen&text=Its%20operating%20pressure%20range%20is,to%20the%20left%20aft%20umbilical.
    gasOx = ct.Oxygen()
    gasOx.TP = 90.15, 140000
    
    #Fuel and ox through low pressure turbopump and heat exchanger
    gasFuel.TP = 533.15, 1000000
    gasOx.TP = 145.9, 2900000 #weighted average of inlet ox temps

    #Create Reactor Infrastructure
    tankFuel = ct.Reservoir(gasFuel)
    tankOx = ct.Reservoir(gasOx)
    CC_exhaust1 = ct.Reservoir(gasExhaust)
    
    CC_reactor = ct.IdealGasReactor(gasExhaust) #create CC with random gas initially filling it
    
    #Connect tanks to CC
    flow_controller_fuel = ct.MassFlowController(upstream=tankFuel,downstream=CC_reactor,mdot=mdot_f)
    flow_controller_ox = ct.MassFlowController(upstream=tankOx,downstream=CC_reactor,mdot=mdot_ox)

    #Connect CC to exhaust reservoir
    press_controller_ox = ct.PressureController(upstream=CC_reactor, downstream=CC_exhaust1, master=flow_controller_ox, K=1e-5)
    press_controller_fuel = ct.PressureController(upstream=CC_reactor, downstream=CC_exhaust1, master=flow_controller_fuel, K=1e-5)
    #Add pressure controllers (fuel + ox)

    CC_reactorNet = ct.ReactorNet([CC_reactor])

    #React!
    t = 0
    del_t = 0.5
    state = ct.SolutionArray(gasExhaust, extra=['t'])

    while t < 1:
        t = t + del_t
        CC_reactorNet.advance(t)
        
    # Extract the final state of the reactor
    state.append(CC_reactor.thermo.state, t=t)
    
    print("Final CC Composition: ", gasExhaust.report())
    
    return(state)
