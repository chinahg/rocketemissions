#!/opt/conda/envs/lae2020/bin/python3
from importlib.abc import Loader
import sys
import os
import ruamel.yaml

sys.path.insert(0,"/home/chinahg/GCresearch/cantera/build/python")
sys.path.insert(1,"/home/chinahg/GCresearch/rocketemissions")


import cantera as ct
ct.add_directory('/user/chinahg')
ct.__file__

import numpy as np
import time
import math as math
#import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize

from combustion_chamber import *
from shocks import *
from nozzle import *

import yaml 
loader=Loader

if __name__ == '__main__':

    stream = open("foo.yaml", 'r')
    dictionary = yaml.safe_load(stream)
#    for key, value in dictionary.items():
#        print (key + " : " + str(value))
print(dictionary['altitude']) #prints value assigned to variable altitude in yaml file

#############################################################################
### LAUNCH VEHICLE STATS ###
#############################################################################

#Fuel and Oxidizer
fuel_comp = 'H2'
ox_comp = 'O2'
mdot_f = 67.35 #[kg/s] https://s3.amazonaws.com/www.tomnoyes.com/shuttle/SSME.jpg
mdot_ox = 404.79 #[kg/s] LOx https://s3.amazonaws.com/www.tomnoyes.com/shuttle/SSME.jpg
mdotCC = mdot_f + mdot_ox
rho_ox = 1141.7 #[kg/m3] LOx
rho_f = 71 #[kg/m3] LH2
#rho_CCavg = (rho_ox+rho_f)/2

#Combustion Chamber Geometry
r_CC = 0.450596/2 #Combustion chamber radius [m]
A_CC = 3.1415*r_CC**2 #Cross sectional area of CC [m2]
L_CC = 0.37338 #Length of combustion chamber [m]
V_CC = A_CC*L_CC #Volume of CC [m3] (assuming cylinder, not tapered)

#Thermodynamic Conditions
P_CC0 = 2.07*10**7 #Chamber pressure [Pa]

#Environmental Conditions
h = 18000 #altitude [m]
P_atm = 7565 #change to fxn of altitude [Pa]


#############################################################################
### COMBUSTION ###
#############################################################################
state = combustion_chamber(P_CC0, V_CC, mdot_ox, mdot_f)

#############################################################################
### CALL NOZZLE REACTOR FUNCTION ###
#############################################################################
# Combustion Chamber Exit to Throat # 
#############################################################################

#PFR (lagrangian, no surface rxns)
#input parameters are those at station 1 (CC1)
n = len(state.T)-1
T_Noz1 = state.T[n]
P_Noz1 = state.P[n] #OR known pressure at CC exit
rho_Noz1 = state.density[n]
comp_Noz1 = state.X[n]
mdot_Noz = mdot_f+mdot_ox

#Nozzle Geometry
L_Noz = 3.0734 #[m]
A_throat = 0.0599 #[m] PLACEHOLDER
A_exit = 4.1317 #[m] PLACEHOLDER

#Call nozzle function
Noz_states = nozzle_react(T_Noz1, P_Noz1, comp_Noz1, A_throat, A_exit, L_Noz, mdot_ox, mdot_f)

#############################################################################
### SHOCKS/EXPANSION ###
#############################################################################

### 1-2 ###
M1 = 3 #need to update
#n = len(Noz.states)
P1 = state.P[0] #Noz_states.P[n]
T1 = state.T[0] #Noz_states.T[n]
P2 = P_atm #Pa

gasPlume4 = shock_calc(M1, P1, T1, P2)
#update with gasCC composition

print(gasPlume4.report())
