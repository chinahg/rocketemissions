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

from combustion_chamber import *

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
rho_CCavg = (rho_ox+rho_f)/2

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
### RUN CC FXN ###
#############################################################################
state = combustion_chamber(P_CC0, mdot_ox, mdot_f)