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
import matplotlib.pyplot as plt
import scipy as sp
import scipy.optimize

from combustion_chamber import *
from shocks import *
from nozzle import *

import yaml 
loader=Loader

#############################################################################
### LAUNCH VEHICLE STATS ###
#############################################################################

#Fuel and Oxidizer
fuel_comp = 'H2'
ox_comp = 'O2'
mdot_f = 67.35 #[kg/s] https://s3.amazonaws.com/www.tomnoyes.com/shuttle/SSME.jpg
mdot_ox = 404.79 #[kg/s] LOx https://s3.amazonaws.com/www.tomnoyes.com/shuttle/SSME.jpg
mdotCC = mdot_f + mdot_ox #Assume constant mdot (no throttle) for now
#rho_ox = 1141.7 #[kg/m3] LOx
#rho_f = 71 #[kg/m3] LH2

#Combustion Chamber Geometry
r_CC = 0.450596/2 #Combustion chamber radius [m]
A_CC = 3.1415*r_CC**2 #Cross sectional area of CC [m2]
L_CC = 0.37338 #Length of combustion chamber [m]
V_CC = A_CC*L_CC #Volume of CC [m3] (assuming cylinder, not tapered)

#Thermodynamic Conditions
P_CC0 = 2.07*10**7 #Chamber pressure [Pa]

#Environmental Conditions
### USER INPUT HERE FOR ALTITUDE ######################################
h = 20000 #altitude [m]
#######################################################################

if h > 25000:
    T_atm = -131.21 + 0.00299*h + 273.14 #[K]
    P_atm = 2.488*(T_atm/216.6)/1000 #[Pa]
    
elif 11000 < h < 25000:
    T_atm = -56.46 + 273.14 #[K]
    P_atm = 22.65*math.e**(1.73-0.000157*h)

else:
    print("ERROR: OUT OF ALTITUDE RANGE")

#Calculate composition of O2, N2, and H2O
X_N2 = 0.78084
X_O2 = 0.2095

P_N2 = P_atm*X_N2
P_O2 = P_atm*X_O2

A = np.array([[1,-P_atm],[1,0]])
B = np.array([0,P_atm-P_N2-P_O2])
S = np.linalg.solve(A,B)

X_H2O = float(S[1])
P_H2O = float(S[0])

#SAVE VARIABLES FOR DOCUMENTATION
dictionary = [{'Altitude [m]' : [h]},
{'Temperature [K]' : [T_atm]}, {'Pressure [Pa]' : [P_atm]}, 
{'Nitrogen Mole Fraction' : [X_N2]}, {'Oxygen Mole Fraction' : [X_O2]},
{'Water Mole Fraction' : [X_H2O]}]

#Creating YAML file to save conditions
with open(str(h)+"_altitude.yaml", 'w') as file:
    documents = yaml.dump(dictionary, file)

#############################################################################
### COMBUSTION ###
#############################################################################
state = combustion_chamber(P_CC0, V_CC, mdot_ox, mdot_f)
n = len(state.T) - 1

#SAVE NEW STATES TO YAML
#Load YAML file to append new data
stream = open(str(h)+"_altitude.yaml", 'r')
dictionary = yaml.safe_load(stream)
dictionary.append({'Combustion Chamber Mechanism':['gri30']})
dictionary.append({'Combustion Chamber Exit Temperature [K]':[float(state.T[n])]})
dictionary.append({'Combustion Chamber Exit Pressure [Pa]':[float(state.P[n])]})

stream2 = open('gri30.yaml', 'r')
gri30 = yaml.safe_load(stream2)
gri30_species = gri30['phases']
gri30_species = gri30_species[0]['species']

i=0
while i < 53:
    if state.X[n][i] != 0:
        dictionary.append({str(gri30_species[i]):[float(state.X[n][i])]})
    i = i+1

#Save new dictionary to YAML file
with open(str(h)+"_altitude.yaml", 'w') as file:
    documents = yaml.dump(dictionary, file)

#PLOT TEMPERATURE
plt.figure()
plt.plot(state.t,state.T, label = 'Combustion Chamber Temperature')
plt.xlabel("Time [s]")
plt.ylabel("Temperature [K]")
plt.savefig("rockettests/altitude10km/CC_T.png")

#PLOT PRESSURE
plt.figure()
plt.plot(state.t,state.P, label = 'Combustion Chamber Pressure')
plt.xlabel("Time [s]")
plt.ylabel("Pressure [Pa]")
plt.savefig("rockettests/altitude10km/CC_P.png")

#PLOT MOLE FRACTION
plt.figure()

plt.plot(state.t,state("O2").X, label="O2") #O2
plt.plot(state.t,state("H2").X, label="H2") #H2
plt.plot(state.t,state("H2O").X, label="H2O") #H2
plt.plot(state.t,state("OH").X, label="OH") #H2

plt.xlabel("Time [s]")
plt.ylabel("Mole Fraction")
plt.legend()
plt.savefig("rockettests/altitude10km/CC_X.png")

#PLOT ENTHALPY
plt.figure()

plt.plot(state.t,state.enthalpy, label="enthalpy in reactor")
plt.plot(state.t,state.eox, label="enthalpy ox")
plt.plot(state.t,state.efuel, label="enthalpy fuel")
plt.xlabel("Time [s]")
plt.ylabel("enthalpy [J/mole]")
plt.legend()
plt.savefig("rockettests/altitude10km/CC_H.png")

#############################################################################
### CALL NOZZLE REACTOR FUNCTION ###
#############################################################################
# Combustion Chamber Exit to Exit # 
#############################################################################

#PFR (no surface rxns)
#input parameters are those at station 1 (CC1)
n = len(state.T)-1
T_Noz1 = state.T[n]
P_Noz1 = state.P[n]
rho_Noz1 = state.density[n]
comp_Noz1 = state.X[n]
mdot_Noz = mdot_f+mdot_ox

#Nozzle Geometry
L_Noz = 3.0734 #[m]
A_throat = 0.0599 #[m]
A_exit = 4.1317 #[m]

#Call nozzle function
Noz_states = nozzle(T_Noz1, P_Noz1, comp_Noz1, A_throat, A_exit, L_Noz, mdot_ox, mdot_f)
n = len(Noz_states.T)-1

#SAVE NEW STATES TO YAML
#Load YAML file to append new data
stream = open(str(h)+"_altitude.yaml", 'r')
dictionary = yaml.safe_load(stream)
dictionary.append({'Nozzle Mechanism':['gri30']})
dictionary.append({'Nozzle Exit Temperature [K]':[float(Noz_states.T[n])]})
dictionary.append({'Nozzle Exit Pressure [Pa]':[float(Noz_states.P[n])]})

i=0
while i < 53:
    if Noz_states.X[n][i] != 0:
        dictionary.append({str(gri30_species[i]):[float(Noz_states.X[n][i])]})
    i = i+1

#Save new dictionary to YAML file
with open(str(h)+"_altitude.yaml", 'w') as file:
    documents = yaml.dump(dictionary, file)

#PLOT AREA(X)
drdx = 0.5
dAdx = 3.1415*2*drdx

A = A_throat+dAdx*Noz_states.x

plt.figure()
L1 = plt.plot(Noz_states.x, A, color='r', label='A', lw=2)
plt.xlabel('distance (m)')
plt.ylabel('Area')
plt.savefig("rockettests/altitude10km/nozzle_area.png")

#PLOT MOLE FRACTIONS
plt.figure()
fig, ax = plt.subplots()
plt.ylabel('Mole Fraction')
ax.plot(Noz_states.x, Noz_states('H2O').X, 'k--', label='H2O')
ax.plot(Noz_states.x, Noz_states('OH').X, 'k:', label='OH')
ax.plot(Noz_states.x, Noz_states('H2').X, 'k', label='H2')
ax.plot(Noz_states.x, Noz_states('O2').X, 'k-.', label='O2')

legend = ax.legend(loc='center right', shadow=True, fontsize='x-large')
plt.savefig("rockettests/altitude10km/nozzle_X.png")

#PLOT TEMPERATURE
plt.figure()
L1 = plt.plot(Noz_states.x, Noz_states.T, color='r', label='T', lw=2)
plt.xlabel('distance (m)')
plt.ylabel('Temperature (K)')

plt.legend(L1, [line.get_label() for line in L1], loc='lower right')
plt.savefig("rockettests/altitude10km/nozzle_T.png")

#PLOT PRESSURE
plt.figure()
L1 = plt.plot(Noz_states.x, Noz_states.P, color='r', label='P', lw=2)
plt.xlabel('distance (m)')
plt.ylabel('Pressure (Pa)')

plt.legend(L1, [line.get_label() for line in L1], loc='lower right')
plt.savefig("rockettests/altitude10km/nozzle_P.png")

#############################################################################
### SHOCKS/EXPANSION ###
#############################################################################

### 1-2 ###
n = len(Noz_states.T)-1
u = mdot_Noz/(Noz_states.density[n]*A[n]) #velocity at exit of nozzle (m/s)
gamma = 1.1

P1 = Noz_states.P[n]
T1 = Noz_states.T[n]
M1 = u/math.sqrt(gamma*P1/Noz_states.density[n])
P2 = P_atm #Pa
print("T1: \n", T1)
print("\nNOZZLE EXIT VELOCITY: ",u, "m/s\n", "M1 = ", M1)

gasPlume4 = shock_calc(M1, P1, T1, P2) #returns gas object

print("\nINITIAL PLUME CONDITIONS: ", gasPlume4.report())

#SAVE NEW STATES TO YAML
#Load YAML file to append new data
stream = open(str(h)+"_altitude.yaml", 'r')
dictionary = yaml.safe_load(stream)
dictionary.append({'Shocks Mechanism':['NONE']})
dictionary.append({'Shocks Exit Temperature [K]':[float(gasPlume4.T)]})
dictionary.append({'Shocks Exit Pressure [Pa]':[float(gasPlume4.P)]})
dictionary.append({'Shocks Exit Velocity [m/s]':[u]})

#same as nozzle, no chemistry in shocks
i=0
while i < 53:
    if Noz_states.X[n][i] != 0:
        dictionary.append({str(gri30_species[i]):[float(Noz_states.X[n][i])]})
    i = i+1

#Save new dictionary to YAML file
with open(str(h)+"_altitude.yaml", 'w') as file:
    documents = yaml.dump(dictionary, file)

""" #PLOT MOLE FRACTIONS
plt.figure()
fig, ax = plt.subplots()
plt.ylabel('Mole Fraction')
ax.plot(gasPlume4.x, gasPlume4('H2O').X, 'k--', label='H2O')
ax.plot(gasPlume4.x, gasPlume4('OH').X, 'k:', label='OH')
ax.plot(gasPlume4.x, gasPlume4('H2').X, 'k', label='H2')
ax.plot(gasPlume4.x, gasPlume4('O2').X, 'k-.', label='O2')

legend = ax.legend(loc='center right', shadow=True, fontsize='x-large')
plt.savefig("rockettests/altitude10km/shocks_X.png") """