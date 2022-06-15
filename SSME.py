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

import h5py
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
loader=yaml.SafeLoader

altitudes = [16000, 20000, 24000, 28000, 32000, 36000, 40000]
n_species = 53

#Define variables to save for each altitude
ambient_T = 0
ambient_P = 0
ambient_X = np.zeros(n_species)

results_T = np.zeros(3)
results_P = np.zeros(3)
results_X = np.zeros((3,n_species))
g=0

while g<len(altitudes):

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
    h = altitudes[g] #altitude [m]
    #######################################################################

    if h > 25000:
        T_atm = -131.21 + 0.00299*h + 273.14 #[K]
        P_atm = 2.488*(T_atm/216.6)*100 #[Pa]
        ambient_T = T_atm
        ambient_P = P_atm
        
    elif 11000 < h < 25000:
        T_atm = -56.46 + 273.14 #[K]
        P_atm = (22.65*math.e**(1.73-0.000157*h))*1000
        ambient_T = T_atm
        ambient_P = P_atm

    else:
        print("ERROR: OUT OF ALTITUDE RANGE")

    #Calculate composition of O2, and N2
    X_N2 = 0.78084
    X_O2 = 0.2095

    P_N2 = P_atm*X_N2
    P_O2 = P_atm*X_O2

    ambient_X = [X_N2, X_O2]

    #SAVE VARIABLES FOR DOCUMENTATION
    dictionary = [{'Altitude [m]' : [h],
    'Temperature [K]' : [T_atm], 'Pressure [Pa]' : [P_atm], 
    'Nitrogen Mole Fraction' : [X_N2], 'Oxygen Mole Fraction' : [X_O2]}]

    #Creating YAML file to save conditions
    with open("rockettests/"+str(h)+"m/"+str(h)+"_altitude.yaml", 'w') as file:
        documents = yaml.dump(dictionary, file)

    #############################################################################
    ### COMBUSTION ###
    #############################################################################
    state = combustion_chamber(P_CC0, V_CC, mdot_ox, mdot_f)
    n = len(state.T) - 1

    #save CC state
    results_T[0] = state.T[n]
    results_P[0] = state.P[n]
    results_X[0,:] = state.X[n]

    #SAVE NEW STATES TO YAML
    #Load YAML file to append new data
    stream = open("rockettests/"+str(h)+"m/"+str(h)+"_altitude.yaml", 'r')
    dictionary = yaml.safe_load(stream)
    dictionary.append({'Combustion Chamber Mechanism':['gri30']})
    dictionary.append({'Combustion Chamber Exit Temperature [K]':[float(state.T[n])]})
    dictionary.append({'Combustion Chamber Exit Pressure [Pa]':[float(state.P[n])]})

    stream2 = open('gri30.yaml', 'r')
    gri30 = yaml.safe_load(stream2)
    gri30_species = gri30['phases']
    gri30_species = gri30_species[0]['species']

    i=0
    while i < n_species:
        if state.X[n][i] != 0:
            dictionary.append({str(gri30_species[i]):[float(state.X[n][i])]})
        i = i+1

    #Save new dictionary to YAML file
    with open("rockettests/"+str(h)+"m/"+str(h)+"_altitude.yaml", 'w') as file:
        documents = yaml.dump(dictionary, file)

    #PLOT TEMPERATURE
    plt.figure()
    plt.plot(state.t,state.T, label = 'Combustion Chamber Temperature')
    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [K]")
    plt.savefig("rockettests/"+str(h)+"m/CC_T.png")

    #PLOT PRESSURE
    plt.figure()
    plt.plot(state.t,state.P, label = 'Combustion Chamber Pressure')
    plt.xlabel("Time [s]")
    plt.ylabel("Pressure [Pa]")
    plt.savefig("rockettests/"+str(h)+"m/CC_P.png")

    #PLOT MOLE FRACTION
    plt.figure()

    plt.plot(state.t,state("O2").X, label="O2") #O2
    plt.plot(state.t,state("H2").X, label="H2") #H2
    plt.plot(state.t,state("H2O").X, label="H2O") #H2
    plt.plot(state.t,state("OH").X, label="OH") #H2

    plt.xlabel("Time [s]")
    plt.ylabel("Mole Fraction")
    plt.legend()
    plt.savefig("rockettests/"+str(h)+"m/CC_X.png")

    #PLOT ENTHALPY
    plt.figure()

    plt.plot(state.t,state.enthalpy, label="enthalpy in reactor")
    plt.plot(state.t,state.eox, label="enthalpy ox")
    plt.plot(state.t,state.efuel, label="enthalpy fuel")
    plt.xlabel("Time [s]")
    plt.ylabel("enthalpy [J/mole]")
    plt.legend()
    plt.savefig("rockettests/"+str(h)+"m/CC_H.png")

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

    #save Nozzle final state
    results_T[1] = Noz_states.T[n]
    results_P[1] = Noz_states.P[n]
    results_X[1,:] = Noz_states.X[n]

    #SAVE NEW STATES TO YAML
    #Load YAML file to append new data
    stream = open("rockettests/"+str(h)+"m/"+str(h)+"_altitude.yaml", 'r')
    dictionary = yaml.safe_load(stream)
    dictionary.append({'Nozzle Mechanism':['gri30']})
    dictionary.append({'Nozzle Exit Temperature [K]':[float(Noz_states.T[n])]})
    dictionary.append({'Nozzle Exit Pressure [Pa]':[float(Noz_states.P[n])]})

    i=0
    while i < n_species:
        if Noz_states.X[n][i] != 0:
            dictionary.append({str(gri30_species[i]):[float(Noz_states.X[n][i])]})
        i = i+1

    #Save new dictionary to YAML file
    with open("rockettests/"+str(h)+"m/"+str(h)+"_altitude.yaml", 'w') as file:
        documents = yaml.dump(dictionary, file)

    #PLOT AREA(X)
    drdx = 0.5
    dAdx = 3.1415*2*drdx

    A = A_throat+dAdx*Noz_states.x

    plt.figure()
    L1 = plt.plot(Noz_states.x, A, color='r', label='A', lw=2)
    plt.xlabel('distance (m)')
    plt.ylabel('Area')
    plt.savefig("rockettests/"+str(h)+"m/nozzle_area.png")

    #PLOT MOLE FRACTIONS
    plt.figure()
    fig, ax = plt.subplots()
    plt.ylabel('Mole Fraction')
    ax.plot(Noz_states.x, Noz_states('H2O').X, 'k--', label='H2O')
    ax.plot(Noz_states.x, Noz_states('OH').X, 'k:', label='OH')
    ax.plot(Noz_states.x, Noz_states('H2').X, 'k', label='H2')
    ax.plot(Noz_states.x, Noz_states('O2').X, 'k-.', label='O2')

    legend = ax.legend(loc='center right', shadow=True, fontsize='x-large')
    plt.savefig("rockettests/"+str(h)+"m/nozzle_X.png")

    #PLOT TEMPERATURE VS ISENTROPIC
    plt.figure()
    L1 = plt.plot(Noz_states.x, Noz_states.T, color='r', label='T', lw=2)
    plt.xlabel('distance (m)')
    plt.ylabel('Temperature (K)')

    plt.legend(L1, [line.get_label() for line in L1], loc='lower right')
    plt.savefig("rockettests/"+str(h)+"m/nozzle_T.png")

    #PLOT PRESSURE
    plt.figure()
    L1 = plt.plot(Noz_states.x, Noz_states.P, color='r', label='P', lw=2)
    plt.xlabel('distance (m)')
    plt.ylabel('Pressure (Pa)')

    plt.legend(L1, [line.get_label() for line in L1], loc='lower right')
    plt.savefig("rockettests/"+str(h)+"m/nozzle_P.png")

    #PLOT PRESSURE VS ISENTROPIC

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

    gasPlume4 = shock_calc(M1, P1, T1, P2) #returns gas object
    #n = (so can plot intermediate states)

    #save shock final state
    results_T[2] = gasPlume4.T
    results_P[2] = gasPlume4.P
    results_X[2,:] = gasPlume4.X

    #SAVE NEW STATES TO YAML
    #Load YAML file to append new data
    stream = open("rockettests/"+str(h)+"m/"+str(h)+"_altitude.yaml", 'r')
    dictionary = yaml.safe_load(stream)
    dictionary.append({'Shocks Mechanism':['NONE']})
    dictionary.append({'Shocks Exit Temperature [K]':[float(gasPlume4.T)]})
    dictionary.append({'Shocks Exit Pressure [Pa]':[float(gasPlume4.P)]})
    dictionary.append({'Shocks Exit Velocity [m/s]':[float(u)]})

    #same as nozzle, no chemistry in shocks
    i=0
    while i < n_species:
        if Noz_states.X[n][i] != 0:
            dictionary.append({str(gri30_species[i]):[float(Noz_states.X[n][i])]})
        i = i+1

    #Save new dictionary to YAML file
    with open("rockettests/"+str(h)+"m/"+str(h)+"_altitude.yaml", 'w') as file:
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

    #Save states for altitude
    with h5py.File('plot_data.h5', 'a') as hdf:
        G1 = hdf.create_group(str(altitudes[g])+'m')
        G1.create_dataset('T_a', data = ambient_T)
        G1.create_dataset('P_a', data = ambient_P)
        G1.create_dataset('X_a', data = ambient_X)

        G1.create_dataset('T', data = results_T)
        G1.create_dataset('P', data = results_P)
        G1.create_dataset('X', data = results_X)

    g = g+1

### PLOTTING ###

#Load data from all altitudes
g=0
while g<len(altitudes):
    h = altitudes[g]
    with open('rockettests/'+str(h)+'m/'+str(h)+'_altitude.yaml') as f:
        globals()['dictionary_'+str(h)] = yaml.load(f, Loader=yaml.SafeLoader)
    g = g+1

#PLOT TEMPERATURES (CC, NOZZLE, SHOCKS, EXIT)
plt.ylabel('Fluid Temperature [K]')
"""
T_16000 = [dictionary_16000[2]['Combustion Chamber Exit Temperature [K]'],
            dictionary_16000[13]['Nozzle Exit Temperature [K]'],
            dictionary_16000[24]['Shocks Exit Temperature [K]']]

T_20000 = [dictionary_20000[2]['Combustion Chamber Exit Temperature [K]'],
            dictionary_20000[13]['Nozzle Exit Temperature [K]'],
            dictionary_20000[24]['Shocks Exit Temperature [K]']]

T_24000 = [dictionary_24000[2]['Combustion Chamber Exit Temperature [K]'],
            dictionary_24000[13]['Nozzle Exit Temperature [K]'],
            dictionary_24000[24]['Shocks Exit Temperature [K]']]

T_28000 = [dictionary_28000[2]['Combustion Chamber Exit Temperature [K]'],
            dictionary_28000[13]['Nozzle Exit Temperature [K]'],
            dictionary_28000[24]['Shocks Exit Temperature [K]']]

T_32000 = [dictionary_32000[2]['Combustion Chamber Exit Temperature [K]'],
            dictionary_32000[13]['Nozzle Exit Temperature [K]'],
            dictionary_32000[24]['Shocks Exit Temperature [K]']]

T_36000 = [dictionary_36000[2]['Combustion Chamber Exit Temperature [K]'],
            dictionary_36000[13]['Nozzle Exit Temperature [K]'],
            dictionary_36000[24]['Shocks Exit Temperature [K]']]

T_40000 = [dictionary_40000[2]['Combustion Chamber Exit Temperature [K]'],
            dictionary_40000[13]['Nozzle Exit Temperature [K]'],
            dictionary_40000[24]['Shocks Exit Temperature [K]']]
"""

#PLOT TEMPERATURE AT EACH STATION FOR EACH ALT
steps = [1,2,3]
fig, ax = plt.subplots()
ax.plot(steps, T_16000, label='16 km')
ax.plot(steps, T_20000, label='20 km')
ax.plot(steps, T_24000, label='24 km')
ax.plot(steps, T_28000, label='28 km')
ax.plot(steps, T_32000, label='32 km')
ax.plot(steps, T_36000, label='36 km')
ax.plot(steps, T_40000, label='40 km')

legend = ax.legend(loc='center right', shadow=True, fontsize='large')
plt.savefig("rockettests/multi_alt/agg_CCtemp.png")

#PLOT PLUME INITIAL TEMPERATURE FOR EACH ALT
fig, ax = plt.subplots()
plt.xlabel('Fluid Temperature [K]')
T_shock_all = [T_16000[2],
                T_20000[2],
                T_24000[2],
                T_28000[2],
                T_32000[2],
                T_36000[2],
                T_40000[2]]

ax.plot(altitudes, T_shock_all)
plt.savefig("rockettests/multi_alt/agg_shocktemp.png")

plt.close('all')