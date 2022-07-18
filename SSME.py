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
import yaml
import shutil

from combustion_chamber import *
from shocks import *
from nozzle import *

#For defining h5py groups
Gx = ["G1"]
#delete old files
dir = 'rockettests'
for f in os.listdir(dir):
    shutil.rmtree(os.path.join(dir, f))

if os.path.exists("plot_data.h5"):
    os.remove("plot_data.h5")

loader=yaml.Loader

#altitudes = [16000]
upper = 40000 #[m]
lower = 16000 #[m]
space = int((upper-lower)/250)+1 #250

altitudes = np.linspace(16000, 40000, space, dtype = int)

w = 0
while w < len(altitudes):
    original_umask = os.umask(0)
    newpath = r'/home/chinahg/GCresearch/rocketemissions/rockettests/'+str(altitudes[w])+'m'

    if not os.path.exists(newpath):
        os.mkdir(newpath, 0o755)
    w = w + 1

n_species = 53

#Define variables to save for each altitude
ambient_T = 0
ambient_P = 0
ambient_X = np.zeros(n_species)

results_T = np.zeros(3)
results_P = np.zeros(3)
results_X = np.zeros((3,n_species))
results_u = np.zeros(3)

g = 0

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
    #P_atm = 101325*(1 - 2.25577*(10**-5) * h)**5.25588
     
    if h >= 25000:
        T_atm = (-131.21 + 0.00299*h) + 273.14 #[K]
        P_atm = (2.488 * ((T_atm)/ 216.6)**-11.388)*1000 #[Pa]
        ambient_T = T_atm
        ambient_P = P_atm
        
    elif 11000 < h < 25000:
        T_atm = -56.46 + 273.14 #[K]
        P_atm = (22.65*10**(1.73-0.000157*h))*1000
        ambient_T = T_atm
        ambient_P = P_atm

    else:
        print("ERROR: OUT OF ALTITUDE RANGE")

    #calculate exit velocity
    #Nozzle Geometry
    Area_ratio = 19.8
    A_throat = 0.0599 #[m]
    A_exit = A_throat*Area_ratio #[m]

    F = 2188080 #[N] thrust at 104.5% RPL
    P_e = 13798.5 #[Pa] exit pressure of nozzle
    u_e = (F-((P_e-ambient_P)*A_exit))/(mdot_f+mdot_ox)
    print(u_e)

    #Calculate composition of O2, and N2
    X_N2 = 0.78084
    X_O2 = 0.2095

    ambient_X = [X_N2, X_O2]

    #SAVE VARIABLES FOR DOCUMENTATION
    dictionary = [{'Altitude [m]' : [float(h)],
    'Temperature [K]' : [float(T_atm)], 'Pressure [Pa]' : [float(P_atm)], 
    'Nitrogen Mole Fraction' : [X_N2], 'Oxygen Mole Fraction' : [X_O2]}]

    #Creating YAML file to save conditions
    with open(r'rockettests/'+str(h)+'m/'+str(h)+'_altitude.yaml', 'w') as file:
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
    results_u[0] = (mdot_f+mdot_ox)/(state.density[n]*A_CC)

    #SAVE NEW STATES TO YAML
    #Load YAML file to append new data
    stream = open("rockettests/"+str(h)+"m/"+str(h)+"_altitude.yaml", 'r')
    dictionary = yaml.unsafe_load(stream)
    dictionary.append({'Combustion Chamber Mechanism':['gri30']})
    dictionary.append({'Combustion Chamber Exit Temperature [K]':[float(state.T[n])]})
    dictionary.append({'Combustion Chamber Exit Pressure [Pa]':[float(state.P[n])]})

    stream2 = open('gri30.yaml', 'r')
    gri30 = yaml.unsafe_load(stream2)
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

    L_Noz = (math.sqrt(1/3.1415 *math.sqrt(A_throat*Area_ratio)) - (1/3.1415 *math.sqrt(A_throat)))**2 #[m]
    print(L_Noz)
    #Call nozzle function
    Noz_states = nozzle(T_Noz1, P_Noz1, comp_Noz1, A_throat, A_exit, L_Noz, mdot_ox, mdot_f)
    n = len(Noz_states.T)-1

    #save Nozzle final state
    results_T[1] = Noz_states.T[n]
    results_P[1] = Noz_states.P[n]
    results_X[1,:] = Noz_states.X[n]
    results_u[1] = u_e #[m/s] from F = m_dot*u_e (rocket equation)

    #SAVE NEW STATES TO YAML
    #Load YAML file to append new data
    stream = open("rockettests/"+str(h)+"m/"+str(h)+"_altitude.yaml", 'r')
    dictionary = yaml.unsafe_load(stream)
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

    #Area function
    A_throat = 0.0599
    r_in = math.sqrt(A_throat)/3.1415
    dAdx = np.zeros(len(Noz_states.x)-1)
    A = np.zeros(len(Noz_states.x)-1)
    t = 0
    for t in range(n):
        A[t] = 3.1415*(np.sqrt(t)+r_in)**2

    #############################################################################
    ### SHOCKS/EXPANSION ###
    #############################################################################

    ### 1-2 ###
    n = len(Noz_states.T)-1
    
    u = results_u[1] #velocity at exit of nozzle (m/s)
    gamma = 1.1
    print(u)
    P1 = P_e #Noz_states.P[n] PLACEHOLDER, need design exit pressure or design alt
    T1 = Noz_states.T[n]
    M1 = u/math.sqrt(gamma*P1/Noz_states.density[n])
    P2 = P_atm #Pa

    gasPlume4 = shock_calc(M1, P1, T1, P2) #returns gas object
    #n = (so can plot intermediate states)
    
    #save shock final state
    R = 8.3145 #universal gas constant
    results_T[2] = gasPlume4.T[0]
    results_P[2] = gasPlume4.P[0]
    results_X[2,:] = gasPlume4.X[0]
    results_u[2] = gasPlume4.M4[0]*math.sqrt(gamma*gasPlume4.P[0]/gasPlume4.density[0])

    #SAVE NEW STATES TO YAML
    #Load YAML file to append new data
    stream = open("rockettests/"+str(h)+"m/"+str(h)+"_altitude.yaml", 'r')
    dictionary = yaml.unsafe_load(stream)
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

   
    #Save states for altitude
    with h5py.File('plot_data.h5', 'a') as hdf:
        Gx[g] = hdf.create_group(str(altitudes[g])+'m')
        Gx[g].create_dataset('T_a', data = ambient_T)
        Gx[g].create_dataset('P_a', data = ambient_P)
        Gx[g].create_dataset('X_a', data = ambient_X)

        Gx[g].create_dataset('T', data = results_T)
        Gx[g].create_dataset('P', data = results_P)
        Gx[g].create_dataset('X', data = results_X)
        Gx[g].create_dataset('u', data = results_u)

    g = g+1
    Gx += ["G" + str(g)]

print("done! :)")