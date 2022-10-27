import sys
import os
from warnings import catch_warnings

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

def shock_calc(M1, P1, T1, P_a):
    print("M1 = ", M1, "\n")
    print("P1 = ", P1, "\n")
    print("T1 = ", T1, "\n")
    print("P_a = ", P_a, "\n")

    gamma = 1.1

    #1-2 expansion fan (underexpanded)
    nuM1 = 0
    nuM2 = nuM1 + math.sqrt((gamma+1)/(gamma-1)) *math.atan(math.sqrt(((gamma-1)/(gamma+1))*(M1**2 -1))) - math.atan(math.sqrt(M1**2 -1)) #7.56 #Prandtl-Meyer function
    M2 = sp.optimize.newton(Prandtl_Meyer_Mach, 1.5, args=(nuM2, gamma)) #Prandtl-Meyer
    #P2P1 = P2/P1 #((1+(gamma-1)*0.5*(M1**2))/(1+(gamma-1)*0.5*(M2**2)))**(gamma/(gamma-1))
    P2 = P_a
    T2T02 = (1 + (gamma-1)/2 * M2**2)**-1 #T4/T04 from isentropic relations
    T01T1 = (1 + (gamma-1)/2 * M1**2) #T03/T3 from isentropic relations
    P1P01 = (T01T1**-1)**(gamma/(gamma-1))
    P01 = (P1P01/P1)**-1
    P2P02 = (T2T02)**(gamma/(gamma-1))
    P02 = (P2P02/P2)**-1
    T02T01 = (P02/P01)**((gamma-1)/gamma)
    T2 = T1*T2T02*T02T01*T01T1
    
    print("After exp fan")
    print("M2 = ", M2, "\n")
    print("P2 = ", P2, "\n")
    print("T2 = ", T2, "\n")

    #2-3 shock 1
    P3 = P_a
    P3P2 = P3/P2 #p2/p1
    M2n = math.sqrt((P3P2 + (gamma-1)/(gamma+1))*((gamma+1)/(2*gamma)))
    M3n = math.sqrt(((gamma-1)*M2n**2 +2)/(2*gamma*M2n**2 - (gamma-1))) #normal shock equation
#   
    beta2rad = math.asin(M2n/M2) #1st shock wave angle [rad]
    beta2deg = beta2rad * 180/math.pi #rad to deg

    #theta-beta-mach relation
    theta2rad = math.atan(2*(math.tan(beta2rad)**-1)*((M2n**2 - 1)/(M2**2 * (gamma+math.cos(2*beta2rad)) +2))) #1st shock wave deflection angle [deg]
    theta2deg = theta2rad*180/math.pi


    M3 = M3n/math.sin(beta2rad-theta2rad) #Mach number after 1st shock
    T3T2 = P3P2*((2*gamma*M2n**2)/(gamma+1))# -(gamma-1))*((gamma-1)*M2n**2 +2))/((gamma+1)**2*M2n**2) #1.29 #T2/T1 from normal shock
    T3 = T2*T3T2  

    print("After shock")
    print("M3 = ", M3, "\n")
    print("P3 = ", P3, "\n")
    print("T3 = ", T3, "\n")
   #### 2-3 ###
   #beta2deg = beta2rad*180/math.pi
   #M2n = M2*math.sin(beta2rad)
   #M3n = math.sqrt(((gamma-1)*M2n**2 +2)/(2*gamma*M2n**2 - (gamma-1))) #normal shock table

   #P3P2 = (2*gamma*M2n**2 -(gamma-1))/(gamma+1) #normal shock table
   #T3T2 = ((2*gamma*M2n**2 -(gamma-1))*((gamma-1)*M2n**2 +2))/((gamma+1)**2*M2n**2) #1.2 #normal shock table
   #M3 = abs(M3n/math.sin(beta2rad-theta1rad))
   #
   #P3 = P3P2*P2
   #T3 = T3T2*T2
    
    #Set new gas state
    gasPlume4 = ct.Solution('h2o2.yaml')
    gasPlume4.TP = T3, P3

    state = ct.SolutionArray(gasPlume4, extra=['M4'])
    state.append(gasPlume4.state,M4=M3)
    
    #print("M2:",M2)
    #print("P2 [Pa]:",P2)
    #print("T2 [K]:",T2)
    #print("Theta2 [deg]:",theta1deg, "\n")
    #
    #print("M3:",M3)
    #print("P3 [Pa]:",P3)
    #print("T3 [K]:",T3)
    #print("Theta3 [deg]:",theta1deg, "\n")
    #
    #print("M4:",M4)
    #print("P4 [Pa]:",P4)
    #print("T4 [K]:",T4)
    #print("Theta4 [deg]:",theta4deg)
        
    return(state)

def Prandtl_Meyer_Mach(x, nuM4, gamma):
    #returns mach number given the turning angle
    return(math.sqrt((gamma+1)/(gamma-1)) * math.atan(math.sqrt((gamma-1)*(x**2 - 1)/(gamma+1))) - math.atan(math.sqrt(x**2-1)) - nuM4)

def TBM_Beta(x, theta1rad, M2, gamma):
    #returns beta given theta and mach number
    return ( (((2*math.cos(x)/math.sin(x))*(M2**2*math.sin(x)**2-1))/(M2**2*(gamma+math.cos(2*x))+2)) - math.tan(theta1rad) )