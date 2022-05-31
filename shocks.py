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

def shock_calc(M1, P1, T1, P2):
    gamma = 1.1
    P2P1 = P2/P1 #p2/p1
    M1n = M1
    M2n = math.sqrt(((gamma-1)*M1n**2 +2)/(2*gamma*M1n**2 - (gamma-1))) #normal shock equation
    print("M1n = ", M1n)
    print("M2n = ", M2n)

    beta1rad = math.asin(M1n/M1) #1st shock wave angle [rad]
    beta1deg = beta1rad * 180/math.pi #rad to deg
    
    #theta-beta-mach relation
    theta1rad = math.atan(2*(math.tan(beta1rad)**-1)*((M1n**2 - 1)/(M1**2 * (gamma+math.cos(2*beta1rad)) +2))) #1st shock wave deflection angle [deg]
    theta1deg = theta1rad*180/math.pi
    
    M2 = M2n/math.sin(beta1rad-theta1rad) #Mach number after 1st shock
    T2T1 = ((2*gamma*M1n**2 -(gamma-1))*((gamma-1)*M1n**2 +2))/((gamma+1)**2*M1n**2) #1.29 #T2/T1 from normal shock
    T2 = T1*T2T1
    
    beta2rad = sp.optimize.newton(TBM_Beta, 0.25, args=(theta1rad,M2,gamma))

    ### 2-3 ###
    beta2deg = beta2rad*180/math.pi
    M2n = M2*math.sin(beta2rad)
    M3n = math.sqrt(((gamma-1)*M2n**2 +2)/(2*gamma*M2n**2 - (gamma-1))) #normal shock table
    P3P2 = (2*gamma*M2n**2 -(gamma-1))/(gamma+1) #1.86 #normal shock table
    T3T2 = ((2*gamma*M2n**2 -(gamma-1))*((gamma-1)*M2n**2 +2))/((gamma+1)**2*M2n**2) #1.2 #normal shock table
    M3 = M3n/math.sin(beta2rad-theta1rad)
    P3 = P3P2*P2
    T3 = T3T2*T2
    
    ### 3-4 ###
    P4 = P2
    theta4deg = theta1deg #check why, make4 sure should be in deg not rad
    theta4rad = theta1rad
    nuM3 = math.sqrt((gamma+1)/(gamma-1)) *math.atan(math.sqrt(((gamma-1)/(gamma+1))*(M3**2 -1))) - math.atan(math.sqrt(M3**2 -1)) #7.56 #Prandtl-Meyer function
    nuM4 = nuM3 + theta4rad
    
    M4 = sp.optimize.newton(Prandtl_Meyer_Mach, 1.5, args=(nuM4, gamma)) #Prandtl-Meyer TABLE ****************** USER INPUT DURING RUN, NEED TABLE
    T4T04 = (1 + (gamma-1)/2 * M4**2)**-1 #T4/T04 from isentropic relations
    T03T3 = (1 + (gamma-1)/2 * M3**2) #T03/T3 from isentropic relations
    P3P03 = (T03T3**-1)**(gamma/(gamma-1))
    P03 = (P3P03/P3)**-1
    P4P04 = (T4T04)**(gamma/(gamma-1))
    P04 = (P4P04/P4)**-1
    T04T03 = (P04/P03)**((gamma-1)/gamma)
    T4 = T3*T4T04*T04T03*T03T3
    
    #Set new gas state
    gasPlume4 = ct.Solution('gri30.yaml')
    gasPlume4.TP = T4, P4
    
    print("M2:",M2)
    print("P2 [Pa]:",P2)
    print("T2 [K]:",T2)
    print("Theta2 [deg]:",theta1deg, "\n")
    
    print("M3:",M3)
    print("P3 [Pa]:",P3)
    print("T3 [K]:",T3)
    print("Theta3 [deg]:",theta1deg, "\n")
    
    print("M4:",M4)
    print("P4 [Pa]:",P4)
    print("T4 [K]:",T4)
    print("Theta4 [deg]:",theta4deg)
        
    return(gasPlume4)

def Prandtl_Meyer_Mach(x, nuM4, gamma):
    #returns mach number given the turning angle
    return(math.sqrt((gamma+1)/(gamma-1)) * math.atan(math.sqrt((gamma-1)*(x**2 - 1)/(gamma+1))) - math.atan(math.sqrt(x**2-1)) - nuM4)

def TBM_Beta(x, theta1rad, M2, gamma):
    #returns beta given theta and mach number
    return ( (((2*math.cos(x)/math.sin(x))*(M2**2*math.sin(x)**2-1))/(M2**2*(gamma+math.cos(2*x))+2)) - math.tan(theta1rad) )