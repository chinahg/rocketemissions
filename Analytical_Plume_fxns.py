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
from scipy import integrate
import matplotlib.pyplot as plt

#ZEF
def dBdx(t,B): #(x,B_star)
   beta_s = 0.154
   U_s = abs((-1 + np.sqrt(1 - 4/(-np.pi*B**2))/2))
   fB = beta_s*(U_s/(U_s + 1))
   return fB

def diffusion(nx, nr, xmax, r_ZFE, r_ZEF, nT, T_a, T0, flag):
   """
   Returns the velocity field and distance for 1D linear convection
   """
   if flag == "ZEF":
      r = r_ZEF
   elif flag == "ZFE":
      r = r_ZFE
   else:
      return(print("Invalid domain specified: must be ZFE or ZEF of type string"))

   # Increments
   dr = np.zeros(nx)
   dx = xmax/(nx-1)

   for i in range(nx):
      dr[i] = r[i,np.size(r[i,:])-1]/(nr-1)
   

   # Initialise data structures
   T = np.zeros((nr,nx))

   # Boundary conditions
   T[0,:] = T[nr-1,:] = T_a
    
   # Initial conditions
   for i in range(1,nr-1):
      rmax = r[i,np.size(r[i,:])-1]
      if(r[i-1, rmax] + dr > (rmax/2)-1.45 and r[i-1, rmax] + dr < (rmax/2)+1.45):
         T[i,0] = T0
      else:
         T[i,0] = T_a

   # Loop
   for n in range(0,nx-1):
      for i in range(0,nr-1):
         T[i,n+1] = T[i,n] + nT*(dx/dr**2.0)*(T[i+1,n]-2.0*T[i,n]+T[i-1,n])

   return T, r

def plot_diffusion(T,r,nx,xmax,title):
   """
   Plots the 1D velocity field
   """
   import matplotlib.pyplot as plt
   import matplotlib.cm as cm
   plt.figure()
   colour=iter(cm.rainbow(np.linspace(0,xmax,nx)))
   for i in range(0,nx):
      c=next(colour)
      plt.plot(r-100,T[:,i],c=c)
   plt.xlabel('r (m)')
   plt.ylabel('T (K)')
   plt.title(title)
   plt.show()

def compute_ZFE(u_ZFE, c_ZFE, u_ZFE0, c_ZFE0, r_ZFE, x_ZFE, beta_ZFE, lmbda_ZFE, D, nx_ZFE, nr, xmax_ZFE, A_c, A_u, U_a):
   #ZFE
   b = np.zeros(nx_ZFE) #plume radius
   R = np.zeros(nx_ZFE) #core radius

   sigma_ZFEu = np.zeros(nx_ZFE)
   sigma_ZFEc = np.zeros(nx_ZFE)

   A_u_sq = np.zeros(nx_ZFE)
   A_c_sq = np.zeros(nx_ZFE)

   #u[x,r], c[x,r]
   for i in range(0,nx_ZFE):
       b[i] = beta_ZFE*x_ZFE[i] + (D/2)
       r_ZFE[i,:] = np.linspace(0,b[i]+50,nr)
       R[i] = x_ZFE[i]*(-D/2)/xmax_ZFE + (D/2)  #core width fxn (y = mx + b)

       A_u_sq[i] = R[i]*u_ZFE0 #rectangular area
       A_c_sq[i] = R[i]*c_ZFE0

       sigma_ZFEu[i] = 2*(A_u-A_u_sq[i])/(u_ZFE0*np.sqrt(2*np.pi)) #WHERE FROM????
       sigma_ZFEc[i] = 2.35*(A_c-A_c_sq[i])/(c_ZFE0*np.sqrt(2*np.pi))

       for j in range(0,nr):
           if r_ZFE[i,j] < R[i]:
               u_ZFE[i,j] = u_ZFE0
               c_ZFE[i,j] = c_ZFE0
           else:
               u_ZFE[i,j] = u_ZFE0*np.exp((-(r_ZFE[i,j]-R[i])**2)/(2*sigma_ZFEu[i]**2))
               c_ZFE[i,j] = c_ZFE0*np.exp((-(r_ZFE[i,j]-R[i])**2)/(lmbda_ZFE**2 * sigma_ZFEc[i]**2))
   return u_ZFE, c_ZFE, r_ZFE

def compute_ZEF(U0, C0, U_a, r0, nx_ZFE, nx_ZEF, r_ZEF, xmax_ZFE, xmax_ZEF, nr, A_c, A_u):
   A0 = np.pi*(r0)**2
   Q0 = (U0+U_a)*A0
   B0 = r0 #initial plume radius
   M_e = np.pi*B0**2 *U0*(U0-U_a)
   lm_star = np.sqrt(np.sqrt(M_e)/U_a)
   B_star = np.zeros(nx_ZEF)
   U_star = np.zeros(nx_ZEF)
   U_star[0] = (U0-U_a)/U_a

   u_ZEF = np.zeros((nx_ZEF,nr))
   c_ZEF = np.zeros((nx_ZEF,nr))

   B_star[0] = np.sqrt((np.pi*(U_star[0]**2 + U_star[0]))**(-1)) #estimated initial plume width
   xmax = xmax_ZFE + xmax_ZEF
   x = np.linspace(0, xmax, nx_ZEF)
   x1 = np.linspace(xmax_ZFE,xmax_ZEF,nx_ZEF)

   #Solve for B*(x)
   sol =  sp.integrate.solve_ivp(dBdx,[xmax_ZFE,xmax_ZEF], [B_star[0]], method='RK45', t_eval=x1)

   for i in range(0,nx_ZEF):
       U_star[i] = abs((-1 + np.sqrt(1 - 4*(1/(-np.pi*sol.y[0,i]**2))))/2)

   #reformat and save specific solution variables
   B_star = sol.y[0,:]
   x_ZEF = sol.t

   #Compute B(x), U(x), C(x)
   B = np.zeros(nx_ZEF)
   U = np.zeros(nx_ZEF)
   Sc = np.zeros(nx_ZEF)
   C_m = np.zeros(nx_ZEF)
   k = 0.174
   
   for i in range(0,nx_ZEF):
       B[i] = B_star[i]*lm_star
       U[i] = U_star[i]*U_a + U_a
       
       Sc[i] = k*(x1[i]-x1[0])/lm_star +1
       C_m[i] = C0/Sc[i]
       r_ZEF[i,:] = np.linspace(0,B[i]+50,nr)
   C_m[0] = C0

   #Compute u(x,r) and c(x,r)
   miu_u = 0
   miu_c = 0

   #i = x, j = r
   sigma_ZEFu = np.zeros(nx_ZEF)
   sigma_ZEFc = np.zeros(nx_ZEF)

   u_ZEF0 = np.zeros(nx_ZEF)
   c_ZEF0 = np.zeros(nx_ZEF)
   u_ZEF0[0] = U0
   c_ZEF0[0] = C0

   for i in range(nx_ZEF-1):
       sigma_ZEFu[i] = (A_u)/(U[i]*np.sqrt(2*np.pi)) #WHERE FROM????
       sigma_ZEFc[i] = (A_c)/(C_m[i]*np.sqrt(2*np.pi))
       for j in range(nr):
           u_ZEF[i,j] = U[i]*np.exp((-(r_ZEF[i,j]-miu_u)**2)/(2*sigma_ZEFu[i]**2))
           c_ZEF[i,j] = C_m[i]*np.exp((-(r_ZEF[i,j]-miu_c)**2)/(2*sigma_ZEFc[i]**2))

#           if j == 0:
#                u_ZEF0[i+1] = u_ZEF[i,0]
#                c_ZEF0[i+1] = c_ZEF[i,0]
   print(sigma_ZEFc)
   return u_ZEF, c_ZEF, r_ZEF

def solve_reaction(X, T, P, x_1, x_2, u_1, u_2, gas, j, n_species):

    reactor = ct.IdealGasConstPressureReactor(gas)

    #FOR NO REACTIONS
    #gas.set_multiplier(0)

    states = ct.SolutionArray(gas)
    X_new = np.zeros_like(X)

    for i in range(1,np.size(X[:, 1])): #index through all "r"s

        try:
            gas.TPX = T[i], P, X[i, :] #X[r value, species index]
        except:
            print(i)
            print(gas.report())

        reactor.syncState()
        reactorNet = ct.ReactorNet([reactor])
        t_final = u

        reactorNet.advance(t_final, apply_limit=false)
        
        states.append(reactor.thermo.state)

        #save new composition for that r position
        X_new[i, :] = 10**6 * reactor.thermo.X #mole fraction to ppm state.X[len,:] 10^9 #kmol/m^3s, assume 1 m^3, to ppm #rates for specific y (i) and all species UNITS

    return X_new, states #for all y and all species [s,n_species]