#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 14:24:18 2024

@author: pierre
"""


import sys
sys.path.append('/Users/pierre/Dropbox/OnGoingResearch/Fluid_and_earthquakes/Softwares/DemystiFicatioN_release/wrapper_python')

from DemystiFicatioN import *
import numpy as np
from matplotlib import pyplot as plt
from scipy.special import gamma




# Variable 
fault_length1 = 10.0


# Variable smaller faults
x_min = -8.0
x_max = 8.0
y_min = -4.0
y_max =4.0
x_num = 40
y_num = 20
L_small = 0.1


# Create new problem
nb_fault = x_num*y_num+1
#nb_fault  = 1
#x_num = 0
#y_num = 0
problem = DemystiFicatioN(nb_fault)

#------------------------------------------------------------------------------
# Generate rough fault
#------------------------------------------------------------------------------

##Function that generate fractal (from Tatsuhiko Saito)
def von_Karman(epsi,ax,kap,seedXi,dx,lambda_min):
#
#           ver 1.0   2004.12.27
#           T. Saito, modified by P. Romanet
#                                                 
# 

    # Lenght of the array
    nx=len(seedXi)
    print('{0}'.format(nx))
    # Create wavenumber
    kx = 2*np.pi/(nx*dx)*np.hstack((np.arange(0,nx/2),np.arange( -nx/2,0)));

    # Uniformely distributed phase between [0:2pi]
    pang=2*np.pi*seedXi

    # Von Karman (Sato et al., 2012, p 24)
    yy_c=2*np.pi**(1/2)*gamma(kap+1/2)*epsi**2*ax/(gamma(kap)*(1+(kx*ax)**2)**(kap+0.5))
    # # yy_c1= (2*pi)^3*(1e-4)^2./kx.^3;

    # Realisation of random media (Sato et al., 2012, p 21)
    Xi_kx=np.sqrt(yy_c)*np.exp(1j*pang)*np.sqrt(nx*dx)      # [m*m]

    # Remove mean of the signal
    Xi_kx[0] = 0.0

    # Remove small perturvation
    # Calculate max wavenumber
    kmax = 2*np.pi/lambda_min

    # Remove wavenumber that are higher than max wavenumber
    index = np.where(abs(kx)>kmax)
    Xi_kx[index] = 0

    # ---   Spectrum to Fluctuation
    Xi=np.fft.ifft(Xi_kx)/dx
    Xi = np.real(Xi)


    # Calculate RMS
    alpha_check = np.sqrt(dx*sum(Xi**2)/((nx*dx)))/(nx*dx);
    ratio = epsi/(2*np.sqrt(2)*np.pi*ax)/alpha_check;
    
    # Print some information
    print('Theoretical alpha = {0}, Calculated alpha = {1}'.format(epsi/(2*np.sqrt(2)*np.pi*ax),alpha_check))
    print('Ratio= {0}'.format(ratio))



    return Xi








#------------------------------------------------------------------------------
# Hyperparameter
#------------------------------------------------------------------------------
Dinj  = 7000.0 # original one
#seed = 5  # original one
seed = 5


# Simulation name
problem.hyperparameter.simulation_name = "LSBB_with_800faults_noLoading_D{0}_seed_{1}".format(Dinj,seed)
# Frature mode
problem.hyperparameter.fracture_mode = "modeII"
# Way of calculating the static kernel
problem.hyperparameter.static_kernel = "hmatrix"
# Type of rate and state friction law
problem.hyperparameter.friction_law = "RateStateAging"
# Tolerance for the solver
problem.hyperparameter.tol_solver = 1e-6
# Frequence of writting output file 
problem.hyperparameter.freq_writing_file = 1000# Original parameter
# Number of time step after which the simulation stops     
problem.hyperparameter.max_it = 100000
# Save every stride_time time step
problem.hyperparameter.stride_time = 10 # Original parameter


# Accuracy of Fast Multipole Method (See code by Greengard) 
problem.hyperparameter.iprec = 4

# Permeability_coupling
problem.hyperparameter.permeability_coupling = 0


#------------------------------------------------------------------------------
# Loading and material
#------------------------------------------------------------------------------
problem.material_loading.mu = 32.04e9;
problem.material_loading.cp = 5000.0;
problem.material_loading.cs = 3.464e3;
 
# Uniform stress rate loading
problem.material_loading.Sigma31_dot = 0.0
problem.material_loading.Sigma32_dot = 0.0
problem.material_loading.Sigma11_dot = 0.0
problem.material_loading.Sigma12_dot = 0.0
problem.material_loading.Sigma22_dot = 0.0
 
# Uniform traction rate loading
problem.material_loading.normal_loading_dot = 0.0
problem.material_loading.shear_loading_dot = 0.0000

problem.material_loading.rock_comp = 0.1e-8
problem.material_loading.fluid_comp = 0.9e-8
problem.material_loading.fluid_density = 1000.0
problem.material_loading.porosity = 0.1
problem.material_loading.dyn_viscosity = 1e-3


 
#------------------------------------------------------------------------------ 
# Fault properties 1
#------------------------------------------------------------------------------
# Geometry
problem.faults[0].nb_element = 1999;
 
# Friction
problem.faults[0].a = 0.006*np.ones(problem.faults[0].nb_element)
problem.faults[0].b = 0.005*np.ones(problem.faults[0].nb_element)

# Make VW patches
# problem.faults[0].a[0:1000] = 0.005
# problem.faults[0].b[0:1000] = 0.007
# problem.faults[0].a[-1000:] = 0.005
# problem.faults[0].b[-1000:] = 0.007


print('size: {0}'.format(np.size(problem.faults[0].a[0:1000])))
print('size: {0}'.format(np.size(problem.faults[0].a[-1000:])))

problem.faults[0].Dc = 0.000001*np.ones(problem.faults[0].nb_element)
problem.faults[0].f0 = 0.6*np.ones(problem.faults[0].nb_element)
problem.faults[0].V0 = 1e-6*np.ones(problem.faults[0].nb_element)

# Stresses
problem.faults[0].sigmaN[:] = -5.0e6;
problem.faults[0].normal_loading_dot[:] = 0.0;
problem.faults[0].shear_loading_dot[:] = 0.00;

# Initial parameters (V and theta)
problem.faults[0].V = 1e-12*np.ones(problem.faults[0].nb_element)
problem.faults[0].theta = np.log(problem.faults[0].V0[0]/problem.faults[0].V[0])*np.ones(problem.faults[0].nb_element);
problem.faults[0].P = 0.0*np.ones(problem.faults[0].nb_element);


# Permeability
problem.faults[0].permeability = 1e-15*np.ones(problem.faults[0].nb_element)

# Injection
problem.faults[0].nb_source = 1
problem.faults[0].Q =  1.3416e-6
problem.faults[0].t_injection_beg = 10.0*365.25*86400.0
problem.faults[0].t_injection_end = 10.0*365.25*86400.0+Dinj
problem.faults[0].index_injection = np.floor(problem.faults[0].nb_element/2-200).astype(int)



# # Get estimate for Lnuc
a = problem.faults[0].a[0]
b = problem.faults[0].b[0]
Lb = problem.material_loading.mu*problem.faults[0].Dc[0]/(np.abs(problem.faults[0].sigmaN[0])*b);
            
# # Compute Lnuc from Viesca [2016]
Lnuc = 2*Lb / (np.pi*(1-a/b)**2)
print('Lb (large fault) = {0}'.format(Lb))
print('Lnuc (large fault) = {0}'.format(Lnuc))
 

#------------------------------------------------------------------------------
# Create fault geometry 1
#------------------------------------------------------------------------------
# Amplitude to wavelength ratio
alpha1 = 1e-3 

# Parameter for generation of rough geometry
aa=10000000         # [m], big value make the Von Karman PSD tends to the one in Dunham et al., 2011b (autocorrelation distance)
epsi= 2*np.sqrt(2)*np.pi*aa*alpha1   # Theoretical value of epsilon based on alpha and aa to match Dunham et al., 2011b

# Order of the von Karman
kap=1 # For this value, the Von Karman PSDF asymtotically converge toward P(k)=k^(-3) which is self-similar   

# Calculate lengthscale
ds = fault_length1/problem.faults[0].nb_element


# Minimum wavelength of the roughness
lambda_min = 40*ds;



# Generate geometry
height = np.zeros((1,problem.faults[0].nb_element));
# Random number for the phase
#np.random.seed(5)
np.random.seed(seed)

seedXi=np.squeeze(np.random.rand(1,problem.faults[0].nb_element+1));


height_rough = von_Karman(epsi,aa,kap,seedXi,ds,lambda_min);



problem.faults[0].node[0,:] = np.linspace(-fault_length1/2,fault_length1/2,problem.faults[0].nb_element+1)
problem.faults[0].node[1,:] = height_rough

# Record all the nodes to avoid overlapping
nodex = problem.faults[0].node[0,:]
nodey = problem.faults[0].node[1,:]

#------------------------------------------------------------------------------ 
# Create smaller faults
#------------------------------------------------------------------------------
fault_id = 1
# For every faults
for x_center in np.linspace(x_min,x_max,x_num):
    for y_center in np.linspace(y_min,y_max,y_num):
    
        # Geometry
        problem.faults[fault_id].nb_element = 2;
        problem.faults[fault_id].nb_element = 1;

        # Friction
        problem.faults[fault_id].a = 0.005*np.ones(problem.faults[fault_id].nb_element)
        problem.faults[fault_id].b = 0.007*np.ones(problem.faults[fault_id].nb_element)
        problem.faults[fault_id].Dc = 0.000000004*np.ones(problem.faults[fault_id].nb_element)
        problem.faults[fault_id].f0 = 0.6*np.ones(problem.faults[fault_id].nb_element)
        problem.faults[fault_id].V0 = 1e-6*np.ones(problem.faults[fault_id].nb_element)
    
        # Stresses
        #print(fault_id)
        problem.faults[fault_id].sigmaN = -5.0e6*np.ones(problem.faults[fault_id].nb_element);
        problem.faults[fault_id].normal_loading_dot = 0.0*np.ones(problem.faults[fault_id].nb_element);
        problem.faults[fault_id].shear_loading_dot = 0.000*np.ones(problem.faults[fault_id].nb_element);
    
        # Initial parameters (V and theta)
        problem.faults[fault_id].V = 1e-12*np.ones(problem.faults[fault_id].nb_element)
        problem.faults[fault_id].theta = np.log(problem.faults[0].V0[0]/problem.faults[fault_id].V[0])*np.ones(problem.faults[fault_id].nb_element)
        problem.faults[fault_id].P = np.zeros(problem.faults[fault_id].nb_element);
    
    
        # Permeability
        problem.faults[fault_id].permeability = np.ones(problem.faults[fault_id].nb_element)
        
        
    
        # Injection
        problem.faults[fault_id].nb_source = 0
    
        
        # Create fault geometry
        # Select random length
        L_r = L_small
        
        # Select random orientation
        angle_r = (np.random.rand()*2.0-1)*30/180*np.pi
        
        
        # Find closest x point on the main fault 
        min_id = np.argmin(np.abs(problem.faults[0].node[0,:]- x_center))
        
        
        problem.faults[fault_id].node[0,:] = np.cos(angle_r)*np.linspace(-L_r/2,L_r/2,problem.faults[fault_id].nb_element+1) + x_center
        problem.faults[fault_id].node[1,:] = np.sin(angle_r)*np.linspace(-L_r/2,L_r/2,problem.faults[fault_id].nb_element+1) + y_center + problem.faults[0].node[1,min_id]
        
        
        
    



        fault_id += 1








# # To remove 
# del problem.faults[0]
# problem._nb_fault = x_num*y_num
# nb_fault = x_num*y_num
# for fault_id in range(1,nb_fault):
#       problem.hyperparameter.nb_fault = nb_fault

#------------------------------------------------------------------------------ 
# Check fault geometry
#------------------------------------------------------------------------------
if nb_fault>=2:
    # # Get estimate for Lnuc
    a = 0.005
    b = 0.007
    Lb = problem.material_loading.mu*problem.faults[1].Dc[0]/(np.abs(problem.faults[1].sigmaN[0])*b);
            
    # # Compute Lnuc from Viesca [2016]
    Lnuc = 2*Lb / (np.pi*(1-a/b)**2)
    print('Lb (smaller fault) = {0}'.format(Lb))
    print('Lnuc (smaller fault) = {0}'.format(Lnuc))




# Plot geometry
plt.figure(1)
plt.plot(problem.faults[0].node[0,:],problem.faults[0].node[1,:],'k',linewidth=2)
for fault_id in range(1,nb_fault):
    plt.plot(problem.faults[fault_id].node[0,:],problem.faults[fault_id].node[1,:])

# plt.xlim(-50,50)
# plt.ylim(-50,50)

plt.xlabel('Distance (m)')
plt.ylabel('Distance (m)')


plt.savefig('geometry_LSBB.pdf')


ax = plt.gca()
# ax.set_aspect('equal', adjustable='box')

plt.show()


#------------------------------------------------------------------------------
# Write input file
#------------------------------------------------------------------------------
path_problem = "/Users/pierre/Dropbox/OnGoingResearch/Fluid_and_earthquakes/Softwares/DemystiFicatioN_release/problems/{0}/".format(problem.hyperparameter.simulation_name)

problem.write(path_problem)


 





 