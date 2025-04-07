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


nb_fault = 10
problem = DemystiFicatioN(nb_fault)

ds = 0.1


#------------------------------------------------------------------------------
# Hyperparameter
#------------------------------------------------------------------------------
# Simulation name
problem.hyperparameter.simulation_name = "network_rate_strengthening"
# Frature mode
problem.hyperparameter.fracture_mode = "modeII"
# Way of calculating the static kernel
problem.hyperparameter.static_kernel = "hmatrix"
# Type of rate and state friction law
problem.hyperparameter.friction_law = "RateStateAging"
# Tolerance for the solver
problem.hyperparameter.tol_solver = 1e-6
# Frequence of writting output file 
problem.hyperparameter.freq_writing_file = 1000
# Number of time step after which the simulation stops     
problem.hyperparameter.max_it = 100000
# Save every stride_time time step
problem.hyperparameter.stride_time = 1
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
problem.material_loading.shear_loading_dot = 0.001

problem.material_loading.rock_comp = 0.1e-8
problem.material_loading.fluid_comp = 0.9e-8
problem.material_loading.fluid_density = 1000.0
problem.material_loading.porosity = 0.1
problem.material_loading.dyn_viscosity = 1e-3


#------------------------------------------------------------------------------
# Create fault geometry 1
#------------------------------------------------------------------------------
angle = np.pi/180*0
L_fault = 10.0

# Geometry
problem.faults[0].nb_element = np.floor(L_fault/ds).astype(int);

# Create nodes of the fault 
problem.faults[0].node[0,:] = np.linspace(-L_fault/2,L_fault/2,problem.faults[0].nb_element+1)
problem.faults[0].node[1,:] = np.zeros(np.size(problem.faults[0].node[0,:]))


#------------------------------------------------------------------------------
# Create fault geometry 2
#------------------------------------------------------------------------------
angle = np.pi/180*60
L_fault = 10.0
offset_x = -3.

# Geometry
problem.faults[1].nb_element = np.floor(L_fault/ds).astype(int);

# Create nodes of the fault 
x = np.linspace(-L_fault/2,L_fault/2,problem.faults[1].nb_element+1)
y = np.zeros(np.size(problem.faults[1].node[0,:]))
problem.faults[1].node[0,:] = np.cos(angle)*(x-offset_x) -np.sin(angle)*y-4.0
problem.faults[1].node[1,:] = np.sin(angle)*(x-offset_x) + np.cos(angle)*y


#------------------------------------------------------------------------------
# Create fault geometry 3
#------------------------------------------------------------------------------
angle = np.pi/180*60
L_fault = 9.0
offset_x = -3.0


# Geometry
problem.faults[2].nb_element = np.floor(L_fault/ds).astype(int);

# Create nodes of the fault 
# Create nodes of the fault 
x = np.linspace(-L_fault/2,L_fault/2,problem.faults[2].nb_element+1)
y = np.zeros(np.size(problem.faults[2].node[0,:]))
problem.faults[2].node[0,:] = np.cos(angle)*(x-offset_x) -np.sin(angle)*y+4.0
problem.faults[2].node[1,:] = np.sin(angle)*(x-offset_x) + np.cos(angle)*y


#------------------------------------------------------------------------------
# Create fault geometry 4
#------------------------------------------------------------------------------
angle = np.pi/180*0
L_fault = 13.0

height1 = 5.0
dy1 = ds*np.round(height1/(ds*np.sin(np.pi/180*60)))*np.sin(np.pi/180*60)
dx = np.round(np.round(height1/(ds*np.sin(np.pi/180*60)))*np.cos(np.pi/180*60)/ds)*ds-np.round(height1/(ds*np.sin(np.pi/180*60)))*np.cos(np.pi/180*60)




# Geometry
problem.faults[3].nb_element = np.floor(L_fault/ds).astype(int);

# Create nodes of the fault 
problem.faults[3].node[0,:] = np.linspace(-L_fault/2,L_fault/2,problem.faults[3].nb_element+1)+dx+3
problem.faults[3].node[1,:] = np.zeros(np.size(problem.faults[3].node[0,:]))+dy1

#------------------------------------------------------------------------------
# Create fault geometry 5
#------------------------------------------------------------------------------
angle = np.pi/180*60
L_fault = 12.0
offset_x = -4.0


# Geometry
problem.faults[4].nb_element = np.floor(L_fault/ds).astype(int);

# Create nodes of the fault 
# Create nodes of the fault 
x = np.linspace(-L_fault/2,L_fault/2,problem.faults[4].nb_element+1)
y = np.zeros(np.size(problem.faults[4].node[0,:]))
problem.faults[4].node[0,:] = np.cos(angle)*(x-offset_x) -np.sin(angle)*y-2.5
problem.faults[4].node[1,:] = np.sin(angle)*(x-offset_x) + np.cos(angle)*y+dy1

#------------------------------------------------------------------------------
# Create fault geometry 6
#------------------------------------------------------------------------------
angle = np.pi/180*60
L_fault = 7.0
offset_x = -2.0


# Geometry
problem.faults[5].nb_element = np.floor(L_fault/ds).astype(int);

# Create nodes of the fault 
# Create nodes of the fault 
x = np.linspace(-L_fault/2,L_fault/2,problem.faults[5].nb_element+1)
y = np.zeros(np.size(problem.faults[5].node[0,:]))
problem.faults[5].node[0,:] = np.cos(angle)*(x-offset_x) -np.sin(angle)*y+1
problem.faults[5].node[1,:] = np.sin(angle)*(x-offset_x) + np.cos(angle)*y+dy1

#------------------------------------------------------------------------------
# Create fault geometry 7
#------------------------------------------------------------------------------
angle = np.pi/180*60
L_fault = 12.0
offset_x = -3.0


# Geometry
problem.faults[6].nb_element = np.floor(L_fault/ds).astype(int);

# Create nodes of the fault 
# Create nodes of the fault 
x = np.linspace(-L_fault/2,L_fault/2,problem.faults[6].nb_element+1)
y = np.zeros(np.size(problem.faults[6].node[0,:]))
problem.faults[6].node[0,:] = np.cos(angle)*(x-offset_x) -np.sin(angle)*y+4
problem.faults[6].node[1,:] = np.sin(angle)*(x-offset_x) + np.cos(angle)*y+dy1

#------------------------------------------------------------------------------
# Create fault geometry 8
#------------------------------------------------------------------------------
angle = np.pi/180*0
L_fault = 8.0

height1 = 8.0
dy = ds*np.round(height1/(ds*np.sin(np.pi/180*60)))*np.sin(np.pi/180*60)
dx = np.round(np.round(height1/(ds*np.sin(np.pi/180*60)))*np.cos(np.pi/180*60)/ds)*ds-np.round(height1/(ds*np.sin(np.pi/180*60)))*np.cos(np.pi/180*60)

# Geometry
problem.faults[7].nb_element = np.floor(L_fault/ds).astype(int);

# Create nodes of the fault 
problem.faults[7].node[0,:] = np.linspace(-L_fault/2,L_fault/2,problem.faults[7].nb_element+1)+L_fault/2+dx
problem.faults[7].node[1,:] = np.zeros(np.size(problem.faults[7].node[0,:]))+dy


#------------------------------------------------------------------------------
# Create fault geometry 10
#------------------------------------------------------------------------------
angle = np.pi/180*60
L_fault = 10.0
offset_x = -4.0


# Geometry
problem.faults[9].nb_element = np.floor(L_fault/ds).astype(int);

# Create nodes of the fault 
# Create nodes of the fault 
x = np.linspace(-L_fault/2,L_fault/2,problem.faults[9].nb_element+1)
y = np.zeros(np.size(problem.faults[9].node[0,:]))
problem.faults[9].node[0,:] = np.cos(angle)*(x-offset_x) -np.sin(angle)*y+1.5
problem.faults[9].node[1,:] = np.sin(angle)*(x-offset_x) + np.cos(angle)*y+dy

#------------------------------------------------------------------------------
# Create fault geometry 9
#------------------------------------------------------------------------------



angle = np.pi/180*0
L_fault = 11

height1 = 12.0
dy = ds*np.round(height1/(ds*np.sin(np.pi/180*60)))*np.sin(np.pi/180*60)
dx = np.round(np.round(height1/(ds*np.sin(np.pi/180*60)))*np.cos(np.pi/180*60)/ds)*ds-np.round(height1/(ds*np.sin(np.pi/180*60)))*np.cos(np.pi/180*60)



# Geometry
problem.faults[8].nb_element = np.floor(L_fault/ds).astype(int);

# Create nodes of the fault 
problem.faults[8].node[0,:] = np.linspace(-L_fault/2,L_fault/2,problem.faults[8].nb_element+1)+4+ds/2
problem.faults[8].node[1,:] = np.zeros(np.size(problem.faults[8].node[0,:]))+dy





#------------------------------------------------------------------------------ 
# Generate same properties for all faults
#------------------------------------------------------------------------------


for fault_id in range(nb_fault):
    
   
 
    # Friction
    problem.faults[fault_id].a = 0.005*np.ones(problem.faults[fault_id].nb_element)
    problem.faults[fault_id].b = 0.004*np.ones(problem.faults[fault_id].nb_element)
    problem.faults[fault_id].Dc = 0.000001*np.ones(problem.faults[fault_id].nb_element)
    problem.faults[fault_id].f0 = 0.6*np.ones(problem.faults[fault_id].nb_element)
    problem.faults[fault_id].V0 = 1e-6*np.ones(problem.faults[fault_id].nb_element)

    # Stresses
    problem.faults[fault_id].sigmaN[:] = -5.0e6;
    problem.faults[fault_id].normal_loading_dot[:] = 0.0;
    problem.faults[fault_id].shear_loading_dot[:] = 0.001;

    # Initial parameters (V and theta)
    problem.faults[fault_id].V = 1e-12*np.ones(problem.faults[fault_id].nb_element)
    problem.faults[fault_id].theta = np.log(problem.faults[fault_id].V0[fault_id]/problem.faults[fault_id].V[0])*np.ones(problem.faults[fault_id].nb_element);
    problem.faults[fault_id].P = 0.0*np.ones(problem.faults[fault_id].nb_element);


    # Permeability
    problem.faults[fault_id].permeability = 1e-15*np.ones(problem.faults[fault_id].nb_element)

    # No injection at this stage
    problem.faults[fault_id].nb_source = 0
    problem.faults[fault_id].Q =  1e-4
    problem.faults[fault_id].t_injection_beg = 10.0*365.25*86400.0
    problem.faults[fault_id].t_injection_end = 10.0*365.25*86400.0+2*1000
    problem.faults[fault_id].index_injection = np.floor(problem.faults[0].nb_element/2-200).astype(int)
    
    # Variable permeability
    # Permeability variation
    problem.faults[fault_id].Lk = 0.0005*np.ones(problem.faults[fault_id].nb_element)
    problem.faults[fault_id].Tk = 1e8*np.ones(problem.faults[fault_id].nb_element)
    problem.faults[fault_id].kmin = 8e-16*np.ones(problem.faults[fault_id].nb_element)
    problem.faults[fault_id].kmax = 1e-13**np.ones(problem.faults[fault_id].nb_element)
    problem.faults[fault_id].Snk = 5.0e6*np.ones(problem.faults[fault_id].nb_element)



# # Get estimate for Lnuc
a = problem.faults[0].a[0]
b = problem.faults[0].b[0]
Lb = problem.material_loading.mu*problem.faults[0].Dc[0]/(np.abs(problem.faults[0].sigmaN[0])*b);
            
# # Compute Lnuc from Viesca [2016]
Lnuc = 2*Lb / (np.pi*(1-a/b)**2)
print('Lb (large fault) = {0}'.format(Lb))
print('Lnuc (large fault) = {0}'.format(Lnuc))
 


#------------------------------------------------------------------------------ 
# Injection on fault 9
#------------------------------------------------------------------------------
# No injection at this stage
problem.faults[9].nb_source = 1
problem.faults[9].Q =  2e-06
problem.faults[9].t_injection_beg = 10.0*365.25*86400.0
problem.faults[9].t_injection_end = 10.0*365.25*86400.0+2000
problem.faults[9].index_injection = 80

 # Variable permeability
 # Permeability variation
problem.faults[9].Lk = 0.0005*np.ones(problem.faults[fault_id].nb_element)
problem.faults[9].Tk = 1e8*np.ones(problem.faults[fault_id].nb_element)
problem.faults[9].kmin = 8e-16*np.ones(problem.faults[fault_id].nb_element)
problem.faults[9].kmax = 1e-13**np.ones(problem.faults[fault_id].nb_element)
problem.faults[9].Snk = 5.0e6*np.ones(problem.faults[fault_id].nb_element)
       

#------------------------------------------------------------------------------ 
# Check fault geometry
#------------------------------------------------------------------------------
# if nb_fault>=2:
    # # # Get estimate for Lnuc
    # a = 0.005
    # b = 0.007
    # Lb = problem.material_loading.mu*problem.faults[1].Dc[0]/(np.abs(problem.faults[1].sigmaN[0])*b);
            
    # # # Compute Lnuc from Viesca [2016]
    # Lnuc = 2*Lb / (np.pi*(1-a/b)**2)
    # print('Lb (smaller fault) = {0}'.format(Lb))
    # print('Lnuc (smaller fault) = {0}'.format(Lnuc))




# Plot geometry
plt.figure(1)
for fault_id in range(0,nb_fault):
    plt.plot(problem.faults[fault_id].node[0,:],problem.faults[fault_id].node[1,:],'+')

# plt.xlim(-50,50)
# plt.ylim(-50,50)

plt.xlabel('Distance (m)')
plt.ylabel('Distance (m)')


plt.savefig('geometry_network.pdf')


ax = plt.gca()
# ax.set_aspect('equal', adjustable='box')

plt.show()


#------------------------------------------------------------------------------
# Write input file
#------------------------------------------------------------------------------
path_problem = "/Users/pierre/Dropbox/OnGoingResearch/Fluid_and_earthquakes/Softwares/DemystiFicatioN__release/problems/{0}/".format(problem.hyperparameter.simulation_name)

problem.write(path_problem)


 





 