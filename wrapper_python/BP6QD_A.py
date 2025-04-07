from DemystiFicatioN import *
import numpy as np
from matplotlib import pyplot as plt
# Variable 
nb_fault = 1
fault_length = 40000.0

# Create new problem
problem = DemystiFicatioN(nb_fault)


# Geometry
problem.faults[0].nb_element = 4000

#------------------------------------------------------------------------------
# Hyperparameter
#------------------------------------------------------------------------------
# Simulation name
problem.hyperparameter.simulation_name = "BP6QD_A_{0}".format(problem.faults[0].nb_element)
# Frature mode
problem.hyperparameter.fracture_mode = "modeIII"
# Way of calculating the static kernel
problem.hyperparameter.static_kernel = "hmatrix"
# Type of rate and state friction law
problem.hyperparameter.friction_law = "RateStateAging"
# Tolerance for the solver
problem.hyperparameter.tol_solver = 1e-6
# Frequence of writting output file 
problem.hyperparameter.freq_writing_file = 1000
# Number of time step after which the simulation stops     
problem.hyperparameter.max_it = 30000
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
problem.material_loading.cp = 1.0;
problem.material_loading.cs = 3.464e3;
 
# Uniform stress rate loading
problem.material_loading.Sigma31_dot = 0.0
problem.material_loading.Sigma32_dot = 0.0
problem.material_loading.Sigma11_dot = 0.0
problem.material_loading.Sigma12_dot = 0.0
problem.material_loading.Sigma22_dot = 0.0
 
# Uniform traction rate loading
problem.material_loading.normal_loading_dot = 0.0
problem.material_loading.shear_loading_dot = 0.00

problem.material_loading.rock_comp = 0.1e-8
problem.material_loading.fluid_comp = 0.9e-8
problem.material_loading.fluid_density = 1000.0
problem.material_loading.porosity = 0.1
problem.material_loading.dyn_viscosity = 1e-3
 
 
#------------------------------------------------------------------------------ 
# Fault properties
#------------------------------------------------------------------------------

 
# Friction
problem.faults[0].a = 0.007*np.ones(problem.faults[0].nb_element)
problem.faults[0].b = 0.005*np.ones(problem.faults[0].nb_element)
problem.faults[0].Dc = 0.004*np.ones(problem.faults[0].nb_element)
problem.faults[0].f0 = 0.6*np.ones(problem.faults[0].nb_element)
problem.faults[0].V0 = 1e-6*np.ones(problem.faults[0].nb_element)

# Stresses
problem.faults[0].sigmaN[:] = -50.0e6;
problem.faults[0].normal_loading_dot[:] = 0.0;
problem.faults[0].shear_loading_dot[:] = 0.00;

# Initial parameters (V and theta)
sigmaT = 29.2e6;
problem.faults[0].V = 1e-12*np.ones(problem.faults[0].nb_element)
problem.faults[0].theta[:] = problem.faults[0].Dc[0]/problem.faults[0].V0[0]*np.exp(problem.faults[0].a[0]/problem.faults[0].b[0]*np.log(2.0*problem.faults[0].V0[0]/problem.faults[0].V[0]*np.sinh(sigmaT/(problem.faults[0].a[0]*np.abs(problem.faults[0].sigmaN[0]))))-problem.faults[0].f0[0]/problem.faults[0].b[0]);
problem.faults[0].theta[:] = np.log(problem.faults[0].V0[0]*problem.faults[0].theta[0]/problem.faults[0].Dc[0]);
problem.faults[0].P = 0.0*np.ones(problem.faults[0].nb_element);


# Permeability
problem.faults[0].permeability = 1e-13*np.ones(problem.faults[0].nb_element)


 # Variable permeability
 # Permeability variation
problem.faults[0].Lk = 0.0005*np.ones(problem.faults[0].nb_element)
problem.faults[0].Tk = 1e8*np.ones(problem.faults[0].nb_element)
problem.faults[0].kmin = 8e-16*np.ones(problem.faults[0].nb_element)
problem.faults[0].kmax = 1e-13**np.ones(problem.faults[0].nb_element)
problem.faults[0].Snk = 5.0e6*np.ones(problem.faults[0].nb_element)

# Injection
problem.faults[0].nb_source = 1
problem.faults[0].Q =  1.25e-6
problem.faults[0].t_injection_beg = 0.0
problem.faults[0].t_injection_end = 100.0*86400.0
#problem.faults[0].index_injection = np.floor(problem.faults[0].nb_element/2).astype(int)





#------------------------------------------------------------------------------
# Create fault geometry 
#------------------------------------------------------------------------------
ds = fault_length/problem.faults[0].nb_element
problem.faults[0].node[0,:] = np.linspace(-fault_length/2+ds/2,fault_length/2+ds/2 ,problem.faults[0].nb_element+1)
problem.faults[0].node[1,:] = np.zeros(problem.faults[0].nb_element+1)

element = (problem.faults[0].node[:,1:]+problem.faults[0].node[:,0:-1])/2


# Find index of injection
idx = np.argmin(np.abs(element[0,:]))
# Be careful that fortran indices and python indices are not the same
problem.faults[0].index_injection = idx+1


#------------------------------------------------------------------------------
# Write input file
#------------------------------------------------------------------------------
path_problem = "/Users/pierre/Dropbox/OnGoingResearch/Fluid_and_earthquakes/Softwares/DemystiFicatioN_release/problems/{0}/".format(problem.hyperparameter.simulation_name)

problem.write(path_problem)













 