# Imports
import numpy as np
import h5py
import os


#------------------------------------------------------------------------------------------------------------------------------------------------------------
# FAULT PROPERTIES
#------------------------------------------------------------------------------------------------------------------------------------------------------------
class Fault:
    # Function that is called when creating an object
    def __init__(self):
        # Define the number of elements of the fault
        self._nb_element = 1
        self._nb_source = 0
        
    # Define the property for nb_element
    @property
    def nb_source(self):
        return self._nb_source
    @nb_source.setter
    def nb_source(self,value):
        # Check if the given value is ok
        if isinstance(value, int):
            if value>=1:
                # Update all the other array if the new value is changing 
                self._nb_source = value
                self.Q = np.zeros(value)
                self.t_injection_beg = np.zeros(value)
                self.t_injection_end = np.zeros(value)
                self.index_injection = np.zeros(value,dtype=int)
                    
        else:
            raise ValueError("The given value must be an integer")
        
    
    
    # Define the property for nb_element
    @property
    def nb_element(self):
        return self._nb_element
    @nb_element.setter
    def nb_element(self,value):
        # Check if the given value is ok
        if isinstance(value, int) or isinstance(value,np.int64):
            if value>=1:
                # Update all the other array if the new value is changing 
                if self._nb_element != value:
                    # Update the value
                    self._nb_element = value
                    self._nb_node = value + 1
                    
                    ## Reinitialise all the other array
                    # Friction
                    self.a = np.zeros(value)
                    self.b = np.zeros(value)
                    self.Dc = np.zeros(value)
                    self.f0 = np.zeros(value)
                    self.V0 = np.zeros(value)
                    self.node = np.zeros((2,value+1))
                    
                    # Stresses
                    self.sigmaN = np.zeros(value)
                    
                    # Initial parameters
                    self.V = np.zeros(value)
                    self.theta = np.zeros(value)
                    self.P = np.zeros(value)
                    
                    # Normal and shear loading
                    self.normal_loading_dot = np.zeros(value)
                    self.shear_loading_dot = np.zeros(value)
                    
                    # Permeability
                    self.permeability = np.ones(value)
                    
                    # Permeability variation
                    self.Lk = np.zeros(value)
                    self.Tk = np.zeros(value)
                    self.kmin = np.zeros(value)
                    self.kmax = np.zeros(value)
                    self.Snk = np.zeros(value)
                    
                    
            else:
                raise ValueError("The number of elements must be greater than 0")
        else:
            raise ValueError("The given value must be an integer")
    
    
    
    
    # Define the node
    @property
    def nb_node(self):
        return self.nb_element+1
    @nb_node.setter
    def nb_node(self,value):
        self.nb_element = value - 1
        self._nb_node = value
    @property
    def node(self):
        return self._node
    @node.setter
    def node(self,value):
        self._node = value

    # Friction parameters
    @property
    def a(self):
        return self._a
    @a.setter
    def a(self,value):
        self._a = value
    @property
    def b(self):
        return self._b
    @b.setter
    def b(self,value):
        self._b = value
    @property
    def Dc(self):
        return self._Dc
    @Dc.setter
    def Dc(self,value):
        self._Dc = value
    @property
    def f0(self):
        return self._f0
    @f0.setter
    def f0(self,value):
        self._f0 = value
    @property
    def V0(self):
        return self._V0
    @V0.setter
    def V0(self,value):
        self._V0 = value
    
    
    # Normal Traction on the fault
    @property
    def sigmaN(self):
        return self._sigmaN
    @sigmaN.setter
    def sigmaN(self,value):
        self._sigmaN = value
        
        
    # V and theta initial on the fault
    @property
    def V(self):
        return self._V
    @V.setter
    def V(self,value):
        self._V= value
    @property
    def theta(self):
        return self._theta
    @theta.setter
    def theta(self,value):
        self._theta= value
    @property
    def P(self):
        return self._P
    @P.setter
    def P(self,value):
        self._P= value
        
    # V and theta initial on the fault
    @property
    def normal_loading_dot(self):
        return self._normal_loading_dot
    @normal_loading_dot.setter
    def normal_loading_dot(self,value):
        self._normal_loading_dot= value
    @property
    def shear_loading_dot(self):
        return self._shear_loading_dot
    @shear_loading_dot.setter
    def shear_loading_dot(self,value):
        self._shear_loading_dot= value
        
     
      # g1.create_dataset('density', data=self.fluid.density)
        # g1.create_dataset('dyn_viscosity', data=self.fluid.dyn_viscosity)
        # g1.create_dataset('compressibility', data=self.fluid.compressibility)
        # g1.create_dataset('ini_pressure', data=self.fluid.P0)

        
        # # # Rock
        # g1 = hf.create_group('rock')
        # g1.create_dataset('compressibility', data=self.rock.compressibility)
        # g1.create_dataset('porosity', data=self.rock.porosity)
        # g1.create_dataset('permeability', data=self.rock.permeability)   
     
    # Fluid and rock properties
    @property
    def permeability(self):
        return self._permeability
    @permeability.setter
    def permeability(self,value):
        self._permeability= value

         
    @property
    def Lk(self):
        return self._Lk
    @Lk.setter
    def Lk(self,value):
        self._Lk= value
             
    @property
    def Tk(self):
        return self._Tk
    @Tk.setter
    def Tk(self,value):
        self._Tk= value
        
    @property
    def kmin(self):
        return self._kmin
    @kmin.setter
    def kmin(self,value):
        self._kmin= value
                     
    @property
    def kmax(self):
        return self._kmax
    @kmax.setter
    def kmax(self,value):
        self._kmax= value
                         
    @property
    def Snk(self):
        return self._Snk
    @Snk.setter
    def Snk(self,value):
        self._Snk= value
     
        
#------------------------------------------------------------------------------------------------------------------------------------------------------------
# HYPERPARAMETERS
#------------------------------------------------------------------------------------------------------------------------------------------------------------
class Hyperparameter:
    # Function that is called when creating an object
    def __init__(self):
        # For all the simulations
        self.simulation_name = 'test'
        self.fracture_mode = 'modeII'
        self.static_kernel = 'direct'
        self.friction_law = 'RateStateAging'
        self.tol_solver = 1e-4
        self.freq_writing_file = 10       
        self.max_it = 10
        self.stride_time = 1
        self.permeability_coupling = 0

        
        # Specific to Fast Multipole Method
        self.iprec = 4
        
        # Number of faults
        self.nb_fault = 1
        
        
        
        
        
        
#------------------------------------------------------------------------------------------------------------------------------------------------------------
# MATERIAL AND LOADING
#------------------------------------------------------------------------------------------------------------------------------------------------------------
class Material_Loading:
    # Function that is called when creating an object
    def __init__(self):
        # Material properties
        self.mu = 30e9;
        self.cp = 5000.0;
        self.cs = 3000.0;
        
        # Uniform stress rate loading
        self.Sigma31_dot = 0.0
        self.Sigma32_dot = 0.0
        self.Sigma11_dot = 0.0
        self.Sigma12_dot = 0.0
        self.Sigma22_dot = 0.0
        
        # Uniform traction rate loading
        self.normal_loading_dot = 0.0
        self.shear_loading_dot = 0.001
        
        # Fluid diffuction
        self.rock_comp = 1e-9
        self.fluid_comp = 1e-8
        self.fluid_density = 1000.0
        self.porosity = 0.1
        self.dyn_viscosity = 1e-3
#------------------------------------------------------------------------------------------------------------------------------------------------------------
# SIMULATION PARAMETERS
#------------------------------------------------------------------------------------------------------------------------------------------------------------
class DemystiFicatioN:
    # Function that is called when creating an object
    def __init__(self,nb_flt=1):
        self.hyperparameter = Hyperparameter()
        self.hyperparameter.nb_fault = nb_flt
        self.faults = [Fault() for i in range(0,self.hyperparameter.nb_fault)]
        self.material_loading = Material_Loading()
        
        
    @property
    def nb_fault(self):
        return self.hyperparameter.nb_fault
    @nb_fault.setter
    def nb_fault(self,value):
        self.hyperparameter.nb_fault = value
        self.faults = [Fault() for i in range(0,self._nb_fault)]
   
        

    def initiate_values(self):

        # Initiate the volume of the element
        self.grid.surface_el = 1.0



#------------------------------------------------------------------------------------------------------------------------------------------------------------
# SIMULATION PARAMETERS
    def write(self,path_problem):
        
        
        #python program to check if a directory exists
        if not os.path.exists(path_problem):
            # Create a new directory because it does not exist
            os.makedirs(path_problem)
        
        
        # Create the HDF5 file
        hf = h5py.File('{0}{1}'.format(path_problem,'config'), 'w')

        ## Create subgroups
        # Hyperparameters
        g1 = hf.create_group('hyperparameter')
        g1.create_dataset('tol_solver', data=self.hyperparameter.tol_solver)
        g1.create_dataset('freq_writing_file', data=self.hyperparameter.freq_writing_file)
        g1.create_dataset('max_it', data=self.hyperparameter.max_it)
        g1.create_dataset('stride_time', data=self.hyperparameter.stride_time)
        g1.create_dataset('iprec', data=self.hyperparameter.iprec)
        g1.create_dataset('nb_fault', data=self.hyperparameter.nb_fault)
        g1.create_dataset('simulation_name', data=self.hyperparameter.simulation_name)
        g1.create_dataset('fracture_mode', data=self.hyperparameter.fracture_mode)
        g1.create_dataset('static_kernel', data=self.hyperparameter.static_kernel)
        g1.create_dataset('friction_law', data=self.hyperparameter.friction_law)
        g1.create_dataset('permeability_coupling', data=self.hyperparameter.permeability_coupling )

        
        # Create subgroup for each fault
        for i in  range(0,self.nb_fault):
            # Make name for the fault
            name_subgroup = 'fault'+str(i+1)
            
            # Create a group for each fault
            g1 = hf.create_group(name_subgroup)
            g1.create_dataset('nb_element', data=self.faults[i].nb_element)
            g1.create_dataset('a', data=self.faults[i].a)
            g1.create_dataset('b', data=self.faults[i].b)
            g1.create_dataset('Dc', data=self.faults[i].Dc)
            g1.create_dataset('f0', data=self.faults[i].f0)
            g1.create_dataset('V0', data=self.faults[i].V0)
            g1.create_dataset('sigmaN', data=self.faults[i].sigmaN)
            g1.create_dataset('V', data=self.faults[i].V)
            g1.create_dataset('theta', data=self.faults[i].theta)
            g1.create_dataset('P', data=self.faults[i].P)

            # The transpose is needed because of C vs fortran convention row major vs column major
            g1.create_dataset('node', data=np.transpose(self.faults[i].node),shape=(self.faults[i].nb_node,2))
            g1.create_dataset('normal_loading_dot', data=self.faults[i].normal_loading_dot)
            g1.create_dataset('shear_loading_dot', data=self.faults[i].shear_loading_dot )
            
            # Permeability
            g1.create_dataset('permeability', data=self.faults[i].permeability )
            
            # Permeability change
            g1.create_dataset('Lk', data=self.faults[i].Lk )
            g1.create_dataset('Tk', data=self.faults[i].Tk )
            g1.create_dataset('kmin', data=self.faults[i].kmin)
            g1.create_dataset('kmax', data=self.faults[i].kmax)
            g1.create_dataset('Snk', data=self.faults[i].Snk)

            # Source
            g1.create_dataset('nb_source', data=self.faults[i].nb_source)
            if self.faults[i].nb_source>=1:
                g1.create_dataset('Q', data=self.faults[i].Q)
                g1.create_dataset('t_injection_beg', data=self.faults[i].t_injection_beg)
                g1.create_dataset('t_injection_end', data=self.faults[i].t_injection_end)
                g1.create_dataset('index_injection', data=self.faults[i].index_injection)



        
        # Material and loading
        g1 = hf.create_group('material_loading')
        g1.create_dataset('mu', data=self.material_loading.mu)
        g1.create_dataset('cp', data=self.material_loading.cp)
        g1.create_dataset('cs', data=self.material_loading.cs)
        g1.create_dataset('Sigma31_dot', data=self.material_loading.Sigma31_dot)
        g1.create_dataset('Sigma32_dot', data=self.material_loading.Sigma32_dot)
        g1.create_dataset('Sigma11_dot', data=self.material_loading.Sigma11_dot)
        g1.create_dataset('Sigma12_dot', data=self.material_loading.Sigma12_dot)
        g1.create_dataset('Sigma22_dot', data=self.material_loading.Sigma22_dot)

        # Fluid diffution
        g1.create_dataset('rock_comp', data=self.material_loading.rock_comp)
        g1.create_dataset('fluid_comp', data=self.material_loading.fluid_comp)
        g1.create_dataset('dyn_viscosity', data=self.material_loading.dyn_viscosity)
        g1.create_dataset('porosity', data=self.material_loading.porosity)
        g1.create_dataset('fluid_density', data=self.material_loading.fluid_density)


        # Close the file
        hf.close() 
        
        
        
        

#------------------------------------------------------------------------------------------------------------------------------------------------------------
# SIMULATION PARAMETERS
    # def remesh(self,path_problem):
        
        
        
        
        
        
        
        
        
#------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read parameters
    def read(self,path_problem,problem_name):

        # Open the file to read
        hf = h5py.File('{0}{1}'.format(path_problem,problem_name), 'r')

        # Read from the file
        # Hyperparameters
        self.hyperparameter.freq_writing_file = hf.get('hyperparameter/freq_writing_file')[()]
        self.hyperparameter.tol_solver = hf.get('hyperparameter/tol_solver')[()]
        self.hyperparameter.max_it = hf.get('hyperparameter/max_it')[()]
        self.hyperparameter.stride_time = hf.get('hyperparameter/stride_time')[()]
         
        # Grid
        self.grid.nbx = hf.get('grid/nbx')[()]
        self.grid.dx = hf.get('grid/dx')[()]

        # Fluid
        self.fluid.density = hf.get('fluid/density')[()]
        self.fluid.dyn_viscosity = hf.get('fluid/dyn_viscosity')[()]
        self.fluid.compressibility = hf.get('fluid/compressibility')[()]
        self.fluid.P0 = hf.get('fluid/ini_pressure')[()]

        
        # # Rock
        self.rock.compressibility = hf.get('rock/compressibility')[()]
        self.rock.porosity = hf.get('rock/porosity')[()]
        self.rock.permeability = hf.get('rock/permeability')[()]

        # # Source
        self.source.nb_source = hf.get('source/nb_source')[()]
        self.source.index_position_x = hf.get('source/index_position_x')[()]
        self.source.t_beg = hf.get('source/t_beg')[()]
        self.source.t_end = hf.get('source/t_end')[()]
        self.source.injection_rate = hf.get('source/injection_rate')[()]

        
        # Close the file
        hf.close()
        
        
        
def load_data(data_path):
    # Open the file to read
    hf = h5py.File(data_path, 'r')

    # Read from the file
    # Hyperparameters
    time = hf.get('time')[()]
    data = hf.get('P')[()]
    data = np.asfortranarray(data)
    
    
    return time,data
        
        
        
        
        