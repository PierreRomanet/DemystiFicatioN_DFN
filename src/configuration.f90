!
!
!
!
!
!
!
module configuration
use variables
use hdf5
!
implicit none
!
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Main routine that read config, data, and then allocate everything
!
subroutine config()
!
! Get path to config file and geometry file
call get_location_config()
!
! Read config file
call read_config()

!
! Show parameter
call show_variables()


end subroutine 




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read location of config.in and geometry.in from command line argument
!
subroutine get_location_config()
integer         :: narg
!
!   Get number of argument
narg = command_argument_count()
!
!   If there is one argument, get it
if(narg==1) then
        call get_command_argument(1,config_file)
else
    print*,'Error: please specify a directory containing config'
    stop
end if
!
!   Build path to config_file
config_file = trim(config_file)//'config'
!
end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!          
!          Read initial parameter to allocate arrays
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Read the main configuration file
!
subroutine read_config()
integer(HID_T)  :: file, dset
integer :: hdferr
integer(HSIZE_T), DIMENSION(1) :: dims 
integer(HSIZE_T), DIMENSION(2) :: dims_node
integer :: fault_id,node_count

!
! Open the HDF5 file
CALL h5open_f(hdferr) ! Initialise fortran interface
CALL h5fopen_f( "./"//trim(config_file), H5F_ACC_RDWR_F , file, hdferr) ! Open the file 
!





dims(1) = 1
CALL h5dopen_f (file, "hyperparameter/tol_solver", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,tol_solver,dims, hdferr)
CALL h5dopen_f (file, "hyperparameter/max_it", dset, hdferr)
CALL h5dread_f(dset, H5T_STD_I64LE,max_it,dims, hdferr)
CALL h5dopen_f(file, "hyperparameter/stride_time", dset, hdferr)
CALL h5dread_f(dset, H5T_STD_I64LE ,stride_time,dims, hdferr)
CALL h5dopen_f(file, "hyperparameter/iprec", dset, hdferr)
CALL h5dread_f(dset, H5T_STD_I64LE,iprec,dims, hdferr)
CALL h5dopen_f(file, "hyperparameter/nb_fault", dset, hdferr)
CALL h5dread_f(dset, H5T_STD_I64LE,nb_fault,dims, hdferr)
CALL h5dopen_f(file, "hyperparameter/freq_writing_file", dset, hdferr)
CALL h5dread_f(dset, H5T_STD_I64LE,freq_writing_file,dims, hdferr)
CALL h5dopen_f(file, "hyperparameter/permeability_coupling", dset, hdferr)
CALL h5dread_f(dset, H5T_STD_I64LE,permeability_coupling,dims, hdferr)
CALL h5dclose_f(dset , hdferr)


! Read all the hyperparameters (there is a big bug is reading text before other parameters - change the variable)
call  read_vl_string(file, "hyperparameter/fracture_mode", fracture_mode )
call  read_vl_string(file, "hyperparameter/static_kernel",  static_kernel )
call  read_vl_string(file, "hyperparameter/simulation_name",  simulation_name )
call  read_vl_string(file, "hyperparameter/friction_law",  fric_law )



!#########################################################################################
! Find number of nodes and number of elements
allocate(faults(nb_fault))

! For each fault 
nb_node = 0
do fault_id=1,nb_fault

    ! Read number of element
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/nb_element", dset, hdferr)
    CALL h5dread_f(dset, H5T_STD_I64LE,faults(fault_id)%nb_element,dims, hdferr)
    
    ! Update the number of nodes
    nb_node = nb_node + faults(fault_id)%nb_element + 1
    
    ! Read number of injection source
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/nb_source", dset, hdferr)
    CALL h5dread_f(dset, H5T_STD_I64LE,faults(fault_id)%nb_source,dims, hdferr)
    
    
enddo

! Update the total number of elements (including fake element with 0 slip in between faults)
nb_element = nb_node-1


! 
!#########################################################################################



! Read material and loading
CALL h5dopen_f (file, "material_loading/Sigma11_dot", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,Sigma11_dot,dims, hdferr)
CALL h5dopen_f (file, "material_loading/Sigma12_dot", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,Sigma12_dot,dims, hdferr)
CALL h5dopen_f (file, "material_loading/Sigma22_dot", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,Sigma22_dot,dims, hdferr)
CALL h5dopen_f (file, "material_loading/Sigma31_dot", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,Sigma31_dot,dims, hdferr)
CALL h5dopen_f (file, "material_loading/Sigma32_dot", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,Sigma32_dot,dims, hdferr)
CALL h5dopen_f (file, "material_loading/cp", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,cp,dims, hdferr)
CALL h5dopen_f (file, "material_loading/cs", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,cs,dims, hdferr)
CALL h5dopen_f (file, "material_loading/mu", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,mu,dims, hdferr)
CALL h5dopen_f (file, "material_loading/porosity", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,porosity,dims, hdferr)
CALL h5dopen_f (file, "material_loading/fluid_comp", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,fluid_comp,dims, hdferr)
CALL h5dopen_f (file, "material_loading/rock_comp", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,rock_comp,dims, hdferr)
CALL h5dopen_f (file, "material_loading/fluid_density", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,fluid_density,dims, hdferr)
CALL h5dopen_f (file, "material_loading/dyn_viscosity", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,dyn_viscosity,dims, hdferr)

! allocate all the array 
! Geometry
allocate(element(2,nb_element))
allocate(node(2,nb_node))
allocate(tangent(2,nb_element))
allocate(normal(2,nb_element))
allocate(ds(nb_element))
allocate(mapping_fault(nb_fault,2))


! Problem 
! allocate(Y(4*nb_element))
! allocate(Y_dot(4*nb_element))
length_Y=3*nb_element
allocate(Y(3*nb_element)) ! Without permeability change
allocate(Y_dot(3*nb_element))


allocate(slip(nb_element))

! Rate and state parameters
allocate(a(nb_element))
allocate(b(nb_element))
allocate(Dc(nb_element))
allocate(f0(nb_element))
allocate(V0(nb_element))
! a, V0, to avoid division by 0 in-between elements, it is initialised to 1
V0 = 1._dp
a = 1._dp
b = 1._dp
Dc = 1._dp
f0 = 1._dp

! Stresses
allocate(normal_loading_dot(nb_element))
allocate(shear_loading_dot(nb_element))
allocate(SigmaN_dot(nb_element))           ! Derivative of the normal traction on the faults
allocate(SigmaT_dot(nb_element))           ! Derivative of the tangential traction on the faults



! Allocate permeability
allocate(permeability(nb_element))           ! Derivative of the tangential traction on the faults
allocate(permeability_star(nb_element))           ! Derivative of the tangential traction on the faults
allocate(permeability_x(1:nb_element+1))

! Variable permeability
allocate(Lk(nb_element)) 
allocate(Tk(nb_element))  
allocate(Snk(nb_element))  
allocate(kmin(nb_element))  
allocate(kmax(nb_element))  
! Snk, Tk, and Lk appear as inverse (1/Snk), to avoid division by 0 in-between elements, it is initialised to 1
Snk = 1._dp
Lk = 1._dp
Tk = 1._dp

! allocate
allocate(P(nb_element))
allocate(P_dot(nb_element))

! Initialise pointers
! Pointer for faults 
V          => Y(1:nb_element)
theta      => Y(nb_element+1:2*nb_element)
SigmaN     => Y(2*nb_element+1:3*nb_element)
!permeability_star => Y(3*nb_element+1:4*nb_element)


node_count = 0
do fault_id=1,nb_fault

! Mapping of the fault
    mapping_fault(fault_id,1) = node_count+1 
    mapping_fault(fault_id,2) = node_count+faults(fault_id)%nb_element
!   Nodes
    faults(fault_id)%node       => node(1 : 2 , mapping_fault(fault_id,1) : mapping_fault(fault_id,2)+1)
!   Slip, V and theta
    faults(fault_id)%slip       => slip(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%V          => Y(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%theta      => Y(nb_element+mapping_fault(fault_id,1) : nb_element+mapping_fault(fault_id,2))
!   Rate and state parameters
    faults(fault_id)%a          => a(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%b          => b(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%Dc         => Dc(mapping_fault(fault_id,1): mapping_fault(fault_id,2))
    faults(fault_id)%V0         => V0(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%f0         => f0(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
!   Tractions
    faults(fault_id)%SigmaN     => SigmaN(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%SigmaN_dot => SigmaN_dot(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    !faults(fault_id)%SigmaT     => SigmaT(node_count+1 : node_count+faults(fault_id)%nb_element)
    faults(fault_id)%SigmaT_dot => SigmaT_dot(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
!   Tangent and normal
    faults(fault_id)%tangent    => tangent(1:2,mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%normal     => normal(1:2,mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
!   ds
    faults(fault_id)%ds         => ds(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
! Pore Pressure
    faults(fault_id)%P          => P(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%P_dot      => P_dot(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
! Background loading
    faults(fault_id)%normal_loading_dot  => normal_loading_dot(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%shear_loading_dot   => shear_loading_dot(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
! Permeability
    faults(fault_id)%permeability  => permeability(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))

    faults(fault_id)%permeability_star  => permeability_star(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%permeability_x(0:faults(fault_id)%nb_element)  => &
                                      permeability_x(node_count : mapping_fault(fault_id,2))

! Variable Permeability
    faults(fault_id)%Tk  => Tk(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%Lk  => Lk(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%Snk  => Snk(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%kmin  => kmin(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))
    faults(fault_id)%kmax  => kmax(mapping_fault(fault_id,1) : mapping_fault(fault_id,2))

    
!   Node count
    node_count = node_count + faults(fault_id)%nb_element + 1 
    
!   Allocate sources
    if (faults(fault_id)%nb_source>=1) then
        allocate(faults(fault_id)%Q(faults(fault_id)%nb_source))
        allocate(faults(fault_id)%t_injection_beg(faults(fault_id)%nb_source))
        allocate(faults(fault_id)%t_injection_end(faults(fault_id)%nb_source))
        allocate(faults(fault_id)%index_injection(faults(fault_id)%nb_source))
    endif
end do

dims(1) = nb_element
! For each fault read initial parameters
do fault_id=1,nb_fault
    ! Read all the values of friction
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/a", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%a,dims, hdferr)  ! Read the data using the default properties.
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/b", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%b,dims, hdferr)
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/Dc", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%Dc,dims, hdferr)
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/f0", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%f0,dims, hdferr)
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/V0", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%V0,dims, hdferr)
        
    !
    ! Read all the values of loading 
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/normal_loading_dot", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%normal_loading_dot,dims, hdferr)
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/shear_loading_dot", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%shear_loading_dot,dims, hdferr)
    ! 
    

    ! Read all the values of initial_parameters
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/V", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%V,dims, hdferr)
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/theta", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%theta,dims, hdferr)
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/sigmaN", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%sigmaN,dims, hdferr)    
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/P", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%P,dims, hdferr)    
    !
    ! Read permeability
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/permeability", dset, hdferr) ! Only part where it is assign to permeability_star
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%permeability,dims, hdferr)  
    !
    ! Read variable permeability
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/kmax", dset, hdferr) ! 
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%kmax,dims, hdferr)  
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/kmin", dset, hdferr) ! 
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%kmin,dims, hdferr) 
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/Tk", dset, hdferr) ! 
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%Tk,dims, hdferr) 
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/Lk", dset, hdferr) ! 
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%Lk,dims, hdferr) 
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/Snk", dset, hdferr) ! 
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%Snk,dims, hdferr)  
    !
    ! Read nodes
    dims_node(1) = 2
    dims_node(2) = nb_node
    CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/node", dset, hdferr)
    CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%node,dims_node, hdferr)  
   
    !
    ! Read source parameters
    if(faults(fault_id)%nb_source>=1) then
        dims = faults(fault_id)%nb_source
        CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/Q", dset, hdferr)
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%Q,dims, hdferr)    
        CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/t_injection_beg", dset, hdferr)
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%t_injection_beg,dims, hdferr)    
        CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/t_injection_end", dset, hdferr)
        CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,faults(fault_id)%t_injection_end,dims, hdferr)    
        CALL h5dopen_f (file, "fault"//trim(int2str(fault_id))//"/index_injection", dset, hdferr)
        CALL h5dread_f(dset, H5T_STD_I64LE,faults(fault_id)%index_injection,dims, hdferr)    
    endif
  
enddo


! Close and release resources.
CALL h5dclose_f(dset , hdferr)
CALL h5fclose_f(file , hdferr)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!          
!          Read variable length string 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine read_vl_string(file, dataset, string)
INTEGER(HSIZE_T), PARAMETER       :: dim0      = 1
INTEGER(HSIZE_T), PARAMETER       :: sdim      = 200
INTEGER(HID_T)                    :: file, filetype, space, dset ! Handles
INTEGER                           :: hdferr
INTEGER(HSIZE_T), DIMENSION(1:1)  :: dims = (/dim0/)
INTEGER(HSIZE_T), DIMENSION(1:2)  :: maxdims
CHARACTER(LEN=sdim), DIMENSION(1) :: rdata ! Read buffer
character(len=*)                  :: string, dataset
INTEGER(HSIZE_T), DIMENSION(2)    :: data_dims = (/sdim,dim0/)
INTEGER(SIZE_T), DIMENSION(1)     :: str_len = 200



! Open the dataset
CALL h5dopen_f (file, trim(dataset), dset, hdferr) ! open the dataset
!
! Get the datatype.
CALL H5Dget_type_f(dset, filetype, hdferr) 
!
! Get dataspace and allocate memory for read buffer.
CALL H5Dget_space_f(dset, space, hdferr)
CALL H5Sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
!
! Read the data.
CALL h5dread_vl_f(dset, filetype, rdata, data_dims, str_len, hdferr, space)
string = trim(rdata(1))


  CALL h5dclose_f(dset , hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL H5Tclose_f(filetype, hdferr)
  
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!          
!   Convert an integer to string
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Read the main configuration file
!
character(len=20) function int2str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (int2str, *) k
    int2str = adjustl(int2str)
end function int2str

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!       Print every parameter at screen
!
subroutine show_variables()
integer                             :: k,l
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%% Modeling %%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*
print*,"You are modeling ",nb_fault,"fault(s) with a ",trim(fric_law)," friction law in ",trim(fracture_mode)
print*
print*,"You chose to calculate the kernel with: ",trim(static_kernel) 
print*
print*,"%%%%%%%%%%%%%% Material properties %%%%%%%%%%%%%"
print*
print*,"mu = ", mu
print*,"Cs = ", cs
print*,"Cp = ", cp
print*
print*,"%%%%%%%%%%%%%% Fluid diffusion %%%%%%%%%%%%%%%"
print*
print*,"Fluid compressibility = ", fluid_comp
print*,"Fluid density = ", fluid_density
print*,"Dynamic viscosity = ", dyn_viscosity
print*,"Rock compressibility = ", rock_comp
print*,"Porosity = ", porosity
print*
do k=1,nb_fault
    print*,'faults(k)%nb_source',faults(k)%nb_source
    do l=1,faults(k)%nb_source
        print*
        print*,'Fault ',k, ', Source ',l
        print*,"index_injection = ", faults(k)%index_injection(l)
        print*,"t_injection_beg = ", faults(k)%t_injection_beg(l)
        print*,"t_injection_end = ", faults(k)%t_injection_end(l)
        print*,"Q = ", faults(k)%Q(l)
    end do
end do
! print*,"%%%%%%%%%%%%%%% Plate loading %%%%%%%%%%%%%%%%%%"
! print*
! if(fracture_mode == 'modeIII') then
!     print*,"Sigma31_dot = ", Sigma31_dot
!     print*,"Sigma32_dot = ", Sigma32_dot
! 
! elseif(fracture_mode == 'modeII') then
!     print*,"Sigma11_dot = ", Sigma11_dot
!     print*,"Sigma22_dot = ", Sigma22_dot
!     print*,"Sigma12_dot = ", Sigma12_dot
! 
! end if
! print*
print*,"%%%%%%%%%%%%%%% Hyperparameters %%%%%%%%%%%%%%"
print*
print*,"Tolerance :", tol_solver
print*,"Tolerance for FMM :",iprec
print*,"Frequence writing file :", freq_writing_file
print*,"Stride time :", stride_time

print*
print*,"%%%%%%%%%%%%%%% Faults Parameters %%%%%%%%%%%%%%%"
! do k=1,nb_fault
!     print*,'Fault ',k
!     print*
!     print*,'Nb of elements = ', faults(k)%nb_element
!     print*
!     print*
! end do
print*,'Please press ENTER'
!read(*,*)

end subroutine
end module
