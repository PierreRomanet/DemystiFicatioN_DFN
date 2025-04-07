module read_input

use variables
use hdf5

implicit none
!
!
!
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Main routine that read config, data, and then allocate everything
!
subroutine config()
!
! Get problem path
call get_problem_path()
!
! Read size of the arrays
call read_hyperparameters()
!
! Allocate list
call allocate_arrays()
!
! Read 
call read_ini()
!
! Show all the variables
call show_variables()
end subroutine 


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Get location of the problem_path
!
subroutine get_problem_path()
integer         :: narg,stat1
!
!   Get number of argument
narg = command_argument_count()
!
!   If there is one argument, get it
if(narg==1) then
        call get_command_argument(1,path_problem)
else
    print*,'Error: please specify to the diffusion problem'
    stop
end if

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Allocate all the arrays
!
subroutine allocate_arrays()
integer :: i

! Wells
allocate(index_position_x(nb_source))
allocate(t_beg(nb_source))
allocate(t_end(nb_source))
allocate(injection_rate(nb_source))
allocate(injection_pressure(nb_source))
allocate(constant_pressure(nb_source))
! For each well 
! do i=1,nb_source
! 	allocate(wells(i)%t_beg(wells(i)%nb_injection_phase))
! 	allocate(wells(i)%t_end(wells(i)%nb_injection_phase))
! 	allocate(wells(i)%injection_rate(wells(i)%nb_injection_phase))
! 	allocate(wells(i)%injection_pressure(wells(i)%nb_injection_phase))
! 	allocate(wells(i)%constant_pressure(wells(i)%nb_injection_phase))
! enddo

! Pressure
allocate(P0(nbx,nby))

! Rock
allocate(permeability(nbx,nby)) 

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Read all the given nodes
!
subroutine read_hyperparameters()
integer(HID_T)  :: file, dset
integer :: hdferr
integer(HSIZE_T), DIMENSION(1) :: dims 

! Calculate dims
dims(1) = 1
	
! Open the file
CALL h5open_f(hdferr) ! Initialise fortran interface
CALL h5fopen_f( trim(path_problem), H5F_ACC_RDWR_F , file, hdferr) ! Open the file 
print*,path_problem
! Read the data
CALL h5dopen_f (file, "/source/nb_source", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,nb_source,dims, hdferr)
CALL h5dopen_f (file, "/grid/nbx", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_INTEGER,nbx,dims, hdferr)
CALL h5dopen_f (file, "/grid/nby", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_INTEGER,nby,dims, hdferr)

! Close and release resources.
CALL h5dclose_f(dset , hdferr)
CALL h5fclose_f(file , hdferr)
!
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Read all the given nodes
!
subroutine read_ini()
integer(HID_T)  :: file, dset
integer :: hdferr
integer(HSIZE_T), DIMENSION(2) :: dims
integer(HSIZE_T), DIMENSION(1) :: dim

! 
! Calculate dims
dims(1) = nbx
dims(2) = nby
dim=1

! 
CALL h5open_f(hdferr) ! Initialise fortran interface
CALL h5fopen_f( trim(path_problem), H5F_ACC_RDWR_F , file, hdferr) ! Open the file 
! 
! Hyperparameters
! Read all the values concerning the hyperparameters
CALL h5dopen_f (file, "hyperparameter/tol_solver", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,tol_solver,dim, hdferr)  ! Read the data using the default properties.
CALL h5dopen_f (file, "hyperparameter/freq_writing_file", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_INTEGER,freq_writing_file,dim, hdferr)
CALL h5dopen_f (file, "hyperparameter/stride_time", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_INTEGER,stride_time,dim, hdferr)
CALL h5dopen_f (file, "hyperparameter/max_it", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_INTEGER,max_it,dim, hdferr)
!
! Read all the values concerning the fluid
CALL h5dopen_f (file, "fluid/compressibility", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,fluid_compressibility,dim, hdferr)  ! Read the data using the default properties.
CALL h5dopen_f (file, "fluid/density", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,density,dim, hdferr)
CALL h5dopen_f (file, "fluid/dyn_viscosity", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,dyn_viscosity,dim, hdferr)
CALL h5dopen_f (file, "fluid/ini_pressure", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,P0,dims, hdferr)
!
! Read all the values concerning the grid
CALL h5dopen_f (file, "grid/dx", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,dx,dim, hdferr)
CALL h5dopen_f (file, "grid/dy", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,dy,dim, hdferr)
! 
! Read all the values concerning the rock
CALL h5dopen_f (file, "rock/compressibility", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,rock_compressibility,dim, hdferr)
CALL h5dopen_f (file, "rock/porosity", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,porosity,dim, hdferr)
CALL h5dopen_f (file, "rock/permeability", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,permeability,dims, hdferr)
! 
! Read all the values concerning the source
dim = nb_source
CALL h5dopen_f (file, "source/index_position_x", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_INTEGER,index_position_x,dim, hdferr)
CALL h5dopen_f (file, "source/index_position_y", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_INTEGER,index_position_y,dim, hdferr)
CALL h5dopen_f (file, "source/constant_pressure", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,t_beg,dim, hdferr)
CALL h5dopen_f (file, "source/constant_pressure", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,t_end,dim, hdferr)
CALL h5dopen_f (file, "source/injection_pressure", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,injection_pressure,dim, hdferr)
CALL h5dopen_f (file, "source/injection_rate", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_DOUBLE,injection_rate,dim, hdferr)
CALL h5dopen_f (file, "source/constant_pressure", dset, hdferr)
CALL h5dread_f(dset, H5T_NATIVE_INTEGER,constant_pressure,dim, hdferr)
! ! 
! Close and release resources.
CALL h5dclose_f(dset , hdferr)
CALL h5fclose_f(file , hdferr)

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Print all the parameters
!
subroutine show_variables()
use omp_lib
integer :: nb_thread_mkl,nb_thread_omp,i

!

print*,' //////////////////////////////////////////////////'
print*,' //                                              //'
print*,' //                 2DiffUnion                   //'
print*,' //               Pierre Romanet                 //'
print*,' //                                              //'
print*,' //////////////////////////////////////////////////'
print*,''
! nb_thread_omp = omp_get_max_threads()
! nb_thread_mkl = mkl_get_max_threads()
! print*,'Thread MKL = ', nb_thread_mkl
! print*,'Thread OMP = ', nb_thread_omp
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%% Hyperparameters %%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*
print*,"max_it = ", max_it
print*,"freq_writing_file = ", freq_writing_file
print*,"stride_time = ", stride_time
print*,"tol_solver = ", tol_solver
print*
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%%   Grid   %%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*
print*,"nbx, nby = ", nbx
print*,"dx, dy = ", dx
print*
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%%   Fluid  %%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*
print*,"Compressibility = ",fluid_compressibility
print*,"Dynamic viscosity = ",dyn_viscosity
print*,"Density = ",density
print*,"Initial pressure = ",P0(1,1)
print*
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%%   Rock   %%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*
print*,"Compressibility = ",rock_compressibility
print*,"Porosity = ",porosity
print*,"Permeability = ",permeability(1,1)
print*
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%%  Source  %%%%%%%%%%%%%%%%%%%%"
print*,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
print*,'nb_source = ',nb_source
do i=1,nb_source
print*,'------------------------------'
print*,'Source ',i
print*,"index_position_x = ", index_position_x(i)
print*,"injection_pressure = ", injection_pressure(i)
print*,"injection_rate = ", injection_rate(i)
print*,"constant_pressure = ", constant_pressure(i)
enddo

end subroutine




end module