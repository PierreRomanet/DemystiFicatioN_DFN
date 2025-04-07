module solver


use variables
use omp_lib
use FVM
use ode_bs
use hdf5


implicit none
private

! Solver
real(kind=dp), dimension(:), allocatable    :: yscal         ! Tmp array for solver (scaling)

! Saving
real(kind=dp), dimension(:,:,:), allocatable           :: P_save            ! Tmp array for solver (scaling)
real(kind=dp), dimension(:), allocatable           :: time_save            ! Tmp array for solver (scaling)
integer                                            :: it_count_save        ! Number of calls to bsstep

! Time
real(kind=dp)                                   :: t_try, t_did, t_next
integer											:: time_id


!   Public 
public :: solve, calculate_derivative
!
!
!
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Main routine that is using an adaptive-time-step-solver 
!
subroutine solve()
!
! Allocate arrays
call allocate_arrays()
!
! Initialisation
it_count_save = 0
t_try = 1.0_dp
time = 0.0_dp
!
! Enter loop 
do time_id=1,max_it
	!
	! Print tome information
	print*,'time',time
	
	! save data
	call save_data()
	!
	! Calculate the derivative
    call calculate_derivative(len_P,P1D, dPdt1D)
    !
    ! Make the scaling
     yscal = tol_solver
     yscal = yscal + dabs(P1D)+dabs(dPdt1D*t_try)
    !
    !   Call the solver
    call bsstep(P1D,dPdt1D,len_P,time,t_try,tol_solver,yscal,t_did,t_next) 
    
    
     
    print*,'---------------------------------------------'
    !
    ! Set new time for next step
    t_try = t_next
	!
	! Save variable
    if(it_count_save==freq_writing_file) then
        call write_data()
        it_count_save = 0
    end if     
    !
end do
! Deallocate arrays.
call deallocate_arrays()

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Allocate all the arrays
!
subroutine allocate_arrays()
allocate(time_save(freq_writing_file))
allocate(P_save(freq_writing_file,nbx,nby))
allocate(yscal(len_P))
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Deallocate all the arrays
!
subroutine deallocate_arrays()
deallocate(time_save)
deallocate(P_save)
deallocate(yscal)
end subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   save data 
!
subroutine save_data()   
   if (modulo(time_id-1,stride_time) == 0) then 

        it_count_save = it_count_save + 1
        P_save(it_count_save,:,:) = P(1:nbx,1:nby)
        time_save(it_count_save) = time

    end if
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Write data
!
subroutine write_data()
character(len=30)   :: num2str
integer :: error
integer(HSIZE_T) :: data_dims_array(2), data_dims_time(1)
integer(HID_T)   :: file_id, dspace_id_P,dspace_id_time, dset_id_P,dset_id_time

if (it_count_save.gt.0) then

	! Initialize Fortran interface
	CALL h5open_f(error)   
	!
	! Create a new file
	write(num2str,'(I7.7)') time_id
	!CALL h5fcreate_f(trim(path_problem), H5F_ACC_TRUNC_F, file_id, error)
	CALL h5fcreate_f('solution', H5F_ACC_TRUNC_F, file_id, error)
	!
	! Specify ranks and dims
	
	data_dims_array(1) = freq_writing_file
	data_dims_array(2) = nbx
	data_dims_time = freq_writing_file
	!
	! Create dataspace (the dataset is next) "dspace_id" is returned
	CALL h5screate_simple_f(3, data_dims_array, dspace_id_P, error)
	CALL h5screate_simple_f(1, data_dims_time, dspace_id_time, error)

	
	! Create dataset with default properties "dset_id" is returned
	CALL h5dcreate_f(file_id, 'P', H5T_NATIVE_DOUBLE, dspace_id_P, dset_id_P, error)
	CALL h5dcreate_f(file_id, 'time', H5T_NATIVE_DOUBLE, dspace_id_time, dset_id_time, error)

                     
	! Write dataset 
	CALL h5dwrite_f(dset_id_P, H5T_NATIVE_DOUBLE, P_save(1:it_count_save,1:nbx,1:nby), data_dims_array, error)
	CALL h5dwrite_f(dset_id_time, H5T_NATIVE_DOUBLE, time_save(1:it_count_save), data_dims_time, error)
                      
	! Close access to dataset 1
	CALL h5dclose_f(dset_id_P, error)
	CALL h5dclose_f(dset_id_time, error)

	! Close access to data space 1
	CALL h5sclose_f(dspace_id_P, error)
	CALL h5sclose_f(dspace_id_time, error)
	
	! Close access to file 1
	CALL h5fclose_f(file_id, error)
	
	! Close fortran interface
	call h5close_f(error)
	
    
end if
end subroutine



end module