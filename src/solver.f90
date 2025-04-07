!
!       ODE solver using bsstep from numerical recipes
!
module solver
use variables
use fluid_diffusion, only:compute_P
!use RateState
use ode_bs
use hdf5

!
!   Variables
implicit none
private
real(kind=dp), dimension(:,:), allocatable  :: V_save, theta_save,P_save         ! Tmp arrays to save data
real(kind=dp), dimension(:,:), allocatable  :: permeability_save   ! Tmp arrays to save data
real(kind=dp), dimension(:,:), allocatable  :: tractionN_save                 ! Tmp arrays to save data
real(kind=dp), dimension(:,:), allocatable  :: tractionT_save                 ! Tmp arrays to save data
real(kind=dp), dimension(:), allocatable    ::  time_save               ! Tmp array for solver (scaling)


integer                                     :: istop                          ! Integer for stop criteria
integer                                     :: it_count, it_count_save        ! Number of calls to bsstep
integer                                     :: k                              ! Integer for for-loop  


! To REMOVE
real(kind=dp), dimension(:), allocatable     :: P_temp, yscal


real(kind=dp)                                   :: t_try, t_did, t_next, time_temp 

!
!   Public 
public :: solve
! 
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Solver routine 
!
subroutine solve()
integer(kind=dp)                                :: fault_id, source_id

!
!   Allocate temporary data
call allocate_tmp_data()
!
!   Initializing variables
istop = 0
t_try = 100._dp
it_count = 0
it_count_save = 0

! Update time_temp
time_temp = time 

!   Main loop of the solver
do while(istop/=1)                  
    it_count = it_count + 1
    !!
    
   
    if (t_try> 10._dp) then
     ! Be careful not to skip injection
     do fault_id=1,nb_fault
         ! For each source 
         do source_id=1,faults(fault_id)%nb_source
         ! If time of injection
            if((time<faults(fault_id)%t_injection_beg(source_id) )&
                 .and.(time+t_try> faults(fault_id)%t_injection_beg(source_id))) then
                        t_try = 10._dp
            endif
!             
         enddo
     enddo 
     endif


     ! Check if there is no problem of convergence for the given time step
     call compute_P(t_try,P_temp,P_dot,ier)
!     print*,'permeability'
!     print*,permeability
!     print*,'permeability_x'
! print*,permeability_x
! stop
    
     do while(ier==1)
     print*,'did not converge'     
    t_try = t_try/2._dp
     call compute_P(t_try,P_temp,P_dot,ier)
    enddo




    time_temp = 0._dp ! Initialised time to 0
    
    
   !   Do the scaling     
   call fric_law_ptr(t_try,length_Y,Y, Y_dot)   
   yscal = abs(Y)+abs(Y_dot*t_try)+ tol_solver ! KEEP THIS SCALING IS SUPER IMPORTANT

   !call bsstep(Y,Y_dot,length_Y,time_temp,t_try,tol_solver,yscal,t_did,t_next,fric_law_ptr)
   !call bsstep(Y,Y_dot,length_Y,time_temp,t_try,tol_solver,yscal,t_did,t_next)

   ! print*,'permeability',permeability/porosity/(fluid_comp+rock_comp)/dyn_viscosity
  !  print*,'ds',ds
!    stop
! !     stop
   call rkqs(Y,Y_dot,length_Y,time_temp,t_try,tol_solver,yscal,t_did,t_next,fric_law_ptr)

      
      ! Set new time for next step
      print*,'t_try',t_try
      print*,'t_next',t_next
      t_try = t_next
      if (t_try>24*60*60)then
               t_try=24*60*60
    end if
        
        
        
    !  if (it_count == 14268) then
!               print*,'Y(1:nb_element)'
!               print*,Y(2000:nb_element) 
!               print*,'-----------------------------------------------------------------------'
!               print*,'Y(2*nb_element+2001) '
!               print*,Y(2*nb_element+2001:3*nb_element) 
!               
!               
!               print*,'-----------------------------------------------------------------------'
!               print*,'Y_dot(1:nb_element)'
!               print*,Y_dot(2000:nb_element) 
!               print*,'-----------------------------------------------------------------------'
!               print*,'Y_dot(2*nb_element:3nb_element)'
!               print*,Y_dot(2*nb_element+2001:3*nb_element) 
!     endif
!         
!    if (it_count == 14774) then
!               print*,'Y(1:nb_element)'
!               print*,Y(2000:nb_element) 
!               print*,'-----------------------------------------------------------------------'
!               print*,'Y(2*nb_element+2001) '
!               print*,Y(2*nb_element+2001:3*nb_element) 
!               
!               
!               print*,'-----------------------------------------------------------------------'
!               print*,'Y_dot(1:nb_element)'
!               print*,Y_dot(2000:nb_element) 
!               print*,'-----------------------------------------------------------------------'
!               print*,'Y_dot(2*nb_element:3nb_element)'
!               print*,Y_dot(2*nb_element+2001:3*nb_element) 
!     endif
!         
!        if (it_count == 9416) then
!               print*,'Y(1:nb_element)'
!               print*,Y(2000:nb_element) 
!               print*,'-----------------------------------------------------------------------'
!               print*,'Y(2*nb_element+2001) '
!               print*,Y(2*nb_element+2001:3*nb_element) 
!               
!               
!               print*,'-----------------------------------------------------------------------'
!               print*,'Y_dot(1:nb_element)'
!               print*,Y_dot(2000:nb_element) 
!               print*,'-----------------------------------------------------------------------'
!               print*,'Y_dot(2*nb_element:3nb_element)'
!               print*,Y_dot(2*nb_element+2001:3*nb_element) 
!     endif
        
        
    ! Update pressure
    call compute_P(t_did,P,P_dot,ier)

    ! Update time
    time = time+t_did
    !
    ! Save data if it_count-1=stride_time
    call save_data()
    !
    !    Print information
    call print_inf()
    ! print*,'t_next = ',t_next/(30*24*60*60)
    !
    !   Check stop criteria 
    call check_stop()
    !
    !   Write data if it_count=freq_writing_file
    if(it_count_save==freq_writing_file) then
        call write_data()
        it_count_save = 0
    end if     
    !
    ! If value are non meaningfull, especially if traction is extensive, stop
    call check_errors()
    !
end do
!
!   Write data
call write_data()
!
!
!   Deallocate temporary data 
call deallocate_tmp_data()
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   save data 
!
subroutine save_data()   
   if (modulo(it_count-1,stride_time) == 0) then 
        it_count_save = it_count_save + 1
        V_save(it_count_save,:) = Y(1:nb_element) 
        theta_save(it_count_save,:) = Dc/V0*exp(Y(nb_element+1:2*nb_element))
        P_save(it_count_save,:)  = P(1:nb_element )
        tractionN_save(it_count_save,:) = SigmaN(1:nb_element)
        tractionT_save(it_count_save,:) = SigmaT_dot(1:nb_element)  
        time_save(it_count_save) = time 
        permeability_save(it_count_save,:) = permeability(1:nb_element)
    end if
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Write data in a file
!
subroutine write_data()       
character(len=30)   :: num2str
integer :: error
integer(HSIZE_T) :: data_dims_array(2), data_dims_time(1)
integer(HID_T)   :: file_id, dspace_id_V,dspace_id_time, dset_id_V,dset_id_time
integer(HID_T)   :: dspace_id_theta,dset_id_theta
integer(HID_T)   :: dspace_id_P,dset_id_P
integer(HID_T)   :: dspace_id_tractionTel,dset_id_tractionTel
integer(HID_T)   :: dspace_id_tractionNel, dset_id_tractionNel
integer(HID_T)   :: dspace_id_permeability, dset_id_permeability





if (it_count_save.gt.0) then

    ! Initialize Fortran interface
    CALL h5open_f(error)   
    !
    ! Create a new file
    write(num2str,'(I7.7)') it_count
    CALL h5fcreate_f('problems/'//trim(simulation_name)//'/output'//trim(num2str), H5F_ACC_TRUNC_F, file_id, error)
    !
    ! Specify ranks and dims
    
    data_dims_array(1) = freq_writing_file
    data_dims_array(2) = nb_element
    data_dims_time = freq_writing_file
    !
    ! Create dataspace (the dataset is next) "dspace_id" is returned
    CALL h5screate_simple_f(2, data_dims_array, dspace_id_V, error)
    CALL h5screate_simple_f(2, data_dims_array, dspace_id_theta, error)
    CALL h5screate_simple_f(2, data_dims_array, dspace_id_P, error)
    CALL h5screate_simple_f(2, data_dims_array, dspace_id_tractionTel, error)
    CALL h5screate_simple_f(2, data_dims_array, dspace_id_tractionNel, error)
    CALL h5screate_simple_f(2, data_dims_array, dspace_id_permeability, error)

    CALL h5screate_simple_f(1, data_dims_time, dspace_id_time, error)

    
    ! Create dataset with default properties "dset_id" is returned
    CALL h5dcreate_f(file_id, 'V', H5T_NATIVE_DOUBLE, dspace_id_V, dset_id_V, error)
    CALL h5dcreate_f(file_id, 'theta', H5T_NATIVE_DOUBLE, dspace_id_theta, dset_id_theta, error)
    CALL h5dcreate_f(file_id, 'P', H5T_NATIVE_DOUBLE, dspace_id_P, dset_id_P, error)
    CALL h5dcreate_f(file_id, 'tractionTel', H5T_NATIVE_DOUBLE, dspace_id_tractionTel, dset_id_tractionTel, error)
    CALL h5dcreate_f(file_id, 'tractionNel', H5T_NATIVE_DOUBLE, dspace_id_tractionNel, dset_id_tractionNel, error)
    CALL h5dcreate_f(file_id, 'permeability', H5T_NATIVE_DOUBLE, dspace_id_permeability, dset_id_permeability, error)
    CALL h5dcreate_f(file_id, 'time', H5T_NATIVE_DOUBLE, dspace_id_time, dset_id_time, error)

                     
    ! Write dataset 
    CALL h5dwrite_f(dset_id_V, H5T_NATIVE_DOUBLE, V_save(1:it_count_save,1:nb_element), data_dims_array, error)
    CALL h5dwrite_f(dset_id_theta, H5T_NATIVE_DOUBLE, theta_save(1:it_count_save,1:nb_element), data_dims_array, error)
    CALL h5dwrite_f(dset_id_P, H5T_NATIVE_DOUBLE, P_save(1:it_count_save,1:nb_element), data_dims_array, error)
    CALL h5dwrite_f(dset_id_tractionTel, H5T_NATIVE_DOUBLE, tractionT_save(1:it_count_save,1:nb_element), data_dims_array, error)
    CALL h5dwrite_f(dset_id_tractionNel, H5T_NATIVE_DOUBLE, tractionN_save(1:it_count_save,1:nb_element), data_dims_array, error)
    CALL h5dwrite_f(dset_id_permeability, H5T_NATIVE_DOUBLE, permeability_save(1:it_count_save,1:nb_element), data_dims_array, error)
    CALL h5dwrite_f(dset_id_time, H5T_NATIVE_DOUBLE, time_save(1:it_count_save), data_dims_time, error)
                      
    ! Close access to dataset 1
    CALL h5dclose_f(dset_id_V, error)
    CALL h5dclose_f(dset_id_theta, error)
    CALL h5dclose_f(dset_id_P, error)
    CALL h5dclose_f(dset_id_tractionTel, error)
    CALL h5dclose_f(dset_id_tractionNel, error)
    CALL h5dclose_f(dset_id_permeability, error)

    CALL h5dclose_f(dset_id_time, error)

    ! Close access to data space 1
    CALL h5sclose_f(dspace_id_V, error)
    CALL h5sclose_f(dspace_id_theta, error)
    CALL h5sclose_f(dspace_id_P, error)
    CALL h5sclose_f(dspace_id_tractionTel, error)
    CALL h5sclose_f(dspace_id_tractionNel, error)
    CALL h5sclose_f(dspace_id_permeability, error)
    CALL h5sclose_f(dspace_id_time, error)
    
    ! Close access to file 1
    CALL h5fclose_f(file_id, error)
    
    ! Close fortran interface
    call h5close_f(error) 
    
    
    
end if


end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   
!   Check stop criteria, istop=1 if it stops.
!
subroutine check_stop()
! if(stop_criteria==0) then
!     if(maxval(Y(1:nb_element))>cut_vel) istop = 1
! else if (stop_criteria==1 ) then
!     if(it_count>=max_it) then
!         istop = 1
!     endif
! else if (stop_criteria==2) then
!     if(time>=final_time) istop = 1
! end if
if(it_count>=max_it) then
         istop = 1
endif
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Allocate local variables
!
subroutine allocate_tmp_data()
allocate(yscal(4*nb_element))
allocate(time_save(freq_writing_file))
allocate(V_save(freq_writing_file,nb_element))
allocate(theta_save(freq_writing_file,nb_element))
allocate(P_save(freq_writing_file,nb_element))
allocate(tractionN_save(freq_writing_file,nb_element))
allocate(tractionT_save(freq_writing_file,nb_element))
allocate(permeability_save(freq_writing_file,nb_element))


! TO REMOVE
allocate(P_temp(nb_element))


end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Deallocate local variables
!
subroutine deallocate_tmp_data()
deallocate(yscal)
deallocate(time_save)
deallocate(V_save)
deallocate(theta_save)
deallocate(P_save)
deallocate(tractionN_save)
deallocate(permeability_save)

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Print information during running
!
subroutine print_inf()
!write(*,'(I7.7,F16.5)',advance='no') it_count,time
print*,it_count,time/86400, ' days'
print*,'t_next',t_try
 
do k=1,1!min(nb_fault,5)
   ! print*,'Fault number ', k
    !write(*,'(ES3.5)',advance='no')maxval(faults(k)%V)
    print*,'max V = ',maxval(faults(k)%V)
    print*,'min V = ',minval(faults(k)%V)
    print*,'min s = ', minval(faults(k)%sigmaN)
    print*,'max s = ', maxval(faults(k)%sigmaN)
    print*,'min theta = ', minval(faults(k)%theta)
    print*,'max theta = ', maxval(faults(k)%theta)
    print*,'max P = ',maxval(faults(k)%P(1:faults(k)%nb_element))
    print*,'min P = ',minval(faults(k)%P(1:faults(k)%nb_element))
    print*,'max k = ',maxval(faults(k)%permeability)
    print*,'min k = ',minval(faults(k)%permeability)
    print*,'max k* = ',maxval(faults(k)%permeability_star)
    print*,'min k* = ',minval(faults(k)%permeability_star)
!     print*,'max Pdot = ',maxval(faults(k)%P_dot)
!     print*,'min Pdot = ',minval(faults(k)%P_dot)
end do

 print*,'max V = ',maxval(V)
print*,'max theta = ',maxval(theta)
 

print*, ''
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Print information during running
!
subroutine check_errors()
if (any(sigmaN>0)) then
    print*,'Simulation has stoped because normal traction became extensional'
    stop
end if



end subroutine
end module
