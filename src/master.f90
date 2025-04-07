!
!
!
!
!
!
program earthquake_cycle
use configuration
use initialisation
use solver
! use omp_lib
! real(kind=dp)           :: t_beg,t_end




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%% Reading data %%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print*, "Reading configuration file" 
call config()
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%% Initialisation %%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!print*, "Initialising arrays"
call initialise()
!
!call show_variables()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%% Running cycles %%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! t_beg = omp_get_wtime()
call solve()
! t_end = omp_get_wtime()
! open(unit=1,file='./problems/'//trim(simulation_name)//'/stats',status='UNKNOWN')
!     write(1,*) max_it, nb_fault, nb_element, t_end-t_beg, static_kernel, &
!                tol_solver, iprec, tol_interp
! close(1)
! print*,'Total Time',t_end-t_beg



end program
