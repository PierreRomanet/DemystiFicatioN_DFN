!
!
!
!
!
module initialisation
use variables
implicit none
!
!
!
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Main subroutine 
!
subroutine initialise()
integer :: Status
!
!   Allocate all the arrays 
call allocate_ini()
!
!   Create Initial values 
call ini_value()
!

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   allocate arrays 
!   
subroutine allocate_ini()
!
!   ALLOCATE COMMON ARRAYS
!
allocate(permeability_x(0:nbx))
!
! Allocate P and dPdt
len_P = (nbx+2)
allocate(P1D(len_P))
allocate(dPdt1D(len_P))
!
! Make 1D pointer
P(0:nbx+1) => P1d(1:len_P) 
dPdt(0:nbx+1) => dPdt1d(1:len_P) 

P = 0._dp
P(1:nbx) = P0
dPdt = 0._dp


end subroutine
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   initialise the value for permeability
!
subroutine ini_value()
! Initialise the permeability to 0
permeability_x = 0._dp

! Make the geometric average
permeability_x(1:nbx-1) = sqrt(permeability(1:nbx-1)*permeability(2:nbx))


end subroutine


end module 