module FVM

use variables

implicit none
private



!   Public 
public :: calculate_derivative
!
!
!
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Calculate the derivative of P using Finite Volume Method
!
subroutine calculate_derivative(nb,P1D_,dPdt1D_)
real(kind=dp), dimension(nb), target :: P1d_, dPdt1d_
real(kind=dp), dimension(:,:), pointer :: P_
real(kind=dp), dimension(:,:), pointer :: dPdt_
integer :: source_id
integer :: nb

! Assign pointer
P_(0:nbx+1) => P1D_(1:len_P)
dPdt_(0:nbx+1) => dPdt1D_(1:len_P)
 
! Compute derivative 
dPdt_ = 0.0_dp   

dPdt_(1:nbx) = 1._dp/(porosity*dyn_viscosity*(fluid_compressibility+rock_compressibility)) * &
            ((permeability_x(1:nbx)*(P(2:nbx+1)-P(1:nbx))  &
             -permeability_x(0:nbx-1)*(P(1:nbx)-P(0:nbx-1))/dx**2) 

! Add Sources of fluid
! do source_id=1,nb_source
! 	!If injection rate is chosen 
!     if(constant_pressure(source_id)==0) then
!         dPdt(index_position_x+1,index_position_y+1,index_position_z+1) = 1!&
!        ! dPdt(index_position_x+1,index_position_y+1,index_position_z+1) &
!        !     +injection_rate/(porosity*dyn_viscosity*(rock_compressibility+fluid_compressibility))
!     endif
! 	! If constant pressure is chosen
!     if(constant_pressure(source_id)==1) then
!         dPdt(index_position_x+1,index_position_y+1,index_position_z+1) = 0.0
!     endif
! 
! end do
end subroutine
end module