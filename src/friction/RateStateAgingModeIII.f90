!
!   Rate and state law with ageing law
!   This module provides variables and subroutines for rate and state law 
!
!
module RateStateAgingModeIII
use variables
use fluid_diffusion, only: compute_P

use omp_lib
implicit none
private
public:: ode_RateStateAgingModeIII
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!       ODEs to supply to the solver
!
subroutine ode_RateStateAgingModeIII(dt,nv,Y1,Y1_dot)
integer :: nv
real(kind=dp), intent(in)                    ::dt
integer(kind=dp)                             :: fault_id,source_id, crossing_id, idx, id_P
integer(kind=dp)                             :: fault_id1,fault_id2, node_id1, node_id2,id_ini1,id_ini2
real(kind=dp)                                :: permeability_x_crossing_left, permeability_x_crossing_right
real(kind=dp), dimension(:),INTENT(IN)       :: Y1
real(kind=dp), dimension(:),INTENT(OUT)      :: Y1_dot
real(kind=dp), dimension(nb_element)         :: P_temp
real(kind=dp), dimension(nb_element)         :: B1, C1, D1 ! coefficient for regularized version 
!real(kind=dp)           :: t_beg,t_end
!character(len=20) :: int2str
ier = 0                                                    
!
! Calculate kernel
call kernel_ptr(nb_element,node,Y1(1:nb_element),iprec,SigmaT_dot,SigmaN_dot,ier,mu)


! Calculate the pressure 
call compute_P(dt,P_temp,P_dot,ier)



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PERMEABILITY EVOLUTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Update permeability_star
! if (permeability_coupling == 1) then
!     Y1_dot(3*nb_element+1:4*nb_element) = -Y1(1:nb_element) &
!                                            /Lk(1:nb_element)*(Y1(3*nb_element+1:4*nb_element)-kmax(1:nb_element)) &
!                                         -1._dp/Tk(1:nb_element)*(Y1(3*nb_element+1:4*nb_element)-kmin(1:nb_element))
! else
!     Y1_dot(3*nb_element+1:4*nb_element) = 0._dp
! endif
! 
! 
! if(isnan(sum(Y1_dot(3*nb_element+1:4*nb_element) ))) then
! print*,'Permeability has NaN'
! stop
! endif




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! NORMAL TRACTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y1_dot(2*nb_element+1:3*nb_element) = normal_loading_dot(1:nb_element)



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! VELOCITY
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate B and C
B1(1:nb_element) = Y1(1:nb_element)/(2._dp*V0(1:nb_element))
C1(1:nb_element) = exp((f0(1:nb_element)+b(1:nb_element)*Y1(nb_element+1:2*nb_element) )/ &
                   a(1:nb_element))
D1(1:nb_element) = 1._dp / C1(1:nb_element)**2 +B1(1:nb_element)**2
D1(1:nb_element) = 1._dp / sqrt(D1(1:nb_element))

Y1_dot(1:nb_element) = (shear_loading_dot(1:nb_element) + SigmaT_dot(1:nb_element)             &
      + (Y1_dot(2*nb_element+1:3*nb_element)+P_dot(1:nb_element))                           & 
      * a(1:nb_element)*asinh(B1(1:nb_element)*C1(1:nb_element))                                   &
      + B1(1:nb_element)*D1(1:nb_element)                                                          &
      * (Y1(2*nb_element+1:3*nb_element)+P_temp(1:nb_element))*b(1:nb_element)* &
      ( V0(1:nb_element) /Dc(1:nb_element) *exp(-Y1(nb_element+1:2*nb_element))  & 
                                  - abs(Y1(1:nb_element))/Dc(1:nb_element) ))                   &
      / (mu/(2._dp*cs)-D1(1:nb_element)*(Y1(2*nb_element+1:3*nb_element)+P_temp(1:nb_element) ) &
                                       *a(1:nb_element)/(2._dp*V0(1:nb_element)))




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THETA
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y1_dot(nb_element+1:2*nb_element) = V0(1:nb_element) /Dc(1:nb_element) *exp(-Y1(nb_element+1:2*nb_element))  & 
                                  - abs(Y1(1:nb_element))/Dc(1:nb_element) 





!
!   Make sure that Y1_dot between faults is 0
idx = 0
do fault_id=1,nb_fault
    idx = idx + faults(fault_id)%nb_element + 1
    if (idx<=nb_element) then
        Y1_dot(idx) = 0
        Y1_dot(idx+nb_element) = 0
        Y1_dot(idx+2*nb_element) = 0
    end if
end do


! print*,Y1(1:2*nb_element)
! print*,Y1_dot(1:2*nb_element)
! stop


end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module
