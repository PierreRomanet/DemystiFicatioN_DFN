!
!   Rate and state law with ageing law
!   This module provides variables and subroutines for rate and state law 
!
!
module RateStateAgingModeII
use variables
use fluid_diffusion, only: compute_P
use omp_lib
implicit none
private
public:: ode_RateStateAgingModeII



contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!       ODEs to supply to the solver
!
subroutine ode_RateStateAgingModeII(dt,nv,Y1,Y1_dot)
integer(kind=dp)                             :: idx
integer                               :: k, ier,nv
real(kind=dp), dimension(:),INTENT(IN)  :: Y1
real(kind=dp), dimension(:),INTENT(OUT)  :: Y1_dot
real(kind=dp), dimension(nb_element)         :: B1, C1, D1 
real(kind=dp), dimension(nb_element)         :: P_temp
real(kind=dp), intent(in)                    :: dt
integer(kind=dp) :: i
ier = 0                                                  
!
! Calculate elastic stress
call kernel_ptr(nb_element,node,Y1(1:nb_element),iprec,SigmaT_dot,SigmaN_dot,ier,mu)
!
! Calculate the pressure 
call compute_P(dt,P_temp,P_dot,ier)
if (ier==1) then 
print*,'Did not converge (inside RateState)'
stop
endif


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PERMEABILITY EVOLUTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Update permeability_star
! if (permeability_coupling == 1) then
!     Y1_dot(3*nb_element+1:4*nb_element) = -abs(Y1(1:nb_element)) &
!                                            /Lk(1:nb_element)*(Y1(3*nb_element+1:4*nb_element)-kmax(1:nb_element)) &
!                                         -1._dp/Tk(1:nb_element)*(Y1(3*nb_element+1:4*nb_element)-kmin(1:nb_element))
! else
!     Y1_dot(3*nb_element+1:4*nb_element) = 0._dp
! endif




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! NORMAL TRACTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y1_dot(2*nb_element+1:3*nb_element) =  SigmaN_dot(1:nb_element)             &
                                    +  normal_loading_dot(1:nb_element) 




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
                                  - abs(Y1(1:nb_element))/Dc(1:nb_element) ))                   & ! Check abs here 
      / (mu/(2._dp*cs)-D1(1:nb_element)*(Y1(2*nb_element+1:3*nb_element)+P_temp(1:nb_element) ) &
                                       *a(1:nb_element)/(2._dp*V0(1:nb_element)))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THETA
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y1_dot(nb_element+1:2*nb_element) = V0(1:nb_element) /Dc(1:nb_element) *exp(-Y1(nb_element+1:2*nb_element))  & 
                                  - abs(Y1(1:nb_element))/Dc(1:nb_element) 


! Make sure that Y1_dot between faults is 0
idx = 0
do k=1,nb_fault
    idx = idx + faults(k)%nb_element + 1
    if (idx<=nb_element) then
        Y1_dot(idx) = 0
        Y1_dot(idx+nb_element) = 0
        Y1_dot(idx+2*nb_element) = 0
    end if
end do

end subroutine
end module
