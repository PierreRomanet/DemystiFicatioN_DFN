!     
!
!     Direct implementation of the tension T3 for mode III rupture.
!     We used the scheme from Tada and Madariaga 2011.
!     Author: Romanet Pierre (romanet@ipgp.fr)
!
subroutine directII(nelement,node,slip,iprec,stressT,stressN,ier,mu)
use variables, only: cp,cs, normal_left, normal_right, normal, element
implicit none
! Variables
integer, intent(in) :: nelement, iprec, ier
! 
real(kind=8), intent(in), dimension(2,nelement+1) :: node
double precision, dimension(nelement) :: slip
double precision, intent(out), dimension(nelement) :: stressN
double precision, intent(out), dimension(nelement) :: stressT
!
real(kind=8) :: mu
real(kind=8) :: r_left , I1_n_left,  I2_n_left,  I1_t_left,  I2_t_left,  gam_1_left,  gam_2_left
real(kind=8) :: r_right, I1_n_right, I2_n_right, I1_t_right, I2_t_right, gam_1_right, gam_2_right
integer :: k, l
!
real(kind=8), parameter :: pi=4.d0*atan(1.d0)





!
!!$OMP PARALLEL 
!!$OMP DO PRIVATE(k,r_left,r_right,I1_t_left,I2_t_left,I1_t_right,I2_t_right,gam_1_left,gam_2_left,gam_1_right,gam_2_right)
! For each value of element on the fault =>
do k=1,nelement
     stressT(k) = 0.d0 ! Initialization
     stressN(k) = 0.d0 ! Initialization
!
!    Compute the interaction due to a node on the fault => 
     do l=1,nelement

        !Calcul of r
        r_left = sqrt((element(1,k)-node(1,l))**2      &
           + (element(2,k)-node(2,l))**2)
        r_right = sqrt((element(1,k)-node(1,l+1))**2      &
           + (element(2,k)-node(2,l+1))**2)
        !print*,'r:',r_left,r_right
        ! Calculation of gam_1 gam_2 (see formulation of Taku Tada 1995 Thesis)
        gam_1_left  = (element(1,k)-node(1,l))/r_left
        gam_2_left  = (element(2,k)-node(2,l))/r_left
        gam_1_right = (element(1,k)-node(1,l+1))/r_right
        gam_2_right = (element(2,k)-node(2,l+1))/r_right

        !print*,'gam:',gam_1_left,gam_2_left,gam_1_right,gam_2_right
        !Calculation of I shear  
        I1_t_left  = (4.d0*normal(1,k)*normal(2,k)*gam_1_left*gam_2_left + & 
                    (normal(2,k)**2-normal(1,k)**2)*(gam_2_left**2-gam_1_left**2))*gam_1_left/r_left 
        I2_t_left  = (4.d0*normal(1,k)*normal(2,k)*gam_1_left*gam_2_left + &
                    (normal(2,k)**2-normal(1,k)**2)*(gam_2_left**2-gam_1_left**2))*gam_2_left/r_left 
        I1_t_right = (4.d0*normal(1,k)*normal(2,k)*gam_1_right*gam_2_right+ &
                    (normal(2,k)**2-normal(1,k)**2)*(gam_2_right**2-gam_1_right**2))*gam_1_right/r_right 
        I2_t_right = (4.d0*normal(1,k)*normal(2,k)*gam_1_right*gam_2_right+ &
                    (normal(2,k)**2-normal(1,k)**2)*(gam_2_right**2-gam_1_right**2))*gam_2_right/r_right
        !print*,'I:',I1_t_left,I2_t_left,I1_t_right,I2_t_right
        ! Calculation of I normal
        I1_n_left  = gam_1_left/r_left   - gam_2_left/r_left*(2.d0*normal(1,k)*normal(2,k)*(gam_2_left**2-gam_1_left**2)  &
                                         -  (normal(2,k)**2-normal(1,k)**2)*2.d0*gam_1_left*gam_2_left)
        I2_n_left  = gam_2_left/r_left   + gam_1_left/r_left*(2.d0*normal(1,k)*normal(2,k)*(gam_2_left**2-gam_1_left**2)  &
                                         -  (normal(2,k)**2-normal(1,k)**2)*2.d0*gam_1_left*gam_2_left)
        I1_n_right = gam_1_right/r_right - gam_2_right/r_right*(2.d0*normal(1,k)*normal(2,k)*(gam_2_right**2-gam_1_right**2) &
                                         -(normal(2,k)**2-normal(1,k)**2)*2.d0*gam_1_right*gam_2_right)
        I2_n_right = gam_2_right/r_right + gam_1_right/r_right*(2.d0*normal(1,k)*normal(2,k)*(gam_2_right**2-gam_1_right**2) &
                                         -(normal(2,k)**2-normal(1,k)**2)*2.d0*gam_1_right*gam_2_right)
                                                      
                                         
        stressT(k) = stressT(k) + slip(l)*(normal_left(2,l)*I1_t_left-normal_right(2,l)*I1_t_right &
                                          -normal_left(1,l)*I2_t_left+normal_right(1,l)*I2_t_right)
        
        stressN(k) = stressN(k) + slip(l)*(normal_left(1,l)*I1_n_left-normal_right(1,l)*I1_n_right &
                                          +normal_left(2,l)*I2_n_left-normal_right(2,l)*I2_n_right)
                                          
                 
                                          
!          print*,'stressT(k)',k,l,stressT(k)
! 	   print*,' ' 
! 	   print*,'i,j',k,l                                  
! 	   print*,'I1l',I1_n_left
! 	   print*,'I1r',I1_n_right
! 	   print*,'I2l',I2_n_left
! 	   print*,'I2r',I2_n_right
! 	   print*,'n2l',normal_left(2,l)
! 	   print*,'n2r',normal_right(2,l)
! !       
!        
!        print*,I2_n_left
!        print*,'1',normal_left(2,l)*I2_n_left-normal_right(2,l)*I2_n_right
!        print*,'2',normal_left(1,l)*I1_n_left-normal_right(1,l)*I1_n_right
     end do 
end do

!!$OMP END DO
!!$OMP END PARALLEL 
    
stressN = mu/(pi)*(1.d0-cs**2/cp**2)*stressN
stressT = mu/(pi)*(1.d0-cs**2/cp**2)*stressT



end subroutine directII
