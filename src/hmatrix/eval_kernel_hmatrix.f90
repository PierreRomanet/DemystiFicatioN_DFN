module eval_kernel_hmatrix
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
FUNCTION eval_kernel_ModeIII(ind_row,ind_col,Ynode,cnod_row,nbnod_row,id_row,id_col)
use variables, only: node_left, node_right, tangent, element
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    !
! FUNCTION NAME: eval_kernel                                         !
!                                                                    !
! INPUTS:  Xnode    : array storing the coordinates of the point x   !
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
 
  REAL(kind=8),DIMENSION(2)           :: Ynode
  INTEGER                             :: ind_row,ind_col,nbnod_row
  INTEGER                             :: id_row,id_col
  REAL(kind=8),DIMENSION(2,nbnod_row) :: cnod_row
  !Output
  COMPLEX(kind=8):: eval_kernel_ModeIII

  !Local Variables
  REAL(kind=8):: r1,r2,I1,I2
 
 
  r1 = dsqrt((node_left(1,id_col+ind_col-1)-element(1,ind_row+id_row-1))**2     &
     +      (node_left(2,id_col+ind_col-1)-element(2,ind_row+id_row-1))**2)
  r2 = dsqrt((node_right(1,id_col+ind_col-1)-element(1,ind_row+id_row-1))**2      &
     +      (node_right(2,id_col+ind_col-1)-element(2,ind_row+id_row-1))**2)
  ! 
  !     Calcul of I
  I1 = tangent(1,ind_row+id_row-1)/r1**2 * (element(1,ind_row+id_row-1)-node_left(1,id_col+ind_col-1))  &
     + tangent(2,ind_row+id_row-1)/r1**2 * (element(2,ind_row+id_row-1)-node_left(2,id_col+ind_col-1)) 
  I2 = tangent(1,ind_row+id_row-1)/r2**2 * (element(1,ind_row+id_row-1)-node_right(1,id_col+ind_col-1)) &
     + tangent(2,ind_row+id_row-1)/r2**2 * (element(2,ind_row+id_row-1)-node_right(2,id_col+ind_col-1)) 

  eval_kernel_ModeIII = CMPLX(-I1+I2,0.0d0,8)
  
  
  RETURN

END FUNCTION eval_kernel_ModeIII
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
FUNCTION eval_kernel_ModeII_T(ind_row,ind_col,Ynode,cnod_row,nbnod_row,id_row,id_col)
use variables, only: node_left, node_right, normal, element, normal_left, normal_right
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    !
! FUNCTION NAME: eval_kernel                                         !
!                                                                    !
! INPUTS:  Xnode    : array storing the coordinates of the point x   !
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
 
  REAL(kind=8),DIMENSION(2)           :: Ynode
  INTEGER                             :: ind_row,ind_col,nbnod_row
  INTEGER                             :: id_row,id_col
  REAL(kind=8),DIMENSION(2,nbnod_row) :: cnod_row
  integer                             :: k,l
  !Output
  COMPLEX(kind=8):: eval_kernel_ModeII_T

  real(kind=8) :: r_left , I1_t_left,  I2_t_left,  gam_1_left,  gam_2_left
  real(kind=8) :: r_right, I1_t_right, I2_t_right, gam_1_right, gam_2_right
 
  ! Calculation of index
  k = ind_row+id_row-1
  l = id_col+ind_col-1
 
  !Calcul of r
  r_left = sqrt((element(1,k)-node_left(1,l))**2      &
         + (element(2,k)-node_left(2,l))**2)
  r_right = sqrt((element(1,k)-node_right(1,l))**2      &
         + (element(2,k)-node_right(2,l))**2)
        
! Calculation of gam_1 gam_2 (see formulation of Taku Tada 1995 Thesis)
gam_1_left  = (element(1,k)-node_left(1,l))/r_left
gam_2_left  = (element(2,k)-node_left(2,l))/r_left
gam_1_right = (element(1,k)-node_right(1,l))/r_right
gam_2_right = (element(2,k)-node_right(2,l))/r_right

        
!Calculation of I shear  
I1_t_left  = (4.d0*normal(1,k)*normal(2,k)*gam_1_left*gam_2_left + & 
             (normal(2,k)**2-normal(1,k)**2)*(gam_2_left**2-gam_1_left**2))*gam_1_left/r_left 
I2_t_left  = (4.d0*normal(1,k)*normal(2,k)*gam_1_left*gam_2_left + &
              (normal(2,k)**2-normal(1,k)**2)*(gam_2_left**2-gam_1_left**2))*gam_2_left/r_left 
I1_t_right = (4.d0*normal(1,k)*normal(2,k)*gam_1_right*gam_2_right+ &
            (normal(2,k)**2-normal(1,k)**2)*(gam_2_right**2-gam_1_right**2))*gam_1_right/r_right 
I2_t_right = (4.d0*normal(1,k)*normal(2,k)*gam_1_right*gam_2_right+ &
            (normal(2,k)**2-normal(1,k)**2)*(gam_2_right**2-gam_1_right**2))*gam_2_right/r_right

!print*,stressT(k) 
eval_kernel_ModeII_T = CMPLX( normal_left(2,l)*I1_t_left-normal_right(2,l)*I1_t_right     &
                            - normal_left(1,l)*I2_t_left+normal_right(1,l)*I2_t_right ,0.0d0,8)


  RETURN

END FUNCTION eval_kernel_ModeII_T
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
FUNCTION eval_kernel_ModeII_N(ind_row,ind_col,Ynode,cnod_row,nbnod_row,id_row,id_col)
use variables, only: node_left, node_right, normal, tangent, element, normal_left, normal_right
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    !
! FUNCTION NAME: eval_kernel                                         !
!                                                                    !
! INPUTS:  Xnode    : array storing the coordinates of the point x   !
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
 
  REAL(kind=8),DIMENSION(2)           :: Ynode
  INTEGER                             :: ind_row,ind_col,nbnod_row
  INTEGER                             :: id_row,id_col
  REAL(kind=8),DIMENSION(2,nbnod_row) :: cnod_row
  integer                             :: k,l
  !Output
  COMPLEX(kind=8):: eval_kernel_ModeII_N

 real(kind=8) :: r_left , I1_n_left,  I2_n_left, gam_1_left,  gam_2_left
 real(kind=8) :: r_right, I1_n_right, I2_n_right, gam_1_right, gam_2_right

  
  ! Calculation of index
  k = ind_row+id_row-1
  l = id_col+ind_col-1
  
 !Calcul of r
r_left = sqrt((element(1,k)-node_left(1,l))**2      &
     + (element(2,k)-node_left(2,l))**2)
r_right = sqrt((element(1,k)-node_right(1,l))**2      &
     + (element(2,k)-node_right(2,l))**2)

! Calculation of gam_1 gam_2 (see formulation of Taku Tada 1995 Thesis)
gam_1_left  = (element(1,k)-node_left(1,l))/r_left
gam_2_left  = (element(2,k)-node_left(2,l))/r_left
gam_1_right = (element(1,k)-node_right(1,l))/r_right
gam_2_right = (element(2,k)-node_right(2,l))/r_right

        
! Calculation of I normal
I1_n_left  = gam_1_left/r_left   - gam_2_left/r_left*(2.d0*normal(1,k)*normal(2,k)*(gam_2_left**2-gam_1_left**2)  &
                                 -  (normal(2,k)**2-normal(1,k)**2)*2.d0*gam_1_left*gam_2_left)
I2_n_left  = gam_2_left/r_left   + gam_1_left/r_left*(2.d0*normal(1,k)*normal(2,k)*(gam_2_left**2-gam_1_left**2)  &
                                 -  (normal(2,k)**2-normal(1,k)**2)*2.d0*gam_1_left*gam_2_left)
I1_n_right = gam_1_right/r_right - gam_2_right/r_right*(2.d0*normal(1,k)*normal(2,k)*(gam_2_right**2-gam_1_right**2) &
                                 -(normal(2,k)**2-normal(1,k)**2)*2.d0*gam_1_right*gam_2_right)
I2_n_right = gam_2_right/r_right + gam_1_right/r_right*(2.d0*normal(1,k)*normal(2,k)*(gam_2_right**2-gam_1_right**2) &
                                 -(normal(2,k)**2-normal(1,k)**2)*2.d0*gam_1_right*gam_2_right)

 
!         print*,' ' 
!         print*,'i,j',k,l                                  
!         print*,'I1l',I1_n_left
!         print*,'I1r',I1_n_right
!         print*,'I2l',I2_n_left
!         print*,'I2r',I2_n_right
!         print*,'n2l',normal_left(2,l)
!         print*,'n2r',normal_right(2,l)
 
       !print*,stressT(k) 
 
  eval_kernel_ModeII_N = CMPLX( normal_left(1,l)*I1_n_left - normal_right(1,l)*I1_n_right  &
                              + normal_left(2,l)*I2_n_left - normal_right(2,l)*I2_n_right ,0.0d0,8)
!  
  RETURN

END FUNCTION eval_kernel_ModeII_N
end module
