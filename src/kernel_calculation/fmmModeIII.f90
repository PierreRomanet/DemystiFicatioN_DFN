!   
!   Ando scheme used with fast multipole method from Greengard
!
!   Author: Romanet Pierre   (romanet@ipgp.fr)
!
!
!
!
      subroutine fmmIII(nelement,node,slip,iprec,stressT,stressN,ier,mu)
!
use omp_lib
use variables, only: element, tangent
        implicit none
!        
        integer ier, iprec, nelement
        integer ifcharge, ifdipole, ifpot, ifgrad, ifhess
        integer ifpottarg, ifgradtarg, ifhesstarg        
!
        real(kind=8) :: node(2,nelement+1)
        real *8 :: dipstr
        real *8 :: dipvec(2,1)
!       
        real *8 :: pot
        real *8 :: grad(2,1)
        real *8 :: hess
!
        real *8 :: pottarg
        real *8 :: gradtarg(2,nelement)
        real *8 :: hesstarg
!
        real *8 :: slip(nelement)
        real *8 :: coeff(nelement+1) 
        
        real *8 :: mu
        real *8 :: stressN(nelement)
        real *8 :: stressT(nelement)
!
        real(kind=8), parameter :: pi=4*atan(1.d0)
        integer  :: k
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!    Compute the coefficient for Ando-scheme (Constant slip in a box)
      do k=2,nelement
      coeff(k) = mu/(2*pi)*(slip(k) - slip(k-1))
      end do
      coeff(1) = mu/(2*pi)*slip(1)
      coeff(nelement+1) = -mu/(2*pi)*slip(nelement)
!
!    Set parameters for fmm
!     Compute only the gradient at target points
      ifcharge = 1        
      ifdipole = 0
      ifpot = 0
      ifgrad = 0
      ifhess = 0
      ifpottarg = 0
      ifgradtarg = 1
      ifhesstarg = 0
!
!
!    Call fast multipole methods with no dipoles and no hessian 
      call rfmm2dparttarg(ier,iprec,             &
                       nelement+1,node,          &
                       ifcharge,coeff,           &
                       ifdipole,dipstr,dipvec,   &
                       ifpot,pot,                &
                       ifgrad,grad,              &
                       ifhess,hess,              &
                       nelement,element,         &
                       ifpottarg,pottarg,        &
                       ifgradtarg,gradtarg,      &
                       ifhesstarg,hesstarg)
!
!     Scalar product with t
      do k=1,nelement
      stressT(k) = - tangent(1,k)*gradtarg(1,k) - tangent(2,k)*gradtarg(2,k)
      end do

!
      end subroutine 
