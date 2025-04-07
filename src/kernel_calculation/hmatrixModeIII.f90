!   
!   Ando scheme used with fast multipole method from Greengard
!
!   Author: Romanet Pierre   (romanet@ipgp.fr)
!
!
!
!
subroutine hmatIII(nelement,node,slip,iprec,stressT,stressN,ier,mu)
use variables, only: dimA_T, dimB_T, dimF_T,                          & 
                     A_matrices_T, B_matrices_T, F_matrices_T,       &
                     nbclus, ind_clusters, nbblclus, blcluster_tree, &
                     info_blclustersT, mat_corr,dp
use calculation_hmat
use omp_lib
!
        implicit none
!        
integer(kind=4)                             ::nelement
real(kind=8)                                :: node(2,nelement+1)
real(kind=dp), dimension(nelement)          :: slip
complex(kind=dp), dimension(nelement)       :: slip_h
real(kind=dp)                               :: mu
real(kind=dp)                               :: stressN(nelement)
real(kind=dp)                               :: stressT(nelement)
complex(kind=dp)                            :: stress_h(nelement)
!
COMPLEX(kind=8)           :: alpha, beta
!
real(kind=8), parameter :: pi=4*atan(1.d0)
integer  :: k, ier, iprec
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!   Set parameter for Hgemv routine (matrix-vector product)
alpha = cmplx(1.0d0,0.d0,8)
beta  = cmplx(0.0d0,0.0d0,8)
!
!   Convert into H-matrix representation (reorder)
stress_h = cmplx(0.d0,0.0d0,8)
do k=1,nelement
    slip_h(k) = cmplx(slip(mat_corr(k)),0.0d0,8)
end do
!
call Hgemv_right(nbclus,ind_clusters,nbblclus,                &
                blcluster_tree,info_blclustersT,              & 
                alpha,beta,1,nelement,nelement,               &
                slip_h,stress_h,dimA_T,dimB_T,dimF_T,         &
                 A_matrices_T,B_matrices_T,F_matrices_T)
!
!    Convert into real representation (reorder)
do k=1,nelement
    slip(mat_corr(k)) = real(slip_h(k),8)
    stressT(mat_corr(k)) = mu/(2.0d0*pi)*real(stress_h(k),8)
end do

! print*,stressT
! stop
end subroutine 
