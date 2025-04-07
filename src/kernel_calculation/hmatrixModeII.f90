!   
!   Ando scheme used with fast multipole method from Greengard
!
!   Author: Romanet Pierre   (romanet@ipgp.fr)
!
!
!
!
subroutine hmatII(nelement,node,slip,iprec,stressT,stressN,ier,mu)
use variables, only: dimA_T, dimB_T, dimF_T, A_matrices_T, B_matrices_T, F_matrices_T, &
                     dimA_N, dimB_N, dimF_N, A_matrices_N, B_matrices_N, F_matrices_N, &
                     nbclus, ind_clusters, nbblclus, blcluster_tree,                   &
                     info_blclustersT, info_blclustersN, mat_corr,dp,                  &
                     cp, cs
use calculation_hmat
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
complex(kind=dp)                            :: stress_temp_T(nelement)
complex(kind=dp)                            :: stress_temp_N(nelement)
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
!   Initialise parameters
stress_temp_T = cmplx(0.d0,0.0d0,8)
stress_temp_N = cmplx(0.d0,0.0d0,8)
!
!   Convert into H-matrix representation (reorder)
do k=1,nelement
    slip_h(k) = cmplx(slip(mat_corr(k)),0.0d0,8)
end do
!
! Call H matrix multiplication for sigmaT
call Hgemv_right(nbclus,ind_clusters,nbblclus,                &
                blcluster_tree,info_blclustersT,              & 
                alpha,beta,1,nelement,nelement,               &
                slip_h,stress_temp_T,dimA_T,dimB_T,dimF_T,    &
                 A_matrices_T,B_matrices_T,F_matrices_T)
!
! Call H matrix multiplication for sigmaN
call Hgemv_right(nbclus,ind_clusters,nbblclus,                &
                blcluster_tree,info_blclustersN,              & 
                alpha,beta,1,nelement,nelement,               &
                slip_h,stress_temp_N,dimA_N,dimB_N,dimF_N,    &
                A_matrices_N,B_matrices_N,F_matrices_N)
!
!    Convert into real representation (reorder)
do k=1,nelement
    slip(mat_corr(k)) = real(slip_h(k),8)
    stressT(mat_corr(k)) = mu/pi*(1.d0-cs**2/cp**2)*real(stress_temp_T(k),8)
    stressN(mat_corr(k)) = mu/pi*(1.d0-cs**2/cp**2)*real(stress_temp_N(k),8)
end do
 !print*,stressN
! print*,stressT
! stop

end subroutine 
