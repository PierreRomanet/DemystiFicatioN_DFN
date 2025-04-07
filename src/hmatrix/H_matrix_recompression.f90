module H_matrix_recompression
use get_info_matrices
implicit none
private
public :: compress_H_matrix



contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
subroutine compress_H_matrix(nbclus,ind_clusters,nbblclus,&
       blcluster_tree,info_blclusters,dimA,A_matrices,dimB,B_matrices,&
       eps_ACA,stor_lowrank)
  ! Variables
  integer                                                   :: nbclus, nbblclus
  integer(kind=8)                                           :: dimA, dimB
  complex(kind=8), dimension(:), allocatable, intent(INOUT) :: A_matrices
  complex(kind=8), dimension(:), allocatable, intent(INOUT)  :: B_matrices
  integer(kind=8),dimension(nbblclus,6),intent(inout)       :: info_blclusters
  integer(kind=4), dimension(2*nbclus), intent(inout)       :: ind_clusters
  integer(kind=4),dimension(:,:),allocatable                :: blcluster_tree
  integer(kind=8),intent(inout)                             :: stor_lowrank
  real(kind=8)                                              ::  eps_ACA
  real                                              :: comp_rate
  
  !****Step 3: RECOMPRESSION of the H-MATRIX
  PRINT*,'RECOMPRESSION of the MATRIX'

  !Initialize the value of variables stor_lowrank
  stor_lowrank=INT(0,8)

  !Recompression of the H-matrix
  CALL recompression_low_rank_blocks(nbclus,ind_clusters,nbblclus,&
       blcluster_tree,info_blclusters,dimA,A_matrices,dimB,B_matrices,&
       eps_ACA,stor_lowrank)




  ! Compute the new maximum low-rank TODO
!   max_low_rank=MAXVAL(INT(info_blclusters(:,5),4))
! 
!   Compute the new compression rate TODO
   

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE recompression_low_rank_blocks(nbclus,ind_clusters,nbblclus,&
     blcluster_tree,info_blclusters,dimA,A_matrices,dimB,B_matrices,&
     eps,stor_lowrank)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Inputs
  INTEGER(kind=4),INTENT(IN):: nbclus, nbblclus
  INTEGER(kind=8),INTENT(IN):: dimA, dimB

  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(IN):: blcluster_tree

  REAL(kind=8),INTENT(IN):: eps
  
  !Input/Outputs
  INTEGER(kind=8),INTENT(INOUT):: stor_lowrank
  
  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(INOUT):: info_blclusters

  COMPLEX(kind=8),DIMENSION(dimA),INTENT(INOUT):: A_matrices
  COMPLEX(kind=8),DIMENSION(dimB),INTENT(INOUT):: B_matrices

  !Local Variables
  INTEGER(kind=4):: ind_block
  INTEGER(kind=4):: nbrow, nbcol, inrank, outrank


  INTEGER(kind=8):: startA, endA, startB, endB


  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Loop on the block-clusters
  DO ind_block=1,nbblclus

     !Check if the current block is an admissible leaf-block
     IF(blcluster_tree(ind_block,1).EQ.1) THEN

        !Compute the number of rows of the current block
        nbrow=cmp_nbrow(ind_block,nbclus,ind_clusters,&
             nbblclus,blcluster_tree)

        !Compute the number of columns of the current block
        nbcol=cmp_nbcol(ind_block,nbclus,ind_clusters,&
             nbblclus,blcluster_tree)

        !Compute the rank of the current block
        inrank=cmp_rank(ind_block,nbblclus,info_blclusters)

        !Compute the start position of matrices A and B in arrays
        !A_matrices and B_matrices respectively
        startA=cmp_startA(ind_block,nbblclus,info_blclusters)
        startB=cmp_startB(ind_block,nbblclus,info_blclusters)

        !Compute the end position of matrices A and B in arrays 
        !A_matrices and B_matrices respectively
        endA=startA+INT(nbrow*inrank-1,8)
        endB=startB+INT(nbcol*inrank-1,8)

        !Initialize variable outrank
        outrank=INT(0,4)
     
        !Compression of matrices A1 and B1
        CALL compression_low_rank(nbrow,nbcol,inrank,&
             A_matrices(startA:endA),B_matrices(startB:endB),&
             eps,outrank)

        !Store the informations about matrix A in info_blclusters
        info_blclusters(ind_block,2)=startA+INT(nbrow*outrank-1,8)

        !Store the informations about matrix B in info_blclusters
        info_blclusters(ind_block,4)=startB+INT(nbcol*outrank-1,8)

        !Store the informations about the low-rank of the current 
        !block
        info_blclusters(ind_block,5)=INT(outrank,8)

        !Update the value of variable stor_lowrank
        stor_lowrank=stor_lowrank+INT((nbrow+nbcol)*outrank,8)

     END IF

  END DO

  RETURN

END SUBROUTINE recompression_low_rank_blocks
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE compression_low_rank(nbrow,nbcol,inrank,vectA,vectB,&
     eps,outrank)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Inputs
  INTEGER(kind=4),INTENT(IN):: nbrow, nbcol, inrank

  REAL(kind=8),INTENT(IN):: eps

  !Input/Outputs
  INTEGER(kind=4),INTENT(INOUT):: outrank

  COMPLEX(kind=8),DIMENSION(nbrow*inrank),INTENT(INOUT):: vectA
  COMPLEX(kind=8),DIMENSION(nbcol*inrank),INTENT(INOUT):: vectB

  !Local Variables
  INTEGER(kind=4):: ind_row, ind_col, endA, endB, posR

  COMPLEX(kind=8),DIMENSION(inrank):: rowRa, rowRb
  
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:):: tauA, tauB
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:):: vectU, vectV, vectF
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:):: vectQa, vectQb, vectR

  !Function Variables
  COMPLEX(kind=8),EXTERNAL:: ZDOTU

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  IF(inrank.LT.MIN(nbrow,nbcol)) THEN
  
     !Allocate arrays tauA and tauB
     ALLOCATE(tauA(MIN(nbrow,inrank)),tauB(MIN(nbcol,inrank)))

     !Compute QR-decomposition of matrix A, i.e. A=Qa*Ra
     CALL new_compute_QR(nbrow,inrank,vectA,tauA)

     !Compute QR-decomposition of matrix B, i.e. B=Qb*Rb
     CALL new_compute_QR(nbcol,inrank,vectB,tauB)

     !Allocate array vectR
     ALLOCATE(vectR(inrank*inrank))

     !Initialize matrix R=Ra*Rb**(T)
     vectR=DCMPLX(0.d0)

     !Compute end positions of matrices A and B
     endA=nbrow*inrank
     endB=nbcol*inrank

     !Compute matrix R=Ra*Rb**(T)
     DO ind_col=1,inrank

        rowRb=DCMPLX(0.d0)
        rowRb(ind_col:inrank)=&
             vectB((nbcol+1)*ind_col-nbcol:endB:nbcol)

        posR=(ind_col-1)*inrank

        DO ind_row=1,inrank

           rowRa=DCMPLX(0.d0)
           rowRa(ind_row:inrank)=&
                vectA((nbrow+1)*ind_row-nbrow:endA:nbrow)

           vectR(posR+ind_row)=ZDOTU(inrank,rowRa,1,rowRb,1)

        END DO
    
     END DO

     !Allocate arrays vectUr and vectVr
     ALLOCATE(vectU(inrank*inrank))
     ALLOCATE(vectV(inrank*inrank))

     !Initialize arrays vectU and vectV
     vectU=DCMPLX(0.d0)
     vectV=DCMPLX(0.d0)

     !Initialize variable new_rank
     outrank=INT(0,4)

     !Compute the SVD of (small) matrix R
     CALL SVD_factorization(inrank,inrank,vectR,inrank,inrank,&
          vectU,inrank,inrank,vectV,eps,outrank)

     !Allocate array vectQa
     ALLOCATE(vectQa(nbrow*inrank))

     !Copy array vectA into array vectQa
     vectQa=vectA

     !Re-initialize array vectA
     vectA=DCMPLX(0.d0)

     !Assemble the leading inrank columns of matrix Qa
     CALL assemble_matrix_Q(nbrow,inrank,vectQa,tauA)

     !Compute A <- Qa*U
     CALL ZGEMM('N','N',nbrow,outrank,inrank,DCMPLX(1.d0),vectQa,&
          nbrow,vectU,inrank,DCMPLX(0.d0),vectA(1:nbrow*outrank),&
          nbrow)

     !Deallocate arrays tauA, vectU and vectQa
     DEALLOCATE(tauA,vectU,vectQa)

     !Allocate array vectQb
     ALLOCATE(vectQb(nbcol*inrank))

     !Copy array vectB into array vectQb
     vectQb=vectB

     !Re-initialize array vectA
     vectB=DCMPLX(0.d0)

     !Assemble the leading inrank columns of matrix Qb
     CALL assemble_matrix_Q(nbcol,inrank,vectQb,tauB)

     !Compute B <- Qb*V
     CALL ZGEMM('N','N',nbcol,outrank,inrank,DCMPLX(1.d0),vectQb,&
          nbcol,vectV,inrank,DCMPLX(0.d0),vectB(1:nbcol*outrank),&
          nbcol)

     !Deallocate arrays tauB, vectV and vectQb
     DEALLOCATE(tauB,vectV,vectQb)

  ELSE

     !Allocate array vectF
     ALLOCATE(vectF(nbrow*nbcol))

     !Initialize array vectF
     vectF=DCMPLX(0.d0)

     !Compute F <- A*B**(T)
     CALL ZGEMM('N','T',nbrow,nbcol,inrank,DCMPLX(1.d0),vectA,&
          nbrow,vectB,nbcol,DCMPLX(0.d0),vectF,nbrow)

     !Re-initialize arrays vectA and vectB
     vectA=DCMPLX(0.d0)
     vectB=DCMPLX(0.d0)

     !Initialize variable new_rank
     outrank=INT(0,4)

     !Compute full-pivoted ACA of matrix F
     CALL ACA_full_pivot(vectF,nbrow,nbcol,vectA,INT(nbrow*inrank,8),&
          vectB,INT(nbcol*inrank,8),eps,outrank)

  END IF

  RETURN

END SUBROUTINE compression_low_rank
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE ACA_full_pivot(vectF,nbrowF,nbcolF,vectA,dimA,vectB,dimB,&
     eps_ACA,rank)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: ACA_full_pivot   
!
! GOAL: Compute the low-rank approximation of a matrix F using the  
!       ACA (Adaptive Cross Approximation) with full pivoting  
!                                                                    
! INPUTS: nbrowF  : number of rows of matrix F
!         nbcolF  : number of columns of matrix F
!         vectF   : matrix F to approximate (stored column-wise)
!         dimA    : dimension of array vectA
!         dimB    : dimension of array vectB
!         eps_ACA : ACA threshold          
!
! INPUT/OUTPUTS: vectA : array where the matrix A is stored
!                vectB : array where the matrix B is stored
!
! OUTPUT:        rank  : low-rank of matrices A and B  
!
! BIBLIOGRAPHY: L.Grasedyck, "Adaptive Recompression of H-Matrices
!               for BEM", 2005 
! 
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: VIII/19/2015                                                   
!                                                                    
! MODIFIED: 
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Inputs
  INTEGER(kind=4),INTENT(IN):: nbrowF, nbcolF
  INTEGER(kind=8),INTENT(IN):: dimA, dimB

  REAL(kind=8),INTENT(IN):: eps_ACA

  !Output
  INTEGER(kind=4),INTENT(OUT):: rank

  !Input/Outputs
  COMPLEX(kind=8),DIMENSION(nbrowF*nbcolF),INTENT(INOUT):: vectF

  COMPLEX(kind=8),DIMENSION(dimA),INTENT(INOUT):: vectA
  COMPLEX(kind=8),DIMENSION(dimB),INTENT(INOUT):: vectB

  !Local Variables
  INTEGER(kind=4):: ind, ind_i, ind_j

  REAL(kind=8):: maxR, norm_froR, norm_froF

  COMPLEX(kind=8):: gamma

  !Function Variables
  INTEGER,EXTERNAL:: IZAMAX

  REAL(kind=8),EXTERNAL:: DZASUM, ZLANGE

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Norm of the input matrix
  norm_froF=DSQRT(SUM(CDABS(vectF)**2,nbrowF*nbcolF))

  DO rank=1,MIN(nbrowF,nbcolF)

     !Compute the maximum element of the error matrix
     maxR=MAXVAL(CDABS(vectF));

     IF(maxR.EQ.0.d0) THEN
        RETURN
     END IF

     !Compute the index of the maximum element of the error matrix
     ind=IZAMAX(nbrowF*nbcolF,vectF,1)

     !Find the indices of the pivot element
     IF(MOD(ind,nbrowF).EQ.0) THEN
        ind_j=ind/nbrowF
     ELSE
        ind_j=ind/nbrowF+1
     END IF
     ind_i=ind-(ind_j-1)*nbrowF

     !Normalising constant
     gamma=DCMPLX(1.d0)/vectF(ind)

     !New column
     vectA((rank-1)*nbrowF+1:rank*nbrowF)=gamma*&
          vectF((ind_j-1)*nbrowF+1:ind_j*nbrowF)

     !New row
     vectB((rank-1)*nbcolF+1:rank*nbcolF)=&
          vectF(ind_i:nbrowF*nbcolF:nbrowF)

     !New error matrix
     CALL ZGEMM('N','N',nbrowF,nbcolF,1,DCMPLX(-1.d0),&
          vectA((rank-1)*nbrowF+1:rank*nbrowF),nbrowF,&
          vectB((rank-1)*nbcolF+1:rank*nbcolF),1,&
          DCMPLX(1.d0),vectF,nbrowF)

     !Norm of the new error matrix
     norm_froR=DSQRT(SUM(CDABS(vectF)**2,nbrowF*nbcolF))
     
     !Stopping criterion
     IF(norm_froR.LE.eps_ACA*norm_froF) THEN
        EXIT
     END IF

  END DO

  RETURN

END SUBROUTINE ACA_full_pivot
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE assemble_matrix_Q(m1,k,vectF,tau)

  IMPLICIT NONE
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: m1, k

  COMPLEX(kind=8),DIMENSION(MIN(m1,k)),INTENT(IN):: tau

  !Input/Output
  
  COMPLEX(kind=8),DIMENSION(m1*k),INTENT(INOUT):: vectF

  !Local Variables
  INTEGER(kind=4):: info, lwork, lwmax

  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:):: work

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!
  
  !Set the parameter lwmax
  lwmax=1000

  !Allocate array work
  ALLOCATE(work(lwmax))

  !Query the optimal workspace
  lwork=-1
  
  !Compute the optimal size of the work array
  CALL ZUNGQR(m1,k,k,vectF,m1,tau,work,lwork,info)
 
  !The optimal value of lwork is returned as the first entry of array
  !work
  lwork=INT(work(1))

  !Deallocate array work
  DEALLOCATE(work)

  !Allocate array work
  ALLOCATE(work(lwork))
  
  !////////////////////////////////////////////////////////////////!
 
  !Compute first k columns of matrix Q such that F=Q*R
  CALL ZUNGQR(m1,k,k,vectF,m1,tau,work,lwork,info)
 
  !Deallocate array work
  DEALLOCATE(work)

  RETURN

END SUBROUTINE assemble_matrix_Q
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE SVD_factorization(nbrowF,nbcolF,vectF,nbrowA,nbcolA,&
     vectA,nbrowB,nbcolB,vectB,eps,rank)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: new_SVD_compression
!
! GOAL: the subroutine provides the low-rank approximation of a given
!       matrix F, using its Singular Values Decomposition (SVD) in
!       according with a given accuracy.  
!  
! INPUTS:  vectF   : array storing matrix F
!          nbrowF  : number of rows of matrix F
!          nbcolF  : number of columns of matrix F
!          nbrowA  : number of rows of matrix A
!          nbcolA  : number of columns of matrix A  
!          nbrowB  : number of rows of matrix B
!          nbcolB  : number of columns of matrix B
!          eps     : given accuracy to perform H-matrix arithmetic  
!                                                                    
! INPUT/OUTPUTS: vectA, vectB  : arrays storing matrices A and B
!                                such that F=A*B**(T)
!                rank          : low-rank of matrices A and B
!
! BIBLIOGRAPHY:
! 1. S.Borm, L.Grasedyck, W.Hackbusch, "Introduction to hierarchical
!    matrices with applications", 2003.
! 2. A.K.Saibaba, S.Ambikasaran, J.Yue Li, P.K. Kitanidis, E.F. Darve
!    "Application of Hierarchical Matrices to Linear Inverse Problems
!    in Geostatistics", 2012.  
!                                                                 
! AUTHOR: Luca Desiderio
!                                                                    
! DATE: I/28/2015
!                                                                    
! MODIFIED: V/8/2015 (Luca Desiderio)
!                                                                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Inputs
  INTEGER(kind=4),INTENT(IN):: nbrowF, nbcolF
  INTEGER(kind=4),INTENT(IN):: nbrowA, nbcolA, nbrowB, nbcolB

  REAL(kind=8),INTENT(IN):: eps

  COMPLEX(kind=8),DIMENSION(nbrowF*nbcolF),INTENT(IN):: vectF

  !Input/Outputs
  INTEGER(kind=4),INTENT(INOUT):: rank

  COMPLEX(kind=8),DIMENSION(nbrowA*nbcolA),INTENT(INOUT):: vectA
  COMPLEX(kind=8),DIMENSION(nbrowB*nbcolB),INTENT(INOUT):: vectB

  !Local Variables
  INTEGER(kind=4):: ind_sv
  
  REAL(kind=8),ALLOCATABLE,DIMENSION(:):: vectS

  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:):: vectU, vectVH

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Nullify the value of variable rank
  rank=0
  
  !Allocate arrays vectU, vectVH and vectS
  ALLOCATE(vectU(nbrowF*nbrowF),vectVH(nbcolF*nbcolF))
  ALLOCATE(vectS(MIN(nbrowF,nbcolF)))
  
  !Compute SVD of matrix F
  CALL new_compute_SVD(vectF,nbrowF,nbcolF,vectU,vectVH,vectS)

  !Compute the low-rank of the approximation
  DO ind_sv=1,MIN(nbrowF,nbcolF)

     !Update the rank
     rank=rank+1
     
     !Check if the allocated memory is enough
     !IF((rank.GT.nbcolA).OR.(rank.GT.nbcolB)) THEN

     !   rank=-1
     !   RETURN

     !END IF

     !Compute the current column of matrix A
     vectA((ind_sv-1)*nbrowA+1:ind_sv*nbrowA)=&
          DCMPLX(vectS(ind_sv))*&
          vectU((ind_sv-1)*nbrowA+1:ind_sv*nbrowA)
     
     !Compute the current row of matrix B
     vectB((ind_sv-1)*nbrowB+1:ind_sv*nbrowB)=&
          vectVH(ind_sv:nbcolF*nbcolF:nbcolF)

     !Check if the current singular value is significant
     IF(vectS(ind_sv).LE.eps*vectS(1)) THEN
        EXIT
     END IF

  END DO

  !Deallocate arrays vectU, vectVH and vectS
  DEALLOCATE(vectS,vectU,vectVH)
  
  RETURN

END SUBROUTINE SVD_factorization
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE new_compute_QR(nbrow,nbcol,vectF,tau)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: new_compute_QR
!
! GOAL: the subroutine forms the QR factorization of an arbitrary
!       rectangular complex m by n matrix F.
!  
! INPUTS:  nbrow : number of rows of matrix F
!          nbcol : number of columns of matrix F  
!                                                                    
! INPUT/OUTPUTS: vectF : array storing F overwritten by Q and R
!                        factors  
!                tau   : array containing further details of the
!                        unitary matrix Q.  
!
! BIBLIOGRAPHY: V.Comincioli, "Metodi numerici e statistici per le
!               Scienze Applicate", 2004.
!                                                                 
! AUTHOR: Luca Desiderio
!                                                                    
! DATE: I/27/2015
!                                                                    
! MODIFIED: V/8/2015 (Luca Desiderio)
!                                                                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  IMPLICIT NONE
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: nbrow, nbcol

  !Input/Output
  COMPLEX(kind=8),DIMENSION(MIN(nbrow,nbcol)):: tau
  COMPLEX(kind=8),DIMENSION(nbrow*nbcol),INTENT(INOUT):: vectF

  !Local Variables
  INTEGER(kind=4):: info, lwork, lwmax

  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:):: work
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!! 
  
  !Set the parameter lwmax
  lwmax=1000

  !Allocate array work
  ALLOCATE(work(lwmax))

  !Query the optimal workspace
  lwork=-1
  
  !Compute the optimal size of the work array
  CALL ZGEQRF(nbrow,nbcol,vectF,nbrow,tau,work,lwork,info)
  
  !The optimal value of lwork is returned as the first entry of array
  !work
  lwork=INT(work(1))

  !Deallocate array work
  DEALLOCATE(work)
  
  !Allocate array work
  ALLOCATE(work(lwork))
  
  !////////////////////////////////////////////////////////////////!
  
  !QR-decomposition
  CALL ZGEQRF(nbrow,nbcol,vectF,nbrow,tau,work,lwork,info)
  
  !Check for convergence
  IF(info.GT.0) THEN
     PRINT*,'ERROR, see compute_QR'
  END IF

  !Deallocate array work
  DEALLOCATE(work)
  
  RETURN

END SUBROUTINE new_compute_QR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE new_compute_SVD(vectF,nbrowF,nbcolF,vectU,vectVH,vectS)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: compute_SVD
!
! GOAL: the subroutine provides the Singular Values Decomposition
!       (SVD) of a given matrix m-by-n F, i.e. F=U*S*V**(H), where
!       U and V are unitary matrices and S is a diagonal matrix,
!       whose values are the singular values of F.  
!  
! INPUTS:  nbrowF : number of rows of matrix F
!          nbcolF : number of columns of matrix F
!          vectF  : array storing matrix F    
!                                                                    
! INPUT/OUTPUTS: vectU  : array storing the m-by-m unitary matrix U
!                vectVH : array storing the n-by-n unitary matrix VH
!                vectS  : array containing the singular values of F,
!                         stored so that S(i)>=S(i+1)
!
! BIBLIOGRAPHY: V.Comincioli, "Metodi numerici e statistici per le
!               Scienze Applicate", 2004.
!                                                                 
! AUTHOR: Luca Desiderio
!                                                                    
! DATE: I/27/2015
!                                                                    
! MODIFIED: V/9/2015 (Luca Desiderio)
!                                                                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Inputs
  INTEGER(kind=4),INTENT(IN):: nbrowF, nbcolF

  COMPLEX(kind=8),DIMENSION(nbrowF*nbcolF),INTENT(IN):: vectF

  !Input/Outputs
  REAL(kind=8),DIMENSION(MIN(nbrowF,nbcolF)),INTENT(INOUT):: vectS

  COMPLEX(kind=8),DIMENSION(nbrowF*nbrowF),INTENT(INOUT):: vectU
  COMPLEX(kind=8),DIMENSION(nbcolF*nbcolF),INTENT(INOUT):: vectVH

  !Local Variables
  INTEGER(kind=4):: info, lwork, lwmax

  REAL(kind=8),DIMENSION(5*MIN(nbrowF,nbcolF)):: rwork
  
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:):: work
 
  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!
  
  !Set the parameter lwmax
  lwmax=1000

  !Allocate array work
  ALLOCATE(work(lwmax))

  !Query the optimal workspace
  lwork=-1
  
  !Compute the optimal size of the work array
  CALL ZGESVD('A','A',nbrowF,nbcolF,vectF,nbrowF,vectS,vectU,nbrowF,&
       vectVH,nbcolF,work,lwork,rwork,info)
 
  !The optimal value of lwork is returned as the first entry of array
  !work
  lwork=INT(work(1))

  !Deallocate array work
  DEALLOCATE(work)

  !Allocate array work
  ALLOCATE(work(lwork))
  
  !////////////////////////////////////////////////////////////////!
 
  !Compute SVD of the matrix F
  CALL ZGESVD('A','A',nbrowF,nbcolF,vectF,nbrowF,vectS,vectU,nbrowF,&
       vectVH,nbcolF,work,lwork,rwork,info)
 
  !Deallocate array work
  DEALLOCATE(work)
 
  RETURN

END SUBROUTINE new_compute_SVD
end module 
