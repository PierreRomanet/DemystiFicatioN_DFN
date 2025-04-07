module calculation_hmat
contains
RECURSIVE SUBROUTINE Hgemv_right(nbclus,ind_clusters,nbblclus,&
     blcluster_tree,info_blclusters,alpha,beta,indH,nbrowH,nbcolH,&
     vectx,vecty,dimA,dimB,dimF,A_matrices,B_matrices,F_matrices)
use get_info_matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: Hgemv_right
!                                                                    
! GOAL: the subroutine performs y <- alpha*H*x+beta*y when H is a
!       H-matrix. alpha and beta are scalars, x and y are vectors.
!                                                                    
! INPUTS:  nbclus           : number of clusters
!          nbblclus         : number of block clusters
!          ind_clusters     : array storing the cluster tree
!          blcluster_tree   : matrix storing the block cluster tree
!          info_blclusters  : informations about block clusters
!          alpha            : scalar
!          beta             : scalar
!          indH             : index of the current block
!          nbrowH           : first dimension of current block
!          nbcolH           : second dimension of current block  
!          vectx            : array storing vector x
!          vecty            : array storing vector y  
!          dimA             : dimension of array A_matrices
!          dimB             : dimension of array B_matrices
!          dimF             : dimension of array F_matrices
!          A_matrices       : array storing A matrices
!          B_matrices       : array storing B matrices
!          F_matrices       : array storing F matrices  
!                                                                    
! INPUT/OUTPUTS: vecty      : array storing vector y
!
! BIBLIOGRAPHY: S.Borm, L.Grasedyck, W.Hackbusch, "Introduction to
!               hierarchical matrices with applications", 2003.
!                                                                    
! AUTHOR: Luca Desiderio
!                                                                    
! DATE: I/29/2015
!                                                                    
! MODIFIED: III/16/2015 (Luca Desiderio)
!                                                                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: nbclus, nbblclus
  INTEGER(kind=4),INTENT(IN):: indH, nbrowH, nbcolH
  INTEGER(kind=8),INTENT(IN):: dimA, dimB, dimF

  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(IN):: blcluster_tree
  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(IN):: info_blclusters

  COMPLEX(kind=8),INTENT(IN):: alpha, beta

  COMPLEX(kind=8),DIMENSION(nbcolH),INTENT(IN):: vectx

  COMPLEX(kind=8),DIMENSION(dimA),INTENT(IN):: A_matrices
  COMPLEX(kind=8),DIMENSION(dimB),INTENT(IN):: B_matrices
  COMPLEX(kind=8),DIMENSION(dimF),INTENT(IN):: F_matrices

  !Input/Output
  COMPLEX(kind=8),DIMENSION(nbrowH),INTENT(INOUT):: vecty

  !Local Variables
  INTEGER(kind=4):: nbrowH11, nbcolH11, rankH
  INTEGER(kind=8):: startA, endA, startB, endB, startF, endF
  
  INTEGER(kind=4),DIMENSION(4):: sonsH

  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:):: vectz

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!
  
  SELECT CASE(blcluster_tree(indH,1))
     
  !The current block is not a leaf (splitted matrix) 
  CASE(2)
     
     !The indices of the block sons are stored in array sonsH
     sonsH=INT(info_blclusters(indH,1:4),4)
     
     !Number of row of block H11
     nbrowH11=cmp_nbrow(sonsH(1),nbclus,ind_clusters,nbblclus,&
          blcluster_tree)
      
     !Number of column of block H11
     nbcolH11=cmp_nbcol(sonsH(1),nbclus,ind_clusters,nbblclus,&
          blcluster_tree)

     !Compute y1 <- alpha*(H11*x1)+beta*y1
     CALL Hgemv_right(nbclus,ind_clusters,nbblclus,blcluster_tree,&
          info_blclusters,alpha,beta,sonsH(1),nbrowH11,nbcolH11,&
          vectx(1:nbcolH11),vecty(1:nbrowH11),dimA,dimB,dimF,&
          A_matrices,B_matrices,F_matrices)

     !Compute y1 <- alpha*(H12*x2)+y1
     CALL Hgemv_right(nbclus,ind_clusters,nbblclus,blcluster_tree,&
          info_blclusters,alpha,DCMPLX(1.d0),sonsH(2),&
          nbrowH11,nbcolH-nbcolH11,vectx(nbcolH11+1:nbcolH),&
          vecty(1:nbrowH11),dimA,dimB,dimF,A_matrices,B_matrices,&
          F_matrices)
     
     !Compute y2 <- alpha*(H21*x1)+beta*y2
     CALL Hgemv_right(nbclus,ind_clusters,nbblclus,blcluster_tree,&
          info_blclusters,alpha,beta,sonsH(3),nbrowH-nbrowH11,&
          nbcolH11,vectx(1:nbcolH11),vecty(nbrowH11+1:nbrowH),&
          dimA,dimB,dimF,A_matrices,B_matrices,F_matrices)
     
     !Compute y2 <- alpha*(H22*x2)+y2
     CALL Hgemv_right(nbclus,ind_clusters,nbblclus,blcluster_tree,&
          info_blclusters,alpha,DCMPLX(1.d0),sonsH(4),&
          nbrowH-nbrowH11,nbcolH-nbcolH11,vectx(nbcolH11+1:nbcolH),&
          vecty(nbrowH11+1:nbrowH),dimA,dimB,dimF,A_matrices,&
          B_matrices,F_matrices)

  !The current block is an admissible leaf (low-rank matrix)
  CASE(1)

     !Compute the rank of the block
     rankH=cmp_rank(indH,nbblclus,info_blclusters)

     !Compute the start and the end positions of matrix BT
     startB=cmp_startB(indH,nbblclus,info_blclusters)
     endB=cmp_endB(indH,nbblclus,info_blclusters)
    
     !Allocate array vectz
     ALLOCATE(vectz(rankH))
     
     !Perform z <- (BT*x)+z
     CALL ZGEMV('T',nbcolH,rankH,DCMPLX(1.d0),B_matrices(startB:endB),&
          nbcolH,vectx,1,DCMPLX(0.d0),vectz,1)

     !Compute the start and the end positions of matrix A
     startA=cmp_startA(indH,nbblclus,info_blclusters)
     endA=cmp_endA(indH,nbblclus,info_blclusters)
      
     !Perform y <- alpha*(A*z)+beta*y
     CALL ZGEMV('N',nbrowH,rankH,alpha,A_matrices(startA:endA),&
          nbrowH,vectz,1,beta,vecty,1)

     !Deallocate array vectz
     DEALLOCATE(vectz)

  !The current block is a non-admissible leaf (full-matrix)    
  CASE(0)

     !Compute the start and the end position of matrix F
     startF=cmp_startF(indH,nbblclus,info_blclusters)
     endF=cmp_endF(indH,nbblclus,info_blclusters)
      
     !Perform y <- alpha*(F*x)+beta*y
     CALL ZGEMV('N',nbrowH,nbcolH,alpha,F_matrices(startF:endF),&
          nbrowH,vectx,1,beta,vecty,1)
    
  END SELECT
  
  RETURN   
  
END SUBROUTINE Hgemv_right
end module