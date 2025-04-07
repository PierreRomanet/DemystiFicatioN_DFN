module get_info_matrices
implicit none
contains


FUNCTION cmp_nbrow(index,nbclus,ind_clusters,nbblclus,blcluster_tree)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbclus, nbblclus

  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(IN):: blcluster_tree

  !Output
  INTEGER(kind=4):: cmp_nbrow

  !Local Variable
  INTEGER(kind=4):: clus_row

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Compute the cluster row
  clus_row=2*blcluster_tree(index,2)

  !Compute the number of rows
  cmp_nbrow=ind_clusters(clus_row)-ind_clusters(clus_row-1)+1

  cmp_nbrow=cmp_nbrow  

  RETURN

END FUNCTION cmp_nbrow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_nbcol(index,nbclus,ind_clusters,nbblclus,blcluster_tree)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbclus, nbblclus

  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(IN):: blcluster_tree

  !Output
  INTEGER(kind=4):: cmp_nbcol

  !Local Variable
  INTEGER(kind=4):: clus_col

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Compute the cluster column
  clus_col=2*blcluster_tree(index,3)

  !Compute the number of columns
  cmp_nbcol=ind_clusters(clus_col)-ind_clusters(clus_col-1)+1

  cmp_nbcol=cmp_nbcol 

  RETURN

END FUNCTION cmp_nbcol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_rank(index,nbblclus,info_blclusters)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbblclus

  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(IN):: info_blclusters

  !Output
  INTEGER(kind=4):: cmp_rank

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  cmp_rank=INT(info_blclusters(index,5),4)

  RETURN

END FUNCTION cmp_rank

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_max_rank(index,nbblclus,info_blclusters)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbblclus

  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(IN):: info_blclusters

  !Output
  INTEGER(kind=4):: cmp_max_rank

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  cmp_max_rank=INT(info_blclusters(index,6),4)

  RETURN

END FUNCTION cmp_max_rank

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_startA(index,nbblclus,info_blclusters)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbblclus

  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(IN):: info_blclusters

  !Output
  INTEGER(kind=8):: cmp_startA

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  cmp_startA=info_blclusters(index,1)

  RETURN

END FUNCTION cmp_startA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_endA(index,nbblclus,info_blclusters)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbblclus

  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(IN):: info_blclusters

  !Output
  INTEGER(kind=8):: cmp_endA

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  cmp_endA=info_blclusters(index,2)

  RETURN

END FUNCTION cmp_endA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_startB(index,nbblclus,info_blclusters)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbblclus

  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(IN):: info_blclusters

  !Output
  INTEGER(kind=8):: cmp_startB

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  cmp_startB=info_blclusters(index,3)

  RETURN

END FUNCTION cmp_startB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_endB(index,nbblclus,info_blclusters)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbblclus

  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(IN):: info_blclusters

  !Output
  INTEGER(kind=8):: cmp_endB

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  cmp_endB=info_blclusters(index,4)

  RETURN

END FUNCTION cmp_endB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_startF(index,nbblclus,info_blclusters)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbblclus

  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(IN):: info_blclusters

  !Output
  INTEGER(kind=8):: cmp_startF

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  cmp_startF=info_blclusters(index,1)

  RETURN

END FUNCTION cmp_startF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_endF(index,nbblclus,info_blclusters)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbblclus

  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(IN):: info_blclusters

  !Output
  INTEGER(kind=8):: cmp_endF

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  cmp_endF=info_blclusters(index,2)

  RETURN

END FUNCTION cmp_endF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_startRHS(index,nbclus,nbblclus,ind_clusters,&
     blcluster_tree)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbclus, nbblclus

  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(IN):: blcluster_tree

  !Output
  INTEGER(kind=8):: cmp_startRHS

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

   
  cmp_startRHS=INT((ind_clusters(2*blcluster_tree(index,2)-1)-1)+1,8)
 
  RETURN

END FUNCTION cmp_startRHS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_endRHS(index,nbclus,nbblclus,ind_clusters,blcluster_tree)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: index, nbclus, nbblclus

  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(IN):: blcluster_tree

  !Output
  INTEGER(kind=8):: cmp_endRHS

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !3D elastodynamics  
  cmp_endRHS=INT(ind_clusters(2*blcluster_tree(index,2)),8)

  RETURN

END FUNCTION cmp_endRHS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cmp_stor_lowrank(nbclus,ind_clusters,nbblclus,blcluster_tree,&
     info_blclusters)

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: nbclus, nbblclus

  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(IN):: blcluster_tree
  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(IN):: info_blclusters

  !Output
  INTEGER(kind=8):: cmp_stor_lowrank

  !Local Variables
  INTEGER(kind=4):: ind_block, nbrow, nbcol, rank
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Initialize the output
  cmp_stor_lowrank=INT(0,8)

  DO ind_block=1,nbblclus

     !The current block is an admissible leaf
     IF(blcluster_tree(ind_block,1).EQ.1) THEN

        !Compute the number of rows of the current block
        nbrow=cmp_nbrow(ind_block,nbclus,ind_clusters,nbblclus,&
             blcluster_tree)
        
        !Compute the number of columns of the current block
        nbcol=cmp_nbcol(ind_block,nbclus,ind_clusters,nbblclus,&
             blcluster_tree)
        
        !Compute the rank of the current block
        rank=cmp_rank(ind_block,nbblclus,info_blclusters)

        !Update the value of the output
        cmp_stor_lowrank=cmp_stor_lowrank+&
             INT(nbrow*rank+rank*nbcol,8)
 
     END IF

  END DO
  

  RETURN

END FUNCTION cmp_stor_lowrank

end module
