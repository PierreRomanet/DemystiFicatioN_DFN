module block_cluster_tree
implicit none
private
public :: build_block_cluster_tree


contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
subroutine build_block_cluster_tree(nbclus,ind_clusters,sons_clusters,&
                nbnod,cnod2,nbblclus,blcluster_tree,info_blclusters,eta)
  
  
integer                                                         :: nbnod, nbclus, nbblclus         
real(kind=8)  , dimension(2,nbnod)                              :: cnod2
integer(kind=4), dimension(2*nbclus), intent(inout)             :: ind_clusters
integer(kind=4), dimension(:), allocatable, intent(inout)       :: sons_clusters  
integer(kind=8),dimension(:,:), allocatable ,intent(inout)      :: info_blclusters
integer(kind=4),dimension(:,:), allocatable                     :: blcluster_tree

real(kind=8)                                         :: eta

INTEGER(kind=4)                                      :: AllocateStatus
INTEGER(kind=4)                                      :: ind_cur
  

nbblclus=INT(1,4)

  !Computation of the number of block-clusters in the Block Cluster Tree
  CALL parameters_blcluster_tree(1,1,nbclus,ind_clusters,sons_clusters,&
       nbnod,cnod2,eta,nbblclus)

  !Stop if nbblclus is unrepresentable as INTEGER(kind=4)
  IF(nbblclus.GT.2147483647) STOP "*** See parameter nbblclus ***"
  
  !Allocation of the variable blcluster_tree, containing the block
  !clusters (type, cluster row and cluster column), and of the array
  !info_blclusters, containing the informations about the block clusters
  ALLOCATE(blcluster_tree(nbblclus,3),info_blclusters(nbblclus,6),&
       STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  info_blclusters=INT(0,8)!Inizialize the array info_blclusters
  ind_cur=INT(1,4) !Initialize the index for create_blcluster_tree  
  
  !Construction of the Block Cluster Tree
  CALL create_blcluster_tree(1,1,nbclus,ind_clusters,sons_clusters,&
       nbnod,cnod2,nbblclus,blcluster_tree,info_blclusters,eta,&
       ind_cur)
  
  DEALLOCATE(sons_clusters)
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
RECURSIVE SUBROUTINE create_blcluster_tree(clus_row,clus_col,nbclus,&
     ind_clusters,sons_clusters,nbnod,cnod,nbblclus,blcluster_tree,&
     sons_blclusters,eta,ind_rec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: create_blcluster_tree                             
!                                                                    
! INPUTS:   clus_row      : index of row cluster                     
!           clus_col      : index of column cluster                  
!           nbclus        : number of clusters in the tree           
!           ind_clusters  : array storing the cluster tree           
!           sons_clusters : array storing the sons of the  cluster   
!                           (the value 0 is assigned leaf-clusters)  
!           nbnod         : number of nodes in the problem           
!           cnod          : coordinates of the nodes                 
!           nbblclus      : number of block clusters                   
!           eta           : parameter to check if the block          
!                           admits a low-rank representation           
!                                                                    
! INPUTS/OUTPUTS: blcluster_tree  : matrix storing the informations  
!                                   on Block Cluster Tree            
!                 sons_blclustree : matrix storing the informations  
!                                   on block-cluster affiliation       
!                 ind_rec         : index to recursion               
!                                                                    
! GOAL: create the Block Cluster Tree                                
                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: clus_row, clus_col
  INTEGER(kind=4),INTENT(IN):: nbclus, nbnod, nbblclus

  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: sons_clusters

  REAL(kind=8),INTENT(IN):: eta

  REAL(kind=8),DIMENSION(2,nbnod),INTENT(IN):: cnod
  
  !Input/Output
  INTEGER(kind=4),INTENT(INOUT):: ind_rec

  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(INOUT):: blcluster_tree
  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(INOUT):: sons_blclusters

  !Local Variables
  INTEGER(kind=4):: ind_row, ind_col, ind_sons, temp_ind_rec
  LOGICAL:: par_adm

  !!!!!!!!!!!!!!!!!!!!!!! BODY of the FUNCTION !!!!!!!!!!!!!!!!!!!!!!!
  
  !Check the admissibility condition
  par_adm=admissibility_condition(ind_clusters(2*clus_row-1),&
       ind_clusters(2*clus_row),ind_clusters(2*clus_col-1),&
       ind_clusters(2*clus_col),cnod,nbnod,eta)
  
  SELECT CASE(par_adm)

  CASE(.FALSE.) !Non-admissible cluster-pair

     !Both of the two clusters are non-leaf cluster
     IF((sons_clusters(2*clus_row).NE.0).AND.&
          ((sons_clusters(2*clus_col).NE.0))) THEN

        !Add the block-cluster to the Block Cluster Tree, inserting
        !its nature, its cluster_row and its cluster_column into the
        !array blcluster_tree
        blcluster_tree(ind_rec,1)=2 !Flag for splitted blocks
        blcluster_tree(ind_rec,2)=clus_row
        blcluster_tree(ind_rec,3)=clus_col

        !Initialize the flag for the sons
        ind_sons=0
        !Fix the index of the block cluster
        temp_ind_rec=ind_rec
        
        DO ind_row=-1,0
           DO ind_col=-1,0
              !Updating the flag for the sons
              ind_sons=ind_sons+1
              !Updating the flag for recursion
              ind_rec=ind_rec+1
              !Add the son to sons_blclusters
              sons_blclusters(temp_ind_rec,ind_sons)=INT(ind_rec,8) 
              !Recursive procedure on the sons
              CALL create_blcluster_tree(&
                   sons_clusters(2*clus_row+ind_row),&
                   sons_clusters(2*clus_col+ind_col),nbclus,&
                   ind_clusters,sons_clusters,nbnod,cnod,nbblclus,&
                   blcluster_tree,sons_blclusters,eta,ind_rec)     
           END DO
        END DO

        !One of the two clusters is a leaf cluster  
     ELSE IF ((sons_clusters(2*clus_row).EQ.0).OR.&
          ((sons_clusters(2*clus_col).EQ.0))) THEN
        
        !Add the block-cluster to the Block Cluster Tree, inserting
        !its nature, its cluster_row and its cluster_column into the
        !array blcluster_tree
        blcluster_tree(ind_rec,1)=0 !Flag for non-admissible leaf
        blcluster_tree(ind_rec,2)=clus_row
        blcluster_tree(ind_rec,3)=clus_col

     END IF
     
  CASE(.TRUE.) !Admissible cluster-pair
     
     !Add the block-cluster to the Block Cluster Tree, inserting
     !its nature, its cluster_row and its cluster_column into the
     !array blcluster_tree
     blcluster_tree(ind_rec,1)=1 !Flag for admissible leaf block

     blcluster_tree(ind_rec,2)=clus_row
     blcluster_tree(ind_rec,3)=clus_col
     
  END SELECT

  RETURN
  
END SUBROUTINE create_blcluster_tree
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
RECURSIVE SUBROUTINE parameters_blcluster_tree(clus_row,clus_col,&
     nbclus,ind_clusters,sons_clusters,nbnod,cnod,eta,nbblclus)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: create_blcluster_tree                             
!                                                                    
! INPUTS:   clus_row      : index of row cluster                     
!           clus_col      : index of column cluster                  
!           nbclus        : number of clusters in the tree           
!           ind_clusters  : array storing the cluster tree           
!           sons_clusters : array storing the sons of the  cluster   
!                           (the value 0 is assigned leaf-clusters)  
!           nbnod         : number of nodes in the problem           
!           cnod          : coordinates of the nodes                 
!           eta           : parameter to check if the block          
!                           admits a low-rank representation           
!                                                                    
! INPUT/OUTPUTS: nbblclus : number of block-clusters in the Block    
!                           Cluster Tree                             
!                                                                    
! GOAL: evaluation of the number of block-clusters in the Block      
!       Cluster Tree                                                 
!                                                                    
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/17/2014                                                  
!                                                                    
! MODIFIED: VIII/28/2015 (Luca Desiderio)                                                          
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: clus_row, clus_col
  INTEGER(kind=4),INTENT(IN):: nbclus, nbnod

  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: sons_clusters

  REAL(kind=8),INTENT(IN):: eta

  REAL(kind=8),DIMENSION(2,nbnod),INTENT(IN):: cnod

  !Input/Output
  INTEGER(kind=4),INTENT(INOUT):: nbblclus

  !Local Variables
  INTEGER(kind=4):: ind_row, ind_col
  
  LOGICAL:: par_adm

  !!!!!!!!!!!!!!!!!!!!!!! BODY of the FUNCTION !!!!!!!!!!!!!!!!!!!!!!!

  !Check the admissibility condition
  par_adm=admissibility_condition(ind_clusters(2*clus_row-1),&
       ind_clusters(2*clus_row),ind_clusters(2*clus_col-1),&
       ind_clusters(2*clus_col),cnod,nbnod,eta)

  SELECT CASE(par_adm)

  CASE(.FALSE.) !Non-admissible cluster-pair

     !Both of the two clusters are non-leaf cluster
     IF((sons_clusters(2*clus_row).NE.0).AND.&
          ((sons_clusters(2*clus_col).NE.0))) THEN

        DO ind_row=-1,0
           DO ind_col=-1,0

              !Updating the number of block-cluster in the
              !Block Cluster Tree
              nbblclus=nbblclus+1
             
              !Recursive procedure on the sons
              CALL parameters_blcluster_tree(&
                   sons_clusters(2*clus_row+ind_row),&
                   sons_clusters(2*clus_col+ind_col),nbclus,&
                   ind_clusters,sons_clusters,nbnod,cnod,eta,nbblclus)
             
           END DO
        END DO   
     END IF
  END SELECT

  RETURN

END SUBROUTINE parameters_blcluster_tree
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
LOGICAL FUNCTION admissibility_condition(first_ind_clus_row,&
     last_ind_clus_row,first_ind_clus_col,last_ind_clus_col,cnod,&
     nbnod,eta)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! FUNCTION NAME: admissibility_condition                             
!                                                                    
! INPUTS:  first_ind_clus_ror : index of first node in the row       
!                               cluster                              
!          last_ind_clus_row  : index of last node in the row        
!                               cluster                              
!          first_ind_clus_col : index of first node in the column    
!                               cluster                               
!          last_ind_clus_col  : index of last node in the column     
!                               cluster                              
!          nbnod              : number of nodes                      
!          cnod               : coordinates of the nodes             
!          eta                : parameter to check the admissibility 
!                                                                    
! OUTPUTS: admissibility_condition: parameter whose value is .TRUE.  
!                                   if the pair is admissible, or    
!                                   .FALSE. otherwise                
!                                                                    
! GOAL: decide if a pair of clusters is admissible                   
!                                                                    
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/11/2014                                                  
!                                                                    
! MODIFIED: VIII/28/2015 (Luca Desiderio)                                                        
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: first_ind_clus_row, last_ind_clus_row
  INTEGER(kind=4),INTENT(IN):: first_ind_clus_col, last_ind_clus_col
  INTEGER(kind=4),INTENT(IN):: nbnod

  REAL(kind=8),DIMENSION(2,nbnod),INTENT(IN):: cnod

  !Local Variables
  REAL(kind=8):: eta
  REAL(kind=8):: dist, diam_row, diam_col

  !!!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the FUNCTION !!!!!!!!!!!!!!!!!!!

  !Compute the diameters of the two clusters
  diam_row=diameter(first_ind_clus_row,last_ind_clus_row,cnod,nbnod)
  diam_col=diameter(first_ind_clus_col,last_ind_clus_col,cnod,nbnod)

  !Compute the distance between the two clusters
  dist=distance(first_ind_clus_row,last_ind_clus_row,&
       first_ind_clus_col,last_ind_clus_col,&
       cnod,nbnod)

  !Check the admissibility
  IF(eta*dist.GE.MAX(diam_row,diam_col)) THEN
     admissibility_condition=.TRUE.      !admissible clusters
  ELSE
     admissibility_condition=.FALSE.     !non-admissible clusters
  END IF
  
  RETURN

END FUNCTION admissibility_condition
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
FUNCTION distance(first_ind_clus_row,last_ind_clus_row,&
     first_ind_clus_col,last_ind_clus_col,cnod,nbnod)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! FUNCTION NAME: distance                                            
!                                                                    
! INPUTS:  first_ind_clus_ror : index of first node in the row       
!                               cluster                              
!          last_ind_clus_row  : index of last node in the row        
!                               cluster                              
!          first_ind_clus_col : index of first node in the column    
!                               cluster                               
!          last_ind_clus_col  : index of last node in the column     
!                               cluster                              
!          nbnod              : number of nodes                      
!          cnod               : coordinates of the nodes             
!                                                                    
! OUTPUTS: distance           : distance between the two clusters    
!                                                                    
! GOAL: compute the distance between two clusters                    
!                                                                    
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/11/2014                                                  
!                                                                    
! MODIFIED: VIII/28/2015 (Luca Desiderio)                                                        
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: first_ind_clus_row, last_ind_clus_row
  INTEGER(kind=4),INTENT(IN):: first_ind_clus_col, last_ind_clus_col
  INTEGER(kind=4),INTENT(IN):: nbnod

  REAL(kind=8),DIMENSION(2,nbnod),INTENT(IN):: cnod

  !Output
  REAL(kind=8):: distance 

  !Local Variables
  INTEGER(kind=4):: ind_dir

  !!!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the FUNCTION !!!!!!!!!!!!!!!!!!! 

  !Initialize the value of distance
  distance=0.d0

  DO ind_dir=1,2
     distance=distance+MAX(0.d0,&
      MINVAL(cnod(ind_dir,first_ind_clus_row:last_ind_clus_row))-&
      MAXVAL(cnod(ind_dir,first_ind_clus_col:last_ind_clus_col)))**2+&
      MAX(0.d0,MINVAL(cnod(ind_dir,first_ind_clus_col:last_ind_clus_col))-&
      MAXVAL(cnod(ind_dir,first_ind_clus_row:last_ind_clus_row)))**2
  END DO

  distance=SQRT(distance)
  
  RETURN

END FUNCTION distance
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
FUNCTION diameter(first_ind_clus,last_ind_clus,cnod,nbnod)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! FUNCTION NAME: diameter                                            
!                                                                    
! INPUTS:  first_ind_clus  : index of first node in the cluster      
!          last_ind_clus   : index of last node in the cluster       
!          nbnod           : number of nodes                         
!          cnod            : coordinates of the nodes                
!                                                                    
! OUTPUTS: diameter        : diameter of the cluster                 
!                                                                    
! GOAL: compute the diameter of a given cluster                      
!                                                                    
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/11/2014                                                  
!                                                                    
! MODIFIED: VIII/28/2015 (Luca Desiderio)                                                      
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: first_ind_clus, last_ind_clus
  INTEGER(kind=4),INTENT(IN):: nbnod

  REAL(kind=8),DIMENSION(2,nbnod),INTENT(IN):: cnod

  !Output
  REAL(kind=8):: diameter

  !Local Variables
  INTEGER(kind=4):: ind_dir

  !!!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the FUNCTION !!!!!!!!!!!!!!!!!!! 

  !Initialize the value of diameter
  diameter=0.d0

  DO ind_dir=1,2
     diameter=diameter+(&
          MAXVAL(cnod(ind_dir,first_ind_clus:last_ind_clus))-&
       MINVAL(cnod(ind_dir,first_ind_clus:last_ind_clus)))**2
  END DO

  diameter=SQRT(diameter)
  
  RETURN

END FUNCTION diameter

end module
