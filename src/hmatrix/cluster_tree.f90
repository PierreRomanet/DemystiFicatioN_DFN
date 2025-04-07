module cluster_tree
implicit none
private
public :: build_cluster_tree 
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
subroutine build_cluster_tree (nbnod,nbclus,ind_clusters,sons_clusters,&
                               cnod2,nleaf,mat_corr)

!Local variables
INTEGER(kind=4)::AllocateStatus
INTEGER depth,ind_cur
INTEGER ind
INTEGER(kind=4),DIMENSION(:),ALLOCATABLE                 :: sons_clusters
integer(kind=4), dimension(:), allocatable, intent(inout):: ind_clusters
integer(kind=4)                                          :: nbnod, nbclus
real(kind=8)  , dimension(2,nbnod)                       :: cnod2 
integer(kind=4)                                          :: nleaf

COMPLEX(kind=8),DIMENSION(:),ALLOCATABLE::vect_RHS2 
INTEGER(kind=4), dimension(:)  , allocatable       ::  mat_corr

ALLOCATE(vect_RHS2(nbnod),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
  vect_RHS2=CMPLX(0.D0,0.D0,8)!Initialize array vect_RHS2

  nbclus=INT(1,4)!Initialize the number of clusters in the tree

  !print*,nleaf

  ALLOCATE(mat_corr(nbnod))
  DO ind=1,nbnod
     mat_corr(ind)=ind
  END DO
  
  !Computation of the number of clusters in the Cluster Tree
  CALL parameters_cluster_tree(1,nbnod,nbnod,cnod2,vect_RHS2,&
       nleaf,nbclus,mat_corr)

  !Stop if nbclus is unrepresentable as INTEGER(kind=4)
  IF(nbclus.GT.2147483647) STOP "*** See parameter nbclus ***"

  !Allocation of the array ind_clusters, needed to represent the Cluster
  !Tree and containing for each cluster the first and the last nodes
  !(REMARK: each cluster has at most two sons), and of the array
  !sons_clusters, containing the sons of each cluster
  ALLOCATE(ind_clusters(2*nbclus),sons_clusters(2*nbclus),&
       STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  depth=INT(0,4)   !Initialize the depth of the tree
  ind_cur=INT(1,4) !Initialize the index for create_cluster_tree
  
  !Construction of the Cluster Tree
  CALL create_cluster_tree(1,nbnod,nbclus,ind_clusters,sons_clusters,&
       nbnod,cnod2,vect_RHS2,nleaf,ind_cur,0,depth,mat_corr)
  DEALLOCATE(vect_RHS2)





end subroutine 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
RECURSIVE SUBROUTINE create_cluster_tree(first_ind_clus,&
     last_ind_clus,nbclus,ind_clusters,sons_clusters,nbnod,cnod,&
     vect_RHS,nleaf,ind_rec,lev,depth,mat_corr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: create_cluster_tree                               
!                                                                    
! INPUTS:   first_ind_clus  : index of first node in the cluster     
!           last_ind_clus   : index of last node in the cluster      
!           nbclus          : number of clusters in the tree         
!           nbnod           : number of nodes in the problem         
!           nleaf           : cardinality of the leaf-clusters         
!           lev             : level of the current cluster       
!                                                                    
! INPUT/OUTPUTS: ind_clusters  : array storing the cluster tree      
!                sons_clusters : array storing the sons of the       
!                                cluster (the value 0 is assigned    
!                                leaf-clusters)                      
!                cnod          : coordinates of the nodes            
!                vect_RHS      : RHS of the problem                  
!                ind_rec       : index for recursion                  
!                depth         : depth of the tree               
!                                                                    
! GOAL: construction of the Cluster Tree                             
!                                                                    
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/5/2014                                                   
!                                                                    
! MODIFIED: VIII/28/2015 (Luca Desiderio)                             
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: nbclus, nbnod, nleaf, lev
  INTEGER(kind=4),INTENT(IN):: first_ind_clus, last_ind_clus

  !Input/Output
  INTEGER(kind=4),INTENT(INOUT):: ind_rec, depth
  INTEGER(kind=4),DIMENSION(nbnod),INTENT(INOUT)::mat_corr
  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(INOUT):: ind_clusters
  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(INOUT):: sons_clusters

  REAL(kind=8),DIMENSION(2,nbnod),INTENT(INOUT):: cnod 

  COMPLEX(kind=8),DIMENSION(nbnod),INTENT(INOUT):: vect_RHS

  !Local Variables
  INTEGER(kind=4):: temp_ind_rec 
  
  INTEGER(kind=4),DIMENSION(4):: ind_sons
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Update the depth of the tree
  IF(lev.GT.depth) THEN
     depth=lev
  END IF
  
  !Add the cluster to the tree, inserting its first and its last node
  !into the array ind_clusters
  ind_clusters(2*ind_rec-1)=first_ind_clus
  ind_clusters(2*ind_rec)=last_ind_clus
   
  !Check if the cluster is a leaf (if the number of nodes is smaller
  !than the input parameter nleaf)
  IF((last_ind_clus-first_ind_clus+1).LE.nleaf) THEN 

     sons_clusters(2*ind_rec-1)=0
     sons_clusters(2*ind_rec)=0

     RETURN

  ELSE !For each non-leaf cluster

     !Cluster affiliation
     CALL cluster_affiliation(first_ind_clus,last_ind_clus,&
          ind_sons,cnod,vect_RHS,nbnod,.FALSE.,mat_corr)

     !Fix the index of the cluster
     temp_ind_rec=ind_rec

     !Updating of the parameter ind_rec
     ind_rec=ind_rec+1

     !Add the first son to the array sons_clusters
     sons_clusters(2*temp_ind_rec-1)=ind_rec
     
     !Recursive procedure on the first son
     CALL create_cluster_tree(ind_sons(1),ind_sons(2),nbclus,&
          ind_clusters,sons_clusters,nbnod,cnod,vect_RHS,nleaf,&
          ind_rec,lev+1,depth,mat_corr)

     !Updating of the parameter ind_rec
     ind_rec=ind_rec+1

     !Add the second son to the array sons_clusters
     sons_clusters(2*temp_ind_rec)=ind_rec
     
     !Recursive procedure on the second son
     CALL create_cluster_tree(ind_sons(3),ind_sons(4),nbclus,&
     ind_clusters,sons_clusters,nbnod,cnod,vect_RHS,nleaf,&
          ind_rec,lev+1,depth,mat_corr)

  END IF

  RETURN
  
END SUBROUTINE create_cluster_tree
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
RECURSIVE SUBROUTINE parameters_cluster_tree(first_ind_clus,&
     last_ind_clus,nbnod,cnod,vect_RHS,nleaf,nbclus,mat_corr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: max_parameters_cluster_tree                       
!                                                                    
! INPUTS:   first_ind_clus   : index of first node in the cluster    
!           last_ind_clus    : index of last node in the cluster     
!           nleaf            : cardinality of the leaf-clusters      
!           nbnod            : number of nodes                         
!                                                                    
! INPUT/OUTPUTS:  cnod     : coordinates of the nodes                
!                 vect_RHS : RHS of the problem                      
!                 nbclus   : number of the cluster in the tree         
!                                                                    
! GOAL: compute the number of recursive subdivisions of the set of   
!       discretization nodes necessary to build the cluster tree,    
!       using an algorithm based on geometric bisection.                                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Input

  INTEGER(kind=4),INTENT(IN):: first_ind_clus, last_ind_clus
  INTEGER(kind=4),INTENT(IN):: nbnod, nleaf

  !Input/Output
  INTEGER(kind=4),INTENT(INOUT):: nbclus
  INTEGER(kind=4),DIMENSION(nbnod),INTENT(INOUT):: mat_corr
  
  REAL(kind=8),DIMENSION(2,nbnod),INTENT(INOUT):: cnod

  COMPLEX(kind=8),DIMENSION(nbnod),INTENT(INOUT):: vect_RHS

  !Local Variables
  INTEGER(kind=4),DIMENSION(4):: ind_sons

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  
  !Check if the cluster is a leaf (if the number of nodes is smaller
  !than the input parameter nleaf)
  IF((last_ind_clus-first_ind_clus+1).LE.nleaf) THEN 
     
     RETURN
     
  ELSE !For each non-leaf cluster

     !Cluster affiliation
     CALL cluster_affiliation(first_ind_clus,last_ind_clus,&
          ind_sons,cnod,vect_RHS,nbnod,.TRUE.,mat_corr)

     !Update the number of cluster in the Tree
     nbclus=nbclus+2
    
     !Recursive procedure on the first son
     CALL parameters_cluster_tree(ind_sons(1),ind_sons(2),&
          nbnod,cnod,vect_RHS,nleaf,nbclus,mat_corr)

     !Recursive procedure on the second son
     CALL parameters_cluster_tree(ind_sons(3),ind_sons(4),&
          nbnod,cnod,vect_RHS,nleaf,nbclus,mat_corr)
     
  END IF

  RETURN
  
END SUBROUTINE parameters_cluster_tree
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE  cluster_affiliation(first_ind_clus,last_ind_clus,ind_sons,&
     cnod,vect_RHS,nbnod,flag,mat_corr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: cluster_affiliation                               
!                                                                    
! INPUTS:  first_ind_clus  : index of first node in the cluster      
!          last_ind_clus   : index of last node in the cluster       
!          nbnod           : number of nodes                         
!          flag            : parameter whose values are .TRUE. (the  
!                            unknowns have to be reorganized) or     
!                            .FALSE. (the unknowns do not have to be 
!                            reorganized)                              
!                                                                    
! OUTPUTS: ind_sons        : array where the first and the last      
!                            nodes of each son are stored            
!                                                                    
! INPUT/OUTPUTS: vect_RHS  : RHS of the problem                      
!                cnod     : coordinates of the nodes                 
!                                                                    
! GOAL: compute the two sons of a given cluster                      
!                                                                    
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/5/2014                                                   
!                                                                    
! MODIFIED: VIII/28/2015 (Luca Desiderio)                                              
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Input
  INTEGER(kind=4),INTENT(IN):: nbnod
  INTEGER(kind=4),INTENT(IN):: first_ind_clus, last_ind_clus

  LOGICAL,INTENT(IN):: flag

  !Output
  INTEGER(kind=4),DIMENSION(4),INTENT(OUT):: ind_sons

  !Input/Output
  INTEGER(kind=4),DIMENSION(nbnod),INTENT(INOUT):: mat_corr
  
  REAL(kind=8),DIMENSION(2,nbnod),INTENT(INOUT):: cnod

  COMPLEX(kind=8),DIMENSION(nbnod),INTENT(INOUT):: vect_RHS

  !Local Variables
  INTEGER(kind=4):: max_dir

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Compute the direction of maximum expansion of the cluster
  CALL compute_maximum_direction(first_ind_clus,last_ind_clus,cnod,&
       nbnod,max_dir)

  !Find the two sons
  CALL compute_cluster_affiliation(first_ind_clus,last_ind_clus,&
       max_dir,ind_sons,cnod,vect_RHS,nbnod,flag,mat_corr)

  RETURN

END SUBROUTINE cluster_affiliation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE  compute_maximum_direction(first_ind_clus,last_ind_clus,&
     cnod,nbnod,max_dir)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: compute_maximum_direction                         
!                                                                    
! INPUTS:  first_ind_clus  : index of first node in the cluster      
!          last_ind_clus   : index of last node in the cluster       
!          nbnod           : number of nodes                         
!          cnod            : coordinates of the nodes                
!                                                                    
! OUTPUTS: max_dir         : direction of maximum expansion          
!                                                                    
! GOAL: compute the direction of maximum expansion of the bounding   
!       box enclosing a given cluster                                
!                                                                    
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/5/2014                                                   
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
  INTEGER(kind=4),INTENT(OUT):: max_dir

  !Local Variables
  INTEGER(kind=4):: ind_dir
  
  REAL(kind=8),DIMENSION(2):: boun_box

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Bounding Box enclosing the cluster

  boun_box=(/MAXVAL(cnod(1,first_ind_clus:last_ind_clus))-&
       MINVAL(cnod(1,first_ind_clus:last_ind_clus)),&
       MAXVAL(cnod(2,first_ind_clus:last_ind_clus))-&
       MINVAL(cnod(2,first_ind_clus:last_ind_clus))/)

  !Initialize the variable max_dir
  max_dir=1
  
  !Find the direction of maximum expasion
  
  DO ind_dir=1,2
     IF(boun_box(ind_dir).GE.boun_box(max_dir)) THEN
        max_dir=ind_dir
     END IF
  END DO

  RETURN
  
END SUBROUTINE compute_maximum_direction
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE  compute_cluster_affiliation(first_ind_clus,last_ind_clus,&
     max_dir,ind_sons,cnod,vect_RHs,nbnod,flag,mat_corr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: compute_cluster_affiliation                       
!                                                                    
! INPUTS: first_ind_clus   : index of first node in the cluster      
!         last_ind_clus    : index of last node in the cluster       
!         max_dir          : maximum direction of expansion          
!         nbnod            : number of nodes                         
!         flag             : parameter whose values are .TRUE. (the  
!                            unknowns have to be reorganized) or     
!                            .FALSE. (the unknowns do not have to be 
!                            reorganized)                              
!                                                                    
! INPUT/OUTPUTS:  ind_sons : array where the first and the last      
!                            nodes of each son are stored            
!                 cnod     : coordinates of the nodes                
!                 vect_RHS : RHS of the problem                      
!                                                                    
! GOAL: compute the two sons of a given cluster                      
!                                                                    
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/5/2014                                                   
!                                                                    
! MODIFIED: - VIII/28/2015 (Luca Desiderio)                           
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  
  INTEGER(kind=4),INTENT(IN):: nbnod
  INTEGER(kind=4),INTENT(IN):: first_ind_clus, last_ind_clus
  INTEGER(kind=4),INTENT(IN):: max_dir

  LOGICAL,INTENT(IN):: flag

  !Input/Output
  INTEGER(kind=4),DIMENSION(4),INTENT(INOUT):: ind_sons
  INTEGER(kind=4),DIMENSION(nbnod),INTENT(INOUT)::mat_corr
  
  REAL(kind=8),DIMENSION(2,nbnod),INTENT(INOUT):: cnod

  COMPLEX(kind=8),DIMENSION(nbnod),INTENT(INOUT):: vect_RHS
  
  !Local Variables
  INTEGER(kind=4):: ind_aff
  INTEGER(kind=4):: count1, count2
  INTEGER(kind=4):: dim_son1, dim_son2

  INTEGER,DIMENSION(:),ALLOCATABLE:: son1,son2

  REAL(kind=8):: min_bd, max_bd
  
  LOGICAL,DIMENSION(last_ind_clus-first_ind_clus+1):: affiliation

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Extremes of the bounding box in the direction of maximum expansion

  min_bd=MINVAL(cnod(max_dir,first_ind_clus:last_ind_clus))
  max_bd=MAXVAL(cnod(max_dir,first_ind_clus:last_ind_clus))

  !Check in which sons the nodes are

  DO ind_aff=1,last_ind_clus-first_ind_clus+1
     IF(cnod(max_dir,first_ind_clus+ind_aff-1).GE.&
          (min_bd+(max_bd-min_bd)/2)) THEN
        affiliation(ind_aff)=.TRUE.
     ELSE
        affiliation(ind_aff)=.FALSE.
     END IF
  END DO

  !Find the dimension of the two sons
  
  dim_son1=COUNT(affiliation)
  dim_son2=(last_ind_clus-first_ind_clus+1)-dim_son1
  
  ALLOCATE(son1(dim_son1),son2(dim_son2))

  !Determine the cluster affiliation
  
  count1=1
  count2=1
  
  DO ind_aff=1,last_ind_clus-first_ind_clus+1
     IF(affiliation(ind_aff).EQV.(.TRUE.)) THEN !node in the first son
        son1(count1)=first_ind_clus+ind_aff-1
        count1=count1+1
     ELSE !node in the second son
        son2(count2)=first_ind_clus+ind_aff-1
        count2=count2+1
     END IF
  END DO

  SELECT CASE (flag)
     
  CASE(.TRUE.)

     !Ordering from the minimum dof to the maximum dof
     CALL ordering_vector(son1,dim_son1)
     CALL ordering_vector(son2,dim_son2)
 
     !Reordering of the unknowns and the RHS vector
     CALL reordering_unknowns(son1,son2,dim_son1,dim_son2,ind_sons,&
          cnod,vect_RHS,nbnod,mat_corr)

     RETURN
     
  CASE(.FALSE.)

     IF(son1(1).LT.son2(1)) THEN !the smallest index is in son1
        ind_sons(1)=son1(1)
        ind_sons(2)=son1(dim_son1)
        ind_sons(3)=son2(1)
        ind_sons(4)=son2(dim_son2)
     ELSE !the smallest index is in son2
        ind_sons(1)=son2(1)
        ind_sons(2)=son2(dim_son2)
        ind_sons(3)=son1(1)
        ind_sons(4)=son1(dim_son1)
     END IF
     
     RETURN

  END SELECT
  
END SUBROUTINE compute_cluster_affiliation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE ordering_vector(vect,len_vect)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: ordering_vector                                   
!                                                                    
! INPUTS:         len_vect    : length of the vector to order        
!                                                                    
! INPUT/OUTPUTS:  vect           : vector to order                   
!                                                                    
! GOAL: order a vector from minimum to maximum                       
!                                                                    
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/5/2014                                                   
!                                                                    
! MODIFIED: VIII/28/2015 (Luca Desiderio)                  
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Input
  INTEGER(kind=4),INTENT(IN):: len_vect
  
  !Input/Output
  INTEGER(kind=4),DIMENSION(len_vect),INTENT(INOUT):: vect

  !Local Variables
  INTEGER(kind=4):: i, j
  INTEGER(kind=4):: temp

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  DO i=1,len_vect
     !REMARK: the loop over i is not redundant because it is necessary
     !to move the minimum at the beginning of the array vect.
     DO j=1,len_vect-1
        IF(vect(j)>vect(j+1)) THEN
           temp=vect(j)
           vect(j)=vect(j+1)
           vect(j+1)=temp
        ELSE IF (vect(j).EQ.vect(j+1)) THEN
           print*, 'Error: vect(j)=vect(j+1)'
           STOP
        END IF
     END DO
  END DO

  RETURN
  
END SUBROUTINE ordering_vector
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE reordering_unknowns(son1,son2,dim_son1,dim_son2,ind_sons,&
     cnod,vect_RHS,nbnod,mat_corr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: reordering_unknowns                               
!                                                                    
! INPUTS: son1       : first son                                     
!         son2       : second son                                    
!         dim_son1   : dimension of the first son                    
!         dim_son2   : dimension of the second son                   
!         nbnod      : number of the nodes in the whole problem      
!                                                                    
! INPUT/OUTPUTS: ind_sons   : array where the first and the last     
!                             nodes of each son are stored           
!                 cnod       : coordinates of the nodes              
!                 vect_RHS   : RHS of the problem                    
!                                                                    
! GOAL: reorder the unknowns and the RHS vector                      
!                                                                    
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/5/2014                                                   
!                                                                    
! MODIFIED: VIII/28/2015 (Luca Desiderio)                                              
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: nbnod
  INTEGER(kind=4),INTENT(IN):: dim_son1, dim_son2

  INTEGER(kind=4),DIMENSION(dim_son1),INTENT(INOUT):: son1
  INTEGER(kind=4),DIMENSION(dim_son2),INTENT(INOUT):: son2

  !Input/Output

  INTEGER(kind=4),DIMENSION(4),INTENT(INOUT):: ind_sons
  INTEGER(kind=4),DIMENSION(nbnod),INTENT(INOUT)::mat_corr
  
  REAL(kind=8),DIMENSION(2,nbnod),INTENT(INOUT):: cnod

  COMPLEX(kind=8),DIMENSION(nbnod),INTENT(INOUT):: vect_RHS
  
  !Local Variables
  INTEGER(kind=4):: ind
  
  INTEGER(kind=4),DIMENSION(nbnod)::copy_mat_corr

  REAL(kind=8),DIMENSION(2,nbnod):: copy_cnod

  COMPLEX(kind=8),DIMENSION(nbnod):: copy_vect_RHS

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  copy_mat_corr=mat_corr
  
  !Copy of the matrix containing the coordinates of the nodes
  copy_cnod=cnod

  !Copy of the array containing the RHS vector
  copy_vect_RHS=vect_RHS

  IF(son1(1).LE.son2(1)) THEN !the smallest index is in son1

   mat_corr(son1(1):son1(1)+dim_son1-1)=copy_mat_corr(son1)
   mat_corr(son1(1)+dim_son1:son1(1)+dim_son1+dim_son2-1)=&
        copy_mat_corr(son2)
     
   !Reorder the coordinates of the unknowns  
   cnod(:,son1(1):son1(1)+dim_son1-1)=copy_cnod(:,son1)
   cnod(:,son1(1)+dim_son1:son1(1)+dim_son1+dim_son2-1)=&
        copy_cnod(:,son2)

   !Reordering the RHS vector (a modifier)
  
   vect_RHS(son1(1):son1(1)+dim_son1-1)=copy_vect_RHS(son1)
   vect_RHS(son1(1)+dim_son1:son1(1)+dim_son1+dim_son2-1)=&
        copy_vect_RHS(son2)  

   !Extrapolation of the informations about the sons
   son1=(/(ind,ind=son1(1),son1(1)+dim_son1-1)/)
   son2=(/(ind,ind=son1(1)+dim_son1,son1(1)+dim_son1+dim_son2-1)/)

   ind_sons=(/son1(1),son1(dim_son1),son2(1),son2(dim_son2)/)
     
  ELSE !the smallest index is in son2

   mat_corr(son2(1):son2(1)+dim_son2-1)=copy_mat_corr(son2)
   mat_corr(son2(1)+dim_son2:son2(1)+dim_son1+dim_son2-1)=&
        copy_mat_corr(son1)
     
   !Reorder the unknowns
   cnod(:,son2(1):son2(1)+dim_son2-1)=copy_cnod(:,son2)
   cnod(:,son2(1)+dim_son2:son2(1)+dim_son1+dim_son2-1)=&
        copy_cnod(:,son1)

   !Reordering the RHS vector (a modifier)
 
   vect_RHS(son2(1):son2(1)+dim_son2-1)=copy_vect_RHS(son2)
   vect_RHS(son2(1)+dim_son2:son2(1)+dim_son1+dim_son2-1)=&
        copy_vect_RHS(son1)

   !Extrapolation of the informations about the sons
   son2=(/(ind,ind=son2(1),son2(1)+dim_son2-1)/)
   son1=(/(ind,ind=son2(1)+dim_son2,son2(1)+dim_son1+dim_son2-1)/)

   ind_sons=(/son2(1),son2(dim_son2),son1(1),son1(dim_son1)/)
     
  END IF

  RETURN

END SUBROUTINE reordering_unknowns
end module cluster_tree
