module H_matrix_construction
use ACA_algorithm
use get_info_matrices
implicit none
private
public :: build_H_matrix



contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
subroutine build_H_matrix(nbnod,cnod2,nbclus,ind_clusters,nbblclus,&
       blcluster_tree,info_blclusters,&
       stor_lowrank,eps_ACA,dimA,A_matrices,dimB,B_matrices,&
       dimF,F_matrices,eval_kernel)
                
interface
    function eval_kernel(ind_row,ind_col,Ynode,cnod_row,nbnod_row,clus_row,clus_col)
        ! Input
        REAL(kind=8),DIMENSION(2)           :: Ynode
        INTEGER                             :: ind_row,ind_col,nbnod_row, clus_row,clus_col
        REAL(kind=8),DIMENSION(2,nbnod_row) :: cnod_row
        !Output
        COMPLEX(kind=8):: eval_kernel
    end function
end interface         
       

  
  ! Variables
  integer                                                   :: nbnod, nbclus, nbblclus
  real(kind=8)  , dimension(2,nbnod)                        :: cnod2
  integer(kind=8)                                           :: dimA, dimB, dimF   
  complex(kind=8), dimension(:), allocatable, intent(INOUT) :: A_matrices
  complex(kind=8), dimension(:), allocatable, intent(INOUT) :: B_matrices
  complex(kind=8), dimension(:), allocatable, intent(INOUT) :: F_matrices   
  integer(kind=8),dimension(nbblclus,6),intent(inout)       :: info_blclusters
  integer(kind=4), dimension(2*nbclus), intent(inout)       :: ind_clusters
  integer(kind=4),dimension(:,:),allocatable                :: blcluster_tree
  integer(kind=8),intent(inout)                             :: stor_lowrank
  real(kind=8)                                              ::  eps_ACA
  
  
  
  ! Local variables
  integer(kind=4),dimension(:),allocatable              :: flag_ACA 
  integer(kind=4)                                       :: AllocateStatus
  real                                                  :: comp_rate
  integer(kind=4)                                       :: max_low_rank 


  !Estimation of the maximum low-rank (linear regression)
  max_low_rank=INT(1.1685d0*nbnod**0.516242d0,4) ! A MODIFIER PAR PIERRE
 
  max_low_rank=max_low_rank+INT(2*max_low_rank,4)
  print*,'max low rank before',max_low_rank
  max_low_rank = INT(20,4)
  print*,'max low rank after',max_low_rank
  !Allocate array flag_ACA
  ALLOCATE(flag_ACA(nbblclus))
  !Initialize array flag_ACA (vectorial partial-pivot ACA)
  flag_ACA=INT(1,4)
  
  print*,'Starting compression of the matrix'
   
  !****Step 1: PARAMETERS of the H-MATRIX
  dimA=INT(0,8) !Inizialize the dimension of array A_matrices
  dimB=INT(0,8) !Inizialize the dimension of array B_matrices
  dimF=INT(0,8) !Inizialize the dimension of array F_matrices

  !Initialize the value of variables stor_lowrank
  stor_lowrank=INT(0,8)

  print*,'parameters_H_matrix'
  print*,'eps_ACA',eps_ACA
  !Estimation of the dimensions of the arrays to store the matrices A,
  !B for low rank approximations (M=A*B**(T)) of admissible leaf-blocks
  CALL parameters_H_matrix(nbnod,cnod2,nbclus,ind_clusters,nbblclus,&
       blcluster_tree,info_blclusters,dimA,dimB,dimF,max_low_rank,&
       stor_lowrank,eps_ACA,eval_kernel)
  
  !Compute the new maximum low-rank
  max_low_rank=MAXVAL(INT(info_blclusters(:,5),4))
  print*,'max_low_rank',max_low_rank 
  !Compute the compression rate
  comp_rate=REAL(stor_lowrank+dimF,4)/REAL(nbnod,4)
  comp_rate=comp_rate/REAL(nbnod,4)
  print*,'Compression rate =',comp_rate
  
  !****Step 2: COMPUTATION of the H-MATRIX
  
  ALLOCATE(A_matrices(dimA),B_matrices(dimB),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory A and B ***"
  
  ALLOCATE(F_matrices(dimF),STAT=AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory F ***"

  A_matrices=DCMPLX(0.d0) !Initialization of the array A_matrices
  B_matrices=DCMPLX(0.d0) !Initialization of the array B_matrices
  F_matrices=DCMPLX(0.d0) !Initialization of the array F_matrices
  print*,'compute_H_matrix'
  !Computation of the admissible and non-admissible blocks
  CALL compute_H_matrix(nbnod,cnod2,nbclus,ind_clusters,nbblclus,&
       blcluster_tree,info_blclusters,dimA,A_matrices,dimB,B_matrices,&
       dimF,F_matrices,eps_ACA,eval_kernel)
 

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE parameters_H_matrix(nbnod,cnod,nbclus,ind_clusters,nbblclus,&
     blcluster_tree,info_blclusters,dimA,dimB,dimF,max_low_rank,&
     stor_lowrank,eps_ACA,eval_kernel)
                
interface
    function eval_kernel(ind_row,ind_col,Ynode,cnod_row,nbnod_row,clus_row,clus_col)
        ! Input
        REAL(kind=8),DIMENSION(2)               :: Ynode
        INTEGER                                 :: ind_row,ind_col,nbnod_row, clus_row,clus_col
        REAL(kind=8),DIMENSION(2,nbnod_row)     :: cnod_row
        !Output
        COMPLEX(kind=8):: eval_kernel
    end function
end interface         
 
 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: parameters_H_matrix
!                                                                    
! INPUTS:  blcluster_tree : matrix storing the informations on Block 
!                           Cluster Tree                             
!          nbblclus       : number of block clusters                 
!          ind_clusters   : array storing the cluster tree           
!          nbclus         : number of clusters in the tree           
!          max_low_rank   : maximum of the ranks of low-rank
!                           approximations (fixed parameter)
!          ind_start      : index to start the investigation of the
!                           block clusters  
!                                                                    
! INPUT/OUTPUTS: dimA       : expected dimension of the array to store
!                             the matrices A
!                dimB       : expected dimension of the array to store
!                             the matrices B
!                                                                    
! GOAL: evaluate the dimension of the arrays to store the matrices A 
!       and B that realize the low-rank representation of admissible  
!       leaf-blocks (M=A*B**(T)).               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: nbnod, nbclus, nbblclus 
  INTEGER(kind=4),INTENT(IN):: max_low_rank 
  
  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(IN):: blcluster_tree
 
  REAL(kind=8),INTENT(IN):: eps_ACA
 
  REAL(kind=8),DIMENSION(2,nbnod),INTENT(IN):: cnod
  
  !Input/Output
  INTEGER(kind=8),INTENT(INOUT):: dimA, dimB, dimF, stor_lowrank
 
  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(INOUT):: info_blclusters

  !Local Variables
  INTEGER(kind=4):: ind_rec, nbrow, nbcol, clus_row, clus_col, id_row, id_col
  INTEGER(kind=4):: AllocateStatus 
  logical        :: flag_ACA
  INTEGER(kind=4):: low_rank
  INTEGER(kind=8):: dimA_block, dimB_block

  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:):: A_matrix, B_matrix


 !PIERRE: A et B peuvent etre reels
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

 ! print*,nbblclus
  
  !Loop on the block-clusters
  DO ind_rec=1,nbblclus
    ! print*,'ind_block',ind_rec
     !Check if the current block is an admissible leaf-block
     IF(blcluster_tree(ind_rec,1).EQ.1) THEN
     
        !Cluster row and cluster column of the current block
        clus_row=2*blcluster_tree(ind_rec,2)
        clus_col=2*blcluster_tree(ind_rec,3)
 
        !Number of rows of the current block
        nbrow=cmp_nbrow(ind_rec,nbclus,ind_clusters,nbblclus,&
             blcluster_tree)
        
        !Number of columns of the current block
        nbcol=cmp_nbcol(ind_rec,nbclus,ind_clusters,nbblclus,&
             blcluster_tree)

        !Compute the maximum memory requirement for matrix A
        dimA_block=INT(nbrow*max_low_rank,8)

        !Compute the maximum memory requirement for matrix B
        dimB_block=INT(nbcol*max_low_rank,8)

        !Allocate arrays A_matrix and B_matrix
        ALLOCATE(A_matrix(dimA_block),STAT=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory A_matrix ***"
        ALLOCATE(B_matrix(dimB_block),STAT=AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory B_matrix ***"

        !Initialize the value of low_rank
        low_rank=INT(0,4)

       
        !Initialize arrays A_matrix and B_matrix
        A_matrix=DCMPLX(0.d0)
        B_matrix=DCMPLX(0.d0)
 
        id_row = ind_clusters(clus_row-1)
        id_col = ind_clusters(clus_col-1)

        CALL compute_full_ACA(cnod(:,ind_clusters(clus_row-1):&
             ind_clusters(clus_row)),nbrow,&  
             cnod(:,ind_clusters(clus_col-1):&
             ind_clusters(clus_col)),nbcol,&  
             A_matrix,dimA_block,B_matrix,dimB_block,eps_ACA,&
             low_rank,INT(max_low_rank,4),flag_ACA,ind_rec,id_row,id_col,eval_kernel)

        

        !Deallocate arrays A_matrix and B_matrix
        DEALLOCATE(A_matrix,B_matrix)

        stor_lowrank=stor_lowrank+INT((nbrow+nbcol)*low_rank,8)

        !Update the values of dimA and dimB
        dimA=dimA+INT(4.0d0*nbrow*low_rank,8)
        dimB=dimB+INT(4.0d0*nbcol*low_rank,8)

        !Store the informations about the low-rank of the current 
        !block
        info_blclusters(ind_rec,5)=INT(low_rank,8)
        
        !print*,ind_rec,low_rank,blcluster_tree(ind_rec,1),nbrow,nbcol
        !pause
        !Store the informations about the maximum low-rank of the  
        !current block
        info_blclusters(ind_rec,6)=INT(low_rank,8)

     ELSE IF (blcluster_tree(ind_rec,1).EQ.0) THEN !Full rank

        !Number of rows in the matrix
        nbrow=cmp_nbrow(ind_rec,nbclus,ind_clusters,nbblclus,&
             blcluster_tree)
        
        !Number of columns in the matrix
        nbcol=cmp_nbcol(ind_rec,nbclus,ind_clusters,nbblclus,&
             blcluster_tree)

        !Update the value of dimF
        dimF=dimF+INT(nbrow*nbcol,8)

     END IF
     
  END DO
   
  RETURN

END SUBROUTINE parameters_H_matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE compute_H_matrix(nbnod,cnod,nbclus,ind_clusters,nbblclus,&
     blcluster_tree,info_blclusters,dimA,A_matrices,dimB,B_matrices,&
     dimF,F_matrices,eps_ACA,eval_kernel)
                
interface
    function eval_kernel(ind_row,ind_col,Ynode,cnod_row,nbnod_row,clus_row,clus_col)
        ! Input
        REAL(kind=8),DIMENSION(2)               :: Ynode
        INTEGER                                 :: ind_row,ind_col,nbnod_row, clus_row,clus_col
        REAL(kind=8),DIMENSION(2,nbnod_row)     :: cnod_row
        !Output
        COMPLEX(kind=8):: eval_kernel
    end function
end interface       

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: parameters_H_matrix
!                                                                    
! INPUTS:  blcluster_tree : matrix storing the informations on Block 
!                           Cluster Tree                             
!          nbblclus       : number of block clusters                 
!          ind_clusters   : array storing the cluster tree           
!          nbclus         : number of clusters in the tree           
!          max_low_rank   : maximum of the ranks of low-rank
!                           approximations (fixed parameter)
!          ind_start      : index to start the investigation of the
!                           block clusters  
!                                                                    
! INPUT/OUTPUTS: dimA       : expected dimension of the array to store
!                             the matrices A
!                dimB       : expected dimension of the array to store
!                             the matrices B
!                                                                    
! GOAL: evaluate the dimension of the arrays to store the matrices A 
!       and B that realize the low-rank representation of admissible  
!       leaf-blocks (M=A*B**(T)).
                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: nbnod, nbclus, nbblclus 
  
  INTEGER(kind=8),INTENT(IN):: dimA, dimB, dimF

 
  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(IN):: blcluster_tree
 
  REAL(kind=8),INTENT(IN):: eps_ACA

  REAL(kind=8),DIMENSION(2,nbnod),INTENT(IN):: cnod
  
  !Input/Output
  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(INOUT):: info_blclusters

  COMPLEX(kind=8),DIMENSION(dimA),INTENT(INOUT):: A_matrices
  COMPLEX(kind=8),DIMENSION(dimB),INTENT(INOUT):: B_matrices
  COMPLEX(kind=8),DIMENSION(dimF),INTENT(INOUT):: F_matrices

  !Local Variables
  INTEGER(kind=4)                              :: ind_block 
  INTEGER(kind=8)                              :: posA, posB, posF
  logical                                      :: flag_ACA
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  posA=INT(1,8)!Inizialize the values of posA 
  posB=INT(1,8)!Inizialize the values of posB
  posF=INT(0,8)!Inizialize the values of posF 

 
  !Loop on the block-clusters
  DO ind_block=1,nbblclus
     
     !Check if the current block is an admissible leaf-block
     IF(blcluster_tree(ind_block,1).EQ.1) THEN
 
        CALL compute_low_rank_approx(nbnod,cnod,nbclus,ind_clusters,&
             nbblclus,blcluster_tree,info_blclusters,ind_block,dimA,&
             A_matrices,posA,dimB,B_matrices,posB,eps_ACA,&
             flag_ACA,eval_kernel)
 
     ELSE IF (blcluster_tree(ind_block,1).EQ.0) THEN

        CALL compute_full_matrices(nbnod,cnod,nbclus,ind_clusters,&
             nbblclus,blcluster_tree,info_blclusters,ind_block,dimF,&
             F_matrices,posF,eval_kernel)
 
     END IF
     
  END DO
  
  RETURN

END SUBROUTINE compute_H_matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE compute_low_rank_approx(nbnod,cnod,nbclus,ind_clusters,&
     nbblclus,blcluster_tree,info_blclusters,ind_block,&
     dimA,A_matrices,pos_A,dimB,B_matrices,pos_B,eps_ACA,&
     flag_ACA,eval_kernel)
                
interface
    function eval_kernel(ind_row,ind_col,Ynode,cnod_row,nbnod_row,id_row,id_col)
        ! Input
        REAL(kind=8),DIMENSION(2)            :: Ynode
        INTEGER                              :: ind_row,ind_col,nbnod_row, id_row,id_col
        REAL(kind=8),DIMENSION(2,nbnod_row)  :: cnod_row
        !Output
        COMPLEX(kind=8):: eval_kernel
    end function
end interface         
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: compute_low_rank_approx
!                                                                    
! INPUTS:  nbnod          : number of nodes 
!          nbclus         : number of clusters
!          nbblclus       : number of block clusters                 
!          dimA           : dimension of array A_matrices           
!          dimB           : dimension of array B_matrices           
!          cnod           : coordinates of the nodes
!          ind_clusters   : array storing the cluster tree
!          blcluster_tree : matrix storing the block cluster tree
!          eps_ACA        : ACA threshold
!
! OUTPUT:  flag_memory    : flag to check if the memory allocated is
!                           sufficient  
!                                                                    
! INPUT/OUTPUTS: ind_cur         : index of the current block cluster
!                max_low_rank    : mamimun low-rank of the
!                                  approximations
!                pos_A           : starting position to store
!                                  A matrices
!                pos_B           : starting position to store
!                                  B matrices
!                info_blclusters : informations about block clusters
!                A_matrices      : array storing A matrices
!                B_matrices      : array storing B matrices
!                                                                    
! GOAL: evaluate the dimension of the arrays to store the matrices A 
!       and B that realize the low-rank representation of admissible  
!       leaf-blocks (M=A*B**(T)).
!
! REMARK: Each row of the array info_blclusters correspons to a block
!         cluster. If the block cluster is a SPLITTED BLOCK, the
!         corresponding row of info_blclusters has the following
!         structure:
!
!         * * @ @ 0 
!
!         where * * are the sons of the row cluster and @ @ are the 
!         sons of column cluster.
!  
!         If the bluck cluster is an ADMISSIBLE LEAF BLOCK, the
!         corresponding row of info_blclusters has the following
!         structure:
!
!         *s *e @s @e k
!
!         where *s and *e represent the starting and ending position
!         of matrix A in array A_matrices, @s and @e represent the
!         starting and ending position of matrix B in array B_matrices
!         and k the low-rank of the block.
!
!         If the bluck cluster is a NON-ADMISSIBLE LEAF BLOCK, the
!         corresponding row of info_blclusters has the following
!         structure:
!
!         *s *e 0 0 0
!
!         where *s and *e represent the starting and ending position
!         of the block in array full_matrices.                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: nbnod, nbclus, nbblclus ,ind_block
  INTEGER(kind=4)           :: id_col, id_row
  INTEGER(kind=8),INTENT(IN):: dimA, dimB
 
  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(IN):: blcluster_tree
 
  REAL(kind=8),INTENT(IN):: eps_ACA

  REAL(kind=8),DIMENSION(2,nbnod),INTENT(IN):: cnod
 

  !Input/Output
  INTEGER(kind=8),INTENT(INOUT):: pos_A, pos_B
 
  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(INOUT):: info_blclusters
  
  COMPLEX(kind=8),DIMENSION(dimA),INTENT(INOUT):: A_matrices
  COMPLEX(kind=8),DIMENSION(dimB),INTENT(INOUT):: B_matrices

  !Local Variables
  INTEGER(kind=4):: clus_row, clus_col, nbrow, nbcol
  INTEGER(kind=4):: rank, low_rank
  INTEGER(kind=4):: ind_rec
  INTEGER(kind=8):: dimA_block, dimB_block

  !Input/Output
  LOGICAL,INTENT(INOUT):: flag_ACA
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  ind_rec=0
  
  !Cluster row and cluster column of the current block
  clus_row=2*blcluster_tree(ind_block,2)
  clus_col=2*blcluster_tree(ind_block,3)

  !Number of rows of the current block
  nbrow=cmp_nbrow(ind_block,nbclus,ind_clusters,nbblclus,blcluster_tree)
        
  !Number of columns of the current block
  nbcol=cmp_nbcol(ind_block,nbclus,ind_clusters,nbblclus,blcluster_tree)

  !Low-rank of the current block
  rank=cmp_rank(ind_block,nbblclus,info_blclusters)

  !print*,'rank',rank,ind_block
  
  !Compute the memory requirement for the low_rank approximation of the
  !current block
  dimA_block=INT(nbrow*rank,8)
  dimB_block=INT(nbcol*rank,8)

  !Initialize the value of low_rank (necessary to call subroutine
  !ACA_3D_vect)
  low_rank=INT(0,4)
  
  ! Calculate id_row, id_col
  id_row = ind_clusters(clus_row-1)
  id_col = ind_clusters(clus_col-1)
  
  
  !Adaptive Cross Approximation of the current block
  CALL compute_full_ACA(cnod(:,ind_clusters(clus_row-1):&
       ind_clusters(clus_row)),nbrow,&  
       cnod(:,ind_clusters(clus_col-1):&
       ind_clusters(clus_col)),nbcol,& 
       A_matrices(pos_A:pos_A+dimA_block-1),dimA_block,&
       B_matrices(pos_B:pos_B+dimB_block-1),dimB_block,&
       eps_ACA,low_rank,INT(rank,4),flag_ACA,ind_rec,id_row,id_col,eval_kernel)
  !print*,'after',low_rank
  !Store the informations about matrix A in info_blclusters
  info_blclusters(ind_block,1)=pos_A
  info_blclusters(ind_block,2)=pos_A+INT(nbrow*low_rank-1,8)

  !Store the informations about matrix B in info_blclusters
  info_blclusters(ind_block,3)=pos_B
  info_blclusters(ind_block,4)=pos_B+INT(nbcol*low_rank-1,8)
       
  !Update the value of pos_A and pos_B
  pos_A=pos_A+INT(4.0d0*nbrow*low_rank,8)
  pos_B=pos_B+INT(4.0d0*nbcol*low_rank,8)

  RETURN

END SUBROUTINE compute_low_rank_approx
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE compute_full_matrices(nbnod,cnod,nbclus,ind_clusters,&
     nbblclus,blcluster_tree,info_blclusters,ind_block,dimF,&
     F_matrices,pos_F,eval_kernel)
                
interface
    function eval_kernel(ind_row,ind_col,Ynode,cnod_row,nbnod_row,id_row,id_col)
        ! Input
        REAL(kind=8),DIMENSION(2)            :: Ynode
        INTEGER                              :: ind_row,ind_col,nbnod_row, id_row,id_col
        REAL(kind=8),DIMENSION(2,nbnod_row)  :: cnod_row
        !Output
        COMPLEX(kind=8):: eval_kernel
    end function
end interface         
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: compute_full_matrices
!                                                                    
! INPUTS:  nbnod          : number of nodes 
!          nbclus         : number of clusters
!          nbblclus       : number of block clusters                 
!          dimF           : dimension of array F_matrices           
!          cnod           : coordinates of the nodes
!          ind_clusters   : array storing the cluster tree
!          blcluster_tree : matrix storing the block cluster tree
!          wave_num       : wave number
!          alpha          : parameter to avoid the singularity of
!                           Helmholtz kerne
!
! INPUT/OUTPUTS: info_blclusters : informations about block clusters
!                F_matrices      : array storing the full matrices
!                                                                    
! GOAL: compute the full matrices to approximate the non-admissible
!       leaf blocks  
!
! REMARK: Each row of the array info_blclusters correspons to a block
!         cluster. If the block cluster is a SPLITTED BLOCK, the
!         corresponding row of info_blclusters has the following
!         structure:
!
!         * * @ @ 0 
!
!         where * * are the sons of the row cluster and @ @ are the 
!         sons of column cluster.
!  
!         If the bluck cluster is an ADMISSIBLE LEAF BLOCK, the
!         corresponding row of info_blclusters has the following
!         structure:
!
!         *s *e @s @e k
!
!         where *s and *e represent the starting and ending position
!         of matrix A in array A_matrices, @s and @e represent the
!         starting and ending position of matrix B in array B_matrices
!         and k the low-rank of the block.
!
!         If the bluck cluster is a NON-ADMISSIBLE LEAF BLOCK, the
!         corresponding row of info_blclusters has the following
!         structure:
!
!         *s *e 0 0 0
!
!         where *s and *e represent the starting and ending position
!         of the block in array full_matrices.  
!                                                                    
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: I/13/2015                                                  
!                                                                    
! MODIFIED: Stephanie Chaillat (VII/09/2015)
!           Luca Desiderio (II/9/2016)  
!                                                                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: nbnod, nbclus, nbblclus 
  INTEGER(kind=4),INTENT(IN)::  ind_block
  
  INTEGER(kind=8),INTENT(IN):: dimF

  INTEGER(kind=4),DIMENSION(2*nbclus),INTENT(IN):: ind_clusters
  INTEGER(kind=4),DIMENSION(nbblclus,3),INTENT(IN):: blcluster_tree
 
  REAL(kind=8),DIMENSION(2,nbnod),INTENT(IN):: cnod
   
  !Input/Output
  INTEGER(kind=8),INTENT(INOUT):: pos_F
  
  INTEGER(kind=8),DIMENSION(nbblclus,6),INTENT(INOUT):: info_blclusters
  
  COMPLEX(kind=8),DIMENSION(dimF),INTENT(INOUT):: F_matrices

  !Local Variables
  INTEGER(kind=4):: ind_row, ind_col, id_row, id_col
  INTEGER(kind=4):: clus_row, clus_col, nbnod_row, nbnod_col

  REAL(kind=8),DIMENSION(2):: Y 
  REAL(kind=8),DIMENSION(:,:),ALLOCATABLE::cnod_row
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Cluster row and cluster column of the current block
  clus_row=2*blcluster_tree(ind_block,2)
  clus_col=2*blcluster_tree(ind_block,3)
  
  !Number of rows of the current block
  nbnod_row=ind_clusters(clus_row)-ind_clusters(clus_row-1)+1

  !Number of columns of the current block
  nbnod_col=ind_clusters(clus_col)-ind_clusters(clus_col-1)+1

  !Store the starting position of the current full matrix
  info_blclusters(ind_block,1)=pos_F+INT(1,8)

  !The current non-admissible leaf is stored col-wise in array
  !F_matrices

  ALLOCATE(cnod_row(2,nbnod_row))
  cnod_row(1,1:nbnod_row)=cnod(1,ind_clusters(clus_row-1)-1+1:ind_clusters(clus_row-1)-1+nbnod_row)
  cnod_row(2,1:nbnod_row)=cnod(2,ind_clusters(clus_row-1)-1+1:ind_clusters(clus_row-1)-1+nbnod_row)
  
  ! Calculate id_row, id _col
  id_row = ind_clusters(clus_row-1)
  id_col = ind_clusters(clus_col-1)
  
  DO ind_col=1,nbnod_col
     DO ind_row=1,nbnod_row
 
        Y(1:2)=cnod(1:2,ind_clusters(clus_col-1)-1+ind_col)

        
        F_matrices(pos_F+INT((ind_col-1)*nbnod_row+&
             (ind_row-1)+1,8))=eval_kernel(ind_row,ind_col,Y,cnod_row,nbnod_row,id_row,id_col)
 
     END DO
  END DO
  !Update the value of pos_F
  pos_F=pos_F+INT(nbnod_row*nbnod_col,8)

  !Store the ending position of the current full matrix
  info_blclusters(ind_block,2)=pos_F

  RETURN
  
END SUBROUTINE compute_full_matrices
 

end module 
