module initialisation_hmat
use variables, only: tangent, normal, node_left, node_right, node, nb_element,element
use variables, only: normal_right, normal_left, fracture_mode
use cluster_tree
use block_cluster_tree
use H_matrix_construction
use H_matrix_recompression
implicit none
private
public :: initiate_hmat, backward_reorder
 
contains
subroutine initiate_hmat(nbnod,cnod2,nbclus,ind_clusters,       &
                nbblclus, blcluster_tree,info_blclusters,              &
                dimA,dimB,dimF, A_matrices,B_matrices,F_matrices,      &
                eta, nleaf,eps_ACA, mat_corr,eval_kernel)
                
interface
    function eval_kernel(ind_row,ind_col,Ynode,cnod_row,nbnod_row,id_row,id_col)
        ! Input
        REAL(kind=8),DIMENSION(2)              :: Ynode
        INTEGER                                :: ind_row,ind_col,nbnod_row, id_row,id_col
        REAL(kind=8),DIMENSION(2,nbnod_row)    :: cnod_row
        !Output
        COMPLEX(kind=8):: eval_kernel
    end function
end interface
      
                
integer                                              :: nbnod, nbclus, nbblclus 
integer(kind=8)                                      :: dimA, dimB, dimF  
complex(kind=8), dimension(:), allocatable, intent(INOUT)      :: A_matrices
complex(kind=8), dimension(:), allocatable, intent(INOUT)      :: B_matrices
complex(kind=8), dimension(:), allocatable, intent(INOUT)      :: F_matrices              
real(kind=8)   , dimension(2,nbnod)                            :: cnod2

integer(kind=4), dimension(:), allocatable, intent(inout):: mat_corr
integer(kind=4), dimension(:), allocatable, intent(inout):: ind_clusters
integer(kind=4), dimension(:), allocatable               :: sons_clusters  
integer(kind=8),dimension(:,:), allocatable,intent(inout):: info_blclusters
integer(kind=4),dimension(:,:),allocatable               :: blcluster_tree
integer(kind=8)                                          :: stor_lowrank
 
real(kind=8)                                          ::  eps_ACA
real(kind=8)                                          ::  eta
integer(kind=4)                                       ::  nleaf


! Internal variable
integer                                               :: k
real												:: comp_rate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 BUILD the CLUSTER TREE                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call build_cluster_tree(nbnod,nbclus,ind_clusters,sons_clusters,&
       cnod2,nleaf,mat_corr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            REORDER IN H MATRIX FORMAT                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

call forward_reorder(nb_element,node_left,mat_corr)   
call forward_reorder(nb_element,node_right,mat_corr)  
call forward_reorder(nb_element,tangent,mat_corr)  
call forward_reorder(nb_element,normal,mat_corr)  
if (fracture_mode=='modeII') then    
    call forward_reorder(nb_element,normal_left,mat_corr)  
    call forward_reorder(nb_element,normal_right,mat_corr)  
end if  
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            BUILD the BLOCK CLUSTER TREE                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call build_block_cluster_tree(nbclus,ind_clusters,sons_clusters,&
       nbnod,cnod2,nbblclus,blcluster_tree,info_blclusters,eta)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              H-MATRIX REPRESENTATION                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call build_H_matrix(nbnod,cnod2,nbclus,ind_clusters,nbblclus,&
       blcluster_tree,info_blclusters,&
       stor_lowrank,eps_ACA,dimA,A_matrices,dimB,B_matrices,&
       dimF,F_matrices,eval_kernel)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              H-MATRIX RECOMPRESSION                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call compress_H_matrix(nbclus,ind_clusters,nbblclus,&
       blcluster_tree,info_blclusters,dimA,A_matrices,dimB,B_matrices,&
       eps_ACA,stor_lowrank)
comp_rate=REAL(stor_lowrank+dimF)/REAL(nbnod) 
   comp_rate=comp_rate /REAL(nbnod)

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Reorder in forward manner (reorder ini )
!
subroutine forward_reorder(nb,list,mat_corr)
integer                            :: nb
real(kind=8), dimension(2,nb)      :: list, list_temp
integer(kind=4),dimension(nb)      :: mat_corr
!
integer                            :: k
!
! Initiate temp_list
list_temp = list
! Reordering
do k=1,nb
    list(:,k) = list_temp(:,mat_corr(k))
end do
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
subroutine backward_reorder(nb,list,mat_corr)
integer                            :: nb
real(kind=8), dimension(2,nb)      :: list, list_temp
integer(kind=4),dimension(nb)      :: mat_corr
!
integer                            :: k
!
! Initiate temp_list
list_temp = list
! Reordering 
do k=1,nb
    list(:,mat_corr(k)) = list_temp(:,k)
end do
end subroutine
end module
