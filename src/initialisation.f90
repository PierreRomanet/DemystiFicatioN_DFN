!
!
!
!
!
module initialisation
use variables
use RateStateAgingModeIII, only: ode_RateStateAgingModeIII
use RateStateAgingModeII,  only: ode_RateStateAgingModeII
use initialisation_hmat
use eval_kernel_hmatrix
implicit none
external :: directII, directIII, hmatII, hmatIII, fmmIII

!
!
!
!
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   Main subroutine 
!
subroutine initialise()
!
!   Check if output directory exist 
call mkdir_simulation_name()
!
!   copy input files in the output directory
call cp_inputFile()
!
!   Allocate all the arrays
call allocate_ini()
!
!   Create Initial values for each faults  
call ini_value()
!
!   Find node where faults are crossing
call find_fault_crossing()
 
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Make a directory which name is the simulation name
!
subroutine mkdir_simulation_name()
character(len=1000)              :: command_line
logical                         :: dir_exist
!
command_line = 'mkdir ./problems/'//trim(simulation_name)
inquire(file='./problems/'//trim(simulation_name),exist=dir_exist)
if(.NOT. dir_exist) then
    call system(command_line)
end if
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    Copy config.in, data_mesh.in, config_fault.in in the output_dir
!
subroutine cp_inputFile()
character(len=1000)   :: command_line
!
command_line = 'cp '//trim(config_file)//' '//'./problems/'//trim(simulation_name)
print*,trim(command_line)
call system(trim(command_line))
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   allocate arrays and pointer needed by fault-type and solver 
!   
subroutine allocate_ini()

! MODE II
print*,'fracture_mode',fracture_mode
if (fracture_mode == 'modeII') then
 
    ! Normal between elements (specific to mode II)
    allocate(normal_left(2,nb_element))
    allocate(normal_right(2,nb_element))
    
    
! MODE III  
elseif(fracture_mode == 'modeIII') then

else
    print*,'Please chose a fracture mode: "modeII" or "modeIII"'
    stop
end if
!
!
!   Allocate pointers 
! Pointer to fric_law
if((fric_law == 'RateStateAging').and.(fracture_mode == 'modeIII')) then
    fric_law_ptr => ode_RateStateAgingModeIII
elseif((fric_law == 'RateStateAging').and.(fracture_mode == 'modeII')) then
    fric_law_ptr => ode_RateStateAgingModeII
else
    print*,'Problem in config.in: wrong friction law, please choose between "RateStateAgeing" or "RateStateSlip"'
    print*,'Maybe the friction with the associated mode of fracture you wanted is not available yet'
    stop
end if
!
! MODE III
if((static_kernel == 'fmm').and.(fracture_mode == 'modeIII')) then
    kernel_ptr => fmmIII
elseif((static_kernel == 'direct').and.(fracture_mode == 'modeIII')) then
    kernel_ptr => directIII
elseif((static_kernel == 'hmatrix').and.(fracture_mode == 'modeIII')) then
    kernel_ptr => hmatIII
    ! Allocate nodes attach to element
    allocate(node_left(2,nb_element))
    allocate(node_right(2,nb_element))
!    
! MODE II   
elseif((static_kernel == 'direct').and.(fracture_mode == 'modeII')) then
    kernel_ptr => directII
elseif((static_kernel == 'hmatrix').and.(fracture_mode == 'modeII')) then
    kernel_ptr => hmatII
    allocate(node_left(2,nb_element))
    allocate(node_right(2,nb_element))
  
!    
! IF NOT AVAILABLE     
else
    print*,'Problem in config.in: wrong static kernel or mode of fracture'
    stop
end if



end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
subroutine ini_value()
integer                                   :: k
!
! Initialise value by reading HMAT files
!
! Calculate element, tangent, normal, background loadings
do k=1,nb_element
  element(1,k) = (node(1,k+1) + node(1,k))*0.5_dp
  element(2,k) = (node(2,k+1) + node(2,k))*0.5_dp
  ds(k) = sqrt((node(1,k+1) - node(1,k))**2+(node(2,k+1) - node(2,k))**2)
  tangent(1,k) = (node(1,k+1) - node(1,k))/ds(k)
  tangent(2,k) = (node(2,k+1) - node(2,k))/ds(k)
  normal(1,k) = -tangent(2,k)
  normal(2,k) = tangent(1,k)
end do
!
! Change variable for theta (now define over -infty to + infty)
!theta = log(theta*V0/Dc)

!
if (fracture_mode == 'modeIII') then
    ! Background shear stress rate projected on the fault
    shear_loading_dot(1:nb_element) = shear_loading_dot(1:nb_element) &
                                    - tangent(2,1:nb_element)*Sigma31_dot  &
                                    + tangent(1,1:nb_element)*Sigma32_dot

    
    ! Normal stress from a global stress state
    SigmaN(1:nb_element) = SigmaN(1:nb_element) +                                        & ! Specific normal stress distribution on a fault              
                                          normal(1,1:nb_element)*                        & ! + global stress state
     (normal(1,1:nb_element)*Sigma11_dot+normal(2,1:nb_element)*Sigma12_dot)             &
                                        + normal(2,1:nb_element)*                        &
     (normal(1,1:nb_element)*Sigma12_dot+normal(2,1:nb_element)*Sigma22_dot) 
    
elseif (fracture_mode == 'modeII') then
    ! Normal stress from a global stress state
    SigmaN(1:nb_element) = SigmaN(1:nb_element) +                                        & ! Specific normal stress distribution on a fault           
                                          normal(1,1:nb_element)*                        & ! + global stress state 
     (normal(1,1:nb_element)*Sigma11+normal(2,1:nb_element)*Sigma12)             &
                                        + normal(2,1:nb_element)*                        &
     (normal(1,1:nb_element)*Sigma12+normal(2,1:nb_element)*Sigma22) 
    ! 
    ! Background normal stress rate projected on the fault
    normal_loading_dot(1:nb_element) = normal_loading_dot(1:nb_element) &
                                        +normal(1,1:nb_element) *                       & 
     (normal(1,1:nb_element)*Sigma11_dot+normal(2,1:nb_element)*Sigma12_dot)     &
                                        + normal(2,1:nb_element)*                        &
     (normal(1,1:nb_element)*Sigma12_dot+normal(2,1:nb_element)*Sigma22_dot)      
     
    ! Background shear stress rate projected on the fault
    shear_loading_dot(1:nb_element) = shear_loading_dot(1:nb_element) &
                                    +tangent(1,1:nb_element)*                       &
     (normal(1,1:nb_element)*Sigma11_dot+normal(2,1:nb_element)*Sigma12_dot)     &
                                    + tangent(2,1:nb_element)*                       &
     (normal(1,1:nb_element)*Sigma12_dot+normal(2,1:nb_element)*Sigma22_dot) 



    
    ! Calculate normal_left and normal_right
    normal_left(:,:) = normal(:,:)
    normal_right(:,:) = normal(:,:)
    
else
    print*,'Check fracture mode'
    stop
end if

!
!
! Initialisation specific for H-matrix
if (static_kernel == 'hmatrix') then
    ! Initialise node left and node right
    node_left(:,1:nb_element)  =  node(:,1:nb_element)
    node_right(:,1:nb_element) =  node(:,2:nb_element+1)
   if (fracture_mode == 'modeIII') then
       call initiate_hmat(nb_element,element,nbclus,ind_clusters,         &
                   nbblclus, blcluster_tree,info_blclustersT,             &
                   dimA_T,dimB_T,dimF_T,                                  &
                   A_matrices_T,B_matrices_T,F_matrices_T, &
                   eta, nleaf, eps_ACA, mat_corr,eval_kernel_ModeIII)   
                
                
   elseif (fracture_mode == 'modeII') then
       call initiate_hmat(nb_element,element,nbclus,ind_clusters,         &
                   nbblclus, blcluster_tree,info_blclustersN,             &
                   dimA_N,dimB_N,dimF_N,                                  &
                   A_matrices_N,B_matrices_N,F_matrices_N,                &
                   eta, nleaf, eps_ACA, mat_corr,eval_kernel_ModeII_N)  
               
       !   REORDER IN NORMAL FORMAT      
       call backward_reorder(nb_element,tangent,mat_corr)
       call backward_reorder(nb_element,normal,mat_corr)
       call backward_reorder(nb_element,element,mat_corr)
       call backward_reorder(nb_element,node_left,mat_corr)
       call backward_reorder(nb_element,node_right,mat_corr)
       call backward_reorder(nb_element,normal_left,mat_corr)
       call backward_reorder(nb_element,normal_right,mat_corr)
      
       deallocate(mat_corr,ind_clusters,blcluster_tree)
    
       ! Create H matrix for tangential traction on the fault
       call initiate_hmat(nb_element,element,nbclus,ind_clusters,         &
                   nbblclus, blcluster_tree,info_blclustersT,             &
                   dimA_T,dimB_T,dimF_T,                                  &
                   A_matrices_T,B_matrices_T,F_matrices_T,                &
                   eta, nleaf, eps_ACA, mat_corr,eval_kernel_ModeII_T) 
  
   end if
end if



!
! Initialise permeability_x (for fluid diffusion)
permeability_x = 0._dp

! Initialise permeability based on pore pressure, sigmaN, Sk and kmin
if (permeability_coupling ==1)  then
    permeability_star = kmin + (permeability-kmin) * exp(abs(sigmaN+P)/abs(Snk))
else 
    permeability_star = 0._dp
endif


! Make the harmonic average
permeability_x(2:nb_element) = 2._dp*permeability(1:nb_element-1)*permeability(2:nb_element)/ &
                            (permeability(1:nb_element-1)+permeability(2:nb_element))





end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
subroutine find_fault_crossing()
real(kind=dp)                                      :: dist
integer                  :: fault_id1,fault_id2,node_id1,node_id2,crossing_id
!
! Find total number of crossing
print*,'Find fault '
!print*,'min(minval(faults(fault_id1)%ds),minval(faults(fault_id2)%ds))',min(minval(faults(1)%ds),minval(faults(2)%ds))

nb_crossing =0
do fault_id1 =1,nb_fault
    do fault_id2=fault_id1+1,nb_fault
        !
        ! for each pair of node
        do node_id1=1,faults(fault_id1)%nb_element+1
            do node_id2=1,faults(fault_id2)%nb_element+1
                ! Calculate the distance between the two nodes
                dist = sqrt((faults(fault_id1)%node(1,node_id1)-faults(fault_id2)%node(1,node_id2))**2 & 
                           +(faults(fault_id1)%node(2,node_id1)-faults(fault_id2)%node(2,node_id2))**2)
            
                ! If distance < min(ds)*01
                if (dist<0.01*min(minval(faults(fault_id1)%ds),minval(faults(fault_id2)%ds))) then
                    print*,'dist',dist
                    !
                    ! Update value of nb_crossing
                    nb_crossing = nb_crossing+1
                    !
                endif
            
            end do
        end do
    end do
end do
!
! Allocate info about crossing
allocate(crossing(nb_crossing,6))

crossing_id =0
do fault_id1 =1,nb_fault
    do fault_id2=fault_id1+1,nb_fault
        !
        ! for each pair of node
        do node_id1=1,faults(fault_id1)%nb_element+1
            do node_id2=1,faults(fault_id2)%nb_element+1
                ! Calculate the distance between the two nodes
                dist = sqrt((faults(fault_id1)%node(1,node_id1)-faults(fault_id2)%node(1,node_id2))**2 & 
                           +(faults(fault_id1)%node(2,node_id1)-faults(fault_id2)%node(2,node_id2))**2)
            
                ! If distance < min(ds)*01
                if (dist<0.01*min(minval(faults(fault_id1)%ds),minval(faults(fault_id2)%ds))) then
                       crossing_id = crossing_id + 1
                       crossing(crossing_id,1) = fault_id1
                       crossing(crossing_id,2) = fault_id2
                       crossing(crossing_id,3) = node_id1 + mapping_fault(fault_id1,1)-1 
                       crossing(crossing_id,4) = node_id2 + mapping_fault(fault_id2,1)-1 
                       crossing(crossing_id,5) = node_id1
                       crossing(crossing_id,6) = node_id2
                endif
            
            end do
        end do
    end do
end do

print*,'nb_crossing',nb_crossing

end subroutine

end module 
