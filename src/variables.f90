module variables
implicit none
!
!  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  INTERFACE FOR POINTERS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abstract interface
    subroutine friction_law(x,nv, Y, Y_dot)
    real(kind=8), INTENT(IN)                      :: x
    real(kind=8), dimension(:), intent(in)        :: Y
    real(kind=8), dimension(:), intent(out)       :: Y_dot
    integer :: nv
    end subroutine
    
    
    subroutine kernel(nelement,node,slip,iprec,stressT,stressN,ier,mu)
    integer                                        :: nelement
    real(kind=8), dimension(2,nelement)            :: node
    real(kind=8), dimension(nelement)              :: slip
    integer                                        :: iprec
    real(kind=8), dimension(nelement)              :: stressT
    real(kind=8), dimension(nelement)              :: stressN
    integer                                        :: ier
    real(kind=8)                                   :: mu
    end subroutine
    
    
    function eval_kernel_hmat(ind_row,ind_col,Ynode,cnod_row,nbnod_row,id_row,id_col)
        integer                                       :: ind_row, ind_col
        integer                                       :: nbnod_row
        real(kind=8),dimension(2)                     :: ynode
        real(kind=8),dimension(2,nbnod_row)           :: cnod_row
        integer                                       :: id_row, id_col
        complex(kind=8):: eval_kernel_hmat
    end function
end interface
!
!   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PARAMETERS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integer, parameter                                 :: dp=kind(1.d0)
real(kind=dp), parameter                           :: pi=3.141592653589793238_dp
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! HYPERPARAMETERS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
character(len=200)                                 :: simulation_name              ! Read in config.in
character(len=200)                                 :: config_file = ""             !"config.in"
character(len=200)                                 :: path_problem = ""             !"config.in"
character(len=200)                                 :: fracture_mode                ! Read in config.in
character(len=200)                                 :: fric_law          ! Read in config.in
character(len=200)                                 :: static_kernel                ! Read in config.in
integer                                            :: max_it
integer                                            :: stride_time                  ! Save time each stride_time
integer                                            :: freq_writing_file            ! Read in config.in
integer                                            :: permeability_coupling        



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! MESHING VARIABLES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integer                                            :: nb_element                   ! Calculated in interpolater
integer                                            :: nb_node                      ! Calculated in interpolater
integer                                            :: nb_crossing                      ! Calculated in interpolater
integer(kind=dp), dimension(:,:), allocatable      :: crossing
integer(kind=dp), dimension(:,:), allocatable      :: mapping_fault                 ! mapping(fault_id,1) = first_index; mapping(fault_id,2) = last_index



real(kind=dp), dimension(:,:), allocatable, target :: element                      ! Allocated and calculated in initialisation
real(kind=dp), dimension(:), allocatable, target   :: ds                   ! Allocated and calculated in initialisation
real(kind=dp), dimension(:,:), allocatable, target :: node                         ! Allocated in initialisation
real(kind=dp), dimension(:,:), allocatable, target :: tangent                      ! Allocated in initialisation
real(kind=dp), dimension(:,:), allocatable, target :: normal                       ! Allocated in initialisation
!
!
!   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROBLEM VARIABLES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integer                                            :: nb_fault                     ! Read in data.in
real(kind=dp)                                      :: time
integer                                            :: length_Y
real(kind=dp), dimension(:), allocatable, target   :: Y            ! Big array composed of values of V and theta (Mode III) or V, theta and sigmaN (Mode II), allocated in initialisation
real(kind=dp), dimension(:), allocatable, target   :: Y_dot        ! Derivative of Y
real(kind=dp), dimension(:), pointer               :: V
real(kind=dp), dimension(:), pointer               :: theta
real(kind=dp)                                      :: mu                   ! Shear modulus
real(kind=dp), dimension(:), allocatable, target   :: slip                 ! Slip on the fault
real(kind=dp)                                      :: Sigma32_dot     ! Loading stress on the faults
real(kind=dp)                                      :: Sigma31_dot      ! Loading stress on the faults
real(kind=dp)                                      :: Sigma12_dot      ! Loading stress on the faults
real(kind=dp)                                      :: Sigma11_dot      ! Loading stress on the faults
real(kind=dp)                                      :: Sigma22_dot     ! Loading stress on the faults
real(kind=dp)                                      :: Sigma32     ! Loading stress on the faults
real(kind=dp)                                      :: Sigma31      ! Loading stress on the faults
real(kind=dp)                                      :: Sigma12     ! Loading stress on the faults
real(kind=dp)                                      :: Sigma11     ! Loading stress on the faults
real(kind=dp)                                      :: Sigma22    ! Loading stress on the faults

real(kind=dp), dimension(:), allocatable, target   :: normal_loading_dot! Background Loading (calculated in initialisation)
real(kind=dp), dimension(:), allocatable, target   :: shear_loading_dot! Background Loading (calculated in initialisation)
real(kind=dp), dimension(:), pointer               :: SigmaN               ! Total Normal traction on the faults   
real(kind=dp), dimension(:), allocatable,target    :: SigmaN_dot           ! Derivative of the normal traction on the faults
!real(kind=dp), dimension(:), allocatable, target   :: SigmaT               ! Tangential traction on the faults
real(kind=dp), dimension(:), allocatable, target   :: SigmaT_dot           ! Derivative of the tangential traction on the faults
! real(kind=dp)                                      :: SigmaN_yield         ! Maximum normal stress on the fault
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FAST-MULTIPOLE-METHOD VARIABLES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integer                                            :: ier
integer                                            :: iprec
!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! H-MATRICES VARIABLES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integer(kind=4)                                    ::  nleaf=100
real(kind=8)                                       ::  eta=3.d0
REAL(kind=8)                                       ::  eps_ACA=1.000e-8 
INTEGER(kind=4), dimension(:)  , allocatable       ::  mat_corr
REAL(kind=dp)   , dimension(:,:), allocatable      ::  node_left, node_right
REAL(kind=dp)   , dimension(:,:), allocatable      ::  normal_right,normal_left
INTEGER(kind=4)                                    ::  nbclus
INTEGER(kind=4),DIMENSION(:),ALLOCATABLE           ::  ind_clusters
INTEGER(kind=4)                                    ::  nbblclus
INTEGER(kind=4),DIMENSION(:,:),ALLOCATABLE         ::  blcluster_tree
INTEGER(kind=8),ALLOCATABLE,DIMENSION(:,:)         ::  info_blclustersT
INTEGER(kind=8),ALLOCATABLE,DIMENSION(:,:)         ::  info_blclustersN
INTEGER(kind=8)                                    ::  dimA_T, dimB_T, dimF_T
INTEGER(kind=8)                                    ::  dimA_N, dimB_N, dimF_N
COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE         ::  A_matrices_T
COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE         ::  B_matrices_T
COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE         ::  F_matrices_T
COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE         ::  A_matrices_N
COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE         ::  B_matrices_N
COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE         ::  F_matrices_N
!



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fluids
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real(kind=dp)                                      :: dyn_viscosity
real(kind=dp)                                      :: fluid_density              
real(kind=dp), dimension(:), allocatable, target                        :: P  
real(kind=dp), dimension(:), allocatable, target                        :: P_dot  


real(kind=dp), dimension(:), pointer               :: E, E_dot ! Hydraulic opening
real(kind=dp)                                      :: fluid_comp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Rock
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real(kind=dp)                                      :: porosity                
real(kind=dp)                                      :: rock_comp
real(kind=dp), dimension(:), allocatable, target   :: permeability 
 
real(kind=dp), dimension(:), pointer   :: permeability_star 
real(kind=dp), dimension(:), allocatable, target   :: permeability_x 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! All the parameters for variable permeability
real(kind=dp), dimension(:), allocatable, target   :: Lk, Tk, kmin, kmax, Snk ! (From Zhu et al., 2000)




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SOLVER VARIABLES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
real(kind=dp)                                      :: tol_solver
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! RATE AND STATE VARIABLES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
real(kind=dp), dimension(:), allocatable, target   :: V0, f0
real(kind=dp), dimension(:), allocatable, target   :: a, b
real(kind=dp), dimension(:), allocatable, target   :: Dc
real(kind=dp)                                      :: cs
real(kind=dp)                                      :: cp
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! !   INTERNAL VARIABLES: Pointers
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
procedure(kernel), pointer                         :: kernel_ptr   => null()   ! Allocated in initialisation
procedure(friction_law), pointer                   :: fric_law_ptr => null()   ! Allocated in initialisation
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! INTERNAL VARIABLES: Type
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
type fault
    integer                                         :: nb_element           
    real(kind=dp), dimension(:,:), pointer          :: node             => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:,:), pointer          :: tangent          => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:,:), pointer          :: normal           => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: ds               => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: slip             => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: V                => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: theta            => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: a                => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: b                => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: Dc               => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: V0               => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: f0               => null()    ! TODO NOT USED YET: Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: SigmaN           => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: SigmaN_dot       => null()    ! TODO NOT USED YET: Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: SigmaT           => null()    ! TODO NOT USED YET: Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: SigmaT_dot       => null()    ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: P                => null()   ! Allocated in initialisation
    real(kind=dp), dimension(:), pointer            :: P_dot                => null()   ! Allocated in initialisation

    real(kind=dp), dimension(:), pointer            :: normal_loading_dot => null()
    real(kind=dp), dimension(:), pointer            :: shear_loading_dot  => null()
    real(kind=dp), dimension(:), pointer            :: permeability     => null()
    real(kind=dp), dimension(:), pointer            :: permeability_x   => null()
    real(kind=dp), dimension(:), pointer            :: permeability_star => null()

    ! Variable permeability
    real(kind=dp), dimension(:), pointer            :: Lk => null()
    real(kind=dp), dimension(:), pointer            :: Tk => null()
    real(kind=dp), dimension(:), pointer            :: Snk => null()
    real(kind=dp), dimension(:), pointer            :: kmin => null()
    real(kind=dp), dimension(:), pointer            :: kmax => null()


    
    ! Source of fluid
    integer                                         :: nb_source                                                          
    real(kind=dp), dimension(:), allocatable        :: Q  
    real(kind=dp), dimension(:), allocatable        :: t_injection_beg
    real(kind=dp), dimension(:), allocatable        :: t_injection_end
    integer(kind=dp), dimension(:), allocatable     :: index_injection
    
end type
!
!   INTERNAL VARIABLES: type(fault) faults
type(fault), dimension(:), allocatable              :: faults                   ! Allocated in input
!
end module variables