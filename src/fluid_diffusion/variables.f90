module variables
implicit none
integer, parameter                                 :: dp=kind(1.d0)
!
!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Problem variable
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Meshing
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integer           	                               :: nbx
integer         								   :: len_P
real(kind=dp)                                      :: dx                       
                        

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Wells
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integer                                           :: nb_source
integer,dimension(:), allocatable             :: index_position_x,index_position_y
real(kind=dp), dimension(:), allocatable      :: t_beg 
real(kind=dp), dimension(:), allocatable      :: t_end
real(kind=dp), dimension(:), allocatable      :: injection_rate
real(kind=dp), dimension(:), allocatable      :: injection_pressure
integer(kind=dp), dimension(:), allocatable   :: constant_pressure
! type well
!     integer                                       :: nb_injection_phase=1
!     integer                                  	  :: index_position_x,index_position_y,index_position_z
! 	real(kind=dp), dimension(:), allocatable      :: t_beg 
! 	real(kind=dp), dimension(:), allocatable      :: t_end
! 	real(kind=dp), dimension(:), allocatable      :: injection_rate
! 	real(kind=dp), dimension(:), allocatable      :: injection_pressure
! 	integer(kind=dp), dimension(:), allocatable   :: constant_pressure
! end type
! 
! type(well), dimension(:), allocatable             ::    wells           ! Allocated in input


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Fluids
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real(kind=dp)                                      :: dyn_viscosity
real(kind=dp)                                      :: density              
real(kind=dp), dimension(:), pointer               :: P => null(), dPdt => null()
real(kind=dp), dimension(:), allocatable           :: P0      
real(kind=dp)									   :: fluid_compressibility

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Rock
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real(kind=dp)                                      :: porosity                
real(kind=dp)									   :: rock_compressibility
real(kind=dp), dimension(:,:), allocatable       :: permeability               
real(kind=dp), dimension(:,:), allocatable       :: permeability_x 



end module