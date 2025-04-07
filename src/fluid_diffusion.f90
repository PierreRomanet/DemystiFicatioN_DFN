module fluid_diffusion
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FINITE VOLUME METHOD 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   
!   Compute Pressure P for new time using implicit scheme
!   TODO: change ds that is wrong for the FVM (it is the size of element, not the distance between the two points of FVM)
!   global input: permeablity, P
use variables
private
public :: compute_P 





contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   
!   Compute Pressure P for new time using implicit scheme
!   It is using conjugate gradient method with a sparse matrix representation
!
subroutine compute_P(dt,P_temp,P_dot_temp,ier)
real(kind=dp)      :: dt
integer(kind=dp)   :: fault_id, source_id,iteration_count
integer            :: ier
real(kind=dp), dimension(nb_element) :: pk,rk,rk1,xk,P_temp, P_dot_temp       ! Tmp array for implicit solver
integer                                     :: k                              ! Integer for for-loop  

! Variable related to fluid diffusion
real(kind=dp)                               :: alpha_k, beta_k,max_imp, errmax      ! Tmp array for implicit solver
real(kind=dp), dimension(nb_element)  :: b1,diag, Pscal                          ! Tmp array for implicit solver



! Initialisation
ier = 0 ! At the beginning there is no problem


! Preconditioner
call compute_diag(dt,diag)

!######### Compute the vector b
do k=1,nb_element
! Compute vector b1
b1(k) = P(k)/diag(k)
enddo


! Add injection point        
do fault_id=1,nb_fault
    ! For each source 
    do source_id=1,faults(fault_id)%nb_source
        ! If time of injection
        if ((time+dt>=faults(fault_id)%t_injection_beg(source_id) ).and.(time+dt< faults(fault_id)%t_injection_end(source_id))) then
            b1(mapping_fault(fault_id,1)-1 + faults(fault_id)%index_injection) = & 
                                   b1(mapping_fault(fault_id,1)-1 + faults(fault_id)%index_injection) &
                                 + (dt*faults(fault_id)%Q/(porosity*(fluid_comp+rock_comp)) &
                                 /ds(mapping_fault(fault_id,1)-1 + faults(fault_id)%index_injection))/&
                                 diag(mapping_fault(fault_id,1)-1 + faults(fault_id)%index_injection)
        end if
    enddo
enddo


!######### Compute the permeability based on permeability_star
if (permeability_coupling == 1) then
permeability = (permeability_star - kmin)*exp(-abs(sigmaN+P)/abs(Snk)) + kmin
end if

!######### Compute the new permeability_x
! Make the harmonic average
permeability_x(2:nb_element) = 2._dp*permeability(1:nb_element-1)*permeability(2:nb_element)/ &
                            (permeability(1:nb_element-1)+permeability(2:nb_element))


! 
!######### Conjugate gradient
! Initialisation
xk = P
rk = b1-compute_MatrixVector(xk,dt,diag)
!print*,'rk ini', dsqrt(sum(rk**2))/nb_element
pk = rk
iteration_count = 0


! Tentative of scaling for the error
Pscal = abs(P) +dt*abs(Pdot)+ tol_solver
Pscal = 1._dp
errmax=maxval(abs(rk(:)/Pscal(:)))/tol_solver !Evaluate accuracy.
! print*,'errmax',errmax
! print*,'Pscal',Pscal

! print*,'dt',dt
! do while ((dsqrt(sum(rk**2))/nb_element>=tol_solver))
do while (errmax >= 1.0)

! print*,'rk',dsqrt(sum(rk**2))/nb_element

! Count the number of iteration
iteration_count = iteration_count + 1

alpha_k = sum(rk**2)/sum(pk*compute_MatrixVector(pk,dt,diag))

xk = xk + alpha_k*pk
rk1 = rk - alpha_k*compute_MatrixVector(pk,dt,diag)

! Update pk
beta_k = sum(rk1**2)/sum(rk**2)
pk = rk1 + beta_k*pk


rk = rk1

! if (maxval(abs(alpha_k*pl)-tol_solver*abs(xk))<=0._dp) then
!     return
!     endif

if(iteration_count>200)then
  ier = 1
    print*,' did no converge (compute_P)','dt = ',dt
  return
!print*,'dsqrt(sum(rk**2))/nb_element',dsqrt(sum(rk**2))/nb_element
end if


errmax=maxval(abs(rk(:)/Pscal(:)))/tol_solver !Evaluate accuracy.
!print*,'maxval(P)',minval(xk),maxval(xk)

end do



!######### Calculate new P_dot at new time
if (dt<1e-9) then
P_dot = 0._dp
else
P_dot = (xk-P)/dt
endif




! Update P
P_temp = xk

if (maxval(abs(P_temp))>huge(1._dp)) then
 ier = 1
 end if

end subroutine





!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   
!   Compute Pressure P for new time using implicit scheme
!   Function that return the result of the matrix vector multiplication with A 
!
function compute_MatrixVector(vector,dt,diag)
real(kind=dp), dimension(nb_element)    ::  compute_MatrixVector,vector
real(kind=dp)                           :: dt,A1,A2,A3
integer(kind=dp)      :: crossing_id, fault_id1, fault_id2,node_id1,node_id2
real(kind=dp)                                :: permeability_x_crossing_left, permeability_x_crossing_right
real(kind=dp), dimension(nb_element)  :: diag   ! Tmp array for implicit solver
integer                               :: k

!-----------------------------------------------------------------------------------------
! Without coupling between faults
!-----------------------------------------------------------------------------------------
! Fist and last componants
A2 = 1._dp+(permeability_x(2))*dt &
                /(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(1)**2)
A3 = -(permeability_x(2))*dt   &
           /(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(1)**2)
             
   
compute_MatrixVector(1) = (A2*vector(1)+A3*vector(2))/diag(1)
           
A1 = -(permeability_x(nb_element))*dt   &
                               /(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(nb_element)**2)
A2 = 1._dp+(permeability_x(nb_element))*dt &
                                  /(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(nb_element)**2)

compute_MatrixVector(nb_element) = (A1*vector(nb_element-1)+A2*vector(nb_element))/diag(nb_element)

! For each element on the diagonal
do k=2,nb_element-1
    A1 = -(permeability_x(k))*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)
    
    
    A2 = 1._dp+(permeability_x(k)+permeability_x(k+1))   &
                    *dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k)**2)


    A3 = -(permeability_x(k+1))*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k+1)**2)
    
    
    compute_MatrixVector(k) = (A1*vector(k-1)+A2*vector(k)+A3*vector(k+1))/diag(k)
    
end do 



!-----------------------------------------------------------------------------------------
! Add the coupling between the fault network
!-----------------------------------------------------------------------------------------
!######### Add Coupling value when faults are crossing
! For each time faults are crossing each other
do crossing_id=1,nb_crossing
    ! Only for readability
    fault_id1 = crossing(crossing_id,1)
    fault_id2 = crossing(crossing_id,2)
    node_id1 = crossing(crossing_id,3)
    node_id2 = crossing(crossing_id,4)
    
    ! Update pressure gradient for the 4 points at the triple junction
    !#!# element1(node_id1-1) 
    ! Calculate permeability at the interface
    permeability_x_crossing_left = 2._dp*permeability(node_id1-1) * permeability(node_id2-1)/&
                                        (permeability(node_id1-1) + permeability(node_id2-1))
    permeability_x_crossing_right = 2._dp*permeability(node_id1-1) * permeability(node_id2)/&
                                        (permeability(node_id1-1) + permeability(node_id2))
                                        
    A1 = -permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)
    A3 =   -permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)   
    A2 = permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2) +  &
         permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)
    
    
    
    
    compute_MatrixVector(node_id1-1) = compute_MatrixVector(node_id1-1) + &
                 (A1*vector(node_id2-1) +A2*vector(node_id1-1)+A3*vector(node_id2))/diag(node_id1-1) 
    
    !#!# element1(node_id1)
    ! Calculate permeability at the interface
    permeability_x_crossing_left = 2._dp*permeability(node_id1) * permeability(node_id2-1)/&
                                        (permeability(node_id1) + permeability(node_id2-1))
    permeability_x_crossing_right = 2._dp*permeability(node_id1) * permeability(node_id2)/&
                                        (permeability(node_id1) + permeability(node_id2))
    !       
    A1 = -permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)
    A3 =   -permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)   
    A2 = permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2) +  &
               permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)

     compute_MatrixVector(node_id1) = compute_MatrixVector(node_id1) + &
                 (A1*vector(node_id2-1) +A2*vector(node_id1)+A3*vector(node_id2))/(node_id1) 




    !#!# element2(node_id2-1)
    ! Calculate permeability at the interface
    permeability_x_crossing_left = 2._dp*permeability(node_id1-1) * permeability(node_id2-1)/&
                                        (permeability(node_id1-1) + permeability(node_id2-1))
    permeability_x_crossing_right = 2._dp*permeability(node_id1) * permeability(node_id2-1)/&
                                        (permeability(node_id1) + permeability(node_id2-1))
                                         
                                     
    
    !            
    A1 = -permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)
    A3 =   -permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)   
    A2 = permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2) +  &
               permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)
           
     compute_MatrixVector(node_id2-1) = compute_MatrixVector(node_id2-1) + &
                 (A1*vector(node_id1-1) +A2*vector(node_id2-1)+A3*vector(node_id1))/diag(node_id2-1) 




    !#!# element2(node_id2)
    ! Calculate permeability at the interface
    permeability_x_crossing_left = 2._dp*permeability(node_id1-1) * permeability(node_id2)/&
                                        (permeability(node_id1-1) + permeability(node_id2))
    permeability_x_crossing_right = 2._dp*permeability(node_id1) * permeability(node_id2)/&
                                        (permeability(node_id1) + permeability(node_id2))
    !        
    A1 = -permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)
    A3 =   -permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)   
    A2 = permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2) +  &
        permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k-1)**2)
    
     compute_MatrixVector(node_id2) = compute_MatrixVector(node_id2) + &
                 (A1*vector(node_id1-1) +A2*vector(node_id2)+A3*vector(node_id1))/diag(node_id2)
    

    
end do


end function




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   
!   Compute Pressure P for new time using implicit scheme
!   Function that return the result of the matrix vector multiplication with A 
!
subroutine compute_diag(dt,diag)
real(kind=dp)                           :: dt
integer(kind=dp)      :: crossing_id, fault_id1, fault_id2,node_id1,node_id2
real(kind=dp)                                :: permeability_x_crossing_left, permeability_x_crossing_right
real(kind=dp), dimension(nb_element)  :: diag   ! Tmp array for implicit solver
integer                               :: k

!


! compute_MatrixVector(:) = 0._dp
diag = 0._dp

! Fist and last componants
diag(1) = 1._dp+(permeability_x(2))*dt &
                /(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(1)**2)
           

diag(nb_element) = 1._dp+(permeability_x(nb_element))*dt &
                                  /(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(nb_element)**2)


! For each element on the diagonal
do k=2,nb_element-1  
    diag(k) = 1._dp+(permeability_x(k)+permeability_x(k+1))   &
                    *dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(k)**2)   
end do 




!######### Add Coupling value when faults are crossing
! For each time faults are crossing each other
do crossing_id=1,nb_crossing
    ! Only for readability
    fault_id1 = crossing(crossing_id,1)
    fault_id2 = crossing(crossing_id,2)
    node_id1 = crossing(crossing_id,3)
    node_id2 = crossing(crossing_id,4)
    
    ! Update pressure gradient for the 4 points at the triple junction
    !#!# element1(node_id1-1) 
    ! Calculate permeability at the interface
    permeability_x_crossing_left = 2._dp*permeability(node_id1-1) * permeability(node_id2-1) &
                                   /(permeability(node_id1-1)+permeability(node_id2-1))
    permeability_x_crossing_right = 2._dp*permeability(node_id1-1) * permeability(node_id2) &
                                  /(permeability(node_id1-1) + permeability(node_id2))
               
   
    diag(node_id1-1) =  diag(node_id1-1) & 
       + permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(1)**2) +  &
         permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(1)**2)
    
    
    !#!# element1(node_id1)
    ! Calculate permeability at the interface
    permeability_x_crossing_left = 2._dp*permeability(node_id1) * permeability(node_id2-1) &
                                         /(permeability(node_id1) + permeability(node_id2-1))
    permeability_x_crossing_right = 2._dp*permeability(node_id1) * permeability(node_id2)  &
                                          /(permeability(node_id1) + permeability(node_id2))
    !       
   
    diag(node_id1) = diag(node_id1) & 
             + permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(1)**2) +  &
               permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(1)**2)



    !#!# element2(node_id2-1)
    ! Calculate permeability at the interface
    permeability_x_crossing_left = 2._dp*permeability(node_id1-1) * permeability(node_id2-1)&
                                        /(permeability(node_id1-1) + permeability(node_id2-1))
    permeability_x_crossing_right = 2._dp*permeability(node_id1) *  permeability(node_id2-1)&
                                       / (permeability(node_id1) +  permeability(node_id2-1))
                                         
                                     
    
    !            
    diag(node_id2-1) = diag(node_id2-1) & 
              +permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(1)**2) +  &
               permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(1)**2)
           
    

    !#!# element2(node_id2)
    ! Calculate permeability at the interface
    permeability_x_crossing_left = 2._dp*permeability(node_id1-1) * permeability(node_id2)&
                                        /(permeability(node_id1-1) + permeability(node_id2))
    permeability_x_crossing_right = 2._dp*permeability(node_id1) * permeability(node_id2)&
                                        /(permeability(node_id1) + permeability(node_id2))
    !        
    diag(node_id2) = diag(node_id2) & 
       + permeability_x_crossing_left*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(1)**2) +  &
        permeability_x_crossing_right*dt/(porosity*dyn_viscosity*(fluid_comp+rock_comp)*ds(1)**2)
    

end do
end subroutine

end module