module ACA_algorithm 
use get_info_matrices
implicit none
private
public :: ACA,compute_full_ACA


contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE ACA(cnod_row,nbnod_row,cnod_col,nbnod_col,&
     A_matrices,dimA,B_matrices,dimB,eps_ACA,&
     low_rank,max_low_rank,flag_ACA,ind_rec1,id_row,id_col,eval_kernel)
                

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: ACA                          
!                                                                    
! INPUTS: nbnod_row  : number of dof in cluster row
!         nbnod_col  : number of dof in cluster column
!         dimA       : dimension of array A_matrices           
!         dimB       : dimension of array B_matrices           
!         eps_ACA    : ACA threshold                          
!         cnod_row   : coordinates of points corresponding to dofs
!                      in cluster row          
!         cnod_col   : coordinates of points corresponding to dofs
!                      in cluster column          
!
! OUTPUTS: A_matrices : array where are stored the matrices A
!          B_matrices : array where are stored the matrices B
!          low_rank   : low-rank of the approximation  
!                                                                    
! INPUTS/OUTPUTS: flag_ACA : flag whose value is TRUE if the memory
!                            is enough to the approximation, FALSE
!                            otherwise.  
!                                                                    
! GOAL: Compute the low-rank approximation of an admissible block
!       using ACA (Adaptive Cross Approximation)  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

abstract interface
    function eval_kernel(ind_row,ind_col,Ynode,cnod_row,nbnod_row,id_row,id_col)
        ! Input
        REAL(kind=8),DIMENSION(2)            :: Ynode
        INTEGER                              :: ind_row,ind_col,nbnod_row, id_row,id_col
        REAL(kind=8),DIMENSION(2,nbnod_row)  :: cnod_row
        !Output
        COMPLEX(kind=8):: eval_kernel
    end function
end interface        
 

  !Input
  INTEGER(kind=4),INTENT(IN):: nbnod_row, nbnod_col 
  INTEGER(kind=4),INTENT(IN):: max_low_rank,ind_rec1
  INTEGER(kind=4),INTENT(IN):: id_row, id_col

  
  INTEGER(kind=8),INTENT(IN):: dimA, dimB
 
  REAL(kind=8),INTENT(IN):: eps_ACA
 
  REAL(kind=8),DIMENSION(2,nbnod_row),INTENT(IN):: cnod_row
  REAL(kind=8),DIMENSION(2,nbnod_col),INTENT(IN):: cnod_col
  
  !Input/Output
  LOGICAL,INTENT(INOUT):: flag_ACA
  !Output
  INTEGER,INTENT(OUT):: low_rank
  COMPLEX(kind=8),DIMENSION(dimA),INTENT(OUT):: A_matrices
  COMPLEX(kind=8),DIMENSION(dimB),INTENT(OUT):: B_matrices

  !Local Variables
  INTEGER(kind=4):: ind_rec 
  
  INTEGER(kind=4),DIMENSION(:),ALLOCATABLE:: ind_row_picked, ind_row
  INTEGER(kind=4),DIMENSION(:),ALLOCATABLE:: ind_col_picked, ind_col

  REAL(kind=8):: norm_app
 
  !Function Variables
  
  REAL(kind=8),EXTERNAL:: DZNRM2, ZLANGE

  
  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!
 
 
  ALLOCATE(ind_row_picked(nbnod_row+1),ind_row(nbnod_row+1))
  ALLOCATE(ind_col_picked(nbnod_col+1),ind_col(nbnod_col+1))
   
  !Indices of the columns picked up from the block
  ind_col_picked=0
  !Indices of the rows picked up from the block
  ind_row_picked=0

  !Row indices to search for maximum in error matrix
  ind_row=(/(ind_rec, ind_rec=1,nbnod_row)/)
  !Column indices to search for maximum in error matrix
  ind_col=(/(ind_rec, ind_rec=1,nbnod_col)/)

!!$  print*,'size mat',nbnod_row,nbnod_col
!!$  ind_file=100+int(rand(0)*10)
!!$  OPEN(ind_file,FORM='FORMATTED',ACTION='WRITE')
!!$  ALLOCATE(mat_ACA(nbnod_row,nbnod_col))
!!$  ind=0
!!$  DO ind1=1,nbnod_row
!!$     DO ind2=1,nbnod_col
!!$        ind=ind+1
!!$
!!$        
!!$        mat_ACA(ind1,ind2)=1/sqrt((cnod_row(1,ind1)-cnod_col(1,ind2))**2+ &
!!$             (cnod_row(2,ind1)-cnod_col(2,ind2))**2)
!!$       
!!$        WRITE(ind_file,*)  REAL(mat_ACA(ind1,ind2))
!!$
!!$         IF(sqrt((cnod_row(1,ind1)-cnod_col(1,ind2))**2+ &
!!$             (cnod_row(2,ind1)-cnod_col(2,ind2))**2)<0.00000000001)THEN
!!$           print*,'zero row',ind1,cnod_row(:,ind1)
!!$           print*,'zero col',ind2,cnod_col(:,ind2)
!!$        END IF
!!$     END DO
!!$  END DO
!!$  print*,'numero fichier',ind_file, ind
!!$ DEALLOCATE(mat_ACA)
!!$ CLOSE(ind_file)
  
  !Initialize the norm of approximate of the matrix to test error
  norm_app=0.d0

  !Initialize the 1st row index
  ind_row_picked(1)=1
 
  !Nullify the first 1st row index to search for maximum in error
  !matrix
  ind_row(1)=0
 
  
  DO ind_rec=1,MIN(nbnod_row,nbnod_col)
  
     !Check the status of memory
     IF(ind_rec.GT.max_low_rank) THEN
        print*, ind_rec, max_low_rank
        STOP 'Memory is not enough, see subroutine ACA   ****'
     END IF
     
     !Compute the ind_rec-th row of the approximate error matrix 
     CALL compute_Rik(cnod_row,nbnod_row,cnod_col,nbnod_col,&
          ind_row_picked(ind_rec),A_matrices,dimA,B_matrices,&
          dimB,ind_rec,id_row,id_col,eval_kernel)
                


!!$     IF(ind_rec1==15)then        
!!$        ! print*,'B',B_matrices(1+(ind_rec-1)*nbnod_col:ind_rec*nbnod_col)
!!$        print*,'A',A_matrices(1+(ind_rec-1)*nbnod_row:ind_rec*nbnod_row)
!!$        print*,'row',cnod_row(1,1)
!!$        print*,'col',cnod_col(1,ind_col)
!!$        print*,'kernel',eval_kernel(1, cnod_col(:,ind_col))
!!$        pause
!!$     end IF
     !Find the ind_rec-th column index
     
     ind_col_picked(ind_rec)=find_max_ACA(nbnod_col,ind_col,&
          B_matrices(1+(ind_rec-1)*nbnod_col:ind_rec*nbnod_col))
     
     !Nullify the ind_rec-th column index to search for maximum in
     !error matrix
     ind_col(ind_col_picked(ind_rec))=0

     !Terminate if the maximum of the row is 0
     IF(CDABS(B_matrices(ind_col_picked(ind_rec)+(ind_rec-1)*&
          nbnod_col)).LE.10.0**(-14)) THEN
        B_matrices(1+(ind_rec-1)*nbnod_col:ind_rec*nbnod_col)=&
             (0.d0,0.d0)
        low_rank=ind_rec-1
        flag_ACA=.TRUE.
       
        EXIT
     END IF
    
     !Normalize the ind_rec-th row
     CALL normalize_Rik(nbnod_col,ind_col_picked(ind_rec),&
          B_matrices(1+(ind_rec-1)*nbnod_col:ind_rec*nbnod_col))

!!$     IF(ind_rec1==15)then
!!$         print*,'B',B_matrices(1+(ind_rec-1)*nbnod_col:ind_rec*nbnod_col)
!!$        !print*,'A',A_matrices(1+(ind_rec-1)*nbnod_row:ind_rec*nbnod_row)
!!$     end IF
     !Compute the ind_rec-th column of the approximate error matrix
     CALL compute_Rjk(cnod_row,nbnod_row,cnod_col,nbnod_col,&
          ind_col_picked(ind_rec),A_matrices,&
          dimA,B_matrices,dimB,ind_rec,id_row,id_col,eval_kernel)
      
     !Compute the norm of approximate of the matrix
     CALL compute_norm_ACA(ind_rec,A_matrices,dimA,B_matrices,&
          dimB,nbnod_row,nbnod_col,norm_app)
 
     !Check the convergence
!!$     IF(ind_rec1==15)then
!!$   !     print*,'B',B_matrices(1+(ind_rec-1)*nbnod_col:ind_rec*nbnod_col)
!!$   !     print*,'A',A_matrices(1+(ind_rec-1)*nbnod_row:ind_rec*nbnod_row)
!!$        print*,DZNRM2(nbnod_row,A_matrices(1+(ind_rec-1)*nbnod_row:&
!!$             ind_rec*nbnod_row),1)*DZNRM2(nbnod_col,B_matrices(1+&
!!$             (ind_rec-1)*nbnod_col:ind_rec*nbnod_col),1)
!!$        print*,'eps',eps_ACA*DSQRT(norm_app)
!!$       ! pause
!!$     end IF
     IF(DZNRM2(nbnod_row,A_matrices(1+(ind_rec-1)*nbnod_row:&
          ind_rec*nbnod_row),1)*DZNRM2(nbnod_col,B_matrices(1+&
          (ind_rec-1)*nbnod_col:ind_rec*nbnod_col),1).LE.eps_ACA*&
          DSQRT(norm_app)) THEN
        low_rank=ind_rec
        flag_ACA=.TRUE.
      !  print*,'low 2',low_rank
        EXIT
     END IF

     !Check the status of memory
     IF(((ind_rec+1)*nbnod_row.GT.dimA).OR.&
          ((ind_rec+1)*nbnod_col.GT.dimB)) THEN
        low_rank=ind_rec
       ! print*,'low 3',low_rank
        flag_ACA=.FALSE.
        EXIT
     END IF

     !Find the next row index
     ind_row_picked(ind_rec+1)=find_max_ACA(nbnod_row,ind_row,&
          A_matrices(1+(ind_rec-1)*nbnod_row:ind_rec*nbnod_row))
     
     !Nullify the next row index to search for maximum in error
     !matrix
     ind_row(ind_row_picked(ind_rec+1))=0
     
  END DO

  DEALLOCATE(ind_row_picked,ind_row)
  DEALLOCATE(ind_col_picked,ind_col)

  IF(low_rank.EQ.0)THEN
     !print*,'Full rank'
     !TODO pause
     low_rank=ind_rec
  END IF

  print*,'flag_ACA',flag_ACA
  print*,'low_rank',low_rank
!!$  IF(ind_rec1==15)then
!!$     pause
!!$  end IF
END SUBROUTINE ACA
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE compute_Rik (cnod_row,nbnod_row,cnod_col,nbnod_col,&
     ind_row,A_matrices,dimA,B_matrices,dimB,ind_rec,id_row,id_col,eval_kernel)
                
       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: compute_Rik                             
!                                                                    
! INPUTS: nbnod_row  : number of dof in cluster row
!         nbnod_col  : number of dof in cluster column
!         dimA       : dimension of array A_matrices           
!         dimB       : dimension of array B_matrices           
!         ind_row    : index of the current row  
!         ind_rec    : index of recursion in ACA           
               
!         cnod_row   : coordinates of points corresponding to dofs
!                      in cluster row          
!         cnod_col   : coordinates of points corresponding to dofs
!                      in cluster column          
!         A_matrices : array where are stored the matrices A
!                                                                    
! INPUTS/OUTPUTS: B_matrices : array where are stored the matrices B
!                                                                    
! GOAL: Compute a row of the approximate error matrix
!                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

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


  !Input
  INTEGER(kind=4),INTENT(IN):: nbnod_row, nbnod_col
  INTEGER(kind=8),INTENT(IN):: dimA, dimB
  INTEGER(kind=4),INTENT(IN):: ind_row, ind_rec
  INTEGER(kind=4),INTENT(IN):: id_row, id_col

 ! REAL(kind=8):: wave_num, alpha

  REAL(kind=8),DIMENSION(2,nbnod_row),INTENT(IN):: cnod_row
  REAL(kind=8),DIMENSION(2,nbnod_col),INTENT(IN):: cnod_col

  COMPLEX(kind=8),DIMENSION(dimA),INTENT(IN):: A_matrices

  !Input/Output
  COMPLEX(kind=8),DIMENSION(dimB),INTENT(INOUT):: B_matrices

  !Local Variables
  INTEGER(kind=4):: ind_col
  
  COMPLEX(kind=8),DIMENSION(ind_rec-1):: a,  b

  
  COMPLEX(kind=8),EXTERNAL:: ZDOTU

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!
  
  !Compute the current row of matrix A
  a=A_matrices(ind_row:(ind_rec-1)*nbnod_row:nbnod_row)
   
  !Loop on the number of columns
  DO ind_col=1,nbnod_col
     
     !Evaluation of the Helmholtz kernel function
     B_matrices(ind_col+(ind_rec-1)*nbnod_col)=&
          eval_kernel(ind_row,ind_col,cnod_col(:,ind_col),cnod_row,nbnod_row,id_row, id_col)

     !Compute the current column of matrix B
     b=B_matrices(ind_col:(ind_rec-1)*nbnod_col:nbnod_col)

     !Update the current row of the approximate error matrix
     B_matrices(ind_col+(ind_rec-1)*nbnod_col)=&
          B_matrices(ind_col+(ind_rec-1)*nbnod_col)-&
          ZDOTU(ind_rec-1,a,1,b,1)

  END DO
  
  RETURN

END SUBROUTINE compute_Rik
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE normalize_Rik(dim_Rik,ind_max,Rik)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: normalize_Rik                             
!                                                                    
! INPUTS: dim_Rik  : dimension of the vector Rik
!         ind_max  : index of the maximum element of Rik
!                                                                    
! INPUTS/OUTPUTS: Rik : current row of the approximate error matrix
!                                                                    
! GOAL: Normalize a row of the approximate error matrix
!
! BIBLIOGRAPHY: K.Zhao, "The Adaptive Cross Approximation Algorithm
!               for Accelerated Method of Moments Computations of   
!               EMC Problems" 
! 
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/23/2014                                                  
!                                                                    
! MODIFIED: VIII/29/2015 (Luca Desiderio)                             
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: dim_Rik, ind_max

  !Input/Output
  COMPLEX(kind=8),DIMENSION(dim_Rik),INTENT(INOUT):: Rik

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  CALL ZSCAL(dim_Rik,CMPLX(1.d0,0.d0,8)/Rik(ind_max),Rik,1)

  RETURN

END SUBROUTINE normalize_Rik
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE compute_Rjk (cnod_row,nbnod_row,cnod_col,nbnod_col,&
     ind_col,A_matrices,dimA,B_matrices,dimB,ind_rec,id_row,id_col,eval_kernel)
                
     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: compute_Rjk                             
!                                                                    
! INPUTS: nbnod_row  : number of dof in cluster row
!         nbnod_col  : number of dof in cluster column
!         dimA       : dimension of array A_matrices           
!         dimB       : dimension of array B_matrices           
!         ind_col    : index of the current column
!         ind_rec    : index of recursion in ACA           
!         wave_num   : wave number                 
!         alpha      : parameter to avoid the singulary of Helmholtz
!                      kernel                   
!         cnod_row   : coordinates of points corresponding to dofs
!                      in cluster row          
!         cnod_col   : coordinates of points corresponding to dofs
!                      in cluster column          
!         B_matrices : array where are stored the matrices B
!
!                                                                    
! INPUTS/OUTPUTS: A_matrices : array where are stored the matrices A
!                                                                    
! GOAL: Compute a column of the approximate error matrix
!
! BIBLIOGRAPHY: K.Zhao, "The Adaptive Cross Approximation Algorithm
!               for Accelerated Method of Moments Computations of   
!               EMC Problems" 
! 
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/24/2014                                                  
!                                                                    
! MODIFIED: VIII/29/2015 (Luca Desiderio)                             
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

interface
    function eval_kernel(ind_row,ind_col,Ynode,cnod_row,nbnod_row,id_row,id_col)
        ! Input
        REAL(kind=8),DIMENSION(2)            :: Ynode
        INTEGER                              :: ind_row,ind_col,nbnod_row,id_row,id_col
        REAL(kind=8),DIMENSION(2,nbnod_row)  :: cnod_row
        !Output
        COMPLEX(kind=8):: eval_kernel
    end function
end interface   

  !Input
  INTEGER(kind=4),INTENT(IN):: nbnod_row, nbnod_col
  INTEGER(kind=8),INTENT(IN):: dimA, dimB
  INTEGER(kind=4),INTENT(IN):: ind_col, ind_rec
  INTEGER(kind=4)           :: id_col, id_row
 
  REAL(kind=8),DIMENSION(2,nbnod_row),INTENT(IN):: cnod_row
  REAL(kind=8),DIMENSION(2,nbnod_col),INTENT(IN):: cnod_col

  COMPLEX(kind=8),DIMENSION(dimB),INTENT(IN):: B_matrices

  !Input/Output
  COMPLEX(kind=8),DIMENSION(dimA),INTENT(INOUT):: A_matrices

  !Local Variables
  INTEGER(kind=4):: ind_row
  
  COMPLEX(kind=8),DIMENSION(ind_rec-1):: a, b

  
  COMPLEX(kind=8),EXTERNAL:: ZDOTU

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!
  
  !Compute the current row of matrix B
  b=B_matrices(ind_col:(ind_rec-1)*nbnod_col:nbnod_col)
  
  !Loop on the number of columns of the block
  DO ind_row=1,nbnod_row
     
     !Evaluate of the Helmholtz kernel function
     A_matrices(ind_row+(ind_rec-1)*nbnod_row)=&
          eval_kernel(ind_row,ind_col,cnod_col(:,ind_col),cnod_row,nbnod_row,id_row,id_col)
     !Compute the current column of matrix A
     a=A_matrices(ind_row:(ind_rec-1)*nbnod_row:nbnod_row)

     !Update the current row of the approximate error matrix
     A_matrices(ind_row+(ind_rec-1)*nbnod_row)=&
          A_matrices(ind_row+(ind_rec-1)*nbnod_row)-&
          ZDOTU(ind_rec-1,a,1,b,1)

  END DO
  
  RETURN
  
END SUBROUTINE compute_Rjk
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
SUBROUTINE compute_norm_ACA(ind_rec,A_matrices,dimA,B_matrices,dimB,&
     nbnod_row,nbnod_col,norm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: compute_norm_ACA                             
!                                                                    
! INPUTS: nbnod_row  : number of dof in cluster row
!         nbnod_col  : number of dof in cluster column
!         dimA       : dimension of array A_matrices           
!         dimB       : dimension of array B_matrices           
!         ind_rec    : index of recursion in ACA           
!         A_matrices : array where are stored the matrices A
!         B_matrices : array where are stored the matrices B
!                                                                    
! INPUTS/OUTPUTS: norm : array where are stored the matrices B
!                                                                    
! GOAL: Compute the norm of approximate of a block
!
! BIBLIOGRAPHY: K.Zhao, "The Adaptive Cross Approximation Algorithm
!               for Accelerated Method of Moments Computations of   
!               EMC Problems" 
! 
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/25/2014                                                  
!                                                                    
! MODIFIED: VIII/29/2015 (Luca Desiderio)                             
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: ind_rec
  INTEGER(kind=8),INTENT(IN):: dimA, dimB
  INTEGER(kind=4),INTENT(IN):: nbnod_row, nbnod_col

  COMPLEX(kind=8),DIMENSION(dimA),INTENT(IN):: A_matrices
  COMPLEX(kind=8),DIMENSION(dimB),INTENT(IN):: B_matrices

  !Input/Output
  REAL(kind=8),INTENT(INOUT):: norm

  !Local variables
  INTEGER(kind=4):: ind_loop
  
  COMPLEX(kind=8),DIMENSION(nbnod_row):: Rjk, a
  COMPLEX(kind=8),DIMENSION(nbnod_col):: Rik, b

  !Function Variables
  REAL(kind=8),EXTERNAL:: DZNRM2
  COMPLEX(kind=8),EXTERNAL:: ZDOTU

  !!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the SUBROUTINE !!!!!!!!!!!!!!!!!!

  !Compute the current row of matrix B
  Rik=B_matrices(1+(ind_rec-1)*nbnod_col:ind_rec*nbnod_col)

  !Compute the current column of matrix A
  Rjk=A_matrices(1+(ind_rec-1)*nbnod_row:ind_rec*nbnod_row)
  
  !First updating of norm
  norm=norm+(DZNRM2(nbnod_row,Rjk,1)**2)*(DZNRM2(nbnod_col,Rik,1)**2)

  DO ind_loop=1,ind_rec-1

     !Compute the ind_loop-th column of matrix A 
     a=A_matrices(1+(ind_loop-1)*nbnod_row:ind_loop*nbnod_row)

     !Compute the ind_loop-th row of matrix B 
     b=B_matrices(1+(ind_loop-1)*nbnod_col:ind_loop*nbnod_col)

     !Update of norm
     norm=norm+2.d0*CDABS(ZDOTU(nbnod_row,a,1,Rjk,1))*&
          CDABS(ZDOTU(nbnod_col,Rik,1,b,1))

  END DO

  RETURN

END SUBROUTINE compute_norm_ACA
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
FUNCTION find_max_ACA(dim_vect,vect_ind,vect)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! FUNCTION NAME: find_max_ACA                             
!                                                                    
! INPUTS: dim_vect : dimension of the vector vect
!         vect     : array whose maximum is searched  
!         vect_ind : array contained the useless (.eq.0) and useful
!                    (.ne.0) indices  
!                                                                    
! GOAL: Compute the maximum absolute value of the elements of an array,
!       ignoring certain indices.
! 
! AUTHOR: Luca Desiderio                                             
!                                                                    
! DATE: XII/23/2014                                                  
!                                                                    
! MODIFIED: VIII/29/2015 (Luca Desiderio)                             
!                                                                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

  !Input
  INTEGER(kind=4),INTENT(IN):: dim_vect

  INTEGER(kind=4),DIMENSION(dim_vect),INTENT(IN):: vect_ind

  COMPLEX(kind=8),DIMENSION(dim_vect),INTENT(IN):: vect

  !Output
  INTEGER(kind=4):: find_max_ACA

  !Local Variables
  INTEGER(kind=4):: ind_loop

  COMPLEX(kind=8),DIMENSION(dim_vect):: vect_copy

  !Function Variables
  INTEGER(kind=4),EXTERNAL:: IZAMAX

  !!!!!!!!!!!!!!!!!!!!!!!!!!! BODY of the FUNCTION !!!!!!!!!!!!!!!!!!!

  !Copy of the array vect
  CALL ZCOPY(dim_vect,vect,1,vect_copy,1)

  DO ind_loop=1,dim_vect
     IF(vect_ind(ind_loop).EQ.0) THEN
        !Nullify the useless elements
        vect_copy(ind_loop)=(0.d0,0.d0)
     ELSE
        !Compute the absolute value of the useful elements
        vect_copy(ind_loop)=CDABS(vect_copy(ind_loop))
     END IF
  END DO
 
  !Find the position of the maximum
  find_max_ACA=IZAMAX(dim_vect,vect_copy,1)
 
  RETURN

END FUNCTION find_max_ACA
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
subroutine compute_full_ACA(cnod_row,nbnod_row,cnod_col,nbnod_col,&
     A_matrices,dimA,B_matrices,dimB,eps_ACA,&
     low_rank,max_low_rank,flag_ACA,ind_rec1,id_row,id_col,eval_kernel)
                

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    
! SUBROUTINE NAME: compute_full_ACA                          
!                                                                    
! INPUTS: nbnod_row  : number of dof in cluster row
!         nbnod_col  : number of dof in cluster column
!         dimA       : dimension of array A_matrices           
!         dimB       : dimension of array B_matrices           
!         eps_ACA    : ACA threshold                          
!         cnod_row   : coordinates of points corresponding to dofs
!                      in cluster row          
!         cnod_col   : coordinates of points corresponding to dofs
!                      in cluster column          
!
! OUTPUTS: A_matrices : array where are stored the matrices A
!          B_matrices : array where are stored the matrices B
!          low_rank   : low-rank of the approximation  
!                                                                    
! INPUTS/OUTPUTS: flag_ACA : flag whose value is TRUE if the memory
!                            is enough to the approximation, FALSE
!                            otherwise.  
!                                                                    
! GOAL: Compute the low-rank approximation of an admissible block
!       using full ACA (Adaptive Cross Approximation)  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!

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
 

  !Input
  INTEGER(kind=4),INTENT(IN):: nbnod_row, nbnod_col 
  INTEGER(kind=4),INTENT(IN):: max_low_rank,ind_rec1
  INTEGER(kind=4),INTENT(IN):: id_row, id_col
  

  
  INTEGER(kind=8),INTENT(IN):: dimA, dimB
 
  REAL(kind=8),INTENT(IN):: eps_ACA
 
  REAL(kind=8),DIMENSION(2,nbnod_row),INTENT(IN):: cnod_row
  REAL(kind=8),DIMENSION(2,nbnod_col),INTENT(IN):: cnod_col
  
  !Input/Output
  LOGICAL,INTENT(INOUT):: flag_ACA
  !Output
  INTEGER,INTENT(OUT):: low_rank
  COMPLEX(kind=8),DIMENSION(dimA),INTENT(OUT):: A_matrices
  COMPLEX(kind=8),DIMENSION(dimB),INTENT(OUT):: B_matrices

  !Local Variables
  INTEGER(kind=4):: ind_rec 
  
  INTEGER(kind=4),DIMENSION(2):: ind_pivot
  INTEGER(kind=4) :: ind_row,ind_col

  REAL(kind=8):: norm_A, norm_R
  COMPLEX(kind=8),DIMENSION(nbnod_row,nbnod_col) :: residual_matrix

 
  !Function Variables
  
  REAL(kind=8),EXTERNAL:: ZLANGE


  ! Calculate the full matrix
  DO ind_row=1,nbnod_row
    DO ind_col=1,nbnod_col
        residual_matrix(ind_row,ind_col) = eval_kernel(ind_row,ind_col,cnod_col(:,ind_col),cnod_row,nbnod_row,id_row,id_col)
    END DO
  END DO  
   
   
  
  ! Compute the Frobenius norm of the full matrix (ZLANGE from lapack)
  norm_A= zlange('F', nbnod_row, nbnod_col, residual_matrix, nbnod_row,1)
  
  ! Calculate low rank approximation  
  DO ind_rec=1,MIN(nbnod_row,nbnod_col)
     
     !Check the status of memory
     IF(ind_rec.GT.max_low_rank) THEN
        print*, ind_rec, max_low_rank
        STOP 'Memory is not enough, see subroutine ACA   ****'
     END IF
     
       
 
     !Find the maximum index in the matrix
     ind_pivot = maxloc(abs(residual_matrix))
    
     !Terminate if the pivot is 0
     IF(ABS(residual_matrix(ind_pivot(1),ind_pivot(2))).LE.10.0**(-14)) THEN
        B_matrices(1+(ind_rec-1)*nbnod_col:ind_rec*nbnod_col)=&
             (0.d0,0.d0)
        low_rank=ind_rec-1
        flag_ACA=.TRUE.
        EXIT
     END IF
     
     
     
     ! Update A
     A_matrices(1+(ind_rec-1)*nbnod_row:nbnod_row+(ind_rec-1)*nbnod_row) = residual_matrix(:,ind_pivot(2))
     ! Update B
     B_matrices(1+(ind_rec-1)*nbnod_col:nbnod_col+(ind_rec-1)*nbnod_col) = residual_matrix(ind_pivot(1),:) &
                                                                           /residual_matrix(ind_pivot(1),ind_pivot(2))

      ! update residual_matrix
      DO ind_row=1,nbnod_row
        DO ind_col=1,nbnod_col
          residual_matrix(ind_row,ind_col) = & 
                        residual_matrix(ind_row,ind_col) - &
                         A_matrices(ind_row+(ind_rec-1)*nbnod_row)*B_matrices(ind_col+(ind_rec-1)*nbnod_col)
        END DO
      END DO  
      
     ! Calculate norm of residual matrix
     norm_R= zlange('F', nbnod_row, nbnod_col, residual_matrix, nbnod_row, 1)
     !print*,'norm_R',norm_R
     ! Check if the wanted accuracy is reached 
     IF(norm_R.LE.eps_ACA*norm_A) THEN
        low_rank=ind_rec
        flag_ACA=.TRUE.
        EXIT
     END IF

     !Check the status of memory
     IF(((ind_rec+1)*nbnod_row.GT.dimA).OR.&
          ((ind_rec+1)*nbnod_col.GT.dimB)) THEN
        low_rank=ind_rec
       ! print*,'low 3',low_rank
        flag_ACA=.FALSE.
        EXIT
     END IF     
  END DO

  IF(low_rank.EQ.0)THEN
     !print*,'Full rank'
     !TODO pause
     low_rank=ind_rec
  END IF
END SUBROUTINE 
end module
