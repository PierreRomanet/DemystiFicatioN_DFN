! Bulirsch-Stoer ODE solver from Numerical Recipes
module ode_bs
  use variables, only: fric_law_ptr
  private
  public :: bsstep,rkqs
contains

      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
    USE nrtype; USE nrutil, ONLY : arth,assert_eq,cumsum,iminloc,nrerror,&
        outerdiff,outerprod,upper_triangle
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: y
    REAL(dp), DIMENSION(:), INTENT(IN) :: dydx,yscal
    REAL(dp), INTENT(INOUT) :: x
    REAL(dp), INTENT(IN) :: htry,eps
    REAL(dp), INTENT(OUT) :: hdid,hnext
    INTEGER :: nv

    INTERFACE
        SUBROUTINE derivs(x,nv,y,dydx)
        USE nrtype
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: x
        REAL(dp), DIMENSION(:), INTENT(IN) :: y
        REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx
        INTEGER :: nv
        END SUBROUTINE derivs
    END INTERFACE
    INTEGER(I4B), PARAMETER :: IMAX=9, KMAXX=IMAX-1
    REAL(dp), PARAMETER :: SAFE1=0.25_dp,SAFE2=0.7_dp,REDMAX=1.0e-5_dp,&
        REDMIN=0.7_dp,TINY=1.0e-30_dp,SCALMX=0.1_dp
    INTEGER(I4B) :: k,km,ndum
    INTEGER(I4B), DIMENSION(IMAX) :: nseq = (/ 2,4,6,8,10,12,14,16,18 /)
    INTEGER(I4B), SAVE :: kopt,kmax
    REAL(dp), DIMENSION(KMAXX,KMAXX), SAVE :: alf
    REAL(dp), DIMENSION(KMAXX) :: err
    REAL(dp), DIMENSION(IMAX), SAVE :: a
    REAL(dp), SAVE :: epsold = -1.0_dp,xnew
    REAL(dp) :: eps1,errmax,fact,h,red,scale,wrkmin,xest
    REAL(dp), DIMENSION(size(y)) :: yerr,ysav,yseq
    LOGICAL(LGT) :: reduct
    LOGICAL(LGT), SAVE :: first=.true.
    ndum=assert_eq(size(y),size(dydx),size(yscal),'bsstep')
    if (eps /= epsold) then
        hnext=-1.0e29_dp
        xnew=-1.0e29_dp
        eps1=SAFE1*eps
        a(:)=cumsum(nseq,1)
        where (upper_triangle(KMAXX,KMAXX)) alf=eps1** &
            (outerdiff(a(2:),a(2:))/outerprod(arth( &
            3.0_dp,2.0_dp,KMAXX),(a(2:)-a(1)+1.0_dp)))
        epsold=eps
        do kopt=2,KMAXX-1
            if (a(kopt+1) > a(kopt)*alf(kopt-1,kopt)) exit
        end do
        kmax=kopt
    end if
    h=htry
    ysav(:)=y(:)
    if (h /= hnext .or. x /= xnew) then
        first=.true.
        kopt=kmax
    end if
    reduct=.false.
    main_loop: do
        do k=1,kmax
            xnew=x+h
            if (xnew == x) call nrerror('step size underflow in bsstep')
            call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
            xest=(h/nseq(k))**2
            call pzextr(k,xest,yseq,y,yerr)
            if (k /= 1) then
                errmax=maxval(abs(yerr(:)/yscal(:)))
                errmax=max(TINY,errmax)/eps
                !print*,'errmax',errmax,'loc',maxloc(abs(yerr(:)/yscal(:)))

                km=k-1
                err(km)=(errmax/SAFE1)**(1.0_dp/(2*km+1))
            end if
            if (k /= 1 .and. (k >= kopt-1 .or. first)) then
                if (errmax < 1.0) exit main_loop
                if (k == kmax .or. k == kopt+1) then
                    red=SAFE2/err(km)
                    exit
                else if (k == kopt) then
                    if (alf(kopt-1,kopt) < err(km)) then
                        red=1.0_dp/err(km)
                        exit
                    end if
                else if (kopt == kmax) then
                    if (alf(km,kmax-1) < err(km)) then
                        red=alf(km,kmax-1)*SAFE2/err(km)
                        exit
                    end if
                else if (alf(km,kopt) < err(km)) then
                    red=alf(km,kopt-1)/err(km)
                    exit
                end if
            end if
        end do
        red=max(min(red,REDMIN),REDMAX)
        h=h*red
        reduct=.true.
    end do main_loop
    x=xnew
    hdid=h
    first=.false.
    kopt=1+iminloc(a(2:km+1)*max(err(1:km),SCALMX))
    scale=max(err(kopt-1),SCALMX)
    wrkmin=scale*a(kopt)
    hnext=h/scale
    if (kopt >= k .and. kopt /= kmax .and. .not. reduct) then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if (a(kopt+1)*fact <= wrkmin) then
            hnext=h/fact
            kopt=kopt+1
        end if
    end if
    END SUBROUTINE bsstep
    
! !#########################################################################################
! 
!     
!     
SUBROUTINE rkqs(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs) 
USE nrtype
USE nrutil, ONLY : assert_eq,nrerror
IMPLICIT NONE
REAL(dp), DIMENSION(:), INTENT(INOUT) :: y 
REAL(dp), DIMENSION(:), INTENT(IN) :: dydx,yscal 
REAL(dp), INTENT(INOUT) :: x
REAL(dp), INTENT(IN) :: htry,eps
REAL(dp), INTENT(OUT) :: hdid,hnext
INTEGER ::nv
INTERFACE
SUBROUTINE derivs(x,nv,y,dydx)
USE nrtype
IMPLICIT NONE
REAL(dp), INTENT(IN) :: x
REAL(dp), DIMENSION(:), INTENT(IN) :: y 
REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx 
INTEGER ::nv
END SUBROUTINE derivs
END INTERFACE
!Fifth order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and adjust stepsize. Input are the dependent variable vector y and its derivative dydx at the starting value of the independent variable x. Also input are the stepsize to be attempted htry, the required accuracy eps, and the vector yscal against which the error is scaled. y, dydx, and yscal are all of the same length. On output, y and x are replaced by their new values, hdid is the stepsize that was actually accomplished, and hnext is the estimated next stepsize. derivs is the user-supplied subroutine that computes the right-hand-side derivatives.
INTEGER(I4B) :: ndum
REAL(dp) :: errmax,h,htemp,xnew
REAL(dp), DIMENSION(size(y)) :: yerr,ytemp
REAL(dp), PARAMETER :: SAFETY=0.9_dp,PGROW=-0.2_dp,PSHRNK=-0.25_dp,&
ERRCON=1.89e-4
! The value ERRCON equals (5/SAFETY)**(1/PGROW), see use below.

ndum=assert_eq(size(y),size(dydx),size(yscal),'rkqs')

h=htry ! Set stepsize to the initial trial value. 
do
    call rkck(y,dydx,nv,x,h,ytemp,yerr,derivs) !Take a step. 
    errmax=maxval(abs(yerr(:)/yscal(:)))/eps !Evaluate accuracy.
   ! print*,'errmax',errmax,'loc',maxloc(abs(yerr(:)/yscal(:)))
    if (errmax <= 1.0) exit ! Step succeeded. 
    htemp=SAFETY*h*(errmax**PSHRNK) ! Truncation error too large, reduce stepsize. 
    h=sign(max(abs(htemp),0.1_dp*abs(h)),h) ! No more than a factor of 10. 
    xnew=x+h
    if (xnew == x) call nrerror('stepsize underflow in rkqs')
end do
if (errmax > ERRCON) then
    hnext=SAFETY*h*(errmax**PGROW)
 else
    hnext=5.0_dp*h 
end if
hdid=h
x=x+h
y(:)=ytemp(:)
END SUBROUTINE rkqs
!#########################################################################################

SUBROUTINE rkck(y,dydx,nv,x,h,yout,yerr,derivs) 
USE nrtype 
USE nrutil, ONLY : assert_eq 
IMPLICIT NONE
REAL(dp), DIMENSION(:), INTENT(IN) :: y,dydx 
REAL(dp), INTENT(IN) :: x,h
REAL(dp), DIMENSION(:), INTENT(OUT) :: yout,yerr 
INTEGER :: nv
INTERFACE
SUBROUTINE derivs(x,nv,y,dydx)
USE nrtype
IMPLICIT NONE
REAL(dp), INTENT(IN) :: x
REAL(dp), DIMENSION(:), INTENT(IN) :: y 
REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx 
INTEGER :: nv
END SUBROUTINE derivs
END INTERFACE
!Given values for N variables y and their derivatives dydx known at x, use the fifth or- der Cash-Karp Runge-Kutta method to advance the solution over an interval h and return
!the incremented variables as yout. Also return an estimate of the local truncation er- ror in yout using the embedded fourth order method. The user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.

INTEGER(I4B) :: ndum
REAL(dp), DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp 
REAL(dp), PARAMETER :: A2=0.2_dp,A3=0.3_dp,A4=0.6_dp,A5=1.0_dp,&
A6=0.875_dp,B21=0.2_dp,B31=3.0_dp/40.0_dp,B32=9.0_dp/40.0_dp,& 
B41=0.3_dp,B42=-0.9_dp,B43=1.2_dp,B51=-11.0_dp/54.0_dp,&
B52=2.5_dp,B53=-70.0_dp/27.0_dp,B54=35.0_dp/27.0_dp,&
B61=1631.0_dp/55296.0_dp,B62=175.0_dp/512.0_dp,&
B63=575.0_dp/13824.0_dp,B64=44275.0_dp/110592.0_dp,&
B65=253.0_dp/4096.0_dp,C1=37.0_dp/378.0_dp,&
C3=250.0_dp/621.0_dp,C4=125.0_dp/594.0_dp,& 
C6=512.0_dp/1771.0_dp,DC1=C1-2825.0_dp/27648.0_dp,&
DC3=C3-18575.0_dp/48384.0_dp,DC4=C4-13525.0_dp/55296.0_dp,&
DC5=-277.0_dp/14336.0_dp,DC6=C6-0.25_dp

ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),'rkck')
ytemp=y+B21*h*dydx
call derivs(x+A2*h,nv,ytemp,ak2)
ytemp=y+h*(B31*dydx+B32*ak2)
call derivs(x+A3*h,nv,ytemp,ak3)
ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
call derivs(x+A4*h,nv,ytemp,ak4)
ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
call derivs(x+A5*h,nv,ytemp,ak5) 
ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
call derivs(x+A6*h,nv,ytemp,ak6) ! Sixth step. 
yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6) !Accumulate increments with proper weights. 
yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
!Estimate error as difference between fourth and fifth order methods.
END SUBROUTINE rkck
!#########################################################################################
SUBROUTINE mmid(y,dydx,nv,xs,htot,nstep,yout,derivs) 
USE nrtype; USE nrutil, ONLY : assert_eq,swap 
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: nstep
REAL(DP), INTENT(IN) :: xs,htot
REAL(DP), DIMENSION(:), INTENT(IN) :: y,dydx
REAL(DP), DIMENSION(:), INTENT(OUT) :: yout 
INTEGER :: nv
INTERFACE 
    SUBROUTINE derivs(x,nv,y,dydx)
    USE nrtype
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: y 
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx 
    INTEGER :: nv
    END SUBROUTINE derivs
END INTERFACE
! Modified midpoint step. Dependent variable vector y and its derivative vector dydx are input at xs. Also input is htot, the total step to be taken, and nstep, the number of substeps to be used. The output is returned as yout, which need not be a distinct array from y; if it is distinct, however, then y and dydx are returned undamaged. y, dydx, and yout must all have the same length.
INTEGER(I4B) :: n,ndum
REAL(DP) :: h,h2,x
REAL(DP), DIMENSION(size(y)) :: ym,yn 
ndum=assert_eq(size(y),size(dydx),size(yout),'mmid')
h=htot/nstep ! Stepsize this trip.
ym=y
yn=y+h*dydx ! First step.
x=xs+h
call derivs(x,nv,yn,yout)  ! Will use yout for temporary storage of derivatives.
h2=2.0_dp*h
do n=2,nstep ! General step.
    call swap(ym,yn) 
    yn=yn+h2*yout
    x=x+h
    call derivs(x,nv,yn,yout) 
end do
yout=0.5_dp*(ym+yn+h*yout)  ! Last step
END SUBROUTINE mmid



!#########################################################################################
SUBROUTINE pzextr(iest,xest,yest,yz,dy)
USE nrtype; USE nrutil, ONLY : assert_eq,nrerror 
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: iest
REAL(dp), INTENT(IN) :: xest
REAL(dp), DIMENSION(:), INTENT(IN) :: yest 
REAL(dp), DIMENSION(:), INTENT(OUT) :: yz,dy
! Use polynomial extrapolation to evaluate N functions at x = 0 by fitting a polynomial
! to a sequence of estimates with progressively smaller values x = xest, and 
! corresponding function vectors yest. This call is number iest in the sequence of calls.
!  Extrapolated function values are output as yz, and their estimated error is output as dy. 
! yest, yz, and dy are arrays of length N.
INTEGER(I4B), PARAMETER :: IEST_MAX=16 
INTEGER(I4B) :: j,nv
INTEGER(I4B), SAVE :: nvold=-1 
REAL(dp) :: delta,f1,f2
REAL(dp), DIMENSION(size(yz)) :: d,tmp,q
REAL(dp), DIMENSION(IEST_MAX), SAVE :: x
REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: qcol 
nv=assert_eq(size(yz),size(yest),size(dy),'pzextr') 
if (iest > IEST_MAX) call &
    nrerror('pzextr: probable misuse, too much extrapolation') 
if (nv /= nvold) then        !Set up internal storage.
    if (allocated(qcol)) deallocate(qcol)   
    allocate(qcol(nv,IEST_MAX))
    nvold=nv
end if
x(iest)=xest  ! Save current independent variable. 
dy(:)=yest(:) 
yz(:)=yest(:)
if (iest == 1) then ! Store first estimate in first column.
    qcol(:,1)=yest(:) 
else
    d(:)=yest(:) 
    do j=1,iest-1
        delta=1.0_dp/(x(iest-j)-xest)
        f1=xest*delta
        f2=x(iest-j)*delta
        q(:)=qcol(:,j) !Propagate tableau 1 diagonal more.
        qcol(:,j)=dy(:) 
        tmp(:)=d(:)-q(:) 
        dy(:)=f1*tmp(:) 
        d(:)=f2*tmp(:) 
        yz(:)=yz(:)+dy(:)
    end do
    qcol(:,iest)=dy(:) 
end if
END SUBROUTINE pzextr
end module ode_bs

