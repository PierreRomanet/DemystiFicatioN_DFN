c     
c
c     Direct implementation of the tension T3 for mode III rupture.
c     We used the scheme from Tada and Madariaga 2011.
c     Author: Romanet Pierre (romanet@ipgp.fr)
c
      subroutine directIII(nelement,node,slip,iprec,stressT,
     &                     stressN,ier,mu)
      implicit none
c     Variables
      integer, intent(in) :: nelement, iprec, ier
       
c
      real(kind=8), intent(in), dimension(2,nelement+1) :: node
      double precision, dimension(nelement) :: slip
      double precision, intent(out), dimension(nelement) :: stressN
      double precision, intent(out), dimension(nelement) :: stressT
c
      real(kind=8), dimension(2,nelement) :: element
      real(kind=8), dimension(2,nelement) :: t
      real(kind=8) :: r1, r2, I1, I2, n_t, mu
      integer :: k, l
c
c
      real(kind=8), parameter :: pi=4*atan(1.d0)
c
c     Compute the vector t = n1e2 - n2e1 (tangential unit vector)      
         do k=1,nelement
         t(1,k) = node(1,k+1)-node(1,k)
         t(2,k) = node(2,k+1)-node(2,k)
c       
         ! Normalization
         n_t = sqrt(t(1,k)**2+t(2,k)**2)
         t(1,k) = t(1,k)/n_t
         t(2,k) = t(2,k)/n_t 
         ! Center of elements
         element(1,k) = node(1,k) + 0.5*(node(1,k+1)-node(1,k))
         element(2,k) = node(2,k) + 0.5*(node(2,k+1)-node(2,k))
         end do

c
c 
!$OMP PARALLEL 
!$OMP DO PRIVATE(k,r1,r2,I1,I2)
c     For each value of stress =>
      do k=1,nelement
      stressT(k) = 0.d0 ! Initialization
c 
c     For each value on the fault =>
      do l=1,nelement+1
c
c     Calcul of r1 and r2  
      r1 = sqrt((element(1,k)-node(1,l))**2
     &   + (element(2,k)-node(2,l))**2)
c
c     Calcul of I
      I1 = t(1,k)/r1**2 * (element(1,k)-node(1,l))
     &   + t(2,k)/r1**2 * (element(2,k)-node(2,l)) 
      if ((l-1.ge.1).and.(l.le.nelement)) then
      stressT(k) = stressT(k) - (slip(l)-slip(l-1))*I1
      else if (l.eq.1) then
      stressT(k) = stressT(k) - slip(l)*I1 
      else
      stressT(k) = stressT(k) + slip(l-1)*I1
      end if
c
      end do 
      end do
c
!$OMP END DO
!$OMP END PARALLEL 
      stressT = stressT*mu/(2*pi)
      
      end subroutine directIII
