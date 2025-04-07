cc Copyright (C) 2011: Leslie Greengard, Zydrunas Gimbutas, Vladimir Rokhlin
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date: 2011-04-15 22:14:18 -0400 (Fri, 15 Apr 2011) $
c    $Revision: 1834 $
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       this is the end of the debugging code and the beginning
c       of the actual logic subroutines for the FMM in R^2
c
c       Fortran 95 version
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine d2tstrcr(ier,z,n,nbox,
     1    nboxes,iz,laddr,nlev,center,size,
     $    ztarg,ntarg,iztarg,w,lw,lused777)
        implicit real *8 (a-h,o-z)
        integer iz(*),iztarg(*),w(*),laddr(2,*)
        real *8 z(2,*),ztarg(2,*),center(2),corners(2,4)
        real *8, allocatable :: z0(:,:),z0targ(:,:)
ccc        save
c
c        this subroutine constructs the logical structure for the 
c        fully adaptive FMM in two dimensions and stores it in the
c        array w in the form of a link-list. after that, the user 
c        can obtain the information about various boxes and lists 
c        in it by calling the entries d2tgetb, d2tgetl, d2tlinfo
c        of this subroutine (see).
c
c              note on the list conventions. 
c
c    list 1 of the box ibox - the list of all boxes with which the
c           box ibox interacts directly, including the boxes on the 
c           same level as ibox, boxes on the finer levels, and boxes 
c           on the coarser levels. obviously, list 1 is empty for any
c           box that is not childless.
c
c    list 2 of the box ibox - the list of all boxes with which the
c           box ibox interacts in the regular multipole fashion, i.e. 
c           boxes on the same level as ibox that are separated from it
c           but whose daddies are not separated from the daddy of ibox.
c
c    list 3 of the box ibox - for a childless ibox, the list of all 
c           boxes on the levels finer than that of ibox, which are 
c           separated from ibox, and whose daddys are not separated 
c           from ibox. for a box with children, list 3 is empty.
c           
c    list 4 is dual to list 3, i.e. jbox is on the list 4 of ibox if 
c           and only if ibox is on the list 3 of jbox. 
c
c    list 5 of the box ibox - the list of all boxes at the same level 
c           that are adjacent to the box ibox - the list of colleagues
c
c                            input parameters:
c
c  z - the user-specified points in the space
c  n - the number of elements in array z
c  nbox - the maximum number of points in a box on the finest level
c  lw - the amount of memory in the array w (in integer elements)
c
c                            output parameters:
c
c  ier - error return code
c    ier=0   means successful execution
c    ier=32  means that the amount lw of space in array w
c                 is insufficient
c    ier=16 means that the subroutine attempted to construct more 
c        than 199 levels of subdivision; indicates bad trouble.
c    ier=64  means that the amount lw of space in array w
c                 is severely insufficient
c  nboxes - the total number of boxes created
c  iz - the integer array addressing the particles in 
c         all boxes. 
c       explanation: for a box ibox, the particles living in
c         it are:
c
c         (z(1,j),z(2,j)),(z(1,j+1),z(2,j+1)),
c         (z(1,j+2),z(2,j+2)), . . . 
c         (z(1,j+nj-1),z(2,j+nj-1)),
c         (z(1,j+nj),z(2,j+nj)),
c
c         with j=boxes(9,ibox), and nj=boxes(10,ibox)
c        
c  iztarg - the integer array addressing the targets in 
c         all boxes. 
c
c         (ztarg(1,j),ztag(2,j)),
c         (ztarg(1,j+1),ztarg(2,j+1)),
c         (ztarg(1,j+2),ztarg(2,j+2)), . . . 
c         (ztarg(1,j+nj-1),ztarg(2,j+nj-1)),
c         (ztarg(1,j+nj),ztarg(2,j+nj)),
c
c         with j=boxes(11,ibox), and nj=boxes(12,ibox)
c
c  laddr - an integer array dimensioned (2,nlev), describing the
c         numbers of boxes on various levels of sybdivision, so that
c         the first box on level (i-1) has sequence number laddr(1,i),
c         and there are laddr(2,i) boxes on level i-1
c  nlev - the maximum level number on which any boxes have 
c         been created. the maximum number possible is 200. 
c         it is recommended that the array laddr above be 
c         dimensioned at least (2,200), in case the user underestimates
c         the number of levels required.
c  center - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c  w - the array containing all tables describing boxes, lists, etc. 
c         it is a link-list (for the most part), and can only be accessed
c         via the entries d2tgetb, d2tgetl, d2tlinfo, of this  subroutine 
c         (see below). the first lused 777 integer locations of 
c         this array should not be altered between the call to this
c         entry and subsequent calls to the entries d2tgetb, d2tgetl,
c         d2tlinfo, of this  subroutine 
c        
c  lused777 - the amount of space in the array w (in integer words)
c        that is occupied by various tables on exit from this 
c        subroutine. this space should not be altered between the
c        call to this entry and subsequent calls to entries d2tgetb,
c        d2tgetl, d2tlinfo, of this  subroutine (see below).
c  
c        . . . construct the quad-tree structure for the user-specified 
c              set of points
c
c       size of real *8 must not exceed the size of two integers
c
        if( n .lt. 1 ) then
c       number of particles less than one, abort
        ier=128
        return
        endif
c
        ier=0
c
        ninire=2
c
        iptr=1
        lptr=500
c
        iiwork=iptr+lptr
        liwork=n+ntarg+4
c     
        iboxes=iiwork+liwork
        lboxes=lw-n-5
        maxboxes=lboxes/15-1
c
c
c        if the memory is insufficient - bomb
c
        if( lw .lt. 12*(n+ntarg) ) then
        ier=64
ccc        call prinf('in d2tstrcr before d2tallb, ier=*',ier,1)
        return
        endif
c
c	 initialize the sorting index 
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,n
	iz(i)=i
	enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,ntarg
	iztarg(i)=i
	enddo
C$OMP END PARALLEL DO
c
        allocate( z0(2,n) )
        if( ntarg .gt. 0 ) then
        allocate( z0targ(2,ntarg) )
        else
        allocate( z0targ(2,1) )
        endif
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,n
	z0(1,i)=z(1,i)
	z0(2,i)=z(2,i)
	enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,ntarg
	z0targ(1,i)=ztarg(1,i)
	z0targ(2,i)=ztarg(2,i)
	enddo
C$OMP END PARALLEL DO
c
        t1=second()
C$        t1=omp_get_wtime()
        ifempty=0
        minlevel=0
        maxlevel=100
        call d2tallbem(ier,z0,n,nbox,w(iboxes),maxboxes,
     1    nboxes,iz,laddr,nlev,center,size,w(iiwork),
     $     ifempty,minlevel,maxlevel,z0targ,ntarg,iztarg)
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('time in d2tallbem=*',t2-t1,1)
c
ccc        call d2tallb(ier,z,n,nbox,w(iboxes),maxboxes,
ccc     1    nboxes,iz,laddr,nlev,center,size,w(iiwork) )
c
c        if the memory is insufficient - bomb
c
        if(ier .eq. 0) goto 1100
ccc           call prinf('in d2tstrcr after d2tallb, ier=*',ier,1)
        if(ier .eq. 4) ier=32
        return
 1100 continue
c
c       compress the array w
c   
        nn=nboxes*15
        do 1200 i=1,nn
        w(iiwork+i-1)=w(iboxes+i-1)
 1200 continue
        iboxes=iiwork
 1300 continue
        lboxes=nboxes*15+16
c
c       ... align array for real *16 storage
        lboxes=lboxes+4-mod(lboxes,4)
c
c       construct the centers and the corners for all boxes
c       in the quad-tree
c
        icenters=iboxes+lboxes
        lcenters=(nboxes*2+2)*ninire
c
        icorners=icenters+lcenters
        lcorners=(nboxes*8+2)*ninire
c
        iwlists=icorners+lcorners
        lwlists=lw-iwlists-6
c
        call prinf('lused: ccenters(k)=*', (lcenters+lcorners)/1000,1)
        call prinf('lused(k)=*', (iwlists)/1000,1)
        call d2tcentc(center,size,w(iboxes),nboxes,
     1      w(icenters),w(icorners) )
c
c       now, construct all lists for all boxes
c
        
ccc           call prinf('before d2tlsts, lwlists=*',lwlists,1)
        t1=second()
C$        t1=omp_get_wtime()
        call d2tlsts(ier,w(iboxes),nboxes,w(icorners),
     1        w(iwlists),lwlists,lused)
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('time in d2tlsts=*',t2-t1,1)
c
        lused777=lused+iwlists
ccc        call prinf('after d2tlsts, ier=*',ier,1)
        call prinf('lused(k)=*', (lused777)/1000,1)
c
c       store all pointers
c
        w(1)=nboxes
        w(2)=iboxes
        w(3)=icorners
        w(4)=icenters
        w(5)=iwlists
        w(6)=lused777
c
        w(7)=n
        w(8)=nbox
        w(9)=nlev
        w(10)=ier
        w(11)=0
        w(12)=ifempty
        w(13)=minlevel
        w(14)=maxlevel
        do i=1,200
        w(100+2*i-2)=laddr(1,i)
        w(100+2*i-1)=laddr(2,i)
        enddo
c
        return
        end
c
c
c
c
        subroutine d2tstrcrem(ier,z,n,nbox,
     1    nboxes,iz,laddr,nlev,center,size,
     $    ztarg,ntarg,iztarg,w,lw,lused777,
     $    ifempty,minlevel,maxlevel)
        implicit real *8 (a-h,o-z)
        integer iz(*),iztarg(*),w(*),laddr(2,*)
        real *8 z(2,*),ztarg(2,*),center(2),corners(2,4)
        real *8, allocatable :: z0(:,:),z0targ(:,:)
ccc        save
c
c        this subroutine constructs the logical structure for the 
c        fully adaptive FMM in two dimensions and stores it in the
c        array w in the form of a link-list. after that, the user 
c        can obtain the information about various boxes and lists 
c        in it by calling the entries d2tgetb, d2tgetl, d2tlinfo
c        of this subroutine (see).
c
c              note on the list conventions. 
c
c    list 1 of the box ibox - the list of all boxes with which the
c           box ibox interacts directly, including the boxes on the 
c           same level as ibox, boxes on the finer levels, and boxes 
c           on the coarser levels. obviously, list 1 is empty for any
c           box that is not childless.
c
c    list 2 of the box ibox - the list of all boxes with which the
c           box ibox interacts in the regular multipole fashion, i.e. 
c           boxes on the same level as ibox that are separated from it
c           but whose daddies are not separated from the daddy of ibox.
c
c    list 3 of the box ibox - for a childless ibox, the list of all 
c           boxes on the levels finer than that of ibox, which are 
c           separated from ibox, and whose daddys are not separated 
c           from ibox. for a box with children, list 3 is empty.
c           
c    list 4 is dual to list 3, i.e. jbox is on the list 4 of ibox if 
c           and only if ibox is on the list 3 of jbox. 
c
c    list 5 of the box ibox - the list of all boxes at the same level 
c           that are adjacent to the box ibox - the list of colleagues
c
c                            input parameters:
c
c  z - the user-specified points in the space
c  n - the number of elements in array z
c  nbox - the maximum number of points in a box on the finest level
c  lw - the amount of memory in the array w (in integer elements)
c
c  ifempty - ifempty=0 - remove empty boxes, ifempty=1 - keep empty boxes 
c  minlevel - minimum level of refinement
c  maxlevel - maximum level of refinement
c
c                            output parameters:
c
c  ier - error return code
c    ier=0   means successful execution
c    ier=32  means that the amount lw of space in array w
c                 is insufficient
c    ier=16 means that the subroutine attempted to construct more 
c        than 197 levels of subdivision; indicates bad trouble.
c    ier=64  means that the amount lw of space in array w
c                 is severely insufficient
c  nboxes - the total number of boxes created
c  iz - the integer array addressing the particles in 
c         all boxes. 
c       explanation: for a box ibox, the particles living in
c         it are:
c
c         (z(1,j),z(2,j)),(z(1,j+1),z(2,j+1)),
c         (z(1,j+2),z(2,j+2)), . . . 
c         (z(1,j+nj-1),z(2,j+nj-1)),
c         (z(1,j+nj),z(2,j+nj)),
c
c         with j=boxes(9,ibox), and nj=boxes(10,ibox)
c        
c  iztarg - the integer array addressing the targets in 
c         all boxes. 
c
c         (ztarg(1,j),ztag(2,j)),
c         (ztarg(1,j+1),ztarg(2,j+1)),
c         (ztarg(1,j+2),ztarg(2,j+2)), . . . 
c         (ztarg(1,j+nj-1),ztarg(2,j+nj-1)),
c         (ztarg(1,j+nj),ztarg(2,j+nj)),
c
c         with j=boxes(11,ibox), and nj=boxes(12,ibox)
c
c  laddr - an integer array dimensioned (2,nlev), describing the
c         numbers of boxes on various levels of sybdivision, so that
c         the first box on level (i-1) has sequence number laddr(1,i),
c         and there are laddr(2,i) boxes on level i-1
c  nlev - the maximum level number on which any boxes have 
c         been created. the maximim number possible is 200. 
c         it is recommended that the array laddr above be 
c         dimensioned at least (2,200), in case the user underestimates
c         the number of levels required.
c  center - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c  w - the array containing all tables describing boxes, lists, etc. 
c         it is a link-list (for the most part), and can only be accessed
c         via the entries d2tgetb, d2tgetl, d2tlinfo, of this  subroutine 
c         (see below). the first lused 777 integer locations of 
c         this array should not be altered between the call to this
c         entry and subsequent calls to the entries d2tgetb, d2tgetl,
c         d2tlinfo, of this  subroutine 
c        
c  lused777 - the amount of space in the array w (in integer words)
c        that is occupied by various tables on exit from this 
c        subroutine. this space should not be altered between the
c        call to this entry and subsequent calls to entries d2tgetb,
c        d2tgetl, d2tlinfo, of this  subroutine (see below).
c  
c        . . . construct the quad-tree structure for the user-specified 
c              set of points
c
c       size of real *8 must not exceed the size of two integers
c
        if( n .lt. 1 ) then
c       number of particles less than one, abort
        ier=128
        return
        endif
c
        ier=0
c
        ninire=2
c
        iptr=1
        lptr=500
c
        iiwork=iptr+lptr
        liwork=n+ntarg+4
c     
        iboxes=iiwork+liwork
        lboxes=lw-n-5
        maxboxes=lboxes/15-1
c
c
c        if the memory is insufficient - bomb
c
        if( lw .lt. 12*(n+ntarg) ) then
        ier=64
ccc        call prinf('in d2tstrcr before d2tallb, ier=*',ier,1)
        return
        endif
c
c	 initialize the sorting index 
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,n
	iz(i)=i
	enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,ntarg
	iztarg(i)=i
	enddo
C$OMP END PARALLEL DO
c
        allocate( z0(2,n) )
        if( ntarg .gt. 0 ) then
        allocate( z0targ(2,ntarg) )
        else
        allocate( z0targ(2,1) )
        endif
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,n
	z0(1,i)=z(1,i)
	z0(2,i)=z(2,i)
	enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)        
        do i=1,ntarg
	z0targ(1,i)=ztarg(1,i)
	z0targ(2,i)=ztarg(2,i)
	enddo
C$OMP END PARALLEL DO
c
        t1=second()
C$        t1=omp_get_wtime()
c        ifempty=0
c        minlevel=0
c        maxlevel=100
        call d2tallbem(ier,z0,n,nbox,w(iboxes),maxboxes,
     1    nboxes,iz,laddr,nlev,center,size,w(iiwork),
     $     ifempty,minlevel,maxlevel,z0targ,ntarg,iztarg)
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('time in d2tallbem=*',t2-t1,1)
c
ccc        call d2tallb(ier,z,n,nbox,w(iboxes),maxboxes,
ccc     1    nboxes,iz,laddr,nlev,center,size,w(iiwork) )
c
c        if the memory is insufficient - bomb
c
        if(ier .eq. 0) goto 1100
ccc           call prinf('in d2tstrcr after d2tallb, ier=*',ier,1)
        if(ier .eq. 4) ier=32
        return
 1100 continue
c
c       compress the array w
c   
        nn=nboxes*15
        do 1200 i=1,nn
        w(iiwork+i-1)=w(iboxes+i-1)
 1200 continue
        iboxes=iiwork
 1300 continue
        lboxes=nboxes*15+16
c
c       ... align array for real *16 storage
        lboxes=lboxes+4-mod(lboxes,4)
c
c       construct the centers and the corners for all boxes
c       in the quad-tree
c
        icenters=iboxes+lboxes
        lcenters=(nboxes*2+2)*ninire
c
        icorners=icenters+lcenters
        lcorners=(nboxes*8+2)*ninire
c
        iwlists=icorners+lcorners
        lwlists=lw-iwlists-6
c
        call prinf('lused: ccenters(k)=*', (lcenters+lcorners)/1000,1)
        call prinf('lused(k)=*', (iwlists)/1000,1)
        call d2tcentc(center,size,w(iboxes),nboxes,
     1      w(icenters),w(icorners) )
c
c       now, construct all lists for all boxes
c
        
ccc           call prinf('before d2tlsts, lwlists=*',lwlists,1)
        t1=second()
C$        t1=omp_get_wtime()
        call d2tlsts(ier,w(iboxes),nboxes,w(icorners),
     1        w(iwlists),lwlists,lused)
        t2=second()
C$        t2=omp_get_wtime()
        call prin2('time in d2tlsts=*',t2-t1,1)
c
        lused777=lused+iwlists
ccc        call prinf('after d2tlsts, ier=*',ier,1)
        call prinf('lused(k)=*', (lused777)/1000,1)
c
c       store all pointers
c
        w(1)=nboxes
        w(2)=iboxes
        w(3)=icorners
        w(4)=icenters
        w(5)=iwlists
        w(6)=lused777
c
        w(7)=n
        w(8)=nbox
        w(9)=nlev
        w(10)=ier
        w(11)=0
        w(12)=ifempty
        w(13)=minlevel
        w(14)=maxlevel
        do i=1,200
        w(100+2*i-2)=laddr(1,i)
        w(100+2*i-1)=laddr(2,i)
        enddo
c
        return
        end
c
c
c
c
        subroutine d2tnkids(box,nkids)
        implicit real *8 (a-h,o-z)
        integer box(15)
c       
        nkids=0
        do ikid=1,4
        if( box(4+ikid) .ne. 0 ) nkids=nkids+1
        enddo
        return
        end
c
c
c
c
c
        subroutine d2tprint(w,lw)
        implicit real *8 (a-h,o-z)
        integer w(*),box(15)
        real *8 center0(2),corners0(2,4)
c
c       ... retrieve stored quad-tree parameters 
c
        ibox=1
        call d2tgetb(ier,ibox,box,center0,corners0,w)
c
        call prin2('after d2tstrcr, center0=*',center0,3)
c
        size0=corners0(1,3)-corners0(1,1)
        call prin2('after d2tstrcr, size0=*',size0,1)
c
        nlev0=w(9)
        call prinf('after d2tstrcr, nlev0=*',nlev0,1)
c
        nbox0=w(8)
        call prinf('after d2tstrcr, nbox0=*',nbox0,1)        
c
        call prinf('after d2tstrcr, laddr0=*',w(100),2*(nlev0+1))
c
        return
        end
c
c
c
        subroutine d2trestore(nboxes,laddr,nlev,center,size,w,lw)
        implicit real *8 (a-h,o-z)
        integer w(*),box(15)
        real *8 center(2),center0(2),corners0(2,4)
        integer laddr(2,200)
c
c       ... retrieve stored quad-tree parameters 
c
        nboxes=w(1)
c
        do i=1,nlev+1
        laddr(1,i)=w(100+2*i-2)
        laddr(2,i)=w(100+2*i-1)
        enddo
c
        nlev=w(9)
c
        ibox=1
        call d2tgetb(ier,ibox,box,center0,corners0,w)
c
        center(1)=center0(1)
        center(2)=center0(2)
c
        size=corners0(1,3)-corners0(1,1)
c
        return
        end
c
c
c
c
        subroutine d2tgetb(ier,ibox,box,center,corners,w)
        implicit real *8 (a-h,o-z)
        integer w(*),box(15),nums(*)
        real *8 center(2),corners(2,4)
c
c        this entry returns to the user the characteristics of
c        user-specified box  ibox.  
c
c                     input parameters:
c
c  ibox - the box number for which the information is desired
c  w - storage area as created by the entry d2tstrcr (see above)
c
c                     output parameters:
c
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=4 means that ibox is either greater than the number of boxes
c           in the structure or less than 1.
c  box - an integer array dimensioned box(15). its elements describe 
c        the box number ibox, as follows:
c
c       1. level - the level of subdivision on which this box 
c             was constructed; 
c       2, 3  - the coordinates of this box among  all
c             boxes on this level
c       4 - the daddy of this box, identified by it address
c             in array boxes
c       5,6,7,8 - the  list of children of this box 
c             (eight of them, and the child is identified by its address
c             in the array boxes; if a box has only one child, only the
c             first of the four child entries is non-zero, etc.)
c       9 - the location in the array iz of the particles 
c             living in this box
c       10 - the number of particles living in this box
c       11 - the location in the array iztarg of the targets
c             living in this box
c       12 - the number of targets living in this box
c       13 - source box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       14 - target box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       15 - reserved for future use
c  center - the center of the box number ibox 
c  corners - the corners of the box number ibox 
c
c       . . . return to the user all information about the box ibox
c 
        nboxes=w(1)
        iboxes=w(2)
        icorners=w(3)
        icenters=w(4)
        iwlists=w(5)
c 
        ier=0
        if( (ibox .ge. 1)  .and. (ibox .le. nboxes) ) goto 2100
        ier=4
        return
 2100 continue
c
        ibox0=iboxes+(ibox-1)*15-1
        do 2200 i=1,15
        box(i)=w(ibox0+i)
 2200 continue
c
c      return to the user the center and the corners of the box ibox
c
        call d2tcpcc(w(icenters),w(icorners),ibox,center,corners) 
c
        return
c
c
c
c
         entry d2tgetl(ier,ibox,itype,list,nlist,w)
c
c  ibox - the box number for which the information is desired
c  itype - the type of the desired list for the box ibox
c  w - storage area as created by the entry d2tstrcr (see above)
c
c                     output parameters:
c
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=4 means that the list  itype  for the box  ibox  is empty
c  list - the list  itype  for the box  ibox 
c  nlist - the number of elements in array  list
c
c       return to the user the list number itype for the box ibox
c
        iwlists=w(5)
c
        call d2tlinkretr(ier,itype,ibox,list,nlist,w(iwlists),lused)
        return
c
c
c
c
        entry d2tlinfo(w,lused77,nums)
c
c       this entry returns to the user some of the information
c       about the storage area. namely, it returns the total 
c       amount lused777 of memory utilized in the array w (in integer 
c       locations), and the integer array nums containing the numbers
c       of elements in each of the lists
c
        iwlists=w(5)
c       
        call d2tlinkinfo(w(iwlists),lused77,nums)
        return
        end
c
c
c
c
c
        subroutine d2tlsts(ier,boxes,nboxes,corners,w,lw,lused)
        implicit real *8 (a-h,o-z)
        integer boxes(15,*),collkids(50000),w(*),
     1      dadcolls(2000),list5(20000),stack(60000)
        real *8 corners(2,4,*)
ccc        save
c
c        this subroutine constructs all lists for all boxes 
c        and stores them in the storage area w in the form
c        of a link list. the resulting data can be accessed 
c        by calls to various entries of the subroutine d2tlinkinit (see).
c
c                          input parameters:
c
c  boxes - an integer array dimensioned (15,nboxes), as created by 
c        the subroutine d2tallb (see).  each 15-element column
c         describes one box, as follows:
c
c
c       1. level - the level of subdivision on which this box 
c             was constructed; 
c       2, 3  - the coordinates of this box among  all
c             boxes on this level
c       4 - the daddy of this box, identified by it address
c             in array boxes
c       5,6,7,8 - the  list of children of this box 
c             (eight of them, and the child is identified by its address
c             in the array boxes; if a box has only one child, only the
c             first of the four child entries is non-zero, etc.)
c       9 - the location in the array iz of the particles 
c             living in this box
c       10 - the number of particles living in this box
c       11 - the location in the array iztarg of the targets
c             living in this box
c       12 - the number of targets living in this box
c       13 - source box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       14 - target box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       15 - reserved for future use
c
c  nboxes - the total number of boxes created
c  corners - the array of corners of all the boxes to in array boxes
c  lw - the total amount of storage in array w (in integer words)
c 
c              output parameters:
c
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=32 means that the amount lw of space in array w is
c           insufficient. it is a fatal error.
c  w - storage area containing all lists for all boxes in 
c        the form of link-lists, accessible by the subroutine 
c        d2tlinkretr (see).
c  lused - the amount of space in the array w (in integer words)
c        that is occupied by various tables on exit from this 
c        subroutine. this space should not be altered between the
c        call to this subroutine and subsequent calls to the 
c        entries d2tlinkretr, etc. of the subroutine d2tlinkinit (see).
c
c       . . . initialize the storage-retrieval routine for all 
c             boxes
c
ccc         call prinf('in d2tlsts, nboxes=*',nboxes,1)
c       
         lused=0
         lused2=0
c
        ier=0
        ntypes=5
        call d2tlinkinit(jer,nboxes,ntypes,w,lw)
cccc         call prinf('in d2tlsts after d2tlinkinit, ier=*',ier,1)
c
c        construct lists 5,2 for all boxes
c
        do 2000 ibox=2,nboxes
c
c       find this guy's daddy
c
        idad=boxes(4,ibox)
c
c       find daddy's collegues, including daddy himself
c
        dadcolls(1)=idad
        itype5=5
        itype2=2
        call d2tlinkretr(jer,itype5,idad,dadcolls(2),ncolls,w,lused)
        ncolls=ncolls+1
c
c        find the children of the daddy's collegues
c
        nkids=0
        do 1600 i=1,ncolls
        icoll=dadcolls(i)
        do 1400 j=1,4
        kid=boxes(4+j,icoll)
        if(kid .le. 0) goto 1600
        if(kid .eq. ibox) goto 1400
        nkids=nkids+1
        collkids(nkids)=kid
 1400 continue
 1600 continue
c
c       sort the kids of the daddy's collegues into the 
c       lists 2, 5 of the box ibox
c
        nlist1=1
        do 1800 i=1,nkids
c
c       check if this kid is touching the box ibox
c
        kid=collkids(i)
ccc        call d2tifint(corners(1,1,kid),corners(1,1,ibox),ifinter)
        call d2tifint2(boxes(1,kid),boxes(1,ibox),ifinter)
c
        if(ifinter .eq. 1)
     1    call d2tlinkstor(ier,itype5,ibox,kid,nlist1,w,lused)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
        if(ifinter .eq. 0)
     1    call d2tlinkstor(ier,itype2,ibox,kid,nlist1,w,lused)
 1800 continue
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
 2000 continue
c
ccc        call prinf('constructed lists 5,2; lused=*',lused,1)
c
c       now, construct lists 1, 3
c
        do 3000 i=1,nboxes
c
c       if this box has kids - its lists 1, 3 are empty;
c       do not construct them
c
        if(boxes(5,i) .gt. 0) goto 3000
c
c       do not construct lists 1, 3 for the main box
c
        if(boxes(1,i) .eq. 0) goto 3000
c
        call d2tlinkretr(jer,itype5,i,list5,nlist,w,lused)  
c
        if(jer .eq. 4) goto 3000
c
        do 2200 j=1,nlist
        jbox=list5(j)
        call d2tlst31(ier,i,jbox,boxes,nboxes,
     1    corners,w,stack,lused)
c
c        if storage capacity has been exceeded - bomb
c
        if(ier .eq. 32) return        
 2200 continue
 3000 continue
c
ccc        call prinf('constructed lists 1,3, lused=*',lused,1)
c
        if( 1 .eq. 2 ) then
c
c       copy all elements of lists 1, 2, 3, and 5 while skipping list 4
c       this is not needed, d2tlst31 is skipping list 4 anyway
c
        ntypes5=5
        call d2tlinkinit(jer,nboxes,ntypes5,w(lused+1),lw-(lused+5))
        do 3600 ibox=1,nboxes
        do 2400 itype=1,5
c
        call d2tlinkretr(jer,itype,ibox,list5,nlist,w,lused)
        if(jer .eq. 4) goto 2400
        call d2tlinkstor(ier,itype,ibox,list5,nlist,w(lused+1),lused2)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return
 2400 continue
 3600 continue
c
c       compress array w
c
        do 4000 i=1,lused2
        w(i)=w(lused+i)
 4000 continue
        lused=lused2
c
        endif
c
c        finally, construct the lists 4 for all boxes 
c        that need them
c
        itype3=3
        itype4=4
        nlist1=1
        do 4400 ibox=1,nboxes
c
        call d2tlinkretr(jer,itype3,ibox,list5,nlist,w,lused)
        if(jer .eq. 4) goto 4400
        do 4200 j=1,nlist
        call d2tlinkstor(ier,itype4,list5(j),ibox,nlist1,w,lused2)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
 4200 continue
 4400 continue
        lused=lused2

ccc         call prinf('exiting d2tlsts, lused=*',lused,1)
        return
        end        
c
c
c
c
c
        subroutine d2tlst31(ier,ibox,jbox0,boxes,nboxes,
     1    corners,w,stack,lused)
        implicit real *8 (a-h,o-z)
        integer w(*)
        real *8 corners(2,4,*)
        integer boxes(15,*),stack(3,*)
        data itype1/1/,itype2/2/,itype3/3/,itype4/4/,itype5/5/,
     1      nlist1/1/
ccc        save
c
c       this subroutine constructs all elements of lists 1 and 3 
c       resulting from the subdivision of one element of list 5
c       of the box ibox. all these elements of lists 1, 3 are 
c       stored in the link-lists by the subroutine linstro (see)
c
c        input parameters:
c
c  ibox - the box whose lists are being constructed
c  
c  jbox0 - the element of list 5 of the box ibox that is being
c          subdivided
c  boxes - the array boxes as created by the subroutine d2tallb (see)
c  nboxes - the number of boxes in array boxes
c  corners - the array of corners of all the boxes to in array boxes
c  w - the storage area formatted by the subroutine d2tlinkinit (see) 
c          to be used to store the elements of lists 1, 3 constructed
c          by this subroutine. obviously, by this time, it contains
c          plenty of other lists.
c  
c                      output parameters:
c
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=32 means that the amount lw of space in array w is
c           insufficient. it is a fatal error.
c  w - the augmented storage area, containing all the boxes just 
c          created, in addition to whatever had been stored previously
c  lused - the total length of array w (in integer words) used 
c          on exit from this subroutine
c
c                      work arrays:
c
c  stack - must be at least 600 integer locations long
c
c       . . . starting with the initial element of list 5, 
c             subdivide the boxes recursively and store 
c             the pieces where they belong
c
c        . . . initialize the process
c
        jbox=jbox0
        istack=1
        stack(1,1)=1
        stack(2,1)=jbox
c
        nsons=0
        do 1200 j=5,8
c
        if(boxes(j,jbox) .gt. 0) nsons=nsons+1
 1200 continue
c
        stack(3,1)=nsons
c
c       . . . move up and down the stack, generating the elements 
c             of lists 1, 3 for the box jbox, as appropriate
c
        do 5000 ijk=1,1 000 000 000
c
c       if this box is separated from ibox - store it in list 3;
c       enter this fact in the daddy's table; pass control
c       to the daddy
c
ccc        call d2tifint(corners(1,1,ibox),corners(1,1,jbox),ifinter)
        call d2tifint2(boxes(1,ibox),boxes(1,jbox),ifinter)
c
        if(ifinter .eq. 1) goto 2000
        call d2tlinkstor(ier,itype3,ibox,jbox,nlist1,w,lused)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
        istack=istack-1
        stack(3,istack)=stack(3,istack)-1
        jbox=stack(2,istack)
        goto 5000
 2000 continue
c
c       this box is not separated from ibox. if it is childless 
c       - enter it in list 1; enter this fact in the daddy's table; 
c       pass control to the daddy
c       
        if(boxes(5,jbox) .ne. 0) goto 3000
        call d2tlinkstor(ier,itype1,ibox,jbox,nlist1,w,lused)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
c
c       . . . entered jbox in the list1 of ibox; if jbox
c             is on the finer level than ibox - enter ibox
c             in the list 1 of jbox
c
        if(boxes(1,jbox) .eq. boxes(1,ibox)) goto 2400
        call d2tlinkstor(ier,itype1,jbox,ibox,nlist1,w,lused)
c
c        if storage capacity has been exceeed - bomb
c
        if(ier .eq. 32) return        
 2400 continue
c
c       if we have processed the whole box jbox0, get out
c       of the subroutine
c
        if(jbox .eq. jbox0) return
c
        istack=istack-1
        stack(3,istack)=stack(3,istack)-1
        jbox=stack(2,istack)
        goto 5000
 3000 continue
c
c       this box is not separated from ibox, and has children. if 
c       the number of unprocessed sons of this box is zero 
c       - pass control to his daddy
c
        nsons=stack(3,istack)
        if(nsons .ne. 0) goto 4000
c
        if(jbox .eq. jbox0) return
c
        istack=istack-1
        stack(3,istack)=stack(3,istack)-1
        jbox=stack(2,istack)
        goto 5000
 4000 continue
c
c       this box is not separated from ibox; it has sons, and
c       not all of them have been processed. construct the stack
c       element for the appropriate son, and pass the control
c       to him. 
c
        jbox=boxes(4+nsons,jbox)
        istack=istack+1
c
        nsons=0
        do 4600 j=5,8
c
        if(boxes(j,jbox) .gt. 0) nsons=nsons+1
 4600 continue
c
        stack(1,istack)=istack
        stack(2,istack)=jbox
        stack(3,istack)=nsons
c
 5000 continue
c
        return
        end
c
c
c
c
c
        subroutine d2tifint(c1,c2,ifinter)
        implicit real *8 (a-h,o-z)
        real *8 c1(2,4),c2(2,4)
ccc        save
c
c        this subroutine determines if two boxes in the square 
c        intersect or touch.
c
c                input parameters:
c
c  c1 - the four corners of the first box
c  c2 - the four corners of the second box
c
c                output parametes:
c 
c  ifinter - the indicator.
c      ifinter=1 means that the boxes intersect
c      ifinter=0 means that the boxes do not intersect
c
c       . . . find the maximum and minimum coordinates
c             for both boxes
c
        xmin1=c1(1,1)
        ymin1=c1(2,1)
c
        xmax1=c1(1,1)
        ymax1=c1(2,1)
c
        xmin2=c2(1,1)
        ymin2=c2(2,1)
c
        xmax2=c2(1,1)
        ymax2=c2(2,1)
c
c        xmin1=1.0d50
c        ymin1=1.0d50
c
c        xmax1=-1.0d50
c        ymax1=-1.0d50
c
c        xmin2=1.0d50
c        ymin2=1.0d50
c
c        xmax2=-1.0d50
c        ymax2=-1.0d50
c
        do 1200 i=1,4
c
        if(xmin1 .gt. c1(1,i)) xmin1=c1(1,i)
        if(ymin1 .gt. c1(2,i)) ymin1=c1(2,i)
c
        if(xmax1 .lt. c1(1,i)) xmax1=c1(1,i)
        if(ymax1 .lt. c1(2,i)) ymax1=c1(2,i)
c
c
        if(xmin2 .gt. c2(1,i)) xmin2=c2(1,i)
        if(ymin2 .gt. c2(2,i)) ymin2=c2(2,i)
c
        if(xmax2 .lt. c2(1,i)) xmax2=c2(1,i)
        if(ymax2 .lt. c2(2,i)) ymax2=c2(2,i)
c
 1200 continue
c        
c        decide if the boxes intersect
c
        eps=xmax1-xmin1
        if(eps .gt. xmax2-xmin2) eps=xmax2-xmin2
        if(eps .gt. ymax2-ymin2) eps=ymax2-ymin2
c
        if(eps .gt. ymax1-ymin1) eps=ymax1-ymin1
c
        eps=eps/10000
c
        ifinter=1
        if(xmin1 .gt. xmax2+eps) ifinter=0
        if(xmin2 .gt. xmax1+eps) ifinter=0
c
        if(ymin1 .gt. ymax2+eps) ifinter=0
        if(ymin2 .gt. ymax1+eps) ifinter=0
c
        return
        end
c
c
c
c
c
        subroutine d2tifint2(box1,box2,ifinter)
        implicit real *8 (a-h,o-z)
        integer box1(15),box2(15)
c
c        this subroutine determines if two boxes in the square 
c        intersect or touch.
c
c                input parameters:
c
c  box1 - the integer array describing the first box
c  box1 - the integer array describing the second box
c
c       integer arrays are dimensioned (15), as produced by d2tallb
c
c                output parametes:
c 
c  ifinter - the indicator.
c      ifinter=1 means that the boxes intersect
c      ifinter=0 means that the boxes do not intersect
c
c
        ifinter=1
c
        do i=1,2
        level1=box1(1)
        level2=box2(1)
        ip1=box1(i+1)-1
        ip2=box2(i+1)-1
	if( (ip1+1)*2**(level2-level1) .lt. (ip2  ) ) ifinter=0
	if( (ip1  )*2**(level2-level1) .gt. (ip2+1) ) ifinter=0
        if( ifinter .eq. 0 ) return
        enddo
c
        return
        end
c
c
c
c
c
        subroutine d2tcpcc(centers,corners,ibox,center,corner)
        implicit real *8 (a-h,o-z)
        real *8 centers(2,1),corners(2,4,1),center(2),corner(2,4)
ccc        save
c
        center(1)=centers(1,ibox)        
        center(2)=centers(2,ibox)        
c
        do 1200 i=1,4
        corner(1,i)=corners(1,i,ibox)        
        corner(2,i)=corners(2,i,ibox)        
c
 1200 continue
c
        return
        end

c
c
c
c
c
        subroutine d2tallbem(ier,z,n,nbox,boxes,maxboxes,
     1    nboxes,iz,laddr,nlev,center0,size,iwork,
     $    ifempty,minlevel,maxlevel,ztarg,ntarg,iztarg)
        implicit real *8 (a-h,o-z)
        integer boxes(15,*),iz(*),iztarg(*),laddr(2,*),iwork(*),
     1      iisons(4),jjsons(4)
        integer, allocatable :: is(:,:), ns(:,:)
        integer, allocatable :: istarg(:,:), nstarg(:,:)
        integer, allocatable :: son_idx(:)        
        real *8 z(2,1),ztarg(2,1),center0(2),center(2)
        data jjsons/1,1,2,2/,
     1      iisons/1,2,1,2/
ccc        save
c
c        this subroutine constructs a quad-tree corresponding
c        to the user-specified collection of points in the plane
c
c              input parameters:
c
c  z - the set of points in the plane
c  n - the number of elements in z
c
c  ztarg - the user-specified targets in the space
c  ntarg - the number of targets in array ztarg
c
c  nbox - the maximum number of points permitted in a box on 
c        the finest level. in other words, a box will be further
c        subdivided if it contains more than nbox points.
c  maxboxes - the maximum total number of boxes the subroutine 
c        is permitted to create. if the points z are such that 
c        more boxes are needed, the error return code ier is
c        set to 4, and the execution of the subroutine is
c        terminated.
c  
c              output parameters:
c
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=4 means that the subroutine attempted to create more 
c        than maxboxes boxes
c    ier=16 means that the subroutine attempted to construct more 
c        than 197 levels of subdivision.
c  boxes - an integer array dimensioned (15,nboxes). each 15-element
c        column describes one box, as follows:
c
c       1. level - the level of subdivision on which this box 
c             was constructed; 
c       2, 3  - the coordinates of this box among  all
c             boxes on this level
c       4 - the daddy of this box, identified by it address
c             in array boxes
c       5,6,7,8 - the  list of children of this box 
c             (eight of them, and the child is identified by its address
c             in the array boxes; if a box has only one child, only the
c             first of the four child entries is non-zero, etc.)
c       9 - the location in the array iz of the particles 
c             living in this box
c       10 - the number of particles living in this box
c       11 - the location in the array iztarg of the targets
c             living in this box
c       12 - the number of targets living in this box
c       13 - source box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       14 - target box type: 0 - empty, 1 - leaf node, 2 - sub-divided
c       15 - reserved for future use
c
c    important warning: the array boxes has to be dimensioned 
c                       at least (15,maxboxes)!! otherwise, 
c                       the subroutine is likely to bomb, since
c                       it assumes that it has that much space!!!!
c  nboxes - the total number of boxes created
c  iz - the integer array addressing the particles in 
c         all boxes. 
c       explanation: for a box ibox, the particles living in
c         it are:
c
c         (z(1,j),z(2,j)),(z(1,j+1),z(2,j+1)),
c         (z(1,j+2),z(2,j+2)), . . . 
c         (z(1,j+nj-1),z(2,j+nj-1)),
c         (z(1,j+nj),z(2,j+nj)),
c
c         with j=boxes(9,ibox), and nj=boxes(10,ibox)
c        
c  iztarg - the integer array addressing the targets in 
c         all boxes. 
c
c         (ztarg(1,j),ztag(2,j)),
c         (ztarg(1,j+1),ztarg(2,j+1)),
c         (ztarg(1,j+2),ztarg(2,j+2)), . . . 
c         (ztarg(1,j+nj-1),ztarg(2,j+nj-1)),
c         (ztarg(1,j+nj),ztarg(2,j+nj)),
c
c         with j=boxes(11,ibox), and nj=boxes(12,ibox)
c
c  laddr - an integer array dimensioned (2,numlev), containing
c         the map of array boxes, as follows:
c       laddr(1,i) is the location in array boxes of the information
c         pertaining to level=i-1
c       laddr(2,i) is the number of boxes created on the level i-1
c
c  nlev - the maximum level number on which any boxes have 
c         been created. the maximim number possible is 200. 
c         it is recommended that the array laddr above be 
c         dimensioned at least (2,200), in case the user underestimates
c         the number of levels required.
c  center0 - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c
c                  work arrays:
c
c  iwork - must be at least n+2 integer*4 elements long.
c
c
c      ... allocate work arrays
c       
        allocate( is(4,maxboxes),ns(4,maxboxes),
     $      istarg(4,maxboxes),nstarg(4,maxboxes),son_idx(maxboxes))

        do i=1,maxboxes
        ns(1,i)=0
        ns(2,i)=0
        ns(3,i)=0
        ns(4,i)=0
        nstarg(1,i)=0
        nstarg(2,i)=0
        nstarg(3,i)=0
        nstarg(4,i)=0
        enddo

c        . . . find the main box containing the whole picture
c
cccc          call prinf('in d2tallb, maxboxes=*',maxboxes,1)
c
        ier=0
        xmin=z(1,1)
        xmax=z(1,1)
        ymin=z(2,1)
        ymax=z(2,1)
c
c        xmin=1.0d50
c        xmax=-xmin
c        ymin=1.0d50
c        ymax=-ymin
c
        do 1100 i=1,n
        if(z(1,i) .lt. xmin) xmin=z(1,i)
        if(z(1,i) .gt. xmax) xmax=z(1,i)
        if(z(2,i) .lt. ymin) ymin=z(2,i)
        if(z(2,i) .gt. ymax) ymax=z(2,i)
 1100 continue
        do 1150 i=1,ntarg
        if(ztarg(1,i) .lt. xmin) xmin=ztarg(1,i)
        if(ztarg(1,i) .gt. xmax) xmax=ztarg(1,i)
        if(ztarg(2,i) .lt. ymin) ymin=ztarg(2,i)
        if(ztarg(2,i) .gt. ymax) ymax=ztarg(2,i)
 1150 continue
        size=xmax-xmin
        sizey=ymax-ymin
        if(sizey .gt. size) size=sizey
c
        center0(1)=(xmin+xmax)/2
        center0(2)=(ymin+ymax)/2
c
ccc         call prin2('in d2tallb, center0=*',center0,3)
ccc         call prin2('in d2tallb, size=*',size,1)
c
        boxes(1,1)=0
        boxes(2,1)=1
        boxes(3,1)=1
        boxes(4,1)=0
        boxes(5,1)=0     
        boxes(6,1)=0     
        boxes(7,1)=0     
        boxes(8,1)=0     
        boxes(9,1)=1
        boxes(10,1)=n
        boxes(11,1)=1
        boxes(12,1)=ntarg
        if( n .le. 0 ) boxes(13,1)=0
        if( ntarg .le. 0 ) boxes(14,1)=0
        if( n .gt. 0 ) boxes(13,1)=1
        if( ntarg .gt. 0 ) boxes(14,1)=1
        boxes(15,1)=0
c
        laddr(1,1)=1
        laddr(2,1)=1
c
        do 1200 i=1,n 
        iz(i)=i
 1200 continue
c
        do 1250 i=1,ntarg
        iztarg(i)=i
 1250 continue
c
c       recursively (one level after another) subdivide all 
c       boxes till none are left with more than nbox particles
c
        maxson=maxboxes
c
ccc         call prinf('in d2tallb, maxson=*',maxson,1) 
c
        maxlev=198
        if( maxlevel .le. maxlev ) maxlev=maxlevel

c
        ison=1
        nlev=0
cccc         call prinf('in d2tallb, nbox=*',nbox,1) 
cccc         call prinf('in d2tallb, n=*',n,1) 
        do 3000 level=0,maxlev-1
cccc          call prinf('in d2tallb, level=*',level,1) 
        laddr(1,level+2)=laddr(1,level+1)+laddr(2,level+1)
        nlevson=0
        idad0=laddr(1,level+1)
        idad1=idad0+laddr(2,level+1)-1
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idad,numpdad,numtdad)
C$OMP$PRIVATE(ii,jj,iiz,nz,iiztarg,nztarg,center)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 2000 idad=idad0,idad1
c
c       subdivide the box number idad (if needed)
c
        numpdad=boxes(10,idad)
        numtdad=boxes(12,idad)
c
c       ... refine on sources only
ccc        if(numpdad .le. nbox .and. level .ge. minlevel ) goto 2000
c
c       ... not a leaf node on sources
ccc        if( numpdad .gt. nbox ) then
ccc        if( boxes(13,idad) .eq. 1 ) boxes(13,idad)=2
ccc        if( boxes(14,idad) .eq. 1 ) boxes(14,idad)=2
ccc        endif
c
c       ... refine on targets only
cccc        if(numtdad .le. nbox .and. level .ge. minlevel ) goto 2000
c
cccc        if( numtdad .gt. nbox ) then
cccc        if( boxes(13,idad) .eq. 1 ) boxes(13,idad)=2
cccc        if( boxes(14,idad) .eq. 1 ) boxes(14,idad)=2
cccc        endif
c
c       ... refine on both sources and targets
        if(numpdad .le. nbox .and. numtdad .le. nbox .and.
     $     level .ge. minlevel ) goto 2000
c
c       
c       ... not a leaf node on sources or targets
        if( numpdad .gt. nbox .or. numtdad .gt. nbox ) then
        if( boxes(13,idad) .eq. 1 ) boxes(13,idad)=2
        if( boxes(14,idad) .eq. 1 ) boxes(14,idad)=2
        endif
c
c
        ii=boxes(2,idad)
        jj=boxes(3,idad)
        call d2tcentf(center0,size,level,ii,jj,center)
c
        iiz=boxes(9,idad)
        nz=boxes(10,idad)
c
        iiztarg=boxes(11,idad)
        nztarg=boxes(12,idad)
c
cccc        call prinf('before d2tsepa1, nz=*',nz,1)
c
        call d2tsepa1(center,z(1,iiz),iz(iiz),nz,iwork(iiz),
     1    is(1,idad),ns(1,idad))
c
        call d2tsepa1(center,ztarg(1,iiztarg),iztarg(iiztarg),
     $     nztarg,iwork(n+iiztarg),istarg(1,idad),nstarg(1,idad))
c
cccc        call prinf('after d2tsepa1, is=*',is,4)
cccc        call prinf('after d2tsepa1, ns=*',ns,4)
c
 2000   continue
C$OMP END PARALLEL DO
c
        ison_cnt=ison
        do 2010 idad=idad0,idad1 
c
        son_idx(idad)=0
        do 1610 i=1,4
        if(ns(i,idad) .eq. 0 .and. nstarg(i,idad) .eq. 0 .and. 
     $     ifempty .ne. 1) goto 1610
        if( son_idx(idad) .eq. 0 ) son_idx(idad)=ison_cnt
        nlevson=nlevson+1
        nlev=level+1
        ison_cnt=ison_cnt+1
c
c       . . . if the user-provided array boxes is too
c             short - bomb out
c
        if(ison_cnt .le. maxson) goto 1410
        ier=4
        return
 1410 continue
c
 1610   continue
 2010   continue
c
ccc        write(*,*) level, ison, ison_cnt
c
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idad)
C$OMP$PRIVATE(ii,jj,iiz,nz,iiztarg,nztarg)
C$OMP$PRIVATE(idadson,ison,i,lll,iison,jjson)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 2020 idad=idad0,idad1 
c
        ii=boxes(2,idad)
        jj=boxes(3,idad)
c
        iiz=boxes(9,idad)
        nz=boxes(10,idad)
c
        iiztarg=boxes(11,idad)
        nztarg=boxes(12,idad)
c
c       store in array boxes the sons obtained by the routine 
c       d2tsepa1
c
        if(son_idx(idad) .eq. 0 ) goto 2020
        ison=son_idx(idad)
         idadson=5
        do 1600 i=1,4
        if(ns(i,idad) .eq. 0 .and. nstarg(i,idad) .eq. 0 .and. 
     $     ifempty .ne. 1) goto 1600
        ison=ison+1
c
c        store in array boxes all information about this son
c
        do 1500 lll=5,8
        boxes(lll,ison)=0
 1500 continue
c
        boxes(1,ison)=level+1
        iison=(ii-1)*2+iisons(i)
        jjson=(jj-1)*2+jjsons(i)
        boxes(2,ison)=iison
        boxes(3,ison)=jjson
        boxes(4,ison)=idad        
c
        boxes(9,ison)=is(i,idad)+iiz-1
        boxes(10,ison)=ns(i,idad)
c
        boxes(11,ison)=istarg(i,idad)+iiztarg-1
        boxes(12,ison)=nstarg(i,idad)
c
        if( ns(i,idad) .le. 0 ) boxes(13,ison)=0
        if( nstarg(i,idad) .le. 0 ) boxes(14,ison)=0
        if( ns(i,idad) .gt. 0 ) boxes(13,ison)=1
        if( nstarg(i,idad) .gt. 0 ) boxes(14,ison)=1
        boxes(15,ison)=0
c
        boxes(idadson,idad)=ison
        idadson=idadson+1
 1600 continue
 2020 continue
C$OMP END PARALLEL DO

        ison=ison_cnt
        nboxes=ison
        laddr(2,level+2)=nlevson
         if(nlevson .eq. 0) goto 4000
         level1=level
 3000 continue
        if( level1 .ge. 197 ) ier=16
 4000 continue
        nboxes=ison
        return
        end
c
c
c
c
c
        subroutine d2tcentc(center0,size,boxes,nboxes,
     1      centers,corners)
        implicit real *8 (a-h,o-z)
        integer boxes(15,*)
        real *8 centers(2,*),corners(2,4,*),center(2),center0(2)
ccc        save
c
c       this subroutine produces arrays of centers and 
c       corners for all boxes in the quad-tree structure.
c
c              input parameters:
c
c  center0 - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c  boxes - an integer array dimensioned (15,nboxes), as produced 
c        by the subroutine d2tallb (see)
c        column describes one box, as follows:
c  nboxes - the total number of boxes created
c
c  
c              output parameters:
c
c  centers - the centers of all boxes in the array boxes
c  corners - the corners of all boxes in the array boxes
c
c       . . . construct the corners for all boxes
c
        x00=center0(1)-size/2
        y00=center0(2)-size/2
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP$PRIVATE(i,level,side,side2,ii,jj,center)   
        do 1400 i=1,nboxes
        level=boxes(1,i)
        side=size/2**level
        side2=side/2
        ii=boxes(2,i)
        jj=boxes(3,i)
        center(1)=x00+(ii-1)*side+side2
        center(2)=y00+(jj-1)*side+side2        
c
        centers(1,i)=center(1)
        centers(2,i)=center(2)
c
        corners(1,1,i)=center(1)-side/2
        corners(1,2,i)=center(1)+side/2
        corners(1,3,i)=center(1)+side/2
        corners(1,4,i)=center(1)-side/2
c
        corners(2,1,i)=center(2)-side/2
        corners(2,2,i)=center(2)-side/2
        corners(2,3,i)=center(2)+side/2
        corners(2,4,i)=center(2)+side/2
 1400 continue
C$OMP END PARALLEL DO
         return
         end
c
c
c
c
c
        subroutine d2tcentf(center0,size,level,i,j,center) 
        implicit real *8 (a-h,o-z)
        real *8 center(2),center0(2)
ccc        data level0/-1/
ccc        save
c
c       this subroutine finds the center of the box 
c       number (i,j) on the level level. note that the 
c       box on level 0 is assumed to have the center 
c       center0, and the side size
c
ccc        if(level .eq. level0) goto 1200
        side=size/2**level
        side2=side/2
        x0=center0(1)-size/2
        y0=center0(2)-size/2
        level0=level
 1200 continue
        center(1)=x0+(i-1)*side+side2
        center(2)=y0+(j-1)*side+side2
        return
        end

c
c
c
c
c
        subroutine d2tsepa1(cent,z,iz,n,iwork,
     1    is,ns)
        implicit real *8 (a-h,o-z)
        real *8 cent(2),z(2,*)
        integer iz(*),iwork(*),is(*),ns(*)
ccc        save
c
c        this subroutine reorders the particles in a box,
c        so that each of the children occupies a contigious 
c        chunk of array iz
c
c        note that we are using a standard numbering convention 
c        for the children:
c
c
c             3,4   
c                                   
c             1,2
c
c 
c                        input parameters:
c
c  cent - the center of the box to be subdivided
c  z - the list of all points in the box to be subdivided
c  iz - the integer array specifying the transposition already
c       applied to the points z, before the subdivision of 
c       the box into children
c  n - the total number of points in array z
c  
c                        output parameters:
c
c  iz - the integer array specifying the transposition already
c       applied to the points z, after the subdivision of 
c       the box into children
c  is - an integer array of length 8 containing the locations 
c       of the sons in array iz
c  ns - an integer array of length 8 containig the numbers of 
c       elements in the sons
c
c                        work arrays:
c
c  iwork - must be n integer elements long
c
c        . . . separate all particles in this box in x
c
        n1=0
        n2=0
        n3=0
        n4=0
c
        n12=0
        n34=0
c
        itype=2
        thresh=cent(2)
        call d2tsepa0ipt(z,iz,n,itype,thresh,iwork,n12)
        n34=n-n12
c
c       at this point, the contents of sons number 1,2 are in
c       the part of array iz with numbers 1,2,...n12
c       the contents of sons number 3,4  are in
c       the part of array iz with numbers n34+1,n34+2,...,n
c
c        . . . separate the boxes 1, 2 and boxes 3, 4
c
        itype=1
        thresh=cent(1)
        if(n12 .ne. 0) 
     1    call d2tsepa0ipt(z,iz,n12,itype,thresh,iwork,n1)
        n2=n12-n1
c
        if(n34 .ne. 0) 
     1    call d2tsepa0ipt
     $     (z(1,n12+1),iz(n12+1),n34,itype,thresh,iwork(n12+1),n3)
        n4=n34-n3
c
c
c       store the information about the sonnies in appropriate arrays
c
        is(1)=1
        ns(1)=n1
c
        is(2)=is(1)+ns(1)
        ns(2)=n2
c
        is(3)=is(2)+ns(2)
        ns(3)=n3
c
        is(4)=is(3)+ns(3)
        ns(4)=n4
c
c
cccc        call prinf('is as created*',is,4)
cccc        call prinf('ns as created*',ns,4)
        return
        end
c
c
c
c
c
        subroutine d2tsepa0(z,iz,n,itype,thresh,iwork,n1)
        implicit real *8 (a-h,o-z)
        real *8 z(2,*)
        integer iz(*),iwork(*)
ccc        save
c
c       subdivide the points in this box, horizontally
c       or vertically, depending on the parameter itype
c
cccc        call prin2('in d2tsepa0,thresh=*',thresh,1)
c
        i1=0
        i2=0
        do 1400 i=1,n
        j=iz(i)
        if(z(itype,j) .gt. thresh) goto 1200
        i1=i1+1
        iz(i1)=j
        goto 1400
c
 1200 continue
        i2=i2+1
        iwork(i2)=j
 1400 continue
c
        do 1600 i=1,i2
        iz(i1+i)=iwork(i)
 1600 continue
        n1=i1
        return
        end
c
c
c
c
c
        subroutine d2tsepa0ip(z,iz,n,itype,thresh,iwork,n1)
        implicit real *8 (a-h,o-z)
        real *8 z(2,*)
        integer iz(*),iwork(*)
ccc        save
c
c       subdivide the points in this box, horizontally
c       or vertically, depending on the parameter itype
c
c       in place sorting algorithm, slightly slower
c
cccc        call prin2('in d2tsepa0,thresh=*',thresh,1)
c
        i1=1
        i2=n

        do while ( i1 .le. i2 )

        do while ( (z(itype,iz(i1)) .le. thresh) ) 
        i1=i1+1
        if( i1 .gt. i2 ) exit
        enddo

        do while ( (z(itype,iz(i2)) .gt. thresh) ) 
        i2=i2-1
        if( i1 .gt. i2 ) exit
        enddo

        if( i1 .lt. i2 ) then
        i3=iz(i2)
        iz(i2)=iz(i1)
        iz(i1)=i3
        i1=i1+1
        endif

        enddo

        n1=i1-1

c        do i=1,n1
c        if( z(itype,iz(i)) .gt. thresh ) then
c        write(*,*) i1,i2,n1
c        write(*,*) i,n, '+', z(itype,iz(i))
c        pause
c        endif
c        enddo

c        do i=n1+1,n
c        if( z(itype,iz(i)) .le. thresh ) then
c        write(*,*) i1,i2,n1
c        write(*,*) i,n, '-', z(itype,iz(i))
c        pause
c        endif
c        enddo

        return
        end
c
c
c
c
c
        subroutine d2tsepa0ipt(z,iz,n,itype,thresh,iwork,n1)
        implicit real *8 (a-h,o-z)
        real *8 z(2,*)
        integer iz(*),iwork(*)
ccc        save
c
c       subdivide the points in this box, horizontally
c       or vertically, depending on the parameter itype
c
c       in place sorting algorithm, sort particles into boxes directly
c       without indirect addressing, this is much faster due to better
c       data access
c
cccc        call prin2('in d2tsepa0,thresh=*',thresh,1)
c
        i1=1
        i2=n

        do while ( i1 .le. i2 )

        do while ( (z(itype,i1) .le. thresh) ) 
        i1=i1+1
        if( i1 .gt. i2 ) exit
        enddo

        do while ( (z(itype,i2) .gt. thresh) ) 
        i2=i2-1
        if( i1 .gt. i2 ) exit
        enddo

        if( i1 .lt. i2 ) then
        i3=iz(i2)
        iz(i2)=iz(i1)
        iz(i1)=i3
        d=z(1,i2)
        z(1,i2)=z(1,i1)
        z(1,i1)=d
        d=z(2,i2)
        z(2,i2)=z(2,i1)
        z(2,i1)=d
        i1=i1+1
        endif

        enddo

        n1=i1-1

c        do i=1,n1
c        if( z(itype,iz(i)) .gt. thresh ) then
c        write(*,*) i1,i2,n1
c        write(*,*) i,n, '+', z(itype,iz(i))
c        pause
c        endif
c        enddo

c        do i=n1+1,n
c        if( z(itype,iz(i)) .le. thresh ) then
c        write(*,*) i1,i2,n1
c        write(*,*) i,n, '-', z(itype,iz(i))
c        pause
c        endif
c        enddo

        return
        end
c
c
c
c
c
        subroutine d2tlinkinit(ier,nboxes,ntypes,w,lw)
        integer w(*),nums(*),inums(20)
        data iilistad/1/,iilists/2/,inumele/3/,inboxes/4/,
     1      intypes/5/,ilw/6/,iltot/7/,
     2      inums/11,12,13,14,15,16,17,18,19,20,21,22,23,24,
     3            25,26,27,28,29,30/
ccc        save
c
c        this is the initialization entry point for the link-list
c        storage-retrieval facility. it formats the array w, to
c        be used by the entries d2tlinkstor, d2tlinkretr, d2tlinkrem,
c        d2tlinkinfo below. 
c
c                     input parameters:
c
c  nboxes - the number of boxes for which the various types of lists
c        will be stored
c  ntypes - the number of types of lists that will be stored
c  lw - the total amount of space in the area w to be used for storage
c        (in integer locations)
c
c                     output parameters:
c
c  ier - error return code;
c    ier=0 means successful execution
c    ier=1024 means that the amount of space in array w is grossly 
c        insufficient
c  w - the formatted area for storage
c
c
c       . . . allocate memory for the storage facility
c
        ier=0
c
        ilistadd=32
        nlistadd=nboxes*ntypes+10
c
        ilists=ilistadd+nlistadd
        numele=0
        ltot=ilists+numele*3+10
c
         if(ltot+100 .lt. lw) goto 1200
         ier=1024
         return
 1200 continue
c
         do 1400 i=1,20
         w(inums(i))=0
 1400 continue
c
        w(iilistad)=ilistadd
        w(iilists)=ilists
        w(inumele)=numele
        w(inboxes)=nboxes
        w(intypes)=ntypes      
        w(ilw)=lw
        w(iltot)=ltot
c
        call d2tlinkini0(w(ilistadd),nboxes,ntypes)
        return
c
c
c
c
        entry d2tlinkstor(ier,itype,ibox,list,nlist,w,lused)
c
c       this entry stores dynamically a list of positive numbers 
c       in the storage array w.
c
c                      input parameters:
c
c  itype - the type of the elements being stored
c  ibox - the box to which these elements corresponds
c  list - the list of positive integer elements to be stored
c  nlist - the number of elements in the array list
c  w - the storage area used by this subroutine; must be first 
c          formatted by the entry d2tlinkinit of this subroutine
c          (see above)
c
c                      output parameters:
c
c  ier - error return code; 
c         ier=0 means successful execution
c         ier=32 means that the storage area w would be exceeded
c                 by this storage request
c  w - the storage area used by this subroutine
c  lused - the number of integer elements used in array
c          w after this call.
c
c       . . . if this storage request exceeds the available memory - bomb
c
        ier=0
        if(w(iilists)+w(inumele)*3+nlist*3 .lt. w(ilw) ) goto 2200
        ier=32
        return
 2200 continue
c
c       store the user-specified list in array w
c
c       Old versions of gfortran <=4.3 will naively optimize w(inumele)
c       away if w(inumele)=0, at optimization levels greater than 1.
c       We use a temporary variable to fix this nasty compiler bug.
c
        itemp = w(inumele)
        call d2tlinksto0(itype,ibox,list,nlist,w(w(iilistad)),
     1      w(inboxes),w(w(iilists)),itemp,w(inums(itype)) )
        w(inumele)=itemp
c
ccc        call d2tlinksto0(itype,ibox,list,nlist,w(w(iilistad)),
ccc     1      w(inboxes),w(w(iilists)),w(inumele),w(inums(itype)) )
c
c       augment the amount of storage used 
c        
        lused=w(iilists)+w(inumele)*3+10
        return
c
c
c
c
        entry d2tlinkretr(ier,itype,ibox,list,nlist,w,lused)
c
c       this entry retrieves from the storage area  w 
c       a list of positive numbers that has been stored there
c       by the entry d2tlinkstor (see above).
c
c                      input parameters:
c
c  itype - the type of the elements to be retrieved
c  ibox - the box to which these elements correspond
c  w - the storage area from which the information is to be 
c          retrieved
c
c                      output parameters:
c
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4 means that no elements are present of the type
c                  itype and the box ibox
c  list - the list of positive integer elements retrieved 
c  nlist - the number of elements in the array list
c  lused - the number of integer elements used in array
c          w after this call.
c
        call d2tlinkret0(ier,itype,ibox,w(w(iilistad)),
     1      w(w(iilists)),list,w(inboxes),nlist)
c
        lused=w(iilists)+w(inumele)*3+10
c
        return
c
c
c
c
        entry d2tlinkrem(ier,itype,ibox,list,nlist,w,lused)
c
c       this entry deletes elements in array lists corresponding 
c       the user-specified list  list. actually, it does not 
c       delete anything, but rather marks these elements for
c       future destruction by making them negative.
c
c                      input parameters:
c
c  itype - the type of the elements to be destroyed
c  ibox - the box to which these elements correspond
c  list - the list of positive integer elements to be destroyed
c  nlist - the number of elements in the array list
c  w - the storage area from which the information is to be 
c          retrieved
c
c                      output parameters:
c
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4*k means that k of the elements on the user-specified
c                   list  list were not present
c          ier=22 means that  no elements whatsoever were present 
c                   for the type  itype   and the box  ibox
c  w - the storage area from which the information is to be 
c          retrieved
c  lused - the number of integer elements used in array
c          w both before and after this call.
c
c       mark for destruction the user-specified elements 
c
        call d2tlinkrem0(ier,itype,ibox,list,nlist,w(w(iilistad)),
     1      w(inboxes),w(w(iilists)),w(inums(itype)) )
c
        lused=w(iilists)+w(inumele)*3+10
        return
c
c
c
c
        entry d2tlinkinfo(w,lused,nums)
c
c       this entry returns to the user some of the information
c       about the storage area. namely, it returns the total 
c       amount lused of memory utilized in the array w (in integer 
c       locations), and the integer array nums containing the numbers
c       of elements in each of the lists
c
        lused=w(iilists)+w(inumele)*3+10
        ntypes7=w(intypes)
         call prinf('in d2tlinkinfo, lused=*',lused,1)
         call prinf('in d2tlinkinfo, ntypes7=*',ntypes7,1)
         call prinf('in w(inumele)=*',w(inumele),1)
        do 6200 i=1,ntypes7
        nums(i)=w(inums(i))
 6200 continue
        return
        end
c
c
c
c
c
        subroutine d2tlinksto0(itype,ibox,list,nlist,listaddr,
     1      nboxes,lists,numele,numtype)
        integer listaddr(nboxes,*),lists(3,*),list(*)
ccc        save
c
c       this entry stores dynamically a list of positive numbers 
c       in the storage array lists, while entering the information
c       about this event in the array listaddr.
c
c                      input parameters:
c
c  itype - the type of the elements being stored
c  ibox - the box to which these elements corresponds
c  list - the list of positive integer elements to be stored
c  nlist - the number of elements in the array list
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c  nboxes - the total number of boxes indexing the elements 
c           in array lists
c  lists - the main storage area used by this subroutine
c  numele - the number of elements stored in array lists on entry 
c             to this subroutinec

c                      output parameters:
c
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c  lists - the main storage area used by this subroutine
c  numele - the number of elements stored in array lists on exit
c             from this subroutine
c  numtype - the total number of elements in all lists of the 
c           type itype AFTER this call
c
c .........................................................................
c
c        interpretation of the entries in arrays lists, listaddr:
c
c  lists(1,i) - the location in array lists of the preceding
c        element in the list with (box,type) combination as the 
c        user-supplied (ibox,itype)
c     lists(1,i) .leq. 0 means that this is the first element 
c        of its type,
c  lists(2,i) - the box number on the interaction list.
c
c
c  listaddr(ibox,itype) - the address in the array lists of the last element
c        of this type for this box;
c     listaddr(ibox,itype) .leq. 0 means that there are no elements 
c        in the list of this type for this box.
c
c       . . . store the user-supplied list elements in the array lists,
c             and enter information about this change in array listaddr
c
        ilast=listaddr(ibox,itype)
        do 1200 i=1,nlist
        numele=numele+1
        numtype=numtype+1
        lists(1,numele)=ilast
        lists(2,numele)=list(i)
c
        ilast=numele
 1200 continue
        listaddr(ibox,itype)=ilast
c
        return
c
c
c
c
        entry d2tlinkret0(ier,itype,ibox,listaddr,lists,list,
     1      nboxes,nlist)
c
c       this entry retrieves from the main storage array lists
c       a list of positive numbers that has been stored there
c       by the entry d2tlinksto0 (see above).
c
c                      input parameters:
c
c  itype - the type of the elements being to be retrieved
c  ibox - the box to which these element corresponds
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c  nboxes - the total number of boxes indexing the elements 
c           in array lists
c  lists - the main storage area used by this subroutine
c
c                      output parameters:
c
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4 means that no elements are present of the type
c                  itype and the box ibox
c  list - the list of positive integer elements to be stored
c  nlist - the number of elements in the array list
c
c
c       . . . retrieve and store in array list the list of the 
c             type itype for the box ibox
c
        ier=0
        ilast=listaddr(ibox,itype)
        if(ilast .gt. 0) goto 2200
        nlist=0
        ier=4
        return
 2200 continue
c
        nlist=0
        do 2400 i=1,1 000 000 000
c
        if(lists(2,ilast) .le. 0) goto 2300       
        nlist=nlist+1
        list(nlist)=lists(2,ilast)
 2300 continue
        ilast=lists(1,ilast)
        if(ilast .le. 0) goto 2600
 2400 continue
 2600 continue
c
        if(nlist .gt. 0) goto 2650
        ier=4
        return
 2650 continue
c
c        flip the retrieved array
c
        if(nlist .eq. 1) return
        do 2700 i=1,nlist/2
        j=list(i)
        list(i)=list(nlist-i+1)
        list(nlist-i+1)=j
 2700 continue
        return
c
c
c
c
        entry d2tlinkini0(listaddr,nboxes,ntypes)
c
c       this subroutine initializes the array listaddr to be used 
c       later by other entries of this subroutine
c
c                        input parameters:
c
c  nboxes - the number of boxes for which the various types of lists
c        will be stored
c  ntypes - the number of types of lists that will be stored
c
c                        output parameters:
c
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c
c       . . . initialize the array listaddr
c
        do 3000 k=1,ntypes
        do 2800 i=1,nboxes
        listaddr(i,k)=-1
 2800 continue
 3000 continue
        return
c
c
c
c
        entry d2tlinkrem0(ier,itype,ibox,list,nlist,listaddr,
     1      nboxes,lists,numtype)
c
c       this entry deletes elements in array lists corresponding 
c       the user-specified list  list. actually, it does not 
c       delete anything, but rather marks these elements for
c       future destruction by making them negative.
c
c                      input parameters:
c
c  itype - the type of the elements to be destroyed
c  ibox - the box to which these elements correspond
c  list - the list of positive integer elements to be destroyed
c  nlist - the number of elements in the array list
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c  nboxes - the total number of boxes indexing the elements 
c           in array lists
c  lists - the main storage area used by this subroutine
c
c                      output parameters:
c
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4*k means that k of the elements on the user-specified
c                   list  list were not present
c          ier=22 means that  no elements whatsoever were present 
c                   for the type  itype   and the box  ibox
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call 
c             to the entry d2tlinkini0 of this subroutine (see below).
c  lists - the main storage area used by this subroutine
c
c       . . . mark for destruction the elements of type (itype,ibox)
c             that are in the list  list.
c
        ier=0
        do 4000 i=1,nlist 
        ilast=listaddr(ibox,itype)
        if(ilast. gt. 0) goto 3200
        ier=22
        return
 3200 continue
c
        iffound=0
        do 3600 j=1,1 000 000 000
c
        if(ilast .le. 0) goto 3800
        if(lists(2,ilast) .ne. list(i)) goto 3400
        lists(2,ilast)=-lists(2,ilast)
        numtype=numtype-1
        iffound=1
 3400 continue
        ilast=lists(1,ilast)
 3600 continue
 3800 continue
         if(iffound .eq. 0) ier=ier+4
c
 4000 continue
        return
        end
c
c
c
c
c
        subroutine d2tgetbbox(n,z,center,size,corners)
        implicit real *8 (a-h,o-z)
        real *8 z(2,*),center(2),corners(2,4)
c
c       this subroutine finds the center and 
c       corners for the top level box in the quad-tree structure
c
c
        xmin=z(1,1)
        xmax=z(1,1)
        ymin=z(2,1)
        ymax=z(2,1)
c
c        xmin=1.0d50
c        xmax=-xmin
c        ymin=1.0d50
c        ymax=-ymin
c
        do 1100 i=1,n
        if(z(1,i) .lt. xmin) xmin=z(1,i)
        if(z(1,i) .gt. xmax) xmax=z(1,i)
        if(z(2,i) .lt. ymin) ymin=z(2,i)
        if(z(2,i) .gt. ymax) ymax=z(2,i)
 1100 continue
        size=xmax-xmin
        sizey=ymax-ymin
        if(sizey .gt. size) size=sizey
c
c
        center(1)=(xmin+xmax)/2
        center(2)=(ymin+ymax)/2
c
c
        corners(1,1)=center(1)-size/2
        corners(1,2)=corners(1,1)+size
        corners(1,3)=corners(1,1)+size
        corners(1,4)=corners(1,1)
c
c
        corners(2,1)=center(2)-size/2
        corners(2,2)=corners(2,1)
        corners(2,3)=corners(2,1)+size
        corners(2,4)=corners(2,3)+size
c
        return
        end

cc Copyright (C) 2010: Leslie Greengard, Zydrunas Gimbutas, Vladimir Rokhlin
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the FMM tree plotting routines in R^2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine plot_points2d(iw,z,n)
        implicit real *8 (a-h,o-z)
        real *8 z(2,*)
c
        do 1200 i=1,n
        write(iw,1000) z(1,i),z(2,i)
 1200   continue
c
 1000   format(6(1x,e11.5))
c
        return
        end
c
c
c
c
c
        subroutine plot_box2d(iw,center,size)
        implicit real *8 (a-h,o-z)
        real *8 center(2)
c
        write(iw,1000) 
     $     center(1)-size/2,center(2)-size/2
        write(iw,1000) 
     $     center(1)-size/2,center(2)+size/2
        write(iw,1000) 
     $     center(1)+size/2,center(2)+size/2
        write(iw,1000) 
     $     center(1)+size/2,center(2)-size/2
        write(iw,1000) 
     $     center(1)-size/2,center(2)-size/2
        write(iw,1200)
c
 1000   format(6(1x,e11.5))
 1200   format(80a1)
        return
        end 
c
c
c
c
c
        subroutine plot_label2d(iw,center,size,itag,label)
        implicit real *8 (a-h,o-z)
        real *8 center(2)
c        
        if( label .ge. 1 .and. label .lt. 10 )
     $     write(iw,1010) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
        if( label .ge. 10 .and. label .lt. 100 )
     $     write(iw,1020) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
        if( label .ge. 100 .and. label .lt. 1000 )
     $     write(iw,1030) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
        if( label .ge. 1000 .and. label .lt. 10000 )
     $     write(iw,1040) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
        if( label .ge. 10000 .and. label .lt. 100000 )
     $     write(iw,1050) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
        if( label .ge. 100000 )
     $     write(iw,*) 'set label ', itag, ' "', label, 
     $     '" at ', center(1), ', ', center(2), ' center'  
c
 1010   format(a11,i7,a3,i1,a5,e13.7,a3,e13.7,a7)
 1020   format(a11,i7,a3,i2,a5,e13.7,a3,e13.7,a7)
 1030   format(a11,i7,a3,i3,a5,e13.7,a3,e13.7,a7)
 1040   format(a11,i7,a3,i4,a5,e13.7,a3,e13.7,a7)
 1050   format(a11,i7,a3,i5,a5,e13.7,a3,e13.7,a7)
 1060   format(a11,i7,a3,i6,a5,e13.7,a3,e13.7,a7)

        return
        end
cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date: 2011-02-22 17:34:23 -0500 (Tue, 22 Feb 2011) $
c    $Revision: 1670 $
c
c
c     Computation of  Bessel functions via recurrence
c
c**********************************************************************
      subroutine jfuns2d(ier,nterms,z,scale,fjs,ifder,fjder,
     1	      lwfjs,iscale,ntop)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c PURPOSE:
c
c	This subroutine evaluates the first NTERMS  Bessel 
c	functions and if required, their derivatives.
c	It incorporates a scaling parameter SCALE so that
c       
c		fjs_n(z)=j_n(z)/SCALE^n
c		fjder_n(z)=\frac{\partial fjs_n(z)}{\partial z}
c
c	NOTE: The scaling parameter SCALE is meant to be used when
c             abs(z) < 1, in which case we recommend setting
c	      SCALE = abs(z). This prevents the fjs_n from 
c             underflowing too rapidly.
c	      Otherwise, set SCALE=1.
c	      Do not set SCALE = abs(z) if z could take on the 
c             value zero. 
c             In an FMM, when forming an expansion from a collection of
c             sources, set SCALE = min( abs(k*r), 1)
c             where k is the Helmholtz parameter and r is the box dimension
c             at the relevant level.
c
c INPUT:
c
c    nterms (integer): order of expansion of output array fjs 
c    z     (complex *16): argument of the  Bessel functions
c    scale    (real *8) : scaling factor (discussed above)
c    ifder  (integer): flag indicating whether to calculate "fjder"
c		          0	NO
c		          1	YES
c    lwfjs  (integer): upper limit of input arrays 
c                         fjs(0:1) and iscale(0:1)
c    iscale (integer): integer workspace used to keep track of 
c                         internal scaling
c
c OUTPUT:
c
c    ier    (integer): error return code 
c                         ier=0 normal return;
c                         ier=8 insufficient array dimension lwfjs
c    fjs   (complex *16): array of scaled Bessel functions.
c    fjder (complex *16): array of derivs of scaled Bessel functions.
c    ntop  (integer) : highest index in arrays fjs that is nonzero
c
c       NOTE, that fjs and fjder arrays must be at least (nterms+2)
c       complex *16 elements long.
c
c
      integer iscale(0:lwfjs)
      complex *16 wavek,fjs(0:lwfjs),fjder(0:lwfjs)
      complex *16 z,zinv,com,fj0,fj1,zscale,ztmp
c
      complex *16 psi,zmul,zsn,zmulinv
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
c
      data upbound/1.0d+32/, upbound2/1.0d+40/, upbound2inv/1.0d-40/
      data tiny/1.0d-200/,done/1.0d0/,zero/0.0d0/
c
c ... Initializing ...
c
      ier=0
c
c       set to asymptotic values if argument is sufficiently small
c
      if (abs(z).lt.tiny) then
         fjs(0) = done
         do i = 1, nterms
            fjs(i) = zero
	 enddo
c
	 if (ifder.eq.1) then
	    do i=0,nterms
	       fjder(i)=zero
	    enddo
	    fjder(1)=done/(2*scale)
	 endif
c
         RETURN
      endif
c
c ... Step 1: recursion up to find ntop, starting from nterms
c
      ntop=0
      zinv=done/z
      fjs(nterms)=done
      fjs(nterms-1)=zero
c
      do 1200 i=nterms,lwfjs
         dcoef=2*i
         ztmp=dcoef*zinv*fjs(i)-fjs(i-1)
         fjs(i+1)=ztmp
c
         dd = dreal(ztmp)**2 + dimag(ztmp)**2
         if (dd .gt. upbound2) then
            ntop=i+1
            goto 1300
         endif
 1200 continue
 1300 continue
      if (ntop.eq.0) then
         ier=8
         return
      endif
c
c ... Step 2: Recursion back down to generate the unscaled jfuns:
c             if magnitude exceeds UPBOUND2, rescale and continue the 
c	      recursion (saving the order at which rescaling occurred 
c	      in array iscale.
c
      do i=0,ntop
         iscale(i)=0
      enddo
c
      fjs(ntop)=zero
      fjs(ntop-1)=done
      do 2200 i=ntop-1,1,-1
	 dcoef=2*i
         ztmp=dcoef*zinv*fjs(i)-fjs(i+1)
         fjs(i-1)=ztmp
c
         dd = dreal(ztmp)**2 + dimag(ztmp)**2
         if (dd.gt.UPBOUND2) then
            fjs(i) = fjs(i)*UPBOUND2inv
            fjs(i-1) = fjs(i-1)*UPBOUND2inv
            iscale(i) = 1
         endif
 2200 continue
c
c ...  Step 3: go back up to the top and make sure that all
c              Bessel functions are scaled by the same factor
c              (i.e. the net total of times rescaling was invoked
c              on the way down in the previous loop).
c              At the same time, add scaling to fjs array.
c
      ncntr=0
      scalinv=done/scale
      sctot = 1.0d0
      do i=1,ntop
         sctot = sctot*scalinv
         if(iscale(i-1).eq.1) sctot=sctot*UPBOUND2inv
         fjs(i)=fjs(i)*sctot
      enddo
c
c ... Determine the normalization parameter:
c
c     From Abramowitz and Stegun (9.1.47) and (9.1.48), Euler's identity
c
        psi = 0d0
c
        if (dimag(z) .lt. 0) zmul = +ima
        if (dimag(z) .ge. 0) zmul = -ima
        zsn = zmul**(mod(ntop,4))
c
        zmulinv=1/zmul
        do i = ntop,1,-1
           psi = scale*psi+fjs(i)*zsn
           zsn = zsn*zmulinv
        enddo
        psi = 2*psi*scale+fjs(0)
c
        if (dimag(z) .lt. 0) zscale = cdexp(+ima*z) / psi
        if (dimag(z) .ge. 0) zscale = cdexp(-ima*z) / psi
c
c
c ... Scale the jfuns by zscale:
c
      ztmp=zscale
      do i=0,nterms
         fjs(i)=fjs(i)*ztmp
      enddo
c
c ... Finally, calculate the derivatives if desired:
c
      if (ifder.eq.1) then
         fjs(nterms+1)=fjs(nterms+1)*ztmp
c
         fjder(0)=-fjs(1)*scale
         do i=1,nterms
            dc1=0.5d0
            dc2=done-dc1
            dc1=dc1*scalinv
            dc2=dc2*scale
            fjder(i)=dc1*fjs(i-1)-dc2*fjs(i+1)
         enddo
      endif
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Printing routines
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C
        SUBROUTINE PRINI(IP1,IQ1)
        CHARACTER *1 MES(1), AA(1)
         save
        REAL *4 A(1)
        REAL *8 A2(1)
        REAL *8 A4(1)
ccc        INTEGER *4 IA(1)
        INTEGER IA(1)
        INTEGER *4 IA1(1)
        INTEGER *2 IA2(1)
        data IP/0/,IQ/0/
        IP=IP1
        IQ=IQ1
        RETURN

C
C
C
C
C
        ENTRY PRIN(MES,A,N)
        CALL  MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
 1200 FORMAT(6(2X,E11.5))
         RETURN
C
C
C
C
        ENTRY PRIN2(MES,A2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
        RETURN
C
C
C
C
        ENTRY PRINQ(MES,A4,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1500)(A4(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1500)(A4(J),J=1,N)
 1500 FORMAT(6(2X,e11.5))
        RETURN
C
C
C
C
        ENTRY PRINF(MES,IA,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
C
C
C
C
        ENTRY PRINF1(MES,IA1,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA1(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA1(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY PRINF2(MES,IA2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY PRINA(MES,AA,N)
        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
        END
c
c
c
c
c
        SUBROUTINE MESSPR(MES,IP,IQ)
        CHARACTER *1 MES(1),AST
        DATA AST/'*'/
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
         RETURN
         END
c
c
c
c
c
        subroutine msgmerge(a,b,c)
        character *1 a(1),b(1),c(1),ast
        data ast/'*'/
c
        do 1200 i=1,1000
c
        if(a(i) .eq. ast) goto 1400
        c(i)=a(i)       
        iadd=i
 1200 continue
c
 1400 continue
c
        do 1800 i=1,1000
c
        c(iadd+i)=b(i)
        if(b(i) .eq. ast) return
 1800 continue
        return
        end
c
cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c       
c     This file contains the main FMM routines and some related
c     subroutines for evaluating Cauchy sums due to
c     point charges and dipoles.  (FORTRAN 90 VERSION)
c
c       cfmm2d does the Cauchy type sums of complex valued
c       charges and dipoles - dipole vectors have no meaning in cfmm2d.
c
c cfmm2d: charge and dipstr are complex valued, z are complex numbers,
c       generalized Cauchy sums.
c
c        Note, that the complex valued logarithm is a multi-valued
c        function, so the potential values have to be interpreted
c        carefully, if charges are specified.  For example, only the
c        real part of potential is meaningful for real valued charges.
c        The gradients and hessians are valid for arbitrary complex charges.
c
c \phi(z_i) = \sum_{j\ne i} charge_j *\log(z_i-z_j) + dipstr_j *(1/(z_i-z_j))
c
c        In this routine, we define the gradient as the first
c        derivative with respect to z, and the hessian as the second
c        derivative with respect to z.
c
c \gradient \phi(z_i) = \frac{\partial \phi(z_i)}{\partial z}
c \hessian  \phi(z_i) = \frac{\partial^2 \phi(z_i)}{\partial z^2}
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine cfmm2dparttarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c       
c       Generalized Cauchy FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets.
c
c       We use log(z) for the Green's function.
c       Self-interactions are not included.
c   
c cfmm2d: charge and dipstr are complex valued, z are complex numbers.
c
c        Note, that the complex valued logarithm is a multi-valued
c        function, so the potential values have to be interpreted
c        carefully, if charges are specified.  For example, only the
c        real part of potential is meaningful for real valued charges.
c        The gradients and hessians are valid for arbitrary complex charges.
c
c \phi(z_i) = \sum_{j\ne i} charge_j \log(z_i-z_j) + dipstr_j \frac{1}{z_i-z_j}
c
c        In this routine, we define the gradient as the first
c        derivative with respect to z, and the hessian as the second
c        derivative with respect to z.
c
c \gradient \phi(z_i) = \frac{\partial \phi(z_i)}{\partial z}
c \hessian  \phi(z_i) = \frac{\partial^2 \phi(z_i)}{\partial z^2}
c
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine cfmm2dparttargmain.
c
c       INPUT PARAMETERS:
c
c       iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c                  5 => tolerance =.5d-15
c
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: complex *16 (nsource): dipole strengths
c       dipvec: real *8 (2,nsource): dipole orientation vectors. 
c
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       ifgrad:  gradient flag (1=compute gradient, otherwise no)
c       ifhess:  hessian flag (1=compute hessian, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (2,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       ifgradtarg:  target gradient flag 
c                   (1=compute gradient, otherwise no)
c       ihesstarg:  target hessian flag 
c                   (1=compute hessian, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       ier   =  error return code
c                ier=0     =>  normal execution
c                ier=4     =>  cannot allocate tree workspace
c                ier=8     =>  cannot allocate bulk FMM  workspace
c                ier=16    =>  cannot allocate mpole expansion
c                              workspace in FMM
c
c       pot: complex *16 (nsource): potential at source locations
c       grad: complex *16 (nsource): gradient  at source locations
c       hess: complex *16 (nsource): hessian at source locations
c       pottarg: complex *16 (ntarget): potential at target locations 
c       gradtarg: complex *16 (ntarget): gradient  at target locations 
c       hesstarg: complex *16 (ntarget): hessian at target locations
c
c
c
        real *8 source(2,*)
        complex *16 charge(*)
        complex *16 dipstr(*)
        complex *16 ima
        complex *16 pot(*)
        complex *16 grad(*)
        complex *16 hess(*)
c
        real *8 target(2,*)
        complex *16 pottarg(*)
        complex *16 gradtarg(*)        
        complex *16 hesstarg(*)
c
        real *8 timeinfo(10)
c       
        real *8 center(2)
c       
c
c     Note: various arrays dimensioned here to 200.
c     That allows for 200 levels of refinement, which is 
c     more than enough for any non-pathological case.
c
        integer laddr(2,1000)     ! Change to 1000
        real *8 scale(0:1000)
        real *8 bsize(0:1000)
        integer nterms(0:1000)
c       
        complex *16 ptemp,gtemp(2),htemp(3)
c       
        integer box(15)
        real *8 center0(2),corners0(2,4)
c       
        integer box1(15)
        real *8 center1(2),corners1(2,4)
c       
        real *8, allocatable :: w(:)
        real *8, allocatable :: wlists(:)
        real *8, allocatable :: wrmlexp(:)
c
        data ima/(0.0d0,1.0d0)/
c       
        ier=0
        lused7=0
c       
        done=1
        pi=4*atan(done)
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c       
        ifprint=1
c
c     set fmm tolerance based on iprec flag.
c       
        if( iprec .eq. -2 ) epsfmm=.5d-0 
        if( iprec .eq. -1 ) epsfmm=.5d-1
        if( iprec .eq. 0 ) epsfmm=.5d-2
        if( iprec .eq. 1 ) epsfmm=.5d-3
        if( iprec .eq. 2 ) epsfmm=.5d-6
        if( iprec .eq. 3 ) epsfmm=.5d-9
        if( iprec .eq. 4 ) epsfmm=.5d-12
        if( iprec .eq. 5 ) epsfmm=.5d-15
        if( iprec .eq. 6 ) epsfmm=0
c       
        if(ifprint .eq. 1) call prin2('epsfmm=*',epsfmm,1)
c
c     set criterion for box subdivision (number of sources per box)
c
        if( iprec .eq. -2 ) nbox=3
        if( iprec .eq. -1 ) nbox=5
        if( iprec .eq. 0 ) nbox=8
        if( iprec .eq. 1 ) nbox=10
        if( iprec .eq. 2 ) nbox=15
        if( iprec .eq. 3 ) nbox=20
        if( iprec .eq. 4 ) nbox=25
        if( iprec .eq. 5 ) nbox=45
        if( iprec .eq. 6 ) nbox=nsource+ntarget
c
c
        if(ifprint .eq. 1) call prinf('nbox=*',nbox,1)
c
c
c     create quad-tree data structure
c
        t1=second()
C$        t1=omp_get_wtime()
        ntot = 100*(nsource+ntarget)+10000
        do ii = 1,10
           allocate (wlists(ntot))
           call lfmm2dparttree(ier,iprec,
     $        nsource,source,ntarget,target,
     $        nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $        nboxes,laddr,nlev,center,size,
     $        wlists,ntot,lused7)
           if (ier.ne.0) then
              deallocate(wlists)
              ntot = ntot*1.5
              call prinf(' increasing allocation, ntot is *',ntot,1)
           else
             goto 1200
           endif
        enddo
1200    continue
        if (ier.ne.0) then
           call prinf(' exceeded max allocation, ntot is *',ntot,1)
           ier = 4
           return
        endif
        t2=second()
C$        t2=omp_get_wtime()
c
c     lused7 is counter that steps through workspace,
c     keeping track of total memory used.
c
        lused7=1
c
c
c       ... prepare data structures 
c
        do i = 0,nlev
        scale(i) = 1.0d0
        boxsize = abs((size/2.0**i))
ccc        if (boxsize .lt. 1) scale(i) = boxsize
        scale(i) = boxsize
        enddo
c       
c
c       carve up workspace further
c
c     isourcesort is pointer for sorted source coordinates
c     itargetsort is pointer for sorted target locations
c     ichargesort is pointer for sorted charge densities
c     idipvecsort is pointer for sorted dipole orientation vectors
c     idipstrsort is pointer for sorted dipole densities
c
        isourcesort = lused7 + 5
        lsourcesort = 2*nsource
        itargetsort = isourcesort+lsourcesort
        ltargetsort = 2*ntarget
        ichargesort = itargetsort+ltargetsort
        lchargesort = 2*nsource
        idipvecsort = ichargesort+lchargesort
        if (ifdipole.eq.1) then
          ldipvec = 2*nsource
          ldipstr = 2*nsource
        else
          ldipvec = 2
          ldipstr = 2
        endif
        idipstrsort = idipvecsort + ldipvec
        lused7 = idipstrsort + ldipstr
c
c       ... allocate the potential and gradient arrays
c
        ipot = lused7
        lpot = 2*nsource
        lused7=lused7+lpot
c       
        igrad = lused7
        if( ifgrad .eq. 1) then
        lgrad = 2*nsource
        else
        lgrad=4
        endif
        lused7=lused7+lgrad
c      
        ihess = lused7
        if( ifhess .eq. 1) then
        lhess = 2*nsource
        else
        lhess=6
        endif
        lused7=lused7+lhess
c      
        ipottarg = lused7
        lpottarg = 2*ntarget
        lused7=lused7+lpottarg
c       
        igradtarg = lused7
        if( ifgradtarg .eq. 1) then
        lgradtarg = 2*ntarget
        else
        lgradtarg=4
        endif
        lused7=lused7+lgradtarg
c      
        ihesstarg = lused7
        if( ifhesstarg .eq. 1) then
        lhesstarg = 2*ntarget
        else
        lhesstarg=6
        endif
        lused7=lused7+lhesstarg
c      
        if(ifprint .eq. 1) call prinf(' lused7 is *',lused7,1)
c
c       based on FMM tolerance, compute expansion lengths nterms(i)
c      
        nmax = 0

        do i = 0,nlev
           bsize(i)=size/2.0d0**i
           call l2dterms(epsfmm, nterms(i), ier)
           if (nterms(i).gt. nmax .and. i.ge. 2) nmax = nterms(i)
        enddo
c
        if (ifprint.eq.1) 
     $     call prin2('in lfmm2dpart, bsize(0)=*',
     $     abs(bsize(0)),1)
c
        if (ifprint.eq.1) call prin2('bsize=*',bsize,nlev+1)
        if (ifprint.eq.1) call prinf('nterms=*',nterms,nlev+1)
c
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c
        iiaddr = lused7 
        imptemp = iiaddr + 2*nboxes
        lmptemp = (2*nmax+1)*2
        lused7 = imptemp + lmptemp
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate bulk FMM workspace,
     1                  lused7 is *',lused7,1)
           ier = 8
           return
        endif
c
c     reorder sources, targets so that each box holds
c     contiguous list of source/target numbers.

c
        call c2dreorder(nsource,source,ifcharge,charge,wlists(iisource),
     $     ifdipole,dipstr,
     1     w(isourcesort),w(ichargesort),w(idipstrsort)) 
c       
        call l2dreordertarg(ntarget,target,wlists(iitarget),
     $     w(itargetsort))
c
        if(ifprint .eq. 1) then
        call prinf('finished reordering=*',ier,1)
        call prinf('ier=*',ier,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('nlev=*',nlev,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('lused7=*',lused7,1)
        endif
c
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
        call l2dmpalloc(wlists(iwlists),w(iiaddr),nboxes,lmptot,nterms)
c
        if(ifprint .eq. 1) call prinf(' lmptot is *',lmptot,1)
c       
        irmlexp = 1
        lused7 = irmlexp + lmptot 
        if(ifprint .eq. 1) call prinf(' lused7 is *',lused7,1)
        allocate(wrmlexp(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate mpole expansion workspace,
     1                  lused7 is *',lused7,1)
           ier = 16
           return
        endif
c
c       
ccc        do i=lused7+1,lused7+1+100
ccc        w(i)=777
ccc        enddo
c
c     Memory allocation is complete. 
c     Call main fmm routine. There are, unfortunately, a lot
c     of parameters here. ifevalfar and ifevalloc determine
c     whether far gradient and local gradients (respectively) are to 
c     be evaluated. Setting both to 1 means that both will be
c     computed (which is the normal scenario).
c
        ifevalfar=1
        ifevalloc=1
c
        t1=second()
C$        t1=omp_get_wtime()
        call cfmm2dparttargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,w(isourcesort),wlists(iisource),
     $     ifcharge,w(ichargesort),
     $     ifdipole,w(idipstrsort),
     $     ifpot,w(ipot),ifgrad,w(igrad),ifhess,w(ihess),
     $     ntarget,w(itargetsort),wlists(iitarget),
     $     ifpottarg,w(ipottarg),ifgradtarg,w(igradtarg),
     $     ifhesstarg,w(ihesstarg),
     $     epsfmm,w(iiaddr),wrmlexp(irmlexp),w(imptemp),lmptemp,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists(iwlists),lwlists)
        t2=second()
C$        t2=omp_get_wtime()
        if( ifprint .eq. 1 ) call prin2('time in fmm main=*',t2-t1,1)
c
c       parameter ier from targmain routine is currently meaningless, reset to 0
        if( ier .ne. 0 ) ier = 0
c
        if(ifprint .eq. 1) call prinf('lwlists=*',lwlists,1)
        if(ifprint .eq. 1) then
        call prinf('lused total=*',lused7,1)
        call prinf('lused total(k)=*',lused7/1000,1)
        call prinf('lused total(M)=*',lused7/1000000,1)
        endif
c       
        if(ifprint .eq. 1) 
     $     call prin2('memory / point = *',(lused7)/dble(nsource),1)
c       
ccc        call prin2('after w=*', w(1+lused7-100), 2*100)
c
        if(ifpot .eq. 1) 
     $     call l2dpsort(nsource,wlists(iisource),w(ipot),pot)
        if(ifgrad .eq. 1) 
     $     call l2dpsort(nsource,wlists(iisource),w(igrad),grad)
        if(ifhess .eq. 1) 
     $     call l2dpsort(nsource,wlists(iisource),w(ihess),hess)
c
        if(ifpottarg .eq. 1 )
     $     call l2dpsort(ntarget,wlists(iitarget),w(ipottarg),pottarg)
        if(ifgradtarg .eq. 1) 
     $     call l2dpsort(ntarget,wlists(iitarget),w(igradtarg),gradtarg)
        if(ifhesstarg .eq. 1) 
     $     call l2dpsort(ntarget,wlists(iitarget),w(ihesstarg),hesstarg)
c       
        return
        end
c
c
c
c
c
        subroutine cfmm2dparttargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,sourcesort,isource,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,ntarget,
     $     targetsort,itarget,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg,
     $     epsfmm,iaddr,rmlexp,mptemp,lmptemp,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists,lwlists)
        implicit real *8 (a-h,o-z)
        real *8 sourcesort(2,*)
        integer isource(*)
        complex *16 chargesort(*)
        complex *16 dipstrsort(*)
        complex *16 ima
        complex *16 pot(*)
        complex *16 grad(*)
        complex *16 hess(*)
        real *8 targetsort(2,*)
        integer itarget(*)
        complex *16 pottarg(*)
        complex *16 gradtarg(*)
        complex *16 hesstarg(*)
        real *8 wlists(*)
        integer iaddr(2,nboxes)
        real *8 rmlexp(*)
        complex *16 mptemp(lmptemp)
        real *8 timeinfo(10)
        real *8 center(3)
        integer laddr(2,200)
        real *8 scale(0:200)
        real *8 bsize(0:200)
        integer nterms(0:200)
        integer list(10 000)
        complex *16 ptemp,gtemp,htemp
        integer box(15)
        real *8 center0(2),corners0(2,4)
        integer box1(15)
        real *8 center1(2),corners1(2,4)
        integer nterms_eval(4,0:200)
c
        data ima/(0.0d0,1.0d0)/
ccc        save
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
        ifprint=1
c
c
c       ... set the potential and gradient to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( ifgrad .eq. 1) then
           grad(i)=0
        endif
        if( ifhess .eq. 1) then
           hess(i)=0
        endif
        enddo
C$OMP END PARALLEL DO
c       
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( ifgradtarg .eq. 1) then
           gradtarg(i)=0
        endif
        if( ifhesstarg .eq. 1) then
           hesstarg(i)=0
        endif
        enddo
C$OMP END PARALLEL DO
c
        do i=1,10
        timeinfo(i)=0
        enddo
c
c
        if( ifevalfar .eq. 0 ) goto 8000
c       
c
        if (ifprint .ge. 2) 
     $     call prinf('nterms_eval=*',nterms_eval,4*(nlev+1))
c
c       ... set all multipole and local expansions to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP$PRIVATE(ibox,box,center0,corners0,level,ier)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do ibox = 1,nboxes
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
        if( level .ge. 2 ) then
        call l2dzero(rmlexp(iaddr(1,ibox)),nterms(level))
        call l2dzero(rmlexp(iaddr(2,ibox)),nterms(level))
        endif
        enddo
C$OMP END PARALLEL DO
c
c
        if(ifprint .ge. 1) 
     $     call prinf('=== STEP 1 (form mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,gtemp,htemp,cd) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 1200 ibox=1,nboxes
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        level=box(1)
        if( level .lt. 2 ) goto 1200
c
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,15)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
c        ipts=box(9)
c        npts=box(10)
c        call prinf('ipts=*',ipts,1)
c        call prinf('npts=*',npts,1)
        npts=box(10)
        if (ifprint .ge. 2) then
           call prinf('npts=*',npts,1)
           call prinf('isource=*',isource(box(9)),box(10))
        endif
        endif
c
c       ... prune all sourceless boxes
c
        if( box(10) .eq. 0 ) goto 1200
c
        if (nkids .eq. 0) then
c
c       ... form multipole expansions
c
	    radius = (corners0(1,1) - center0(1))**2
	    radius = radius + (corners0(2,1) - center0(2))**2
	    radius = sqrt(radius)
c
            call l2dzero(rmlexp(iaddr(1,ibox)),nterms(level))
            if_use_trunc = 0

            if( ifcharge .eq. 1 ) then
            call l2dformmp(ier,scale(level),sourcesort(1,box(9)),
     1  	chargesort(box(9)),npts,center0,nterms(level),
     2          rmlexp(iaddr(1,ibox)))        
            endif
c 
cc               call prin2('after formmp, rmlexp=*',
cc     $            rmlexp(iaddr(1,ibox)),2*(2*nterms(level)+1))
c
            if (ifdipole .eq. 1 ) then
               call l2dzero(mptemp,nterms(level))
               call l2dformmp_dp(ier,scale(level),
     $           sourcesort(1,box(9)),
     1           dipstrsort(box(9)),
     $           npts,center0,nterms(level),
     2           mptemp)
              call l2dadd(mptemp,rmlexp(iaddr(1,ibox)),nterms(level))
            endif
         endif
c
 1200    continue
C$OMP END PARALLEL DO
 1300    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(1)=t2-t1
c       
         if(ifprint .ge. 1)
     $      call prinf('=== STEP 2 (form lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 2, adaptive part, form local expansions, 
c           or evaluate the potentials and gradients directly
c 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,itype,list,nlist)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect3,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,gtemp,htemp,cd,ilist,npts) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
         do 3251 ibox=1,nboxes
c
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
c
ccc         if( box(10) .eq. 0 ) goto 3251
c
         itype=4
         call d2tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list3=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
ccc         if( box(10) .eq. 0 ) nlist=0
c
c
c       ... note that lists 3 and 4 are dual
c
c       ... form local expansions for all boxes in list 3
c       ... if target is childless, evaluate directly (if cheaper)
c        
         do 3250 ilist=1,nlist
            jbox=list(ilist)
            call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
c        
            npts=box1(10)            
            if( npts .eq. 0 ) goto 3250
c
            level0=box(1)
            level1=box1(1)
c
c            ifdirect3 = 0
c
c            if( box1(10) .lt. (nterms(level1)+1)/2 .and.
c     $          box(10) .lt. (nterms(level1)+1)/2 ) ifdirect3 = 1
c
            ifdirect3 = 0
c
            if( ifdirect3 .eq. 0 ) then
               npts=box1(10)
               if_use_trunc = 0
c
               if( ifcharge .eq. 1 ) then
               call l2dformta_add(ier,scale(level0),
     $            sourcesort(1,box1(9)),
     1            chargesort(box1(9)),npts,center0,nterms(level0),
     2            rmlexp(iaddr(2,ibox)))
               endif
c
               if( ifdipole .eq. 1 ) then
               call l2dformta_dp_add(ier,scale(level0),
     1            sourcesort(1,box1(9)),dipstrsort(box1(9)),
     2            npts,center0,
     3            nterms(level0),rmlexp(iaddr(2,ibox)))
               endif
c
            else

            call cfmm2dpart_direct(box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,
     $         ifpot,pot,ifgrad,grad,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $         ifhesstarg,hesstarg)
c
            endif
 3250    continue
c
 3251    continue
C$OMP END PARALLEL DO
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(2)=t2-t1
c
c
        if(ifprint .ge. 1)
     $      call prinf('=== STEPS 3,4,5 ====*',i,0)
        ifprune_list2 = 1
        if (ifpot.eq.1) ifprune_list2 = 0
        if (ifgrad.eq.1) ifprune_list2 = 0
        if (ifhess.eq.1) ifprune_list2 = 0
        call lfmm2d_list2
     $     (bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,
     $     ifprune_list2)
c
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (eval mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 6, adaptive part, evaluate multipole expansions, 
c           or evaluate the potentials and gradients directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,itype,list,nlist)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect4,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,gtemp,htemp,cd,ilist,level) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
         do 3252 ibox=1,nboxes
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
c
c       ... prune all sourceless boxes
c
ccc         if( box(10) .eq. 0 ) goto 3252
c
         itype=3
         call d2tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list4=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
ccc         if( box(10) .eq. 0 ) nlist=0
c
c       ... note that lists 3 and 4 are dual
c
c       ... evaluate multipole expansions for all boxes in list 4 
c       ... if source is childless, evaluate directly (if cheaper)
c
         do ilist=1,nlist
            jbox=list(ilist)
            call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
c
            level=box1(1)
c
c            ifdirect4 = 0
c
c            if (box1(10) .lt. (nterms(level)+1)/2 .and.
c     $         box(10) .lt. (nterms(level)+1)/2 ) ifdirect4 = 1
c
            ifdirect4 = 0
c
            if (ifdirect4 .eq. 0) then

            if( box(10) .gt. 0 ) 
     $         call c2dmpevalall(scale(level),center1,
     $         rmlexp(iaddr(1,jbox)),nterms(level),
     $         sourcesort(1,box(9)),box(10),
     $         ifpot,pot(box(9)),
     $         ifgrad,grad(box(9)),
     $         ifhess,hess(box(9)))

            if( box(12) .gt. 0 ) 
     $         call c2dmpevalall(scale(level),center1,
     $         rmlexp(iaddr(1,jbox)),nterms(level),
     $         targetsort(1,box(11)),box(12),
     $         ifpottarg,pottarg(box(11)),
     $         ifgradtarg,gradtarg(box(11)),
     $         ifhesstarg,hesstarg(box(11)))

            else
            
            call cfmm2dpart_direct(box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,
     $         ifpot,pot,ifgrad,grad,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $         ifhesstarg,hesstarg)

            endif
        enddo
c
 3252   continue
C$OMP END PARALLEL DO
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(6)=t2-t1
c

        if(ifprint .ge. 1)
     $     call prinf('=== STEP 7 (eval lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 7, evaluate local expansions
c       and all gradients directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,ier)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6201 ibox=1,nboxes
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,15)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(10)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(9)),box(10))
            endif
        endif
c
        if (nkids .eq. 0) then
c
c       ... evaluate local expansions
c       
        level=box(1)
        npts=box(10)
c       
        if (level .ge. 2) then

        if( box(10) .gt. 0 )
     $     call c2dtaevalall(scale(level),center0,
     $     rmlexp(iaddr(2,ibox)),nterms(level),
     $     sourcesort(1,box(9)),box(10),
     $     ifpot,pot(box(9)),ifgrad,grad(box(9)),
     $     ifhess,hess(box(9)))

        if( box(12) .gt. 0 )
     $     call c2dtaevalall(scale(level),center0,
     $     rmlexp(iaddr(2,ibox)),nterms(level),
     $     targetsort(1,box(11)),box(12),
     $     ifpottarg,pottarg(box(11)),ifgradtarg,gradtarg(box(11)),
     $     ifhesstarg,hesstarg(box(11)))

        endif
c
        endif
c
 6201   continue
C$OMP END PARALLEL DO
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(7)=t2-t1
c
c
 8000   continue
c
c
        if( ifevalloc .eq. 0 ) goto 9000
c 
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 8 (direct) =====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 8, evaluate direct interactions 
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,nkids,list,nlist,npts)
C$OMP$PRIVATE(jbox,box1,center1,corners1)
C$OMP$PRIVATE(ier,ilist,itype) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6202 ibox=1,nboxes
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,15)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(10)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(9)),box(10))
            endif
        endif
c
c
        if (nkids .eq. 0) then
c
c       ... evaluate self interactions
c
        call cfmm2dpart_direct_self_sym(box,sourcesort,
     $     ifcharge,chargesort,ifdipole,dipstrsort,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     targetsort,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d2tgetl(ier,ibox,itype,list,nlist,wlists)
        if (ifprint .ge. 2) call prinf('list1=*',list,nlist)
c
c       ... for all pairs in list #1, 
c       evaluate the potentials and gradients directly
c
            do 6203 ilist=1,nlist
               jbox=list(ilist)
               call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
         if( box1(10) .eq. 0 ) goto 6203
c    
            call cfmm2dpart_direct(box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,
     $         ifpot,pot,ifgrad,grad,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $         ifhesstarg,hesstarg)
c
 6203       continue
        endif
c
 6202   continue
C$OMP END PARALLEL DO
c
ccc        call prin2('inside fmm, pot=*',pot,2*nsource)
ccc        call prin2('inside fmm, grad=*',grad,2*nsource)
ccc        call prin2('inside fmm, hess=*',hess,2*nsource)
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
 9000   continue
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
        if (ifprint .ge. 1) call prin2('timeinfo=*',timeinfo,8)
c       
        d=0
        do i=1,8
        d=d+timeinfo(i)
        enddo
c       
        if (ifprint .ge. 1) call prin2('sum(timeinfo)=*',d,1)
c
        if (ifprint .ge. 1) then
        call prinf('nboxes=*',nboxes,1)
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        endif
c       
        return
        end
c
c
c
c
c
        subroutine cfmm2dpart_direct_self(box,
     $     source,ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        integer box(15),box1(15)
c
        real *8 source(2,*)
        complex *16 charge(*),dipstr(*)
        real *8 target(2,*)
c
        complex *16 pot(*),grad(*),hess(*)
        complex *16 pottarg(*),gradtarg(*),hesstarg(*)
        complex *16 ptemp,gtemp,htemp
c
c
        do 6160 j=box(9),box(9)+box(10)-1
        do 6150 i=box(9),box(9)+box(10)-1
            if (i .eq. j) goto 6150
            call cpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge(i),ifdipole,dipstr(i),
     1         source(1,j),ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
            if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
            if (ifgrad .eq. 1) then
               grad(j)=grad(j)+gtemp
            endif
            if (ifhess .eq. 1) then
               hess(j)=hess(j)+htemp
            endif
 6150   continue
 6160   continue
        do j=box(11),box(11)+box(12)-1
        do i=box( 9),box( 9)+box(10)-1
            call cpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge(i),ifdipole,dipstr(i),
     1         target(1,j),ifpottarg,ptemp,ifgradtarg,gtemp,
     $         ifhesstarg,htemp)
            if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
            if (ifgradtarg .eq. 1) then
               gradtarg(j)=gradtarg(j)+gtemp
            endif
            if (ifhesstarg .eq. 1) then
               hesstarg(j)=hesstarg(j)+htemp
            endif
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine cfmm2dpart_direct_self_sym(box,
     $     source,ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        integer box(15),box1(15)
c
        real *8 source(2,*)
        complex *16 charge(*),dipstr(*)
        real *8 target(2,*)
c
        complex *16 pot(*),grad(*),hess(*)
        complex *16 pottarg(*),gradtarg(*),hesstarg(*)
        complex *16 ptemp,gtemp,htemp
        complex *16 ptemp1,gtemp1,htemp1
        complex *16 ptemp2,gtemp2,htemp2
c
c
        do 6160 j=box(9),box(9)+box(10)-1
        do 6150 i=box(9),box(9)+box(10)-1
            if (i .ge. j) goto 6150
c            call cpotgrad2d_sdp_sym(source(1,i),source(1,j),
c     $         ifcharge,charge(i),charge(j),
c     $         ifdipole,dipstr(i),dipstr(j),
c     1         ifpot,ptemp1,ptemp2,
c     $         ifgrad,gtemp1,gtemp2,ifhess,htemp1,htemp2)
c            if (ifpot .eq. 1) then 
c               pot(i)=pot(i)+ptemp1
c               pot(j)=pot(j)+ptemp2
c            endif
c            if (ifgrad .eq. 1) then
c               grad(i)=grad(i)+gtemp1
c               grad(j)=grad(j)+gtemp2
c            endif
c            if (ifhess .eq. 1) then
c               hess(i)=hess(i)+htemp1
c               hess(j)=hess(j)+htemp2
c            endif
            call cpotgrad2d_sdp_sym_add(source(1,i),source(1,j),
     $         ifcharge,charge(i),charge(j),
     $         ifdipole,dipstr(i),dipstr(j),
     1         ifpot,pot(i),pot(j),
     $         ifgrad,grad(i),grad(j),ifhess,hess(i),hess(j))
 6150   continue
 6160   continue
        do j=box(11),box(11)+box(12)-1
        do i=box( 9),box( 9)+box(10)-1
            call cpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge(i),ifdipole,dipstr(i),
     1         target(1,j),ifpottarg,ptemp,ifgradtarg,gtemp,
     $         ifhesstarg,htemp)
            if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
            if (ifgradtarg .eq. 1) then
               gradtarg(j)=gradtarg(j)+gtemp
            endif
            if (ifhesstarg .eq. 1) then
               hesstarg(j)=hesstarg(j)+htemp
            endif
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine cfmm2dpart_direct(box,box1,
     $     source,ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        integer box(15),box1(15)
c
        real *8 source(2,*)
        complex *16 charge(*),dipstr(*)
        real *8 target(2,*)
c
        complex *16 pot(*),grad(*),hess(*)
        complex *16 pottarg(*),gradtarg(*),hesstarg(*)
        complex *16 ptemp,gtemp,htemp
c
c
        do j=box1(9),box1(9)+box1(10)-1
        call cpotgrad2dall_sdp
     $     (source(1,box(9)),box(10),ifcharge,charge(box(9)),
     $     ifdipole,dipstr(box(9)),source(1,j),
     1     ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
        if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
        if (ifgrad .eq. 1) then
          grad(j)=grad(j)+gtemp
        endif
        if (ifhess .eq. 1) then
          hess(j)=hess(j)+htemp
        endif
        enddo
c
        do j=box1(11),box1(11)+box1(12)-1
        call cpotgrad2dall_sdp
     $     (source(1,box(9)),box(10),ifcharge,charge(box(9)),
     $     ifdipole,dipstr(box(9)),target(1,j),
     $     ifpottarg,ptemp,ifgradtarg,gtemp,ifhesstarg,htemp)
        if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
        if (ifgradtarg .eq. 1) then
          gradtarg(j)=gradtarg(j)+gtemp
        endif
        if (ifhesstarg .eq. 1) then
          hesstarg(j)=hesstarg(j)+htemp
        endif
        enddo
c       
        return
        end
c
c
c
c
c
        subroutine c2dpartdirect(nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
c       Generalized Cauchy interactions in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets via direct O(N^2) algorithm.
c
c       We use log(z) for the Green's function.
c       Self-interactions are not included.
c   
c c2d: charge and dipstr are complex valued, z are complex numbers.
c
c        Note, that the complex valued logarithm is a multi-valued
c        function, so the potential values have to be interpreted
c        carefully, if charges are specified.  For example, only the
c        real part of potential is meaningful for real valued charges.
c        The gradients and hessians are valid for arbitrary complex charges.
c
c \phi(z_i) = \sum_{j\ne i} charge_j \log(z_i-z_j) + dipstr_j \frac{1}{z_i-z_j}
c
c        In this routine, we define the gradient as the first
c        derivative with respect to z, and the hessian as the second
c        derivative with respect to z.
c
c \gradient \phi(z_i) = \frac{\partial \phi(z_i)}{\partial z}
c \hessian  \phi(z_i) = \frac{\partial^2 \phi(z_i)}{\partial z^2}
c
c       INPUT PARAMETERS:
c
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: complex *16 (nsource): dipole strengths
c
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       ifgrad:  gradient flag (1=compute gradient, otherwise no)
c       ifhess:  hessian flag (1=compute hessian, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (2,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       ifgradtarg:  target gradient flag 
c                   (1=compute gradient, otherwise no)
c       ihesstarg:  target hessian flag 
c                   (1=compute hessian, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       pot: complex *16 (nsource): potential at source locations
c       grad: complex *16 (nsource): gradient  at source locations
c       hess: complex *16 (nsource): hessian at source locations
c       pottarg: complex *16 (ntarget): potential at target locations 
c       gradtarg: complex *16 (ntarget): gradient  at target locations 
c       hesstarg: complex *16 (ntarget): hessian at target locations
c
c
c
        real *8 source(2,*)
        complex *16 charge(*),dipstr(*)
        real *8 target(2,*)
c
        complex *16 pot(*),grad(*),hess(*)
        complex *16 pottarg(*),gradtarg(*),hesstarg(*)
        complex *16 ptemp,gtemp,htemp
c
c
        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( ifgrad .eq. 1) then
           grad(i)=0
           grad(i)=0
        endif
        if( ifhess .eq. 1) then
           hess(i)=0
           hess(i)=0
           hess(i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( ifgradtarg .eq. 1) then
           gradtarg(i)=0
           gradtarg(i)=0
        endif
        if( ifhesstarg .eq. 1) then
           hesstarg(i)=0
           hesstarg(i)=0
           hesstarg(i)=0
        endif
        enddo
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 .or. ifhess .eq. 1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,gtemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 6160 j=1,nsource
        do 6150 i=1,nsource
            if (i .eq. j) goto 6150
            call cpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge(i),ifdipole,dipstr(i),
     1         source(1,j),ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
            if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
            if (ifgrad .eq. 1) then
               grad(j)=grad(j)+gtemp
            endif
            if (ifhess .eq. 1) then
               hess(j)=hess(j)+htemp
            endif
 6150   continue
 6160   continue
C$OMP END PARALLEL DO
        endif
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 
     $      .or. ifhesstarg .eq. 1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,gtemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do j=1,ntarget
        do i=1,nsource
            call cpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge(i),ifdipole,dipstr(i),
     1         target(1,j),ifpottarg,ptemp,ifgradtarg,gtemp,
     $         ifhesstarg,htemp)
            if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
            if (ifgradtarg .eq. 1) then
               gradtarg(j)=gradtarg(j)+gtemp
            endif
            if (ifhesstarg .eq. 1) then
               hesstarg(j)=hesstarg(j)+htemp
            endif
        enddo
        enddo
C$OMP END PARALLEL DO
        endif
c
        return
        end
c
c
c
c
c
        subroutine c2dreorder(nsource,source,
     $     ifcharge,charge,isource,ifdipole,
     1     dipstr,sourcesort,chargesort,dipstrsort) 
        implicit real *8 (a-h,o-z)
        real *8 source(2,*),sourcesort(2,*)
        integer isource(*)
        complex *16 charge(*),chargesort(*),dipstr(*),dipstrsort(*)
c       
ccc        call prinf('nsource=*',nsource,1)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i = 1,nsource
        sourcesort(1,i) = source(1,isource(i))
        sourcesort(2,i) = source(2,isource(i))
        if( ifcharge .ge. 1 ) then
        chargesort(i) = charge(isource(i))
        endif
        if (ifdipole .ge. 1) then
        dipstrsort(i) = dipstr(isource(i))
        endif
        enddo
C$OMP END PARALLEL DO
        return
        end
c
c
c
c
c
cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c       
c     This file contains the main FMM routines and some related
c     subroutines for evaluating real-valued Laplace potentials and
c     gradients due to point charges and dipoles.  (FORTRAN 90 VERSION)
c
c lfmm2d: charge and dipstr are complex valued, x \in R^2
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                    + dipstr_j (dipvec_j \dot \grad_j log|x_i-x_j|)
c
c or 
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                   + dipstr_j (dipvec_j \dot (x_i-x_j)) * (-1/|x_i-x_j|^2)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     rfmm2dpart - Laplace FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     rfmm2dpartself - Laplace FMM in R^2: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     rfmm2dparttarg - Laplace FMM in R^2: evaluate all pairwise
c         particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c     r2dpartdirect - Laplace interactions in R^2:  evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c         interactions with targets via direct O(N^2) algorithm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the Laplace particle FMM in R^2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine rfmm2dpart(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c              
c       Laplace FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c
c       We use log(r) for the Green's function.
c       Self-interactions are not included.
c   
c rfmm2d: charge and dipstr are real valued, x \in R^2
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                    + dipstr_j (dipvec_j \dot \grad_j log|x_i-x_j|)
c
c or 
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                   + dipstr_j (dipvec_j \dot (x_i-x_j)) * (-1/|x_i-x_j|^2)
c
c       The main FMM routine permits both evaluation at sources
c       and at a collection of targets. 
c       This subroutine is used to simplify the user interface 
c       (by setting the number of targets to zero) and calling the more 
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c  
        real *8 source(2,*)
        real *8 charge(*)
        real *8 dipstr(*)
        real *8 dipvec(2,*)
        real *8 pot(*)
        real *8 grad(2,*)
        real *8 hess(3,*)
c
        real *8 target(2,1)
        real *8 pottarg(1)
        real *8 gradtarg(2,1)
        real *8 hesstarg(3,1)
c       
        ntarget=0
        ifpottarg=0
        ifgradtarg=0
        ifhesstarg=0
c
        call rfmm2dparttarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
c
c
c
        subroutine rfmm2dpartself(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess)
        implicit real *8 (a-h,o-z)
c              
c       Laplace FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c
c       We use log(r) for the Green's function.
c       Self-interactions are not included.
c   
c rfmm2d: charge and dipstr are real valued, x \in R^2
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                    + dipstr_j (dipvec_j \dot \grad_j log|x_i-x_j|)
c
c or 
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                   + dipstr_j (dipvec_j \dot (x_i-x_j)) * (-1/|x_i-x_j|^2)
c
c       The main FMM routine permits both evaluation at sources
c       and at a collection of targets. 
c       This subroutine is used to simplify the user interface 
c       (by setting the number of targets to zero) and calling the more 
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c  
        real *8 source(2,*)
        real *8 charge(*)
        real *8 dipstr(*)
        real *8 dipvec(2,*)
        real *8 pot(*)
        real *8 grad(2,*)
        real *8 hess(3,*)
c
        real *8 target(2,1)
        real *8 pottarg(1)
        real *8 gradtarg(2,1)
        real *8 hesstarg(3,1)
c       
        ntarget=0
        ifpottarg=0
        ifgradtarg=0
        ifhesstarg=0
c
        call rfmm2dparttarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
c
c
c
        subroutine rfmm2dparttarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c       
c       Laplace FMM in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets.
c
c       We use log(r) for the Green's function.
c       Self-interactions are not included.
c   
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine cfmm2dparttarg
c
c rfmm2d: charge and dipstr are real valued, x \in R^2
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                    + dipstr_j (dipvec_j \dot \grad_j log|x_i-x_j|)
c
c or 
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                   + dipstr_j (dipvec_j \dot (x_i-x_j)) * (-1/|x_i-x_j|^2)
c
c       INPUT PARAMETERS:
c
c       iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c                  5 => tolerance =.5d-15
c
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: real *8 (nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: real *8 (nsource): dipole strengths
c       dipvec: real *8 (2,nsource): dipole orientation vectors. 
c
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       ifgrad:  gradient flag (1=compute gradient, otherwise no)
c       ifhess:  hessian flag (1=compute hessian, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (2,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       ifgradtarg:  target gradient flag 
c                   (1=compute gradient, otherwise no)
c       ihesstarg:  target hessian flag 
c                   (1=compute hessian, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       ier   =  error return code
c                ier=0     =>  normal execution
c                ier=4     =>  cannot allocate tree workspace
c                ier=8     =>  cannot allocate bulk FMM  workspace
c                ier=16    =>  cannot allocate mpole expansion
c                              workspace in FMM
c
c       pot: real *8 (nsource): potential at source locations
c       grad: real *8 (2,nsource): gradient  at source locations
c       hess: real *8 (3,nsource): hessian at source locations
c       pottarg: real *8 (ntarget): potential at target locations 
c       gradtarg: real *8 (2,ntarget): gradient  at target locations 
c       hesstarg: real *8 (3,ntarget): hessian at target locations
c
c
c
        real *8 source(2,*)
        real *8 charge(*)
        real *8 dipstr(*)
        real *8 dipvec(2,*)
        real *8 pot(*)
        real *8 grad(2,*)
        real *8 hess(3,*)
c
        real *8 target(2,*)
        real *8 pottarg(*)
        real *8 gradtarg(2,*)
        real *8 hesstarg(3,*)
c       
        real *8 timeinfo(10)
c
        complex *16 ima
c
        complex *16, allocatable :: charge1(:)
        complex *16, allocatable :: dipstr1(:)
        complex *16, allocatable :: pot1(:)
        complex *16, allocatable :: grad1(:)
        complex *16, allocatable :: hess1(:)
        complex *16, allocatable :: pottarg1(:)
        complex *16, allocatable :: gradtarg1(:)
        complex *16, allocatable :: hesstarg1(:)
c
        data ima/(0.0d0,1.0d0)/
c
c
        if( ifcharge .eq. 1 ) then
        allocate(charge1(nsource))
        else
        allocate(charge1(1))
        endif
c
        if( ifdipole .eq. 1 ) then
        allocate(dipstr1(nsource))
        else
        allocate(dipstr1(1))
        endif
c
        if( ifpot .eq. 1 ) then
        allocate(pot1(nsource))
        else
        allocate(pot1(1))
        endif
c
        if( ifgrad .eq. 1 ) then
        allocate(grad1(nsource))
        else
        allocate(grad1(1))
        endif
c
        if( ifhess .eq. 1 ) then
        allocate(hess1(nsource))
        else
        allocate(hess1(1))
        endif
c
        if( ifpottarg .eq. 1 ) then
        allocate(pottarg1(ntarget))
        else
        allocate(pottarg1(1))
        endif
c
        if( ifgradtarg .eq. 1 ) then
        allocate(gradtarg1(ntarget))
        else
        allocate(gradtarg1(1))
        endif
c
        if( ifhesstarg .eq. 1 ) then
        allocate(hesstarg1(ntarget))
        else
        allocate(hesstarg1(1))
        endif
c
        if( ifcharge .eq. 1 ) then
        do i=1,nsource
        charge1(i)=charge(i)
        enddo
        endif
c
        if( ifdipole .eq. 1 ) then
        do i=1,nsource
        dipstr1(i)=-dipstr(i)*(dipvec(1,i)+ima*dipvec(2,i))
        enddo
        endif
c
        call cfmm2dparttarg(ier,iprec,
     $     nsource,source,
     $     ifcharge,charge1,ifdipole,dipstr1,
     $     ifpot,pot1,ifgrad,grad1,ifhess,hess1,
     $     ntarget,target,ifpottarg,pottarg1,ifgradtarg,gradtarg1,
     $     ifhesstarg,hesstarg1)
c
        do i=1,nsource
        if( ifpot .eq. 1 ) then
        pot(i)=dble(pot1(i))
        endif
        if( ifgrad .eq. 1 ) then
        grad(1,i)=dble(grad1(i))
        grad(2,i)=-imag(grad1(i))
        endif
        if( ifhess .eq. 1 ) then
        hess(1,i)=dble(hess1(i))
        hess(2,i)=-imag(hess1(i))
        hess(3,i)=-dble(hess1(i))
        endif
        enddo
c
        do i=1,ntarget
        if( ifpottarg .eq. 1 ) then
        pottarg(i)=dble(pottarg1(i))
        endif
        if( ifgradtarg .eq. 1 ) then
        gradtarg(1,i)=dble(gradtarg1(i))
        gradtarg(2,i)=-imag(gradtarg1(i))
        endif
        if( ifhesstarg .eq. 1 ) then
        hesstarg(1,i)=dble(hesstarg1(i))
        hesstarg(2,i)=-imag(hesstarg1(i))
        hesstarg(3,i)=-dble(hesstarg1(i))
        endif
        enddo
c
        return
        end
c
c
c
c
c
        subroutine r2dpartdirect(nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,ifgrad,grad,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
c       Laplace interactions in R^2: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets via direct O(N^2) algorithm.
c
c       We use log(r) for the Green's function.
c       Self-interactions are not included.
c   
c r2d: charge and dipstr are real valued, x \in R^2
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                    + dipstr_j (dipvec_j \dot \grad_j log|x_i-x_j|)
c
c or 
c
c \phi(x_i) = \sum_{j\ne i}   charge_j \log |x_i-x_j|  
c                   + dipstr_j (dipvec_j \dot (x_i-x_j)) * (-1/|x_i-x_j|^2)
c
c       INPUT PARAMETERS:
c
c       nsource: integer:  number of sources
c       source: real *8 (2,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: real *8 (nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: real *8 (nsource): dipole strengths
c       dipvec: real *8 (2,nsource): dipole orientation vectors. 
c
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       ifgrad:  gradient flag (1=compute gradient, otherwise no)
c       ifhess:  hessian flag (1=compute hessian, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (2,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       ifgradtarg:  target gradient flag 
c                   (1=compute gradient, otherwise no)
c       ihesstarg:  target hessian flag 
c                   (1=compute hessian, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       pot: real *8 (nsource): potential at source locations
c       grad: real *8 (2,nsource): gradient  at source locations
c       hess: real *8 (3,nsource): hessian at source locations
c       pottarg: real *8 (ntarget): potential at target locations 
c       gradtarg: real *8 (2,ntarget): gradient  at target locations 
c       hesstarg: real *8 (3,ntarget): hessian at target locations
c
c
c
        real *8 source(2,*)
        real *8 charge(*),dipstr(*),dipvec(2,*)
        real *8 target(2,*)
c
        real *8 pot(*),grad(2,*),hess(3,*)
        real *8 pottarg(*),gradtarg(2,*),hesstarg(3,*)
        real *8 ptemp,gtemp(2),htemp(3)
c
c
        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( ifgrad .eq. 1) then
           grad(1,i)=0
           grad(2,i)=0
        endif
        if( ifhess .eq. 1) then
           hess(1,i)=0
           hess(2,i)=0
           hess(3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( ifgradtarg .eq. 1) then
           gradtarg(1,i)=0
           gradtarg(2,i)=0
        endif
        if( ifhesstarg .eq. 1) then
           hesstarg(1,i)=0
           hesstarg(2,i)=0
           hesstarg(3,i)=0
        endif
        enddo
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 .or. ifhess .eq. 1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,gtemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=1,nsource
        do i=1,nsource
            if (i .eq. j) cycle
            call rpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge(i),ifdipole,dipstr(i),dipvec(1,i),
     1         source(1,j),ifpot,ptemp,ifgrad,gtemp,ifhess,htemp)
            if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
            if (ifgrad .eq. 1) then
               grad(1,j)=grad(1,j)+gtemp(1)
               grad(2,j)=grad(2,j)+gtemp(2)
            endif
            if (ifhess .eq. 1) then
               hess(1,j)=hess(1,j)+htemp(1)
               hess(2,j)=hess(2,j)+htemp(2)
               hess(3,j)=hess(3,j)+htemp(3)
            endif
        enddo
        enddo
C$OMP END PARALLEL DO
        endif
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 
     $      .or. ifhesstarg .eq. 1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,gtemp,htemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=1,ntarget
        do i=1,nsource
            call rpotgrad2d_sdp(source(1,i),
     $         ifcharge,charge(i),ifdipole,dipstr(i),dipvec(1,i),
     1         target(1,j),ifpottarg,ptemp,ifgradtarg,gtemp,
     $         ifhesstarg,htemp)
            if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
            if (ifgradtarg .eq. 1) then
               gradtarg(1,j)=gradtarg(1,j)+gtemp(1)
               gradtarg(2,j)=gradtarg(2,j)+gtemp(2)
            endif
            if (ifhesstarg .eq. 1) then
               hesstarg(1,j)=hesstarg(1,j)+htemp(1)
               hesstarg(2,j)=hesstarg(2,j)+htemp(2)
               hesstarg(3,j)=hesstarg(3,j)+htemp(3)
            endif
        enddo
        enddo
C$OMP END PARALLEL DO
        endif
c
        return
        end
c
c
c
c
c
cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c       
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the routines for Helmholtz FMM in R^2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine lfmm2dparttree(ier,iprec,
     $     nsource,source,ntarget,target,
     $     nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $     nboxes,laddr,nlev,center,size,
     $     w,lw,lused7)
        implicit real *8 (a-h,o-z)
c       
c       Helmholtz FMM in R^2: build the quad-tree
c
c     INPUT PARAMETERS:
c
c       nsource: integer:  number of sources
c       source: real *8 (2,n):  source locations
c       w: real *8 (lw): workspace
c       lw:  length of workspace
c
c     OUTPUT PARAMETERS:
c
c       ier   =  error return code
c       lused7 = the amount of workspace used
c
c
        real *8 source(2,*),target(2,*)
c       
        real *8 center(2)
c       
        integer laddr(2,200)
        integer box(15)
        real *8 center0(2),corners0(2,4)
c       
        integer box1(15)
        real *8 center1(2),corners1(2,4)
c
        real *8 w(*)
c       
        ier=0
c       
        done=1
        pi=4*atan(done)
c       
        lused7=0
        ifprint=0
c       
        iisource=1
        lused7=lused7+nsource
        if (ifprint.eq.1) call prinf('lused7=*',lused7,1)
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c
        iitarget=iisource+nsource
        lused7=lused7+ntarget
        if (ifprint.eq.1) call prinf('lused7=*',lused7,1)
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c
        iwlists=iisource+lused7+10
c
c       ... construct the adaptive FMM quad-tree structure
c
c        call d2tstrcr(ier,source,nsource,nbox,
c     $     nboxes,w(iisource),laddr,nlev,center,size,
c     $     target,ntarget,w(iitarget),w(iwlists),lw-lused7,lused)
        ifempty=0
        minlevel=0
        maxlevel=100 !TODO was before 30
        call d2tstrcrem(ier,source,nsource,nbox,
     $     nboxes,w(iisource),laddr,nlev,center,size,
     $     target,ntarget,w(iitarget),w(iwlists),lw-lused7,lused,
     $     ifempty,minlevel,maxlevel)
c
        if( ier .ne. 0 ) return
c
        lwlists=lused
        lused7=lused7+lwlists
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c       
c       ... optional, plot the oct-tree in gnuplot compatible format
c
        ifplot = 0
        if (ifplot .eq. 1 .and. nsource .lt. 10000 ) then
c
c       ... plot the boxes
c
        iw=51
        call plot_box2d(iw,center,size)
c       
        itag=1
        iw=52
        call plot_label2d(iw,center,size,itag,itag)
c
        iw=60
        call plot_points2d(iw,source,nsource)
c       
        iw=63
        call plot_points2d(iw,target,ntarget)
c       
        do ibox=1,nboxes
           call d2tgetb(ier,ibox,box,center0,corners0,w(iwlists))
           level=box(1)
           size0=size/2**level
c
           iw=61
           call plot_box2d(iw,center0,size0)
c       
           itag=ibox
           iw=62
           call plot_label2d(iw,center0,size0,itag,itag)
        enddo  
c      
        endif
c
        return
        end
c
c
c
c
c
        subroutine lfmm2d_list2
     $     (bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,
     $     ifprune_list2)
        implicit real *8 (a-h,o-z)
c
        integer iaddr(2,*),laddr(2,*),nterms(0:*)
        real *8 rmlexp(*),scale(0:*)
        integer itable(-3:3,-3:3)
c
        integer list(10 000)
c
        integer box(15)
        real *8 bsize(0:200)
        real *8 center0(2),corners0(2,4)
c       
        integer box1(15)
        real *8 center1(2),corners1(2,4)
c
        real *8 wlists(*)
        complex *16 mptemp(lmptemp)
        complex *16 ptemp,ftemp(2),htemp(3)
c
        real *8 timeinfo(10)
c
        real *8, allocatable :: carray(:,:)
c               
c
        ldc = 100
        allocate( carray(0:ldc,0:ldc) )
        call l2d_init_carray(carray,ldc)
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
        ifprint=1
c
         if (ifprint .ge. 1) 
     $     call prinf('=== STEP 3 (merge mp) ===*',i,0)
         t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 3, merge all multipole expansions
c       
ccc         do 2200 ibox=nboxes,1,-1
         do 2300 ilev=nlev,3,-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,ftemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
         do 2200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
         call d2tgetb(ier,ibox,box,center0,corners0,wlists)
         call d2tnkids(box,nkids)
c
c       ... prune all sourceless boxes
c
         if( box(10) .eq. 0 ) goto 2200
c
         if (nkids .ne. 0) then
c
         level0=box(1)
         if( level0 .ge. 2 ) then
ccc         if (level0 .ge. 0) then
            radius = (corners0(1,1) - center0(1))**2
            radius = radius + (corners0(2,1) - center0(2))**2
            radius = sqrt(radius)
c       
            if( ifprint .ge. 2 ) then
               call prin2('radius=*',radius,1)
               call prinf('ibox=*',ibox,1)
               call prinf('box=*',box,15)
               call prinf('nkids=*',nkids,1)
            endif
c
c       ... merge multipole expansions of the kids 
c
            call l2dzero(rmlexp(iaddr(1,ibox)),nterms(level0))
            if (ifprint .ge. 2) then
               call prin2('center0=*',center0,2)
            endif
c
            do 2100 i = 1,4
               jbox = box(4+i)
               if (jbox.eq.0) goto 2100
               call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
               if (ifprint .ge. 2) then
               call prinf('jbox=*',jbox,1)
               call prin2('center1=*',center1,2)
               endif
               level1=box1(1)
               if( nterms(level0)+nterms(level1) .gt. 95 ) then
               call l2dmpmp(scale(level1),center1,
     1            rmlexp(iaddr(1,jbox)),nterms(level1),scale(level0),
     1            center0,mptemp,nterms(level0))
               else
               call l2dmpmp_carray(scale(level1),center1,
     1            rmlexp(iaddr(1,jbox)),nterms(level1),scale(level0),
     1            center0,mptemp,nterms(level0),carray,ldc)
               endif
               call l2dadd(mptemp,rmlexp(iaddr(1,ibox)),
     1            nterms(level0))
 2100       continue
            if (ifprint .ge. 2) then
            call prinf('=============*',x,0)
            endif
c       ... mark the local expansion of all kids and the parent
c
            endif
         endif
 2200    continue
C$OMP END PARALLEL DO
 2300    continue
c
c
ccc        call prinf('=== UPWARD PASS COMPLETE ===*',i,0)
c
c------------------------------------------------------------
c      DEBUGGING SEGMENT - once all multipole expansions are merged
c      to top level, one can compare it to the direct formation of the
c      expansion at the top level from the source locations.
c
ccc        call prinm(rmlexp(iaddr(1,1)),nterms(0))
c
ccc        call h2dformmp(ier,scale(0),source,charge,n,
ccc     1  	center,nterms(0),mptemp)
c
ccc        call prinm(mptemp,nterms(0))
c
ccc        call h2dmperr(rmlexp(iaddr(1,1)),mptemp,nterms(0),d)
ccc        call prin2('error in upward pass=*',d,1)
c
ccc        pause
ccc        stop
c      END DEBUGGING SEGMENT
c------------------------------------------------------------
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(3)=t2-t1
c
        if (ifprint .ge. 1) 
     $     call prinf('=== STEP 4 (mp to lo) ===*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 4, convert multipole expansions into the local ones
c
cc        call prinf('laddr=*',laddr,2*(nlev+1))
cc        call prin2('bsize=*',bsize,(nlev+1))
cc        do 4200 ibox=1,nboxes
ccc        ntops=0

        call l2dterms_list2(epsfmm, itable, ier)
c        call prinf('itable=*',itable,7*7)
c
        do 4300 ilev=3,nlev+1
c        t3=second()
cC$        t3=omp_get_wtime()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,list,nlists,nlist,itype)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect2,radius)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,ftemp,htemp,cd,ilist)
C$OMP$PRIVATE(if_use_trunc,nterms_trunc,ii,jj) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1)
        do 4200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,15)
        endif
        level0=box(1)
        if (level0 .ge. 2) then
c
c       ... retrieve list #2
c
           itype=2
           call d2tgetl(ier,ibox,itype,list,nlist,wlists)
           if (ifprint .ge. 2) then
              call prinf('list2=*',list,nlist)
           endif
c
c       ... prune all sourceless boxes
c
ccc           if( box(10) .eq. 0 ) nlist=0

c       ... for all pairs in list #2, apply the translation operator 
c
           do 4150 ilist=1,nlist
              jbox=list(ilist)
              call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
              if( box1(10) .eq. 0 ) goto 4150
              if (jbox.eq.0) goto 4150
              if ((box(12).eq.0).and.(ifprune_list2.eq.1))
     $           goto 4150
c              radius = (corners1(1,1) - center1(1))**2
c              radius = radius + (corners1(2,1) - center1(2))**2
c              radius = sqrt(radius)
c
c       ... convert multipole expansions for all boxes in list 2 in local exp
c       ... if source is childless, evaluate directly (if cheaper)
c

              level1=box1(1)
c              ifdirect2 = 0
c              if( box1(10) .lt. (nterms(level1)+1)/2 .and. 
c     $            box(10) .lt. (nterms(level1)+1)/2  ) 
c     $               ifdirect2 = 1        
c
              ifdirect2 = 0
c
              if_use_trunc = 0
c
              if (ifdirect2 .eq. 0) then
              if( if_use_trunc .eq. 0) then

                 call l2dzero(mptemp,nterms(level1))
                 if( nterms(level0)+nterms(level1) .gt. 95 ) then
                 call l2dmploc(scale(level1),center1,
     1              rmlexp(iaddr(1,jbox)),nterms(level1),
     $              scale(level0),
     1              center0,mptemp,nterms(level0))
                 else
                 call l2dmploc_carray(scale(level1),center1,
     1              rmlexp(iaddr(1,jbox)),nterms(level1),
     $              scale(level0),
     1              center0,mptemp,nterms(level0),
     $              carray,ldc)
                 endif
                 call l2dadd(mptemp,rmlexp(iaddr(2,ibox)),
     1              nterms(level0))
c
c              call l2dmploc_add(scale(level1),center1,
c     $           rmlexp(iaddr(1,jbox)),nterms(level1),
c     $           scale(level0),center0,rmlexp(iaddr(2,ibox)),
c     $           nterms(level0))

              else

              ii=box1(2)-box(2)
              jj=box1(3)-box(3)
              nterms_trunc=itable(ii,jj)
              nterms_trunc=min(nterms(level0),nterms_trunc)
              nterms_trunc=min(nterms(level1),nterms_trunc)

                 call l2dzero(mptemp,nterms_trunc)
                 if( nterms_trunc+nterms_trunc .gt. 95 ) then
                 call l2dmploc(scale(level1),center1,
     1              rmlexp(iaddr(1,jbox)),nterms_trunc,
     $              scale(level0),
     1              center0,mptemp,nterms_trunc)
                 else
                 call l2dmploc_carray(scale(level1),center1,
     1              rmlexp(iaddr(1,jbox)),nterms_trunc,
     $              scale(level0),
     1              center0,mptemp,nterms_trunc,
     $              carray,ldc)
                 endif
                 call l2dadd(mptemp,rmlexp(iaddr(2,ibox)),
     1              nterms_trunc)

c              call l2dmploc_add_trunc(scale(level1),center1,
c     $           rmlexp(iaddr(1,jbox)),nterms(level1),nterms_trunc,
c     $           scale(level0),center0,rmlexp(iaddr(2,ibox)),
c     $           nterms(level0))

c              call l2dmploc_add(scale(level1),center1,
c     $           rmlexp(iaddr(1,jbox)),nterms_trunc,
c     $           scale(level0),center0,rmlexp(iaddr(2,ibox)),
c     $           nterms_trunc)

              endif
              endif

 4150       continue
        endif
 4200   continue
C$OMP END PARALLEL DO
c        t4=second()
cC$        t4=omp_get_wtime()
c        write(*,*) 'level ', ilev, ' time in list2:', t4-t3
ccc        write(*,*) 'time in list2:', second()-t1
ccc        write(*,*) 'ntops:', ntops
ccc        write(*,*) 'speed:', ntops/(second()-t1)
 4300   continue
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(4)=t2-t1
c       
        if (ifprint .ge. 1) 
     $     call prinf('=== STEP 5 (split lo) ===*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 5, split all local expansions
c
ccc        do 5200 ibox=1,nboxes
        do 5300 ilev=3,nlev
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(mptemp,lused,ier,i,j,ptemp,ftemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 5200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        call d2tnkids(box,nkids)
c       
        if (nkids .ne. 0) then
            level0=box(1)
            if (level0 .ge. 2) then
               if (ifprint .ge. 2) then
                  call prinf('ibox=*',ibox,1)
                  call prinf('box=*',box,15)
                  call prinf('nkids=*',nkids,1)
                  call prin2('center0=*',center0,2)
               endif
c
c       ... split local expansion of the parent box
c
               do 5100 i = 1,4
	          jbox = box(4+i)
	          if (jbox.eq.0) goto 5100
                  call d2tgetb(ier,jbox,box1,center1,corners1,wlists)
                  radius = (corners1(1,1) - center1(1))**2
                  radius = radius + (corners1(2,1) - center1(2))**2
                  radius = sqrt(radius)
                  if (ifprint .ge. 2) then
                     call prinf('jbox=*',jbox,1)
                     call prin2('radius=*',radius,1)
                     call prin2('center1=*',center1,2)
                  endif
                  level1=box1(1)
                  if( nterms(level0)+nterms(level1) .gt. 95 ) then
                  call l2dlocloc(scale(level0),center0,
     1               rmlexp(iaddr(2,ibox)),nterms(level0),
     1               scale(level1),center1,mptemp,nterms(level1))
                  else
                  call l2dlocloc_carray(scale(level0),center0,
     1               rmlexp(iaddr(2,ibox)),nterms(level0),
     1               scale(level1),center1,mptemp,nterms(level1),
     1               carray,ldc)
                  endif
                  call l2dadd(mptemp,rmlexp(iaddr(2,jbox)),
     1   	       nterms(level1))
 5100          continue
               if (ifprint .ge. 2) call prinf('=============*',x,0)
            endif
        endif
c
        if (nkids .ne. 0) then
            level=box(1)
            if (level .ge. 2) then
               if( ifprint .ge. 2 ) then
                  call prinf('ibox=*',ibox,1)
                  call prinf('box=*',box,15)
                  call prinf('nkids=*',nkids,1)
               endif
            endif
        endif
 5200   continue
C$OMP END PARALLEL DO
 5300   continue
c       
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(5)=t2-t1
c
        return
        end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the auxiliary routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine l2dpsort(n,isource,psort,pot)
        implicit real *8 (a-h,o-z)
        integer isource(*)
        complex *16 pot(*),psort(*)
c
ccc        call prinf('isource=*',isource,n)
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,n
        pot(isource(i))=psort(i)
        enddo
C$OMP END PARALLEL DO
c
        return
        end
c
c
c
c
c
        subroutine l2dfsort(n,isource,fldsort,fld)
        implicit real *8 (a-h,o-z)
        integer isource(*)
        complex *16 fld(2,*),fldsort(2,*)
c        
ccc        call prinf('isource=*',isource,n)
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,n
        fld(1,isource(i))=fldsort(1,i)
        fld(2,isource(i))=fldsort(2,i)
        enddo
C$OMP END PARALLEL DO
c
        return
        end
c
c
c
c
c
        subroutine l2dhsort(n,isource,hesssort,hess)
        implicit real *8 (a-h,o-z)
        integer isource(*)
        complex *16 hess(3,*),hesssort(3,*)
c        
ccc        call prinf('isource=*',isource,n)
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,n
        hess(1,isource(i))=hesssort(1,i)
        hess(2,isource(i))=hesssort(2,i)
        hess(3,isource(i))=hesssort(3,i)
        enddo
C$OMP END PARALLEL DO
c
        return
        end
c
c
c
c
c
        subroutine l2dreorder(nsource,source,
     $     ifcharge,charge,isource,ifdipole,
     1     dipstr,dipvec,sourcesort,chargesort,dipvecsort,dipstrsort) 
        implicit real *8 (a-h,o-z)
        real *8 source(2,*),sourcesort(2,*)
        integer isource(*)
        real *8 dipvec(2,*),dipvecsort(2,*)
        complex *16 charge(*),chargesort(*),dipstr(*),dipstrsort(*)
c       
ccc        call prinf('nsource=*',nsource,1)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i = 1,nsource
        sourcesort(1,i) = source(1,isource(i))
        sourcesort(2,i) = source(2,isource(i))
        if( ifcharge .ge. 1 ) then
        chargesort(i) = charge(isource(i))
        endif
        if (ifdipole .ge. 1) then
        dipstrsort(i) = dipstr(isource(i))
        dipvecsort(1,i) = dipvec(1,isource(i))
        dipvecsort(2,i) = dipvec(2,isource(i))
        endif
        enddo
C$OMP END PARALLEL DO
        return
        end
c
c
c
c
c
        subroutine l2dreordertarg(ntarget,target,itarget,targetsort)
        implicit real *8 (a-h,o-z)
        real *8 target(2,*),targetsort(2,*)
        integer itarget(*)
c       
ccc        call prinf('ntarget=*',ntarget,1)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i = 1,ntarget
        targetsort(1,i) = target(1,itarget(i))
        targetsort(2,i) = target(2,itarget(i))
        enddo
C$OMP END PARALLEL DO
        return
        end
c
c
c
c
c
        subroutine l2dzero(mpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       ... set multipole to zero
c
        complex *16 mpole(0:nterms)
c       
        do n=0,nterms
        mpole(n)=0
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l2dadd(mpole,mpole2,nterms)
        implicit real *8 (a-h,o-z)
        complex *16 mpole(0:nterms)
        complex *16 mpole2(0:nterms)
c       
        do n=0,nterms
        mpole2(n)=mpole2(n)+mpole(n)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l2dmpalloc(wlists,iaddr,nboxes,lmptot,nterms)
        implicit real *8 (a-h,o-z)
        integer box(20)
        integer nterms(0:*)
        integer iaddr(2,nboxes)
        real *8 center0(2),corners0(2,4)
        real *8 wlists(*)
c
c       ... construct pointer array iaddr for addressing multipole and
c       local expansion
c
        iptr=1
        do ibox=1,nboxes
        call d2tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
c
c       ... first, allocate memory for the multipole expansion
c       
        iaddr(1,ibox)=iptr
        iptr=iptr+(nterms(level)+1)*2
c
c       ... then, allocate memory for the local expansion
c       
        iaddr(2,ibox)=iptr
        iptr=iptr+(nterms(level)+1)*2
c       
        enddo
        lmptot = iptr
        return
        end
c
c
c
c
c
        subroutine l2d_init_carray(carray,ldc)
        implicit real *8 (a-h,o-z)
        real *8 carray(0:ldc,0:ldc)

        do l = 0,ldc
        carray(l,0) = 1.0d0
        enddo
        do m=1,ldc
        carray(m,m) = 1.0d0
        do l=m+1,ldc
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
        enddo
        enddo
c
        return
        end
        
cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date: 2010-12-28 14:28:30 -0500 (Tue, 28 Dec 2010) $
c    $Revision: 1571 $
c
c
c-----------------------------------------------------------------------------
c
c      l2dterms - determine number of terms in mpole expansions 
c
c      l2dterms_list2 - build the number of terms table for all boxes 
c           in list 2
c
c      l2dterms_list2w - build the number of terms table for all boxes 
c           in list 2, worst case error in multipole to local translation 
c
c      l2dterms_list2ew - build the number of terms table for all boxes
c           in extended list 2, worst case error in multipole to local
c           translation
c
c-----------------------------------------------------------------------------
c
c
c
      subroutine l2dterms(eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions.
c
c     The method is based on examining the decay of \rho^n / r^{n+1}
c     for \rho a worst case source and r a worst case target.
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      ier = 0
c
      ntmax = 1000
c       
      z1 = 1.5d0
      do i = 0,ntmax
         hfun(i) = 1.0d0/(z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          return
        endif
c
      enddo
      return
      end
c
c
c
c
c
      subroutine l2dterms_far(eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk=0. 
c
c     The method is based on examining the decay of h_n * j_n.
c
c     This routine assumes slightly larger separation of boxes: the
c     first unit box is located at the origin and the second box is
c     located at (3,0).
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      ier = 0
c
      ntmax = 1000
c
      z1 = 2.5d0
      do i = 0,ntmax
         hfun(i) = 1.0d0/(z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          return
        endif
c
      enddo
      return
      end
c
c
c
c
c
      subroutine l2dterms_list2(eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in list 2
c
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      integer nterms_table(2:3,0:3)
      integer itable(-3:3,-3:3)
c
      ier = 0
c
      do 1800 ii=2,3
      do 1800 jj=0,3
c
        dx=ii
        dy=jj
c       
        if( dx .gt. 0 ) dx=dx-.5
        if( dy .gt. 0 ) dy=dy-.5
c
        rr=sqrt(dx*dx+dy*dy)
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(3.0d0)/2*5,1)
c
      ntmax = 1000
c
      z1 = rr
      do i = 0,ntmax
         hfun(i) = 1.0d0/( z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          goto 1600
        endif
      enddo
 1600   continue
c
        nterms_table(ii,jj)=nterms
c
 1800   continue
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in list 2
c
        do i=-3,3
        do j=-3,3
        itable(i,j)=0
        enddo
        enddo
c
        do 2200 i=-3,3
        do 2200 j=-3,3
c
        if( abs(i) .gt. 1 ) then
        itable(i,j)=nterms_table(abs(i),abs(j))
        else if( abs(j) .gt. 1) then
        itable(i,j)=nterms_table(abs(j),abs(i))
        endif
c
 2200   continue
c
      return
      end
c
c
c
c
c
      subroutine l2dterms_list2w(eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in list 2
c     Estimate worst case multipole to local translation operator errors.
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      integer nterms_table(2:3,0:3)
      integer itable(-3:3,-3:3)
c
      ier = 0
c
      do 1800 ii=2,3
      do 1800 jj=0,3
c
        dx=ii
        dy=jj
c       
c        if( dx .gt. 0 ) dx=dx-.5
c        if( dy .gt. 0 ) dy=dy-.5
c
        rr=sqrt(dx*dx+dy*dy)
        rr=rr-sqrt(2.0d0)/2
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(3.0d0)/2*5,1)
c
      ntmax = 1000
c
      z1 = rr
      do i = 0,ntmax
         hfun(i) = 1.0d0/( z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          goto 1600
        endif
      enddo
 1600   continue
c
        nterms_table(ii,jj)=nterms
c
 1800   continue
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in list 2
c
        do i=-3,3
        do j=-3,3
        itable(i,j)=0
        enddo
        enddo
c
        do 2200 i=-3,3
        do 2200 j=-3,3
c
        if( abs(i) .gt. 1 ) then
        itable(i,j)=nterms_table(abs(i),abs(j))
        else if( abs(j) .gt. 1) then
        itable(i,j)=nterms_table(abs(j),abs(i))
        endif
c
 2200   continue
c
      return
      end
c
c
c
c
c
      subroutine l2dterms_list2e(eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in extended list 2
c
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      integer nterms_table(2:7,0:7)
      integer itable(-7:7,-7:7)
c
      ier = 0
c
      do 1800 ii=2,7
      do 1800 jj=0,7
c
        dx=ii
        dy=jj
c       
        if( dx .gt. 0 ) dx=dx-.5
        if( dy .gt. 0 ) dy=dy-.5
c
        rr=sqrt(dx*dx+dy*dy)
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(2.0d0)/2*5,1)
c
      ntmax = 1000
c
      z1 = rr
      do i = 0,ntmax
         hfun(i) = 1.0d0/( z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          goto 1600
        endif
      enddo
 1600   continue
c
        nterms_table(ii,jj)=nterms
c
 1800   continue
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in extended list 2
c
        do i=-7,7
        do j=-7,7
        itable(i,j)=0
        enddo
        enddo
c
        do 2200 i=-7,7
        do 2200 j=-7,7
c
        if( abs(i) .gt. 2 ) then
        itable(i,j)=nterms_table(abs(i),abs(j))
        else if( abs(j) .gt. 2) then
        itable(i,j)=nterms_table(abs(j),abs(i))
        endif
c
 2200   continue
c
      return
      end
c
c
c
c
c
      subroutine l2dterms_list2ew(eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in extended list 2
c     Estimate worst case multipole to local translation operator errors.
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      integer nterms_table(2:7,0:7)
      integer itable(-7:7,-7:7)
c
      ier = 0
c
      do 1800 ii=2,7
      do 1800 jj=0,7
c
        dx=ii
        dy=jj
c       
c        if( dx .gt. 0 ) dx=dx-.5
c        if( dy .gt. 0 ) dy=dy-.5
c
        rr=sqrt(dx*dx+dy*dy)
        rr=rr-sqrt(2.0d0)/2
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(2.0d0)/2*5,1)
c
      ntmax = 1000
c
      z1 = rr
      do i = 0,ntmax
         hfun(i) = 1.0d0/( z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          goto 1600
        endif
      enddo
 1600   continue
c
        nterms_table(ii,jj)=nterms
c
 1800   continue
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in extended list 2
c
        do i=-7,7
        do j=-7,7
        itable(i,j)=0
        enddo
        enddo
c
        do 2200 i=-7,7
        do 2200 j=-7,7
c
        if( abs(i) .gt. 2 ) then
        itable(i,j)=nterms_table(abs(i),abs(j))
        else if( abs(j) .gt. 2) then
        itable(i,j)=nterms_table(abs(j),abs(i))
        endif
c
 2200   continue
c
      return
      end
c
c
c
c
c
      subroutine l2dterms_eval(itype, eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      ier = 0
      ntmax = 1000
c
      z1 = 1.5d0
      do i = 0,ntmax
         hfun(i) = 1.0d0/(z1**(i+1))
      enddo
c
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
        z2 = dsqrt(2d0)/2.d0
c
c       corners included
        if( itype .eq. 1 ) z2 = dsqrt(2d0)/2.d0
c       edges included, no corners
        if( itype .eq. 2 ) z2 = dsqrt(1d0)/2.d0
c       center only
        if( itype .eq. 3 ) z2 = 1.0d0/2.d0
c       center only, small interior sphere
        if( itype .eq. 4 ) z2 = 0.8d0/2.d0
c
c
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          return
        endif
c
      enddo
      return
      end

*8cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c
c      This file contains the basic subroutines for 
c      forming and evaluating multipole (partial wave) expansions
c      in two dimensions.
c
c      Since log(z) is a multivalued complex function, we use
c      the real part Re(log(z)) = log(abs(z)) in all computations.
c
c      All multipole and local expansions are properly scaled 
c
c-----------------------------------------------------------------------
c
c      L2DFORMMP: creates multipole expansion (outgoing) due to 
c                 a collection of charge sources.
c
c      L2DFORMTA: creates local expansion due to 
c                 a collection of charge sources.
c
c      L2DFORMMP_DP: creates multipole expansion (outgoing) due to 
c                 a collection of dipole sources.
c
c      L2DFORMTA_DP: creates local expansion due to 
c                 a collection of dipole sources.
c
c-----------------------------------------------------------------------
c
c      Multipole and local translation operators
c
c      L2DMPMP:     converts multipole expansion to a multipole expansion.
c      L2DMPLOC:     converts multipole expansion to a local expansion.
c      L2DLOCLOC:     converts local expansion to a local expansion.
c
c-----------------------------------------------------------------------
c
c      Multipole and local translation operators with precomputed storage
c
c      L2DMPMP_CARRAY:   converts multipole expansion to a multipole expansion.
c      L2DMPLOC_CARRAY:   converts multipole expansion to a local expansion.
c      L2DLOCLOC_CARRAY:   converts local expansion to a local expansion.
c
c-----------------------------------------------------------------------
c       
c      Complex valued Cauchy sums
c       __ depreciated functions __
c
c      L2DMPEVALALL: computes potential and grad(potential)
c                 due to a multipole expansion for a collection of targets
c      L2DMPEVAL: computes potential and grad(potential)
c                 due to a multipole expansion.
c
c      L2DTAEVALALL: computes potential and grad(potential) 
c                  due to local expansion for a collection of targets
c      L2DTAEVAL: computes potential and grad(potential) 
c                  due to local expansion.
c
c      LPOTGRAD2DALL:  direct calculation for a collection of charge sources
c      LPOTGRAD2D : direct calculation for a single charge source
c
c      LPOTGRAD2DALL_DP:  direct calculation for a collection of dipoles
c      LPOTGRAD2D_DP : direct calculation for a single dipole
c
c
c      LPOTGRAD2DALL_SDP:  direct calculation for 
c                 a collection of charges and dipoles
c      LPOTGRAD2D_SDP : direct calculation for a single charge and a dipole
c
c-----------------------------------------------------------------------
c
c      Complex valued Cauchy sums
c
c      C2DMPEVALALL: computes potential and grad(potential)
c                 due to a multipole expansion for a collection of targets
c      C2DMPEVAL: computes potential and grad(potential)
c                 due to a multipole expansion.
c
c      C2DTAEVALALL: computes potential and grad(potential) 
c                  due to local expansion for a collection of targets
c      C2DTAEVAL: computes potential and grad(potential) 
c                  due to local expansion.
c
c      CPOTGRAD2DALL_SDP:  direct calculation for 
c                 a collection of charges and dipoles
c      CPOTGRAD2D_SDP : direct calculation for a single charge and a dipole
c
c      CPOTGRAD2D_SDP_SYM : direct calculation for
c      a pair of charges and a dipoles, uses symmetries
c
c-----------------------------------------------------------------------
c
c      Complex valued Laplace potentials
c
c      RCPOTGRAD2DALL_SDP:  direct calculation for 
c                 a collection of charges and dipoles
c      RCPOTGRAD2D_SDP : direct calculation for a single charge and a dipole
c
c-----------------------------------------------------------------------
c
c      Real valued Laplace potentials
c
c      RPOTGRAD2DALL_SDP:  direct calculation for 
c                 a collection of charges and dipoles
c      RPOTGRAD2D_SDP : direct calculation for a single charge and a dipole
c
c
c
c**********************************************************************
      subroutine l2dmpeval(rscale,center,mpole,nterms,ztarg,
     1      pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n / z^n  + mpole_0 log(abs(z))
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot,grad(2),hess(3),mpole(0:nterms)
        real *8 center(2),ztarg(2),zdiff(2)
        complex *16 ztemp1, ztemp2
ccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000)
c
        complex *16 ima,z,cd
        data ima/(0.0d0,1.0d0)/
c
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z=dcmplx(zdiff(1),zdiff(2))
c
c
        nmax = nterms + 3
ccc        allocate( z0pow(0:nmax) )
        ztemp1=rscale/z
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
        pot=mpole(0)*log(abs(z))
        do n=1,nterms
        cd=mpole(n)*z0pow(n)
        pot=pot+cd
        enddo
c
c
        if( ifgrad .eq. 1 ) then

        rinv=1/rscale
        cd=mpole(0)*z0pow(1)
        grad(1)=cd
        grad(2)=cd*ima
        do n=1,nterms
        cd=-mpole(n)*z0pow(n+1)*n
        grad(1)=grad(1)+cd
        grad(2)=grad(2)+cd*ima 
        enddo
        grad(1)=grad(1)*rinv
        grad(2)=grad(2)*rinv

        endif
c
c
        if( ifhess .eq. 1 ) then

        rinv2=1/rscale**2
        cd=-mpole(0)*z0pow(2)
        hess(1)=cd*(1)
        hess(2)=cd*(ima)
        hess(3)=cd*(-1)
        do n=1,nterms
        cd=mpole(n)*z0pow(n+2)*n*(n+1)
        hess(1)=hess(1)+cd*(1) 
        hess(2)=hess(2)+cd*(ima)
        hess(3)=hess(3)+cd*(-1) 
        enddo
        hess(1)=hess(1)*rinv2
        hess(2)=hess(2)*rinv2
        hess(3)=hess(3)*rinv2

        endif
c
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine l2dtaeval(rscale,center,mpole,nterms,ztarg,
     1      pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n * z^n  + mpole_0 
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the local expansion
c     ztarg  :    target location
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot,grad(2),hess(3),mpole(0:nterms)
        real *8 center(2),ztarg(2),zdiff(2)
        complex *16 ztemp1, ztemp2
ccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000)
c
        complex *16 ima,z,cd
        data ima/(0.0d0,1.0d0)/
c
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z=dcmplx(zdiff(1),zdiff(2))
c
c
c
        nmax = nterms
ccc        allocate( z0pow(0:nmax) )
        ztemp1=z/rscale
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
c
        pot=mpole(0)
        do n=1,nterms
        cd=mpole(n)*z0pow(n)
        pot=pot+cd
        enddo
c
c
        if( ifgrad .eq. 1 ) then

        grad(1)=0
        grad(2)=0

        rinv=1/rscale
        do n=1,nterms
        cd=mpole(n)*z0pow(n-1)*n
        grad(1)=grad(1)+cd
        grad(2)=grad(2)+cd*ima
        enddo
        grad(1)=grad(1)*rinv
        grad(2)=grad(2)*rinv

        endif
c
c
        if( ifhess .eq. 1 ) then

        hess(1)=0
        hess(2)=0
        hess(3)=0

        rinv2=1/rscale**2
        do n=2,nterms
        cd=mpole(n)*z0pow(n-2)*n*(n-1)
        hess(1)=hess(1)+cd*(1) 
        hess(2)=hess(2)+cd*(ima)
        hess(3)=hess(3)+cd*(-1)
        enddo
        hess(1)=hess(1)*rinv2
        hess(2)=hess(2)*rinv2
        hess(3)=hess(3)*rinv2

        endif
c
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine l2dmpevalall(rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n / z^n  + mpole_0 log(abs(z))
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot(1),grad(2,1),hess(3,1),mpole(0:nterms)
        complex *16 potloc,gradloc(2),hessloc(3)
        real *8 center(2),ztarg(2,ntarg),zdiff(2)
c
        complex *16 ima,z
        data ima/(0.0d0,1.0d0)/
c
        do i = 1,ntarg
        call l2dmpeval(rscale,center,mpole,nterms,ztarg(1,i),
     1     potloc,ifgrad,gradloc,ifhess,hessloc)
        if (ifpot.eq.1) pot(i) = pot(i) + potloc
        if (ifgrad.eq.1) then
        grad(1,i) = grad(1,i) + gradloc(1)
        grad(2,i) = grad(2,i) + gradloc(2)
        endif
        if (ifhess.eq.1) then
        hess(1,i) = hess(1,i) + hessloc(1)
        hess(2,i) = hess(2,i) + hessloc(2)
        hess(3,i) = hess(3,i) + hessloc(3)
        endif
        enddo
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine l2dtaevalall(rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n * z^n  + mpole_0 
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the local expansion
c     ztarg  :    target location
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot(1),grad(2,1),hess(3,1),mpole(0:nterms)
        complex *16 potloc,gradloc(2),hessloc(3)
        real *8 center(2),ztarg(2,ntarg),zdiff(2)
c
        complex *16 ima,z
        data ima/(0.0d0,1.0d0)/
c
        do i = 1,ntarg
        call l2dtaeval(rscale,center,mpole,nterms,ztarg(1,i),
     1     potloc,ifgrad,gradloc,ifhess,hessloc)
        if (ifpot.eq.1) pot(i) = pot(i) + potloc
        if (ifgrad.eq.1) then
        grad(1,i) = grad(1,i) + gradloc(1)
        grad(2,i) = grad(2,i) + gradloc(2)
        endif
        if (ifhess.eq.1) then
        hess(1,i) = hess(1,i) + hessloc(1)
        hess(2,i) = hess(2,i) + hessloc(2)
        hess(3,i) = hess(3,i) + hessloc(3)
        endif
        enddo
c
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformmp(ier,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_0  =  sum charge_j 
c                  j  
c
c     mpole_n  = -sum charge_j 1/n (z)^n /rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),charge(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do n=0,nterms
        mpole(n)=0
        enddo

        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
        mpole(0)=mpole(0)+charge(j)
c
        ztemp1=z0/rscale
        ztemp2=ztemp1
        do n=1,nterms
        mpole(n)=mpole(n)-charge(j)*ztemp1/n
        ztemp1=ztemp1*ztemp2
        enddo

        enddo

        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformta(ier,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_0  =  sum charge_j log(abs(z)) 
c                  j  
c
c     mpole_n  = -sum charge_j 1/n (1/z)^n *rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),charge(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do n=0,nterms
        mpole(n)=0
        enddo

        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
        mpole(0)=mpole(0)+charge(j)*log(abs(-z0))
c
        zinv=rscale/z0
        ztemp1=zinv
        do n=1,nterms
        mpole(n)=mpole(n)-charge(j)*ztemp1/n
        ztemp1=ztemp1*zinv
        enddo

        enddo

        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformta_add(ier,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_0  =  sum charge_j log(abs(z)) 
c                  j  
c
c     mpole_n  = -sum charge_j 1/n (1/z)^n *rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),charge(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
        mpole(0)=mpole(0)+charge(j)*log(abs(-z0))
c
        zinv=rscale/z0
        ztemp1=zinv
        do n=1,nterms
        mpole(n)=mpole(n)-charge(j)*ztemp1/n
        ztemp1=ztemp1*zinv
        enddo

        enddo

        return
        end
c
c
c
c
c
        subroutine l2dmpmp(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a multipole expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted multipole expansion
C           center2 = center of shifted multipole expansion
C           nterms2 = order of shifted multipole expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted multipole expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8, allocatable :: carray(:,:)
        complex *16, allocatable :: z0pow(:), z0powm(:)
ccc        complex *16, allocatable :: z0pow1(:), z0pow2(:)
        complex *16 z0pow1(0:1000), z0pow2(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
cc        done=1
cc        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = max(nterms1,nterms2)
        allocate( carray(0:nmax,0:nmax) )
c
        do l = 0,nmax
        carray(l,0) = 1.0d0
        enddo
        do m=1,nmax
        carray(m,m) = 1.0d0
        do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
        enddo
        enddo
c
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0pow1(0:nmax) )
        ztemp1=z0/rscale1
        ztemp2=ztemp1
        z0pow1(0)=1
        do i=1,nmax
        z0pow1(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
ccc        allocate( z0pow2(0:nmax) )
        ztemp1=z0/rscale2
        ztemp2=ztemp1
        z0pow2(0)=1
        do i=1,nmax
        z0pow2(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        jexp(0) = hexp(0)
        do i = 1,nterms2
        jexp(i) = jexp(i) - hexp(0)*z0pow2(i)/i
        do j = 1,min(i,nterms1)
        jexp(i) = jexp(i) +
     $     hexp(j)*z0pow2(i)/z0pow1(j)*carray(i-1,j-1)
        enddo
        enddo
c
ccc        call prin2('jexp=*',jexp,2*(nterms2+1))
c
        return
        end
c
c
c
c
c
        subroutine l2dmploc(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8, allocatable :: carray(:,:)
ccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000), z0powm(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = nterms1+nterms2
        allocate( carray(0:nmax,0:nmax) )
c
        do l = 0,nmax
        carray(l,0) = 1.0d0
        enddo
        do m=1,nmax
        carray(m,m) = 1.0d0
        do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
        enddo
        enddo
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0pow(0:nmax) )
        ztemp1=1/z0
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2        
        enddo
c
ccc        allocate( z0powm(0:nmax) )
        ztemp1=1/(-z0)
        ztemp2=ztemp1
        z0powm(0)=1
        do i=1,nmax
        z0powm(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        jexp(0) = hexp(0)*log(abs(-z0))
        do j = 1,nterms1
        jexp(0) = jexp(0) + hexp(j)*rscale1**j*z0powm(j)
        enddo

        do i = 1,nterms2
        jexp(i) = jexp(i) - hexp(0)/i
        do j = 1,nterms1
        jexp(i) = jexp(i) + hexp(j)*rscale1**j
     $     *z0powm(j)*carray(i+j-1,j-1)
        enddo
        jexp(i) = jexp(i)*z0pow(i)
        enddo
c
        do i = 1,nterms2
        jexp(i) = jexp(i) *rscale2**i
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l2dmploc_add(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,cd,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8, allocatable :: carray(:,:)
        complex *16, allocatable :: z0pow(:), z0powm(:)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = nterms1+nterms2
        allocate( carray(0:nmax,0:nmax) )
c
        do l = 0,nmax
        carray(l,0) = 1.0d0
        enddo
        do m=1,nmax
        carray(m,m) = 1.0d0
        do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
        enddo
        enddo
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
        allocate( z0pow(0:nmax) )
        ztemp1=1/z0
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
        allocate( z0powm(0:nmax) )
        ztemp1=1/(-z0)
        ztemp2=ztemp1
        z0powm(0)=1
        do i=1,nmax
        z0powm(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
ccc        do i = 0,nterms2
ccc        jexp(i) = 0
ccc        enddo
c
        jexp(0) = jexp(0) + hexp(0)*log(abs(-z0))
        do j = 1,nterms1
        jexp(0) = jexp(0) + hexp(j)*rscale1**j*z0powm(j)
        enddo

        do i = 1,nterms2
        cd = - hexp(0)/i
        do j = 1,nterms1
        cd = cd + hexp(j)*rscale1**j*z0powm(j)*carray(i+j-1,j-1)
        enddo
        jexp(i) = jexp(i) + cd*z0pow(i)*rscale2**i
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l2dlocloc(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts local expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original local expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original local expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8, allocatable :: carray(:,:)
        complex *16, allocatable :: z0pow(:), z0powm(:)
ccc        complex *16, allocatable :: z0powm1(:), z0powm2(:)
        complex *16 z0powm1(0:1000), z0powm2(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = max(nterms1,nterms2)
        allocate( carray(0:nmax,0:nmax) )
c
        do l = 0,nmax
        carray(l,0) = 1.0d0
        enddo
        do m=1,nmax
        carray(m,m) = 1.0d0
        do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
        enddo
        enddo
c
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0powm1(0:nmax) )
        ztemp1=(-z0)/rscale1
        ztemp2=ztemp1
        z0powm1(0)=1
        do i=1,nmax
        z0powm1(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
ccc        allocate( z0powm2(0:nmax) )
        ztemp1=(-z0)/rscale2
        ztemp2=ztemp1
        z0powm2(0)=1
        do i=1,nmax
        z0powm2(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        do i = 0,nterms2
        do j = i,nterms1
        jexp(i) = jexp(i) + hexp(j)
     $     *z0powm1(j)/z0powm2(i)*carray(j,i)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l2dmpmp_carray(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2,
     $     carray,ldc)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a multipole expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted multipole expansion
C           center2 = center of shifted multipole expansion
C           nterms2 = order of shifted multipole expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted multipole expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8 carray(0:ldc,0:ldc)
        complex *16, allocatable :: z0pow(:), z0powm(:)
ccc        complex *16, allocatable :: z0pow1(:), z0pow2(:)
        complex *16 z0pow1(0:1000), z0pow2(0:1000)
        real *8 rfactors(0:1000)
        complex *16 hexp1(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
cc        done=1
cc        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = max(nterms1,nterms2)
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0pow1(0:nmax) )
        ztemp1=1/(z0/rscale1)
        ztemp2=ztemp1
        z0pow1(0)=1
        do i=1,nmax
        z0pow1(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
ccc        allocate( z0pow2(0:nmax) )
        ztemp1=z0/rscale2
        ztemp2=ztemp1
        z0pow2(0)=1
        do i=1,nmax
        z0pow2(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
        rfactors(0)=1
        do i=1,max(nterms1,nterms2)
        rfactors(i)=rfactors(i-1)*(rscale1/rscale2)
        enddo
c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        do i = 0,nterms1
        hexp1(i) = hexp(i)*z0pow1(i)
        enddo
c
        jexp(0) = hexp(0)
        do i = 1,nterms2
        jexp(i) = jexp(i) - hexp1(0)/i
        do j = 1,min(i,nterms1)
        jexp(i) = jexp(i) + hexp1(j)*carray(i-1,j-1)
        enddo
        jexp(i)=jexp(i)*z0pow2(i)
        enddo
c
ccc        call prin2('jexp=*',jexp,2*(nterms2+1))
c
        return
        end
c
c
c
c
c
        subroutine l2dmploc_carray(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2,
     $     carray,ldc)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2,ztemp3
        real *8 carray(0:ldc,0:ldc)
cccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000), z0powm(0:1000)
        complex *16 hexp1(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
c        done=1
c        pi=4*atan(done)
c
ccc        nmax = nterms1+nterms2
        nmax = max(nterms1,nterms2)
c
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0pow(0:nmax) )
ccc        allocate( z0powm(0:nmax) )
        ztemp1=1/z0
        z0pow(0)=1
        z0powm(0)=1
        ztemp2=+ztemp1*rscale2
        ztemp3=-ztemp1*rscale1
        do i=1,nmax
        z0pow(i)=ztemp2
        z0powm(i)=ztemp3
        ztemp2=+ztemp2*ztemp1*rscale2
        ztemp3=-ztemp3*ztemp1*rscale1
        enddo

c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        do i = 0,nterms1
        hexp1(i) = hexp(i)*z0powm(i)
        enddo
c
        jexp(0) = hexp1(0)*log(abs(-z0))
        do j = 1,nterms1
        jexp(0) = jexp(0) + hexp1(j)
        enddo

        do i = 1,nterms2
        jexp(i) = jexp(i) - hexp1(0)/i
        do j = 1,nterms1
        jexp(i) = jexp(i) + hexp1(j)*carray(i+j-1,j-1)
        enddo
        jexp(i) = jexp(i)*z0pow(i) 
        enddo
c
ccc        call prin2('jexp=*',jexp,2*(nterms2+1))
c
        return
        end
c
c
c
c
c
        subroutine l2dlocloc_carray(
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2,
     $     carray,ldc)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts local expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           rscale1 = scaling parameter for original local expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original local expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,z0,ima, zmul,zinv,ztemp1,ztemp2
        real *8 carray(0:ldc,0:ldc)
ccc        complex *16, allocatable :: z0powm1(:), z0powm2(:)
        complex *16 z0powm1(0:ldc), z0powm2(0:ldc)
        complex *16 hexp1(0:1000)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        nmax = max(nterms1,nterms2)
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(-zdiff(1),-zdiff(2))
c
c
ccc        allocate( z0powm1(0:nmax) )
        ztemp1=(-z0)/rscale1
        ztemp2=ztemp1
        z0powm1(0)=1
        do i=1,nmax
        z0powm1(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
ccc        allocate( z0powm2(0:nmax) )
        ztemp1=1/((-z0)/rscale2)
        ztemp2=ztemp1
        z0powm2(0)=1
        do i=1,nmax
        z0powm2(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
        do i = 0,nterms2
        jexp(i) = 0
        enddo
c
        do i = 0,nterms1
        hexp1(i) = hexp(i)*z0powm1(i)
        enddo
c
        do i = 0,nterms2
        do j = i,nterms1
        jexp(i) = jexp(i) + hexp1(j)*carray(j,i)
        enddo
        jexp(i) = jexp(i)*z0powm2(i)
        enddo
c
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformmp_dp(ier,rscale,source,dipstr,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_0  =  0
c                 
c
c     mpole_n  =  sum dipstr_j (z_0)^(n-1)/z^n /rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),dipstr(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do n=0,nterms
        mpole(n)=0
        enddo

        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
c
        zmul=z0/rscale
        ztemp1=1/rscale
        do n=1,nterms
        mpole(n)=mpole(n)+dipstr(j)*ztemp1
        ztemp1=ztemp1*zmul
        enddo

        enddo

        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformta_dp(ier,rscale,source,dipstr,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_n  = -sum dipstr_j (z_0)^(n-1)/z^n /rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),dipstr(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do n=0,nterms
        mpole(n)=0
        enddo

        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
        zinv=rscale/z0
        ztemp1=zinv/rscale
        do n=0,nterms
        mpole(n)=mpole(n)-dipstr(j)*ztemp1
        ztemp1=ztemp1*zinv
        enddo

        enddo

        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine l2dformta_dp_add(ier,rscale,source,dipstr,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_n  = -sum dipstr_j (z_0)^(n-1)/z^n /rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),dipstr(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c
        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
        zinv=rscale/z0
        ztemp1=zinv/rscale
        do n=0,nterms
        mpole(n)=mpole(n)-dipstr(j)*ztemp1
        ztemp1=ztemp1*zinv
        enddo

        enddo

        return
        end
c
c
c
c
c
c**********************************************************************
c
c       Direct evaluation of Cauchy type sums, obsoleted version
c
c**********************************************************************
      subroutine lpotgrad2dall(ifgrad,ifhess,sources,charge,ns,
     1                   target,pot,grad,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     charges at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = log(abs(z))
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      complex *16 charge(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
         call lpotgrad2d(ifgrad,ifhess,sources(1,i),charge(i),target,
     1        potloc,gradloc,hessloc)
         pot = pot + potloc
         if (ifgrad.eq.1) then
         grad(1) = grad(1) + gradloc(1)
         grad(2) = grad(2) + gradloc(2)
         endif
         if (ifhess.eq.1) then
         hess(1) = hess(1) + hessloc(1)
         hess(2) = hess(2) + hessloc(2)
         hess(3) = hess(3) + hessloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotgrad2d(ifgrad,ifhess,source,charge,target,
     1                pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD
c     and Hessian HESS at the target point TARGET, due to a charge at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = log(abs(z))
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      real *8 source(2),target(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 charge
      complex *16 z, cd, zk, ima
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      pot=charge*log(abs(z))
c
      if (ifgrad.eq.1) then
         cd = charge/z
         grad(1) = cd
         grad(2) = cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = -charge/z**2
         hess(1) = cd
         hess(2) = cd*ima
         hess(3) = -cd
      endif
c
      return
      end
c
c
c
c**********************************************************************
      subroutine lpotgrad2dall_dp(ifgrad,ifhess,sources,dipstr,ns,
     1                   target,pot,grad,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of dipoles
c     at SOURCE(2,ns).  The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     dipstr        : dipole strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      complex *16 dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
         call lpotgrad2d_dp(ifgrad,ifhess,sources(1,i),dipstr(i),target,
     1        potloc,gradloc,hessloc)
         pot = pot + potloc
         if (ifgrad.eq.1) then
         grad(1) = grad(1) + gradloc(1)
         grad(2) = grad(2) + gradloc(2)
         endif
         if (ifhess.eq.1) then
         hess(1) = hess(1) + hessloc(1)
         hess(2) = hess(2) + hessloc(2)
         hess(3) = hess(3) + hessloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotgrad2d_dp(ifgrad,ifhess,source,dipstr,target,
     1                pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD
c     and Hessian HESS at the target point TARGET, due to a dipole at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     dipstr    : dipole strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      real *8 source(2),target(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 dipstr
      complex *16 z, cd, zk, ima
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      pot=dipstr/z
c
      if (ifgrad.eq.1) then
         cd = -dipstr/z**2
         grad(1) = cd
         grad(2) = cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = 2*dipstr/z**3
         hess(1) = cd
         hess(2) = cd*ima
         hess(3) = -cd
      endif
c
      return
      end
c
c
c
c**********************************************************************
      subroutine lpotgrad2dall_sdp(sources,ns,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of charges
c     and dipoles at SOURCE(2,ns).  The scaling is that required of the
c     delta function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot         : flag for computing potential
c	                 	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      complex *16 charge(ns),dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      if (ifpot.eq.1) pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
c        call lpotgrad2d_sdp(sources(1,i),
c     $     ifcharge,charge(i),ifdipole,dipstr(i),
c     1     target,ifpot,potloc,ifgrad,gradloc,ifhess,hessloc)
c         if (ifpot.eq.1) pot = pot + potloc
c         if (ifgrad.eq.1) then
c         grad(1) = grad(1) + gradloc(1)
c         grad(2) = grad(2) + gradloc(2)
c         endif
c         if (ifhess.eq.1) then
c         hess(1) = hess(1) + hessloc(1)
c         hess(2) = hess(2) + hessloc(2)
c         hess(3) = hess(3) + hessloc(3)
c         endif
        call lpotgrad2d_sdp_add(sources(1,i),
     $     ifcharge,charge(i),ifdipole,dipstr(i),
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotgrad2dall_sdp_add(sources,ns,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of charges
c     and dipoles at SOURCE(2,ns).  The scaling is that required of the
c     delta function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot         : flag for computing potential
c	                 	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      complex *16 charge(ns),dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
c      if (ifpot.eq.1) pot = 0.0d0
c      if (ifgrad.eq.1) then
c         grad(1) = 0.0d0
c         grad(2) = 0.0d0
c      endif
c      if (ifhess.eq.1) then
c         hess(1) = 0.0d0
c         hess(2) = 0.0d0
c         hess(3) = 0.0d0
c      endif
c
      do i = 1,ns
c        call lpotgrad2d_sdp(sources(1,i),
c     $     ifcharge,charge(i),ifdipole,dipstr(i),
c     1     target,ifpot,potloc,ifgrad,gradloc,ifhess,hessloc)
c         if (ifpot.eq.1) pot = pot + potloc
c         if (ifgrad.eq.1) then
c         grad(1) = grad(1) + gradloc(1)
c         grad(2) = grad(2) + gradloc(2)
c         endif
c         if (ifhess.eq.1) then
c         hess(1) = hess(1) + hessloc(1)
c         hess(2) = hess(2) + hessloc(2)
c         hess(3) = hess(3) + hessloc(3)
c         endif
        call lpotgrad2d_sdp_add(sources(1,i),
     $     ifcharge,charge(i),ifdipole,dipstr(i),
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine lpotgrad2d_sdp(source,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      if (ifpot.eq.1) then
         pot = 0
      endif
c
      if (ifgrad.eq.1) then
         grad(1) = 0
         grad(2) = 0
      endif
c
      if (ifhess.eq.1) then
         hess(1) = 0
         hess(2) = 0
         hess(3) = 0
      endif

        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=charge*log(abs(z))
c
      if (ifgrad.eq.1) then
         cd = charge*zinv
         grad(1) = cd
         grad(2) = cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = -charge*zinv2
         hess(1) = cd
         hess(2) = cd*ima
         hess(3) = -cd
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) pot=pot+dipstr*zinv
c
      if (ifgrad.eq.1) then
         cd = -dipstr*zinv2
         grad(1) = grad(1)+cd
         grad(2) = grad(2)+cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = 2*dipstr/z**3
         hess(1) = hess(1)+cd
         hess(2) = hess(2)+cd*ima
         hess(3) = hess(3)-cd
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
      subroutine lpotgrad2d_sdp_add(source,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = (d/dx, d/dy)
c		hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2)
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=pot+charge*log(abs(z))
c
      if (ifgrad.eq.1) then
         cd = charge*zinv
         grad(1) = grad(1)+cd
         grad(2) = grad(2)+cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = -charge*zinv2
         hess(1) = hess(1)+cd
         hess(2) = hess(2)+cd*ima
         hess(3) = hess(3)-cd
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) pot=pot+dipstr*zinv
c
      if (ifgrad.eq.1) then
         cd = -dipstr*zinv2
         grad(1) = grad(1)+cd
         grad(2) = grad(2)+cd*ima
      endif
c
      if (ifhess.eq.1) then
         cd = 2*dipstr*zinv*zinv2
         hess(1) = hess(1)+cd
         hess(2) = hess(2)+cd*ima
         hess(3) = hess(3)-cd
      endif
      
      endif
c
c
      return
      end
c
c
c
c
c
c**********************************************************************
c
c     Multipole and local expansion evaluation routines for Cauchy sums
c
c**********************************************************************
      subroutine c2dmpeval(rscale,center,mpole,nterms,ztarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n / z^n  + mpole_0 log(abs(z)), if ifpot = 1.
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifpot  :   flag controlling evaluation of potential:
c                   ifpot = 0, do not compute potential.
c                   ifpot = 1, compute potential.
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg (if requested)
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot,grad,hess,mpole(0:nterms)
        real *8 center(2),ztarg(2),zdiff(2)
        complex *16 ztemp1, ztemp2
ccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000)
c
        complex *16 ima,z,cd
        data ima/(0.0d0,1.0d0)/
c
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z=dcmplx(zdiff(1),zdiff(2))
c
c
        nmax = nterms + 3
ccc        allocate( z0pow(0:nmax) )
        ztemp1=rscale/z
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
        pot=mpole(0)*log(abs(z))
        do n=1,nterms
        pot=pot+mpole(n)*z0pow(n)
        enddo
c
c
        if( ifgrad .eq. 1 ) then

        rinv=1/rscale
        grad=mpole(0)*z0pow(1)
        do n=1,nterms
        grad=grad-mpole(n)*z0pow(n+1)*n
        enddo
        grad=grad*rinv

        endif
c
c
        if( ifhess .eq. 1 ) then

        rinv2=1/rscale**2
        hess=-mpole(0)*z0pow(2)
        do n=1,nterms
        hess=hess+mpole(n)*z0pow(n+2)*n*(n+1)
        enddo
        hess=hess*rinv2

        endif
c
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine c2dtaeval(rscale,center,mpole,nterms,ztarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n * z^n  + mpole_0, if ifpot = 1.
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the local expansion
c     ztarg  :    target location
c     ifpot  :   flag controlling evaluation of potential:
c                   ifpot = 0, do not compute potential.
c                   ifpot = 1, compute potential.
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg (if requested)
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot,grad,hess,mpole(0:nterms)
        real *8 center(2),ztarg(2),zdiff(2)
        complex *16 ztemp1, ztemp2
ccc        complex *16, allocatable :: z0pow(:), z0powm(:)
        complex *16 z0pow(0:1000)
c
        complex *16 ima,z,cd
        data ima/(0.0d0,1.0d0)/
c
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z=dcmplx(zdiff(1),zdiff(2))
c
c
c
        nmax = nterms
ccc        allocate( z0pow(0:nmax) )
        ztemp1=z/rscale
        ztemp2=ztemp1
        z0pow(0)=1
        do i=1,nmax
        z0pow(i)=ztemp1
        ztemp1=ztemp1*ztemp2
        enddo
c
c
c
        pot=mpole(0)
        do n=1,nterms
        pot=pot+mpole(n)*z0pow(n)
        enddo
c
c
        if( ifgrad .eq. 1 ) then

        grad=0

        rinv=1/rscale
        do n=1,nterms
        grad=grad+mpole(n)*z0pow(n-1)*n
        enddo
        grad=grad*rinv

        endif
c
c
        if( ifhess .eq. 1 ) then

        hess=0

        rinv2=1/rscale**2
        do n=2,nterms
        hess=hess+mpole(n)*z0pow(n-2)*n*(n-1)
        enddo
        hess=hess*rinv2
        endif
c
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine c2dmpevalall(rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n / z^n  + mpole_0 log(abs(z)) if ifgrad = 1.
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifpot  :   flag controlling evaluation of potential:
c                   ifpot = 0, do not compute potential.
c                   ifpot = 1, compute potential.
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg (if requested)
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot(1),grad(1),hess(1),mpole(0:nterms)
        complex *16 potloc,gradloc,hessloc
        real *8 center(2),ztarg(2,ntarg),zdiff(2)
c
        complex *16 ima,z
        data ima/(0.0d0,1.0d0)/
c
        do i = 1,ntarg
        call c2dmpeval(rscale,center,mpole,nterms,ztarg(1,i),
     1     ifpot,potloc,ifgrad,gradloc,ifhess,hessloc)
        if (ifpot.eq.1) pot(i) = pot(i) + potloc
        if (ifgrad.eq.1) then
        grad(i) = grad(i) + gradloc
        endif
        if (ifhess.eq.1) then
        hess(i) = hess(i) + hessloc
        endif
        enddo
c
        return
        end
c
c
c
c**********************************************************************
      subroutine c2dtaevalall(rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n * z^n  + mpole_0 if ifpot = 1.
c                n=1  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the local expansion
c     ztarg  :    target location
c     ifpot  :   flag controlling evaluation of potential:
c                   ifpot = 0, do not compute potential.
c                   ifpot = 1, compute potential.
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg (if requested)
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 pot(1),grad(1),hess(1),mpole(0:nterms)
        complex *16 potloc,gradloc,hessloc
        real *8 center(2),ztarg(2,ntarg),zdiff(2)
c
        complex *16 ima,z
        data ima/(0.0d0,1.0d0)/
c
        do i = 1,ntarg
        call c2dtaeval(rscale,center,mpole,nterms,ztarg(1,i),
     1     ifpot,potloc,ifgrad,gradloc,ifhess,hessloc)
        if (ifpot.eq.1) pot(i) = pot(i) + potloc
        if (ifgrad.eq.1) then
        grad(i) = grad(i) + gradloc
        endif
        if (ifhess.eq.1) then
        hess(i) = hess(i) + hessloc
        endif
        enddo
c
        return
        end
c
c
c
c
c
c**********************************************************************
c
c     Cauchy sums, direct evaluation routines
c
c**********************************************************************
      subroutine cpotgrad2dall_sdp(sources,ns,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     charges and dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot         : flag for computing potential
c	                 	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad,hess,potloc,gradloc,hessloc
      complex *16 h0,h1,cd,eye,z
      complex *16 charge(ns),dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      if (ifpot.eq.1) pot = 0.0d0
      if (ifgrad.eq.1) then
         grad = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess = 0.0d0
      endif
c
      do i = 1,ns
        call cpotgrad2d_sdp_add(sources(1,i),
     $     ifcharge,charge(i),ifdipole,dipstr(i),
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine cpotgrad2dall_sdp_add(sources,ns,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     charges and dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot         : flag for computing potential
c	                 	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 pot,grad,hess,potloc,gradloc,hessloc
      complex *16 h0,h1,cd,eye,z
      complex *16 charge(ns),dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      do i = 1,ns
        call cpotgrad2d_sdp_add(sources(1,i),
     $     ifcharge,charge(i),ifdipole,dipstr(i),
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine cpotgrad2d_sdp(source,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2)
      complex *16 pot,grad,hess
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      if (ifpot.eq.1) then
         pot = 0
      endif
c
      if (ifgrad.eq.1) then
         grad = 0
      endif
c
      if (ifhess.eq.1) then
         hess = 0
      endif

        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=charge*log(abs(z))
c
      if (ifgrad.eq.1) then
         grad = charge*zinv
      endif
c
      if (ifhess.eq.1) then
         hess = -charge*zinv2
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) pot=pot+dipstr*zinv
c
      if (ifgrad.eq.1) then
         grad = grad-dipstr*zinv2
      endif
c
      if (ifhess.eq.1) then
         hess = hess+2*dipstr*zinv2*zinv
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
      subroutine cpotgrad2d_sdp_add(source,
     $     ifcharge,charge,ifdipole,dipstr,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2)
      complex *16 pot,grad,hess
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=pot+charge*log(abs(z))
c
      if (ifgrad.eq.1) then
         grad = grad+charge*zinv
      endif
c
      if (ifhess.eq.1) then
         hess = hess-charge*zinv2
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) pot=pot+dipstr*zinv
c
      if (ifgrad.eq.1) then
         grad = grad-dipstr*zinv2
      endif
c
      if (ifhess.eq.1) then
         hess = hess+2*dipstr*zinv2*zinv
      endif
      
      endif
c
c
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine cpotgrad2d_sdp_sym(source1,source2,
     $     ifcharge,charge1,charge2,ifdipole,dipstr1,dipstr2,
     1     ifpot,pot1,pot2,ifgrad,grad1,grad2,ifhess,hess1,hess2)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source1(2),source2(2)
      complex *16 pot1,grad1,hess1
      complex *16 pot2,grad2,hess2
      complex *16 charge1,dipstr1
      complex *16 charge2,dipstr2
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=source2(1)-source1(1)
      ydiff=source2(2)-source1(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      if (ifpot.eq.1) then
         pot1 = 0
         pot2 = 0
      endif
c
      if (ifgrad.eq.1) then
         grad1 = 0
         grad2 = 0
      endif
c
      if (ifhess.eq.1) then
         hess1 = 0
         hess2 = 0
      endif

        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) then
        cd = log(abs(z))
        pot2 = charge1*cd
        pot1 = charge2*cd
      endif
c
      if (ifgrad.eq.1) then
         cd = zinv
         grad2 = +charge1*cd
         grad1 = -charge2*cd
      endif
c
      if (ifhess.eq.1) then
         cd = zinv2
         hess2 = -charge1*cd
         hess1 = -charge2*cd
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd = zinv
         pot2 = pot2+dipstr1*cd
         pot1 = pot1-dipstr2*cd
      endif
c
      if (ifgrad.eq.1) then
         cd = zinv2
         grad2 = grad2-dipstr1*cd
         grad1 = grad1-dipstr2*cd
      endif
c
      if (ifhess.eq.1) then
         cd = 2*zinv2*zinv
         hess2 = hess2+dipstr1*cd
         hess1 = hess1-dipstr2*cd
      endif
      
      endif
c
c
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine cpotgrad2d_sdp_sym_add(source1,source2,
     $     ifcharge,charge1,charge2,ifdipole,dipstr1,dipstr2,
     1     ifpot,pot1,pot2,ifgrad,grad1,grad2,ifhess,hess1,hess2)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = charge log(abs(z)) + dipstr 1/z
c		grad = gradient = d/dz
c		hess = Hessian = d^2/dz^2
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source1(2),source2(2)
      complex *16 pot1,grad1,hess1
      complex *16 pot2,grad2,hess2
      complex *16 charge1,dipstr1
      complex *16 charge2,dipstr2
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=source2(1)-source1(1)
      ydiff=source2(2)-source1(2)
ccc      rr=xdiff*xdiff+ydiff*ydiff
ccc      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
        zinv=1/z
        zinv2=zinv*zinv

      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) then
        cd = log(abs(z))
        pot2 = pot2+charge1*cd
        pot1 = pot1+charge2*cd
      endif
c
      if (ifgrad.eq.1) then
         cd = zinv
         grad2 = grad2+charge1*cd
         grad1 = grad1-charge2*cd
      endif
c
      if (ifhess.eq.1) then
         cd = zinv2
         hess2 = hess2-charge1*cd
         hess1 = hess1-charge2*cd
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd = zinv
         pot2 = pot2+dipstr1*cd
         pot1 = pot1-dipstr2*cd
      endif
c
      if (ifgrad.eq.1) then
         cd = zinv2
         grad2 = grad2-dipstr1*cd
         grad1 = grad1-dipstr2*cd
      endif
c
      if (ifhess.eq.1) then
         cd = 2*zinv2*zinv
         hess2 = hess2+dipstr1*cd
         hess1 = hess1-dipstr2*cd
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
c
c       Complex valued Laplace particle direct evaluation routines
c
c**********************************************************************
      subroutine rcpotgrad2d_sdp(source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c     pot = charge log(r) + dipstr (dipvec \dot \grad log(r) )
c     grad = gradient = (d/dx, d/dy) 
c     hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2) 
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     dipstr    : dipole strength
c     dipvec    : dipole orientation vector
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2),dipvec(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
c
      z=dcmplx(xdiff,ydiff)
c
      if (ifpot.eq.1) then
         pot = 0
      endif
c
      if (ifgrad.eq.1) then
         grad(1) = 0
         grad(2) = 0
      endif
c
      if (ifhess.eq.1) then
         hess(1) = 0
         hess(2) = 0
         hess(3) = 0
      endif

c
      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=charge*log(r)
c
      if (ifgrad.eq.1) then
         cd = charge/rr
         grad(1) = xdiff*cd
         grad(2) = ydiff*cd
      endif
c
      if (ifhess.eq.1) then
         cd = charge/rr**2
         hess(1) = cd*(rr-2*xdiff*xdiff)
         hess(2) = cd*(  -2*xdiff*ydiff)
         hess(3) = cd*(rr-2*ydiff*ydiff)
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd=dipstr/rr
         pot=pot-cd*(xdiff*dipvec(1)+ydiff*dipvec(2))
      endif
c
      if (ifgrad.eq.1) then
         cd = dipstr/rr**2
         derxx=+rr-2*xdiff*xdiff
         derxy=   -2*xdiff*ydiff
         deryy=+rr-2*ydiff*ydiff
         grad(1) = grad(1)-cd*(derxx*dipvec(1)+derxy*dipvec(2))
         grad(2) = grad(2)-cd*(derxy*dipvec(1)+deryy*dipvec(2))
      endif
c
      if (ifhess.eq.1) then
         cd = dipstr/rr**3
         derxxx=-6*xdiff*rr+8*xdiff**3
         derxxy=-2*ydiff*rr+8*xdiff**2*ydiff
         derxyy=-2*xdiff*rr+8*xdiff*ydiff**2
         deryyy=-6*ydiff*rr+8*ydiff**3
         hess(1) = hess(1)-cd*(derxxx*dipvec(1)+derxxy*dipvec(2))
         hess(2) = hess(2)-cd*(derxxy*dipvec(1)+derxyy*dipvec(2))
         hess(3) = hess(3)-cd*(derxyy*dipvec(1)+deryyy*dipvec(2))
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
      subroutine rcpotgrad2d_sdp_add(source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c     pot = charge log(r) + dipstr (dipvec \dot \grad log(r) )
c     grad = gradient = (d/dx, d/dy) 
c     hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2) 
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     dipstr    : dipole strength
c     dipvec    : dipole orientation vector
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2),dipvec(2)
      complex *16 pot,grad(2),hess(3)
      complex *16 charge,dipstr
      complex *16 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
c
c
      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=pot+charge*log(r)
c
      if (ifgrad.eq.1) then
         cd = charge/rr
         grad(1) = grad(1) + xdiff*cd
         grad(2) = grad(2) + ydiff*cd
      endif
c
      if (ifhess.eq.1) then
         cd = charge/rr**2
         hess(1) = hess(1) + cd*(rr-2*xdiff*xdiff)
         hess(2) = hess(2) + cd*(  -2*xdiff*ydiff)
         hess(3) = hess(3) + cd*(rr-2*ydiff*ydiff)
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd=dipstr/rr
         pot=pot-cd*(xdiff*dipvec(1)+ydiff*dipvec(2))
      endif
c
      if (ifgrad.eq.1) then
         cd = dipstr/rr**2
         derxx=+rr-2*xdiff*xdiff
         derxy=   -2*xdiff*ydiff
         deryy=+rr-2*ydiff*ydiff
         grad(1) = grad(1)-cd*(derxx*dipvec(1)+derxy*dipvec(2))
         grad(2) = grad(2)-cd*(derxy*dipvec(1)+deryy*dipvec(2))
      endif
c
      if (ifhess.eq.1) then
         cd = dipstr/rr**3
         derxxx=-6*xdiff*rr+8*xdiff**3
         derxxy=-2*ydiff*rr+8*xdiff**2*ydiff
         derxyy=-2*xdiff*rr+8*xdiff*ydiff**2
         deryyy=-6*ydiff*rr+8*ydiff**3
         hess(1) = hess(1)-cd*(derxxx*dipvec(1)+derxxy*dipvec(2))
         hess(2) = hess(2)-cd*(derxxy*dipvec(1)+derxyy*dipvec(2))
         hess(3) = hess(3)-cd*(derxyy*dipvec(1)+deryyy*dipvec(2))
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
c
c       Real valued Laplace particle direct evaluation routines
c
c**********************************************************************
      subroutine rpotgrad2d_sdp(source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c     pot = charge log(r) + dipstr (dipvec \dot \grad log(r) )
c     grad = gradient = (d/dx, d/dy) 
c     hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2) 
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     dipstr    : dipole strength
c     dipvec    : dipole orientation vector
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2),dipvec(2)
      real *8 pot,grad(2),hess(3)
      real *8 charge,dipstr
      real *8 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
c
      if (ifpot.eq.1) then
         pot = 0
      endif
c
      if (ifgrad.eq.1) then
         grad(1) = 0
         grad(2) = 0
      endif
c
      if (ifhess.eq.1) then
         hess(1) = 0
         hess(2) = 0
         hess(3) = 0
      endif

c
      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=charge*log(r)
c
      if (ifgrad.eq.1) then
         cd = charge/rr
         grad(1) = xdiff*cd
         grad(2) = ydiff*cd
      endif
c
      if (ifhess.eq.1) then
         cd = charge/rr**2
         hess(1) = cd*(rr-2*xdiff*xdiff)
         hess(2) = cd*(  -2*xdiff*ydiff)
         hess(3) = cd*(rr-2*ydiff*ydiff)
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd=dipstr/rr
         pot=pot-cd*(xdiff*dipvec(1)+ydiff*dipvec(2))
      endif
c
      if (ifgrad.eq.1) then
         cd = dipstr/rr**2
         derxx=+rr-2*xdiff*xdiff
         derxy=   -2*xdiff*ydiff
         deryy=+rr-2*ydiff*ydiff
         grad(1) = grad(1)-cd*(derxx*dipvec(1)+derxy*dipvec(2))
         grad(2) = grad(2)-cd*(derxy*dipvec(1)+deryy*dipvec(2))
      endif
c
      if (ifhess.eq.1) then
         cd = dipstr/rr**3
         derxxx=-6*xdiff*rr+8*xdiff**3
         derxxy=-2*ydiff*rr+8*xdiff**2*ydiff
         derxyy=-2*xdiff*rr+8*xdiff*ydiff**2
         deryyy=-6*ydiff*rr+8*ydiff**3
         hess(1) = hess(1)-cd*(derxxx*dipvec(1)+derxxy*dipvec(2))
         hess(2) = hess(2)-cd*(derxxy*dipvec(1)+derxyy*dipvec(2))
         hess(3) = hess(3)-cd*(derxyy*dipvec(1)+deryyy*dipvec(2))
      endif
      
      endif
c
c
      return
      end
c
c
c
c**********************************************************************
      subroutine rpotgrad2d_sdp_add(source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c     pot = charge log(r) + dipstr (dipvec \dot \grad log(r) )
c     grad = gradient = (d/dx, d/dy) 
c     hess = Hessian = (d^2/dx^2, d^2/dxdy, d^2/dy^2) 
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing potential
c	                 	ifpot = 0 -> don't compute 
c		                ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     dipstr    : dipole strength
c     dipvec    : dipole orientation vector
c     target    : location of the target
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(2),target(2),dipvec(2)
      real *8 pot,grad(2),hess(3)
      real *8 charge,dipstr
      real *8 z, cd, zk, ima, zinv, zinv2
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
c
c
      if( ifcharge .eq. 1 ) then
c
      if (ifpot.eq.1) pot=pot+charge*log(r)
c
      if (ifgrad.eq.1) then
         cd = charge/rr
         grad(1) = grad(1) + xdiff*cd
         grad(2) = grad(2) + ydiff*cd
      endif
c
      if (ifhess.eq.1) then
         cd = charge/rr**2
         hess(1) = hess(1) + cd*(rr-2*xdiff*xdiff)
         hess(2) = hess(2) + cd*(  -2*xdiff*ydiff)
         hess(3) = hess(3) + cd*(rr-2*ydiff*ydiff)
      endif
c
      endif
c
c
c
      if( ifdipole .eq. 1 ) then

      if (ifpot.eq.1) then
         cd=dipstr/rr
         pot=pot-cd*(xdiff*dipvec(1)+ydiff*dipvec(2))
      endif
c
      if (ifgrad.eq.1) then
         cd = dipstr/rr**2
         derxx=+rr-2*xdiff*xdiff
         derxy=   -2*xdiff*ydiff
         deryy=+rr-2*ydiff*ydiff
         grad(1) = grad(1)-cd*(derxx*dipvec(1)+derxy*dipvec(2))
         grad(2) = grad(2)-cd*(derxy*dipvec(1)+deryy*dipvec(2))
      endif
c
      if (ifhess.eq.1) then
         cd = dipstr/rr**3
         derxxx=-6*xdiff*rr+8*xdiff**3
         derxxy=-2*ydiff*rr+8*xdiff**2*ydiff
         derxyy=-2*xdiff*rr+8*xdiff*ydiff**2
         deryyy=-6*ydiff*rr+8*ydiff**3
         hess(1) = hess(1)-cd*(derxxx*dipvec(1)+derxxy*dipvec(2))
         hess(2) = hess(2)-cd*(derxxy*dipvec(1)+derxyy*dipvec(2))
         hess(3) = hess(3)-cd*(derxyy*dipvec(1)+deryyy*dipvec(2))
      endif
      
      endif
c
c
      return
      end
c
c
c
