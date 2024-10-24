      program main
      implicit none
      double precision :: x(100),y(100),w(100),h(100)
      double precision :: r,xt,yt,d,ht,wt
      integer :: i,j,n,m
      double precision :: xc,yc

      call random_seed()
      xc=10.0;yc=0.0
      x(1)=0.0;y(1)=0.0;w(1)=1.0;h(1)=2.0
      x(2)=20.0;y(2)=0.0;w(2)=1.7;h(2)=3.5
      n=2
      do while (n.lt.100)
! First generate the coordinate
 10      call random_number(r)
         xt=20.0*r
         call random_number(r)
         yt=10.0*r-5.0
         d=sqrt((xt-xc)**2.0+(yt-yc)**2.0)
         if (d.ge.10.0) goto 10 ! out of range

! Second, generate the height and width
         call random_number(r)
         ht=0.5+1.2*r  ! height from 0.5 to 1.7 eV
         call random_number(r)
         wt=0.5+0.5*r  ! width from 0.5 to 1.0

! Then judge whether too close
         do i=1,n
            d=sqrt((xt-x(i))**2.0+(yt-y(i))**2.0)
            if (d.lt.(wt+w(i))/1.5) goto 10
         end do
         n=n+1
         x(n)=xt;y(n)=yt;w(n)=wt;h(n)=ht
      end do
      open(1,file='PES_info')
      do i=1,100
         write(1,'(4f16.8)') x(i),y(i),h(i),w(i)
      end do
      
      stop
      end program
