      program main
      implicit none
      double precision :: x(1000),y(1000),e(1000),xt,yt,et
      integer :: n_c,n_f ! represent coarse and fine, respectively
      integer :: i,j,k,m,p,char

      open(1,file='length')
      read(1,*) n_c
      close(1)

      open(2,file='coarse_min.dat')
      m=1;p=1
      read(2,*) char,x(m),y(m),e(m)

      do while (p.lt.n_c)
         p=p+1
         read(2,*) char,xt,yt,et
         do i=1,m
            if (abs(xt-x(i)).le.0.001.and.abs(yt-y(i)).le.0.001) goto 10
         end do
         m=m+1
         x(m)=xt;y(m)=yt;e(m)=et
 10   end do

      print *,p,m
      
      open(3,file='fine_min.dat')
      do i=1,m
         write(3,*) x(i),y(i),e(i)
      end do
      
      stop
      end
