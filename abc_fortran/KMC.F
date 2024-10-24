      program main
      implicit none
      double precision :: bar(1000,1000),r,tmin,tbar
      double precision :: T=1000.0,kb=1.38d-23,ev=1.6d-19,k0=1.0d13
      double precision :: time,prob(1000),sum_prob
      integer :: p(1000,1000),num(1000)
      integer :: n_c,n_f ! represent coarse and fine, respectively
      integer :: i,j,k,m,char
      integer :: init,final

      open(1,file='length')
      read(1,*) n_f
      close(1)

      call random_seed()
      p=0;num=0
      open(2,file='barMatrix.dat')
      do i=1,n_f
         read(2,*) char,tmin,m
         num(i)=m
         do j=1,num(i)
            read(2,*) p(i,j),r,tbar
            bar(i,j)=tbar-tmin
         end do
         read(2,*)
      end do

      open(3,file='results')
      do k=1,11

      time=0.0
      do m=1,20
      init=1;final=2
      i=init
      do while (i.ne.final)
         sum_prob=0.0;prob=0.0
         do j=1,num(i)
            prob(j)=k0*exp(-bar(i,j)*ev/kb/T)
            sum_prob=sum_prob+prob(j)
         end do
         call random_number(r)
         do j=1,num(i)
            prob(j)=prob(j)/sum_prob
            if (r.lt.prob(j)) goto 10
            r=r-prob(j)
         end do
 10      i=p(i,j)
         time=time+1.0/sum_prob
      end do
      end do
      write(3,*) 1.0/T, log(time/20.0)
      T=T-50.0
      end do

      stop
      end
