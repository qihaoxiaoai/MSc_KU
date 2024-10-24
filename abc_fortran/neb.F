      subroutine neb
      use constants
      use variables
      IMPLICIT DOUBLE PRECISION(A-H,O-Z),INTEGER(I-N)

      l=0;judge=.true.
      do while (judge .and. l.le.50000)
         l=l+1
         judge=.false.
         do i=1,19
            t_img(i,1)=r_img(i+1,1)-r_img(i-1,1)
            t_img(i,1)=t_img(i,1)/sqrt((r_img(i+1,1)-r_img(i-1,1))**2.+(r_img(i+1,2)-r_img(i-1,2))**2.)
            t_img(i,2)=r_img(i+1,2)-r_img(i-1,2)
            t_img(i,2)=t_img(i,2)/sqrt((r_img(i+1,1)-r_img(i-1,1))**2.+(r_img(i+1,2)-r_img(i-1,2))**2.)
            x(1,2)=r_img(i,1);y(1,2)=r_img(i,2)
            call calc_force
            g_img(i,1)=-fx(1);g_img(i,2)=-fy(1)
            g_per(i,1)=g_img(i,1)-(g_img(i,1)*t_img(i,1)+g_img(i,2)*t_img(i,2))*t_img(i,1)
            g_per(i,2)=g_img(i,2)-(g_img(i,1)*t_img(i,1)+g_img(i,2)*t_img(i,2))*t_img(i,2)
            if (sqrt(g_per(i,1)**2.+g_per(i,2)**2.).ge.torr) judge=.true.
            r_img(i,1)=r_img(i,1)-g_per(i,1)*delta
            r_img(i,2)=r_img(i,2)-g_per(i,2)*delta
         end do
      end do

      ensad=-1000.0
      do i=0,20
         x(1,2)=r_img(i,1);y(1,2)=r_img(i,2)
         call calc_force
         ensad=max(ensad,etot)
      end do
      print *, "converged after",l,"steps",ensad
      end
