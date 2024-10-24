      subroutine add_gaus
      use constants
      use variables
      use gaus

      IMPLICIT DOUBLE PRECISION(A-H,O-Z),INTEGER(I-N)


      do ig=1,ng
         delg=wg
         sum=0.0d0
         xx=x(1,2)-r0gx(1,ig)
         yy=y(1,2)-r0gy(1,ig)
         sum=sum+xx*xx+yy*yy
         delg=delg*(dexp(-xx**2/2.0d0/dgg2) &
         *dexp(-yy**2/2.0d0/dgg2))
         sum=sum/3.0d0
 !        sum=sum/dble(na)/3.0d0
         if(sum.gt.4.0d0*dgg2) delg=0.0d0
!         print *, ig, delg,xx,yy
!         print *,'before',etot,etot+delg
         etot=etot+delg

         delfx=delg/dgg2*xx !*dexp(-xx**2/2.0d0/dgg2)*dexp(-yy**2/2.0d0/dgg2)*(dexp(dexp(-xx**2/2.0d0/dgg2)*dexp(-yy**2/2.0d0/dgg2))-1.0)
         delfy=delg/dgg2*yy !*dexp(-xx**2/2.0d0/dgg2)*dexp(-yy**2/2.0d0/dgg2)*(dexp(dexp(-xx**2/2.0d0/dgg2)*dexp(-yy**2/2.0d0/dgg2))-1.0)
               
         fx(1)=fx(1)+delfx
         fy(1)=fy(1)+delfy
      enddo


      end
