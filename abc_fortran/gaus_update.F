      subroutine gaus_update
      use constants
      use variables
      use gaus

      IMPLICIT DOUBLE PRECISION(A-H,O-Z),INTEGER(I-N)


      print *, ngaus,ng

      if(ng.lt.ngaus) then ! add gaussian
!         do i=1,na
!            r0gx(i,ng+1)=x(i,2)+dble(nx(i))*hmat(1,1)
!            r0gy(i,ng+1)=y(i,2)+dble(ny(i))*hmat(2,2)
!            r0gz(i,ng+1)=z(i,2)+dble(nz(i))*hmat(3,3)
!         enddo        
            r0gx(1,ng+1)=x(1,2)
            r0gy(1,ng+1)=y(1,2)

      endif



      end
