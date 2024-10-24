      subroutine velocity
      use constants
      use variables

      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z),INTEGER(I-N)
      double precision ::  r1,r2,r


      open(222,file='vxvy_init')
      read(222,*) vx(1),vy(1)
      close(222)
      x(1,2)=x(1,2)-vx(1)*dt
      y(1,2)=y(1,2)-vy(1)*dt
      end


