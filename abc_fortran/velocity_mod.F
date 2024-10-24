      subroutine velocity_mod
      use constants
      use variables

      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z),INTEGER(I-N)

!      call random_seed()

      x(1,1)=x(1,2)
      y(1,1)=y(1,2)

      call random_number(r1)
      vx(1)=0.001*(2.0*r1-1.0)
      call random_number(r2)
      vy(1)=0.001*(2.0*r2-1.0)

      x(1,2)=x(1,2)-vx(1)*dt
      y(1,2)=y(1,2)-vy(1)*dt


      end



