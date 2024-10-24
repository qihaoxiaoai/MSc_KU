      subroutine gaus_setup
      use gaus
      use constants
      use variables

      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z),INTEGER(I-N)


      allocate(r0gx(na,ngaus*2),r0gy(na,ngaus*2),r0gz(na,ngaus*2))


      end
