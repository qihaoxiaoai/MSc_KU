      subroutine initpos
      use constants
      use variables
           
      IMPLICIT DOUBLE PRECISION(A-H,O-Z),INTEGER(I-N)
      integer fanfan
      integer ind(2050)
      real r

      call random_seed()
      id=0
      na=iix*iiy
      allocate(x(na,3),y(na,3),z(na,3))
      allocate(fx(na),fy(na),fz(na))
      allocate(fxo(na),fyo(na),fzo(na))
      allocate(px(na),py(na),pz(na),vx(na),vy(na),vz(na))

      open(222,file='xy_pos') 
      read(222,*) x(1,2),y(1,2)
      x(1,1)=x(1,2);y(1,1)=y(1,2)
      close(222)

      end
