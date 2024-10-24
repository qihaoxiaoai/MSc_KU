!-----------------------------------------------------------
!     Simulation code Ising 2010
!-----------------------------------------------------------

      module constants
      integer :: na ! number of atoms
      integer :: nfix=0,ncalc,nn
      integer :: iix=1,iiy=1,iiz=10 ! Unit cell repetition
      double precision :: J = 1.0d0 ! interaction strength
      double precision :: wm = 1.1d3 ! inertia of atom
      double precision :: dt = 0.01d0 !time step
      double precision :: nrelax=10000 ! simulation steps
      double precision :: t=0.01d0 !temperature
      double precision :: am = 1.6726485d-27 ! Mass of proton [kg]
      double precision :: bm = 1.6749543d-27 ! Mass of neutron [kg]
      double precision :: rs, rs2 ! Bookkeeping radius
      double precision :: rc, rc2 ! Cutoff-radius
      double precision :: bl=1.3806503d-23, blt !Boltzman constant  
      double precision :: pi=3.14159265d0,u(6)
      integer :: mr = 521232427 !48828215
      integer :: jjj ! system's time for random seed
      integer, dimension(:), allocatable :: seed
      double precision :: ang=1.0d-10 ![m/A]
      double precision :: ev=1.60217733d-19 ![J/eV] 
      character(1) :: num(0:9)
      character(10) :: number
      character(1) :: char
      character(9) :: char9
      end module

      module variables
      integer :: mm1,mm2,mm3,ii(10),iwrite
      double precision :: ss1,ss2,ss3
      double precision :: a(3,3),b(3,3)
      double precision, allocatable :: x(:,:),y(:,:),z(:,:) ! atom position
      double precision, allocatable :: negx(:,:),negy(:,:),negz(:,:) ! negtive start position
      double precision, allocatable :: xt(:),yt(:),zt(:),it(:)
      double precision, allocatable :: xsad(:),ysad(:),zsad(:),pot_atom(:)
      double precision, allocatable :: fx(:),fy(:),fz(:) ! force on atom
      double precision, allocatable :: fxo(:),fyo(:),fzo(:)
      double precision, allocatable :: px(:),py(:),pz(:) ! kinetic moment
      double precision, allocatable :: vx(:),vy(:),vz(:) ! velocity
      integer, allocatable :: ni(:),nij(:,:),nij1(:,:) ! BK param.
      integer, allocatable :: nij2(:,:),nij3(:,:) ! BK param.
      integer, allocatable :: nx(:),ny(:),nz(:),ns(:)
      integer, allocatable :: nxc(:),nyc(:),nzc(:),ifix(:)
      double precision :: tots,tots1,avsp,fluc,etot,aven,fluce ! total spin and average spin
      double precision :: temperature ! temperature
      double precision :: stress(6)
      double precision, allocatable :: stress_i(:,:)
      
      ! for neb use
      double precision :: r_img(0:20,2),t_img(20,2)
      double precision :: g_img(20,2),g_per(20,2)
      double precision :: delta=0.001,torr=1.0E-5
      double precision :: basin_x(100),basin_y(100),ebas(100)
      logical :: judge=.true.

      ! for creating PES
      integer :: n_pes=100
      double precision :: pes_x(100),pes_y(100),pes_a(100),pes_r(100)

      ! for neighboring matrix
      double precision :: sad(100,100),ensad
      integer :: con(100,100),neighb(100,100)
      end module

      module eam
      integer :: mesh=500000,m
      double precision :: pairmin,pairmax,densmin,densmax
      double precision :: embedmin,embedmax
      double precision, allocatable :: pairx(:),densx(:),embedx(:)
      double precision, allocatable:: pairy(:),densy(:),embedy(:)
      double precision, allocatable :: paird(:),pair2d(:)
      double precision, allocatable :: densd(:),dens2d(:)
      double precision, allocatable :: embedd(:),embed2d(:)
      double precision, allocatable :: alpx(:),aldx(:),alex(:)
      double precision, allocatable :: alpy(:),aldy(:),aley(:)
      double precision, allocatable :: pairc(:,:)
      double precision :: df(2)
      integer :: iopt(2)
      end module

      module gaus
      integer :: ngaus=0,ng
      double precision :: dgg2=0.4d0,wg=0.05d0,wg_1=0.0,dgg2_1=1.0   
      double precision, allocatable :: r0gx(:,:),r0gy(:,:),r0gz(:,:)
      end module

      program main
      use constants
      use variables
      use eam
      use gaus


      IMPLICIT DOUBLE PRECISION(A-H,O-Z),INTEGER(I-N)

      print '("-----------------------------------------------------")'
      print '("|                   Ising program                    ")'
      print '("-----------------------------------------------------")' 

      num(0)='0'
      num(1)='1'
      num(2)='2'
      num(3)='3'
      num(4)='4'
      num(5)='5'
      num(6)='6'
      num(7)='7'
      num(8)='8'
      num(9)='9'
      

      call random_seed()
      call initpos  ! initial angle distribution
      call gaus_setup
      allocate(xsad(na),ysad(na),zsad(na))
      open(1,file='gaus_wh')
      read(1,*) dgg2,wg
      close(1)


      ncalc=na

      open(55,file='conMatrix.dat')
      open(66,file='barMatrix.dat')
      iwrite=0
!      write(55,*) iwrite,tots

      open(1,file='length')
      read(1,*) n_min
      close(1)

      n_pes=100
      open(1,file='PES_info')
      do i=1,n_pes
         read(1,*) pes_x(i),pes_y(i),pes_a(i),pes_r(i)
      end do
      close(1)

      open(1,file='fine_min.dat')
      do i=1,n_min
         read(1,*) basin_x(i),basin_y(i),ebas(i)
      end do

      con=0;neighb=0
      do i=1,n_min
         do j=1,i
            dis=(basin_x(i)-basin_x(j))**2.0+(basin_y(i)-basin_y(j))**2.0
            if (dis.le.16.0 .and. i.ne.j) then
               con(i,j)=1;con(j,i)=1
! First create images config
               r_img(0,1)=basin_x(i);r_img(0,2)=basin_y(i)
               r_img(20,1)=basin_x(j);r_img(20,2)=basin_y(j)
               do k=1,19
                  r_img(k,1)=r_img(0,1)+dble(k)/20.0*(r_img(20,1)-r_img(0,1))
                  r_img(k,2)=r_img(0,2)+dble(k)/20.0*(r_img(20,2)-r_img(0,2))
               end do
! Then call neb subroutine
               call neb
               sad(i,j)=ensad;sad(j,i)=ensad
            end if
         end do
      end do

      do i=1,n_min
         ij=0
         do j=1,n_min
            if (con(i,j).eq.1) then
               ij=ij+1
               neighb(i,ij)=j
            end if
         end do
         write(55,*) i,basin_x(i),basin_y(j),ij
         write(66,*) i,ebas(i),ij
         do ijk=1,ij
            write(55,*) neighb(i,ijk),basin_x(neighb(i,ijk)),basin_y(neighb(i,ijk))
            write(66,*) neighb(i,ijk),ebas(neighb(i,ijk)),sad(i,neighb(i,ijk))
         end do
         write(55,*)
         write(66,*)
      end do

      stop

      do i=4,4 !n_min
         do j=1,1 !1,i
            if (con(i,j).eq.1) then
! First create images config
               r_img(0,1)=basin_x(i);r_img(0,2)=basin_y(i)
               r_img(20,1)=basin_x(j);r_img(20,2)=basin_y(j)
               do k=1,19
                  r_img(k,1)=r_img(0,1)+dble(k)/20.0*(r_img(20,1)-r_img(0,1))
                  r_img(k,2)=r_img(0,2)+dble(k)/20.0*(r_img(20,2)-r_img(0,2))
               end do
! Then call neb subroutine
               call neb
            end if
         end do
      end do

      stop
            


      do ij=1,1
         x(1,2)=pes_x(ij);y(1,2)=pes_y(ij)
         x(1,1)=x(1,2);y(1,1)=y(1,2)

      do ig=0,ngaus
         pot_m=-1000000000.0d0
         do irelax=1,nrelax

            call calc_force
            enpot0=etot
            if(enpot0.gt.pot_m) then
               pot_m=enpot0
               xsad(1)=x(1,2)
               ysad(1)=y(1,2)
            end if
            ng=ig

            if (ng.gt.0) call add_gaus
!            print *, 'after',etot
            call verlet

            if(mod(irelax,50).eq.0) then
               iwrite=irelax
!               call write_pos
               write(55,'(i10,3d16.8)') irelax,x(1,2),y(1,2),etot
             end if
         enddo
         write(60,'(i10,3d16.8)') ng,etot,enpot0,pot_m
         if ((etot-enpot0).le.1.0d-5) then
            write(65,*) x(1,2),y(1,2)
            write(66,'(i5,6d16.8)') ng,enpot0,pot_m,x(1,2),y(1,2),xsad(1),ysad(1)
         end if
         write(70,'(i10,4d16.8)') ng,x(1,2),y(1,2),xsad(1),ysad(1)
         ng=ig
         if (abs(x(1,2)-19.97).le.0.05 .and. abs(y(1,2)).le.0.03) then
!            write(85,*) pot_m,1.0
            stop
         end if
!         if (abs(x(1,2)+1.97).le.0.05) then
!            write(85,*) pot_m,2.0
!            stop
!         end if
         if (ng.lt.ngaus) call gaus_update
         call velocity_mod
       enddo
!       write(1,*) x(1,2),y(1,2),etot
       end do
!       close(1)
       open(1,file='mesh')
       do i=1,100
          do j=1,100
             x(1,2)=-2.0+dble(i-1)*0.25
             y(1,2)=-12.5+dble(j-1)*0.25
             call calc_force
             write(1,*) x(1,2),y(1,2),etot
          end do
          write(1,*)
       end do
          
      print *, "ok",na
      end
