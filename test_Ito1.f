      program main
          use ito1
          implicit none

      real, dimension(:,:), allocatable :: R
      real, dimension(:,:,:), allocatable :: phi
      integer, dimension(:), allocatable :: iflags
      real, dimension(:,:,:), allocatable :: rs
      real, dimension(:), allocatable :: pgamma
      real, dimension(:,:), allocatable :: p
      integer :: i
      real, dimension(:,:,:,:), allocatable :: dphi
      real, dimension(:,:,:), allocatable :: phi1, phi2
      real, dimension(:), allocatable :: a, b, c, d
      real :: v
      real, dimension(:,:), allocatable :: eps
      real, dimension(:,:), allocatable :: sigma
      real :: lam, mu, kA
      real, dimension(:), allocatable :: x
      real, parameter :: pi = 4.0*atan(1.0)
      character(len=128):: fname
      real, dimension(:,:), allocatable :: rc
      real, dimension(:,:), allocatable :: tshear
      real, dimension(:,:), allocatable :: rshear

      allocate(R(3,3))
      allocate(phi(3,3,3))
      allocate(iflags(12))
      allocate(rs(12,3,3))
      allocate(pgamma(24))
      allocate(p(3,3))
      allocate(dphi(3,3,3,3))
      allocate(phi1(3,3,3), phi2(3,3,3))
      allocate(a(3), b(3), c(3), d(3))
      allocate(eps(3,3), sigma(3,3))
      allocate(x(3))
      allocate(rc(3,3))
      allocate(tshear(3,3))
      allocate(rshear(3,3))

      print *, " euler z1x2z3 " 
      call euler_zxz(R, 0.0, 45.0, 90.0)
      write(6, "(3F10.4)") R
      
      lam = 3.0
      mu = 1.0

      kA = (lam+mu)/8.0/pi/mu/(lam+2.0*mu)
      print *, " kelvin potential "
      call kelvin_potential(phi, kA, 3.0, 1.0, 
     1       0.0, 0.0, 0.0,
     2       1.0, 1.0, 1.0)
      write(6, "(3F10.4)") phi
      print *, " check kelvin potential "
      call kelvin_potential(phi, kA, 3.0, 1.0, 
     1       1.0, 1.0, 1.0,
     2       0.0, 0.0, 0.0)
      write(6, "(3F10.4)") phi
      print *, " check kelvin potential "
      call kelvin_potential(phi, kA, 3.0, 1.0, 
     1       1.0, 0.0, 0.0,
     2       0.0, 0.0, 0.0)
      write(6, "(3F10.4)") phi

      iflags = (/1, -1, -1,
     1          -1, -1, -1,
     2          -1, -1, -1,
     3          -1, -1, -1/)

      rs(1,:,:) = reshape((/1, 0, 0,
     1              0, 1, 0,
     2              0, 0, 1/), (/3,3/))

      do i = 2, 12
        rs(i,:,:) = rs(1,:,:)
      end do

      pgamma = 0
      pgamma(1) = 1

      rs(1,:,:) = R

      print *, " plastic strain "
      call pstrain(p, pgamma, rs, iflags, r)
      write(6, "(3F10.4)") p

      print *, " kelvin potential gradient "
      call  dkelvin_potential(dphi, kA, 3.0, 1.0, 
     1       0.0, 0.0, 0.0,
     2       1.0, 1.0, 1.0)
      write(6, "(27F10.4)") dphi

      print *, " d kelvin potential  "
      call kelvin_potential(phi1, kA, 3.0, 1.0, 
     1       0.0, 0.0, 0.0,
     2       1.0-1e-3, 1.0, 1.0)
      call kelvin_potential(phi2, kA, 3.0, 1.0, 
     1       0.0, 0.0, 0.0,
     2       1.0+1e-3, 1.0, 1.0)
      write(6, "(27F10.4)") (phi2-phi1)/2e-3

      call kelvin_potential(phi1, kA, 3.0, 1.0, 
     1       0.0, 0.0, 0.0,
     2       1.0, 1.0-1e-3, 1.0)
      call kelvin_potential(phi2, kA, 3.0, 1.0, 
     1       0.0, 0.0, 0.0,
     2       1.0, 1.0+1e-3, 1.0)
      write(6, "(27F10.4)") (phi2-phi1)/2e-3

      call kelvin_potential(phi1, kA, 3.0, 1.0, 
     1       0.0, 0.0, 0.0,
     2       1.0, 1.0, 1.0-1e-3)
      call kelvin_potential(phi2, kA, 3.0, 1.0, 
     1       0.0, 0.0, 0.0,
     2       1.0, 1.0, 1.0+1e-3)
      write(6, "(27F10.4)") (phi2-phi1)/2e-3

      print *, " tet volume "
      v = 0
      a = (/1.0, 0.0, 0.0/)
      b = (/0.0, 1.0, 0.0/)
      c = (/0.0, 0.0, 1.0/)
      d = (/0.0, 0.0, 0.0/)
      v = orient3dfast(a, b, c, d) 
      write(6, "(F10.4)") v

      print *, " triangle area "
      v = triangle_area(sqrt(2.0), sqrt(2.0), sqrt(2.0))
      write(6, "(F10.4)") v

      print *, " tri area "
      v = triarea(a, b, c)
      write(6, "(F10.4)") v

      print *, " face normal "
      call facenormal(a, b, c, d, 1)
      write(6, "(3F10.4)") d 
      
      print *, " body force " 
      lam = 3.0
      mu = 1.0
      x(1) = 0.3
      x(2) = 0.4
      x(3) = 0.5
      kA = (lam+mu)/8.0/pi/mu/(lam+2.0*mu)
      call body_force(eps, kA, lam, mu, 
     1       x(1), x(2), x(3),
     2       a, b, c, d, 
     3       pgamma, rs, iflags, r)
      write(6, "(3F10.4)") eps

      write(fname, '(a128)') "cube1x1x1"
      print *, " init ", fname
      call init(fname, 1)
      
      sigma = 0
      sigma(2,3) = 1
      sigma(3,2) = 1
      !call gshear(eps, sigma, lam, mu, 10, 1)     
      print *, " sorting "
      rs = 0   
      call sort4plane(pgamma, rs, iflags, rc, sigma, 5, 1)     
      write(6, '(4F10.4)') pgamma(1:4)
      write(6, *) iflags
      print *, " slip system " 
      write(6, "(9F10.4)") rs(1,1:3,1:3)
      write(6, "(9F10.4)") rs(4,1:3,1:3)
 
      pgamma(5:8) = 0.1
      call pstrain(p, pgamma(5:8), rs, iflags, rc)
      print *, " single crystal " 
      write(6, '(3F10.4)') 2*mu*p
      print *, " its plastic strain "
      write(6, '(3F10.4)') p


      pgamma = 0.1
      call gshear(tshear, rshear, pgamma, sigma, lam, mu, 5, 1)
      print *, " poly crystal "
      write(6, '(3F10.4)') tshear
      print *, "  stress "
      write(6, '(3F10.4)') -tshear-2*mu*p
      print *, " homogeneous stress "
      write(6, '(3F10.4)') sigma
      print *, "  stress "
      write(6, '(3F10.4)') rshear
      

      end program
