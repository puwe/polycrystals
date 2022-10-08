      module ito1
          implicit none

          real, dimension(:,:), allocatable :: node
          integer, dimension(:,:), allocatable :: ele
          integer, dimension(:,:), allocatable :: face
          integer, dimension(:,:), allocatable :: edge

          real, dimension(:,:), allocatable :: aele

          integer, dimension(:,:), allocatable :: fcc_slipnor 
          integer, dimension(:,:), allocatable :: fcc_primary
          integer, dimension(:,:), allocatable :: fcc_conjugate
          integer, dimension(:,:), allocatable :: fcc_crossglide
          integer, dimension(:,:), allocatable :: fcc_critical 
          
          interface cross_product
              module procedure cross_productReal
              module procedure cross_productInteger
          end interface
      contains
          subroutine init(f, num)
          character(len=128), intent(in) :: f
          integer, intent(in) :: num

          character(len=256) :: fname
          integer :: n, d, a, b
          integer :: en, ed, ea
          integer :: fn, fb
          integer :: egn, egb
          integer :: aen, aed, aea

          integer :: i, m 

          !allocate(fname(256))
          write(fname, '(a, a, i0, a)') f, '.', num, '.node'
          fname = adjustl(fname)
          !print *, fname
          open(1, file=fname, status="old",
     1      action="read")
          read(1, *) n, d, a, b
          allocate(node(n, d+a+b))
          node = 0
          do i = 1, n
            read(1, *) m, node(i, :)
          end do
          close(1)
          
          write(fname, '(a, a, i0, a)') f, '.', num, '.ele'
          fname = adjustl(fname)
          open(2, file=fname, status="old",
     1      action="read")
          read(2, *) en, ed, ea
          allocate(ele(en, ed+ea))
          ele = 0
          do i = 1, en
            read(2, *) m, ele(i, :)
          end do
          close(2)

          write(fname, '(a, a, i0, a)') f, '.', num, '.face'
          fname = adjustl(fname)
          open(3, file=fname, status="old",
     1      action="read")
          read(3, *) fn, fb
          ! print *, fn, fb
          allocate(face(fn, 3+fb))
          face = 0
          do i = 1, fn
            read(3, *) m, face(i, :)
            ! print *, m, face(i, :)
          end do
          close(3)
          
          write(fname, '(a, a, i0, a)') f, '.', num, '.edge'
          fname = adjustl(fname)
          open(4, file=fname, status="old",
     1       action="read")
          read(4, *) egn, egb
          allocate(edge(egn, 2+egb))
          edge = 0
          do i = 1, egn
            read(4, *) m, edge(i, :)
          end do
          close(4)

          
          write(fname, '(a, a, i0, a)') 'Ito', '.', num, '.angle'  
          fname = adjustl(fname)
          open(5, file=fname, status='old',
     1      action='read')
          read(5, *) aen, aed, aea
          allocate(aele(aen, aed+aea))
          aele = 0
          do i = 1, aen
            read(5, *) m, aele(i, :)
            print *, m, aele(i,:)
          end do
          close(5)

              allocate(fcc_slipnor(3, 4))
              allocate(fcc_primary(3, 3))
              allocate(fcc_conjugate(3, 3))
              allocate(fcc_crossglide(3, 3))
              allocate(fcc_critical(3, 3))


          
              fcc_slipnor = reshape((/
     1                  1, 1, -1, 
     2                  1, -1, 1,
     3                  1, -1, -1,
     4                  1,  1,  1 
     5                 /), shape(fcc_slipnor))
                  
              fcc_primary = reshape((/
     1                   1, 0, 1,
     2                  0, -1, -1,
     3                  -1, 1, 0
     4                  /), shape(fcc_primary))

              fcc_conjugate = reshape((/
     1                  -1, -1, 0,
     2                  0, 1, 1,
     3                  1, 0, -1
     4                  /), shape(fcc_conjugate))

              fcc_crossglide = reshape((/
     1                  -1, 0, -1,
     2                  1, 1, 0,
     3                  0, -1, 1
     4                  /), shape(fcc_crossglide))

              fcc_critical = reshape((/
     1                  -1, 0, 1,
     2                  0, 1, -1,
     3                  1, -1, 0
     4                  /), shape(fcc_critical))
 
          !deallocate(fname)
          end subroutine

          ! axb = c
          subroutine cross_productInteger(c, b, a)
          integer, dimension(:), intent(in) :: b, a
          integer, dimension(:), intent(out) :: c

          c(1) = a(2)*b(3) - a(3)*b(2)
          c(2) = a(3)*b(1) - a(1)*b(3)
          c(3) = a(1)*b(2) - a(2)*b(1)

          end subroutine


          subroutine slip_system(slipsys, slipdir, slipnor)
          integer, dimension(:), intent(in) :: slipdir, slipnor
          real, dimension(:,:), intent(out) :: slipsys

          integer, dimension(3) :: sliplin
          real :: er, es, et
          er = 0
          es = 0
          et = 0

          call cross_product(sliplin, slipdir, slipnor)

          er = sqrt(real(dot_product(slipnor, slipnor)))
          es = sqrt(real(dot_product(slipdir, slipdir)))
          et = sqrt(real(dot_product(sliplin, sliplin)))

          slipsys(:,1) = slipnor/er
          slipsys(:,2) = slipdir/es
          slipsys(:,3) = sliplin/et 

          end subroutine



          subroutine euler_zxz(R, alpha1, beta, alpha2)
          real, intent(in) :: alpha1, beta, alpha2
          real, dimension(:,:), intent(out) :: R 

          real :: a1, b2, a3
          real, parameter :: pi = 4.0*atan(1.0)
          
          a1 = alpha1*pi/180
          b2 = beta*pi/180
          a3 = alpha2*pi/180

          R(1,1) = cos(a1)*cos(a3) - cos(b2)*sin(a1)*sin(a3)
          R(2,1) = cos(a3)*sin(a1) + cos(a1)*cos(b2)*sin(a3)
          R(3,1) = sin(b2)*sin(a3)

          R(1,2) = -cos(a1)*sin(a3) - cos(b2)*cos(a3)*sin(a1)
          R(2,2) = cos(a1)*cos(b2)*cos(a3) - sin(a1)*sin(a3)
          R(3,2) = cos(a3)*sin(b2)

          R(1,3) = sin(a1)*sin(b2)
          R(2,3) = -cos(a1)*sin(b2)
          R(3,3) = cos(b2)

          end subroutine

          subroutine kelvin_potential(phi, A, lam, mu, 
     1       x1, x2, x3,
     2       y1, y2, y3)
          real, intent(in) :: A, lam, mu
          real, intent(in) :: x1, x2, x3
          real, intent(in) :: y1, y2, y3
          real, dimension(:,:,:), intent(out) :: phi

          integer :: i, j, k
          real :: r1, r2
          real :: c1, c2
          real :: d1, d2
          real, dimension(:), allocatable :: x
          real, dimension(:), allocatable :: y
          allocate(x(3))
          allocate(y(3))
          x = 0
          y = 0

          r1 = 0
          r2 = 0
          c1 = 0
          c2 = 0
          d1 = 0
          d2 = 0

          r2 = (x1-y1)**2 + (x2-y2)**2 + (x3-y3)**2
          r1 = sqrt(r2)

          if (r1 .le. 0) then
              print *, " singularity at "
              write(6, "(F10.3)") x1, x3, x3
              call abort
          end if

          c1 = -6.0*mu*A
          c2 = 2.0*mu**2*A/(lam+mu)
          
          x = (/x1, x2, x3/)
          y = (/y1, y2, y3/)


          do i = 1, 3
            do j = 1, 3
              do k = 1, 3

                phi(i, j, k) = c1/r1**5*
     &              (x(i)-y(i))*(x(j)-y(j))*(x(k)-y(k))
                if (i .eq. j) then
                    phi(i, j, k) = phi(i, j, k) + c2/r1**3*(x(k)-y(k))
                end if
                if (i .eq. k) then
                    phi(i, j, k) = phi(i, j, k) - c2/r1**3*(x(j)-y(j))
                end if
                if (j .eq. k) then
                    phi(i, j, k) = phi(i, j, k) - c2/r1**3*(x(i)-y(i))
                end if

              end do
            end do
          end do


          deallocate(x)
          deallocate(y)
          end subroutine

          subroutine pstrain(p, pgamma, rs, iflags, rc)
          real, dimension(:), intent(in) :: pgamma
          real, dimension(:,:,:), intent(in) :: rs
          integer, dimension(:), intent(in) :: iflags
          real, dimension(:,:), intent(in) :: rc

          real, dimension(:,:), intent(out) :: p

          integer, dimension(:), allocatable :: nslipsys
          real, dimension(:,:), allocatable :: r
          real, dimension(:), allocatable :: en, em, el
          integer :: n
          integer :: i,j,k,l

          n = 0
          do i = 1, 12
            if (iflags(i) >= 0) then
                n = n + 1
            end if
          end do

          allocate(nslipsys(n))
          allocate(r(3,3))
          allocate(en(3), em(3), el(3))

          nslipsys = 0
          r = 0
          en = 0
          em = 0
          el = 0

          n = 0
          do i = 1, 12
            if (iflags(i) >= 0) then
                n = n + 1
                ! index from 0
                nslipsys(n) = i-1
            end if
          end do

          p = 0
          do j = 1, n
            r = rs(j,1:3,1:3)
            en = r(:,1)
            em = r(:,2)
            el = r(:,3)

            do k = 1, 3
              do l = 1, 3                
                p(k, l) = p(k, l) + pgamma(j)*
     &           0.5*(em(k)*en(l) + em(l)*en(k)) 
              end do
            end do
          end do

          r = matmul(matmul(transpose(rc), p), rc)
          p = r

          deallocate(nslipsys)
          deallocate(r)
          deallocate(en, em, el)
          end subroutine

          subroutine dkelvin_potential(dphi, A, lam, mu, 
     1       x1, x2, x3,
     2       y1, y2, y3)
          real, intent(in) :: A, lam, mu
          real, intent(in) :: x1, x2, x3
          real, intent(in) :: y1, y2, y3
          real, dimension(:,:,:,:), intent(out) :: dphi

          integer :: i, j, k, l
          real :: r1, r2
          real :: c1, c2
          real :: d1, d2
          real, dimension(:), allocatable :: x
          real, dimension(:), allocatable :: y
          real :: dphi1
          allocate(x(3))
          allocate(y(3))
          x = 0
          y = 0

          r1 = 0
          r2 = 0
          c1 = 0
          c2 = 0
          d1 = 0
          d2 = 0

          r2 = (x1-y1)**2 + (x2-y2)**2 + (x3-y3)**2
          r1 = sqrt(r2)

          if (r1 .le. 0) then
              print *, " singularity at "
              write(6, "(F10.3)") x1, x3, x3
              call abort
          end if

          c1 = -6.0*mu*A
          c2 = 2.0*mu**2*A/(lam+mu)
          d1 = 30*mu*A
          d2 = -6.0*mu**2*A/(lam+mu)
          
          x = (/x1, x2, x3/)
          y = (/y1, y2, y3/)

          dphi1 = 0
          do i = 1, 3
            do j = 1, 3
              do k = 1, 3
                do l = 1, 3
                  dphi1 = d1/r1**7*
     1              (x(i)-y(i))*(x(j)-y(j))*(x(k)-y(k))*(x(l)-y(l))
                  if ( i == j ) then
                    dphi1 = dphi1 + d2/r1**5*(x(k)-y(k))*(x(l)-y(l))
                    if ( k == l ) then
                        dphi1 = dphi1 + c2/r1**3
                    end if
                  end if
                  if ( i == k ) then
                    dphi1 = dphi1 - d2/r1**5*(x(j)-y(j))*(x(l)-y(l))
                    if ( j == l ) then
                        dphi1 = dphi1 - c2/r1**3
                    end if
                  end if
                  if ( j == k ) then
                    dphi1 = dphi1 - d2/r1**5*(x(i)-y(i))*(x(l)-y(l))
                    if ( i == l ) then
                        dphi1 = dphi1 - c2/r1**3
                    end if
                  end if
                  if ( i == l ) then
                    dphi1 = dphi1 + c1/r1**5*(x(j)-y(j))*(x(k)-y(k))
                  end if
                  if ( j == l ) then
                    dphi1 = dphi1 + c1/r1**5*(x(i)-y(i))*(x(k)-y(k))
                  end if
                  if ( k == l ) then
                    dphi1 = dphi1 + c1/r1**5*(x(i)-y(i))*(x(j)-y(j))
                  end if

                  dphi(i,j,k,l) = -dphi1

                end do
              end do
            end do
          end do

          deallocate(x)
          deallocate(y)
          end subroutine

          ! predicates.cxx
          real function orient3dfast(pa, pb, pc, pd) result(vtet)
          real, dimension(:), intent(in) :: pa, pb, pc, pd
          real :: v

          real :: adx, bdx, cdx
          real :: ady, bdy, cdy
          real :: adz, bdz, cdz

          adx = pa(1) - pd(1)
          bdx = pb(1) - pd(1)
          cdx = pc(1) - pd(1)

          ady = pa(2) - pd(2)
          bdy = pb(2) - pd(2)
          cdy = pc(2) - pd(2)

          adz = pa(3) - pd(3)
          bdz = pb(3) - pd(3)
          cdz = pc(3) - pd(3)

          v = adx * (bdy * cdz - bdz * cdy)
     1      + bdx * (cdy * adz - cdz * ady)
     2      + cdx * (ady * bdz - adz * bdy)

          vtet = v/6.0 

          end function

          real function triangle_area(a, b, c) result(atri)
          real, intent(in) :: a, b, c
          real :: p, s

          p = (a + b + c)/2.0
          s = sqrt(p*(p-a)*(p-b)*(p-c))

          atri = s
          end function

          ! axb = c
          subroutine cross_productReal(c, b, a)
          real, dimension(:), intent(in) :: b, a
          real, dimension(:), intent(out) :: c

          c(1) = a(2)*b(3) - a(3)*b(2)
          c(2) = a(3)*b(1) - a(1)*b(3)
          c(3) = a(1)*b(2) - a(2)*b(1)

          end subroutine

          ! tetgen.cxx
          real function triarea(pa, pb, pc) result(area)
          real, dimension(:), intent(in) :: pa, pb, pc
          real, dimension(:,:), allocatable :: a
          allocate(a(3,3))
          a = 0

          a(1,1) = pb(1) - pa(1)
          a(2,1) = pb(2) - pa(2)
          a(3,1) = pb(3) - pa(3)

          a(1,2) = pc(1) - pa(1)
          a(2,2) = pc(2) - pa(2)
          a(3,2) = pc(3) - pa(3)

          call cross_product(a(:,3), a(:,2), a(:,1))
          
          area = 0.5*sqrt(dot_product(a(:,3),a(:,3)))

          deallocate(a) 
          end function

          ! tetgen.cxx
          subroutine facenormal(pa, pb, pc, nor, pivot) 
          real, dimension(:), intent(in) :: pa, pb, pc
          integer, intent(in) :: pivot
          real, dimension(:), intent(out) :: nor

          real, dimension(:), allocatable :: v1, v2, v3
          real, dimension(:), allocatable :: pv1, pv2
          real :: l1, l2, l3

          allocate(v1(3),v2(3),v3(3))
          v1 = 0
          v2 = 0
          v3 = 0
          l1 = 0
          l2 = 0
          l3 = 0

          v1(1) = pb(1) - pa(1)
          v1(2) = pb(2) - pa(2)
          v1(3) = pb(3) - pa(3)

          v2(1) = pa(1) - pc(1)
          v2(2) = pa(2) - pc(2)
          v2(3) = pa(3) - pc(3)

          if ( pivot >= 0 ) then
              v3(1) = pc(1) - pb(1)
              v3(2) = pc(2) - pb(2)
              v3(3) = pc(3) - pb(3)

              l1 = dot_product(v1, v1)
              l2 = dot_product(v2, v2)
              l3 = dot_product(v3, v3)
              if ( l1 < l2 ) then
                  if ( l2 < l3 ) then
                      pv1 = v1
                      pv2 = v2
                  else
                      pv1 = v3
                      pv2 = v1

                  end if
              else
                  if ( l1 < l3 ) then
                      pv1 = v1
                      pv2 = v2
                  else
                      pv1 = v2
                      pv2 = v3
                  end if
              end if
          else 
              pv1 = v1
              pv2 = v2
          end if

          call cross_product(nor, pv1, pv2)

          nor = nor/sqrt(dot_product(nor, nor))

          deallocate(v1, v2, v3)
          end subroutine

          subroutine body_force(eps, A, lam, mu, 
     1       x1, x2, x3,
     2       pa, pb, pc, pd, 
     3       pgamma, rs, iflags, rc)
          real, intent(in) :: A, lam, mu
          real, intent(in) :: x1, x2, x3
          real, dimension(:), intent(in) :: pa, pb, pc, pd
          real, dimension(:), intent(in) :: pgamma
          real, dimension(:,:,:), intent(in) :: rs
          integer, dimension(:), intent(in) :: iflags
          real, dimension(:,:), intent(in) :: rc
          real, dimension(:,:), intent(out) :: eps
          
          real, dimension(:,:), allocatable :: p
          real, dimension(:,:,:), allocatable :: phi
          real, dimension(:,:,:,:), allocatable :: dphi

          real :: y1, y2, y3
          real, dimension(:), allocatable :: nor, pg
          real :: tria, vtet
          integer :: pivot, i, j
          real, dimension(:,:), allocatable :: eps1, eps2, eps3, eps4
          real, dimension(:,:), allocatable :: eps0

          allocate(p(3,3))
          allocate(phi(3,3,3))
          allocate(dphi(3,3,3,3))
          allocate(nor(3), pg(3))
          allocate(eps1(3,3))
          allocate(eps2(3,3))
          allocate(eps3(3,3))
          allocate(eps4(3,3))
          allocate(eps0(3,3))

          p = 0
          phi = 0
          dphi = 0
          nor = 0
          tria = 0
          vtet = 0
          pivot = 1

          eps1 = 0
          eps2 = 0
          eps3 = 0
          eps4 = 0
          eps0 = 0

          y1 = 0
          y2 = 0
          y3 = 0
          
          y1 = (pa(1) + pb(1) + pc(1) + pd(1))/4.0
          y2 = (pa(2) + pb(2) + pc(2) + pd(2))/4.0
          y3 = (pa(3) + pb(3) + pc(3) + pd(3))/4.0
          
          call pstrain(p, pgamma, rs, iflags, rc)
          
          call kelvin_potential(phi, A, lam, mu, 
     1       x1, x2, x3,
     2       y1, y2, y3)

          call dkelvin_potential(dphi, A, lam, mu, 
     1       x1, x2, x3,
     2       y1, y2, y3)
          ! 1->2->3
          call facenormal(pa, pb, pc, nor, pivot) 
          tria = triarea(pa, pb, pc) 
          pg = matmul(p, nor)*tria
          do i = 1, 3
            do j = 1, 3
            eps1(i,j) = dot_product(phi(i,j,1:3), pg)
            end do
          end do
          
          ! 4->3->2
          call facenormal(pd, pc, pb, nor, pivot)
          tria = triarea(pd, pc, pb)
          pg = matmul(p, nor)*tria
          do i = 1, 3
            do j = 1, 3
              eps2(i,j) = dot_product(phi(i,j,1:3), pg)
            end do
          end do
          
          ! 1->3->4
          call facenormal(pa, pc, pd, nor, pivot)
          tria = triarea(pa, pc, pd)
          pg = matmul(p, nor)*tria
          do i = 1, 3
            do j = 1, 3
              eps3(i,j) = dot_product(phi(i,j,1:3), pg)
            end do
          end do
          
          ! 4->2->1
          call facenormal(pd, pb, pa, nor, pivot)
          tria = triarea(pd, pb, pa)
          pg = matmul(p, nor)*tria
          do i = 1, 3
            do j = 1, 3
              eps4(i,j) = dot_product(phi(i,j,1:3), pg)
            end do
          end do

          vtet = orient3dfast(pa, pb, pc, pd)    
          eps0 = 0
          do i = 1, 3
            do j = 1, 3
              eps0(i,j) = eps0(i,j) + sum(dphi(i,j,1:3,1:3)*p)
            end do
          end do
          eps0 = eps0*vtet
          
          eps = eps1+eps2+eps3+eps4-eps0
          
          deallocate(p)
          deallocate(phi)
          deallocate(dphi)
          deallocate(nor,pg)
          deallocate(eps1,eps2,eps3,eps4,eps0)

          end subroutine

          integer function slipsys_enum(nslipdir, nslipnor) result(n)
          integer, intent(in) :: nslipdir, nslipnor
          integer :: nor 
          if ((nslipnor .lt. 1) .or. (nslipnor .gt. 4)) then
              print *, " maximum number of normals (4) ", nslipnor
              call abort
          end if
          if ((nslipdir .lt. 1) .or. (nslipdir .gt. 3)) then
              print *, " maixmum number of directions (3) ", nslipdir
              call abort
          end if
          nor = 3*(nslipnor-1)+(nslipdir-1)+1
          n = nor
          end function

          subroutine sort4plane(pshear, rs, iflags, rc, sigma, 
     1                            nrc, nc) 
          integer, intent(in) :: nc
          integer, intent(in) :: nrc
          real, dimension(:,:), intent(in) :: sigma
          integer, dimension(:), intent(out) :: iflags
          real, dimension(:,:,:), intent(out) :: rs
          real, dimension(:,:), intent(out) :: rc
          real, dimension(:), intent(out) :: pshear

          integer :: n, d
          integer :: en, ed
          integer :: aen, aed
          real :: v
          
          
          real, dimension(:,:), allocatable :: lsig
          real, dimension(:), allocatable :: lp
          real, dimension(:), allocatable :: pgamma
          integer, dimension(:), allocatable :: sd, sn
          !real, dimension(:), allocatable :: pshear
          integer :: i, j, k, l, m
          
          allocate(lsig(3,3)) 
          allocate(lp(3))
          !allocate(pshear(4))
          allocate(sd(3), sn(3))

           
          lsig = 0
          lp = 0
          pshear = 0
          sd = 0
          sn = 0

          en = size(ele, dim=1)
          ed = size(ele, dim=2)
          n = size(node, dim=1)
          d = size(node, dim=2)
          aen = size(aele, dim=1)
          aed = size(aele, dim=2)

          if (nc > en) then
              print *, " number of elements: ", en
              call abort
          end if

          if (nrc > aen) then
              print *, " number of orientations: ", aen
              call abort
          end if

          allocate(pgamma(4)) 
          pgamma = 0
          iflags = -1
          rs = 0
          rc = 0
    
          call euler_zxz(rc, aele(nrc,1), aele(nrc,2), aele(nrc,3))
          lsig = matmul(rc, matmul(sigma, transpose(rc)))
          do j = 1, 4 
            lp = matmul(lsig, fcc_slipnor(1:3, j))
            if (j == 1) then
                do l = 1, 3
                  pgamma(l) = dot_product(fcc_primary(1:3, l), lp)
                  rs(j, l, 1) = pgamma(l) 
                end do
                pgamma(4) = maxval(abs(pgamma(1:3)))
                do l = 1, 3
                  if( abs(pgamma(l)) == pgamma(4) ) then
                      i = slipsys_enum(l, j)
                      iflags(i) = 1
                      exit
                  end if
                end do
                pshear(j) = pgamma(4) 
            end if
            if (j == 2) then
                do l = 1, 3
                  pgamma(l) = dot_product(fcc_conjugate(1:3, l), lp)
                  rs(j, l, 1) = pgamma(l) 
                end do
                pgamma(4) = maxval(abs(pgamma(1:3)))
                do l = 1, 3
                  if( abs(pgamma(l)) == pgamma(4) ) then
                      i = slipsys_enum(l, j)
                      iflags(i) = 1
                      exit
                  end if
                end do
                pshear(j) = pgamma(4)
            end if
            if (j == 3) then
                do l = 1, 3
                  pgamma(l) = dot_product(fcc_crossglide(1:3, l), lp)
                  rs(j, l, 1) = pgamma(l) 
                end do
                pgamma(4) = maxval(abs(pgamma(1:3)))
                do l = 1, 3
                  if( abs(pgamma(l)) == pgamma(4) ) then
                      i = slipsys_enum(l, j)
                      iflags(i) = 1
                      exit
                  end if
                end do
                pshear(j) = pgamma(4)
            end if
            if (j == 4) then
                do l = 1, 3
                  pgamma(l) = dot_product(fcc_critical(1:3, l), lp)
                  rs(j, l, 1) = pgamma(l) 
                end do
                pgamma(4) = maxval(abs(pgamma(1:3)))
                do l = 1, 3
                  if( abs(pgamma(l)) == pgamma(4) ) then
                      i = slipsys_enum(l, j)
                      iflags(i) = 1
                      exit
                  end if
                end do
                pshear(j) = pgamma(4)
            end if
          end do
 
          ! v = minval(pshear)
          v = pshear(1)
          k = 1
          do j = 1, 4
              if ( pshear(j) .lt. v ) then
                  k = j
                  v = pshear(j)
              end if
          end do
          ! index of min
          ! j = minloc(pshear)
          do l = 1, 3
            i = slipsys_enum(l, k)
            ! print *, i
            iflags(i) = -1
          end do
          ! v = maxval(pshear)
          v = pshear(1)
          m = 1
          do j = 1, 4
            if ( pshear(j) .gt. v) then
                m = j
                v = pshear(j)
            end if
          end do
          ! index of max
          ! j = maxloc(pshear)
          v = minval(rs(m, 1:3, 1))
          do l = 1, 3
            i = slipsys_enum(l, m)
            if ( rs(m, l, 1) .eq. v ) then
              iflags(i) = -1
            else
              iflags(i) = 1
            end if
          end do
          !write(6, *) iflags
          !write(6, '(4F10.4)') pshear
          rs = 0
          m = 0
          do j = 1, 4
            sn = fcc_slipnor(1:3, j)
            do l = 1, 3
              i = slipsys_enum(l, j)
              if ( iflags(i) > 0 ) then
                m = m + 1                 
                if ( j .eq. 1) then
                  sd = fcc_primary(1:3, l)
                end if
                if ( j .eq. 2) then
                  sd = fcc_conjugate(1:3, l)
                end if
                if ( j .eq. 3) then
                  sd = fcc_crossglide(1:3, l)
                end if
                if ( j .eq. 4) then
                  sd = fcc_critical(1:3, l)
                end if
                call slip_system(rs(m, 1:3, 1:3), sd, sn)
              end if
            end do
          end do

          deallocate(lsig)
          deallocate(lp, sd, sn)
          deallocate(pgamma)

          end subroutine
      
          subroutine gshear(tshear, rshear, pgamma, sigma, lam, mu, 
     1                      nrc, nc)
          real, intent(in) :: lam, mu
          integer, intent(in) :: nrc, nc
          real, dimension(:,:), intent(in) :: sigma
          real, dimension(:), intent(in) :: pgamma

          real, dimension(:,:), intent(out) :: tshear 
          real, dimension(:,:), intent(out) :: rshear

          real, parameter :: pi = 4.0*atan(1.0)
          
          integer :: n, d
          integer :: en, ed
          integer :: aen, aed

          real :: kA
          real :: x1, x2, x3 
          real, dimension(:), allocatable :: pa, pb, pc, pd
          real, dimension(:,:), allocatable :: eps
          real, dimension(:,:), allocatable :: eps1

          integer, dimension(:), allocatable :: iflags
          real, dimension(:,:,:), allocatable :: rs
          real, dimension(:,:), allocatable :: rc
          real, dimension(:), allocatable :: pshear

          !integer :: nrc
          integer :: i, j, k, l, m
          allocate(iflags(12))
          allocate(rs(4, 3, 3))
          allocate(pa(3), pb(3), pc(3), pd(3))
          allocate(eps(3,3))
          allocate(rc(3,3))
          allocate(eps1(3,3)) 
          allocate(pshear(4))

          pa = 0
          pb = 0
          pc = 0
          pd = 0
          eps = 0
          kA = 0
          kA = (lam+mu)/8.0/pi/mu/(lam+2.0*mu)
          iflags = -1
          rs = 0
          rc = 0
          eps1 = 0

          en = size(ele, dim=1)
          ed = size(ele, dim=2)
          n = size(node, dim=1)
          d = size(node, dim=2)
          aen = size(aele, dim=1)
          aed = size(aele, dim=2)

          
          do j = 1, n
            if (j .eq. ele(nc, 1)) then
              pa = node(j, 1:3)
            end if
            if (j .eq. ele(nc, 2)) then
              pb = node(j, 1:3)
            end if
            if (j .eq. ele(nc, 3)) then
              pc = node(j, 1:3)
            end if
            if (j .eq. ele(nc, 4)) then
              pd = node(j, 1:3) 
            end if
          end do
          x1 = 0.25*(pa(1) + pb(1) + pc(1) + pd(1))
          x2 = 0.25*(pa(2) + pb(2) + pc(2) + pd(2))
          x3 = 0.25*(pa(3) + pb(3) + pc(3) + pd(3))

          tshear = 0
          do k = 1, en
            if (k .ne. nc) then
              do l = 1, n
                if (l .eq. ele(k, 1)) then
                  pa = node(l, 1:3)
                end if
                if (l .eq. ele(k, 2)) then
                  pb = node(l, 1:3)
                end if
                if (l .eq. ele(k, 3)) then
                  pc = node(l, 1:3)
                end if
                if (l .eq. ele(k, 4)) then
                  pd = node(l, 1:3)
                end if
              end do

              ! 
              m = mod(k, aen)

              call sort4plane(pshear, rs, iflags, rc, sigma,
     1                        m, k) 

              call body_force(eps, kA, lam, mu, 
     1                        x1, x2, x3,
     2                        pa, pb, pc, pd, 
     3                        pgamma(4*k-3:4*k), rs, iflags, rc)
              tshear = tshear + eps
            end if
          end do

          call sort4plane(pshear, rs, iflags, rc, sigma,
     1                    nrc, nc)
          call pstrain(eps1, pgamma(4*nc-3:4*nc), 
     1                 rs, iflags, rc)

          rshear = -tshear - 2*mu*eps1
      
          deallocate(iflags)
          deallocate(rs)
          deallocate(eps)
          deallocate(pa, pb, pc, pd)
          deallocate(rc) 
          deallocate(eps1)
          deallocate(pshear)
          end subroutine


      end module
