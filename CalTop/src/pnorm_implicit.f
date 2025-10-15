

      subroutine pnorm_implicit_c3d4(co,kon,ipkon,lakon,ne,mi,
     &     xstiff, v, lam, design, penal,
     &     nea, neb, list, ilist, djdrho)

c--- Implicit term for dJ/drho_e using ce = penal * rho^(penal-1)
c    dJ/drho|imp = - ce * (lambda_e^T * K0_e * u_e)
c    C3D4 only; K0_e u_e is assembled via B^T D B u * (xsj*weight)
c
c    Inputs:
c      co(3,nk)            : nodal coordinates
c      kon(*)              : connectivity vector
c      ipkon(ne)           : element pointer into kon
c      lakon(ne)           : element type string; selects C3D4
c      ne                  : number of elements
c      mi(*)               : CalculiX size array (mi(1)=#GP max, mi(2)=dofs)
c      xstiff(27,mi(1),ne) : material D in CalculiX packed layout
c      v(0:mi(2),nk)       : primal nodal field (displacements)
c      lam(0:mi(2),nk)     : adjoint nodal field (same layout as v)
c      design(ne)          : element densities (clamped to [0,1] here)
c      penal               : SIMP exponent
c      nea,neb             : element range (can pass 1,ne)
c      list, ilist(*)      : optional selection; if list=1 use ilist(k)
c
c    Output:
c      djdrho_impl(ne)          : accumulates the implicit contribution

      implicit none

c--- mesh
      integer ne, kon(*), ipkon(*), mi(*), nea, neb, list, ilist(*)
      character*8 lakon(*), lakonl
      integer i,k,j,idx,nope,iflag,jj
      integer konl(4)

c--- geometry
      real*8 co(3,*)

c--- fields
      real*8 v(0:mi(2),*), lam(0:mi(2),*)

c--- material (D) store
      real*8 xstiff(27,mi(1),*)

c--- design / params
      real*8 design(*), penal, heav

c--- output
      real*8 djdrho(*)

c--- locals
      real*8 xl(3,4), shp(4,4), xsj, xi, et, ze, weight
      real*8 B(6,12), ue(12), le(12), k0u(12), dotlam
      real*8 eps(6), sig(6), rho, ce, rho_eff
      integer a, m, n, m1

c--- Gauss rule (C3D4, 1 point)
      include 'gauss.f'

c--- init
      iflag = 3
      nope  = 4

c--- element loop
      do k = nea, neb
         if (list.eq.1) then
            i = ilist(k)
         else
            i = k
         endif

         if (ipkon(i).lt.0) cycle
         lakonl = lakon(i)
         if (lakonl(4:4).ne.'4') cycle  ! C3D4 only

c------ gather connectivity & coords
         idx = ipkon(i)
         do j=1,4
            konl(j) = kon(idx+j)
            xl(1,j) = co(1,konl(j))
            xl(2,j) = co(2,konl(j))
            xl(3,j) = co(3,konl(j))
         enddo

c------ local u and lambda (3 dofs per node, stacked)
         m = 0
         do j=1,4
            ue(m+1) = v(1,konl(j))
            ue(m+2) = v(2,konl(j))
            ue(m+3) = v(3,konl(j))
            le(m+1) = lam(1,konl(j))
            le(m+2) = lam(2,konl(j))
            le(m+3) = lam(3,konl(j))
            m = m + 3
         enddo

c------ single Gauss point for C3D4
         jj     = 1
         xi     = gauss3d4(1,jj)
         et     = gauss3d4(2,jj)
         ze     = gauss3d4(3,jj)
         weight = weight3d4(jj)

c------ shape & jacobian
         call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
         if (xsj.le.0.d0) cycle

c------ build B (6x12) from global derivatives dN/dx = shp(1,*), etc.
c       Voigt order: [exx eyy ezz exy exz eyz]
         do a=1,12
            k0u(a) = 0.d0
         enddo

         call build_B_c3d4(shp,B)

c------ strain eps = B * u
         do m=1,6
            eps(m) = 0.d0
            do n=1,12
               eps(m) = eps(m) + B(m,n)*ue(n)
            enddo
         enddo

c------ stress sig = D * eps  (use xstiff mapping like in your RHS code)
         call mult_D_vec(sig, eps, xstiff(1,jj,i))

c------ internal nodal forces k0u += B^T * sig * vol
         do n=1,12
            do m=1,6
               k0u(n) = k0u(n) + B(m,n)*sig(m)
            enddo
            k0u(n) = k0u(n) * (xsj*weight)
         enddo

c------ lambda^T * (K0 u)
         dotlam = 0.d0
         do n=1,12
            dotlam = dotlam + le(n)*k0u(n)
         enddo

c------ ce = penal * rho^(penal-1)   (same logic as e_c3d_se.f)
         rho = design(i)
         if (rho .lt. 0.d0) rho = 0.d0
         if (rho .gt. 1.d0 ) rho = 1.d0
         rho_eff = dmax1(rho, 1e-06)

         if (rho .gt. 1e-06) then
            heav = 1.d0
         else
            heav = 0.d0
         endif

         !if (rho.le.0.d0) then
         !   ce = 0.d0
         !else
         !   ce = penal * rho**(penal-1.d0)
         !endif

         ce = penal * rho_eff**(penal-1.d0) * heav

c------ accumulate implicit sensitivity
         djdrho(i) = djdrho(i) - ce * dotlam

      enddo

      return
      end


c======================================================================
c  Build B for a 4-node tet from CalculiX shape derivatives
c  shp(1,j)=dNj/dx, shp(2,j)=dNj/dy, shp(3,j)=dNj/dz
c======================================================================
      subroutine build_B_c3d4(shp,B)
      implicit none
      real*8 shp(4,4), B(6,12)
      integer j, col

      do j=1,12
         B(1,j)=0.d0; B(2,j)=0.d0; B(3,j)=0.d0
         B(4,j)=0.d0; B(5,j)=0.d0; B(6,j)=0.d0
      enddo

      do j=1,4
         col = 3*(j-1)
c        exx
         B(1,col+1) = shp(1,j)
c        eyy
         B(2,col+2) = shp(2,j)
c        ezz
         B(3,col+3) = shp(3,j)
c        exy
         B(4,col+1) = shp(2,j)
         B(4,col+2) = shp(1,j)
c        exz
         B(5,col+1) = shp(3,j)
         B(5,col+3) = shp(1,j)
c        eyz
         B(6,col+2) = shp(3,j)
         B(6,col+3) = shp(2,j)
      enddo

      return
      end


c======================================================================
c  sig = D * eps using CalculiX xstiff packed layout (same as your RHS)
c  xD(1..27) holds upper-triangular 6x6 mapping in CCX’s order.
c======================================================================
      subroutine mult_D_vec(sig,eps,xD)
      implicit none
      real*8 sig(6), eps(6), xD(27)

c  Unrolled like in your RHS (ptv computation), matching CCX layout
      sig(1)= xD( 1)*eps(1) + xD( 2)*eps(2) + xD( 4)*eps(3)
     &      + xD( 7)*eps(4) + xD(11)*eps(5) + xD(16)*eps(6)

      sig(2)= xD( 2)*eps(1) + xD( 3)*eps(2) + xD( 5)*eps(3)
     &      + xD( 8)*eps(4) + xD(12)*eps(5) + xD(17)*eps(6)

      sig(3)= xD( 4)*eps(1) + xD( 5)*eps(2) + xD( 6)*eps(3)
     &      + xD( 9)*eps(4) + xD(13)*eps(5) + xD(18)*eps(6)

      sig(4)= xD( 7)*eps(1) + xD( 8)*eps(2) + xD( 9)*eps(3)
     &      + xD(10)*eps(4) + xD(14)*eps(5) + xD(19)*eps(6)

      sig(5)= xD(11)*eps(1) + xD(12)*eps(2) + xD(13)*eps(3)
     &      + xD(14)*eps(4) + xD(15)*eps(5) + xD(20)*eps(6)

      sig(6)= xD(16)*eps(1) + xD(17)*eps(2) + xD(18)*eps(3)
     &      + xD(19)*eps(4) + xD(20)*eps(5) + xD(21)*eps(6)

      return
      end
