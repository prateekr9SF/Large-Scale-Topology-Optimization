

      subroutine pnorm_implicit(co,kon,ipkon,lakon,ne,mi,
     &     xstiff, v, lam, design, penal, pexp, relax, sig0,
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
c      pexp                : the exponential term for Pnorm calculations
c      relax               : stress relaxation factor
c      sig0                : minimum stress allowable
c      nea,neb             : element range (can pass 1,ne)
c      list, ilist(*)      : optional selection; if list=1 use ilist(k)

c
c    Output:
c      djdrho_impl(ne)          : accumulates the implicit contribution

      implicit none

c--- mesh
      integer ne, kon(*), ipkon(*), mi(*), nea, neb, list, ilist(*)
      character*8 lakon(*), lakonl
      integer i,k,j,idx,nope,mint3d, iflag,jj
      integer konl(26)

c--- geometry
      real*8 co(3,*)

c--- fields
      real*8 v(0:mi(2),*), lam(0:mi(2),*)

c--- material (D) store
      real*8 xstiff(27,mi(1),*)
!     THIS IS NOT THE CONST. MATRIX! 
!     xstiff(1,jj,i) is E and xstiff(2,jj,i) is nu

c--- design / params
      real*8 design(*), penal, heav, pexp, relax
      real*8 sig0

c--- output
      real*8 djdrho(*)

c--- locals
      real*8 xl(3,26), shp(4,26), xsj, xi, et, ze, weight
      real*8 B(6,30), ue(30), le(30), k0u(30), dotlam
      real*8 eps(6), sig(6), rho, ce, rho_eff,sigP(6)
      integer a, m, n, m1
      real*8 vm,vm2,sige, vm_sum
      real*8 Epen
      real*8 E
      real*8 expli1
      real*8 expli2
c--- Gauss rule (C3D4, 1 point)
      include 'gauss.f'

c--- init
      iflag = 3

c--- element loop
      do k = nea, neb
         if (list.eq.1) then
            i = ilist(k)
         else
            i = k
         endif

         if (ipkon(i).lt.0) cycle
         lakonl = lakon(i)
         if (lakonl(4:4).eq.'4') then
            nope = 4
            mint3d = 1
         elseif (lakonl(4:5).eq.'10') then
            nope = 10
            mint3d = 4
         else
            cycle
         endif

c------ gather connectivity & coords
         idx = ipkon(i)
         do j=1,nope
            konl(j) = kon(idx+j)
            xl(1,j) = co(1,konl(j))
            xl(2,j) = co(2,konl(j))
            xl(3,j) = co(3,konl(j))
         enddo

c------ local u and lambda (3 dofs per node, stacked)
         m = 0
         do j=1,nope
            ue(m+1) = v(1,konl(j))
            ue(m+2) = v(2,konl(j))
            ue(m+3) = v(3,konl(j))
            le(m+1) = lam(1,konl(j))
            le(m+2) = lam(2,konl(j))
            le(m+3) = lam(3,konl(j))
            m = m + 3
         enddo

c------ gaussian integration for C3D4&D10
         vm_sum = 0.d0
         do jj=1,mint3d
            if(lakonl(4:5).eq.'10') then
               xi=gauss3d5(1,jj); et=gauss3d5(2,jj); ze=gauss3d5(3,jj)
               weight=weight3d5(jj)
               call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(lakonl(4:4).eq.'4') then
               xi=gauss3d4(1,jj); et=gauss3d4(2,jj); ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
               call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            endif
            xsj = dabs(xsj)
            do a=1,3*nope
               k0u(a) = 0.d0
            enddo
            call build_B_tet(shp,B,nope)
c------ strain eps = B * u
            do m=1,6
               eps(m) = 0.d0
               do n=1,3*nope
                  eps(m) = eps(m) + B(m,n)*ue(n)
               enddo
            enddo
            rho = design(i)

            Epen = (rho**(penal-1))*xstiff(1,jj,i)*penal
            E    = (rho**(penal))*xstiff(1,jj,i)
c------ stress sig = D * eps  (use xstiff mapping like in your RHS code)

            call mult_D_vec(sig, eps, E,xstiff(2,jj,1))
            call mult_D_vec(sigP, eps, Epen,xstiff(2,jj,1))
c------ calculate effective vm stress
            vm2 = (sig(1)-sig(2))**2
            vm2 = vm2 + (sig(2)-sig(3))**2
            vm2 = vm2 + (sig(3)-sig(1))**2
            vm2 = 0.5d0*vm2
            vm2 = vm2 + 3.d0* (sig(4))**2 
            vm2 = vm2 + 3.d0* (sig(5))**2 
            vm2 = vm2 + 3.d0* (sig(6))**2 
            vm_sum = vm_sum+dsqrt(max(vm2,0.d0))
c------ internal nodal forces k0u += B^T * sig * vol
            do n=1,3*nope
               do m=1,6
                  k0u(n) = k0u(n) + B(m,n)*sigP(m)
               enddo
               k0u(n) = k0u(n) * (xsj*weight)
            enddo
         enddo  !-------END of GP loop
c------ element-average vm and sige
         vm = vm_sum / dble(mint3d)
         rho = design(i)
         if (rho .lt. 0.d0) rho = 0.d0
         if (rho .gt. 1.d0) rho = 1.d0
         rho_eff = dmax1(rho, 1e-06)
         sige = vm/sig0 + relax - relax/rho_eff
         if (sige .lt. 0.d0) sige = 0.d0

c------ lambda^T * (K0 u)
         dotlam = 0.d0
         do n=1,3*nope
            dotlam = dotlam + le(n)*k0u(n)
         enddo

c------ accumulate implicit sensitivity
         djdrho(i) = djdrho(i)-dotlam !Implicit
         expli1 = (sige**(pexp-1))*relax/(rho_eff**2.d0) 
         expli2 = (sige**(pexp-1))*vm*penal/(rho_eff*sig0)
         djdrho(i) = (djdrho(i)+expli1+expli2)
         !write(*,*),"IMP1",dotlam
         !now we have qbar*dkdrho*q
      enddo

      return
      end


c======================================================================
c  Build B for a 4-node tet from CalculiX shape derivatives
c  shp(1,j)=dNj/dx, shp(2,j)=dNj/dy, shp(3,j)=dNj/dz
c======================================================================
      subroutine build_B_tet(shp,B,nope)
      implicit none
      integer nope
      real*8 shp(4,26), B(6,30)
      integer j, col

      do j=1,3*nope
         B(1,j)=0.d0; B(2,j)=0.d0; B(3,j)=0.d0
         B(4,j)=0.d0; B(5,j)=0.d0; B(6,j)=0.d0
      enddo

      do j=1,nope
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
         B(6,col+1) = shp(3,j)
         B(6,col+3) = shp(1,j)
c        eyz
         B(5,col+2) = shp(3,j)
         B(5,col+3) = shp(2,j)
      enddo

      return
      end


c======================================================================
c  sig = D * eps using CalculiX xstiff packed layout (same as your RHS)
c  xD(1..27) holds upper-triangular 6x6 mapping in CCX’s order.
c======================================================================
      subroutine mult_D_vec(sig,eps,E,nu)
      implicit none
      real*8 sig(6), eps(6), E, nu,al,um
      real*8 xD(7)

      um = E/(2.d0*(1.d0+nu))                  ! G
      al = nu*E/((1.d0+nu)*(1.d0-2.d0*nu))     ! lambda
      xD( 1)=(al+2.d0*um)  ! C11
      xD( 2)= al           ! C12
      xD( 3)= (al+2.d0*um)  ! C22
      xD( 4)= al           ! C13
      xD( 5)= al           ! C23
      xD( 6)= (al+2.d0*um)  ! C33
      xD(7)= um           ! C44 (τ12/ε12_tensorial)

c  Unrolled like in your RHS (ptv computation), matching CCX layout
      sig(1)= xD( 1)*eps(1) + xD( 2)*eps(2) + xD( 4)*eps(3)

      sig(2)= xD( 2)*eps(1) + xD( 3)*eps(2) + xD( 5)*eps(3)

      sig(3)= xD( 4)*eps(1) + xD( 5)*eps(2) + xD( 6)*eps(3)

      sig(4)= xD(7)*eps(4) 

      sig(5)= xD(7)*eps(5)

      sig(6)= xD(7)*eps(6)

      return
      end
