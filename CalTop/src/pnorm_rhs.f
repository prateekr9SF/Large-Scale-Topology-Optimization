      subroutine pnorm_rhs(co,kon,ipkon,lakon,ne,
     &     stx,xstiff,mi,rhs,alpha,p,design,penal,
     &     sig0,eps_relax,rho_min,
     &     nea,neb,list,ilist)

c  Adjoint RHS for epsilon-relaxed stress p-norm, C3D4 only
c  J = ( sum w * phi^p )^(1/p),  phi = vm/(rho^p*sig0) + eps - eps/rho_eff
c  dJ/dsigma = alpha * phi^(p-1) * (1/(rho^p*sig0)) * dvm/dsigma
c  RHS = ∫ B^T [ D * (dJ/dsigma) ] dV  (assemble in Voigt: ptv = D*g)

      implicit none

c--- mesh / topology
      integer ne,kon(*),ipkon(*),mi(*),nea,neb,list
      integer ilist(*)
      character*8 lakon(*),lakonl

c--- geometry
      real*8 co(3,*)

c--- stresses and material stiffness at IPs
      real*8 stx(6,mi(1),*), xstiff(27,mi(1),*)

c--- output RHS (nodal forces: first index 0..mi(2))
      real*8 rhs(0:mi(2),*)

c--- parameters
      real*8 alpha,p,design(*),penal,sig0,eps_relax,rho_min

c--- locals
      integer i,k,jj,indexe,nope,iflag
      integer m,j
      integer konl(4)
      real*8 xl(3,4),shp(4,4),xsj,xi,et,ze,weight
      real*8 sx,sy,sz,txy,txz,tyz,vm2,vm,invvm
      real*8 rho,rho_eff,rho_p,phi,fac
      real*8 g(6),ptv(6)

c--- NEW (Voigt assembly)
      real*8 B(6,12), rhs_loc(12)

c--- Gauss rules
      include 'gauss.f'

c--- init
      iflag  = 3
      nope   = 4

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

c-------- gather nodal connectivity and coords
         indexe = ipkon(i)
         do j=1,4
            konl(j) = kon(indexe+j)
            xl(1,j) = co(1,konl(j))
            xl(2,j) = co(2,konl(j))
            xl(3,j) = co(3,konl(j))
         enddo

c-------- element density (filtered)
         rho = design(i)
         if (rho .lt. 0.d0) rho = 0.d0
         if (rho .gt. 1.d0) rho = 1.d0
         rho_eff = dmax1(rho, rho_min)
         rho_p   = rho_eff**penal

c-------- single Gauss point for C3D4
         jj     = 1
         xi     = gauss3d4(1,jj)
         et     = gauss3d4(2,jj)
         ze     = gauss3d4(3,jj)
         weight = weight3d4(jj)

c-------- shapes & jacobian
         call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
         if (xsj.le.0.d0) cycle

c-------- stresses at this IP (Cauchy in Voigt order)
         sx  = stx(1,jj,i)
         sy  = stx(2,jj,i)
         sz  = stx(3,jj,i)
         txy = stx(4,jj,i)
         txz = stx(5,jj,i)
         tyz = stx(6,jj,i)

c-------- von Mises (consistent with CCX)
         vm2 = (sx-sy)*(sx-sy) + (sy-sz)*(sy-sz) + (sz-sx)*(sz-sx)
         vm2 = 0.5d0*vm2 + 3.d0*(txy*txy + txz*txz + tyz*tyz)
         if (vm2 .le. 0.d0) cycle
         vm  = dsqrt(vm2)
         invvm = 1.d0/vm

c-------- epsilon-relaxed phi
         phi = vm/(rho_p*sig0) + eps_relax - eps_relax/rho_eff
         if (phi .le. 0.d0) cycle

c-------- dJ/dsigma factor (no extra p — alpha = J^(1-p))
         fac = alpha * (phi**(p-1)) / (rho_p*sig0)

c-------- dvm/dsigma in Voigt, scale by fac
         g(1) = fac * ((2.d0*sx - sy - sz) * 0.5d0 * invvm)
         g(2) = fac * ((2.d0*sy - sx - sz) * 0.5d0 * invvm)
         g(3) = fac * ((2.d0*sz - sx - sy) * 0.5d0 * invvm)
         g(4) = fac * (3.d0*txy * invvm)
         g(5) = fac * (3.d0*txz * invvm)
         g(6) = fac * (3.d0*tyz * invvm)


c-------- ptv = D * g  (same mapping you already use)
         ptv(1)= g(1)*xstiff( 1,jj,i) + g(2)*xstiff( 2,jj,i)
     &          + g(3)*xstiff( 4,jj,i) + g(4)*xstiff( 7,jj,i)
     &          + g(5)*xstiff(11,jj,i) + g(6)*xstiff(16,jj,i)
         ptv(2)= g(1)*xstiff( 2,jj,i) + g(2)*xstiff( 3,jj,i)
     &          + g(3)*xstiff( 5,jj,i) + g(4)*xstiff( 8,jj,i)
     &          + g(5)*xstiff(12,jj,i) + g(6)*xstiff(17,jj,i)
         ptv(3)= g(1)*xstiff( 4,jj,i) + g(2)*xstiff( 5,jj,i)
     &          + g(3)*xstiff( 6,jj,i) + g(4)*xstiff( 9,jj,i)
     &          + g(5)*xstiff(13,jj,i) + g(6)*xstiff(18,jj,i)
         ptv(4)= g(1)*xstiff( 7,jj,i) + g(2)*xstiff( 8,jj,i)
     &          + g(3)*xstiff( 9,jj,i) + g(4)*xstiff(10,jj,i)
     &          + g(5)*xstiff(14,jj,i) + g(6)*xstiff(19,jj,i)
         ptv(5)= g(1)*xstiff(11,jj,i) + g(2)*xstiff(12,jj,i)
     &          + g(3)*xstiff(13,jj,i) + g(4)*xstiff(14,jj,i)
     &          + g(5)*xstiff(15,jj,i) + g(6)*xstiff(20,jj,i)
         ptv(6)= g(1)*xstiff(16,jj,i) + g(2)*xstiff(17,jj,i)
     &          + g(3)*xstiff(18,jj,i) + g(4)*xstiff(19,jj,i)
     &          + g(5)*xstiff(20,jj,i) + g(6)*xstiff(21,jj,i)

c-------- NEW: build B and assemble rhs_loc = (xsj*weight) * B^T * ptv
         call build_B_c3d42(shp,B)

         do m=1,12
            rhs_loc(m) = 0.d0
         enddo

         do m=1,12
            do j=1,6
               rhs_loc(m) = rhs_loc(m) + B(j,m) * ptv(j)
            enddo
            rhs_loc(m) = rhs_loc(m) * (xsj * weight)
         enddo

c-------- scatter rhs_loc(1..12) into nodal rhs(1..3,konl(1..4))
         do j=1,4
            rhs(1,konl(j)) = rhs(1,konl(j)) + rhs_loc(3*(j-1)+1)
            rhs(2,konl(j)) = rhs(2,konl(j)) + rhs_loc(3*(j-1)+2)
            rhs(3,konl(j)) = rhs(3,konl(j)) + rhs_loc(3*(j-1)+3)
         enddo

      enddo  ! elements

      return
      end

c======================================================================
c  Build B for a 4-node tet from CalculiX shape derivatives
c  shp(1,j)=dNj/dx, shp(2,j)=dNj/dy, shp(3,j)=dNj/dz
c  >>> Engineering shear convention (exy=γxy, etc.) <<<
c======================================================================
      subroutine build_B_c3d42(shp,B)
      implicit none
      real*8 shp(4,4), B(6,12)
      integer j, col

      do j=1,12
         B(1,j)=0.d0; B(2,j)=0.d0; B(3,j)=0.d0
         B(4,j)=0.d0; B(5,j)=0.d0; B(6,j)=0.d0
      enddo

      do j=1,4
         col = 3*(j-1)
c        normal strains
         B(1,col+1) = shp(1,j)            ! exx
         B(2,col+2) = shp(2,j)            ! eyy
         B(3,col+3) = shp(3,j)            ! ezz
c        shear (engineering) strains
         B(4,col+1) = shp(2,j)            ! exy = dNy/dx + dNx/dy
         B(4,col+2) = shp(1,j)
         B(5,col+1) = shp(3,j)            ! exz
         B(5,col+3) = shp(1,j)
         B(6,col+2) = shp(3,j)            ! eyz
         B(6,col+3) = shp(2,j)
      enddo

      return
      end

