      subroutine pnorm_explicit(co,kon,ipkon,lakon,ne,
     &     stx,mi,design,penal,sig0,eps_relax,rho_min,
     &     alpha,pexp,nea,neb,list,ilist,djdrho)

c  Explicit term for dJ/drho_e (epsilon-relaxed, unnormalized p-norm)
c
c  J = ( sum_g w_g * phi_g^p )^(1/p),  alpha = J^(1-p)
c  phi = vm/(rho_eff^q * sig0) + eps - eps/rho_eff
c  dphi/drho = H(rho-rho_min) * [- q*vm/(sig0*rho_eff^(q+1)) + eps/rho_eff^2]
c  dJ/drho|exp = sum_g alpha * phi^(p-1) * dphi/drho * w_g
c
c  C3D4 only (one Gauss point). stx holds Cauchy stresses in Voigt order.

      implicit none

c--- mesh / topology
      integer ne,kon(*),ipkon(*),mi(*),nea,neb,list
      integer ilist(*)
      character*8 lakon(*),lakonl

c--- geometry
      real*8 co(3,*)

c--- stresses at IPs
      real*8 stx(6,mi(1),*)

c--- design / params
      real*8 design(*),penal,sig0,eps_relax,rho_min
      real*8 alpha,pexp

c--- output: per-element explicit sensitivity
      real*8 djdrho(*)

c--- locals
      integer i,k,j,indexe,nope,mint3d,iflag,jj
      integer konl(4)
      real*8 xl(3,4),shp(4,4),xsj,xi,et,ze,weight
      real*8 rho,rho_eff,rho_p,phi,vm2,vm
      real*8 sx,sy,sz,txy,txz,tyz
      real*8 dphi,heav

c--- Gauss rules
      include 'gauss.f'

c--- init (C3D4)
      iflag  = 3
      nope   = 4
      mint3d = 1

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

c-------- connectivity & coords
         indexe = ipkon(i)
         do j=1,4
            konl(j) = kon(indexe+j)
         enddo
         do j=1,4
            xl(1,j) = co(1,konl(j))
            xl(2,j) = co(2,konl(j))
            xl(3,j) = co(3,konl(j))
         enddo

c-------- element density and effective density
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
         !print *, 'gp=', xi, et, ze, ' w=', weight

         call shape4tet(xi,et,ze,xl,xsj,shp,iflag)

c-------- stresses (Cauchy, Voigt order)
         sx  = stx(1,jj,i)
         sy  = stx(2,jj,i)
         sz  = stx(3,jj,i)
         txy = stx(4,jj,i)
         txz = stx(5,jj,i)
         tyz = stx(6,jj,i)

c-------- von Mises
         vm2 = (sx-sy)*(sx-sy) + (sy-sz)*(sy-sz) + (sz-sx)*(sz-sx)
         vm2 = 0.5d0*vm2 + 3.d0*(txy*txy + txz*txz + tyz*tyz)
         if (vm2 .le. 0.d0) cycle
         vm  = dsqrt(vm2)

c-------- phi (same as in resultsmech)
         phi = vm/(rho_p*sig0) + eps_relax - eps_relax/rho_eff
         if (phi .le. 0.d0) cycle

c-------- piecewise derivative of rho_eff = max(rho, rho_min)
         if (rho .gt. rho_min) then
            heav = 1.d0
         else
            heav = 0.d0
         endif

c-------- dphi/drho (explicit)
         dphi = heav * ( - penal*vm/(sig0 * rho_eff**(penal+1)) +
     &                   eps_relax/(rho_eff*rho_eff) )

c-------- accumulate explicit term for element i
         djdrho(i) = djdrho(i) + alpha * (phi**(pexp-1)) * dphi *
     &               xsj * weight

      enddo

      return
      end
