      subroutine pnorm_rhs(co,kon,ipkon,lakon,ne,
     &     stx,xstiff,mi,rhs,alpha,p,design,penal,
     &     sig0,eps_relax,rho_min,
     &     nea,neb,list,ilist)

c  Adjoint RHS for epsilon-relaxed stress p-norm, C3D4 only
c  J = ( sum w * phi^p )^(1/p),  phi = vm/(rho^p*sig0) + eps - eps/rho_eff
c  dJ/dsigma = alpha * p * phi^(p-1) * (1/(rho^p*sig0)) * dvm/dsigma
c  RHS = ∫ B^T [ D * (dJ/dsigma) ] dV

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
      integer i,k,jj,indexe,nope,mint3d,iflag
      integer m1,m2,m3,j
      integer konl(4)
      real*8 xl(3,4),shp(4,4),xsj,xi,et,ze,weight
      real*8 sx,sy,sz,txy,txz,tyz,vm2,vm
      real*8 rho,rho_eff,rho_p,phi,fac,invvm
      real*8 g(6),ptv(6),pt(3,3)

c--- Gauss rules
      include 'gauss.f'

c--- init
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

c-------- gather nodal connectivity and coords
         indexe = ipkon(i)
         do j=1,4
            konl(j) = kon(indexe+j)
         enddo
         do j=1,4
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
         !print *, 'gp= ', xi
         
c-------- shapes & jacobian
         call shape4tet(xi,et,ze,xl,xsj,shp,iflag)

c-------- stresses at this IP (Cauchy in Voigt order)
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
         invvm = 1.d0/vm

c-------- epsilon-relaxed phi (same as in resultsmech)
         phi = vm/(rho_p*sig0) + eps_relax - eps_relax/rho_eff
         if (phi .le. 0.d0) cycle

c-------- dJ/dsigma factor
         fac = alpha * (phi**(p-1)) / (rho_p*sig0)

c-------- dvm/dsigma in 3D Voigt, then scale by fac
         g(1) = fac * ((2.d0*sx - sy - sz) * 0.5d0 * invvm)
         g(2) = fac * ((2.d0*sy - sx - sz) * 0.5d0 * invvm)
         g(3) = fac * ((2.d0*sz - sx - sy) * 0.5d0 * invvm)
         g(4) = fac * (3.d0*txy * invvm)
         g(5) = fac * (3.d0*txz * invvm)
         g(6) = fac * (3.d0*tyz * invvm)

c-------- ptv = D * g (use increments to avoid long lines)
         ptv(1)=0.d0
         ptv(1)=ptv(1)+g(1)*xstiff( 1,jj,i)
         ptv(1)=ptv(1)+g(2)*xstiff( 2,jj,i)
         ptv(1)=ptv(1)+g(3)*xstiff( 4,jj,i)
         ptv(1)=ptv(1)+g(4)*xstiff( 7,jj,i)
         ptv(1)=ptv(1)+g(5)*xstiff(11,jj,i)
         ptv(1)=ptv(1)+g(6)*xstiff(16,jj,i)

         ptv(2)=0.d0
         ptv(2)=ptv(2)+g(1)*xstiff( 2,jj,i)
         ptv(2)=ptv(2)+g(2)*xstiff( 3,jj,i)
         ptv(2)=ptv(2)+g(3)*xstiff( 5,jj,i)
         ptv(2)=ptv(2)+g(4)*xstiff( 8,jj,i)
         ptv(2)=ptv(2)+g(5)*xstiff(12,jj,i)
         ptv(2)=ptv(2)+g(6)*xstiff(17,jj,i)

         ptv(3)=0.d0
         ptv(3)=ptv(3)+g(1)*xstiff( 4,jj,i)
         ptv(3)=ptv(3)+g(2)*xstiff( 5,jj,i)
         ptv(3)=ptv(3)+g(3)*xstiff( 6,jj,i)
         ptv(3)=ptv(3)+g(4)*xstiff( 9,jj,i)
         ptv(3)=ptv(3)+g(5)*xstiff(13,jj,i)
         ptv(3)=ptv(3)+g(6)*xstiff(18,jj,i)

         ptv(4)=0.d0
         ptv(4)=ptv(4)+g(1)*xstiff( 7,jj,i)
         ptv(4)=ptv(4)+g(2)*xstiff( 8,jj,i)
         ptv(4)=ptv(4)+g(3)*xstiff( 9,jj,i)
         ptv(4)=ptv(4)+g(4)*xstiff(10,jj,i)
         ptv(4)=ptv(4)+g(5)*xstiff(14,jj,i)
         ptv(4)=ptv(4)+g(6)*xstiff(19,jj,i)

         ptv(5)=0.d0
         ptv(5)=ptv(5)+g(1)*xstiff(11,jj,i)
         ptv(5)=ptv(5)+g(2)*xstiff(12,jj,i)
         ptv(5)=ptv(5)+g(3)*xstiff(13,jj,i)
         ptv(5)=ptv(5)+g(4)*xstiff(14,jj,i)
         ptv(5)=ptv(5)+g(5)*xstiff(15,jj,i)
         ptv(5)=ptv(5)+g(6)*xstiff(20,jj,i)

         ptv(6)=0.d0
         ptv(6)=ptv(6)+g(1)*xstiff(16,jj,i)
         ptv(6)=ptv(6)+g(2)*xstiff(17,jj,i)
         ptv(6)=ptv(6)+g(3)*xstiff(18,jj,i)
         ptv(6)=ptv(6)+g(4)*xstiff(19,jj,i)
         ptv(6)=ptv(6)+g(5)*xstiff(20,jj,i)
         ptv(6)=ptv(6)+g(6)*xstiff(21,jj,i)

c-------- Voigt -> tensor for assembly
         pt(1,1)=ptv(1)
         pt(2,2)=ptv(2)
         pt(3,3)=ptv(3)
         pt(2,1)=ptv(4)
         pt(3,1)=ptv(5)
         pt(3,2)=ptv(6)
         pt(1,2)=pt(2,1)
         pt(1,3)=pt(3,1)
         pt(2,3)=pt(3,2)

c-------- assemble RHS: r += ∫ B^T P dV
         do m1=1,nope
            do m2=1,3
               do m3=1,3
                  rhs(m2,konl(m1)) = rhs(m2,konl(m1)) +
     &               xsj * pt(m2,m3) * shp(m3,m1) * weight
               enddo
            enddo
         enddo

      enddo  ! elements

      return
      end
