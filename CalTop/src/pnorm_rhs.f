c=======================================================================
c  Assemble adjoint RHS for ε-relaxed stress p-norm
c
c  J = ( Σ_g  w_g * φ_g^p )^(1/p) ,   w_g = xsj * weight
c
c  φ = max( 0 ,  vm/(ρ^p * σ0) + ε - ε/ρ )
c   with ρ = clamp(design_e, 0..1), ρ_eff = max(ρ, ρ_min), ρ^p = ρ_eff^penal
c
c  For φ>0:
c    dJ/dσ = α * φ^(p-1) * dφ/dσ
c    dφ/dσ = (1/(ρ^p * σ0)) * d(vm)/dσ
c  where α = J^(1-p)   (NOTE: no 1/Σw normalization)
c
c  We assemble   RHS = ∫ B^T * [ D * (dJ/dσ) ] dV
c  Implementation follows CalculiX internal-force stencil: σ → P = D*g
c=======================================================================
      subroutine pnorm_rhs(co,kon,ipkon,lakon,ne,
     &     stx,xstiff,mi,rhs,alpha,pexp,design,penal,
     &     sig0,eps_relax,rho_min,
     &     nea,neb,list,ilist)

      implicit none

c--- mesh / topology
      integer            ne,kon(*),ipkon(*),mi(*),nea,neb,list
      integer            ilist(*)
      character*8        lakon(*),lakonl

c--- geometry
      real*8             co(3,*)

c--- integration-point (already computed in resultsmech)
      real*8             stx(6,mi(1),*), xstiff(27,mi(1),*)

c--- output RHS (same layout as mechanical force array)
      real*8             rhs(0:mi(2),*)

c--- controls / parameters
      real*8             alpha, pexp, design(*), penal, sig0
      real*8             eps_relax, rho_min

c--- working vars
      integer            i,j,k,jj,nope,mint3d,indexe
      integer            m1,m2,m3,iflag
      real*8             xl(3,26),shp(4,26),xsj
      real*8             xi,et,ze,weight
      real*8             sx,sy,sz,txy,txz,tyz,vm2,vm
      real*8             rho, rho_eff, rho_p, phi, p, fac
      real*8             g(6),ptv(6),pt(3,3)
      integer            konl(26)

c--- shape helpers
      real*8             gs(8,4),a

c--- Gauss rules
      include 'gauss.f'

c-----------------------------------------------------------------------
c  Init
c-----------------------------------------------------------------------
      iflag = 3
      p     = pexp

c-----------------------------------------------------------------------
c  Element loop over [nea..neb]
c-----------------------------------------------------------------------
      do k=nea,neb

         if (list.eq.1) then
            i = ilist(k)
         else
            i = k
         endif

         if (ipkon(i).lt.0) cycle
         lakonl = lakon(i)

c        structural solids only
         if ( (lakonl(1:1).eq.'F') .or.
     &        (lakonl(1:7).eq.'DCOUP3D') ) cycle


c-------- #nodes per element
         indexe = ipkon(i)
         if     (lakonl(1:5).eq.'C3D8I') then
            nope = 11
         elseif (lakonl(4:5).eq.'20') then
            nope = 20
         elseif (lakonl(4:4).eq.'8') then
            nope = 8
         elseif (lakonl(4:5).eq.'10') then
            nope = 10
         elseif (lakonl(4:4).eq.'4') then
            nope = 4
         elseif (lakonl(4:5).eq.'15') then
            nope = 15
         elseif (lakonl(4:4).eq.'6') then
            nope = 6
         else
            cycle
         endif

c-------- #gauss points per element (solid variants)
         if     (lakonl(4:5).eq.'8R') then
            mint3d = 1
         elseif ( (lakonl(4:4).eq.'8') .or.
     &            (lakonl(4:6).eq.'20R') ) then

            mint3d = 8
         elseif (lakonl(4:4).eq.'2') then
            mint3d = 27
         elseif (lakonl(4:5).eq.'10') then
            mint3d = 4
         elseif (lakonl(4:4).eq.'4') then
            mint3d = 1
         elseif (lakonl(4:5).eq.'15') then
            mint3d = 9
         elseif (lakonl(4:4).eq.'6') then
            mint3d = 2
         else
            cycle
         endif

c-------- gather nodal coords
         do j=1,nope
            konl(j) = kon(indexe+j)
            xl(1,j) = co(1,konl(j))
            xl(2,j) = co(2,konl(j))
            xl(3,j) = co(3,konl(j))
         enddo

c-------- element density (filtered is already in design())
         rho = design(i)
         if (rho .lt. 0.d0) rho = 0.d0
         if (rho .gt. 1.d0) rho = 1.d0
         rho_eff = dmax1(rho, rho_min)
         rho_p   = rho_eff**penal

c--------------------------------------------------------------------
c        Gauss loop
c--------------------------------------------------------------------
         do jj=1,mint3d

c           Gauss point & weight
            if     (lakonl(4:5).eq.'8R') then
               xi=gauss3d1(1,jj); et=gauss3d1(2,jj); ze=gauss3d1(3,jj)
               weight=weight3d1(jj)
            elseif ( (lakonl(4:4).eq.'8') .or.
     &               (lakonl(4:6).eq.'20R') ) then

               xi=gauss3d2(1,jj); et=gauss3d2(2,jj); ze=gauss3d2(3,jj)
               weight=weight3d2(jj)
            elseif (lakonl(4:4).eq.'2') then
               xi=gauss3d3(1,jj); et=gauss3d3(2,jj); ze=gauss3d3(3,jj)
               weight=weight3d3(jj)
            elseif (lakonl(4:5).eq.'10') then
               xi=gauss3d5(1,jj); et=gauss3d5(2,jj); ze=gauss3d5(3,jj)
               weight=weight3d5(jj)
            elseif (lakonl(4:4).eq.'4') then
               xi=gauss3d4(1,jj); et=gauss3d4(2,jj); ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
            elseif (lakonl(4:5).eq.'15') then
               xi=gauss3d8(1,jj); et=gauss3d8(2,jj); ze=gauss3d8(3,jj)
               weight=weight3d8(jj)
            elseif (lakonl(4:4).eq.'6') then
               xi=gauss3d7(1,jj); et=gauss3d7(2,jj); ze=gauss3d7(3,jj)
               weight=weight3d7(jj)
            endif

c           Shapes & Jacobian
            if     (lakonl(1:5).eq.'C3D8R') then
               call shape8hr(xl,xsj,shp,gs,a)
            elseif (lakonl(1:5).eq.'C3D8I') then
               call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
            elseif (nope.eq.20) then
               call shape20h(xi,et,ze,xl,xsj,shp,iflag)
            elseif (nope.eq.8) then
               call shape8h(xi,et,ze,xl,xsj,shp,iflag)
            elseif (nope.eq.10) then
               call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif (nope.eq.4) then
               call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif (nope.eq.15) then
               call shape15w(xi,et,ze,xl,xsj,shp,iflag)
            else
               call shape6w(xi,et,ze,xl,xsj,shp,iflag)
            endif

c----------- stresses at this IP
            sx  = stx(1,jj,i)
            sy  = stx(2,jj,i)
            sz  = stx(3,jj,i)
            txy = stx(4,jj,i)
            txz = stx(5,jj,i)
            tyz = stx(6,jj,i)

c----------- von Mises
            vm2 = (sx-sy)*(sx-sy) + (sy-sz)*(sy-sz) + (sz-sx)*(sz-sx)
            vm2 = 0.5d0*vm2 + 3.d0*(txy*txy + txz*txz + tyz*tyz)
            if (vm2.le.0.d0) cycle
            vm  = dsqrt(vm2)

c----------- ε-relaxed φ
            phi = vm/(rho_p*sig0) + eps_relax -
     &             eps_relax/dmax1(rho,1.d-16)

            if (phi.le.0.d0) cycle

c----------- dJ/dσ = α * φ^(p-1) * (1/(ρ^p σ0)) * dvm/dσ   (Voigt)
            fac = alpha * (phi**(p-1)) / (rho_p*sig0)

            g(1) =  fac*((2.d0*sx - sy - sz)/(2.d0*vm))
            g(2) =  fac*((2.d0*sy - sx - sz)/(2.d0*vm))
            g(3) =  fac*((2.d0*sz - sx - sy)/(2.d0*vm))
            g(4) =  fac*(3.d0*txy/vm)
            g(5) =  fac*(3.d0*txz/vm)
            g(6) =  fac*(3.d0*tyz/vm)

c----------- P = D * g  (Voigt); D entries use xstiff(1..21)
            ptv(1)= g(1)*xstiff(1 ,jj,i)+g(2)*xstiff(2 ,jj,i)+
     &              g(3)*xstiff(4 ,jj,i)+g(4)*xstiff(7 ,jj,i)+
     &              g(5)*xstiff(11,jj,i)+g(6)*xstiff(16,jj,i)
            ptv(2)= g(1)*xstiff(2 ,jj,i)+g(2)*xstiff(3 ,jj,i)+
     &              g(3)*xstiff(5 ,jj,i)+g(4)*xstiff(8 ,jj,i)+
     &              g(5)*xstiff(12,jj,i)+g(6)*xstiff(17,jj,i)
            ptv(3)= g(1)*xstiff(4 ,jj,i)+g(2)*xstiff(5 ,jj,i)+
     &              g(3)*xstiff(6 ,jj,i)+g(4)*xstiff(9 ,jj,i)+
     &              g(5)*xstiff(13,jj,i)+g(6)*xstiff(18,jj,i)
            ptv(4)= g(1)*xstiff(7 ,jj,i)+g(2)*xstiff(8 ,jj,i)+
     &              g(3)*xstiff(9 ,jj,i)+g(4)*xstiff(10,jj,i)+
     &              g(5)*xstiff(14,jj,i)+g(6)*xstiff(19,jj,i)
            ptv(5)= g(1)*xstiff(11,jj,i)+g(2)*xstiff(12,jj,i)+
     &              g(3)*xstiff(13,jj,i)+g(4)*xstiff(14,jj,i)+
     &              g(5)*xstiff(15,jj,i)+g(6)*xstiff(20,jj,i)
            ptv(6)= g(1)*xstiff(16,jj,i)+g(2)*xstiff(17,jj,i)+
     &              g(3)*xstiff(18,jj,i)+g(4)*xstiff(19,jj,i)+
     &              g(5)*xstiff(20,jj,i)+g(6)*xstiff(21,jj,i)

c----------- Voigt → tensor for assembly
            pt(1,1)=ptv(1); pt(2,2)=ptv(2); pt(3,3)=ptv(3)
            pt(2,1)=ptv(4); pt(3,1)=ptv(5); pt(3,2)=ptv(6)
            pt(1,2)=pt(2,1); pt(1,3)=pt(3,1); pt(2,3)=pt(3,2)

c----------- assemble RHS
            do m1=1,nope
               do m2=1,3
                  do m3=1,3
                     rhs(m2,konl(m1)) = rhs(m2,konl(m1)) +
     &                    xsj * pt(m2,m3) * shp(m3,m1) * weight
                  enddo
               enddo
            enddo

         enddo   ! gauss
      enddo      ! element

      return
      end
