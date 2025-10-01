c=======================================================================
c  Assemble adjoint RHS for stress p-norm objective
c
c  J = ( sum_g H_e * w_g * vm_g^p / sum_g H_e * w_g )^(1/p)
c
c  Discrete derivative at Gauss point:
c      dJ/dsigma = alpha * H_e * vm^(p-1) * d(vm)/d(sigma)
c  where alpha = J^(1-p) / (sum_g H_e * w_g)
c
c  We then assemble:  RHS = ∫ B^T * D^T * (dJ/dsigma) dV
c  (D is symmetric => D^T = D; CalculiX assembly uses the same pattern
c  as internal forces: σ replaced by P = D * (dJ/dσ) in the same
c  3×3 mapping and multiplied by ∂N/∂x and jacobian/weight.)
c=======================================================================
      subroutine pnorm_rhs(co,kon,ipkon,lakon,ne,
     &     stx,xstiff,mi,rhs,alpha,pexp,design,
     &     nea,neb,list,ilist)

      implicit none

c--- mesh / topology
      integer            ne,kon(*),ipkon(*),mi(*),nea,neb,list
      integer            ilist(*)
      character*8        lakon(*),lakonl

c--- geometry
      real*8             co(3,*)

c--- integration-point data (already computed in resultsmech)
c    stx: Cauchy stresses at IPs (Voigt: 1..6)
c    xstiff: 21 independent entries of elasticity matrix D per IP
      real*8             stx(6,mi(1),*),xstiff(27,mi(1),*)

c--- output RHS (same shape as mechanical force array)
      real*8             rhs(0:mi(2),*)

c--- p-norm controls passed in
      real*8             alpha,pexp,design(*)

c --- element weightin
      real*8 dens, dens_eff, rho_min, wexp      

c--- working vars
      integer            i,j,k,jj,kk,nope,mint3d,indexe
      integer            m1,m2,m3,ki,kl
      integer            iflag,imat,istrainfree
      real*8             xl(3,26),shp(4,26),xsj
      real*8             xi,et,ze,weight
      real*8             sx,sy,sz,txy,txz,tyz,vm2,vm
      real*8             g(6),ptvec(6),pt(3,3)
      real*8             Hloc, p, epsvm
      integer            konl(26)

c--- gauss-bookkeeping and (rarely used here) layer vars
      integer            nlayer,mint2d,nopes,ilayer
      real*8             tlayer(4),dlayer(4),xlayer(mi(3),4),xsj2(3)
      real*8             xl2(3,8),xs2(3,7),shp2(7,8)

c--- shape helper (hourglass/incompatible, etc.)
      real*8             gs(8,4),a

c--- pull in Gauss rules
      include 'gauss.f'

c-----------------------------------------------------------------------
c  Initialize
c-----------------------------------------------------------------------
      iflag = 3
      p     = pexp
      epsvm = 1.d-16

c  Zero RHS once (caller may pass a pre-zeroed array; doing it here is safe)
      do i=1,mi(2)
         do j=1,mi(1)
c           nothing: mi(1) is #gp per element; rhs is nodal-sized.
         enddo
      enddo

c-----------------------------------------------------------------------
c  Element loop over the thread slice [nea..neb]
c-----------------------------------------------------------------------
      do k=nea,neb

         if(list.eq.1) then
            i = ilist(k)
         else
            i = k
         endif

         if(ipkon(i).lt.0) cycle   ! skip inactive
         lakonl = lakon(i)

c        Only structural elements (skip fluids/user/fluid-coupled)
         if((lakonl(1:1).eq.'F').or.(lakonl(1:7).eq.'DCOUP3D')) cycle

c        Reactivated element handling (not used here, but keep parity)
         if( (1).eq.0 ) then
            istrainfree=0
         endif

c-------- element connectivity and type -------------------------------
         indexe = ipkon(i)

c        #nodes per element
         if(lakonl(1:5).eq.'C3D8I') then
            nope = 11
         elseif(lakonl(4:5).eq.'20') then
            nope = 20
         elseif(lakonl(4:4).eq.'8') then
            nope = 8
         elseif(lakonl(4:5).eq.'10') then
            nope = 10
         elseif(lakonl(4:4).eq.'4') then
            nope = 4
         elseif(lakonl(4:5).eq.'15') then
            nope = 15
         elseif(lakonl(4:4).eq.'6') then
            nope = 6
         else
            cycle
         endif

c        #gauss points per element (solids only here)
         if(lakonl(4:5).eq.'8R') then
            mint3d=1
         elseif(lakonl(4:7).eq.'20RB') then
c           (beam variants omitted in RHS assembly)
            mint3d=50
         elseif( (lakonl(4:4).eq.'8') .or.
     &           (lakonl(4:6).eq.'20R') ) then
            mint3d=8
         elseif(lakonl(4:4).eq.'2') then
            mint3d=27
         elseif(lakonl(4:5).eq.'10') then
            mint3d=4
         elseif(lakonl(4:4).eq.'4') then
            mint3d=1
         elseif(lakonl(4:5).eq.'15') then
            mint3d=9
         elseif(lakonl(4:4).eq.'6') then
            mint3d=2
         else
            cycle
         endif

c        gather nodal coords for this element
         do j=1,nope
            konl(j) = kon(indexe+j)
            xl(1,j) = co(1,konl(j))
            xl(2,j) = co(2,konl(j))
            xl(3,j) = co(3,konl(j))
         enddo

c        SIMP-style objective weight
         dens = design(i)
         if (dens .lt. 0.d0) dens = 0.d0
         if (dens .gt. 1.d0) dens = 1.d0

c        choose the same paramters used  when accumulating the p-norm
c        rho_min: density floor (set to 0.d0 if we want voids)
c        wexp: objective's density exponent (usually 1.d0)
         rho_min = 1.d-3
         wexp = 1.d0

         dens_eff = rho_min + (1.d0 - rho_min) * dens
         Hloc     = dens_eff**wexp


c--------------------------------------------------------------------
c        Gauss loop
c--------------------------------------------------------------------
         do jj=1,mint3d

c           Gauss point and weight
            if(lakonl(4:5).eq.'8R') then
               xi=gauss3d1(1,jj)
               et=gauss3d1(2,jj)
               ze=gauss3d1(3,jj)
               weight=weight3d1(jj)
            elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) then
               xi=gauss3d2(1,jj)
               et=gauss3d2(2,jj)
               ze=gauss3d2(3,jj)
               weight=weight3d2(jj)
            elseif(lakonl(4:4).eq.'2') then
               xi=gauss3d3(1,jj)
               et=gauss3d3(2,jj)
               ze=gauss3d3(3,jj)
               weight=weight3d3(jj)
            elseif(lakonl(4:5).eq.'10') then
               xi=gauss3d5(1,jj)
               et=gauss3d5(2,jj)
               ze=gauss3d5(3,jj)
               weight=weight3d5(jj)
            elseif(lakonl(4:4).eq.'4') then
               xi=gauss3d4(1,jj)
               et=gauss3d4(2,jj)
               ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
            elseif(lakonl(4:5).eq.'15') then
               xi=gauss3d8(1,jj)
               et=gauss3d8(2,jj)
               ze=gauss3d8(3,jj)
               weight=weight3d8(jj)
            elseif(lakonl(4:4).eq.'6') then
               xi=gauss3d7(1,jj)
               et=gauss3d7(2,jj)
               ze=gauss3d7(3,jj)
               weight=weight3d7(jj)
            endif

c           Shape derivatives and Jacobian at (xi,et,ze)
            if(lakonl(1:5).eq.'C3D8R') then
               call shape8hr(xl,xsj,shp,gs,a)
            elseif(lakonl(1:5).eq.'C3D8I') then
               call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.20) then
               call shape20h(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.8) then
               call shape8h(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.10) then
               call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.4) then
               call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.15) then
               call shape15w(xi,et,ze,xl,xsj,shp,iflag)
            else
               call shape6w(xi,et,ze,xl,xsj,shp,iflag)
            endif

c----------- Von Mises and its sigma-gradient (no weight here!) ------
c           Pull Cauchy stress at this IP (saved earlier in stx)
            sx  = stx(1,jj,i)
            sy  = stx(2,jj,i)
            sz  = stx(3,jj,i)
            txy = stx(4,jj,i)
            txz = stx(5,jj,i)
            tyz = stx(6,jj,i)

c           vm (same formula used when accumulating in resultsmech)
            vm2 = (sx-sy)*(sx-sy) + (sy-sz)*(sy-sz) + (sz-sx)*(sz-sx)
            vm2 = 0.5d0*vm2 + 3.d0*(txy*txy + txz*txz + tyz*tyz)
            if(vm2.le.0.d0) cycle
            vm  = dsqrt(vm2)
            if(vm.le.epsvm) cycle

c           dJ/dsigma (Voigt) at this IP (no jacobian/weight here):
c             g = alpha * H_e * vm^(p-1) * d(vm)/d(sigma)
c             with d(vm)/d(sigma) =
c               [ (2σx-σy-σz)/(2vm), (2σy-σx-σz)/(2vm), (2σz-σx-σy)/(2vm),
c                 3τxy/vm, 3τxz/vm, 3τyz/vm ]
            g(1) = alpha * Hloc * (vm**(p-1)) * ( (2.d0*sx - sy - sz)/(2.d0*vm) )
            g(2) = alpha * Hloc * (vm**(p-1)) * ( (2.d0*sy - sx - sz)/(2.d0*vm) )
            g(3) = alpha * Hloc * (vm**(p-1)) * ( (2.d0*sz - sx - sy)/(2.d0*vm) )
            g(4) = alpha * Hloc * (vm**(p-1)) * ( 3.d0*txy / vm )
            g(5) = alpha * Hloc * (vm**(p-1)) * ( 3.d0*txz / vm )
            g(6) = alpha * Hloc * (vm**(p-1)) * ( 3.d0*tyz / vm )

c----------- P = D * g  (Voigt).  D is symmetric, stored in xstiff(1:21)
c           Voigt mapping (same as in resultsmech: mattyp=3 branch)
            ptvec(1)= g(1)*xstiff(1 ,jj,i)+g(2)*xstiff(2 ,jj,i)+
     &                g(3)*xstiff(4 ,jj,i)+g(4)*xstiff(7 ,jj,i)+
     &                g(5)*xstiff(11,jj,i)+g(6)*xstiff(16,jj,i)
            ptvec(2)= g(1)*xstiff(2 ,jj,i)+g(2)*xstiff(3 ,jj,i)+
     &                g(3)*xstiff(5 ,jj,i)+g(4)*xstiff(8 ,jj,i)+
     &                g(5)*xstiff(12,jj,i)+g(6)*xstiff(17,jj,i)
            ptvec(3)= g(1)*xstiff(4 ,jj,i)+g(2)*xstiff(5 ,jj,i)+
     &                g(3)*xstiff(6 ,jj,i)+g(4)*xstiff(9 ,jj,i)+
     &                g(5)*xstiff(13,jj,i)+g(6)*xstiff(18,jj,i)
            ptvec(4)= g(1)*xstiff(7 ,jj,i)+g(2)*xstiff(8 ,jj,i)+
     &                g(3)*xstiff(9 ,jj,i)+g(4)*xstiff(10,jj,i)+
     &                g(5)*xstiff(14,jj,i)+g(6)*xstiff(19,jj,i)
            ptvec(5)= g(1)*xstiff(11,jj,i)+g(2)*xstiff(12,jj,i)+
     &                g(3)*xstiff(13,jj,i)+g(4)*xstiff(14,jj,i)+
     &                g(5)*xstiff(15,jj,i)+g(6)*xstiff(20,jj,i)
            ptvec(6)= g(1)*xstiff(16,jj,i)+g(2)*xstiff(17,jj,i)+
     &                g(3)*xstiff(18,jj,i)+g(4)*xstiff(19,jj,i)+
     &                g(5)*xstiff(20,jj,i)+g(6)*xstiff(21,jj,i)

c----------- Map Voigt ptvec() into 3x3 tensor for CalculiX assembly
            pt(1,1)=ptvec(1)
            pt(2,2)=ptvec(2)
            pt(3,3)=ptvec(3)
            pt(2,1)=ptvec(4)
            pt(3,1)=ptvec(5)
            pt(3,2)=ptvec(6)
            pt(1,2)=pt(2,1)
            pt(1,3)=pt(3,1)
            pt(2,3)=pt(3,2)

c----------- Assemble: rhs += ∫ (∂N/∂x)^T : P  dV
c           (identical stencil to internal force with σ→P)
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
