!
!     CalTop: stress P-norm value for geometrically nonlinear (NLGEOM)
!     analysis, using the von Mises stress of the CAUCHY stress.
!
      subroutine stresspnorm_nl(co,kon,ipkon,lakon,ne,v,
     &  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
     &  iprestr,eme,iperturb,fn,iout,qa,vold,nmethod,
     &  veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &  xstateini,xstiff,xstate,npmat_,matname,mi,ielas,icmd,
     &  ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,
     &  springarea,reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,
     &  ikin,nal,ne0,thicke,emeini,pslavsurf,
     &  pmastsurf,mortar,clearini,nea,neb,ielprop,prop,kscale,
     &  list,ilist, design, penal, sig0, eps_relax, rho_min, pexp)
!
!     Same argument list as stresspnorm (linear). Differences:
!       - stress is the penalized Cauchy stress of a total Lagrangian
!         St. Venant-Kirchhoff model (see cauchy_vm_nl.f)
!       - stx is NOT written (it keeps the Cauchy stress from
!         resultsmech for post-processing)
!     Per element:  vm_e  = (1/mint3d) sum_gp vm_gp
!                   phi_e = vm_e/sig0 + eps_relax - eps_relax/rho_eff
!     Returns qa(3) = sum_e phi_e**pexp
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*80 amat,matname(*)
!
      integer kon(*),konl(26),nea,neb,mi(*),
     &  nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),
     &  ielorien(mi(3),*),ntmat_,ipkon(*),ne0,iflag,kscale,
     &  istep,iinc,ne,mattyp,ithermal(2),iprestr,i,j,k,jj,
     &  nener,indexe,nope,norien,iperturb(*),iout,
     &  nal,icmd,ihyper,nmethod,kode,imat,mint3d,iorien,ielas,
     &  istiff,ncmat_,nstate_,ikin,ielprop(*),
     &  nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,calcul_fn,
     &  calcul_cauchy,calcul_qa,mortar,list,ilist(*),m
!
      real*8 co(3,*),v(0:mi(2),*),shp(4,26),stiini(6,mi(1),*),
     &  stx(6,mi(1),*),xl(3,26),vl(0:mi(2),26),prop(*),
     &  elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),
     &  alcon(0:6,ntmat_,*),vini(0:mi(2),*),
     &  alzero(*),orab(7,*),elas(21),rho,fn(0:mi(2),*),
     &  t0(*),t1(*),prestr(6,mi(1),*),eme(6,mi(1),*),
     &  vold(0:mi(2),*),veold(0:mi(2),*),
     &  springarea(2,*),elconloc(21),eth(6),
     &  ener(mi(1),*),eei(6,mi(1),*),enerini(mi(1),*),
     &  xi,et,ze,t0l,t1l,
     &  xsj,qa(*),dtime,weight,pgauss(3),time,ttime,
     &  plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mi(1),*),xstate(nstate_,mi(1),*),plconloc(802),
     &  xstateini(nstate_,mi(1),*),reltime,
     &  thicke(mi(3),*),emeini(6,mi(1),*),clearini(3,9,*),
     &  pslavsurf(3,*),pmastsurf(6,*)
!
      real*8 design(*), penal
      real*8 vm, vm_gp, vm_sum, g_sump, pexp
      real*8 rho_e, rho_min, rho_eff, rho_p, eps_relax, sig0, phi
      real*8 sig(6), dum(1), epen, phi_max
!
      intent(in) co,kon,ipkon,lakon,ne,v,
     &  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielorien,norien,orab,ntmat_,t0,t1,ithermal,
     &  iprestr,iperturb,iout,vold,nmethod,
     &  veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &  xstateini,xstate,npmat_,matname,mi,ielas,icmd,
     &  ncmat_,nstate_,stiini,vini,enerini,istep,iinc,
     &  springarea,reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,
     &  ikin,ne0,thicke,pslavsurf,
     &  pmastsurf,mortar,clearini,nea,neb,ielprop,prop,kscale,
     &  list,ilist,design, penal, sig0, eps_relax, rho_min, pexp
!
      intent(inout) nal,qa,fn,xstiff,ener,eme,eei,stx,ielmat,prestr,
     &  emeini
!
      include "gauss.f"
!
      iflag=3
      t0l=0.d0
      t1l=0.d0
      nal=0
      qa(3)=-1.d0
      qa(4)=0.d0
      g_sump=0.d0
      phi_max=0.d0
!
      do m=nea,neb
!
         if(list.eq.1) then
            i=ilist(m)
         else
            i=m
         endif
         if(ipkon(i).lt.0) cycle
!
         lakonl=lakon(i)
!
         if(lakonl(7:8).ne.'LC') then
            imat=ielmat(1,i)
            amat=matname(imat)
            if(norien.gt.0) then
               iorien=ielorien(1,i)
            else
               iorien=0
            endif
            if(nelcon(1,imat).lt.0) then
               ihyper=1
            else
               ihyper=0
            endif
         endif
!
         indexe=ipkon(i)
!
!        Set number of nodes and integration pts for C3D4 and C3D10
         if(lakonl(4:5).eq.'10') then
            nope=10
            mint3d=4
         elseif(lakonl(4:4).eq.'4') then
            nope=4
            mint3d=1
         else
            cycle
         endif
!
!        Gather coordinates and displacements
         do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl(j))
               vl(k,j)=v(k,konl(j))
            enddo
         enddo
!
!        Density / SIMP
         rho_e=design(i)
         if (rho_e .lt. 0.d0) rho_e = 0.d0
         if (rho_e .gt. 1.d0) rho_e = 1.d0
         rho_eff=max(rho_e, rho_min)
         rho_p=rho_eff**penal
!
         vm_sum=0.d0
         do jj=1,mint3d
!
            if(lakonl(4:5).eq.'10') then
               xi=gauss3d5(1,jj)
               et=gauss3d5(2,jj)
               ze=gauss3d5(3,jj)
               weight=weight3d5(jj)
               call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            else
               xi=gauss3d4(1,jj)
               et=gauss3d4(2,jj)
               ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
               call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            endif
!
!           material data (isotropic, kode=2: elconloc(1)=E, (2)=nu)
            kode=nelcon(1,imat)
            istiff=0
            call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,
     &           nalcon,imat,amat,iorien,pgauss,orab,ntmat_,
     &           elas,rho,i,ithermal,alzero,mattyp,t0l,t1l,ihyper,
     &           istiff,elconloc,eth,kode,plicon,nplicon,
     &           plkcon,nplkcon,npmat_,plconloc,mi(1),dtime,jj,
     &           xstiff,ncmat_)
            if(kode.ne.2) then
               write(*,*) '*ERROR in stresspnorm_nl: only isotropic'
               write(*,*) '       elastic materials are supported'
               call exit(201)
            endif
            epen=elconloc(1)*rho_p
!
!           penalized Cauchy stress and its von Mises value
            call cauchy_vm_nl(shp,vl,mi,nope,epen,elconloc(2),0,
     &           sig,vm_gp,dum,dum)
            vm_sum=vm_sum+vm_gp
         enddo
!
         vm=vm_sum/dble(mint3d)
         phi=(vm/sig0)+eps_relax-(eps_relax/rho_eff)
         if (phi .lt. 0.d0) phi = 0.d0
         g_sump=g_sump+(phi**pexp)
         if(phi.gt.phi_max) phi_max=phi
      enddo
!
      qa(3)=g_sump
      qa(4)=phi_max
!
      return
      end

!
!     CalTop: adjoint RHS d(sum phi^p)/du (without the p and the
!     P^(1-p) prefactor) for the NLGEOM stress P-norm (Cauchy stress).
!
      subroutine pnorm_rhs_nl(co,kon,ipkon,lakon,ne,v,
     &  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
     &  iprestr,eme,iperturb,fn,iout,qa,vold,nmethod,
     &  veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &  xstateini,xstiff,xstate,npmat_,matname,mi,ielas,icmd,
     &  ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,
     &  springarea,reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,
     &  ikin,nal,ne0,thicke,emeini,pslavsurf,
     &  pmastsurf,mortar,clearini,nea,neb,ielprop,prop,kscale,
     &  list,ilist, rhs, design, penal, sig0, eps_relax,
     &  rho_min, pexp, alpha)
!
!     Same argument list as pnorm_rhs (linear). Differences:
!       - stress is the penalized Cauchy stress of a total Lagrangian
!         St. Venant-Kirchhoff model (see cauchy_vm_nl.f)
!       - stx is NOT written (it keeps the Cauchy stress from
!         resultsmech for post-processing)
!     Per element (same definitions as stresspnorm_nl):
!        rhs_e = phi_e**(pexp-1)/sig0 * (1/mint3d) sum_gp dvm_gp/du_e
!     scattered into rhs(1..3,node). Matches pnorm_rhs for the linear
!     case: the caller multiplies lambda by P^(1-p) afterwards.
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*80 amat,matname(*)
!
      integer kon(*),konl(26),nea,neb,mi(*),
     &  nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),
     &  ielorien(mi(3),*),ntmat_,ipkon(*),ne0,iflag,kscale,
     &  istep,iinc,ne,mattyp,ithermal(2),iprestr,i,j,k,jj,
     &  nener,indexe,nope,norien,iperturb(*),iout,
     &  nal,icmd,ihyper,nmethod,kode,imat,mint3d,iorien,ielas,
     &  istiff,ncmat_,nstate_,ikin,ielprop(*),
     &  nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,calcul_fn,
     &  calcul_cauchy,calcul_qa,mortar,list,ilist(*),m
!
      real*8 co(3,*),v(0:mi(2),*),shp(4,26),stiini(6,mi(1),*),
     &  stx(6,mi(1),*),xl(3,26),vl(0:mi(2),26),prop(*),
     &  elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),
     &  alcon(0:6,ntmat_,*),vini(0:mi(2),*),
     &  alzero(*),orab(7,*),elas(21),rho,fn(0:mi(2),*),
     &  t0(*),t1(*),prestr(6,mi(1),*),eme(6,mi(1),*),
     &  vold(0:mi(2),*),veold(0:mi(2),*),
     &  springarea(2,*),elconloc(21),eth(6),
     &  ener(mi(1),*),eei(6,mi(1),*),enerini(mi(1),*),
     &  xi,et,ze,t0l,t1l,
     &  xsj,qa(*),dtime,weight,pgauss(3),time,ttime,
     &  plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mi(1),*),xstate(nstate_,mi(1),*),plconloc(802),
     &  xstateini(nstate_,mi(1),*),reltime,
     &  thicke(mi(3),*),emeini(6,mi(1),*),clearini(3,9,*),
     &  pslavsurf(3,*),pmastsurf(6,*)
!
      real*8 design(*), penal
      real*8 vm, vm_gp, vm_sum, g_sump, pexp
      real*8 rho_e, rho_min, rho_eff, rho_p, eps_relax, sig0, phi
      real*8 sig(6), dum(1), epen, alpha, coeff
      real*8 rhs(0:mi(2),*), dvm(30), rhs_e(30), fgp(30)
!
      intent(in) co,kon,ipkon,lakon,ne,v,
     &  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielorien,norien,orab,ntmat_,t0,t1,ithermal,
     &  iprestr,iperturb,iout,vold,nmethod,
     &  veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &  xstateini,xstate,npmat_,matname,mi,ielas,icmd,
     &  ncmat_,nstate_,stiini,vini,enerini,istep,iinc,
     &  springarea,reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,
     &  ikin,ne0,thicke,pslavsurf,
     &  pmastsurf,mortar,clearini,nea,neb,ielprop,prop,kscale,
     &  list,ilist,design, penal, sig0, eps_relax, rho_min, pexp,
     &  alpha
!
      intent(inout) nal,qa,fn,xstiff,ener,eme,eei,stx,ielmat,prestr,
     &  emeini,rhs
!
      include "gauss.f"
!
      iflag=3
      t0l=0.d0
      t1l=0.d0
!
      do m=nea,neb
!
         if(list.eq.1) then
            i=ilist(m)
         else
            i=m
         endif
         if(ipkon(i).lt.0) cycle
!
         lakonl=lakon(i)
!
         if(lakonl(7:8).ne.'LC') then
            imat=ielmat(1,i)
            amat=matname(imat)
            if(norien.gt.0) then
               iorien=ielorien(1,i)
            else
               iorien=0
            endif
            if(nelcon(1,imat).lt.0) then
               ihyper=1
            else
               ihyper=0
            endif
         endif
!
         indexe=ipkon(i)
!
!        Set number of nodes and integration pts for C3D4 and C3D10
         if(lakonl(4:5).eq.'10') then
            nope=10
            mint3d=4
         elseif(lakonl(4:4).eq.'4') then
            nope=4
            mint3d=1
         else
            cycle
         endif
!
!        Gather coordinates and displacements
         do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl(j))
               vl(k,j)=v(k,konl(j))
            enddo
         enddo
!
!        Density / SIMP
         rho_e=design(i)
         if (rho_e .lt. 0.d0) rho_e = 0.d0
         if (rho_e .gt. 1.d0) rho_e = 1.d0
         rho_eff=max(rho_e, rho_min)
         rho_p=rho_eff**penal
!
         vm_sum=0.d0
         do k=1,3*nope
            rhs_e(k)=0.d0
         enddo
         do jj=1,mint3d
!
            if(lakonl(4:5).eq.'10') then
               xi=gauss3d5(1,jj)
               et=gauss3d5(2,jj)
               ze=gauss3d5(3,jj)
               weight=weight3d5(jj)
               call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            else
               xi=gauss3d4(1,jj)
               et=gauss3d4(2,jj)
               ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
               call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            endif
!
!           material data (isotropic, kode=2: elconloc(1)=E, (2)=nu)
            kode=nelcon(1,imat)
            istiff=0
            call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,
     &           nalcon,imat,amat,iorien,pgauss,orab,ntmat_,
     &           elas,rho,i,ithermal,alzero,mattyp,t0l,t1l,ihyper,
     &           istiff,elconloc,eth,kode,plicon,nplicon,
     &           plkcon,nplkcon,npmat_,plconloc,mi(1),dtime,jj,
     &           xstiff,ncmat_)
            if(kode.ne.2) then
               write(*,*) '*ERROR in pnorm_rhs_nl: only isotropic'
               write(*,*) '       elastic materials are supported'
               call exit(201)
            endif
            epen=elconloc(1)*rho_p
!
!           penalized Cauchy stress and its von Mises value
            call cauchy_vm_nl(shp,vl,mi,nope,epen,elconloc(2),2,
     &           sig,vm_gp,fgp,dvm)
            if(vm_gp.gt.0.d0) then
               vm_sum=vm_sum+vm_gp
               do k=1,3*nope
                  rhs_e(k)=rhs_e(k)+dvm(k)
               enddo
            endif
         enddo
!
         vm=vm_sum/dble(mint3d)
         phi=(vm/sig0)+eps_relax-(eps_relax/rho_eff)
         if (phi .lt. 0.d0) phi = 0.d0
!
!        coeff = phi_e^(p-1)/sig0 * (1/mint3d)
         coeff=(phi**(pexp-1))/(sig0*dble(mint3d))
!
!        scatter to global rhs(1..3,node)
         do j=1,nope
            do k=1,3
               rhs(k,konl(j))=rhs(k,konl(j))+coeff*rhs_e(3*(j-1)+k)
            enddo
         enddo
      enddo
!
      return
      end

!
!     CalTop: dJ/drho for the NLGEOM stress P-norm (Cauchy stress),
!     implicit (adjoint) + explicit parts, before the P^(1-p) prefactor.
!
      subroutine pnorm_implicit_nl(co,kon,ipkon,lakon,ne,mi,
     &     xstiff, v, lam, design, penal, pexp, relax, sig0,
     &     nea, neb, list, ilist, djdrho)
!
!     Same argument list as pnorm_implicit (linear). Per element:
!
!       f_int,e = sum_gp B_L^T S dV       (S: penalized PK2, total
!                                          Lagrangian internal force)
!       implicit = -(penal/rho) * lambda_e . f_int,e
!                  (dR/drho = d f_int/drho = (penal/rho) f_int, SIMP)
!       explicit = phi^(p-1) relax/rho^2 + phi^(p-1) vm penal/(rho sig0)
!                  (penalized Cauchy stress scales with rho^penal at
!                   fixed u, exactly as in the linear routine)
!
!     xstiff(1,jj,i) = E, xstiff(2,jj,i) = nu (isotropic, as stored by
!     resultsmech/linel).
!
      implicit none
!
      integer ne, kon(*), ipkon(*), mi(*), nea, neb, list, ilist(*)
      character*8 lakon(*), lakonl
      integer i,k,j,idx,nope,mint3d,iflag,jj,n
      integer konl(26)
!
      real*8 co(3,*), v(0:mi(2),*), lam(0:mi(2),*)
      real*8 xstiff(27,mi(1),*)
      real*8 design(*), penal, pexp, relax, sig0, djdrho(*)
!
      real*8 xl(3,26), vl(0:mi(2),26), shp(4,26), xsj, xi, et, ze,
     &  weight
      real*8 le(30), fint(30), fgp(30), sig(6), dum(1)
      real*8 rho, rho_eff, epen, vm, vm_gp, vm_sum, sige, dotlam
      real*8 expli1, expli2
!
      include 'gauss.f'
!
      iflag = 3
!
      do k = nea, neb
         if (list.eq.1) then
            i = ilist(k)
         else
            i = k
         endif
!
         if (ipkon(i).lt.0) cycle
         lakonl = lakon(i)
         if (lakonl(4:5).eq.'10') then
            nope = 10
            mint3d = 4
         elseif (lakonl(4:4).eq.'4') then
            nope = 4
            mint3d = 1
         else
            cycle
         endif
!
!------ gather coordinates, displacement and adjoint
         idx = ipkon(i)
         do j=1,nope
            konl(j) = kon(idx+j)
            xl(1,j) = co(1,konl(j))
            xl(2,j) = co(2,konl(j))
            xl(3,j) = co(3,konl(j))
            vl(1,j) = v(1,konl(j))
            vl(2,j) = v(2,konl(j))
            vl(3,j) = v(3,konl(j))
            le(3*(j-1)+1) = lam(1,konl(j))
            le(3*(j-1)+2) = lam(2,konl(j))
            le(3*(j-1)+3) = lam(3,konl(j))
         enddo
!
         rho = design(i)
!
!------ Gauss integration: f_int (penalized) and element vm
         do n=1,3*nope
            fint(n) = 0.d0
         enddo
         vm_sum = 0.d0
         do jj=1,mint3d
            if(lakonl(4:5).eq.'10') then
               xi=gauss3d5(1,jj); et=gauss3d5(2,jj); ze=gauss3d5(3,jj)
               weight=weight3d5(jj)
               call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            else
               xi=gauss3d4(1,jj); et=gauss3d4(2,jj); ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
               call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            endif
            xsj = dabs(xsj)
!
            epen = (rho**penal)*xstiff(1,jj,i)
            call cauchy_vm_nl(shp,vl,mi,nope,epen,xstiff(2,jj,i),1,
     &           sig,vm_gp,fgp,dum)
            vm_sum = vm_sum + vm_gp
            do n=1,3*nope
               fint(n) = fint(n) + fgp(n)*(xsj*weight)
            enddo
         enddo
!
!------ element-average vm and relaxed stress (as in stresspnorm_nl)
         vm = vm_sum / dble(mint3d)
         if (rho .lt. 0.d0) rho = 0.d0
         if (rho .gt. 1.d0) rho = 1.d0
         rho_eff = dmax1(rho, 1.d-6)
         sige = vm/sig0 + relax - relax/rho_eff
         if (sige .lt. 0.d0) sige = 0.d0
!
!------ implicit: -(penal/rho) lambda . f_int
         dotlam = 0.d0
         do n=1,3*nope
            dotlam = dotlam + le(n)*fint(n)
         enddo
         dotlam = (penal/rho_eff)*dotlam
!
         djdrho(i) = djdrho(i) - dotlam
         expli1 = (sige**(pexp-1))*relax/(rho_eff**2.d0)
         expli2 = (sige**(pexp-1))*vm*penal/(rho_eff*sig0)
         djdrho(i) = djdrho(i) + expli1 + expli2
      enddo
!
      return
      end

!
!     CalTop: Gauss-point kinematics for the NLGEOM stress P-norm
!
      subroutine cauchy_vm_nl(shp,vl,mi,nope,e,un,igrad,
     &     sig,vm,fgp,dvm)
!
!     Total Lagrangian, St. Venant-Kirchhoff, isotropic:
!
!        F     = I + dU/dX                (vl, shp from shape*tet, iflag=3)
!        E     = 1/2 (F^T F - I)          (Green-Lagrange)
!        S     = lambda tr(E) I + 2 mu E  (PK2, e = penalized modulus)
!        sigma = J^-1 F S F^T             (Cauchy)
!        vm    = von Mises of sigma
!
!     igrad = 0 : sig, vm
!     igrad = 1 : + fgp(3*nope)  = dE/du : S  (= (F S g_a)_k; times dV
!                                              gives the internal force)
!     igrad = 2 : + dvm(3*nope)  = d vm / d u_e
!
!     For node a, dof k (g_a = dN_a/dX, gx_a = F^-T g_a, h_a = F S g_a):
!        d sigma/d u_ak = -gx_a(k) sigma + (e_k h_a^T + h_a e_k^T)/J
!                         + F dS_ak F^T / J
!        dS_ak = lambda tr(dE) I + 2 mu dE,
!        dE_ij = 1/2 (F_ki g_aj + F_kj g_ai)
!     and d vm = (V sigma) . d sigma / vm (V: von Mises selector,
!     shear weight 3, tensorial shear components).
!
!     sig(6) is returned in CCX order: 11,22,33,12,13,23.
!
      implicit none
!
      integer mi(*),nope,igrad,a,k,i,j,m
      real*8 shp(4,26),vl(0:mi(2),26),e,un,sig(6),vm,
     &  fgp(*),dvm(*)
      real*8 xf(3,3),xe(3,3),s(3,3),t(3,3),fs(3,3),xfi(3,3),
     &  xj,um,al,tre,vm2,vs(3,3),g(3),gx(3),h(3),de(3,3),ds(3,3),
     &  dt(3,3),tmp(3,3),trde,dot
!
!     deformation gradient F = I + sum_a u_a (x) g_a
      do i=1,3
         do j=1,3
            xf(i,j)=0.d0
            do a=1,nope
               xf(i,j)=xf(i,j)+vl(i,a)*shp(j,a)
            enddo
         enddo
         xf(i,i)=xf(i,i)+1.d0
      enddo
!
!     Green-Lagrange strain
      do i=1,3
         do j=1,3
            xe(i,j)=0.d0
            do m=1,3
               xe(i,j)=xe(i,j)+xf(m,i)*xf(m,j)
            enddo
            xe(i,j)=0.5d0*xe(i,j)
         enddo
         xe(i,i)=xe(i,i)-0.5d0
      enddo
!
!     PK2 (SVK)
      um=e/(2.d0*(1.d0+un))
      al=un*e/((1.d0+un)*(1.d0-2.d0*un))
      tre=xe(1,1)+xe(2,2)+xe(3,3)
      do i=1,3
         do j=1,3
            s(i,j)=2.d0*um*xe(i,j)
         enddo
         s(i,i)=s(i,i)+al*tre
      enddo
!
!     J and Cauchy stress
      xj=xf(1,1)*(xf(2,2)*xf(3,3)-xf(2,3)*xf(3,2))
     &  -xf(1,2)*(xf(2,1)*xf(3,3)-xf(2,3)*xf(3,1))
     &  +xf(1,3)*(xf(2,1)*xf(3,2)-xf(2,2)*xf(3,1))
!
!     fs = F S
      do i=1,3
         do j=1,3
            fs(i,j)=0.d0
            do m=1,3
               fs(i,j)=fs(i,j)+xf(i,m)*s(m,j)
            enddo
         enddo
      enddo
!     t = F S F^T / J
      do i=1,3
         do j=1,3
            t(i,j)=0.d0
            do m=1,3
               t(i,j)=t(i,j)+fs(i,m)*xf(j,m)
            enddo
            t(i,j)=t(i,j)/xj
         enddo
      enddo
!
      sig(1)=t(1,1)
      sig(2)=t(2,2)
      sig(3)=t(3,3)
      sig(4)=t(1,2)
      sig(5)=t(1,3)
      sig(6)=t(2,3)
!
      vm2=0.5d0*((t(1,1)-t(2,2))**2+(t(2,2)-t(3,3))**2
     &          +(t(3,3)-t(1,1))**2)
     &   +3.d0*(t(1,2)**2+t(1,3)**2+t(2,3)**2)
      vm=dsqrt(max(vm2,0.d0))
!
      if(igrad.lt.1) return
!
!     internal-force density: fgp(ak) = (F S g_a)_k
      do a=1,nope
         do k=1,3
            fgp(3*(a-1)+k)=fs(k,1)*shp(1,a)+fs(k,2)*shp(2,a)
     &                    +fs(k,3)*shp(3,a)
         enddo
      enddo
!
      if(igrad.lt.2) return
!
      do m=1,3*nope
         dvm(m)=0.d0
      enddo
      if(vm.le.0.d0) return
!
!     V*sigma as a symmetric tensor (off-diagonal weight 3), so that
!     sum_ij vs(i,j)*dt(i,j) = (V sigma)_voigt . dsigma_voigt
      vs(1,1)=t(1,1)-0.5d0*(t(2,2)+t(3,3))
      vs(2,2)=t(2,2)-0.5d0*(t(1,1)+t(3,3))
      vs(3,3)=t(3,3)-0.5d0*(t(1,1)+t(2,2))
      vs(1,2)=1.5d0*t(1,2)
      vs(2,1)=vs(1,2)
      vs(1,3)=1.5d0*t(1,3)
      vs(3,1)=vs(1,3)
      vs(2,3)=1.5d0*t(2,3)
      vs(3,2)=vs(2,3)
!
!     F^-1 (cofactor / J)
      xfi(1,1)=(xf(2,2)*xf(3,3)-xf(2,3)*xf(3,2))/xj
      xfi(1,2)=(xf(1,3)*xf(3,2)-xf(1,2)*xf(3,3))/xj
      xfi(1,3)=(xf(1,2)*xf(2,3)-xf(1,3)*xf(2,2))/xj
      xfi(2,1)=(xf(2,3)*xf(3,1)-xf(2,1)*xf(3,3))/xj
      xfi(2,2)=(xf(1,1)*xf(3,3)-xf(1,3)*xf(3,1))/xj
      xfi(2,3)=(xf(1,3)*xf(2,1)-xf(1,1)*xf(2,3))/xj
      xfi(3,1)=(xf(2,1)*xf(3,2)-xf(2,2)*xf(3,1))/xj
      xfi(3,2)=(xf(1,2)*xf(3,1)-xf(1,1)*xf(3,2))/xj
      xfi(3,3)=(xf(1,1)*xf(2,2)-xf(1,2)*xf(2,1))/xj
!
      do a=1,nope
         do i=1,3
            g(i)=shp(i,a)
         enddo
!        gx = F^-T g  (gx_i = sum_j xfi(j,i) g_j)
         do i=1,3
            gx(i)=xfi(1,i)*g(1)+xfi(2,i)*g(2)+xfi(3,i)*g(3)
         enddo
!        h = F S g
         do i=1,3
            h(i)=fs(i,1)*g(1)+fs(i,2)*g(2)+fs(i,3)*g(3)
         enddo
!
         do k=1,3
!           dE_ij = 1/2 (F_ki g_j + F_kj g_i)
            do i=1,3
               do j=1,3
                  de(i,j)=0.5d0*(xf(k,i)*g(j)+xf(k,j)*g(i))
               enddo
            enddo
            trde=de(1,1)+de(2,2)+de(3,3)
            do i=1,3
               do j=1,3
                  ds(i,j)=2.d0*um*de(i,j)
               enddo
               ds(i,i)=ds(i,i)+al*trde
            enddo
!           tmp = F dS
            do i=1,3
               do j=1,3
                  tmp(i,j)=0.d0
                  do m=1,3
                     tmp(i,j)=tmp(i,j)+xf(i,m)*ds(m,j)
                  enddo
               enddo
            enddo
!           dt = -gx_k sigma + (e_k h^T + h e_k^T)/J + F dS F^T / J
            do i=1,3
               do j=1,3
                  dt(i,j)=-gx(k)*t(i,j)
                  do m=1,3
                     dt(i,j)=dt(i,j)+tmp(i,m)*xf(j,m)/xj
                  enddo
               enddo
            enddo
            do j=1,3
               dt(k,j)=dt(k,j)+h(j)/xj
               dt(j,k)=dt(j,k)+h(j)/xj
            enddo
!
            dot=0.d0
            do i=1,3
               do j=1,3
                  dot=dot+vs(i,j)*dt(i,j)
               enddo
            enddo
            dvm(3*(a-1)+k)=dot/vm
         enddo
      enddo
!
      return
      end
