!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine stresspnorm(co,kon,ipkon,lakon,ne,v,
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
!

!     calculates the relaxed P-norm-aggregated stress based on 
!     Duysinx, P and Sigmund O, New developments in handling stress
!     constraints in optimal material distribution
!     7th AIAA symposium on multidisciplinary analysis and optimization
      
      implicit none
!
      integer cauchy
!
      character*8 lakon(*),lakonl
      character*80 amat,matname(*)
!
      
      integer kon(*),konl(26),nea,neb,mi(*),mint2d,nopes,
     &  nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),
     &  ielorien(mi(3),*),ntmat_,ipkon(*),ne0,iflag,null,kscale,
     &  istep,iinc,mt,ne,mattyp,ithermal(2),iprestr,i,j,k,m1,m2,jj,
     &  i1,m3,m4,kk,nener,indexe,nope,norien,iperturb(*),iout,
     &  nal,icmd,ihyper,nmethod,kode,imat,mint3d,iorien,ielas,
     &  istiff,ncmat_,nstate_,ikin,ilayer,nlayer,ki,kl,ielprop(*),
     &  nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,calcul_fn,
     &  calcul_cauchy,calcul_qa,nopered,mortar,jfaces,igauss,
     &  istrainfree,nlgeom_undo,list,ilist(*),m
!
      real*8 co(3,*),v(0:mi(2),*),shp(4,26),stiini(6,mi(1),*),
     &  stx(6,mi(1),*),xl(3,26),vl(0:mi(2),26),stre(6),prop(*),
     &  elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),xs2(3,7),
     &  alcon(0:6,ntmat_,*),vini(0:mi(2),*),thickness,
     &  alzero(*),orab(7,*),elas(21),rho,fn(0:mi(2),*),
     &  fnl(3,10),skl(3,3),beta(6),q(0:mi(2),26),xl2(3,8),
     &  vkl(0:3,3),t0(*),t1(*),prestr(6,mi(1),*),eme(6,mi(1),*),
     &  ckl(3,3),vold(0:mi(2),*),eloc(9),veold(0:mi(2),*),
     &  springarea(2,*),elconloc(21),eth(6),xkl(3,3),voldl(0:mi(2),26),
     &  xikl(3,3),ener(mi(1),*),emec(6),eei(6,mi(1),*),enerini(mi(1),*),
     &  emec0(6),vel(1:3,26),veoldl(0:mi(2),26),xsj2(3),shp2(7,8),
     &  e,un,al,um,am1,xi,et,ze,tt,exx,eyy,ezz,exy,exz,eyz,
     &  xsj,qa(*),vj,t0l,t1l,dtime,weight,pgauss(3),vij,time,ttime,
     &  plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &  xstiff(27,mi(1),*),xstate(nstate_,mi(1),*),plconloc(802),
     &  vokl(3,3),xstateini(nstate_,mi(1),*),vikl(3,3),
     &  gs(8,4),a,reltime,tlayer(4),dlayer(4),xlayer(mi(3),4),
     &  thicke(mi(3),*),emeini(6,mi(1),*),clearini(3,9,*),
     &  pslavsurf(3,*),pmastsurf(6,*)
!     
!     Element density and penalization vector
      real*8 design(*), penal

!     p-norm variables 
      real*8 sx, sy, sz, txy, txz, tyz, vm, vm2, wgt
      real*8 g_sump, g_vol, pexp
      real*8 rho_e, rho_min, rho_eff, rho_p, eps_relax, sig0, phi
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
      null=0
!
      mt=mi(2)+1
      nal=0
      qa(3)=-1.d0
      qa(4)=0.d0

      ! --- INIT global p-norm accumulators ---
      g_sump = 0.d0
      g_vol  = 0.d0
      !pexp = 4.d0     ! choose your global p




!     Loop over all elements starting from nea...neb
      do m=nea,neb
!
         if(list.eq.1) then
            i=ilist(m)
         else
            i=m
         endif
!
!
         lakonl=lakon(i)
!
!

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

!        Set number of nodes for C3D4 element
         if(lakonl(4:4).eq.'4') then
            nope=4
         endif

!        Set number integration points for C3D4
         if(lakonl(4:4).eq.'4') then
            mint3d=1  ! One integration point per each C3D4 element
         endif

!        Gather nodal displacements
         do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl(j))
               vl(k,j)=v(k,konl(j))
               voldl(k,j)=vold(k,konl(j))
            enddo
         enddo
!

!     Begin loop over all integrations points per element
         do jj=1,mint3d
            ! Get integration weight for this element
            if(lakonl(4:4).eq.'4') then
               xi=gauss3d4(1,jj)
               et=gauss3d4(2,jj)
               ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
            endif
            ! Compute shape func values
            if(nope.eq.4) then
               call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            endif

!
!                 vkl(m2,m3) contains the derivative of the m2-
!                 component of the displacement with respect to
!                 direction m3
!
            do m2=1,3
               do m3=1,3
                  vkl(m2,m3)=0.d0
               enddo
            enddo
!
            do m1=1,nope
               do m2=1,3
                  do m3=1,3
                     vkl(m2,m3)=vkl(m2,m3)+shp(m3,m1)*vl(m2,m1)
                  enddo
c                  write(*,*) 'vnoeie',i,konl(m1),(vkl(m2,k),k=1,3)
               enddo
            enddo
!
!
            kode=nelcon(1,imat)
!
!           calculating the strain
!
!           attention! exy,exz and eyz are engineering strains!
!           Small-strain (engineering shears)
            exx=vkl(1,1)
            eyy=vkl(2,2)
            ezz=vkl(3,3)
            exy=vkl(1,2)+vkl(2,1)
            exz=vkl(1,3)+vkl(3,1)
            eyz=vkl(2,3)+vkl(3,2)
!
!
!              storing the local strains
!              Store tenorial components
            if(iperturb(1).ne.-1) then
               eloc(1)=exx
               eloc(2)=eyy
               eloc(3)=ezz
               eloc(4)=exy/2.d0
               eloc(5)=exz/2.d0
               eloc(6)=eyz/2.d0
            endif
!
!
!           material data; for linear elastic materials
!           this includes the calculation of the stiffness
!           matrix
!
            istiff=0
!     
            call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,
     &           nalcon,imat,amat,iorien,pgauss,orab,ntmat_,
     &           elas,rho,i,ithermal,alzero,mattyp,t0l,t1l,ihyper,
     &           istiff,elconloc,eth,kode,plicon,nplicon,
     &           plkcon,nplkcon,npmat_,plconloc,mi(1),dtime,jj,
     &           xstiff,ncmat_)
!
!           determining the mechanical strain
!
            if(ithermal(1).ne.0) then
               do m1=1,6
                  emec(m1)=eloc(m1)-eth(m1)
               enddo
            else
               do m1=1,6
                  emec(m1)=eloc(m1)
               enddo
            endif
            if(kode.le.-100) then
               do m1=1,6
                  emec0(m1)=emeini(m1,jj,i)
               enddo
            endif
!
!           subtracting the initial strains
!
            if(iprestr.eq.2) then
               if(istrainfree==0) then
                  do m1=1,6
                     emec(m1)=emec(m1)-prestr(m1,jj,i)
                  enddo
               else
                  do m1=1,6
                     prestr(m1,jj,i)=emec(m1)
                     emeini(m1,jj,i)=emec(m1)
                     eme(m1,jj,i)=emec(m1)
                     emec(m1)=0.d0
                  enddo
               endif
            endif

! --- SIMP penalization for linear isotropic (mattyp == 1) ---
            if (mattyp .eq. 1) then
!              write(*,*) 'Scaling the C matrix!'
              rho_e   = design(i)
              if (rho_e .lt. 0.d0) rho_e = 0.d0
              if (rho_e .gt. 1.d0) rho_e = 1.d0
              rho_eff = max(rho_e, rho_min)
!              Based on the definition of effective von Misses 
!              stress in Duysinx and Sigmnd, the penalized
!              rho cancels out in the stress term with
!              relaxation. Pass rho_p = 1 below
!             rho_p   = rho_eff**penal
               rho_p = 1.d0
!              Note: Since we are not penalizing the stress here
!              von Misses must be penalized in ccx_2.15.c 
!              when writing to .vtu file
              ! Pass rho_p to linel.f through mechmodel.f
!              write(*,*) 'Done!'
            endif
!
!           calculating the local stiffness and stress
!           Constitutive law
            nlgeom_undo=0
            call mechmodel(elconloc,elas,emec,kode,emec0,ithermal,
     &           icmd,beta,stre,xkl,ckl,vj,xikl,vij,
     &           plconloc,xstate,xstateini,ielas,
     &           amat,t1l,dtime,time,ttime,i,jj,nstate_,mi(1),
     &           iorien,pgauss,orab,eloc,mattyp,qa(3),istep,iinc,
     &           ipkon,nmethod,iperturb,qa(4),nlgeom_undo)
!
            if(((nmethod.ne.4).or.(iperturb(1).ne.0)).and.
     &         (nmethod.ne.5).and.(icmd.ne.3)) then
!              No scaling of stress matrix since effective von Misses
!              with relaxation only depends on xstiff for rhp_p = 1
              do m1=1,21
                  xstiff(m1,jj,i)=elas(m1)
               enddo
            endif
!
!
            if((iout.gt.0).or.(iout.eq.-2).or.(kode.le.-100)) then
               eei(1,jj,i)=eloc(1)
               eei(2,jj,i)=eloc(2)
               eei(3,jj,i)=eloc(3)
               eei(4,jj,i)=eloc(4)
               eei(5,jj,i)=eloc(5)
               eei(6,jj,i)=eloc(6)
            endif

!           Write stress into stx (integration-point storage)
            skl(1,1)=stre(1)
            skl(2,2)=stre(2)
            skl(3,3)=stre(3)
            skl(2,1)=stre(4)
            skl(3,1)=stre(5)
            skl(3,2)=stre(6)
!
            stx(1,jj,i)=skl(1,1)
            stx(2,jj,i)=skl(2,2) 
            stx(3,jj,i)=skl(3,3)
            stx(4,jj,i)=skl(2,1)
            stx(5,jj,i)=skl(3,1)
            stx(6,jj,i)=skl(3,2)
!
            skl(1,2)=skl(2,1)
            skl(1,3)=skl(3,1)
            skl(2,3)=skl(3,2)
!
!
!           calculation of the Cauchy stresses (skip for liner CalTop)
            if((calcul_cauchy.eq.1).and.(nlgeom_undo.eq.0)) then

               write(*,*), "Skipping Cauchy stress eval"
           endif
!--------------------------------------------------------------!
!                       BEGIN P-NORM EVAL                          
!--------------------------------------------------------------!
! --- Read element stress values
            sx  = stx(1,jj,i) 
            sy  = stx(2,jj,i)
            sz  = stx(3,jj,i)
            txy = stx(4,jj,i)
            txz = stx(5,jj,i)
            tyz = stx(6,jj,i)

            !print *, 'sx = ', sx
            !print *, 'sy = ', sy

!  --- von Mises (without penalization)
            vm2 = (sx-sy)*(sx-sy) + (sy-sz)*(sy-sz) + (sz-sx)*(sz-sx)
            vm2 = 0.5d0*vm2 + 3.d0*(txy*txy + txz*txz + tyz*tyz)
            vm  = dsqrt(vm2)

!            write(*,*), 'Currrent vm:', vm

!  --- filtered design alread in [0,1] (clamp defenseively)  ---
            rho_e = design(i)
            ! (optional clamp, safe if design may drift)
            if (rho_e .lt. 0.d0) rho_e = 0.d0
            if (rho_e .gt. 1.d0) rho_e = 1.d0

            rho_eff = max(rho_e, rho_min)
            rho_p = rho_eff**penal

! --- Duysinx-Sigmund effective von Misses stress measure for an element
            phi = (vm/ (sig0)) + eps_relax - (eps_relax/rho_eff)
! --- Consider effective von Misses stress (phi) > 0 only
            if (phi .lt. 0.d0) phi = 0.d0

!           write(*,*), 'Phi: ', phi

! --- With effective von Misses stress calculated, raise to pexp
! --- and sum over all elements
! --- P-norm aggregation for this element
            g_sump = g_sump + (phi**pexp)
!
         enddo  ! <--- end of integration over element Gauss points

      enddo   ! <--- end of loop over all elements
!

!

      qa(3) = g_sump


! ------------------------
      return
      end
