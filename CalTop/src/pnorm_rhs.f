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
       subroutine pnorm_rhs(co,kon,ipkon,lakon,ne,v,
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
      real*8 rho_e, rho_min, rho_eff, rho_p, eps_relax, sig0, phi, alpha

!     variables for accumulating adjoint rhs
      real*8 g(6), ptv(6), invvm, fac
      real*8 B(6,12), rhs_loc(12)
      real*8 rhs(0:mi(2),*)
      real*8 C(6,6)

!     work arrays for Me*u_e path (C3D4 block)
      integer ia, ib
      real*8 Bten(6,12)
      real*8 Vec(6,6)
      real*8 uel(12)
      real*8 s(6), t(6)
      real*8 coeff
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
     &  list,ilist,design, penal, sig0, eps_relax, rho_min, pexp, alpha
!
      intent(inout) nal,qa,fn,xstiff,ener,eme,eei,stx,ielmat,prestr,
     &  emeini, rhs
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

!     
! --- Begin loop over all elements starting from nea
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


!        Set orientation flags
         if(lakonl(7:8).ne.'LC') then
!
            imat=ielmat(1,i)
            amat=matname(imat)

            if(norien.gt.0) then
               iorien=ielorien(1,i)
            else
               iorien=0
            endif
!
            if(nelcon(1,imat).lt.0) then
               ihyper=1
            else
               ihyper=0
            endif
         endif


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
!           vkl(m2,m3) contains the derivative of the m2-
!           component of the displacement with respect to
!           direction m3

!           Initialize vkl to zero
            do m2=1,3
               do m3=1,3
                  vkl(m2,m3)=0.d0
               enddo
            enddo

!           Compute vkl
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
!           storing the local strains
!           Store tenorial components
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

! ---       SIMP penalization for linear isotropic (mattyp == 1) ---
            if (mattyp .eq. 1) then
              rho_e   = design(i)
              if (rho_e .lt. 0.d0) rho_e = 0.d0
              if (rho_e .gt. 1.d0) rho_e = 1.d0
              rho_eff = max(rho_e, rho_min)
              rho_p   = rho_eff**penal

! ---         Do not scale using rho_p anywhere since the adjoint RHS 
!             depends on the rho_p = 1 state
!              TODO: Remove above conditional
              rho_p = 1.d0
            endif
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

!               Build full 6X6 C (voigt, tensoral shear)
                if (mattyp.eq.1) then
!               Isotropic
                  e  = elas(1)
                  un = elas(2)

                  um = e/(2.d0*(1.d0+un))                  ! G
                  al = un*e/((1.d0+un)*(1.d0-2.d0*un))     ! lambda

!                 Pack upper-triangular Voigt (1..21) with tensorial shear:
!                 1:C11  2:C12  3:C22  4:C13  5:C23  6:C33
!                 7:C14  8:C24  9:C34 10:C44 11:C15 12:C25
!                 13:C35 14:C45 15:C55 16:C16 17:C26 18:C36
!                 19:C46 20:C56 21:C66
                  xstiff( 1,jj,i)=(al+2.d0*um)  ! C11
                  xstiff( 2,jj,i)= al           ! C12
                  xstiff( 3,jj,i)= (al+2.d0*um)  ! C22
                  xstiff( 4,jj,i)= al           ! C13
                  xstiff( 5,jj,i)= al           ! C23
                  xstiff( 6,jj,i)= (al+2.d0*um)  ! C33

                  xstiff( 7,jj,i)=0.d0                ! C14
                  xstiff( 8,jj,i)=0.d0                ! C24
                  xstiff( 9,jj,i)=0.d0                ! C34
                  xstiff(10,jj,i)= um           ! C44 (τ12/ε12_tensorial)
                  xstiff(11,jj,i)=0.d0                ! C15
                  xstiff(12,jj,i)=0.d0                ! C25
                  xstiff(13,jj,i)=0.d0                ! C35
                  xstiff(14,jj,i)=0.d0                ! C45
                  xstiff(15,jj,i)= um           ! C55 (τ13/ε13_tensorial)
                  xstiff(16,jj,i)=0.d0                ! C16
                  xstiff(17,jj,i)=0.d0                ! C26
                  xstiff(18,jj,i)=0.d0                ! C36
                  xstiff(19,jj,i)=0.d0                ! C46
                  xstiff(20,jj,i)=0.d0                ! C56
                  xstiff(21,jj,i)= um           ! C66 (τ23/ε23_tensorial)
               endif
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

!           calculation of the Cauchy stresses (skip for linear CalTop)
            if((calcul_cauchy.eq.1).and.(nlgeom_undo.eq.0)) then
            
               write(*,*), "Skipping Cauchy stress eval"
            endif
!--------------------------------------------------------------!
!                       BEGIN P-NORM RHS EVAL                          
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

!  --- von Mises 
            vm2 = (sx-sy)*(sx-sy) + (sy-sz)*(sy-sz) + (sz-sx)*(sz-sx)
            vm2 = 0.5d0*vm2 + 3.d0*(txy*txy + txz*txz + tyz*tyz)
            vm  = dsqrt(vm2)

!  --- filtered design alread in [0,1] (clamp defenseively)  ---
            rho_e = design(i)
            ! (optional clamp, safe if design may drift)
            if (rho_e .lt. 0.d0) rho_e = 0.d0
            if (rho_e .gt. 1.d0) rho_e = 1.d0

            rho_eff = max(rho_e, rho_min)
            rho_p = rho_eff**penal

! --- Duysinx-Sigmund effective von Misses stress measure for this element
            phi = vm/ (sig0) + eps_relax - eps_relax/rho_eff
            if (phi .lt. 0.d0) phi = 0.d0

! --- With effective von Misses stress calculated, raise to pexp
! --- and sum over all elements
!            g_sump = g_sump + (phi**pexp)
            !g_vol  = g_vol  + wgt  <-- valid only for p-mean


c----- RHS accumulation (minimal: C3D4 only) --------------------------
            if (nope.eq.4) then

c           skip if vm or phi yields zero gradient
              if (vm.gt.0.d0) then

!  ---        evaluate 
                invvm = 1.d0/(sig0 * vm)
                coeff   = (phi**(pexp-1)) * invvm

! ---- unpack xstiff(:,jj,i) -> full symmetric C(6,6) in tensorial Voigt ----
               
               C = 0.d0
               C(1,1)=xstiff( 1,jj,i);
               C(2,1)=xstiff( 2,jj,i); 
               C(2,2)=xstiff( 3,jj,i);
               C(3,1)=xstiff( 4,jj,i); 
               C(3,2)=xstiff( 5,jj,i); 
               C(3,3)=xstiff( 6,jj,i);
               C(4,1)=xstiff( 7,jj,i); 
               C(4,2)=xstiff( 8,jj,i); 
               C(4,3)=xstiff( 9,jj,i); 
               C(4,4)=xstiff(10,jj,i);
               C(5,1)=xstiff(11,jj,i);
               C(5,2)=xstiff(12,jj,i); 
               C(5,3)=xstiff(13,jj,i); 
               C(5,4)=xstiff(14,jj,i); 
               C(5,5)=xstiff(15,jj,i);
               C(6,1)=xstiff(16,jj,i); 
               C(6,2)=xstiff(17,jj,i); 
               C(6,3)=xstiff(18,jj,i); 
               C(6,4)=xstiff(19,jj,i); 
               C(6,5)=xstiff(20,jj,i); 
               C(6,6)=xstiff(21,jj,i);
! mirror upper triangle
               C(1,2)=C(2,1); 
               C(1,3)=C(3,1); 
               C(1,4)=C(4,1); 
               C(1,5)=C(5,1); 
               C(1,6)=C(6,1);
               C(2,3)=C(3,2);
               C(2,4)=C(4,2); 
               C(2,5)=C(5,2); 
               C(2,6)=C(6,2);
               C(3,4)=C(4,3);
               C(3,5)=C(5,3);
               C(3,6)=C(6,3);
               C(4,5)=C(5,4);
               C(4,6)=C(6,4);
               C(5,6)=C(6,5);


c             Build B for C3D4 (engineering shear), from shp(derivs)
c             shp(1,j)=dNj/dx, shp(2,j)=dNj/dy, shp(3,j)=dNj/dz
               do m1=1,12
                  B(1,m1)=0.d0; B(2,m1)=0.d0; B(3,m1)=0.d0
                  B(4,m1)=0.d0; B(5,m1)=0.d0; B(6,m1)=0.d0
               enddo

               do m1=1,4
                  m2 = 3*(m1-1)
                  B(1,m2+1) = shp(1,m1)
                  B(2,m2+2) = shp(2,m1)
                  B(3,m2+3) = shp(3,m1)
                  B(4,m2+1) = shp(2,m1)   ! exy engineering
                  B(4,m2+2) = shp(1,m1)
                  B(6,m2+1) = shp(3,m1)   ! exz engineering
                  B(6,m2+3) = shp(1,m1)
                  B(5,m2+2) = shp(3,m1)   ! eyz engineering
                  B(5,m2+3) = shp(2,m1)
               enddo

! ---          Convert to tensorial shear: Bten = R * Beng, R=diag(1,1,1,1/2,1/2,1/2)
               do m1=1,12
                  Bten(1,m1)=B(1,m1)
                  Bten(2,m1)=B(2,m1)
                  Bten(3,m1)=B(3,m1)
                  Bten(4,m1)=B(4,m1)
                  Bten(5,m1)=B(5,m1)
                  Bten(6,m1)=B(6,m1)
               enddo
               !0.5d0*
!  ---         von-Misses selecter Vec in tensorial Voigt (3D)
               Vec = 0.d0
               Vec(1,1)=1.d0
               Vec(1,2)=-0.5d0
               Vec(1,3)=-0.5d0
               Vec(2,1)=-0.5d0
               Vec(2,2)=1.d0
               Vec(2,3)=-0.5d0
               Vec(3,1)=-0.5d0
               Vec(3,2)=-0.5d0
               Vec(3,3)=1.d0
               Vec(4,4)=3.d0
               Vec(5,5)=3.d0
               Vec(6,6)=3.d0

! ---          Gather u_e (12 X1) from vl(:,:)
               do m1=1,4
                  uel(3*(m1-1)+1) = vl(1,m1)
                  uel(3*(m1-1)+2) = vl(2,m1)
                  uel(3*(m1-1)+3) = vl(3,m1)
               enddo
! ---         s = C * (Bten * u_e) (tensorial stress)
               do ia=1,6
                  s(ia)=0.d0
                  do ib=1,12
                     s(ia)=s(ia)+C(ia,1)*Bten(1,ib)*uel(ib)
     &               +C(ia,2)*Bten(2,ib)*uel(ib)
     &               +C(ia,3)*Bten(3,ib)*uel(ib)
     &               +C(ia,4)*Bten(4,ib)*uel(ib)
     &               +C(ia,5)*Bten(5,ib)*uel(ib)
     &               +C(ia,6)*Bten(6,ib)*uel(ib)
                  enddo
               enddo

! ---          t = V * s
               do ia=1,6
                  t(ia)=0.d0
                  do ib=1,6
                     t(ia)=t(ia)+Vec(ia,ib)*s(ib)
                  enddo
               enddo

!  ---         p  = C * t 
               do ia=1,6
                  ptv(ia)=0.d0
                  do ib=1,6
                     ptv(ia)=ptv(ia)+C(ia,ib)*t(ib)
                  enddo
               enddo

! ---          Build Me0 u_e
               do m1=1,12
                  rhs_loc(m1)=0.d0
                  do ia = 1,6
                     rhs_loc(m1)=rhs_loc(m1) + Bten(ia,m1)*ptv(ia)
                  enddo
               enddo

! ---          Multiply with scalar coefficient
               do m1=1,12
                  rhs_loc(m1) = coeff * rhs_loc(m1)
!                  rhs_loc(m1) = coeff 
               enddo

! ---          Scatter to global rhs(1..3, node)
                do m1=1,4   ! Loop over all DOFs.
                  rhs(1,konl(m1)) = rhs(1,konl(m1))
     &                             + rhs_loc(3*(m1-1)+1)
                  rhs(2,konl(m1)) = rhs(2,konl(m1)) 
     &                             + rhs_loc(3*(m1-1)+2)
                  rhs(3,konl(m1)) = rhs(3,konl(m1)) 
     &                             + rhs_loc(3*(m1-1)+3)
                enddo
              endif
            endif ! End if nope .eq. 4 condition
c----------------------------------------------------------------------

!
         enddo  ! <--- end of integration over element Gauss points
      enddo ! <--- end of loop over all elements

! ------------------------
      return
      end
