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
!        q contains the nodal forces per element; initialisation of q
!
         if((iperturb(1).ge.2).or.((iperturb(1).le.0).and.(iout.lt.1))) 
     &      then
            do m1=1,nope
               do m2=0,mi(2)
                  q(m2,m1)=fn(m2,konl(m1))
               enddo
            enddo
         endif
!

         do jj=1,mint3d
            if(lakonl(4:5).eq.'8R') then
               xi=gauss3d1(1,jj)
               et=gauss3d1(2,jj)
               ze=gauss3d1(3,jj)
               weight=weight3d1(jj)
            elseif(lakonl(4:7).eq.'20RB') then
               if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
                  xi=gauss3d13(1,jj)
                  et=gauss3d13(2,jj)
                  ze=gauss3d13(3,jj)
                  weight=weight3d13(jj)
               else
                  call beamintscheme(lakonl,mint3d,ielprop(i),prop,
     &                 jj,xi,et,ze,weight)
               endif
            elseif((lakonl(4:4).eq.'8').or.
     &             (lakonl(4:6).eq.'20R'))
     &        then
               if(lakonl(7:8).ne.'LC') then
                  xi=gauss3d2(1,jj)
                  et=gauss3d2(2,jj)
                  ze=gauss3d2(3,jj)
                  weight=weight3d2(jj)
               else
                  kl=mod(jj,8)
                  if(kl.eq.0) kl=8
!
                  xi=gauss3d2(1,kl)
                  et=gauss3d2(2,kl)
                  ze=gauss3d2(3,kl)
                  weight=weight3d2(kl)
!
                  ki=mod(jj,4)
                  if(ki.eq.0) ki=4
!
                  if(kl.eq.1) then
                     ilayer=ilayer+1
                     if(ilayer.gt.1) then
                        do k=1,4
                           dlayer(k)=dlayer(k)+xlayer(ilayer-1,k)
                        enddo
                     endif
                  endif
                  ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &                 tlayer(ki)-1.d0
                  weight=weight*xlayer(ilayer,ki)/tlayer(ki)
!
!                 material and orientation
!
                  imat=ielmat(ilayer,i)
                  amat=matname(imat)
                  if(norien.gt.0) then
                     iorien=ielorien(ilayer,i)
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
!            elseif(lakonl(4:4).eq.'2') then
!               xi=gauss3d3(1,jj)
!               et=gauss3d3(2,jj)
!               ze=gauss3d3(3,jj)
!               weight=weight3d3(jj)
!            elseif(lakonl(4:5).eq.'10') then
!               xi=gauss3d5(1,jj)
!               et=gauss3d5(2,jj)
!               ze=gauss3d5(3,jj)
!               weight=weight3d5(jj)
            elseif(lakonl(4:4).eq.'4') then
               xi=gauss3d4(1,jj)
               et=gauss3d4(2,jj)
               ze=gauss3d4(3,jj)
               weight=weight3d4(jj)
            elseif(lakonl(4:5).eq.'15') then
               if(lakonl(7:8).ne.'LC') then
                  xi=gauss3d8(1,jj)
                  et=gauss3d8(2,jj)
                  ze=gauss3d8(3,jj)
                  weight=weight3d8(jj)
               else
                  kl=mod(jj,6)
                  if(kl.eq.0) kl=6
!
                  xi=gauss3d10(1,kl)
                  et=gauss3d10(2,kl)
                  ze=gauss3d10(3,kl)
                  weight=weight3d10(kl)
!
                  ki=mod(jj,3)
                  if(ki.eq.0) ki=3
!
                  if(kl.eq.1) then
                     ilayer=ilayer+1
                     if(ilayer.gt.1) then
                        do k=1,3
                           dlayer(k)=dlayer(k)+xlayer(ilayer-1,k)
                        enddo
                     endif
                  endif
                  ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &                 tlayer(ki)-1.d0
                  weight=weight*xlayer(ilayer,ki)/tlayer(ki)
!
!                 material and orientation
!
                  imat=ielmat(ilayer,i)
                  amat=matname(imat)
                  if(norien.gt.0) then
                     iorien=ielorien(ilayer,i)
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
            elseif(lakonl(4:4).eq.'6') then
               xi=gauss3d7(1,jj)
               et=gauss3d7(2,jj)
               ze=gauss3d7(3,jj)
               weight=weight3d7(jj)
            endif
!
c     Bernhardi start
            if(lakonl(1:5).eq.'C3D8R') then
               call shape8hr(xl,xsj,shp,gs,a)
            elseif(lakonl(1:5).eq.'C3D8I') then
               call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.20) then
c     Bernhardi end
               if(lakonl(7:7).eq.'A') then
                  call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
               elseif((lakonl(7:7).eq.'E').or.
     &                (lakonl(7:7).eq.'S')) then
                  call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
               else
                  call shape20h(xi,et,ze,xl,xsj,shp,iflag)
               endif
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
!           for frequency analysis or buckling with preload the
!           strains are calculated with respect to the deformed
!           configuration
!           for a linear iteration within a nonlinear increment:
!           the tangent matrix is calculated at strain at the end
!           of the previous increment
!
            if((iperturb(1).eq.1).or.(iperturb(1).eq.-1))then
               do m2=1,3
                  do m3=1,3
                     vokl(m2,m3)=0.d0
                  enddo
               enddo
!
               do m1=1,nope
                  do m2=1,3
                     do m3=1,3
                        vokl(m2,m3)=vokl(m2,m3)+
     &                       shp(m3,m1)*voldl(m2,m1)
                     enddo
                  enddo
               enddo
            endif
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
            if(iperturb(2).eq.1) then
!     
!                 Lagrangian strain
!     
               exx=exx+(vkl(1,1)**2+vkl(2,1)**2+vkl(3,1)**2)/2.d0
               eyy=eyy+(vkl(1,2)**2+vkl(2,2)**2+vkl(3,2)**2)/2.d0
               ezz=ezz+(vkl(1,3)**2+vkl(2,3)**2+vkl(3,3)**2)/2.d0
               exy=exy+vkl(1,1)*vkl(1,2)+vkl(2,1)*vkl(2,2)+
     &              vkl(3,1)*vkl(3,2)
               exz=exz+vkl(1,1)*vkl(1,3)+vkl(2,1)*vkl(2,3)+
     &              vkl(3,1)*vkl(3,3)
               eyz=eyz+vkl(1,2)*vkl(1,3)+vkl(2,2)*vkl(2,3)+
     &              vkl(3,2)*vkl(3,3)
!
!           for frequency analysis or buckling with preload the
!           strains are calculated with respect to the deformed
!           configuration
!
            elseif(iperturb(1).eq.1) then
               exx=exx+vokl(1,1)*vkl(1,1)+vokl(2,1)*vkl(2,1)+
     &              vokl(3,1)*vkl(3,1)
               eyy=eyy+vokl(1,2)*vkl(1,2)+vokl(2,2)*vkl(2,2)+
     &              vokl(3,2)*vkl(3,2)
               ezz=ezz+vokl(1,3)*vkl(1,3)+vokl(2,3)*vkl(2,3)+
     &              vokl(3,3)*vkl(3,3)
               exy=exy+vokl(1,1)*vkl(1,2)+vokl(1,2)*vkl(1,1)+
     &              vokl(2,1)*vkl(2,2)+vokl(2,2)*vkl(2,1)+
     &              vokl(3,1)*vkl(3,2)+vokl(3,2)*vkl(3,1)
               exz=exz+vokl(1,1)*vkl(1,3)+vokl(1,3)*vkl(1,1)+
     &              vokl(2,1)*vkl(2,3)+vokl(2,3)*vkl(2,1)+
     &              vokl(3,1)*vkl(3,3)+vokl(3,3)*vkl(3,1)
               eyz=eyz+vokl(1,2)*vkl(1,3)+vokl(1,3)*vkl(1,2)+
     &              vokl(2,2)*vkl(2,3)+vokl(2,3)*vkl(2,2)+
     &              vokl(3,2)*vkl(3,3)+vokl(3,3)*vkl(3,2)
            endif
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
            else
!
!              linear iteration within a nonlinear increment:
!
               eloc(1)=vokl(1,1)+
     &           (vokl(1,1)**2+vokl(2,1)**2+vokl(3,1)**2)/2.d0
               eloc(2)=vokl(2,2)+
     &           (vokl(1,2)**2+vokl(2,2)**2+vokl(3,2)**2)/2.d0
               eloc(3)=vokl(3,3)+
     &           (vokl(1,3)**2+vokl(2,3)**2+vokl(3,3)**2)/2.d0
               eloc(4)=(vokl(1,2)+vokl(2,1)+vokl(1,1)*vokl(1,2)+
     &              vokl(2,1)*vokl(2,2)+vokl(3,1)*vokl(3,2))/2.d0
               eloc(5)=(vokl(1,3)+vokl(3,1)+vokl(1,1)*vokl(1,3)+
     &              vokl(2,1)*vokl(2,3)+vokl(3,1)*vokl(3,3))/2.d0
               eloc(6)=(vokl(2,3)+vokl(3,2)+vokl(1,2)*vokl(1,3)+
     &              vokl(2,2)*vokl(2,3)+vokl(3,2)*vokl(3,3))/2.d0
            endif
!
!                 calculating the deformation gradient (needed to
!                 convert the element stiffness matrix from spatial
!                 coordinates to material coordinates
!                 deformation plasticity)
!
            if((kode.eq.-50).or.(kode.le.-100)) then
!
!                    calculating the deformation gradient
!
c     Bernhardi start
               xkl(1,1)=vkl(1,1)+1.0d0
               xkl(2,2)=vkl(2,2)+1.0d0
               xkl(3,3)=vkl(3,3)+1.0d0
c     Bernhardi end
               xkl(1,2)=vkl(1,2)
               xkl(1,3)=vkl(1,3)
               xkl(2,3)=vkl(2,3)
               xkl(2,1)=vkl(2,1)
               xkl(3,1)=vkl(3,1)
               xkl(3,2)=vkl(3,2)
!
!                    calculating the Jacobian
!
               vj=xkl(1,1)*(xkl(2,2)*xkl(3,3)-xkl(2,3)*xkl(3,2))
     &              -xkl(1,2)*(xkl(2,1)*xkl(3,3)-xkl(2,3)*xkl(3,1))
     &              +xkl(1,3)*(xkl(2,1)*xkl(3,2)-xkl(2,2)*xkl(3,1))
!
!              inversion of the deformation gradient (only for
!              deformation plasticity)
!
               if(kode.eq.-50) then
!
                  ckl(1,1)=(xkl(2,2)*xkl(3,3)-xkl(2,3)*xkl(3,2))/vj
                  ckl(2,2)=(xkl(1,1)*xkl(3,3)-xkl(1,3)*xkl(3,1))/vj
                  ckl(3,3)=(xkl(1,1)*xkl(2,2)-xkl(1,2)*xkl(2,1))/vj
                  ckl(1,2)=(xkl(1,3)*xkl(3,2)-xkl(1,2)*xkl(3,3))/vj
                  ckl(1,3)=(xkl(1,2)*xkl(2,3)-xkl(2,2)*xkl(1,3))/vj
                  ckl(2,3)=(xkl(2,1)*xkl(1,3)-xkl(1,1)*xkl(2,3))/vj
                  ckl(2,1)=(xkl(3,1)*xkl(2,3)-xkl(2,1)*xkl(3,3))/vj
                  ckl(3,1)=(xkl(2,1)*xkl(3,2)-xkl(2,2)*xkl(3,1))/vj
                  ckl(3,2)=(xkl(3,1)*xkl(1,2)-xkl(1,1)*xkl(3,2))/vj
!
!                 converting the Lagrangian strain into Eulerian
!                 strain (only for deformation plasticity)
!
                  cauchy=0
                  call str2mat(eloc,ckl,vj,cauchy)
               endif
!                        
            endif
!
!                 calculating fields for incremental plasticity
!
            if(kode.le.-100) then
!
!              calculating the deformation gradient at the
!              start of the increment
!
!              calculating the displacement gradient at the
!              start of the increment
!
               do m2=1,3
                  do m3=1,3
                     vikl(m2,m3)=0.d0
                  enddo
               enddo
!
               do m1=1,nope
                  do m2=1,3
                     do m3=1,3
                        vikl(m2,m3)=vikl(m2,m3)
     &                       +shp(m3,m1)*vini(m2,konl(m1))
                     enddo
                  enddo
               enddo
!
!              calculating the deformation gradient of the old
!              fields
!
               xikl(1,1)=vikl(1,1)+1.d0
               xikl(2,2)=vikl(2,2)+1.d0
               xikl(3,3)=vikl(3,3)+1.d0
               xikl(1,2)=vikl(1,2)
               xikl(1,3)=vikl(1,3)
               xikl(2,3)=vikl(2,3)
               xikl(2,1)=vikl(2,1)
               xikl(3,1)=vikl(3,1)
               xikl(3,2)=vikl(3,2)
!
!              calculating the Jacobian
!
               vij=xikl(1,1)*(xikl(2,2)*xikl(3,3)
     &              -xikl(2,3)*xikl(3,2))
     &              -xikl(1,2)*(xikl(2,1)*xikl(3,3)
     &              -xikl(2,3)*xikl(3,1))
     &              +xikl(1,3)*(xikl(2,1)*xikl(3,2)
     &              -xikl(2,2)*xikl(3,1))
!
!              stresses at the start of the increment
!
               do m1=1,6
                  stre(m1)=stiini(m1,jj,i)
               enddo
!
            endif
!
!                 prestress values
!
            if(iprestr.eq.1) then
               do kk=1,6
                  beta(kk)=-prestr(kk,jj,i)
               enddo
            else
               do kk=1,6
                  beta(kk)=0.d0
               enddo
            endif
!
            if(ithermal(1).ge.1) then
!
!              calculating the temperature difference in
!              the integration point
!
               t0l=0.d0
               t1l=0.d0
               if(ithermal(1).eq.1) then
                  if((lakonl(4:5).eq.'8 ').or.
     &               (lakonl(4:5).eq.'8I')) then
                     do i1=1,8
                        t0l=t0l+t0(konl(i1))/8.d0
                        t1l=t1l+t1(konl(i1))/8.d0
                     enddo
                  elseif(lakonl(4:6).eq.'20 ') then
                     nopered=20
                     call lintemp(t0,t1,konl,nopered,jj,t0l,t1l)
                  elseif(lakonl(4:6).eq.'10T') then
                     call linscal10(t0,konl,t0l,null,shp)
                     call linscal10(t1,konl,t1l,null,shp)
                  else
                     do i1=1,nope
                        t0l=t0l+shp(4,i1)*t0(konl(i1))
                        t1l=t1l+shp(4,i1)*t1(konl(i1))
                     enddo
                  endif
               elseif(ithermal(1).ge.2) then
                  if((lakonl(4:5).eq.'8 ').or.
     &               (lakonl(4:5).eq.'8I')) then
                     do i1=1,8
                        t0l=t0l+t0(konl(i1))/8.d0
                        t1l=t1l+vold(0,konl(i1))/8.d0
                     enddo
                  elseif(lakonl(4:6).eq.'20 ') then
                     nopered=20
                     call lintemp_th(t0,vold,konl,nopered,jj,t0l,t1l,mi)
                  elseif(lakonl(4:6).eq.'10T') then
                     call linscal10(t0,konl,t0l,null,shp)
                     call linscal10(vold,konl,t1l,mi(2),shp)
                  else
                     do i1=1,nope
                        t0l=t0l+shp(4,i1)*t0(konl(i1))
                        t1l=t1l+shp(4,i1)*vold(0,konl(i1))
                     enddo
                  endif
               endif
               tt=t1l-t0l
            endif
!
!                 calculating the coordinates of the integration point
!                 for material orientation purposes (for cylindrical
!                 coordinate systems)
!
            if((iorien.gt.0).or.(kode.le.-100)) then
               do j=1,3
                  pgauss(j)=0.d0
                  do i1=1,nope
                     pgauss(j)=pgauss(j)+shp(4,i1)*co(j,konl(i1))
                  enddo
               enddo
            endif
!
!                 material data; for linear elastic materials
!                 this includes the calculation of the stiffness
!                 matrix
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
            if((iperturb(1).eq.-1).and.(nlgeom_undo.eq.0)) then
!
!                    if the forced displacements were changed at
!                    the start of a nonlinear step, the nodal
!                    forces due do this displacements are 
!                    calculated in a purely linear way, and
!                    the first iteration is purely linear in order
!                    to allow the displacements to redistribute
!                    in a quasi-static way (only applies to
!                    quasi-static analyses (*STATIC))
               write(*,*) 'In manual Hooks law loop'
!
               eloc(1)=exx-vokl(1,1)
               eloc(2)=eyy-vokl(2,2)
               eloc(3)=ezz-vokl(3,3)
               eloc(4)=exy-(vokl(1,2)+vokl(2,1))
               eloc(5)=exz-(vokl(1,3)+vokl(3,1))
               eloc(6)=eyz-(vokl(2,3)+vokl(3,2))
!
               if(mattyp.eq.1) then
                  e=elas(1)
                  un=elas(2)
                  um=e/(1.d0+un)
                  al=un*um/(1.d0-2.d0*un)
                  um=um/2.d0
                  am1=al*(eloc(1)+eloc(2)+eloc(3))
                  stre(1)=am1+2.d0*um*eloc(1)
                  stre(2)=am1+2.d0*um*eloc(2)
                  stre(3)=am1+2.d0*um*eloc(3)
                  stre(4)=um*eloc(4)
                  stre(5)=um*eloc(5)
                  stre(6)=um*eloc(6)
               elseif(mattyp.eq.2) then
                  stre(1)=eloc(1)*elas(1)+eloc(2)*elas(2)
     &                 +eloc(3)*elas(4)
                  stre(2)=eloc(1)*elas(2)+eloc(2)*elas(3)
     &                 +eloc(3)*elas(5)
                  stre(3)=eloc(1)*elas(4)+eloc(2)*elas(5)
     &                 +eloc(3)*elas(6)
                  stre(4)=eloc(4)*elas(7)
                  stre(5)=eloc(5)*elas(8)
                  stre(6)=eloc(6)*elas(9)
               elseif(mattyp.eq.3) then
                  stre(1)=eloc(1)*elas(1)+eloc(2)*elas(2)+
     &                 eloc(3)*elas(4)+eloc(4)*elas(7)+
     &                 eloc(5)*elas(11)+eloc(6)*elas(16)
                  stre(2)=eloc(1)*elas(2)+eloc(2)*elas(3)+
     &                 eloc(3)*elas(5)+eloc(4)*elas(8)+
     &                 eloc(5)*elas(12)+eloc(6)*elas(17)
                  stre(3)=eloc(1)*elas(4)+eloc(2)*elas(5)+
     &                 eloc(3)*elas(6)+eloc(4)*elas(9)+
     &                 eloc(5)*elas(13)+eloc(6)*elas(18)
                  stre(4)=eloc(1)*elas(7)+eloc(2)*elas(8)+
     &                 eloc(3)*elas(9)+eloc(4)*elas(10)+
     &                 eloc(5)*elas(14)+eloc(6)*elas(19)
                  stre(5)=eloc(1)*elas(11)+eloc(2)*elas(12)+
     &                 eloc(3)*elas(13)+eloc(4)*elas(14)+
     &                 eloc(5)*elas(15)+eloc(6)*elas(20)
                  stre(6)=eloc(1)*elas(16)+eloc(2)*elas(17)+
     &                 eloc(3)*elas(18)+eloc(4)*elas(19)+
     &                 eloc(5)*elas(20)+eloc(6)*elas(21)
               endif
            endif
! 
!           updating the internal energy and mechanical strain
!           for user materials (kode<=-100) the mechanical strain has to
!           be updated at the end of each increment (also if no output
!           is requested), since it is input to the umat routine
!
            if((iout.gt.0).or.(iout.eq.-2).or.(kode.le.-100).or.
     &         ((nmethod.eq.4).and.
     &          ((iperturb(1).gt.1).and.(nlgeom_undo.eq.0)).and.
     &          (ithermal(1).le.1))) then
               if(ithermal(1).eq.0) then
                  do m1=1,6
                     eth(m1)=0.d0
                  enddo
               endif
               if(nener.eq.1) then
                  ener(jj,i)=enerini(jj,i)+
     &                 ((eloc(1)-eth(1)-emeini(1,jj,i))*
     &                  (stre(1)+stiini(1,jj,i))+
     &                  (eloc(2)-eth(2)-emeini(2,jj,i))*
     &                  (stre(2)+stiini(2,jj,i))+
     &                  (eloc(3)-eth(3)-emeini(3,jj,i))*
     &                  (stre(3)+stiini(3,jj,i)))/2.d0+
     &         (eloc(4)-eth(4)-emeini(4,jj,i))*(stre(4)+stiini(4,jj,i))+
     &         (eloc(5)-eth(5)-emeini(5,jj,i))*(stre(5)+stiini(5,jj,i))+
     &         (eloc(6)-eth(6)-emeini(6,jj,i))*(stre(6)+stiini(6,jj,i))
               endif
               eme(1,jj,i)=eloc(1)-eth(1)
               eme(2,jj,i)=eloc(2)-eth(2)
               eme(3,jj,i)=eloc(3)-eth(3)
               eme(4,jj,i)=eloc(4)-eth(4)
               eme(5,jj,i)=eloc(5)-eth(5)
               eme(6,jj,i)=eloc(6)-eth(6)
            endif
!
            if((iout.gt.0).or.(iout.eq.-2).or.(kode.le.-100)) then
!
               eei(1,jj,i)=eloc(1)
               eei(2,jj,i)=eloc(2)
               eei(3,jj,i)=eloc(3)
               eei(4,jj,i)=eloc(4)
               eei(5,jj,i)=eloc(5)
               eei(6,jj,i)=eloc(6)
            endif
!
!     updating the kinetic energy
!
            if(ikin.eq.1) then
!               
               call materialdata_rho(rhcon,nrhcon,imat,rho,t1l,
     &              ntmat_,ithermal)
               do m1=1,3
                  vel(m1,1)=0.d0
                  do i1= 1,nope
                     vel(m1,1)=vel(m1,1)+shp(4,i1)*veoldl(m1,i1)
                  enddo
               enddo
               ener(jj,i+ne)=rho*(vel(1,1)*vel(1,1)+
     &              vel(2,1)*vel(2,1)+ vel(3,1)*vel(3,1))/2.d0
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
!                 calculation of the nodal forces
!
            if(calcul_fn.eq.1)then
!
!                    calculating fn using skl
!
               do m1=1,nope
                  do m2=1,3
!
!                          linear elastic part
!                           
                     do m3=1,3
                        fn(m2,konl(m1))=fn(m2,konl(m1))+
     &                       xsj*skl(m2,m3)*shp(m3,m1)*weight
                     enddo
!
!                          nonlinear geometric part
!
                     if((iperturb(2).eq.1).and.(nlgeom_undo.eq.0)) then
                        do m3=1,3
                           do m4=1,3
                              fn(m2,konl(m1))=fn(m2,konl(m1))+
     &                             xsj*skl(m4,m3)*weight*
     &                             (vkl(m2,m4)*shp(m3,m1)+
     &                             vkl(m2,m3)*shp(m4,m1))/2.d0
                           enddo
                        enddo
                     endif
!
                  enddo
               enddo
c     Bernhardi start
               if(lakonl(1:5).eq.'C3D8R') then
                  call hgforce (fn,elas,a,gs,vl,mi,konl)
               endif
c     Bernhardi end
            endif
!
!           calculation of the Cauchy stresses
!
            if((calcul_cauchy.eq.1).and.(nlgeom_undo.eq.0)) then

               write(*,*), "Checking for Bernhardi start..."
!
!              changing the displacement gradients into
!              deformation gradients
!
c               if(kode.ne.-50) then
               if((kode.ne.-50).and.(kode.gt.-100)) then
                  write(*,*), "I am in kode.neq -50"
c     Bernhardi start
                  xkl(1,1)=vkl(1,1)+1.0d0
                  xkl(2,2)=vkl(2,2)+1.0d0
                  xkl(3,3)=vkl(3,3)+1.0d0
c     Bernhardi end   
                  xkl(1,2)=vkl(1,2)
                  xkl(1,3)=vkl(1,3)
                  xkl(2,3)=vkl(2,3)
                  xkl(2,1)=vkl(2,1)
                  xkl(3,1)=vkl(3,1)
                  xkl(3,2)=vkl(3,2)
!
                  vj=xkl(1,1)*(xkl(2,2)*xkl(3,3)-xkl(2,3)*xkl(3,2))
     &                 -xkl(1,2)*(xkl(2,1)*xkl(3,3)-xkl(2,3)*xkl(3,1))
     &                 +xkl(1,3)*(xkl(2,1)*xkl(3,2)-xkl(2,2)*xkl(3,1))
               endif
!
               do m1=1,3
                  do m2=1,m1
                     ckl(m1,m2)=0.d0
                     do m3=1,3
                        do m4=1,3
                           ckl(m1,m2)=ckl(m1,m2)+
     &                          skl(m3,m4)*xkl(m1,m3)*xkl(m2,m4)
                        enddo
                     enddo
                     ckl(m1,m2)=ckl(m1,m2)/vj
                  enddo
               enddo
!
               stx(1,jj,i)=ckl(1,1) 
               stx(2,jj,i)=ckl(2,2) 
               stx(3,jj,i)=ckl(3,3)
               stx(4,jj,i)=ckl(2,1)
               stx(5,jj,i)=ckl(3,1)
               stx(6,jj,i)=ckl(3,2)
            endif
! --- p-norm accumulation (AFTER the Cauchy block) ------
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
         enddo  ! <--- end of jj=1, mint3d
!
!        q contains the contributions to the nodal force in the nodes
!        belonging to the element at stake from other elements (elements
!        already treated). These contributions have to be
!        subtracted to get the contributions attributable to the element
!        at stake only
!
         if(calcul_qa.eq.1) then
            do m1=1,nope
               do m2=1,3
                  qa(1)=qa(1)+dabs(fn(m2,konl(m1))-q(m2,m1))
               enddo
            enddo
            nal=nal+3*nope
         endif
      enddo
!

!

      qa(3) = g_sump
      !qa(4) = g_vol

! ------------------------
      return
      end
