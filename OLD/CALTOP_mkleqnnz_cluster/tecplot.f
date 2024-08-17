
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
      subroutine tecplot(set,nset,istartset,iendset,ialset,nprint,
     &  prlab,prset,v,t1,fn,ipkon,lakon,stx,eei,xstate,ener,
     &  mi,nstate_,ithermal,co,kon,qfx,ttime,trab,inotr,ntrans,
     &  orab,ielorien,norien,nk,ne,inum,filab,vold,ikin,ielmat,thicke,
     &  eme,islavsurf,mortar,time,ielprop,prop,veold,orname,
     &  nelemload,nload,sideload,xload)

      implicit none

      logical force

      character*1 cflag
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 noset,elset,orname(*)
      character*81 set(*),prset(*)
      character*87 filab(*)

      integer nset,istartset(*),iendset(*),ialset(*),nprint,ipkon(*),
     &  mi(*),nstate_,ii,jj,iset,l,limit,node,ipos,ithermal,ielem,
     &  nelem,kon(*),inotr(2,*),ntrans,ielorien(mi(3),*),norien,nk,ne,
     &  inum(*),nfield,ikin,nodes,ne0,nope,mt,ielmat(mi(3),*),iface,
     &  jfaces,mortar,islavsurf(2,*),ielprop(*),nload,
     &  nelemload(2,*)

      real*8 v(0:mi(2),*),t1(*),fn(0:mi(2),*),stx(6,mi(1),*),bhetot,
     &  eei(6,mi(1),*),xstate(nstate_,mi(1),*),ener(mi(1),*),energytot,
     &  volumetot,co(3,*),qfx(3,mi(1),*),rftot(0:3),ttime,time,
     &  trab(7,*),orab(7,*),vold(0:mi(2),*),enerkintot,thicke(mi(3),*),
     &  eme(6,mi(1),*),prop(*),veold(0:mi(2),*),xload(2,*)

      intent(in) set,nset,istartset,iendset,ialset,nprint,
     &  prlab,prset,v,t1,fn,ipkon,lakon,stx,eei,xstate,ener,
     &  mi,nstate_,ithermal,co,kon,qfx,ttime,trab,inotr,ntrans,
     &  orab,ielorien,norien,nk,ne,inum,filab,vold,ikin,ielmat,thicke,
     &  eme,islavsurf,mortar,time,ielprop,prop,veold,orname,
     &  nelemload,nload,sideload,xload

      mt=mi(2)+1

      !Open the .vtk file for writing
      open(unit=10, file='output.vtk', status='unknown')


      !     Write VTK file header
      write(10,*) '# vtk DataFile Version 3.0'
      write(10,*) 'CalculiX Output Data'
      write(10,*) 'ASCII'
      write(10,*) 'DATASET UNSTRUCTURED_GRID'

      !     Write nodal coordinates
      write(10,*) 'POINTS', nk, 'float'
      do node=1, nk
         write(10, '(3F12.5)') co(1,node), co(2,node), co(3,node)
      enddo

      do ielem=1, ne
        select case (lakon(ielem)(1:2))
        case ('C3D8')  ! 8-node quad element
            print*, 'Eight node quad element'
        end select
      end do

      close(10)

      return 
      end







