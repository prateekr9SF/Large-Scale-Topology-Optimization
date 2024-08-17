
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
     &  nelemload(2,*), j

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
      open(unit=10, file='elastic_field.dat', status='unknown')

      !Write Tecplot 360 file header
      write(10,*) 'TITLE = "Simulation Results"'
      write(10, '(A)') 'VARIABLES = "X", "Y", "Z",' //
     &                '"DispX", "DispY", "DispZ"'

     

    ! Inspect lakon to identify element type
      select case (lakon(1))
      case ('C3D8')  ! 8-node quad element
       write(10, '(A)') 'ZONE STRANDID=1,'
       write(10, '(A)') ' SOLUTIONTIME=1,'
       write(10, '(A, I0)') ' NODES=', nk
       write(10, '(A, I0)') ' ELEMENTS=', ne
       write(10, '(A)') ', DATAPACKING=POINT,'
       write(10, '(A)') ' ZONETYPE=FEQUADRILATERAL'
      

      case ('C3D4')  ! 4-node tetrahedral element
       write(10, '(A)') 'ZONE STRANDID=1,'
       write(10, '(A)') ' SOLUTIONTIME=1,'
       write(10, '(A, I0)') ' NODES=', nk
       write(10, '(A, I0)') ' ELEMENTS=', ne
       write(10, '(A)') ', DATAPACKING=POINT,'
       write(10, '(A)') ' ZONETYPE=TETRAHEDRON'
      end select

      !     Write nodal coordinates
      do node = 1, nk
        write(10, '(6F12.5)') co(1,node), co(2,node), co(3,node),
     &        v(1,node), v(2,node), v(3,node)
      end do



      ! Write element connectivity based on the element type
      ! Write element connectivity based on the element type
      do ielem = 1, ne
        select case (lakon(ielem)(1:4))
        case ('C3D8')  ! 8-node hexahedral element
           do j = 1, 8
              write(10, '(I6)', advance='no') kon(ipkon(ielem) + j)
           end do
           write(10,*)  ! Newline after each element's nodes
        case ('C3D4')  ! 4-node tetrahedral element
           do j = 1, 4
              write(10, '(I6)', advance='no') kon(ipkon(ielem) + j)
           end do
           write(10,*)  ! Newline after each element's nodes
        case ('C3D6')  ! 6-node wedge element
           do j = 1, 6
              write(10, '(I6)', advance='no') kon(ipkon(ielem) + j)
           end do
           write(10,*)  ! Newline after each element's nodes
        case ('C3D20')  ! 20-node quadratic hexahedral element
           do j = 1, 20
              write(10, '(I6)', advance='no') kon(ipkon(ielem) + j)
           end do
           write(10,*)  ! Newline after each element's nodes
        ! Add more cases as needed for other element types
        case default
           ! Handle unknown or unsupported element types
           write(10,*) 'Unknown element type: ', lakon(ielem)(1:4)
        end select
      end do
          

      close(10)

      return 
      end







