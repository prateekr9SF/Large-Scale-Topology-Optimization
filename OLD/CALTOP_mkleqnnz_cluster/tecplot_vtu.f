
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 2022-2025 Prateek Ranjan
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
     
      subroutine tecplot_vtu(set,nset,istartset,iendset,ialset,nprint,
     &  prlab,prset,v,t1,fn,ipkon,lakon,stx,eei,xstate,ener,
     &  mi,nstate_,ithermal,co,kon,qfx,ttime,trab,inotr,ntrans,
     &  orab,ielorien,norien,nk,ne,inum,filab,vold,ikin,ielmat,thicke,
     &  eme,islavsurf,mortar,time,ielprop,prop,veold,orname,
     &  nelemload,nload,sideload,xload)

      implicit none

      character*1 cflag
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 sideload(*)
      character*80 noset,elset,orname(*)
      character*81 set(*),prset(*)
      character*87 filab(*)
      character*200 header

      character*200 subheader_1
      character*200 subheader_2
      character*200 subheader_3
      character*200 subheader_4
      character*200 subheader_5



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

       open(unit=20, file='ipkon_debug.txt', status='unknown')

       do ielem = 1, ne
           ! Skip invalid elements
           if (ipkon(ielem) == -1) then
               write(20,*) 'WARNING: ipkon(', ielem, ') = -1'
               cycle
           end if
       
           ! Write 4-node connectivity matrix properly formatted
           do j = 1, 4
               write(20, '(I6)', advance='no') kon(ipkon(ielem) + j)
           end do
           write(20,*)  ! Move to new line for next element
       end do
       
       close(20)
       
       

! Open the .vtu file for writing
       open(unit=10, file='elastic_Field.vtu', status='unknown')

! Write VTU XML header
     
       header = '<?xml version="1.0"?>' // 
     & '<VTKFile type="UnstructuredGrid" version="0.1" ' // 
     & 'byte_order="LittleEndian">'
       write(10,'(A)') header

       write(10,'(A)') '  <UnstructuredGrid>'
       
      
       write(10,'(A)', advance='no') '    <Piece NumberOfPoints="'
       write(10,'(I0)', advance='no') nk
       write(10,'(A)', advance='no') '" NumberOfCells="'
       write(10,'(I0)', advance='no') ne
       write(10,'(A)') '" >'
       
       
       

! Write nodal coordinates
      write(10,'(A)') '      <Points>'
      subheader_1 = '        <DataArray type="Float64" ' //
     &            'NumberOfComponents="3" ' //
     &            'format="ascii">'
      write(10,'(A)') subheader_1
 
      do node = 1, nk
        write(10,'(3F12.5)') co(1,node), co(2,node), co(3,node)
      end do
      write(10,'(A)') '        </DataArray>'
      write(10,'(A)') '      </Points>'


! Write connectivity
      write(10,'(A)') '      <Cells>'
      subheader_2 = '        <DataArray type="Int32" ' //
     &                      'Name="connectivity" ' //
     &                      'format="ascii">'
      write(10,'(A)') subheader_2
    
      do ielem = 1, ne
        do j = 1, 4  ! TODO: Add conditionals for element type
          write(10,'(I6)', advance='no') kon(ipkon(ielem) + j)
        end do
        write(10,*)
      end do
      write(10,'(A)') '        </DataArray>'


! Write offsets
      subheader_3 = '        <DataArray type="Int32" ' //
     &                      'Name="offsets" ' //
     &                      'format="ascii">'
      write(10,'(A)') subheader_3

      do ielem = 1, ne
        write(10,'(I0)') ielem * 4
      end do
      write(10,'(A)') '        </DataArray>'

! Write cell types
      subheader_4 = '        <DataArray type="UInt8" ' //
     &                      'Name="types" ' //
     &                      'format="ascii">'
      write(10,'(A)') subheader_4

      do ielem = 1, ne
        write(10,'(I6)') 10
      end do
      write(10,'(A)') '        </DataArray>'
      write(10,'(A)') '      </Cells>'


! Write nodal displacements

      write(10,'(A)') '      <PointData Scalars="Displacement">'

      subheader_5 = '        <DataArray type="Float64" ' //
     &                      'Name="Displacement" ' //
     &                       'NumberOfComponents="3" ' //                   
     &                      'format="ascii">'
      write(10,'(A)') subheader_5
      
      do node = 1, nk
        write(10,'(3F12.8)') v(1,node), v(2,node), v(3,node)
      end do
      write(10,'(A)') '        </DataArray>'
      write(10,'(A)') '       </PointData>'

      ! Close XML tags
      write(10,'(A)') '    </Piece>'
      write(10,'(A)') '  </UnstructuredGrid>'
      write(10,'(A)') '</VTKFile>'

      close(10)

      return
      end