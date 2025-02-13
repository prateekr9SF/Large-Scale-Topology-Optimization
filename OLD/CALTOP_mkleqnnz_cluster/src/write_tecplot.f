      subroutine write_vtk(co, v, stx, nk, ne, ipkon, kon, lakon, mi)
      implicit none

      ! Parameters
      character*20 filename
      

      integer nk, ne, ipkon(*), kon(*), mi(*)
      character*8 lakon(*)
      real*8 co(3,*), v(0:mi(2),*), stx(6,mi(1),*)

      ! Local variables
      integer i, elem, node
      integer n1, n2, n3, n4
      character*2 eltype
      integer mt, unit

      mt = mi(2) + 1
      unit = 10

      filename = 'elastic_field.vtk'  ! Fixed file name

      ! Open the file for writing
      open(unit=unit, file=filename, status='unknown')

      ! Write VTK header
      write(unit, '(A)') '# vtk DataFile Version 3.0'
      write(unit, '(A)') 'VTK output from Fortran'
      write(unit, '(A)') 'ASCII'
      write(unit, '(A)') 'DATASET UNSTRUCTURED_GRID'

      ! Write nodal coordinates
      write(unit, '(A,I10)') 'POINTS ', nk, ' double'
      do node = 1, nk
         write(unit, '(3(1X,E14.7))') co(1,node), co(2,node), co(3,node)
      end do

      ! Write element connectivity
      write(unit, '(A,I10)') 'CELLS ', ne, ', ', 4*ne
      do elem = 1, ne
         eltype = lakon(elem)(1:2)  ! Extract the element type

         if (eltype .eq. 'C3') then
            ! For 3-node triangular elements
            n1 = kon(ipkon(elem))
            n2 = kon(ipkon(elem) + 1)
            n3 = kon(ipkon(elem) + 2)
            write(unit, '(I10, 3I10)') 3, n1-1, n2-1, n3-1  ! 0-based indexing

         else if (eltype .eq. 'C4') then
            ! For 4-node quadrilateral elements
            n1 = kon(ipkon(elem))
            n2 = kon(ipkon(elem) + 1)
            n3 = kon(ipkon(elem) + 2)
            n4 = kon(ipkon(elem) + 3)
            write(unit, '(I10, 4I10)') 4, n1-1, n2-1, n3-1, n4-1
         end if
      end do

      ! Write cell types
      write(unit, '(A,I10)') 'CELL_TYPES ', ne
      do elem = 1, ne
         eltype = lakon(elem)(1:2)

         if (eltype .eq. 'C3') then
            write(unit, '(I10)') 5   ! VTK type 5: VTK_TRIANGLE
         else if (eltype .eq. 'C4') then
            write(unit, '(I10)') 9   ! VTK type 9: VTK_QUAD
         end if
      end do

      ! Write point data (displacements and stresses)
      write(unit, '(A,I10)') 'POINT_DATA ', nk
      write(unit, '(A)') 'VECTORS Displacement double'
      do node = 1, nk
         write(unit, '(3(1X,E14.7))') v(0,node), v(1,node), v(2,node)
      end do

      ! Writing the stress tensor components individually to avoid line continuation
    !  write(unit, '(A)') 'TENSORS Stress double'
    !  do node = 1, nk
    !    write(unit, '(1X,E14.7,1X,E14.7,1X,E14.7)') stx(1,node), stx(2,node), stx(3,node)
    !    write(unit, '(1X,E14.7,1X,E14.7,1X,E14.7)') stx(4,node), stx(5,node), stx(6,node)
    !  end do

      ! Close the file
      close(unit)

      return
      end
