      subroutine write_tecplot(co, v, stx, nk, ne, ipkon, kon, lakon, mi)
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

      ! Properly initialize the filename
      filename = 'elastic_field.dat'  ! Fixed file name for Tecplot

      ! Open the file for writing
      open(unit=unit, file=filename, status='unknown')

      ! Write Tecplot header
      write(unit, '(A)') 'TITLE = "Elastic Field Data"'
      write(unit, '(A)') 'VARIABLES = "X", "Y", "Z", "DisplX", "DisplY", "DisplZ"'

    ! Write zone information based on the element type
      if (lakon(1)(1:2) .eq. 'C3') then
        write(unit, '(A, I10, A, I10, A)') 'ZONE N = ', nk, ', E = ', ne, ', DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE'
      else if (lakon(1)(1:2) .eq. 'C4') then
        write(unit, '(A, I10, A, I10, A)') 'ZONE N = ', nk, ', E = ', ne, ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL'
      else if (lakon(1)(1:2) .eq. 'C10') then
        write(unit, '(A, I10, A, I10, A)') 'ZONE N = ', nk, ', E = ', ne, ', DATAPACKING=BLOCK, ZONETYPE=FETETRAHEDRON'
      endif

    ! Write node coordinates (X, Y, Z)
      write(unit, '(A)') 'NODAL DATA:'
      do node = 1, nk
        write(unit, '(3(1X,E14.7))') co(1,node), co(2,node), co(3,node)
      end do

      ! Write displacements (DisplX, DisplY, DisplZ)
      write(unit, '(A)') 'DISPLACEMENTS:'
      do node = 1, nk
        write(unit, '(3(1X,E14.7))') v(0,node), v(1,node), v(2,node)
      end do

      ! Write element connectivity
      write(unit, '(A)') 'ELEMENT CONNECTIVITY:'
      do elem = 1, ne
        eltype = lakon(elem)(1:2)  ! Extract the element type

        if (eltype .eq. 'C3') then
            ! For 3-node triangular elements
            n1 = kon(ipkon(elem))
            n2 = kon(ipkon(elem) + 1)
            n3 = kon(ipkon(elem) + 2)
            write(unit, '(I10, 3I10)') n1, n2, n3  ! 1-based indexing for Tecplot

        else if (eltype .eq. 'C4') then
            ! For 4-node quadrilateral elements
            n1 = kon(ipkon(elem))
            n2 = kon(ipkon(elem) + 1)
            n3 = kon(ipkon(elem) + 2)
            n4 = kon(ipkon(elem) + 3)
            write(unit, '(I10, 4I10)') n1, n2, n3, n4  ! 1-based indexing for Tecplot

        else if (eltype .eq. 'C10') then
            ! For 4-node tetrahedral elements
            n1 = kon(ipkon(elem))
            n2 = kon(ipkon(elem) + 1)
            n3 = kon(ipkon(elem) + 2)
            n4 = kon(ipkon(elem) + 3)
            write(unit, '(I10, 4I10)') n1, n2, n3, n4  ! 1-based indexing for Tecplot

        endif
      end do

      ! Close the file
      close(unit)

      return
      end