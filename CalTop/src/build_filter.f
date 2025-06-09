      subroutine mafillsm_filter2_full(ne, ttime, time,
     &  ne0, nea, neb, elCentroid, rmin, filternnz,
     &  FilterMatrixs, rowFilters, colFilters,
     &  filternnzElems, elarr, fnnzassumed)

      implicit none

      integer ne0, filternnz, ne, nea, neb, i, j, ii, fnnzassumed
      integer elarr(*), filternnzElems(*)
      integer rowFilters(fnnzassumed,*), colFilters(fnnzassumed,*)
      integer rowval, colval, dummy1

      real*8 FilterMatrixs(fnnzassumed,*), ttime, time, rmin
      real*8 elCentroid(3,*), xi, yi, zi, xj, yj, zj
      real*8 d, rmind

      dummy1 = fnnzassumed / 3
      filternnz = 0

c-- Construct upper triangle of filter matrix
      do ii = nea, neb
        i = elarr(ii) + 1
        xi = elCentroid(1, i)
        yi = elCentroid(2, i)
        zi = elCentroid(3, i)
        filternnzElems(i) = 0

        do j = i, ne
          xj = elCentroid(1, j)
          yj = elCentroid(2, j)
          zj = elCentroid(3, j)

          d = sqrt((xj - xi)**2 + (yj - yi)**2 + (zj - zi)**2)
          rmind = rmin - d

          if (rmind .gt. 0.d0) then
            filternnzElems(i) = filternnzElems(i) + 1
            if (filternnzElems(i) .gt. dummy1) exit

            rowFilters(filternnzElems(i), i) = i
            colFilters(filternnzElems(i), i) = j
            FilterMatrixs(filternnzElems(i), i) = rmind**2
            filternnz = filternnz + 1
          endif
        enddo
      enddo

c-- Expand symmetrically: H_ji = H_ij
      do i = 1, ne
        do j = 1, filternnzElems(i)
          rowval = rowFilters(j, i)
          colval = colFilters(j, i)
          if (colval .gt. i) then
            filternnzElems(colval) = filternnzElems(colval) + 1
            if (filternnzElems(colval) .gt. dummy1) cycle
            rowFilters(filternnzElems(colval), colval) = colval
            colFilters(filternnzElems(colval), colval) = i
            FilterMatrixs(filternnzElems(colval), colval) =
     &        FilterMatrixs(j, i)
          endif
        enddo
      enddo

      return
      end
