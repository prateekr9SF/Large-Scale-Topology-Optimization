!
!     mafillsm_filter2 - Constructs a sparse filter matrix for density filtering
!     in topology optimization based on element centroid distances.
!
!     For each element i, neighbors j within a filter radius rmin are stored
!     in sparse row format: (i,j) with weight = rmin - distance(i,j)
!
!     The filter matrix is asymmetric and unnormalized.
!
      subroutine mafillsm_filter2(ne, ttime, time,
     &  ne0, nea, neb, elCentroid,
     &  rmin, filternnz,
     &  FilterMatrixs, rowFilters, colFilters,
     &  filternnzElems, elarr, fnnzassumed)

      implicit none

!==== INPUTS ====
      integer ne, ne0, nea, neb, fnnzassumed
      real*8 ttime, time, rmin
      real*8 elCentroid(3, *)              ! Element centroids (x, y, z)
      integer elarr(*)                     ! Element renumbering array

!==== IN/OUT ====
      integer filternnz                    ! Total number of nonzeros (global counter)
      real*8 FilterMatrixs(fnnzassumed, *) ! Filter weights
      integer rowFilters(fnnzassumed, *)   ! Row indices for sparse matrix
      integer colFilters(fnnzassumed, *)   ! Column indices for sparse matrix
      integer filternnzElems(*)            ! Local nnz count per row

!==== LOCALS ====
      integer i, j, ii, row_idx, dummy1
      real*8 xi, yi, zi, xj, yj, zj
      real*8 dist_sq, rmind, weight

! Conservative cap on the number of neighbors per row
      dummy1 = fnnzassumed / 3

      

! Loop over element range assigned to this thread
      do ii = nea, neb

        i = elarr(ii) + 1    ! Convert from 0-based to 1-based indexing

        ! Coordinates of element i
        xi = elCentroid(1, i)
        yi = elCentroid(2, i)
        zi = elCentroid(3, i)

        
       
        filternnzElems(i) = 0
        filternnz = 0

        ! Only search neighbors j >= i
        do j = i, ne

          ! Coordinates of element j
          xj = elCentroid(1, j)
          yj = elCentroid(2, j)
          zj = elCentroid(3, j)

          dist_sq = (xj - xi)**2 + (yj - yi)**2 + (zj - zi)**2
          rmind = rmin*rmin - dist_sq

          if (rmind .GE. 0.d0) then
            filternnz = filternnz + 1
            filternnzElems(i) = filternnzElems(i) + 1
            
            weight = rmin - dist_sq**0.5
            rowFilters(filternnz, i) = i
            colFilters(filternnz, i) = j
            FilterMatrixs(filternnz, i) = weight
          end if

          if (filternnzElems(i) .EQ. dummy1) then
          EXIT
          end if

        end do

      end do

      return
      end

