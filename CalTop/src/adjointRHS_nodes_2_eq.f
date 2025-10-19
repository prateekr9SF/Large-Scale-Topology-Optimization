c======================================================================
c  Build adjoint RHS in equation space from a nodal adjoint RHS
c  Linear static, mechanical only, SPCs only (no MPCs, no thermal)
c======================================================================
      subroutine adjrhs_scatter_linstatic_nompc(
     &     nk, neq, mi, nactdof,
     &     rhs_nodal,           ! (0:mi(2), 1..nk) use only 1..mi(2)
     &     rhs_eq,              ! (1..neq)
     &     nboun, nodeboun, ndirboun)

      implicit none
c----- inputs
      integer           nk, neq, nboun
      integer           mi(*)
      integer           nactdof(0:mi(2),*)
      integer           nodeboun(*), ndirboun(*)
      real*8            rhs_nodal(0:mi(2),*)
c----- outputs
      real*8            rhs_eq(*)

c----- locals
      integer           i, j, eq

c===== 0) zero equation-space RHS
      do eq = 1, neq
         rhs_eq(eq) = 0.d0
      enddo

c===== 1) enforce homogeneous adjoint on primal SPC DOFs: rhs_nodal = 0 there
      do i = 1, nboun
         if (ndirboun(i) .ge. 1 .and. ndirboun(i) .le. mi(2)) then
            rhs_nodal(ndirboun(i), nodeboun(i)) = 0.d0
         endif
      enddo

c===== 2) scatter nodal → equation space via nactdof
c      (only active DOFs appear in the equation system)
      do i = 1, nk
         do j = 1, mi(2)
            eq = nactdof(j,i)
            if (eq .gt. 0) then
               rhs_eq(eq) = rhs_nodal(j,i)
            endif
         enddo
      enddo

      return
      end
