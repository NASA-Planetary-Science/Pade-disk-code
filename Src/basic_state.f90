!----------------------------------------------------------------------------------85

! This module is used to store a basic state so fluctuation statistics can be computed
! if desired.  The arrays are allocated in subroutine initialize and given values
! in subroutine store_basic_state (for a non-restart).

! For a restart from a type 1 restart file, these arrays are given values in
! subroutine read_restart.  This is now deprecated and it is reccomended that
! the basic state be specified for each run, restart or not.
module basic_state
   implicit none
   logical :: have_basic_state
   real(8), allocatable, dimension(:, :) :: rho_basic, uphi_basic
   ! These are used to compute the fluctuation vorticity.
   real(8), allocatable, dimension(:, :) :: rho_basic_r_space, uphi_basic_r_space

   ! Basic state angular momentum and internal energy.  These will be allocated
   ! and used if the smoothing is to be done on the fluctuation field relative
   ! to the basic state:
   real(8), allocatable, dimension(:, :) :: amom_basic, eint_basic

contains
   
!----------------------------------------------------------------------------------85

subroutine store_basic_state(filter_relative_to_basic_state)

! Store the field currently in q as the basic state.

! This assumes that the field currently in q is axisymmetric so we can just pick
! the sphi index of the processor.

! This routine is typically called only at the initialization of a run, before
! the addition of any perturbations.  For a restart, the basic state is read from
! the restart file.

use q_array
use transposes_of_q_and_qdot
use grid
use dof_indices
use partition_data
use total_allocated_words, only: n_words_allocated
implicit none

logical, intent(in) :: filter_relative_to_basic_state

! local:
integer :: ir, iz, idof

! First executable:
allocate(rho_basic (sr:er, nz))
allocate(uphi_basic(sr:er, nz))
allocate(amom_basic(sr:er, nz), eint_basic(sr:er, nz))

! Used in vorticity computation to avoid a large memory stride.
allocate(uphi_basic_r_space(sz_r:ez_r, nr))
allocate(rho_basic_r_space (sz_r:ez_r, nr))
n_words_allocated = n_words_allocated + 4*mr*nz + mz_r*nr

do iz = 1, nz
   do ir = sr, er
      rho_basic (ir, iz) = q(ir, sphi, iz, irho)
      uphi_basic(ir, iz) = q(ir, sphi, iz, amom) / q(ir, sphi, iz, irho) / rgrid(ir)
      amom_basic(ir, iz) = q(ir, sphi, iz, amom)
      eint_basic(ir, iz) = q(ir, sphi, iz, ener)      
   end do
end do

if (nr .ne. 1) then
   call transpose_z_to_r (ndof, q, q_r_space)
   do ir = 1, nr
      do iz = sz_r, ez_r
         ! This is needed for the vorticity computation:
         uphi_basic_r_space(iz, ir) = q_r_space(sphi, iz, amom, ir) / &
              q_r_space(sphi, iz, irho, ir) / rgrid(ir)
         rho_basic_r_space(iz, ir) = q_r_space(sphi, iz, irho, ir)
      end do
   end do
end if

have_basic_state = .true.

end subroutine store_basic_state

!----------------------------------------------------------------------------------85

end module basic_state

!----------------------------------------------------------------------------------85
