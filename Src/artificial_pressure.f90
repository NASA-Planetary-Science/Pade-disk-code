!----------------------------------------------------------------------------------85

! The artificial pressure treatment is the same as the artificial bulk viscosity
! treatment suggested by Cook and Cabot (though our implementation may be a little
! different).

! p_art = - beta_Delta div u, where beta_Delta is the bulk viscosity.

!----------------------------------------------------------------------------------85

module artificial_pressure_module
   logical :: have_dilatation
   logical :: apply_artificial_pressure
   real(8) :: C_ap ! order unity constant in front
   ! These are allocated in subroutine activate_artificial_pressure which is
   ! located in activate_routines.f90
   real(8), allocatable, dimension(:, :, :) :: dil

   ! Location of the maximum of the eigenvalue lambda_max_ap.
   integer :: ir_max_ap, iz_max_ap, iphi_max_ap
   ! Bulk viscosity at the above location.
   real(8) :: beta_Delta_max   
end module artificial_pressure_module

!----------------------------------------------------------------------------------85

subroutine get_artificial_pressure (get_lambda_max, q, p_art, lambda_max_ap_global)

! Computes the artificial pressure using the length scale suggested by Mani et al.
! lambda_max is the returned "eigenvalue" to aid in time-step selection.

! Note: Information on the location (ir_max_ap, iz_max_ap, iphi_max_ap) of the maximum
! eigenvalue is in module artificial pressure as is the bulk viscosity at this
! location (beta_Delta_max).

use grid
use partition_data
use artificial_pressure_module
use transposes_of_q_and_qdot
use dof_indices
use math_constants
use pade_filter, only: smooth_field
#ifdef mpi_code
use mpi
#endif
implicit none
logical,                                        intent(in)  :: get_lambda_max
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in)  :: q
real(8), dimension(sr:er, sphi:ephi, nz),       intent(out) :: p_art
real(8),                                        intent(out) :: lambda_max_ap_global

! Locals:
real(8), dimension(sr:er, sphi:ephi, nz) :: beta_Delta ! Artificial bulk viscosity.

real(8) :: eps_filter_now, lambda_now, lambda_max
integer :: iz, iphi, ir
#ifdef mpi_code
integer :: ier, status(mpi_status_size), maxloc(3)
#endif
real(8), dimension(2) :: lambda_max_in, lambda_max_out

! This subroutine is below.  If we have the strain-tensor or the velocity gradient tensor, we will use them.
! The dilatation is sitting "artificial_pressure_module".
call dilatation(q) 

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! The array beta_Delta is the artificial bulk viscosity:
         ! On physical grounds, the bulk viscosity must always be positive:
         ! Option 1:
         beta_Delta(ir,iphi,iz) = C_ap * q(ir,iphi,iz,irho) * l_grid_squared(ir,iz) * abs(dil(ir,iphi,iz))
         ! Options 2 and 3:
         !beta_Delta(ir,iphi,iz) = -C_ap * q(ir,iphi,iz,irho) * l_grid_squared(ir,iz) * min(0.d0, dil(ir,iphi,iz))
      end do
   end do
end do
!call tecplot_scalar_in_horizontal_plane('bD_b4_smth', beta_delta, 1, 1.d0, 0.d0, 0)
eps_filter_now = 0.2d0
call smooth_field(1, beta_Delta, eps_filter_now)

!call tecplot_scalar_in_horizontal_plane('beta_Delta', beta_delta, 1, 1.d0, 0.d0, 0)
!call tecplot_scalar_in_horizontal_plane('dil_______', dil,        1, 1.d0, 0.d0, 0)
!stop


lambda_max = 0.0d0
ir_max_ap   = sr
iz_max_ap   = 1
iphi_max_ap = sphi
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! One other pi factor comes later:
         lambda_now  = pi * beta_Delta(ir,iphi,iz) / q(ir,iphi,iz,irho) / l_grid_squared(ir,iz)

         if (lambda_now .gt. lambda_max) then
            lambda_max = lambda_now
            ir_max_ap   = ir
            iz_max_ap   = iz
            iphi_max_ap = iphi
         end if
      end do
   end do
end do

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! Options 1 and 2:
         p_art(ir,iphi,iz) = -beta_Delta(ir,iphi,iz) * dil(ir,iphi,iz)
         ! Option 3
         !p_art(ir,iphi,iz) = -beta_Delta(ir,iphi,iz) * min(dil(ir,iphi,iz), 0.d0)
      end do
   end do
end do

beta_Delta_max = beta_Delta(ir_max_ap, iphi_max_ap, iz_max_ap)

if (get_lambda_max) then
#ifdef mpi_code
   lambda_max_in(1) = lambda_max
   lambda_max_in(2) = my_node
   maxloc(1) = ir_max_ap
   maxloc(2) = iz_max_ap
   maxloc(3) = iphi_max_ap
   ! The rank of the processor having the maximum will be int(lambda_max_out(2)).
   call mpi_allreduce(lambda_max_in, lambda_max_out, 1, mpi_2double_precision, mpi_maxloc, &
        mpi_comm_world, ier)
   lambda_max_ap_global = lambda_max_out(1)

   ! The processor with the max broadcasts the location of the maximum to all because
   ! they might be calling a routine to plot the plane with the max.
   call mpi_bcast(maxloc, 3, mpi_integer, int(lambda_max_out(2)), mpi_comm_world, ier)
   ir_max_ap   = maxloc(1)
   iz_max_ap   = maxloc(2)
   iphi_max_ap = maxloc(3)
   
   if (int(lambda_max_out(2)) .eq. my_node) then
      ! 1 is the tag
      call mpi_send(beta_Delta_max, 1, mpi_double,  0, 1, mpi_comm_world, ier)
   end if

   if (my_node .eq. 0) then
      call mpi_recv(beta_Delta_max, 1, mpi_double,  int(lambda_max_out(2)), 1, mpi_comm_world, status, ier)
   end if
#else
   lambda_max_ap_global = lambda_max
#endif
end if ! first substep

end subroutine get_artificial_pressure
   
!----------------------------------------------------------------------------------85

subroutine dilatation(q)

! Dilatation computation for artificial pressure (arising from bulk viscosity) term.

use grid
use partition_data
use artificial_pressure_module
use viscous
use dof_indices
use fargo_or_plotting_shift
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in) :: q

! Locals:
integer :: iz, iphi, ir

real(8), dimension(sr:er, sphi:ephi, nz) :: F, dF ! z space

real(8), dimension(sphi:ephi, sz_r:ez_r, nr)   :: F_r_space,   dF_r_space  
real(8), dimension(sr:er, sz_phi:ez_phi, nphi) :: F_phi_space, dF_phi_space

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: routine dilatation: first executable'
   print *, ' have_dilatation = ', have_dilatation
end if
#endif

! If we have the velocity gradient tensor, we can simply take its trace:
if (have_velocity_gradient_tensor) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            dil(ir,iphi,iz) = G(ir,iphi,iz,z_comp,z_comp) + G(ir,iphi,iz,r_comp,r_comp) + &
                 G(ir,iphi,iz,p_comp,p_comp)
         end do
      end do
   end do
   have_dilatation = .true.
   return
end if

! If we have the strain tensor, we can simply take its trace:
if (have_strain_tensor) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            dil(ir,iphi,iz) = T_zz(ir,iphi,iz) + T_rr(ir,iphi,iz) + T_pp(ir,iphi,iz)
         end do
      end do
   end do
   have_dilatation = .true.
   return
end if

! We don't have the dilatation so we need to compute it:
dil = 0.d0

if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
   !print *, ' debug: About to take z der for dil enter anything'
   !read(5, *)
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er   
            F(ir,iphi,iz) = q(ir,iphi,iz,zmom) / q(ir,iphi,iz,irho)
            ! print *, ' debug: iz = ', iz, ' F = ', F(ir,iphi,iz)
         end do
      end do
   end do
   call pade_diff_z(nbundle_z, F, dF)
   dil = dil + dF
end if

if (nr .ne. 1) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            F(ir,iphi,iz) = q(ir,iphi,iz,rmom) / q(ir,iphi,iz,irho) * rgrid(ir)
         end do
      end do
   end do   
   call transpose_z_to_r (1, F, F_r_space)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, F_r_space, dF_r_space)   
   call transpose_r_to_z (1, dF_r_space, dF)
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            dil(ir,iphi,iz) = dil(ir,iphi,iz) + dF(ir,iphi,iz) / rgrid(ir)
         end do
      end do
   end do

   ! Add Fargo extra operator if needed:
   if ((add_fargo_extra_operator_now) .and. (nphi .ne. 1)) then
      call transpose_z_to_phi(1, F, F_phi_space)
      call pade_diff_periodic(nbundle_phi, nphi, dphi, F_phi_space, dF_phi_space)   
      call transpose_phi_to_z (1, dF_phi_space, dF)
      do iz = 1, nz
         do iphi = sphi, ephi
            do ir = sr, er
               dil(ir,iphi,iz) = dil(ir,iphi,iz) - fargo_factor(ir)*dF(ir,iphi,iz)/rgrid(ir)
            end do
         end do
      end do
   end if
end if ! r derivatives

if (nphi .ne. 1) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            F(ir,iphi,iz) = q(ir,iphi,iz,amom) / q(ir,iphi,iz,irho) / rgrid(ir)
         end do
      end do
   end do
   call transpose_z_to_phi(1, F, F_phi_space)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, F_phi_space, dF_phi_space)   
   call transpose_phi_to_z (1, dF_phi_space, dF)
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            dil(ir,iphi,iz) = dil(ir,iphi,iz) + dF(ir,iphi,iz) / rgrid(ir)
         end do
      end do
   end do
end if

have_dilatation = .true.

end subroutine dilatation

!----------------------------------------------------------------------------------85

subroutine max_abs_of_array(array, array_max_abs, ir_max, iz_max, iphi_max)

! Obtains the maximum absolute value of array "array" in z space and its grid
! location over all processors. "array_max_abs" and (ir_max, iz_max, iphi_max) will
! be correct on all processors.

#ifdef mpi_code
use mpi
#endif
use partition_data
use grid
implicit none
real(8), dimension(sr:er, sphi:ephi, nz), intent(in) :: array
real(8), intent(out)                                 :: array_max_abs
integer, intent(out)                                 :: ir_max, iz_max, iphi_max

! Local:
integer :: iz, iphi, ir, ier, maxloc(3)
real(8) :: abs_now, my_max, array_max_in(2), array_max_out(2)

my_max   = 0.d0
ir_max   = 0
iz_max   = 0
iphi_max = 0
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         abs_now = abs(array(ir, iphi, iz))
         if (abs_now .gt. my_max) then
            my_max   = abs_now
            ir_max   = ir
            iz_max   = iz
            iphi_max = iphi
         end if
      end do
   end do
end do

#ifdef mpi_code
   array_max_in(1) = my_max
   array_max_in(2) = my_node
   ! The rank of the processor having the maximum will be int(lambda_max_out(2)).
   call mpi_allreduce(array_max_in, array_max_out, 1, mpi_2double_precision, mpi_maxloc, &
        mpi_comm_world, ier)
   array_max_abs = array_max_out(1)

   ! The processor with the max broadcasts the location of the maximum to all.
   maxloc(1) = ir_max
   maxloc(2) = iz_max
   maxloc(3) = iphi_max   
   call mpi_bcast(maxloc, 3, mpi_integer, int(array_max_out(2)), mpi_comm_world, ier)
   ir_max   = maxloc(1)
   iz_max   = maxloc(2)
   iphi_max = maxloc(3)
#else
   array_max_abs = my_max
#endif

end subroutine max_abs_of_array

!----------------------------------------------------------------------------------85

