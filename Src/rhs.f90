!----------------------------------------------------------------------------------85

subroutine rhs(q, qdot, substep_number, t, t0, lambda_max, lambda_max_extra_operator, &
               unphysical, iphi_bad)

! RHS of ODE system for governing equations for mass, momentum, and internal energy.
! t and t0 are needed to compute fargo_factor (if needed) which has (t - t0).

use grid
use partition_data
use dof_indices
use gravity
use boundary_condition_data
use logical_units, only: lun_lambda
use math_constants
use sponge
use rotating_frame
#ifdef mpi_code
   use mpi
#endif
use viscous ! Laminar terms and LES model
use transposes_of_q_and_qdot
use fargo_or_plotting_shift
! Location of maximum eigenvalue lambda_max_ap_global and bulk viscosity at
! this location.
use artificial_pressure_module, only: apply_artificial_pressure, have_dilatation, &
     ir_max_ap, iz_max_ap, iphi_max_ap, beta_Delta_max
use pade_filter
use thermal_parameters

implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q, qdot
integer, intent(in)                            :: substep_number
real(8), intent(in)                            :: t          ! needed to compute fargo_factor which has (t-t0)
real(8), intent(in )                           :: t0         ! needed to compute fargo_factor
real(8), intent(out)                           :: lambda_max ! maximum eigenvalue/pi from all the terms
                                                             ! in the first sub-step.
real(8), intent(out)                           :: lambda_max_extra_operator
logical, intent(out)                           :: unphysical
integer, intent(out)                           :: iphi_bad

! Local:

! Eigenvalues for time step determination.
real(8) :: lambda_max_viscous, lambda_max_non_euler
real(8) :: lambda_max_ap_global, &
           lambda_max_extra_operator_global, lambda_max_viscous_global
real(8) :: lambda_max_euler, lambda_max_euler_fargo, lambda_max_euler_non_fargo
real(8) :: lambda_max_without_fargo, lambda_max_with_fargo

integer :: iz, iphi, ir, ier
logical :: worth_doing_fargo, dominated_by_phi, get_lambda_max
real(8) :: eint_target ! For Newtonian cooling

! These should be global.
real(8), dimension(sr:er, sphi:ephi, nz) :: pressure, p_art

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' rhs has been called'
      print *, ' nz = ', nz
      print *, ' suppress_z_derivatives_when_nz_not_1 = ', suppress_z_derivatives_when_nz_not_1
   end if
#endif

! "extra_operator" refers to the extra operator for Fargo resulting from the chain rule.
! See Shariff & Wray, ApJ, Supplement Series.   
lambda_max_extra_operator        = 0.d0 ! Within each processor
lambda_max_extra_operator_global = 0.d0 ! Global over all processors
lambda_max_euler_fargo           = 0.d0
lambda_max_euler_non_fargo       = 0.d0

! These flags are used to possibly avoid calculating the dilatation in an expensive way.
! We are beginning a new substep so all old quantities should be recomputed if needed.   
have_dilatation               = .false.   
have_velocity_gradient_tensor = .false.
have_strain_tensor            = .false.

qdot = 0.0d0


if (substep_number .eq. 1) then
   if (apply_fargo) then
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, about to call get_lambda_euler_fargo'
#endif
call get_lambda_euler_fargo(q, lambda_max_euler_non_fargo, lambda_max_euler_fargo)
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, returned from get_lambda_euler_fargo'
#endif
   else
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, about to call get_lambda_euler_non_fargo'
#endif   
   call get_lambda_euler_non_fargo(q, lambda_max_euler_non_fargo, dominated_by_phi)
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, returned from get_lambda_euler_non_fargo'
#endif
      lambda_max_euler_fargo = 0.d0
   end if
   add_fargo_extra_operator_now = .false.
else
   ! Note: apply_fargo_this_step will have been determined in the first sub-step.
   if (apply_fargo_this_step) then
      add_fargo_extra_operator_now = apply_fargo_extra_operator
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, doing loop for fargo_factor'
#endif      
      do ir = 1, nr
         ! The division by "r" is because the r derivative terms have an 1/r in front: 
         fargo_factor_over_r(ir) = - (t - t0) * d_Omega_bar_dr(ir) / rgrid(ir)
         fargo_factor       (ir) = - (t - t0) * d_Omega_bar_dr(ir)
      end do
   end if
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, finished loop for fargo_factor'
#endif   
end if

! With Fargo factor set we can now compute derivatives with the extra operator needed:

! --------------------------Artificial pressure ----------------------
if (apply_artificial_pressure) then
   ! This subroutine is in artificial_pressure.f90
#ifdef debug_print
   if (my_node .eq. 0) print *, ' routine rhs, node 0: Calling get_artificial_pressure'
#endif
   ! Adding in of artificial pressure happens later. This routine is in artificial_pressure.f90,
   ! It is only in the first-substep that we obtain lambda_max_ap_global (the eigenvalue for
   ! the artificial pressure term.)
   get_lambda_max = (substep_number .eq. 1)

   !print *, ' debug: About to call get_artificial_pressure'
   !print *, ' substep # = ', substep_number
   !read(5, *)
   
   call get_artificial_pressure (get_lambda_max, q, p_art, lambda_max_ap_global)
   if ((my_node .eq. 0) .and. get_lambda_max) then
      ! The items in the last two lines are in module artificial_pressure
      print *, ' lambda_max_ap_global = ', lambda_max_ap_global, &
           ' at ir = ', ir_max_ap, ' iz = ', iz_max_ap, ' iphi = ', iphi_max_ap, &
           ' beta_Delta_max = ', beta_Delta_max
   end if
else
   lambda_max_ap_global = 0.d0
end if

! Calculate pressure:
if (isothermal) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            pressure(ir, iphi, iz) = q(ir, iphi, iz, irho) * ci_squared_initial(ir, iz)
         end do
      end do
   end do
else
   ! Note: gm1 is assigned in set_up_thermal_parameters:
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            pressure(ir, iphi, iz) = q(ir, iphi, iz, ener) * gm1 ! From internal energy.
         end do
      end do
   end do
end if

! Note: Artificial pressure is what we are calling the artificial bulk viscosity treatment.
if (apply_artificial_pressure) pressure = pressure + p_art

! Check for negative density, pressure, and internal energy (for the non-isothermal case):
call check_for_unphysical(q, pressure, substep_number, unphysical, iphi_bad)
if (unphysical) then
   ! Debug
   if (my_node .eq. 0) then
      print *, ' In rhs.  t = ', t
      print *, ' substep_number = ', substep_number
      print *, ' unphysical = ', unphysical
      print *, ' returning'
    end if
    return
end if

! We do the viscous terms first because if we compute the strain rate or velocity gradient
! tensor it can be used to compute the dilatation for the artificial pressure.
if (apply_viscosity) then
   ! This is located in viscous_terms.f90
   call add_viscous_terms_to_qdot(t, q, qdot, lambda_max_viscous) ! In viscous_terms.f90
#ifdef debug_print
   if (my_node .eq. 0) print *, ' routine rhs, node 0: Returned from add_viscous_terms_to_qdot'
#endif
   if (substep_number .eq. 1) then
#ifdef mpi_code
      ! Note: Only proc. 0 will have lambda_max_viscous_global
      call mpi_reduce(lambda_max_viscous, lambda_max_viscous_global, 1, mpi_double, mpi_max, &
           0, mpi_comm_world, ier)
#else
      lambda_max_viscous_global = lambda_max_viscous
#endif
   end if
else
   lambda_max_viscous_global = 0.d0
end if

if (apply_viscosity .and. (my_node .eq. 0) .and. (substep_number .eq. 1)) then
   print *, ' lambda_max_viscous_global = ', lambda_max_viscous_global
end if


if (substep_number .eq. 1) then
   if (my_node .eq. 0) then
      ! Eigenvalues lambda for the purpose of satisfying the cfl constraint
      ! lammba_max * dt = cfl
      lambda_max_non_euler = max(lambda_max_viscous_global, lambda_max_ap_global)
      lambda_max_without_fargo = max(lambda_max_euler_non_fargo, lambda_max_non_euler)
      lambda_max_with_fargo    = max(lambda_max_euler_fargo,     lambda_max_non_euler)

      worth_doing_fargo = (lambda_max_with_fargo .lt. lambda_max_without_fargo)

      ! See if we should suggest Fargo to the user:
      if ((nphi .ne. 1) .and. (.not. apply_fargo)) then
         if (worth_doing_fargo) then
            print *, ' subroutine rhs: We suggest that you use Fargo'
            print *, ' lambda_max_with_fargo    = ', lambda_max_with_fargo
            print *, ' lambda_max_without_fargo = ', lambda_max_without_fargo

            print *, ' lambda_max_euler_non_fargo = ', lambda_max_euler_non_fargo
            print *, ' lambda_max_euler_fargo     = ', lambda_max_euler_fargo            
            print *, ' lambda_max_non_euler       = ', lambda_max_non_euler
         end if
      end if

      ! See if we should do Fargo this step:
      if (apply_fargo .and. worth_doing_fargo) then
         apply_fargo_this_step = .true.
         lambda_max = lambda_max_with_fargo
      else
         apply_fargo_this_step = .false.
         lambda_max = lambda_max_without_fargo
      end if
   end if ! my_node = 0
#ifdef mpi_code
   call mpi_bcast(apply_fargo_this_step, 1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(lambda_max,            1, mpi_double,  0, mpi_comm_world, ier)
#endif
end if ! first_substep


if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' routine rhs, node 0: Calling add_z_derivatives'
#endif
   call add_z_derivatives(q, pressure, qdot)
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: about to do if check on gravity flag in subroutine rhs'
#endif

! Gravity terms in the momentum equations:   
if (gravity_flag) call add_gravity_terms(q, qdot) ! In gravity.f90

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: in rhs about to add centrifugal term'
#endif

! Centrifugal term in radial momentum eq.:
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         qdot(ir,iphi,iz,rmom) = qdot(ir,iphi,iz,rmom) + &
                  q(ir,iphi,iz,amom)**2/q(ir,iphi,iz,irho)/rgrid(ir)**3
      end do
   end do
end do

! Rotating frame terms if required:
if (apply_rotating_frame) call add_coriolis_and_centrifugal_terms(q, qdot)

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: in rhs finished adding centrifugal term'
#endif

if (nr .ne. 1) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: in rhs about to call add_radial_derivatives'
#endif      
   call add_radial_derivatives(q, pressure, qdot)
end if

if (nphi .ne. 1) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: about to call add_phi_derivatives'
#endif
   call add_phi_derivatives(q, qdot, pressure, lambda_max_extra_operator)
#ifdef mpi_code
   ! This reduction takes place to proc. 0:
   call mpi_reduce(lambda_max_extra_operator, lambda_max_extra_operator_global, 1, mpi_double, mpi_max, &
        0, mpi_comm_world, ier)
#else
   lambda_max_extra_operator_global = lambda_max_extra_operator
#endif
end if

! Newton's law of cooling if desired:
if ((.not. isothermal) .and. (apply_newtonian_cooling)) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            eint_target = q(ir, iphi, iz, irho) * ci_squared_initial(ir, iz) / gm1
            qdot(ir,iphi,iz,ener) = qdot(ir,iphi,iz,ener) + &
                     (eint_target - q(ir,iphi,iz,ener)) / tau_newtonian_cooling
         end do
      end do
   end do
   lambda_max = max(lambda_max, 1.d0/tau_newtonian_cooling)
end if


! The dilatation, etc. will now pertain to the middle of a substep.
have_dilatation               = .false.   
have_velocity_gradient_tensor = .false.
have_strain_tensor            = .false.

! Sync. for safety:
#ifdef mpi_code
call mpi_barrier(mpi_comm_world, ier)
#endif

end subroutine rhs

!----------------------------------------------------------------------------------85
