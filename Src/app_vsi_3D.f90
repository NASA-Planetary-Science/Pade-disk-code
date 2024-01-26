!----------------------------------------------------------------------------------85

subroutine app_vsi_3D

use precision ! This defines xp
use dof_indices, only: irho, amom, zmom, rmom, ener
use grid
use boundary_condition_types 
use gravity, only: gz, gr
use q_array, only: q
use logical_units
use control_parameters, only: restart_file, save_file
use total_allocated_words, only: n_words_allocated
use math_constants, only: pi
use partition_data
use viscosity_types
use basic_state
use cpu_timing_module
use activate_routines
use set_up_routines
use gravity_types
#ifdef mpi_code
   use mpi
#endif
implicit none

external rhs ! name of rhs subroutine passed to rk4.
logical :: restart
integer :: istep_of_restart_file, file_version
real(xp) :: t_of_restart_file

integer :: istep, nsteps, istep0
real(xp) :: t, cfl, dt, dt_previous, dt_min
real(xp) :: H0, rho0, c0, r0
real(xp) :: zmax_over_H0, H0_over_r0, rsize_over_H0, lambda_r, V_Kep0, &
     Omega0, T_orbital0, phi_max_over_pi, rmin_over_H0

! Pade filter stuff:
logical :: apply_pade_filter, filter_relative_to_basic_state =.true.
character(3) :: eps_or_tau
real(xp) :: eps_filter, tau_filter

! For thermal set-up:
real(xp) :: gamma      = 1.4_xp
logical :: isothermal, apply_newtonian_cooling
!> Read from input file:
real(xp) :: tau_cooling_over_t_Kepler_at_mid_radius
!> Passed to set_up_thermal_parameters
real(xp) :: tau_newtonian_cooling

! For gravity activation:
real(xp) :: GM
real(xp) :: gz_uniform ! unused
logical :: gravity_flag
integer :: gravity_type

real(xp) :: pexp, qexp, c_exponent, H_exponent ! exponents for rho, temp, etc.
! Formula temps.:
real(xp) :: Omega, rho_temporary, term1, term2, term3, R_spherical, arg

integer :: iz, ir, iphi, ier
logical :: output_profiles, output_profiles_now
integer :: tecplot_interval, profiles_interval, fluctuation_ke_interval, save_interval
integer :: rms_interval
integer :: n_waves_in_r

! Local functions of r:
real(xp), allocatable, dimension(:) :: rho_midplane, Omega_K, ci, H

! To satisfy numerical vertical hydrostatic balance:
real(xp), allocatable, dimension(:, :, :) :: pressure, dpdz
logical :: perturb, wavy_perturbation

logical :: apply_fargo_trick, apply_fargo_correction, integer_shifts = .false., &
     suppress_z_derivatives_when_nz_not_1_arg = .false., periodic_z_arg = .false., &
     use_supplied_dt = .false.

! To satisfy numerical centrifugal balance:
real(xp), allocatable, dimension(:, :, :) :: p_r_space, dpdr_r_space, dpdr
real(xp) :: uphi_squared

logical :: apply_artificial_pressure, get_lambda_max_ap
real(xp) :: lambda_max_ap
logical :: use_rsize_for_domain, use_Manger_p
real(xp) :: C_ap

! Viscosity stuff, e.g., LES, etc.
logical :: apply_viscosity
integer :: viscosity_type
real(xp) :: C_Smag, nu_molecular, nu_b_molecular, Pr_molecular
real(xp) :: C_DDSV ! Coefficient for dilatation-dependent shear (will be read from input file).

integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced
real(8) :: c_sound_rmin, c_sound_rmax ! Note used

! Fluctuation kinetic energy:
real(xp) :: afke_z_middle, afke_r_middle, afke_phi_middle, &
           afke_z_upper, afke_r_upper, afke_phi_upper, &
           afke_total_middle, afke_total_upper

logical :: name_using_step

logical, parameter :: stretched_r = .false., stretched_z = .false.
real(xp), parameter :: r0_unif = 0._xp, z0 = 0._xp
integer, parameter :: nr_u = 0, nz_u = 0

character(7) :: sponge_type
logical :: apply_sponge
real(xp) :: rho1, rho2, d1, d2, tau_decay
integer :: n_decay_steps

character(25) :: filename_given
logical :: unphysical_flag
integer :: iphi_bad, istatus

! Arguments for tecplot output calls.
integer :: iphi_plot, iz_plot_midplane, iz_plot_one_and_half_H0
real(xp), parameter :: L_scale=1.0d0, T_scale=1.0d0, rho_scale=1.0d0

logical :: plot_pert = .false., plot_curl_rho_u = .true.

! For computing Reynolds phi averages from which Favre phi-time averaged means and stresses
! can be computed using TOOLS/phi_time_Favre_averages.f90.  This is valid only in a stationary
! state.
logical :: output_phi_Reynolds_averages
real(xp) :: time_interval_for_phi_Reynolds_averages
real(xp) :: next_t_for_phi_Reynolds_averages
logical, parameter :: average_in_z = .false. ! For Reynolds averages

real(xp), allocatable, dimension(:, :) :: ci_squared_initial
real(xp), allocatable, dimension(:)    :: d_ci_dr_inner, d_ci_dr_outer

! For call to rk4:
logical :: hit_target, target_met
real(8) :: t_target
! For irrelevant BC:
real(xp), parameter :: dummy = 0.d0

! For outputting many z planes if required:
logical :: plot_many_horiz_planes_only
integer :: nz_planes, nphi_planes
real(8), allocatable, dimension(:) :: z_list
integer, allocatable, dimension(:) :: iz_list, iphi_list

logical :: perturbation, compute_curl_rho_u

namelist /vsi_3D_input/ &
     plot_many_horiz_planes_only, &
     restart, file_version, isothermal, nr, nz, nphi, phi_max_over_pi, &
     cfl, dt_min, nsteps, &
     hit_target, t_target, &
     tecplot_interval, iphi_plot, output_profiles, profiles_interval, &
     output_phi_Reynolds_averages, time_interval_for_phi_Reynolds_averages, &
     rms_interval, &
     fluctuation_ke_interval, save_interval, &
     perturb, wavy_perturbation, &
     zmax_over_H0, &
     use_rsize_for_domain, rsize_over_H0, rmin_over_H0, &
     n_waves_in_r, &
     use_Manger_p, &
     apply_newtonian_cooling, &
     tau_cooling_over_t_Kepler_at_mid_radius, &
     apply_pade_filter, eps_or_tau, eps_filter, &
     tau_filter, apply_artificial_pressure, C_ap, apply_viscosity, viscosity_type, &
     C_DDSV, name_using_step, &
     rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, apply_fargo_trick, &
     apply_fargo_correction, &
     apply_sponge, sponge_type, rho1, rho2, d1, d2, tau_decay, n_decay_steps

if (my_node .eq. 0) then
   print *, ' node 0: First executable of vsi_3D'
   print *, ' my_node = ', my_node, ' num_nodes = ', num_nodes   
end if

! Read namelist input:
if (my_node .eq. 0) then
   print *, ' node 0: about to open and read namelist file for app_vsi_3D'      
   open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', &
        status = 'old')
   read (lun_general_purpose, nml = vsi_3D_input)
   close (lun_general_purpose)
end if

#ifdef mpi_code
   call mpi_bcast(plot_many_horiz_planes_only, 1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(restart,                   1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(file_version,              1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(isothermal,                1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(nr,                        1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nz,                        1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nphi,                      1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(phi_max_over_pi,           1, my_mpi_real,  0, mpi_comm_world, ier)
   
   call mpi_bcast(cfl,                       1, my_mpi_real,  0, mpi_comm_world, ier)
   call mpi_bcast(dt_min,                    1, my_mpi_real,  0, mpi_comm_world, ier)   
   call mpi_bcast(nsteps,                    1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(hit_target,                1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(t_target,                    1, my_mpi_real, 0, mpi_comm_world, ier)
   
   call mpi_bcast(tecplot_interval,          1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(iphi_plot,                 1, mpi_integer, 0, mpi_comm_world, ier)   

   call mpi_bcast(output_profiles,           1, mpi_logical, 0, mpi_comm_world, ier)   
   call mpi_bcast(profiles_interval,         1, mpi_integer, 0, mpi_comm_world, ier)

   call mpi_bcast(output_phi_Reynolds_averages, 1, mpi_logical, 0, mpi_comm_world, ier)   
   call mpi_bcast(time_interval_for_phi_Reynolds_averages, 1, my_mpi_real, &
        0, mpi_comm_world, ier)   
   
   call mpi_bcast(fluctuation_ke_interval,   1, mpi_integer, 0, mpi_comm_world, ier)   
   call mpi_bcast(save_interval,             1, mpi_integer, 0, mpi_comm_world, ier)
   
   call mpi_bcast(perturb,                   1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(wavy_perturbation,         1, mpi_logical, 0, mpi_comm_world, ier)
   
   call mpi_bcast(zmax_over_H0,              1, my_mpi_real,  0, mpi_comm_world, ier)
   call mpi_bcast(use_rsize_for_domain,      1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(rsize_over_H0,             1, my_mpi_real,  0, mpi_comm_world, ier)
   call mpi_bcast(rmin_over_H0,              1, my_mpi_real,  0, mpi_comm_world, ier)      
   call mpi_bcast(n_waves_in_r,              1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(use_Manger_p,              1, mpi_logical, 0, mpi_comm_world, ier)

   call mpi_bcast(apply_newtonian_cooling,                 1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(tau_cooling_over_t_Kepler_at_mid_radius, 1, my_mpi_real, 0, mpi_comm_world, ier)   

   call mpi_bcast(apply_pade_filter,         1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(eps_or_tau,                3, mpi_character, 0, mpi_comm_world, ier)
   call mpi_bcast(eps_filter,                1, my_mpi_real,  0, mpi_comm_world, ier)   
   call mpi_bcast(tau_filter,                1, my_mpi_real,  0, mpi_comm_world, ier)
   call mpi_bcast(apply_artificial_pressure, 1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(C_ap,                      1, my_mpi_real,  0, mpi_comm_world, ier)
   call mpi_bcast(apply_viscosity,           1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(viscosity_type,            1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(C_DDSV,                    1, my_mpi_real,  0, mpi_comm_world, ier)   
   call mpi_bcast(name_using_step,           1, mpi_logical, 0, mpi_comm_world, ier)

   call mpi_bcast(rmin_BC,                   1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(rmax_BC,                   1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(zmin_BC,                   1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(zmax_BC,                   1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(ibalanced,                 1, mpi_integer, 0, mpi_comm_world, ier)

   call mpi_bcast(apply_fargo_trick,         1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(apply_fargo_correction,    1, mpi_logical, 0, mpi_comm_world, ier)

   call mpi_bcast(apply_sponge,              1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(sponge_type,               7, mpi_character, 0, mpi_comm_world, ier)   
   call mpi_bcast(rho1,                      1, my_mpi_real,  0, mpi_comm_world, ier)
   call mpi_bcast(rho2,                      1, my_mpi_real,  0, mpi_comm_world, ier)
   call mpi_bcast(d1,                        1, my_mpi_real,  0, mpi_comm_world, ier)
   call mpi_bcast(d2,                        1, my_mpi_real,  0, mpi_comm_world, ier)
   call mpi_bcast(tau_decay,                 1, my_mpi_real,  0, mpi_comm_world, ier)
   call mpi_bcast(n_decay_steps,             1, mpi_integer,  0, mpi_comm_world, ier)   
#endif

if (my_node .eq. 1) then
   print *,' Output from node 1 to see if it has input variables'      
   print *,' plot_many_horiz_planes_only = ', plot_many_horiz_planes_only 
   print *,' restart                   = ', restart
   print *,' file_version              = ', file_version
   print *,' isothermal                = ', isothermal
   print *,' nr, nz, nphi              = ', nr, nz, nphi
   print *,' phi_max_over_pi           = ', phi_max_over_pi
   print *,' cfl                       = ', cfl
   print *,' dt_min                    = ', dt_min
   print *,' nsteps                    = ', nsteps
   print *,' hit_target                = ', hit_target
   print *,' t_target                  = ', t_target
   print *,' tecplot_interval          = ', tecplot_interval
   print *,' iphi_plot                 = ', iphi_plot
   print *,' profiles_interval         = ', profiles_interval
   print *,' fluctuation_ke_interval   = ', fluctuation_ke_interval
   print *,' save_interval             = ', save_interval
   print *,''
   print *,' perturb                   = ', perturb
   print *,' wavy_perturbation         = ', wavy_perturbation
   print *,''
   print *,' zmax_over_H0              = ', zmax_over_H0
   print *,''
   print *,' use_rsize_for_domain      = ',  use_rsize_for_domain
   print *,' rsize_over_H0             = ', rsize_over_H0
   print *,' n_waves_in_r              = ', n_waves_in_r
   print *,' use_Manger_p              = ', use_Manger_p
   print *,''
   print *,' apply_pade_filter         = ', apply_pade_filter
   print *,' eps_or_tau                = ', eps_or_tau
   print *,' eps_filter                = ', eps_filter
   print *,' tau_filter                = ', tau_filter
   print *,''
   print *,' apply_artificial_pressure = ', apply_artificial_pressure
   print *,' C_ap                      = ', C_ap
   print *,''
   print *,' apply_viscosity           = ', apply_viscosity
   print *,' viscosity_type            = ', viscosity_type
   print *,' C_DDSV                    = ', C_DDSV
   print *,''
   print *,' name_using_step           = ', name_using_step
   print *,''
   print *,' rmin_BC                   = ', rmin_BC
   print *,' rmax_BC                   = ', rmax_BC
   print *,' zmin_BC                   = ', zmin_BC
   print *,' zmax_BC                   = ', zmax_BC
   print *,' ibalanced                 = ', ibalanced
   print *,''
   print *,' apply_fargo_trick         = ', apply_fargo_trick
   print *,' apply_fargo_correction    = ', apply_fargo_correction
   print *,''
   print *,' apply_sponge              = ', apply_sponge
   print *,' rho1                      = ', rho1
   print *,' rho2                      = ', rho2
   print *,' tau_decay                 = ', tau_decay
   print *,' n_decay_steps             = ', n_decay_steps
end if

! The basic state follows Nelson et al. (2013).

! We can set three things to unity:
!GM    = 1.0_xp
pi         = 4.0*ATAN(1.0_xp)
!> Orbital period at mid-radius of the domain:   
T_orbital0 = 1.0_xp
H0         = 1.0_xp
rho0       = 1.0_xp

! This is a parameter:
H0_over_r0 = 0.10_xp
r0 = H0 / H0_over_r0

! Consequence of above:
Omega0 = 2._xp * pi / T_orbital0
GM     = Omega0**2 * r0**3
V_Kep0 = SQRT(GM/r0)
c0     = H0_over_r0 * V_Kep0

! Gravity set-up parameters:
gravity_flag = .true.

! Exponents:
if (use_Manger_p) then
   pexp = -2._xp/3._xp  ! for density.  Manger
else
   pexp = -3._xp/2._xp  ! for density Nelson
end if
qexp = -1.0_xp  ! for temperature
c_exponent = 0.5_xp*qexp
H_exponent = (qexp + 3._xp) / 2._xp

zmin = - zmax_over_H0 * H0
zmax =   zmax_over_H0 * H0

! Most amplified mode according to Matt:
lambda_r = pi * abs(qexp) * H0_over_r0 * H0

! Domain size in r:
if (.not. use_rsize_for_domain) then 
   rsize_over_H0 = n_waves_in_r * lambda_r
end if

!rmin   = r0 - 0.5_xp * rsize_over_H0
! This is to allow a computational domain that is not centered at r0:
rmin   = rmin_over_H0
rmax   = rmin + rsize_over_H0

phi_min = 0.0_xp
phi_max = phi_max_over_pi * pi

! These calls should be in every application subroutine:

call set_up_domain_mesh_and_partition(rmin, rmax, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
     suppress_z_derivatives_when_nz_not_1_arg, periodic_z_arg, &
     stretched_r, stretched_z, r0_unif, nr_u, z0, nz_u)

if (my_node .eq. 0) print *, ' app_vsi_3D: Returned from set-up_domain_mesh_and_partition'

! Note: We compute the basic state even for a restart run in order to save it for
! basic state subtraction for diagnostics.

! Functions of r for the basic state:
allocate(rho_midplane(nr), Omega_K(nr), ci(nr), H(nr))
do ir = 1, nr
   rho_midplane(ir) = rho0 * (rgrid(ir) / r0)**pexp
   Omega_K     (ir) = SQRT(GM/rgrid(ir)**3)
   ci          (ir) = c0 * (rgrid(ir)/r0)**c_exponent
   H           (ir) = H0 * (rgrid(ir)/r0)**H_exponent
end do

if (my_node .eq. 0) print *, ' app_vsi_3D: Finished with functions of r for basic state'

allocate(ci_squared_initial(nr, nz))
do iz = 1, nz
   do ir = 1, nr
      ci_squared_initial(ir, iz) = ci(ir)**2
   end do
end do

if (my_node .eq. 0) print *, ' app_vsi_3D: Finished with with ci_squared initial'

tau_newtonian_cooling = tau_cooling_over_t_Kepler_at_mid_radius * T_orbital0
call set_up_thermal_parameters(gamma, isothermal, apply_newtonian_cooling, &
           tau_newtonian_cooling, ci_squared_initial)

if (my_node .eq. 0) print *, ' node 0: Returned from set_up_thermal_parameters'

! Needed for isothermal centrifugally balanced non-reflective BC:
if (isothermal) then
   allocate(d_ci_dr_inner(nz), d_ci_dr_outer(nz))
   do iz = 1, nz
      d_ci_dr_inner(iz) = c0 * c_exponent / r0 * (rgrid(1 )/r0)**(c_exponent - 1._xp)
      d_ci_dr_outer(iz) = c0 * c_exponent / r0 * (rgrid(nr)/r0)**(c_exponent - 1._xp)
   end do
end if

call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, &
     d_ci_dr_inner, d_ci_dr_outer, c_sound_rmin, c_sound_rmax, isothermal)

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: Returned from set_up_boundary_conditions'
#endif   
   
gravity_type = mass_at_origin
gz_uniform   = 0._xp   
call activate_gravity(gravity_flag, gravity_type, GM, gz_uniform)   
if (my_node .eq. 0) print *, ' node 0: Returned from activate_gravity'

if (apply_fargo_trick) then
   ! The arguments are Boolean:
   call activate_fargo(integer_shifts, apply_fargo_correction)
   if (my_node .eq. 0) print *, ' node 0: Returned from activate_fargo'
end if

call activate_pade_filter(apply_pade_filter, eps_or_tau, eps_filter, &
     tau_filter, filter_relative_to_basic_state)
if (my_node .eq. 0) print *, ' node 0: Returned from activate_pade_filter'
if (apply_artificial_pressure) then
   if (my_node .eq. 0) print *, ' node 0: Calling activate_artificial_pressure with C_ap = ', C_ap
   if (my_node .eq. 1) print *, ' node 1: Calling activate_artificial_pressure with C_ap = ', C_ap      
   call activate_artificial_pressure (C_ap)
end if

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: app_vsi_3D: Returned from activate_artificial_pressure'
#endif

if (apply_viscosity) then
   C_Smag     = 0.22_xp**2
   nu_molecular   = 0._xp
   nu_b_molecular = 0._xp ! bulk viscosity
   Pr_molecular = 0._xp
   ! subroutine activate_viscosity(apply_viscosity_arg, viscosity_type_arg, isothermal_arg, &
   ! nu_molecular_arg, nu_b_molecular_arg, Pr_molecular_arg, gamma_arg, C_DDSV_arg, C_Smag_arg)
   
   call activate_viscosity(apply_viscosity, viscosity_type, isothermal, &
     nu_molecular, nu_b_molecular, Pr_molecular, gamma, C_DDSV, C_Smag)
end if

! ------------------------------------------------------
! Compute basic state/initial condition and put it in q:
! ------------------------------------------------------
! NOTE: This is now done for both fresh and restart runs.

! For an non-restart run q will be the initial condition.
! For BOTH restart and non-restart runs, part of q will
! be used to store the basic state.

! We compute the basic state so it can be used to output fluctuations.
! Previously we used to store the basic state

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: app_vsi_3D: About to allocate pressure'
#endif

allocate(pressure(sr:er, sphi:ephi, nz), dpdz(sr:er, sphi:ephi, nz))

! For numerical centrifugal balance:
allocate(p_r_space   (sphi:ephi, sz_r:ez_r, nr))
allocate(dpdr_r_space(sphi:ephi, sz_r:ez_r, nr))
allocate(dpdr        (sr:er, sphi:ephi, nz))

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         R_spherical = SQRT(rgrid(ir)**2 + zgrid(iz)**2)
         arg         = GM / ci(ir)**2 * (1._xp/R_spherical - 1._xp/rgrid(ir))
         rho_temporary          = rho_midplane(ir) * EXP(arg)
         pressure(ir, iphi, iz) = rho_temporary * ci(ir)**2
      end do
   end do
end do

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: app_vsi_3D: Finished calculating pressure'
#endif

! Initial condition/Basic state.
! For a restart run we need to go through the initial condition set-up
! in order to set-up the basic state needed for the sponge, and other things.   
call pade_diff_z(mr*mphi, pressure, dpdz)

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: app_vsi_3D: Returned from pade_diff_z'
#endif   

! Set rho to satisfy numerical vertical hydrostatic balance:
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! Required since gz = 0 at the midplane:
         if (gz(ir, iz) .ne. 0._xp) then
            q(ir, iphi, iz, irho) = dpdz(ir, iphi, iz)/gz(ir, iz)
         else
            q(ir, iphi, iz, irho) = rho_midplane(ir)
         end if
         ! Correct the pressure:
         pressure(ir, iphi, iz) = q(ir,iphi,iz,irho) * ci(ir)**2
      end do
   end do
end do

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: app_vsi_3D: Finished correcting the pressure'
#endif   

! Set uphi to satisfy numerical centrifugal balance:
call transpose_z_to_r (1, pressure, p_r_space)
call pade_diff_bundle(mphi*mz_r, nr, Ji_r, p_r_space, dpdr_r_space)
call transpose_r_to_z (1, dpdr_r_space, dpdr)   
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         q(ir, iphi, iz, rmom) = 0._xp
         q(ir, iphi, iz, zmom) = 0.0_xp

         ! term1 = (p + qexp) * (H(ir)/rgrid(ir))**2
         ! term2 = 1._xp + qexp
         ! R_spherical = SQRT(rgrid(ir)**2 + zgrid(iz)**2)
         ! term3 = - qexp*rgrid(ir) / R_spherical
         ! Omega = Omega_K(ir) * (term1 + term2 + term3)**(0.5_xp)
         ! uphi  = Omega * rgrid(ir)

         ! Look at NRR-3 notes:
         uphi_squared = rgrid(ir)*(1._xp/q(ir,iphi,iz,irho)*dpdr(ir,iphi,iz) - gr(ir,iz))
         q(ir, iphi, iz, amom) = q(ir, iphi, iz, irho) * sqrt(uphi_squared) * rgrid(ir)

         ! This is provisional and will be set below for the non-isothermal case:
         q(ir, iphi, iz, ener) = 0._xp
      end do
   end do
end do

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: app_vsi_3D: Finished setting uphi to satisfy numerical balance'
#endif   

if (.not. isothermal) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, iz, ener) = pressure(ir, iphi, iz) / (gamma - 1._xp)
         end do
      end do
   end do
end if
! Finished with basic state

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: app_vsi_3D: Finished with basic state'
#endif   

! This is now read from the input file.   
!iphi_plot    = 1
if (nphi .eq. 1) then
   ! For an axisymmetric run we don't want to plot a horizontal plane:
   iz_plot_midplane = 0
else
   iz_plot_midplane = iz_mid
end if


! Output basic state vorticity:
istep = 999
t = 0.d0
plot_curl_rho_u = .false.
plot_pert = .false.
call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, &
     iz_plot_midplane, L_scale, T_scale, rho_scale, t, istep)

! Store what is in q as the basic state.  This is done for both
! fresh and restart runs.
call store_basic_state(filter_relative_to_basic_state)

deallocate(pressure, dpdz, p_r_space, dpdr_r_space, dpdr)

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: app_vsi_3D: Finished deallocating pressure'
if (my_node .eq. 1) print *, ' node 0: app_vsi_3D: Finished deallocating pressure'
#endif

if (restart) then
   call read_restart (file_version, istep_of_restart_file, t_of_restart_file)
else
   ! Assign initial condition using the basic state:
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, iz, irho) = rho_basic(ir, iz)
            q(ir, iphi, iz, rmom) = 0._xp
            q(ir, iphi, iz, zmom) = 0._xp                        
            q(ir, iphi, iz, amom) = rho_basic(ir, iz) * uphi_basic(ir, iz) * rgrid(ir)

            ! Set energy in case we run adiabatically:
            ! q(ir, iphi, iz, ener) = pressure(ir,iphi,iz)/gm1
            q(ir, iphi, iz, ener) = 0._xp
         end do
      end do
   end do
end if
! Finished with initial condition.

if (apply_sponge) then
   call activate_sponge(sponge_type, rho1, rho2, rho_basic, d1, d2, tau_decay, n_decay_steps)   
end if

! Needed for the isothermal or Newtonian cooling options.  ci_squared_initial sits in
! module thermal_parameters.
if (isothermal .or. apply_newtonian_cooling) then
   do iz = 1, nz
      do ir = 1, nr
         ci_squared_initial(ir, iz) = ci(ir)**2
      end do
   end do
end if

if (.not. restart) then
   istep0 = 0
   t      = 0.0_xp
else
   istep0 = istep_of_restart_file
   t      = t_of_restart_file
end if

if (perturb) then
   if (wavy_perturbation) then
      if (my_node .eq. 0) print *, ' app_vsi_3D: Calling add_wavy_perturbation'
      call add_wavy_perturbation(c0, H0_over_r0, H0, lambda_r)
      if (my_node .eq. 0) print *, ' app_vsi_3D: Returned from wavy_perturbation'
   else
      call add_random_perturbation(c0)
      if (my_node .eq. 0) print *, ' app_vsi_3D: Returned from add_random_perturbation'      
   end if
end if

! These routines are currently only serial.
#ifndef mpi_code
   ! call output_gravity_profile_at_mid_radius
   ! call vertical_profiles(ir_mid, 1, 'midr', t) ! mid ir
   ! call radial_profiles  (iz_mid, 1, 'midz', t) ! midplane
   ! call radial_profiles  (1,      1, 'botz', t) ! bottom
   ! call radial_profiles  (nz,     1, 'topz', t) ! top
#endif

call sanity_check

if (my_node .eq. 0) print *, ' node 0: in vsi_3D Mb allocated = ', float(8*n_words_allocated)/1.d6

if (my_node .eq. 0) then
   print *, ' app_vsi_3D: Calling tecplot_vort_and_dil_in_merid_horiz_planes'
end if

!plot_curl_rho_u = .true.
!plot_pert = .true.
!call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot_midplane, L_scale, T_scale, rho_scale, t, istep0)

if (plot_many_horiz_planes_only) then
   nz_planes   = 21
   nphi_planes = 0   
   allocate(z_list(nz_planes), iz_list(nz_planes))
   ! For some reason I have to allocate this even if I don't use it:
   allocate(iphi_list(nphi_planes))   
   z_list = (/ -1.0d0, -0.9d0, -0.8d0, -0.7d0, -0.6d0, -0.5d0, -0.4d0, -0.3d0, -0.2d0, -0.1d0, 0.d0, &
               0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, 0.9d0, 1.0d0/)
   call make_iz_list(nz_planes, z_list, iz_list)
   call tecplot_many_horizontal_planes(q, nz_planes, iz_list, z_list, T_scale, t, istep0)

   perturbation = .true.
   compute_curl_rho_u = .false.
   call tecplot_vort_and_dil_in_many_merid_horiz_planes(perturbation, &
     compute_curl_rho_u, q, nphi_planes, iphi_list, nz_planes, iz_list, &
     z_list, L_scale, T_scale, rho_scale, t, istep0)
   
   call terminate_with_no_save(0)
end if

call many_tecplot(istep0, iphi_plot, iz_plot_midplane, t, &
                  L_scale, T_scale, rho_scale)

if (output_phi_Reynolds_averages) then
   call phi_Reynolds_averages(istep0, t, average_in_z)
   next_t_for_phi_Reynolds_averages = t + time_interval_for_phi_Reynolds_averages
end if

call history_write_ave_fluctuation_ke(t)
call history_write_rms(t, output_profiles)

! Assume normal return.  Will be changed if we hit a snag:
istatus = 0
do istep = istep0 + 1, istep0 + nsteps
   call rk4(rhs, q, cfl, t, dt, use_supplied_dt, unphysical_flag, iphi_bad, &
                       hit_target, t_target, target_met)
   call output_conservation_diagnostics(t, q)
   if (hit_target .and. target_met) then
      istatus = 2
      exit
   end if
   
   if (my_node .eq. 0) then
      print *, ' routine app_vsi_3D, node 0: finished istep = ', istep, ' t = ', t
   end if

   !call print_max(q)

   if (mod(istep - istep0, tecplot_interval) .eq. 0) &
        call many_tecplot(istep, iphi_plot, iz_plot_midplane, t, &
                        L_scale, T_scale, rho_scale)
  
   if ((output_phi_Reynolds_averages) .and. (t .ge. next_t_for_phi_Reynolds_averages)) then
      call phi_Reynolds_averages(istep, t, average_in_z)
      next_t_for_phi_Reynolds_averages = t + time_interval_for_phi_Reynolds_averages
   end if   

   if (mod(istep - istep0, fluctuation_ke_interval) .eq. 0) then
      call history_write_ave_fluctuation_ke(t)
      call history_write_rms(t, output_profiles)
   end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node = 0, vsi_3D: about to check for a save'
#endif

   if (mod(istep - istep0, save_interval) .eq. 0) then
      call write_save_file(istep, t, filename_given)
   end if

   ! Detect too small time step:
   if (istep .gt. istep0+1) then
      if (dt .lt. dt_min) then
         if (my_node .eq. 0) then
            print *, ' my_node = ', my_node
            print *, ' dt < dt_min'
            print *, ' dt          = ', dt
            print *, ' dt_min      = ', dt_min
            print *, ' dt_previous = ', dt_previous
         end if
      end if
   end if

   if ((dt .lt. dt_min) .or. (unphysical_flag)) then
      istatus = 1 ! Flag for non-normal return
      exit
   end if
   dt_previous = dt
end do

if (my_node .eq. 0) then
   print *, ' ***** average stepping time per step = ', total_cpu_for_stepping/nsteps
end if

close(lun_history(1))
call many_tecplot(istep, iphi_plot, iz_plot_midplane, t, &
                  L_scale, T_scale, rho_scale)
! istatus:
! 0 : Normal return. Finished nsteps
! 1 : Non-normal return.  Either unphysical or dt < dt_min
! 2 : hit_target = .true. and target_met = .true.  We set this to non-zero
!     to prevent an auto restart by the pbs script.

! Need the istep - 1 since a fortran do loop increments the counter at the end
call terminate_with_save(istatus, istep-1, t)

end subroutine app_vsi_3D

!----------------------------------------------------------------------------------85

subroutine add_wavy_perturbation(c0, H0_over_r0, H0, lambda_r)

use precision
use q_array
use grid
use dof_indices
use partition_data
#ifdef mpi_code
   use mpi, only: mpi_comm_world
#endif

! ifort compatibility library so we can use the same random number generator as gfortran.   
#ifdef ifort
   use IFPORT 
#endif
   
implicit none
real(xp), intent(in) :: c0, H0_over_r0, H0, lambda_r
integer, parameter :: nkr = 11, nkz = 12, nkphi_dim = 12
integer :: nkphi
real(xp) :: kr(nkr), kz(nkz), kphi(nkphi_dim)
! Local:
integer :: ir, iz, iphi, ikr, ikz, ikphi, ier
real(xp) :: pi, eps, lambda_z

real(xp) :: lambda_r_modulation, kr_modulation, lambda_z_modulation, kz_modulation, &
     r_modulation, z_modulation, modulation
real(xp) :: rho, uz, ur, uphi
real(xp) :: uz_prime, ur_prime, uphi_prime
real(xp) :: phase_ur, phase_uz, phase_uphi
complex(8) :: exp_arg, exp_fac, ic = CMPLX(0.0_xp, 1.0_xp)

complex(8), allocatable, dimension(:,:,:) :: A_ur, A_uz, A_uphi

if (my_node .eq. 0) print *, ' add_wavy_perturbation has been called'

if (nphi .eq. 1) then
   nkphi = 1
else
   nkphi = nkphi_dim
end if

allocate(A_ur(nkphi,nkz,nkr), A_uz(nkphi,nkz,nkr), A_uphi(nkphi,nkz,nkr))

! Establish seed for random number generator.  Note that each processor will have the same
! seed and hence the same phase:
call srand(987654)
pi = 4.0_xp * atan(1.0_xp)
!eps = 1.e-3_xp*c0
eps = 0.001_xp*c0

! The purpose of the modulation is to make normal velocities = 0 at boundaries.
! The modulation is a half sin.
lambda_r_modulation = 2._xp * (rgrid(nr) - rgrid(1))
kr_modulation       = 2._xp * pi / lambda_r_modulation
lambda_z_modulation = 2._xp * (zgrid(nz) - zgrid(1))
kz_modulation       = 2._xp * pi / lambda_z_modulation

! Note: nkr = 11
kr(1) = 2._xp * pi / lambda_r
do ikr = 2, 8
   kr(ikr) = ikr * kr(1)
end do
! Sub-harmonics:
kr(9)  = kr(1)/2._xp
kr(10) = kr(1)/3._xp
kr(11) = kr(1)/4._xp

lambda_z = zgrid(nz) - zgrid(1)
kz(1) = 2._xp * pi / lambda_z
do ikz = 2, nkz
   kz(ikz) = ikz * kz(1)
end do

do ikphi = 1, nkphi
   kphi(ikphi) = ikphi - 1
end do

if (my_node .eq. 0) then
   print *, ' lambda_r = ', lambda_r
   print *, ' # of grid points per lambda_r = ', nr * lambda_r / (rmax - rmin)
   print *, ' # of waves in r domain = ', (rgrid(nr) - rgrid(1)) / lambda_r
end if

#ifdef mpi_code
   call mpi_barrier(mpi_comm_world, ier)
#endif

! Note that this will create jumps between processors.

! Get complex amplitudes for each mode:
do ikr = 1, nkr
   do ikz = 1, nkz
      do ikphi = 1, nkphi
         ! Zero argument means next random number in the sequence.
         phase_ur   = rand(0) * 2._xp * pi                  
         phase_uz   = rand(0) * 2._xp * pi
         phase_uphi = rand(0) * 2._xp * pi

         A_uz  (ikphi, ikz, ikr) = cos(phase_uz)   + ic*sin(phase_uz)
         A_ur  (ikphi, ikz, ikr) = cos(phase_ur)   + ic*sin(phase_ur)
         A_uphi(ikphi, ikz, ikr) = cos(phase_uphi) + ic*sin(phase_uphi)
      end do
   end do
end do
   
do iz = 1, nz
   ! print *, ' in point loop iz = ', iz
   ! print *, ' sphi = ', sphi, ' ephi = ', ephi
   z_modulation = sin(kz_modulation*(zgrid(iz) - zgrid(1)))
   do iphi = sphi, ephi
      do ir = sr, er
         ! r_modulation = sin(kr_modulation*(rgrid(ir) - rgrid(1)))
         ! modulation = r_modulation * z_modulation
         modulation = z_modulation
         
         rho  = q(ir, iphi, iz, irho)
         uz   = q(ir, iphi, iz, zmom) / rho
         ur   = q(ir, iphi, iz, rmom) / rho
         uphi = q(ir, iphi, iz, amom) / rho / rgrid(ir)

         ! Sum over modes:
         uz_prime   = 0.0_xp
         ur_prime   = 0.0_xp
         uphi_prime = 0.0_xp
         do ikr = 1, nkr
            do ikz = 1, nkz
               do ikphi = 1, nkphi
                  exp_arg = ic * (kr(ikr)*rgrid(ir) + kz(ikz)*zgrid(iz) + &
                       kphi(ikphi)*phi_grid(iphi))
                  exp_fac = EXP(exp_arg)
                  uz_prime   = uz_prime   + A_uz  (ikphi,ikz,ikr) * exp_fac
                  ur_prime   = ur_prime   + A_ur  (ikphi,ikz,ikr) * exp_fac
                  uphi_prime = uphi_prime + A_uphi(ikphi,ikz,ikr) * exp_fac
                  ! print *, ' phase_uz = ', phase_uz
                  ! print *, ' A_uz = ', A_uz, ' pi = ', pi
               end do
            end do
         end do ! sum over modes
         uz_prime   = eps * modulation * uz_prime
         ur_prime   = eps * modulation * ur_prime
         uphi_prime = eps * modulation * uphi_prime

         q(ir, iphi, iz, zmom) = rho*(uz   + uz_prime)
         q(ir, iphi, iz, rmom) = rho*(ur   + ur_prime)
         q(ir, iphi, iz, amom) = rho*(uphi + uphi_prime)*rgrid(ir)
      end do
   end do
end do

end subroutine add_wavy_perturbation

!----------------------------------------------------------------------------------85

subroutine add_random_perturbation(c0)

use precision, only:xp
use q_array
use grid
use dof_indices
use partition_data
use math_constants
#ifdef mpi_code
   use mpi, only: mpi_comm_world
#endif

! ifort compatibility library so we can use the same random number generator as gfortran.   
#ifdef ifort
   use IFPORT 
#endif
   
implicit none
real(xp), intent(in) :: c0
integer :: ir, iz, iphi, ier
real(xp) :: eps

real(xp) :: lambda_r_modulation, kr_modulation, lambda_z_modulation, kz_modulation, &
     r_modulation, z_modulation, modulation
real(xp) :: rho, uz, ur, uphi
real(xp) :: uz_prime, ur_prime, uphi_prime

if (my_node .eq. 0) print *, ' random_perturbation has been called'

! Establish seed for random number generator:
call srand(987654 + 789*my_node)
eps = 1.d-3*c0

! The purpose of the modulation is to make normal velocities = 0 at boundaries.
! The modulation is a half sin.
lambda_r_modulation = 2._xp * (rgrid(nr) - rgrid(1))
kr_modulation       = 2._xp * pi / lambda_r_modulation
lambda_z_modulation = 2._xp * (zgrid(nz) - zgrid(1))
kz_modulation       = 2._xp * pi / lambda_z_modulation

#ifdef mpi_code
   call mpi_barrier(mpi_comm_world, ier)
#endif

do iz = 1, nz, 2
   ! print *, ' in point loop iz = ', iz
   ! print *, ' sphi = ', sphi, ' ephi = ', ephi
   z_modulation = sin(kz_modulation*(zgrid(iz) - zgrid(1)))
   do iphi = sphi, ephi, 2
      do ir = sr, er, 2
         r_modulation = sin(kr_modulation*(rgrid(ir) - rgrid(1)))
         modulation = r_modulation * z_modulation
         modulation = z_modulation
         modulation = 1._xp
         
         rho  = q(ir, iphi, iz, irho)
         uz   = q(ir, iphi, iz, zmom) / rho
         ur   = q(ir, iphi, iz, rmom) / rho
         uphi = q(ir, iphi, iz, amom) / rho / rgrid(ir)

         uz_prime   = eps*rand(0)*modulation
         ur_prime   = eps*rand(0)*modulation
         uphi_prime = eps*rand(0)*modulation

         q(ir, iphi, iz, zmom) = rho*(uz   + uz_prime)
         q(ir, iphi, iz, rmom) = rho*(ur   + ur_prime)
         q(ir, iphi, iz, amom) = rho*(uphi + uphi_prime)*rgrid(ir)
      end do
   end do
end do

end subroutine add_random_perturbation

!----------------------------------------------------------------------------------85

subroutine many_tecplot(istep, iphi_plot, iz_plot_midplane, t, &
                        L_scale, T_scale, rho_scale)

use precision
use q_array
use grid
implicit none
integer, intent(in) :: istep, iphi_plot, iz_plot_midplane
real(8), intent(in) :: t, L_scale, T_scale, rho_scale

! Local:
logical, parameter :: plot_curl_rho_u = .false.
logical :: plot_pert

call tecplot_meridional_plane (q, iphi_plot, L_scale, T_scale, t, istep)
if (nphi .ne. 1) call tecplot_horizontal_plane (q, iz_plot_midplane, T_scale, t, istep)

plot_pert = .true.
call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, &
     iz_plot_midplane, L_scale, T_scale, rho_scale, t, istep)
call tecplot_fluc_vel_meridional_plane(q, iphi_plot, L_scale, T_scale, t, istep)

end subroutine many_tecplot

!----------------------------------------------------------------------------------85





