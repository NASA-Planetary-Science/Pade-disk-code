!----------------------------------------------------------------------------------85

subroutine app_hydrostatic_test

! Hydrostatic in z with only one grid point in r and phi.  There is no azimuthal
! velocity.  NOTE: This is an interactive and serial routine.

! To run without an acoustic pulse set eps = 0.d0 in the input file.

!use hydrostatic_parameters
use grid, only: nz, nr, nphi, ndof, rmin, rmax, zmin, zmax, phi_min, phi_max, zgrid, rgrid
use q_array
use boundary_condition_types
use boundary_condition_routines
use gravity_types
use logical_units
use set_up_routines
use activate_routines
use basic_state
implicit none

! Locals:
external rhs
real(8), allocatable, dimension(:) :: rho_bar, p_bar

! For domain (read from namelist):
real(8) :: zmax_over_H

! For passing to thermal parameters:
real(8) :: gamma, tau_newtonian_cooling
logical :: isothermal, apply_newtonian_cooling
real(8), allocatable, dimension(:, :) :: ci_squared_initial

! For local use:
real(8) :: gm1 ! gamma - 1

! For subroutine activate_gravity:
logical :: gravity_flag
integer :: gravity_type
real(8) :: gz_uniform
real(8) :: GM

logical :: hit_target, target_met
integer :: istep, i_output_step, n_output_steps
real(8) :: dt_output, t_target
real(8) :: cfl, t, dt, Omega, r0
integer :: ir, iz, ir_mid

real(8) :: eps, sigma, z0

! Isothermal and adiabatic sound speeds:
real(8) :: ci, ca

real(8) :: rho0, H
logical :: suppress_z_derivatives_when_nz_not_1 = .false., periodic_z = .false., &
     use_supplied_dt = .false.

! For restart if needed:
logical :: restart
integer :: istep_of_restart_file, file_version
real(8) :: t_of_restart_file

integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced
real(8), allocatable :: d_ci_dr_inner(:), d_ci_dr_outer(:)
real(8) :: c_sound_rmin, c_sound_rmax

character(3) :: zmin_BC_string, zmax_BC_string

logical :: apply_pade_filter, filter_relative_to_basic_state = .true.
character(3) :: eps_or_tau
real(8) :: eps_filter, tau_filter

logical :: apply_artificial_pressure
real(8) :: vke, dt_old, C_ap

logical, parameter :: stretched_r = .false., stretched_z = .false.
real(8), parameter :: r0_unif = 0.d0, z0_unif = 0.d0
integer, parameter :: nr_u = 0, nz_u = 0
logical :: unphysical
integer :: iphi_bad
! For irrelevant arguments: 
real(8) :: dummy

namelist /hydrostatic_test_input/ gravity_type, isothermal, zmin_BC_string, zmax_BC_string, &
     nz, cfl, dt_output, n_output_steps, apply_pade_filter, eps_or_tau, eps_filter, tau_filter, &
     eps, sigma, z0, zmax_over_H, apply_artificial_pressure, C_ap

open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', &
      status = 'old')
read (lun_general_purpose, nml = hydrostatic_test_input)
close (lun_general_purpose)

print *, ' Some parameters from the namelist input file:'
print *, ' gravity_type = ', gravity_type
print *, ' isothermal   = ', isothermal
print *, ' zmin_BC_string = ', zmin_BC_string
print *, ' zmax_BC_string = ', zmax_BC_string

print *, ' nz = ', nz
print *, ' pulse amplitude, eps = ', eps

print *, ' gravity types are'
print *, ' mass_at_origin = 1'
print *, ' thin_disk_mass_at_origin = 2'
print *, ' thin_disk_mass_at_origin_no_radial = 3'
print *, ' uniform_gz = 4'

print *, ' Enter anything to continue'
read(5, *)

allocate(rho_bar(nz), p_bar(nz))
! In case we run isothermally:
allocate(ci_squared_initial(1, nz))

! Not needed but just to be clean:
allocate(d_ci_dr_inner(nz), d_ci_dr_outer(nz))

gamma      = 1.4d0
gm1        = gamma - 1.d0

! This is the radius at which the calculation will be performed.
r0 = 1.d0

! For units can set three things to unity:

if (gravity_type .eq. uniform_gz) then
   gz_uniform = -1.d0 ! Uniform gravity pointing down
   ci         =  1.d0 ! Isothermal sound speed:
   rho0       =  1.d0 ! Midplane density

   ! Scale height:
   H = ci**2 / abs(gz_uniform)

   zmin = 0.d0
   zmax = zmax_over_H*H

   print *, ' zmin = ', zmin
   print *, ' zmax = ', zmax
   print *, ' Enter anything to continue'
   read(5, *)
else if (gravity_type .eq. thin_disk_mass_at_origin_no_radial) then
   GM   = 1.d0
   ci   = 1.d0
   rho0 = 1.d0

   Omega = SQRT(GM/r0**3)

   ! Scale height:
   H = ci/Omega

   zmin = -zmax_over_H*H
   zmax =  zmax_over_H*H   
end if

ca    = SQRT(gamma) * ci   ! adiabatic sound speed

! Needed for the isothermal option.  ci_squared_initial is passed to
! set_up_thermal_parameters.
if (isothermal) then
   do iz = 1, nz
      do ir = 1, nr
         ci_squared_initial(ir, iz) = ci**2
      end do
   end do
end if

rmin    =  r0
rmax    =  r0
phi_min = 0.0d0
phi_max = 0.0d0
nr      = 1
nphi    = 1

gravity_flag = .true.
restart = .false.
file_version = 2

! Boundary conditions:
rmin_BC   = null
rmax_BC   = null

if (zmin_BC_string .eq. 'znm') then
   zmin_BC = zero_normal_momentum
else if (zmin_BC_string .eq. 'nrf') then
   zmin_BC = non_reflective
else
   print *, ' Unrecognized zmin_BC_string = ', zmin_BC_string
end if

if (zmax_BC_string .eq. 'znm') then
   zmax_BC = zero_normal_momentum
else if (zmin_BC_string .eq. 'nrf') then
   zmax_BC = non_reflective
else
   print *, ' Unrecognized zmax_BC_string = ', zmax_BC_string
end if

ibalanced = 1

! These three set-up calls should be in every application subroutine:
call set_up_domain_mesh_and_partition(rmin, rmax, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
     suppress_z_derivatives_when_nz_not_1, periodic_z, &
     stretched_r, stretched_z, r0_unif, nr_u, z0_unif, nz_u)
apply_newtonian_cooling = .false.
tau_newtonian_cooling   = 0.d0
call set_up_thermal_parameters(gamma, isothermal, apply_newtonian_cooling, &
                                     tau_newtonian_cooling, ci_squared_initial)

call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, &
     d_ci_dr_inner, d_ci_dr_outer, c_sound_rmin, c_sound_rmax, isothermal)

call activate_gravity(gravity_flag, gravity_type, GM, gz_uniform)

if (apply_pade_filter) then
   call activate_pade_filter(apply_pade_filter, eps_or_tau, eps_filter, &
        tau_filter, filter_relative_to_basic_state)
end if

if (apply_artificial_pressure) then
   call activate_artificial_pressure (C_ap)
end if

if (restart) then
   call read_restart(file_version, istep_of_restart_file, t_of_restart_file)
end if
   
#ifdef debug_print
   print *, ' Finished initializating ci_squared_initial'
#endif

call initial_condition

#ifdef debug_print
   print *, ' Returned from initial condition'
#endif
      
call store_basic_state(filter_relative_to_basic_state)   
call output_gravity_profile_at_mid_radius   

#ifdef debug_print
   print *, ' Finished hydrostatic initial condition'
#endif

t = 0.0d0
call output
open (unit = 101, file = 'vke.dat', form = 'formatted', status = 'unknown')

hit_target = .true.
t_target   = t + dt_output
target_met = .false.
istep      = 0
do i_output_step = 1, n_output_steps
   do while (.not. target_met)
      call rk4(rhs, q, cfl, t, dt, use_supplied_dt, unphysical, iphi_bad, &
           hit_target, t_target, target_met)
      if (unphysical) call terminate_with_no_save(1)   
      istep = istep + 1
      print *, ' finished rk4 istep = ', istep, ' dt = ', dt, ' t = ', t
      print *, ' target_met = ', target_met, ' t_target = ', t_target

      call vertical_kinetic_energy (q, vke)
      write (101, "(2(1x, e12.5))") t, vke
   end do

   ! We just met the target time:
   call output   

   t_target = t_target + dt_output ! Next target
   target_met = .false.
end do

close(101)
call terminate_with_no_save(0)

contains

!----------------------------------------------------------------------------------85

subroutine output

use dof_indices
implicit none
character(40), dimension(5) :: filename
integer       :: iz, ifile
real(8)       :: rho, uz, p, ci

write (filename(1), "('rho_t=',   i3.3, f0.4, '.dat')") int(t), t - int(t)
write (filename(2), "('uz_t=',    i3.3, f0.4, '.dat')") int(t), t - int(t)
write (filename(3), "('p_t=',     i3.3, f0.4, '.dat')") int(t), t - int(t)
write (filename(4), "('rhop_t=',  i3.3, f0.4, '.dat')") int(t), t - int(t)

do ifile = 1, 4
   open (unit = 10+ifile, file = filename(ifile), form = 'formatted', status = 'unknown')
end do

do iz = 1, nz
   rho       = q(1,1, iz, irho)
   uz        = q(1,1, iz, zmom) / rho

   if (isothermal) then
      p = rho * ci**2
   else
      p = q(1,1, iz, ener) * gm1
   end if
   
   write (11, "(2(1x, e13.5e3))") zgrid(iz), rho
   write (12, "(2(1x, e13.5e3))") zgrid(iz), uz
   write (13, "(2(1x, e13.5e3))") zgrid(iz), p
   write (14, "(2(1x, e13.5e3))") zgrid(iz), rho - rho_bar(iz)     

   if (q(1, 1, iz, amom) .ne. 0.d0) then
      print *, ' angular momentum .ne. 0 for iz = ', iz
      stop
   end if

   if (q(1, 1, iz, rmom) .ne. 0.d0) then
      print *, ' radial momentum .ne. 0 for iz = ', iz
      stop
   end if
end do
do ifile = 1, 4
   close(10+ifile)
end do
end subroutine output

!~~~~~~~~~~~~~~~
subroutine initial_condition

use dof_indices, only: irho, rmom, amom, zmom, ener
use grid
use gravity
use thermal_parameters
implicit none

! Local:
integer :: ir, iphi, iz
real(8) :: rho_provisional, c2, N2, pi
real(8), dimension(nz) :: dpdz, u_pert, theta, d_theta_dz, dpdz_exact

real(8) :: dpdz_err, rho_err, c2_err

real(8) :: z_shift

! Relative density perturbation in the notation of Liepmann and Roshko:
real(8):: s_tilde

! For quantities after acoustic perturbation has been added:
real(8) :: rho1, p1, uz1

if (gravity_type .eq. uniform_gz) then
   do iz = 1, nz
      rho_provisional = rho0 * EXP(-zgrid(iz) / H)
      p_bar(iz)       = rho_provisional * ci**2
      dpdz_exact(iz)  = -p_bar(iz)/ H
   end do
else if (gravity_type .eq. thin_disk_mass_at_origin_no_radial) then
   do iz = 1, nz
      rho_provisional = rho0 * EXP(-zgrid(iz)**2 / (2.d0*H**2)  )
      p_bar(iz)       = rho_provisional * ci**2
      dpdz_exact(iz)  = p_bar(iz) * (-2.d0 * zgrid(iz) / (2.d0*H**2)   )
   end do
end if

call pade_diff_z(nr*nphi, p_bar, dpdz)

! Set rho to satisfy numerical hydrostatic balance:
do iz = 1, nz
   ! Required since gz = 0 at the midplane for the thin disk gravity since
   ! it equals Omega^2 *z:
   if (gz(1, iz) .ne. 0.d0) then
      rho_bar(iz) = dpdz(iz)/gz(1, iz)
   else
      rho_bar(iz) = rho0
   end if
end do

open (unit = 1, file = 'dp_dz_error.dat', form = 'formatted', status = 'unknown')
do iz = 1, nz
   dpdz_err = ABS(dpdz(iz) - dpdz_exact(iz))
   write (1, "(2(1x, e12.5))") zgrid(iz), dpdz_err
end do
close (1)

do iz = 1, nz
   do iphi = 1, nphi
      do ir = 1, nr
         z_shift = zgrid(iz) - z0
         s_tilde = eps*exp(-z_shift**2/sigma**2)
         if (isothermal) then
            uz1 = ci * s_tilde
         else
            uz1 = ca * s_tilde
         end if
         rho1 = rho_bar(iz) * (1.d0 + s_tilde)
         p1   =   p_bar(iz) * (1.d0 + s_tilde)
         q(ir, iphi, iz, irho) = rho1
         q(ir, iphi, iz, rmom) = 0.d0
         q(ir, iphi, iz, amom) = 0.d0
         q(ir, iphi, iz, zmom) = rho1 * uz1

         ! In case we are adiabatic:
         q(ir, iphi, iz, ener) = p1 / gm1
      end do
   end do
end do

print *, ' Finished with initial condition'
print *, ' Enter anything to continue'
read(5, *)

! Calculate the Brunt-Vaisala frequency:
do iz = 1, nz
   c2 = gamma * p_bar(iz) / rho_bar(iz)
   ! This is proportional to the potential temperature:
   theta(iz) = c2 * p_bar(iz)**(gamma/gm1)
end do

print *,' calling pade_diff_bundle'
call pade_diff_z(1, theta, d_theta_dz)
print *, ' returned from pade_diff_bundle'

open (unit = 1, file = 'Brunt_Vaisala_Squared.dat', form = 'formatted', &
     status = 'unknown')
do iz = 1, nz
   ! Square of the Brunt-Vaisala frequency:
   N2 = gz(1, iz) / theta(iz) * d_theta_dz(iz)
   write (1, "(2(1x, e12.5))") zgrid(iz), N2
end do
close (1)

return
end subroutine initial_condition
!~~~~~~~

end subroutine app_hydrostatic_test

!----------------------------------------------------------------------------------85

subroutine vertical_kinetic_energy (q, vke)

use partition_data
use grid
use dof_indices
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: vke

! Local:
real(8) :: integrand_1, integrand_2
integer :: iz

vke = 0.0d0
do iz = 1, nz-1
   integrand_1 = q(1, 1, iz,   zmom)**2 / q(1, 1, iz,   irho)*Jz(iz)
   integrand_2 = q(1, 1, iz+1, zmom)**2 / q(1, 1, iz+1, irho)*Jz(iz+1)
   vke = 0.5d0 * (integrand_1 + integrand_2)
end do
vke = 0.5d0*vke

end subroutine vertical_kinetic_energy

!----------------------------------------------------------------------------------85

subroutine output_gravity_profile_at_mid_radius

! This is a serial routine.

use grid
use logical_units
use gravity
implicit none

! Local:
integer :: iz

open (unit = lun_profile(1), file = 'gz_gr_vs_z.dat', form = 'formatted', &
     status = 'unknown')
do iz = 1, nz
   write (lun_profile(1), 1) zgrid(iz), gz(ir_mid, iz), gr(ir_mid, iz)
   1 format (3(1x, e12.5))
end do

close (lun_profile(1))
end

!----------------------------------------------------------------------------------85




