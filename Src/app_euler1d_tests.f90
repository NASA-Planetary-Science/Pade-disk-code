!----------------------------------------------------------------------------------85

module acoustic_pulse_parameters
   implicit none
   real(8) :: p0, rho0, z0, eps, sigma
end module acoustic_pulse_parameters

module shock_tube_parameters
   implicit none
   ! Ms is the shock Mach number.  It is needed for the exact solution and is
   ! obtained by initializing subroutine exact_shock_tube.
   real(8) :: rhoL, pL, uL, rhoR, pR_Sod, uR, z0_shock_tube_actual, Ms

   ! Computed in exact_shock_tube for passing to function compatibility during
   ! initialization of the exact solution.
   real(8) :: aL, aR   
end module shock_tube_parameters

module problem_type_module
   integer, parameter :: acoustic = 1, shock_tube = 2, density_waves = 3
   integer :: problem_type
end module problem_type_module

!----------------------------------------------------------------------------------85

subroutine app_euler1d_tests

use activate_routines
use set_up_routines
use boundary_condition_types
! This is something defined only in this application.  See above.
use problem_type_module
use logical_units
use boundary_condition_types
use boundary_condition_routines
use q_array, only: q

! Use defined.
use acoustic_pulse_parameters, only: p0, rho0, z0, eps, sigma
use grid, only: zgrid
! Debug:
use artificial_pressure_module
implicit none

! Locals:
external rhs
real(8) :: t, cfl, dt, dt_output, t_target
integer :: n_output_steps, i_output_step, istep
logical :: hit_target, target_met

! For mesh and domain:
real(8) :: rmin, rmax, zmin, zmax, phi_min, phi_max
integer :: nr, nz, nphi
logical :: suppress_z_derivatives_when_nz_not_1, periodic_z

! For thermal parameters:
real(8) :: gamma
logical :: isothermal
real(8), allocatable, dimension(:, :) :: ci_squared_initial
logical :: apply_newtonian_cooling
real(8) :: tau_newtonian_cooling

real(8) :: gm1

! For subroutine rk4
logical :: use_supplied_dt = .false.

! Pade filter stuff:
logical :: apply_pade_filter
character(3) :: eps_or_tau
real(8) :: eps_filter, tau_filter

! For restart:
logical :: restart
integer :: istep_of_restart_file
real(8) :: t_of_restart_file

real(8), allocatable, dimension(:) :: f, fp

! For set_up_boundary_conditions
integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced
real(8), allocatable :: d_ci_dr_inner(:), d_ci_dr_outer(:)
real(8) :: c_sound_rmin, c_sound_rmax

integer :: iz
real(8) :: twopi

logical, parameter :: stretched_r = .false., stretched_z = .false.
real(8), parameter :: r0 = 0.d0, z0_unif = 0.d0
integer, parameter :: nr_u = 0, nz_u = 0
logical :: unphysical
integer :: iphi_bad

! For initializing the exact solution to the shock-tube problem:
logical :: initialize
! Used for irrelevant argumemts:
real(8) :: dummy = 0.d0

! Debug
logical, parameter :: get_lambda_max_ap = .true.
real(8) :: lambda_max_ap

logical :: filter_relative_to_basic_state = .false.

namelist /euler1d_tests_input/ problem_type, nz, n_output_steps, dt_output, cfl, &
     isothermal, periodic_z, apply_pade_filter, eps_filter, apply_artificial_pressure, &
     C_ap 

open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', status = 'old')
read (lun_general_purpose, nml = euler1d_tests_input)
close (lun_general_purpose)

twopi = 4.0d0 * atan(1.0d0)
allocate(ci_squared_initial(1, nz))

! For set-up domain mesh and partition:
rmin    = 1.0d0
rmax    = 1.0d0
phi_min = 0.0d0
phi_max = 0.0d0
nr      = 1
nphi    = 1
suppress_z_derivatives_when_nz_not_1 = .false.

if (problem_type .eq. acoustic) then
   zmin    = 0.0d0
   zmax    = 4.0d0
else if (problem_type .eq. shock_tube) then
   zmin = 0.d0
   zmax = 2.d0
else if (problem_type .eq. density_waves) then
   zmin = -5.d0
   zmax =  5.d0
end if

! For set_up_thermal_parameters:
gamma      = 1.4d0
gm1        = gamma - 1.d0

! For set-up_boundary_conditions:
rmin_BC = null
rmax_BC = null

if (problem_type .eq. acoustic) then
   apply_artificial_pressure = .false.
   apply_pade_filter         = .false.
   if (periodic_z) then
      zmin_BC = periodic
      zmax_BC = periodic
   else
      zmin_BC = non_reflective
      zmax_BC = non_reflective
      ibalanced = 1   
   end if
   eps = 0.005d0   ! amplitude of pulse
   z0  = 2.0d0     ! origin of pulse
   sigma = 0.04d0 * (zmax - zmin)  ! width of pulse
   ! These make the sound speed equal to unity:
   if (isothermal) then
      p0 = 1.0d0
   else
      p0 = 1.d0 / gamma
   end if
   rho0 = 1.d0
else if ((problem_type .eq. shock_tube) .or. (problem_type .eq. density_waves)) then
   isothermal = .false.
   periodic_z = .false.
   zmin_BC = z_dirichlet
   zmax_BC = z_dirichlet
   print *, ' isothermal  = ', isothermal
   print *, ' periodic_z  = ', periodic_z
   print *, ' zmin_BC = ', zmin_BC
   print *, ' zmax_BC = ', zmax_BC
end if

! These calls should be in every application subroutine:
call set_up_domain_mesh_and_partition(rmin, rmax, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
     suppress_z_derivatives_when_nz_not_1, periodic_z, &
     stretched_r, stretched_z, r0, nr_u, z0_unif, nz_u)
! finite_sphere = .false., R_sphere = 0.0d0
if (isothermal) then

   ! Not used here:
   allocate(d_ci_dr_inner(nz), d_ci_dr_outer(nz))
   c_sound_rmin = 1.0d0
   c_sound_rmax = 1.d0
   do iz = 1, nz
      ci_squared_initial(1, iz) = 1.0d0
   end do
end if

apply_newtonian_cooling = .false.
tau_newtonian_cooling   = 0.d0
call set_up_thermal_parameters(gamma, isothermal, apply_newtonian_cooling, &
     tau_newtonian_cooling, ci_squared_initial)

print *, ' zmin_BC = ', zmin_BC, ' zmax_BC = ', zmax_BC

allocate(d_ci_dr_inner(1), d_ci_dr_outer(1))
call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, &
     d_ci_dr_inner, d_ci_dr_outer, c_sound_rmin, c_sound_rmax, isothermal)

if (apply_pade_filter) then
   eps_or_tau = 'eps'
   call activate_pade_filter(apply_pade_filter, eps_or_tau, eps_filter, &
        tau_filter, filter_relative_to_basic_state)
   print *, ' apply_pade_filter = ', apply_pade_filter
   print *, ' eps_filter = ', eps
end if

if (apply_artificial_pressure) then
   call activate_artificial_pressure(C_ap)
   print *, ' apply_artificial_pressure = ', apply_artificial_pressure
   print *, ' C_ap = ', C_ap
end if

if (problem_type .eq. acoustic) then
   call acoustic_pulse_initial_condition(isothermal, gamma, q)
else if (problem_type .eq. shock_tube) then
   call shock_tube_initial_condition(q)
   ! This call must come after the previous one because this routine needs
   ! parameters from the shock tube initial condition:
   initialize = .true.
   call exact_shock_tube(initialize, dummy, dummy, dummy, dummy, dummy)
else if (problem_type .eq. density_waves) then
   call density_wave_initial_condition(q)
end if

t     = 0.0d0
istep = 0
if (problem_type .eq. acoustic) then
   call acoustic_pulse_output(gamma, gm1, q, t)
else if (problem_type .eq. shock_tube) then
   call shock_tube_output(q, t, istep)
else if (problem_type .eq. density_waves) then
   call density_wave_output(q, t, istep)   
end if   

hit_target = .true.
t_target   = t + dt_output
target_met = .false.
do i_output_step = 1, n_output_steps
   do while (.not. target_met)
      call rk4(rhs, q, cfl, t, dt, use_supplied_dt, unphysical, iphi_bad, &
           hit_target, t_target, target_met)
      if (unphysical) call terminate_with_no_save(1)   
      istep = istep + 1
      print *, ' finished rk4 istep = ', istep, ' dt = ', dt, ' t = ', t
      print *, ' target_met = ', target_met, ' t_target = ', t_target
   end do
   
   if (problem_type .eq. acoustic) then
      call acoustic_pulse_output(gamma, gm1, q, t)
   else if (problem_type .eq. shock_tube) then
      print *, ' calling shock_tube_output'
      call shock_tube_output(q, t, istep)
   else if (problem_type .eq. density_waves) then
      call density_wave_output(q, t, istep)   
   end if   

   t_target = t_target + dt_output
   target_met = .false.   
end do

call terminate_with_no_save(0)

end subroutine app_euler1d_tests

!----------------------------------------------------------------------------------85

subroutine acoustic_pulse_initial_condition(isothermal, gamma, q)

use dof_indices, only: irho, rmom, amom, zmom, ener
use grid
implicit none

logical :: isothermal
real(8), intent(in)  :: gamma
real(8), intent(out) :: q(nr, nphi, nz, ndof)

! Local:
integer :: ir, iphi, iz

real(8) :: rho, rho_uz, eint

do iz = 1, nz
   do iphi = 1, nphi
      do ir = 1, nr
         call acoustic_pulse(isothermal, gamma, zgrid(iz), rho, rho_uz, eint)
         q(ir, iphi, iz, irho) = rho
         q(ir, iphi, iz, rmom) = 0.d0
         q(ir, iphi, iz, amom) = 0.d0
         q(ir, iphi, iz, zmom) = rho_uz
         q(ir, iphi, iz, ener) = eint
      end do
   end do
end do

end subroutine acoustic_pulse_initial_condition

!----------------------------------------------------------------------------------85

subroutine acoustic_pulse(locally_isothermal, gamma, z, rho, rho_uz, eint)

! Gaussian acoustic pulse in "z" in a medium at rest.
! eps      : relative amplitude in density.
! sigma    : width parameter.
! z0       : center of the pulse.
! p0, rho0 : basic state pressure and density.
! gamma    : ratio of specific heats.
! z        : value of z at which rho, rho_u and energy are output.

use acoustic_pulse_parameters, only: p0, rho0, z0, eps, sigma
implicit none
logical :: locally_isothermal
real(8), intent(in)  :: gamma, z
real(8), intent(out) :: rho, rho_uz, eint

! Locals:
real(8) :: a0, z_shift, s_tilde, p, uz

z_shift = z - z0
! s_tilde is the notation of Lieppman and Roshko.
s_tilde = eps*exp(-z_shift**2/sigma**2)

if (locally_isothermal) then
   a0      = sqrt(p0/rho0) ! speed of sound in the medium
   rho     = rho0*(1.d0 + s_tilde)
   p       = p0*(1.d0 + s_tilde)
   uz      = a0*s_tilde   
else
   a0      = sqrt(gamma*p0/rho0) ! speed of sound in the medium
   rho     = rho0*(1.d0 + s_tilde)
   p       = p0*(1.d0 + gamma*s_tilde)
   uz      = a0*s_tilde
end if

rho_uz  = rho * uz
eint  = p/(gamma - 1) ! inconsequential for the isothermal case

return
end subroutine acoustic_pulse

!----------------------------------------------------------------------------------85

subroutine acoustic_pulse_output(gamma, gm1, q, t)

use grid
use dof_indices
use acoustic_pulse_parameters, only: p0, rho0
implicit none

! Arguments:
real(8), intent(in) :: gamma, gm1
real(8), intent(in) :: q(nr, nphi, nz, ndof)
real(8) :: t

! Local:
character(40) :: filename
integer       :: iz
real(8)       :: rho_prime, uz, p, p_prime
real(8)       :: pressure_at_point ! function called

write (filename, "('rhop_u_p_t=', f5.3, '.dat')") t
open (unit = 1, file = filename, form = 'formatted', status = 'unknown')
do iz = 1, nz
   rho_prime = q(1,1, iz, irho) - rho0
   uz        = q(1,1, iz, zmom) / q(1,1, iz, irho)
   p = pressure_at_point(q, 1, 1, iz)
   p_prime = p - p0
   write (1, "(4(1x, e13.5e3))") zgrid(iz), rho_prime, uz, p_prime

   if (q(1, 1, iz, amom) .ne. 0.d0) then
      print *, ' angular momentum .ne. 0 for iz = ', iz
      print *, ' q(1, 1, iz, amom) = ', q(1, 1, iz, amom) 
      stop
   end if

   if (q(1, 1, iz, rmom) .ne. 0.d0) then
      print *, ' radial momentum .ne. 0 for iz = ', iz
      stop
   end if
end do
close(1)
end subroutine acoustic_pulse_output

!----------------------------------------------------------------------------------85

subroutine shock_tube_initial_condition(q)

use thermal_parameters
use shock_tube_parameters
use grid
use dof_indices
! We set-up Dirichlet BC here:
use boundary_condition_data
implicit none
real(8), intent(out) :: q(nr, nphi, nz, ndof)

! Local:
real(8) :: rho, u, p
integer :: iz, idof

! rhoL = 1.d0
! pL   = 1.d0
! uL   = 0.d0

! rhoR = 0.125d0
! pR_Sod   = 0.1d0
! uR   = 0.d0

! Hofner:
rhoL = 8.d0
pL   = 10.d0/gamma
uL   = 0.d0

rhoR   = 1.d0
pR_Sod = 1.d0/gamma
uR     = 0.d0

! This is a little to the left because the difference in
! speed in the expansion wave.
z0_shock_tube_actual = 0.8d0

do iz = 1, nz
   if (zgrid(iz) .lt. z0_shock_tube_actual) then
      rho = rhoL
      p   = pL
      u   = uL
   else
      rho = rhoR
      p   = pR_Sod
      u   = uR
   end if

   q(1, 1, iz, irho) = rho
   q(1, 1, iz, zmom) = rho*u   
   q(1, 1, iz, ener) = p/gm1   
   q(1, 1, iz, rmom) = 0.d0
   q(1, 1, iz, amom) = 0.d0      
end do

do idof = 1, ndof
   qLeft_z_BC (idof)  = q(1, 1, 1,  idof)
   qRight_z_BC(idof)  = q(1, 1, nz, idof)
end do
   
end

!----------------------------------------------------------------------------------85

subroutine density_wave_initial_condition(q)

use thermal_parameters
use shock_tube_parameters
use grid
use dof_indices
! We set-up Dirichlet BC here:
use boundary_condition_data
implicit none
real(8), intent(out) :: q(nr, nphi, nz, ndof)

! Local:
real(8) :: rho, u, p
integer :: iz, idof

do iz = 1, nz
   if (zgrid(iz).lt. -4.d0) then
      rho = 3.857143
      p   = 10.33333
      u   = 2.629369
   else
      rho = 1.d0 + 0.2d0*sin(5.d0*zgrid(iz))
      p   = 1.d0
      u   = 0.d0
   end if

   q(1, 1, iz, irho) = rho
   q(1, 1, iz, zmom) = rho*u   
   q(1, 1, iz, ener) = p/gm1   
   q(1, 1, iz, rmom) = 0.d0
   q(1, 1, iz, amom) = 0.d0      
end do

do idof = 1, ndof
   qLeft_z_BC (idof)  = q(1, 1, 1,  idof)
   qRight_z_BC(idof)  = q(1, 1, nz, idof)
end do
   
end

!----------------------------------------------------------------------------------85

subroutine density_wave_output(q, t, istep)

use thermal_parameters
use dof_indices
use grid
implicit none

! Arguments:
real(8), intent(in) :: q(nr, nphi, nz, ndof)
real(8), intent(in) :: t
integer, intent(in) :: istep

! Locals:
integer :: iz
real(8) :: rho, u, p
character(40) :: filename

write (filename, 1) int(t), t - int(t), istep
1 format ('rho_t_', i4.4, f0.4, '_step_', i6.6, '.dat')
open(unit = 1, file = filename, form = 'formatted', status = 'unknown')

write (filename, 2) int(t), t - int(t), istep
2 format ('u_t_', i4.4, f0.4, '_step_', i6.6, '.dat')
open(unit = 2, file = filename, form = 'formatted', status = 'unknown')

write (filename, 3) int(t), t - int(t), istep
3 format ('p_t_', i4.4, f0.4, '_step_', i6.6, '.dat')
open(unit = 3, file = filename, form = 'formatted', status = 'unknown')

do iz = 1, nz
   rho = q(1, 1, iz, irho)
   u   = q(1, 1, iz, zmom) / rho
   p   = q(1, 1, iz, ener) * gm1
   write(1, "(2(1x, e12.5))") zgrid(iz), rho
   write(2, "(2(1x, e12.5))") zgrid(iz), u
   write(3, "(2(1x, e12.5))") zgrid(iz), p
end do

close(1); close(2); close(3)
return
end subroutine density_wave_output

!----------------------------------------------------------------------------------85

subroutine shock_tube_output(q, t, istep)

use thermal_parameters
use dof_indices
use grid
implicit none

! Arguments:
real(8), intent(in) :: q(nr, nphi, nz, ndof)
real(8), intent(in) :: t
integer, intent(in) :: istep

! Locals:
integer :: iz
real(8) :: rho, u, p, rho_ex, u_ex, p_ex
logical, parameter :: initialize = .false.
character(40) :: filename

print *, ' in shock tube output'
print *, ' t = ', ' istep = ', istep
read(5, *)

write (filename, 1) int(t), t - int(t), istep
1 format ('rho_t_', i4.4, f0.4, '_step_', i6.6, '.dat')
open(unit = 1, file = filename, form = 'formatted', status = 'unknown')

write (filename, 2) int(t), t - int(t), istep
2 format ('u_t_', i4.4, f0.4, '_step_', i6.6, '.dat')
open(unit = 2, file = filename, form = 'formatted', status = 'unknown')

write (filename, 3) int(t), t - int(t), istep
3 format ('p_t_', i4.4, f0.4, '_step_', i6.6, '.dat')
open(unit = 3, file = filename, form = 'formatted', status = 'unknown')

do iz = 1, nz
   call exact_shock_tube(initialize, zgrid(iz), t, rho_ex, u_ex, p_ex)
   rho = q(1, 1, iz, irho)
   u   = q(1, 1, iz, zmom) / rho
   p   = q(1, 1, iz, ener) * gm1
   write(1, "(3(1x, e12.5))") zgrid(iz), rho, rho_ex
   write(2, "(3(1x, e12.5))") zgrid(iz), u,   u_ex
   write(3, "(3(1x, e12.5))") zgrid(iz), p,   p_ex   
end do

close(1); close(2); close(3)
return
end subroutine shock_tube_output

!-----------------------------------------------------------------------------80

subroutine exact_shock_tube(initialize, z, t, rho, u, p)

use thermal_parameters
use shock_tube_parameters
implicit none
logical, intent(in)  :: initialize
real(8), intent(in)  :: z, t
real(8), intent(out) :: u, p, rho

! Local:
real(8), parameter :: tol = 1.d-6
real(8) :: Us
real(8) :: p1_over_pR, rhoR_over_rho1
real(8) :: u1, p1, rho1, u2, p2, rho2
real(8) :: z1, z2, z3, z4
real(8) :: a, exponent, a2

! Debugging:
integer :: i
real(8) :: Ms1, Ms2, del

! Passed function name:
real(8) :: compatibility
external compatibility

if (initialize) then
   aL = SQRT(gamma*pL/rhoL)
   aR = SQRT(gamma*pR_Sod/rhoR)
   print *, ' aL = ', aL, ' aR = ', aR, ' rhoL = ', rhoL
   read(5, *)

   ! Plot in case we cannot find a root:
   Ms1 = 0.5d0
   Ms2 = 10.d0
   del = (Ms2 - Ms1) / 99.d0
   open(unit = 1, file = 'compat.dat', form = 'formatted', status = 'unknown')
   do i = 1, 100
      Ms = Ms1 + (i-1)*del
      write(1, "(2(1x, e12.5))") Ms, compatibility(Ms)
   end do
   close(1)
   print *, ' calling zeroin'
   call zeroin(0.5d0, 10.d0, compatibility, tol, Ms)
   print *, ' Shock Mach number = ', Ms
   return
end if

! Shock speed:
Us = Ms * aR

p1_over_pR     = 2.d0*gamma/(gamma + 1.d0)*Ms**2 - (gamma - 1.d0)/(gamma + 1.d0)
rhoR_over_rho1 = 2.d0/(gamma + 1.d0)/Ms**2       + (gamma - 1.d0)/(gamma + 1.d0)

p1   = p1_over_pR * pR_Sod
rho1 = rhoR / rhoR_over_rho1
u1   = 2.d0/(gamma + 1.d0) * (Ms - 1.d0/Ms)

! Across the contact surface:
u2   = u1
p2   = p1
rho2 = rhoL*(p2/pL)**(1.d0/gamma)
a2   = SQRT(gamma*p2/rho2)

z1 = z0_shock_tube_actual - aL*t
z2 = z0_shock_tube_actual + (u2 - a2)*t
z3 = z0_shock_tube_actual + u2*t
z4 = z0_shock_tube_actual + Us*t

if (z .lt. z1) then
   rho = rhoL
   u   = uL
   p   = pL
else if ((z .ge. z1) .and. (z .lt. z2)) then
   ! Expansion fan:
   u        = 2.d0/(gamma + 1.d0) * (aL + (z-z0_shock_tube_actual)/t)
   a        = aL - (gamma - 1.d0)*u/2
   exponent = 2.*gamma / (gamma - 1.d0)
   p = pL*(a/aL)**exponent
   rho = gamma*p/a**2
else if ((z .ge. z2) .and. (z .lt. z3)) then
   ! Region 2:
   u   = u2
   p   = p2
   rho = rho2
else if ((z .ge. z3) .and. (z .lt. z4)) then
   ! Region 1:
   u   = u1
   p   = p1
   rho = rho1
else if (z .ge. z4) then
   u   = uR
   p   = pR_Sod
   rho = rhoR
else
   print *, ' exact_shock_tube: Invalid location z = ', z
end if

end subroutine exact_shock_tube

!-----------------------------------------------------------------------------80

real(8) function compatibility(Ms_var)

! Compatibility condition whose root determines the shock Mach number.
! We use the name Ms_var to avoid conflict with Ms in
! module initial_condition_info

! From Susanne Hofner class notes.

use thermal_parameters
use shock_tube_parameters, only: aL, pR_Sod, pL
implicit none
real(8) :: Ms_var

! Locals:
real(8) :: lhs, rhs, fac, paren, bracket, curly, exponent

lhs      = Ms_var - 1.d0/Ms_var

fac      = aL*(gamma + 1.d0)/(gamma - 1.d0)
paren    = 2.d0*gamma/(gamma + 1.)*Ms_var**2 - (gamma - 1.d0)/(gamma + 1.d0)
bracket  = pR_Sod/pL*paren
exponent = (gamma - 1.d0) / (2.d0*gamma) 
curly   = 1.d0 - bracket**exponent
rhs     = fac*curly

compatibility = rhs - lhs

end function compatibility

!----------------------------------------------------------------------------------85
