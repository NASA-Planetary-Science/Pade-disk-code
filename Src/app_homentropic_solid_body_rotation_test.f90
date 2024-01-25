!----------------------------------------------------------------------------------85

module parameters_for_homentropic_solid_body_rotation_test
   implicit none
   real(8), parameter :: gamma = 1.4d0, gm1 = gamma - 1.d0
   ! We can set three quantities to unity.
   real(8), parameter :: Omega_SBR = 1.d0, r_outer_SBR = 1.d0, rho_outer_SBR = 1.d0, &
        Mach_outer_SBR = 0.30d0, c_outer = Omega_SBR * r_outer_SBR / Mach_outer_SBR, &
        p_outer = rho_outer_SBR * c_outer**2 / gamma, &
        entropy_constant = p_outer / rho_outer_SBR**gamma
end module parameters_for_homentropic_solid_body_rotation_test

!----------------------------------------------------------------------------------85

subroutine app_homentropic_solid_body_rotation_test

use set_up_routines
use grid
use q_array
use parameters_for_homentropic_solid_body_rotation_test
use boundary_condition_types
use boundary_condition_routines

!Debug
use dof_indices
implicit none

! Local:
external rhs
real(8) :: t, cfl, dt
integer :: istep, nsteps, output_interval

! For thermal set-up:
logical :: isothermal = .false., apply_newtonian_cooling
logical :: restart = .false., &
     suppress_z_derivatives_when_nz_not_1_arg = .true., periodic_z_arg = .false., &
     use_supplied_dt = .false.

! For set_up_boundary_conditions
integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced
real(8), allocatable :: d_ci_dr_inner(:), d_ci_dr_outer(:)
real(8) :: c_sound_rmin, c_sound_rmax

logical, parameter :: stretched_r = .false., stretched_z = .false.
real(8), parameter :: r0 = 0.d0, z0 = 0.d0
integer, parameter :: nr_u = 0, nz_u = 0

real(8) :: tau_newtonian_cooling
real(8), allocatable, dimension(:,:) :: ci_squared_initial
logical :: unphysical, hit_target, target_met
real(8) :: t_target
integer :: iphi_bad
! For irrelevant arguments:
real(8) :: dummy = 0.d0

rmin = 0.5d0
rmax = 1.0d0
zmin = 0.0d0
zmax = 0.0d0
phi_min = 0.d0
phi_max = 0.d0

nsteps = 200
cfl    = 0.5d0
output_interval = 50

nr   = 64
nz   = 32
nphi = 1
zmin = 0.0d0
zmax = 0.0d0


rmin_BC = null
rmax_BC = null
zmin_BC = null
zmax_BC = null
ibalanced = null

! These three set_up calls should be in every application subroutine:
call set_up_domain_mesh_and_partition(rmin, rmax, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
     suppress_z_derivatives_when_nz_not_1_arg, periodic_z_arg, &
     stretched_r, stretched_z, r0, nr_u, z0, nz_u)

! Note: I am leaving out the last two arguments "ci_squared_initial" and
! "tau_newtonian_cooling" which are optional and not needed here.
tau_newtonian_cooling   = 0.d0
apply_newtonian_cooling = .false.
call set_up_thermal_parameters(gamma, isothermal, apply_newtonian_cooling, tau_newtonian_cooling, &
     ci_squared_initial)

call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, &
     d_ci_dr_inner, d_ci_dr_outer, c_sound_rmin, c_sound_rmax, isothermal)

!if (restart) then
!   call read_restart (istep_of_restart_file, t_of_restart_file)
!end if

print *, ' calling initial_condition_for_HSBR'
call initial_condition_for_HSBR(gamma, gm1, q)
print *, ' returned from initial_condition_for_HSBR'

istep = 0
t = 0.0d0

print *, ' Calling output_for_HSBR'
call output_for_HSBR(istep, t, q)
print *, ' Returned from output_for_HSBR'

hit_target = .false.
t_target   = dummy
do istep = 1, nsteps
   ! call euler (q, cfl, t)
   call rk4(rhs, q, cfl, t, dt, use_supplied_dt, unphysical, iphi_bad, &
                       hit_target, t_target, target_met)
   print *, ' finished istep = ', istep, ' t = ', t

   !if (mod(istep, filter_interval) .eq. 0) then
   !   call filter_r (q, filter_param)
   !end if
   
   if (mod(istep, output_interval) .eq. 0) then
      call output_for_HSBR(istep, t, q)
   end if
end do

end subroutine app_homentropic_solid_body_rotation_test

!----------------------------------------------------------------------------------85

subroutine initial_condition_for_HSBR(gamma, gm1, q)

use grid
use dof_indices
use partition_data
implicit none
real(8) :: gamma, gm1
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
integer :: ir, iz, iphi, idof
real(8) :: fac, rho, uz, ur, uphi, pressure, energy

do iphi = sphi, ephi
   do iz = 1, nz
      do ir = sr, er
         call HSBR_flow (zgrid(iz), rgrid(ir), uz, ur, uphi, rho, pressure, energy)
         q(ir, iphi, iz, irho) = rho
         q(ir, iphi, iz, rmom) = 0.0d0
         q(ir, iphi, iz, zmom) = 0.0d0
         q(ir, iphi, iz, amom) = rho * uphi * rgrid(ir)
         q(ir, iphi, iz, ener) = energy
         ! print *, ' IC routine: ir = ', ir, ' iz = ', iz, ' iphi = ', iphi
         ! print *, ' rho    = ', q(ir, iphi, iz, irho)
         ! print *, ' amom   = ', q(ir, iphi, iz, amom)
         ! print *, ' energy = ', q(ir, iphi, iz, ener)
      end do
      ! print *, ' enter anything to continue'
      ! read (5, *)
   end do
end do

end subroutine initial_condition_for_HSBR

!----------------------------------------------------------------------------------85

subroutine HSBR_flow(z, r, uz, ur, uphi, rho, pressure, energy)
! Homentropic solid body rotation flow.
use parameters_for_homentropic_solid_body_rotation_test
implicit none
real(8) :: z, r, uz, ur, uphi, rho, pressure, energy

! Local:
real(8) :: fac

uz       = 0.0d0
ur       = 0.0d0
fac      = 1.d0 + gm1 / c_outer**2 * Omega_SBR**2 / 2.d0 * (r**2 - r_outer_SBR**2)
rho      = fac**(1.d0/gm1) * rho_outer_SBR
uphi     = Omega_SBR * r
pressure = entropy_constant * rho**gamma
energy   = pressure/gm1

end subroutine HSBR_flow

!----------------------------------------------------------------------------------85

subroutine output_for_HSBR(istep, t, q)

use grid
use logical_units
use dof_indices
use partition_data

implicit none
integer :: istep
real(8) :: t
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
integer, parameter :: nfiles = 3
character*40, dimension(5) :: filename
integer :: ir, iz, iphi, i
real(8) :: rho, ur, uz, uphi

real(8) rho_exact, uphi_exact, p_exact, e_exact, uz_exact, ur_exact

write (filename(1), "('ur',          i4.4, '.dat')") istep
write (filename(2), "('uphi_',       i4.4, '.dat')") istep
write (filename(3), "('rho_',        i4.4, '.dat')") istep

do i = 1, nfiles
   open (unit = lun_profile(i), file = filename(i), form = 'formatted', status = 'unknown', access = 'direct', recl = 13*3+1)
end do

! Select the mdiplane:
if ((nz .eq. 1) .or. suppress_z_derivatives_when_nz_not_1) then
   iz = 1
else
   iz = (nz - 1)/2
end if
iphi = 1
do ir = sr, er
   call HSBR_flow(zgrid(iz), rgrid(ir), uz_exact, ur_exact, uphi_exact, rho_exact, p_exact, e_exact)
   rho  = q(ir, iphi, iz, irho)
   uphi = q(ir, iphi, iz, amom) / rgrid(ir) / rho
   ur   = q(ir, iphi, iz, rmom) / rho
   write (lun_profile(1), 4, rec = ir) rgrid(ir), ur,   0.0d0,      char(10)
   write (lun_profile(2), 4, rec = ir) rgrid(ir), uphi, uphi_exact, char(10)
   write (lun_profile(3), 4, rec = ir) rgrid(ir), rho,  rho_exact,  char(10)
   4  format (3(1x, e12.5), a1)
end do

do i = 1, nfiles
   close (lun_profile(i))
end do

end subroutine output_for_HSBR

!----------------------------------------------------------------------------------85

