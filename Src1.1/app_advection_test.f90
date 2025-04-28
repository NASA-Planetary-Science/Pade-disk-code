!----------------------------------------------------------------------------------85

module advection_test_params
real(8), parameter :: rho0 = 1.d0, eps_rho = 0.5d0, uz_advec = 1.d0, Omega_advec = 1.d0
integer, parameter :: kphi = 2, kz = 2
end module advection_test_params

!----------------------------------------------------------------------------------85

subroutine app_advection_test

use grid, only: rgrid, zgrid, phi_grid
use dof_indices
use boundary_condition_types 
use q_array, only: q
use logical_units
use control_parameters, only: restart_file, save_file
use total_allocated_words, only: n_words_allocated
use math_constants, only: pi
use partition_data
use set_up_routines
use activate_routines
use advection_test_params
#ifdef mpi_code
   use mpi
#endif
implicit none

logical :: restart
integer :: istep_of_restart_file
integer, parameter :: file_version = 2
real(8) :: t_of_restart_file

integer :: nr, nz, nphi
real(8) :: zmin, zmax, rmin, rmax, phi_min, phi_max, dz, dphi
integer :: istep, nsteps, istep0
real(8) :: t, cfl, dt_used, tmax_over_pi, tmax, lambda_z, lambda_phi, lambda_max

! Pade filter stuff:
logical :: apply_pade_filter, filter_relative_to_basic_state = .false.
character(3) :: eps_or_tau
real(8) :: eps_filter, tau_filter
integer :: filtering_interval = 1

! For thermal set-up:
real(8), parameter :: gamma = 1.4d0
logical :: isothermal = .true.
logical :: apply_newtonian_cooling = .false.
real(8) :: Pr_molecular = 0.7d0 ! Not needed for isothermal.
real(8) :: tau_newtonian_cooling = 0.d0 ! dummy

integer :: iz, ir, iphi, ier

logical :: perturb, suppress_z_derivatives_when_nz_not_1_arg = .false., periodic_z_arg = .true., &
    use_supplied_dt = .false., name_using_step = .false., output_profiles
real(8) :: dt_supplied = 0.d0

! For passing to subroutine set_up_boundary_conditions.
integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced
real(8), allocatable, dimension(:) :: d_ci_dr_inner, d_ci_dr_outer

real(8), parameter :: c_sound = 1.d0
real(8) :: pressure

! Stretched mesh stuff.
logical, parameter :: stretched_r = .false., stretched_z = .false.
real(8), parameter :: r0 = 0.d0, z0 = 0.d0
integer, parameter :: nr_u = 0, nz_u = 0

character(25) :: filename_given
logical :: unphysical
integer :: iphi_bad

real(8), allocatable, dimension(:,:) :: ci_squared_initial

logical :: hit_target = .false., target_met
real(8) :: t_target, err_L2, rho_func

namelist /advec_input/ restart, nr, nz, nphi, tmax_over_pi, cfl, apply_pade_filter, eps_filter

print *, ' first executable in app_advection_test', ' my_node = ', my_node

pi = 4.d0*atan(1.d0)

if (my_node .eq. 0) then
   print *, ' node 0: about to open and read namelist file for advection test'      
   open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', status = 'old')
   read (lun_general_purpose, nml = advec_input)
   close (lun_general_purpose)

   print *, ' node 0: read namelist input file'
   print *, ' node 0: nr = ', nr, ' nz = ', nz, ' nphi = ', nphi
end if

#ifdef mpi_code
   call mpi_bcast(restart,                   1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(nr,                        1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nz,                        1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nphi,                      1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(tmax_over_pi,              1, mpi_double,  0, mpi_comm_world, ier)   
   call mpi_bcast(cfl,                       1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(apply_pade_filter,         1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(eps_filter,                1, mpi_double,  0, mpi_comm_world, ier)   
#endif

if (my_node .eq. 0) then
   print *, ' my_node = ', my_node
   print *, ' file_version = ', file_version
   print *, ' isothermal   = ', isothermal
end if

pi = 4.0*ATAN(1.0d0)
! For thermal set up and initial density profile:
eps_or_tau = 'eps'

allocate(d_ci_dr_inner(nz), d_ci_dr_outer(nz))

phi_min = 0.0d0
phi_max = 2.0*pi
   
zmin = 0.0d0
zmax = 2.d0*pi

rmin = 1.d0
rmax = 2.d0

rmin_BC = null
rmax_BC = null
zmin_BC = periodic
zmax_BC = periodic
ibalanced = 0
periodic_z_arg = .true. ! For mesh set-up.

use_supplied_dt = .true.

dz = 2.d0 * pi / nz
lambda_z = uz_advec / dz

dphi = 2.0 * pi / nphi
lambda_phi = Omega_advec / dphi

lambda_max = max(lambda_z, lambda_phi)
dt = cfl / (pi * lambda_max)

tmax = tmax_over_pi * pi
nsteps = tmax / dt + 2

hit_target = .true.
t_target   = tmax

! These calls should be in every application subroutine:
call set_up_domain_mesh_and_partition(rmin, rmax, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
     suppress_z_derivatives_when_nz_not_1_arg, periodic_z_arg, &
     stretched_r, stretched_z, r0, nr_u, z0, nz_u)

! Needed for the isothermal option.  ci_squared_initial sits in
! module thermal_parameters.
if (isothermal) then
   allocate(ci_squared_initial(nr,nz))
   do iz = 1, nz
      do ir = 1, nr ! Note each processor has all of r.
         ci_squared_initial(ir, iz) = c_sound**2
      end do
   end do
end if

call set_up_thermal_parameters(gamma, isothermal, apply_newtonian_cooling, tau_newtonian_cooling, &
     ci_squared_initial)


print *, ' all nodes.  my_node = ', my_node, ' advection_test: returned from set_up_thermal_parameters'
call mpi_barrier(mpi_comm_world, ier)

! Note: c_sound_i and c_sound_o are used only for viscous_wall BC at rmin and rmax
! and the non-isothermal case.
call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, &
     d_ci_dr_inner, d_ci_dr_outer, c_sound, c_sound, isothermal)

print *, ' all nodes, my_node = ', my_node, ' app_advection_test: returned from set_up_boundary_conditions'
call mpi_barrier(mpi_comm_world, ier)

call activate_pade_filter(apply_pade_filter, eps_or_tau, eps_filter, &
     tau_filter, filter_relative_to_basic_state, filtering_interval)

! Initial condition:
! ~~~~~~~~~~~~~~~~~~
if (.not. restart) then
   t = 0.d0
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, iz, irho) = rho_func(t, phi_grid(iphi), zgrid(iz))
            ! print *, ' In main rho = ', q(ir, iphi, iz, irho)
            q(ir, iphi, iz, zmom) = q(ir, iphi, iz, irho) * uz_advec
            q(ir, iphi, iz, rmom) = 0.0d0
            q(ir, iphi, iz, amom) = q(ir, iphi, iz, irho) * Omega_advec * rgrid(ir)**2
            ! In case we are running non-isothermal:
            pressure = rho0*c_sound**2.d0
            q(ir, iphi, iz, ener) = pressure / (gamma - 1.d0)
         end do
      end do
   end do
end if

print *, ' all nodes, my_node = ', my_node, ' finished initial condition'
call mpi_barrier(mpi_comm_world, ier)

if (.not. restart) then
   istep0 = 0
   t      = 0.0d0
   print *, ' all nodes, my_node = ', my_node, ' after not restart'   
else
   call read_restart (file_version, istep_of_restart_file, t_of_restart_file)   
   istep0 = istep_of_restart_file
   t      = t_of_restart_file
end if

print *, ' all nodes, my_node = ', my_node, ' before ir = 1'
call mpi_barrier(mpi_comm_world, ier)

ir = 1
call tecplot_cylinder(q, ir, t, istep)

if (my_node .eq. 0.d0) then
   open(unit = lun_general_purpose, file = 'err_L2vs_t.dat', form = 'formatted', status = 'unknown')
end if

print *, ' my_node = ', my_node, ' calling error'
call error(t, q, err_L2)               
if (my_node .eq. 0) then
   write(lun_general_purpose, "(2(1x, e12.5))") t, err_L2
end if

do istep = istep0 + 1, istep0 + nsteps   
   call rk4(istep, q, cfl, t, dt_used, use_supplied_dt, dt_supplied, unphysical, iphi_bad, &
        hit_target, t_target, target_met)
   call error(t, q, err_L2)                 
   if (my_node .eq. 0) then
      print *, ' node 0: finished istep = ', istep, ' t = ', t, ' L2 error = ', err_L2
      write(lun_general_purpose, "(2(1x, e12.5))") t, err_L2
   end if

   if (target_met) then
      print *, ' target met'
      exit
   end if
end do

call tecplot_cylinder(q, 1, t, istep)

call error(t, q, err_L2)               
if (my_node .eq. 0) then
   write(lun_general_purpose, "(2(1x, e12.5))") t, err_L2
   print *, ' node 0: finished istep = ', istep, ' t = ', t, ' L2 error = ', err_L2
   close(lun_general_purpose)
end if

! Need the minus 1 since a fortran do loop increments the counter at the end of the
! loop.
call terminate_with_save(0, istep-1, t)

end subroutine app_advection_test

!----------------------------------------------------------------------------------85

real(8) function rho_func(t, phi, z)

use advection_test_params
implicit none
real(8), intent(in) :: t, phi, z

! Local:
real(8) :: arg_phi, arg_z

arg_phi = phi - Omega_advec*t
arg_z   = z   - uz_advec*t
rho_func = rho0 + eps_rho*sin(kphi*arg_phi)*sin(kz*arg_z)

!print *, ' t = ', t, ' phi = ', phi, ' z = ', z
!print *, ' arg_phi = ', arg_phi, ' arg_z = ', arg_z
!print *, ' kphi = ', kphi, ' kz = ', kz
!read(5, *)

end function rho_func

!----------------------------------------------------------------------------------85

subroutine fix_velocity(q)

use advection_test_params
use grid
use partition_data
use dof_indices
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(inout) :: q

! Local:
integer :: iz, iphi, ir

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         q(ir, iphi, iz, zmom) = q(ir, iphi, iz, irho) * uz_advec
         q(ir, iphi, iz, rmom) = 0.d0
         q(ir, iphi, iz, amom) = q(ir, iphi, iz, irho) * Omega_advec*rgrid(ir) * rgrid(ir) 
      end do
   end do
end do

return

end subroutine fix_velocity

!----------------------------------------------------------------------------------85

subroutine error(t, q, err_L2)

#ifdef mpi_code
   use mpi
#endif

use advection_test_params
use grid
use partition_data
use dof_indices

implicit none
real(8), intent(in) :: t
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in) :: q
real(8), intent(out) :: err_L2

! Local:
integer :: iz, ir, iphi, ier
real(8) :: sum_of_squares, my_sum_of_squares, rho_func

if (my_node .eq. 0) then
   print *, ' In subroutine error'
end if
   

ir = 1
my_sum_of_squares = 0.d0

if ( (ir .ge. sr) .and. (ir .le. er) ) then
   do iz = 1, nz
      do iphi = sphi, ephi
         my_sum_of_squares = my_sum_of_squares + (q(ir, iphi, iz, irho) - rho_func(t, phi_grid(iphi), zgrid(iz)))**2
      end do
   end do
end if

print *, ' my_node = ', my_node, ' my_sum_of_squares = ', my_sum_of_squares

#ifdef mpi_code
call mpi_reduce(my_sum_of_squares, sum_of_squares, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ier)
print *, ' returned from mpi_reduce'
#else
sum_of_squares = my_sum_of_squares
#endif

if (my_node .eq. 0) then
   err_L2 = sqrt(sum_of_squares / (nphi * nz))
   print *, ' In subroutine error, err_L2 = ', err_L2, ' sum_of_squares = ', sum_of_squares
end if

end subroutine error

!----------------------------------------------------------------------------------85
