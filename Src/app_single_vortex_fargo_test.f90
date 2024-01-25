!----------------------------------------------------------------------------------85

subroutine app_single_vortex_fargo_test

use dof_indices
use grid, only: ndof, nz, nr, nphi, rmin, rmax, zmin, zmax, zgrid, rgrid, &
     Ji_r, phi_min, phi_max

use boundary_condition_types
use boundary_condition_routines
use gravity, only: gz
use q_array, only: q
use transposes_of_q_and_qdot
use logical_units
use total_allocated_words, only: n_words_allocated
use partition_data
use cpu_timing_module
use set_up_routines
use activate_routines
use gravity_types
use basic_state
#ifdef mpi_code
   use mpi
#endif
implicit none

external rhs
real(8) :: pi   

! For restart:
logical :: restart
integer :: istep_of_restart_file, file_version
real(8) :: t_of_restart_file

integer :: istep, nsteps, istep1
real(8) :: t, cfl, dt, dt_previous

! Things set to unity:
real(8) :: T_orbital0, r0, Sigma0, Omega0

! Temps:
real(8) :: V_Kep0, uphi2, ci0

! Exponents:
real(8) :: q_Sigma, q_T, q_ci

! The only free parameter beside the exponents:
real(8) :: Keplerian_Mach

! For passing to subroutine set-up_thermal_parameters
real(8) :: gamma = 1.4d0
logical :: isothermal ! will be read from input file.

! For subroutine activate_gravity
logical :: gravity_flag
real(8) :: GM
real(8) :: gz_uniform ! Unused
integer :: gravity_type

! Field smoothing
logical :: apply_pade_filter = .true., filter_relative_to_basic_state = .true.
character(3) :: eps_or_tau = 'tau'
real(8) :: eps_filter, tau_filter

integer :: iz, ir, iphi, ier
integer :: tecplot_interval, profiles_interval, save_interval

! Functions of r:
real(8), allocatable, dimension(:) :: Sigma, uK, ci, Pbar

! Used to satisfy numerical radial balance:
real(8), allocatable, dimension(:) :: dPdr

logical :: apply_fargo_trick, integer_shifts, apply_fargo_correction

logical :: perturb, plot_curl_rho_u = .false.

logical :: plot_perturbation_field1 = .false., plot_vorticity_over_rho1 = .true., &
           plot_perturbation_field2 = .true.,  plot_vorticity_over_rho2 = .false.

logical :: suppress_z_derivatives_when_nz_not_1, periodic_z, use_supplied_dt

! For passing to set_up_boundary_conditions
integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced
real(8), allocatable, dimension(:) :: d_ci_dr_inner, d_ci_dr_outer
real(8) :: c_sound_rmin, c_sound_rmax

logical, parameter :: stretched_r = .false., stretched_z = .false.
real(8), parameter :: r0_unif = 0.d0, z0 = 0.d0
integer, parameter :: nr_u = 0, nz_u = 0
logical :: unphysical
integer :: iphi_bad
logical :: apply_newtonian_cooling
real(8) :: tau_newtonian_cooling
real(8), allocatable, dimension(:,:) :: ci_squared_initial
! For irrelevant arguments:
real(8) :: dummy = 0.d0
logical :: apply_artificial_pressure
real(8) :: C_ap ! Coefficient for artificial pressure

logical :: hit_target = .false., target_met
real(8) :: t_target

namelist /single_vortex_input/ restart, file_version, nr, nphi, nz, cfl, nsteps, tecplot_interval, &
     tau_filter, apply_artificial_pressure, C_ap, apply_fargo_trick, integer_shifts, &
     apply_fargo_correction, use_supplied_dt, dt, isothermal

if (my_node .eq. 0) print *, ' Running vertically_integrated_disk'

pi = 4.0d0 * atan(1.0d0)

! Read namelist input:
if (my_node .eq. 0) then
   print *, ' node 0: about to open and read namelist file for VID'      
   open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', &
        status = 'old')
   read (lun_general_purpose, nml = single_vortex_input)
   close (lun_general_purpose)

   print *, ' node 0: read namelist input file'
   print *, ' node 0: nr = ', nr, ' nphi = ', nphi
   print *, ' isothermal = ', isothermal
end if

#ifdef mpi_code
   call mpi_bcast(restart,                1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(file_version,           1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nr,                     1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nphi,                   1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nz,                     1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(cfl,                    1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(nsteps,                 1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(tecplot_interval,       1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(tau_filter,             1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(apply_fargo_trick,      1, mpi_logical,  0, mpi_comm_world, ier)
   call mpi_bcast(integer_shifts,         1, mpi_logical,  0, mpi_comm_world, ier)   
   call mpi_bcast(apply_fargo_correction, 1, mpi_logical,  0, mpi_comm_world, ier)   
   call mpi_bcast(use_supplied_dt,        1, mpi_logical,  0, mpi_comm_world, ier)
   call mpi_bcast(dt,                     1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(isothermal,             1, mpi_logical,  0, mpi_comm_world, ier)
#endif
   
if (my_node .eq. 0) print *, ' node 0: tau_filter = ', tau_filter

! Trick to get a parallel partitioning for an (r, phi) grid.
#ifdef mpi_code
   nz = 2  ! arrays will be artificially dimensioned larger and some work will be twice.
   suppress_z_derivatives_when_nz_not_1 = .true.
   zmin = 0.0
   zmax = 1.0
#else
   nz = 1
   suppress_z_derivatives_when_nz_not_1 = .false.
   zmin = 0.0
   zmax = 0.0
#endif

periodic_z = .false.
   
! We can set three things to unity:
T_orbital0   = 1.0d0
r0           = 1.0d0
Sigma0       = 1.0d0 ! Surface density

! The vertically integrated equations (so-called) look the same as the
! non-integrated ones with the identification:
! rho ---> Sigma
! p   ---> P_bar, vertically integrated pressure
! Hence in the coding below (and in the rest of the code) we sometimes
! use the notation on the left.

! Consequence of above:
Omega0 = 2.d0 * pi / T_orbital0
GM     = Omega0**2 * r0**3
V_Kep0 = SQRT(GM/r0)

! Gravity set-up parameters:
gravity_flag = .true.

! Exponents:
q_Sigma = -1.5d0  ! for surface density
q_T     = -0.5d0  ! for temperature or equivalently c_i^2
q_ci    = q_T * 0.5d0 ! for the isothermal sound-speed

! This is a parameter:
Keplerian_Mach = 30.0d0

! Therefore:
ci0 = V_Kep0/Keplerian_Mach

! Domain
rmin   = 0.5d0
rmax   = 1.5d0
phi_min = 0.0d0
phi_max = 2.d0 * pi

rmin_BC   = zero_normal_momentum
rmax_BC   = zero_normal_momentum
zmin_BC   = null
zmax_BC   = null
ibalanced = 0

! These three set_up calls should be in every application subroutine:
call set_up_domain_mesh_and_partition(rmin, rmax, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
     suppress_z_derivatives_when_nz_not_1, periodic_z, &
     stretched_r, stretched_z, r0_unif, nr_u, z0, nz_u)

allocate(Sigma(nr), uK(nr), ci(nr), Pbar(nr), dPdr(nr))

! Functions of r for the basic state:
do ir = 1, nr
   Sigma(ir) = Sigma0 * (rgrid(ir) / r0)**q_Sigma
   uK(ir)    = SQRT(GM/rgrid(ir))
   ci(ir)    = ci0 * (rgrid(ir)/r0)**q_ci
   ! Vertically integrated pressure from isothermal EOS:
   Pbar(ir) = Sigma(ir) * ci(ir)**2
end do

! Needed for the isothermal option and/or for isothermal wall bc.
! ci_squared_initial sits in module thermal_parameters.
allocate(ci_squared_initial(nr, nz))
do iz = 1, nz
   do ir = 1, nr
      ci_squared_initial(ir, iz) = ci(ir)**2
   end do
end do

apply_newtonian_cooling = .false.
tau_newtonian_cooling   = 0.d0
call set_up_thermal_parameters(gamma, isothermal, apply_newtonian_cooling, &
     tau_newtonian_cooling, ci_squared_initial)

allocate(d_ci_dr_inner(nz), d_ci_dr_outer(nz))
call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, &
     d_ci_dr_inner, d_ci_dr_outer, c_sound_rmin, c_sound_rmax, isothermal)

gravity_flag = .true.
gravity_type = mass_at_origin
gz_uniform   = 0.d0   ! Unused
call activate_gravity(gravity_flag, gravity_type, GM, gz_uniform)

call activate_pade_filter(apply_pade_filter, eps_or_tau, eps_filter, &
     tau_filter, filter_relative_to_basic_state)

if (apply_artificial_pressure) then
   if (my_node .eq. 0) print *, ' node 0: Calling activate_artificial_pressure with C_ap = ', C_ap
   if (my_node .eq. 1) print *, ' node 1: Calling activate_artificial_pressure with C_ap = ', C_ap      
   call activate_artificial_pressure (C_ap)
end if

if (apply_fargo_trick) then
   call activate_fargo(integer_shifts, apply_fargo_correction)
end if



! Basic state:
call pade_diff_bundle(1, nr, Ji_r, Pbar, dPdr)
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         q(ir, iphi, iz, irho) = Sigma(ir)
         ! Satisfy radial balance:
         uphi2 = uK(ir)**2 + rgrid(ir)/Sigma(ir) * dPdr(ir) 
         q(ir, iphi, iz, amom) = Sigma(ir) * SQRT(uphi2) * rgrid(ir)
         !q(ir, iphi, iz, amom) = 0.0d0 ! For testing pure vortex field
         q(ir, iphi, iz, rmom) = 0.0d0
         q(ir, iphi, iz, zmom) = 0.0d0
         if (.not. isothermal) q(ir, iphi, iz, ener) = Pbar(ir)/(gamma - 1.d0)         
      end do
   end do
end do

call store_basic_state(filter_relative_to_basic_state)

! Initial condition:
if (.not. restart) then
   ! Start with basic state then add vortex
   call add_vortex (GM)
end if

! These arrays are not needed anymore:
deallocate(Pbar)

if (.not. restart) then
   istep1 = 1
   t      = 0.0d0
else
   call read_restart (file_version, istep_of_restart_file, t_of_restart_file)
   istep1 = istep_of_restart_file + 1
   t      = t_of_restart_file
end if

call tecplot_horizontal_plane (q, 1, T_orbital0, t, istep)
call tecplot_vort_z_and_dil(q, t, plot_perturbation_field1, plot_vorticity_over_rho1, plot_curl_rho_u, 'style1')
call tecplot_vort_z_and_dil(q, t, plot_perturbation_field2, plot_vorticity_over_rho2, plot_curl_rho_u, 'style2')
call tecplot_baroclinic_term (q, t)

call sanity_check

if (my_node .eq. 0) print *, ' node 0: Mb allocated = ', float(8*n_words_allocated)/1.d6

apply_pade_filter = .true.
do istep = istep1, istep1 + nsteps - 1
   call rk4(rhs, q, cfl, t, dt, use_supplied_dt, unphysical, iphi_bad, &
                       hit_target, t_target, target_met)
   if (unphysical) call terminate_with_save(1, istep, t)   
   ! call euler(q, cfl, t, dt, use_supplied_dt)
   
   if (my_node .eq. 0) print *, ' node 0: finished istep = ', istep, ' t = ', t

   
   if (mod(istep, tecplot_interval ) .eq. 0) then
      call tecplot_horizontal_plane (q, 1, T_orbital0, t, istep)
      call tecplot_vort_z_and_dil(q, t, plot_perturbation_field1, plot_vorticity_over_rho1, plot_curl_rho_u, 'style1')
      call tecplot_vort_z_and_dil(q, t, plot_perturbation_field2, plot_vorticity_over_rho2, plot_curl_rho_u, 'style2')
      call tecplot_baroclinic_term (q, t)      
   end if

   ! Detect blowing-up:
   if (istep .gt. 1) then
      if (dt .lt. dt_previous/20.d0) then
         print *, ' dt decreased more than a factor of 20 in one step'
         print *, ' dt          = ', dt
         print *, ' dt_previous = ', dt_previous
         call terminate_with_save(1, istep, t)
      end if
   end if
   dt_previous = dt

   if (use_supplied_dt .and. (cfl .lt. 1.0d-4)) then
      print *, ' cfl has fallen to cfl = ', cfl
      call terminate_with_save(1, istep, t)
   end if
end do

100 continue
if (my_node .eq. 0) then
   print *, ' average stepping time per step = ', total_cpu_for_stepping/nsteps
end if

! Need the minus 1 since a fortran do loop increments the counter at the end.
call terminate_with_save(0, istep-1, t)

end subroutine app_single_vortex_fargo_test

!----------------------------------------------------------------------------------85

subroutine add_vortex (GM)

use grid
use dof_indices
use q_array
use math_constants, only: pi
use partition_data
implicit none
real(8) :: GM

integer :: iz, iphi, ir
real(8) :: x, y, xcen, ycen, dx, dy, s, theta, u_theta, arg, ux, uy, r, phi, &
     peak_vorticity, Keplerian_vorticity, rmid, Gamma, core_radius, Sigma, ur, &
     uphi
real(8) :: s_max, s_cut_off, lambda, kr, func ! debug
! Function called:
real(8) :: Melander

rmid = 0.5d0 * (rmin + rmax)
Keplerian_vorticity = 0.5d0 * SQRT(GM / rmid**3)
if (my_node .eq. 0) print *, ' Keplerian vorticity at rmid = ', Keplerian_vorticity

! Desired peak vorticity of the vortex:
peak_vorticity = -0.10d0 * Keplerian_vorticity
if (my_node .eq. 0) print *, ' peak vorticity of vortex = ', peak_vorticity

! Vortex center:
xcen = rmid
ycen = 0.0d0

! Parameters for 1/4 cosine cut-off to ensure zero perturbation at the
! boundaries:
! s_max = rmax - rmid
s_max = 0.8d0*(rmax - rmid)
s_cut_off   = 0.50d0*s_max
lambda      = 4.0d0*(s_max - s_cut_off)
kr          = 2.d0*pi/lambda

! Desired core radius:
core_radius = 0.40d0*s_max

! Resulting vortex circulation:
Gamma = peak_vorticity * pi * core_radius**2

print *, ' s_max = ', s_max

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         r   = rgrid(ir)
         phi = phi_grid(iphi)
         x = r * cos(phi)
         y = r * sin(phi)
         dx = x - xcen
         dy = y - ycen

         ! Polar coordinates of point relative to vortex center:
         s = SQRT(dx**2 + dy**2)
         theta = ATAN2(dy, dx)

         arg = -s**2/core_radius**2
         u_theta = Gamma / (2.d0*pi*s) * (1.d0 - exp(arg))

         if (s .lt. s_cut_off) then
            u_theta = u_theta
            ! print *, ' s<s_cut u_theta = ', u_theta
         else if (s .lt. s_max) then
            ! func = cos(kr*(s - s_cut_off))
            ! u_theta = u_theta*func
            u_theta = u_theta * Melander(s, s_cut_off, s_max)
            ! print *, ' s>s_cut: s = ', s, ' u_theta = ', u_theta, ' func = ', func
            ! print *, ' arg/pi = ', kr*(s - s_cut_off)/pi
         else
            u_theta = 0.0d0
            ! print *, ' s>s_max: s = ', s, ' u_theta = ', u_theta
         end if

         ! Test of solid body rotation.
         ! u_theta = s ! debug

         ux = -u_theta * sin(theta)
         uy =  u_theta * cos(theta)

         ur   =  ux*cos(phi) + uy*sin(phi)
         uphi = -ux*sin(phi) + uy*cos(phi)

         Sigma = q(ir, iphi, iz, irho)
         q(ir, iphi, iz, rmom) = q(ir, iphi, iz, rmom) +  Sigma*ur
         q(ir, iphi, iz, amom) = q(ir, iphi, iz, amom) +  Sigma*uphi*r         
      end do
   end do
end do

end subroutine add_vortex

!----------------------------------------------------------------------------------85









