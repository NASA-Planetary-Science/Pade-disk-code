!----------------------------------------------------------------------------------85

subroutine app_vortex_pair  ! Vortex pair in vertically integrated disk.

use set_up_routines
use dof_indices
use grid, only: ndof, nz, nr, nphi, rmin, rmax, zmin, zmax, zgrid, rgrid, &
     Ji_r, phi_min, phi_max
use boundary_condition_types
use gravity, only: gz
use q_array, only: q
use transposes_of_q_and_qdot
use logical_units
use total_allocated_words, only: n_words_allocated
use partition_data
use cpu_timing_module
use activate_routines
use gravity_types
use basic_state
#ifdef mpi_code
   use mpi
#endif
implicit none
external rhs
character(8) :: vortex_type   
   
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

! For passing to subroutine activate_gravity
logical :: gravity_flag
integer :: gravity_type
real(8) :: GM
real(8) :: gz_uniform ! Unused

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

logical :: perturb, plot_perturbation_field = .false., &
     plot_vorticity_over_rho = .true., &
     suppress_z_derivatives_when_nz_not_1, periodic_z, &
     use_supplied_dt, apply_artificial_pressure

real(8) :: C_ap

integer :: itype
real(8) :: Keplerian_vorticity, peak_vorticity, Gamma_Gaussian, sigma_Gaussian, &
     V_shear, u_theta_max_fac, u_theta_max

! Parameters for Seligman vortex:
real(8) :: k_delta, Gamma_Seligman ! For Seligman vortex

! For passing to subroutine set_up_boundary_conditions.
integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced
real(8), allocatable, dimension(:) :: d_ci_dr_inner, d_ci_dr_outer
real(8) :: c_sound_rmin, c_sound_rmax

integer :: n_wedges

logical :: apply_plotting_shift = .true., plot_curl_rho_u = .false.

logical, parameter :: stretched_r = .false., stretched_z = .false.
real(8), parameter :: r0_unif = 0.d0, z0 = 0.d0
integer, parameter :: nr_u = 0, nz_u = 0
logical :: unphysical
integer :: iphi_bad
! For irrelevant parameters:
real(8), parameter :: dummy = 0.d0

! For passing to subroutine initialize thermal set-up:
real(8) :: gamma = 1.4d0
logical :: isothermal ! will be read from input file.
real(8), allocatable, dimension(:,:) :: ci_squared_initial
logical :: apply_newtonian_cooling
real(8) :: tau_newtonian_cooling

logical :: hit_target = .false., target_met
real(8) :: t_target = 0.d0

namelist /vp_input/ restart, nr, nphi, nz, cfl, nsteps, tecplot_interval, &
     tau_filter, apply_fargo_trick, integer_shifts, &
     apply_fargo_correction, use_supplied_dt, dt, isothermal, &
     apply_artificial_pressure, C_ap

if (my_node .eq. 0) then
   print *, ' Running subroutine vortex_pair'
end if

vortex_type = 'Gaussian' !'Seligman' or 'Gaussian'
pi = 4.0d0 * atan(1.0d0)

! Read namelist input:
if (my_node .eq. 0) then
   print *, ' node 0: about to open and read namelist file for vortex_pair'      
   open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', &
        status = 'old')
   read (lun_general_purpose, nml = vp_input)
   close (lun_general_purpose)

   print *, ' node 0: read namelist input file'
   print *, ' node 0: nr = ', nr, ' nphi = ', nphi
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
   call mpi_bcast(apply_artificial_pressure, 1, mpi_logical,  0, mpi_comm_world, ier)
   call mpi_bcast(C_ap,                   1, mpi_double, 0, mpi_comm_world, ier)            
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
r0      = 1.0d0
Sigma0  = 1.0d0 ! Surface density
Omega0  = 1.0d0

ci0     = 0.05d0
k_delta = Omega0 / ci0

! The vertically integrated equations (so-called) look the same as the
! non-integrated ones with the identification:
! rho ---> Sigma
! p   ---> P_bar, vertically integrated pressure
! Hence in the coding below (and in the rest of the code) we sometimes
! use the notation on the left.

if (vortex_type .eq. 'Gaussian') then
   u_theta_max    = 0.025d0
   sigma_Gaussian = 0.10d0/k_delta

   u_theta_max_fac = (1.d0 - exp(-1.121**2)) / (2.d0 * pi * 1.121)
   Gamma_Gaussian = -u_theta_max * sigma_Gaussian / u_theta_max_fac
   Gamma_Gaussian = -0.04/30.d0
   print *, ' Gamma_Gaussian = ', Gamma_Gaussian
else if (vortex_type .eq. 'Seligman') then
   Gamma_Seligman = -0.01d0 / 30.d0
end if

T_orbital0  = 2.d0 * pi / Omega0
GM          = Omega0**2 * r0**3
V_Kep0      = SQRT(GM/r0)

print *, ' Max shear velocity = ', 1.5d0*Omega0*0.02d0
print *, ' Keplerian Mach number = ', V_Kep0/ci0


! Exponents:
q_Sigma = 0.0d0  ! for surface density
q_T     = 0.0d0  ! for temperature or equivalently c_i^2
q_ci    = 0.0d0  ! for the isothermal sound-speed

! Domain
n_wedges = 79 ! number of wedges around 2 pi.  A wedge is the computational domain.
rmin     = 0.96d0
rmax     = 1.04d0
phi_min  = 0.0d0
phi_max  = 2.d0*pi/n_wedges  ! approx 0.04 radians

print *, ' Omega0       = ', Omega0
print *, ' GM           = ', GM

rmin_BC   = zero_normal_momentum
rmax_BC   = zero_normal_momentum
zmin_BC   = null
zmax_BC   = null
ibalanced = 0

! These three set-up calls should be in every application subroutine:
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

call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, &
     d_ci_dr_inner, d_ci_dr_outer, c_sound_rmin, c_sound_rmax, isothermal)

gravity_flag = .true.
gravity_type = mass_at_origin
gz_uniform   = 0.d0   
call activate_gravity(gravity_flag, gravity_type, GM, gz_uniform)   

call activate_pade_filter(apply_pade_filter, eps_or_tau, eps_filter, &
     tau_filter, filter_relative_to_basic_state)

if (apply_fargo_trick) then
   call activate_fargo(integer_shifts, apply_fargo_correction)
end if

! This is also optional.
if (apply_plotting_shift) then
   call activate_plotting_shift(Omega0, t)
end if

if (apply_artificial_pressure) then
   print *, ' About to call activate_artificial_pressure'
   call activate_artificial_pressure (C_ap)
   print *, ' Returned from activate artificial pressure'
end if

print *, ' about to differentiate Pbar'
call pade_diff_bundle(1, nr, Ji_r, Pbar, dPdr)

! Initial condition.  For restart, we still calculate the basic state because
! we are filtering relative to the basic state:
! Start with basic state then add vortex:
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         q(ir, iphi, iz, irho) = Sigma(ir)
         ! Satisfy radial balance:
         uphi2 = uK(ir)**2 + rgrid(ir)/Sigma(ir) * dPdr(ir) 
         q(ir, iphi, iz, amom) = Sigma(ir) * SQRT(uphi2) * rgrid(ir)
         ! q(ir, iphi, iz, amom) = 0.0d0 ! For testing pure vortex field
         q(ir, iphi, iz, rmom) = 0.0d0
         q(ir, iphi, iz, zmom) = 0.0d0
         if (.not. isothermal) q(ir, iphi, iz, ener) = Pbar(ir)/(gamma - 1.d0)
      end do
   end do
end do

call store_basic_state(filter_relative_to_basic_state)

if (.not. restart) then
   call add_vortex_pair (vortex_type, Gamma_Gaussian, sigma_Gaussian, GM, Gamma_Seligman, &
        k_delta, n_wedges)
end if

if (restart) then
   call read_restart (file_version, istep_of_restart_file, t_of_restart_file)
   istep1 = istep_of_restart_file + 1
   t      = t_of_restart_file
else
   istep1 = 1
   t      = 0.0d0
end if

! These arrays are not needed anymore:
deallocate(Pbar)

! call tecplot_horizontal_plane (q, 1, T_orbital0, t, istep)
plot_perturbation_field = .false.
plot_vorticity_over_rho = .true.
call tecplot_vort_z_and_dil(q, t, plot_perturbation_field, plot_vorticity_over_rho, plot_curl_rho_u, 'style3')
! plot_perturbation_field = .false.
! plot_vorticity_over_rho = .false.
! call tecplot_vort_and_dil(q, t, plot_perturbation_field, plot_vorticity_over_rho, plot_curl_rho_u, 'style2')
! call tecplot_baroclinic_term (q, t)

call sanity_check

if (my_node .eq. 0) print *, ' node 0: Mb allocated = ', float(8*n_words_allocated)/1.d6

do istep = istep1, istep1 + nsteps - 1
   call rk4(rhs, q, cfl, t, dt, use_supplied_dt, unphysical, iphi_bad, &
                  hit_target, t_target, target_met)
   ! call euler(q, cfl, t, dt, use_supplied_dt)
   
   if (my_node .eq. 0) print *, ' node 0: finished istep = ', istep, ' t = ', t
   
   if (mod(istep, tecplot_interval ) .eq. 0) then
      if (apply_plotting_shift) call plotting_shift(q, t)      
      ! call tecplot_horizontal_plane (q, 1, T_orbital0, t, istep)
      plot_perturbation_field = .false.
      plot_vorticity_over_rho = .true.
      call tecplot_vort_z_and_dil(q, t, plot_perturbation_field, plot_vorticity_over_rho, plot_curl_rho_u, 'style3')
      ! plot_perturbation_field = .false.
      ! plot_vorticity_over_rho = .false.
      ! call tecplot_vort_and_dil(q, t, plot_perturbation_field, plot_vorticity_over_rho, plot_curl_rho_u, 'style2')

      ! if (apply_fargo_trick) call fargo_informational_output (t)      
      ! call tecplot_baroclinic_term (q, t)      
   end if

   if (apply_fargo_trick .and. (istep .eq. istep1)) call fargo_informational_output (t)         

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

! Need the minus 1 since a fortran do loop increments the counter at the end.

if (my_node .eq. 0) then
   print *, ' average stepping time per step = ', total_cpu_for_stepping/nsteps
end if

call terminate_with_save(0, istep-1, t)

end subroutine app_vortex_pair

!----------------------------------------------------------------------------------85

subroutine add_vortex_pair (vortex_type, Gamma_Gaussian, sigma_Gaussian, GM, &
                         Gamma_Seligman, k_delta, n_wedges)

use grid
use dof_indices
use q_array
use math_constants, only: pi
use partition_data
implicit none
integer, parameter :: n_points = 400
real(8), dimension(n_points) :: svec, u_theta_vec

character(8) :: vortex_type
real(8)      :: Gamma_Gaussian, sigma_Gaussian, GM, Gamma_Seligman, k_delta
integer      :: n_wedges

integer :: iz, iphi, ir
integer :: i_wedge, ivort

real(8), dimension(2, n_wedges) :: rcen, phi_cen, xcen, ycen

! Melander cut-off parameters:
real(8) :: s_cut_off, s_max_Melander

real(8) :: x, y, dx, dy, s, smid, theta, u_theta, uphi, ur, ux, uy, Sigma, r, phi, vort

! Function called:
real(8) :: u_theta_func

real(8) :: delta, phi_mid, delta_phi_wedge

! For plotting u_theta profile:
real(8) :: smin_plot, smax_plot, ds
integer :: i, n

! Two-vortex case:
s_cut_off      = 0.006d0
s_max_Melander = 0.012d0
smin_plot = 0.00001d0
smax_plot = 0.02d0

delta_phi_wedge = 2.0d0 * pi / n_wedges
phi_mid         = 0.5d0 * delta_phi_wedge
delta = 0.05d0 / k_delta
print *, ' delta = ', delta

!!$do i_wedge = 1, n_wedges 
!!$   rcen   (1, i_wedge)  = 1.0d0   + 0.0025d0
!!$   phi_cen(1, i_wedge)  = phi_mid + 0.0025d0 + (i_wedge - 1)*delta_phi_wedge
!!$   rcen   (2, i_wedge)  = 1.0d0   - 0.0025d0
!!$   phi_cen(2, i_wedge)  = phi_mid - 0.0025d0 + (i_wedge - 1)*delta_phi_wedge
!!$end do
!!$s_cut_off      = 0.010d0
!!$s_max_Melander = 0.02d0 - 0.0025d0



do i_wedge = 1, n_wedges 
   rcen   (1, i_wedge)  = 1.0d0   + 0.0050d0
   phi_cen(1, i_wedge)  = phi_mid + 0.0100d0 + (i_wedge - 1)*delta_phi_wedge
   rcen   (2, i_wedge)  = 1.0d0   - 0.0050d0
   phi_cen(2, i_wedge)  = phi_mid - 0.0100d0 + (i_wedge - 1)*delta_phi_wedge
end do

do i_wedge = 1, n_wedges
   do ivort = 1, 2
      xcen(ivort, i_wedge) = rcen(ivort, i_wedge)*cos(phi_cen(ivort, i_wedge))
      ycen(ivort, i_wedge) = rcen(ivort, i_wedge)*sin(phi_cen(ivort, i_wedge))
   end do
end do

! Profile plots:
ds   = (smax_plot - smin_plot) / (n_points - 1)

open (unit = 1, file = 'u_theta.dat', form = 'formatted', status = 'unknown')
do i = 1, n_points
   svec(i) = smin_plot + (i - 1)*ds
   u_theta_vec(i) = u_theta_func(vortex_type, pi, Gamma_Gaussian, sigma_Gaussian, &
        Gamma_Seligman, k_delta, delta, svec(i), s_cut_off, s_max_Melander)               
   write (1, "(2(1x, e12.5))") svec(i), abs(u_theta_vec(i))
end do
close (1)

open (unit = 1, file = 'vort.dat', form = 'formatted', status = 'unknown')
do i = 1, n_points - 1
   smid = 0.5d0*(svec(i) + svec(i+1))
   vort = 1.d0/smid * (svec(i+1)*u_theta_vec(i+1) - svec(i)*u_theta_vec(i)) / ds
   write (1, "(2(1x, e12.5))") smid, vort
end do
close (1)


do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         r   = rgrid(ir)
         phi = phi_grid(iphi)
         x = r * cos(phi)
         y = r * sin(phi)

         ux = 0.0d0
         uy = 0.0d0
         do ivort = 1, 2
            do i_wedge = 1, n_wedges
               dx = x - xcen(ivort, i_wedge)
               dy = y - ycen(ivort, i_wedge)
               
               ! Polar coordinates of point relative to vortex center:
               s = SQRT(dx**2 + dy**2)
               theta = atan2(dy, dx)

               u_theta = u_theta_func(vortex_type, pi, Gamma_Gaussian, sigma_Gaussian, &
                    Gamma_Seligman, k_delta, delta, s, s_cut_off, s_max_Melander)               

               ! Test of solid body rotation.
               ! u_theta = s ! debug
               ux = ux - u_theta * sin(theta)
               uy = uy + u_theta * cos(theta)            
            end do
         end do

         ur   =  ux*cos(phi) + uy*sin(phi)
         uphi = -ux*sin(phi) + uy*cos(phi)

         Sigma = q(ir, iphi, iz, irho)
         q(ir, iphi, iz, rmom) = q(ir, iphi, iz, rmom) +  Sigma*ur
         q(ir, iphi, iz, amom) = q(ir, iphi, iz, amom) +  Sigma*uphi*r         
      end do
   end do
end do

end subroutine add_vortex_pair

!----------------------------------------------------------------------------------85

real(8) function u_theta_func(vortex_type, pi, Gamma_Gaussian, sigma_Gaussian, &
     Gamma_Seligman, k_delta, delta, s, s_cut_off, s_max_Melander)

implicit none
character(8) :: vortex_type
real(8)      :: pi, Gamma_Gaussian, sigma_Gaussian, Gamma_Seligman, k_delta, delta, s, &
     s_cut_off, s_max_Melander

! Called functons:
real(8) :: seligman_u_theta, Melander3, gaussian_u_theta

if (vortex_type .eq. 'Seligman') then
!   u_theta_func = seligman_u_theta (pi, Gamma_Seligman, k_delta, delta, s) * &
!        Melander3(s, s_cut_off, s_max_Melander)
   u_theta_func = seligman_u_theta (pi, Gamma_Seligman, k_delta, delta, s)
else if (vortex_type .eq. 'Gaussian') then
   u_theta_func = gaussian_u_theta (pi, Gamma_Gaussian, sigma_Gaussian, s) * &
        Melander3(s, s_cut_off, s_max_Melander)
!     u_theta_func = gaussian_u_theta (pi, Gamma_Gaussian, sigma_Gaussian, s)
else
   ! To avoid "May be unitialized warning" from the compiler.
   u_theta_func = 0.0d0
   print *, ' Unknown type'
   call terminate_with_no_save(1)
end if

return
end function u_theta_func

!----------------------------------------------------------------------------------85

real(8) function gaussian_u_theta(pi, Gamma, sigma_Gaussian, s)

real(8) :: s, pi, Gamma, sigma_Gaussian

! Locals:
real(8) :: arg

arg = -s**2/sigma_Gaussian**2
gaussian_u_theta = Gamma / (2.d0*pi*s) * (1.d0 - exp(arg))

return
end

!----------------------------------------------------------------------------------85

real(8) function seligman_u_theta (pi, Gamma, k_delta, delta, s)

! delta = core radius, s = distance from vortex center.
real(8) :: pi, Gamma, k_delta, delta, s

! Local:
integer :: want, status
real(8) :: k_delta_s, k_delta_delta, I1_s, K1_s, I1_delta, K1_delta

print *, ' Seligman vortex currently commented out'
stop

want = 3 ! compute both

k_delta_s = k_delta * s
!!!call DBI1K1 (k_delta_s, I1_s, K1_s, want, status)
if (status .ne. 0) then
   print *, ' DBI1K1 returned status = ', status
   print *, ' want = ', want
   call terminate_with_no_save(1)
end if

k_delta_delta = k_delta * delta
!!!call DBI1K1 (k_delta_delta, I1_delta, K1_delta, want, status)
if (status .ne. 0) then
   print *, ' DBI1K1 returned status = ', status
   print *, ' want = ', want
   call terminate_with_no_save(1)
end if

if (s .le. delta) then
   seligman_u_theta = Gamma/(pi * delta) * K1_delta * I1_s * exp(-100.d0 * s)
else
   seligman_u_theta = Gamma/(pi * delta) * I1_delta * K1_s * exp(-100.d0 * s)
end if

end function seligman_u_theta

!----------------------------------------------------------------------------------85








