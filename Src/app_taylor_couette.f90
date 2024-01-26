!----------------------------------------------------------------------------------85

module taylor_couette_params
! C1 and C2 are constants for the exact solution in the case when the cylindrical
! walls are moving axially. Subscript i: inner wall, o: outer wall.
! rho0 : initial uniform density.
real(8) :: Omega_i, Omega_o, Ui, uz_i, uz_o, r_i, r_o, eta, A_tc, B_tc, C1, C2, &
     rho0, c_sound_i, gamma, Re_d
real(8) :: eps_pert
integer :: itype, isubtype
integer, parameter :: marcus = 1, dong = 2, axial = 3, meyer = 4, moser = 5, &
     check_viscous = 6, dong_counter_rotating = 7, laminar_moser = 8
! Sub-types:
integer, parameter :: moser_meyer = 1, moser_DS_narrow = 2, moser_DS_wide = 3 
end module taylor_couette_params

!----------------------------------------------------------------------------------85

subroutine app_taylor_couette

use grid, only: rgrid, zmin, zmax, iz_mid
use taylor_couette_params
use dof_indices
use boundary_condition_types 
use q_array, only: q
use logical_units
use control_parameters, only: restart_file, save_file
use total_allocated_words, only: n_words_allocated
use math_constants, only: pi
use partition_data
use partition_data_for_alan, only: ng1, ng2
use viscosity_types
use set_up_routines
use activate_routines
use basic_state
#ifdef mpi_code
   use mpi
#endif
implicit none

! Function called:
external rhs
real(8) :: tc_density

logical :: restart
integer :: istep_of_restart_file, file_version
real(8) :: t_of_restart_file

integer :: nr, nz, nphi
integer :: istep, nsteps, istep0
real(8) :: t, cfl, dt, dt_previous

! Pade filter stuff:
logical :: apply_pade_filter, filter_relative_to_basic_state = .false.
character(3) :: eps_or_tau
real(8) :: eps_filter, tau_filter

! For thermal set-up:
logical :: isothermal = .true.
logical :: apply_newtonian_cooling = .false.
real(8) :: Pr_molecular = 0.7d0 ! Not needed for isothermal.
real(8) :: tau_newtonian_cooling = 0.d0 ! dummy

! For gravity set-up
real(8) :: GM
integer :: gravity_flag=0, i_thin, i_no_radial

integer :: iz, ir, iphi, ier
integer :: tecplot_interval, profiles_interval, save_interval, history_interval

logical :: perturb, suppress_z_derivatives_when_nz_not_1_arg = .false., periodic_z_arg = .true., &
     use_supplied_dt = .false., name_using_step = .false., output_profiles

! For passing to subroutine set_up_boundary_conditions.
integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced
real(8), allocatable, dimension(:) :: d_ci_dr_inner, d_ci_dr_outer
real(8) :: c_sound_rmin, c_sound_rmax

real(8) :: uz_rms, ur_rms

real(8) :: Mach_i, d, Re, nu_molecular, nu_b_molecular, mu_molecular, Uo, &
     bracket_i, rho, p
real(8) :: c_sound, c_sound_o
real(8) :: kr, eps_profile, uz_prime
real(8) :: phi_min, phi_max, zmax_over_d
real(8) :: L_scale = 1.0d0, T_scale = 1.0d0

real(8) :: Re_Dong, Re_c
logical :: apply_viscosity
integer :: viscosity_type
real(8) :: C_DDSV ! unused
real(8) :: C_Smag

logical, parameter :: stretched_r = .false., stretched_z = .false.
real(8), parameter :: r0 = 0.d0, z0 = 0.d0
integer, parameter :: nr_u = 0, nz_u = 0

character(25) :: filename_given
logical :: unphysical
integer :: iphi_bad
logical, parameter :: average_in_z = .true. ! For Reynolds averages

! Logical units:
integer :: lun_torque
real(8) :: det ! determinant for laminar solution

! z domain size for Dong counter-rotating:
real(8) :: Lz_over_pi
real(8) :: pressure
real(8), allocatable, dimension(:,:) :: ci_squared_initial

logical :: hit_target = .false., target_met
real(8) :: t_target

namelist /tc_input/ itype, isubtype, Re_Dong, restart, file_version, isothermal, nr, nz, nphi, &
     Lz_over_pi, cfl, nsteps, tecplot_interval, profiles_interval, save_interval, history_interval,&
     perturb, apply_pade_filter, eps_or_tau, eps_filter, tau_filter

lun_torque = lun_history(1)

pi = 4.d0*atan(1.d0)

if (my_node .eq. 0) then
   print *, ' node 0: First executable of taylor_couette'
   print *, ' my_node = ', my_node, ' num_nodes = ', num_nodes
end if

#ifdef mpi_code
   call mpi_barrier(mpi_comm_world, ier)
#endif

! Read namelist input:
if (my_node .eq. 0) then
   print *, ' node 0: about to open and read namelist file for Taylor-Couette Flow Test Case'      
   open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', status = 'old')
   read (lun_general_purpose, nml = tc_input)
   close (lun_general_purpose)

   print *, ' node 0: read namelist input file'
   print *, ' node 0: nr = ', nr, ' nz = ', nz, ' nphi = ', nphi
end if

#ifdef mpi_code
   call mpi_bcast(itype,                     1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(isubtype,                  1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(Re_Dong,                   1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(restart,                   1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(file_version,              1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(isothermal,                1, mpi_logical, 0, mpi_comm_world, ier)      
   call mpi_bcast(nr,                        1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nz,                        1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nphi,                      1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(Lz_over_pi,                1, mpi_double,  0, mpi_comm_world, ier)   
   call mpi_bcast(cfl,                       1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(nsteps,                    1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(tecplot_interval,          1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(profiles_interval,         1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(save_interval,             1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(history_interval,          1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(perturb,                   1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(apply_pade_filter,         1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(eps_or_tau,                3, mpi_character, 0, mpi_comm_world, ier)
   call mpi_bcast(eps_filter,                1, mpi_double,  0, mpi_comm_world, ier)   
   call mpi_bcast(tau_filter,                1, mpi_double,  0, mpi_comm_world, ier)
#endif

if (my_node .eq. 0) then
   print *, ' my_node = ', my_node
   print *, ' tau_filter = ', tau_filter
   print *, ' file_version = ', file_version
   print *, ' isothermal   = ', isothermal
end if

pi = 4.0*ATAN(1.0d0)
! For thermal set up and initial density profile:
gamma = 1.4d0
eps_pert = 0.01d0

allocate(d_ci_dr_inner(nz), d_ci_dr_outer(nz))

if (itype .eq. marcus) then
   ! Marcus (1984), Sec. 3: Axisymmetric Taylor-vortex flow.  His Figure 2.
   nphi = 1
   phi_min = 0.0d0
   phi_max = 2.d0*pi
   
   ! Can set three quantities to unity:
   rho0  = 1.0d0
   d     = 1.0d0
   Ui    = 1.0d0

   eta = 0.875d0
   Re_c = 118.16
   Re  = 1.179d0 * Re_c
   !Re = 0.8d0*Re_c
   zmax_over_d = 2.5d0*d  
   
   r_o = d/(1.d0 - eta)
   r_i = r_o - d   

   ! For passing to boundary conditions:
   uz_i = 0.0d0
   uz_o = 0.0d0
   Omega_i = Ui / r_i
   Omega_o = 0.0d0   

   ! This is a parameter.  Mach_i = Omega_i * r_i / c_sound_i
   Mach_i = 0.10d0

   zmin = 0.0d0
   zmax = zmax_over_d * d

   ! Temperature BC for non-isothermal BC:
   c_sound_i = Ui / Mach_i
   c_sound_o = c_sound_i

   A_tc = -Omega_i * eta**2 / (1.d0 - eta**2)
   B_tc =  Omega_i * r_i**2 / (1.d0 - eta**2)

   ! nu_molecular is passed to activate_viscosity
   nu_molecular = Ui*d/Re

#ifndef mpi_code
   print *, ' A_tc = ', A_tc, ' B_tc = ', B_tc
#endif
else if (itype .eq. dong) then
   ! Dong (2007) standard turbulent Taylor-Couette flow with
   ! only the inner cylinder rotating.
   phi_min = 0.0d0
   phi_max = 2.d0*pi   

   ! Can set three quantities to unity:
   rho0 = 1.0d0
   d    = 1.0d0
   Ui   = 1.0d0
   
   ! For passing to viscous wall BC:
   uz_i = 0.d0
   uz_o = 0.d0

   ! Domain:
   R_i = 1.d0
   R_o = R_i + d
   eta = R_i / R_o
   
   zmin = 0.d0
   zmax = 2.d0*pi

   ! Speeds:
   Uo = 0.d0

   Omega_i = Ui / r_i ! Passed to viscous wall BC
   Omega_o = 0.0d0     ! Passed to viscous wall BC
   
   Mach_i = 0.1d0

   ! Passed to viscous wall BC:
   c_sound_i = Ui / Mach_i
   c_sound_o = c_sound_i   

   A_tc = -Omega_i * eta**2 / (1.d0 - eta**2)
   B_tc =  Omega_i * r_i**2 / (1.d0 - eta**2)

   ! nu is passed to activate_viscosity
   nu_molecular = Ui*d/Re_Dong
else if (itype .eq. axial) then
   ! Axial wall motion:

   ! No rotational motion:
   Omega_i = 0.0d0
   Omega_o = 0.0d0
   
   rho0  = 1.0d0
   uz_o  = 1.0d0
   r_i   = 1.0d0

   r_o       = 1.5d0
   uz_i      = 0.0d0
   c_sound_i = 10.0d0 ! This makes the Mach # = 0.10
   c_sound_o = 10.0d0
   nu_molecular = 1.d0/20.d0
   mu_molecular = rho0 * nu_molecular

   zmin = 0.0d0
   zmax = 2.5d0

   C1 = uz_o / log(r_o/r_i)
   C2 = -C1 * log(r_i)

   ! Amplitude of sinusoidal perturbation on laminar profile:
   eps_profile = 0.10d0
   kr  = 2.d0 * pi / (r_o - r_i)

   print *, ' itype = ', itype
   print *, ' r_i = ', r_i, ' r_o = ', r_o
   print *, ' zmin = ', zmin, ' zmax = ', zmax
else if (itype .eq. meyer) then
   ! Meyer-Spasche and Keller (Narrow gap)
   nphi = 1
   phi_min = 0.0d0
   phi_max = 0.0d0

   ! Setting three things:
   Ui   = 1.d0  ! cm/s
   rho0 = 1.0d0
   d    = 0.1d0 ! cm

   zmin = 0.d0
   zmax_over_d = 2.007d0
   zmax = zmax_over_d * d

   r_i = 1.9d0
   r_o = r_i + d

   Re_d = 40.d0
   
   ! For passing to boundary conditions:
   uz_i = 0.0d0
   uz_o = 0.0d0

   omega_i = Ui / r_i
   omega_o = 0.d0

   ! This is a parameter.  Mach_i = Omega_i * r_i / c_sound_i
   Mach_i = 0.10d0

   c_sound_i = Ui / Mach_i
   c_sound_o = c_sound_i

   eta  = r_i/r_o
   A_tc = -Omega_i * eta**2 / (1.d0 - eta**2)
   B_tc =  Omega_i * r_i**2 / (1.d0 - eta**2)

   ! nu_molecular is passed to activate_viscosity
   nu_molecular = Ui*d/Re_d
else if (itype .eq. moser) then
   nphi = 1
   phi_min = 0.d0
   phi_max = 0.d0

   ! Setting three things to unity:
   Ui   = 1.0d0
   rho0 = 1.0d0
   d    = 1.0d0

   ! Three parameters:
   if (isubtype .eq. moser_meyer) then
      eta = 5.d0 / 6.d0
      zmax_over_d = 1.05d0
      Re_d = 400.d0
   else if (isubtype .eq. moser_DS_narrow) then
      eta = 0.95d0
      zmax_over_d = 2.009d0
      Re_d = 195.d0
   else if (isubtype .eq. moser_DS_wide) then      
      eta = 0.5d0
      zmax_over_d = 1.988d0
      Re_d = 78.8d0
   end if
   
   zmin = 0.d0
   zmax = zmax_over_d * d

   ! This assumes that d = 1:
   r_i = 1.d0 / (1.d0/eta - 1.d0)
   r_o = r_i + d
   
   ! For passing to boundary conditions:
   uz_i = 0.0d0
   uz_o = 0.0d0

   Omega_i = Ui / r_i
   Omega_o = 0.d0

   ! This is a parameter.  Mach_i = Omega_i * r_i / c_sound_i
   Mach_i = 0.05d0

   c_sound_i = Ui / Mach_i
   c_sound_o = c_sound_i

   eta  = r_i/r_o
   A_tc = -Omega_i * eta**2 / (1.d0 - eta**2)
   B_tc =  Omega_i * r_i**2 / (1.d0 - eta**2)

   ! nu_molecular is passed to activate_viscosity
   nu_molecular = Ui*d/Re_d
else if (itype  .eq. laminar_moser) then
   nphi = 1
   phi_min = 0.d0
   phi_max = 0.d0

   ! Setting three things to unity:
   Ui   = 1.0d0
   rho0 = 1.0d0
   d    = 1.0d0

   ! Three parameters:
   if (isubtype .eq. moser_meyer) then
      eta = 5.d0 / 6.d0
      zmax_over_d = 1.05d0
      Re_d = 20.d0
   else if (isubtype .eq. moser_DS_narrow) then
      eta = 0.95d0
      zmax_over_d = 2.009d0
      Re_d = 20.d0
   else if (isubtype .eq. moser_DS_wide) then      
      eta = 0.5d0
      zmax_over_d = 1.988d0
      Re_d = 20.d0
   end if
   
   zmin = 0.d0
   zmax = zmax_over_d * d

   ! This assumes that d = 1:
   r_i = 1.d0 / (1.d0/eta - 1.d0)
   r_o = r_i + d
   
   ! For passing to boundary conditions:
   uz_i = 0.0d0
   uz_o = 0.0d0

   Omega_i = Ui / r_i
   Omega_o = 0.d0

   ! This is a parameter.  Mach_i = Omega_i * r_i / c_sound_i
   Mach_i = 0.10d0

   c_sound_i = Ui / Mach_i
   c_sound_o = c_sound_i

   eta  = r_i/r_o
   A_tc = -Omega_i * eta**2 / (1.d0 - eta**2)
   B_tc =  Omega_i * r_i**2 / (1.d0 - eta**2)

   ! nu_molecular is passed to activate_viscosity
   nu_molecular = Ui*d/Re_d   
else if (itype .eq. check_viscous) then
   zmin = 0.d0
   zmax = 1.d0
   r_i  = 1.d0
   r_o  = 2.d0
   phi_min = 0.d0
   phi_max = 2.d0*pi
   nu_molecular = 1.d0
   print *, ' check viscous: nphi = ', nphi
else if (itype .eq. dong_counter_rotating) then
   ! Can set three things to unity.
   rho0 = 1.d0
   d    = 1.d0
   Ui   = 1.d0

   Uo = -Ui

   ! Domain:
   ! Radius ratio is 1/2:
   r_i = 1.d0
   r_o = r_i + d ! = 2   
   zmin = 0.d0
   ! Lz/d = 2*pi
   zmax = 2.d0 * pi
   phi_min = 0.d0
   phi_max = 2.d0*pi

   ! This is a parameter.  Mach_i = Omega_i * r_i / c_sound_i
   Mach_i    = 0.10d0
   c_sound_i = Ui / Mach_i
   c_sound_o = c_sound_i
   
   ! For passing to boundary conditions:
   Omega_i = Ui / r_i
   Omega_o = Uo / r_o   
   uz_i = 0.0d0 ! No axial motion of the cylinders
   uz_o = 0.0d0
   
   nu_molecular = Ui*d/Re_Dong

   ! Coeffs of the laminar flow uphi = Ar + B/r
   det  = r_i/r_o - r_o/r_i
   A_tc = (Ui/r_o - Uo/r_i) / det
   B_tc = r_i*Uo - r_o*Ui
end if

! Gravity set-up parameters:
gravity_flag = 0
i_thin       = 0
i_no_radial  = 0

rmin_BC = viscous_wall
rmax_BC = viscous_wall
zmin_BC = periodic
zmax_BC = periodic
ibalanced = 0
periodic_z_arg = .true. ! For mesh set-up.

! These calls should be in every application subroutine:
call set_up_domain_mesh_and_partition(r_i, r_o, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
     suppress_z_derivatives_when_nz_not_1_arg, periodic_z_arg, &
     stretched_r, stretched_z, r0, nr_u, z0, nz_u)

! Needed for the isothermal option.  ci_squared_initial sits in
! module thermal_parameters.
if (isothermal) then
   allocate(ci_squared_initial(nr,nz))
   do iz = 1, nz
      do ir = 1, nr ! Note each processor has all of r.
         ci_squared_initial(ir, iz) = c_sound_i**2
      end do
   end do
end if

call set_up_thermal_parameters(gamma, isothermal, apply_newtonian_cooling, tau_newtonian_cooling, &
     ci_squared_initial)

if (my_node .eq. 0) then
   print *, ' my_node = 0, app_taylor_couette: returned from set_up_thermal_parameters'
end if

! Note: c_sound_i and c_sound_o are used only for viscous_wall BC at rmin and rmax
! and the non-isothermal case.
call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, &
     d_ci_dr_inner, d_ci_dr_outer, c_sound_i, c_sound_o, isothermal)

if (my_node .eq. 0) then
   print *, ' my_node = 0, app_taylor_couette: returned from set_up_boundary_conditions'
end if

! I am now doing an isothermal calculation.
! This gives c_sound_o based on a homentropic initial condition:
!call tc_laminar_profiles(r_o, rho, p, c_sound_o)

! This call is required for viscous wall BC:
call set_up_viscous_wall_conditions(Omega_i, Omega_o, uz_i, uz_o, c_sound_i, c_sound_o)

call activate_pade_filter(apply_pade_filter, eps_or_tau, eps_filter, &
     tau_filter, filter_relative_to_basic_state)

apply_viscosity = .true.
viscosity_type  = molecular
! isothermal has been set in the declarations
! nu_molecular has been set above
nu_b_molecular = 0.d0 ! Bulk kinematic viscosity
! Pr_molecular has been set in the declarations
! gamma has been set above.
C_DDSV = 0.0d0 ! Unused
C_Smag = 0.d0 ! unused
print *, ' nu_molecular = ', nu_molecular
call activate_viscosity(apply_viscosity, viscosity_type, isothermal, &
     nu_molecular, nu_b_molecular, Pr_molecular, gamma, C_DDSV, C_Smag)

if (itype .eq. check_viscous) then
   ! Currently this works in serial mode only
   call check_viscous_with_phi
   stop
end if

! Initial condition/basic state:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         if (itype .ne. axial) then
            !call tc_laminar_profiles(rgrid(ir), rho, p, c_sound)
            q(ir, iphi, iz, irho) = rho0
            q(ir, iphi, iz, zmom) = 0.0d0
            q(ir, iphi, iz, rmom) = 0.0d0
            q(ir, iphi, iz, amom) = rho0 * rgrid(ir) * (A_tc*rgrid(ir) + B_tc/rgrid(ir))
            ! In case we are running non-isothermal:
            pressure = rho0*c_sound_i**2.d0
            q(ir, iphi, iz, ener) = pressure / (gamma - 1.d0)
            !print *, ' ir = ', ir, ' iphi = ', iphi, ' iz = ', iz
            !print *, ' rho0 = ', rho0
            !print *, ' q(ir, iphi, iz, amom) = ', q(ir, iphi, iz, amom)
            !read(5, *)
         else if (itype .eq. axial) then
            q(ir, iphi, iz, irho) = rho0
            uz_prime = eps_profile*sin(kr*rgrid(ir))
            q(ir, iphi, iz, zmom) = rho0 * (C1*log(rgrid(ir)) + C2 + uz_prime) 
            q(ir, iphi, iz, rmom) = 0.0d0
            q(ir, iphi, iz, amom) = 0.0d0
         end if
      end do
   end do
end do
call store_basic_state(filter_relative_to_basic_state)

if (restart) then
   ! file_version (1 or 2) is read from the input file.
   call read_restart (file_version, istep_of_restart_file, t_of_restart_file)
end if

if (.not. restart) then
   istep0 = 0
   t      = 0.0d0
else
   istep0 = istep_of_restart_file
   t      = t_of_restart_file
end if

if (perturb) call add_perturbation_tc

call tecplot_meridional_plane(q, 1, 1.0d0, 1.0d0, t, istep0)
call output_torque(t, lun_torque)

call sanity_check
call laminar_profile_output

do istep = istep0 + 1, istep0 + nsteps   
   !print *, ' Before rk4 istep = ', istep
   !read(5, *)
   call rk4(rhs, q, cfl, t, dt, use_supplied_dt, unphysical, iphi_bad, &
                  hit_target, t_target, target_met)
   !print *, ' After rk4 istep = ', istep
   !read(5, *)
   ! if (unphysical) call terminate_with_save(1, istep, t)
   ! call print_max(q)
   if (my_node .eq. 0) then
      print *, ' node 0: finished istep = ', istep, ' t = ', t
   end if

   if (mod(istep - istep0, tecplot_interval) .eq. 0) then
      if (my_node .eq. 0) print *, ' about to call tecplot_meridional_plane'
      call tecplot_meridional_plane (q, 1, L_scale, T_scale, t, istep)
      if (nphi .ne. 1) call tecplot_horizontal_plane (q, iz_mid, T_scale, t, istep)      
      ! call phi_Favre_averaged_statistics(istep, t)

      ! This will perform phi and z averages:
      !call phi_Reynolds_averages(istep, t, average_in_z)
   end if

   if (mod(istep - istep0, history_interval) .eq. 0) then
      call output_torque(t, lun_torque)
   end if

   if (mod(istep - istep0, save_interval) .eq. 0) then
      call write_save_file(istep, t, filename_given)
   end if

   ! Detect blowing-up:
   if (istep .gt. istep0+1) then
      if (dt .lt. dt_previous/20.d0) then
         print *, ' my_node = ', my_node
         print *, ' dt decreased more than a factor of 20'
         print *, ' dt          = ', dt
         print *, ' dt_previous = ', dt_previous
         call tecplot_meridional_plane (q, 1, L_scale, T_scale, t, istep)         
         call terminate_with_save(1, istep, t)
      end if
   end if
   dt_previous = dt
end do

100 continue

call output_uphi_versus_r(t, lun_general_purpose, iz_mid)
! Need the minus 1 since a fortran do loop increments the counter at the end of the
! loop.
call terminate_with_save(0, istep-1, t)

end subroutine app_taylor_couette

!----------------------------------------------------------------------------------85

subroutine laminar_profile_output

use partition_data
use grid
use logical_units
use taylor_couette_params
implicit none

! Local:
integer :: ir
real(8) :: uphi

if (my_node .eq. 0) then
   open(unit = lun_general_purpose, file = 'laminar_uphi.dat', form = 'formatted', &
        status = 'unknown')
   
   do ir = 1, nr
      uphi = (A_tc*rgrid(ir) + B_tc/rgrid(ir))
      write(lun_general_purpose, "(2(1x, e12.5))") rgrid(ir), uphi
   end do
   close(unit = lun_general_purpose)
end if

end subroutine laminar_profile_output

!----------------------------------------------------------------------------------85

subroutine add_perturbation_tc

use q_array
use grid
use dof_indices
use partition_data
use taylor_couette_params
#ifdef mpi_code
   use mpi, only: mpi_comm_world
#endif

! ifort compatibility library so we can use the same random number generator as gfortran.   
#ifdef ifort
   use IFPORT 
#endif
   
implicit none
integer, parameter :: nkr = 7, nkz = 9, nkphi_dim = 8
integer :: nkphi
real(8) :: kr(nkr), kz(nkz), kphi(nkphi_dim)
! Local:
integer :: ir, iz, iphi, ikr, ikz, ikphi, ier
real(8) :: pi, eps, lambda_z

real(8) :: lambda_r_modulation, kr_modulation, lambda_z_modulation, kz_modulation, &
     r_modulation, z_modulation, modulation
real(8) :: rho, uz, ur, uphi
real(8) :: uz_prime, ur_prime, uphi_prime
real(8) :: phase_ur, phase_uz, phase_uphi
complex(8) :: A_ur, A_uz, A_uphi, exp_arg, exp_fac, ic = CMPLX(0.0d0, 1.0d0)

if (my_node .eq. 0) print *, ' add_perturbation_tc has been called'

if (nphi .eq. 1) then
   nkphi = 1
else
   nkphi = nkphi_dim
end if

call srand(987654)
pi = 4.0d0 * atan(1.0d0)

! The purpose of the modulation is to make normal velocities = 0 at radial boundaries.
! The modulation is a half sin.
lambda_r_modulation = 2.d0 * (rgrid(nr) - rgrid(1))
kr_modulation       = 2.d0 * pi / lambda_r_modulation

do ikr = 1, nkr
   kr(ikr) = 2.d0 * pi / (rmax - rmin) * ikr
end do

lambda_z = nz*dz(1)
do ikz = 1, nkz
   kz(ikz) = ikz * 2.d0 * pi / lambda_z
end do

do ikphi = 1, nkphi
   kphi(ikphi) = ikphi
end do

#ifdef mpi_code
   call mpi_barrier(mpi_comm_world, ier)
#endif

do iz = 1, nz
   ! print *, ' in point loop iz = ', iz
   ! print *, ' sphi = ', sphi, ' ephi = ', ephi
   do iphi = sphi, ephi
      do ir = sr, er
         ! if (my_node .eq. 0) print *, ' ir = ', ir, ' iphi = ', iphi, ' iz = ', iz      
         r_modulation = sin(kr_modulation*(rgrid(ir) - rgrid(1)))
         ! modulation = r_modulation * z_modulation
         modulation = r_modulation
         
         rho  = q(ir, iphi, iz, irho)
         uz   = q(ir, iphi, iz, zmom) / rho
         ur   = q(ir, iphi, iz, rmom) / rho
         uphi = q(ir, iphi, iz, amom) / rho / rgrid(ir)

         ! Sum over modes:
         uz_prime   = 0.0d0
         ur_prime   = 0.0d0
         uphi_prime = 0.0d0
         ! print *, ' beginning sum over modes'
         do ikr = 1, nkr
            do ikz = 1, nkz
               do ikphi = 1, nkphi
                  phase_ur   = rand(0) * 2.d0 * pi                  
                  phase_uz   = rand(0) * 2.d0 * pi
                  phase_uphi = rand(0) * 2.d0 * pi

                  ! print *, ' phase_uz = ', phase_uz
                  ! print *, ' phase_ur = ', phase_ur
                  ! print *, ' phase_uphi = ', phase_uphi
                  A_uz   = cos(phase_uz)   + ic*sin(phase_uz)
                  A_ur   = cos(phase_ur)   + ic*sin(phase_ur)
                  A_uphi = cos(phase_uphi) + ic*sin(phase_uphi)
                  exp_arg = ic * (kr(ikr)*rgrid(ir) + kz(ikz)*zgrid(iz) + &
                       kphi(ikphi)*phi_grid(iphi))
                  exp_fac = EXP(exp_arg)
                  uz_prime   = uz_prime   + A_uz   * exp_fac
                  ur_prime   = ur_prime   + A_ur   * exp_fac
                  uphi_prime = uphi_prime + A_uphi * exp_fac
                  ! print *, ' phase_uz = ', phase_uz
                  ! print *, ' A_uz = ', A_uz, ' pi = ', pi
               end do
            end do
         end do ! sum over modes
         ! print *, ' finished sum over modes'
         uz_prime   = eps_pert * modulation * uz_prime
         ur_prime   = eps_pert * modulation * ur_prime
         uphi_prime = eps_pert * modulation * uphi_prime

         !print *, ' eps_pert = ', eps_pert
         !print *, ' uz_prime = ', uz_prime
         !print *, ' ur_prime = ', ur_prime
         !print *, ' uphi_prime = ', uphi_prime

         q(ir, iphi, iz, zmom) = rho*(uz   + uz_prime)
         q(ir, iphi, iz, rmom) = rho*(ur   + ur_prime)
         q(ir, iphi, iz, amom) = rho*(uphi + uphi_prime)*rgrid(ir)

         ! print *, ' ir = ', ir, ' iz = ', iz, ' iphi = ', iphi
         ! print *, ' q(zmom) = ', q(ir, iphi, iz, zmom)
         ! print *, ' q(rmom) = ', q(ir, iphi, iz, rmom)
         ! print *, ' q(amom) = ', q(ir, iphi, iz, amom)
         ! read (5, *)
      end do
   end do
end do

end subroutine add_perturbation_tc

!----------------------------------------------------------------------------------85

subroutine output_uphi_versus_r(t, lun, iz)

! Output the profile at the given iz.

use taylor_couette_params
use dof_indices
use grid
use q_array
use transposes_of_q_and_qdot
use partition_data

implicit none
real(8) :: t
integer :: lun, iz

! Local:
integer :: ir
character(40) :: filename
real(8) :: A, B

call transpose_z_to_r(ndof, q, q_r_space)

if ((iz .ge. sz_r) .and. (iz .le. ez_r)) then
   write(filename, "('uphi_iz_', i4.4, '_t_', i6.6, f0.4, '.dat')") iz, int(t), t - int(t)
   open (unit = lun, file = filename, form = 'formatted', status = 'unknown')

   do ir = 1, nr
      write(lun, "(3(1x,e12.5))") rgrid(ir), q(ir, 1, iz, amom)/q(ir, 1, iz, irho)/rgrid(ir), &
           A_tc*rgrid(ir) + B_tc/rgrid(ir)
   end do
   close(lun)
end if

end subroutine output_uphi_versus_r

!----------------------------------------------------------------------------------85

subroutine output_uz_versus_r(t, lun)

use taylor_couette_params
use dof_indices
use grid
use q_array
implicit none
real(8) :: t
integer :: lun

! Local:
integer :: ir
character(40) :: filename
real(8) :: A, B

iz_mid = nz/2

write(filename, "('uz_t_', i6.6, f0.4, '.dat')") int(t), t - int(t)
open (unit = lun, file = filename, form = 'formatted', status = 'unknown')

do ir = 1, nr
   write(lun, "(3(1x,e12.5))") rgrid(ir), q(ir, 1, iz_mid, zmom)/q(ir, 1, iz_mid, irho), &
        C1*log(rgrid(ir)) + C2
end do
close(lun)
end subroutine output_uz_versus_r

!----------------------------------------------------------------------------------85

subroutine output_torque(t, lun_torque)

! Calculates and outputs the torque coeffocient at the inner and outer walls.

! A Fargo correction is not applied since we assume that this calculation is for
! diagnostic purposes and a "shift-back" has been performed."

! See notes of Dec. 4, 2018.

#ifdef mpi_code
   use mpi
#endif
use partition_data
use q_array
use dof_indices
use grid
use viscous, only: nu_molecular
use math_constants
use taylor_couette_params
use logical_units

implicit none
real(8) :: t
integer :: lun_torque  ! The logical unit into which we will output

! Local:
real(8), dimension(sr:er, sphi:ephi, nz) :: uphi_over_r
! These are dimensioned in z-space but will have data in other spaces since pencils in
! all the spaces have the same size:
real(8), dimension(sr:er, sphi:ephi, nz) :: transposed, deriv
real(8), dimension(sr:er, sphi:ephi, nz) :: deriv_in_z_space
real(8) :: torque_outer

integer :: ir, iphi, iz, ier
real(8) :: normalization, mu_molecular, my_stress_sum, stress_sum, torque, &
     inner_torque_coeff, outer_torque_coeff, average_stress

! Dong's normalization:
if ((itype .eq. dong) .or. (itype .eq. dong_counter_rotating)) then
   normalization = 0.5d0 * rho0 * Ui**2 * pi * r_i**2 * (zmax - zmin)
else if (itype .eq. meyer) then
   normalization = 1.d0 / (omega_i * r_i**2)
else if (itype .eq. moser) then
   !if (isubtype .eq. moser_meyer) then
   !   normalization = (zmax - zmin)
   !else
      normalization = (zmax - zmin) * rho0 * nu_molecular**2
   !end if
else if (itype .eq. marcus) then  
   normalization = (zmax - zmin) ! To get a torque per unit length.
end if

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         uphi_over_r(ir,iphi,iz) = q(ir,iphi,iz,amom) / q(ir,iphi,iz,irho) / rgrid(ir)**2         
      end do
   end do
end do

call transpose_z_to_r(1, uphi_over_r, transposed)
call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
call transpose_r_to_z(1, deriv, deriv_in_z_space)

! -----------
! Inner wall:
! -----------
if (sr .eq. 1) then ! If I have the inner wall then
   ir = 1
   my_stress_sum = 0.0d0
   do iz = 1, nz
      do iphi = sphi, ephi
         mu_molecular = nu_molecular * q(ir,iphi,iz,irho)
         ! Note the factor of r since we differenitated u_phi/r
         my_stress_sum = my_stress_sum + mu_molecular*rgrid(ir)*deriv_in_z_space(ir,iphi,iz)
      end do
   end do
else
   ! I don't have the inner wall:
   my_stress_sum = 0.d0
end if

#ifdef mpi_code
call mpi_reduce(my_stress_sum, stress_sum, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ier)
#else
stress_sum = my_stress_sum
#endif
! Node 0 now has stress_sum

if (my_node .eq. 0) then
   ! torque = average stress * area * radius.
   ! Torque of the wall on the fluid:
   torque = -stress_sum/(nphi*nz) * 2.d0*pi*r_i*(zmax - zmin) * r_i
   inner_torque_coeff = torque / normalization
end if

! -----------
! Outer wall:
! -----------
if (er .eq. nr) then ! If I have the outer wall then
   ir = nr
   my_stress_sum = 0.0d0
   do iz = 1, nz
      do iphi = sphi, ephi
         mu_molecular = nu_molecular * q(ir,iphi,iz,irho)
         ! Note the factor of r:
         my_stress_sum = my_stress_sum + mu_molecular*rgrid(ir)*deriv_in_z_space(ir,iphi,iz)
      end do
   end do
else
   ! I don't have the outer wall:
   my_stress_sum = 0.d0
end if

#ifdef mpi_code
! Reduced result sits in rank 0.
call mpi_reduce(my_stress_sum, stress_sum, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ier)
#else
stress_sum = my_stress_sum
#endif
if (my_node .eq. 0) then
   torque = stress_sum/(nphi*nz) * 2.d0*pi*r_o*(zmax - zmin) * r_o
   outer_torque_coeff = torque / normalization
end if

if (my_node .eq. 0) then
   open(unit = lun_torque, file = 'torque.his',    form = 'formatted', &
        status = 'unknown', access = 'append')
   write(lun_torque, "(3(1x, e16.9))") t, inner_torque_coeff, -outer_torque_coeff
   ! So we don't lose the history:
   close(unit = lun_torque)
end if

end subroutine output_torque

!----------------------------------------------------------------------------------85

subroutine check_viscous

! Checks the calculation of the viscous term for an analytical velocity field.

use dof_indices
use partition_data
use grid
use q_array
use viscous
use math_constants
implicit none

! Local:
real(8) :: kr, kz, lambda_max, rho0, ur, uz, up, laplacian1, laplacian2, laplacian, &
     Fr_exact, krr, kzz, kz_fundamental, Fz_exact, Fp_exact
real(8) :: sin_r, cos_r, sin_z, cos_z, bracket, r_term, z_term, r, z
integer :: iz, iphi, ir

kz_fundamental = 2.d0 * pi / (zmax - zmin)
kr   = 4.d0
kz   = 2.d0*kz_fundamental
rho0 = 1.0d0
! Check viscous Fr for a test velocity field:
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
        krr = kr*rgrid(ir)
        kzz = kz*zgrid(iz)
        ur = -kz * sin(krr) * cos(kzz)
        uz = sin(kzz)/rgrid(ir) * (krr*cos(krr) + sin(krr))
        up = sin(krr)*sin(kzz);
        q(ir, iphi, iz, irho) = rho0
        q(ir, iphi, iz, zmom) = rho0 * uz
        q(ir, iphi, iz, rmom) = rho0 * ur
        q(ir, iphi, iz, amom) = rho0 * up * rgrid(ir)
      end do
   end do
end do

call strain_tensor(q)
call molecular_viscous_stress_and_heat_term(q, lambda_max)
call viscous_force_and_heat_term

open(unit = 1, file = 'viscous_forces.dat',    form = 'formatted', status = 'unknown')
iz   = 10
iphi = 1

do ir = 1, nr
   r = rgrid(ir)
   z = zgrid(iz)

   Fr_exact = (kz*sin(kr*r)*cos(kz*z))/r**2+kz**3*sin(kr*r)*cos(kz*z)+ &
        kr**2*kz*sin(kr*r)*cos(kz*z)-(kr*kz*cos(kr*r)*cos(kz*z))/r
   Fr_exact = Fr_exact * nu_molecular

   Fz_exact = (-(kz**2*sin(kr*r)*sin(kz*z))/r)-(2*kr**2*sin(kr*r)*sin(kz*z))/r+( &
        sin(kr*r)*sin(kz*z))/r**3-(kr*cos(kr*r)*sin(kz*z))/r**2-kr*kz** &
        2*cos(kr*r)*sin(kz*z)-kr**3*cos(kr*r)*sin(kz*z)
   Fz_exact = Fz_exact * nu_molecular

   Fp_exact = &
      (-(sin(kr*r)*sin(kz*z))/r**2)-kz**2*sin(kr*r)*sin(kz*z)-kr**2*sin( &
      kr*r)*sin(kz*z)+(kr*cos(kr*r)*sin(kz*z))/r
   Fp_exact = Fp_exact * nu_molecular
      
   write(1, "(7(1x, e12.5))") rgrid(ir), Fr(ir, iphi, iz), Fr_exact, Fz(ir, iphi, iz), Fz_exact, &
        Fp(ir, iphi, iz), Fp_exact 
end do

end subroutine check_viscous

!----------------------------------------------------------------------------------85

subroutine check_viscous_with_phi

use dof_indices
use partition_data
use grid
use q_array
use viscous
use math_constants
implicit none

! Local:
real(8) :: kr, kz, kp, lambda_max, rho0, ur, uz, up, laplacian1, laplacian2, laplacian, &
     Fr_exact, Fz_exact, Fp_exact, kr_fundamental, kz_fundamental
real(8) :: r, z, p
integer :: iz, iphi, ir

kr_fundamental = 2.d0 * pi / (rmax - rmin)
kz_fundamental = 2.d0 * pi / (zmax - zmin)

kr   = kr_fundamental
kz   = kz_fundamental
kp   = 2.d0
rho0 = 1.0d0
nu_molecular = 1.d0
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
        r = rgrid(ir)
        z = zgrid(iz)
        p = phi_grid(iphi)

        ur = (2*kp*sin(kp*p)*cos(kr*r)*sin(kz*z))/r-3*kz*cos(kp*p)*sin(kr*r)*sin(kz*z)
        uz = (kp*cos(kp*p)*sin(kr*r)*cos(kz*z))/r-(3*cos(kp*p)*sin(kr*r)*cos(kz*z))/r-&
                3*kr*cos(kp*p)*cos(kr*r)*cos(kz*z)
        up = kz*sin(kp*p)*sin(kr*r)*sin(kz*z)-2*kr*cos(kp*p)*sin(kr*r)*sin(kz*z)

        q(ir, iphi, iz, irho) = rho0
        q(ir, iphi, iz, zmom) = rho0 * uz
        q(ir, iphi, iz, rmom) = rho0 * ur
        q(ir, iphi, iz, amom) = rho0 * up * rgrid(ir)
      end do
   end do
end do

call strain_tensor(q)
call molecular_viscous_stress_and_heat_term(q, lambda_max)
call viscous_force_and_heat_term

open(unit = 1, file = 'viscous_forces.dat',    form = 'formatted', status = 'unknown')
iz   = 10
iphi = 7

print *, ' nphi = ', nphi

do ir = 1, nr
   r = rgrid(ir)
   z = zgrid(iz)
   p = phi_grid(iphi)

   print *, ' r, z, p, kr, kz, kp = ', r, z, p, kr, kz, kp
   Fr_exact = (r*((4*kp*kr*sin(kp*p)*sin(kr*r)*sin(kz*z))/r**2+ &
        3*kr**2*kz*cos(kp*p)*sin(kr*r)*sin(kz*z)- &
        (2*kp*kr**2*sin(kp*p)*cos(kr*r)*sin(kz*z))/r+ &
        (4*kp*sin(kp*p)*cos(kr*r)*sin(kz*z))/r**3)-&
        (2*kp*kr*sin(kp*p)*sin(kr*r)*sin(kz*z))/r-&
        (2*kp*sin(kp*p)*cos(kr*r)*sin(kz*z))/r**2-&
        3*kr*kz*cos(kp*p)*cos(kr*r)*sin(kz*z))/r-&
        (2*(2*kp*kr*sin(kp*p)*sin(kr*r)*sin(kz*z)+&
        kp*kz*cos(kp*p)*sin(kr*r)*sin(kz*z)))/r**2+&
        (3*kp**2*kz*cos(kp*p)*sin(kr*r)*sin(kz*z)-&
        (2*kp**3*sin(kp*p)*cos(kr*r)*sin(kz*z))/r)/r**2-&
        ((2*kp*sin(kp*p)*cos(kr*r)*sin(kz*z))/r-&
        3*kz*cos(kp*p)*sin(kr*r)*sin(kz*z))/r**2+&
        3*kz**3*cos(kp*p)*sin(kr*r)*sin(kz*z)-&
        (2*kp*kz**2*sin(kp*p)*cos(kr*r)*sin(kz*z))/r

   Fz_exact = (r*((-(kp*kr**2*cos(kp*p)*sin(kr*r)*cos(kz*z))/r)+&
          (3*kr**2*cos(kp*p)*sin(kr*r)*cos(kz*z))/r+&
          (2*kp*cos(kp*p)*sin(kr*r)*cos(kz*z))/r**3-&
          (6*cos(kp*p)*sin(kr*r)*cos(kz*z))/r**3-&
          (2*kp*kr*cos(kp*p)*cos(kr*r)*cos(kz*z))/r**2+&
          (6*kr*cos(kp*p)*cos(kr*r)*cos(kz*z))/r**2+&
          3*kr**3*cos(kp*p)*cos(kr*r)*cos(kz*z))-&
          (kp*cos(kp*p)*sin(kr*r)*cos(kz*z))/r**2+&
          (3*cos(kp*p)*sin(kr*r)*cos(kz*z))/r**2+&
          3*kr**2*cos(kp*p)*sin(kr*r)*cos(kz*z)+&
          (kp*kr*cos(kp*p)*cos(kr*r)*cos(kz*z))/r-&
          (3*kr*cos(kp*p)*cos(kr*r)*cos(kz*z))/r)/r+&
          ((-(kp**3*cos(kp*p)*sin(kr*r)*cos(kz*z))/r)+&
          (3*kp**2*cos(kp*p)*sin(kr*r)*cos(kz*z))/r+&
          3*kp**2*kr*cos(kp*p)*cos(kr*r)*cos(kz*z))/r**2-&
          (kp*kz**2*cos(kp*p)*sin(kr*r)*cos(kz*z))/r+&
          (3*kz**2*cos(kp*p)*sin(kr*r)*cos(kz*z))/r+&
          3*kr*kz**2*cos(kp*p)*cos(kr*r)*cos(kz*z)

   Fp_exact = (r*(2*kr**3*cos(kp*p)*sin(kr*r)*sin(kz*z)-&
        kr**2*kz*sin(kp*p)*sin(kr*r)*sin(kz*z))+&
        kr*kz*sin(kp*p)*cos(kr*r)*sin(kz*z)-&
        2*kr**2*cos(kp*p)*cos(kr*r)*sin(kz*z))/r+&
        (2*kp**2*kr*cos(kp*p)*sin(kr*r)*sin(kz*z)-&
        kp**2*kz*sin(kp*p)*sin(kr*r)*sin(kz*z))/r**2+&
        (2*(3*kp*kz*sin(kp*p)*sin(kr*r)*sin(kz*z)+&
        (2*kp**2*cos(kp*p)*cos(kr*r)*sin(kz*z))/r))/r**2+&
        (2*kr*cos(kp*p)*sin(kr*r)*sin(kz*z)-&
        kz*sin(kp*p)*sin(kr*r)*sin(kz*z))/r**2-&
        kz**3*sin(kp*p)*sin(kr*r)*sin(kz*z)+&
        2*kr*kz**2*cos(kp*p)*sin(kr*r)*sin(kz*z)

   write(1, "(7(1x, e12.5))") rgrid(ir), Fr(ir, iphi, iz), Fr_exact, Fz(ir, iphi, iz), Fz_exact, &
        Fp(ir, iphi, iz), Fp_exact 
end do

end subroutine check_viscous_with_phi

!----------------------------------------------------------------------------------85

!!$subroutine tc_laminar_profiles(r, rho, p, c_sound)
!!$
!!$! Unused.
!!$! Density field for homentropic Taylor-Couette flow.
!!$
!!$use taylor_couette_params
!!$implicit none
!!$real(8) :: r, rho, p, c_sound
!!$
!!$! Local:
!!$real(8) :: bracket, stuff, ratio, gm1, K
!!$
!!$gm1 = gamma - 1
!!$
!!$bracket = 0.5d0*A_tc**2*(r**2 - r_i**2) - 0.5d0*B_tc**2*(r**(-2) - r_i**(-2)) + &
!!$           2.d0*A_tc*B_tc*log(r/r_i)
!!$stuff = gm1/c_sound_i**2*bracket + 1.d0
!!$ratio = stuff**(1.d0/gm1)
!!$rho   = ratio*rho_i
!!$
!!$K = c_sound_i**2 / gamma / rho_i**gm1
!!$p = K * rho**gamma
!!$c_sound = sqrt(gamma*p/rho)
!!$
!!$end subroutine tc_laminar_profiles

!----------------------------------------------------------------------------------85
