!----------------------------------------------------------------------------------85

subroutine pade_diff_test

! Tests Pade differentiation routines.  This subroutine is strictly serial.

use grid
use math_constants

#ifdef mpi_code
   use mpi
#endif
implicit none

! Functions of r:
real(8), allocatable, dimension(:,:,:) :: f 

logical, parameter :: stretched_r = .false., stretched_z = .false.
real(8), parameter :: r0_unif = 0.d0, z0 = 0.d0
integer, parameter :: nr_u = 0, nz_u = 0

nr   = 64
nz   = 64
nphi = 64

rmin = 1.d0
rmax = 2.d0
zmin = -1.d0
zmax =  1.d0

phi_min = 0.d0
phi_max = 2.d0*pi

if (my_node .eq. 0) then
   print *, ' node 0: First executable of pade_diff_test'
   print *, ' my_node = ', my_node, ' num_nodes = ', num_nodes   
end if

allocate(f(sr:er, sphi:ephi, nz)
call set_up_domain_mesh_and_partition(rmin, rmax, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
     suppress_z_derivatives_when_nz_not_1_arg, periodic_z_arg, &
     stretched_r, stretched_z, r0_unif, nr_u, z0, nz_u)               

dz = 

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         f(ir, iphi, iz) = rho_temporary * ci(ir)**2
      end do
   end do
end do

! Initial condition:
if (.not. restart) then
   call pade_diff_z(mr*mphi, pressure, dpdz)

   ! Set rho to satisfy numerical vertical hydrostatic balance:
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            ! Required since gz = 0 at the midplane:
            if (gz(ir, iz) .ne. 0.d0) then
               q(ir, iphi, iz, irho) = dpdz(ir, iphi, iz)/gz(ir, iz)
            else
               q(ir, iphi, iz, irho) = rho_midplane(ir)
            end if
            ! Correct the pressure:
            pressure(ir, iphi, iz) = q(ir,iphi,iz,irho) * ci(ir)**2
         end do
      end do
   end do

   ! Set uphi to satisfy numerical centrifugal balance:
   call transpose_z_to_r (1, pressure, p_r_space)
   call pade_diff_bundle(mphi*mz_r, nr, Ji_r, p_r_space, dpdr_r_space)
   call transpose_r_to_z (1, dpdr_r_space, dpdr)   
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, iz, rmom) = 0.d0
            q(ir, iphi, iz, zmom) = 0.0d0

            ! term1 = (p + qexp) * (H(ir)/rgrid(ir))**2
            ! term2 = 1.d0 + qexp
            ! R_spherical = SQRT(rgrid(ir)**2 + zgrid(iz)**2)
            ! term3 = - qexp*rgrid(ir) / R_spherical
            ! Omega = Omega_K(ir) * (term1 + term2 + term3)**(0.5d0)
            ! uphi  = Omega * rgrid(ir)

            ! Look at NRR-3 notes:
            uphi_squared = rgrid(ir)*(1.d0/q(ir,iphi,iz,irho)*dpdr(ir,iphi,iz) - gr(ir,iz))
            q(ir, iphi, iz, amom) = q(ir, iphi, iz, irho) * sqrt(uphi_squared) * rgrid(ir)

            ! Set energy in case we run adiabatically:
            ! q(ir, iphi, iz, ener) = pressure(ir,iphi,iz)/gm1 + 0.5d0 * q(ir, iphi, iz, irho) * uphi**2
            q(ir, iphi, iz, ener) = 0.d0
         end do
      end do
   end do
end if
! Finished with basic state

deallocate(pressure, dpdz, p_r_space, dpdr_r_space, dpdr)

! Store what is in q as the basic state:
call store_basic_state(filter_relative_to_basic_state)

if (restart) then
   call read_restart (file_version, istep_of_restart_file, t_of_restart_file)
else
   ! Assign initial condition using the basic state:
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, iz, irho) = rho_basic(ir, iz)
            q(ir, iphi, iz, rmom) = 0.d0
            q(ir, iphi, iz, zmom) = 0.d0                        
            q(ir, iphi, iz, amom) = rho_basic(ir, iz) * uphi_basic(ir, iz) * rgrid(ir)

            ! Set energy in case we run adiabatically:
            ! q(ir, iphi, iz, ener) = pressure(ir,iphi,iz)/gm1
            q(ir, iphi, iz, ener) = 0.d0
         end do
      end do
   end do
end if
! Finished with initial condition.

if (apply_sponge) then
   call activate_sponge(sponge_type, rho1, rho2, rho_basic, d1, d2, tau_decay)   
end if

! Needed for the isothermal option.  ci_squared_initial sits in
! module thermal_parameters.
if (isothermal) then
   do iz = 1, nz
      do ir = 1, nr
         ci_squared_initial(ir, iz) = ci(ir)**2
      end do
   end do
end if

if (.not. restart) then
   istep0 = 0
   t      = 0.0d0
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
   print *, ' calling tecplot_fluc_vel_meridional_plane'
end if

iphi_plot    = 1
iz_plot      = iz_mid
!iz_plot      = 192
L_scale      = 1.0d0
T_scale      = 1.0d0
rho_scale    = 1.0d0
plot_curl_rho_u = .true.
call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot, &
     L_scale, T_scale, rho_scale, t, istep0)
plot_curl_rho_u = .false.
call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot, &
     L_scale, T_scale, rho_scale, t, istep0)
call tecplot_fluc_vel_meridional_plane         (q, iphi_plot,          L_scale, T_scale, t, istep0)
call tecplot_meridional_plane                  (q, iphi_plot,          L_scale, T_scale, t, istep0)
if (nphi .ne. 1) call tecplot_horizontal_plane (q,            iz_plot,          T_scale, t, istep0)
call phi_Favre_averaged_statistics(istep0, t)

if (output_phi_Reynolds_averages) then
   call phi_Reynolds_averages(istep0, t, average_in_z)
   next_t_for_phi_Reynolds_averages = t + time_interval_for_phi_Reynolds_averages
end if

call write_ave_fluctuation_ke(t)

! Time stepping loop:
do istep = istep0 + 1, istep0 + nsteps
   call rk4(rhs, q, cfl, t, dt, use_supplied_dt, unphysical_flag, iphi_bad, &
        hit_target, t_target, target_met)
   call output_conservation_diagnostics(t, q)
   
   if (my_node .eq. 0) then
      print *, ' routine app_vsi_3D, node 0: finished istep = ', istep, ' t = ', t
   end if

   !call print_max(q)

   if (mod(istep - istep0, tecplot_interval) .eq. 0) then
      call tecplot_meridional_plane (q, 1, 1.0d0, 1.0d0, t, istep)
      if (nphi .ne. 1) call tecplot_horizontal_plane (q, iz_plot, 1.0d0, t, istep)

      plot_curl_rho_u = .true.      
      call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot, &
           L_scale, T_scale, rho_scale, t, istep)
      plot_curl_rho_u = .false.            
      call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot, &
           L_scale, T_scale, rho_scale, t, istep)      
      call tecplot_fluc_vel_meridional_plane(q, 1, 1.0d0, 1.0d0, t, istep)
      
      ! Debug:
      ! Plot the artificial pressure and dilatation in the plane in which the eigenvalue lambda_max_ap
      ! is the largest:
      !if (apply_artificial_pressure) then
      !   get_lambda_max_ap = .true.
      !   if (my_node .eq. 0) print *, ' calling get_artificial_pressure' 
      !   call get_artificial_pressure(get_lambda_max_ap, q, lambda_max_ap)
      !   if (my_node .eq. 0) print *, ' calling tecplot_p_art_in_meridional_plane'
      !   call tecplot_scalar_in_meridional_plane('p_art', p_art, iphi_max_ap, 1.0d0, 1.0d0, &
      !        1.0d0, t, istep)
      !   call tecplot_scalar_in_meridional_plane('dilat', dil,   iphi_max_ap, 1.0d0, 1.0d0, &
      !                                            1.0d0, t, istep)         
      !end if
   end if

   if ((output_phi_Reynolds_averages) .and. (t .ge. next_t_for_phi_Reynolds_averages)) then
      call phi_Reynolds_averages(istep, t, average_in_z)
      next_t_for_phi_Reynolds_averages = t + time_interval_for_phi_Reynolds_averages
   end if   

   if (mod(istep - istep0, fluctuation_ke_interval) .eq. 0) then
      call write_ave_fluctuation_ke(t)
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
      call tecplot_meridional_plane (q, 1, 1.0d0, 1.0d0, t, istep)
      if (nphi .ne. 1) call tecplot_horizontal_plane (q, iz_plot, 1.0d0, t, istep)

      call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot, &
           L_scale, T_scale, rho_scale, t, istep)      
      call tecplot_fluc_vel_meridional_plane(q, 1, 1.0d0, 1.0d0, t, istep)
      close(lun_history(1))
      ! 1 = istatus : Abnormal return.
      call terminate_with_save(1, istep, t)
   end if
   
   dt_previous = dt
end do

if (my_node .eq. 0) then
   print *, ' ***** average stepping time per step = ', total_cpu_for_stepping/nsteps
end if

close(lun_history(1))
! 0 = istatus: Normal return
! Need the minus 1 since a fortran do loop increments the counter at the end
call terminate_with_save(0, istep-1, t)

end subroutine vsi_3D

!----------------------------------------------------------------------------------85

subroutine add_wavy_perturbation(c0, H0_over_r0, H0, lambda_r)

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
real(8), intent(in) :: c0, H0_over_r0, H0, lambda_r
integer, parameter :: nkr = 11, nkz = 12, nkphi_dim = 12
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
complex(8) :: exp_arg, exp_fac, ic = CMPLX(0.0d0, 1.0d0)

complex(8), allocatable, dimension(:,:,:) :: A_ur, A_uz, A_uphi

if (my_node .eq. 0) print *, ' modified add_perturbation_3D has been called'

if (nphi .eq. 1) then
   nkphi = 1
else
   nkphi = nkphi_dim
end if

allocate(A_ur(nkphi,nkz,nkr), A_uz(nkphi,nkz,nkr), A_uphi(nkphi,nkz,nkr))

! Establish seed for random number generator different for each processor:
call srand(987654 + my_node*769)
pi = 4.0d0 * atan(1.0d0)
eps = 1.d-3*c0

! The purpose of the modulation is to make normal velocities = 0 at boundaries.
! The modulation is a half sin.
lambda_r_modulation = 2.d0 * (rgrid(nr) - rgrid(1))
kr_modulation       = 2.d0 * pi / lambda_r_modulation
lambda_z_modulation = 2.d0 * (zgrid(nz) - zgrid(1))
kz_modulation       = 2.d0 * pi / lambda_z_modulation

! Note: nkr = 11
kr(1) = 2.d0 * pi / lambda_r
do ikr = 2, 8
   kr(ikr) = ikr * kr(1)
end do
! Sub-harmonics:
kr(9)  = kr(1)/2.d0
kr(10) = kr(1)/3.d0
kr(11) = kr(1)/4.d0

lambda_z = zgrid(nz) - zgrid(1)
kz(1) = 2.d0 * pi / lambda_z
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

! This entire section in incorrect for what you want to do.


! Get complex amplitudes for each mode:
do ikr = 1, nkr
   do ikz = 1, nkz
      do ikphi = 1, nkphi
         ! Zero argument means next random number in the sequence.
         phase_ur   = rand(0) * 2.d0 * pi                  
         phase_uz   = rand(0) * 2.d0 * pi
         phase_uphi = rand(0) * 2.d0 * pi

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
         uz_prime   = 0.0d0
         ur_prime   = 0.0d0
         uphi_prime = 0.0d0
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
real(8), intent(in) :: c0
integer :: ir, iz, iphi, ier
real(8) :: eps

real(8) :: lambda_r_modulation, kr_modulation, lambda_z_modulation, kz_modulation, &
     r_modulation, z_modulation, modulation
real(8) :: rho, uz, ur, uphi
real(8) :: uz_prime, ur_prime, uphi_prime

if (my_node .eq. 0) print *, ' random_perturbation has been called'

! Establish seed for random number generator:
call srand(987654 + 789*my_node)
eps = 1.d-3*c0

! The purpose of the modulation is to make normal velocities = 0 at boundaries.
! The modulation is a half sin.
lambda_r_modulation = 2.d0 * (rgrid(nr) - rgrid(1))
kr_modulation       = 2.d0 * pi / lambda_r_modulation
lambda_z_modulation = 2.d0 * (zgrid(nz) - zgrid(1))
kz_modulation       = 2.d0 * pi / lambda_z_modulation

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
         modulation = 1.d0
         
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





