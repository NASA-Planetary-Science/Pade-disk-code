!----------------------------------------------------------------------------------85

module activate_routines

!> This module contains routines that the user may call in the set-up portion of
!> his/her application to activate various features.  <br>
!> These are <br>
!>    FARGO for time step reduction <br>
!>    Artificial pressure (also known as artificial bulk viscosity) for
!>    shock-capturing <br>
!>    Pade filtering to filter 2-delta waves (or more) <br>
!>    Plotting shift in which the output will follow a feature rotating around
!>    the azimuth at Keplerian velocity.

contains

! Activation routines for:
! gravity, artificial_pressure, pade_filter, fargo, plotting_shift

!----------------------------------------------------------------------------------85

subroutine activate_fargo(integer_shifts_arg, apply_fargo_extra_operator_arg)

!> Issue a call to this optional routines in the set-up portion of your user app to
!> activate the FARGO method for reducing the time step if it is dominated by
!> Keplerian advection due to a central mass.

use grid
use partition_data
use fargo_or_plotting_shift
use total_allocated_words, only: n_words_allocated
use transposes_of_q_and_qdot
use math_constants
use gravity, only: gravity_flag, GM
implicit none
!> Set this to .true. if you want the (discontinuous) integer shifting of the
!> original Masset FARGO method.  Set this to .false. if you want real-valued
!> continuous shifting of the corrected scheme.
logical, intent(in) :: integer_shifts_arg
!> Set this to .true. if you want the extra chain rule terms to be included.
!> Set this to .false. if you don't want those terms.  You can use this
!> option to save time during approach to a statistically stationary state.
logical, intent(in) :: apply_fargo_extra_operator_arg

! We should not apply Fargo if nphi = 1:
if (nphi .eq. 1) then
   if (my_node .eq. 0) then
      print *, ' Message from subroutine activate_fargo:'
      print *, ' Cannot apply the Fargo trick when nphi = 1'
      print *, ' Returning without activating Fargo'
   end if
   return
end if

apply_fargo                = .true.
integer_shifts             = integer_shifts_arg   
apply_fargo_extra_operator = apply_fargo_extra_operator_arg   

allocate(Omega_fargo(nr), uphi_fargo_subtract(nr), fargo_factor(nr), fargo_factor_over_r(nr))
allocate(nshift(nr), inew_phi(nr, nphi))
n_words_allocated = n_words_allocated + 3*nr + nr/2 + nr*nphi/2

! The second test is there because set_up_plotting_shift could have allocated it.
if ( (.not. integer_shifts) .and. (.not. allocated(complex_shift_factor)) ) then
   allocate(complex_shift_factor(nphi/2+1, sr:er))
   n_words_allocated = n_words_allocated + 2*(nphi/2+1)*mr
end if

if (.not. have_Omega_bar) call setup_Omega_bar_for_fargo(gravity_flag, GM)

if (.not. fft_has_been_initialized) then
#ifdef fftw
   call make_plans_for_fftw (q_phi_space)
   fft_has_been_initialized = .true.
#else
   call initialize_rogallo_fft_for_fargo
   fft_has_been_initialized = .true.
#endif
end if

twopi_over_Delta_phi_domain = 2.d0 * pi / Delta_phi_domain

end subroutine activate_fargo

!----------------------------------------------------------------------------------85

subroutine activate_plotting_shift(Omega0_plotting_shift_arg, t_of_previous_shift_arg)

!> A call to this routine is optional.
!> Set-up for Fargo trick for phi advection.

use grid
use partition_data
use fargo_or_plotting_shift
use total_allocated_words, only: n_words_allocated
use transposes_of_q_and_qdot
use math_constants
implicit none
logical :: apply_plotting_shift_arg
real(8) :: Omega0_plotting_shift_arg, t_of_previous_shift_arg

apply_plotting_shift  = .true.
Omega0_plotting_shift = Omega0_plotting_shift_arg
t_of_previous_shift   = t_of_previous_shift_arg

if (.not. allocated(complex_shift_factor) ) then
   allocate(complex_shift_factor(nphi/2+1, sr:er))
   n_words_allocated = n_words_allocated + 2*(nphi/2+1)*mr
end if

if (.not. fft_has_been_initialized) then
#ifdef fftw
   call make_plans_for_fftw (q_phi_space)
   fft_has_been_initialized = .true.
#else
   call initialize_rogallo_fft_for_fargo
   fft_has_been_initialized = .true.
#endif
end if

twopi_over_Delta_phi_domain = 2.d0 * pi / Delta_phi_domain

!print *, ' In set_up_plotting_shift'
!print *, ' twopi_over_Delta_phi_domain = ', twopi_over_Delta_phi_domain
!print *, ' enter anything to continue'
!read (5, *)

end subroutine activate_plotting_shift

!----------------------------------------------------------------------------------85

subroutine activate_gravity(gravity_flag_in, gravity_type_in, GM_in, gz_uniform_in)

!> This is an optional routine.  The defult is no gravity.
!> This routine is called to activate or deactivate gravity via
!> the logical variable gravity_flag_in.
!>
!> Arguments:
!> gravity_flag_in : 
!>
!> gravity_type_in (the integer value for these types are provided in
!> module gravity_types).
!>    mass_at_origin : g : - GM / R^2 Rhat where R is the spherical radius and
!>                     Rhat is the unit vector along the spherical radius.
!>    thin_disk_mass_at_origin : g = -GM/r^2 * (rhat + (z/r) zhat) where rhat and zhat are
!>                               unit vectors.
!>    thin_disk_mass_at_origin_no_radial : Similar above where the radial component is
!>                                         set to zero.
!>    uniform_gz : A uniform gravity in the z direction.  Make it < 0 if you want it to
!>                 point in the minus z direction.
!>
!> GM_in : Value of GM for mass_at_origin or thin_disk_mass_at_origin
!>
!> gz_uniform_in : Value of gz_uniform for the uniform_gz option.

use grid
use gravity_types
use gravity
use total_allocated_words, only: n_words_allocated
use partition_data
use math_constants
use fargo_or_plotting_shift
implicit none
! Arguments:
logical, intent(in) :: gravity_flag_in !! Set it to .true. to activate gravity.
integer, intent(in) :: gravity_type_in !! Type of gravity selected from module gravity_types
real(8), intent(in) :: GM_in           !! Value of GM for the first three gravity types
real(8), intent(in) :: gz_uniform_in   !! Value gz for the uniform_gz type. < 0 if gravity is in negative z

! Local:
integer iz, ir
real(8) :: R_spherical,g_mag, cos_theta, sin_theta

#ifdef debug_print
   if (my_node .eq. 0) print *, ' activate_gravity has been called'
#endif
   
gravity_flag = gravity_flag_in
gravity_type = gravity_type_in   
GM           = GM_in
gz_uniform   = gz_uniform_in

allocate (gr(nr, nz), gz(nr, nz))
n_words_allocated = n_words_allocated + 2*nr*nz

! We do this because we would like to suggest Fargo to the user
! if we determine that it might be advantageous.
if(.not. have_Omega_bar) then
   ! Note: If gravity_flag = .false. then this routine sets Omega_bar to zero so it can
   ! be used by the time-step determination routine.
   call setup_Omega_bar_for_fargo(gravity_flag, GM)
end if

if (.not. gravity_flag) then
   ! Zero gravity:
   do iz = 1, nz
      do ir = 1, nr
         gr(ir, iz) = 0.0d0
         gz(ir, iz) = 0.0d0
      end do
   end do
else if (gravity_type .eq. mass_at_origin) then
   do iz = 1, nz
      do ir = 1, nr
         R_spherical = SQRT(rgrid(ir)**2 + zgrid(iz)**2)
         g_mag       = GM / R_spherical**2
         cos_theta   = rgrid(ir) / R_spherical
         sin_theta   = zgrid(iz) / R_spherical
         gr(ir, iz)  = - g_mag * cos_theta
         gz(ir, iz)  = - g_mag * sin_theta
      end do
   end do  
else if (gravity_type .eq. thin_disk_mass_at_origin) then
   do iz = 1, nz
      do ir = 1, nr
         gr(ir, iz) = -GM/rgrid(ir)**2
         gz(ir, iz) = gr(ir, iz) * (zgrid(iz) / rgrid(ir))
      end do
   end do
else if (gravity_type .eq. thin_disk_mass_at_origin_no_radial) then
   do iz = 1, nz
      do ir = 1, nr
         gr(ir, iz) = -GM/rgrid(ir)**2         
         gz(ir, iz) = gr(ir, iz) * (zgrid(iz) / rgrid(ir))
         ! Zero out the radial component:
         gr(ir, iz) = 0.d0         
      end do
   end do
else if (gravity_type .eq. uniform_gz) then
   do iz = 1, nz
      do ir = 1, nr
         gr(ir, iz) = 0.d0
         gz(ir, iz) = - abs(gz_uniform)
      end do
   end do
else
   if (my_node .eq. 0) then
      print *, ' my_node = 0.  subroutine activate_gravity'
      print *, ' Unrecognized gravity_type = ', gravity_type
   end if
   call terminate_with_no_save(1)
end if

! Keplerian velocity in case needed for Fargo option:
if (gravity_flag) then
   if ((gravity_type .eq. mass_at_origin) .or. (gravity_type .eq. thin_disk_mass_at_origin)) then
      allocate (u_Kepler(nr), Omega_Kepler(nr))
      n_words_allocated = n_words_allocated + 2*nr
      do ir = 1, nr
         Omega_Kepler(ir) = SQRT((GM / rgrid(ir)**3))
         u_Kepler    (ir) = Omega_Kepler(ir) * rgrid(ir)
      end do
   end if
end if

end subroutine activate_gravity

!----------------------------------------------------------------------------------85

subroutine activate_artificial_pressure (C_ap_arg)

use partition_data
use artificial_pressure_module
use grid
use total_allocated_words, only: n_words_allocated

real(8), intent(in) :: C_ap_arg  ! constant for artificial pressure term

! Locals:
integer :: ir, iphi, iz 

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: activate_artificial_pressure has been called'
#endif

apply_artificial_pressure = .true.
C_ap = C_ap_arg

if (.not. allocated(dil)) then
   allocate(dil(sr:er, sphi:ephi, nz))
   n_words_allocated = n_words_allocated + mr*mphi*nz
end if

if (.not. partitioned) then
   print *, ' You need to partition the grid before calling activate_artificial_pressure'
   call terminate_with_no_save(1)
end if

if (.not. gridded) then
   print *, ' You need to generate the grid before calling activate_artificial_pressure'
   call terminate_with_no_save(1)
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: About to return from activate artificial pressure'
#endif
   
end subroutine activate_artificial_pressure

!----------------------------------------------------------------------------------85

subroutine activate_pade_filter(apply_pade_filter_arg, eps_or_tau_arg, &
     eps_filter_arg, tau_filter_arg, filter_relative_to_basic_state_arg)

use partition_data
use pade_filter
use logical_units
implicit none
logical,      intent(in) :: apply_pade_filter_arg
character(3), intent(in) :: eps_or_tau_arg
real(8),      intent(in) :: eps_filter_arg, tau_filter_arg
logical,      intent(in) :: filter_relative_to_basic_state_arg

apply_pade_filter = apply_pade_filter_arg
eps_or_tau = eps_or_tau_arg
eps_filter = eps_filter_arg
tau_filter = tau_filter_arg
filter_relative_to_basic_state = filter_relative_to_basic_state_arg

if (my_node .eq. 0) then
   print *, ' In activate_pade_filter eps_or_tau = ', eps_or_tau
   print *, '    eps_filter = ', eps_filter, ' tau_filter = ', tau_filter
end if

end subroutine activate_pade_filter

!----------------------------------------------------------------------------------85

end module activate_routines

!----------------------------------------------------------------------------------85
