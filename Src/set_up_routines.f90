!----------------------------------------------------------------------------------85

module set_up_routines

contains

!----------------------------------------------------------------------------------85

subroutine set_up_thermal_parameters(gamma_arg, isothermal_arg, apply_newtonian_cooling_arg, &
     tau_newtonian_cooling_arg, ci_squared_initial_arg)

use grid
use thermal_parameters
use logical_units
use partition_data
use total_allocated_words, only: n_words_allocated
implicit none
!> gamma_arg : Ratio of specific heats
real(8), intent(in) :: gamma_arg
!> Boolean to run locally isothermal using the isothermal speed of sound, ci0(ir), below
logical, intent(in) :: isothermal_arg

!> Boolean to apply Newtonian cooling.
logical, intent(in) :: apply_newtonian_cooling_arg

!> Time scale parameter for Newtonian cooling.
real(8), intent(in) :: tau_newtonian_cooling_arg

!> Isothermal speed of sound squared used to set the pressure for the isothermal option, or,
!> to set the target pressure for the Newtonian cooling option
!> Ignore if isothermal = .false. or apply_newtonian_cooling = .false. or not present.
real(8), intent(in), dimension(nr, nz) :: ci_squared_initial_arg

! Local:
integer :: ir, iz

gamma                     = gamma_arg
isothermal                = isothermal_arg

apply_newtonian_cooling = apply_newtonian_cooling_arg
apply_newtonian_cooling = .false.

! This if check is necessary since tau_newtonian_cooling is an optional argument:
if (apply_newtonian_cooling) then
   tau_newtonian_cooling = tau_newtonian_cooling_arg
end if

gm1 = gamma - 1.d0

! Note: This is for the whole field.  I am having each processor have
! all of r.  I am not sure why I did this.  I'll try to fix it later.
! May be you did it so it could be used in all spaces.
if (isothermal .or. apply_newtonian_cooling) then
   allocate (ci_squared_initial(nr, nz))
   n_words_allocated = n_words_allocated + nr*nz
   ci_squared_initial = ci_squared_initial_arg
end if
   
#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' my_node = 0: set_up_thermal_parameters: Finished allocating ci_squared_initial'
      print *, ' # Bytes total = ', n_words_allocated/1e6*8, ' Mb'
   end if
#endif

if (my_node .eq. 0) then
   write(6, *) ' In subroutine set_up_thermal_parameters:'
   write(6, *) ' isothermal = ', isothermal
end if   

end subroutine set_up_thermal_parameters

!----------------------------------------------------------------------------------85

subroutine set_up_boundary_conditions(rmin_BC_arg, rmax_BC_arg, zmin_BC_arg, &
     zmax_BC_arg, ibalanced_arg, d_ci_dr_inner_arg, d_ci_dr_outer_arg, &
     c_sound_rmin_arg, c_sound_rmax_arg, isothermal_arg)

use grid
use partition_data
use boundary_condition_data
use boundary_condition_types
use logical_units
use total_allocated_words, only: n_words_allocated
use thermal_parameters
implicit none
integer, intent(in) :: rmin_BC_arg, rmax_BC_arg, zmin_BC_arg, zmax_BC_arg, ibalanced_arg
!> These need to be specified only if you want balanced non-reflective boundary conditions in
!> the radial direction.
real(8), intent(in), dimension(nz) :: d_ci_dr_inner_arg, d_ci_dr_outer_arg

! These are proxies for specifying the temperature at rmin and rmax when we are running
! adiabatically (non-isothermlly).
real(8) :: c_sound_rmin_arg, c_sound_rmax_arg
logical :: isothermal_arg

! Local:
integer :: iz

rmin_BC = rmin_BC_arg
rmax_BC = rmax_BC_arg
zmin_BC = zmin_BC_arg
zmax_BC = zmax_BC_arg
ibalanced = ibalanced_arg

if ((rmin_BC .eq. viscous_wall) .and. (.not. isothermal_arg)) c_sound_rmin = c_sound_rmin_arg
if ((rmax_BC .eq. viscous_wall) .and. (.not. isothermal_arg)) c_sound_rmax = c_sound_rmax_arg

if (zmin_BC .eq. z_dirichlet) then
   allocate(qLeft_z_BC(ndof))
end if

if (zmax_BC .eq. z_dirichlet) then
   allocate(qRight_z_BC(ndof))
end if

if ((ibalanced .eq. 0) .or. (ibalanced .eq. 1)) then
   if (my_node .eq. 0) write(6, *) ' ibalanced = ', ibalanced
else
   if (my_node .eq. 0) then
      print *, ' ibalanced = ', ibalanced, ' must be 0 or 1.  Aborting'
   end if
   call terminate_with_no_save(1)
end if

if (periodic_z .and. (zmin_BC .ne. periodic)) then
   write (6, 3) zmin_BC
3  format (' Changed zmin_BC from ', i2, ' to periodic')
   zmin_BC = periodic
end if

if (periodic_z .and. (zmin_BC .ne. periodic)) then
   write(6, 4) zmax_BC
4  format (' Changed zmax_BC from ', i2, ' to periodic')
   zmax_BC = periodic
end if

! This is for boundary conditions.
if (isothermal) then
   allocate(d_ci_dr_inner(nz), d_ci_dr_outer(nz))
   n_words_allocated = n_words_allocated + 2*nz

   do iz = 1, nz
      d_ci_dr_inner(iz) = d_ci_dr_inner_arg(iz)
      d_ci_dr_outer(iz) = d_ci_dr_outer_arg(iz)
   end do
end if

if (my_node .eq. 0) then
   write(6, 1) rmin_BC, rmax_BC, zmin_BC, zmax_BC
   1 format (' rmin_BC = ', i2,/, &
             ' rmax_BC = ', i2,/, &
             ' zmin_BC = ', i2,/, &
             ' zmax_BC = ', i2)

   write(6, 2) null, non_reflective, zero_normal_momentum, &
        Cassen_Moosman_BC, periodic, viscous_wall
   2 format (' Codes are:             ',     /, &
             ' null                 = ', i2, /, &
             ' non_reflective       = ', i2, /, &
             ' zero_normal_momentum = ', i2, /, &
             ' Cassen-Moosman       = ', i2, /, &
             ' periodic             = ', i2, /, &
             ' viscous_wall         = ', i2)
end if

if (my_node .eq. 0) print *, ' node 0: About to return from set-up_boundary_conditions'

end subroutine set_up_boundary_conditions

!----------------------------------------------------------------------------------85

subroutine set_up_viscous_wall_conditions(Omega_rmin_arg, Omega_rmax_arg, &
     uz_rmin_arg, uz_rmax_arg, c_sound_rmin_arg, c_sound_rmax_arg)

! This routine should be called if rmin_BC or rmax_BC = viscous wall.
! The speed of sound at the wall is used as a proxy for temperature.

use boundary_condition_data
implicit none
real(8) :: Omega_rmin_arg, Omega_rmax_arg, uz_rmin_arg, uz_rmax_arg, &
     c_sound_rmin_arg, c_sound_rmax_arg

Omega_rmin   = Omega_rmin_arg
Omega_rmax   = Omega_rmax_arg
uz_rmin      = uz_rmin_arg
uz_rmax      = uz_rmax_arg
c_sound_rmin = c_sound_rmin_arg
c_sound_rmax = c_sound_rmax_arg

specify_viscous_wall_conditions_was_called = .true.

return
end subroutine set_up_viscous_wall_conditions

!----------------------------------------------------------------------------------85

end module set_up_routines

!----------------------------------------------------------------------------------85
