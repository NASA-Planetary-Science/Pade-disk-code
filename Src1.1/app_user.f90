!----------------------------------------------------------------------------------85

subroutine user_application

use activate_routines
use set_up_routines
use dof_indices, only: irho, amom, zmom, rmom, ener
! Flow field array.  Actually a "z" pencil of the flowfield array.
use q_array, only: q
#ifdef mpi_code
   use mpi
#endif

! Starting and ending indices in phi and r for pencil oriented along z.
! All the basic work of the code is done with "z" pencils.
! my_node : Rank of this processor.
use partition_data, only: sphi, ephi, sr, er, my_node

! rgrid(ir) is the r coordinate and can be used below to intialize the angular
! momentum = rho * uphi * r
use grid, only: rgrid

implicit none

! Indices for looping over the grid and looping in time:
integer :: ir, iz, iphi, istep

! Previous step completed. = 0 for non-restart run:
integer :: istep0

! Grid size and number of steps:
integer :: nr, nz, nphi, nsteps

! Time:
real(8) :: t

! Logical unit for namelist input:
integer, parameter :: lun_input = 1

! Flag for restarting:
logical :: restart

! Time and step read from restart file:
real(8) :: t_of_restart_file
integer :: istep_of_restart_file

! Temporary angular velocity variable.  Could be used to initialize angular
! momentum:
real(8) :: uphi

! ier = Error flag for mpi calls.
! istatus = Completion status of the run for passing to the subroutine terminate_with_save
integer :: ier, istatus

! File version for the save/restart file."2" is the newer version where
! we don't store the basic state.
integer, parameter :: file_version = 2

namelist /user_input/ restart, nr, nz, nphi, nsteps

! Read namelist input:
if (my_node .eq. 0) then
   print *, ' node 0: about to open and read namelist file for app_user'      
   open (unit = lun_input, file = 'input_file', form = 'formatted', &
        status = 'old')
   read (lun_input, nml = user_input)
   close (lun_input)
end if

#ifdef mpi_code
   call mpi_bcast(restart, 1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(nr,      1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nz,      1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nphi,    1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nsteps,  1, mpi_integer, 0, mpi_comm_world, ier)   
#endif   

! Set-up routines are located in set_up_routines.f90.  All thre must be called.

!call set_up_domain_mesh_and_partition(rmin, rmax, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
!     suppress_z_derivatives_when_nz_not_1_arg, periodic_z_arg, &
!     stretched_r, stretched_z, r0_unif, nr_u, z0, nz_u)

!call set_up_thermal_parameters(gamma, isothermal, apply_newtonian_cooling, &
!           tau_newtonian_cooling, ci_squared_initial)

!call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, &
!     d_ci_dr_inner, d_ci_dr_outer, c_sound_rmin, c_sound_rmax, isothermal)

! Optional activate routines.  They are located in activate_routines.f90

!call activate_gravity(gravity_flag_in, gravity_type_in, GM_in, gz_uniform_in)

!call activate_fargo(integer_shifts_arg, apply_fargo_extra_operator_arg)

!call activate_artificial_pressure (C_ap_arg)

!call activate_pade_filter(apply_pade_filter_arg, eps_or_tau_arg, &
!     eps_filter_arg, tau_filter_arg, filter_relative_to_basic_state_arg, filtering_interval_arg)

!call activate_viscosity(apply_viscosity_arg, viscosity_type_arg, isothermal_arg, &
!     nu_molecular_arg, nu_b_molecular_arg, Pr_molecular_arg, gamma_arg, C_DDSV_arg, C_Smag_arg)

!call activate_plotting_shift(Omega0_plotting_shift_arg, t_of_previous_shift_arg)

! Initial condition:
if (.not. restart) then   
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, iz, irho) = 0.d0         
            q(ir, iphi, iz, rmom) = 0.d0
            q(ir, iphi, iz, zmom) = 0.d0
            ! Angular velocity:
            uphi = 0.d0
            q(ir, iphi, iz, amom) = q(ir, iphi, iz, irho) * uphi * rgrid(ir)
            ! Internal energy = pressure / (gamma - 1)
            q(ir, iphi, iz, ener) = 0.d0
         end do
      end do
   end do
else
   call read_restart(file_version, istep_of_restart_file, t_of_restart_file)
end if

if (.not. restart) then
   istep0 = 0
   t      = 0.0
else
   istep0 = istep_of_restart_file
   t      = t_of_restart_file
end if

! Assume normal completion:
istatus = 0
do istep = istep0 + 1, istep0 + nsteps
   !call rk4(istep, q, cfl, t, dt_used, use_supplied_dt, dt_supplied, unphysical, iphi_bad, &
   !         hit_target, t_target, target_met)
end do

! istatus:
! 0 : Normal completion. Finished nsteps
! 1 : Non-normal return.  Either unphysical or dt < dt_min
! 2 : hit_target = .true. and target_met = .true.  We set this to non-zero
!     to prevent an auto restart by the pbs script.

! Need the istep - 1 since a fortran do loop increments the counter at the end
! istatus is written to return_status.dat
call terminate_with_save(istatus, istep-1, t)

end subroutine user_application

!----------------------------------------------------------------------------------85
