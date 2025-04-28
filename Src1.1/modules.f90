!----------------------------------------------------------------------------------85

module precision
#ifdef mpi_code
   use mpi
   integer, parameter :: my_mpi_real = mpi_double_precision
#endif
   integer, parameter :: xp = 8   
end module precision

module total_allocated_words
   integer :: n_words_allocated
end module total_allocated_words

module control_parameters
   implicit none
   integer :: istep, nsteps
   real(8) :: t, cfl
   integer, parameter :: restart_file = 1, save_file = 2
end module control_parameters   

module cpu_timing_module
   integer :: transpose_count, start_count_for_routine, end_count_for_routine
   real(8) :: total_cpu_for_stepping, total_cpu_for_transposes
end module cpu_timing_module

! Use these to improve code readability.
module dof_indices
   implicit none
   integer, parameter :: irho = 1, rmom = 2, zmom = 3, amom = 4, ener = 5
   ! For second-rank symmetric tensors:
   integer, parameter :: zz = 1, zr = 2, phiz = 3, rr = 4, phir = 5, phiphi = 6, pz = 3, pr = 5, pp = 6
   ! For spatial direction indices:
   integer, parameter :: z_comp = 1, r_comp = 2, phi_comp = 3, p_comp = 3

   ! For primitive variables (with irho = 1):
   integer, parameter :: iur = 2, iuz = 3, iuphi = 4, ip = 5
end module dof_indices

module grid
   !> Contains data pertaining to the grid:
   implicit none
   !> Flag to indicate that the grid has been set-up 
   logical :: gridded
   !> Flag to indicate that the z direction is periodic   
   logical :: periodic_z
   !> Number of flow degrees of freedom per grid point   
   integer, parameter :: ndof = 5
!> Number of grid points in the specified direction   
   integer :: nr, nz, nphi        

   !> Widths of each grid cell:
   real(8), dimension(:), allocatable :: dr, dz, r_dphi 
   real(8) :: dphi, dz_periodic
   
   !> Grid point locations
   real(8), dimension(:), allocatable :: rgrid, zgrid, phi_grid

   !> Jacobians
   real(8), dimension(:), allocatable :: Jr, Jz

   !> Inverse Jacobians evaluated at grid points:
   real(8), dimension(:), allocatable :: Ji_r, Ji_z
   real(8) :: dphi_inv

   !> Domain specification
   real(8) :: rmin, rmax, zmin, zmax, phi_min, phi_max, Delta_phi_domain

   !> For time step computation due to Yoshizawa model.  These arrays are dimensioned (sr:er,nz) and
   !> is assigned values in subroutine make_grid
   real(8), dimension(:,:), allocatable :: min_grid_size, l_grid_squared

   !> Grid spacing vector for Vreman SGS model:
   real(8), allocatable, dimension(:, :, :) :: Delta_vector
   
   !> Weights for trapezoidal rule integration:
   real(8), dimension(:), allocatable :: trap_weight_z, trap_weight_r

   !> This flag is needed to run an nz = 1 case, i.e., (r, phi) plane case in
   !> parallel.  We fake it by using nz = 2 and suppress z derivatives and smoothing.
   logical :: suppress_z_derivatives_when_nz_not_1
   integer :: nz_actual

   !> Useful for outputting the midplane and mid-radius.
   integer :: iz_mid, ir_mid
end module grid   

module logical_units
   implicit none
   integer, parameter :: lun_general_purpose = 1, lun_tecplot = 2
   ! This file will be open throughout the run to output the maximum
   ! time stepping eigenvalues.  So do not use this unit for anything
   ! else.
   integer, parameter :: lun_lambda = 3
   integer, parameter, dimension(9) :: lun_profile = [11,12,13,14,15,16,17,18,19]
   integer, parameter, dimension(9) :: lun_history = [21,22,23,24,25,26,27,28,29]
   integer, parameter :: lun_general_purpose_1 = 30, lun_general_purpose_2 = 31
   integer, parameter :: lun_probe = 32 ! Open throughout a run if probe is applied.
   ! History of |dilatation|_max in subroutine attack_dilatation.
   integer, parameter :: lun_dil_abs_max = 33
end module logical_units

module math_constants
   implicit none
   !> Given its value in domain_and_mesh_partition.f90
   real(8) :: pi
   real(8), parameter :: athird = 1.d0/3.d0, asixth = 1.d0/6.d0
end module math_constants

!module pade_filter
!   logical :: apply_pade_filter, filter_relative_to_basic_state
!   ! Three characters ('tau' or 'eps').
!   ! 'eps' : use the specified eps_filter
!   ! 'tau' : determine eps_filter from tau.
!   character(3) :: eps_or_tau
!   real(8) eps_filter, tau_filter
!end module pade_filter   

module partition_data
   logical :: partitioned
   integer :: num_nodes, my_node
   ! Range of indices held by processors.
   integer :: sr, sz_r, sphi, sz_phi
   integer :: er, ez_r, ephi, ez_phi
   integer :: mr, mz_r, mphi, mz_phi
   integer :: nbundle_z, nbundle_r, nbundle_phi
end module partition_data

module partition_data_for_alan
   ! These should equal the variables in the line directly above.
   integer :: mx, mz_x, my, mz_y
   integer :: mpi_comm_xz, mpi_comm_yz
   integer :: nz_for_alan ! Since I don't want to use grid in Alan's routines.
   integer :: ng1, ng2
   integer :: nx, ny
end module partition_data_for_alan

module physical_constants
   implicit none
   real(8), parameter :: Gconst = 6.673d-8 ! cm^3 g^-1 s^-2

   ! Stefan-Boltzmann:
   real(8), parameter :: sigma_SB = 5.6705d-5 ! erg s^-1 cm^-2 K^-4
   ! Boltzmann:
   real(8), parameter :: kB       = 1.3807d-16 ! erg K^-1
   real(8), parameter :: Msolar   = 1.989d33 ! gm
   real(8), parameter :: Rsolar   = 69.57d9  ! cm
   real(8), parameter :: AU       = 1.496d13 ! cm
   real(8), parameter :: parsec    = 3.0857d18 ! cm
   real(8), parameter :: kilometer = 1000.d0*100.d0 ! cm 
   real(8), parameter :: kilometer_per_sec = 1000.d0*100.d0 ! cm/sec
   real(8), parameter :: Rgas     = 3.5871d7  ! cm^2 s^-2 K^-1
   real(8), parameter :: year     = 3.16d7  ! s
   ! Mass of atomic hydrogen:
   real(8), parameter :: mH       = 1.6735d-24 ! gm
   real(8), parameter :: mH2      = 2.d0 * mH  ! gm
end module physical_constants

! Conservative variable array and pressure.
module q_array
   implicit none
   real(8), allocatable, dimension(:, :, :, :) :: q
   real(8), allocatable, dimension(:, :, :   ) :: pressure   
end module q_array

! For summation-by-parts z derivative using the Carpenter et al. scheme.
! Not currently used.
module sbp_z_derivative_module
   logical :: use_sbp_z_derivative
   real(8), allocatable, dimension(:)    :: a_sbp, b_sbp, c_sbp
   real(8), allocatable, dimension(:, :) :: Q_top_sbp, Q_bot_sbp
end module sbp_z_derivative_module

! Similarly, this module was created to avoid dynamic storage allocations and
! deallocations.  These are allocated in...
module transposes_of_q_and_qdot
   real(8), allocatable, dimension(:, :, :, :) :: q_r_space,   qdot_r_space
   real(8), allocatable, dimension(:, :, :, :) :: q_phi_space, qdot_phi_space
   real(8), allocatable, dimension(:, :, :   ) :: p_phi_space, p_r_space   
end module transposes_of_q_and_qdot

module transpose_buffer
   ! allocate(trans_buffer(sx:sx+mx-1, sy:sy+my-1, nz, nvar))

   ! Declaration in z_to_x and x_to_z:
   ! b(my,mz_x,nvar,mx,0:ng1-1)  mz_x*ng1 = nz

   ! Declaration in z_to_y and y_to_z
   ! b(mx,mz_y,nvar,my,0:ng1-1)  mz_y*ng2 = nz

   real(8), allocatable, dimension(:, :, :, :) :: trans_buffer
end module transpose_buffer
   
! This module was created to avoid dynamic storage allocations and deallocations
! everytime rk4 is called.  These arrays are allocated in subroutine initialize.
module rk4_arrays
   ! Accumulator, next argument to rhs, rhs array.
   real(8), allocatable, dimension(:, :, :, :) :: q_accum, q_next_arg, qdot
end module rk4_arrays
   
module thermal_parameters
   real(8) :: gamma, gm1
   logical :: isothermal
   ! Initial isothermal speed of sound squared.  Dimensions are (nr, nz).
   ! Isothermality is imposed by setting p = rho * ci_squared_initial(r, z)
   real(8), allocatable, dimension(:, :) :: ci_squared_initial

   logical :: apply_newtonian_cooling
   real(8) :: tau_newtonian_cooling
end module thermal_parameters

! This is for plotting purposes and is created in subroutine make_grid.
module xy_coordinates_of_grid
   real(8), allocatable, dimension(:, :) :: xgrid, ygrid
end module xy_coordinates_of_grid

!----------------------------------------------------------------------------------85
