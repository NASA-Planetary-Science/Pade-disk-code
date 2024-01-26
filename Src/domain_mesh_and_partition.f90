!----------------------------------------------------------------------------------85>

subroutine set_up_domain_mesh_and_partition(rmin_arg, rmax_arg, zmin_arg, zmax_arg,&
     phi_min_arg, phi_max_arg, nr_arg, nz_arg, nphi_arg, &
     suppress_z_derivatives_when_nz_not_1_arg, periodic_z_arg, stretched_r_arg, &
     stretched_z_arg, r0_arg, nr_u_arg, z0_arg, nz_u_arg)

use grid
use partition_data
use math_constants
use total_allocated_words, only: n_words_allocated
use q_array
use rk4_arrays
use transposes_of_q_and_qdot
use transpose_buffer
use basic_state
use cpu_timing_module
use fargo_or_plotting_shift
use pade_filter
use logical_units
use viscous
use thermal_parameters ! To set defaults
use boundary_condition_data
use boundary_condition_types
use stretched_mesh
use sponge ! To set defaults
use gravity ! To set-up defaults
use partition_data_for_alan

implicit none
real(8) :: rmin_arg, rmax_arg, zmin_arg, zmax_arg, phi_min_arg, phi_max_arg
integer :: nr_arg, nz_arg, nphi_arg
logical :: suppress_z_derivatives_when_nz_not_1_arg, periodic_z_arg
logical :: stretched_r_arg, stretched_z_arg
real(8) :: r0_arg, z0_arg
integer :: nr_u_arg, nz_u_arg

if (my_node .eq. 0) then
   print *, ' my_node=0: subroutine set_up_domain_mesh_and_partition called'
   print *, '      nr_arg = ', nr_arg, ' nz_arg = ', nz_arg, ' nphi_arg = ', nphi_arg
   print *, '      rmin_arg = ', rmin_arg, ' rmax_arg = ', rmax_arg
   print *, '      zmin_arg = ', zmin_arg, ' zmax_arg = ', zmax_arg
   print *, '      phi_min_arg = ', phi_min_arg, ' phi_max_arg = ', phi_max_arg
   print *, '      suppress_z_derivatives_when_nz_not_1_arg = ', suppress_z_derivatives_when_nz_not_1_arg
   print *, '      periodic_z_arg  = ', periodic_z_arg

   print *, '      stretched_r_arg = ', stretched_r_arg
   print *, '      stretched_z_arg = ', stretched_z_arg
   if (stretched_r_arg) then
      print *, '   Stretching parameters for r-direction:'
      print *, '      r0_arg   = ', r0_arg
      print *, '      nr_u_arg = ', nr_u_arg
   end if

   if (stretched_z_arg) then
      print *, ' Stretching parameters for z-direction'
      print *, '      z0_arg   = ', z0_arg
      print *, '      nz_u_arg = ', nz_u_arg
   end if

   ! Check that none of the domain arguments is Nan:
   if (isnan(rmin_arg) .or. isnan(rmax_arg) .or. isnan(zmin_arg) .or. isnan(zmax_arg) &
        .or. isnan(phi_min_arg) .or. isnan(phi_max_arg)) then
      if (my_node .eq. 0) then
         print *, ' node 0, subroutine set_up_domain_mesh_and_partition'
         print *, ' One of these is Nan:'
         print *, ' rmin_arg    = ', rmin_arg,    ' rmax_arg    = ', rmax_arg
         print *, ' zmin_arg    = ', zmin_arg,    ' zmax_arg    = ', zmax_arg
         print *, ' phi_min_arg = ', phi_min_arg, ' phi_max_arg = ', phi_max_arg
      end if
      call terminate_with_no_save(1)
   end if
end if

rmin    = rmin_arg
rmax    = rmax_arg
zmin    = zmin_arg
zmax    = zmax_arg
phi_min = phi_min_arg
phi_max = phi_max_arg
nr   = nr_arg
nz   = nz_arg
nphi = nphi_arg
suppress_z_derivatives_when_nz_not_1 = suppress_z_derivatives_when_nz_not_1_arg
periodic_z = periodic_z_arg

stretched_r = stretched_r_arg
stretched_z = stretched_z_arg
r0   = r0_arg
z0   = z0_arg
nr_u = nr_u_arg
nz_u = nz_u_arg

if (suppress_z_derivatives_when_nz_not_1) then
   nz_actual = 1
else
   nz_actual = nz
end if

if (my_node .eq. 0) then
   print *, ' nr = ', nr, ' nz = ', nr, ' nphi = ', nphi
end if

! This sits is in module control_parameters and keeps track of how much memory
! has been allocated:
n_words_allocated = 0

! Assign this variable of module math_constants:
pi     = 4.0d0 * DATAN(1.0d0)

! Partitioning for parallel processing:
#ifdef mpi_code
   if (my_node .eq. 0) print *, ' node = 0: calling partition'
   call partition
   if (my_node .eq. 0) print *, ' node = 0: calling make_communicators_for_alans_transposes'
   call make_communicators_for_alans_transposes ! in transpose_routines.f90
   if (my_node .eq. 0) print *, ' node = 0: returned from make_communicators_for_alans_transposes'   
#else
   ! For non-parallel code:
   my_node = 0
   sr     = 1
   sz_r   = 1
   sphi   = 1
   sz_phi = 1
   
   er     = nr
   ez_r   = nz
   ephi   = nphi
   ez_phi = nz

   mr     = nr
   mz_r   = nz
   mphi   = nphi
   mz_phi = nz
   
   partitioned = .true.
#endif

#ifdef mpi_code
   if (my_node .eq. 0) then
      print *, ' mpi run with num_nodes = ', num_nodes
      print *, ' cpu layout: ng1 = ', ng1, ' ng2 = ', ng2
   end if
#endif   

! Bundle sizes for differentiating in each direction:   
nbundle_z   = mr   * mphi
nbundle_r   = mphi * mz_r
nbundle_phi = mr * mz_phi 
   
if (my_node .eq. 0) print *, ' node = 0: About to allocate arrays'
if (my_node .eq. 0) print *, ' node = 0: nr = ', nr, ' nphi = ', nphi, ' nz = ', nz, ' ndof = ', ndof

! q array.  In z space
allocate (q(sr:er, sphi:ephi, nz, ndof))
allocate (pressure(sr:er, sphi:ephi, nz))
n_words_allocated = n_words_allocated + mr*mphi*nz*(ndof + 1)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' my_node = 0: set_up_domain_mesh_and_partition: Finished allocating q'
#endif

! rk4 arrays.  These are also in z space:
allocate (q_accum   (sr:er, sphi:ephi, nz, ndof))
allocate (q_next_arg(sr:er, sphi:ephi, nz, ndof))
allocate (qdot      (sr:er, sphi:ephi, nz, ndof))
n_words_allocated = n_words_allocated + 3*mr*mphi*nz*ndof
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node = 0: set_up_domain_mesh_and_partition: Finished allocating rk4 arrays'
#endif

! r-space transposes:
allocate (q_r_space   (sphi:ephi, sz_r:ez_r, ndof, nr))
allocate (qdot_r_space(sphi:ephi, sz_r:ez_r, ndof, nr))
allocate (p_r_space   (sphi:ephi, sz_r:ez_r,       nr))
n_words_allocated = n_words_allocated + (2*ndof + 1)*mphi*nz**nr

! We allocate these even when nphi = 1 because we pass them as work arrays for
! vorticity and dilatation computation:
! phi-space transposes:
allocate (q_phi_space   (sr:er, sz_phi:ez_phi, ndof, nphi))
allocate (qdot_phi_space(sr:er, sz_phi:ez_phi, ndof, nphi))
allocate (p_phi_space   (sr:er, sz_phi:ez_phi,       nphi))
n_words_allocated = n_words_allocated + (2*ndof + 1)*mr*mz_phi*nphi
#ifdef debug_print
   if (my_node .eq. 0) print *, ' my_node = 0: subroutine set_up_domain_mesh_and_partition: Finished allocating transposes'
#endif

#ifdef mpi_code
   allocate(trans_buffer(sr:er, sphi:ephi, nz, ndof))
   n_words_allocated = n_words_allocated + mr*mphi*nz*ndof 
#endif

#ifdef debug_print
   if (my_node .eq. 0) print *, ' my_node = 0: subroutine set_up_domain_mesh_and_partition: Finished allocating arrays'
#endif

! Allocate and Generate the mesh.  This routine is in mesh.f90
call make_grid
#ifdef debug_print
   if (my_node .eq. 0) print *, ' my_node = 0: subroutine set_up_domain_mesh_and_partition: Returned from make_grid'
#endif

! This routine is in flux_diff.f90:
call set_up_pade_coefficients      

if (my_node .eq. 0) then
   print *, ' my_node = 0: subroutine set_up_domain_mesh_and_partition: Returned from set_up_pade...'
   print *, ' about to call get_xy_coordinates_of_grid'
end if

! This routine is in plotting_output.f90
if (nphi .ne. 1) then
   call get_xy_coordinates_of_grid ! for plotting purposes
   if (my_node .eq. 0) then
      print *, ' my_node = 0: set_up_domain_mesh_and_partition: Returned from get_xy_coordinates_of_grid'
   end if
end if

total_cpu_for_stepping   = 0.0d0

#ifdef transpose_timing
   total_cpu_for_transposes = 0.d0
#endif

! Set-up some defaults which can be over-ridden by calls to other set-up routines:
apply_fargo                  = .false.
apply_plotting_shift         = .false.
fft_has_been_initialized     = .false.
apply_pade_filter            = .false.
apply_viscosity              = .false.
specify_viscous_wall_conditions_was_called = .false.
apply_sponge                 = .false.
gravity_flag                 = .false.

! These are the defaults for Newtonian cooling which can be over-ridden by
! calling activate_newtonian_cooling(apply_newtonian_cooling, tau_newtonian_cooling)
apply_newtonian_cooling = .false.
tau_newtonian_cooling   = 0.d0

rmin_BC   = null
rmax_BC   = null
zmin_BC   = null
zmax_BC   = null

! Sanity check:
if ((nphi .ne. 1) .and. (dphi .eq. 0.d0)) then
   if (my_node .eq. 0) then
      print *, ' nphi .ne. 1 but dphi = 0'
      print *, ' You likely forgot to set phi_min and phi_max'
   end if
   call terminate_with_no_save(1)
end if

end subroutine set_up_domain_mesh_and_partition

!----------------------------------------------------------------------------------85

#ifdef mpi_code
subroutine partition

! Partition the mesh into processors.  Run the tool "part_tool" beforehand to make
! sure you get a solution for the partitioning problem.

! Look at notes of 12/29/2017.

! Alan      Me
! x    <--> r
! y    <--> phi
! z    <--> z

use grid
use partition_data
use partition_data_for_alan
use mpi
implicit none
integer :: ier

! Code copied from Alan's stellarbox code.
ng1 = 1
Do while(ng1 <= num_nodes)
   ng2 = num_nodes/ng1
   If (ng1*ng2 == num_nodes) then      ! num_nodes is divisible by these factors
      mr     = nr/ng1
      mphi   = nphi/ng2
      mz_r   = nz/ng1
      mz_phi = nz/ng2
      ! Ask Alan why the first solution (lowest ng1) is the most efficient.
      If (mr*ng1==nr .and. mphi*ng2==nphi .and. mz_r*ng1==nz .and. mz_phi*ng2==nz) go to 100
   End If
   ng1 = ng1 + 1
End Do

if (my_node .eq. 0) then
   print *, ' No solution found for partitioning'
   print *, ' nr = ', nr, ' nz = ', nz, ' nphi = ', nphi, ' num_nodes = ', num_nodes
end if
call mpi_abort(mpi_comm_world, 99, ier)

100 continue

! The left-hand side guys is what Alan's transpose routines need in their notation:
mx   = mr
my   = mphi
mz_x = mz_r
mz_y = mz_phi
nz_for_alan = nz

nx = nr
ny = nphi

partitioned = .true.
end subroutine partition
#endif

!----------------------------------------------------------------------------------85
