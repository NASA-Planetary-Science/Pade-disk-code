!----------------------------------------------------------------------------------85

subroutine stop_if_unphysical(q)

! Stops the execution if:
! The density if negative or The internal energy is negative in a non-isothermal
! run.

use dof_indices
use grid
use thermal_parameters
use partition_data
#ifdef mpi_code
   use mpi
#endif
implicit none
real(8), intent(in), dimension (sr:er, sphi:ephi, nz, ndof) :: q

! Local
integer :: ier, ir, iphi, iz
logical :: my_abort ! this processor is requesting abort
logical :: all_abort    ! flag for everyone to abort

my_abort  = .false.
all_abort = .false.
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         if (q(ir, iphi, iz, irho) .le. 0.d0) then
            print *, ' In stop_if_unphysical, density is <= 0'
            print *, ' iz = ', iz, ' ir = ', ir, ' iphi = ', iphi
            my_abort = .true.
            go to 100
         end if
      end do
   end do
end do

if (.not. isothermal) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            if (q(ir, iphi, iz, ener) .le. 0.d0) then
               print *, ' In stop_if_unphysical, internal energy is <= 0'
               print *, ' iz = ', iz, ' ir = ', ir, ' iphi = ', iphi
               my_abort = .true.
               go to 100
            end if
         end do
      end do
   end do
end if

100 continue

#ifdef mpi_code
   call mpi_allreduce(my_abort, all_abort, 1, mpi_logical, mpi_lor, mpi_comm_world, ier)
#else
   all_abort = my_abort
#endif
   
if (all_abort) call terminate_with_no_save(1)

end subroutine stop_if_unphysical

!----------------------------------------------------------------------------------85

subroutine check_for_unphysical(q, pressure, substep_number, unphysical, iphi_bad)

use dof_indices
use grid
use thermal_parameters
use partition_data
use logical_units
#ifdef mpi_code
   use mpi
#endif
implicit none
real(8), intent(in), dimension (sr:er, sphi:ephi, nz, ndof) :: q
real(8), intent(in), dimension (sr:er, sphi:ephi, nz      ) :: pressure
integer, intent(in)  :: substep_number
logical, intent(out) :: unphysical
integer, intent(out) :: iphi_bad

! Local
integer :: ier, ir, iphi, iz, my_bad_iphi
logical :: i_have_unphysical ! this processor

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: subroutine check_for_unphysical has been called'
   print *, ' substep_number = ', substep_number
end if
#endif

i_have_unphysical = .false.
unphysical        = .false.
my_bad_iphi       = 0
iphi_bad          = 0

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! Check for negative density:
         if (q(ir, iphi, iz, irho) .lt. 0.0d0) then
            print *, ' density is negative at sub-step = ', substep_number
            print *, ' ir = ', ir, ' iz = ', iz, ' iphi = ', iphi
            i_have_unphysical = .true.
            my_bad_iphi = iphi
            goto 100
         end if

         ! Check for negative pressure:
         !if (pressure(ir, iphi, iz) .lt. 0.d0) then
         !   print *, ' pressure is negative at sub-step = ', substep_number
         !   print *, ' ir = ', ir, ' iz = ', iz, ' iphi = ', iphi
         !   i_have_unphysical = .true.
         !   my_bad_iphi = iphi
         !   goto 100
         !end if
      end do   
   end do
end do

if (.not. isothermal) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            if (q(ir, iphi, iz, ener) .lt. 0.0d0) then
               print *, ' internal energy is negative at sub-step = ', substep_number
               print *, ' ir = ', ir, ' iz = ', iz, ' iphi = ', iphi
               i_have_unphysical = .true.
               my_bad_iphi = iphi
               goto 100
            end if
         end do   
      end do
   end do
end if

100 continue

#ifdef mpi_code
   call mpi_allreduce(i_have_unphysical, unphysical, 1, &
        mpi_logical, mpi_lor, mpi_comm_world, ier)
   call mpi_allreduce(my_bad_iphi, iphi_bad, 1, &
                      mpi_integer, mpi_max, mpi_comm_world, ier)   
#else
   unphysical = i_have_unphysical
   iphi_bad   = my_bad_iphi
#endif

end subroutine check_for_unphysical

!----------------------------------------------------------------------------------85

subroutine stop_if_nan(q)

! Stops the execution if any q variable is Nan.

use dof_indices
use grid
use thermal_parameters
use partition_data
#ifdef mpi_code
   use mpi
#endif
implicit none
real(8), intent(in), dimension (sr:er, sphi:ephi, nz, ndof) :: q

! Local
integer :: ier, ir, iphi, iz
logical :: my_abort ! this processor is requesting abort
logical :: all_abort    ! flag for everyone to abort

my_abort  = .false.
all_abort = .false.
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         if (isnan(q(ir, iphi, iz, irho))) then
            print *, ' In stop_if_nan, density is Nan'
            print *, ' iz = ', iz, ' ir = ', ir, ' iphi = ', iphi
            my_abort = .true.
            go to 100
         end if
      
         if (isnan(q(ir, iphi, iz, rmom))) then
            print *, ' In stop_if_nan, r-momentum is Nan'
            print *, ' iz = ', iz, ' ir = ', ir, ' iphi = ', iphi
            my_abort = .true.
            go to 100
         end if

         if (isnan(q(ir, iphi, iz, zmom))) then
            print *, ' In stop_if_nan, z-momentum is Nan'
            print *, ' iz = ', iz, ' ir = ', ir, ' iphi = ', iphi
            my_abort = .true.
            go to 100
         end if

         if (isnan(q(ir, iphi, iz, amom))) then
            print *, ' In stop_if_nan, z-momentum is Nan'
            print *, ' iz = ', iz, ' ir = ', ir, ' iphi = ', iphi
            my_abort = .true.
            go to 100
         end if

         if (.not. isothermal) then
            if (isnan(q(ir, iphi, iz, ener))) then
               print *, ' In stop_if_nan, z-momentum is Nan'
               print *, ' iz = ', iz, ' ir = ', ir, ' iphi = ', iphi
               my_abort = .true.
               go to 100
            end if
         end if
      end do
   end do
end do

100 continue

#ifdef mpi_code
   call mpi_allreduce(my_abort, all_abort, 1, mpi_logical, mpi_lor, mpi_comm_world, ier)
#else
   all_abort = my_abort
#endif
   
if (all_abort) call terminate_with_no_save(1)

end subroutine stop_if_nan

!----------------------------------------------------------------------------------85

subroutine sanity_check

use partition_data
use boundary_condition_types
use boundary_condition_data
use thermal_parameters
use grid
use fargo_or_plotting_shift
implicit none

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: sanity_check'
#endif

if (isothermal) then
   if (isnan(ci_squared_initial(1, 1))) then
      if (my_node .eq. 0) then
         print *, ' In subroutine sanity_check'         
         print *, ' isothermal = ', isothermal
         print *, ' but ci_squared_initial is Nan'
      end if
      call terminate_with_no_save(1)
   end if
end if

if (isothermal) then
   if (ci_squared_initial(1, 1) .eq. 0.0d0) then
      if (my_node .eq. 0) then
         print *, ' In subroutine sanity_check'
         print *, ' isothermal = ', isothermal
         print *, ' but ci_squared_initial is zero'
      end if
      call terminate_with_no_save(1)
   end if
end if

if (rmin_BC .eq. viscous_wall) then
   if (.not. specify_viscous_wall_conditions_was_called) then
      if (my_node .eq. 0) then
         print *, ' In subroutine sanity_check'         
         print *, ' rmin_BC = viscous_wall but subroutine specify_viscous_wall_conditions'
         print *, ' was not called.'
      end if
      call terminate_with_no_save(1)
   end if
end if

if (rmax_BC .eq. viscous_wall) then
   if (.not. specify_viscous_wall_conditions_was_called) then
      if (my_node .eq. 0) then
         print *, ' In subroutine sanity_check'         
         print *, ' rmax_BC = viscous_wall but subroutine specify_viscous_wall_conditions'
         print *, ' was not called.'
      end if
      call terminate_with_no_save(1)
   end if
end if

if ((nphi .eq. 1) .and. (apply_fargo)) then
   apply_fargo = .false.
   if (my_node .eq. 0) then
         print *, ' In subroutine sanity_check'      
         print *, ' Since nphi = 1, apply_fargo was set to false' 
   end if
end if

end subroutine sanity_check

!----------------------------------------------------------------------------------85

subroutine terminate_with_save(istatus, istep, t)

#ifdef mpi_code
   use mpi
#endif
use logical_units
use partition_data
implicit none
integer, intent(in) :: istatus ! 0 = normal, 1 = abnormal
integer, intent(in) :: istep
real(8), intent(in) :: t

! Local:
character(25) :: filename_given ! output from write_save_file
integer :: ier

call write_save_file(istep, t, filename_given)

if (my_node .eq. 0) then
   open(unit = lun_general_purpose_1, file = 'return_status.dat', form = 'formatted', &
        status = 'unknown')
   write(lun_general_purpose_1, "(i1)") istatus
   close(lun_general_purpose_1)

   open(unit = lun_general_purpose_1, file = 'name_of_save_file_at_termination.dat', &
        form = 'formatted', &
        status = 'unknown')
   write(lun_general_purpose_1, *) filename_given
   close(lun_general_purpose_1)   
end if

#ifdef mpi_code
   ! This is to make sure that all processors finish what they are doing before
   ! we wrap up.
   call mpi_barrier(mpi_comm_world, ier)
   call mpi_finalize(ier)
#endif

stop
end subroutine terminate_with_save

!----------------------------------------------------------------------------------85

subroutine terminate_with_no_save(istatus)

#ifdef mpi_code
   use mpi
#endif
use logical_units
use partition_data
implicit none
integer, intent(in) :: istatus ! 0 = normal, 1 = abnormal

integer :: ier

if (my_node .eq. 0) then
   open(unit = lun_general_purpose_1, file = 'return_status.dat', form = 'formatted', &
        status = 'unknown')
   write(lun_general_purpose_1, "(i1)") istatus
   close(lun_general_purpose_1)
end if

#ifdef mpi_code
   ! This is to make sure that all processors finish what they are doing before
   ! we wrap up.
   call mpi_barrier(mpi_comm_world, ier)
   call mpi_finalize(ier)
#endif

stop
end subroutine terminate_with_no_save

!----------------------------------------------------------------------------------85

subroutine read_namelist_run_type (lun_general_purpose, i_run_type)

#ifdef mpi_code
   use mpi
#endif
use partition_data
implicit none
integer :: lun_general_purpose, i_run_type, ier

namelist /run_type/ i_run_type

! Read namelist input:

if (my_node .eq. 0) then   
   open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', status = 'old')
   read (lun_general_purpose, nml = run_type)
   close (lun_general_purpose)
end if

#ifdef mpi_code
   call mpi_bcast(i_run_type, 1, mpi_integer, 0, mpi_comm_world, ier)
#endif

if (my_node .eq. 0) then
   print *, ' node 0: i_run_type = ', i_run_type, ' in input_file'
end if
   
end subroutine read_namelist_run_type

!----------------------------------------------------------------------------------85

subroutine print_max(q, t, var_abs_max)

! In rank 0 only, obtains and prints (to stdout) the max values of
! rho, |ur|, |uz|, |uphi| and |eint| in the domain.

#ifdef mpi_code
use mpi
#endif

use grid
use partition_data
use physical_constants
use dof_indices
use logical_units
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in) :: q
real(8),                                        intent(in) :: t
real(8), dimension(5) :: var_abs_max

! Local:
integer :: iz, iphi, ir, ivar, idir
real(8), dimension(5) :: var
integer, dimension(5, 3) :: loc

#ifdef mpi_code
integer :: status(mpi_status_size)
real(8), dimension(5,    0: num_nodes-1) :: var_abs_max_array
integer, dimension(5, 3, 0: num_nodes-1) :: loc_array
integer :: ier, inode
#endif

!print *, ' my_node = ', my_node, ' print_max has been called'

do ivar = 1, 5
   var_abs_max(ivar) = 0.0d0
   loc(ivar, 1) = 0
   loc(ivar, 2) = 0
   loc(ivar, 3) = 0
end do

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         var(1) = q(ir, iphi, iz, irho)
         var(2) = q(ir, iphi, iz, rmom) / var(1)
         var(3) = q(ir, iphi, iz, zmom) / var(1)
         var(4) = q(ir, iphi, iz, amom) / var(1) / rgrid(ir)
         var(5) = q(ir, iphi, iz, ener)

         do ivar = 1, 5
            if (abs(var(ivar)) .gt. var_abs_max(ivar)) then
               var_abs_max(ivar) = abs(var(ivar))
               loc(ivar, 1) = ir
               loc(ivar, 2) = iz
               loc(ivar, 3) = iphi

               !print *, ' my_node = ', my_node
               !print *, ' ivar = ', ivar
               !print *, ' var_abs_max(ivar) = ', var_abs_max(ivar)
               !print *, ' loc(ivar, 1) = ', loc(ivar, 1)
               !print *, ' loc(ivar, 2) = ', loc(ivar, 2)
               !print *, ' loc(ivar, 3) = ', loc(ivar, 3)
            end if
         end do
      end do
   end do
end do

#ifdef mpi_code
call mpi_send(var_abs_max,  5, mpi_double,  0,  my_node,           mpi_comm_world, ier)
call mpi_send(loc,         15, mpi_integer, 0,  num_nodes+my_node, mpi_comm_world, ier)

if (my_node .eq. 0) then
   do inode = 0, num_nodes-1
      call mpi_recv(var_abs_max_array(1,    inode),    5, mpi_double,  inode, inode,           mpi_comm_world, status, ier)
      call mpi_recv(loc_array        (1, 1, inode),   15, mpi_integer, inode, num_nodes+inode, mpi_comm_world, status, ier)
   end do

   do ivar = 1, 5
      var_abs_max(ivar) = 0.0d0
      do inode = 0, num_nodes-1
         if (var_abs_max_array(ivar, inode) .gt. var_abs_max(ivar)) then
            var_abs_max(ivar) = var_abs_max_array(ivar, inode)
            do idir = 1, 3
               loc(ivar, idir) = loc_array(ivar, idir, inode)
            end do
         end if
      end do      
   end do
end if
#endif

if (my_node .eq. 0) then
   print *, ' Output from subroutine print_max, my_node   = ', my_node, ' maxima are:'
   print *, ' rho_max     = ', var_abs_max(1), ' at ir = ', loc(1,1), ' iz = ', loc(1,2), ' iphi = ', loc(1, 3)
   print *, ' |ur|_max    = ', var_abs_max(2), ' at ir = ', loc(2,1), ' iz = ', loc(2,2), ' iphi = ', loc(2, 3)
   print *, ' |uz|_max    = ', var_abs_max(3), ' at ir = ', loc(3,1), ' iz = ', loc(3,2), ' iphi = ', loc(3, 3)
   print *, ' |uphi|_max  = ', var_abs_max(4), ' at ir = ', loc(4,1), ' iz = ', loc(4,2), ' iphi = ', loc(4, 3)
   print *, ' |eint|_max  = ', var_abs_max(5), ' at ir = ', loc(5,1), ' iz = ', loc(5,2), ' iphi = ', loc(5, 3)            
   ! write(lun_history(2), "(2(1x, e12.5))")      t/year, var_abs_max(1)
   ! write(lun_history(3), "(1x, e12.5, 1x, i7)") t/year, loc(1,1)    
end if

end subroutine print_max

!----------------------------------------------------------------------------------85
