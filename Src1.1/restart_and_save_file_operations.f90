!----------------------------------------------------------------------------------85

subroutine write_save_file(istep, t, filename_given)

use q_array
use grid
use partition_data
use logical_units, only: lun_general_purpose
use basic_state
implicit none
integer, intent(in) :: istep
real(8), intent(in) :: t
character(25), intent(out) :: filename_given ! name given to the save file

! Local:
integer       :: ir, iz, iphi, idof

#ifdef mpi_code
   if (my_node .eq. 0) print *, ' node = 0: calling mpi_write_save_file'
   !call mpi_write_save_file (istep, t, filename_given)
   call mpi_write_save_file_version2(istep, t, filename_given)   
   if (my_node .eq. 0) print *, ' node = 0: returned from mpi_write_save_file'   
#else
! Name of save file for a serial (non-mpi)  run.
write (filename_given, "('save_version2_', i7.7)") istep
open (unit = lun_general_purpose, file = filename_given, form = 'unformatted', &
         status = 'unknown')

write (lun_general_purpose) istep, t
write (lun_general_purpose) nr, nz, nphi
write (lun_general_purpose) rmin, rmax, zmin, zmax, phi_max
write (lun_general_purpose) have_basic_state

write (lun_general_purpose) &
     ((((q(ir,iphi,iz,idof),ir=1,nr),iphi=1,nphi),iz=1,nz),idof=1,ndof)

!if (have_basic_state) then
!   write (lun_general_purpose) ((rho_basic (ir, iz),ir=1,nr),iz=1,nz)
!   write (lun_general_purpose) ((uphi_basic(ir, iz),ir=1,nr),iz=1,nz)   
!end if
close (lun_general_purpose)

print *, ' Wrote unformatted file = ', filename_given, ' for istep = ', istep, ' t = ', t
#endif

!if (my_node .eq. 0) print *, ' node 0: About to return from write_save_file'

end subroutine write_save_file

!----------------------------------------------------------------------------------85

#ifdef mpi_code
subroutine mpi_write_save_file (istep, t, filename_given)

use mpi
use partition_data
use grid
use q_array
use basic_state
implicit none
integer, intent(in) :: istep
real(8), intent(in) :: t
character(25), intent(out) :: filename_given

integer :: len, new_ier
character(132) :: string

! Local
integer(kind = mpi_offset_kind)     :: zero_offset = 0
integer, dimension(mpi_status_size) :: status
integer :: file_handle, ier

write (filename_given, "('mpi_save_', i7.7)") istep
call mpi_file_open(mpi_comm_world, filename_given, mpi_mode_create + mpi_mode_wronly, mpi_info_null, &
     file_handle, ier)

if ((ier .ne. mpi_success) .and. (my_node .eq. 0)) then
   print *, ' node 0: error opening mpi_save. ier = ', ier
   call mpi_error_string(ier, string, len, ier)
   print *, ' len = ', len
   print *, ' string = ', string
end if   

call mpi_file_set_view(file_handle, zero_offset, mpi_byte, mpi_byte, "native", mpi_info_null, ier)

call mpi_barrier(mpi_comm_world, ier)

! Only node 0 writes the header portion:
if (my_node .eq. 0) then
   call mpi_file_write_shared(file_handle, num_nodes,  1, mpi_integer, status, ier)
   call mpi_file_write_shared(file_handle, ndof,       1, mpi_integer, status, ier)
   call mpi_file_write_shared(file_handle, istep,      1, mpi_integer, status, ier)
   call mpi_file_write_shared(file_handle, t,          1, mpi_double,  status, ier)

   print *, ' wrote following to the mpi_save file'
   print *, ' num_nodes = ', num_nodes
   print *, ' ndof      = ', ndof
   print *, ' istep     = ', istep
   print *, ' t         = ', t
      
   call mpi_file_write_shared(file_handle, nr,      1, mpi_integer, status, ier)
   call mpi_file_write_shared(file_handle, nz,      1, mpi_integer, status, ier)
   call mpi_file_write_shared(file_handle, nphi,    1, mpi_integer, status, ier)

   call mpi_file_write_shared(file_handle, rmin,    1, mpi_double, status, ier)
   call mpi_file_write_shared(file_handle, rmax,    1, mpi_double, status, ier)
   call mpi_file_write_shared(file_handle, zmin,    1, mpi_double, status, ier)
   call mpi_file_write_shared(file_handle, zmax,    1, mpi_double, status, ier)
   call mpi_file_write_shared(file_handle, phi_max, 1, mpi_double, status, ier)

   call mpi_file_write_shared(file_handle, have_basic_state, 1, mpi_logical,  status, ier)
   !print *, ' node 0 finished with header'
end if ! node 0

if (have_basic_state) then
   call mpi_file_write_ordered(file_handle, rho_basic,          mr*nz,   mpi_double,  status, ier)
   !print *, ' my_node = ', my_node, ' finished writing rho_basic ier = ', ier
   !print *, ' status = ', status
   call MPI_Error_string(ier, string, len, new_ier)
   !print *, ' string = ', string
   call mpi_file_write_ordered(file_handle, uphi_basic,         mr*nz,   mpi_double,  status, ier)
   !print *, ' my_node = ', my_node, ' finished writing uphi_basic ier = ', ier

   ! Note: uphi_basic_r_space is needed for vorticity computation:
   call mpi_file_write_ordered(file_handle, uphi_basic_r_space, mz_r*nr, mpi_double,  status, ier)
   !print *, ' my_node = ', my_node, ' finished writing uphi_basic_r_space ier = ', ier
end if

call mpi_file_write_ordered(file_handle, q, mr*mphi*nz*ndof, mpi_double, status, ier)

if ((ier .ne. mpi_success) .and. (my_node .eq. 0)) then
   print *, ' node 0: error writing pencil to restart file. ier = ', ier
   call mpi_error_string(ier, string, len, ier)
   print *, ' len = ', len
   print *, ' error string = ', string
end if

call mpi_file_close(file_handle, ier)
if (my_node .eq. 0) print *, ' node 0: finished calling mpi_file_close for mpi_save'
end subroutine mpi_write_save_file
#endif

!----------------------------------------------------------------------------------85

subroutine read_restart(version, istep, t)

use q_array
use grid
use logical_units, only: lun_general_purpose
use basic_state
use partition_data
use total_allocated_words, only: n_words_allocated
implicit none
integer, intent(in)  :: version  ! 1 or 2
integer, intent(out) :: istep
real(8), intent(out) :: t

! Local:
integer :: ir, iz, iphi, idof

if (my_node .eq. 0) print *, ' my_node = 0: read_restart has been called'

#ifdef mpi_code
if (version == 2) then
   call mpi_read_restart_version2(istep, t)
else
   call mpi_read_restart(istep, t)
end if
#else
open (unit = lun_general_purpose, file = 'restart', form = 'unformatted', &
     status = 'old')
read (lun_general_purpose) istep, t
read (lun_general_purpose) nr, nz, nphi
read (lun_general_purpose) rmin, rmax, zmin, zmax, phi_max
read (lun_general_purpose) have_basic_state

! read (lun_general_purpose) (rgrid   (ir  ), ir   = 1, nr)
! read (lun_general_purpose) (zgrid   (iz  ), iz   = 1, nz)
! read (lun_general_purpose) (phi_grid(iphi), iphi = 1, nphi)

! read (lun_general_purpose) (Ji_r(ir), ir = 1, nr)
! read (lun_general_purpose) (Ji_z(iz), iz = 1, nz)
! read (lun_general_purpose) dphi, dphi_inv

read (lun_general_purpose) &
     ((((q(ir,iphi,iz,idof),ir=1,nr),iphi=1,nphi),iz=1,nz),idof=1,ndof)

! For serial we assume that we always have version 2:
!if (have_basic_state) then
!   read (lun_general_purpose) ((rho_basic (ir, iz),ir=1,nr),iz=1,nz)
!   read (lun_general_purpose) ((uphi_basic(ir, iz),ir=1,nr),iz=1,nz)   
!end if

close (lun_general_purpose)
#endif

end subroutine read_restart

!----------------------------------------------------------------------------------85

#ifdef mpi_code
subroutine mpi_read_restart(istep, t)

use mpi
use partition_data
use grid
use q_array
use basic_state
use total_allocated_words
implicit none
integer, intent(out) :: istep
real(8), intent(out) :: t
integer :: len
character(80) :: string

! Local
integer(kind = mpi_offset_kind)     :: zero_offset = 0
integer, dimension(mpi_status_size) :: status
integer :: file_handle, ier

integer :: num_nodes_read, ndof_read, nr_read, nz_read, nphi_read

call mpi_file_open(mpi_comm_world, 'mpi_restart', mpi_mode_rdonly, mpi_info_null, &
     file_handle, ier)
call my_mpi_error (my_node, ier, 'opening mpi_restart ')
call mpi_file_set_view(file_handle, zero_offset, mpi_byte, mpi_byte, "native", mpi_info_null, ier)
call my_mpi_error (my_node, ier, 'set_view            ')

! Only node 0 reads the header portion:
if (my_node .eq. 0) then
   call mpi_file_read_shared(file_handle, num_nodes_read,  1, mpi_integer, status, ier)
   call mpi_file_read_shared(file_handle, ndof_read,       1, mpi_integer, status, ier)
   call mpi_file_read_shared(file_handle, istep,           1, mpi_integer, status, ier)
   call mpi_file_read_shared(file_handle, t,               1, mpi_double,  status, ier)

   print *, ' node 0: Read the following from restart file:'
   print *, '    num_nodes_read = ', num_nodes_read
   print *, '    ndof_read      = ', ndof_read
   print *, '    istep          = ', istep
   print *, '    t              = ', t
   
   if (num_nodes_read .ne. num_nodes) then
      print *, ' node 0: Number of nodes in restart = ', num_nodes_read
      print *, '       does not match current value = ', num_nodes
      call mpi_abort(mpi_comm_world, 99, ier)
      stop
   end if   

   if (ndof_read .ne. ndof) then
      print *, ' node 0: ndof in restart      = ', ndof_read
      print *, ' does not match current value = ', ndof
      call mpi_abort(mpi_comm_world, 99, ier)
      stop
   end if   
   
   call mpi_file_read_shared(file_handle, nr_read,      1, mpi_integer, status, ier)
   call mpi_file_read_shared(file_handle, nz_read,      1, mpi_integer, status, ier)
   call mpi_file_read_shared(file_handle, nphi_read,    1, mpi_integer, status, ier)
   if ((nr_read .ne. nr) .or. (nz_read .ne. nz) .or. (nphi_read .ne. nphi)) then
      print *, ' node 0: Grid sizes in input_file and restart file do not match'
      print *, ' Values in restart file are:'
      print *, ' nr = ', nr_read, ' nz = ', nz_read, ' nphi = ', nphi_read
      call mpi_abort(mpi_comm_world, 99, ier)
      stop
   end if
   
   print *, ' nr_read = ', nr_read, ' nz_read = ', nz_read, ' nphi_read = ', nphi_read

   call mpi_file_read_shared(file_handle, rmin,    1, mpi_double, status, ier)
   call mpi_file_read_shared(file_handle, rmax,    1, mpi_double, status, ier)
   call mpi_file_read_shared(file_handle, zmin,    1, mpi_double, status, ier)
   call mpi_file_read_shared(file_handle, zmax,    1, mpi_double, status, ier)
   call mpi_file_read_shared(file_handle, phi_max, 1, mpi_double, status, ier)
   print *, ' rmin = ', rmin, ' rmax = ', rmax, ' zmin = ', zmin, ' zmax = ', zmax
   print *, ' phi_max = ', phi_max
   call mpi_file_read_shared(file_handle, have_basic_state, 1, mpi_logical,  status, ier)
   print *, ' have_basic_state = ', have_basic_state
end if ! node 0

call mpi_bcast(istep, 1, mpi_integer, 0, mpi_comm_world, ier)
call mpi_bcast(t,     1, mpi_double,  0, mpi_comm_world, ier)
      
call mpi_bcast(nr,    1, mpi_integer,  0, mpi_comm_world, ier)
call mpi_bcast(nz,    1, mpi_integer,  0, mpi_comm_world, ier)
call mpi_bcast(nphi,  1, mpi_integer,  0, mpi_comm_world, ier)

call mpi_bcast(rmin,    1, mpi_double,  0, mpi_comm_world, ier)
call mpi_bcast(rmax,    1, mpi_double,  0, mpi_comm_world, ier)
call mpi_bcast(zmin,    1, mpi_double,  0, mpi_comm_world, ier)
call mpi_bcast(zmax,    1, mpi_double,  0, mpi_comm_world, ier)
call mpi_bcast(phi_max, 1, mpi_double,  0, mpi_comm_world, ier)

call mpi_bcast(have_basic_state, 1, mpi_logical,  0, mpi_comm_world, ier)
call mpi_bcast(rgrid(1),    nr,   mpi_double,  0, mpi_comm_world, ier)
call mpi_bcast(zgrid(1),    nz,   mpi_double,  0, mpi_comm_world, ier)
call mpi_bcast(phi_grid(1), nphi, mpi_double,  0, mpi_comm_world, ier)

call mpi_bcast(Ji_r(1), nr, mpi_double,  0, mpi_comm_world, ier)
call mpi_bcast(Ji_z(1), nz, mpi_double,  0, mpi_comm_world, ier)

call mpi_bcast(dphi,     1, mpi_double,  0, mpi_comm_world, ier)
call mpi_bcast(dphi_inv, 1, mpi_double,  0, mpi_comm_world, ier)     

if (have_basic_state) then
   ! These allocations might have been done in basic_state.f90
   if (.not.(allocated(rho_basic)))          allocate(rho_basic (sr:er, nz))
   if (.not.(allocated(uphi_basic)))         allocate(uphi_basic(sr:er, nz))
   if (.not.(allocated(uphi_basic_r_space))) allocate(uphi_basic_r_space(sz_r:ez_r, nr)) ! needed for vorticity computation
   n_words_allocated = n_words_allocated + 2*mr*nz + mz_r*nr
   
   call mpi_file_read_ordered(file_handle, rho_basic,          mr*nz,   mpi_double,  status, ier)
   call mpi_file_read_ordered(file_handle, uphi_basic,         mr*nz,   mpi_double,  status, ier)
   ! Needed for vorticity computation:
   call mpi_file_read_ordered(file_handle, uphi_basic_r_space, mz_r*nr, mpi_double,  status, ier)      
end if

! Read the flow field:
call mpi_file_read_ordered(file_handle, q, mr*mphi*nz*ndof, mpi_double, status, ier)

if ((ier .ne. mpi_success) .and. (my_node .eq. 0)) then
   print *, ' node 0: error reading pencil from restart file. ier = ', ier
   call mpi_error_string(ier, string, len, ier)
   print *, ' len = ', len
   print *, ' error string = ', string
end if

call mpi_file_close(file_handle, ier)
if (my_node .eq. 0) print *, ' node 0: finished calling mpi_file_close for restart file'
end subroutine mpi_read_restart

!----------------------------------------------------------------------------------85

subroutine my_mpi_error (my_node, ier, where20)

use mpi
implicit none
integer :: my_node, ier
character(20) :: where20

! Local:
character(80) :: string
integer :: ierror, len

if (ier .ne. mpi_success) then
   if (my_node .eq. 0) then
      print *, ' error from node = ', my_node, ' and possibly other nodes, Location = ', where20
      call mpi_error_string(ier, string, len, ierror)
      print *, ' len = ', len, ' string = ', string
   end if
end if

end subroutine my_mpi_error

#endif

!----------------------------------------------------------------------------------85

#ifdef mpi_code
subroutine mpi_write_save_file_version2(istep, t, filename_given)

use mpi
use partition_data
use partition_data_for_alan
use grid
use q_array
use basic_state
use transposes_of_q_and_qdot
implicit none
integer, intent(in) :: istep
real(8), intent(in) :: t
character(25), intent(out) :: filename_given

integer :: len, new_ier
character(132) :: string

! Local
integer(kind = mpi_offset_kind)     :: zero_offset = 0
integer, dimension(mpi_status_size) :: status
integer :: file_handle, ier

! For reference.
! q           (sr:er, sphi:ephi, nz, ndof)           : in z space
! q_r_space   (sphi:ephi, sz_r:ez_r, ndof, nr)
! q_phi_space (sr:er, sz_phi:ez_phi, ndof, nphi)

! Temps. containing a single variable:
real(8), dimension(sr:er, sphi:ephi, nz)       :: temp_z
real(8), dimension(sphi:ephi, sz_r:ez_r, nr)   :: temp_r
real(8), dimension(sr:er, sz_phi:ez_phi, nphi) :: temp_phi

real(8), dimension(nr, sz_r:ez_r) :: Qsave

! Loop indices:
integer :: i, j, k, m

integer :: ng_r, ng_phi

! Alan's notation:
ng_r   = ng1
ng_phi = ng2

write (filename_given, "('mpi_save_version2_', i7.7)") istep
call mpi_file_open(mpi_comm_world, filename_given, mpi_mode_create + mpi_mode_wronly, mpi_info_null, &
     file_handle, ier)

If (ier /= MPI_SUCCESS) then
   write(6,*) "Can't MPI_FILE_OPEN the save file "//filename_given
   Call MPI_ABORT(MPI_COMM_WORLD, 99, ier)
End If

call mpi_file_set_view(file_handle, zero_offset, mpi_byte, mpi_byte, "native", mpi_info_null, ier)

have_basic_state = .false.
! Only node 0 writes the header portion:
if (my_node .eq. 0) then
   ! 13 items in header:
   call mpi_file_write_shared(file_handle, num_nodes,  1, mpi_integer, status, ier)
   call mpi_file_write_shared(file_handle, ndof,       1, mpi_integer, status, ier)
   call mpi_file_write_shared(file_handle, istep,      1, mpi_integer, status, ier)
   call mpi_file_write_shared(file_handle, t,          1, mpi_double,  status, ier)
   call mpi_file_write_shared(file_handle, nr,      1, mpi_integer, status, ier)
   call mpi_file_write_shared(file_handle, nz,      1, mpi_integer, status, ier)
   call mpi_file_write_shared(file_handle, nphi,    1, mpi_integer, status, ier)
   call mpi_file_write_shared(file_handle, rmin,    1, mpi_double, status, ier)
   call mpi_file_write_shared(file_handle, rmax,    1, mpi_double, status, ier)
   call mpi_file_write_shared(file_handle, zmin,    1, mpi_double, status, ier)
   call mpi_file_write_shared(file_handle, zmax,    1, mpi_double, status, ier)
   call mpi_file_write_shared(file_handle, phi_max, 1, mpi_double, status, ier)
   call mpi_file_write_shared(file_handle, have_basic_state, 1, mpi_logical,  status, ier)
   print *, ' wrote following to the mpi_save_version2 file'
   print *, ' num_nodes = ', num_nodes
   print *, ' ndof      = ', ndof
   print *, ' istep     = ', istep
   print *, ' t         = ', t   
end if ! node 0

! Go to phi-space to rearrange to phi sequential across processors
Call transpose_z_to_phi(ndof, q, q_phi_space)

Do m = 1, ndof
   Do j = 1, nphi
      Do k = sz_phi, ez_phi
         Do i = sr, er
            ! change from y sequential within each processor to y sequential across processors (after y_to_z)   
            temp_phi(i,k,j) = q_phi_space(i,k,m,1+mod(j-1,mphi)*ng_phi+(j-1)/mphi)    
         End Do
      End Do
   End Do

   Call transpose_phi_to_z(1, temp_phi, temp_z)     ! go to z-space
   Call transpose_z_to_r  (1, temp_z,   temp_r)     ! go to x-space for output

   do j = sphi, ephi
      do k = sz_r, ez_r
         do i = 1, nr
            Qsave(i,k) = temp_r(j,k,i)
         end do
      end do
      call mpi_file_write_ordered(file_handle, qsave, nr*mz_r, mpi_double_precision, status, ier)
      call mpi_barrier(mpi_comm_world, ier)
   end do
end do

Call MPI_FILE_CLOSE(file_handle, ier)
if (my_node .eq. 0) then
   print *, ' called mpi_file_close'
end if

end subroutine mpi_write_save_file_version2
#endif

!----------------------------------------------------------------------------------85

#ifdef mpi_code
subroutine mpi_read_restart_version2(istep, t)

! This routine is from Alan Wray's stellar box code and when used for restart allows:
!    (a) A change in the number of processors (This must have a partitioning solution.
!        If not, the partitioning routine should have signaled an error.)
!    (b) A doubling of nr
!    (c) A doubling of nphi
!    (d) A doubling of the phi domain size, phi_max
!    Note: A doubling of nz is not allowed.

use mpi
use partition_data
use partition_data_for_alan
use grid
use q_array
use transposes_of_q_and_qdot
use basic_state
implicit none
integer, intent(out) :: istep
real(8), intent(out) :: t
integer :: len
character(80) :: string

! Local
integer(kind = mpi_offset_kind)     :: zero_offset = 0
integer, dimension(mpi_status_size) :: status
integer :: file_handle, ier

integer :: num_nodes_read, ndof_read, nr_read, nz_read, nphi_read
real(8) :: rmin_read, rmax_read, zmin_read, zmax_read, phi_max_read
integer :: i, j, k, m, n
integer :: ng_r, ng_phi

! For reference:
! q          (sr:er, sphi:ephi, nz, ndof)           : in z space
! q_r_space  (sphi:ephi, sz_r:ez_r, ndof, nr)
! q_phi_space(sr:er, sz_phi:ez_phi, ndof, nphi)

real(8), dimension(nr, sz_r:ez_r)              :: Qres
real(8), dimension(sr:er,     sphi:ephi, nz)   :: temp_z
real(8), dimension(sphi:ephi, sz_r:ez_r, nr)   :: temp_r
real(8), dimension(sr:er, sz_phi:ez_phi, nphi) :: temp_phi

! This does not have any meaning anymore since we don't store the basic state
! in the save file.
logical :: have_basic_state_read

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' first executable in mpi_read_restart_version2'
   end if
#endif

! Alan's notation:
ng_r   = ng1
ng_phi = ng2

! Open restart file

call mpi_file_open(mpi_comm_world, 'mpi_restart_version2', mpi_mode_rdonly, mpi_info_null, &
     file_handle, ier)

if (ier/=MPI_SUCCESS) then
   write(6,*) "Can't MPI_FILE_OPEN the restart file mpi_restart_version2"
   call MPI_ABORT(MPI_COMM_WORLD, 99, ier)
end if

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' mpi_read_restart_version2: about to call mpi_file_set_view'
   end if
#endif

! Create a "view" for the restart file
call mpi_file_set_view(file_handle, zero_offset, mpi_byte, mpi_byte, "native", mpi_info_null, ier)

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' mpi_read_restart_version2: after mpi_file_set_view'
   end if
#endif

! Only node 0 reads the header portion of the restart file:
if (my_node .eq. 0) then
   ! 13 items in header:
   call mpi_file_read_shared(file_handle, num_nodes_read,        1, mpi_integer, status, ier)
   call mpi_file_read_shared(file_handle, ndof_read,             1, mpi_integer, status, ier)
   call mpi_file_read_shared(file_handle, istep,                 1, mpi_integer, status, ier)
   call mpi_file_read_shared(file_handle, t,                     1, mpi_double,  status, ier)
   call mpi_file_read_shared(file_handle, nr_read,               1, mpi_integer, status, ier)
   call mpi_file_read_shared(file_handle, nz_read,               1, mpi_integer, status, ier)
   call mpi_file_read_shared(file_handle, nphi_read,             1, mpi_integer, status, ier)
   call mpi_file_read_shared(file_handle, rmin_read,             1, mpi_double,  status, ier)
   call mpi_file_read_shared(file_handle, rmax_read,             1, mpi_double,  status, ier)
   call mpi_file_read_shared(file_handle, zmin_read,             1, mpi_double,  status, ier)
   call mpi_file_read_shared(file_handle, zmax_read,             1, mpi_double,  status, ier)
   call mpi_file_read_shared(file_handle, phi_max_read,          1, mpi_double,  status, ier)
   ! This should be removed at some point since we don't anymore store the basic state
   ! in the save file.
   call mpi_file_read_shared(file_handle, have_basic_state_read, 1, mpi_logical, status, ier)
   print *, ' node 0: Read the following from restart file:'
   print *, '    num_nodes read = ', num_nodes_read, ' num_nodes = ', num_nodes
   print *, '    ndof read      = ', ndof_read, ' ndof = ', ndof
   print *, '    istep          = ', istep
   print *, '    t              = ', t
   print *, ' nr read = ', nr_read, ' nz read = ', nz_read, ' nphi read = ', nphi_read
   print *, ' nr      = ', nr,      ' nz      = ', nz,      ' nphi      = ', nphi
   print *, ' rmin read = ', rmin_read, ' rmax read = ', rmax_read, &
            ' zmin read = ', zmin_read, ' zmax read = ', zmax_read, &
            ' phi_max read = ', phi_max_read   
   print *, ' have_basic_state = ', have_basic_state
end if ! node 0

call mpi_bcast(ndof_read,      1, mpi_integer, 0, mpi_comm_world, ier)
call mpi_bcast(nr_read,        1, mpi_integer, 0, mpi_comm_world, ier)
call mpi_bcast(nz_read,        1, mpi_integer, 0, mpi_comm_world, ier)
call mpi_bcast(nphi_read,      1, mpi_integer, 0, mpi_comm_world, ier)
call mpi_bcast(phi_max_read,   1, mpi_double,  0, mpi_comm_world, ier)

! Check to see that the changes have been kept to what is allowed:
if (nz_read .ne. nz) then
   if (my_node .eq. 0) then
      print *, ' subroutine mpi_read_restart_version2'
      print *, ' A change in nz is not allowed'
      print *, ' nz in restart file = ', nz, ' nz = ', nz
   end if
   call terminate_with_no_save(1)
end if

if ((nr_read .ne. nr) .and. (nr .ne. 2*nr_read)) then
   if (my_node .eq. 0) then
      print *, ' subroutine mpi_read_restart_version2'
      print *, ' You can only double nr'
      print *, ' nr in restart file = ', nr_read, ' nr = ', nr
   end if
   call terminate_with_no_save(1)
end if

if ((nphi_read .ne. nphi) .and. (nphi .ne. 2*nphi_read)) then
   if (my_node .eq. 0) then
      print *, ' subroutine mpi_read_restart_version2'
      print *, ' You can only double nphi'
      print *, ' nphi in restart file = ', nphi_read, ' nphi = ', nphi
   end if
   call terminate_with_no_save(1)
end if

! Read the data
! Variables with _r on the end of their names are the restart-file values.
! These may differ in certain ways from the values input for the current run.
! E.g., nr may = 2*nr_r, doubling the mesh in the x direction; or,
! Lphi may be = 2*Lphi_r, meaning that the domain size is doubled in the phi direction.
! ndof_r will differ from ndof if we are restarting from a non-MHD run but are
! doing MHD in the current run.
! Only certain combinations of these changes are allowed, for simplicity and usefulness.
! I did not include here the possible changes in nz and Lz as they are quite messy.

Do m = 1,ndof_read
   Do j = sphi,sphi+mphi-1
      ! a block of nr_r*mz_r, but don't go past the end of the restart file's z's      
      n = nr_read*max(0,min(mz_r,nz_read-mod(my_node,ng_r)*mz_r))
      ! don't read past the end of the restart file's y's; must read even
      ! if count is 0 according to the rules for ORDERED I/O
      If (my_node / ng_r+(j-sphi)*ng_phi >= nphi_read) n = 0
      ! Read the data in node order      
      call MPI_FILE_READ_ORDERED(file_handle, Qres, n, MPI_DOUBLE_PRECISION, status, ier)

      ! What is allowed: (a) Doubing nr (b) Doubling nphi
      Do k = sz_r, ez_r
         If (nr_read == nr) then
            Do i = 1, nr
               temp_r(j,k,i) = Qres(i,k)   
            End Do
         !Else if (Lx==2*Lx_r) then
         !   doubling the r-domain at the same resolution.
         !   I don't allow this because my "r" is never periodic.
         !   Do i = 1,nr_r
         !      temp_r(j,k,i)      = Qres(i,k)
         !      temp_r(j,k,i+nr_r) = Qres(i,k)
         !   End Do
         Else  ! linearly interpolate onto doubled r-mesh
            Do i = 1, nr_read - 1
               temp_r(j,k,2*i-1) = Qres(i,k)
               temp_r(j,k,2*i  ) = 0.5d0*(Qres(i,k) + Qres(i+1,k))
            End Do
            temp_r(j,k,2*nr_read-1) = Qres(nr_read,k)
            temp_r(j,k,2*nr_read  ) = 0.5d0*(Qres(nr_read,k) + Qres(1,k))
         End If
      End Do
   End Do
#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' node 0, mpi_read_restart_version2: calling transpose_r_to_z for idof = ', m
   end if
#endif
   
   Call transpose_r_to_z  (1, temp_r, temp_z  )     ! go to z-space

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' node 0, mpi_read_restart_version2: returned from transpose_r_to_z for idof = ', m
   end if
#endif
   
   Call transpose_z_to_phi(1, temp_z, temp_phi)     ! go to phi-space

#ifdef debug_print
   !if (my_node .eq. 0) then
      print *, ' node = ', my_node, ' mpi_read_restart_version2: returned from transpose_z_to_phi for idof = ', m
   !end if
#endif   
   
   Do j = 1, nphi
      Do k = sz_phi, ez_phi
         Do i = sr, er
            ! change from y sequential across processors to y sequential within each processor            
            q_phi_space(i,k,m,1+mod(j-1,mphi)*ng_phi+(j-1)/mphi) = temp_phi(i,k,j)    
         End Do
      End Do
   End Do
   If (nphi == 2*nphi_read) then
      If (phi_max == 2*phi_max_read) then      ! doubling the phi-domain at the same resolution
         Do j = 1, nphi_read
            Do k = sz_phi, ez_phi
               Do i = sr, er
                  temp_phi(i,k,j          ) = q_phi_space(i,k,m,j)
                  temp_phi(i,k,j+nphi_read) = q_phi_space(i,k,m,j)
               End Do
            End Do
         End Do
      Else  ! linearly interpolate onto doubled phi mesh
         Do j = 1, nphi_read - 1
            Do k = sz_phi, ez_phi
               Do i = sr, er
                  temp_phi(i,k,2*j-1) = q_phi_space(i,k,m,j)
                  temp_phi(i,k,2*j  ) = 0.5d0*(q_phi_space(i,k,m,j) + q_phi_space(i,k,m,j+1))
               End Do
            End Do
         End Do
         Do k = sz_phi, ez_phi
            Do i = sr, er
               temp_phi(i,k,2*nphi_read-1) = q_phi_space(i,k,m,nphi_read)
               temp_phi(i,k,2*nphi_read  ) = 0.5d0*(q_phi_space(i,k,m,nphi_read) + q_phi_space(i,k,m,1))
            End Do
         End Do
      End If

      ! q_phi_space(sr:er, sz_phi:ez_phi, ndof, nphi)     
      Do j = 1,nphi
         Do k = sz_phi, ez_phi
            Do i = sr, er
               q_phi_space(i,k,m,j) = temp_phi(i,k,j)
            End Do
         End Do
      End Do
   End If
End Do

#ifdef debug_print
!   if (my_node .eq. 0) then
      print *, ' node = ', my_node, ' mpi_read_restart_version2: calling transpose_phi_to_z'
!   end if
#endif

call transpose_phi_to_z(ndof, q_phi_space, q)

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' node 0, mpi_read_restart_version2: returned from transpose_phi_to_z'
   end if
#endif

call mpi_bcast(istep,            1, mpi_integer, 0, mpi_comm_world, ier)
call mpi_bcast(t,                1, mpi_double,  0, mpi_comm_world, ier)
call mpi_bcast(have_basic_state, 1, mpi_logical, 0, mpi_comm_world, ier)

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' node 0, mpi_read_restart_version2: returning'
   end if
#endif

end subroutine mpi_read_restart_version2
#endif

!----------------------------------------------------------------------------------85



