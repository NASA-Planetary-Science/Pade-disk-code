use HDF5
implicit none
include 'mpif.h'
integer, parameter :: nx = 512, ny = 512
integer :: ier

integer(HID_T) :: file_id, dspace_id
integer(HID_T) :: dset_id_for_x, dset_id_for_y, dset_id_for_f
integer(HID_T) :: fa_plist_id   ! file acccess  property list id
integer(HID_T) :: dx_plist_id   ! data transfer property list id

integer, PARAMETER :: rank = 2
integer(HSIZE_T), DIMENSION(2) :: dims_ds, dims_slab, offset, count
REAL(8), allocatable, DIMENSION(:, :) :: x, y, f
integer :: ix, iy, my, sy, ey
real(8) :: dx, dy
real(8) :: pi

integer :: info
integer :: mpi_size, mpi_rank

info = MPI_INFO_NULL
call MPI_INIT(ier)
call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ier)
call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ier)

my = ny/mpi_size     ! number of points each processor has in the direction
sy = mpi_rank*my + 1 ! starting index in x for this processor 
ey = (mpi_rank+1)*my ! ending   index in x for this processor

allocate(x(nx, sy:ey), y(nx, sy:ey), f(nx, sy:ey))
count(1) = nx
count(2) = my

! Create the function which is to be our test dataset:
pi = 4.0d0*atan(1.d0)
dx = 2.d0*pi / (nx - 1)
dy = 2.d0*pi / (ny - 1)

! Generate the function portion contained in this processor.
do iy = sy, ey
   do ix = 1, nx
      x(ix, iy) = (ix - 1)*dx
      y(ix, iy) = (iy - 1)*dy         
      f(ix, iy) = sin(x(iy,ix))
   end do
end do
   
! Initialize FORTRAN interface of HDF5.
call h5open_f(ier)
if (ier .ne. 0) then
   print *, ' h5open_f returned ier = ', ier
   stop
end if

! Set up file access property list w/parallel I/O access
call h5Pcreate_f(H5P_FILE_ACCESS_F, fa_plist_id, ier)
if (ier .ne. 0) then
   print *, ' h5Pcreate_f returned ier = ', ier
   stop
end if
! Stores MPI IO communicator information to the file access property list. 
call h5Pset_fapl_mpio_f(fa_plist_id, MPI_COMM_WORLD, info, ier)
if (ier .ne. 0) then
   print *, ' h5Pset_fapl_mpio_f returned ier = ', ier
   stop
end if

! Create the file:
call h5Fcreate_f("func.h5", H5F_ACC_TRUNC_F, file_id, ier, access_prp = fa_plist_id)
if (ier .ne. 0) then
   print *, ' h5Fcreate_f returned ier = ', ier
   stop
end if

! Done with fa_plist_id:
call h5Pclose_f(fa_plist_id, ier)
if (ier .ne. 0) then
   print *, ' h5Pclose_f returned ier = ', ier
   stop
end if

! Create the data space for the whole dataset (covering all the processors):
dims_ds(1) = nx
dims_ds(2) = ny
call h5Screate_simple_f(rank, dims_ds, dspace_id, ier)
if (ier .ne. 0) then
   print *, ' h5Screate_simple_f returned ier = ', ier
   stop
end if

! Create the datasets with default properties.
call h5Dcreate_f(file_id, "x", H5T_NATIVE_DOUBLE, dspace_id, dset_id_for_x, ier)
if (ier .ne. 0) then
   print *, ' h5Dcreate_f returned ier = ', ier
   stop
end if
call h5Dcreate_f(file_id, "y", H5T_NATIVE_DOUBLE, dspace_id, dset_id_for_y, ier)
if (ier .ne. 0) then
   print *, ' h5Dcreate_f returned ier = ', ier
   stop
end if
call h5Dcreate_f(file_id, "f", H5T_NATIVE_DOUBLE, dspace_id, dset_id_for_f, ier)
if (ier .ne. 0) then
   print *, ' h5Dcreate_f returned ier = ', ier
   stop
end if

! Set up data transfer property list w/collective MPI-IO
call h5Pcreate_f(H5P_DATASET_XFER_F, dx_plist_id, ier)
if (ier .ne. 0) then
   print *, ' h5Pcreate_f returned ier = ', ier
   stop
end if
! Sets data transfer mode (collective versus independent): 
call h5Pset_dxpl_mpio_f(dx_plist_id, H5FD_MPIO_COLLECTIVE_F, ier)
if (ier .ne. 0) then
   print *, ' h5Pset_dxpl_mpio_f returned ier = ', ier
   stop
end if

! Select hyperslab within the full dataspace:
offset(1) = 0                 ! hdf5 counts from 0
offset(2) = mpi_rank*count(2) ! hdf5 counts from 0
count(1) = nx
count(2) = my
CALL h5Sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, ier)
if (ier .ne. 0) then
   print *, ' h5Sselect_hyperslab_f returned ier = ', ier
   stop
end if

! Write the datasets collectively.
dims_slab(1) = nx
dims_slab(1) = my
call h5Dwrite_f(dset_id_for_x, H5T_NATIVE_DOUBLE, x, dims_slab, ier, xfer_prp = dx_plist_id)
if (ier .ne. 0) then
   print *, ' h5Dwrite_f returned ier = ', ier
   stop
end if
call h5Dwrite_f(dset_id_for_y, H5T_NATIVE_DOUBLE, y, dims_slab, ier, xfer_prp = dx_plist_id)
if (ier .ne. 0) then
   print *, ' h5Dwrite_f returned ier = ', ier
   stop
end if
call h5Dwrite_f(dset_id_for_f, H5T_NATIVE_DOUBLE, f, dims_slab, ier, xfer_prp = dx_plist_id)
if (ier .ne. 0) then
   print *, ' h5Dwrite_f returned ier = ', ier
   stop
end if

! Deallocate data buffers.
DEALLOCATE(x, y, f)

! Close resources.
call h5Dclose_f(dset_id_for_x, ier)
call h5Dclose_f(dset_id_for_y, ier)
call h5Dclose_f(dset_id_for_f, ier)
call h5Sclose_f(dspace_id,     ier)
call h5Pclose_f(dx_plist_id,   ier)
call h5Fclose_f(file_id,       ier)
! Attempt to remove the data file.  Remove the line if the compiler
! does not support it.
!call unlink(filename)

call h5close_f(ier)
call MPI_FINALIZE(ier)

end


   
