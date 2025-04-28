!
! Number of processes is assumed to be 1 or multiples of 2 (1,2,4,6,8)
!

     PROGRAM DATASET_BY_COL

     USE HDF5 ! This module contains all necessary modules 
        
     IMPLICIT NONE

     include 'mpif.h'
     integer, parameter :: nx = 512, ny = 512
     integer :: ix, iy, my, sy, ey
     real(8) :: dx, dy
     real(8) :: pi
     CHARACTER(LEN=10), PARAMETER :: filename = "sds_col.h5"  ! File name
     CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name

     INTEGER(HID_T) :: file_id       ! File identifier 
     INTEGER(HID_T) :: dset_id_for_x, dset_id_for_y, dset_id_for_f  ! Dataset identifiers
     INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
     INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
     INTEGER(HID_T) :: plist_id      ! Property list identifier 

     INTEGER(HSIZE_T), DIMENSION(2) :: dimsf = (/nx,ny/) ! Dataset dimensions.
     INTEGER(HSIZE_T), DIMENSION(2) :: dimsfi = (/nx,ny/)

     INTEGER(HSIZE_T), DIMENSION(2) :: count  
     INTEGER(HSSIZE_T), DIMENSION(2) :: offset 
     real(8), ALLOCATABLE, dimension(:,:) :: x, y, f  ! Data to write
     INTEGER :: rank = 2 ! Dataset rank 

     INTEGER :: error, error_n  ! Error flags

     ! For MPI definitions and calls.
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank

     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL
     CALL MPI_INIT(mpierror)
     CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
     CALL MPI_COMM_RANK(comm, mpi_rank, mpierror)

     pi = 4.d0 * atan(1.0d0)

     ! Initialize data buffer with data.
     my = ny / mpi_size     
     sy = mpi_rank*my + 1
     ey = (mpi_rank + 1)*my
     ALLOCATE(x(nx, sy:ey), y(nx, sy:ey), f(nx, sy:ey))
     dx = 2.d0*pi / (nx - 1)
     dy = 2.d0*pi / (ny - 1)
     do iy = sy, ey
        do ix = 1, ny
           x(ix, iy) = (ix - 1)*dx
           y(ix, iy) = (iy - 1)*dy
           f(ix, iy) = sin(3.d0*x(ix, iy))*cos(y(ix,iy))
        end do
     end do
     
     ! Initialize FORTRAN predefined datatypes
     CALL h5open_f(error) 

     ! Setup file access property list with parallel I/O access.
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)

     ! Create the file collectively. 
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
     CALL h5pclose_f(plist_id, error)

     ! Create the data space for the dataset. 
     CALL h5screate_simple_f(rank, dimsf, filespace, error)

     ! Create the dataset with default properties.
     CALL h5dcreate_f(file_id, 'x', H5T_NATIVE_DOUBLE, filespace, dset_id_for_x, error)
     CALL h5dcreate_f(file_id, 'y', H5T_NATIVE_DOUBLE, filespace, dset_id_for_y, error)
     CALL h5dcreate_f(file_id, 'f', H5T_NATIVE_DOUBLE, filespace, dset_id_for_f, error)     
     CALL h5sclose_f(filespace, error)
     !
     ! Each process defines dataset in memory and writes it to the hyperslab
     ! in the file. 
     !
     dimsf(1) = nx
     dimsf(2) = my
     
     count(1) = nx
     count(2) = my 
     offset(1) = 0
     offset(2) = mpi_rank * my 
     CALL h5screate_simple_f(rank, count, memspace, error) 

     ! Select hyperslab in the file.
     CALL h5dget_space_f(dset_id_for_x, filespace, error)
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)


     ! Create property list for collective dataset write
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     
     ! Write the dataset collectively. 
     CALL h5dwrite_f(dset_id_for_x, H5T_NATIVE_DOUBLE, x, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
     CALL h5dwrite_f(dset_id_for_y, H5T_NATIVE_DOUBLE, y, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
     CALL h5dwrite_f(dset_id_for_f, H5T_NATIVE_DOUBLE, f, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

     DEALLOCATE(x, y, f)

     ! Close dataspaces.
     CALL h5sclose_f(filespace, error)
     CALL h5sclose_f(memspace, error)

     ! Close the dataset and property list.
     CALL h5dclose_f(dset_id_for_x, error)
     CALL h5dclose_f(dset_id_for_y, error)
     CALL h5dclose_f(dset_id_for_f, error)     
     CALL h5pclose_f(plist_id, error)

     ! Close the file.
     CALL h5fclose_f(file_id, error)

     ! Close FORTRAN predefined datatypes.
     CALL h5close_f(error)

     CALL MPI_FINALIZE(mpierror)

     END PROGRAM DATASET_BY_COL
