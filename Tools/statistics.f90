!**********************************************************************************85

! Computes phi and time averaged statistics.

#ifdef mpi_code
   use mpi
#endif
implicit none
integer :: nfiles, nr, nz, nphi, ier, ifile
character(40), allocatable, dimension(:) :: filename

#ifdef mpi_code
   call mpi_barrier(mpi_comm_world, ier)
#endif

! Read namelist input:
if (my_node .eq. 0) then
   print *, ' node 0: about to open and read namelist file for app_statistics'      
   open (unit = lun_general_purpose, file = 'list_of_files', form = 'formatted', &
        status = 'old')
   read(lun_general_purpose, *) nfiles
   allocate(filename(nfiles)
   do ifile = 1, nfiles
      read(lun_general_purpose, *) filename(ifile)
   end do   
   close (lun_general_purpose)
end if

#ifdef mpi_code
   call mpi_bcast(nfiles, 1, mpi_integer, 0, mpi_comm_world, ier)

#endif

end program

!**********************************************************************************85
