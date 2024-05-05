!**********************************************************************************85

! Performs a time average of profile files.
! The program will ask you for a file containing a list of profile files.
! The list file should consist of

! # of files    # points in each profile file
! file1
! file2
! etc.

implicit none

real(8), dimension(:), allocatable :: x, v, ave
integer :: i, npnts, ifile, nfiles, nchar
character(40) :: list_file, output_file, profile_file

integer, parameter :: lun_list = 1, lun_profile = 2, lunout = 3

write(6, "(' enter name of file containing list of profile files to average--->', $)")
read (5, *) list_file

open(unit = lun_list, file = list_file, form = 'formatted', status = 'unknown')
read(lun_list, *) nfiles, npnts
print *, ' nfiles = ', nfiles, ' nchar = ', nchar, ' npnts = ', npnts
allocate(x(npnts), v(npnts), ave(npnts))
allocate(profile_file_array(nchar))

ave = 0.d0
do ifile = 1, nfiles
   read(lun_list, *) profile_file
   print *, ' opening profile file = ', profile_file
   open(unit = lun_profile, file = profile_file, form = 'formatted', status = 'unknown')

   do i = 1, npnts
      read(lun_profile, *) x(i), v(i)
      ave(i) = ave(i) + v(i)
   end do
   close(lun_profile)
end do

write(6, "(' enter name of output file containing the time average--->', $)")
read (5, *) output_file
open (unit = lunout, file = output_file, form = 'formatted', status = 'unknown')
do i = 1, npnts
   write(lunout, "(2(1x, e12.5))") x(i), ave(i)/nfiles
end do
close (lunout)

end

!**********************************************************************************85
