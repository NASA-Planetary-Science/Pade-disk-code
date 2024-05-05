!**********************************************************************************85

! This program reads files of the type: phi_z_averages_000011.6485_0024004.dat
! (for instance) and computes phi-z-time averages.

! It was written to post-process Taylor-Couette data.

implicit none
character(80) stats_file
integer :: nfiles, ifile, ir, nr, nz, ivar
integer, parameter :: lun_list_file = 1, lun_Reynolds = 2

! Indices.  These follows the ordering of the quantities in the
! lun_Raynolds file:
integer, parameter :: rho = 1, rho_uphi = 2, rho_ur = 3, rho_uz = 4, &
     rho_uphi_ur = 5, rho_uphi2 = 6, rho_ur2 = 7, &
     rho_uz2 = 8, rho_uphi_uz = 9, rho_ur_uz = 10

real(8) :: value_read

! phi-z-time averages:
real(8), allocatable, dimension(:,:) :: bar

real(8), allocatable, dimension(:) :: rgrid

! Final phi-z-time averages we will output:
real(8), allocatable, dimension(:) :: rho_bar, ur_tilde, uz_tilde, uphi_tilde

! Final stresses of phi-time average of Favre fluctuations we will output
real(8), allocatable, dimension(:, :) :: T

! The above divided by rho (incompressible type stresses):
real(8), allocatable, dimension(:, :) :: F

character(80) :: list_file
character(10) :: var_name

! Last one is the trace:
integer, parameter :: pr = 1, pp = 2, rr = 3, zz = 4, pz = 5, rz = 6, kk = 7

write(6, "(' enter name of file containing list of files to process--->', $)")
read(5, *) list_file

print *, ' Attempting to open file containing list of files'
open(unit = lun_list_file, file = list_file, &
     form = 'formatted', status = 'unknown')

! First add up the Reynolds averages at different times:
do ifile = 1, 1000000
   read(lun_list_file, "(a80)", end = 999) stats_file
   print *, ' Attempting to open file = ', stats_file
   open(unit = lun_Reynolds, file = stats_file, form = 'formatted', status = 'old')
   print *, ' Successful'
   read(lun_Reynolds, *) nr
   print *, ' nr = ', nr

   if (ifile .eq. 1) then
      ! This will be read:
      allocate(bar(nr, 10))

      ! Final output:
      allocate(rgrid(nr), rho_bar(nr), ur_tilde(nr), uz_tilde(nr), uphi_tilde(nr))
      allocate(T(nr, 7))

      ! The above divided by <rho>:
      allocate(F(nr, 7))

      bar = 0.d0
   end if

   do ivar = 1, 10
      read(lun_Reynolds, *) var_name
      print *, ' Reading variable = ', var_name
      do ir = 1, nr
         read(lun_Reynolds, *) rgrid(ir), value_read
         print *, ' ir = ', ir, ' r = ', rgrid(ir), ' value = ', value_read
         bar(ir, ivar) = bar(ir, ivar) + value_read
      end do
   end do
   close(lun_Reynolds)
end do
close(lun_list_file)

999 continue
nfiles = ifile - 1
print *, ' done with reading files.  nfiles = ', nfiles

! Done with reading all the files and doing sums.

! Make the sums into averages:
bar = bar / nfiles

! Now compute Favre averages (see notes of 11/24/20):

! phi-z-time Favre averages
do ir = 1, nr
   rho_bar   (ir) = bar(ir, rho)
   ur_tilde  (ir) = bar(ir, rho_ur  ) / bar(ir, rho)
   uz_tilde  (ir) = bar(ir, rho_uz  ) / bar(ir, rho)
   uphi_tilde(ir) = bar(ir, rho_uphi) / bar(ir, rho)
end do

open(unit = 1, file = 'u_phi_Favre_mean.dat', form = 'formatted', status = 'unknown')
open(unit = 2, file = 'u_r_Favre_mean.dat',   form = 'formatted', status = 'unknown')
open(unit = 3, file = 'u_z_Favre_mean.dat',   form = 'formatted', status = 'unknown')
do ir = 1, nr
   write(1, "(2(1x, e12.5))") rgrid(ir), uphi_tilde(ir)
   write(2, "(2(1x, e12.5))") rgrid(ir), ur_tilde(ir)
   write(3, "(2(1x, e12.5))") rgrid(ir), uz_tilde(ir)   
end do
close(1); close(2); close(3)

! These are <rho u_i" u_j">:
do ir = 1, nr
   T(ir, pr) = bar(ir, rho_uphi_ur) - bar(ir, rho_ur  )*bar(ir, rho_uphi)/bar(ir,rho)
   T(ir, pp) = bar(ir, rho_uphi2  ) - bar(ir, rho_uphi)*bar(ir, rho_uphi)/bar(ir,rho)
   T(ir, rr) = bar(ir, rho_ur2    ) - bar(ir, rho_ur  )*bar(ir, rho_ur  )/bar(ir,rho)
   T(ir, zz) = bar(ir, rho_uz2    ) - bar(ir, rho_uz  )*bar(ir, rho_uz  )/bar(ir,rho)
   T(ir, pz) = bar(ir, rho_uphi_uz) - bar(ir, rho_uphi)*bar(ir, rho_uz  )/bar(ir,rho)
   T(ir, rz) = bar(ir, rho_ur_uz  ) - bar(ir, rho_ur  )*bar(ir, rho_uz  )/bar(ir,rho)
   T(ir, kk) = T(ir, rr) + T(ir, zz) + T(ir, pp)
end do

! Divide the above by bar(ir, rho) to get the incompressible type stresses:
do ivar = 1, 7
   do ir = 1, nr
      F(ir, ivar) = T(ir, pr) / bar(ir, rho)
   end do
end do

open(unit = 1, file = 'Favre_uphi_p_rms.dat', form = 'formatted', status = 'unknown')
do ir = 1, nr
   write(1, "(2(1x, e12.5))") rgrid(ir), sqrt(F(ir, pp))
end do
close(1)

end program

!**********************************************************************************85

