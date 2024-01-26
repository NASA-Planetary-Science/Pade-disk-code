!----------------------------------------------------------------------------------85

! This is needed for FFTW:
use, intrinsic :: iso_c_binding
implicit none

! This is needed for FFTW:
include 'fftw3.f03'

character(80) filename
integer, parameter :: rcomp=1, zcomp=2, pcomp=3
integer, parameter :: lun_list=1, lun_tecplot=2, lun_spectrum=3, lun_hanning=4
real(8) :: time ! Read from header in tecplot file
! For test:
real(8) :: r_prime
! Temporaries for calculating Hanning window:
real(8) :: lambda, rmid, xi, pi
integer :: nr, nphi, nphi_plus_1, ifile, nfiles, ir, iphi, icomp, N, kr

! x, y grid read from tecplot file:
real(8), allocatable, dimension(:,:) :: xgrid, ygrid

real(8), allocatable, dimension(:) :: rgrid
! The vorticity components are all fluctuations from the basic state.
! The three indices are (ir, iphi, icomponent):
real(8), allocatable, dimension(:, :, :) :: omega
character(1) :: dum_char

! For the FFTW library:
integer(8) :: plan

complex(8), allocatable, dimension(:) :: fhat
real(8),    allocatable, dimension(:) :: f, sum_over_components, average_over_phi_and_files, hanning
character(80) :: name_of_list_file

pi = 4.0d0 * atan(1.d0)

write(6, "(' enter name of file containing list of files to process--->', $)")
read(5, *) name_of_list_file

print *, ' Attempting to open list of tecplot files to process'
open(unit = lun_list, file = name_of_list_file, &
     form = 'formatted', status = 'unknown')

read(lun_list, *) nfiles
print *, ' nfiles to process = ', nfiles
print *, ' Enter anything to continue'
read(5, *)

do ifile = 1, nfiles
   read(lun_list, "(a46)") filename
   print *, ' Attempting to open file = ', filename, ' File # = ', ifile
   open(unit = lun_tecplot, file = filename, form = 'formatted', status = 'old')
   print *, ' Successful'

   read(lun_tecplot, 1) time
   1 format (13x, e12.5)
   print *, ' t = ', time
   read(lun_tecplot, 2) dum_char
   2 format(a1)
   ! The plus 1 accounts accounts for the fact that we complete periodicity in the tecplot file:
   read(lun_tecplot, 3) nr, nphi_plus_1
   3 format(7x, i4, 4x, i4)

   nphi = nphi_plus_1 - 1
   print *, ' nr = ', nr, ' nphi = ', nphi

   if (ifile .eq. 1) then
      allocate(xgrid(nr,nphi_plus_1), ygrid(nr,nphi_plus_1), omega(nr,nphi_plus_1,3), rgrid(nr))

      ! For the FFT:
      N = nr - 1 ! Since the FFT counts from 0:
      allocate(f(N))
      allocate(fhat                      (0:N/2))
      allocate(sum_over_components       (0:N/2))
      allocate(average_over_phi_and_files(0:N/2))
      average_over_phi_and_files = 0.d0

      ! Initialize the FFT for real to complex:
      print *, ' calling dfftw_plan_dft_r2c_1d'
      call dfftw_plan_dft_r2c_1d(plan, N, f, fhat, FFTW_ESTIMATE)
      print *, ' returned from dfftw_plan_dft_r2c_1d'
   end if

   ! Read vorticity field in the horizontal plane:
   print *, ' reading vorticity field'
   do iphi = 1, nphi+1
      do ir = 1, nr
         read(lun_tecplot, "(5(1x, e12.5))") xgrid(ir,iphi), ygrid(ir,iphi), &
              omega(ir,iphi,rcomp), omega(ir,iphi,zcomp), omega(ir,iphi,pcomp)
      end do
   end do
   close(lun_tecplot)
   print *, ' Finished reading vorticity field'

   ! Calculate rgrid:
   if (ifile .eq. 1) then
      iphi = 1
      do ir = 1, nr
         rgrid(ir) = sqrt(xgrid(ir,iphi)**2 + ygrid(ir,iphi)**2)
      end do

      ! Hanning window:
      allocate(hanning(nr))
      open(unit = lun_hanning, file = 'hanning.dat', form = 'formatted', status = 'unknown')
      do ir = 1, N
         lambda = rgrid(nr) - rgrid(1)
         rmid = rgrid(1) + 0.5d0*lambda
         xi = rgrid(ir) - rmid
         hanning(ir) = cos(pi*xi/lambda)**2
         write(lun_hanning, "(2(1x, e12.5))")
      end do
      close(lun_hanning)
   end if

   do iphi = 1, nphi
      sum_over_components = 0.d0
      do icomp = 1, 3
         do ir = 1, n
            f(ir) = omega(ir,iphi,icomp)*hanning(ir)      
         end do
         call dfftw_execute(plan)
         fhat(0) = fhat(0)/N ! Normalization
         do kr = 1, N/2
            fhat(kr) = fhat(kr)*2.d0/N
         end do
   
         do kr = 0, N/2-1
            sum_over_components(kr) = sum_over_components(kr) + fhat(kr)*conjg(fhat(kr)) 
         end do
      end do ! sum over components
      average_over_phi_and_files = average_over_phi_and_files + sum_over_components
   end do ! iphi
end do ! loop over files

do kr = 0, N/2 - 1
   average_over_phi_and_files(kr) = average_over_phi_and_files(kr) / (nphi*nfiles)
end do
   
! Test:
!lambda = rgrid(nr) - rgrid(1)
!k      = 2.d0*pi/lambda
!print *, ' lambda = ',lambda
!do ir = 1, N
!   r_prime = rgrid(ir) - rgrid(1)
!   f(ir) = 1.d0 + sin(k*r_prime) + cos(k*r_prime)
!end do

open(unit = lun_spectrum, file = 'enstrophy_spectrum.dat', form = 'formatted', status = 'unknown')
do kr = 1, N/2-1
   write(lun_spectrum, "(e12.5, 1x, e12.5)") kr*2.d0*pi/lambda, average_over_phi_and_files(kr)
end do
close(lun_spectrum)

call dfftw_destroy_plan(plan)

end

!----------------------------------------------------------------------------------85

