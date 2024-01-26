!----------------------------------------------------------------------------------85

! This is needed for FFTW:
use, intrinsic :: iso_c_binding
implicit none

! This is needed for FFTW:
include 'fftw3.f03'

character(80) filename
integer, parameter :: lun_tecplot = 1, lun_list = 2, lun_spec = 3
integer, parameter :: max_files = 1000
integer :: ifile, nfiles, k
real(8) :: t, r_prime, lambda, pi, rmid, x, kH0
integer :: nr, nphi, ir, iphi, irec, N, ir_min, ir_max
real(8), allocatable, dimension(:,:) :: xgrid, ygrid
real(8), allocatable, dimension(:)   :: rgrid
! The velocities are all fluctuations from the basic state:
real(8), allocatable, dimension(:, :) :: ur, uz, uphi_prime

! Irrelevant variables for the velocity spectrum:
real(8) :: rho, rho_prime, uphi, ci

integer :: include_uz, include_ur, include_uphi
character(1) :: dum_char, remove_ave, r_or_p
real(8) :: ur_ave, uz_ave, up_ave

! For the FFTW library:
integer(8) :: plan
real(8),    allocatable, dimension(:) :: u
complex(8), allocatable, dimension(:) :: uhat
real(8), allocatable, dimension(:) :: E, Hann
character(80) :: line

pi = 4.0d0 * atan(1.d0)

write(6, "(' Spectrum wrt r or phi (r,p)--->', $)")
read(5, "(a1)") r_or_p

write(6, "(' enter 3 flags (0/1) for including uz, ur, and uphi--->', $)")
read(5, *) include_uz, include_ur, include_uphi

open(unit = lun_list, file = 'list.dat', form = 'formatted', status = 'old')
nfiles = 0
do ifile = 1, max_files
   read(lun_list, "(a80)", end = 999) filename
   nfiles = nfiles + 1
   print *, ' opening ', filename, ' filename = ', filename, ' nfiles = ', nfiles

   open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
         access = 'direct', recl = 13*9 + 1)

   read(lun_tecplot, 1, rec = 1) t
   1 format (14x, e12.5)
   print *, ' t = ', t
   read(lun_tecplot, 2, rec = 2) dum_char
   2 format(a1)
   print *, ' about to read nr and nphi'
   read(lun_tecplot, 3, rec = 3) line
   3  format(a80)
   print *, ' line = ', line
   read(lun_tecplot, 4, rec = 3) nr, nphi
   4 format(8x, i4, 4x, i4)
   print *, ' nr = ', nr, ' nphi = ', nphi
   print *, ' enter anything to continue'
   read(5, *)
   ! nr = 512
   ! nphi = 1025
   
   ! Remove the data at 2*pi which completes the circle:
   nphi = nphi - 1

   if (ifile .eq. 1) then
      allocate(xgrid(nr,nphi), ygrid(nr,nphi), ur(nr,nphi), uz(nr,nphi), uphi_prime(nr,nphi))
      allocate(rgrid(nr), hann(nr))
   end if
   
   ! Read horizontal plane of data:
   do iphi = 1, nphi
      do ir = 1, nr
         irec = (iphi - 1)*nr + ir + 3
         !read(lun_tecplot, "(9(1x, e12.5))",rec=irec) xgrid(ir,iphi),ygrid(ir,iphi),&
         !     rho, rho_prime, ur(ir,iphi), uphi, uphi_prime(ir,iphi), uz(ir,iphi), ci

         read(lun_tecplot, "(9(1x, e12.5))",rec=irec) xgrid(ir,iphi),ygrid(ir,iphi),&
              rho, rho_prime, ur(ir,iphi), uphi, uphi_prime(ir,iphi), uz(ir,iphi), ci         
      end do
   end do
   close(lun_tecplot)

   ! Calculate rgrid:
   do ir = 1, nr
      rgrid(ir) = sqrt(xgrid(ir,1)**2 + ygrid(ir,2)**2)
   end do

   ! Determine the region that is free of sponge:
   do ir = 1, nr
      if (rgrid(ir) - rgrid(1) .gt. 0.5d0) then
         ir_min = ir
         exit
      end if
   end do

   do ir = nr, 1, -1
      if (rgrid(nr) - rgrid(ir) .gt. 0.5d0) then
         ir_max = ir
         exit
      end if
   end do

   print *, ' ir_min = ', ir_min, ' ir_max = ', ir_max
   print *, ' enter anything to continue'
   read(5, *)

   if (r_or_p .eq. 'r') then
      ! For the FFT:
      N = nr
      open(unit = 1, file = 'hann.dat', form = 'formatted', status = 'unknown')
      do ir = 1, N
         lambda = rgrid(nr) - rgrid(1)
         rmid = rgrid(1) + 0.5d0*lambda
         x = rgrid(ir) - rmid
         hann(ir) = cos(pi*x/lambda)**2
      end do
      close(1)
   else
      N = nphi
   end if

   if (ifile .eq. 1) then
      allocate(u(N))
      allocate(uhat(0: N/2))
      allocate(E(0:N/2))
      E = 0.d0
      ! Initialize the FFT:
      call dfftw_plan_dft_r2c_1d(plan, N, u, uhat, FFTW_ESTIMATE)      
   end if

   if (r_or_p .eq. 'p') then
      ! Spectrum wrt phi:
      ! ~~~~~~~~~~~~~~~~~
      do ir = ir_min, ir_max
         if (include_ur .eq. 1) then
            do iphi = 1, N
               u(iphi) = ur(ir, iphi)
            end do
      
            call dfftw_execute(plan)
            uhat(0) = uhat(0)/N
            do iphi = 1, N/2
               uhat(iphi) = uhat(iphi)*2.d0/N
            end do
   
            do iphi = 0, N/2-1
               E(iphi) = E(iphi) + uhat(iphi)*conjg(uhat(iphi)) 
            end do
         end if

         if (include_uz .eq. 1) then
            do iphi = 1, N
               u(iphi) = uz(ir, iphi)
            end do
            call dfftw_execute(plan)
            uhat(0) = uhat(0)/N
            do iphi = 1, N/2
               uhat(iphi) = uhat(iphi)*2.d0/N
            end do
   
            do iphi = 0, N/2-1
               E(iphi) = E(iphi) + uhat(iphi)*conjg(uhat(iphi)) 
            end do
         end if

         if (include_uphi .eq. 1) then
            do iphi = 1, N
               u(iphi) = uphi_prime(ir, iphi)      
            end do
            call dfftw_execute(plan)
            uhat(0) = uhat(0)/N
            do iphi = 1, N/2
               uhat(iphi) = uhat(iphi)*2.d0/N
            end do
   
            do iphi = 0, N/2-1
               E(iphi) = E(iphi) + uhat(iphi)*conjg(uhat(iphi)) 
            end do
         end if
      end do ! ir
   else if (r_or_p .eq. 'r') then
      ! Spectrum wrt r:
      ! ~~~~~~~~~~~~~~~
      do iphi = 1, nphi
         if (include_ur .eq. 1) then
            do ir = 1, N
               u(ir) = ur(ir, iphi)*hann(ir)
            end do
      
            call dfftw_execute(plan)
            uhat(0) = uhat(0)/N
            do ir = 1, N/2
               uhat(ir) = uhat(ir)*2.d0/N
            end do
   
            do ir = 0, N/2-1
               E(ir) = E(ir) + uhat(ir)*conjg(uhat(ir)) 
            end do
         end if

         if (include_uz .eq. 1) then
            do ir = 1, N
               u(ir) = uz(ir, iphi)*hann(ir)
            end do
            call dfftw_execute(plan)
            uhat(0) = uhat(0)/N
            do ir = 1, N/2
               uhat(ir) = uhat(ir)*2.d0/N
            end do
   
            do ir = 0, N/2-1
               E(ir) = E(ir) + uhat(ir)*conjg(uhat(ir)) 
            end do
         end if

         if (include_uphi .eq. 1) then
            do ir = 1, N
               u(ir) = uphi_prime(ir, iphi)*hann(ir)      
            end do
            call dfftw_execute(plan)
            uhat(0) = uhat(0)/N
            do ir = 1, N/2
               uhat(ir) = uhat(ir)*2.d0/N
            end do
   
            do ir = 0, N/2-1
              E(ir) = E(ir) + uhat(ir)*conjg(uhat(ir)) 
            end do
         end if
      end do
   end if
end do ! ifile

print *, ' max_files exceeded'
stop

999 continue
close(lun_list)

if (r_or_p .eq. 'p') then
   do iphi = 0, N/2 - 1
      E(iphi) = E(iphi) / (ir_max - ir_min + 1) / nfiles
   end do
else
   do ir = 0, N/2 - 1
      E(ir) = E(ir) / nphi / nfiles
   end do
end if

! Write spectrum:
open(unit = lun_spec, file = 'vel_spectrum.dat', form = 'formatted', status = 'unknown')
do k = 1, N/2-1
   if (r_or_p .eq. 'r') then
      kH0 = k * 2.d0*pi / (rgrid(nr) - rgrid(1))
      write(lun_spec, "(2(1x, e12.5))") kH0, E(k)
   else if (r_or_p .eq. 'p') then
      write(lun_spec, "(i4, 1x, e12.5)") k, E(k)
   end if
end do
call dfftw_destroy_plan(plan)
end

!----------------------------------------------------------------------------------85

