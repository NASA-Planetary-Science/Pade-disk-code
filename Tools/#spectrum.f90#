!----------------------------------------------------------------------------------85

! This is needed for FFTW:
use, intrinsic :: iso_c_binding
implicit none

! This is needed for FFTW:
include 'fftw3.f03'

character(80) filename
integer, parameter :: lun_tecplot = 1
real(8) :: t, r_prime, k, lambda, pi, rmid, x
integer :: nr, nz, ir, iz, irec, N
real(8), allocatable, dimension(:) :: rgrid, zgrid
! The velocities are all fluctuations from the basic state:
real(8), allocatable, dimension(:, :) :: rho, ur, uz, uphi
character(1) :: dum_char
integer :: include_uz, include_ur, include_uphi

! For the FFTW library:
integer(8) :: plan
real(8),    allocatable, dimension(:) :: f
complex(8), allocatable, dimension(:) :: fhat
real(8), allocatable, dimension(:) :: S, Hann

pi = 4.0d0 * atan(1.d0)

write(6, "(' enter tecplot file name--->', $)")
read(5, *) filename

write(6, "(' enter 3 flags (0/1) for including uz, ur, and uphi--->', $)")
read(5, *) include_uz, include_ur, include_uphi

open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', access = 'direct', recl = 13*6 + 1)

read(lun_tecplot, 1, rec = 1) t
1 format (13x, e12.5)
print *, ' t = ', t
read(lun_tecplot, 2, rec = 2) dum_char
2 format(a1)
read(lun_tecplot, 3, rec = 3) nr, nz
3 format(7x, i4, 4x, i4)
print *, ' nr = ', nr, ' nz = ', nz

allocate(rgrid(nr), zgrid(nz), rho(nr,nz), ur(nr,nz), uz(nr,nz), uphi(nr,nz))

! For the FFT:
N = nr - 1
allocate(f(N))
allocate(fhat(0: N/2))
allocate(S(0:N/2))

! Initialize the FFT:
print *, ' calling dfftw_plan_dft_r2c_1d'
call dfftw_plan_dft_r2c_1d(plan, N, f, fhat, FFTW_ESTIMATE)
print *, ' returned from dfftw_plan_dft_r2c_1d'

do iz = 1, nz
   do ir = 1, nr
      print *, ' ir = ', ir, ' iz = ', iz
      irec = (iz - 1)*nr + ir + 3
      read(lun_tecplot, "(6(1x, e12.5))", rec = irec) rgrid(ir), zgrid(iz), rho(ir,iz), ur(ir,iz), &
           uz(ir,iz), uphi(ir,iz)
   end do
end do
close(lun_tecplot)

allocate(hann(nr))
open(unit = 1, file = 'hann.dat', form = 'formatted', status = 'unknown')
do ir = 1, N
   lambda = rgrid(nr) - rgrid(1)
   rmid = rgrid(1) + 0.5d0*lambda
   x = rgrid(ir) - rmid
   hann(ir) = cos(pi*x/lambda)**2
end do
close(1)

S = 0.d0
do iz = 1, nz
   if (include_ur .eq. 1) then
      ! ur:
      do ir = 1, N
         f(ir) = ur(ir, iz)*hann(ir)      
      end do
      call dfftw_execute(plan)
      fhat(0) = fhat(0)/N
      do ir = 1, N/2
         fhat(ir) = fhat(ir)*2.d0/N
      end do
   
      do ir = 0, N/2-1
         S(ir) = S(ir) + fhat(ir)*conjg(fhat(ir)) 
      end do
   end if

   if (include_uz .eq. 1) then
      ! Add in uz:
      do ir = 1, N
         f(ir) = uz(ir, iz)*hann(ir)      
      end do
      call dfftw_execute(plan)
      fhat(0) = fhat(0)/N
      do ir = 1, N/2
         fhat(ir) = fhat(ir)*2.d0/N
      end do
   
      do ir = 0, N/2-1
         S(ir) = S(ir) + fhat(ir)*conjg(fhat(ir)) 
      end do
   end if

   if (include_uphi .eq. 1) then
      ! Add in uphi
      do ir = 1, N
         f(ir) = uphi(ir, iz)*hann(ir)      
      end do
      call dfftw_execute(plan)
      fhat(0) = fhat(0)/N
      do ir = 1, N/2
         fhat(ir) = fhat(ir)*2.d0/N
      end do
   
      do ir = 0, N/2-1
         S(ir) = S(ir) + fhat(ir)*conjg(fhat(ir)) 
      end do
   end if
end do

do ir = 0, N/2 - 1
   S(ir) = S(ir) / nz
end do
   
! Test:
!lambda = rgrid(nr) - rgrid(1)
!k      = 2.d0*pi/lambda
!print *, ' lambda = ',lambda
!do ir = 1, N
!   r_prime = rgrid(ir) - rgrid(1)
!   f(ir) = 1.d0 + sin(k*r_prime) + cos(k*r_prime)
!end do

open(unit = 1, file = 'spectrum.dat', form = 'formatted', status = 'unknown')
do ir = 1, N/2-1
   write(1, "(i4, 1x, e12.5)") ir, S(ir)
end do

call dfftw_destroy_plan(plan)

end

!----------------------------------------------------------------------------------85

