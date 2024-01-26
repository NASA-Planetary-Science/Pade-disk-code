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

real(8), allocatable, dimension(:,:) :: rho

! Three velocity components:
real(8), allocatable, dimension(:,:,:) :: u
integer, parameter :: r_comp = 1, phi_comp = 2, z_comp = 3

integer :: icomp

character(1) :: dum_char

! For the FFTW library:
integer(8) :: plan
real(8),    allocatable, dimension(:) :: f
complex(8), allocatable, dimension(:) :: fhat
real(8), allocatable, dimension(:) :: S, Hann

pi = 4.0d0 * atan(1.d0)

write(6, "(' enter tecplot file name for primitive data in a horizontal plane--->', $)")
read(5, *) filename

write(6, "(' enter 3 flags (0/1) for including uz, ur, and uphi--->', $)")
read(5, *) include_uz, include_ur, include_uphi

open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
     access = 'direct', recl = 13*7 + 1)

read(lun_tecplot, 1, rec = 1) t
1  format (' TITLE = "t = ', e12.5)
read(lun_tecplot, 2, rec = 2)
read(lun_tecplot, 3, rec = 3) nr, nphi_plus_1
3 format(' ZONE I=', i4, ',', ' J=', i4)
nphi = nphi_plus_1 - 1

kr_max = nr/2 - 1

allocate(xgrid(nr,nphi), ygrid(nr,nphi), rho(nr,nphi), u(nr,nphi,3))

do iphi = 1, nphi
   do ir = 1, nr
      irec = (iphi - 1)*nr + ir + 3
      read(lun_tecplot, "(6(1x, e12.5))",rec=irec) xgrid(ir,iphi),ygrid(ir,iphi),&
           rho(ir,iphi), u(ir,iphi,r_comp), u(ir,iphi,phi_comp), u(ir,iphi, z_comp)
   end do
end do
close(lun_tecplot)

do icomp
do ir = 1, nr
   rho_bar(ir)  = 0.d0
   ur_(ir)   = 0.d0
   uphi_ave(ir) = 0.d0
   uz_ave(ir)   = 0.d0
   do iphi = 1, nphi
      rho_ave (ir) = rho_ave (ir) + rho(ir,iphi)
      ur_ave  (ir) = ur_ave  (ir) + ur  (ir,iphi)
      uphi_ave(ir) = uphi_ave(ir) + uphi(ir,iphi)
      uz_ave  (ir) = uz_ave  (ir) + uz  (ir,iphi)
   end do

   rho_ave (ir) = rho_ave (ir) / nphi
   ur_ave  (ir) = ur_ave  (ir) / nphi
   uphi_ave(ir) = uphi_ave(ir) / nphi   
   uz_ave  (ir) = uz_ave  (ir) / nphi
end do

! Next compute velocity fluctuations:
do iphi = 1, nphi
   do ir = 1, nr
      u_prime  (ir, iphi, icomp) = ur  (ir, iphi) - ur_ave  (ir)
      uz_prime  (ir, iphi) = uz  (ir, iphi) - uz_ave  (ir)
      uphi_prime(ir, iphi) = uphi(ir, iphi) - uphi_ave(ir)      
   end do
end do

! Next compute the fluctuation ke 
      
! For the FFT:
!ndim_r = nr - 1
allocate(f(nr))
allocate(fhat(0: nr/2))
allocate(Er(0:nr/2), Ez(0:nr/2), Ephi(0:nr/2))

! Initialize the FFT:
print *, ' calling dfftw_plan_dft_r2c_1d'
call dfftw_plan_dft_r2c_1d(plan_r, nr, f, fhat, FFTW_ESTIMATE)
print *, ' returned from dfftw_plan_dft_r2c_1d'

allocate(hann(nr))
open(unit = 1, file = 'hann.dat', form = 'formatted', status = 'unknown')
do ir = 1, nr
   lambda = rgrid(nr) - rgrid(1)
   rmid = rgrid(1) + 0.5d0*lambda
   x = rgrid(ir) - rmid
   hann(ir) = cos(pi*x/lambda)**2
end do
close(1)

! Compute spectrum w.r.t. r:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~

! Radial component:
Er = 0.d0
do iphi = 1, nphi
   do ir = 1, nr
      f(ir) = ur_prime(ir, iphi)*hann(ir)      
   end do
   call dfftw_execute(plan_r)
   fhat(0) = fhat(0)/nr
   do ir = 1, nr/2
      fhat(ir) = fhat(ir)*2.d0/nr
   end do
   
   do kr = 0, nr/2
      Er(kr) = Er(kr) + fhat(kr)*conjg(fhat(kr)) 
   end do
end do
Er = Er / nphi

! z-component:
Ez = 0.d0
do iphi = 1, nphi
   do ir = 1, nr
      f(ir) = uz_prime(ir, iphi)*hann(ir)      
   end do
   call dfftw_execute(plan_r, f, fhat)
   fhat(0) = fhat(0)/nr
   do ir = 1, nr/2
      fhat(ir) = fhat(ir)*2.d0/nr
   end do
   
   do kr = 0, nr/2
      Ez(kr) = Ez(kr) + fhat(kr)*conjg(fhat(kr)) 
   end do
end do
Ez = Ez / nphi

! phi component:
Ephi = 0.d0
do iphi = 1, nphi
   do ir = 1, nr
      f(ir) = uphi_prime(ir, iphi)*hann(ir)      
   end do
   call dfftw_execute(plan_r, f, fhat)
   fhat(0) = fhat(0)/nr
   do ir = 1, nr/2
      fhat(ir) = fhat(ir)*2.d0/nr
   end do
   
   do kr = 0, nr/2
      Ephi(kr) = Ephi(kr) + fhat(kr)*conjg(fhat(kr)) 
   end do
end do
Ez = Ez / nphi





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

