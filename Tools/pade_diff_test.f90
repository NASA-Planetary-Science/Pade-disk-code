! Test of cache-friendly versus cache-unfriendly flux differentiation. 

! Note:  When we use the Shu and Osher trick the number of intervals
! equals the number of points: each grid point lies in the center of
! an interval.

implicit none

! Number of intervals = # of points:
integer, parameter :: nx = 64, ny = 64, nz = 256
real(8), dimension(nx, ny, nz) :: q, F3, Fd3 ! 3d
real(8), dimension(nx, ny, 0:nz) :: F3_in_place
real(8), dimension(nz) :: F1, Fd1, cc_nz ! 1d
real(8), dimension(nx) :: Fx, Fdx ! 1d
real(8), dimension(nz+1) :: cc_nzp1
real(8), dimension(nx+1) :: cc_nxp1
integer :: ix, iy, iz, nbundle
real(8) :: pi, value
real(8) :: dz, z(nz), Ji_z(nz)
real(8) :: dx, x(nx), Ji_x(nz)
real(8) :: exact_deriv(nz)
real(8) :: cpu_secs1, cpu_secs2, cpu_line, cpu_bundle, cpu_line_first_index

print *, ' Calculating derivatives along z'
print *, ' nx = ', nx, ' ny = ', ny, ' nz = ', nz

! Compute q:
! ~~~~~~~~~~

pi = 4.0d0 * atan(1.0d0)

dz = 2.d0*pi / nz
do iz = 1, nz
   z(iz) = 0.5d0*dz + (iz - 1)*dz
   Ji_z(iz) = 1.d0/dz
   value = cos(z(iz))
   exact_deriv(iz) = -sin(z(iz))
   do iy = 1, ny
      do ix = 1, nx
         q(ix, iy, iz) = value
      end do
   end do
end do

! Test routine for differentiation along a single z line:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Precompute an array used in the tridiagonal solver.
call pre_compute_cc (nz+1, 0.25d0, 1.0d0, 3.0d0, cc_nzp1)

call cpu_time(cpu_secs1)
do iy = 1, ny
   do ix = 1, nx
      ! Notional flux calculation:
      do iz = 1, nz
         F1(iz) = q(ix, iy, iz)
      end do
      call flux_diff (nz, F1, Ji_z, cc_nzp1, Fd1)
      do iz = 1, nz
         Fd3(ix, iy, iz) = Fd1(iz)
      end do
   end do
end do
call cpu_time(cpu_secs2)
cpu_line = cpu_secs2 - cpu_secs1

open (unit = 1, file = 'Fd_34_line.dat',  form = 'formatted', status = 'unknown')
open (unit = 2, file = 'err_34_line.dat', form = 'formatted', status = 'unknown')
ix = 3
iy = 4
do iz = 1, nz
   write(1, "(3(1x, e12.5))") z(iz), q(ix, iy, iz), Fd3(ix, iy, iz)
   write(2, "(3(1x, e12.5))") z(iz), ABS(Fd3(ix, iy, iz) - exact_deriv(iz))
end do

close(1)
close(2)
print *, ' wrote Fd_34_line.dat and err_34_line.dat with nz = ', nz

! Straight Pade differentiation (non_Shu):
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
call pre_compute_cc (nz, 0.25d0, 1.0d0, 3.0d0, cc_nz)

call cpu_time(cpu_secs1)
do iy = 1, ny
   do ix = 1, nx
      ! Notional flux calculation:
      do iz = 1, nz
         F1(iz) = q(ix, iy, iz)
      end do
      call pade_diff_precomputed_cc (nz, cc_nz, F1, Fd1)
      do iz = 1, nz
         Fd3(ix, iy, iz) = Fd1(iz) * Ji_z(iz)
      end do
   end do
end do
call cpu_time(cpu_secs2)
cpu_line = cpu_secs2 - cpu_secs1

open (unit = 2, file = 'err_34_straight_Pade.dat', form = 'formatted', status = 'unknown')
ix = 3
iy = 4
do iz = 1, nz
   write(2, "(3(1x, e12.5))") z(iz), ABS(Fd3(ix, iy, iz) - exact_deriv(iz))
end do

close(2)
print *, ' wrote err_34_straight_Pade.dat with nz = ', nz

! Test bundle routine:
! ~~~~~~~~~~~~~~~~~~~~

! Extra (nx, ny, nz) arrays.
!    (1) flux
!    (2) derivative of the flux,
!    (3) H for Shu and Osher trick
!    (4) Derivative of H which is the numerical flux whose difference gives the
!        derivative of the flux

nbundle = nx*ny
print *, ' calling flux_diff_bundle'

call cpu_time(cpu_secs1)
! Notional flux calculation:
do iz = 1, nz
   do iy = 1, ny
      do ix = 1, nx
         F3(ix, iy, iz) = q(ix, iy, iz)
      end do
   end do
end do
call flux_diff_bundle(nbundle, nz, F3(1, 1, 1), Ji_z, Fd3(1, 1, 1))
call cpu_time(cpu_secs2)
cpu_bundle = cpu_secs2 - cpu_secs1

ix = 3
iy = 4

open (unit = 1, file = 'Fd_34_bundle.dat', form = 'formatted', status = 'unknown')
do iz = 1, nz
   write(1, "(3(1x, e12.5))") z(iz), q(ix, iy, iz), Fd3(ix, iy, iz)
end do
close(1)


! Test bundled in-place routine:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Extra (nx, ny, nz) arrays.
!    (1) flux
!    (2) derivative of the flux --- no need
!    (3) H for Shu and Osher trick --- no need
!    (4) Derivative of H which is the numerical flux whose difference gives the
!        derivative of the flux

nbundle = nx*ny

call cpu_time(cpu_secs1)
! Notional flux calculation:
do iz = 1, nz
   do iy = 1, ny
      do ix = 1, nx
         F3_in_place(ix, iy, iz) = q(ix, iy, iz)
      end do
   end do
end do
print *, ' calling flux_diff_bundle_in_place'
call flux_diff_bundle_in_place(nbundle, nz, F3_in_place(1, 1, 0), Ji_z)
call cpu_time(cpu_secs2)
cpu_bundle = cpu_secs2 - cpu_secs1

ix = 3
iy = 4

open (unit = 1, file = 'Fd_34_bundle_in_place.dat', form = 'formatted', status = 'unknown')
do iz = 1, nz
   write(1, "(3(1x, e12.5))") z(iz), q(ix, iy, iz), F3_in_place(ix, iy, iz)
end do
close(1)

! Line by line w.r.t. first index:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dx = 2.d0*pi / (nx - 1)
do ix = 1, nx
   x(ix) = (ix - 1)*dx
   Ji_x(ix) = 1.d0/dx
end do

do iz = 1, nz
   do iy = 1, ny
      do ix = 1, nx
         q(ix, iy, iz) = cos(x(ix))
      end do
   end do
end do

call pre_compute_cc (nx+1, 0.25d0, 1.0d0, 3.0d0, cc_nxp1)

call cpu_time(cpu_secs1)
do iz = 1, nz
   do iy = 1, ny
      ! Notional flux calculation:
      do ix = 1, nx
         Fx(ix) = q(ix, iy, iz)
      end do
      call flux_diff (nx, Fx, Ji_x, cc_nxp1, Fdx)
      do ix = 1, nx
         Fd3(ix, iy, iz) = Fdx(ix)
      end do
   end do
end do
call cpu_time(cpu_secs2)
cpu_line_first_index = cpu_secs2 - cpu_secs1

iy = 3
iz = 4

open (unit = 1, file = 'Fdx_34_line.dat', form = 'formatted', status = 'unknown')

do ix = 1, nx
   write(1, "(3(1x, e12.5))") x(ix), q(ix, iy, iz), Fd3(ix, iy, iz)
end do

close(1)

print *, ' for line routine called many times, cpu secs = ', cpu_line
print *, ' for bundle routine,                 cpu_secs = ', cpu_bundle
print *, '                            ratio line/bundle = ', cpu_line/cpu_bundle
print *, ' for line routine on first index,    cpu secs = ', cpu_line_first_index
print *, '                          ratio x line/bundle = ', cpu_line_first_index/cpu_bundle

stop
end
