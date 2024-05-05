!**********************************************************************************85

implicit none
real(8) :: rmin, r0, rmax, stretch_r, savings_factor_r
real(8) :: zmin, z0, zmax, stretch_z, savings_factor_z
integer :: nr, nz, nr_u, nz_u
logical :: activate_r_stretching, activate_z_stretching

real(8), allocatable, dimension(:) :: rgrid, Jr, Ji_r
real(8), allocatable, dimension(:) :: zgrid, Jz, Ji_z
integer :: i

namelist /input/ activate_r_stretching, activate_z_stretching, &
                 rmin, r0, rmax, nr_u, nr, &
                 zmin, z0, zmax, nz_u, nz

!xmin   = 0.0d0
!x0     = 0.5d0
!xmax   = 1.5d0
!n      = 256
!n_unif = 150

open (unit = 3, file = 'input_for_stretched', form = 'formatted', status = 'old')
read (3, nml = input)
close (3)

print *, ' Read following parameters from file = input_for_stretched'
print *, ' activate_r_stretching = ', activate_r_stretching
print *, ' activate_z_stretching = ', activate_z_stretching

if (activate_r_stretching) then
   print *, '    rmin   = ', rmin
   print *, '    r0     = ', r0
   print *, '    rmax   = ', rmax
   print *, '    nr_u   = ', nr_u
   print *, '    nr     = ', nr
   print *, ' enter anything to continue'
   read(5,*)

   allocate(rgrid(nr), Ji_r(nr), Jr(nr))

   call stretched_mesh_r(rmin, rmax, r0, nr_u, nr, rgrid, Jr, Ji_r, &
        stretch_r, savings_factor_r)

   print *,  ' stretch_r = ', stretch_r, ' savings_factor_r = ', savings_factor_r
   read(5,*)

   open (unit = 1, file = 'rgrid.dat', form = 'formatted', status = 'unknown')
   open (unit = 2, file = 'Jr.dat',    form = 'formatted', status = 'unknown')
   do i = 1, nr
      write(1, "(2(1x, i4, 1x, e12.5))") i, rgrid(i)
      write(2, "(2(1x, i4, 1x, e12.5))") i, Jr(i)   
   end do
   close(1)
   close(2)
   print *, ' Wrote rgrid.dat and Jr.dat'   
end if

if (activate_z_stretching) then
   print *, '    zmin   = ', zmin
   print *, '    z0     = ', z0
   print *, '    zmax   = ', zmax
   print *, '    nz_u   = ', nz_u
   print *, '    nz     = ', nz
   print *, ' enter anything to continue'
   read(5,*)
   
   allocate(zgrid(nz), Ji_z(nz), Jz(nz))

   call stretched_mesh_z(zmin, zmax, z0, nz_u, nz, zgrid, Jz, Ji_z, &
        stretch_z, savings_factor_z)

   print *,  ' stretch_z = ', stretch_z, ' savings_factor_z = ', savings_factor_z
   read(5,*)
   
   open (unit = 1, file = 'zgrid.dat', form = 'formatted', status = 'unknown')
   open (unit = 2, file = 'Jz.dat',    form = 'formatted', status = 'unknown')
   do i = 1, nz
      write(1, "(2(1x, i4, 1x, e12.5))") i, zgrid(i)
      write(2, "(2(1x, i4, 1x, e12.5))") i, Jz(i)   
   end do
   close(1)
   close(2)
   print *, ' Wrote zgrid.dat and Jz.dat'      
end if

end

!----------------------------------------------------------------------------------85

subroutine stretched_mesh_r(rmin, rmax, r0, nr_u, nr, rgrid, Jr, Ji_r, &
     stretch_r, savings_factor_r)

! Creates a geometrically stretched mesh in r.

! Input:
! ~~~~~~
! rmin, rmax : Extent of mesh
! r0         : The region of uniform spacing is [rmin, x0]
! nr_u       : # of points in uniform region
! nr         : Total number of points

! Output:
! ~~~~~~~
! rgrid(n)         : grid points
! Jr(n)            : Jacobian dx/d_xi, where xi is the grid index variable
! Ji_r(n)          : Inverse Jacobian d xi / dx
! stretch_r        : stretching factor
! savings_factor_r : # of grid points if uniform everywhere / actual number of grid points

implicit none
real(8), intent(in) :: rmin, rmax, r0
integer, intent(in) :: nr, nr_u

real(8), dimension(nr), intent(out) :: rgrid, Jr, Ji_r
real(8), intent(out) :: stretch_r, savings_factor_r

! Local:
real(8) :: dr0, alpha, alpha_next
integer :: i, iter
integer, parameter :: maxit = 2000

dr0 = (r0 - rmin) / (nr_u - 1.d0)
savings_factor_r = (rmax - rmin) / dr0 / nr

! Perform Gauss-Seidel itrerations to determine the stretching factor:

alpha = log(1.05d0)
do iter = 1, maxit
   alpha_next = alpha_func(alpha)
   print *, ' iter = ', iter, ' alpha = ', alpha_next
   if (abs(alpha_next - alpha) .lt. 1.d-15) go to 100
   alpha = alpha_next
end do

print *, ' Iterations did not converge'
stop

100 continue

stretch_r = exp(alpha)
print *, ' alpha = ', alpha
print *, ' Stretching factor = ', stretch_r

do i = 1, nr_u
   rgrid(i) = rmin + dr0*(i - 1)
   Jr  (i)  = dr0
   Ji_r(i)  = 1.d0 / dr0
end do

do i = nr_u + 1, nr
   rgrid(i) = r0 + dr0/log(stretch_r) * (stretch_r**(i - nr_u) - 1)
   Jr (i)  = dr0 * stretch_r**(i - nr_u)
   Ji_r(i) = 1.d0 / Ji_r(i) 
end do

contains

!---
real(8) function alpha_func(alpha) ! alpha = ln r where r = stretch ratio
real(8) :: alpha

real(8) :: arg ! Local
arg = 1 + (rmax - r0) / dr0 * alpha
alpha_func = log(arg) / real(nr - nr_u)
end
!---

end subroutine stretched_mesh_r

!----------------------------------------------------------------------------------85

subroutine stretched_mesh_z(zmin, zmax, z0, nz_u, nz, zgrid, Jz, Ji_z, &
     stretch_z, savings_factor_z)

! Creates a geometrically stretched mesh in z, the uniform region being [-z0, z0]
! If n is odd,  nz_u must be odd
! If n is even, nz_u must be even

! Input:
! ~~~~~~
! zmin, zmax : Extent of mesh
! z0         : The region of uniform spacing is [-z0, z0]
! nz_u       : # of points in uniform region including the end-points.
! nz         : Total number of points

! Output:
! ~~~~~~~
! zgrid(n)         : grid points
! Jz(n)            : Jacobian dz/d_xi, where xi is the grid index variable
! Ji_z(n)          : Inverse Jacobian d xi / dz
! stretch_z        : stretching factor
! savings_factor_z : # of grid points if uniform everywhere / actual number of grid points

implicit none
real(8), intent(in) :: zmin, zmax, z0
integer, intent(in) :: nz_u, nz

real(8), dimension(nz), intent(out) :: zgrid, Jz, Ji_z
real(8), intent(out) :: stretch_z, savings_factor_z

! Local:
real(8) :: dz0, alpha, alpha_next
integer :: i, iter, n1, n2
integer, parameter :: maxit = 2000
logical :: even_nz, even_nz_u, odd_nz, odd_nz_u

if (mod(nz, 2) .eq. 0) then
   even_nz = .true.
else
   even_nz = .false.
end if

if (mod(nz_u, 2) .eq. 0) then
   even_nz_u = .true.
else
   even_nz_u = .false.
end if

if (even_nz .and. .not. even_nz_u) then
   print *, ' Both n and nu must be even'
   stop
end if

odd_nz   = .not. even_nz
odd_nz_u = .not. even_nz_u
if (odd_nz .and. .not. odd_nz_u) then
   print *, ' Both n and nu must be odd'
   stop
end if

dz0 = 2.d0*z0 / (nz_u - 1.d0) ! spacing in uniform region.
savings_factor_z = (zmax - zmin) / dz0 / nz

n1 = (nz - nz_u)/2 + 1 ! This is correct for both odd and even
n2 = n1 + nz_u - 1

! Perform Gauss-Seidel iterations to determine the stretching factor in [z0, zmax].

alpha = log(1.05d0)
do iter = 1, maxit
   alpha_next = alpha_func(alpha)
   print *, ' iter = ', iter, ' alpha = ', alpha_next
   if (abs(alpha_next - alpha) .lt. 1.d-15) go to 100
   alpha = alpha_next
end do

print *, ' Iterations did not converge'
stop

100 continue

stretch_z = exp(alpha)
print *, ' alpha = ', alpha
print *, ' Stretching factor in [z0, zmax] = ', stretch_z
 
do i = n1, n2
   zgrid(i) = -z0 + dz0*(i - n1)
   Jz(i)    = dz0
   Ji_z(i)  = 1.d0 / dz0
end do

do i = n2 + 1, nz
   zgrid(i) = z0 + dz0/log(stretch_z) * (stretch_z**(i - n2) - 1.d0)
   Jz  (i)  = dz0 * stretch_z**(i - n2)
   Ji_z(i)  = 1.d0 / Ji_z(i) 
end do

do i = 1, n1 - 1
   zgrid(i) = -zgrid(nz - i + 1)
   Jz   (i) =  Jz   (nz - i + 1)
   Ji_z (i) =  Ji_z (nz - i + 1)
end do

contains

!---
real(8) function alpha_func(alpha) ! alpha = ln r where r = stretch ratio
real(8) :: alpha

real(8) :: arg ! Local
arg = 1.d0 + (zmax - z0) / dz0 * alpha
alpha_func = log(arg) / real(nz - nz_u)
end
!---

end subroutine stretched_mesh_z

!----------------------------------------------------------------------------------85
