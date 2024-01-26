!----------------------------------------------------------------------------------85

module stretched_mesh
logical :: stretched_r, stretched_z
real(8) :: r0, z0
integer :: nr_u, nz_u
end module stretched_mesh

!----------------------------------------------------------------------------------85

subroutine make_grid

! Non cell-centered grid.

use grid
use math_constants, only:pi
use partition_data
use total_allocated_words, only: n_words_allocated
use logical_units
use dof_indices
use stretched_mesh
implicit none

! Local:
integer :: ir, iz, iphi
real(8) :: dr_unif, dz_unif
logical :: planar_case, case1, case2, r_only_case

! For output purposes:
if (mod(nr, 2) .eq. 0) then
   ir_mid = nr/2
else
   ir_mid = (nr + 1)/2
end if

if (mod(nz, 2) .eq. 0) then
   iz_mid = nz/2
else
   iz_mid = (nz + 1)/2
end if

case1 = (nz .eq. 1) .and.                                            (nr .ne. 1) .and. (nphi .ne. 1)
case2 = (nz .ne. 1) .and. suppress_z_derivatives_when_nz_not_1 .and. (nr .ne. 1) .and. (nphi .ne. 1)
planar_case = case1 .or. case2

r_only_case = (nz .eq. 1) .and. (nphi .eq. 1) .and. (nr .ne. 1)

! Grid:
! The nphi + 1 is there so we plot the end-point of thr periodic interval:
allocate (rgrid(nr), zgrid(nz), phi_grid(nphi + 1))
n_words_allocated = n_words_allocated + nr + nz + nphi + 1

! Inverse Jacobians:
allocate (Ji_r(nr), Ji_z(nz))
n_words_allocated = n_words_allocated + nr + nz   

! Jacobians:
allocate (Jr(nr), Jz(nz))
n_words_allocated = n_words_allocated + nr + nz   

! Grid sizes:
allocate (dr(nr), dz(nz), r_dphi(nr))
n_words_allocated = n_words_allocated + nr + nz + nphi

! Minimize grid size at each grid-point (used for time-step determination):
allocate(min_grid_size(sr:er,nz))
n_words_allocated = n_words_allocated + mr*nz

allocate(l_grid_squared(sr:er, nz))
n_words_allocated = n_words_allocated + mr*nz

if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
   if (periodic_z) then
      dz_periodic = (zmax - zmin) / nz
      do iz = 1, nz
         zgrid(iz) = zmin + (iz - 1)*dz_periodic
         Ji_z(iz) = 1.d0 / dz_periodic ! inverse Jacobian
         Jz  (iz) = dz_periodic ! Jacobian
      end do
   else if (stretched_z) then
      call stretched_mesh_z(zmin, zmax, z0, nz_u, nz, zgrid, Jz, Ji_z)
   else 
      dz_unif = (zmax - zmin) / (nz - 1)
      do iz = 1, nz
         zgrid(iz) = zmin + (iz - 1)*dz_unif
         Ji_z(iz) = 1.d0 / dz_unif ! inverse Jacobian
         Jz  (iz) = dz_unif ! Jacobian
      end do
   end if

   do iz = 1, nz - 1
      dz(iz) = zgrid(iz + 1) - zgrid(iz)
   end do
   dz(nz) = zgrid(nz) - zgrid(nz-1)
else
   zmin = 0.0d0
   zmax = 0.0d0
   zgrid(1) = 0.0d0
   dz(1)    = 0.0d0
end if

if (suppress_z_derivatives_when_nz_not_1) then
   do iz = 1, nz
      zgrid(iz) = 0.0d0
      dz(iz)    = 0.0d0
   end do
end if

if (nr .ne. 1) then
   if (stretched_r) then
      call stretched_mesh_r(rmin, rmax, r0, nr_u, nr, rgrid, Jr, Ji_r)
   else
      dr_unif   = (rmax - rmin) / (nr - 1)
      do ir = 1, nr
         rgrid(ir) = rmin + (ir - 1) * dr_unif
         Ji_r(ir) = 1.d0 / dr_unif ! inverse Jacobian
         Jr  (ir) = dr_unif ! inverse Jacobian
      end do
   end if

   do ir = 1, nr - 1
      dr(ir) = rgrid(ir + 1) - rgrid(ir)
   end do
   dr(nr) = rgrid(nr) - rgrid(nr-1)
else
   rgrid(1) = rmin
   dr(1)    = 0.0d0
end if

if (nphi .ne. 1) then
   ! This is used for the Fargo trick when the phi domain is not zero to twopi.
   Delta_phi_domain = phi_max - phi_min
   dphi = Delta_phi_domain / nphi
   do iphi = 1, nphi + 1
      phi_grid(iphi) = phi_min + 0.5d0*dphi + (iphi - 1) * dphi
   end do
   dphi_inv = 1.0d0 / dphi
else
   dphi = 2.d0*pi ! This should really be zero
end if

! We are in subroutine make_grid

! Grid size in the phi direction:
do ir = 1, nr
   r_dphi(ir) = dphi*rgrid(ir)
end do

! min_grid_size  : Minimum grid-size for time-step selection for constant viscosity.
! l_grid_squared : Used for artificial pressure treatment.  Also for Moin et al sgs model.
if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1) .and. (nr .ne. 1) .and. (nphi .ne. 1)) then
   ! 3D case:
   do iz = 1, nz
      do ir = sr, er
         min_grid_size (ir, iz) = MIN(dz(iz), dr(ir), r_dphi(ir))
         l_grid_squared(ir, iz) = (dr(ir)*dz(iz)*r_dphi(ir))**(2.d0/3.d0) ! volume based
      end do
   end do
   ! Grid spacing vector:
   if (my_node .eq. 0) print *, 'Computed min_grid_size for 3D case'
else if (planar_case) then
   ! Planar case (vertically integrated, r, phi):
   do iz = 1, nz
      do ir = sr, er
         min_grid_size (ir, iz) = MIN(dr(ir), r_dphi(ir))
         l_grid_squared(ir, iz) = dr(ir)*r_dphi(ir) ! area based
      end do
   end do
   if (my_node .eq. 0) print *, 'Computed min_grid_size for planar case'
else if ((nz .ne. 1) .and. (nr .ne. 1) .and. (nphi .eq. 1)) then
   ! Axisymmetric case (z, r):
   do iz = 1, nz
      do ir = sr, er
         min_grid_size (ir, iz) = MIN(dz(iz), dr(ir))
         l_grid_squared(ir, iz) = dz(iz)*dr(ir) ! area based         
      end do
   end do
   if (my_node .eq. 0) print *, 'Computed min_grid_size for axisymmetric case'
else if ((nz .ne. 1) .and. (nr .eq. 1) .and. (nphi .eq. 1)) then
   ! Vertical only case:
   do iz = 1, nz
      do ir = sr, er
         min_grid_size (ir, iz) = dz(iz)
         l_grid_squared(ir, iz) = dz(iz)**2                  
      end do
   end do
else if (r_only_case) then
   do iz = 1, nz
      do ir = sr, er
         min_grid_size (ir, iz) = dr(ir)
         l_grid_squared(ir, iz) = dr(ir)*dr(ir)
      end do
   end do
else
   if (my_node .eq. 0) then
      print *, ' subroutine make_grid:'
      print *, ' You are running a new case type for which you need to add coding for'
      print *, ' the minimum grid size in subroutine make_grid'
   end if
   call terminate_with_no_save(1)
end if

! Trapezoidal rule weight for taking averages.  See notes of Nov. 26, 2018
! Also for conservation_diagnostics.f90.
allocate(trap_weight_r(nr))
n_words_allocated = n_words_allocated + nr
if (nr .ne. 1) then
   trap_weight_r(1)  = 0.5d0*dr(1)
   trap_weight_r(nr) = 0.5d0*dr(nr-1)
   do ir = 2, nr-1
      trap_weight_r(ir) = 0.5d0*(dr(ir-1) + dr(ir))
   end do
else
   trap_weight_r(1) = 1.0d0
end if

allocate(trap_weight_z(nz))
n_words_allocated = n_words_allocated + nz
if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
   trap_weight_z(1)  = 0.5d0*dz(1)
   trap_weight_z(nz) = 0.5d0*dz(nz-1)
   do iz = 2, nz-1
      trap_weight_z(iz) = 0.5d0*(dz(iz-1) + dz(iz))
   end do
else
   trap_weight_z(1) = 1.0d0
   if (suppress_z_derivatives_when_nz_not_1) trap_weight_z(2) = 0.0d0
end if

if (my_node .eq. 0) then
   print *, ' subroutine make_grid: Finished setting-up grid:'
   print *, ' nr = ', nr, ' nz = ', nz, ' nphi = ', nphi
end if

gridded = .true.
return
end subroutine make_grid

!----------------------------------------------------------------------------------85

subroutine stretched_mesh_r(rmin, rmax, r0, nr_u, nr, rgrid, Jr, Ji_r)

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

use logical_units
use partition_data
implicit none
real(8), intent(in) :: rmin, rmax, r0
integer, intent(in) :: nr, nr_u

real(8), dimension(nr), intent(out) :: rgrid, Jr, Ji_r

! Local:
real(8) :: dr0, alpha, alpha_next
integer :: i, iter
integer, parameter :: maxit = 2000
real(8) :: stretch_r, savings_factor_r

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' my_node = 0: stretched_mesh_in_r has been called'
   print *, '    rmin = ', rmin, ' rmax = ', rmax
   print *, '    r0 = ', r0
   print *, '    nr_u = ', nr_u, ' nr = ', nr
end if
#endif

dr0 = (r0 - rmin) / (nr_u - 1.d0)
savings_factor_r = (rmax - rmin) / dr0 / nr

! Perform Gauss-Seidel itrerations to determine the stretching factor:

alpha = log(1.05d0)
do iter = 1, maxit
   alpha_next = alpha_func(alpha)
   if (isnan(alpha_next)) then
      print *, ' stretched_mesh_r: alpha_next = ', alpha_next
      call terminate_with_no_save(1)
   end if
   ! if (my_node .eq. 0) print *, ' iter = ', iter, ' alpha = ', alpha_next
   if (abs(alpha_next - alpha) .lt. 1.d-15) go to 100
   alpha = alpha_next
end do

print *, ' Iterations for stretching factor in r did not converge'
call terminate_with_no_save(1)

100 continue

stretch_r = exp(alpha)
if (my_node .eq. 0) then
   print *, ' alpha = ', alpha
   print *, ' Stretching factor = ', stretch_r
   print *, ' Savings factor in r = ', savings_factor_r
end if

do i = 1, nr_u
   rgrid(i) = rmin + dr0*(i - 1)
   Jr  (i)  = dr0
   Ji_r(i)  = 1.d0 / dr0
end do

do i = nr_u + 1, nr
   rgrid(i) = r0 + dr0/log(stretch_r) * (stretch_r**(i - nr_u) - 1)
   Jr (i)  = dr0 * stretch_r**(i - nr_u)
   Ji_r(i) = 1.d0 / Jr(i) 
end do

! Output:
if (my_node .eq. 0) then
   write(6, "(' stretch_r        = ', e12.5)") stretch_r
   write(6, "(' savings_factor_r = ', e12.5)") savings_factor_r
   open (unit = lun_general_purpose_1, file = 'rgrid.dat', form = 'formatted', status = 'unknown')
   open (unit = lun_general_purpose_2, file = 'Jr.dat',    form = 'formatted', status = 'unknown')

   do i = 1, nr
      write(lun_general_purpose_1, "(2(1x, i4, 1x, e12.5))") i, rgrid(i)
      write(lun_general_purpose_2, "(2(1x, i4, 1x, e12.5))") i, Jr(i)   
   end do
   close(lun_general_purpose_1)
   close(lun_general_purpose_2)
end if

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

subroutine stretched_mesh_z(zmin, zmax, z0, nz_u, nz, zgrid, Jz, Ji_z)

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

! Output to stdout:
! stretch_z        : stretching factor
! savings_factor_z : # of grid points if uniform everywhere / actual number of grid points

use partition_data
use logical_units ! for output
implicit none
real(8), intent(in) :: zmin, zmax, z0
integer, intent(in) :: nz_u, nz

real(8), dimension(nz), intent(out) :: zgrid, Jz, Ji_z


! Local:
real(8) :: dz0, alpha, alpha_next
integer :: i, iter, n1, n2
integer, parameter :: maxit = 2000
logical :: even_nz, even_nz_u, odd_nz, odd_nz_u
real(8) :: stretch_z, savings_factor_z

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' my_node = 0: stretched_mesh_in_z has been called'
   print *, '    zmin = ', zmin, ' zmax = ', zmax
   print *, '    z0 = ', z0
   print *, '    nz_u = ', nz_u, ' nz = ', nz
end if
#endif

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
   if (isnan(alpha_next)) then
      if (my_node .eq. 0) print *, ' stretched_mesh_z: alpha_next = ', alpha_next
      call terminate_with_no_save(1)
   end if   
   ! if (my_node .eq. 0) print *, ' iter = ', iter, ' alpha = ', alpha_next
   if (abs(alpha_next - alpha) .lt. 1.d-15) go to 100
   alpha = alpha_next
end do

print *, ' Iterations for stretching factor in z did not converge'
call terminate_with_no_save(1)

100 continue

stretch_z = exp(alpha)
if (my_node .eq. 0) then
   print *, ' alpha = ', alpha
   print *, ' Stretching factor in [z0, zmax] = ', stretch_z
   print *, ' Savings factor in z = ', savings_factor_z
end if
 
do i = n1, n2
   zgrid(i) = -z0 + dz0*(i - n1)
   Jz(i)    = dz0
   Ji_z(i)  = 1.d0 / dz0
end do

do i = n2 + 1, nz
   zgrid(i) = z0 + dz0/log(stretch_z) * (stretch_z**(i - n2) - 1.d0)
   Jz  (i)  = dz0 * stretch_z**(i - n2)
   Ji_z(i)  = 1.d0 / Jz(i) 
end do

do i = 1, n1 - 1
   zgrid(i) = -zgrid(nz - i + 1)
   Jz   (i) =  Jz   (nz - i + 1)
   Ji_z (i) =  Ji_z (nz - i + 1)
end do

! Output:
if (my_node .eq. 0) then
   write(6, "(' stretch_z        = ', e12.5)") stretch_z
   write(6, "(' savings_factor_z = ', e12.5)") savings_factor_z
   open (unit = lun_general_purpose_1, file = 'zgrid.dat', form = 'formatted', status = 'unknown')
   open (unit = lun_general_purpose_2, file = 'Jz.dat',    form = 'formatted', status = 'unknown')

   do i = 1, nz
      write(lun_general_purpose_1, "(2(1x, i4, 1x, e12.5))") i, zgrid(i)
      write(lun_general_purpose_2, "(2(1x, i4, 1x, e12.5))") i, Jz(i)
   end do
   close(lun_general_purpose_1)
   close(lun_general_purpose_2)
end if

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
