!----------------------------------------------------------------------------------85

! See also subroutine activate_gravity in activate_routines.f90 

!----------------------------------------------------------------------------------85

module gravity_types

implicit none
integer :: mass_at_origin = 1, &
     thin_disk_mass_at_origin = 2, &
     thin_disk_mass_at_origin_no_radial = 3, &
     uniform_gz = 4

!>    mass_at_origin : $g = -GM/ R^2 \widehat{R}$ where $R$ is the spherical radius and
!>                     $\widehat{R}$ is the unit vector along the spherical radius.
!>    thin_disk_mass_at_origin : g = -GM/r^2 * (rhat + (z/r) zhat) where rhat and zhat are
!>                               unit vectors.
!>    thin_disk_mass_at_origin_no_radial : Similar above where the radial component is
!>                                         set to zero.
!>    uniform_gz : A uniform vertical gravity with value = gz_uniform.

end module gravity_types

!----------------------------------------------------------------------------------85

module gravity

implicit none
   logical :: gravity_flag
   integer :: gravity_type  ! Will be selected from module gravity_types
   real(8) :: GM, gz_uniform
   ! Radial and vertical gravitational acceleration
   ! Dimensions of gr and gz are (nr, nz).
   real(8), allocatable, dimension(:, :) :: gr, gz
   
   ! Currently used for the Fargo trick.
   real(8), allocatable, dimension(:) :: u_Kepler, Omega_Kepler

!----------------------------------------------------------------------------------85   
   
contains

!----------------------------------------------------------------------------------85

subroutine add_gravity_terms(q, qdot)

!> Adds gravity terms to qdot.  It assumes axisymmetric (r, z) gravity and adds the
!> r and z components of gravity terms to the momentum equation and energy equation.
!> The updating of the energy equation is not necessary for the isothermal case and
!> in the future should be done only if a pre-processor flag is set.

use grid, only: ndof, nz, nphi, nr
use dof_indices, only: irho, zmom, amom, rmom, ener
use partition_data
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in)    :: q
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(inout) :: qdot

! Local:
integer :: iz, iphi, ir
real(8) :: fr, fz, ur, uz

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! Gravitational forces per unit volume:
         fr = q(ir, iphi, iz, irho) * gr(ir, iz)
         fz = q(ir, iphi, iz, irho) * gz(ir, iz)
         ur = q(ir, iphi, iz, rmom) / q(ir, iphi, iz, irho)
         uz = q(ir, iphi, iz, zmom) / q(ir, iphi, iz, irho)
         qdot(ir,iphi,iz,rmom) = qdot(ir,iphi,iz,rmom) + fr 
         qdot(ir,iphi,iz,zmom) = qdot(ir,iphi,iz,zmom) + fz
         qdot(ir,iphi,iz,ener) = qdot(ir,iphi,iz,ener) + fr*ur + fz*uz
      end do
   end do
end do

end subroutine add_gravity_terms

!----------------------------------------------------------------------------------85

end module gravity

!----------------------------------------------------------------------------------85
