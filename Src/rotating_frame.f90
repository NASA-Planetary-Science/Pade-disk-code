!----------------------------------------------------------------------------------85

module rotating_frame

logical :: apply_rotating_frame
real(8) :: Omega_frame

end module rotating_frame

!----------------------------------------------------------------------------------85

subroutine activate_rotating_frame(apply_rotating_frame_arg, Omega_frame_arg)

use rotating_frame
implicit none
logical :: apply_rotating_frame_arg
real(8) :: Omega_frame_arg

apply_rotating_frame = apply_rotating_frame_arg
Omega_frame = Omega_frame_arg

end subroutine activate_rotating_frame

!----------------------------------------------------------------------------------85

subroutine add_coriolis_and_centrifugal_terms(q, qdot)

use grid
use partition_data
use rotating_frame
use dof_indices
implicit none
real(8), intent(in),  dimension (sr:er, sphi:ephi, nz, ndof) :: q
real(8), intent(out), dimension (sr:er, sphi:ephi, nz, ndof) :: qdot

! Local:
integer :: iz, iphi, ir
real(8) :: uphi, ur

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! Centrifugal term:
         qdot(ir, iphi, iz, rmom) = qdot(ir, iphi, iz, rmom) + &
              q(ir, iphi, iz, irho)*Omega_frame**2 * rgrid(ir)

         ! Adding + 2 Omega rho*uphi:
         qdot(ir, iphi, iz, rmom) = qdot(ir, iphi, iz, rmom) + &         
              2.d0*Omega_frame*q(ir, iphi, iz, amom) / rgrid(ir)
         ! Adding - 2 Omega rho*ur * r
         qdot(ir, iphi, iz, amom) = qdot(ir, iphi, iz, amom) - &         
            2.d0*Omega_frame*q(ir, iphi, iz, rmom) * rgrid(ir)         
      end do
   end do
end do

end subroutine add_coriolis_and_centrifugal_terms

!----------------------------------------------------------------------------------85
