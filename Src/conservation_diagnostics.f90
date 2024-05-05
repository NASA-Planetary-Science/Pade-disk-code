!----------------------------------------------------------------------------------85

! The fluxes do not include any pressure which is correct for mass and
! angular momentum only.  See notes of 9/27/2020 for the algebra.

! Works for the axisymmetric case.  Probably will not work for the planar case
! nz = 1.

!----------------------------------------------------------------------------------85

subroutine output_conservation_diagnostics(t, q)

use dof_indices
use grid
use partition_data
use logical_units
use pade_coefficients
implicit none
real(8) :: t
real(8), intent(in), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
real(8) :: mass, ang_mom
real(8) :: left_mass_flux_out, right_mass_flux_out
real(8) :: left_amom_flux_out, right_amom_flux_out
real(8) :: top_mass_flux_out,  bot_mass_flux_out
real(8) :: top_amom_flux_out,  bot_amom_flux_out

call volume_integral(q(sr,  sphi, 1, irho), mass)
call r_boundary_fluxes(irho, q, left_mass_flux_out, right_mass_flux_out)
call z_boundary_fluxes(irho, q, top_mass_flux_out,  bot_mass_flux_out)

call volume_integral(q(sr,  sphi, 1, amom), ang_mom)
call r_boundary_fluxes(amom, q, left_amom_flux_out, right_amom_flux_out)
call z_boundary_fluxes(amom, q, top_amom_flux_out,  bot_amom_flux_out)

if (my_node .eq. 0) then
   open(unit = lun_history(1), file = 'mass_conservation.his', form = 'formatted', &
        status = 'unknown', access = 'append')
   open(unit = lun_history(2), file = 'amom_conservation.his', form = 'formatted', &
        status = 'unknown', access = 'append')

   write(lun_history(1), "(6(1x, e16.9))") t, mass, left_mass_flux_out, &
        right_mass_flux_out, top_mass_flux_out,  bot_mass_flux_out
   write(lun_history(2), "(6(1x, e16.9))") t, ang_mom, left_amom_flux_out, &
        right_amom_flux_out, top_amom_flux_out,  bot_amom_flux_out
   close(lun_history(1))
   close(lun_history(2))
end if

end subroutine output_conservation_diagnostics

!----------------------------------------------------------------------------------85

subroutine r_boundary_fluxes(idof, q, left_flux_out, right_flux_out)

! Note: Only node 0 will have the correct fluxes.

! Computes the r surface fluxes of the degree of freedom "idof" in q.

use grid
use partition_data
use dof_indices
use pade_coefficients
#ifdef mpi_code
use mpi
#endif
implicit none
integer, intent(in) :: idof
real(8), intent(in), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8), intent(out) :: left_flux_out, right_flux_out

! Local:
integer :: iphi, iz, ier
real(8) :: ur
real(8), dimension(sphi:ephi) :: integral_wrt_z_left, integral_wrt_z_right
real(8) :: my_left_flux_out, my_right_flux_out

! Left flux:
if (sr .eq. 1) then ! does this proc. have the left r boundary.
   ! Integrate w.r.t. to z:
   do iphi = sphi, ephi
      integral_wrt_z_left(iphi) = 0.d0
      do iz = 1, nz
         ur = -q(1, iphi, iz, rmom) / q(1, iphi, iz, irho)
         integral_wrt_z_left(iphi) = integral_wrt_z_left(iphi) + &
              conservative_pade_weight_z(iz)*rgrid(1)*ur*q(1, iphi, iz, idof)
      end do
   end do

   ! Integrate w.r.t. phi:
   if (nphi .eq. 1) then
      my_left_flux_out = integral_wrt_z_left(1)
   else
      my_left_flux_out = 0.d0
      do iphi = sphi, ephi
         my_left_flux_out = my_left_flux_out + integral_wrt_z_left(iphi)
      end do
      my_left_flux_out = my_left_flux_out * dphi
   end if
else ! This proc does not have rmin.
   my_left_flux_out = 0.d0
end if

! Right flux:
if (er .eq. nr) then  ! If this proc has the last r surface.
   ! Integrate w.r.t. to z:
   do iphi = sphi, ephi
      integral_wrt_z_right(iphi) = 0.d0
      do iz = 1, nz
         ur = q(nr, iphi, iz, rmom) / q(nr, iphi, iz, irho)
         integral_wrt_z_right(iphi) = integral_wrt_z_right(iphi) + &
              conservative_pade_weight_z(iz)*rgrid(nr)*ur*q(nr, iphi, iz, idof)
      end do
   end do

   ! Integrate w.r.t. phi:
   if (nphi .eq. 1) then
      my_right_flux_out = integral_wrt_z_right(1)
   else
      my_right_flux_out = 0.d0
      do iphi = sphi, ephi
         my_right_flux_out = my_right_flux_out + integral_wrt_z_right(iphi)
      end do
      my_right_flux_out = my_right_flux_out * dphi
   end if
else ! This proc. does not have rmax.
   my_right_flux_out = 0.d0
end if

#ifdef mpi_code
call mpi_reduce(my_left_flux_out,  left_flux_out,  1, mpi_double_precision, &
     mpi_sum, 0, mpi_comm_world, ier)
call mpi_reduce(my_right_flux_out, right_flux_out, 1, mpi_double_precision, &
     mpi_sum, 0, mpi_comm_world, ier)
#else
left_flux_out  = my_left_flux_out
right_flux_out = my_right_flux_out
#endif
      
end subroutine r_boundary_fluxes

!----------------------------------------------------------------------------------85

subroutine z_boundary_fluxes(idof, q, top_flux_out, bot_flux_out)

! z-surface flux of the "q" variable "idof".
! Note: Only node 0 will have the correct "top_flux_out" and "bot_flux_out" 

use grid
use partition_data
use dof_indices
use pade_coefficients
#ifdef mpi_code
use mpi
#endif
implicit none
integer :: idof ! The "q" variable whose fluxes you want.
real(8), intent(in), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8), intent(out) :: top_flux_out, bot_flux_out

! Local:
integer :: ir, iphi, ier
real(8) :: uz
real(8), dimension(sphi:ephi) :: integral_wrt_r_top, integral_wrt_r_bot
real(8) :: my_top_flux_out, my_bot_flux_out

! Integrate with respect to r:
do iphi = sphi, ephi
   integral_wrt_r_top(iphi) = 0.d0
   integral_wrt_r_bot(iphi) = 0.d0   
   do ir = sr, er
      uz = q(ir, iphi, nz, zmom) / q(ir, iphi, nz, irho)
      integral_wrt_r_top(iphi) = integral_wrt_r_top(iphi) + &
           conservative_pade_weight_r(ir)*rgrid(ir)*uz*q(ir, iphi, nz, idof)
      ! Bottom surface:
      uz = -q(ir, iphi, 1, zmom) / q(ir, iphi, 1, irho)       
      integral_wrt_r_bot(iphi) = integral_wrt_r_bot(iphi) + &
           conservative_pade_weight_r(ir)*rgrid(ir)*uz*q(ir, iphi, 1, idof)
   end do
end do

! Integrate with respect to phi:
if (nphi .eq. 1) then
   my_top_flux_out = integral_wrt_r_top(1)
   my_bot_flux_out = integral_wrt_r_bot(1)
else
   my_top_flux_out = 0.d0
   my_bot_flux_out = 0.d0

   do iphi = sphi, ephi
      my_top_flux_out = my_top_flux_out + integral_wrt_r_top(iphi) 
      my_bot_flux_out = my_bot_flux_out + integral_wrt_r_bot(iphi)
   end do

   my_top_flux_out = my_top_flux_out * dphi
   my_bot_flux_out = my_bot_flux_out * dphi
end if

#ifdef mpi_code
call mpi_reduce(my_top_flux_out, top_flux_out, 1, mpi_double_precision, &
     mpi_sum, 0, mpi_comm_world, ier)
call mpi_reduce(my_bot_flux_out, bot_flux_out, 1, mpi_double_precision, &
     mpi_sum, 0, mpi_comm_world, ier)
#else
top_flux_out = my_top_flux_out
bot_flux_out = my_bot_flux_out
#endif
      
end subroutine z_boundary_fluxes

!----------------------------------------------------------------------------------85

subroutine volume_integral(integrand, vol_int)

! Note: Only node 0 will have the correct "vol_int"

! Performs a volume integral of the quantity given in "integrand."
! Note: The Jacobian factor will be inserted here so the "integrand"
! should not include it.

use grid
use partition_data
use pade_coefficients ! Need this for quadrature weights.
#ifdef mpi_code
use mpi
#endif
implicit none
real(8), dimension(sr:er, sphi:ephi, nz), intent(in)  :: integrand
real(8), intent(out) :: vol_int

! Local:
integer :: iphi, ir, iz, ier
real(8), dimension(sr:er, sphi:ephi) :: integral1
real(8), dimension(       sphi:ephi) :: integral2
real(8) :: my_integral

! Integrate w.r.t. z:
do iphi = sphi, ephi
   do ir = sr, er
      integral1(ir, iphi) = 0.d0
      do iz = 1, nz
         integral1(ir, iphi) = integral1(ir, iphi) + &
              conservative_pade_weight_z(iz)*integrand(ir, iphi, iz)
      end do
   end do
end do

! Integrate w.r.t. r:
do iphi = sphi, ephi
   integral2(iphi) = 0.d0
   do ir = sr, er
      integral2(iphi) = integral2(iphi) + &
           conservative_pade_weight_r(ir)*rgrid(ir)*integral1(ir, iphi)
   end do
end do

! Finally, integrate w.r.t. phi:
if (nphi .eq. 1) then
   my_integral = integral2(1)
else
   my_integral = 0.d0
   do iphi = sphi, ephi
      my_integral = my_integral + integral2(iphi)
   end do
   my_integral = my_integral*dphi
end if

#ifdef mpi_code
call mpi_reduce(my_integral, vol_int, 1, mpi_double_precision, &
     mpi_sum, 0, mpi_comm_world, ier)
#else
vol_int = my_integral
#endif

end subroutine volume_integral

!----------------------------------------------------------------------------------85
