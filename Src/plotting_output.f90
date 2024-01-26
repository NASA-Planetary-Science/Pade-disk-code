!----------------------------------------------------------------------------------85

real(8) function c_sound_at_point(q, ir, iphi, iz)

!> Convenience routine used by plotting output routines.
!> Try not to use this inside a do loop for huge arrays.

use grid
use dof_indices
use thermal_parameters
use partition_data
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in) :: q
integer, intent(in) :: ir, iphi, iz

! Local:
real(8) :: pressure

if (isothermal) then
   c_sound_at_point = SQRT(ci_squared_initial(ir, iz))
else
   pressure = q(ir, iphi, iz, ener) * gm1
   c_sound_at_point = SQRT(gamma * pressure / q(ir, iphi, iz, irho))
   if (isnan(c_sound_at_point)) then
      print *, ' c_sound_at_point = ', c_sound_at_point
      print *, ' pressure = ', pressure
      print *, ' rho = ', q(ir, iphi, iz, irho)
   end if
end if

end function c_sound_at_point

!----------------------------------------------------------------------------------85

real(8) function pressure_at_point(q, ir, iphi, iz)

! Convenience routine used by plotting output routines.
! Try not to use this for huge arrays.

use grid
use partition_data
use dof_indices
use thermal_parameters
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
integer :: ir, iphi, iz

! Local:

if (isothermal) then
   pressure_at_point = q(ir, iphi, iz, irho) * ci_squared_initial(ir, iz)
else
   pressure_at_point = q(ir, iphi, iz, ener) * gm1
end if

end function pressure_at_point

!----------------------------------------------------------------------------------85

subroutine vertical_profiles(ir, iphi, location_id, t)

! The location id will be used to compose the file name.

use grid
use q_array
use dof_indices, only: irho, zmom, rmom, amom, ener
use logical_units
use thermal_parameters
implicit none
integer :: ir, iphi
character(4) :: location_id
real(8) :: t

! Local:
character(40), dimension(5) :: file
integer :: iz, i
real(8) :: rho, uphi, uz, ci, p
! Function to be called:
real(8) :: pressure_at_point

write (file(1), "(a4, '_rho_vs_z_t=',   i6.6, f0.4, '.dat')") location_id, int(t), t - int(t)
write (file(2), "(a4, '_Omega_vs_z_t=', i6.6, f0.4, '.dat')") location_id, int(t), t - int(t)
write (file(3), "(a4, '_uz_vs_z_t=',    i6.6, f0.4, '.dat')") location_id, int(t), t - int(t)
write (file(4), "(a4, '_p_vs_z_t=',     i6.6, f0.4, '.dat')") location_id, int(t), t - int(t)
write (file(5), "(a4, '_ci_vs_z_t=',    i6.6, f0.4, '.dat')") location_id, int(t), t - int(t)

do i = 1, 5
   open (lun_profile(i), file = file(i), form = 'formatted', status = 'unknown')
end do

do iz = 1, nz
   rho       = q(ir, iphi, iz, irho)
   uphi      = q(ir, iphi, iz, amom) / rho / rgrid(ir)
   uz        = q(ir, iphi, iz, zmom) / rho

   ! pressure_at_point is a function call:
   p  = pressure_at_point(q, ir, iphi, iz)
   ci = SQRT(p/rho)

   write (lun_profile(1), "(2(1x, e13.5e3))") zgrid(iz), rho
   write (lun_profile(2), "(2(1x, e13.5e3))") zgrid(iz), uphi/rgrid(ir)
   write (lun_profile(3), "(2(1x, e13.5e3))") zgrid(iz), uz
   write (lun_profile(4), "(2(1x, e13.5e3))") zgrid(iz), p
   write (lun_profile(5), "(2(1x, e13.5e3))") zgrid(iz), ci   
end do

do i = 1, 5
   close (lun_profile(i))
end do

print *, ' Wrote vertical profiles ', file(1), file(2), file(3), file(4)

end subroutine vertical_profiles

!----------------------------------------------------------------------------------85

subroutine radial_profiles(iz, iphi, location_id, t)

! The location id will be used to compose the file name.

use grid
use q_array
use dof_indices, only: irho, zmom, rmom, amom, ener
use logical_units
use thermal_parameters
implicit none
integer :: iz, iphi
character(4) :: location_id
real(8) :: t

! Local:
character(40), dimension(8) :: file
integer :: ir, i
real(8) :: rho, uphi, uz, ur, ci, p, cs
! Functions to be called:
real(8) :: pressure_at_point, c_sound_at_point

! The SP forces a sign so there is no blank in the filename.
write (file(1), "(a4, '_rho_vs_r_t=',   f9.3, '.dat')") location_id, t
write (file(2), "(a4, '_Omega_vs_r_t=', f9.3, '.dat')") location_id, t
write (file(3), "(a4, '_uz_vs_r_t=',    f9.3, '.dat')") location_id, t
write (file(4), "(a4, '_p_vs_r_t=',     f9.3, '.dat')") location_id, t
write (file(5), "(a4, '_ci_vs_r_t=',    f9.3, '.dat')") location_id, t
write (file(6), "(a4, '_ur_vs_r_t=',    f9.3, '.dat')") location_id, t
write (file(7), "(a4, '_Mr_vs_r_t=',    f9.3, '.dat')") location_id, t
write (file(8), "(a4, '_Mz_vs_r_t=',    f9.3, '.dat')") location_id, t

do i = 1, 8
   open (lun_profile(i), file = file(i), form = 'formatted', status = 'unknown')
end do

do ir = 1, nr
   rho       = q(ir, iphi, iz, irho)
   uphi      = q(ir, iphi, iz, amom) / rho / rgrid(ir)
   uz        = q(ir, iphi, iz, zmom) / rho
   ur        = q(ir, iphi, iz, rmom) / rho

   ! pressure_at_point is a function call:
   p  = pressure_at_point(q, ir, iphi, iz)
   ci = SQRT(p/rho)
   cs = c_sound_at_point(q, ir, iphi, iz)

   write (lun_profile(1), "(2(1x, e13.5e3))") rgrid(ir), rho
   write (lun_profile(2), "(2(1x, e13.5e3))") rgrid(ir), uphi/rgrid(ir)
   write (lun_profile(3), "(2(1x, e13.5e3))") rgrid(ir), uz
   write (lun_profile(4), "(2(1x, e13.5e3))") rgrid(ir), p
   write (lun_profile(5), "(2(1x, e13.5e3))") rgrid(ir), ci
   write (lun_profile(6), "(2(1x, e13.5e3))") rgrid(ir), ur
   write (lun_profile(7), "(2(1x, e13.5e3))") rgrid(ir), ur/cs
   write (lun_profile(8), "(2(1x, e13.5e3))") rgrid(ir), uz/cs   
end do

do i = 1, 8
   close (lun_profile(i))
end do

print *, ' Wrote radial profiles at t = ', t
!read (5, *)

end subroutine radial_profiles

!----------------------------------------------------------------------------------85

subroutine tecplot_meridional_plane(q, iphi, L_scale, T_scale, t, istep)

!> This subroutine works in either serial or parallel mode.  The direct access
!> and specification of the record number for each data line allows each processor
!> to write to the proper record.
!>
!> L_scale and T_scale are the scalings that will be used to output lengths and velocities.
!> The density and pressure will be left unscaled.

use dof_indices, only: irho, zmom, rmom, amom
use logical_units
use grid
use thermal_parameters
use partition_data
use basic_state
use physical_constants
implicit none
real(8), intent(in), dimension (sr:er, sphi:ephi, nz, ndof) :: q
integer, intent(in) :: iphi
real(8), intent(in) :: t, L_scale, T_scale
integer, intent(in) :: istep

! Local:
character(80) :: filename
integer :: ir, iz, irec
real(8) :: rho, rho_prime, ur, uz, uphi, uphi_prime, uphi_r, p, cs, Mach_r, Mach_z, rho_ur, rho_uz, &
     uphi_basic_scalar, pressure
! Function called:
real(8) :: c_sound_at_point, pressure_at_point
real(8) :: t_scaled, U_scale
logical :: i_have_iphi

t_scaled = t / T_scale
U_scale = L_scale / T_scale

write (filename, "('merid_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')") iphi, &
     int(t_scaled), t_scaled - int(t_scaled), istep

! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
     access = 'direct', recl = 13*13 + 1)

if (my_node .eq. 0) then
   ! For emacs, try just char(10) = new line, i.e., avoid char(13) = CR.
   write (lun_tecplot, 1, rec = 1) t/T_scale, char(10)
   1  format ('TITLE = "t = ', e12.5, '"', a1)
   write (lun_tecplot, 2, rec = 2) char(10)
      2  format('VARIABLES = "r", "z", "rho", "ur", "uz", "uphi", "uphi*r", "Mach_r", "Mach_z", "rho*ur", "rho*uz", "rho_prime", "uphi_prime"', a1)
   write (lun_tecplot, 3, rec = 3) nr, nz, char(10)
   3 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', a1)
end if

! Write only of this processor has the iphi:
i_have_iphi = (iphi .ge. sphi) .and. (iphi .le. ephi)
if (i_have_iphi) then
   do iz = 1, nz
      do ir = sr, er
         rho    = q(ir, iphi, iz, irho)
         ur     = q(ir, iphi, iz, rmom) / rho
         uz     = q(ir, iphi, iz, zmom) / rho
         uphi   = q(ir, iphi, iz, amom) / (rho * rgrid(ir))
         uphi_r = uphi * rgrid(ir)
         cs   = c_sound_at_point(q, ir, iphi, iz)
         Mach_r = ur/cs
         Mach_z = uz/cs
         ! pressure = pressure_at_point(q, ir, iphi, iz)
         rho_prime  = rho  - rho_basic(ir, iz)
         uphi_prime = uphi - uphi_basic(ir, iz)
      
         irec = (iz - 1)*nr + ir + 3

         write(lun_tecplot, "(13(1x, e12.5), a1)", rec = irec) rgrid(ir)/L_scale, &
              zgrid(iz)/L_scale, rho, ur/U_scale, uz/U_scale, uphi/U_scale, &
              uphi_r/U_scale/L_scale, Mach_r, Mach_z, rho*ur, rho*uz, &
              rho_prime, uphi_prime/U_scale, char(10)
      end do
   end do
end if

close(lun_tecplot)

if (my_node .eq. 0) then
   print *, ' rank 0: In subroutine tecplot_meridional_plane wrote tecplot file ', &
        filename
end if

end subroutine tecplot_meridional_plane

!----------------------------------------------------------------------------------85

subroutine tecplot_raw_q_in_meridional_plane (prefix, q, iphi, t)

! Used for debugging.  Can be used for anything with the same structure as "q".

! This subroutine works in either serial or parallel mode.  The direct access
! and specification of the record number for each data line allows each processor
! to write to the proper record.

use dof_indices
use logical_units
use grid
use thermal_parameters
use partition_data
implicit none
character(10) :: prefix
real(8), dimension (sr:er, sphi:ephi, nz, ndof) :: q
integer, intent(in) :: iphi
real(8), intent(in) :: t

! Local:
character(80) :: filename
integer :: ir, iz, idof, irec
logical :: i_have_iphi

write (filename, 4) prefix, int(t), t - int(t)
4 format (a10, '_t_', i6.6, f0.4, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', access = 'direct', recl = 13*7 + 1)

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
   1  format ('TITLE = "t = ', e12.5, '"', a1, a1)
   write (lun_tecplot, 2, rec = 2) char(13), char(10)
   2  format('VARIABLES = "r", "z", "rho", "rmom", "zmom", "amom", "e"', a1, a1)
   write (lun_tecplot, 3, rec = 3) nr, nz, char(13), char(10)
   3 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', a1, a1)
end if

! Write only of this processor has the iphi:
i_have_iphi = (iphi .ge. sphi) .and. (iphi .le. ephi)
if (i_have_iphi) then
   do iz = 1, nz
      do ir = sr, er
         irec = (iz - 1)*nr + ir + 3
         write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) rgrid(ir), zgrid(iz), &
              (q(ir, iphi, iz, idof), idof = 1, 5), char(10)
      end do
   end do
end if

close(lun_tecplot)

if (my_node .eq. 0) then
   print *, ' rank 0: In subroutine tecplot_raw_q_in_meridional_plane wrote'
   print *, ' tecplot file ', filename
end if

end subroutine tecplot_raw_q_in_meridional_plane

!----------------------------------------------------------------------------------85

subroutine tecplot_scalar_in_meridional_plane(prefix5, u, iphi, L_scale, T_scale, &
     u_scale, t, istep)

! Creates a tecplot file for a scalar "u".  The filename will use the 5 character
! variable "prefix5".

! This subroutine works in either serial or parallel mode.  The direct access
! and specification of the record number for each data line allows each processor
! to write to the proper record.

use dof_indices
use logical_units
use grid
use thermal_parameters
use partition_data

implicit none
character(5) :: prefix5
real(8), dimension (sr:er, sphi:ephi, nz) :: u
integer, intent(in) :: iphi, istep
real(8), intent(in) :: L_scale, T_scale, u_scale, t

! Local:
character(80) :: filename
integer :: ir, iz, idof, irec, ier
real(8) :: t_scaled
logical :: i_have_iphi

t_scaled = t/T_scale

write (filename, 4) prefix5, iphi, int(t_scaled), t_scaled - int(t_scaled), istep
4 format (a5, '_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', access = 'direct', recl = 13*3 + 1)

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, iphi, char(13), char(10)
   1  format ('TITLE = "t = ', e12.5, ' iphi = ', i4, '"', a1, a1)
   write (lun_tecplot, 2, rec = 2) prefix5, char(13), char(10)
   2  format('VARIABLES = "r", "z", "', a5, '"', a1, a1)
   write (lun_tecplot, 3, rec = 3) nr, nz, char(13), char(10)
   3 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', a1, a1)
end if

! Write only of this processor has the iphi:
i_have_iphi = (iphi .ge. sphi) .and. (iphi .le. ephi)
if (i_have_iphi) then
   do iz = 1, nz
      do ir = sr, er
         irec = (iz - 1)*nr + ir + 3
         write(lun_tecplot, "(3(1x, e12.5), a1)", rec = irec) rgrid(ir)/L_scale, zgrid(iz)/L_scale, &
              u(ir, iphi, iz)/u_scale, char(10)
      end do
   end do
end if
close(lun_tecplot)

if (my_node .eq. 0) then
   print *, ' rank 0: In subroutine tecplot_scalar_in_meridional_plane'
   print *, ' wrote tecplot file ', filename
end if

end subroutine tecplot_scalar_in_meridional_plane

!----------------------------------------------------------------------------------85

subroutine tecplot_scalar_in_horizontal_plane(prefix10, u, iz, T_scale, t, istep)

! Note: Use iz_mid if you want the mid-plane.

! This subroutine works in either serial or parallel mode.  The direct access
! and specification of the record number for each data line allows each processor
! to write to the proper record.

use dof_indices
use logical_units
use grid
use thermal_parameters
use partition_data
use math_constants, only:pi
use xy_coordinates_of_grid
use basic_state
implicit none
character(10) :: prefix10
real(8), dimension (sr:er, sphi:ephi, nz), intent(in) :: u
integer, intent(in) :: iz, istep
real(8), intent(in) :: t, T_scale

! Local:
character(80) :: filename
integer :: ir, iphi, irec
integer :: iphi_rec

write (filename, 4) prefix10, int(t), t - int(t), istep, iz
4 format (a10, '_horiz_t_', i6.6, f0.4, '_', i7.7, '_iz_', i4.4, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
     access = 'direct', recl = 13*7 + 1)

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
   1  format (' TITLE = "t = ', e12.5, '"', a1, a1)
   write (lun_tecplot, 2, rec = 2) char(13), char(10)
2  format(' VARIABLES = "x", "y", "scalar"', a1, a1)
   write (lun_tecplot, 3, rec = 3) nr, nphi+1, char(13), char(10)
   3 format(' ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', a1, a1)
end if

do iphi = sphi, ephi
   do ir = sr, er
      irec = (iphi - 1)*nr + ir + 3
      write(lun_tecplot, "(3(1x, e12.5), a1)",rec=irec) xgrid(ir,iphi),ygrid(ir,iphi),&
           u(ir,iphi,iz), char(10)
   end do
end do

! Periodic completion in phi:
if (sphi .eq. 1) then
   ! The data for iphi = 1, gets put in the records for nphi + 1
   iphi_rec = nphi + 1
   iphi     = 1
   do ir = sr, er
      irec = (iphi_rec - 1)*nr + ir + 3
      ! Note: I have stored the coordinates for this line in the iphi = 0 slot
      ! since the flow data sits in the same processor:
      write(lun_tecplot, "(3(1x, e12.5), a1)", rec = irec) &
           xgrid(ir,0), ygrid(ir,0), u(ir, iphi, iz), char(10)
   end do
end if

close(lun_tecplot)

if (my_node .eq. 0) then
   print *, ' rank 0: subroutine tecplot_scalar_in_horizontal_plane'
   print *, ' wrote tecplot file ', filename
end if
end subroutine tecplot_scalar_in_horizontal_plane

!----------------------------------------------------------------------------------85

subroutine tecplot_axisymmetric_scalar(prefix, u, L_scale, T_scale, u_scale, t)

! Can be used for anything with the same structure as "u" as declared below.

! This subroutine works in either serial or parallel mode.  The direct access
! and specification of the record number for each data line allows each processor
! to write to the proper record.

use dof_indices
use logical_units
use grid
use thermal_parameters
use partition_data
implicit none
character(10) :: prefix
real(8), dimension (sr:er, nz), intent(in) :: u
real(8),                        intent(in) :: L_scale, T_scale, u_scale, t

! Local:
character(80) :: filename
integer :: ir, iz, idof, irec 
real(8) :: t_scaled
logical :: i_have_iphi
integer, parameter :: iphi = 1

t_scaled = t/T_scale

write (filename, 4) prefix, int(t_scaled), t_scaled - int(t_scaled), iphi
4 format (a10, '_t_', i6.6, f0.4, '_iphi_', i4.4, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', access = 'direct', recl = 13*3 + 1)

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, iphi, char(13), char(10)
   1  format ('TITLE = "t = ', e12.5, ' iphi = ', i4, '"', a1, a1)
   write (lun_tecplot, 2, rec = 2) prefix, char(13), char(10)
   2  format('VARIABLES = "r", "z", "', a10, '"', a1, a1)
   write (lun_tecplot, 3, rec = 3) nr, nz, char(13), char(10)
   3 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', a1, a1)
end if

! Output only of this processor has phi = 1:
i_have_iphi = (iphi .ge. sphi) .and. (iphi .le. ephi)
if (i_have_iphi) then
   do iz = 1, nz
      do ir = sr, er
         irec = (iz - 1)*nr + ir + 3
         write(lun_tecplot, "(3(1x, e12.5), a1)", rec = irec) rgrid(ir)/L_scale, zgrid(iz)/L_scale, &
              u(ir, iz)/u_scale, char(10)
      end do
   end do
end if

close(lun_tecplot)

if (my_node .eq. 0) then
   print *, ' rank 0: subroutine tecplot_axisymmetric_scalar'
   print *, ' wrote tecplot file ', filename
end if

end subroutine tecplot_axisymmetric_scalar

!----------------------------------------------------------------------------------85

subroutine tecplot_horizontal_plane(q, iz, T_scale, t, istep)

! Note: Use iz_mid if you want the mid-plane.

! This subroutine works in either serial or parallel mode.  The direct access
! and specification of the record number for each data line allows each processor
! to write to the proper record.

use dof_indices
use logical_units
use grid
use thermal_parameters
use partition_data
use math_constants, only:pi
use xy_coordinates_of_grid
use basic_state
implicit none
real(8), dimension (sr:er, sphi:ephi, nz, ndof), intent(in) :: q
integer, intent(in) :: iz, istep
real(8), intent(in) :: t, T_scale

! Local:
character(80) :: filename
integer :: ir, iphi, irec
real(8) :: rho, rho_prime, ur, uphi, uphi_prime, uz, ci
! Function called:
real(8) :: pressure_at_point
integer :: iphi_rec

if (my_node .eq. 0) then
   print *, ' rank 0: subroutine tecplot_horizontal_plane has been called'
   print *, ' iz = ', iz
end if

write (filename, 4) int(t), t - int(t), istep, iz
4 format ('horiz_t_', i6.6, f0.4, '_', i7.7, '_iz_', i4.4, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
     access = 'direct', recl = 13*9 + 1)

! Removed char(13) at the beginning.  Previously was char(13), char(10)
if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
   1  format (' TITLE = "t = ', e12.5, '"', a1, a1)
   write (lun_tecplot, 2, rec = 2) char(13), char(10)
2  format(' VARIABLES = "x", "y", "rho", "d_rho", "ur", "uphi", "d_uphi", "uz", "ci"', a1, a1)
   write (lun_tecplot, 3, rec = 3) nr, nphi+1, char(13), char(10)
   3 format(' ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', a1, a1)
end if

do iphi = sphi, ephi
   do ir = sr, er
      rho        = q(ir, iphi, iz, irho)
      ur         = q(ir, iphi, iz, rmom) / rho
      uphi       = q(ir, iphi, iz, amom) / (rho * rgrid(ir))
      uphi_prime = uphi - uphi_basic(ir, iz)
      uz         = q(ir, iphi, iz, zmom) / rho
      rho_prime  = rho - rho_basic(ir, iz)
      ci         =  SQRT(pressure_at_point(q, ir, iphi, iz) / rho)

      irec = (iphi - 1)*nr + ir + 3
      write(lun_tecplot, "(9(1x, e12.5), a1)",rec=irec) xgrid(ir,iphi),ygrid(ir,iphi),&
           rho, rho_prime, ur, uphi, uphi_prime, uz, ci, char(10)
   end do
end do

! Periodic completion in phi:
if (sphi .eq. 1) then
   ! The data for iphi = 1, gets put in the records for nphi + 1
   iphi_rec = nphi + 1
   iphi     = 1
   do ir = sr, er
      rho        = q(ir, iphi, iz, irho)
      ur         = q(ir, iphi, iz, rmom) / rho
      uphi       = q(ir, iphi, iz, amom) / (rho * rgrid(ir))
      uphi_prime = uphi - uphi_basic(ir, iz)
      uz         = q(ir, iphi, iz, zmom) / rho
      rho_prime  = rho - rho_basic(ir, iz)
      ci         =  SQRT(pressure_at_point(q, ir, iphi, iz) / rho)      
         
      irec = (iphi_rec - 1)*nr + ir + 3
      ! Note: I have stored the coordinates for this line in the iphi = 0 slot
      ! since the flow data sits in the same processor:
      write(lun_tecplot, "(9(1x, e12.5), a1)",rec=irec) xgrid(0,iphi),ygrid(0,iphi),&
           rho, rho_prime, ur, uphi, uphi_prime, uz, ci, char(10)
   end do
end if

close(lun_tecplot)

if (my_node .eq. 0) then
   print *, ' rank 0: subroutine tecplot_horizontal_plane'
   print *, ' wrote tecplot file ', filename
end if
end subroutine tecplot_horizontal_plane

!----------------------------------------------------------------------------------85

subroutine make_iz_list(nz_plot, z_list, iz_list)

! From a list of "z" values makes a list of corresponding "iz" values that can be
! given to the subroutine below.

use grid
implicit none
integer, intent(in)  :: nz_plot
real(8), intent(in)  :: z_list (nz_plot)
integer, intent(out) :: iz_list(nz_plot)

! Local:
integer :: iz_plot, iz
real(8) :: dist, min_dist

do iz_plot = 1, nz_plot
   min_dist = 1.d100
   do iz = 1, nz
      dist = abs(z_list(iz_plot) - zgrid(iz)) 
      if (dist .lt. min_dist) then
         min_dist = dist
         iz_list(iz_plot) = iz
      end if
   end do
end do

end subroutine make_iz_list

!----------------------------------------------------------------------------------85

subroutine tecplot_many_horizontal_planes(q, nz_planes, iz_list, z_list, T_scale, &
                                          t, istep)

! This subroutine works in either serial or parallel mode.  The direct access
! and specification of the record number for each data line allows each processor
! to write to the proper record.

use dof_indices
use logical_units
use grid
use thermal_parameters
use partition_data
use math_constants, only:pi
use xy_coordinates_of_grid
use basic_state
implicit none
real(8), dimension (sr:er, sphi:ephi, nz, ndof), intent(in) :: q
integer, intent(in) :: nz_planes, iz_list(nz_planes), istep
! List of z coordinates of the z planes.  This is for file naming only:
real(8), intent(in), dimension(nz_planes) :: z_list
real(8), intent(in) :: t, T_scale

! Local:
character(80) :: filename
integer :: ir, iphi, irec
real(8) :: rho, rho_prime, ur, uphi, uphi_prime, uz, ci
! Function called:
real(8) :: pressure_at_point
integer :: iphi_rec, iplane, iz

if (my_node .eq. 0) then
   print *, ' rank 0: subroutine tecplot_horizontal_plane has been called'
end if

do iplane = 1, nz_planes
   iz     = iz_list(iplane)
   write (filename, 4) int(t), t - int(t), istep, iz
   4 format ('horiz_t_', i6.6, f0.4, '_', i7.7, '_iz_', i4.4, '.tec')
   ! recl is the maximum record length.
   open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
        access = 'direct', recl = 13*9 + 1)

   if (my_node .eq. 0) then
      ! char(10) = new line
      write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
      1  format (' TITLE = "t = ', e12.5, '"', a1, a1)
      write (lun_tecplot, 2, rec = 2) char(13), char(10)
      2  format(' VARIABLES = "x", "y", "rho", "d_rho", "ur", "uphi", "d_uphi", "uz", "ci"', a1, a1)
      write (lun_tecplot, 3, rec = 3) nr, nphi+1, char(13), char(10)
      3 format(' ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', a1, a1)
   end if

   do iphi = sphi, ephi
      do ir = sr, er
         rho        = q(ir, iphi, iz, irho)
         ur         = q(ir, iphi, iz, rmom) / rho
         uphi       = q(ir, iphi, iz, amom) / (rho * rgrid(ir))
         uphi_prime = uphi - uphi_basic(ir, iz)
         uz         = q(ir, iphi, iz, zmom) / rho
         rho_prime  = rho - rho_basic(ir, iz)
         ci         =  SQRT(pressure_at_point(q, ir, iphi, iz) / rho)

         irec = (iphi - 1)*nr + ir + 3
         write(lun_tecplot, "(9(1x, e12.5), a1)",rec=irec) xgrid(ir,iphi),ygrid(ir,iphi),&
              rho, rho_prime, ur, uphi, uphi_prime, uz, ci, char(10)
      end do
   end do

   ! Periodic completion in phi:
   if (sphi .eq. 1) then
      ! The data for iphi = 1, gets put in the records for nphi + 1
      iphi_rec = nphi + 1
      iphi     = 1
      do ir = sr, er
         rho        = q(ir, iphi, iz, irho)
         ur         = q(ir, iphi, iz, rmom) / rho
         uphi       = q(ir, iphi, iz, amom) / (rho * rgrid(ir))
         uphi_prime = uphi - uphi_basic(ir, iz)
         uz         = q(ir, iphi, iz, zmom) / rho
         rho_prime  = rho - rho_basic(ir, iz)
         ci         =  SQRT(pressure_at_point(q, ir, iphi, iz) / rho)      
         
         irec = (iphi_rec - 1)*nr + ir + 3
         ! Note: I have stored the coordinates for this line in the iphi = 0 slot
         ! since the flow data sits in the same processor:
         write(lun_tecplot, "(9(1x, e12.5), a1)",rec=irec) xgrid(0,iphi),ygrid(0,iphi),&
              rho, rho_prime, ur, uphi, uphi_prime, uz, ci, char(10)         
         !write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) &
         !     xgrid(ir,0), ygrid(ir,0), &
         !     rho, ur, uphi, uz, ci, char(10)
      end do
   end if

   close(lun_tecplot)

   if (my_node .eq. 0) then
      print *, ' rank 0: subroutine tecplot_horizontal_plane'
      print *, ' wrote tecplot file ', filename
   end if
end do

end subroutine tecplot_many_horizontal_planes

!----------------------------------------------------------------------------------85

subroutine tecplot_rho_and_u (q, t)

! This subroutine works in either serial or parallel mode.  The direct access
! and specification of the record number for each data line allows each processor
! to write to the proper record.

use dof_indices
use logical_units
use grid
use thermal_parameters
use partition_data
use basic_state
use xy_coordinates_of_grid
implicit none
real(8), intent(in), dimension (sr:er, sphi:ephi, nz, ndof) :: q
real(8), intent(in) :: t

! Local:
character(80) :: filename
integer :: ir, iz, iphi, irec, iphi_rec
real(8) :: rho, ur, uz, uphi

write (filename, 4) int(t), t - int(t)
4 format ('rho_and_u_t_', i6.6, f0.4, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', access = 'direct', recl = 13*7+1)

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
   1  format ('TITLE = "t = ', e12.5, '"', a1, a1)
   write (lun_tecplot, 2, rec = 2) char(13), char(10)
   2  format('VARIABLES = "x", "y", "z" "rho", "ur", "uz", "uphi"', a1, a1)
   ! 2  format('VARIABLES = "ir", "iphi", "iz"', a1, a1)
   write (lun_tecplot, 3, rec = 3) nr, nphi+1, nz, char(13), char(10)
   3 format('ZONE I=', i4, ',', ' J=', i4, ',', ' K=', i4, ' DATAPACKING=POINT', a1, a1)
end if

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         rho  = q(ir, iphi, iz, irho)
         ur   = q(ir, iphi, iz, rmom) / rho
         uz   = q(ir, iphi, iz, zmom) / rho
         uphi = q(ir, iphi, iz, amom) / (rho * rgrid(ir))

         irec = ir + (iphi-1)*nr + (iz-1)*nr*(nphi+1) + 3  

         write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) &
              xgrid(ir,iphi), ygrid(ir,iphi), zgrid(iz), rho, ur, &
              uz, uphi, char(10)

          !write(lun_tecplot, "(4(1x, i3), i9, a1)", rec = irec) &
          !     ir, iphi, iz, my_node, irec, char(10)

          !if (irec .eq. 16388) then
          !   print *, ' irec = ', irec, ' ir = ', ir, ' iphi = ', iphi, ' iz = ', iz, ' my_node = ', my_node
          !end if
      end do
   end do
end do

! Periodic completion in phi:
if (sphi .eq. 1) then
   ! The data for iphi = 1, gets put in the records for nphi + 1
   iphi_rec = nphi + 1
   iphi     = 1
   do iz = 1, nz
      do ir = sr, er
         irec = ir + (iphi_rec-1)*nr + (iz-1)*nr*(nphi+1) + 3
         
         write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) &
              xgrid(ir,0), ygrid(ir,0), zgrid(iz), rho, ur, &
              uz, uphi, char(10)
         
      end do
   end do
end if

close(lun_tecplot)

if (my_node .eq. 0) then
   print *, ' rank 0: subroutine tecplot_rho_and_u wrote tecplot file ', filename
end if

end subroutine tecplot_rho_and_u

!----------------------------------------------------------------------------------85

subroutine get_xy_coordinates_of_grid

! For plotting purposes.

use grid
use xy_coordinates_of_grid
use partition_data
implicit none
integer :: ir, iphi

#ifdef debug_print
   if (my_node .eq. 0) print *, ' my_node = 0: In subroutine get_xy_coordinates of grid'
#endif

! xy coordinates of grid horizontal planes for plotting purposes:

if (sphi .ne. 1) then
   allocate (xgrid(sr:er, sphi:ephi), ygrid(sr:er, sphi:ephi))
else
   ! So we can store the end of the periodic interval for plotting:
   allocate (xgrid(sr:er, 0:ephi), ygrid(sr:er, 0:ephi))
end if

#ifdef debug_print
   print *, ' my_node = ', my_node, ' sphi = ', sphi, ' ephi = ', ephi, ' sr = ', sr, ' er = ', er, &
            ' nphi = ', nphi  
#endif

do iphi = sphi, ephi
   do ir = sr, er
      xgrid(ir, iphi) = rgrid(ir) * cos(phi_grid(iphi))
      ygrid(ir, iphi) = rgrid(ir) * sin(phi_grid(iphi))
   end do
end do

if (sphi .eq. 1) then
   do ir = sr, er
      xgrid(ir, 0) = rgrid(ir) * cos(phi_grid(nphi + 1))
      ygrid(ir, 0) = rgrid(ir) * sin(phi_grid(nphi + 1))
   end do   
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' my_node = 0: About to return from get_xy_coordinates_of_grid'
#endif


end subroutine get_xy_coordinates_of_grid

!----------------------------------------------------------------------------------85

subroutine tecplot_fluc_vel_meridional_plane(q, iphi, L_scale, T_scale, t, istep)

! This subroutine works in either serial or parallel mode.  The direct access
! and specification of the record number for each data line allows each processor
! to write to the proper record.

! L_scale and T_scale are the scalings that will be used to output lengths and velocities.
! The density will be left unscaled.

use dof_indices, only: irho, zmom, rmom, amom
use logical_units
use grid
use thermal_parameters
use partition_data
use basic_state
use physical_constants
implicit none
real(8), intent(in), dimension (sr:er, sphi:ephi, nz, ndof) :: q
integer, intent(in) :: iphi
real(8), intent(in) :: t, L_scale, T_scale
integer, intent(in) :: istep

! Local:
character(90) :: filename
integer :: ir, iz, irec
real(8) :: rho, rho_prime, ur, uz, uphi_prime
! Function called:
real(8) :: t_scaled, U_scale
logical :: i_have_iphi

t_scaled = t / T_scale
U_scale = L_scale / T_scale

write (filename, "('fluc_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')") iphi, &
     int(t_scaled), t_scaled - int(t_scaled), istep

! recl is the maximum record length.  6 variables
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', access = 'direct', recl = 13*6 + 1)

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t/T_scale, char(13), char(10)
   1  format ('TITLE = "t = ', e12.5, '"', a1, a1)
   write (lun_tecplot, 2, rec = 2) char(13), char(10)
   2  format('VARIABLES = "r", "z", "rho", "ur", "uz", "uphi"', a1, a1)
   write (lun_tecplot, 3, rec = 3) nr, nz, char(13), char(10)
   3 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', a1, a1)
end if

! Output only of this processor has the iphi:
i_have_iphi = (iphi .ge. sphi) .and. (iphi .le. ephi)
if (i_have_iphi) then
   do iz = 1, nz
      do ir = sr, er
         rho        = q(ir, iphi, iz, irho)
         rho_prime  = rho - rho_basic(ir, iz)
         ur         = q(ir, iphi, iz, rmom) / rho
         uz         = q(ir, iphi, iz, zmom) / rho
         uphi_prime = q(ir, iphi, iz, amom) / (rho * rgrid(ir)) - uphi_basic(ir, iz)
      
         irec = (iz - 1)*nr + ir + 3

         write(lun_tecplot, "(6(1x, e12.5), a1)", rec = irec) rgrid(ir)/L_scale, &
              zgrid(iz)/L_scale, rho_prime, ur/U_scale, uz/U_scale, uphi_prime/U_scale, &
              char(10)
      end do
   end do
end if
close(lun_tecplot)

if (my_node .eq. 0) then
   print *, ' rank 0: subroutine tecplot_fluc_vel_meridional_plane'
   print *, ' wrote tecplot file ', filename
end if

end subroutine tecplot_fluc_vel_meridional_plane

!----------------------------------------------------------------------------------85

