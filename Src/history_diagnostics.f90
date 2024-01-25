!----------------------------------------------------------------------------------85

subroutine history_write_rms(t, output_profiles)

use grid
use partition_data
use q_array
use dof_indices
use basic_state
use phi_avg_stats
use logical_units
implicit none
real(8) :: t
real(8) :: uz_rms, ur_rms, uphi_rms
logical :: output_profiles

! Local:
real(8), dimension(sr:er, sphi:ephi, nz) :: v
real(8), dimension(nz) :: v_ave
character(17) :: prefix
integer :: ir, iz, iphi
real(8) :: uphi, uphi_prime
logical :: all_reduce = .false.

if (output_profiles) then
   prefix = 'uz_rphi_ave____t_'
   v(:, :, :) = q(:, :, :, zmom) / q(:, :, :, irho)
   call r_phi_average(v, v_ave, all_reduce)
   call output_r_phi_average(t, v_ave, prefix)
end if

if (output_profiles) then
   prefix = 'ur_rphi_ave____t_'
   v(:, :, :) = q(:, :, :, rmom) / q(:, :, :, irho)
   call r_phi_average(v, v_ave, all_reduce)
   call output_r_phi_average(t, v_ave, prefix)
end if

if (output_profiles) then
   prefix = 'uphi_rphi_ave__t_'
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            v(ir, iphi, iz) = q(ir, iphi, iz, rmom) / q(ir, iphi, iz, irho) / rgrid(ir)
         end do
      end do
   end do
   call r_phi_average(v, v_ave, all_reduce)
   call output_r_phi_average(t, v_ave, prefix)
end if

! uz^2:
v(:, :, :) = (q(:, :, :, zmom) / q(:, :, :, irho))**2
call r_phi_average(v, v_ave, all_reduce)
if (output_profiles) then
  prefix = 'uz2_rphi_ave____t_'
  call output_r_phi_average(t, v_ave, prefix)
end if
! Average in z:
if (my_node .eq. 0) then
   uz_rms = SQRT(sum(v_ave, nz) / nz)
end if

! ur^2
v(:, :, :) = (q(:, :, :, rmom) / q(:, :, :, irho))**2
call r_phi_average(v, v_ave, all_reduce)
! Average in z:
if (my_node .eq. 0) then
   ur_rms = SQRT(sum(v_ave, nz) / nz)
end if

! uphi prime ^2
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         uphi = q(ir, iphi, iz, amom) / q(ir, iphi, iz, irho) / rgrid(ir)
         uphi_prime = uphi - uphi_basic(ir, iz)
         v(ir, iphi, iz) = uphi_prime**2
      end do
   end do
end do
call r_phi_average(v, v_ave, all_reduce)
! Average in z:
if (my_node .eq. 0) then
   uphi_rms = SQRT(sum(v_ave, nz) / nz)
end if

if (my_node .eq. 0) then         
   open(unit = lun_history(1), file = 'vel_rms_rzphi.his', form = 'formatted', &
        status = 'unknown', access = 'append')
   write(lun_history(1), "(4(1x, e12.5))") t, ur_rms, uz_rms, uphi_rms
   ! In case the code bombs we have the latest history:
   close(lun_history(1))
end if


end subroutine history_write_rms

!----------------------------------------------------------------------------------85

subroutine history_write_ave_fluctuation_ke(t)

! This routine was written specially as a diagnostic for VSI runs and works only if
! nz - 1 is divisible by 4.  It also assumes a uniform mesh in z.

! afke = average fluctuation kinetic energy.  This is computed in two-parts: the
! middle of the disk and the upper part of the disk.
! fke_z = average of 1/2 rho uz'^2 where uz' = uz - uz basic state
! and similarly for the other components.
! We assume that the basic state has only a phi component.

use grid
use partition_data
use q_array
use dof_indices
use basic_state
use phi_avg_stats
use logical_units
implicit none
real(8), intent(in) :: t

! Local:
real(8) :: afke_z_middle, afke_r_middle, afke_phi_middle, &
           afke_z_upper, afke_r_upper, afke_phi_upper, &
           afke_total_middle, afke_total_upper

real(8), dimension(sr:er, sphi:ephi, nz) :: v
real(8), dimension(nz) :: v_ave
real(8), dimension(4) :: afke_z, afke_r, afke_phi
integer :: ir, iz, iphi, icut, ipart
integer, dimension(5) :: iz_cut
real(8) :: uphi, uphi_prime
! Only processor 0 will end up with the r phi average:
logical :: all_reduce = .false.
integer :: delta_nz

delta_nz = (nz - 1) / 4
do icut = 1, 5
   iz_cut(icut) = 1 + (icut - 1)*delta_nz
end do

! 1/2 rho uz^2:
v(:, :, :) = 0.5d0 * q(:, :, :, zmom)**2 / q(:, :, :, irho)
call r_phi_average(v, v_ave, all_reduce)
! Average in z:
if (my_node .eq. 0) then
   do ipart = 1, 4
      afke_z(ipart) = 0.d0
      do iz = iz_cut(ipart), iz_cut(ipart+1)
         afke_z(ipart) = afke_z(ipart) + v_ave(iz)
      end do
     afke_z(ipart) = afke_z(ipart) / delta_nz      
   end do
end if

! 1/2 rho ur^2
v(:, :, :) = q(:, :, :, rmom)**2 / q(:, :, :, irho)
call r_phi_average(v, v_ave, all_reduce)
! Average in z:
if (my_node .eq. 0) then
   do ipart = 1, 4
      afke_r(ipart) = 0.d0
      do iz = iz_cut(ipart), iz_cut(ipart+1)
         afke_r(ipart) = afke_r(ipart) + v_ave(iz)
      end do
      afke_r(ipart) = afke_r(ipart) / delta_nz      
   end do
end if

! 1/2 rho uphi prime ^2
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         uphi = q(ir, iphi, iz, amom) / q(ir, iphi, iz, irho) / rgrid(ir)
         uphi_prime = uphi - uphi_basic(ir, iz)
         v(ir, iphi, iz) = q(ir, iphi, iz, irho) * uphi_prime**2
      end do
   end do
end do
call r_phi_average(v, v_ave, all_reduce)
! Average in z:
if (my_node .eq. 0) then
   do ipart = 1, 4
      afke_phi(ipart) = 0.d0
      do iz = iz_cut(ipart), iz_cut(ipart+1)
         afke_phi(ipart) = afke_phi(ipart) + v_ave(iz)
      end do
      afke_phi(ipart) = afke_phi(ipart) / delta_nz      
   end do
end if

! Separate the parts into a middle part and an upper part:
afke_z_middle   = afke_z  (2) + afke_z  (3)
afke_r_middle   = afke_r  (2) + afke_r  (3)
afke_phi_middle = afke_phi(2) + afke_phi(3)

afke_z_upper   = afke_z  (1) + afke_z  (4)
afke_r_upper   = afke_r  (1) + afke_r  (4)
afke_phi_upper = afke_phi(1) + afke_phi(4)

afke_total_middle = afke_z_middle + afke_r_middle + afke_phi_middle
afke_total_upper  = afke_z_upper  + afke_r_upper  + afke_phi_upper 

if (my_node .eq. 0) then         
   open(unit = lun_history(1), file = 'ave_fluc_ke_middle_rzphi.his', form = 'formatted', &
        status = 'unknown', access = 'append')
   open(unit = lun_history(2), file = 'ave_fluc_ke_upper_rzphi.his', form = 'formatted', &
        status = 'unknown', access = 'append')
   open(unit = lun_history(3), file = 'ave_fluc_ke_total_mid_and_up.his', form = 'formatted', &
        status = 'unknown', access = 'append')
   open(unit = lun_history(4), file = 'ave_fluc_ke_whole_rzphi.his', form = 'formatted', &
        status = 'unknown', access = 'append')                        
   write(lun_history(1), "(4(1x, e12.5))") t, afke_r_middle, afke_z_middle, afke_phi_middle
   write(lun_history(2), "(4(1x, e12.5))") t, afke_r_upper,  afke_z_upper,  afke_phi_upper
   write(lun_history(3), "(3(1x, e12.5))") t, afke_total_middle,  afke_total_upper
   write(lun_history(4), "(4(1x, e12.5))") t, afke_r_middle + afke_r_upper, afke_z_middle + afke_z_upper, &
                                              afke_phi_middle + afke_phi_upper   
   ! In case the code bombs we have the latest history:
   close(lun_history(1)); close(lun_history(2)); close(lun_history(3)); close(lun_history(4))
end if

end subroutine history_write_ave_fluctuation_ke

!----------------------------------------------------------------------------------85
