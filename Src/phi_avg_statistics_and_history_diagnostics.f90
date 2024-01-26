!----------------------------------------------------------------------------------85

! phi_avg_statistics_and_history_diagnostics.f90

! Note: If the z-direction is periodic then you can perform a further average
!       in z by calling phi_Reynolds_averages with average_in_z = .true.

!----------------------------------------------------------------------------------85

module phi_avg_stats

real(8), allocatable, dimension(:,:) :: rho_bar, ur_tilde, uz_tilde, uphi_tilde
real(8), allocatable, dimension(:,:) :: T_phiphi, T_phir, T_phiz, T_rr, T_rz, T_zz

real(8), allocatable, dimension(:)   :: T_phiphi_zsum, T_phir_zsum, T_phiz_zsum, &
     T_rr_zsum,     T_rz_zsum,   T_zz_zsum

real(8), allocatable, dimension(:)   :: T_phiphi_rsum, T_phir_rsum, T_phiz_rsum, &
     T_rr_rsum,     T_rz_rsum,   T_zz_rsum
end module phi_avg_stats

!----------------------------------------------------------------------------------85

subroutine phi_Reynolds_averages(istep, t, average_in_z)

! This routine computes and outputs those phi averaged statistics that can be used &
! to compute phi-time Favre averages. This processing is done by
! TOOLS/phi_time_averages.f90.  See notes of 11/24/20 for the algebra on how to
! average up the phi Reynolds stresses to get the Favre stresses.

! p.s.: We have to do this because time average of phi Favre averages is not the
! the phi-time average.

#ifdef mpi_code
   use mpi
#endif
use grid
use partition_data
use partition_data_for_alan
use q_array
use transposes_of_q_and_qdot
use dof_indices
use logical_units
implicit none
integer :: istep
real(8) :: t
logical :: average_in_z

! Local:
#ifdef mpi_code
integer :: status(mpi_status_size), ier
#endif
integer :: ir, iz, iphi, irec, ivar
real(8) :: ur, uphi, uz, rho_uphi_temp ! temps
character(80) :: filename

! The ten phi averages we compute:
real(8), dimension(10, sr:er, sz_phi:ez_phi) :: bar

! Indices for the ten quantities:
integer, parameter :: rho = 1, rho_uphi= 2, rho_ur = 3, rho_uz = 4, &
     rho_uphi_ur = 5, rho_uphi2 = 6, rho_ur2 = 7, rho_uz2 = 8, &
     rho_uphi_uz = 9, rho_ur_uz = 10

! These are for performing z averages (on top of phi averages) if required.
real(8), dimension(10, sr:er) :: bar_z_sum
real(8), dimension(10, nr)    :: bar_z_sum_final
real(8), dimension(10)        :: buffer
integer :: i

character(11), parameter, dimension(10) :: &
     names = ['rho        ', &
              'rho_uphi   ', &
              'rho_ur     ', &
              'rho_uz     ', &
              'rho_uphi_ur', &
              'rho_uphi2  ', &
              'rho_ur2    ', &
              'rho_uz2    ', &
              'rho_uphi_uz', &
              'rho_ur_uz  ']

#ifdef debug_print
   if (my_node .eq. 0) print *, ' First executable in phi_Reynolds_averages'
#endif

call transpose_z_to_phi(ndof, q, q_phi_space)
! Note: q_phi_space (sr:er, sz_phi:ez_phi, ndof, nphi)

#ifdef debug_print
   if (my_node .eq. 0) print *, ' phi_Reynolds_averages: Returned from transpose_z_to_phi'
#endif

bar = 0.d0

do iz = sz_phi, ez_phi
   do ir = sr, er
      do iphi = 1, nphi
         rho_uphi_temp = q_phi_space(ir,iz,amom,iphi) / rgrid(ir)
         ur            = q_phi_space(ir,iz,rmom,iphi) / q_phi_space(ir,iz,irho,iphi)
         uphi          = rho_uphi_temp / q_phi_space(ir,iz,irho,iphi)
         uz            = q_phi_space(ir,iz,zmom,iphi) / q_phi_space(ir,iz,irho,iphi)
         
         bar(rho,         ir,iz) = bar(rho,         ir,iz) + q_phi_space(ir,iz,irho,iphi)
         bar(rho_uphi,    ir,iz) = bar(rho_uphi,    ir,iz) + rho_uphi_temp
         bar(rho_ur,      ir,iz) = bar(rho_ur,      ir,iz) + q_phi_space(ir,iz,rmom,iphi)
         bar(rho_uz,      ir,iz) = bar(rho_uz,      ir,iz) + q_phi_space(ir,iz,zmom,iphi)
         bar(rho_uphi_ur, ir,iz) = bar(rho_uphi_ur, ir,iz) + rho_uphi_temp*ur
         bar(rho_uphi2,   ir,iz) = bar(rho_uphi2,   ir,iz) + rho_uphi_temp*uphi
         bar(rho_ur2,     ir,iz) = bar(rho_ur2,     ir,iz) + q_phi_space(ir,iz,rmom,iphi)*ur
         bar(rho_uz2,     ir,iz) = bar(rho_uz2,     ir,iz) + q_phi_space(ir,iz,zmom,iphi)*uz
         bar(rho_uphi_uz, ir,iz) = bar(rho_uphi_uz, ir,iz) + rho_uphi_temp*uz
         bar(rho_ur_uz,   ir,iz) = bar(rho_ur_uz,   ir,iz) + q_phi_space(ir,iz,rmom,iphi)*uz       
      end do
   end do
end do
bar = bar/nphi

! Output:
! ~~~~~~~
! r, z grid:
if (my_node .eq. 0) then
   open(unit = lun_general_purpose, file = 'rz_grid.dat', form = 'formatted', status = 'unknown')
   write(lun_general_purpose, "(i6, i6)") nr, nz
   do ir = 1, nr
      write(lun_general_purpose, "(e16.9)") rgrid(ir)
   end do
   
   do iz = 1, nz
      write(lun_general_purpose, "(e16.9)") zgrid(iz)      
   end do
end if
close(lun_general_purpose)

write (filename, 1) int(t), t - int(t), istep
1 format ('phi_Reynolds_averages_', i6.6, f0.4, '_', i7.7, '.dat')
! recl is the maximum record length.
open (unit = lun_general_purpose, file = filename, form = 'formatted', status = 'unknown', &
      access = 'direct', recl = 17*10 + 1) ! 10 quantities + char(10)

do iz = sz_phi, ez_phi
   do ir = sr, er
      irec = ir + (iz-1)*nr  
      write(lun_general_purpose, "(10(1x, e16.9), a1)", rec = irec) &
      (bar(ivar, ir, iz), ivar = 1, 10), char(10)
   end do
end do
close(lun_general_purpose)

! This is useful if the z-direction is periodic:
if (average_in_z) then
   bar_z_sum = 0.d0
   do iz = sz_phi, ez_phi
      do ir = sr, er
         do ivar = 1, 10
            bar_z_sum(ivar, ir) = bar_z_sum(ivar, ir) + bar(ivar, ir, iz)
         end do
      end do
   end do

#ifdef mpi_code
   ! tag = ir
   do ir = sr, er
      call mpi_send(bar_z_sum(1, ir), 10, mpi_double, 0, ir, mpi_comm_world, ier)
   end do
   
   if (my_node .eq. 0) then
      bar_z_sum_final = 0.d0
      do ir = 1, nr      
         do i = 1, ng2
            !print *, ' received ir = ', ir, ' ier = ', ier
            call mpi_recv(buffer(1), 10, mpi_double, mpi_any_source, ir, &
                 mpi_comm_world, status, ier)
            !print *, ' received ir = ', ir, ' ier = ', ier
            bar_z_sum_final(:, ir) = bar_z_sum_final(:, ir) + buffer(:)
         end do
      end do
      bar_z_sum_final = bar_z_sum_final / nz
   end if
#else
   bar_z_sum_final = bar_z_sum / nz
#endif

   ! Output:
   if (my_node .eq. 0) then
      write(filename, 2) int(t), t - int(t), istep
      2 format ('profiles_of_phi_z_averages_wrt_r', i6.6, f0.4, '_', i7.7, '.dat')
      open (unit = lun_general_purpose, file = filename, form = 'formatted', status = 'unknown')
      write(lun_general_purpose, *) nr
      do ivar = 1, 10
         write(lun_general_purpose, "(a11)") names(ivar)
         do ir = 1, nr
            write(lun_general_purpose, "(2(1x, e12.5))") rgrid(ir), bar_z_sum_final(ivar, ir)
         end do
      end do
      close(lun_general_purpose)
   end if
end if ! average in z

#ifdef debug_print
   if (my_node .eq. 0) print *, ' Returning from phi_Reynolds_averages'
#endif


end subroutine phi_Reynolds_averages

!----------------------------------------------------------------------------------85

subroutine phi_Favre_averaged_statistics(istep, t)

! Obtains various statistics w.r.t. phi averaging.

use grid
use partition_data
use q_array
use transposes_of_q_and_qdot
use dof_indices
use phi_avg_stats
implicit none
integer :: istep
real(8) :: t

! Local:
integer :: ir, iz, iphi
real(8) :: rho, ur_pp, uz_pp, uphi_pp

! Average density and Favre mean velocities:
allocate(rho_bar   (sr:er, sz_phi:ez_phi))
allocate(ur_tilde  (sr:er, sz_phi:ez_phi))
allocate(uz_tilde  (sr:er, sz_phi:ez_phi))
allocate(uphi_tilde(sr:er, sz_phi:ez_phi))

! Reynolds stresses with Favre fluctuations:
allocate (T_phiphi(sr:er, sz_phi:ez_phi))
allocate (T_phir  (sr:er, sz_phi:ez_phi))
allocate (T_phiz  (sr:er, sz_phi:ez_phi))
allocate (T_rr    (sr:er, sz_phi:ez_phi))
allocate (T_rz    (sr:er, sz_phi:ez_phi))
allocate (T_zz    (sr:er, sz_phi:ez_phi))

call transpose_z_to_phi (ndof, q, q_phi_space)
! Note: q_phi_space (sr:er, sz_phi:ez_phi, ndof, nphi)

! Compute phi Favre averages:
rho_bar    = 0.0d0
ur_tilde   = 0.0d0
uz_tilde   = 0.0d0
uphi_tilde = 0.0d0      
do iz = sz_phi, ez_phi
   do ir = sr, er      
      do iphi = 1, nphi
         rho_bar   (ir, iz) = rho_bar   (ir, iz) + q_phi_space(ir, iz, irho, iphi)
         ur_tilde  (ir, iz) = ur_tilde  (ir, iz) + q_phi_space(ir, iz, rmom, iphi)  
         uz_tilde  (ir, iz) = uz_tilde  (ir, iz) + q_phi_space(ir, iz, zmom, iphi)
         uphi_tilde(ir, iz) = uphi_tilde(ir, iz) + q_phi_space(ir, iz, amom, iphi) / rgrid(ir)         
      end do
      rho_bar   (ir, iz) = rho_bar   (ir, iz) / nphi
      ur_tilde  (ir, iz) = ur_tilde  (ir, iz) / nphi / rho_bar(ir, iz)
      uz_tilde  (ir, iz) = uz_tilde  (ir, iz) / nphi / rho_bar(ir, iz)
      uphi_tilde(ir, iz) = uphi_tilde(ir, iz) / nphi / rho_bar(ir, iz)
   end do
end do

! Stresses, e.g. T_phi_r = <rho uphi_pp ur_pp> where _pp is a Favre fluctuation
T_phiphi = 0.0d0
T_phir   = 0.0d0
T_phiz   = 0.0d0
T_rr     = 0.0d0
T_rz     = 0.0d0
T_zz     = 0.0d0
do iz = sz_phi, ez_phi
   do ir = sr, er                
      do iphi = 1, nphi
         ! Favre fluctuations (double prime):
         rho     = q_phi_space(ir,iz,irho,iphi)
         ur_pp   = q_phi_space(ir,iz,rmom,iphi)/rho           - ur_tilde  (ir,iz)
         uz_pp   = q_phi_space(ir,iz,zmom,iphi)/rho           - uz_tilde  (ir,iz)
         uphi_pp = q_phi_space(ir,iz,amom,iphi)/rho/rgrid(ir) - uphi_tilde(ir,iz)

         T_phiphi(ir, iz) = T_phiphi(ir, iz) + rho*uphi_pp*uphi_pp  
         T_phir  (ir, iz) = T_phir  (ir, iz) + rho*uphi_pp*  ur_pp
         T_phiz  (ir, iz) = T_phiz  (ir, iz) + rho*uphi_pp*  uz_pp
         T_rr    (ir, iz) = T_rr    (ir, iz) + rho*ur_pp  *  ur_pp
         T_rz    (ir, iz) = T_rz    (ir, iz) + rho*ur_pp  *  uz_pp
         T_zz    (ir, iz) = T_zz    (ir, iz) + rho*uz_pp  *  uz_pp
      end do

      T_phiphi(ir, iz) = T_phiphi(ir, iz) / nphi
      T_phir  (ir, iz) = T_phir  (ir, iz) / nphi
      T_phiz  (ir, iz) = T_phiz  (ir, iz) / nphi
      T_rr    (ir, iz) = T_rr    (ir, iz) / nphi
      T_rz    (ir, iz) = T_rz    (ir, iz) / nphi
      T_zz    (ir, iz) = T_zz    (ir, iz) / nphi
   end do
end do

call tecplot_phi_avg_stats(t, istep)

deallocate(rho_bar   )
deallocate(ur_tilde  )
deallocate(uz_tilde  )
deallocate(uphi_tilde)

deallocate (T_phiphi)
deallocate (T_phir  )
deallocate (T_phiz  )
deallocate (T_rr    )
deallocate (T_rz    )
deallocate (T_zz    )

end subroutine phi_Favre_averaged_statistics

!----------------------------------------------------------------------------------85

subroutine tecplot_phi_avg_stats(t, istep)

use grid
use xy_coordinates_of_grid
use logical_units, only: lun_tecplot
use partition_data
use transposes_of_q_and_qdot, only: q_r_space, q_phi_space, qdot_r_space, &
     qdot_phi_space
use phi_avg_stats

implicit none
real(8) :: t
integer :: istep

! Locals:
character(80) :: filename
integer :: ir, iz, irec, nz_plot

! Mean density and Favre mean velocities:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write (filename, 1) int(t), t - int(t), istep
1 format ('phi_Favre_means_', i6.6, f0.4, '_', i7.7, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
     access = 'direct', recl = 13*6 + 1) ! recl = 79, 6 quantities + char(10)

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 2, rec = 1) t, char(10)
2  format (' TITLE = "t = ', e12.5, '"', a1)

   write (lun_tecplot, 3, rec = 2) char(10)
3  format(' VARIABLES="r","z","rho_bar","ur_tilde","uz_tilde","uphi_tilde"',a1)
   
   write (lun_tecplot, 4, rec = 3) nr, nz, char(10)
4  format(' ZONE I=',i4,',',' J=',i4,' DATAPACKING=POINT', a1)
end if

do iz = sz_phi, ez_phi
   do ir = sr, er
      irec = ir + (iz-1)*nr + 3  
      write(lun_tecplot, "(6(1x, e12.5), a1)", rec = irec) &
           rgrid(ir), zgrid(iz), &
           rho_bar (ir, iz), ur_tilde(ir, iz), uz_tilde(ir, iz), uphi_tilde(ir, iz),char(10)
   end do
end do

close(lun_tecplot)

! Favre stresses:
! ~~~~~~~~~~~~~~~
write (filename, 5) int(t), t - int(t), istep
5 format ('phi_Favre_stresses_', i6.6, f0.4, '_', i7.7, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
     access = 'direct', recl = 13*8 + 1) ! recl = 105, 8 quantities + char(10)

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 6, rec = 1) t, char(10)
   6  format (' TITLE = "t = ', e12.5, '"', a1)

   write (lun_tecplot, 7, rec = 2) char(10)
   7  format(' VARIABLES="r","z","T_phiphi","T_phir","T_phiz","T_rr","T_rz","T_zz"',a1)
   
   write (lun_tecplot, 8, rec = 3) nr, nz, char(10)
   8  format(' ZONE I=',i4,',',' J=',i4,' DATAPACKING=POINT', a1)
end if

do iz = sz_phi, ez_phi
   do ir = sr, er
      irec = ir + (iz-1)*nr + 3  
      write(lun_tecplot, "(8(1x, e12.5), a1)", rec = irec) &
           rgrid(ir), zgrid(iz), &
           T_phiphi(ir, iz), T_phir  (ir, iz), T_phiz(ir, iz), T_rr(ir, iz), &
           T_rz    (ir, iz), T_zz    (ir, iz), char(10)
   end do
end do
      
close(lun_tecplot)

end subroutine tecplot_phi_avg_stats

!----------------------------------------------------------------------------------85

subroutine r_phi_average(v, v_ave, all_reduce)

! Performs an r-phi average of the quantity v (dimensioned in z space) as a function
! of z and returns it in v_ave.

! If nphi = 1, then the result will be the r average.

! If all_reduce = .true then all processors will end up with the average.  If .false.
! then only processor 0 will have it.

use partition_data
use dof_indices
use grid
use logical_units
use math_constants
#ifdef mpi_code
   use mpi
#endif

real(8), dimension(sr:er, sphi:ephi, nz) :: v
real(8), dimension(nz) :: v_ave
logical :: all_reduce

! Locals:
real(8), dimension(nz) :: my_integral ! Local to the processor
real(8) :: factor

integer :: ir, iz, iphi, ier

if (nphi .eq. 1) then
   ! To ensure we get the r average:
   factor = 1.d0 / (rmax - rmin)
else
   ! To ensure we get the r phi average:
   factor = dphi / (rmax - rmin) / (phi_max - phi_min)
end if

do iz = 1, nz
   my_integral(iz) = 0.0d0
   do iphi = sphi, ephi
      do ir = sr, er
         my_integral(iz) =  my_integral(iz) + v(ir, iphi, iz)*trap_weight_r(ir)
      end do
   end do
end do

#ifdef mpi_code
   if (all_reduce) then
      call mpi_allreduce(my_integral, v_ave, nz, mpi_double_precision, mpi_sum,    mpi_comm_world, ier)
      do iz = 1, nz
         v_ave(iz) = v_ave(iz)*factor
      end do      
   else
      call mpi_reduce   (my_integral, v_ave, nz, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ier)
      if (my_node .eq. 0) then
         do iz = 1, nz
            ! The formula below will reduce to the r average when nphi = 1 because make_grid.f90
            ! sets dphi = 2*pi in this case.
            v_ave(iz) = v_ave(iz)*factor
         end do
      end if
   end if
#else
   do iz = 1, nz
      v_ave(iz) = my_integral(iz)*factor
   end do
#endif

end subroutine r_phi_average

!----------------------------------------------------------------------------------85

subroutine output_r_phi_average(t, v_ave, prefix)

use partition_data
use logical_units
use grid
use math_constants
implicit none
real(8) :: t
real(8), dimension(nz) :: v_ave
character(17) :: prefix

! Local:
character(40) :: filename
integer :: iz

if (my_node .eq. 0) then
   write(filename, "(a17, i6.6, f0.4, '.pro')") prefix, int(t), t - int(t)
   open (unit = lun_profile(1), file = filename, form = 'formatted', status = 'unknown')

   do iz = 1, nz
      write(lun_profile(1), "(2(1x, e12.5))") zgrid(iz), v_ave(iz)
   end do
   close(lun_profile(1))
end if

end

!----------------------------------------------------------------------------------85

! Routines for history diagnostics.

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
   open(unit = lun_history(4), file = 'ave_fluc_ke_whole_rzphi_tot.his', form = 'formatted', &
        status = 'unknown', access = 'append')                        
   write(lun_history(1), "(4(1x, e12.5))") t, afke_r_middle, afke_z_middle, afke_phi_middle
   write(lun_history(2), "(4(1x, e12.5))") t, afke_r_upper,  afke_z_upper,  afke_phi_upper
   write(lun_history(3), "(3(1x, e12.5))") t, afke_total_middle,  afke_total_upper
   write(lun_history(4), "(5(1x, e12.5))") t, afke_r_middle + afke_r_upper, afke_z_middle + afke_z_upper, &
                                           afke_phi_middle + afke_phi_upper, &
                                           afke_total_middle + afke_total_upper   
   ! In case the code bombs we have the latest history:
   close(lun_history(1)); close(lun_history(2)); close(lun_history(3)); close(lun_history(4))
end if

end subroutine history_write_ave_fluctuation_ke

!----------------------------------------------------------------------------------85










   

