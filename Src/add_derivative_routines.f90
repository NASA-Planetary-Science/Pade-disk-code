!----------------------------------------------------------------------------------85

subroutine add_z_derivatives(q, pressure, qdot)

use grid
use partition_data
use dof_indices, only: irho, amom, zmom, rmom, ener
use thermal_parameters
use boundary_condition_routines
use boundary_condition_types
use boundary_condition_data
use gravity
use artificial_pressure_module
use sbp_z_derivative_module

implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8), dimension(sr:er, sphi:ephi, nz      ) :: pressure
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: qdot

! Locals:
! Later see if these can be made static.
real(8), dimension(sr:er, sphi:ephi, nz) :: F, Fd
real(8), dimension(sr:er, sphi:ephi, nz) :: uz, d_uz_dz
integer :: iphi, iz, ir
real(8) :: gravity_term
real(8) :: c_sound, uz1, lambda_max_local, gamma_sound_speed, u_plus_c

! To store boundary flux derivatives for boundary condition modification.
real(8), dimension(sr:er, sphi:ephi, ndof) :: dF_bottom, dF_top

! The order of calculation is necessitated by the indexing order of "q" which
! in turn is dictated by Alan's transpose routines.

#ifdef debug_print
   if (my_node .eq. 0) print *, ' first executable in add_z_derivatives'
#endif

if (isothermal) then
   gamma_sound_speed = 1.0d0
else
   gamma_sound_speed = gamma
end if

! Store uz which is needed below:
uz(:, :, :) = q(:, :, :, zmom) / q(:, :, :, irho)


! Mass equation:
! ~~~~~~~~~~~~~~
F(:, :, :) = q(:, :, :, zmom) ! rho*uz
call pade_diff_z(nbundle_z, F, Fd)

! Store for possible later modification by boundary conditions:
dF_bottom(:, :, irho) = Fd(:, :, 1 )
dF_top   (:, :, irho) = Fd(:, :, nz)

qdot(:, :, :, irho) = qdot(:, :, :, irho) - Fd(:, :, :)

! Angular momentum equation:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
F(:, :, :) = q(:, :, :, amom) * uz(:, :, :)
call pade_diff_z(nbundle_z, F, Fd)

! Store for possible later modification by non-reflective boundary conditions:
dF_bottom(:, :, amom) = Fd(:, :, 1 )
dF_top   (:, :, amom) = Fd(:, :, nz)

qdot(:, :, :, amom) = qdot(:, :, :, amom) - Fd(:, :, :)

! Radial momentum equation:
! ~~~~~~~~~~~~~~~~~~~~~~~~~
F(:, :, :) = q(:, :, :, rmom) * uz(:, :, :)
call pade_diff_z(nbundle_z, F, Fd)
dF_bottom(:, :, rmom) = Fd(:, :, 1 )
dF_top   (:, :, rmom) = Fd(:, :, nz)
qdot(:, :, :, rmom) = qdot(:, :, :, rmom) - Fd(:, :, :)

! Vertical momentum equation:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Note the addition of pressure here:
F(:, :, :) = pressure(:,:,:) + q(:,:,:,zmom) * uz(:,:,:)
call pade_diff_z(nbundle_z, F, Fd)
dF_bottom(:, :, zmom) = Fd(:, :, 1 )
dF_top   (:, :, zmom) = Fd(:, :, nz)
qdot(:, :, :, zmom) = qdot(:, :, :, zmom) - Fd(:, :, :)

! Internal energy equation:
! ~~~~~~~~~~~~~~~~~~~~~~~~~
if (.not. isothermal) then
   ! Advective term: d/dz(e_int uz):
   F(:, :, :) = q(:, :, :, ener)*uz(:, :, :)
   call pade_diff_z(nbundle_z, F, Fd)

   ! Store for possible later modification by boundary conditions:
   dF_bottom(:, :, ener) = Fd(:, :, 1 )
   dF_top   (:, :, ener) = Fd(:, :, nz)
   
   ! For p div_u term:
   call pade_diff_z(nbundle_z, uz, d_uz_dz)

   ! Add to qdot:
   qdot(:, :, :, ener) = qdot(:, :, :, ener) - Fd(:, :, :) - pressure(:,:,:)*d_uz_dz(:,:,:)
end if

! Non-reflective boundary conditions at the top and bottom boundaries.
! qdot is modified.
if  (zmin_BC .eq. non_reflective) then
   iz = 1
   if (isothermal) then
      call z_isothermal_non_reflective_bc(iz, q, uz, dF_bottom, qdot)
   else
      call z_adiabatic_non_reflective_bc (iz, q, uz, dF_bottom, qdot)
   end if
end if

if  (zmax_BC .eq. non_reflective) then
   iz = nz
   if (isothermal) then
      call z_isothermal_non_reflective_bc(iz, q, uz, dF_top, qdot)
   else
      call z_adiabatic_non_reflective_bc (iz, q, uz, dF_top, qdot)      
   end if
end if
! Done with non-reflective BC 
   
#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' returning from add_z_derivatives'
   end if   
#endif

return
end subroutine add_z_derivatives

!----------------------------------------------------------------------------------85

subroutine add_radial_derivatives (q, pressure, qdot)

! Adds radial derivative terms to qdot.
! q and qdot should be in z space.
! We will perform a transpose to r space and then back to z space here.

use grid, only: ndof, nr, nz, Ji_r, rgrid, dphi, Ji_r, dr
use partition_data
use dof_indices, only: irho, rmom, amom, zmom, ener
use thermal_parameters
use transposes_of_q_and_qdot
use artificial_pressure_module
use boundary_condition_types
! So I can call enforce_BC
use boundary_condition_routines
use boundary_condition_data
use thermal_parameters
implicit none

! Note the dimensioning:
! q_r_space(sphi:ephi, sz_r:ez_r, ndof, nr)

real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in)  :: q
real(8), dimension(sr:er, sphi:ephi, nz      ), intent(in)  :: pressure
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(out) :: qdot

! Local:
integer :: ir, iz, iphi, idof

! Again, see later if these can be made static.
real(8), dimension(sphi:ephi, sz_r:ez_r, ndof, nr) :: F, dFdr

! Needed since using the internal energy:
real(8), dimension(sphi:ephi, sz_r:ez_r, nr) :: r_ur, div_u_term

! This is used so often that I decided to make it an array:
real(8), dimension(sphi:ephi, sz_r:ez_r, nr) :: ur

real(8) :: c_sound, lambda_max_local, gamma_sound_speed

#ifdef debug_print
   ! Note: my_node = 0 for the serial code so this works both cases.
   if (my_node .eq. 0) then
      print *, ' add_radial_derivatives has been called'
   end if
#endif

! Perform transposes:
call transpose_z_to_r (ndof, q,        q_r_space   )
call transpose_z_to_r (ndof, qdot,     qdot_r_space)
call transpose_z_to_r (1,    pressure, p_r_space   )

if (isothermal) then
   gamma_sound_speed = 1.0d0
else
   gamma_sound_speed = gamma
end if

! Store ur:
do ir = 1, nr
   do iz = sz_r, ez_r
      do iphi = sphi, ephi
          ur(iphi, iz, ir) = q_r_space(iphi, iz, rmom, ir) / q_r_space(iphi, iz, irho, ir)
      end do
   end do
end do

! Compute fluxes.  On Oct 25, 2018 I put in the pressure into the flux in order to
! implement the non-reflective condition.  This means that I need a p/r term on the rhs.
! This is done below.
do ir = 1, nr
   do iz = sz_r, ez_r
      do iphi = sphi, ephi
         F(iphi,iz,irho,ir) = rgrid(ir)* q_r_space(iphi,iz,rmom,ir)
         F(iphi,iz,amom,ir) = rgrid(ir)* q_r_space(iphi,iz,amom,ir)*ur(iphi,iz,ir)
         F(iphi,iz,rmom,ir) = rgrid(ir)*(q_r_space(iphi,iz,rmom,ir)*ur(iphi,iz,ir) + p_r_space(iphi,iz,ir))
         F(iphi,iz,zmom,ir) = rgrid(ir)* q_r_space(iphi,iz,zmom,ir)*ur(iphi,iz,ir)
      end do
   end do
end do

if (.not. isothermal) then
   ! This is for the internal energy equation:
   do ir = 1, nr
      do iz = sz_r, ez_r
         do iphi = sphi, ephi
            F(iphi, iz, ener,  ir) = q_r_space(iphi,iz,ener,ir) * rgrid(ir) * ur(iphi,iz,ir)
            ! We will differentiate this below:
            r_ur(iphi, iz, ir) = ur(iphi,iz,ir) * rgrid(ir)
         end do
      end do
   end do
end if

! Differentiate the fluxes:
! Note that for the isothermal case we are differentiating one more dof
! than needed.  This can't be helped due to the indexing of F we have chosen.   
call pade_diff_bundle(mphi*mz_r*ndof, nr, Ji_r, F, dFdr)

! Non-reflective boundary conditions at the left and right boundaries.
! dFdr is modified:
if  (rmin_BC .eq. non_reflective) then
   if (isothermal) then
      ir = 1
      ! print *, ' calling r_isothermal_non_reflective_bc'
      call r_isothermal_non_reflective_bc(ir, q_r_space, ur, dFdr)
   end if
end if

if  (rmax_BC .eq. non_reflective) then
   if (isothermal) then
      ir = nr
      call r_isothermal_non_reflective_bc(ir, q_r_space, ur, dFdr)
   end if
end if
! Done with non-reflective BC 
   
   
if (.not.isothermal) then
   ! Calculate the derivative of r*ur which is needed for div u:
   call pade_diff_bundle(mphi*mz_r, nr, Ji_r, r_ur, div_u_term)
end if   

! Note: All the r derivatives are divided by r:
do ir = 1, nr
   do idof = 1, ndof
      do iz = sz_r, ez_r
         do iphi = sphi, ephi
            qdot_r_space(iphi, iz, idof, ir) = qdot_r_space(iphi, iz, idof, ir) - &
                                       dFdr(iphi, iz, idof, ir)/rgrid(ir)
         end do
      end do
   end do
end do

! p/r term since you now have p + rho * ur**2 in the radial momentum equation:
do ir = 1, nr
   do iz = sz_r, ez_r
      do iphi = sphi, ephi
         qdot_r_space(iphi,iz,rmom,ir) = qdot_r_space(iphi,iz,rmom,ir) + p_r_space(iphi,iz,ir)/rgrid(ir)
      end do
   end do
end do

! - p div.u term in internal energy equation
if (.not.isothermal) then
   do ir = 1, nr
      do iz = sz_r, ez_r
         do iphi = sphi, ephi
            ! div_u_term is ddr(r*ur)
            qdot_r_space(iphi, iz, ener, ir) = qdot_r_space(iphi, iz, ener, ir) - &
                                       p_r_space(iphi, iz, ir)/rgrid(ir) * div_u_term(iphi, iz, ir)
         end do
      end do
   end do
end if

! Transpose qdot back to z space:
call transpose_r_to_z (ndof, qdot_r_space, qdot)

return
end subroutine add_radial_derivatives

!----------------------------------------------------------------------------------85

subroutine add_phi_derivatives(q, qdot, pressure, lambda_max_extra_operator)

! Adds phi derivative terms to qdot.
! Both q and qdot should have z-space indexing.  We will perform a transpose to phi
! space here (and then transpose to z-space at the end).

! Note: The extra FARGO terms are treated at the end via a call to
! subroutine add_fargo_extra_operator_terms.

use transposes_of_q_and_qdot
use grid
use partition_data
use dof_indices, only: irho, zmom, rmom, amom, ener
use thermal_parameters
use fargo_or_plotting_shift
implicit none

real(8), dimension(sr:er, sz_phi:ez_phi, ndof, nphi), intent(in)    :: q
real(8), dimension(sr:er, sz_phi:ez_phi, ndof, nphi), intent(inout) :: qdot
real(8), dimension(sr:er, sz_phi:ez_phi, ndof      ), intent(in)    :: pressure
real(8),                                              intent(out)   :: lambda_max_extra_operator

! Local:
integer :: ir, iz, iphi, idof

! Fluxes and their phi derivatives in phi space dimensioning:
! Again, see later if these can be made static.
real(8), dimension(sr:er, sz_phi:ez_phi, ndof, nphi) :: F, dF
real(8) :: c_sound, lambda_max_local, gamma_sound_speed

real(8), dimension(sr:er, sz_phi:ez_phi, nphi) :: uphi

! Residual azimuthal velocity for Fargo trick:
real(8), dimension(sr:er, sz_phi:ez_phi, nphi) :: uphi_prime

! Needed for the internal energy.
real(8), dimension(sr:er, sz_phi:ez_phi, nphi) :: div_u_term

! Artificial pressure related local:
real(8) :: beta_ap

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' add_phi_derivatives has been called'
   end if
#endif

! Transposes:
#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' Calling transpose_z_to_phi for q in add_phi_derivatives'
   end if
#endif
call transpose_z_to_phi(ndof, q,        q_phi_space)
call transpose_z_to_phi(ndof, qdot,     qdot_phi_space)
call transpose_z_to_phi(1,    pressure, p_phi_space)   

if (isothermal) then
   gamma_sound_speed = 1.0d0
else
   gamma_sound_speed = gamma
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' add_phi_derivatives: about to compute pressure'
#endif

! Store uphi:
do iphi = 1, nphi
   do iz = sz_phi, ez_phi
      do ir = sr, er
         uphi(ir,iz,iphi) = q_phi_space(ir,iz,amom,iphi)/q_phi_space(ir,iz,irho,iphi)/rgrid(ir)
      end do
   end do
end do

! Compute fluxes:
if (apply_fargo_this_step) then
   ! Fargo:
   do iphi = 1, nphi
      do iz = sz_phi, ez_phi
         do ir = sr, er
            uphi_prime(ir,iz,iphi) =  uphi(ir,iz,iphi) - uphi_fargo_subtract(ir)

            ! rho*uphi'
            F(ir, iz, irho, iphi) =  q_phi_space(ir, iz, irho, iphi) * uphi_prime(ir,iz,iphi)

            ! rho*uphi*uphi' + p
            F(ir, iz, amom, iphi) =  q_phi_space(ir, iz, amom, iphi)/rgrid(ir) * uphi_prime(ir,iz,iphi) + &
                 p_phi_space(ir, iz, iphi)    ! On Aug. 28, 2023 I noticed that this term was missing.
            ! rho*ur*uphi'
            F(ir, iz, rmom, iphi) =  q_phi_space(ir, iz, rmom, iphi) * uphi_prime(ir,iz,iphi)
            ! rho*uz*uphi'
            F(ir, iz, zmom, iphi) =  q_phi_space(ir, iz, zmom, iphi) * uphi_prime(ir,iz,iphi)
         end do
      end do
   end do

   ! Fargo: Energy equation:
   if (.not. isothermal) then
      do iphi = 1, nphi
         do iz = sz_phi, ez_phi
            do ir = sr, er
               F(ir, iz, ener,  iphi) = q_phi_space(ir, iz, ener, iphi) * uphi_prime(ir,iz,iphi)          
            end do
         end do
      end do
   end if
else
   
#ifdef debug_print
   if (my_node .eq. 0) print *, ' add_phi_derivatives: about to compute non-fargo fluxes'
#endif
   ! No Fargo:

   ! Non-energy fluxes:
   do iphi = 1, nphi
      do iz = sz_phi, ez_phi
         do ir = sr, er
            F(ir, iz, irho, iphi) =  q_phi_space(ir, iz, irho, iphi) * uphi(ir, iz, iphi)
            ! rho uphi^2 + p
            F(ir, iz, amom, iphi) =  q_phi_space(ir, iz, irho, iphi) * uphi(ir, iz, iphi)**2 + &
                                     p_phi_space(ir, iz, iphi)
            F(ir, iz, rmom, iphi) =  q_phi_space(ir, iz, rmom, iphi) * uphi(ir, iz, iphi)
            F(ir, iz, zmom, iphi) =  q_phi_space(ir, iz, zmom, iphi) * uphi(ir, iz, iphi)
         end do
      end do
   end do

   ! Energy fluxes for non-Fargo:
   if (.not. isothermal) then
      do iphi = 1, nphi
         do iz = sz_phi, ez_phi
            do ir = sr, er
               F(ir, iz, ener, iphi) = q_phi_space(ir, iz, ener, iphi) * uphi(ir, iz, iphi)
            end do
         end do
      end do
   end if
end if 

! Differentiate the fluxes in theta.  Note: We have an extra un-needed one for the isothermal case:
call pade_diff_periodic(mr*mz_phi*ndof, nphi, dphi, F(sr, sz_phi, 1, 1), dF(sr, sz_phi, 1, 1))

do iphi = 1, nphi
   do iz = sz_phi, ez_phi
      do ir = sr, er
         qdot_phi_space(ir,iz,irho,iphi) = qdot_phi_space(ir,iz,irho,iphi) - dF(ir,iz,irho,iphi)/rgrid(ir)
         qdot_phi_space(ir,iz,rmom,iphi) = qdot_phi_space(ir,iz,rmom,iphi) - dF(ir,iz,rmom,iphi)/rgrid(ir)
         qdot_phi_space(ir,iz,zmom,iphi) = qdot_phi_space(ir,iz,zmom,iphi) - dF(ir,iz,zmom,iphi)/rgrid(ir)
         qdot_phi_space(ir,iz,amom,iphi) = qdot_phi_space(ir,iz,amom,iphi) - dF(ir,iz,amom,iphi) ! No r
      end do
   end do
end do

if (.not. isothermal) then
   do iphi = 1, nphi
      do iz = sz_phi, ez_phi
         do ir = sr, er
            qdot_phi_space(ir,iz,ener,iphi) = qdot_phi_space(ir,iz,ener,iphi) - dF(ir,iz,ener,iphi)/rgrid(ir)
         end do
      end do
   end do
end if

! - p div u term for internal energy equation:
if (.not.isothermal) then
   ! div_u_term is simply d_uphi_dphi:
   call pade_diff_periodic(mr*mz_phi, nphi, dphi, uphi, div_u_term)
   do iphi = 1, nphi
      do iz = sz_phi, ez_phi
         do ir = sr, er
            qdot_phi_space(ir,iz,ener,iphi) = qdot_phi_space(ir,iz,ener,iphi) - &
                 p_phi_space(ir,iz,iphi)/rgrid(ir)*div_u_term(ir,iz,iphi)
         end do
      end do
   end do
end if                           

if (add_fargo_extra_operator_now) then
   ! This subroutine is located in fargo_and_plotting_shift.f90
   ! Sept. 21, 2023: I now pass p_phi_space here instead of computing it in
   ! subroutine add_fargo_extra_operator_terms:
   call add_fargo_extra_operator_terms (q_phi_space, qdot_phi_space, p_phi_space, &
              lambda_max_extra_operator)
else
   lambda_max_extra_operator = 0.d0
end if
call transpose_phi_to_z(ndof, qdot_phi_space, qdot)

end subroutine add_phi_derivatives

!----------------------------------------------------------------------------------85

subroutine output_midplane_radial_force_balance_diagnostic(q, t, L_scale, T_scale)

! Assumes no Fargo.

use grid
use partition_data
use dof_indices
use thermal_parameters
use transposes_of_q_and_qdot
use gravity
use logical_units
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in)  :: q
real(8) :: t, L_scale, T_scale

! Local:
integer :: ir, iz, iphi, iz1, iz2
real(8), dimension(sphi:ephi, sz_r:ez_r, nr):: p, dpdr
real(8) :: t_scaled, a_scale
real(8) :: pressure_force_1, pressure_force_2, pressure_force, uphi1, uphi2, uphi, &
         centrifugal_force_1, centrifugal_force_2, centrifugal_force, &
         g_force_1, g_force_2, g_force, sum, sum_no_p
character(80) :: filename

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' output_midplane_radial_force_balance_diagnostic has been called'
   end if
#endif

call transpose_z_to_r(ndof, q, q_r_space)

if (isothermal) then
   do ir = 1, nr
      do iz = sz_r, ez_r
         do iphi = sphi, ephi
            p(iphi, iz, ir) = q_r_space(iphi, iz, irho, ir) * ci_squared_initial(ir, iz)
         end do
      end do
   end do
else
   do ir = 1, nr
      do iz = sz_r, ez_r
         do iphi = sphi, ephi
            p(iphi, iz, ir) = q_r_space(iphi, iz, ener, ir) * gm1 ! from internal energy
         end do
      end do
   end do
end if

call pade_diff_bundle(mphi*mz_r, nr, Ji_r, p, dpdr)

t_scaled = t / T_scale
a_scale  = L_scale / T_scale**2 ! force scale per unit mass
iphi = 1

if (mod(nz, 2) .eq. 0) then
!!$   iz1 = nz/2
!!$   iz2 = iz1 + 1
!!$   if ((iz1 .ge. sz_r) .and. (iz1 .le. ez_r) .and. (iz2 .ge. sz_r) .and. (iz2 .le. ez_r)) then
!!$      write(filename, "('radial_balance_t_', i6.6, f0.4, '.dat')") int(t_scaled), t_scaled - int(t_scaled)
!!$      open(unit = lun_general_purpose_1, file = filename, form = 'formatted', status = 'unknown')
!!$      do ir = 1, nr
!!$         pressure_force_1 = -dpdr(iphi, iz1, ir) / q_r_space(iphi, iz1, irho, ir)
!!$         pressure_force_2 = -dpdr(iphi, iz2, ir) / q_r_space(iphi, iz2, irho, ir)
!!$         pressure_force = (pressure_force_1 + pressure_force_2) / 2.d0 / a_scale
!!$
!!$         uphi1 = q_r_space(iphi, iz1, amom, ir) / q_r_space(iphi, iz1, irho, ir) / rgrid(ir)
!!$         uphi2 = q_r_space(iphi, iz2, amom, ir) / q_r_space(iphi, iz2, irho, ir) / rgrid(ir)         
!!$         centrifugal_force_1 = uphi1**2 / rgrid(ir)
!!$         centrifugal_force_2 = uphi2**2 / rgrid(ir)
!!$         centrifugal_force = (centrifugal_force_1 + centrifugal_force_2) / 2.d0 / a_scale
!!$         
!!$         g_force_1 = gr(ir, iz1)
!!$         g_force_2 = gr(ir, iz1)
!!$         g_force = (g_force_1 + g_force_2) / 2.d0 / a_scale
!!$
!!$         sum      = pressure_force + centrifugal_force + g_force
!!$         sum_no_p = centrifugal_force + g_force
!!$         
!!$         write(lun_general_purpose_1, "(5(1x, e12.5))") pressure_force, centrifugal_force, g_force, sum, sum_no_p
   iz = nz / 2
   if ((iz .ge. sz_r) .and. (iz .le. ez_r)) then
      write(filename, "('radial_balance_t_', i6.6, f0.4, '.dat')") int(t_scaled), t_scaled - int(t_scaled)
      open(unit = lun_general_purpose_1, file = filename, form = 'formatted', status = 'unknown')
      do ir = 1, nr
         pressure_force = -dpdr(iphi, iz, ir) / q_r_space(iphi, iz, irho, ir) / a_scale

         uphi = q_r_space(iphi, iz, amom, ir) / q_r_space(iphi, iz, irho, ir) / rgrid(ir)
         centrifugal_force = uphi**2 / rgrid(ir) / a_scale
         
         g_force = gr(ir, iz) / a_scale

         sum = pressure_force + centrifugal_force + g_force
         sum_no_p = centrifugal_force + g_force         
         
         write(lun_general_purpose_1, "(6(1x, e12.5))") rgrid(ir)/L_scale, &
              pressure_force, centrifugal_force, g_force, sum, sum_no_p
      end do
      close(lun_general_purpose_1)
   end if ! I have iz1 and iz2
else
   iz = (nz + 1) / 2
   if ((iz .ge. sz_r) .and. (iz .le. ez_r)) then
      write(filename, "('radial_balance_t_', i6.6, f0.4, '.dat')") int(t_scaled), t_scaled - int(t_scaled)
      open(unit = lun_general_purpose_1, file = filename, form = 'formatted', status = 'unknown')
      do ir = 1, nr
         pressure_force = -dpdr(iphi, iz, ir) / q_r_space(iphi, iz, irho, ir) / a_scale

         uphi = q_r_space(iphi, iz, amom, ir) / q_r_space(iphi, iz, irho, ir) / rgrid(ir)
         centrifugal_force = uphi**2 / rgrid(ir) / a_scale
         
         g_force = gr(ir, iz) / a_scale

         sum = pressure_force + centrifugal_force + g_force
         sum_no_p = centrifugal_force + g_force         
         
         write(lun_general_purpose_1, "(6(1x, e12.5))") rgrid(ir)/L_scale, &
              pressure_force, centrifugal_force, g_force, sum, sum_no_p
      end do
      close(lun_general_purpose_1)
   end if ! I have iz  
end if

return
end subroutine output_midplane_radial_force_balance_diagnostic

!----------------------------------------------------------------------------------85
