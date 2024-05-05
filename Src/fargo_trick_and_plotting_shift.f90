!----------------------------------------------------------------------------------85

module fargo_or_plotting_shift
   logical :: apply_fargo, apply_fargo_extra_operator, integer_shifts, &
        apply_plotting_shift, add_fargo_extra_operator_now, &
        apply_fargo_this_step, have_Omega_bar
   real(8), allocatable, dimension(:)   :: Omega_bar   ! continuous
   real(8), allocatable, dimension(:)   :: Omega_fargo ! not-necessarily continuous
   real(8), allocatable, dimension(:)   :: uphi_fargo_subtract
   integer, allocatable, dimension(:)   :: nshift
   integer, allocatable, dimension(:,:) :: inew_phi ! used for index shifting
   ! Used for chain rule correction terms.  See notes of Feb. 10, 2018.   
   real(8), allocatable, dimension(:)   :: d_Omega_bar_dr
   ! Fargo factor = - dOmegabar/dr (t - t0)
   real(8), allocatable, dimension(:)   :: fargo_factor, fargo_factor_over_r  
   complex(8), allocatable, dimension(:,:) :: complex_shift_factor ! (nphi/2+1, sr:er)
   real(8) :: Omega0_plotting_shift, t_of_previous_shift
   ! Needed because the phi domain size need not be twopi:
   real(8) :: twopi_over_Delta_phi_domain
   logical fft_has_been_initialized
end module fargo_or_plotting_shift

module fftw_for_fargo
   integer(8) :: forward_plan, inverse_plan   
end module fftw_for_fargo

!----------------------------------------------------------------------------------85

subroutine setup_Omega_bar_for_fargo(gravity_flag, GM)

use grid
use fargo_or_plotting_shift
use total_allocated_words, only: n_words_allocated
implicit none
logical, intent(in) :: gravity_flag
real(8), intent(in) :: GM

integer :: ir

allocate(Omega_bar(nr), d_Omega_bar_dr(nr))
n_words_allocated = n_words_allocated + 2*nr

if (gravity_flag) then
   do ir = 1, nr
      Omega_bar(ir)      =          SQRT(GM / rgrid(ir)**3)
      ! Needed for fargo correction terms.
      d_Omega_bar_dr(ir) = -1.5d0 * SQRT(GM / rgrid(ir)**5)
   end do
else
   do ir = 1, nr
      Omega_bar(ir)      = 0.d0
      d_Omega_bar_dr(ir) = 0.d0
   end do
end if

have_Omega_bar = .true.

end subroutine setup_Omega_bar_for_fargo

!----------------------------------------------------------------------------------85

subroutine get_fargo_shifts_and_uphi_fargo_subtract(t, dt)

use grid
use fargo_or_plotting_shift
use partition_data
implicit none
real(8) :: t, dt

! Local:
integer :: ir, iphi, ir_interval
real(8) :: actual_phi_shift
real(8) :: interval_slope(nr)
complex(8) :: ci = CMPLX(0.0, 1.0d0)
integer :: k
real(8), dimension(nr) :: desired_phi_shift, phase_factor

! Note: _bar = continuous; _fargo : not-necessarily so.

if (integer_shifts) then
   do ir = 1, nr
      ! Omega_bar is set-up by subroutine setup_Omega_fargo
      desired_phi_shift(ir) = Omega_bar(ir) * dt
      nshift(ir)            = INT(desired_phi_shift(ir) / dphi + 0.5d0)
      actual_phi_shift      = nshift(ir) * dphi
      Omega_fargo(ir)       = actual_phi_shift / dt
      ! This is what will be subtracted in the phi advection terms.
      uphi_fargo_subtract(ir) = Omega_fargo(ir) * rgrid(ir)  
   end do
else
   do ir = 1, nr
      desired_phi_shift(ir) = Omega_bar(ir) * dt
      Omega_fargo(ir) = Omega_bar(ir)
      uphi_fargo_subtract(ir)  = Omega_fargo(ir) * rgrid(ir)
   end do

   do ir = 1, nr
      phase_factor(ir) = desired_phi_shift(ir) * twopi_over_Delta_phi_domain
   end do

   do ir = sr, er
      do k = 1, nphi/2+1
         complex_shift_factor(k, ir) = EXP(-ci * (k-1) * phase_factor(ir))
      end do
   end do
end if

!!$do ir_interval = 1, nr - 1
!!$   interval_slope(ir) = (Omega_bar(ir + 1) - Omega_bar(ir)) / (rgrid(ir + 1) - rgrid(ir))
!!$end do
!!$
!!$d_Omega_bar_dr(1) = 0.0d0
!!$do ir = 2, nr
!!$   d_Omega_bar_dr(ir) = 0.5d0*interval_slope(ir) + 0.5d0*interval_slope(ir-1)
!!$end do

!!$do ir = 1, nr
!!$   print *, ' ir = ', ir, ' nshift = ', nshift(ir)
!!$end do
!!$read (5, *)

do iphi = 1, nphi
   do ir = 1, nr
      inew_phi(ir, iphi) = mod(iphi + nshift(ir), nphi)
      if (inew_phi(ir, iphi) .eq. 0) inew_phi(ir, iphi) = nphi
   end do
end do

end subroutine get_fargo_shifts_and_uphi_fargo_subtract

!----------------------------------------------------------------------------------85

subroutine apply_fargo_shifts(q)

use grid
use partition_data
use fargo_or_plotting_shift
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

if (integer_shifts) then
   call apply_integer_shifts(q)
else
#ifdef fftw
   ! call apply_real_shifts_using_fftw(q)
   call apply_real_shifts_using_fftw_old(q) ! this can actually be cheaper.
#else
   call apply_real_shifts_using_rogallo_fft(q)
#endif
end if

end subroutine apply_fargo_shifts

!----------------------------------------------------------------------------------85

#ifdef fftw
subroutine make_plans_for_fftw(q_phi_space)

use fftw_for_fargo
use grid
use partition_data
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

! Argument:
real(8), dimension(sr:er, sz_phi:ez_phi, ndof, nphi) :: q_phi_space

! Local:
integer, parameter :: rank = 1
integer :: num_transforms
integer :: stride_i, stride_o, dist_i, dist_o, n_embed_i, n_embed_o

! Make this global at at some point:
complex(8), dimension(sr:er, sz_phi:ez_phi, ndof, nphi/2+1):: q_phi_space_hat

num_transforms = mr*mz_phi*ndof
stride_i       = num_transforms
stride_o       = stride_i
dist_i         = 1
dist_o         = 1
n_embed_i      = nphi
n_embed_o      = nphi

call dfftw_plan_many_dft_r2c(forward_plan, rank, nphi, num_transforms, &
     q_phi_space,     n_embed_i, stride_i, dist_i, &
     q_phi_space_hat, n_embed_o, stride_o, dist_o, FFTW_ESTIMATE)

! Note: The i's and o'e have been reversed compared to above.
call dfftw_plan_many_dft_c2r(inverse_plan, rank, nphi, num_transforms, &
     q_phi_space_hat, n_embed_o, stride_o, dist_o, &
     q_phi_space,     n_embed_i, stride_i, dist_i, FFTW_ESTIMATE)

end subroutine make_plans_for_fftw

!----------------------------------------------------------------------------------85

subroutine apply_real_shifts_using_fftw(q)

! Does a real valued shift in phi.

use fftw_for_fargo
use grid
use partition_data
use fargo_or_plotting_shift
use transposes_of_q_and_qdot
use, intrinsic :: iso_c_binding ! Needed for FFTW stuff.
implicit none
include 'fftw3.f03'
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
integer :: iphi, idof, iz, ir, k
! Make this global at some point or use a global work array.
complex(8), dimension(sr:er, sz_phi:ez_phi, ndof, nphi/2+1):: q_phi_space_hat

call transpose_z_to_phi (ndof, q, q_phi_space)
! q_phi_space (sr:er, sz_phi:ez_phi, ndof, nphi)

call dfftw_execute_dft_r2c(forward_plan, q_phi_space, q_phi_space_hat)

do k = 1, nphi/2+1
   do idof = 1, ndof
      do iz = sz_phi, ez_phi
         do ir = sr, er
            q_phi_space_hat(ir,iz,idof,k) = q_phi_space_hat(ir,iz,idof,k) / nphi * complex_shift_factor(k, ir)
         end do
      end do
   end do
end do

call dfftw_execute_dft_c2r(inverse_plan, q_phi_space_hat, q_phi_space)

call transpose_phi_to_z (ndof, q_phi_space, q)

end subroutine apply_real_shifts_using_fftw

!----------------------------------------------------------------------------------85

subroutine apply_real_shifts_using_fftw_old(q)

! Cache unfriendly but I found this to be cheaper.
! Does a real valued shift in phi.

use grid
use partition_data
use fargo_or_plotting_shift
use transposes_of_q_and_qdot
use, intrinsic :: iso_c_binding ! Needed for FFTW stuff.
implicit none
include 'fftw3.f03'
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
integer :: iphi, idof, iz, ir
! Make these global at some point or use a global work array.
real   (8), dimension(nphi    ) :: f
complex(8), dimension(nphi/2+1) :: fhat(nphi/2 + 1)
integer(8) :: forward_plan, inverse_plan ! 8 byte addresses.
integer :: k

! Debug:
real   (8), dimension(nphi    ) :: f_shifted

call dfftw_plan_dft_r2c_1d(forward_plan, nphi, f, fhat, FFTW_ESTIMATE)
call dfftw_plan_dft_c2r_1d(inverse_plan, nphi, fhat, f, FFTW_ESTIMATE)

call transpose_z_to_phi (ndof, q, q_phi_space)

! q_phi_space (sr:er, sz_phi:ez_phi, ndof, nphi)

! Replace this cache-unfriendly stuff later:
do idof = 1, ndof
   do iz = sz_phi, ez_phi
      do ir = sr, er
         do iphi = 1, nphi
            f(iphi) = q_phi_space(ir, iz, idof, iphi)
         end do

         call dfftw_execute_dft_r2c(forward_plan, f, fhat)         

         do k = 1, nphi/2+1
            fhat(k) = fhat(k) / nphi * complex_shift_factor(k, ir)
         end do

         call dfftw_execute_dft_c2r(inverse_plan, fhat, f_shifted)

         do iphi = 1, nphi
            q_phi_space(ir, iz, idof, iphi) = f_shifted(iphi)
         end do
      end do ! ir
   end do
end do
         
call transpose_phi_to_z (ndof, q_phi_space, q)

end subroutine apply_real_shifts_using_fftw_old

#endif

!----------------------------------------------------------------------------------85

#ifdef rogallo_fft

subroutine initialize_rogallo_fft_for_fargo

use rogallo_fft_for_fargo
use grid
implicit none

allocate (trig_table(nphi))
call rtrig(nphi, trig_table)

end subroutine initialize_rogallo_fft_for_fargo

!----------------------------------------------------------------------------------85

subroutine apply_real_shifts_using_rogallo_fft(q)

! Does a real valued shift in phi.

use rogallo_fft_for_fargo
use grid
use partition_data
use fargo_or_plotting_shift
use transposes_of_q_and_qdot
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
integer :: iphi, idof, iz, ir, k, lot, inc, jump
complex(8) :: fhat ! temporary

if (my_node .eq. 0) print *, ' apply_real_shifts_using_rogallo_fft has been called'

call transpose_z_to_phi (ndof, q, q_phi_space)
! q_phi_space (sr:er, sz_phi:ez_phi, ndof, nphi)

lot  = mr*mz_phi*ndof
inc  = lot
jump = 1
call rfft_dri (nphi, lot, q_phi_space, inc, jump, trig_table, -1)

do iphi = 1, nphi, 2
   k = (iphi+1)/2 ! - 1  ! The k below is actually measured from 1
   do idof = 1, ndof
      do iz = sz_phi, ez_phi
         do ir = sr, er
            fhat = CMPLX(q_phi_space(ir,iz,idof,iphi), q_phi_space(ir,iz,idof,iphi+1)) &
                 / nphi * complex_shift_factor(k, ir) / 2.d0
            q_phi_space(ir,iz,idof,iphi)   = DREAL(fhat)
            q_phi_space(ir,iz,idof,iphi+1) = DIMAG(fhat)
         end do
      end do
   end do
end do

call rfft_dri (nphi, lot, q_phi_space, inc, jump, trig_table, 1)
call transpose_phi_to_z (ndof, q_phi_space, q)

end subroutine apply_real_shifts_using_rogallo_fft

#endif

! End of ifdef section for rogallo fft

!----------------------------------------------------------------------------------85

subroutine apply_integer_shifts(q)

! Does a shift in phi using the array inew_phi(ir, iphi) which was obtained in
! subroutine get_integer_shifts.

use grid
use partition_data
use fargo_or_plotting_shift, only: inew_phi
use transposes_of_q_and_qdot
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
integer :: iphi, idof, iz, ir
! Make this global at some point or use a global work array.
real(8), dimension(sr:er, sz_phi:ez_phi, ndof, nphi) :: q_temp

integer :: inew

if (my_node .eq. 0) print *, ' apply_phi_shift has been called'

call transpose_z_to_phi (ndof, q, q_phi_space)
! q_phi_space (sr:er, sz_phi:ez_phi, ndof, nphi)

do iphi = 1, nphi
   do idof = 1, ndof
      do iz = sz_phi, ez_phi
         do ir = sr, er
            inew = inew_phi(ir, iphi)
            q_temp(ir, iz, idof, inew) = q_phi_space(ir, iz, idof, iphi)
         end do
      end do
   end do
end do

do iphi = 1, nphi
   do idof = 1, ndof
      do iz = sz_phi, ez_phi
         do ir = sr, er
            q_phi_space(ir, iz, idof, iphi) = q_temp(ir, iz, idof, iphi)
         end do
      end do
   end do
end do

call transpose_phi_to_z (ndof, q_phi_space, q)

end subroutine apply_integer_shifts

!----------------------------------------------------------------------------------85

subroutine add_fargo_extra_operator_terms(q, qdot, lambda_extra_operator_max)

! Every radial derivative term needs a corresponding phi derivative.  See notes of
! Feb. 10, 2018.
! Both q and qdot should be in phi-space: (sr:er, sz_phi:ez_phi, ndof, nphi) 

use grid
use partition_data
use dof_indices, only: irho, zmom, rmom, amom, ener
use thermal_parameters
use fargo_or_plotting_shift
implicit none

real(8) :: dt, t, t0
real(8), dimension(sr:er, sz_phi:ez_phi, ndof, nphi), intent(in)  :: q
real(8), dimension(sr:er, sz_phi:ez_phi, ndof, nphi), intent(out) :: qdot
! Maximum eigenvalue.  To be used for time step selection.
real(8),                                intent(out) :: lambda_extra_operator_max  

! Local:
integer :: ir, iz, iphi, idof

! Again, see later if these can be made static.
! ndof ---> ndof+1 for -p div_u term in internal energy equation.
real(8), dimension(sr:er, sz_phi:ez_phi, ndof+1, nphi) :: F, dF
real(8) :: ur, c_sound, lambda_local, gamma_sound_speed

! Pressure associated locals:
real(8), dimension(sr:er, sz_phi:ez_phi, nphi) :: p, dp_dphi

integer :: div_u = ndof + 1 ! Index for div u term required for internal energy equation.

! Note: my_node = 0 for the serial code so this works both cases.
!if (my_node .eq. 0) then
!   print *, ' add_fargo_correction_terms has been called'
!   print *, ' enter anything to continue'
!   read (5, *)
!end if
   
#ifdef debug_print
   ! Note: my_node = 0 for the serial code so this works both cases.
   if (my_node .eq. 0) then
      print *, ' add_fargo_extra_operator_terms has been called'
   end if
#endif

if (isothermal) then
   gamma_sound_speed = 1.0d0
else
#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' gamma = ', gamma
      print *, ' gm1 = ', gm1
   end if
#endif

   gamma_sound_speed = gamma
end if

! Store pressure which is needed below:
if (isothermal) then
   do iphi = 1, nphi
      do iz = sz_phi, ez_phi
         do ir = sr, er
            p(ir, iz, iphi) = q(ir, iz, irho, iphi) * ci_squared_initial(ir, iz)
         end do
      end do
   end do
else
   ! print *, ' calculating pressure for non-isothermal option'
   do iphi = 1, nphi
      do iz = sz_phi, ez_phi
         do ir = sr, er
            p(ir, iz, iphi) = q(ir, iz, ener, iphi) * gm1
         end do
      end do
   end do
end if

! These are the same things as in add_radial_derivatives
lambda_extra_operator_max = 0.0d0
do iphi = 1, nphi
   do iz = sz_phi, ez_phi
      do ir = sr, er
         ur = q(ir, iz, rmom, iphi) / q(ir, iz, irho, iphi)

         F(ir, iz, irho, iphi) =  q(ir, iz, rmom, iphi) * rgrid(ir)
         F(ir, iz, amom, iphi) =  q(ir, iz, amom, iphi) * rgrid(ir) * ur
         F(ir, iz, rmom, iphi) =  rgrid(ir)*(p(ir,iz,iphi) + q(ir, iz, rmom, iphi)*ur)
         F(ir, iz, zmom, iphi) =  q(ir, iz, zmom, iphi) * rgrid(ir) * ur

         ! These are not needed for the isothermal case but does no harm:
         F(ir, iz, ener,  iphi) = q(ir, iz, ener, iphi) * rgrid(ir) * ur
         F(ir, iz, div_u, iphi) = rgrid(ir) * ur         

         c_sound = SQRT(gamma_sound_speed * p(ir, iz, iphi) / q(ir, iz, irho, iphi))
         lambda_local = (ABS(ur) + c_sound) * abs(fargo_factor(ir))
         lambda_extra_operator_max = MAX(lambda_extra_operator_max, lambda_local) 
      end do
   end do
end do

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' in add_fargo_correction_terms: lambda_extra_operator_max = ', lambda_extra_operator_max
   end if
#endif

! Differentiate w.r.t. phi:
call pade_diff_periodic(mr*mz_phi*(ndof+1), nphi, dphi, F(sr, sz_phi, 1, 1), dF(sr, sz_phi, 1, 1))

do iphi = 1, nphi
   do idof = 1, ndof
      do iz = sz_phi, ez_phi
         do ir = sr, er
            qdot(ir, iz, idof, iphi) = qdot(ir, iz, idof, iphi) - &
                                         dF(ir, iz, idof, iphi) * fargo_factor_over_r(ir)
         end do
      end do
   end do
end do

! - p div u term in internal energy equation:
do iphi = 1, nphi
   do iz = sz_phi, ez_phi
      do ir = sr, er
         qdot(ir, iz, ener, iphi) = qdot(ir, iz, ener, iphi) - &
                                    p(ir, iz, iphi) * dF(ir, iz, div_u, iphi) * fargo_factor_over_r(ir)
      end do
   end do
end do

! Remember that at the beginning of a time step, the extra operator is zero.
!print *, ' in add_fargo_extra_operator_terms'
!print *, ' lambda_extra_operator_max = ', lambda_extra_operator_max
!print *, ' enter anything to continue'
!read (5, *)

return
end subroutine add_fargo_extra_operator_terms

!----------------------------------------------------------------------------------85

subroutine fargo_informational_output (t)

use grid
use logical_units
use fargo_or_plotting_shift
implicit none
real(8) :: t

! Locals:
character(80) :: file1, file2
integer :: ir

write (file1, 1) int(t), t - int(t)
1 format ('Omega_fargo_', i6.6, f0.4, '.dat')

write (file2, 2) int(t), t - int(t)
2 format ('nshift_', i6.6, f0.4, '.dat')

open (unit = lun_general_purpose,   file = file1, form = 'formatted', status = 'unknown')
open (unit = lun_general_purpose_2, file = file2, form = 'formatted', status = 'unknown')

do ir = 1, nr
   write (lun_general_purpose,   3) rgrid(ir), Omega_bar(ir), Omega_fargo(ir)
3  format (3(1x, e12.5))
   write (lun_general_purpose_2, 4) rgrid(ir), nshift(ir)
4  format (e12.5, 1x, i6)   
end do
close (lun_general_purpose)
close (lun_general_purpose_2)
print *, ' Wrote Omega_fargo.dat and nshift.dat'

return
end subroutine fargo_informational_output

!----------------------------------------------------------------------------------85

subroutine plotting_shift(q, t)

! To use this subroutine you need to call setup_plotting_shift
! This subroutine must be called by the application subroutine.

use grid
use partition_data
use fargo_or_plotting_shift
use math_constants
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: t

! Local:
integer :: ir, k
complex(8) :: ci = CMPLX(0.0, 1.0d0)
real(8) :: desired_phi_shift, phase_factor

! This is a backward shift:
desired_phi_shift = - Omega0_plotting_shift * (t - t_of_previous_shift)
phase_factor      = desired_phi_shift * twopi_over_Delta_phi_domain

do ir = sr, er
   do k = 1, nphi/2+1
      complex_shift_factor(k, ir) = EXP(-ci * (k-1) * phase_factor)
   end do
end do

#ifdef fftw
   ! call apply_real_shifts_using_fftw(q)
   call apply_real_shifts_using_fftw_old(q) ! this can actually be cheaper than the above
#else
   call apply_real_shifts_using_rogallo_fft(q)
#endif

t_of_previous_shift = t

end subroutine plotting_shift

!----------------------------------------------------------------------------------85

