!----------------------------------------------------------------------------------85

subroutine rk4(rhs, q, cfl, t, dt, use_supplied_dt, unphysical, iphi_bad, &
               hit_target, t_target, target_met)

! Performs one step of RK4 as given in Abramowitz and Stegun,
! Handbook of Mathematical Functions, p. 896, 25.5.10.
! The FARGO trick is used to reduce the time step for advection in the phi-direction.

! If use_supplied = true, we will use the given dt instead of calculating one.
! Then cfl will be the resulting cfl of the step.

use grid
use partition_data
! Work arrays: accumulator, next argument, rhs
use rk4_arrays, only: q_accum, q_next_arg, qdot
use fargo_or_plotting_shift
use pade_filter

use cpu_timing_module
use thermal_parameters
use math_constants
use dof_indices
use sponge
use boundary_condition_data
! So I can call enforce_BC
use boundary_condition_routines
implicit none

! The array being advanced:
real(8), intent(inout), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8), intent(inout) :: cfl
real(8), intent(inout) :: t  ! will be updated
! If use_supplied_dt = false, this will be determined internally and output for
! diagnostic purposes.
real(8), intent(inout) :: dt 
logical, intent(in)    :: use_supplied_dt
logical, intent(out)   :: unphysical
integer, intent(out)   :: iphi_bad
logical, intent(in)    :: hit_target
real(8), intent(in)    :: t_target
logical, intent(out)   :: target_met

! Locals:
integer :: substep_number
real(8) :: lambda_max ! Will be used to determine the time step.
real(8) :: lambda_max_extra_operator ! To print for informational purposes.
real(8) :: eps_filter_now, cfl_given_dt

! For timing:
integer :: start_count_for_step, end_count_for_step, clock_count_rate
real(8) :: cpu_for_this_step, transpose_cpu_this_step

! For clarity, time at the beginning of the step.  For use by Fargo stuff:
real(8) :: t0

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: first executable of rk4'
#endif

! Just to avoid a warning that iphi_bad, an intent_out variables is not
! given a value.  It is, of course, given a value, by subroutine rhs.
iphi_bad = 0

#ifdef transpose_timing
   transpose_count = 0
#endif

call system_clock (count = start_count_for_step)

! Step 1:
!~~~~~~~~
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: about to call rhs in rk step 1'
#endif
t0 = t      ! This is for Fargo.

! Note: In the first sub-step, rhs returns the "eigenvalue" lambda_max
! (and lambda_max_extra_operator for Fargo) which we use to choose the
! time step:
substep_number = 1
! We return iphi_bad so that we can plot the plane with the bad iphi.
call rhs(q, qdot, substep_number, t, t0, lambda_max, lambda_max_extra_operator, unphysical, iphi_bad)
if (unphysical) return
   
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: returned from rhs'
#endif

! Calculate dt or cfl given dt:
if (use_supplied_dt) then
   ! The pi factor is for spectral schemes and is a good upper bound.               
   cfl = pi*lambda_max*dt
   if (my_node .eq. 0) then
      print *, ' Given lambda_max = ', lambda_max, ' & dt = ', dt, ' implies cfl = ', cfl
   end if
else
   dt = cfl/pi/lambda_max
   if (my_node .eq. 0) then
      print *, ' Given lambda_max = ', lambda_max, ' & cfl = ', cfl, ' implies dt = ', dt
   end if
end if

if (hit_target) then
   if (t + dt .gt. t_target) then
      dt = t_target - t
      target_met = .true.
   else
      target_met = .false.
   end if
end if

if (apply_sponge) call add_sponge_term(dt, q, qdot)

q_next_arg = q + 0.5d0    *dt*qdot ! arg for k2
q_accum    = q + 1.d0/6.d0*dt*qdot

call enforce_BC(q_next_arg)

! Now that we have dt we can get the Fargo shifts and uphi_fargo_subtract which will
! be needed in the phi momentum equation.
! Note: apply_fargo_this_step is determined in subroutine rhs at the first substep.
! It will be set to .true. of apply_fargo is .true. AND if it is advantageous to do
! so.  Specifically, if lambda max from the non-Euler terms is > than
! lambda max of Euler terms without fargo, then it is point-less to do Fargo.
if (apply_fargo_this_step) then
   if (my_node .eq. 0) print *, ' Applying FARGO this step'
   ! This subroutine is located in fargo_and_plotting_shift.f90
   call get_fargo_shifts_and_uphi_fargo_subtract(t, dt)
end if

! Step 2:
!~~~~~~~~
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: about to call rhs in rk step 2'
#endif

! gives k2/dt
substep_number = 2
call rhs(q_next_arg, qdot, substep_number, t + 0.5d0*dt, t0, lambda_max, lambda_max_extra_operator, &
         unphysical, iphi_bad)
if (unphysical) return   
if (apply_sponge) call add_sponge_term(dt, q_next_arg, qdot)
q_next_arg = q       +     0.5d0*dt*qdot ! arg for k3
q_accum    = q_accum + 1.d0/3.d0*dt*qdot

call enforce_BC(q_next_arg)

! Step 3:
!~~~~~~~~
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: about to call rhs in rk step 3'
#endif

! gives k3/dt
substep_number = 3   
call rhs(q_next_arg, qdot, substep_number, t + 0.5d0*dt, t0, lambda_max, lambda_max_extra_operator, &
         unphysical, iphi_bad)
if (unphysical) return   
if (apply_sponge) call add_sponge_term(dt, q_next_arg, qdot)   
q_next_arg = q       +           dt*qdot ! arg for k4
q_accum    = q_accum + 1.d0/3.d0*dt*qdot

call enforce_BC(q_next_arg)

! Step 4:
!~~~~~~~~
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: about to call rhs in rk step 4'
#endif

! gives k4/h
substep_number = 4   
call rhs(q_next_arg, qdot, substep_number, t + dt, t0, lambda_max, lambda_max_extra_operator, &
         unphysical, iphi_bad)
if (unphysical) return      
if (apply_sponge) call add_sponge_term(dt, q_next_arg, qdot)
q = q_accum + 1.d0/6.d0*dt*qdot

! Shift flowfield back to original coordinate system:
if (apply_fargo_this_step) then
   call apply_fargo_shifts(q)
   if (my_node .eq. 0) print *, ' In last sub-step, lambda_max_extra_operator = ', &
        lambda_max_extra_operator
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: about to call enforce_BC in rk4'
#endif

call enforce_BC(q)

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: returned from enforce_BC in rk4'
#endif

t = t + dt

if (apply_pade_filter) then
   if (eps_or_tau .eq. 'tau') then
      eps_filter_now = 16.d0 * dt / tau_filter
   else
      eps_filter_now = eps_filter
   end if
   if (my_node .eq. 0) then
      print *, ' node = 0: calling smooth_field with eps_filter = ', eps_filter_now
   end if
   call smooth_field (ndof, q, eps_filter_now)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node = 0: rk4: returned from smooth_field'
#endif   
end if
   
! We don't need to add the Fargo extra operator outside of time-stepping:
add_fargo_extra_operator_now = .false.

!  I am doing this because, currently stuff is being done (without harm) to compute
! terms in th energy equation even for the isothermal case and this is causing the
! q(:, :, :, ener) to grow:
if (isothermal) then
   q(:, :, :, ener) = 0.d0
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node = 0: rk4: returning from rk4'
#endif

if (my_node .eq. 0) then
   call system_clock (count = end_count_for_step)
   call system_clock(count_rate = clock_count_rate)
   cpu_for_this_step = (end_count_for_step - start_count_for_step) / real(clock_count_rate, 8)
   total_cpu_for_stepping = total_cpu_for_stepping + cpu_for_this_step
   print *, ' node = 0: cpu time for this step = ', cpu_for_this_step
   print *, ' node = 0: total cpu secs for stepping = ', total_cpu_for_stepping
end if

#ifdef transpose_timing
   if (my_node .eq. 0) then
      transpose_cpu_this_step = transpose_count / real(clock_count_rate, 8)
      print *, ' cpu fraction for transposes this step = ', &
           transpose_cpu_this_step/cpu_for_this_step
      total_cpu_for_transposes = total_cpu_for_transposes + transpose_cpu_this_step
      print *, ' cpu fraction for transposes till this step = ', &
           total_cpu_for_transposes/total_cpu_for_stepping
    end if
#endif   

end subroutine rk4

!----------------------------------------------------------------------------------85

subroutine euler (q, cfl, t, dt, use_supplied_dt)

! Note that the Fargo extra operator never needs to be applied for the Euler shcheme.

use grid
use partition_data
use rk4_arrays
use fargo_or_plotting_shift
use cpu_timing_module
use fargo_or_plotting_shift
use math_constants
use boundary_condition_data
! So I can call enforce_BC
use boundary_condition_routines
implicit none

! Arguments:
real(8), intent(inout), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: cfl, t, dt
logical :: use_supplied_dt

! Locals:
logical :: first_substep
real(8) :: t0, lambda_max, lambda_max_extra_operator

! For timing:
integer :: start_count_for_step, end_count_for_step, clock_count_rate
real(8) :: cpu_for_this_step

call system_clock (count = start_count_for_step)

t0 = t
first_substep = .true.
call rhs(q, qdot, first_substep, t, t0, lambda_max, lambda_max_extra_operator)

! Calculate dt or cfl given dt:
if (use_supplied_dt) then
   ! The pi factor is for spectral schemes and is a good upper bound.         
   cfl = pi*lambda_max*dt  
   print *, ' Given dt = ', dt, ' implies cfl = ', cfl      
else
   dt = cfl/pi/lambda_max
   print *, ' Given cfl = ', cfl, ' implies dt = ', dt
end if

! Note: We will never use uphi_subtract for the Euler scheme, since the
! the extra operator is non-zero only after the first sub-step.
call get_fargo_shifts_and_uphi_fargo_subtract(t, dt)

q = q + dt*qdot
if (apply_fargo_this_step) then
   call apply_fargo_shifts(q)
end if
call enforce_BC(q)

t = t + dt

call system_clock (count = end_count_for_step)
call system_clock(count_rate = clock_count_rate)
cpu_for_this_step = (end_count_for_step - start_count_for_step) / real(clock_count_rate, 8)
total_cpu_for_stepping = total_cpu_for_stepping + cpu_for_this_step
print *, ' cpu time for this step = ', cpu_for_this_step

end subroutine euler

!----------------------------------------------------------------------------------85

subroutine get_lambda_euler_fargo(q, lambda_max_euler, lambda_max_euler_fargo)

! Here "euler" refers to the inviscid Euler equations.

! lambda is the eigenvalue that goes into time step selection: pi*lambda*dt = cfl.

! lambda for the Euler equation part of the problem we are solving.

! lambda is defined to be lambda = max over the domain of:
! (local max of characteristic speed / grid size)

! Then we have CFL number = pi * lambda * dt

! This subroutine is needed for fargo since we need the time step before
! we begin the step.  Specifically it is needed to get Omega_fargo in
! get_integer_shifts which is then used when we subtract out Omega_fargo.

use grid
use dof_indices, only: irho, zmom, rmom, amom, ener
use partition_data
use fargo_or_plotting_shift
use thermal_parameters
use logical_units
use math_constants
#ifdef mpi_code
   use mpi
#endif
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in) :: q
real(8), intent(out) :: lambda_max_euler, lambda_max_euler_fargo

! Locals:
integer :: iz, iphi, ir, ier
real(8) :: c, p, rho, ur, uz, uphi
real(8) :: lambda_r, lambda_z, lambda_phi, lambda_phi_fargo
real(8) :: lambda_r_max, lambda_z_max, lambda_phi_max, lambda_phi_fargo_max
real(8) :: lambda_r_max_global, lambda_z_max_global, lambda_phi_max_global, &
     lambda_phi_fargo_max_global
real(8) :: uphi_fargo

lambda_r_max         = 0.0d0
lambda_z_max         = 0.0d0
lambda_phi_max       = 0.0d0
lambda_phi_fargo_max = 0.0d0

if (.not. isothermal) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            rho = q(ir, iphi, iz, irho)
            p = q(ir, iphi, iz, ener) * gm1
            c = SQRT(gamma * p / rho)

            ur   = q(ir, iphi, iz, rmom) / rho 
            uz   = q(ir, iphi, iz, zmom) / rho
            uphi = q(ir, iphi, iz, amom) / (rho * rgrid(ir))
            ! Note: Omega_bar(ir) is actually the target.  For integer shifting it will
            ! be different.
            uphi_fargo = uphi - Omega_bar(ir)*rgrid(ir)

            lambda_r   = (ABS(ur)   + c) * Ji_r(ir)
            lambda_z   = (ABS(uz)   + c) * Ji_z(iz)
            ! Without Fargo:
            lambda_phi       = (ABS(uphi)       + c) / (rgrid(ir) * dphi)
            ! With Fargo:
            lambda_phi_fargo = (ABS(uphi_fargo) + c) / (rgrid(ir) * dphi)

            lambda_r_max         = MAX(lambda_r,   lambda_r_max)
            lambda_z_max         = MAX(lambda_z,   lambda_z_max)
            lambda_phi_max       = MAX(lambda_phi, lambda_phi_max)
            lambda_phi_fargo_max = MAX(lambda_phi_fargo, lambda_phi_fargo_max)                      
         end do
      end do
   end do
else
   ! Isothermal:
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            rho = q(ir, iphi, iz, irho)
            c = SQRT(ci_squared_initial(ir, iz))

            ur   = q(ir, iphi, iz, rmom) / rho 
            uz   = q(ir, iphi, iz, zmom) / rho
            uphi = q(ir, iphi, iz, amom) / (rho * rgrid(ir))
            uphi_fargo = uphi - Omega_bar(ir)*rgrid(ir)

            lambda_r   = (ABS(ur)   + c) * Ji_r(ir)
            lambda_z   = (ABS(uz)   + c) * Ji_z(iz)
            lambda_phi       = (ABS(uphi)       + c) / (rgrid(ir) * dphi)            
            lambda_phi_fargo = (ABS(uphi_fargo) + c) / (rgrid(ir) * dphi)            

            lambda_r_max         = MAX(lambda_r,   lambda_r_max)
            lambda_z_max         = MAX(lambda_z,   lambda_z_max)
            lambda_phi_max       = MAX(lambda_phi, lambda_phi_max)
            lambda_phi_fargo_max = MAX(lambda_phi_fargo, lambda_phi_fargo_max)                      
         end do
      end do
   end do   
end if

#ifdef mpi_code
   ! Take max over all the nodes.
   call mpi_reduce(lambda_r_max, lambda_r_max_global, 1, mpi_double, mpi_max, &
        0, mpi_comm_world, ier)
   call mpi_reduce(lambda_z_max, lambda_z_max_global, 1, mpi_double, mpi_max, &
        0, mpi_comm_world, ier)
   call mpi_reduce(lambda_phi_max, lambda_phi_max_global, 1, mpi_double, mpi_max, &
        0, mpi_comm_world, ier)
   call mpi_reduce(lambda_phi_fargo_max, lambda_phi_fargo_max_global, 1, mpi_double, mpi_max, &
        0, mpi_comm_world, ier)
   ! I do this so I don't need two sets of prints for the serial and parallel codes.
   if (my_node .eq. 0) then
      lambda_r_max         = lambda_r_max_global
      lambda_z_max         = lambda_z_max_global
      lambda_phi_max       = lambda_phi_max_global
      lambda_phi_fargo_max = lambda_phi_fargo_max_global
   end if
#endif

! Now the serial and parallel code have the same notation:
   
if (my_node .eq. 0) then
   ! Some of the above lambdas could be irrelevant:
   print *, ' In get_lambda_euler_part:'      
   if (nz   .le. 2) then
      lambda_z_max   = 0.0d0
   else
      print *, '    lambda_z_max           = ', lambda_z_max
   end if
   
   if (nphi .eq. 1) then
      lambda_phi_max       = 0.0d0
      lambda_phi_fargo_max = 0.0d0
   else
      print *, '    lambda_phi_max         = ', lambda_phi_max, &
                    ' lambda_phi_fargo_max = ', lambda_phi_fargo_max          
   end if
   print *, '    lambda_r_max           = ', lambda_r_max
   
   lambda_max_euler       = MAX(lambda_r_max,lambda_z_max,lambda_phi_max      )
   lambda_max_euler_fargo = MAX(lambda_r_max,lambda_z_max,lambda_phi_fargo_max)
   if (nphi .eq. 1) then
      print *, '    lambda_max_euler       = ', lambda_max_euler
   else
      print *, ' apply_fargo = ', apply_fargo
      print *, '    lambda_max_euler       = ', lambda_max_euler, &      
                  ' lambda_max_euler_fargo = ', lambda_max_euler_fargo
      print *, '    reduction ratio        = ', lambda_max_euler_fargo / lambda_max_euler
   end if
end if
   
#ifdef mpi_code
   call mpi_barrier(mpi_comm_world, ier) ! probably not needed.
   call mpi_bcast(lambda_max_euler,       1, mpi_double, 0, mpi_comm_world, ier)
   call mpi_bcast(lambda_max_euler_fargo, 1, mpi_double, 0, mpi_comm_world, ier)   
#endif

end subroutine get_lambda_euler_fargo

!----------------------------------------------------------------------------------85

subroutine get_lambda_euler_non_fargo(q, lambda_max_euler, dominated_by_phi)

! lambda is the eigenvalue that goes into time step selection: pi*lambda*dt = cfl.

! lambda for the Euler equation part of the problem we are solving.

! lambda is defined to be lambda = max over the domain of:
! (local max of characteristic speed / grid size)

! Then we have CFL number = pi * lambda * dt

! This subroutine is needed for fargo since we need the time step before
! we begin the step.  Specifically it is needed to get Omega_fargo in
! get_integer_shifts which is then used when we subtract out Omega_fargo.

use grid
use dof_indices, only: irho, zmom, rmom, amom, ener
use partition_data
use fargo_or_plotting_shift
use thermal_parameters
use logical_units
use math_constants
#ifdef mpi_code
   use mpi
#endif
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in) :: q
real(8), intent(out) :: lambda_max_euler
logical, intent(out) :: dominated_by_phi

! Locals:
integer :: iz, iphi, ir, ier
real(8) :: c, p, rho, ur, uz, uphi
real(8) :: lambda_r, lambda_z, lambda_phi
real(8) :: lambda_r_max, lambda_z_max, lambda_phi_max
real(8) :: lambda_r_max_global, lambda_z_max_global, lambda_phi_max_global

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: First executable of get_lambda_euler_non_fargo'
#endif

lambda_r_max         = 0.0d0
lambda_z_max         = 0.0d0
lambda_phi_max       = 0.0d0

if (.not. isothermal) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            rho = q(ir, iphi, iz, irho)
            p = q(ir, iphi, iz, ener) * gm1
            c = SQRT(gamma * p / rho)

            ur   = q(ir, iphi, iz, rmom) / rho 
            uz   = q(ir, iphi, iz, zmom) / rho
            uphi = q(ir, iphi, iz, amom) / (rho * rgrid(ir))

            lambda_r   = (ABS(ur)   + c) * Ji_r(ir)
            lambda_z   = (ABS(uz)   + c) * Ji_z(iz)
            ! Without Fargo:
            lambda_phi       = (ABS(uphi)       + c) / (rgrid(ir) * dphi)

            lambda_r_max         = MAX(lambda_r,   lambda_r_max)
            lambda_z_max         = MAX(lambda_z,   lambda_z_max)
            lambda_phi_max       = MAX(lambda_phi, lambda_phi_max)
         end do
      end do
   end do
else
   ! Isothermal:
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In get_lambda_euler_non_fargo, starting loop for isothermal'
#endif   
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            rho = q(ir, iphi, iz, irho)
            c = SQRT(ci_squared_initial(ir, iz))

            ur   = q(ir, iphi, iz, rmom) / rho 
            uz   = q(ir, iphi, iz, zmom) / rho
            uphi = q(ir, iphi, iz, amom) / (rho * rgrid(ir))

            lambda_r   = (ABS(ur)   + c) * Ji_r(ir)
            lambda_z   = (ABS(uz)   + c) * Ji_z(iz)
            lambda_phi = (ABS(uphi)       + c) / (rgrid(ir) * dphi)            

            lambda_r_max         = MAX(lambda_r,   lambda_r_max)
            lambda_z_max         = MAX(lambda_z,   lambda_z_max)
            lambda_phi_max       = MAX(lambda_phi, lambda_phi_max)
         end do
      end do
   end do
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In get_lambda_euler_non_fargo, finished loop for isothermal'
#endif   
end if

#ifdef mpi_code
   ! Take max over all the nodes.
   call mpi_reduce(lambda_r_max, lambda_r_max_global, 1, mpi_double, mpi_max, &
        0, mpi_comm_world, ier)
   call mpi_reduce(lambda_z_max, lambda_z_max_global, 1, mpi_double, mpi_max, &
        0, mpi_comm_world, ier)
   call mpi_reduce(lambda_phi_max, lambda_phi_max_global, 1, mpi_double, mpi_max, &
        0, mpi_comm_world, ier)
   ! I do this so I don't need two sets of prints for the serial and parallel codes.
   !if (my_node .eq. 0) then
      lambda_r_max         = lambda_r_max_global
      lambda_z_max         = lambda_z_max_global
      lambda_phi_max       = lambda_phi_max_global
   !end if
#endif

! Now the serial and parallel code have the same notation:
   
if (my_node .eq. 0) then
   ! Some of the above lambdas could be irrelevant:
   print *, ' In get_lambda_euler_part:'      
   if (nz   .le. 2) then
      lambda_z_max   = 0.0d0
   else
      print *, '    lambda_z_max           = ', lambda_z_max
   end if
   
   if (nphi .eq. 1) then
      lambda_phi_max       = 0.0d0
   else
      print *, '    lambda_phi_max         = ', lambda_phi_max
   end if
   print *, '    lambda_r_max           = ', lambda_r_max
   
   lambda_max_euler       = MAX(lambda_r_max,lambda_z_max,lambda_phi_max      )
   print *, '    lambda_max_euler       = ', lambda_max_euler
end if
   
#ifdef mpi_code
   call mpi_barrier(mpi_comm_world, ier) ! probably not needed.
   call mpi_bcast(lambda_max_euler, 1, mpi_double, 0, mpi_comm_world, ier)
#endif

! Boolean which tells the caller that Euler eigenvalue is dominated by the
! phi direction.  This is used for reccomending FARGO to the user if
! appropriate.
if (my_node .eq. 0) dominated_by_phi = lambda_max_euler .eq. lambda_phi_max 

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: Returning from get_lambda_euler_non_fargo'
#endif

end subroutine get_lambda_euler_non_fargo

!----------------------------------------------------------------------------------85

real(8) function Melander3(s, s1, s2)

! Has three parts:
! 1, for s <= s1
! transition for s1 < s < s2
! 0 for s >= s1

use partition_data
implicit none
real(8) :: s, s1, s2

! Local:
real(8) :: r, paren, bracket
real(8), parameter :: kappa = 0.5d0 * exp(2.d0) * log(2.0d0)

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: First executable in Melander3'
   print *, ' s = ', s, ' s1 = ', s1, ' s2 = ', s2
end if
#endif

if (s .le. s1) then
   Melander3 = 1.0d0
else if (s .lt. s2) then
   r        = (s - s1) / (s2 - s1)
   paren    = 1.d0 / (r - 1)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' about to compute bracket and Melander3'
#endif   
   bracket  = - kappa/r * exp(paren)
   Melander3 = 1.d0 - exp(bracket)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' Finished computing bracket and Melander3'
#endif   
else
   Melander3 = 0.0d0
end if

if (Melander3 .lt. 0.d0) then
   print *, ' In Melander3, my_node = ', my_node
   print *, ' s = ', s, ' s1 = ', s1, ' s2 = ', s2, ' Melander3 = ', Melander3
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: Returning from Melander3'
#endif

end function Melander3

!----------------------------------------------------------------------------------85

real(8) function Melander(s, s_cut_off, s_max)

implicit none
real(8) :: s, s_cut_off, s_max

! Local:
real(8) :: r, paren, bracket
real(8), parameter :: kappa = 0.5d0 * exp(2.d0) * log(2.0d0)

r        = (s - s_cut_off) / (s_max - s_cut_off)
paren    = 1.d0 / (r - 1)
bracket  = - kappa/r * exp(paren)
Melander = 1.d0 - exp(bracket)

end function Melander

!**********************************************************************************85
