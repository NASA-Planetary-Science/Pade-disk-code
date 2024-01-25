!----------------------------------------------------------------------------------85

module Cassen_Moosman_infall
   real(8) :: Omega0, t0, Tcloud, rd, c0, Mstar, M0_dot
   real(8), parameter :: m0_Shu = 0.975d0

   ! For deleting infall for small r:
   logical :: apply_rise_in_r
   real(8) :: r_begin_rise_rd, r_end_rise_rd

   logical :: infall_forcing
   real(8) :: vel_forcing_amplitude
   
   ! Boundary values:
   real(8), allocatable, dimension(:,:) :: q_CM_top, q_CM_bot, q_CM_right
   ! Post cooling layer temperature:
   real(8), allocatable, dimension(:) :: T_psl
end module Cassen_Moosman_infall

!----------------------------------------------------------------------------------85

subroutine set_up_Cassen_Moosman_infall (Omega0_arg, t0_arg, Tcloud_arg, GM_out, &
     rd_out, infall_forcing_arg, vel_forcing_Mach_arg)

! This routine must be called before set_up_gravity

! Omega0_arg : Cloud rotation rate [sec^-1]
! t0_arg     : Evolutionary time [yr]
! T0_arg     : Cloud temperature [K]

use partition_data
use physical_constants
use Cassen_Moosman_infall
use logical_units
implicit none
real(8), intent(in)  :: Omega0_arg, t0_arg, Tcloud_arg
real(8), intent(out) :: GM_out, rd_out
logical, intent(in) :: infall_forcing_arg
real(8), intent(in) :: vel_forcing_Mach_arg

Omega0 = Omega0_arg
t0     = t0_arg
Tcloud = Tcloud_arg

c0      = SQRT(Rgas * Tcloud)
M0_dot  = m0_Shu * c0**3 / Gconst
Mstar   = M0_dot * t0
rd      = Omega0**2 * c0 * m0_Shu**3 * t0**3 / 16.d0

GM_out = Gconst*Mstar
rd_out = rd

infall_forcing = infall_forcing_arg

vel_forcing_amplitude = vel_forcing_Mach_arg * c0

if (my_node .eq. 0) then
   print *, ' in set_up_Cassen_Moosman_infall'
   print *, ' rd/AU = ', rd/AU
   print *, ' Mstar/Msolar = ', Mstar/Msolar
   print *, ' c0 = ', c0/100.d0, 'm s^-1'
   print *, ' infall_forcing = ', infall_forcing
   print *, ' vel_forcing_amplitude = ', vel_forcing_amplitude, ' cm/s'
end if

! Defaults for rise.  These can be over-written by calling set_up_Cassen_Moosman_rise
apply_rise_in_r = .false.
r_begin_rise_rd = 0.0d0
r_end_rise_rd   = 0.0d0

end subroutine set_up_Cassen_Moosman_infall

!---------------------------------------------------------------------------------85

subroutine set_up_Cassen_Moosman_rise(apply_rise_in_r_arg, r_begin_rise_rd_arg, r_end_rise_rd_arg)

use Cassen_Moosman_infall
use physical_constants
implicit none
logical :: apply_rise_in_r_arg
real(8) :: r_begin_rise_rd_arg, r_end_rise_rd_arg

apply_rise_in_r = apply_rise_in_r_arg
r_begin_rise_rd = r_begin_rise_rd_arg
r_end_rise_rd   = r_end_rise_rd_arg

end subroutine set_up_Cassen_Moosman_rise

!---------------------------------------------------------------------------------85

subroutine store_Cassen_Moosman_boundary_values

use grid
use total_allocated_words, only: n_words_allocated
use partition_data
use Cassen_Moosman_infall
use dof_indices
use physical_constants
implicit none

! Local:
integer :: ir, iz
real(8) :: uz, ur, uphi, rho, eint 

allocate(q_CM_top(sr:er, ndof), q_CM_bot(sr:er, ndof), q_CM_right(nz, ndof))
n_words_allocated = n_words_allocated + 2*mr*ndof + nz*ndof

do ir = sr, er
   call CM_flow(zgrid(1), rgrid(ir), uz, ur, uphi, rho, eint)
   q_CM_bot(ir, irho) = rho
   q_CM_bot(ir, zmom) = rho * uz
   q_CM_bot(ir, rmom) = rho * ur
   q_CM_bot(ir, amom) = rho * uphi * rgrid(ir)
   q_CM_bot(ir, ener) = eint

   call CM_flow(zgrid(nz), rgrid(ir), uz, ur, uphi, rho, eint)   
   q_CM_top(ir, irho) = rho
   q_CM_top(ir, zmom) = rho * uz
   q_CM_top(ir, rmom) = rho * ur
   q_CM_top(ir, amom) = rho * uphi * rgrid(ir)
   q_CM_top(ir, ener) = eint
end do

ir = nr
do iz = 1, nz
   call CM_flow(zgrid(iz), rgrid(ir), uz, ur, uphi, rho, eint)   
   q_CM_right(iz, irho) = rho
   q_CM_right(iz, zmom) = rho * uz
   q_CM_right(iz, rmom) = rho * ur
   q_CM_right(iz, amom) = rho * uphi * rgrid(ir)
   q_CM_right(iz, ener) = eint
end do

end subroutine store_Cassen_Moosman_boundary_values

!----------------------------------------------------------------------------------85

subroutine apply_Cassen_Moosman_boundary_values (q)

! For the rand function:
#ifdef ifort
   use IFPORT 
#endif

use grid
use total_allocated_words, only: n_words_allocated   
use partition_data
use Cassen_Moosman_infall
use dof_indices
implicit none
real(8), dimension (sr:er, sphi:ephi, nz, ndof) :: q

! Local:
integer :: iphi, ir, iz

if (.not. infall_forcing) then
   do iphi = sphi, ephi
      do ir = sr, er
         q(ir, iphi, 1, irho) = q_CM_bot(ir, irho)
         q(ir, iphi, 1, zmom) = q_CM_bot(ir, zmom)
         q(ir, iphi, 1, rmom) = q_CM_bot(ir, rmom)
         q(ir, iphi, 1, amom) = q_CM_bot(ir, amom)
         q(ir, iphi, 1, ener) = q_CM_bot(ir, ener)

         q(ir, iphi, nz, irho) = q_CM_top(ir, irho)
         q(ir, iphi, nz, zmom) = q_CM_top(ir, zmom)
         q(ir, iphi, nz, rmom) = q_CM_top(ir, rmom)
         q(ir, iphi, nz, amom) = q_CM_top(ir, amom)
         q(ir, iphi, nz, ener) = q_CM_top(ir, ener)
      end do
   end do

   if (er .eq. nr) then
      do iz = 1, nz
         do iphi = sphi, ephi
            q(nr, iphi, iz, irho) = q_CM_right(iz, irho)
            q(nr, iphi, iz, zmom) = q_CM_right(iz, zmom)
            q(nr, iphi, iz, rmom) = q_CM_right(iz, rmom)
            q(nr, iphi, iz, amom) = q_CM_right(iz, amom)
            q(nr, iphi, iz, ener) = q_CM_right(iz, ener)
         end do
      end do
   end if
else ! infall_forcing
   do iphi = sphi, ephi
      do ir = sr, er
         q(ir, iphi, 1, irho) = q_CM_bot(ir, irho)
         q(ir, iphi, 1, zmom) = q_CM_bot(ir, zmom) + q(ir,iphi,1,irho)*(1.d0-2.d0*rand(0))*vel_forcing_amplitude
         q(ir, iphi, 1, rmom) = q_CM_bot(ir, rmom) + q(ir,iphi,1,irho)*(1.d0-2.d0*rand(0))*vel_forcing_amplitude
         q(ir, iphi, 1, amom) = q_CM_bot(ir, amom) + q(ir,iphi,1,irho)*(1.d0-2.d0*rand(0))*vel_forcing_amplitude*rgrid(ir)
         q(ir, iphi, 1, ener) = q_CM_bot(ir, ener)

         q(ir, iphi, nz, irho) = q_CM_top(ir, irho)
         q(ir, iphi, nz, zmom) = q_CM_top(ir, zmom) + q(ir,iphi,nz,irho)*(1.d0-2.d0*rand(0))*vel_forcing_amplitude
         q(ir, iphi, nz, rmom) = q_CM_top(ir, rmom) + q(ir,iphi,nz,irho)*(1.d0-2.d0*rand(0))*vel_forcing_amplitude
         q(ir, iphi, nz, amom) = q_CM_top(ir, amom) + q(ir,iphi,nz,irho)*(1.d0-2.d0*rand(0))*vel_forcing_amplitude*rgrid(ir)
         q(ir, iphi, nz, ener) = q_CM_top(ir, ener)
      end do
   end do

   if (er .eq. nr) then
      do iz = 1, nz
         do iphi = sphi, ephi
            q(nr, iphi, iz, irho) = q_CM_right(iz, irho)
            q(nr, iphi, iz, zmom) = q_CM_right(iz, zmom) + q(nr,iphi,iz,irho)*(1.d0-2.d0*rand(0))*vel_forcing_amplitude
            q(nr, iphi, iz, rmom) = q_CM_right(iz, rmom) + q(nr,iphi,iz,irho)*(1.d0-2.d0*rand(0))*vel_forcing_amplitude
            q(nr, iphi, iz, amom) = q_CM_right(iz, amom) + q(nr,iphi,iz,irho)*(1.d0-2.d0*rand(0))*vel_forcing_amplitude*rgrid(nr)
            q(nr, iphi, iz, ener) = q_CM_right(iz, ener)
         end do
      end do
   end if
end if
   
end subroutine apply_Cassen_Moosman_boundary_values

!----------------------------------------------------------------------------------85

subroutine Cassen_Moosman_initial_condition

use grid
use dof_indices
use partition_data
use q_array
use Cassen_Moosman_infall
implicit none

! Local:
integer :: ir, iz, iphi, idof
real(8) :: uz, ur, uphi, rho, eint

if (my_node .eq. 0) print *, ' Cassen_Moosman_initial_condition has been called'

iphi = sphi
do iz = 1, nz
   do ir = sr, er
      !print *, ' ir = ', ir, ' iz = ', iz, ' iphi = ', iphi
      !print *, ' calling CM_flow'
      call CM_flow(zgrid(iz), rgrid(ir), uz, ur, uphi, rho, eint)
      if (rho .gt. 1.0d-15) rho = 1.0d-15
      q(ir, iphi, iz, irho) = rho
      q(ir, iphi, iz, zmom) = rho * uz
      q(ir, iphi, iz, rmom) = rho * ur
      q(ir, iphi, iz, amom) = rho * uphi * rgrid(ir)
      q(ir, iphi, iz, ener) = eint
   end do
end do

! Copy to the rest of the meridional planes:
do idof = 1, 5
   do iz = 1, nz
      do iphi = sphi+1, ephi
         do ir = sr, er
            q(ir, iphi, iz, idof) = q(ir, sphi, iz, idof)
         end do
      end do
   end do
end do

end subroutine Cassen_Moosman_initial_condition

!----------------------------------------------------------------------------------85

subroutine CM_flow (z, r, uz, ur, uphi, rho, eint)

! Computes the Ulrich-Cassen-Moosman flow at the point (z, r) in cylindrical
! coordinates.  Uses formulas in Chevalier form from Terebey, Shu and Cassen (1984),
! pp. 544--545.

! theta0 = polar angle of the asymptotically radial streamline.

use Cassen_Moosman_infall
use physical_constants
use math_constants, only : pi
use thermal_parameters, only: gm1
use gravity, only: GM
! For rise function:
use grid
implicit none
real(8), intent(in)  :: z, r
real(8), intent(out) :: uz, ur, uphi, rho, eint

! Local:
real(8) :: Rspherical, cost, sint, Bracket4, zeta, fac1, cos_ratio, fac2, fac2m, &
     uRspherical, utheta, sint0, P2, cost0, pressure
real(8) :: gr, gz, Phi_g, KE, PE_abs, relative_error
logical :: in_cavity

! Called function:
real(8) :: rise_function

!print *, ' CM_flow has been called for (r/AU, z/AU) = ', r/AU, z/AU

! Spherical coordinates:
Rspherical = SQRT(z**2 + r**2)
cost = z / Rspherical
sint = r / Rspherical

! Eq. 85 in Tereby et al.:
Bracket4 = (0.5d0 * m0_Shu * c0 * t0)**4
zeta = Omega0**2 * Bracket4 / (Gconst * Mstar * Rspherical)

! Get cos(theta_0) from cos(theta):
!print *, ' calling CM_root for zeta = ', zeta, ' cost = ', cost
call CM_root (zeta, cost, cost0)
!print *, 'returned from CM_solve with cost0 = ', cost0
!read(5, *)

!!$in_cavity = .false.
!!$
!!$if (include_cavity) then
!!$   if (z .gt. 0.0d0) then
!!$      in_cavity = (cost .gt. cos_theta_cavity)
!!$   else
!!$      in_cavity = (cost .gt. cos_theta_cavity)
!!$   end if
!!$end if
!!$
!!$! This will be used to make the density continuous.
!!$if (in_cavity) then
!!$   cost = cost_theta_cavity
!!$   call CM_root (zeta, cost, cost0)
!!$end if

! Velocity in spherical coordinates:
fac1      = SQRT (Gconst * Mstar / Rspherical)
cos_ratio = cost / cost0
fac2      = sqrt(1.d0 + cos_ratio)
fac2m     = sqrt(1.d0 - cos_ratio)

uRspherical  = - fac1 * fac2
utheta       =   fac1 * (cost0 - cost) / sint * fac2

! We choose the positive sign because sin(theta) is positive
! for theta between 0 and pi.
sint0  =   SQRT(1.d0 - cost0**2)
uphi   =   fac1 * sint0 / sint * fac2m

! Convert to cylindrical coordinates:
ur = uRspherical*sint + utheta*cost
uz = uRspherical*cost - utheta*sint

! Density:

fac1 = - M0_dot / (4.d0 * pi * Rspherical**2 * uRspherical)
P2   = 0.5d0 * (3.d0*cost0**2 - 1.d0)
fac2 = 1.d0/(1.d0 + 2.d0*zeta*P2)
rho  = fac1 * fac2

pressure = rho * Rgas * Tcloud
eint     = pressure/gm1

! Verify that the mechanical energy is zero.
KE = 0.5d0 * (ur**2 + uz**2 + uphi**2)
PE_abs = GM/Rspherical
relative_error = ABS(KE - PE_abs) / PE_abs
if (relative_error .gt. 1.d-6) then
   print *, ' In CM_flow, the total mechanical energy is not zero'
   print *, ' KE = ', KE, ' PE = ', PE_abs
   stop
end if

!!$if (in_cavity) then
!!$   ur   = 0.d0
!!$   uz   = 0.d0
!!$   uphi = 0.d0
!!$end if

return
end

!----------------------------------------------------------------------------------85

subroutine CM_root (zeta, cost, cost0)

! Finds the root to the (cubic) equation for the Ulrich-Cassen-Moosman initial
! condition.  Given cost, find cost0.

implicit none
real(8) :: zeta, cost, cost0

real(8) :: cost_c, zeta_c, CM_func
common /CM_func_stuff/ cost_c, zeta_c
external CM_func

! Local:
integer, parameter :: nint_CM = 100
integer :: nroots, int
real(8) :: delx, x1, x2, f1, f2, xa, xb
real(8), parameter :: tol = 1.0d-12

! Set-up common stuff to give to CM_func:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Since the flow is symmetric about the mid-plane, we will find a solution
! on the upper side of the mid-plane and then reflect it if required.
! To make sense of the absolute value,  note that cos(theta) is the projection of
! the particle on the vertical axis.
cost_c = ABS(cost)
zeta_c = zeta

!print *, ' CM_root has been called'
!print *, ' cost = ', cost, ' zeta = ', zeta

! Find an interval containing the solution:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nroots = 0
delx = 1.d0 / nint_CM
x1 = 0.0d0
f1 = CM_func(x1)

do int = 1, nint_CM
   x2 = int * delx
   f2 = CM_func(x2)

   if (f1*f2 .lt. 0) then
      nroots = nroots + 1
      xa = x1
      xb = x2
    end if

    x1 = x2
    f1 = f2
end do

if (nroots .eq. 0) then
   print *, ' Subroutine CM_root:'
   print *, ' No zero crossings found'
   print *, ' zeta = ', zeta, ' cost = ', cost
   print *, ' If this is a point in the mid-plane, then'
   print *, ' I suggest that you choose an even number'
   print *, ' of points to avoid the mid-plane.'
   stop
else if (nroots .gt. 1) then
   print *, ' Subroutine CM_root:'
   print *, ' More than one zero crossing found'
   stop
end if

!print *, ' n zero crossings = ', nroots
!print *, ' xa = ', xa, ' xb = ', xb

call zeroin (xa, xb, CM_func, tol, cost0)
!print *, ' cost0 = ', cost0

! Reflect the solution if required.
! Since a particle now in the lower half originated in the lower
! half, transfer sign of cost to cost0:
cost0 = SIGN (cost0, cost)

!print *, ' returning from CM_root'

return
end

!----------------------------------------------------------------------------------85

real(8) function CM_func (cost0)

!     This is the function we need to zero to get the 
!     co-latitude at infinity for an in-falling particle
!     in the Cassen-Moosman flow.  See Shu & Tereby, p. 544

implicit none
real(8) :: cost0

real(8) :: cost_c, zeta_c
common /CM_func_stuff/ cost_c, zeta_c

! Local:
real(8) :: sin2t0, num, denom

sin2t0 = 1.d0 - cost0**2
num    = cost0 - cost_c
denom  = sin2t0 * cost0
CM_func = zeta_c * denom - num

! print *, ' cost_c = ', cost_c, ' zeta_c = ', zeta_c
! print *, ' rnum = ', rnum
! print *, ' cost0 = ', cost0, ' CM_func = ', CM_func
return
end function CM_func

!----------------------------------------------------------------------------------85

real(8) function rise_function(r, r_begin_rise, r_end_rise)

implicit none
real(8) :: r, r_begin_rise, r_end_rise

! Local:
real(8) :: xi, paren, bracket
real(8), parameter :: kappa = 0.5d0 * exp(2.d0) * log(2.0d0)

if (r .lt. r_begin_rise) then
   rise_function = 0.0d0
else if (r .lt. r_end_rise) then
   xi        = (r - r_begin_rise) / (r_end_rise - r_begin_rise)
   paren    = 1.d0 / (xi - 1)
   bracket  = - kappa/xi * exp(paren)
   rise_function = exp(bracket)
else
   rise_function = 1.0d0
end if

end function rise_function

!----------------------------------------------------------------------------------85

subroutine set_post_shock_layer_temperature(q)

use grid
use partition_data
use artificial_pressure_module
use physical_constants
use Cassen_Moosman_infall
use dof_indices
use thermal_parameters
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
integer :: iphi, ir, iz
real(8) :: max_compressive_dil
real(8), dimension(nr) :: z_shock_bot, z_shock_top
integer, dimension(nr) :: iz_shock_bot, iz_shock_top
real(8) :: rho, pressure

! Detect shock:
do iphi = sphi, ephi
   do ir = sr, er
      max_compressive_dil = 0.0d0
      do iz = 1, nz/2
         if (-dil(ir, iphi, iz) .gt. max_compressive_dil) then
            max_compressive_dil = -dil(ir, iphi, iz)
            z_shock_bot(ir) = zgrid(iz)
            iz_shock_bot(ir) = iz
         end if
      end do

      max_compressive_dil = 0.0d0
      do iz = nz/2, nz
         if (-dil(ir, iphi, iz) .gt. max_compressive_dil) then
            max_compressive_dil = -dil(ir, iphi, iz)
            z_shock_top(ir) = zgrid(iz)
            iz_shock_top(ir) = iz            
         end if
      end do      
   end do
end do

open (unit = 1, file = 'z_shock.dat', form = 'formatted', status = 'unknown')
do ir = 1, nr
   write(1, "(4(1x, e12.5))") rgrid(ir)/AU, z_shock_bot(ir)/AU, z_shock_top(ir)/AU, T_psl(ir)
end do
close(1)
print *, ' wrote z_shock.dat'

do iphi = sphi, ephi
   do ir = sr, er
      iz = iz_shock_top(ir) - 1
      rho = q(ir, iphi, iz, irho)
      pressure = rho * Rgas * T_psl(ir)
      q(ir, iphi, iz, ener) = pressure/gm1

      iz = iz_shock_bot(ir) + 1
      rho = q(ir, iphi, iz, irho)
      pressure = rho * Rgas * T_psl(ir)
      q(ir, iphi, iz, ener) = pressure/gm1      
   end do
end do

end subroutine set_post_shock_layer_temperature

!----------------------------------------------------------------------------------85

subroutine store_post_cooling_layer_temperature

use grid
use gravity
use physical_constants
use math_constants
use Cassen_Moosman_infall
use total_allocated_words, only: n_words_allocated
implicit none
real(8) :: r

! Local:
real(8) :: sin_theta0, cos_theta0, tan_theta0
real(8) :: fac, rho, ur, uz, uphi
real(8) :: ke, pe_abs, rel_diff
integer :: ir

allocate(T_psl(nr))
n_words_allocated = n_words_allocated + nr

do ir = 1, nr
   r = rgrid(ir)
   if (r .lt. rd) then
      fac = SQRT(GM / r)
      sin_theta0 = SQRT(r/rd)
      cos_theta0 = SQRT(1.d0 - r/rd)
      tan_theta0 = sin_theta0 / cos_theta0

      ur   = - fac
      uz   = - fac * cos_theta0
      uphi =   fac * sin_theta0
      rho  = M0_dot / (8.d0 * pi * r**2) / fac * tan_theta0**2

      T_psl(ir) = ( (0.5d0 * rho * abs(uz)**3) / sigma_SB ) ** (0.25d0)

      ! rel_diff should be zero:
      ke     = 0.5d0 * (ur**2 + uphi**2 + uz**2)
      pe_abs = GM/r
      rel_diff = (ke - pe_abs) / ke
      if (rel_diff .gt. 1.0d-6) then
         print *, ' r/rd = ', r/rd
         print *, ' ke = ', 0.5d0 * (ur**2 + uphi**2 + uz**2), ' pe = ', -GM/r 
         stop
      end if
   else
      T_psl(ir) = Tcloud
   end if
end do

open (unit = 1, file = 'T_post_shock_cooling_layer.dat', form = 'formatted', status = 'unknown')
do ir = 1, nr
   write(1, "(2(1x, e12.5))") rgrid(ir), T_psl(ir)
end do
close(1)

end subroutine store_post_cooling_layer_temperature

!----------------------------------------------------------------------------------85

subroutine get_cos_theta_0_max (rmin_infall_at_midplane, cos_theta_0_max)

! Determines the maximum cos(theta0) we can allow (i.e., how close we can get to the
! pole) if we don't want to allow any infall inward of rmin_infall_at_midplane.

use Cassen_Moosman_infall
use physical_constants
implicit none
real(8) :: rmin_infall_at_midplane, cos_theta_0_max

! Local:
real(8) :: z, Rspherical, cost, Bracket4, zeta

z = 0.d0

! Spherical coordinates:
Rspherical = SQRT(z**2 + rmin_infall_at_midplane**2)
cost = z / Rspherical

! Eq. 85 in Tereby et al.:
Bracket4 = (0.5d0 * m0_Shu * c0 * t0)**4
zeta = Omega0**2 * Bracket4 / (Gconst * Mstar * Rspherical)

call CM_root (zeta, cost, cos_theta_0_max)

end subroutine get_cos_theta_0_max

!----------------------------------------------------------------------------------85


