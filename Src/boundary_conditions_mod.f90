!----------------------------------------------------------------------------------85

! Note for non-reflective BC: The flag ibalanced (0 or 1) when set to zero gives
! sets the rate of incoming waves to zero.  When set to 1, we say the rate of
! change of the the incoming waves is set by hydrostatic balance
! in z and centrifugal balance in r.  The flag "ibalanced" sits in
! module boundary_conditions.

!----------------------------------------------------------------------------------85

module boundary_condition_types
   ! Boundary condition types to improve code readability:
   integer, parameter :: null = 0, non_reflective = 1, zero_normal_momentum = 2, &
     Cassen_Moosman_BC = 3, &
     outflow = 5, periodic = 6, viscous_wall = 7, zero_shear_stress = 8, &
     z_dirichlet = 9, hold_basic_state = 10
end module boundary_condition_types

!----------------------------------------------------------------------------------85   

module boundary_conditions
   ! ibalanced is either 0 or 1.  1: Incoming waves maintain force balance.
   !                              0: The rate of change of all incoming waves is zeroed.
   !                                 This can cause loss of hydrostatic (in z) and centrifugal balance
   !                                 (in r).
   integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced
   
   ! These parameters are for viscous wall BC at rmin and rmax:
   ! rotation rate, axial velocity and sound speed squared (proportional to temperature)
   real(8) :: Omega_rmin, Omega_rmax, uz_rmin, uz_rmax, c_sound_rmin, c_sound_rmax
   logical :: specify_viscous_wall_conditions_was_called
   real(8), allocatable :: p_BC_zmin  (:), p_BC_zmax  (:) ! dimension will be nr
   real(8), allocatable :: p_BC_rmin  (:), p_BC_rmax  (:) ! dimension will be nz
   real(8), allocatable :: ci2_BC_rmin(:), ci2_BC_rmax(:) ! dimension will be nz

   real(8), allocatable :: qLeft_z_BC(:), qRight_z_BC(:)

!----------------------------------------------------------------------------------85

contains   

!----------------------------------------------------------------------------------85

subroutine enforce_BC (q)

use grid
use partition_data
use dof_indices
use boundary_condition_types
use thermal_parameters, only: ci_squared_initial, isothermal, gamma, gm1
use basic_state
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
integer :: iphi, ir, iz, idof
real(8) :: rho2_q2, pressure, ur, Mach_r
real(8) :: c_sound_at_point ! called function

real(8) :: uphi_wall, uphi_adjacent, uphi_prime_adjacent
integer :: ira, iza ! For the grid point adjacent to the wall.

! r min BC:
!~~~~~~~~~~
if (sr .eq. 1) then ! This checks if this proc has the boundary.
   ir   = 1
   if (rmin_BC .eq. zero_normal_momentum) then
      do iz = 1, nz
         do iphi = sphi, ephi
            q(ir, iphi, iz, rmom) = 0.0d0
         end do
      end do
   else if (rmin_BC .eq. outflow) then
      do iz = 1, nz
         do iphi = sphi, ephi
            if (q(ir, iphi, iz, rmom) .gt. 0.0d0) q(ir, iphi, iz, rmom) = 0.0d0
         end do
      end do
   else if (rmin_BC .eq. viscous_wall) then
      do iz = 1, nz
         do iphi = sphi, ephi
            q(ir, iphi, iz, rmom) = 0.0d0
            q(ir, iphi, iz, zmom) = uz_rmin
            q(ir, iphi, iz, amom) = q(ir,iphi,iz,irho) * Omega_rmin*rgrid(ir) * rgrid(ir)
            if (.not. isothermal) then
               ! Specified wall temperature (in the form of a sound speed).
               pressure = c_sound_rmin**2 * q(ir,iphi,iz,irho) / gamma
               q(ir,iphi,iz,ener) = pressure / gm1
            end if
         end do
      end do
   else if ((rmin_BC .eq. hold_basic_state) .and. (have_basic_state)) then
      do iz = 1, nz
         do iphi = sphi, ephi
            q(ir, iphi, iz, irho) = rho_basic(ir,iz)
            q(ir, iphi, iz, rmom) = 0.d0
            q(ir, iphi, iz, zmom) = 0.d0
            q(ir, iphi, iz, amom) = amom_basic(ir,iz)
            if (.not. isothermal) then
               q(ir,iphi,iz,ener) = eint_basic(ir, iz)
            end if
         end do
      end do
   end if
end if

! r max BC:
!~~~~~~~~~~
if (er .eq. nr) then
   ir  = nr
   if (rmax_BC .eq. zero_normal_momentum) then
      do iz = 1, nz
         do iphi = sphi, ephi
            q(ir, iphi, iz, rmom) = 0.0d0
         end do
      end do
   else if (rmax_BC .eq. outflow) then
      do iz = 1, nz
         do iphi = sphi, ephi
            if (q(ir, iphi, iz, rmom) .lt. 0.0d0) q(ir, iphi, iz, rmom) = 0.0d0
         end do
      end do      
   else if (rmax_BC .eq. viscous_wall) then
      do iz = 1, nz
         do iphi = sphi, ephi
            q(ir, iphi, iz, rmom) = 0.0d0
            q(ir, iphi, iz, zmom) = uz_rmax
            q(ir, iphi, iz, amom) = q(ir,iphi,iz,irho) * Omega_rmax*rgrid(ir) * rgrid(ir)
            if (.not. isothermal) then
               ! Specified wall temperature (in the form of a sound speed).
               pressure = c_sound_rmax**2 * q(ir,iphi,iz,irho) / gamma
               q(ir,iphi,iz,ener) = pressure / gm1
            end if
         end do
      end do
   else if ((rmax_BC .eq. hold_basic_state) .and. (have_basic_state)) then
      do iz = 1, nz
         do iphi = sphi, ephi
            q(ir, iphi, iz, irho) = rho_basic(ir,iz)
            q(ir, iphi, iz, rmom) = 0.d0
            q(ir, iphi, iz, zmom) = 0.d0
            q(ir, iphi, iz, amom) = amom_basic(ir,iz)
            if (.not. isothermal) then
               q(ir,iphi,iz,ener) = eint_basic(ir, iz)
            end if
         end do
      end do
   end if      
end if

! z min BC:
!~~~~~~~~~~
iz = 1
if (zmin_BC .eq. zero_normal_momentum) then
   do iphi = sphi, ephi
      do ir = sr, er
         q(ir, iphi, iz, zmom) = 0.0d0
      end do
   end do
else if ((zmin_BC .eq. hold_basic_state) .and. (have_basic_state)) then
   do iphi = sphi, ephi
      do ir = sr, er
         q(ir, iphi, iz, irho) = rho_basic(ir, iz)
         q(ir, iphi, iz, rmom) = 0.d0
         q(ir, iphi, iz, zmom) = 0.d0
         q(ir, iphi, iz, amom) = amom_basic(ir, iz)
         if (.not. isothermal) then
            q(ir, iphi, iz, ener) = eint_basic(ir, iz)
         end if
      end do
   end do   
end if

iz = nz
if (zmax_BC .eq. zero_normal_momentum) then
   do iphi = sphi, ephi
      do ir = sr, er
         q(ir, iphi, iz, zmom) = 0.0d0
      end do
   end do
else if ((zmax_BC .eq. hold_basic_state) .and. (have_basic_state)) then
   do iphi = sphi, ephi
      do ir = sr, er
         q(ir, iphi, iz, irho) = rho_basic(ir, iz)
         q(ir, iphi, iz, rmom) = 0.d0
         q(ir, iphi, iz, zmom) = 0.d0
         q(ir, iphi, iz, amom) = amom_basic(ir, iz)
         if (.not. isothermal) then
            q(ir, iphi, iz, ener) = eint_basic(ir, iz)
         end if
      end do
   end do   
end if

if ((rmax_BC .eq. Cassen_Moosman_BC) .and. (zmin_BC .eq. Cassen_Moosman_BC) .and. &
    (zmax_BC .eq. Cassen_Moosman_BC)) then
   call apply_Cassen_Moosman_boundary_values(q)
end if

if (zmin_BC .eq. z_dirichlet) then
   do idof = 1, ndof
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, 1, idof) = qLeft_z_BC(idof)
            ! print *, ' idof = ', idof, ' qL = ', qLeft_z_BC(idof)
         end do
      end do
   end do
end if

if (zmax_BC .eq. z_dirichlet) then
   do idof = 1, ndof
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, nz, idof) = qRight_z_BC(idof)
            ! print *, ' idof = ', idof, ' qR = ', qRight_z_BC(idof)
         end do
      end do
   end do
end if

end subroutine enforce_BC

!----------------------------------------------------------------------------------85

subroutine z_adiabatic_non_reflective_bc(iz, q, uz, dF, qdot)

! Modifies qdot to implement non-reflective BC in the z direction with the assumption
! that incoming waves are hydrostatic and isentropic.
! See notes dated Nov. 7, 2018 and numbered as PZ.  Primitive variables are used
! for the characteristic decomposition.

use thermal_parameters
use partition_data
use dof_indices
use gravity
use grid
implicit none
integer, intent(in) :: iz
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q, qdot
real(8), dimension(sr:er, sphi:ephi, nz)       :: uz
real(8), dimension(sr:er, sphi:ephi, ndof)     :: dF
real(8) :: dt

! Local:
! One sided vertical derivatives of the primitive variables at the boundaries:
real(8), dimension(sr:er, sphi:ephi, ndof) :: d_prim_dz

! Values of primitive variables at the end-points in order to take finite differences:
real(8), dimension(sr:er, sphi:ephi, 1:4,     ndof) :: prim_left
real(8), dimension(sr:er, sphi:ephi, nz-3:nz, ndof) :: prim_right
integer :: iz1, ir, iphi, idof

! Temps:
real(8) :: drho_dz_hydrostatic, dpdz_hydrostatic
real(8) :: cs, uz1, ur, uphi, rho, L1, L2, L3, L4, L5, uz_p_cs, uz_m_cs
integer :: nhat_z ! Direction of unit normal to the z boundary
real(8), dimension(ndof) :: dF_non_reflective, A_dUdz

! Functions called:
real(8) :: c_sound_at_point, pressure_at_point

if (iz .eq. nz) then
   nhat_z = 1
   ! Get primitive variables at last 4 z grid points:
   do iphi = sphi, ephi
      do ir = sr, er
         ! We use iz1 since iz is an input variable:
         do iz1 = nz - 3, nz
            prim_right(ir, iphi, iz1, irho ) = q(ir, iphi, iz1, irho)
            prim_right(ir, iphi, iz1, iuz  ) = q(ir, iphi, iz1, zmom) / q(ir, iphi, iz1, irho)
            prim_right(ir, iphi, iz1, iur  ) = q(ir, iphi, iz1, rmom) / q(ir, iphi, iz1, irho)
            prim_right(ir, iphi, iz1, iuphi) = q(ir, iphi, iz1, amom) / q(ir, iphi, iz1, irho) / rgrid(ir)
            prim_right(ir, iphi, iz1, ip   ) = q(ir, iphi, iz1, ener) * gm1
         end do
      end do
   end do

   ! z finite-difference of primitive variables:
   do idof = 1, 5
      do iphi = sphi, ephi
         do ir = sr, er
            d_prim_dz(ir,iphi,idof) = (11.d0*prim_right(ir,iphi,nz,  idof) - 18.d0*prim_right(ir,iphi,nz-1,idof) + &
                                        9.d0*prim_right(ir,iphi,nz-2,idof) - 2.0d0*prim_right(ir,iphi,nz-3,idof)) / 6.d0 * &
                                        Ji_z(nz)
         end do
      end do
   end do
else if (iz .eq. 1) then
   nhat_z = -1
   ! Get primitive variables at first 4 grid points:
   do iphi = sphi, ephi
      do ir = sr, er
         do iz1 = 1, 4
            prim_left(ir, iphi, iz1, irho ) = q(ir, iphi, iz1, irho)
            prim_left(ir, iphi, iz1, iuz  ) = uz(ir, iphi, iz1)
            prim_left(ir, iphi, iz1, iur  ) = q(ir, iphi, iz1, rmom) / q(ir, iphi, iz1, irho)
            prim_left(ir, iphi, iz1, iuphi) = q(ir, iphi, iz1, amom) / q(ir, iphi, iz1, irho) / rgrid(ir)
            prim_left(ir, iphi, iz1, ip   ) = q(ir, iphi, iz1, ener) * gm1
         end do
      end do
   end do

   ! z-finite_difference of primitive variables:
   do idof = 1, 5 
      do iphi = sphi, ephi
         do ir = sr, er
            d_prim_dz(ir,iphi,idof) = (2.d0*prim_left(ir,iphi,4,idof) -  9.d0*prim_left(ir,iphi,3,idof) + &
                                      18.d0*prim_left(ir,iphi,2,idof) - 11.d0*prim_left(ir,iphi,1,idof)) / 6.d0 * &
                                      Ji_z(1)
         end do
      end do
   end do
else
   nhat_z = 0 ! To avoid "May be uninitialized" warning from the compiler.
   print *, ' In z_adiabatic_primitive_non_reflective_bc invalid iz = ', iz
   call terminate_with_no_save(1)
end if

do iphi = sphi, ephi
   do ir = sr, er
      cs   = c_sound_at_point(q, ir, iphi, iz)
      uz1  = uz(ir, iphi, iz)
      rho  = q(ir, iphi, iz, irho)
      ur   = q(ir, iphi, iz, rmom) / rho
      uphi = q(ir, iphi, iz, amom) / rho / rgrid(ir)
      dpdz_hydrostatic    = rho * gz(ir, iz)
      drho_dz_hydrostatic = dpdz_hydrostatic / cs**2
      uz_m_cs = uz1 - cs
      uz_p_cs = uz1 + cs
      if (uz_m_cs*nhat_z .lt. 0.d0) then
         L1 = ibalanced*uz_m_cs*(                                     + 0.5d0/cs**2*dpdz_hydrostatic     ) ! incoming
      else
         L1 =           uz_m_cs*(-0.5d0*rho/cs*d_prim_dz(ir,iphi,iuz) + 0.5d0/cs**2*d_prim_dz(ir,iphi,ip)) ! outgoing
      end if

      if (uz_p_cs*nhat_z .lt. 0.d0) then
         L2 = ibalanced*uz_p_cs*(                                     + 0.5d0/cs**2*dpdz_hydrostatic     ) ! incoming         
      else
         L2 =           uz_p_cs*(+0.5d0*rho/cs*d_prim_dz(ir,iphi,iuz) + 0.5d0/cs**2*d_prim_dz(ir,iphi,ip)) ! outgoing 
      end if

      if (uz1*nhat_z .lt. 0.d0) then
         L3 = ibalanced*uz1*(drho_dz_hydrostatic      - 1.d0/cs**2*dpdz_hydrostatic     ) ! incoming
         L4 = 0.0d0 ! incoming
         L5 = 0.0d0 ! incoming
      else
         L3 =  uz1 * (d_prim_dz(ir,iphi,irho ) - 1.d0/cs**2*d_prim_dz(ir,iphi,ip)) ! outgoing
         L4 =  uz1 *  d_prim_dz(ir,iphi,iur  )                                     ! outgoing
         L5 =  uz1 *  d_prim_dz(ir,iphi,iuphi)                                     ! outgoing
      end if

      ! S*L:
      A_dUdz(irho ) = L1 + L2 + L3
      A_dUdz(iuz  ) = cs/rho*(-L1 + L2) 
      A_dUdz(iur  ) = L4
      A_dUdz(iuphi) = L5
      A_dUdz(ip   ) = cs**2 * (L1 + L2)

      ! For conservative variables:
      dF_non_reflective(irho) = A_dUdz(irho)
      dF_non_reflective(zmom) =  rho*A_dUdz(iuz  ) + uz1 *A_dUdz(irho)
      dF_non_reflective(rmom) =  rho*A_dUdz(iur  ) + ur  *A_dUdz(irho)
      dF_non_reflective(amom) = (rho*A_dUdz(iuphi) + uphi*A_dUdz(irho)) * rgrid(ir)
      dF_non_reflective(ener) = A_dUdz(ip) / gm1

      ! Remove what you stored in favor of the non-reflective flux derivative:
      do idof = 1, 5
         qdot(ir, iphi, iz, idof) = qdot(ir, iphi, iz, idof) + dF(ir, iphi, idof) - &
              dF_non_reflective(idof)
      end do
   end do
end do

end subroutine z_adiabatic_non_reflective_bc 

!----------------------------------------------------------------------------------85

subroutine z_isothermal_non_reflective_bc(iz, q, uz, dF, qdot)

! Modifies qdot to implement a non-reflective BC in the z direction.
! See notes dated April 26, 2018.
! "Non-reflective BC using conservative variables (isothermal z)."

use thermal_parameters
use partition_data
use dof_indices
use gravity
use grid
implicit none
integer :: iz
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q, qdot
real(8), dimension(sr:er, sphi:ephi, nz)       :: uz
real(8), dimension(sr:er, sphi:ephi, ndof)     :: dF
real(8) :: dt

! Local:
! One sided vertical derivatives of the q variables at the boundaries:
real(8), dimension(sr:er, sphi:ephi, ndof) :: dq_dz
integer :: ir, iphi, idof
! Temps:
real(8) :: drho_dz_hydrostatic, ddz_zmom_hydrostatic, ddz_rmom_hydrostatic, &
           ddz_uphi_r_hydrostatic
real(8) :: ci, uz1, ur, uphi_r, L1, L2, L3, L4, uz_p_ci, uz_m_ci
integer :: nhat_z ! Direction of unit normal to the z boundary
real(8), dimension(ndof) :: dF_non_reflective

real(8) :: Mach_z_incoming, uz_target, zmom_target, Mach_z_outgoing

!print *, ' subroutine z_isothermal_non_reflective_bc called'
!read (5, *)

if (iz .eq. nz) then
   nhat_z = 1
   do idof = 1, 4
      do iphi = sphi, ephi
         do ir = sr, er
            dq_dz(ir,iphi,idof) = (11.d0*q(ir,iphi,nz,  idof) - 18.d0*q(ir,iphi,nz-1,idof) + &
                                    9.d0*q(ir,iphi,nz-2,idof) - 2.0d0*q(ir,iphi,nz-3,idof)) / 6.d0 * &
                                    Ji_z(nz)
         end do
      end do
   end do
else if (iz .eq. 1) then
   nhat_z = -1
   do idof = 1, 4 
      do iphi = sphi, ephi
         do ir = sr, er
            dq_dz(ir,iphi,idof) = (2.d0*q(ir,iphi,4,idof) -  9.d0*q(ir,iphi,3,idof) + &
                                  18.d0*q(ir,iphi,2,idof) - 11.d0*q(ir,iphi,1,idof)) / 6.d0 * &
                                  Ji_z(1)
         end do
      end do
   end do
else
   nhat_z = 0 ! To avoid "May be uninitialized" warning from the compiler.
   print *, ' In z_isothermal_non_reflective, invalid iz = ', iz
   call terminate_with_no_save(1)
end if

do iphi = sphi, ephi
   do ir = sr, er
      ci  = SQRT(ci_squared_initial(ir, iz))
      uz1    = uz(ir, iphi, iz)
      ur     = q(ir, iphi, iz, rmom) / q(ir, iphi, iz, irho)
      uphi_r = q(ir, iphi, iz, amom) / q(ir, iphi, iz, irho)
      
      drho_dz_hydrostatic    = q (ir, iphi, iz, irho) * gz(ir, iz) / ci_squared_initial(ir, iz)
      ddz_zmom_hydrostatic   = uz1    * drho_dz_hydrostatic
      ddz_rmom_hydrostatic   = ur     * drho_dz_hydrostatic
      ddz_uphi_r_hydrostatic = uphi_r * drho_dz_hydrostatic      
      
      uz_m_ci = uz1 - ci
      uz_p_ci = uz1 + ci
      if (uz_m_ci*nhat_z .lt. 0.d0) then
         L1 = ibalanced*uz_m_ci/(2.d0*ci) * (uz_p_ci*drho_dz_hydrostatic - ddz_zmom_hydrostatic) ! incoming
      else
         L1 =           uz_m_ci/(2.d0*ci) * (uz_p_ci*dq_dz(ir,iphi,irho) - dq_dz(ir,iphi,zmom) ) ! outgoing
      end if

      if (uz_p_ci*nhat_z .lt. 0.d0) then
         L2 = ibalanced*uz_p_ci/(2.d0*ci) * (-uz_m_ci*drho_dz_hydrostatic + ddz_zmom_hydrostatic) ! incoming
      else
         L2 =           uz_p_ci/(2.d0*ci) * (-uz_m_ci*dq_dz(ir,iphi,irho) + dq_dz(ir,iphi,zmom) ) ! outgoing
      end if

      if (uz1*nhat_z .lt. 0.d0) then
         L3 = ibalanced*uz1 * (-ur     * drho_dz_hydrostatic + ddz_rmom_hydrostatic  ) ! incoming
         L4 = ibalanced*uz1 * (-uphi_r * drho_dz_hydrostatic + ddz_uphi_r_hydrostatic) ! incoming
      else
         L3 =  uz1 * (    -ur*dq_dz(ir,iphi,irho) + dq_dz(ir,iphi,rmom)) ! outgoing
         L4 =  uz1 * (-uphi_r*dq_dz(ir,iphi,irho) + dq_dz(ir,iphi,amom)) ! outgoing
      end if

      ! S*L:
      dF_non_reflective(irho) =       L1 + L2
      dF_non_reflective(rmom) = ur * (L1 + L2) + L3
      dF_non_reflective(zmom) = (uz1 - ci)*L1 + (uz1 + ci)*L2
      dF_non_reflective(amom) = uphi_r*(L1 + L2) + L4

      ! Remove what you stored in favor of the non-reflective flux derivative:
      do idof = 1, 4
         qdot(ir, iphi, iz, idof) = qdot(ir, iphi, iz, idof) + dF(ir, iphi, idof) - &
              dF_non_reflective(idof)
      end do

!!$      Mach_z_incoming = -uz1/ci * nhat_z
!!$      if ( (Mach_z_incoming .gt. 0.3d0) .and. (.not. isnan(dt)) ) then
!!$         uz_target = -0.3d0 * ci / nhat_z
!!$         zmom_target = q(ir, iphi, iz, irho) * uz_target
!!$         qdot(ir,iphi,iz,zmom) = qdot(ir,iphi,iz,zmom) + &
!!$              (q(ir,iphi,iz,zmom) - zmom_target)/(20.d0 * dt)
!!$         print *, ' dt = ', dt, ' ir = ', ir, ' Mach_z_incoming = ', Mach_z_incoming
!!$         ! read (5, *)
!!$      end if
!!$
!!$      Mach_z_outgoing = uz1/ci * nhat_z
!!$      if ( (Mach_z_outgoing .gt. 0.3d0) .and. (.not. isnan(dt)) ) then
!!$         uz_target = 0.3d0 * ci / nhat_z
!!$         zmom_target = q(ir, iphi, iz, irho) * uz_target
!!$         qdot(ir,iphi,iz,zmom) = qdot(ir,iphi,iz,zmom) + &
!!$              (q(ir,iphi,iz,zmom) - zmom_target)/(20.d0 * dt)
!!$         print *, ' dt = ', dt, ' ir = ', ir, ' Mach_z_outgoing = ', Mach_z_outgoing
!!$         ! read (5, *)
!!$      end if      
   end do
end do

end subroutine z_isothermal_non_reflective_bc

!----------------------------------------------------------------------------------85

subroutine r_isothermal_non_reflective_bc(ir, q, ur, dF)

! Modifies qdot to implement a non-reflective BC in the r direction.  We require
! that incoming waves satisfy centrifigal balance.
! See notes dated Oct 25, 2018.
! "Non-reflective BC in r (isothermal)"

use thermal_parameters
use partition_data
use dof_indices
use gravity
use grid
implicit none
integer :: ir
real(8), dimension(sphi:ephi, sz_r:ez_r, ndof, nr) :: q
real(8), dimension(sphi:ephi, sz_r:ez_r,       nr) :: ur
real(8), dimension(sphi:ephi, sz_r:ez_r, ndof, nr) :: dF

! Local:
! One sided r derivatives of the q variables at the boundaries:
real(8), dimension(sphi:ephi, sz_r:ez_r, ndof) :: dq_dr
integer :: iz, iphi, idof
! Temps:
real(8) :: drho_dr_centrifugal, ddr_zmom_centrifugal, ddr_rmom_centrifugal, &
           ddr_uphi_r_centrifugal, dPdr_centrifugal, d_ci_dr
real(8) :: ci, ur1, uz, uphi_r, L1, L2, L3, L4, ur_p_ci, ur_m_ci
integer :: nhat_r ! Direction of unit normal to the r boundary
real(8), dimension(ndof) :: dF_non_reflective

real(8) :: uphi_squared_over_r, accel

! Debug:
integer :: n_in, iwhich

if (ir .eq. nr) then
   nhat_r = 1
   do idof = 1, 4
      do iz = sz_r, ez_r
         do iphi = sphi, ephi
            dq_dr(iphi,iz,idof) = (11.d0*q(iphi,iz,idof,nr  ) - 18.d0*q(iphi,iz,idof,nr-1) + &
                                    9.d0*q(iphi,iz,idof,nr-2) - 2.0d0*q(iphi,iz,idof,nr-3)) / 6.d0 * &
                                    Ji_r(nr)
         end do
      end do
   end do
else if (ir .eq. 1) then
   nhat_r = -1
   do idof = 1, 4 
      do iz = sz_r, ez_r
         do iphi = sphi, ephi
            dq_dr(iphi,iz,idof) = (2.d0*q(iphi,iz,idof,4) -  9.d0*q(iphi,iz,idof,3) + &
                                  18.d0*q(iphi,iz,idof,2) - 11.d0*q(iphi,iz,idof,1)) / 6.d0 * &
                                  Ji_r(1)
         end do
      end do
   end do
else
   nhat_r = 0 ! To avoid "May be uninitialized" warning from the compiler.
   print *, ' In r_isothermal_non_reflective, invalid ir = ', ir
   call terminate_with_no_save(1)
end if

do iz = sz_r, ez_r
   do iphi = sphi, ephi
      ci  = SQRT(ci_squared_initial(ir, iz))
      ur1    = ur(iphi, iz, ir)
      uz     = q(iphi, iz, zmom, ir) / q(iphi, iz, irho, ir)
      uphi_r = q(iphi, iz, amom, ir) / q(iphi, iz, irho, ir)
      uphi_squared_over_r = uphi_r**2 / rgrid(ir)**3

      accel = gr(ir, iz) + uphi_squared_over_r

      if (ir .eq. 1 ) d_ci_dr = d_ci_dr_inner(iz)
      if (ir .eq. nr) d_ci_dr = d_ci_dr_outer(iz)
      dPdr_centrifugal       = q (iphi, iz, irho, ir) * accel
      drho_dr_centrifugal = (dPdr_centrifugal - 2.d0*q(iphi,iz,irho,ir)*ci*d_ci_dr) / ci_squared_initial(ir,iz)
      ddr_rmom_centrifugal   = ur1    * drho_dr_centrifugal ! assumes d ur/dr = 0
      ddr_zmom_centrifugal   = uz     * drho_dr_centrifugal ! assumes d uz/dr = 0
      ddr_uphi_r_centrifugal = uphi_r * drho_dr_centrifugal ! assumes d (uphi*r)/dr = 0
      
      ur_m_ci = ur1 - ci
      ur_p_ci = ur1 + ci
      n_in = 0
      iwhich = 0
      if (ur_m_ci*nhat_r .lt. 0.d0) then
         n_in = n_in + 1
         iwhich = 1
         ! if (iz .eq. nz/2) print *, ' L1 is incoming'
         L1 = ibalanced*ur_m_ci/(2.d0*ci) * (ur_p_ci*drho_dr_centrifugal - ddr_rmom_centrifugal) ! incoming
      else
         ! if (iz .eq. nz/2) print *, ' L1 is outgoing'
         L1 =           ur_m_ci/(2.d0*ci) * (ur_p_ci*dq_dr(iphi,iz,irho) - dq_dr(iphi,iz,rmom) ) ! outgoing
      end if

      if (ur_p_ci*nhat_r .lt. 0.d0) then
         n_in = n_in + 1
         iwhich = 2
         ! if (iz .eq. nz/2) print *, ' L2 is incoming'         
         L2 = ibalanced*ur_p_ci/(2.d0*ci) * (-ur_m_ci*drho_dr_centrifugal + ddr_rmom_centrifugal) ! incoming
         ! if ((ir .eq. 1) .and. (iz .eq. nz/2)) then
         !   print *, ' ibalanced = ', ibalanced, ' L2 = ', L2
         ! end if
      else
         ! if (iz .eq. nz/2) print *, ' L2 is outgoing'                  
         L2 =           ur_p_ci/(2.d0*ci) * (-ur_m_ci*dq_dr(iphi,iz,irho) + dq_dr(iphi,iz,rmom) ) ! outgoing
      end if

      if (ur1*nhat_r .lt. 0.d0) then
         n_in = n_in + 1
         iwhich = 3
         ! if (iz .eq. nz/2) print *, ' L3 and L4 are incoming'
         L3 = ibalanced*ur1 * (-uz     * drho_dr_centrifugal + ddr_zmom_centrifugal  ) ! incoming
         L4 = ibalanced*ur1 * (-uphi_r * drho_dr_centrifugal + ddr_uphi_r_centrifugal) ! incoming
      else
         ! if (iz .eq. nz/2) print *, ' L3 and L4 are outgoing'                           
         L3 =  ur1 * (    -uz*dq_dr(iphi,iz,irho) + dq_dr(iphi,iz,zmom)) ! outgoing
         L4 =  ur1 * (-uphi_r*dq_dr(iphi,iz,irho) + dq_dr(iphi,iz,amom)) ! outgoing
      end if

      ! if (iz .eq. nz/2) print *, ' ir = ', ir, ' iz = ', iz, ' n_in = ', n_in, ' iwhich = ', iwhich
      ! S*L: (the r factor will be put below):
      dF_non_reflective(irho) =       L1 + L2
      dF_non_reflective(zmom) = uz * (L1 + L2) + L3
      dF_non_reflective(rmom) = (ur_m_ci)*L1 + (ur_p_ci)*L2
      dF_non_reflective(amom) = uphi_r*(L1 + L2) + L4

      ! Replace in favor of the non-reflective flux derivative:
      do idof = 1, 4
         dF(iphi, iz, idof, ir) = dF_non_reflective(idof)*rgrid(ir)
      end do
   end do
end do

end subroutine r_isothermal_non_reflective_bc

!----------------------------------------------------------------------------------85

subroutine r_adiabatic_non_reflective_bc(ir_boundary, q, ur, dF)

! Modifies qdot to implement a non-reflective BC in the r direction.  We require
! that incoming waves satisfy gravity+centrifugal balance.
! See notes dated approximately November 26, 2021.
! titled "Non-reflective BC in r (adiabatic)"

use thermal_parameters
use partition_data
use dof_indices
use gravity
use grid
implicit none
integer, intent(in) :: ir_boundary
real(8), dimension(sphi:ephi, sz_r:ez_r, ndof, nr), intent(in)  :: q
real(8), dimension(sphi:ephi, sz_r:ez_r,       nr), intent(in)  :: ur
real(8), dimension(sphi:ephi, sz_r:ez_r, ndof, nr), intent(out) :: dF

! Local:
integer :: ir, iz, iphi, idof
! Temps:
real(8) :: drho_dr_centrifugal, ddr_zmom_centrifugal, ddr_rmom_centrifugal, &
           ddr_uphi_r_centrifugal, dPdr_centrifugal, d_ci_dr
real(8) :: ci, ur1, uz, uphi_r, ur_p_ci, ur_m_ci
real(8), dimension(5) :: L
integer :: nhat_r ! Direction of unit normal to the r boundary
real(8), dimension(ndof) :: dF_non_reflective

real(8) :: uphi_squared_over_r, accel
integer :: nr_in
real(8) :: coeff1, coeff2

! Primitive variables at the boundaries, uL and uR:
real(8), dimension(sphi:ephi, sz_r:ez_r, ndof,    1:4 ) :: u_L
real(8), dimension(sphi:ephi, sz_r:ez_r, ndof, nr-3:nr) :: u_R

! One sided derivatives of primitive variables:
real(8), dimension(sphi:ephi, sz_r:ez_r, ndof) :: du_dr

! Debug:
integer :: iwhich

if (ir .eq. 1) then
   normal_r_in = 1
else
   normal_r_in = -1
end if

if (ir_boundary .eq. 1) then
   ! Primitive variables and their one-sided differences:   
   do ir = 1, 4
      u_L(:,:, irho,  ir) = q(:,:, irho, ir)
      u_L(:,:, iur,   ir) = q(:,:, rmom, ir) / q(:,:, irho, ir)
      u_L(:,:, iuz,   ir) = q(:,:, zmom, ir) / q(:,:, irho, ir)
      u_L(:,:, iuphi, ir) = q(:,:, amom, ir) / q(:,:, irho, ir) / rgrid(ir)
      u_L(:,:, ip,    ir) = q(:,:, ener, ir) * gm1
   end do
   du_dr(:,:,:) = (2.d0*u_L(:,:,:,4) -  9.d0*u_L(:,:,:,3) + &
                  18.d0*u_L(:,:,:,2) - 11.d0*u_L(:,:,:,1)) / 6.d0 * &
                  Ji_r(1)
else if (ir_boundary .eq. nr) then
   do ir = nr-3, nr
      u_R(:,:, irho,  ir) = q(:,:, irho, ir)
      u_R(:,:, iur,   ir) = q(:,:, rmom, ir) / q(:,:, irho, ir)
      u_R(:,:, iuz,   ir) = q(:,:, zmom, ir) / q(:,:, irho, ir)
      u_R(:,:, iuphi, ir) = q(:,:, amom, ir) / q(:,:, irho, ir) / rgrid(ir)
      u_R(:,:, ip,    ir) = q(:,:, ener, ir) * gm1      
   end do
   du_dr(:,:,:) = (11.d0*u_R(:,:,:,nr)   - 18.d0*u_R(:,:,:,nr-1) + &
                    9.d0*u_R(:,:,:,nr-2) - 2.0d0*u_R(:,:,:,nr-3)) / 6.d0 * &
                    Ji_r(nr)
end if
   
do iphi = sphi, ephi
   do iz = sz_r, ez_r
      if (i_boundary .eq. 1) then
         rho  = uL(iphi, iz, irho,  1)
         ur   = uL(iphi, iz, iur,   1)
         uphi = uL(iphi, iz, iuphi, 1)          
         p    = uL(iphi, iz, ip,    1)
      else 
         rho  = uL(iphi, iz, irho,  nr)
         ur   = uL(iphi, iz, iur,   nr)
         uphi = uL(iphi, iz, iuphi, nr)          
         p    = uL(iphi, iz, ip,    nr)
      end if
      
      c    = sqrt(gamma * p / rho)
      coeff1 = 0.5d0 * rho / c
      coeff2 = 0.5d0 / c**2
      coeff3 = 2.d0*coeff2 ! 1/c^2
      dpdr_balanced = rho * (gr(1, iz) + uphi**2/rgrid(1))         

      ! First characteristic:
         if ((ur - c)*normal_r_in .gt. 0.d0) then
            ! Characteristic wave 1 is incoming.  Use balance for dp/dr
            L(1) = (ur - c) * (-coeff1*du_dr(iphi,iz,iur) + coeff2*dpdr_balanced)            
         else
            ! Outgoing: use one-sided differences:
            L(1) = (ur - c) * (-coeff1*du_dr(iphi,iz,iur) + coeff2*du_dr(iphi,iz,ip))
         end if

         ! Second characteristic:
         if ((ur + c)*normal_r_in .gt. 0.d0) then
            ! Characteristic wave 2 is incoming.  Use balance for dp/dr
            L(2) = (ur + c) * (coeff1*du_dr(iphi,iz,iur) + coeff2*dpdr_balanced)            
         else
            ! Outgoing: use one-sided differences:
            L(2) = (ur + c) * (coeff1*du_dr(iphi,iz,iur) + coeff2*du_dr(iphi,iz,ip))
         end if

         ! Third thru fifth characteristics:
         if (ur*normal_r_in .gt. 0.d0) then
            ! Characteristic wave 3 is incoming.  Use balance for dp/dr
            L(3) = ur * (du_dr(iphi,iz,irho) - coeff3*dpdr_balanced)
            L(4) = ur * du_dr(iphi,iz,iuz  )
            L(5) = ur * du_dr(iphi,iz,iuphi)            
         else
            ! Outgoing: use one-sided differences for dp/dr:
            L(3) = ur * (du_dr(iphi,iz,irho) - coeff3*du_dr(iphi,iz,ip))
            L(4) = 0.d0
            L(5) = 0.d0            
         end if
      end if


   call terminate_with_no_save(1)
end if

do iz = sz_r, ez_r
   do iphi = sphi, ephi
      ci  = SQRT(ci_squared_initial(ir, iz))
      ur1    = ur(iphi, iz, ir)
      uz     = q(iphi, iz, zmom, ir) / q(iphi, iz, irho, ir)
      uphi_r = q(iphi, iz, amom, ir) / q(iphi, iz, irho, ir)
      uphi_squared_over_r = uphi_r**2 / rgrid(ir)**3

      accel = gr(ir, iz) + uphi_squared_over_r

      if (ir .eq. 1 ) d_ci_dr = d_ci_dr_inner(iz)
      if (ir .eq. nr) d_ci_dr = d_ci_dr_outer(iz)
      dPdr_centrifugal       = q (iphi, iz, irho, ir) * accel
      drho_dr_centrifugal = (dPdr_centrifugal - 2.d0*q(iphi,iz,irho,ir)*ci*d_ci_dr) / ci_squared_initial(ir,iz)
      ddr_rmom_centrifugal   = ur1    * drho_dr_centrifugal ! assumes d ur/dr = 0
      ddr_zmom_centrifugal   = uz     * drho_dr_centrifugal ! assumes d uz/dr = 0
      ddr_uphi_r_centrifugal = uphi_r * drho_dr_centrifugal ! assumes d (uphi*r)/dr = 0
      
      ur_m_ci = ur1 - ci
      ur_p_ci = ur1 + ci
      n_in = 0
      iwhich = 0
      if (ur_m_ci*nhat_r .lt. 0.d0) then
         n_in = n_in + 1
         iwhich = 1
         ! if (iz .eq. nz/2) print *, ' L1 is incoming'
         L1 = ibalanced*ur_m_ci/(2.d0*ci) * (ur_p_ci*drho_dr_centrifugal - ddr_rmom_centrifugal) ! incoming
      else
         ! if (iz .eq. nz/2) print *, ' L1 is outgoing'
         L1 =           ur_m_ci/(2.d0*ci) * (ur_p_ci*dq_dr(iphi,iz,irho) - dq_dr(iphi,iz,rmom) ) ! outgoing
      end if

      if (ur_p_ci*nhat_r .lt. 0.d0) then
         n_in = n_in + 1
         iwhich = 2
         ! if (iz .eq. nz/2) print *, ' L2 is incoming'         
         L2 = ibalanced*ur_p_ci/(2.d0*ci) * (-ur_m_ci*drho_dr_centrifugal + ddr_rmom_centrifugal) ! incoming
         ! if ((ir .eq. 1) .and. (iz .eq. nz/2)) then
         !   print *, ' ibalanced = ', ibalanced, ' L2 = ', L2
         ! end if
      else
         ! if (iz .eq. nz/2) print *, ' L2 is outgoing'                  
         L2 =           ur_p_ci/(2.d0*ci) * (-ur_m_ci*dq_dr(iphi,iz,irho) + dq_dr(iphi,iz,rmom) ) ! outgoing
      end if

      if (ur1*nhat_r .lt. 0.d0) then
         n_in = n_in + 1
         iwhich = 3
         ! if (iz .eq. nz/2) print *, ' L3 and L4 are incoming'
         L3 = ibalanced*ur1 * (-uz     * drho_dr_centrifugal + ddr_zmom_centrifugal  ) ! incoming
         L4 = ibalanced*ur1 * (-uphi_r * drho_dr_centrifugal + ddr_uphi_r_centrifugal) ! incoming
      else
         ! if (iz .eq. nz/2) print *, ' L3 and L4 are outgoing'                           
         L3 =  ur1 * (    -uz*dq_dr(iphi,iz,irho) + dq_dr(iphi,iz,zmom)) ! outgoing
         L4 =  ur1 * (-uphi_r*dq_dr(iphi,iz,irho) + dq_dr(iphi,iz,amom)) ! outgoing
      end if

      ! if (iz .eq. nz/2) print *, ' ir = ', ir, ' iz = ', iz, ' n_in = ', n_in, ' iwhich = ', iwhich
      ! S*L: (the r factor will be put below):
      dF_non_reflective(irho) =       L1 + L2
      dF_non_reflective(zmom) = uz * (L1 + L2) + L3
      dF_non_reflective(rmom) = (ur_m_ci)*L1 + (ur_p_ci)*L2
      dF_non_reflective(amom) = uphi_r*(L1 + L2) + L4

      ! Replace in favor of the non-reflective flux derivative:
      do idof = 1, 4
         dF(iphi, iz, idof, ir) = dF_non_reflective(idof)*rgrid(ir)
      end do
   end do
end do

end subroutine r_adiabatic_non_reflective_bc

!----------------------------------------------------------------------------------85

end module boundary_conditions

!----------------------------------------------------------------------------------85
