!----------------------------------------------------------------------------------85

module pade_filter

logical :: apply_pade_filter, filter_relative_to_basic_state
! Three characters ('tau' or 'eps').
! 'eps' : use the specified eps_filter
! 'tau' : determine eps_filter from tau.
character(3) :: eps_or_tau
real(8) eps_filter, tau_filter

contains

!----------------------------------------------------------------------------------85

subroutine smooth_field (nvar, q, eps_filter_now)

! Smooths the q field in all three directions.  Note:  This routine is also used
! to smooth the beta_Delta field for artificial bulk viscosity.  This is a scalar
! field (nvar = 1, nvar_dim = 1).  Note: Ths routine assumes that nvar is also
! "nvar_dim", the nvar dimension of the array that is passed. 

use grid
use partition_data
use transposes_of_q_and_qdot
use dof_indices
use basic_state
use thermal_parameters
implicit none
integer, intent(in) :: nvar
real(8), intent(inout), dimension(sr:er, sphi:ephi, nz, nvar) :: q
real(8), intent(in) :: eps_filter_now

! Local:
integer :: iz, iphi, ir
real(8) :: Q_filter

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' subroutine smooth_field called with nvar = ', nvar
end if
#endif

! Alan's smoothing parameter:
Q_filter = 0.5d0 - 0.25d0*eps_filter_now

if (my_node .eq. 0) then
   print *, ' node 0: In smooth_field eps_filter_now = ', eps_filter_now
end if

! I added the second check since this routine is also called by the artificial_pressure
! module to smooth the bulk viscosity.
if ((filter_relative_to_basic_state) .and. (nvar .ne. 1)) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, iz, irho) = q(ir, iphi, iz, irho) - rho_basic (ir, iz)
            q(ir, iphi, iz, amom) = q(ir, iphi, iz, amom) - amom_basic(ir, iz)
         end do
      end do
   end do

   if (.not. isothermal) then
      do iz = 1, nz
         do iphi = sphi, ephi
            do ir = sr, er
               q(ir, iphi, iz, ener) = q(ir, iphi, iz, ener) - eint_basic (ir, iz)
            end do
         end do
      end do
   end if
end if

if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' smooth_field:  calling smooth_z for z dir'
#endif
   if (periodic_z) then
      call smooth_periodic_z(mr*mphi, nz, nvar, q, Q_filter)
   else
      call smooth_non_periodic(mr*mphi, nz, nvar, q, Ji_z, Q_filter)
   end if
end if

! Smooth in the r direction:
if (nr .ne. 1) then
   call transpose_z_to_r(nvar, q, q_r_space)
   
#ifdef debug_print
   if (my_node .eq. 0) print *, ' smooth_field:  calling smooth_z for r dir'
#endif   
   call smooth_non_periodic(mz_r*mphi*nvar, nr, 1, q_r_space, Ji_r, Q_filter)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' smooth_field:  calling transpose_r_to_z'
#endif   
   call transpose_r_to_z (nvar, q_r_space, q)
end if

if (nphi .ne. 1) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' smooth_field:  calling transpose_z_to_phi'
#endif   
   call transpose_z_to_phi (nvar, q, q_phi_space)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' smooth_field:  calling smooth_periodic'
#endif
   call smooth_periodic(mr*mz_phi*nvar, mr*mz_phi*nvar, nphi, q_phi_space, Q_filter)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node = 0: smooth_field:  calling transpose_phi_to_z'
#endif   
   call transpose_phi_to_z (nvar, q_phi_space, q)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node = 0: smooth_field: returned from transpose_phi_to_z'
#endif   
end if

! Add back basic state if needed:
if ((filter_relative_to_basic_state) .and. (nvar .ne. 1)) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, iz, irho) = q(ir, iphi, iz, irho) + rho_basic (ir, iz)
            q(ir, iphi, iz, amom) = q(ir, iphi, iz, amom) + amom_basic(ir, iz)
         end do
      end do
   end do

   if (.not. isothermal) then
      do iz = 1, nz
         do iphi = sphi, ephi
            do ir = sr, er
               q(ir, iphi, iz, ener) = q(ir, iphi, iz, ener) + eint_basic (ir, iz)
            end do
         end do
      end do
   end if
end if

end subroutine smooth_field

!----------------------------------------------------------------------------------85
   
! Conservative smoothing in z, where conservation means preserving sum(u(k)/rlz(k),k=1..nz).
! The away-from-the-boundaries formula is the same as that used for smoothing in the
! periodic (x,y) directions, namely
!       a*U[i-1] + U[i] + a*U[i+1] = P*(u[i-2]+u[i+2]) + Q*(u[i-1]+u[i+1]) + R*u[i]
! where U is the smoothed u.  Q is the smoothing parameter: Q <= 1/2 (with no smoothing at Q = 1/2).
! The near-boundary formulae are
!       B*U[1]  + C*U[2]    + D*U[3]    = E*u[1]  + F*u[2]    + G*u[3]    + H*u[4]
!       B*U[nz] + C*U[nz-1] + D*U[nz-2] = E*u[nz] + F*u[nz-1] + G*u[nz-2] + H*u[nz-3]
! The boundary values are not changed.
! The smoothing is 4th-order, i.e., it does not change constant, linear, quadratic, or cubic
! polynomials, except at next-to-boundary points.
!
! The coeffcients are derived using the following Mathematica expression.
! Since Mathematica won't allow E as a variable, EE is used in the Mathematica expressions; it
! corresponds to E in the comment paragraph above and in the Fortran below.
!
!eqx = {A + B == q, C + a == q, D + 1 + a == q, a + 1 + a == q,
!       A + EE + P == q, F + Q + P == q, G + R + Q + P == q, 
!       H + Q + R + Q + P == q, P + Q + R + Q + P == q, 
!       a + 1 + a == 2*P + 2*Q + R, a + 2 + 3*a == 4*P + 4*Q + 2*R, 
!       a + 4 + 9*a == 16*P + 10*Q + 4*R, a + 8 + 27*a == 64*P + 28*Q + 8*R,
!       B + C + D == EE + F + G + H, 
!       C + 2*D == F + 2*G + 3*H, A == 1}
! 
! The first four equations define the coefficients of U(k) for the sum over k; the next five
! do the same for u(k).  Since they are made equal at each k, we get q*sum(U/rlz) = q*sum(u/rlz).
! The next four make polynomials of order 0 to 3 unchanged by the smoothing, except at the
! next-to-boundary points.
! The next two equations make polynomials of order 0 to 1 unchanged at the next-to-boundary points.
! The equation 2*P - 2*Q + R == 0, which is not enforced, kills 2-delta-x at the interior
! points, but this case is very pooly conditioned.  See the comment about a_incr in
! module Smoothing_Parameters.
! q is an arbitrary (but non-zero) weighting of all the rhs; it is required to get a solution
! with the diagonal coefficients (other than B) = 1; q depends on Q.
! A is the coefficient of the unchanging boundary values and is taken to be 1, as in tridc.
! Solution in Mathematica is by: solx = Solve[eqx, {A, B, C, D, EE, F, G, H, P, R, a, q}];
! Q is left out of the list in order to get the solution in terms of Q.

Subroutine smooth_non_periodic(mxy, nz, n, u, rlz, Q)
  ! Use Smoothing_Parameters
  use partition_data
  Implicit none
  Integer, intent(in) :: mxy, nz, n
  Real(8), intent(inout) :: u(mxy,nz,n)
  Real(8), intent(in) :: rlz(nz)
  Real(8), intent(in) :: Q
  Integer :: i, k, m
  Real(8) :: u0, um1(mxy), um2(mxy)     ! temps for doing the smoothing in-place

  Real(8) :: a, B, C, D, E, F, G, H, P, R

  a = -0.5d0 + 2.d0*Q
  B = 2*a
  C = 1 + a
  D = a
  E = 0.25d0*(7*a + Q)
  F = 0.25d0*(4 + 7*a - 3*Q)
  G = 0.25d0*(a + 3*Q)
  H = 0.25d0*(a - Q)
  P = 0.25d0*(a - Q)
  R = 0.5d0*(2 + 3*a - 3*Q)

  Do m = 1,n
     Do k = 1,nz
        Do i = 1,mxy
           u(i,k,m) = u(i,k,m)/rlz(k)
        End Do
     End Do

     Do i = 1,mxy
        um2(i) = u(i,1,m)
        um1(i) = u(i,2,m)
     End Do
     Do i = 1,mxy
        u(i,2,m) = E*u(i,1,m) + F*u(i,2,m) + G*u(i,3,m) + H*u(i,4,m)
     End Do
     Do k = 3,nz-2
        Do i = 1,mxy
           u0 = u(i,k,m)
           u(i,k,m) = P*(u(i,k+2,m)+um2(i)) + Q*(u(i,k+1,m)+um1(i)) + R*u(i,k,m)
           um2(i) = um1(i)
           um1(i) = u0
        End Do
     End Do
     Do i = 1,mxy
        u(i,nz-1,m) = E*u(i,nz,m) + F*u(i,nz-1,m) + G*um1(i) + H*um2(i)
     End Do

     call tridc(mxy, nz, 1, a, B, C, D, C, B, u(1,1,m))

     Do k = 1,nz
        Do i = 1,mxy
           u(i,k,m) = rlz(k)*u(i,k,m)
        End Do
     End Do
  End Do

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node = 0: About to return from smooth_non_periodic'
#endif

End Subroutine smooth_non_periodic

!----------------------------------------------------------------------------------85

Subroutine smooth_non_periodic_clean(mxy, nz, ndof, u, rlz, Q)
  ! Use Smoothing_Parameters
  use partition_data
  Implicit none
  Integer, intent(in) :: mxy, nz, ndof
  Real(8), intent(inout) :: u(mxy, nz, ndof)
  Real(8), intent(in) :: rlz(nz)
  Real(8), intent(in) :: Q
  Integer :: i, idof, iz
  Real(8) :: u0, um1(mxy), um2(mxy)     ! temps for doing the smoothing in-place

  Real(8) :: a, B, C, D, E, F, G, H, P, R

  a = -0.5d0 + 2.d0*Q
  B = 2*a
  C = 1 + a
  D = a
  E = 0.25d0*(7*a + Q)
  F = 0.25d0*(4 + 7*a - 3*Q)
  G = 0.25d0*(a + 3*Q)
  H = 0.25d0*(a - Q)
  P = 0.25d0*(a - Q)
  R = 0.5d0*(2 + 3*a - 3*Q)

  Do idof = 1, ndof
     Do iz = 1, nz
        Do i = 1, mxy
           u(i, iz, idof) = u(i, iz, idof) / rlz(iz)
        End Do
     End Do

     Do i = 1, mxy
        um2(i) = u(i, 1, idof)
        um1(i) = u(i, 2, idof)
     End Do
     Do i = 1, mxy
        u(i, 2, idof) = E*u(i, 1, idof) + F*u(i, 2, idof) + G*u(i, 3, idof) + H*u(i, 4, idof)
     End Do
     Do iz = 3, nz - 2
        Do i = 1, mxy
           u0 = u(i, iz, idof)
           u(i, iz, idof) = P*(u(i, iz+2, idof) + um2(i)) + &
                            Q*(u(i, iz+1, idof) + um1(i)) + &
                            R* u(i, iz,   idof)
           um2(i) = um1(i)
           um1(i) = u0
        End Do
     End Do
     Do i = 1,mxy
        u(i,nz-1, idof) = E*u(i,nz, idof) + F*u(i,nz-1,idof) + G*um1(i) + H*um2(i)
     End Do

     call tridc(mxy, nz, 1, a, B, C, D, C, B, u(1,1,idof))

     Do iz = 1, nz
        Do i = 1,mxy
           u(i, iz, idof) = rlz(iz)*u(i, iz, idof)
        End Do
     End Do
  End Do

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node = 0: About to return from smooth_non_periodic'
#endif

End Subroutine smooth_non_periodic_clean

!----------------------------------------------------------------------------------85

subroutine smooth_uz_wrt_z (q, Q_smooth)

use grid
use dof_indices
implicit none
real(8), dimension(nr, nphi, nz, ndof) :: q
real(8) :: Q_smooth

! Local:
real(8) :: uz(nr, nphi, nz)
integer :: iz, iphi, ir

do iz = 1, nz
   do iphi = 1, nphi
      do ir = 1, nr
         uz(ir, iphi, iz) = q(ir, iphi, iz, zmom) / q(ir, iphi, iz, irho)
      end do
   end do
end do

call smooth_non_periodic(nr*nphi, nz, 1, uz(1, 1, 1), Ji_z, Q_smooth)

do iz = 1, nz
   do iphi = 1, nphi
      do ir = 1, nr
         q(ir, iphi, iz, zmom) = q(ir, iphi, iz, irho) * uz(ir, iphi, iz)
      end do
   end do
end do

end subroutine smooth_uz_wrt_z

!----------------------------------------------------------------------------------85

subroutine smooth_periodic(nbundle, nbundle_dim, nphi, u, Q_smooth_periodic)

! Smooth "u" w.r.t. to its second index assumed to be periodic.

  ! Use Smoothing_Parameters
  Implicit none
  Integer, intent(in)    :: nbundle, nbundle_dim, nphi
  Real(8), intent(inout) :: u(nbundle_dim, nphi)
  Real(8), intent(in) :: Q_smooth_periodic
  Integer :: iphi, m
  Real(8) :: u0, um1, um2, u_1, u_2     ! temps for doing the smoothing in-place and periodic
  Real(8) :: aa, C, B
  aa = 2.d0*Q_smooth_periodic-0.5d0; C = 0.25d0*Q_smooth_periodic-0.125d0; B = 2.d0*(Q_smooth_periodic-C)
  Do m = 1, nbundle
     u_1 = u(m, 1)
     u_2 = u(m, 2)
     um2 = u(m, nphi-1)
     um1 = u(m, nphi)
     Do iphi = 1, nphi-2
        u0 = u(m,iphi)
        u(m,iphi) = C*(u(m,iphi+2)+um2) + Q_smooth_periodic*(u(m,iphi+1)+um1) + B*u(m,iphi)
        um2 = um1
        um1 = u0
     End Do
     u0 = u(m, nphi-1)
     u(m,nphi-1) = C*(u_1+um2) + Q_smooth_periodic*(u(m,nphi)+um1) + B*u(m,nphi-1)
     um2 = um1
     um1 = u0
     u(m,nphi) = C*(u_2+um2) + Q_smooth_periodic*(u_1+um1) + B*u(m,nphi)
  End Do
  Call ptrid2nd(nbundle, nbundle_dim, nphi, aa, 1.d0, aa, u)
End Subroutine smooth_periodic

!----------------------------------------------------------------------------------85

subroutine smooth_periodic_z(nbundle, nz, ndof, u, Q)

! Smooth u(nbundle, nz, ndof) w.r.t. to its second index assumed to be periodic.

implicit none
integer, intent(in)    :: nbundle, nz, ndof
real(8), intent(inout) :: u(nbundle, nz, ndof)
Real(8), intent(in) :: Q

! Local:
integer :: idof, m, iz
real(8) :: u0, um1, um2, u_1, u_2     ! temps for doing the smoothing in-place and periodic
real(8) :: aa, C, B
aa = 2.d0*Q - 0.5d0; C = 0.25d0*Q - 0.125d0; B = 2.d0*(Q - C)

do idof = 1, ndof
   do m = 1, nbundle
      u_1 = u(m, 1, idof)
      u_2 = u(m, 2, idof)
      um2 = u(m, nz-1, idof)
      um1 = u(m, nz,   idof)
      do iz = 1, nz-2
         u0 = u(m, iz, idof)
         u(m, iz, idof) = C*(u(m, iz+2, idof) + um2) + Q*(u(m, iz+1, idof) + um1) + B*u(m, iz, idof)
         um2 = um1
         um1 = u0
      end do
      u0 = u(m, nz-1, idof)
      u(m, nz-1, idof) = C*(u_1 + um2) + Q*(u(m, nz, idof) + um1) + B*u(m, nz-1, idof)
      um2 = um1
      um1 = u0
      u(m, nz, idof) = C*(u_2 + um2) + Q*(u_1 + um1) + B*u(m, nz, idof)
   end do
   call ptrid2nd(nbundle, nbundle, nz, aa, 1.d0, aa, u(1,1,idof))
end do
end subroutine smooth_periodic_z

!----------------------------------------------------------------------------------85

end module pade_filter

!----------------------------------------------------------------------------------85
