!----------------------------------------------------------------------------------85

module pade_coefficients

! All the coefficients below assume that the grid-size h = 1.
! The grid size (and variable grid size) is treated via the grid inverse Jacobians
! Ji_z and Ji_r.

! Coefficients of interior scheme.  Currently standard fourth-order Pade.
real(8), parameter :: alpha = 0.25d0, a_over_2 = 0.75d0

! Coefficients of first row of boundary scheme (4.1.4) in Lele's paper.
! This is a third-order scheme.
real(8), parameter :: alpha1 = 2.d0, a1 = -15.d0/6.d0, b1 = 2.0d0, c1 = 0.5d0, &
     d1 = 0.0d0

! Coefficients of the second row.  Currently standard fourth-order Pade.
real(8), parameter :: alpha2 = 0.25d0, a2_over_2 = 0.75d0

! Coefficients of the third row will be given in subroutine set_up_pade_coefficients
real(8) :: alpha3, a3_over_2, b3_over_4

! Quadrature weights for conservation.  See notes paginated CPS.  These weights are
! assigned in subroutine set_up_pade_coefficients below.
real(8) :: w1, w2, w3 ! Provisional Lele weights.
real(8) :: w_final(4), w_interior_final

! Quadrature weights for each grid point:
real(8), allocatable, dimension(:) :: conservative_pade_weight_r, &
                                      conservative_pade_weight_z

end module pade_coefficients

!----------------------------------------------------------------------------------85

subroutine set_up_pade_coefficients

! This routine is called by domain_mesh_and_partition.f90 
! These coefficients are for conservation.  See Sec. 4.2 in Lele and your notes
! paginated CPS.

!use partition_data
use pade_coefficients
use grid
use logical_units
use partition_data
implicit none

! Local:
real(8) :: q, r, s, alpha_hat, numer, denom, a3, b3, sum
integer :: ir, iz, i

! Go over to Lele's notation. b1, c1, and d1 are the rhs coefficients of the
! first row.
q = b1
r = c1
s = d1

alpha_hat = alpha

! Scheme for third row based on conservation (collapsing sum derivative):
! See pg. L-3 ofyour notes for alpha3 (which is denoted alpha'' in the notes.
numer = 7.d0 *(4.d0*alpha_hat - 1.d0)*s +      (40.d0*alpha_hat - 1.d0)*q
denom = 16.d0*(alpha_hat      + 2.d0)*q - 8.d0*(4.d0*alpha_hat  - 1.d0)*s
alpha3 = numer/denom
a3     = 2.d0/3.d0*(alpha3 + 2.d0)
b3     = 1.d0/3.d0*(4.d0*alpha3 - 1.d0)
a3_over_2 = 0.50d0*a3
b3_over_4 = 0.25d0*b3

! "Quadrature" weights:
w1 = (2.d0*alpha_hat + 1.d0) / (2.d0*(q + s))
 
numer = (8.d0*alpha_hat + 7.d0)*s - 6.d0*(2.d0*alpha_hat + 1.d0)*r + (8.d0*alpha_hat + 7.d0)*q
denom = 9.d0 * (q + s)
w2 = numer/denom

numer = 4.d0*(alpha_hat + 2.d0)*q - 2.d0*(4.d0*alpha_hat - 1.d0)*s
w3 = numer/denom ! same denom as for w2

if (my_node .eq. 0) then
   print *, ' Checking satisfaction of equations:'
   print *, ' eq1 = ', w1*b1 - w3*a3_over_2
   print *, ' eq2 = ', w1*c1 + (w2-1.d0)*a_over_2
   print *, ' eq3 = ', w3*a3_over_2 - a_over_2
   print *, ' Lele weights (provisional weights are):'
   print *, ' w1 = ', w1
   print *, ' w2 = ', w2
   print *, ' w3 = ', w3

   print *, ' coefficients of third row'
   print *, ' alpha3    = ', alpha3
   print *, ' a3_over_2 = ', a3_over_2
   print *, ' b3_over_4 = ', b3_over_4
end if

! Get final weights:
denom = w1*a1 - w2*a2_over_2 - w3*b3_over_4
w_final(1) = -(w1 + w2*alpha2) / denom
w_final(2) = -(w1*alpha1 + w2 + w3*alpha3) / denom
w_final(3) = -(w2*alpha2 + w3 + alpha) / denom
w_final(4) = -(w3*alpha3 + 1.d0 + alpha) / denom
! This should be unity:
w_interior_final = -(alpha + 1.d0 + alpha) / denom

if (my_node .eq. 0) then
   write(6, "('Quadrature weights for conservative Pade')")
   print *, ' denominator for final weights = ', denom 
   sum = 0.d0
   do i = 1, 4
      write(6, 1) i, w_final(i)
      1 format(' i = ', i1, ' w_final = ', e12.5)
      sum = sum + w_final(i)
   end do
   write(6, 2) w_interior_final
   2  format(' w_interior_final = ', e12.5)
   write(6, 3) sum
   3 format(' sum of w1 thru w4 = ', e12.5)            
end if

! These are weight arrays to be used for integrations over the domain to
! check for conservation or to compute diagnostics.
if (nr .ne. 1) then
   allocate(conservative_pade_weight_r(nr))
   do i = 1, 4
      conservative_pade_weight_r(i         ) = w_final(i) / Ji_r(i)
      conservative_pade_weight_r(nr - i + 1) = w_final(i) / Ji_r(nr-i+1)
   end do

   do ir = 5, nr - 4
      conservative_pade_weight_r(ir) = w_interior_final / Ji_r(ir)
   end do
end if

if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1))then
   allocate(conservative_pade_weight_z(nz))
   do i = 1, 4
      conservative_pade_weight_z(i         ) = w_final(i) / Ji_z(i)
      conservative_pade_weight_z(nz - i + 1) = w_final(i) / Ji_z(nz-i+1)      
   end do

   do iz = 5, nz - 4
      conservative_pade_weight_z(iz) = w_interior_final / Ji_z(iz)
   end do
end if

! Check that b3_over_4 = 0 as assumed below:
if (b3_over_4 .ne. 0.d0) then
   print *, ' b3_over_4 = ', b3_over_4, ' is not 0'
   call terminate_with_no_save(1)
end if

if (my_node .eq. 0) print *, ' Returning from set_up_pade_coefficients'

#ifdef debug_print
   if (my_node .eq. 0) print *, ' Returning from set_up_pade_coefficients'
#endif

end subroutine set_up_pade_coefficients

!----------------------------------------------------------------------------------85

subroutine pade_diff_z(nbundle, F, Fd)

! Differentiator for the z-direction. Currently very much USED.

use grid
implicit none
integer :: nbundle
real(8), dimension(nbundle, nz) :: F, Fd

if (periodic_z) then
   call pade_diff_periodic(nbundle, nz, dz_periodic, F, Fd)
else
   call pade_diff_bundle  (nbundle, nz, Ji_z,        F, Fd)
end if

end subroutine pade_diff_z

!----------------------------------------------------------------------------------85

subroutine pade_diff_bundle(nbundle, n, Ji, f, fp)

! USED. An the pade scheme is conservative.

! Standard 4th order Pade differentiation with 4th order numerical boundary scheme.
! Ji is the inverse Jacobian: Ji = dxi/dx, where xi is the numerical index.
! Oct. 24, 2017.

use pade_coefficients
implicit none
integer,                        intent(in)  :: nbundle, n
real(8), dimension(n),          intent(in)  :: Ji
real(8), dimension(nbundle, n), intent(in)  :: f
real(8), dimension(nbundle, n), intent(out) :: fp

integer :: i, ib

! fp is the RHS vector and will eventually become the the derivative.
! First and last rows of the matrix.
do ib = 1, nbundle
   fp(ib, 1) =  a1*f(ib, 1) + b1*f(ib, 2  ) + c1*f(ib, 3  ) + d1*f(ib, 4  )
   fp(ib, n) = -a1*f(ib, n) - b1*f(ib, n-1) - c1*f(ib, n-2) - d1*f(ib, n-3)
end do

! Second and second-to-last rows.  Same as interior.
do ib = 1, nbundle
   fp(ib, 2)   = a_over_2 * (f(ib, 3) - f(ib, 1  ))
   fp(ib, n-1) = a_over_2 * (f(ib, n) - f(ib, n-2))
end do

! Third and third-to-last rows.  This assumes that b3_over_4 = 0 which
! is the case for us.
do ib = 1, nbundle
   fp(ib, 3)   = a3_over_2 * (f(ib, 4  ) - f(ib, 2  ))
   fp(ib, n-2) = a3_over_2 * (f(ib, n-1) - f(ib, n-3))
end do

do i = 4, n - 3
   do ib = 1, nbundle
      fp(ib, i) = a_over_2 * (f(ib, i+1) - f(ib, i-1))
   end do
end do

! Call is:
! tridiag_bundle (m, a, b1, c1, bm, am, x)
! Don't confuse with notation in Lele's paper.
! First row is b1, c1 = 1, alpha1
! Last row is  am, bm = alpha1 1

call tridiag_ddx_conservative(nbundle, n, alpha1, alpha, alpha3, fp)

do i = 1, n
   do ib = 1, nbundle
      fp(ib, i) = Ji(i) * fp(ib, i)
   end do
end do

return
end subroutine pade_diff_bundle

!----------------------------------------------------------------------------------85

subroutine pade_diff_periodic(nbundle, n, h, f, fp)

! Standard 4th order Pade differentiation with periodic BC.
! Nov. 17, 2017.

use partition_data
implicit none
integer,                        intent(in)  :: nbundle, n
real(8),                        intent(in)  :: h
real(8), dimension(nbundle, n), intent(in)  :: f
real(8), dimension(nbundle, n), intent(out) :: fp

! Coefficients for standard Pade scheme:
real(8), parameter :: alpha = 0.25d0, a_over_2 = 0.75d0

! Locals:
integer :: i, ib
real(8) :: a_over_2h

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: pade_diff_periodic has been called with nbundle = ', nbundle
#endif

a_over_2h = a_over_2 / h

! This is the RHS which will eventually become the derivative.
do ib = 1, nbundle
   fp(ib, 1) =  a_over_2h * (f(ib, 2) - f(ib, n  ))
   fp(ib, n) =  a_over_2h * (f(ib, 1) - f(ib, n-1))
end do

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: pade_diff_periodic, finished with first loop'
#endif  

do i = 2, n - 1
   do ib = 1, nbundle
      fp(ib, i) = a_over_2h * (f(ib, i+1) - f(ib, i-1))
   end do
end do

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: pade_diff_periodic, finished with RHS calculation'
#endif  

call ptrid2nd(nbundle, nbundle, n, alpha, 1.0d0, alpha, fp)

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: returning from pade_diff_periodic'
#endif

return
end subroutine pade_diff_periodic

!----------------------------------------------------------------------------------85

