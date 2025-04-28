!----------------------------------------------------------------------------------85

module parallel_pade_data

real(8), allocatable, dimension(:) :: a_coeffs_r, b_coeffs_r, c_coeffs_r
logical :: i_have_left_boundary_and_part_of_interior, &
    i_have_only_part_of_interior, i_have_right_boundary_and_part_of_interior, &
    i_have_all_of_r

end module parallel_pade_data

!----------------------------------------------------------------------------------85

subroutine parallel_pade_set_up_data_for_diff_r

use grid
use pade_coefficients
use parallel_pade_data
use partition_data

! Local:
integer :: ir

if ((sr .eq. 1) .and. (er .le. nr - 3)) then
   i_have_left_r_boundary_and_part_of_interior = .true.
else if ((sr .ge. 4) .and. (er .le. nr - 3)) then
   i_have_only_part_of_r_interior = .true.
else if ((sr .eq. nr) .and. (er .lt. nr - 2)) then
   i_have_right_r_boundary_and_part_of_interior = .true.
else if ((sr .eq. 1) .and. (er .eq. 1)) then
   i_have_all_of_r = .true.
else
   print *, ' my_node = ', my_node
   print *, ' None of the coded possibilities fit the range of r I have'
   call terminate_without_save
end if

allocate(a_coeffs_r(nr), b_coeffs_r(nr), c_coeffs_r(nr))

a_coeffs_r(1) = 0.d0 ! Hopefully never used
b_coeffs_r(1) = 1.d0
c_coeffs_r(1) = alpha1

a_coeffs_r(2) = alpha
b_coeffs_r(2) = 1.d0
c_coeffs_r(2) = alpha

a_coeffs_r(3) = alpha3
b_coeffs_r(3) = 1.d0
c_coeffs_r(3) = alpha3

a_coeffs_r(nr) = alpha1
b_coeffs_r(nr) = 1.d0
c_coeffs_r(nr) = 0.d0 ! Hopefully never used

a_coeffs_r(nr-1) = alpha
b_coeffs_r(nr-1) = 1.d0
c_coeffs_r(nr-1) = alpha

a_coeffs_r(nr-2) = alpha3
b_coeffs_r(nr-2) = 1.d0
c_coeffs_r(nr-2) = alpha3

do ir = 4, nr - 3
   a_coeffs_r(ir) = alpha3
   b_coeffs_r(ir) = 1.d0
   c_coeffs_r(ir) = alpha3
end do

end subroutine parallel_pade_set_up_coeffs_for_diff_r

!----------------------------------------------------------------------------------85

subroutine parallel_pade_diff_r(f, fp)

! You need to pass only that portion of the data that this processor has.

use partition_data_for_alan
use parallel_pade_coeffs
use grid
implicit none
real(8), dimension(nr),             intent(in)  :: Ji
real(8), dimension(nz*mphi, sr:er), intent(in)  :: f
real(8), dimension(nz*mphi, sr:er), intent(out) :: fp

! Local:
integer :: ib. nbundle
real(8), dimension(sr: er) :: rhs

nbundle = nz*mphi

! RHS of matrix system:
if (i_have_left_r_boundary_and_part_of_interior) then
   do ib = 1, nbundle
      fp(ib, 1) =  a1*f(ib, 1) + b1*f(ib, 2  ) + c1*f(ib, 3  ) + d1*f(ib, 4  )
      fp(ib, 2) = a_over_2 * (f(ib, 3) - f(ib, 1  ))
      fp(ib, 3) = a3_over_2 * (f(ib, 4  ) - f(ib, 2  ))      
   end if
   do ir = 4, er
      do ib = 1, ib
         fp(ib, ir) = a_over_2 * (f(ib, ir+1) - f(ib, ir-1))
      end do      
   end do
end if

if (i_have_only_part_of_r_interior)
   do ir = sr, er
      do ib = 1, nbundle
         fp(ib, ir) = a_over_2 * (f(ib, i+1) - f(ib, i-1))
      end do
   end do
end if

if (i_have_right_r_boundary_and_part_of_interior) then
   do ib = 1, nbundle
      fp(ib, nr  ) = -a1*f(ib, nr) - b1*f(ib, nr-1) - c1*f(ib, nr-2) - d1*f(ib, nr-3)
      fp(ib, nr-1) = a_over_2 * (f(ib, nr) - f(ib, nr-2))
      fp(ib, nr-2) = a3_over_2 * (f(ib, nr-1) - f(ib, nr-3))      
   end do
   do ir = sr, nr-3
      do ib = 1, nbundle
         fp(ib, ir) = a_over_2 * (f(ib, ir+1) - f(ib, ir-1))   
      end do
   end do
end if

if (i_have_all_of_r) then
   do ib = 1, nbundle
      fp(ib, 1) =  a1*f(ib, 1) + b1*f(ib, 2  ) + c1*f(ib, 3  ) + d1*f(ib, 4  )
      fp(ib, 2) = a_over_2 * (f(ib, 3) - f(ib, 1  ))
      fp(ib, 3) = a3_over_2 * (f(ib, 4  ) - f(ib, 2  ))      
   end if
   do ir = 4, nr - 3
      do ib = 1, ib
         fp(ib, ir) = a_over_2 * (f(ib, ir+1) - f(ib, ir-1))
      end do      
   end d
   do ib = 1, nbundle
      fp(ib, nr  ) = -a1*f(ib, nr) - b1*f(ib, nr-1) - c1*f(ib, nr-2) - d1*f(ib, nr-3)
      fp(ib, nr-1) = a_over_2 * (f(ib, nr) - f(ib, nr-2))
      fp(ib, nr-2) = a3_over_2 * (f(ib, nr-1) - f(ib, nr-3))      
   end do
end if

call PaScaL_TDMA_many_rhs_solve(mp_comm_xz, a_coeffs_r(sr), b_coeffs_r(sr), &
    c_coeffs_r(sr), fp(sr), nbundle, mr)

do ir = sr, er
   do ib = 1, nbundle
      fp(ib, ir) = Ji_r(ir) * fp(ib, ir)
   end do
end do


end subroutine parallel_pade_diff_r

!----------------------------------------------------------------------------------85
