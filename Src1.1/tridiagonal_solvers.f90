!----------------------------------------------------------------------------------85

subroutine tridiag_precomputed_cc (m, cc, a, b1, c1, bm, am, x)

! Tridiagonal solver.  Courtesy of Alan Wray's routine triddx.
! Matrix entries are scalars.  Only first and last rows are allowed to be different.

implicit none
integer, intent(in)    :: m ! # of rows
real(8), intent(in)    :: cc(m)
real(8), intent(in)    :: a, b1, c1, bm, am
real(8), intent(inout) ::  x(m)

! Local:
integer :: j
!real(8) :: cc(m)  Now precomputed

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    (b1 c1 0 . . . . . . . ) x(1)      r(1)
  !         
  !    (a  1  a  0 . . . . .  ) x(2)      r(2)
  !           
  !    (0  a  1  a  0 . . . . ) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !    (. . . .  0   a   1  a ) x(m-1)    r(m-1)
  !                 
  !    (. . . . . .  0   am bm) x(m)      r(m)

! This is precomputed:
!cc(1) = c1/b1
!do j = 2, m-1
!   cc(j) = a/(1.d0 - a*cc(j-1))
!end do

x(1) = x(1)/b1

do j = 2, m-1
   x(j) = (x(j) - a*x(j-1))/(1.d0 - a*cc(j-1))
end do

x(m) = (x(m) - am*x(m-1))/(bm - am*cc(m-1))

Do j = m-1, 1, -1
   x(j) = x(j) - cc(j)*x(j+1)
end Do

end subroutine tridiag_precomputed_cc

!----------------------------------------------------------------------------------85

subroutine tridiag_internal_cc (m, a, b1, c1, bm, am, x)

! Tridiagonal solver.  Courtesy of Alan Wray's routine triddx.
! Matrix entries are scalars: only first and last rows are allowed to be different.

implicit none
integer, intent(in)    :: m ! # of rows
real(8), intent(in)    :: a, b1, c1, bm, am
real(8), intent(inout) ::  x(m)

! Local:
integer :: j
real(8) :: cc(m)

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    (b1 c1 0 . . . . . . . ) x(1)      r(1)
  !         
  !    (a  1  a  0 . . . . .  ) x(2)      r(2)
  !           
  !    (0  a  1  a  0 . . . . ) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !    (. . . .  0   a   1  a ) x(m-1)    r(m-1)
  !                 
  !    (. . . . . .  0   am bm) x(m)      r(m)

cc(1) = c1/b1
do j = 2, m-1
   cc(j) = a/(1.d0 - a*cc(j-1))
end do

x(1) = x(1)/b1

do j = 2, m-1
   x(j) = (x(j) - a*x(j-1))/(1.d0 - a*cc(j-1))
end do

x(m) = (x(m) - am*x(m-1))/(bm - am*cc(m-1))

Do j = m-1, 1, -1
   x(j) = x(j) - cc(j)*x(j+1)
end Do

end subroutine tridiag_internal_cc

!----------------------------------------------------------------------------------85

subroutine tridiag_ddx_conservative(nbundle, m, a1, a, a3, x)

! Tridiagonal solver for the matrix that arises for conservative differentiation.
! See notes of Oct. 27, 2020.

! Note: a1 coefficient in first and last row
!       a  coefficient in second row, second to last row and internal rows
!       a3 coefficient in third row and third to the last row.
!   Solves the following special tridiagonal system; r is passed in x
!
!    (1 a1 0 . . . . . . .        ) x(1)  =   r(1)         
!    (a 1  a  0 . . . . .         ) x(2)  =   r(2)           
!    (0 a3 1  a3  0 . . . .       ) x(3)  =   r(3)
!    (0  0  a  1   a  0 . . .     ) x(4)  =   r(4)
!    (0  0  0  a   1  a   . .     ) x(5)  =   r(5)
!            

!    (                   a3 1 a3 0) x(m-2) =  r(m-2)
!    (. . . .            0  a 1  a) x(m-1) =  r(m-1)
!    (. . . . . . .         0 a1 1) x(m)   =  r(m)

implicit none
integer, intent(in)    :: nbundle, m ! bundle size, # of rows
real(8), intent(in)    :: a1, a, a3
real(8), intent(inout) ::  x(nbundle, m)

! Local:
integer :: j, ib
real(8) :: cc(m)

! cc can be precomputed but it is a small extra computation.
! Eliminate the lower diagonal keeping the
! main diagonal as unity:
cc(1) = a1

cc(2) = a /(1.d0 - a*cc(1))
do ib = 1, nbundle
   x(ib, 2) = (x(ib, 2) - a*x(ib, 1)) / (1.d0 - a*cc(1))
end do

cc(3) = a3/(1.d0 - a3*cc(2))
do ib = 1, nbundle
   x(ib, 3) = (x(ib, 3) - a3*x(ib, 2)) / (1.d0 - a3*cc(2))
end do

do j = 4, m-3
   cc(j) = a/(1.d0 - a*cc(j-1))
end do
do j = 4, m-3
   do ib = 1, nbundle
      x(ib, j) = (x(ib, j) - a*x(ib, j-1))/(1.d0 - a*cc(j-1))
   end do
end do

cc(m-2) = a3/(1.d0 - a3*cc(m-3))
do ib = 1, nbundle
   x(ib, m-2) = (x(ib, m-2) - a3*x(ib, m-3))/(1.d0 - a3*cc(m-3))
end do

cc(m-1) = a/(1.d0 - a*cc(m-2))
do ib = 1, nbundle
   x(ib, m-1) = (x(ib, m-1) - a*x(ib, m-2))/(1.d0 - a*cc(m-2))
end do

! Solve last row:
do ib = 1, nbundle
   x(ib, m) = (x(ib, m) - a1*x(ib, m-1))/(1.d0 - a1*cc(m-1))
end do

! Back solve the rest:
do j = m-1, 1, -1
   do ib = 1, nbundle
      x(ib, j) = x(ib, j) - cc(j)*x(ib, j+1)
   end do
end do

end subroutine tridiag_ddx_conservative

!----------------------------------------------------------------------------------85

subroutine tridiag_bundle(nbundle, m, a, b1, c1, bm, am, x)

! Tridiagonal solver.  Courtesy of Alan Wray's routine triddx.
! Matrix entries are scalars.  Only first and last rows are allowed to be different.

implicit none
integer, intent(in)    :: nbundle, m ! bundle size, # of rows
real(8), intent(in)    :: a, b1, c1, bm, am
real(8), intent(inout) ::  x(nbundle, m)

! Local:
integer :: j, ib
real(8) :: cc(m)

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    (b1 c1 0 . . . . . . . ) x(1)      r(1)
  !         
  !    (a  1  a  0 . . . . .  ) x(2)      r(2)
  !           
  !    (0  a  1  a  0 . . . . ) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !    (. . . .  0   a   1  a ) x(m-1)    r(m-1)
  !                 
  !    (. . . . . .  0   am bm) x(m)      r(m)

! This can be precomputed.  Eliminate th lower diagonal keeping the
! main diagonal as unity:
cc(1) = c1/b1
do j = 2, m-1
   cc(j) = a/(1.d0 - a*cc(j-1))
end do

do ib = 1, nbundle
   x(ib, 1) = x(ib, 1)/b1
end do

do j = 2, m-1
   do ib = 1, nbundle
      x(ib, j) = (x(ib, j) - a*x(ib, j-1))/(1.d0 - a*cc(j-1))
   end do
end do

do ib = 1, nbundle
   x(ib, m) = (x(ib, m) - am*x(ib, m-1))/(bm - am*cc(m-1))
end do

! Back solve:
do j = m-1, 1, -1
   do ib = 1, nbundle
      x(ib, j) = x(ib, j) - cc(j)*x(ib, j+1)
   end do
end Do

end subroutine tridiag_bundle

!----------------------------------------------------------------------------------85

subroutine pre_compute_cc (m, a, b1, c1, cc)

! Precomputes cc for the tridiagonal solver.

implicit none
integer, intent(in) :: m
real(8), intent(in) :: a, b1, c1
real(8), intent(out) :: cc(m)

integer :: j

cc(1) = c1/b1
do j = 2, m-1
   cc(j) = a/(1.d0 - a*cc(j-1))
end do

return
end subroutine pre_compute_cc

!----------------------------------------------------------------------------------85

Subroutine trid_general_bundle(m, n, a, b, c, x)
  Implicit none
  Integer, intent(in) :: m, n
  Real(8), intent(in), dimension(m,n) :: a, b
  Real(8), intent(inout), dimension(m,n) :: c, x

  Integer :: i, j

  !   Solves the following general tridiagonal system; r is passed in x
  !
  !    (b  c  0 . . . . . .) x(1)      r(1)
  !         
  !    (a  b  c  0 . . . . ) x(2)      r(2)
  !           
  !    (0  a  b  c  0 . . .) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !     . . . . . . . . . .
  !
  !    (. . . .  0  a  b  c) x(m-1)    r(m-1)
  !                 
  !    (. . . . . . 0  a  b) x(m)      r(m)

  Do i=1,m
     c(i,1) = c(i,1)/b(i,1)
     x(i,1) = x(i,1)/b(i,1)
  End Do
  Do j = 2,n-1
     Do i=1,m
        c(i,j) = c(i,j)/(b(i,j) - a(i,j)*c(i,j-1))
        x(i,j) = (x(i,j)-a(i,j)*x(i,j-1))/(b(i,j) - a(i,j)*c(i,j-1))
     End Do
  End Do
  Do i=1,m
     x(i,n) = (x(i,n)-a(i,n)*x(i,n-1))/(b(i,n)-a(i,n)*c(i,n-1))
  End Do

  Do j = n-1,1,-1
     Do i=1,m
        x(i,j) = x(i,j) - c(i,j)*x(i,j+1)
     End Do
  End Do
  End Subroutine trid_general_bundle

!----------------------------------------------------------------------------------85

Subroutine trid_unit_diag_bundle(m, n, a, c, x)
  Implicit none
  Integer, intent(in) :: m, n
  Real(8), intent(in), dimension(m,n) :: a
  Real(8), intent(inout), dimension(m,n) :: c, x

  Integer :: i, j

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    (1  c  0 . . . . . .) x(1)      r(1)
  !         
  !    (a  1  c  0 . . . . ) x(2)      r(2)
  !           
  !    (0  a  1  c  0 . . .) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !     . . . . . . . . . .
  !
  !    (. . . .  0  a  1  c) x(m-1)    r(m-1)
  !                 
  !    (. . . . . . 0  a  1) x(m)      r(m)

  Do j = 2,n-1
     Do i=1,m
        c(i,j) = c(i,j)/(1.d0 - a(i,j)*c(i,j-1))
        x(i,j) = (x(i,j)-a(i,j)*x(i,j-1))/(1.d0 - a(i,j)*c(i,j-1))
     End Do
  End Do
  Do i=1,m
     x(i,n) = (x(i,n)-a(i,n)*x(i,n-1))/(1.d0-a(i,n)*c(i,n-1))
  End Do

  Do j = n-1,1,-1
     Do i=1,m
        x(i,j) = x(i,j) - c(i,j)*x(i,j+1)
     End Do
  End Do
  End Subroutine trid_unit_diag_bundle

!----------------------------------------------------------------------------------85

Subroutine trid_unit_diag_middle_index(m, n, nn, a, c, x)
  Implicit none
  Integer, intent(in) :: m, n, nn
  Real(8), intent(in), dimension(m,n,nn) :: a
  Real(8), intent(inout), dimension(m,n,nn) :: c, x

  Integer :: i, j, k

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    (1  c  0 . . . . . .) x(1)      r(1)
  !         
  !    (a  1  c  0 . . . . ) x(2)      r(2)
  !           
  !    (0  a  1  c  0 . . .) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !     . . . . . . . . . .
  !
  !    (. . . .  0  a  1  c) x(m-1)    r(m-1)
  !                 
  !    (. . . . . . 0  a  1) x(m)      r(m)

  Do k = 1,nn
     Do j = 2,n-1
        Do i = 1,m
           c(i,j,k) = c(i,j,k)/(1.d0 - a(i,j,k)*c(i,j-1,k))
           x(i,j,k) = (x(i,j,k)-a(i,j,k)*x(i,j-1,k))/(1.d0 - a(i,j,k)*c(i,j-1,k))
        End Do
     End Do
     Do i = 1,m
        x(i,n,k) = (x(i,n,k)-a(i,n,k)*x(i,n-1,k))/(1.d0-a(i,n,k)*c(i,n-1,k))
     End Do

     Do j = n-1,1,-1
        Do i = 1,m
           x(i,j,k) = x(i,j,k) - c(i,j,k)*x(i,j+1,k)
        End Do
     End Do
  End Do
  End Subroutine trid_unit_diag_middle_index

!----------------------------------------------------------------------------------85

Subroutine trid3sx(m, n, nn, nv, a, c, x)
  Implicit none
  Integer, intent(in) :: m, n, nn, nv
  Real(8), intent(in), dimension(m,nv,nn) :: a
  Real(8), intent(inout), dimension(m,nv,nn) :: c, x

  Integer :: i, j, k

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    (1  c  0 . . . . . .) x(1)      r(1)
  !         
  !    (a  1  c  0 . . . . ) x(2)      r(2)
  !           
  !    (0  a  1  c  0 . . .) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !     . . . . . . . . . .
  !
  !    (. . . .  0  a  1  c) x(m-1)    r(m-1)
  !                 
  !    (. . . . . . 0  a  1) x(m)      r(m)

  Do k = 1,nn
     Do j = 2,n-1
        Do i = 1,m
           c(i,j,k) = c(i,j,k)/(1.d0 - a(i,j,k)*c(i,j-1,k))
           x(i,j,k) = (x(i,j,k)-a(i,j,k)*x(i,j-1,k))/(1.d0 - a(i,j,k)*c(i,j-1,k))
        End Do
     End Do
     Do i = 1,m
        x(i,n,k) = (x(i,n,k)-a(i,n,k)*x(i,n-1,k))/(1.d0-a(i,n,k)*c(i,n-1,k))
     End Do

     Do j = n-1,1,-1
        Do i = 1,m
           x(i,j,k) = x(i,j,k) - c(i,j,k)*x(i,j+1,k)
        End Do
     End Do
  End Do
End Subroutine trid3sx

!----------------------------------------------------------------------------------85

Subroutine tridd(l, m, n, a, b, c, x)
  Implicit none
  Integer, intent(in) :: l, m, n
  Real(8), intent(in) :: a, b, c
  Real(8), intent(inout) ::  x(l,m,n)

  Integer :: i,j,k
  Real(8) :: cc(m)

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    (b  c  0 . . . . . . .) x(1)      r(1)
  !         
  !    (a  1  a  0 . . . . . ) x(2)      r(2)
  !           
  !    (0  a  1  a  0 . . . .) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !    (. . . . .  0  a  1  a) x(m-1)    r(m-1)
  !                 
  !    (. . . . . . . 0  c  b) x(m)      r(m)

  cc(1) = c/b
  Do j = 2,m-1
     cc(j) = a/(1.d0 - a*cc(j-1))
  End Do
  Do k=1,n
     Do i = 1,l
        x(i,1,k) = x(i,1,k)/b
     End Do
  End Do

  Do j = 2,m-1
     Do k=1,n
        Do i = 1,l
           x(i,j,k) = (x(i,j,k) - a*x(i,j-1,k))/(1.d0 - a*cc(j-1))
        End Do
     End Do
  End Do
  Do k=1,n
     Do i = 1,l
        x(i,m,k) = (x(i,m,k) - c*x(i,m-1,k))/(b - c*cc(m-1))
     End Do
  End Do
  Do j = m-1,1,-1
     Do k=1,n
        Do i = 1,l
           x(i,j,k) = x(i,j,k) - cc(j)*x(i,j+1,k)
        End Do
     End Do
  End Do
End Subroutine tridd

!----------------------------------------------------------------------------------85

Subroutine triddv(l, m, n, a, b, c, x)
  Implicit none
  Integer, intent(in) :: l, m, n
  Real(8), intent(in) :: a, b, c
  Real(8), intent(inout) ::  x(l,0:m,n)

  Integer :: i,j,k
  Real(8) :: cc(m-2)

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    x(0): no change
  !
  !    (b  c  0 . . . . . . .) x(1)      r(1)
  !         
  !    (a  1  a  0 . . . . . ) x(2)      r(2)
  !           
  !    (0  a  1  a  0 . . . .) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !    (. . . . .  0  a  1  a) x(m-1)    r(m-1)
  !                 
  !    x(m): no change

  cc(1) = c/b
  Do j = 2,m-2
     cc(j) = a/(1.d0 - a*cc(j-1))
  End Do
  Do k = 1,n
     Do i = 1,l
        x(i,1,k) = x(i,1,k)/b
     End Do
  End Do

  Do j = 2,m-2
     Do k = 1,n
        Do i = 1,l
           x(i,j,k) = (x(i,j,k) - a*x(i,j-1,k))/(1.d0 - a*cc(j-1))
        End Do
     End Do
  End Do
  Do k = 1,n
     Do i = 1,l
        x(i,m-1,k) = (x(i,m-1,k) - c*x(i,m-2,k))/(b - c*cc(m-2))
     End Do
  End Do
  Do j = m-2,1,-1
     Do k = 1,n
        Do i = 1,l
           x(i,j,k) = x(i,j,k) - cc(j)*x(i,j+1,k)
        End Do
     End Do
  End Do
End Subroutine triddv

!----------------------------------------------------------------------------------85

Subroutine triddx(l, m, n, a, b1, c1, bm, am, x)
  Implicit none
  Integer, intent(in) :: l, m, n
  Real(8), intent(in) :: a, b1, c1, bm, am
  Real(8), intent(inout) ::  x(l,m,n)

  Integer :: i,j,k
  Real(8) :: cc(m)

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    (b1 c1 0 . . . . . . . ) x(1)      r(1)
  !         
  !    (a  1  a  0 . . . . .  ) x(2)      r(2)
  !           
  !    (0  a  1  a  0 . . . . ) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !    (. . . .  0   a   1  a ) x(m-1)    r(m-1)
  !                 
  !    (. . . . . .  0   am bm) x(m)      r(m)

  cc(1) = c1/b1
  Do j = 2,m-1
     cc(j) = a/(1.d0 - a*cc(j-1))
  End Do
  Do k=1,n
     Do i = 1,l
        x(i,1,k) = x(i,1,k)/b1
     End Do
  End Do

  Do j = 2,m-1
     Do k=1,n
        Do i = 1,l
           x(i,j,k) = (x(i,j,k) - a*x(i,j-1,k))/(1.d0 - a*cc(j-1))
        End Do
     End Do
  End Do
  Do k=1,n
     Do i = 1,l
        x(i,m,k) = (x(i,m,k) - am*x(i,m-1,k))/(bm - am*cc(m-1))
     End Do
  End Do
  Do j = m-1,1,-1
     Do k=1,n
        Do i = 1,l
           x(i,j,k) = x(i,j,k) - cc(j)*x(i,j+1,k)
        End Do
     End Do
  End Do
End Subroutine triddx

!----------------------------------------------------------------------------------85

Subroutine tridc(l, m, n, a, Bb, Cb, D, Ct, Bt, x)
  Implicit none
  Integer, intent(in) :: l, m, n
  Real(8), intent(in) :: a, Bb, Cb, D, Ct, Bt
  Real(8), intent(inout) ::  x(l,m,n)

  Integer :: i,j,k
  Real(8) :: cc(m)

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    (1  0  . . . . . . . .) x(1)      r(1)
  !         
  !    (Bb Cb D  0 . . . . . ) x(2)      r(2)
  !           
  !    (0  a  1  a  0 . . . .) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !    (. . . . 0  a  1  a  0) x(m-2)    r(m-2)
  !                 
  !    (. . . . .  0  D Ct Bt) x(m-1)    r(m-1)
  !                 
  !    (. . . . . .  0   0  1) x(m)      r(m)

  cc(1) = 0
  cc(2) = D/(Cb - Bb*cc(1))
  Do j = 3,m-2
     cc(j) = a/(1.d0 - a*cc(j-1))
  End Do
  cc(m-1) = Bt/(Ct - D*cc(m-2))

  Do k=1,n
     Do i = 1,l
        x(i,2,k) = (x(i,2,k) - Bb*x(i,1,k))/(Cb - Bb*cc(1))
     End Do
  End Do
  Do j = 3,m-2
     Do k=1,n
        Do i = 1,l
           x(i,j,k) = (x(i,j,k) - a*x(i,j-1,k))/(1.d0 - a*cc(j-1))
        End Do
     End Do
  End Do
  Do k=1,n
     Do i = 1,l
        x(i,m-1,k) = (x(i,m-1,k) - D*x(i,m-2,k))/(Ct - D*cc(m-2))
     End Do
  End Do
  Do j = m-1,2,-1
     Do k=1,n
        Do i = 1,l
           x(i,j,k) = x(i,j,k) - cc(j)*x(i,j+1,k)
        End Do
     End Do
  End Do
End Subroutine tridc

!----------------------------------------------------------------------------------85

Subroutine tridc2(l, m, n, a, b, p, x)
  Implicit none
  Integer, intent(in) :: l, m, n
  Real(8), intent(in) :: a(2), b(2), p
  Real(8), intent(inout) ::  x(l,m,n)

  Integer :: i,j,k
  Real(8) :: cc(m)

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    (1  b1   . . . . . . . . . ) x(1)      r(1)
  !         
  !    (a1  1   p  0  . . . . . . ) x(2)      r(2)
  !           
  !    (0   p   1  p  0 . . . . . ) x(3)      r(3)
  !            
  !     . . . . . . . . . . . . .        =
  !
  !    (. . . . . .  0  p  1  p  0) x(m-2)    r(m-2)
  !                 
  !    (. . . . . . .  0   p  1 a2) x(m-1)    r(m-1)
  !                 
  !    (. . . . . . . . .  0 b2  1) x(m)      r(m)

  cc(1) = b(1)
  cc(2) = p/(1.d0 - a(1)*cc(1))
  Do j = 3,m-2
     cc(j) = p/(1.d0 - p*cc(j-1))
  End Do
  cc(m-1) = a(2)/(1.d0 - p*cc(m-2))

  Do k=1,n
     Do i = 1,l
        x(i,2,k) = (x(i,2,k) - a(1)*x(i,1,k))/(1.d0 - a(1)*cc(1))
     End Do
  End Do
  Do j = 3,m-2
     Do k=1,n
        Do i = 1,l
           x(i,j,k) = (x(i,j,k) - p*x(i,j-1,k))/(1.d0 - p*cc(j-1))
        End Do
     End Do
  End Do
  Do k=1,n
     Do i = 1,l
        x(i,m-1,k) = (x(i,m-1,k) -    p*x(i,m-2,k))/(1.d0 -    p*cc(m-2))
        x(i,m,k)   = (x(i,m,k)   - b(2)*x(i,m-1,k))/(1.d0 - b(2)*cc(m-1))
     End Do
  End Do
  Do j = m-1,1,-1
     Do k=1,n
        Do i = 1,l
           x(i,j,k) = x(i,j,k) - cc(j)*x(i,j+1,k)
        End Do
     End Do
  End Do
End Subroutine tridc2

!----------------------------------------------------------------------------------85

Subroutine ptrid(m, n, a, b, c, x)
  Implicit none
  ! Solves the periodic tridiagonal system a[i]*x[i-1]+b[i]*x[i]+c[i]*x[i+1] = r[i]; r is passed in x
  Integer, intent(in) :: m, n
  Real(8), dimension(m,n), intent(inout) :: a, b, c, x
  Integer :: i, j
  Do j = 2,n
     Do i = 1,m
        b(i,j) =  b(i,j)*b(i,j-1)-a(i,j)*c(i,j-1)
        x(i,j) =  x(i,j)*b(i,j-1)-a(i,j)*x(i,j-1)
        a(i,j) =                 -a(i,j)*a(i,j-1)
        c(i,j) =  c(i,j)*b(i,j-1)
     End Do
  End Do
  Do i = 1,m
     x(i,n) = x(i,n)/(b(i,n)+a(i,n))
     b(i,n) = c(i,n)/(b(i,n)+a(i,n))
  End Do
  Do j = n-1,1,-1
     Do i = 1,m
        x(i,j) = (x(i,j)-c(i,j)*x(i,j+1)-a(i,j)*x(i,n))/b(i,j)
        b(i,j) = -(c(i,j)*b(i,j+1)+a(i,j)*b(i,n))/b(i,j)
     End Do
  End Do
  Do i = 1,m
     x(i,1) = x(i,1)/(1.d0+b(i,1))
  End Do
  Do j = 2,n
     Do i = 1,m
        x(i,j) = x(i,j)-b(i,j)*x(i,1)
     End Do
  End Do
End Subroutine ptrid

!----------------------------------------------------------------------------------85

Subroutine ptrid2nd(l, ldim, m, a, b, c, x)

use partition_data

  Implicit none
  ! Note: ldim is passed for dimensioning purposes only.
  Integer, intent(in) :: l, ldim, m
  Real(8), intent(in) :: a, b, c
  Real(8), intent(inout) ::  x(ldim,m)

  ! Solves the periodic tridiagonal system a*x(i-1)+b*x(i)+c*x(i+1) = r(i); r is passed in x
  Integer :: i, j
  Real(8) :: aa(m),bb(m),cc(m),dd(m)

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: ptrid2nd has been called'
#endif  
  
  aa(1) = a
  bb(1) = b
  cc(1) = c
  Do j = 2,m
     aa(j) = -a*aa(j-1)
     bb(j) = b*bb(j-1)-a*cc(j-1)
     cc(j) = c*bb(j-1)
  End Do
  dd(m) = cc(m)/(bb(m)+aa(m))
  Do j = m-1,1,-1
     dd(j) = -(cc(j)*dd(j+1)+aa(j)*dd(m))/bb(j)
  End Do
  Do j = 2,m
     Do i = 1,l
        x(i,j) = x(i,j)*bb(j-1)-a*x(i,j-1)
     End Do
  End Do
  Do i = 1,l
     x(i,m) = x(i,m)/(bb(m)+aa(m))
  End Do
  Do j = m-1,1,-1
     Do i = 1,l
        x(i,j) = (x(i,j)-cc(j)*x(i,j+1)-aa(j)*x(i,m))/bb(j)
     End Do
  End Do
  Do i = 1,l
     x(i,1) = x(i,1)/(1.d0+dd(1))
  End Do
  Do j = 2,m
     Do i = 1,l
        x(i,j) = x(i,j)-dd(j)*x(i,1)
     End Do
  End Do

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: Returning from ptrid2nd'
#endif  
  
End Subroutine ptrid2nd

!----------------------------------------------------------------------------------85

Subroutine antiptrid2nd(l, lb, m, a, b, c, x)
  Implicit none
  Integer, intent(in) :: l, lb, m
  Real(8), intent(in) :: a, b, c
  Real(8), intent(inout) ::  x(lb,m)

  ! Solves the periodic tridiagonal system a*x(i-1)+b*x(i)+c*x(i+1) = r(i); r is passed in x
  Integer :: i, j
  Real(8) :: aa(m),bb(m),cc(m),dd(m)
  aa(1) = -a
  bb(1) = b
  cc(1) = c
  Do j = 2,m
     aa(j) = -a*aa(j-1)
     bb(j) = b*bb(j-1)-a*cc(j-1)
     cc(j) = c*bb(j-1)
  End Do
  dd(m) = -cc(m)/(bb(m)+aa(m))
  Do j = m-1,1,-1
     dd(j) = -(cc(j)*dd(j+1)+aa(j)*dd(m))/bb(j)
  End Do
  Do j = 2,m
     Do i = 1,l
        x(i,j) = x(i,j)*bb(j-1)-a*x(i,j-1)
     End Do
  End Do
  Do i = 1,l
     x(i,m) = x(i,m)/(bb(m)+aa(m))
  End Do
  Do j = m-1,1,-1
     Do i = 1,l
        x(i,j) = (x(i,j)-cc(j)*x(i,j+1)-aa(j)*x(i,m))/bb(j)
     End Do
  End Do
  Do i = 1,l
     x(i,1) = x(i,1)/(1.d0+dd(1))
  End Do
  Do j = 2,m
     Do i = 1,l
        x(i,j) = x(i,j)-dd(j)*x(i,1)
     End Do
  End Do
End Subroutine antiptrid2nd

!----------------------------------------------------------------------------------85

