!----------------------------------------------------------------------------------85

implicit none
real(8) :: s, s1, s2, smin, smax, ds
integer :: n, i, isel
real(8) :: Melander3

write(6, "('enter 1 to test individual values, 2 to output for plot--->', $)")
read(5, *) isel

if (isel .eq. 1) then
   write(6, "(' enter  s, s1, s2--->', $)")
   read(5,*) s, s1, s2
   print *, ' Melander3 = ', Melander3(s, s1, s2)
else if (isel .eq. 2) then
   write(6, "(' enter range of plot: smin, smax--->', $)")
   read(5,*) smin, smax

   write(6, "(' enter bounds of transition region: s1, s2--->', $)")
   read(5,*) s1, s2

   n  = 500
   ds = (smax - smin) / (n - 1)

   open(unit = 1, file = 'Melander_func.dat', form = 'formatted', status = 'unknown')
   do i = 1, n
      s = smin + (i - 1)*ds
      write(1, "(2(1x, e12.5))") s, Melander3(s, s1, s2)
   end do
   print *, ' wrote Melander_func.dat'
   close(1)
end if

end

!----------------------------------------------------------------------------------85

real(8) function Melander3(s, s1, s2)

! Has three parts:
! 1, for s <= 1
! transition for s1 < s < s2
! 0 for s >= s1

implicit none
real(8) :: s, s1, s2

! Local:
real(8) :: r, paren, bracket
real(8), parameter :: kappa = 0.5d0 * exp(2.d0) * log(2.0d0)

if (s .lt. s1) then
   Melander3 = 1.0d0
else if (s .lt. s2) then
   r        = (s - s1) / (s2 - s1)
   paren    = 1.d0 / (r - 1)
   bracket  = - kappa/r * exp(paren)
   Melander3 = 1.d0 - exp(bracket)
else
   Melander3 = 0.0d0
end if

end function Melander3

!----------------------------------------------------------------------------------85
