implicit none
integer isel
real(8) slope, x1, x2, y1, y2, A

write(6, 1)
1 format(' enter 1 to get a line of desired slope',/, &
         '       2 to get power law by entering two data point',/, &
         ' --->', $)
read(5, *) isel

if (isel .eq. 1) then
   write(6, "(' enter desired slope on log-log--->', $)")
   read(5, *) slope
   write(6, "(' enter x y coordinates of first point--->', $)")
   read(5, *) x1, y1
   write(6, "(' enter x coordinate of second point--->', $)")
   read(5, *) x2

   A  = y1 / x1**slope
   y2 = A*x2**slope

   open(unit = 1, file = 'log_log_line.dat', form = 'formatted', status = 'unknown')
   write(1, "(2(1x, e12.5))") x1, y1
   write(1, "(2(1x, e12.5))") x2, y2
   close(1)
   print *, ' wrote log_log_line.dat'
else if (isel .eq. 2) then
   write(6, "(' enter x y coordinates of first point--->', $)")
   read(5, *) x1, y1
   write(6, "(' enter x y coordinates of second point--->', $)")
   read(5, *) x2, y2
   slope = (log10(y2) - log10(y1)) / (log10(x2) - log10(x1))
   print *, ' slope = ', slope
end if

end

