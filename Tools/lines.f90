implicit none
real(8) :: k0, k1, E0, E1, num, den, slope, A

write(6, "(' enter starting x,y of line--->', $)")
read(5, *) k0, E0

write(6, "(' enter ending x of line--->', $)")
read(5, *) k1

write(6, "(' enter num and den in slope of line--->', $)")
read(5, *) num, den

slope = num/den

! Determine constant:
A = E0 / k0**slope

! Therefore
E1 = A*k1**slope

open(unit = 1, file = 'line.dat', form = 'formatted', status = 'unknown')
write(1, "(2(1x, e12.5))") k0, E0
write(1, "(2(1x, e12.5))") k1, E1
close(1)

end




