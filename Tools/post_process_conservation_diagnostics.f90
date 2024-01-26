!----------------------------------------------------------------------------------85

! Program to post-process conservation diagnostics.

implicit none
integer :: itype, count, i
real(8) :: t

real(8), allocatable, dimension(:) :: t_array
real(8), allocatable, dimension(:) :: vol_int, left_flux_out, right_flux_out, &
     top_flux_out,  bot_flux_out
real(8) :: dt, flux_out_total, left, right, bot, top
character(4), dimension(2) :: suffix
character(19) :: filename

suffix(1) = 'mass'
suffix(2) = 'amom'

write(6, 1)
1 format(' enter 1 or 2 correspoding to the file you want to process:',/, &
         '    1: mass_conservation.his',/, &
         '    2: amom_conservation.his',/, &
         ' --->', $)
read(5, *) itype

if (itype .eq. 1) then
   open(unit = 1, file = 'mass_conservation.his', form = 'formatted', &
        status = 'unknown')
else if (itype .eq. 2) then
      open(unit = 1, file = 'amom_conservation.his', form = 'formatted', &
           status = 'unknown')
end if

! Determine how many lines the file has.
200 continue
read(1, "(1x, e16.9)", end = 100) t
count = count + 1
go to 200

100 continue
print *, ' # lines in file = ', count
close(1)

allocate(t_array       (count))
allocate(vol_int       (count))
allocate(left_flux_out (count))
allocate(right_flux_out(count))
allocate(top_flux_out  (count))
allocate(bot_flux_out  (count))

if (itype .eq. 1) then
   open(unit = 1, file = 'mass_conservation.his', form = 'formatted', &
        status = 'unknown')
else if (itype .eq. 2) then
      open(unit = 1, file = 'amom_conservation.his', form = 'formatted', &
           status = 'unknown')
end if

do i = 1, count
   read(1, "(6(1x, e16.9))") t_array(i), vol_int(i), &
        left_flux_out(i), right_flux_out(i), &
        top_flux_out (i), bot_flux_out  (i)
   print *, ' line = ', i, ' t = ', t_array(i)
end do
close(1)

! Determine the rate of range of the volume integral.
write(filename, "('vol_integr_', a4, '.dat')") suffix(itype)
open(unit = 1, file = filename, form = 'formatted', status = 'unknown')
write(filename, "('out_fluxes_', a4, '.dat')") suffix(itype)
open(unit = 2, file = filename, form = 'formatted', status = 'unknown')
! Time integral of the flux out:
flux_out_total = 0.d0
left  = 0.d0
right = 0.d0
bot   = 0.d0
top   = 0.d0
do i = 2, count - 1
   dt = t_array(i) - t_array(i-1)
   ! Trapezoidal rule for time integration:
   left  = left  + 0.5 * dt * (left_flux_out (i-1) + left_flux_out (i))
   right = right + 0.5 * dt * (right_flux_out(i-1) + right_flux_out(i))
   bot   = bot   + 0.5 * dt * (bot_flux_out  (i-1) + bot_flux_out  (i))
   top   = top   + 0.5 * dt * (top_flux_out  (i-1) + top_flux_out  (i))
   flux_out_total = left + right + bot + top

   ! The third item below should be constant in the exact limit:
   write(1, "(3(1x, e16.9))") t_array(i), vol_int(i), vol_int(i) + flux_out_total
   write(2, "(6(1x, e16.9))") t_array(i), flux_out_total, left, right, bot, top   
end do

close(1)
close(2)
end

!----------------------------------------------------------------------------------85
