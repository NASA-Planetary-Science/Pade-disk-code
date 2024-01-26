! Tool for checking partitioning.

! Look at notes of 12/29/2017.

! Alan      Me
! x    <--> r
! y    <--> phi
! z    <--> z

implicit none

integer :: nr, nphi, nz
integer :: num_nodes
integer :: ng1, ng2   ! number of groups in the cross-plane of pencils
integer :: mr, mphi, mz_r, mz_phi
logical :: found

found = .false.

write (6, 1)
1 format (' enter nr, nz, nphi, num_nodes--->', $)
read (5, *) nr, nz, nphi, num_nodes

ng1 = 1
do while (ng1<=num_nodes)
   ng2 = num_nodes/ng1
   if (ng1*ng2 == num_nodes) then      ! num_nodes is divisible by these factors
      mr     = nr/ng1
      mphi   = nphi/ng2
      mz_r   = nz/ng1
      mz_phi = nz/ng2
      if (mr*ng1==nr .and. mphi*ng2==nphi .and. mz_r*ng1==nz .and. mz_phi*ng2==nz) then
         print *, ' ***** Found a solution *****' 
         print *, ' ng1 = ', ng1, ' mr = ', mr
         print *, ' ng2 = ', ng2, ' mphi = ', mphi
         print *, ' mz_r = ', mz_r, ' mz_phi = ', mz_phi
         found = .true.
      end if
    end if
    ng1 = ng1 + 1
end do

if (.not. found) print *, ' No solution found'
 

stop
end
