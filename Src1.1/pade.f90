!**********************************************************************************85
!> Protoplanetary disk code in cylindrical coordinates (z, r, phi) employing
!> Pade differentiation. <br>
!
!> This is the main program. <br>

!> Degrees of freedom in the "q" array are: <br>
!>        1: Density: rho <br>
!>        2: Radial momentum: rho*ur <br>
!>        3: Vertical momentum: rho*uz <br>
!>        4: Angular momentum: rho*u_phi*r <br>
!>        5: Internal energy per unit volume, e = p/(gamma - 1) <br>

!>        If isothermal = .true. we don't use the energy equation at all and instead <br>
!>        use p = rho * ci_squared_initial(ir, iz) where ci_squared_initial is a <br>
!>        prescribed isothermal sound speed. <br>

! Note: For the Fargo option, we must use the internal energy as explained in the
!       Ap. J. manuscript.

!> For the isothermal option (isothermality = true) we don't need the energy
!> degree of freedom but to keep the coding simple, it is carried anyway without
!> harm.

!> Spatial differencing : 4th order Pade.
!> Time stepping        : 4th order Runge Kutta.

! For reference the pencil directions are:

! Alan      Me
! x    <--> r
! y    <--> phi
! z    <--> z

! Alan's x-space (y, z, q, x) ---> Our r space   (phi, z, q, r)
! Alan's y-space (x, z, q, y) ---> Our phi space (r, z, q, phi)
! Alan's z-space (x, y, z, q) ---> Our z space   (r, phi, z, q)

!> Dimensioning of each pencil:
!> q           (sr:er, sphi:ephi, nz, ndof)           : in z space
!> q_r_space   (sphi:ephi, sz_r:ez_r, ndof, nr)
!> q_phi_space (sr:er, sz_phi:ez_phi, ndof, nphi)

! Things to do:
!    (1) Make ci_squared locally dimensioned to each processor.
!    (2) Make F and Fd global this they are used often.
!
!----------------------------------------------------------------------------------85

! Main program
program pade

#ifdef mpi_code
   use mpi
#endif
use logical_units, only : lun_general_purpose
use partition_data

implicit none
integer :: i_run_type, ier

#ifdef mpi_code
   call mpi_init (ier)
   call mpi_comm_rank (mpi_comm_world, my_node,   ier)
   call mpi_comm_size (mpi_comm_world, num_nodes, ier)
#else
   print *, ' Serial run'
   my_node = 0 ! For a serial run.
#endif

! This subroutine does the broadcast:
call read_namelist_run_type (lun_general_purpose, i_run_type)

if (i_run_type .eq. 0) then
   call user_application
else if (i_run_type .eq. 1) then
   call app_euler1d_tests ! in app_euler1d_tests.f90
else if (i_run_type .eq. 2) then
   call app_hydrostatic_test    ! in app_hydrostatic_test.f90
else if (i_run_type .eq. 3) then
   ! in app_homentropic_solid_body_rotation_test.f90
   call app_homentropic_solid_body_rotation_test
else if (i_run_type .eq. 4) then
   call app_single_vortex_fargo_test
else if (i_run_type .eq. 5) then
   call app_vsi_3D
else if (i_run_type .eq. 6) then
   call app_vortex_pair
else if (i_run_type .eq. 7) then
   call app_taylor_couette
#ifdef advection_test
else if (i_run_type .eq. 8) then
   print *, ' my_node = ', my_node, ' calling app_advcection_test'
   call app_advection_test
#endif
else if (i_run_type .eq. 9) then
   call app_vsi_with_dt_output
else
   if (my_node .eq. 0) print *, ' node 0: Unknown i_run_type = ', i_run_type
   call terminate_with_no_save(1)         
end if

end program pade

!----------------------------------------------------------------------------------85

subroutine read_namelist_run_type (lun_general_purpose, i_run_type)

! Reads the number of the application subroutine to be invoked.

#ifdef mpi_code
   use mpi
#endif
use partition_data
implicit none
integer :: lun_general_purpose, i_run_type, ier

namelist /run_type/ i_run_type

! Read namelist input:

if (my_node .eq. 0) then   
   open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', status = 'old')
   read (lun_general_purpose, nml = run_type)
   close (lun_general_purpose)
end if

#ifdef mpi_code
   call mpi_bcast(i_run_type, 1, mpi_integer, 0, mpi_comm_world, ier)
#endif

if (my_node .eq. 0) then
   print *, ' node 0: i_run_type = ', i_run_type, ' in input_file'
end if
   
end subroutine read_namelist_run_type

!----------------------------------------------------------------------------------85
