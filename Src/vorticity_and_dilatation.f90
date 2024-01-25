!----------------------------------------------------------------------------------85

! subroutine tecplot_vort_and_dil_in_merid_horiz_planes
! subroutine tecplot_vort_z_and_dil  ! Intended for (r, phi) 2D runs.
! subroutine tecplot_vort_and_dil
! subroutine vort_and_dil

!----------------------------------------------------------------------------------85

subroutine tecplot_vort_and_dil_in_merid_horiz_planes(perturbation, compute_curl_rho_u, &
     q, iphi_plot, iz_plot, L_scale, T_scale, rho_scale, t, istep)

! iphi_plot = desired meridional plane.  iz_plot = desired horizontal plane.
! If iphi_plot = 0 we will skip the meridional plane.  If iz_plot = 0 we will skip the horizontal
! plane.

use grid
use xy_coordinates_of_grid
use logical_units, only: lun_tecplot
use partition_data
use transposes_of_q_and_qdot, only: q_r_space, q_phi_space, qdot_r_space, &
     qdot_phi_space
use dof_indices
implicit none
logical :: perturbation ! If .true. the perturbation vorticity will be output.
logical :: compute_curl_rho_u
integer :: iphi_plot, iz_plot, istep
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in)  :: q
real(8) :: t, L_scale, T_scale
real(8) :: rho_scale ! Needed if you compute curl(rho*u).

! Locals:
! Make these global or passed in the argument list.  Note: These are passed as work
! to vort_and_dil which further passes them down to the derivative routines where
! they will be dimensioned differently but have the same size.
real(8), dimension(sr:er, sphi:ephi, nz) :: F, dF

! Vorticity and dilatation (which is in the fourth slot).
real(8), dimension(sr:er, sphi:ephi, nz, 4) :: v 

character(80) :: filename
integer :: ir, iz, iphi, irec, nz_plot, iphi_rec
logical :: i_have_iphi_plot
logical, parameter :: divide_by_rho = .false.
real(8) :: t_scaled
real(8) :: omega_r, omega_z, omega_p, vort_mag, inverse_vorticity_scale

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: tecplot_vort_and_z_in_meridional_plane has been called'
#endif

! In case the user wants a scaled output:
if (compute_curl_rho_u) then
   inverse_vorticity_scale = T_scale / rho_scale  ! rho * u / L = rho / T
else
   inverse_vorticity_scale = T_scale
end if   

! Calculate the vorticity and dilatation.
! q_dot_r_space and q_dot_phi_space are being passed as work space, i.e., for
! the arguments v_r_space and v_phi_space, respectively.
call vort_and_dil(perturbation, divide_by_rho, compute_curl_rho_u, &
        q, v, q_r_space, q_phi_space, qdot_r_space, qdot_phi_space, F, dF)

t_scaled = t/T_scale

! Output meridional plane:
if (iphi_plot .ne. 0) then
   if (suppress_z_derivatives_when_nz_not_1) then
      nz_plot = 1
   else
      nz_plot = nz
   end if

   if (compute_curl_rho_u) then
      if (perturbation) then
         write (filename, 71) iphi_plot, int(t_scaled), t_scaled - int(t_scaled), istep
71       format ('pert_curl_rho_u_and_dil_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
      else
         write (filename, 72) iphi_plot, int(t_scaled), t_scaled - int(t_scaled), istep
72       format ('total_curl_rho_u_and_dil_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')         
      end if
   else
      if (perturbation) then
         write (filename, 73) iphi_plot, int(t_scaled), t_scaled - int(t_scaled), istep
73       format ('pert_vort_and_dil_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
      else
         write (filename, 74) iphi_plot, int(t_scaled), t_scaled - int(t_scaled), istep
74       format ('total_vort_and_dil_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')         
      end if
   end if
   
   open (unit = lun_tecplot, file = filename, form = 'formatted', &
         status = 'unknown', access = 'direct',recl = 13*7 + 1)

   if (my_node .eq. 0) then
      ! char(10) = new line
      write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
      1  format ('TITLE = "t = ', e12.5, '"', 2a1)
      write (lun_tecplot, 2, rec = 2) char(13), char(10)
      2  format('VARIABLES = "r", "z", "wr", "wz", "wphi", "wmag", "dil"', 2a1)
      write (lun_tecplot, 3, rec = 3) nr, nz, char(13), char(10)
      3 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', 2a1)
   end if

   ! Write only of this processor has the iphi:
   i_have_iphi_plot = (iphi_plot .ge. sphi) .and. (iphi_plot .le. ephi)
   if (i_have_iphi_plot) then
      do iz = 1, nz_plot
         do ir = sr, er
            irec = (iz - 1)*nr + ir + 3
            omega_r = v(ir, iphi_plot, iz, r_comp)*inverse_vorticity_scale
            omega_z = v(ir, iphi_plot, iz, z_comp)*inverse_vorticity_scale
            omega_p = v(ir, iphi_plot, iz, p_comp)*inverse_vorticity_scale
            vort_mag = SQRT(omega_r**2 + omega_z**2 + omega_p**2)
            ! Why no char(13) = CR
            write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) &
                 rgrid(ir)/L_scale, zgrid(iz)/L_scale, &
                 omega_r, omega_z, omega_p, vort_mag, &
                 v(ir, iphi_plot, iz, 4)*T_scale, char(10)  ! dil.      
         end do
      end do
   end if
   close(lun_tecplot)
   if (my_node .eq. 0) print *, ' node 0: wrote tecplot file =', filename
end if

! Output horizontal plane if desired and for a non-axisymmetric run:
if ((iz_plot .ne. 0) .and. (nphi .ne. 1)) then
   if (compute_curl_rho_u) then
      if (perturbation) then
         write (filename, 81) iz_plot, int(t_scaled), t_scaled - int(t_scaled), istep
81       format ('pert_curl_rho_u_and_dil_iz_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
      else
         write (filename, 82) iz_plot, int(t_scaled), t_scaled - int(t_scaled), istep
82       format ('total_curl_rho_u_and_dil_iz_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
      end if
   else
      if (perturbation) then
         write (filename, 83) iz_plot, int(t_scaled), t_scaled - int(t_scaled), istep
83       format ('pert_vort_and_dil_iz_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
      else
         write (filename, 84) iz_plot, int(t_scaled), t_scaled - int(t_scaled), istep
84       format ('total_vort_and_dil_iz_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')         
      end if
   end if
   
   open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
        access = 'direct',    recl = 13*7 + 1)

   if (my_node .eq. 0) then
      ! char(10) = new line
      write (lun_tecplot, 11, rec = 1) t, char(13), char(10)
      11  format ('TITLE = "t = ', e12.5, '"', 2a1)
      write (lun_tecplot, 21, rec = 2) char(13), char(10)
      21  format('VARIABLES = "x", "y", "wr", "wz", "wphi", "wmag", "dil"', 2a1)
      write (lun_tecplot, 31, rec = 3) nr, nphi+1, char(13), char(10)
      31 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', 2a1)
   end if

   do iphi = sphi, ephi
      do ir = sr, er
         irec = (iphi - 1)*nr + ir + 3
         omega_r = v(ir, iphi, iz_plot, r_comp)*inverse_vorticity_scale
         omega_z = v(ir, iphi, iz_plot, z_comp)*inverse_vorticity_scale
         omega_p = v(ir, iphi, iz_plot, p_comp)*inverse_vorticity_scale
         vort_mag = SQRT(omega_r**2 + omega_z**2 + omega_p**2)
         ! Why no char(13)
         write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) &
                 xgrid(ir,iphi)/L_scale, ygrid(ir,iphi)/L_scale, &
                 omega_r, omega_z, omega_p, vort_mag, &
                 v(ir, iphi, iz_plot, 4)*T_scale, char(10)  ! dil.      
      end do
   end do

   ! Periodic completion in phi:
   if (sphi .eq. 1) then
      ! The data for iphi = 1, gets put in the records for nphi + 1
      iphi_rec = nphi + 1
      iphi     = 1
      do ir = sr, er
         irec = (iphi_rec - 1)*nr + ir + 3

         omega_r = v(ir, iphi, iz_plot, r_comp)*inverse_vorticity_scale
         omega_z = v(ir, iphi, iz_plot, z_comp)*inverse_vorticity_scale
         omega_p = v(ir, iphi, iz_plot, p_comp)*inverse_vorticity_scale
         vort_mag = SQRT(omega_r**2 + omega_z**2 + omega_p**2)
         
         ! Note: I have stored the coordinates for this line in the iphi = 0 slot
         ! since the flow data sits in the same processor:
         ! Why no char(13)?
         write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) &
                 xgrid(ir,0)/L_scale, ygrid(ir,0)/L_scale, &
                 omega_r, omega_z, omega_p, vort_mag, &
                 v(ir, iphi, iz_plot, 4)*T_scale, char(10)  ! dil.      
      end do         
   end if
   
   close(lun_tecplot)
   if (my_node .eq. 0) print *, ' node 0: wrote tecplot file =', filename
end if

end subroutine tecplot_vort_and_dil_in_merid_horiz_planes

!----------------------------------------------------------------------------------85

subroutine tecplot_vort_and_dil_in_many_merid_horiz_planes(perturbation, &
     compute_curl_rho_u, q, nphi_planes, iphi_list, nz_planes, iz_list, &
     z_list, L_scale, T_scale, rho_scale, t, istep)

! iphi_plot = desired meridional plane.  iz_plot = desired horizontal plane.
! If iphi_plot = 0 we will skip the meridional plane.  If iz_plot = 0 we will skip the horizontal
! plane.

use grid
use xy_coordinates_of_grid
use logical_units, only: lun_tecplot
use partition_data
use transposes_of_q_and_qdot, only: q_r_space, q_phi_space, qdot_r_space, &
     qdot_phi_space
use dof_indices
implicit none
logical, intent(in) :: perturbation ! If .true. the perturbation vorticity will be output.
logical, intent(in) :: compute_curl_rho_u
integer, intent(in) :: nphi_planes, nz_planes, istep
integer, intent(in) :: iphi_list(nphi_planes), iz_list(nz_planes)
! z coordinates of the planes for file naming only:
real(8), intent(in), dimension(nz_planes) :: z_list
real(8), intent(in), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8), intent(in) :: t, L_scale, T_scale
real(8), intent(in) :: rho_scale ! Needed if you compute curl(rho*u).

! Locals:
! Make these global or passed in the argument list.  Note: These are passed as work
! to vort_and_dil which further passes them down to the derivative routines where
! they will be dimensioned differently but have the same size.
real(8), dimension(sr:er, sphi:ephi, nz) :: F, dF

! Vorticity and dilatation (which is in the fourth slot).
real(8), dimension(sr:er, sphi:ephi, nz, 4) :: v 

character(80) :: filename
integer :: ir, iz, iphi, irec, nz_plot, iphi_rec
logical :: i_have_iphi_plot
logical, parameter :: divide_by_rho = .false.
real(8) :: t_scaled
real(8) :: omega_r, omega_z, omega_p, vort_mag, inverse_vorticity_scale
integer :: iphi_plane, iz_plane

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: tecplot_vort_and_z_in_meridional_plane has been called'
#endif

! In case the user wants a scaled output:
if (compute_curl_rho_u) then
   inverse_vorticity_scale = T_scale / rho_scale  ! rho * u / L = rho / T
else
   inverse_vorticity_scale = T_scale
end if   

! Calculate the vorticity and dilatation.
! q_dot_r_space and q_dot_phi_space are being passed as work space, i.e., for
! the arguments v_r_space and v_phi_space, respectively.
call vort_and_dil(perturbation, divide_by_rho, compute_curl_rho_u, &
        q, v, q_r_space, q_phi_space, qdot_r_space, qdot_phi_space, F, dF)

t_scaled = t/T_scale

! Output meridional planes:
do iphi_plane = 1, nphi_planes
   iphi = iphi_list(iphi_plane)
   if (suppress_z_derivatives_when_nz_not_1) then
      nz_plot = 1
   else
      nz_plot = nz
   end if

   if (compute_curl_rho_u) then
      if (perturbation) then
         write (filename, 71) iphi, int(t_scaled), t_scaled - int(t_scaled), istep
         71 format ('pert_curl_rho_u_and_dil_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
      else
         write (filename, 72) iphi, int(t_scaled), t_scaled - int(t_scaled), istep
         72 format ('total_curl_rho_u_and_dil_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')         
      end if
   else
      if (perturbation) then
         write (filename, 73) iphi, int(t_scaled), t_scaled - int(t_scaled), istep
         73  format ('pert_vort_and_dil_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
      else
         write (filename, 74) iphi, int(t_scaled), t_scaled - int(t_scaled), istep
         74 format ('total_vort_and_dil_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')         
      end if
   end if
   
   open (unit = lun_tecplot, file = filename, form = 'formatted', &
         status = 'unknown', access = 'direct',recl = 13*7 + 1)

   if (my_node .eq. 0) then
      ! char(10) = new line
      write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
      1  format ('TITLE = "t = ', e12.5, '"', 2a1)
      write (lun_tecplot, 2, rec = 2) char(13), char(10)
      2  format('VARIABLES = "r", "z", "wr", "wz", "wphi", "wmag", "dil"', 2a1)
      write (lun_tecplot, 3, rec = 3) nr, nz, char(13), char(10)
      3 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', 2a1)
   end if

   ! Write only of this processor has the iphi:
   i_have_iphi_plot = (iphi .ge. sphi) .and. (iphi .le. ephi)
   if (i_have_iphi_plot) then
      do iz = 1, nz_plot
         do ir = sr, er
            irec = (iz - 1)*nr + ir + 3
            omega_r = v(ir, iphi, iz, r_comp)*inverse_vorticity_scale
            omega_z = v(ir, iphi, iz, z_comp)*inverse_vorticity_scale
            omega_p = v(ir, iphi, iz, p_comp)*inverse_vorticity_scale
            vort_mag = SQRT(omega_r**2 + omega_z**2 + omega_p**2)
            ! Why no char(13) = CR
            write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) &
                 rgrid(ir)/L_scale, zgrid(iz)/L_scale, &
                 omega_r, omega_z, omega_p, vort_mag, &
                 v(ir, iphi, iz, 4)*T_scale, char(10)  ! dil.      
         end do
      end do
   end if
   close(lun_tecplot)
   if (my_node .eq. 0) print *, ' node 0: wrote tecplot file =', filename
end do

! Output horizontal plane if desired and for a non-axisymmetric run:
do iz_plane = 1, nz_planes
   iz = iz_list(iz_plane)
   if (compute_curl_rho_u) then
      if (perturbation) then
         write (filename, 81) iz, int(t_scaled), t_scaled - int(t_scaled), istep
         81 format ('pert_curl_rho_u_and_dil_iz_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
      else
         write (filename, 82) iz, int(t_scaled), t_scaled - int(t_scaled), istep
         82 format ('total_curl_rho_u_and_dil_iz_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
      end if
   else
      if (perturbation) then
         write (filename, 83) iz, int(t_scaled), t_scaled - int(t_scaled), istep
         83 format ('pert_vort_and_dil_iz_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
      else
         write (filename, 84) iz, int(t_scaled), t_scaled - int(t_scaled), istep
         84 format ('total_vort_and_dil_iz_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')         
      end if
   end if
   
   open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
        access = 'direct', recl = 13*7 + 1)

   if (my_node .eq. 0) then
      ! char(10) = new line
      write (lun_tecplot, 11, rec = 1) t, char(13), char(10)
      11  format ('TITLE = "t = ', e12.5, '"', 2a1)
      write (lun_tecplot, 21, rec = 2) char(13), char(10)
      21  format('VARIABLES = "x", "y", "wr", "wz", "wphi", "wmag", "dil"', 2a1)
      write (lun_tecplot, 31, rec = 3) nr, nphi+1, char(13), char(10)
      31 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', 2a1)
   end if

   do iphi = sphi, ephi
      do ir = sr, er
         irec = (iphi - 1)*nr + ir + 3
         omega_r = v(ir, iphi, iz, r_comp)*inverse_vorticity_scale
         omega_z = v(ir, iphi, iz, z_comp)*inverse_vorticity_scale
         omega_p = v(ir, iphi, iz, p_comp)*inverse_vorticity_scale
         vort_mag = SQRT(omega_r**2 + omega_z**2 + omega_p**2)
         ! Why no char(13)
         write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) &
                 xgrid(ir,iphi)/L_scale, ygrid(ir,iphi)/L_scale, &
                 omega_r, omega_z, omega_p, vort_mag, &
                 v(ir, iphi, iz, 4)*T_scale, char(10)  ! dil.      
      end do
   end do

   ! Periodic completion in phi:
   if (sphi .eq. 1) then
      ! The data for iphi = 1, gets put in the records for nphi + 1
      iphi_rec = nphi + 1
      iphi     = 1
      do ir = sr, er
         irec = (iphi_rec - 1)*nr + ir + 3

         omega_r = v(ir, iphi, iz, r_comp)*inverse_vorticity_scale
         omega_z = v(ir, iphi, iz, z_comp)*inverse_vorticity_scale
         omega_p = v(ir, iphi, iz, p_comp)*inverse_vorticity_scale
         vort_mag = SQRT(omega_r**2 + omega_z**2 + omega_p**2)
         
         ! Note: I have stored the coordinates for this line in the iphi = 0 slot
         ! since the flow data sits in the same processor:
         ! Why no char(13)?
         write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) &
                 xgrid(ir,0)/L_scale, ygrid(ir,0)/L_scale, &
                 omega_r, omega_z, omega_p, vort_mag, &
                 v(ir, iphi, iz, 4)*T_scale, char(10)  ! dil.      
      end do         
   end if
   
   close(lun_tecplot)
   if (my_node .eq. 0) print *, ' node 0: wrote tecplot file =', filename
end do

end subroutine tecplot_vort_and_dil_in_many_merid_horiz_planes

!----------------------------------------------------------------------------------85

subroutine tecplot_vort_z_and_dil(q, t, perturbation, divide_by_rho, &
                                  compute_curl_rho_u, prefix)

use grid
use xy_coordinates_of_grid
use logical_units, only: lun_tecplot
use partition_data
use transposes_of_q_and_qdot, only: q_r_space, q_phi_space, qdot_r_space, &
     qdot_phi_space
use dof_indices
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in)  :: q
real(8) :: t
logical :: perturbation, divide_by_rho, compute_curl_rho_u
character(6) :: prefix 

! Locals:
! Make these global or passed in the argument list.  Note: These are passed as work
! to vort_and_dil which further passes them down to the derivative routines where
! they will be dimensioned differently but have the same size.
real(8), dimension(sr:er, sphi:ephi, nz) :: F, dF

! Vorticity and dilatation (which is in the fourth slot).
real(8), dimension(sr:er, sphi:ephi, nz, 4) :: v 

character(80) :: filename
integer :: ir, iphi, iz, irec, icomp, iphi_rec, nz_plot

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: tecplot_vort_z_and_dil has been called'
#endif

if (suppress_z_derivatives_when_nz_not_1) then
   nz_plot = 1
else
   nz_plot = nz
end if

! Calculate the vorticity and dilatation.
! q_dot_r_space and q_dot_phi_space are being passed as work space, i.e., for
! the arguments v_r_space and v_phi_space, respectively.
call vort_and_dil(perturbation, divide_by_rho, compute_curl_rho_u, &
     q, v, q_r_space, q_phi_space, qdot_r_space, qdot_phi_space, F, dF)

if (.not. compute_curl_rho_u) then
   write (filename, 7) prefix, int(t), t - int(t)
   7 format (a6, '_vort_and_dil_t_', i6.6, f0.4, '.tec')   
else
   write (filename, 71) prefix, int(t), t - int(t)
   71 format (a6, '_curl_rho_u_and_dil_t_', i6.6, f0.4, '.tec')   
end if

! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
     access = 'direct', recl = 15*5 + 1) ! 3 coordinates + 1 vorticity components + dilatation recl = 76

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
1  format (' TITLE = "t = ', e12.5, '"', a1, a1)

   if (perturbation .and. divide_by_rho) then
      write (lun_tecplot, 2, rec = 2) char(13), char(10)
2     format(' VARIABLES = "x", "y", "z", "wz p/rho", "dil"', a1, a1) ! recl = 47
   else if (perturbation .and. .not. divide_by_rho) then
      write (lun_tecplot, 3, rec = 2) char(13), char(10)
3     format(' VARIABLES = "x", "y", "z", "wz p", "dil"', a1, a1)
   else if (.not. perturbation .and. divide_by_rho) then
      write (lun_tecplot, 4, rec = 2) char(13), char(10)
4     format(' VARIABLES = "x", "y", "z", "wz/rho", "dil"', a1, a1)
   else if (.not. perturbation .and. .not. divide_by_rho) then
      write (lun_tecplot, 5, rec = 2) char(13), char(10)
5     format(' VARIABLES = "x", "y", "z", "wz", "dil"', a1, a1)            
   end if
   
   write (lun_tecplot, 6, rec = 3) nr, nphi+1, nz_plot, char(13), char(10)
6  format(' ZONE I=',i4,',',' J=',i4,',',' K=',i4,',',' DATAPACKING=POINT', a1, a1)
end if

do iz = 1, nz_plot
   do iphi = sphi, ephi
      do ir = sr, er
         ! nphi+1 leaves room for periodic completion in phi:
         irec = ir + (iphi-1)*nr + (iz-1)*nr*(nphi+1) + 3  
         write(lun_tecplot, "(5(1x, e14.7), a1)", rec = irec) &
              xgrid(ir,iphi), ygrid(ir,iphi), zgrid(iz), &
              v(ir, iphi, iz, z_comp), v(ir, iphi, iz, 4), char(10)  ! z-component, dil.
      end do
   end do
end do

! Periodic completion in phi:
if (sphi .eq. 1) then
   ! The data for iphi = 1, gets put in the records for nphi + 1
   iphi_rec = nphi + 1
   iphi     = 1
   do iz = 1, nz_plot
      do ir = sr, er
         irec = ir + (iphi_rec-1)*nr + (iz-1)*nr*(nphi+1) + 3
         write(lun_tecplot, "(5(1x, e14.7), a1)", rec = irec) &
              xgrid(ir,0), ygrid(ir,0), zgrid(iz), &
              v(ir, iphi, iz, z_comp), v(ir, iphi, iz, 4), char(10) ! z-comp, dil.
      end do
   end do
end if
close(lun_tecplot)

if (my_node .eq. 0) print *, ' node 0: Wrote z-vorticity and dilatation tecplot file ', filename

end subroutine tecplot_vort_z_and_dil

!----------------------------------------------------------------------------------85

subroutine tecplot_vort_and_dil(q, t, perturbation, divide_by_rho, &
     compute_curl_rho_u)

use grid
use xy_coordinates_of_grid
use logical_units, only: lun_tecplot
use partition_data
use transposes_of_q_and_qdot, only: q_r_space, q_phi_space, qdot_r_space, &
     qdot_phi_space
use dof_indices
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in)  :: q
real(8) :: t
logical :: perturbation, divide_by_rho, compute_curl_rho_u
!character(6) :: prefix 

! Locals:
! Make these global or passed in the argument list.  Note: These are passed as work
! to vort_and_dil which further passes them down to the derivative routines where
! they will be dimensioned differently but have the same size.
real(8), dimension(sr:er, sphi:ephi, nz) :: F, dF

! Vorticity and dilatation (which is in the fourth slot).
real(8), dimension(sr:er, sphi:ephi, nz, 4) :: v 

character(80) :: filename
integer :: ir, iphi, iz, irec, icomp, iphi_rec, nz_plot
real(8) :: omega_z, omega_r, omega_p, vort_mag

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: tecplot_vort_and_dil has been called'
#endif

if (suppress_z_derivatives_when_nz_not_1) then
   nz_plot = 1
else
   nz_plot = nz
end if

! Calculate the vorticity and dilatation.
! q_dot_r_space and q_dot_phi_space are being passed as work space, i.e., for
! the arguments v_r_space and v_phi_space, respectively.
call vort_and_dil(perturbation, divide_by_rho, compute_curl_rho_u, q, v, q_r_space, q_phi_space, &
     qdot_r_space, qdot_phi_space, F, dF)

if (compute_curl_rho_u) then
   write (filename, 71) int(t), t - int(t)
   71  format ('curl_rho_u_and_dil_t_', i6.6, f0.4, '.tec')   
else
   write (filename, 7) int(t), t - int(t)
   7  format ('vort_and_dil_t_', i6.6, f0.4, '.tec')
end if

! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
     access = 'direct', recl = 15*8 + 1) ! 3 coordinates + 3 vorticity components + mag. + dilatation

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
1  format (' TITLE = "t = ', e12.5, '"', a1, a1)

   if (perturbation .and. divide_by_rho) then
      write (lun_tecplot, 2, rec = 2) char(13), char(10)
2     format(' VARIABLES = "x", "y", "z", "wz p/rho", "wr p/rho", "wphi p/rho", "mag", "dil"', a1, a1) ! recl = 81
   else if (perturbation .and. .not. divide_by_rho) then
      write (lun_tecplot, 3, rec = 2) char(13), char(10)
3     format(' VARIABLES = "x", "y", "z", "wz p", "wr p", "wphi p", "dil", "mag"', a1, a1)
   else if (.not. perturbation .and. divide_by_rho) then
      write (lun_tecplot, 4, rec = 2) char(13), char(10)
4     format(' VARIABLES = "x", "y", "z", "wz/rho", "wr/rho", "wphi/rho", "dil", "mag"', a1, a1)
   else if (.not. perturbation .and. .not. divide_by_rho) then
      write (lun_tecplot, 5, rec = 2) char(13), char(10)
5     format(' VARIABLES = "x", "y", "z", "wz", "wr", "wphi", "dil", "mag"', a1, a1)
   end if
   
   write (lun_tecplot, 6, rec = 3) nr, nphi+1, nz_plot, char(13), char(10)
6  format(' ZONE I=',i4,',',' J=',i4,',',' K=',i4,',',' DATAPACKING=POINT', a1, a1)
end if

do iz = 1, nz_plot
   do iphi = sphi, ephi
      do ir = sr, er
         ! nphi+1 leaves room for periodic completion in phi:
         irec = ir + (iphi-1)*nr + (iz-1)*nr*(nphi+1) + 3
         omega_z  = v(ir, iphi, iz, z_comp  )
         omega_r  = v(ir, iphi, iz, r_comp  )
         omega_p  = v(ir, iphi, iz, phi_comp)
         vort_mag = SQRT(omega_z**2 + omega_r**2 + omega_p**2) 
         write(lun_tecplot, "(8(1x, e14.7), a1)", rec = irec) &
              xgrid(ir,iphi), ygrid(ir,iphi), zgrid(iz), &
              omega_z, omega_r, omega_p, vort_mag, v(ir, iphi, iz, 4), char(10)  ! z, r, phi components.
      end do
   end do
end do

! Periodic completion in phi:
if (sphi .eq. 1) then
   ! The data for iphi = 1, gets put in the records for nphi + 1
   iphi_rec = nphi + 1
   iphi     = 1
   do iz = 1, nz_plot
      do ir = sr, er
         irec = ir + (iphi_rec-1)*nr + (iz-1)*nr*(nphi+1) + 3
         omega_z  = v(ir, iphi, iz, z_comp  )
         omega_r  = v(ir, iphi, iz, r_comp  )
         omega_p  = v(ir, iphi, iz, phi_comp)
         vort_mag = SQRT(omega_z**2 + omega_r**2 + omega_p**2) 
         write(lun_tecplot, "(8(1x, e14.7), a1)", rec = irec) &
              xgrid(ir,0), ygrid(ir,0), zgrid(iz), &
              omega_z, omega_r, omega_p, vort_mag, v(ir, iphi, iz, 4), char(10)  ! z, r, phi components.
      end do
   end do
end if
close(lun_tecplot)

if (my_node .eq. 0) print *, ' node 0: Wrote vort. and dil. tecplot file ', filename

end subroutine tecplot_vort_and_dil

!----------------------------------------------------------------------------------85

subroutine vort_and_dil (perturbation, divide_by_rho, compute_curl_rho_u, &
     q, v, q_r_space, q_phi_space, v_r_space, &
     v_phi_space, F, dF)

! Given the first 4 elements of q (i.e., rho and the 3 momenta) calculate the
! the 3 components of the vorticity vector and the dilatation.
! If divide_by_rho = true we calculate vorticity/rho.

! The determinant form for the curl in cylindrical coordinates is used.
! curl u = | zhat rhat r*phi_hat |
!          | d/dz d/dr d/dphi    |
!          | uz   ur   r*uphi    |  all divided by r.

! dil = ddz uz + 1/r ddr r*ur + 1/r ddphi uphi.

! Work arrays
!    for transposes: q_r_space, q_phi_space, v_r_space, v_phi_space
!    for taking a derivative: F, dF

use grid
use partition_data
use dof_indices
implicit none
! If true both we will  return the dilatation and vorticity relative to the
! basic state.
logical, intent(in) :: perturbation
! If true we will return omega/rho for the vorticity
logical, intent(in) :: divide_by_rho
! If .true. we will return the curl of rho*u for the vorticity.
logical, intent(in) :: compute_curl_rho_u 
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in)  :: q
! The fourth index of "v" is the dilatation.
real(8), dimension(sr:er, sphi:ephi, nz, 4), intent(out) :: v

! Work arrays which are passed to the routine:
real(8), dimension(sphi:ephi, sz_r:ez_r, ndof, nr)   :: q_r_space
real(8), dimension(sr:er, sz_phi:ez_phi, ndof, nphi) :: q_phi_space

real(8), dimension(sphi:ephi, sz_r:ez_r, 4, nr)   :: v_r_space
real(8), dimension(sr:er, sz_phi:ez_phi, 4, nphi) :: v_phi_space

real(8), dimension(sr:er, sphi:ephi, nz) :: F, dF

integer :: icomp, iz, iphi, ir

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: vort_and_dil has been called.  nphi = ', nphi
   print *, ' divide_by_rho  = ', divide_by_rho
   print *, ' perturbation   = ', perturbation
end if
#endif

v = 0.0d0

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: Finished zeroing v'
end if
#endif


#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: About to call add_z_derivatines_to_vort_and_dil'
end if
#endif

if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) &
     call add_z_derivatives_to_vort_and_dil (perturbation, compute_curl_rho_u, &
     q, v, F, dF)

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: About to check if nr .ne. 1'
end if
#endif

if (nr .ne. 1) then
   call transpose_z_to_r (4, q, q_r_space)
   call transpose_z_to_r (4, v, v_r_space)
   call add_r_derivatives_to_vort_and_dil (perturbation, compute_curl_rho_u, &
        q_r_space, v_r_space, F, dF)
   call transpose_r_to_z (4, v_r_space, v)   
end if

if (nphi .ne. 1) then
#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: In vort_and_dil: About to transpose q to phi space'
end if
#endif   
   call transpose_z_to_phi (4, q, q_phi_space)
#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: In vort_and_dil: About to transpose v to phi space'
end if
#endif
   call transpose_z_to_phi (4, v, v_phi_space)
   ! Note: No perturbation argument because the basic state has no
   ! phi derivatives.
#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: In vort_and_dil: About to add_phi_derivatives...'
end if
#endif   
   call add_phi_derivatives_to_vort_and_dil (compute_curl_rho_u, q_phi_space, &
        v_phi_space, F, dF)
   call transpose_phi_to_z (4, v_phi_space, v)   
end if

! Finally, the z and r components need division by r.

do icomp = 1, 2
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            v(ir,iphi,iz,icomp) = v(ir,iphi,iz,icomp) / rgrid(ir)
         end do
      end do
   end do
end do

if (divide_by_rho) then
   do icomp = 1, 3
      do iz = 1, nz
         do iphi = sphi, ephi
            do ir = sr, er
               v(ir,iphi,iz,icomp) = v(ir,iphi,iz,icomp) / q(ir,iphi,iz,irho)
            end do
         end do
      end do
   end do
end if

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: About to return from subroutine vort_sand_dil'
end if
#endif

end subroutine vort_and_dil

!----------------------------------------------------------------------------------85

subroutine add_z_derivatives_to_vort_and_dil (perturbation, compute_curl_rho_u, &
                                              q, v, F, dF)

! q and v must be in z-space.  v = vorticity.
! F and dF are used for work and should be sized to a least a scalar pencil.

use grid
use partition_data
use dof_indices
use basic_state
implicit none
! Arguments:
logical :: perturbation, compute_curl_rho_u
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in)  :: q
real(8), dimension(sr:er, sphi:ephi, nz, 4),    intent(out) :: v

! These are used as work arrays for the thing being differentiated and its derivative:
real(8), dimension(sr:er, sphi:ephi, nz) :: F, dF

! Locals:
integer :: iz, ir, iphi
integer, parameter :: rphi_comp = 3, dilat = 4

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: add z_derivatives_to_vort has been called'
   print *, ' perturbation = ', perturbation, ' compute_curl_rho_u = ', compute_curl_rho_u
end if
#endif

! -d/dz (r uphi) ---> r-component of vorticity:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (.not. perturbation) then
   if (compute_curl_rho_u) then      
      F(:,:,:) = q(:,:,:,amom)  ! r*uphi*rho
   else 
      F(:,:,:) = q(:,:,:,amom) / q(:,:,:,irho) ! r*uphi
   end if
else if (perturbation) then
   if (compute_curl_rho_u) then      
      do iz = 1, nz
         do iphi = sphi, ephi
            do ir = sr, er      
               F(ir,iphi,iz) = q(ir,iphi,iz,amom) - &
                    rho_basic(ir, iz)*uphi_basic(ir, iz)*rgrid(ir)
            end do
         end do
      end do
   else
      do iz = 1, nz
         do iphi = sphi, ephi
            do ir = sr, er      
               F(ir,iphi,iz) = q(ir,iphi,iz,amom) / q(ir,iphi,iz,irho) - &
                    uphi_basic(ir, iz)*rgrid(ir)
            end do
         end do
      end do
   end if
end if

call pade_diff_z(mr*mphi, F, dF)
v(:,:,:,r_comp) = v(:,:,:,r_comp) - dF(:,:,:) 

! d/dz ur ---> rphi component of vorticity:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (compute_curl_rho_u) then
   F(:,:,:) = q(:,:,:,rmom)                 ! rho*ur   
else
   F(:,:,:) = q(:,:,:,rmom) / q(:,:,:,irho) ! ur
end if

call pade_diff_z(mr*mphi, F, dF)
v(:,:,:,rphi_comp) = v(:,:,:,rphi_comp) + dF(:,:,:)

! d/dz uz ---> dilatation:
! ~~~~~~~~~~~~~~~~~~~~~~~~
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! uz
         F(ir,iphi,iz) = q(ir,iphi,iz,zmom) / q(ir,iphi,iz,irho)
      end do
   end do
end do
call pade_diff_z(mr*mphi, F, dF)
v(:,:,:,dilat) = v(:,:,:,dilat) + dF(:,:,:)

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0: about to return from add z_derivatives_to_vort'
end if
#endif

end subroutine add_z_derivatives_to_vort_and_dil

!----------------------------------------------------------------------------------85

subroutine add_r_derivatives_to_vort_and_dil(perturbation, compute_curl_rho_u, &
                                             q, v, F, dF)

! q and v must be in r-space.  v = vorticity and dilatation.

! F and dF are work arrays used for the quantity being differentiated and its
! derivative.  They should each have at least the size of a scalar pencil.

use grid
use partition_data
use dof_indices
use basic_state
implicit none
! Arguments:
real(8), dimension(sphi:ephi, sz_r:ez_r, 4, nr), intent(in)  :: q
real(8), dimension(sphi:ephi, sz_r:ez_r, 4, nr), intent(out) :: v
logical :: perturbation, compute_curl_rho_u

! These are used as work arrays for the thing being differentiated and its derivative:
real(8), dimension(sphi:ephi, sz_r:ez_r, nr) :: F, dF

! Locals:
integer :: iz, ir, iphi
integer, parameter :: rphi_comp = 3 ! phi_hat times r
integer, parameter :: dilat = 4

! d/dr (r uphi) ---> z-component.  Need to divide by r at the end:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (.not. perturbation) then
   if (compute_curl_rho_u) then
      do ir = 1, nr
         do iz = sz_r, ez_r
            do iphi = sphi, ephi
               F(iphi,iz,ir) = q(iphi,iz,amom,ir)
            end do
         end do
      end do      
   else
      do ir = 1, nr
         do iz = sz_r, ez_r
            do iphi = sphi, ephi
               F(iphi,iz,ir) = q(iphi,iz,amom,ir) / q(iphi,iz,irho,ir)
            end do
         end do
      end do
   end if
else if (perturbation) then
   if (compute_curl_rho_u) then
      do ir = 1, nr
         do iz = sz_r, ez_r
            do iphi = sphi, ephi
               F(iphi,iz,ir) = q(iphi,iz,amom,ir) - &
                    rho_basic_r_space(iz,ir)*uphi_basic_r_space(iz,ir)*rgrid(ir)
               !F(iphi,iz,ir) = q(iphi,iz,amom,ir) - &
               !  rho_basic(ir,iz)*uphi_basic(ir,iz)*rgrid(ir)               
            end do
         end do
      end do
   else
      do ir = 1, nr
         do iz = sz_r, ez_r
            do iphi = sphi, ephi
               F(iphi,iz,ir) = q(iphi,iz,amom,ir) / q(iphi,iz,irho,ir) - &
                    uphi_basic_r_space(iz, ir)*rgrid(ir)
               !F(iphi,iz,ir) = q(iphi,iz,amom,ir) / q(iphi,iz,irho,ir) - &
               !                uphi_basic(ir,iz)*rgrid(ir)
               
            end do
         end do
      end do
   end if
end if

call pade_diff_bundle(mphi*mz_r, nr, Ji_r, F, dF)
do ir = 1, nr
   do iz = sz_r, ez_r
      do iphi = sphi, ephi
          v(iphi,iz,z_comp,ir) = v(iphi,iz,z_comp,ir) + dF(iphi,iz,ir)
      end do
   end do
end do   

! -d/dr uz ---> rphi component:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (compute_curl_rho_u) then
   do ir = 1, nr
      do iz = sz_r, ez_r
         do iphi = sphi, ephi
            F(iphi,iz,ir) = q(iphi,iz,zmom,ir)
         end do
      end do
   end do   
else
   do ir = 1, nr
      do iz = sz_r, ez_r
         do iphi = sphi, ephi
            F(iphi,iz,ir) = q(iphi,iz,zmom,ir) / q(iphi,iz,irho,ir)
         end do
      end do
   end do
end if

call pade_diff_bundle(mphi*mz_r, nr, Ji_r, F, dF)
do ir = 1, nr
   do iz = sz_r, ez_r
      do iphi = sphi, ephi
         v(iphi,iz,rphi_comp,ir) = v(iphi,iz,rphi_comp,ir) - dF(iphi,iz,ir)
      end do
   end do
end do

! 1/r ddr r*ur ---> dilatation:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
do ir = 1, nr
   do iz = sz_r, ez_r
      do iphi = sphi, ephi
         F(iphi,iz,ir) = q(iphi,iz,rmom,ir) / q(iphi,iz,irho,ir) * rgrid(ir)
      end do
   end do
end do
call pade_diff_bundle(mphi*mz_r, nr, Ji_r, F, dF)
do ir = 1, nr
   do iz = sz_r, ez_r
      do iphi = sphi, ephi
         v(iphi,iz,dilat,ir) = v(iphi,iz,dilat,ir) + dF(iphi,iz,ir) / rgrid(ir)
      end do
   end do
end do

end subroutine add_r_derivatives_to_vort_and_dil

!----------------------------------------------------------------------------------85

subroutine add_phi_derivatives_to_vort_and_dil (compute_curl_rho_u, q, v, F, dF)

! q and v must be in phi-space.  v = vorticity and dilatation.

! F and dF are work arrays used for the quantity being differentiated and its
! derivative.  They need only have size >= a scalar penci.,

use grid
use partition_data
use dof_indices
use basic_state
implicit none
! Arguments:
logical :: compute_curl_rho_u
real(8), dimension(sr:er, sz_phi:ez_phi, 4, nphi), intent(in)  :: q
real(8), dimension(sr:er, sz_phi:ez_phi, 4, nphi), intent(out) :: v

! These are used as work arrays for the thing being differentiated and its derivative:
real(8), dimension(sr:er, sz_phi:ez_phi, nphi) :: F, dF

! Locals:
integer :: iz, ir, iphi
integer, parameter :: rphi_comp = 3 ! phi_hat times r
integer, parameter :: dilat = 4

! -d/dphi ur ---> z-comp.  Need to divide by "r" at end:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (compute_curl_rho_u) then
   do iphi = 1, nphi
      do iz = sz_phi, ez_phi
         do ir = sr, er
            F(ir,iz,iphi) = q(ir,iz,rmom,iphi) ! rho*ur
         end do
      end do
   end do   
else
   do iphi = 1, nphi
      do iz = sz_phi, ez_phi
         do ir = sr, er
            F(ir,iz,iphi) = q(ir,iz,rmom,iphi) / q(ir,iz,irho,iphi) ! ur
         end do
      end do
   end do
end if

call pade_diff_periodic(mr*mz_phi, nphi, dphi, F, dF)

do iphi = 1, nphi
   do iz = sz_phi, ez_phi
      do ir = sr, er
         v(ir,iz,z_comp,iphi) = v(ir,iz,z_comp,iphi) - dF(ir,iz,iphi)
      end do
   end do
end do

! d/dphi uz ---> r-comp:
!~~~~~~~~~~~~~~~~~~~~~~~
if (compute_curl_rho_u) then
   do iphi = 1, nphi
      do iz = sz_phi, ez_phi
         do ir = sr, er
            F(ir,iz,iphi) = q(ir,iz,zmom,iphi) ! rho*uz
         end do
      end do
   end do   
else
   do iphi = 1, nphi
      do iz = sz_phi, ez_phi
         do ir = sr, er
            F(ir,iz,iphi) = q(ir,iz,zmom,iphi) / q(ir,iz,irho,iphi) ! uz
         end do
      end do
   end do
end if

call pade_diff_periodic(mr*mz_phi, nphi, dphi, F, dF)

do iphi = 1, nphi
   do iz = sz_phi, ez_phi
      do ir = sr, er
         v(ir,iz,r_comp,iphi) = v(ir,iz,r_comp,iphi) + dF(ir,iz,iphi)
      end do
   end do
end do

! 1/r ddphi uphi ---> dilatation: 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
do iphi = 1, nphi
   do iz = sz_phi, ez_phi
      do ir = sr, er
         ! uphi
         F(ir,iz,iphi) = q(ir,iz,amom,iphi) / q(ir,iz,irho,iphi) / rgrid(ir)
      end do
   end do
end do

call pade_diff_periodic(mr*mz_phi, nphi, dphi, F, dF)

do iphi = 1, nphi
   do iz = sz_phi, ez_phi
      do ir = sr, er
         v(ir,iz,dilat,iphi) = v(ir,iz,dilat,iphi) + dF(ir,iz,iphi)/rgrid(ir)
      end do
   end do
end do

end subroutine add_phi_derivatives_to_vort_and_dil

!----------------------------------------------------------------------------------85

