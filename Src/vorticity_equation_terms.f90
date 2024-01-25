!----------------------------------------------------------------------------------85

! As a diagnostic, compute and output for plotting, the terms in the vorticity
! transport equation.

!----------------------------------------------------------------------------------85

module indices_for_vorticity_equation_terms

integer, parameter :: &
z_stretch       = 1, &
z_tilt_from_r   = 2, &
z_tilt_from_phi = 3, &
z_dilatation    = 4, &
z_baroclinic    = 5, &

r_stretch       = 6, &
r_tilt_from_z   = 7, &
r_tilt_from_phi = 8, &
r_dilatation    = 9, &
r_baroclinic    = 10, &
r_geometry      = 11, &
 
phi_stretch     = 12, &
phi_tilt_from_z = 13, &
phi_tilt_from_r = 14, &
phi_dilatation  = 15, &
phi_baroclinic  = 16, &
phi_geometry    = 17, &

z_total   = 18, &
r_total   = 19, &
phi_total = 20, &

n_vort_terms = 20

end module indices_for_vorticity_equation_terms

!----------------------------------------------------------------------------------85

subroutine tecplot_vort_eq_terms_in_merid_horiz_planes(q, iphi_plot, iz_plot, t, istep)

! iphi_plot = desired meridional plane.  iz_plot = desired horizontal plane.
! If iphi_plot = 0 we will skip the meridional plane.  If iz_plot = 0 we will skip
! the horizontal plane.

use grid
use xy_coordinates_of_grid
use logical_units, only: lun_tecplot
use partition_data
use indices_for_vorticity_equation_terms
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in) :: q
integer,                                        intent(in) :: iphi_plot, iz_plot, istep
real(8),                                        intent(in) :: t

! Locals:
real(8), dimension(sr:er, sphi:ephi, nz, n_vort_terms) :: terms
character(80) :: filename
integer :: ir, iz, iphi, irec, nz_plot, iphi_rec, iterm
logical :: i_have_iphi_plot
character(138) :: list

#ifdef debug_print
if (my_node .eq. 0) print *,' node 0: tecplot_vort_eq_terms_in_merid_horiz_planes has been called'
#endif

list = '"r","z","zs","zt_f_r","zt_r_p","zd","zb","rs","rt_f_z","rt_f_phi",'// &
       '"rd","rb","rg","ps","pt_f_z","pt_f_r","pd","pb","pg"' // &
       '"ztot","rtot","ptot"' 
   
call vorticity_equation_terms(q, terms)   

! Output meridional plane:
! ~~~~~~~~~~~~~~~~~~~~~~~~
if (iphi_plot .ne. 0) then
   if (suppress_z_derivatives_when_nz_not_1) then
      nz_plot = 1
   else
      nz_plot = nz
   end if

   write (filename, 1) iphi_plot, int(t), t - int(t), istep
1  format ('vort_eq_terms_iphi_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
   open (unit = lun_tecplot, file = filename, form = 'formatted', &
         status = 'unknown', access = 'direct',recl = 13*(n_vort_terms+2) + 1)
   
   if (my_node .eq. 0) then
      ! char(10) = new line; char(13) = cr
      write (lun_tecplot, 2, rec = 1) t, char(13)
      2  format ('TITLE = "t = ', e12.5, '"', a1)
      write (lun_tecplot, 3, rec = 2) list, char(13)
      3  format('VARIABLES = ', a138, a1)
      write (lun_tecplot, 4, rec = 3) nr, nz_plot, char(13)
      4 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', a1)
   end if

   ! Write only of this processor has the iphi:
   i_have_iphi_plot = (iphi_plot .ge. sphi) .and. (iphi_plot .le. ephi)
   if (i_have_iphi_plot) then
      do iz = 1, nz_plot
         do ir = sr, er
            irec = (iz - 1)*nr + ir + 3
            write(lun_tecplot, "(22(1x, e12.5), a1)", rec = irec) &
                 rgrid(ir), zgrid(iz), (terms(ir,iphi_plot,iz,iterm), iterm = 1, n_vort_terms), &
                 char(10)      
         end do
      end do
   end if
   close(lun_tecplot)
   if (my_node .eq. 0) print *, ' node 0: wrote tecplot file =', filename
end if

! Output horizontal plane if desired and for a non-axisymmetric run:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ((iz_plot .ne. 0) .and. (nphi .ne. 1)) then
   write (filename, 5) iz_plot, int(t), t - int(t), istep
   5 format ('vort_eq_terms_iz_', i4.4, '_t_', i6.6, f0.4, '_', i7.7, '.tec')
   
   open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
        access = 'direct', recl = 13*(n_vort_terms+2) + 1)

   if (my_node .eq. 0) then
      ! char(10) = new line
      write (lun_tecplot, 2, rec = 1) t,          char(10)
      write (lun_tecplot, 3, rec = 2) list,       char(10)
      write (lun_tecplot, 4, rec = 3) nr, nphi+1, char(10)
   end if

   do iphi = sphi, ephi
      do ir = sr, er
         irec = (iphi - 1)*nr + ir + 3
         write(lun_tecplot, "(22(1x, e12.5), a1)", rec = irec) &
              xgrid(ir,iphi), ygrid(ir,iphi), &
              (terms(ir, iphi, iz_plot, iterm), iterm = 1, n_vort_terms), char(10)  
      end do
   end do

   ! Periodic completion in phi:
   if (sphi .eq. 1) then
      ! The data for iphi = 1, gets put in the records for nphi + 1
      iphi_rec = nphi + 1
      iphi     = 1
      do ir = sr, er
         irec = (iphi_rec - 1)*nr + ir + 3

         ! Note: I have stored the coordinates for this line in the iphi = 0 slot
         ! since the flow data sits in the same processor:
         write(lun_tecplot, "(22(1x, e12.5), a1)", rec = irec) &
                 xgrid(ir,0), ygrid(ir,0), &
                 (terms(ir, iphi, iz_plot, iterm), iterm = 1, n_vort_terms), char(10)      
      end do         
   end if
   
   close(lun_tecplot)
   if (my_node .eq. 0) print *, ' node 0: wrote tecplot file =', filename
end if

end subroutine tecplot_vort_eq_terms_in_merid_horiz_planes

!----------------------------------------------------------------------------------85

subroutine vorticity_equation_terms(q, terms)

use grid
use partition_data
! This is work space needed for subroutine vort_and_dil.
use transposes_of_q_and_qdot, only: q_r_space, q_phi_space, qdot_r_space, &
     qdot_phi_space
use dof_indices
use indices_for_vorticity_equation_terms
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof),         intent(in)   :: q
real(8), dimension(sr:er, sphi:ephi, nz, n_vort_terms), intent(out)  :: terms

! Locals:
! ~~~~~~~
integer :: ir, iz, iphi

! Some arguments to vort_and_dil:
logical, parameter :: perturbation       = .false.
logical, parameter :: compute_curl_rho_u = .false.
logical, parameter :: divide_by_rho      = .false.

! Note: These are passed as work
! to vort_and_dil which further passes them down to the derivative routines where
! they will be dimensioned differently but have the same size.
real(8), dimension(sr:er, sphi:ephi, nz) :: F, dF

! Vorticity and dilatation (which is in the fourth slot).
real(8), dimension(sr:er, sphi:ephi, nz, 4) :: v

! Primitive variable and their derivatives:
real(8), dimension(sr:er, sphi:ephi, nz) :: &
                                       uz,              ur,              uphi,     &
     dp_dz,          drho_dz,          duz_dz,          dur_dz,          duphi_dz, &
     dp_dr,          drho_dr,          duz_dr,          dur_dr,          duphi_dr, &
     dp_dphi_over_r, drho_dphi_over_r, duz_dphi_over_r, dur_dphi_over_r, duphi_dphi_over_r

! First executable:
call vort_and_dil(perturbation, divide_by_rho, compute_curl_rho_u, &
     q, v, q_r_space, q_phi_space, qdot_r_space, qdot_phi_space, F, dF)

call prim_vars_and_derivatives(q,      uz,              ur,              uphi,     &
     dp_dz,          drho_dz,          duz_dz,          dur_dz,          duphi_dz, &
     dp_dr,          drho_dr,          duz_dr,          dur_dr,          duphi_dr, &
     dp_dphi_over_r, drho_dphi_over_r, duz_dphi_over_r, dur_dphi_over_r, duphi_dphi_over_r)

! Geometry terms:
! ~~~~~~~~~~~~~~~
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         terms(ir,iphi,iz,phi_geometry) =   v(ir,iphi,iz,phi_comp)*ur  (ir,iphi,iz)/rgrid(ir)
         terms(ir,iphi,iz,r_geometry  ) = - v(ir,iphi,iz,phi_comp)*uphi(ir,iphi,iz)/rgrid(ir)     
      end do
   end do
end do

terms(:,:,:,z_stretch      ) =  v(:,:,:,z_comp)   * duz_dz         (:,:,:)
terms(:,:,:,z_tilt_from_r  ) =  v(:,:,:,r_comp)   * duz_dr         (:,:,:)
terms(:,:,:,z_tilt_from_phi) =  v(:,:,:,phi_comp) * duz_dphi_over_r(:,:,:)
terms(:,:,:,z_dilatation   ) = -v(:,:,:,z_comp  ) * v(:,:,:,4)
terms(:,:,:,z_baroclinic   ) = (drho_dr(:,:,:)*dp_dphi_over_r  (:,:,:) - &
     dp_dr  (:,:,:)*drho_dphi_over_r(:,:,:)) / q(:,:,:,irho)**2
terms(:,:,:,z_total) = terms(:,:,:,z_stretch)       + terms(:,:,:,z_tilt_from_r) + &
                       terms(:,:,:,z_tilt_from_phi) + terms(:,:,:,z_dilatation)  + &
                       terms(:,:,:,z_baroclinic)

terms(:,:,:,r_stretch)       =  v(:,:,:,r_comp)   * dur_dr        (:,:,:)
terms(:,:,:,r_tilt_from_z)   =  v(:,:,:,z_comp)   * dur_dz        (:,:,:)
terms(:,:,:,r_tilt_from_phi) =  v(:,:,:,phi_comp) * dur_dphi_over_r(:,:,:) 
terms(:,:,:,r_dilatation)    = -v(:,:,:,r_comp)   * v(:,:,:,4)
terms(:,:,:,r_baroclinic)    = (drho_dphi_over_r(:,:,:)*dp_dz         (:,:,:) - &
     drho_dz         (:,:,:)*dp_dphi_over_r(:,:,:)) / q(:,:,:,irho)**2
terms(:,:,:,r_total) = terms(:,:,:,r_stretch)       + terms(:,:,:,r_tilt_from_z) + &
                       terms(:,:,:,r_tilt_from_phi) + terms(:,:,:,r_dilatation)  + &
                       terms(:,:,:,r_baroclinic)    + terms(:,:,:,r_geometry)

terms(:,:,:,phi_stretch)     =  v(:,:,:,phi_comp) * duphi_dphi_over_r(:,:,:)
terms(:,:,:,phi_tilt_from_z) =  v(:,:,:,z_comp  ) * duphi_dz         (:,:,:)
terms(:,:,:,phi_tilt_from_r) =  v(:,:,:,r_comp  ) * duphi_dr         (:,:,:)
terms(:,:,:,phi_dilatation)  = -v(:,:,:,phi_comp) * v(:,:,:,4)
terms(:,:,:,phi_baroclinic)  = (drho_dz(:,:,:) *  dp_dr(:,:,:) - &
     dp_dz  (:,:,:) *drho_dr(:,:,:)) / q(:,:,:,irho)**2
terms(:,:,:,phi_total) = terms(:,:,:,phi_stretch)     + terms(:,:,:,phi_tilt_from_z) + &
                         terms(:,:,:,phi_tilt_from_r) + terms(:,:,:,phi_dilatation)  + &
                         terms(:,:,:,phi_baroclinic)  + terms(:,:,:,phi_geometry)

end subroutine vorticity_equation_terms

!----------------------------------------------------------------------------------85

subroutine prim_vars_and_derivatives(q,      uz,              ur,              uphi,     &
           dp_dz,          drho_dz,          duz_dz,          dur_dz,          duphi_dz, &
           dp_dr,          drho_dr,          duz_dr,          dur_dr,          duphi_dr, &
           dp_dphi_over_r, drho_dphi_over_r, duz_dphi_over_r, dur_dphi_over_r, duphi_dphi_over_r)

! This subroutine is for computing terms in the vorticity equation for diagnostic
! purposes and therefore does not apply the Fargo transformation.

use partition_data
use dof_indices
use grid
use thermal_parameters

implicit none
real(8), intent(in),   dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8), intent(out),  dimension(sr:er, sphi:ephi, nz) :: &
                                        uz,              ur,              uphi,     &
      dp_dz,          drho_dz,          duz_dz,          dur_dz,          duphi_dz, &
      dp_dr,          drho_dr,          duz_dr,          dur_dr,          duphi_dr, &
      dp_dphi_over_r, drho_dphi_over_r, duz_dphi_over_r, dur_dphi_over_r, duphi_dphi_over_r

! These are dimensioned in z-space but will have data in other spaces since pencils in
! all the spaces have the same size:
real(8), dimension(sr:er, sphi:ephi, nz) :: transposed, deriv

! Pressure is local:
real(8), dimension(sr:er, sphi:ephi, nz) :: p

integer :: ir, iphi, iz

! Velocities:
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ur  (ir,iphi,iz) = q(ir,iphi,iz,rmom) / q(ir,iphi,iz,irho)
         uz  (ir,iphi,iz) = q(ir,iphi,iz,zmom) / q(ir,iphi,iz,irho)
         uphi(ir,iphi,iz) = q(ir,iphi,iz,amom) / q(ir,iphi,iz,irho) / rgrid(ir)
      end do
   end do
end do

! Pressure:
if (isothermal) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            p(ir,iphi,iz) = q(ir, iphi, iz, irho) * ci_squared_initial(ir, iz)
         end do
      end do
   end do
else
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            p(ir,iphi,iz) = q(ir, iphi, iz, ener) * gm1
         end do
      end do
   end do   
end if

! --------------
! r derivatives:
! --------------
if (nr .ne. 1) then
   call transpose_z_to_r (1, uz, transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, duz_dr)

   call transpose_z_to_r(1, ur, transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, dur_dr)

   call transpose_z_to_r(1, uphi, transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, duphi_dr)

   call transpose_z_to_r(1, p, transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, dp_dr)

   call transpose_z_to_r(1, q(sr,sphi,1,irho), transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, drho_dr)
else
   duz_dr   = 0.d0
   dur_dr   = 0.d0
   duphi_dr = 0.d0
   dp_dr    = 0.d0
   drho_dr  = 0.d0
end if

! ----------------
! phi derivatives:
! ----------------
if (nphi .ne. 1) then
   ! duz/dphi:
   call transpose_z_to_phi(1, uz, transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z (1, deriv, duz_dphi_over_r)

   ! dur/dphi
   call transpose_z_to_phi(1, ur, transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z(1, deriv, dur_dphi_over_r)

   ! duphi/dphi
   call transpose_z_to_phi(1, uphi, transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z (1, deriv, duphi_dphi_over_r)

   ! dp/dr
   call transpose_z_to_phi(1, p, transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z (1, deriv, dp_dphi_over_r)
   
   ! drho_/dphi:
   call transpose_z_to_phi(1, q(sr,sphi,1,irho), transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z (1, deriv, drho_dphi_over_r)

   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            duz_dphi_over_r  (ir,iphi,iz) = duz_dphi_over_r  (ir,iphi,iz) / rgrid(ir)
            dur_dphi_over_r  (ir,iphi,iz) = dur_dphi_over_r  (ir,iphi,iz) / rgrid(ir)
            duphi_dphi_over_r(ir,iphi,iz) = duphi_dphi_over_r(ir,iphi,iz) / rgrid(ir)
            dp_dphi_over_r   (ir,iphi,iz) = dp_dphi_over_r   (ir,iphi,iz) / rgrid(ir)
            drho_dphi_over_r (ir,iphi,iz) = drho_dphi_over_r (ir,iphi,iz) / rgrid(ir)
         end do
      end do
   end do   
else
   duz_dphi_over_r   = 0.d0
   dur_dphi_over_r   = 0.d0
   duphi_dphi_over_r = 0.d0
   dp_dphi_over_r    = 0.d0
   drho_dphi_over_r  = 0.d0
end if

! -------------------------
! Done with phi derivatives
! -------------------------

! --------------
! z derivatives:
! --------------
if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
   call pade_diff_z(nbundle_z, uz,   duz_dz  )
   call pade_diff_z(nbundle_z, ur,   dur_dz  ) 
   call pade_diff_z(nbundle_z, uphi, duphi_dz)
   call pade_diff_z(nbundle_z, p,    dp_dz   )   
   call pade_diff_z(nbundle_z, q(sr,sphi,1,irho), drho_dz)   
else
   duz_dz   = 0.d0
   dur_dz   = 0.d0
   duphi_dz = 0.d0
   dp_dz    = 0.d0
   drho_dz  = 0.d0
end if

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: prim_vars_and_derivatives.  Returning'
#endif

end subroutine prim_vars_and_derivatives

!----------------------------------------------------------------------------------85
