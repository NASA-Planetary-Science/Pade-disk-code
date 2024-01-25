!----------------------------------------------------------------------------------85

subroutine tecplot_baroclinic_term (q, t)

use grid
use xy_coordinates_of_grid
use logical_units, only: lun_tecplot
use partition_data

implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: t

! Local:
real(8), dimension(sr:er, sphi:ephi, nz, 3) :: b ! baroclinic term

character(80) :: filename
integer :: ir, iphi, iz, irec, icomp, iphi_rec, nz_plot

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: tecplot_baroclinic_term has been called'
#endif

if (suppress_z_derivatives_when_nz_not_1) then
   nz_plot = 1
else
   nz_plot = nz
end if

call baroclinic_term (q, b)

write (filename, 5) int(t), t - int(t)
5 format ('baroclinic_term', i6.6, f0.4, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
     access = 'direct', recl = 13*7 + 1) ! 3 coordinates + 3 components + node

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
1  format (' TITLE = "t = ', e12.5, '"', a1, a1)

   write (lun_tecplot, 3, rec = 2) char(13), char(10)
3  format(' VARIABLES = "x", "y", "z", "bz", "br", "bphi", "node"', a1, a1)
   
   write (lun_tecplot, 4, rec = 3) nr, nphi+1, nz_plot, char(13), char(10)
4  format(' ZONE I=',i4,',',' J=',i4,',',' K=',i4,',',' DATAPACKING=POINT', a1, a1)
end if

do iz = 1, nz_plot
   do iphi = sphi, ephi
      do ir = sr, er
         ! nphi+1 leaves room for periodic completion in phi:
         irec = ir + (iphi-1)*nr + (iz-1)*nr*(nphi+1) + 3  
         write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) &
              xgrid(ir,iphi), ygrid(ir,iphi), zgrid(iz), &
              (b(ir, iphi, iz, icomp), icomp = 1, 3), float(my_node), char(10)  ! z, r, phi components.
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
         write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) &
              xgrid(ir,0), ygrid(ir,0), zgrid(iz), &
              (b(ir, iphi, iz, icomp), icomp = 1, 3), float(my_node), char(10)  ! z, r, phi components.
      end do
   end do
end if
close(lun_tecplot)

if (my_node .eq. 0) print *, ' node 0: Wrote baroclinic tecplot file ', filename

end subroutine tecplot_baroclinic_term

!----------------------------------------------------------------------------------85

subroutine baroclinic_term (q, b)

! b = Baroclinic term = (grad rho X grad p) / rho^3.  The 3 components are ordered as
! (z, r, and phi).

use grid
use partition_data
use dof_indices, only: irho, zmom, rmom, amom, ener
use thermal_parameters

implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8), dimension(sr:er, sphi:ephi, nz, 3)    :: b

! Local:
! ~~~~~~
integer, parameter :: nvars = 2
integer, parameter :: z_component = 1, r_component = 2, phi_component = 3
integer, parameter :: rho = 1, p = 2 ! convenient indices
real(8), dimension(sr:er, sphi:ephi, nz, nvars)      :: vars,           ddz_vars, ddr_vars, ddphi_vars
real(8), dimension(sphi:ephi, sz_r:ez_r, nvars, nr)  :: vars_r_space,   ddr_vars_r_space
real(8), dimension(sr:er, sz_phi:ez_phi, ndof, nphi) :: vars_phi_space, ddphi_vars_phi_space
real(8) :: rho_cubed
logical :: have_z_derivatives
integer :: iz, iphi, ir

have_z_derivatives = (nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)

! Calculate rho and p and put them into vars:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         vars(ir, iphi, iz, rho) = q(ir, iphi, iz, irho)
      end do
   end do
end do

if (isothermal) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            vars(ir, iphi, iz, p) = q(ir, iphi, iz, irho) * ci_squared_initial(ir, iz)
         end do
      end do
   end do
else
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            vars(ir, iphi, iz, p) = q(ir, iphi, iz, ener) * gm1
         end do
      end do
   end do
end if

! Now take derivatives of rho and p:
if (have_z_derivatives) then
   call pade_diff_z(nbundle_z, vars(1, 1, 1, 1), ddz_vars(1, 1, 1, 1))
   call pade_diff_z(nbundle_z, vars(1, 1, 1, 2), ddz_vars(1, 1, 1, 2))
end if

if (nr .ne. 1) then
   call transpose_z_to_r    (nvars, vars,   vars_r_space)
   call pade_diff_bundle(mphi*mz_r*nvars, nr, Ji_r, vars_r_space, ddr_vars_r_space)
   call transpose_r_to_z    (nvars, ddr_vars_r_space, ddr_vars)   
end if

if (nphi .ne. 1) then
   call transpose_z_to_phi (nvars, vars, vars_phi_space)
   call pade_diff_periodic(mr*mz_phi*nvars, nphi, dphi, vars_phi_space, ddphi_vars_phi_space)
   call transpose_phi_to_z (nvars, ddphi_vars_phi_space, ddphi_vars)   
end if

! Finally the 3 components of the baroclinic term:

! z-component:
! ~~~~~~~~~~~~
if ((nr .ne. 1) .and. (nphi .ne. 1)) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            rho_cubed = q(ir, iphi, iz, irho)**3
            b(ir, iphi, iz, z_component) = (ddr_vars  (ir, iphi, iz, rho)*  &
                                            ddphi_vars(ir, iphi, iz, p  ) - &
                                            ddr_vars  (ir, iphi, iz, p  )*  &
                                            ddphi_vars(ir, iphi, iz, rho))/ &
                                            rgrid(ir) / rho_cubed
         end do
      end do
   end do   
else
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            b(ir, iphi, iz, z_component) = 0.0d0
         end do
      end do
   end do   
end if

! r-component:
! ~~~~~~~~~~~~
if (have_z_derivatives .and. (nphi .ne. 1)) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            rho_cubed = q(ir, iphi, iz, irho)**3
            b(ir, iphi, iz, r_component) = (ddz_vars  (ir, iphi, iz, p  )*  &
                                            ddphi_vars(ir, iphi, iz, rho) - &
                                            ddz_vars  (ir, iphi, iz, rho)*  &
                                            ddphi_vars(ir, iphi, iz, p  ))/ &
                                            rgrid(ir) / rho_cubed
         end do
      end do
   end do   
else
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            b(ir, iphi, iz, r_component) = 0.0d0
         end do
      end do
   end do   
end if

! phi-component:
! ~~~~~~~~~~~~~~
if (have_z_derivatives .and. (nr .ne. 1)) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            rho_cubed = q(ir, iphi, iz, irho)**3
            b(ir, iphi, iz, phi_component) = (ddz_vars(ir, iphi, iz, rho)*  &
                                              ddr_vars(ir, iphi, iz, p  ) - &
                                              ddz_vars(ir, iphi, iz, p  )*  &
                                              ddr_vars(ir, iphi, iz, rho))/ &
                                              rho_cubed
         end do
      end do
   end do   
else
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            b(ir, iphi, iz, phi_component) = 0.0d0
         end do
      end do
   end do   
end if
   
end subroutine baroclinic_term

!----------------------------------------------------------------------------------85

subroutine z_derivative_of_vars(nvars, vars, ddz_vars)

! Keep this just in case the above z derivatives don't work.

use grid
use partition_data
implicit none
integer :: nvars
real(8), dimension(sr:er, sphi:ephi, nz, nvars) :: vars, ddz_vars

! Locals:
real(8), dimension(sr:er, sphi:ephi, nz) :: F, Fd
integer :: ivar, iz, iphi, ir

do ivar = 1, 2
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            F(ir, iphi, iz) = vars(ir, iphi, iz, ivar)
         end do
      end do
   end do


   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            ddz_vars(ir, iphi, iz, ivar) = Fd(ir, iphi, iz)
         end do
      end do
   end do
end do

end subroutine z_derivative_of_vars

!----------------------------------------------------------------------------------85

