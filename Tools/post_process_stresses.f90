!*****************************************************************************80

implicit none

integer :: nvar, isel
! dataM is data in the middle of the cell:
real(8), allocatable, dimension(:, :, :) :: data, dataW, dataE, dataN, dataS, dataM

real(8) :: rN, rS, rE, rW
real(8) :: dr, dz, rcen, zcen
real(8) :: ratio, at_north, at_south, ur_tilde_theory, ur_tilde_theory2, ur_tilde_theory3, &
     ur_tilde_theory4

! Derivatives:
real(8) :: ddr_r2_Tpr, ddz_r2_Tpz, ddz_upt, ddr_upt_r
real(8) :: upt_r2, duz_dz

! Variable indices for file with stresses:
integer :: iTzz, iTpp, iTrr, iTrz, iTpr, iTpz, iTtrace, irhob, iurt, iuzt, &
           iupt, r_var, z_var

! For file with velocities, etc.
integer :: irho, iur, iuz

integer :: lun_tecplot_in = 1, lun_tecplot_out = 2, lunout = 3, lun_tecplot_out2 = 4
integer :: ir, iz, nr, nz, ivar
real(8) :: uK, right, left, rave, deriv
real(8) :: H0_over_r0

real(8), allocatable, dimension(:) :: Tpr_mid, ur_tilde_mid, rho_bar_mid

! Quantities we will read (phi Reynolds averages at an instant):
real(8) :: bar_read(10)

! phi-time Reynolds averages:
real(8), allocatable, dimension(:,:,:) :: bar

real(8), allocatable, dimension(:) :: rgrid, zgrid

! Final phi-time averages we will output:
!real(8), allocatable, dimension(:, :) :: rho_bar, ur_tilde, uz_tilde, uphi_tilde

! For midplane determination:
integer :: iz_mid, iz_mid1, iz_mid2
logical :: nz_is_odd

real(8) :: pi, T_orbital0, H0, r0, Omega0, GM
real(8) :: ur_tilde, mass_flux, mass_flux_theory, mass_flux_theory2, mass_flux_theory3, &
     mass_flux_theory4

real(8) :: d_rur, rinv_ddr_r_ur, duz_dr
real(8) :: r, rmin, rmax, q, exponent, H

character(80) :: filename

! This is to avoid "uninitialized warnings":
iz_mid1 = 1
iz_mid2 = 2

pi = 4.0*ATAN(1.0d0)
! Orbital period at mid-radius of the domain:   
T_orbital0 = 1.0d0
H0 = 1.d0

! This is a parameter of the simulation:
H0_over_r0 = 0.10d0
r0 = H0 / H0_over_r0

! This is also a paraneter of the simulation:
q = -1.d0

! Consequence of above:
Omega0 = 2.d0 * pi / T_orbital0
GM     = Omega0**2 * r0**3

print *, ' r0 = ', r0
print *, ' Enter anything to continue'
read(5, *)

write(6, 1)
1 format(' Enter',/, &
         ' 1: For stresses file',/, &
         ' 2: For merid velocity perturbation file',/, &
         ' 3: To output H(r)',/, &
         ' --->', $)
read(5, *) isel

go to (100, 200, 300) isel

!*****************************************************************************************************

! Post processes stresses:
! ~~~~~~~~~~~~~~~~~~~~~~~~
100 continue

nvar = 13

iTzz = 1; iTpp = 2; iTrr = 3; iTrz = 4; iTpr = 5; iTpz = 6; iTtrace = 7
irhob = 8; iurt = 9; iuzt = 10; iupt = 11; r_var = 12; z_var = 13; irho = 14; iur = 15; iuz = 16

! Read mean density and Favre mean velocities:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print *, ' opening phi_time_Favre_means.tec'
open (unit = lun_tecplot_in, file = 'phi_time_Favre_means.tec', form = 'formatted', status = 'unknown')
read (lun_tecplot_in, *) ! title
read (lun_tecplot_in, *) ! list of variables
! 2 format(' VARIABLES="r","z","rho_bar","ur_tilde","uz_tilde","uphi_tilde", "rho_ur_bar", "rho_uz_bar"')
read (lun_tecplot_in, "(8x, i4, 4x, i4)") nr, nz
!3 format(' ZONE I=',i4,',',' J=',i4,' DATAPACKING=POINT')
print *, ' nr = ', nr, ' nz = ', nz

! Determine iz for the midplane and calculate required quantities at the midplane:
if (mod(nz, 2) .eq. 1) then
   nz_is_odd = .true.
   iz_mid = (1 + nz)/2
   print *, ' iz_mid = ', iz_mid
else
   nz_is_odd = .false.
   iz_mid1 = nz/2
   iz_mid2 = iz_mid1 + 1
   print *, ' iz_mid1 = ', iz_mid1, ' iz_mid2 = ', iz_mid2
end if
print *, ' enter anything to continue'
read(5, *)

! Ready to allocate:
allocate(data (nr, nz, nvar))
allocate(dataW(nr, nz, nvar))
allocate(dataE(nr, nz, nvar))
allocate(dataN(nr, nz, nvar))
allocate(dataS(nr, nz, nvar))
allocate(dataM(nr, nz, nvar)) ! at midpoint of cell

allocate(bar(10, nr, nz))
allocate(Tpr_mid(nr), ur_tilde_mid(nr), rho_bar_mid(nr))

! Read Favre means
do iz = 1, nz
   do ir = 1, nr
      ! 2 format(' VARIABLES="r","z","rho_bar","ur_tilde","uz_tilde","uphi_tilde", "rho_ur_bar", "rho_uz_bar"')
      read(lun_tecplot_in, "(6(1x, e12.5))") &
           data(ir, iz, r_var), data(ir, iz, z_var), &
           data(ir, iz, irhob), data(ir, iz, iurt), data(ir, iz, iuzt), data(ir, iz, iupt)
   end do
end do
close(lun_tecplot_in)

! Read Favre stresses:
open (unit = lun_tecplot_in, file = 'phi_time_Favre_stresses.tec', form = 'formatted', status = 'unknown')
read (lun_tecplot_in, *) ! title
read (lun_tecplot_in, *) ! Variable list
! 12 format(' VARIABLES="r","z","T_phiphi","T_phir","T_phiz","T_rr","T_rz","T_zz", "Trace"')
read (lun_tecplot_in, *) ! nr, nz
do iz = 1, nz
   do ir = 1, nr
      read(lun_tecplot_in, "(9(1x, e12.5))") &
           data(ir, iz, r_var),  data(ir, iz, z_var), &
           data(ir, iz, iTpp), data(ir, iz, iTpr), data(ir, iz, iTpz), data(ir, iz, iTrr), &
           data(ir, iz, iTrz), data(ir, iz, iTzz), data(ir, iz, iTtrace)
   end do
end do
close(lun_tecplot_in)

do iz = 1, nz - 1
   do ir = 1, nr - 1
      ! Obtain data at north, south, east, and west, and center of each cell:
      do ivar = 1, nvar
         dataW(ir, iz, ivar) = 0.5d0*(data(ir,   iz,   ivar) + data(ir,   iz+1, ivar))
         dataE(ir, iz, ivar) = 0.5d0*(data(ir+1, iz,   ivar) + data(ir+1, iz+1, ivar))
         dataN(ir, iz, ivar) = 0.5d0*(data(ir,   iz+1, ivar) + data(ir+1, iz+1, ivar))
         dataS(ir, iz, ivar) = 0.5d0*(data(ir,   iz,   ivar) + data(ir+1, iz,   ivar))

         dataM(ir, iz, ivar) = 0.25d0*(dataW(ir,iz,ivar)+dataE(ir,iz,ivar)+&
                                       dataN(ir,iz,ivar)+dataS(ir,iz,ivar))
      end do
   end do
end do

! Header of output tecplot file:
open (unit = lun_tecplot_out, file = 'mass_flux_theory.tec', form = 'formatted', status = 'unknown')
write (lun_tecplot_out, 2)
2  format ('TITLE = "ur_tilde"')
write (lun_tecplot_out, 3)
3  format('VARIABLES="r","z","mass_flux_theory", "mass_flux_theory3", "mass_flux"')
write (lun_tecplot_out, 4) nr-1, nz-1
4 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT')

do iz = 1, nz - 1
   do ir = 1, nr - 1
      rcen = dataM(ir, iz, r_var)
      zcen = dataM(ir, iz, z_var)

      dz    = dataN(ir,iz,z_var) - dataS(ir,iz,z_var)
      rW    = dataW(ir,iz,r_var)
      rE    = dataE(ir,iz,r_var)
      rN    = dataN(ir,iz,r_var)
      rS    = dataS(ir,iz,r_var)      
      dr    = rE - rW

      ddr_r2_Tpr   = (rE**2*dataE(ir, iz, iTpr) - rW**2*dataW(ir, iz, iTpr)) / dr
      ddz_r2_Tpz   = (rN**2*dataN(ir, iz, iTpz) - rS**2*dataS(ir, iz, iTpz)) / dz

      ddz_upt = (dataN(ir,iz,iupt) - dataS(ir,iz,iupt)) / dz

      ddr_upt_r = (dataE(ir, iz, iupt)*rE - dataW(ir, iz, iupt)*rW) / dr

      ur_tilde = dataM(ir, iz, iurt)

      uK = sqrt(GM/rcen)

      ratio = 1.d0 / (rcen * dataM(ir, iz, irhob) * ddr_upt_r)
      upt_r2 = dataM(ir, iz, iupt) * rcen**2
      ur_tilde_theory  = ratio*(-ddr_r2_Tpr - ddz_r2_Tpz - rcen**2*dataM(ir,iz,irhob)*data(ir,iz,iuzt)*ddz_upt)
      ! ur_tilde_theory2 = ratio*(-ddr_r2_Tpr - ddz_r2_Tpz)

      ratio = 2.d0/(rcen * dataM(ir, iz, irhob) * uK)
      ur_tilde_theory3 = ratio*(-ddr_r2_Tpr - ddz_r2_Tpz)
      ! ur_tilde_theory4 = ratio*(-ddr_r2_Tpr)      
      
      mass_flux         = 2.d0 * pi* dataM(ir,iz,irhob) * ur_tilde         * rcen
      mass_flux_theory  = 2.d0 * pi* dataM(ir,iz,irhob) * ur_tilde_theory  * rcen
      ! mass_flux_theory2 = 2.d0 * pi* dataM(ir,iz,irhob) * ur_tilde_theory2 * rcen
      mass_flux_theory3 = 2.d0 * pi* dataM(ir,iz,irhob) * ur_tilde_theory3 * rcen
      ! mass_flux_theory4 = 2.d0 * pi* dataM(ir,iz,irhob) * ur_tilde_theory4 * rcen                       
      
      write(lun_tecplot_out, "(5(1x, e12.5))") rcen, zcen, mass_flux_theory, mass_flux_theory3, &
            mass_flux
   end do
end do
close(lun_tecplot_out)

! Determine iz for the midplane and calculate required quantities at the midplane:
!!$if (nz_is_odd) then
!!$   do ir = 1, nr
!!$      Tpr_mid(ir)      = T_pr    (ir, iz_mid)
!!$      ur_tilde_mid(ir) = ur_tilde(ir, iz_mid)
!!$      rho_bar_mid(ir)  = rho_bar (ir, iz_mid)
!!$      print *, ' ir = ', ir, ' Tpr_mid = ', Tpr_mid(ir), ' ur_tilde_mid = ', ur_tilde_mid(ir)
!!$   end do
!!$else
!!$   do ir = 1, nr
!!$      Tpr_mid(ir)      = 0.5d0*(T_pr    (ir, iz_mid1) + T_pr    (ir, iz_mid2))
!!$      ur_tilde_mid(ir) = 0.5d0*(ur_tilde(ir, iz_mid1) + ur_tilde(ir, iz_mid2))
!!$      rho_bar_mid(ir)  = 0.5d0*(rho_bar (ir, iz_mid1) + rho_bar (ir, iz_mid2))
!!$      print *, ' ir = ', ir, ' Tpr_mid = ', Tpr_mid(ir), ' ur_tilde_mid = ', ur_tilde_mid(ir)
!!$      print *, ' iz_mid1 = ', iz_mid1, ' iz_mid2 = ', iz_mid2
!!$      print *, ' ur_tilde(ir, iz_mid1) = ', ur_tilde(ir, iz_mid1)
!!$   end do
!!$end if

! Calculate ur_tilde at the mid-plane predicted from the angular momentum equation:
!!$open(unit = lunout, file = 'ur_tilde_theory.dat', form = 'formatted', status = 'unknown')
!!$do ir = 1, nr-1
!!$   ! Take finite difference:
!!$   right = rgrid(ir+1)**2.d0 * Tpr_mid(ir+1)
!!$   left  = rgrid(ir  )**2.d0 * Tpr_mid(ir)
!!$   dr    = rgrid(ir+1) - rgrid(ir)
!!$   deriv = (right - left) / dr
!!$
!!$   rave  = 0.5d0*(rgrid(ir+1) + rgrid(ir))
!!$   uK = sqrt(GM / rave)
!!$   ur_tilde_theory = - 2.d0*deriv / (rave*rho_bar_mid(ir)*uK)
!!$   write(lunout, "(3(1x, e12.5))") rave, ur_tilde_theory, ur_tilde_mid(ir)
!!$end do
!!$close(lunout)
!!$print *, ' Finished writing ur_tilde_theory.dat'

!*****************************************************************************80

! Post proceess meridional velocity perturbation file.

200 continue

nvar = 5

r_var = 1; z_var = 2; irho = 3; iur = 4; iuz = 5

write(6, "(' enter name of merid file--->', $)")
read(5, *) filename
open(unit = lun_tecplot_in, file = filename, form = 'formatted', status = 'old')

read (lun_tecplot_in, *) ! title
read (lun_tecplot_in, *) ! Variable list
! VARIABLES = "r", "z", "rho", "ur", "uz", "uphi", "uphi*r", "Mach_r", "Mach_z", "rho*ur", "rho*uz"
read (lun_tecplot_in, "(8x, i4, 4x, i4)") nr, nz

print *, ' nr ', nr, ' nz = ', nz
print *, ' enter anything to continue'

! Ready to allocate:
allocate(data (nr, nz, nvar))
allocate(dataW(nr, nz, nvar))
allocate(dataE(nr, nz, nvar))
allocate(dataN(nr, nz, nvar))
allocate(dataS(nr, nz, nvar))
allocate(dataM(nr, nz, nvar)) ! at midpoint of cell

do iz = 1, nz
   do ir = 1, nr
      read(lun_tecplot_in, "(5(1x, e12.5))") &
           data(ir, iz, r_var),  data(ir, iz, z_var), &
           data(ir, iz, irho), data(ir, iz, iur), data(ir, iz, iuz)
   end do
end do
close(lun_tecplot_in)

do iz = 1, nz - 1
   do ir = 1, nr - 1
      ! Obtain data at north, south, east, and west, and center of each cell:
      do ivar = 1, nvar
         dataW(ir, iz, ivar) = 0.5d0*(data(ir,   iz,   ivar) + data(ir,   iz+1, ivar))
         dataE(ir, iz, ivar) = 0.5d0*(data(ir+1, iz,   ivar) + data(ir+1, iz+1, ivar))
         dataN(ir, iz, ivar) = 0.5d0*(data(ir,   iz+1, ivar) + data(ir+1, iz+1, ivar))
         dataS(ir, iz, ivar) = 0.5d0*(data(ir,   iz,   ivar) + data(ir+1, iz,   ivar))

         dataM(ir, iz, ivar) = 0.25d0*(dataW(ir,iz,ivar)+dataE(ir,iz,ivar)+&
                                       dataN(ir,iz,ivar)+dataS(ir,iz,ivar))
      end do
   end do
end do

! Header of output tecplot file:
open (unit = lun_tecplot_out, file = 'vel_deriv.tec', form = 'formatted', status = 'unknown')
write (lun_tecplot_out, 5)
5  format ('TITLE = "du_dz"')
write (lun_tecplot_out, 6)
6  format('VARIABLES="r","z","duz_dz", "1/r*ddr(r ur)", duz_dr')
write (lun_tecplot_out, 7) nr-1, nz-1
7 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT')

do iz = 1, nz - 1
   do ir = 1, nr - 1
      rcen = dataM(ir, iz, r_var)
      zcen = dataM(ir, iz, z_var)

      dz     = dataN(ir,iz,z_var) - dataS(ir,iz,z_var)
      duz_dz = (dataN(ir,iz,iuz  ) - dataS(ir,iz,iuz  )) / dz

      dr    = dataE(ir,iz,r_var) - dataW(ir,iz,r_var)
      d_rur = dataE(ir,iz,r_var)*dataE(ir,iz,iur) - dataW(ir,iz,r_var)*dataW(ir,iz,iur)
      rinv_ddr_r_ur = 1.d0 / rcen * d_rur / dr

      duz_dr = (dataE(ir,iz,iuz  ) - dataW(ir,iz,iuz  )) / dz

      write(lun_tecplot_out, "(5(1x, e12.5))") rcen, zcen, duz_dz, rinv_ddr_r_ur, duz_dr
   end do
end do
close(lun_tecplot_out)      

stop

!*****************************************************************************80

! Output H(r):
! ~~~~~~~~~~~~

300 continue

write(6, "(' enter rmin, rmax--->', $)")
read(5, *) rmin, rmax

exponent = 0.5d0*(3.d0 + q)
nr = 500
dr = (rmax - rmin) / (nr - 1)

open (unit = lun_tecplot_out,  file = 'H_of_r.tec',       form = 'formatted', status = 'unknown')
open (unit = lun_tecplot_out2, file = 'H_minus_of_r.tec', form = 'formatted', status = 'unknown')
write (lun_tecplot_out, 8)
8  format ('TITLE = "H of r"')
write (lun_tecplot_out, 9)
9  format('VARIABLES="r","z"')
write (lun_tecplot_out, 10) nr
10 format('ZONE I=', i4, ' DATAPACKING=POINT')
write (lun_tecplot_out2, 11)
11  format ('TITLE = "H minus of r"')
write (lun_tecplot_out2, 12)
12  format('VARIABLES="r","z"')
write (lun_tecplot_out2, 13) nr
13 format('ZONE I=', i4, ' DATAPACKING=POINT')

do ir = 1, nr
   r = rmin + (ir - 1)*dr
   H = (r/r0)**exponent
   write(lun_tecplot_out,  "(2(1x, e12.5))") r,  H
   write(lun_tecplot_out2, "(2(1x, e12.5))") r, -H
end do
close(lun_tecplot_out)
close(lun_tecplot_out2)
stop

end
!*****************************************************************************80

