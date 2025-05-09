!**********************************************************************************85

! This program reads files of the type: phi_Reynolds_averages_000065.1450_0035000.dat
! (for instance) and computes phi-time averages.

! p.s.: We have to do this because time average of phi Favre averages is not the
! the phi-time average.

! The file contain the list of files must contain just the list of files.

! A file called rz_grid.dat must also be present.  nr and nz will be read from
! it.

implicit none
character(45) filename
integer :: nfiles, ifile, ir, iz, nr, nz, ivar, nr_in_grid, nz_in_grid
integer, parameter :: lun_input = 1, lun_Reynolds = 2, lun_tecplot = 3, lun_grid = 4, &
     lun_nu_t = 7, lun_alpha = 8, lun_Tpr_int = 9

! For Manger check:
integer, parameter :: lun_alpha1 = 10
integer, parameter :: lun_Sigma  = 11, lun_L = 12, lun_vort_z = 13

! Indices:
integer, parameter :: rho = 1, rho_uphi = 2, rho_ur = 3, rho_uz = 4, &
     rho_uphi_ur = 5, rho_uphi2 = 6, rho_ur2 = 7, &
     rho_uz2 = 8, rho_uphi_uz = 9, rho_ur_uz = 10, uphi_index = 11

! Quantities we will read (phi Reynolds averages at an instant):
real(8) :: bar_read(11)

! phi-time Reynolds averages:
real(8), allocatable, dimension(:,:,:) :: bar

real(8), allocatable, dimension(:) :: rgrid, zgrid

! Final phi-time averages we will output:
real(8), allocatable, dimension(:, :) :: rho_bar, ur_tilde, uz_tilde, uphi_tilde

! Final stresses of phi-time average of Favre fluctuations we will output
real(8), allocatable, dimension(:, :) :: T_zz, T_pp, T_rr, T_rz, T_pr, T_pz, T_trace, &
                                         T_pR1 ! For Manger spherical check

! For computing alpha:
real(8), allocatable, dimension(:) :: Sigma, uphi_ave, T_pr_int, T_pR1_int, Omega_ave, &
                                      alpha
real(8) :: Lz, d_Omega_dr, nu_t, nu_t1, dz, Rspherical, alpha1

real(8) :: uz_z, up_z, ur_z, uz_r, up_r, ur_r, ddr_r_ur
real(8) :: nu_t_zz, nu_t_pp, nu_t_rr, nu_t_rz, nu_t_pr, nu_t_pz
real(8) :: S_zz, S_pp, S_rr, S_rz, S_pr, S_pz
real(8) :: S_trace, denom, ddr_r_uphi, S_fac
real(8), allocatable, dimension(:) :: L_Rossby, vort_z

real(8) :: pi, T_orbital0, H0, rho0, H0_over_r0, r0, Omega0, GM, V_Kep0, c0, pexp, qexp, &
     c_exponent, H_exponent
real(8) :: cs, H, P
logical :: use_Manger_p
character(80) :: input_file
integer :: n_yes, nbar, isel, ir_mid
character(2)  :: yorn
real(8), parameter :: gamma = 1.4d0

pi         = 4.0*ATAN(1.0d0)
T_orbital0 = 1.0d0
H0         = 1.0d0
rho0       = 1.0d0

! This is a parameter:
H0_over_r0 = 0.10d0
r0 = H0 / H0_over_r0

! Consequence of above:
Omega0 = 2.d0 * pi / T_orbital0
GM     = Omega0**2 * r0**3
V_Kep0 = SQRT(GM/r0)
c0     = H0_over_r0 * V_Kep0

! Exponents:
use_Manger_p = .false.
if (use_Manger_p) then
   pexp = -2.d0/3.d0  ! for density.  Manger
else
   pexp = -3.d0/2.d0  ! for density Nelson
end if
qexp = -1.0d0  ! for temperature
c_exponent = 0.5d0*qexp
H_exponent = (qexp + 3.d0) / 2.d0

nbar = 10

! Read grid file:
print *, ' Attempting to open rz_grid.dat'
open(unit = lun_grid, file = 'rz_grid.dat', form = 'formatted', status = 'old')
read(lun_grid, "(i6, i6)") nr, nz
print *, ' nr = ', nr, ' nz = ', nz

allocate(L_Rossby(nr), vort_z(nr))

allocate(bar(nbar, nr, nz))
allocate(rgrid(nr), zgrid(nz))

! Final output:
allocate(rho_bar(nr,nz), ur_tilde(nr, nz), uz_tilde(nr, nz), uphi_tilde(nr, nz))
allocate(T_zz(nr,nz), T_pp(nr,nz), T_rr(nr,nz), T_rz(nr,nz), T_pr(nr,nz), T_pz(nr,nz), &
     T_trace(nr, nz))
! For computing alpha:
allocate(Sigma(nr), uphi_ave(nr), T_pr_int(nr), Omega_ave(nr), alpha(nr))

! For Manger spherical:
allocate(T_pR1(nr,nz))
allocate(T_pR1_int(nr))

write(6, "(' enter name of file containing list of files to process--->', $)")
read(5, *) input_file

print *, ' Attempting to open input_file_for_phi_time_Favre_averages'
open(unit = lun_input, file = input_file, &
     form = 'formatted', status = 'unknown')

do ir = 1, nr
   read(lun_grid, "(e16.9)") rgrid(ir)
end do
   
do iz = 1, nz
   read(lun_grid, "(e16.9)") zgrid(iz)      
end do
close(lun_grid)

write(6, 17)
17 format(' enter 1: For analysis, including alpha, using phi-time averages',/, &
          '       2: To approxomate alpha using time average of azimuthal means',/, &
          ' --->', $)
read(5, *) isel

go to (100, 200) isel

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
100 continue

! First add up the Reynolds averages at different times:
bar = 0.d0
n_yes = 0
do ifile = 1, 10000
   read(lun_input, "(a45)", end=999) filename
   yorn = 'y'
   if (yorn .eq. 'y ') then
      n_yes = n_yes + 1
      print *, ' Attempting to open file = ', filename, ' File # = ', ifile
      open(unit = lun_Reynolds, file = filename, form = 'formatted', status = 'old')
      print *, ' Successful'

      do iz = 1, nz
         do ir = 1, nr
            read(lun_Reynolds, "(11(1x, e16.9))") (bar_read(ivar), ivar = 1, nbar)

            ! Add into the sum over times:
            do ivar = 1, nbar
               bar(ivar, ir, iz) = bar(ivar, ir, iz) + bar_read(ivar)
            end do
         end do
      end do
      close(lun_Reynolds)
   end if
end do

999 continue
close(lun_input)

! Done with reading all the files and doing sums.

! Make the sums into averages:
bar = bar / n_yes

! Now compute Favre averages (see notes of 11/24/20):

! phi-time Reynolds averaged density:
rho_bar(:, :) = bar(rho,:,:)
! phi-time Favre mean velocities:
ur_tilde  (:,:) = bar(rho_ur,  :,:) / bar(rho,:,:)
uz_tilde  (:,:) = bar(rho_uz,  :,:) / bar(rho,:,:)
uphi_tilde(:,:) = bar(rho_uphi,:,:) / bar(rho,:,:) 

! Reynolds stresses:
T_pr(:,:) = bar(rho_uphi_ur, :,:) - bar(rho_ur,  :,:)*bar(rho_uphi,:,:)/bar(rho,:,:)
T_pp(:,:) = bar(rho_uphi2,   :,:) - bar(rho_uphi,:,:)*bar(rho_uphi,:,:)/bar(rho,:,:)
T_rr(:,:) = bar(rho_ur2,     :,:) - bar(rho_ur,  :,:)*bar(rho_ur,  :,:)/bar(rho,:,:)
T_zz(:,:) = bar(rho_uz2,     :,:) - bar(rho_uz,  :,:)*bar(rho_uz,  :,:)/bar(rho,:,:)
T_pz(:,:) = bar(rho_uphi_uz, :,:) - bar(rho_uphi,:,:)*bar(rho_uz,  :,:)/bar(rho,:,:)
T_rz(:,:) = bar(rho_ur_uz,   :,:) - bar(rho_ur,  :,:)*bar(rho_uz,  :,:)/bar(rho,:,:)
T_trace = T_rr + T_zz + T_pp

do iz = 1, nz
   do ir = 1, nr
      Rspherical = sqrt(rgrid(ir)**2 + zgrid(iz)**2)
      T_pR1(ir, iz) = (T_pr(ir, iz)*rgrid(ir) + T_pz(ir, iz)*zgrid(iz))/Rspherical
   end do
end do

! Output mean density and Favre mean velocities:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open (unit = lun_tecplot, file = 'phi_time_Favre_means.tec', form = 'formatted', status = 'unknown')
write (lun_tecplot, 1)
1 format(' TITLE = "phi-time Favre means"')
write (lun_tecplot, 2)
2 format(' VARIABLES="r","z","rho_bar","ur_tilde","uz_tilde","uphi_tilde", "rho_ur_bar", "rho_uz_bar"')
write (lun_tecplot, 3) nr, nz
3 format(' ZONE I=',i4,',',' J=',i4,' DATAPACKING=POINT')
do iz = 1, nz
   do ir = 1, nr
      write(lun_tecplot, "(8(1x, e12.5))") &
           rgrid(ir), zgrid(iz), &
           rho_bar (ir, iz), ur_tilde(ir, iz), uz_tilde(ir, iz), uphi_tilde(ir, iz), &
           bar(rho_ur,ir,iz), bar(rho_uz,ir,iz)
   end do
end do
close(lun_tecplot)

! Output stresses:
! ~~~~~~~~~~~~~~~~
open (unit = lun_tecplot, file = 'phi_time_Reynolds_stresses.tec', form = 'formatted', status = 'unknown')
write (lun_tecplot, 11)
11 format (' TITLE = "phi-time Reynolds streses"')
write (lun_tecplot, 12)
12 format(' VARIABLES="r","z","T_phiphi","T_phir","T_phiz","T_rr","T_rz","T_zz", "Trace"')
write (lun_tecplot, 13) nr, nz
13 format(' ZONE I=',i4,',',' J=',i4,' DATAPACKING=POINT')
do iz = 1, nz
   do ir = 1, nr
      write(lun_tecplot, "(9(1x, e12.5))") &
           rgrid(ir), zgrid(iz), &
           T_pp(ir, iz), T_pr(ir, iz), T_pz(ir, iz), T_rr(ir, iz), &
           T_rz(ir, iz), T_zz(ir, iz), T_trace(ir, iz)
   end do
end do
close(lun_tecplot)

! Compute alpha:
T_pr_int  = 0.d0
T_pR1_int = 0.d0  ! To check against Manger and Klahr.
Sigma     = 0.d0 ! Surface density
uphi_ave  = 0.d0 ! Vertical average of uphi_tilde
do iz = 1, nz - 1
   dz = zgrid(iz + 1) - zgrid(iz)
   do ir = 1, nr
      Sigma(ir)    = Sigma(ir)    + 0.5d0*(rho_bar   (ir,iz) + rho_bar   (ir,iz+1))*dz
      uphi_ave(ir) = uphi_ave(ir) + 0.5d0*(uphi_tilde(ir,iz) + uphi_tilde(ir,iz+1))*dz
      T_pr_int(ir)  = T_pr_int (ir) + 0.5d0*(T_pr (ir,iz) + T_pr (ir,iz+1))*dz
      T_pR1_int(ir) = T_pR1_int(ir) + 0.5d0*(T_pR1(ir,iz) + T_pR1(ir,iz+1))*dz
   end do
end do



Lz = zgrid(nz) - zgrid(1)
do ir = 1, nr
   uphi_ave(ir) = uphi_ave(ir) / Lz 
   Omega_ave(ir) = uphi_ave(ir) / rgrid(ir)
end do
ir_mid = nr/2

open(unit = lun_alpha,    file = 'alpha_vs_r.dat',  form = 'formatted', status = 'unknown')
open(unit = lun_Tpr_int,  file = 'Tpr_int.dat',     form = 'formatted', status = 'unknown')
open(unit = lun_alpha1,   file = 'alpha1_vs_r.dat', form = 'formatted', status = 'unknown')
open(unit = lun_Sigma,    file = 'Sigma_norm.dat',       form = 'formatted', status = 'unknown')
open(unit = lun_L,        file = 'L_norm.dat',           form = 'formatted', status = 'unknown')
open(unit = lun_vort_z,   file = 'vort_z_norm.dat',      form = 'formatted', status = 'unknown')
do ir = 2, nr - 1
   d_Omega_dr = (Omega_ave(ir+1) - Omega_ave(ir-1)) / (rgrid(ir+1) - rgrid(ir-1))
   ddr_r_uphi = (uphi_ave(ir+1)*rgrid(ir+1) - uphi_ave(ir-1)*rgrid(ir-1)) / (rgrid(ir+1) - rgrid(ir-1))
   vort_z(ir) = ddr_r_uphi / rgrid(ir)
   ! Check:
   ! d_Omega_dr = -1.5d0 * Omega0/r0 * (rgrid(ir)/r0)**(-2.5d0)
   nu_t       = - T_pr_int(ir) / Sigma(ir) / rgrid(ir) / d_Omega_dr

   cs = c0 * (rgrid(ir)/r0)**c_exponent
   H  = H0 * (rgrid(ir)/r0)**H_exponent   
   write(lun_alpha,   "(2(1x, e12.5))") rgrid(ir), nu_t  / (cs*H)
   write(lun_Tpr_int, "(2(1x, e12.5))") rgrid(ir), T_pr_int(ir)

   P = Sigma(ir) * cs**2
   S_fac = P / Sigma(ir)**gamma
   L_Rossby(ir) = Sigma(ir) / (2.d0 * vort_z(ir)) * S_fac**(2.d0/gamma)

   alpha1 = T_pR1_int(ir) / (Sigma(ir) * cs**2)
   write(lun_alpha1,"(2(1x, e12.5))") rgrid(ir), alpha1
end do

do ir = 1, nr
   write(lun_Sigma,  "(2(1x, e12.5))") rgrid(ir), Sigma(ir)/Sigma(ir_mid)
   write(lun_L,      "(2(1x, e12.5))") rgrid(ir), L_Rossby(ir) / L_Rossby(ir_mid)
   write(lun_vort_z, "(2(1x, e12.5))") rgrid(ir), vort_z(ir) / vort_z(ir_mid)
end do

close(lun_alpha)
close(lun_Tpr_int)
close(lun_alpha1)
close(lun_Sigma)
close(lun_L)
close(lun_vort_z)

! Turbulent viscosity (see notes of Dec. 30, 2020):
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open (unit = lun_tecplot, file = 'turbulent_viscosity.tec', form = 'formatted', status = 'unknown')
write (lun_tecplot, 14)
14 format (' TITLE = "Turbulent Viscosity"')
write (lun_tecplot, 15)
15 format(' VARIABLES="r","z","nuT_pp","nuT_pr","nuT_pz","nuT_rr","nuT_rz","nuT_zz"')
write (lun_tecplot, 16) nr-2, nz-2
16 format(' ZONE I=',i4,',',' J=',i4,' DATAPACKING=POINT')
! Strain-rate tensor:
do iz = 2, nz - 1
   do ir = 2, nr-1
      uz_z = (uz_tilde  (ir, iz+1) - uz_tilde  (ir, iz-1)) / (zgrid(iz+1) - zgrid(iz-1))
      up_z = (uphi_tilde(ir, iz+1) - uphi_tilde(ir, iz-1)) / (zgrid(iz+1) - zgrid(iz-1))
      ur_z = (ur_tilde  (ir, iz+1) - ur_tilde  (ir, iz-1)) / (zgrid(iz+1) - zgrid(iz-1))

      uz_r = (uz_tilde  (ir+1, iz) - uz_tilde  (ir-1, iz)) / (rgrid(ir+1) - rgrid(ir-1))
      up_r = (uphi_tilde(ir+1, iz) - uphi_tilde(ir-1, iz)) / (rgrid(ir+1) - rgrid(ir-1))
      ur_r = (ur_tilde  (ir+1, iz) - ur_tilde  (ir-1, iz)) / (rgrid(ir+1) - rgrid(ir-1))

      ddr_r_ur = (rgrid(ir+1)*ur_tilde(ir+1, iz) - rgrid(ir-1)*ur_tilde(ir-1, iz)) / (rgrid(ir+1) - rgrid(ir-1))     
      S_zz = uz_z
      S_rz = 0.5d0*(uz_r + ur_z)
      S_pz = 0.5d0*up_z
      S_pr = 0.5d0*(-uphi_tilde(ir,iz)/rgrid(ir) + up_r)
      S_rr = ur_r
      S_pp = ur_tilde(ir,iz) / rgrid(ir)
      S_trace = S_rr + S_zz + S_pp

      denom = S_zz - 1.d0/3.d0*S_trace
      nu_t_zz = -(T_zz(ir, iz) - 1.d0/3.d0*T_trace(ir,iz)) / 2.d0 / rho_bar(ir,iz) / denom

      denom = S_rr - 1.d0/3.d0*S_trace
      nu_t_rr = -(T_rr(ir, iz) - 1.d0/3.d0*T_trace(ir,iz)) / 2.d0 / rho_bar(ir,iz) / denom    

      denom = S_pp - 1.d0/3.d0*S_trace
      nu_t_pp = -(T_pp(ir, iz) - 1.d0/3.d0*T_trace(ir,iz)) / 2.d0 / rho_bar(ir,iz) / denom

      nu_t_pr = -T_pr(ir, iz) / 2.d0 / rho_bar(ir,iz) / S_pr
      nu_t_rz = -T_rz(ir, iz) / 2.d0 / rho_bar(ir,iz) / S_rz
      nu_t_pz = -T_pz(ir, iz) / 2.d0 / rho_bar(ir,iz) / S_pz

      write(lun_tecplot, "(8(1x, e12.5))") &
           rgrid(ir), zgrid(iz), &
           nu_t_pp, nu_t_pr, nu_t_pz, nu_t_rr, nu_t_rz, nu_t_zz      
   end do
end do

stop

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
200 continue

print *, ' After statement 200.  nbar = ', nbar

n_yes = 0
alpha = 0.d0
do ifile = 1, 10000
   read(lun_input, "(a45)", end=998) filename
   yorn = 'y'
   if (yorn .eq. 'y ') then
      n_yes = n_yes + 1
      print *, ' Attempting to open file = ', filename, ' File # = ', ifile
      open(unit = lun_Reynolds, file = filename, form = 'formatted', status = 'old')
      print *, ' Successful'

      do iz = 1, nz
         do ir = 1, nr
            read(lun_Reynolds, "(10(1x, e16.9))") (bar(ivar,ir,iz), ivar = 1, nbar)
         end do
      end do
      close(lun_Reynolds)
   end if

   T_pr(:,:) = bar(rho_uphi_ur, :,:) - bar(rho_ur,:,:)*bar(rho_uphi,:,:)/bar(rho,:,:)
   rho_bar(:, :) = bar(rho,:,:)
   uphi_tilde(:,:) = bar(rho_uphi,:,:) / bar(rho,:,:) 

   T_pr_int  = 0.d0
   Sigma     = 0.d0 ! Surface density
   uphi_ave  = 0.d0 ! Vertical average of uphi_tilde
   do iz = 1, nz - 1
      dz = zgrid(iz + 1) - zgrid(iz)
      do ir = 1, nr
         Sigma(ir)    = Sigma(ir)    + 0.5d0*(rho_bar   (ir,iz) + rho_bar   (ir,iz+1))*dz
         uphi_ave(ir) = uphi_ave(ir) + 0.5d0*(uphi_tilde(ir,iz) + uphi_tilde(ir,iz+1))*dz
         T_pr_int(ir) = T_pr_int(ir) + 0.5d0*(T_pr   (ir,iz) + T_pr   (ir,iz+1))*dz
      end do
   end do

   Lz = zgrid(nz) - zgrid(1)
   do ir = 1, nr
      uphi_ave(ir) = uphi_ave(ir) / Lz 
      Omega_ave(ir) = uphi_ave(ir) / rgrid(ir)
   end do
   
   do ir = 2, nr - 1
      d_Omega_dr = (Omega_ave(ir+1) - Omega_ave(ir-1)) / (rgrid(ir+1) - rgrid(ir-1))
      nu_t       = - T_pr_int(ir) / Sigma(ir) / rgrid(ir) / d_Omega_dr

      cs = c0 * (rgrid(ir)/r0)**c_exponent
      H  = H0 * (rgrid(ir)/r0)**H_exponent
      alpha(ir) = alpha(ir) + nu_t / (cs*H)
   end do
end do

998 continue

print *, ' Taking the average of all the alpha s'
alpha = alpha / n_yes

open(unit = lun_alpha,  file = 'alpha1_vs_r.dat',  form = 'formatted', status = 'unknown')
do ir = 2, nr - 1
   write(lun_alpha, "(2(1x, e12.5))") rgrid(ir), alpha(ir)
end do
close(lun_alpha)

close(lun_input)
stop

end program

!**********************************************************************************85

