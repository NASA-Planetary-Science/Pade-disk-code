!***************************************************************************************90

! pdf of vorticity components from a meridional plane file of perturbation vorticity and
! dilatation.

! Perturbation vorticity = d_omega

! Order of components is: r, z, phi, magnitude

integer, parameter :: lun_tecplot = 1, lun_list = 2, n_components = 5
integer, dimension(n_components) :: lunpdf
character(100) :: tec_file, pdf_file, list_file, line

! Suffix for output pdf file to indicate component.  p = phi, m = magnitude 
character(7), dimension(n_components) :: suffix
character(1) :: dum_char
integer :: nr, nz, ir, iz, irec, nbins, ibin, icomp, ifile, nfiles, nsamples
real(8) :: t, check, pi, c, cauchy
real(8), allocatable, dimension(:) :: rgrid, zgrid
real(8), allocatable, dimension(:,:,:) :: d_omega
real(8), dimension(n_components) :: d_omega_min, d_omega_max, bin_width
real(8), allocatable, dimension(:,:) :: pdf
real(8) :: bin_cen

suffix(1) = 'r------'
suffix(2) = 'z------'
suffix(3) = 'p------'
suffix(4) = 'mag----'
suffix(5) = 'p_sgn_z'

do icomp = 1, n_components
   lunpdf(icomp) = 10 + icomp
end do

write(6, 10)
10 format(' enter name of file containing list of pert_vort_and_dil tec files--->', $)
read(5, *) list_file

write(6, 11)
11 format(' enter # of bins--->', $)
read(5, *) nbins

! Last component will be omega_phi*sgn(z)
allocate(pdf(n_components, nbins))

open(unit = lun_list, file = list_file, form = 'formatted', status = 'old')

nsamples = 0
do ifile = 1, 10000
   print *, ' ifile = ', ifile
   
   read(lun_list, "(a100)", end=999) tec_file

   print *, ' Attempting to open ifile = ', ifile, ' tec_file = ', tec_file
   open (unit = lun_tecplot, file = tec_file, form = 'formatted', status = 'old', &
            access = 'direct', recl = 13*7 + 1)

   read(lun_tecplot, 1, rec = 1) t
   1 format(13x, e12.5)
   print *, ' t = ', t
   read(lun_tecplot, 2, rec = 2) dum_char
   2 format(a1)
   print *, ' about to read nr and nz'
   read(lun_tecplot, 3, rec = 3) nr, nz
   3 format(7x, i4, 4x, i4)
   print *, ' nr = ', nr, ' nz = ', nz

   if (ifile .eq. 1) then
      allocate(rgrid(nr), zgrid(nz))
      allocate(d_omega(n_components,nr,nz))
   end if

   ! Read meridional plane of data:
   do iz = 1, nz
      do ir = 1, nr
         irec = (iz - 1)*nr + ir + 3
         read(lun_tecplot, "(7(1x, e12.5))", rec = irec) &
            ! These are all relative to the basic state:
             rgrid(ir), zgrid(iz), (d_omega(icomp,ir,iz), icomp = 1, n_components-1)

         ! omega_phi * sgn(z)
         icomp = n_components
         d_omega(icomp,ir,iz) = d_omega(3,ir,iz) * sign(1.d0, zgrid(iz))
      end do
   end do
   
   close(lun_tecplot)

   if (ifile .eq. 1) then
      ! Determine range of the data:
      do icomp = 1, n_components
         d_omega_min(icomp) =  1.d25
         d_omega_max(icomp) = -1.d25
         do iz = 1, nz
            do ir = 1, nr
               d_omega_min(icomp) = min(d_omega_min(icomp), d_omega(icomp,ir,iz))
               d_omega_max(icomp) = max(d_omega_max(icomp), d_omega(icomp,ir,iz))         
            end do
         end do
         print *, ' icomp = ', icomp, ' min = ', d_omega_min(icomp), ' max = ', d_omega_max(icomp)
      end do
   
      do icomp = 1, n_components
         bin_width(icomp) = (d_omega_max(icomp) - d_omega_min(icomp)) / nbins
      end do
   end if

   do iz = 1, nz
      do ir = 1, nr
         nsamples = nsamples + 1
         do icomp = 1, n_components
            ibin = (d_omega(icomp,ir,iz) - d_omega_min(icomp)) / bin_width(icomp) + 1
            if ((ibin .ge. 1) .and. (ibin .le. nbins)) then
               pdf(icomp, ibin) = pdf(icomp,ibin) + 1
            end if
         end do
      end do
   end do
end do

print *, ' Number of files exceeds maximum allowed (10000)'
stop

999 continue
close(lun_list)
nfiles = ifile - 1

do icomp = 1, n_components
   do ibin = 1, nbins
      pdf(icomp, ibin) = pdf(icomp, ibin) / nsamples / bin_width(icomp)
   end do
end do

! Check that each pdf integrates to unity
do icomp = 1, n_components
   check = 0.d0
   do ibin = 1, nbins
      check = check + pdf(icomp, ibin)*bin_width(icomp)
   end do
   print *, ' icomp = ', icomp, ' check = ', check
end do

do icomp = 1, n_components
   write(pdf_file, "('pdf_', a7, '.dat')") suffix(icomp)   
   open(unit = lunpdf(icomp), file = pdf_file, form = 'formatted', status = 'unknown')    

   do ibin = 1, nbins
      bin_cen = d_omega_min(icomp) + (ibin - 0.5d0) * bin_width(icomp)
      write(lunpdf(icomp), "(5(1x, e12.5))") bin_cen, pdf(icomp, ibin)
   end do
   
   close(lunpdf(icomp))
end do

! Cauchy distribution to compare with pdf of d_omega_phi:
open(unit = 1, file = 'cauchy.dat', form = 'formatted', status = 'unknown')
icomp = 3
pi    = 4.d0*atan(1.d0)
c     = 4.35d0
do ibin = 1, nbins
   bin_cen = d_omega_min(icomp) + (ibin - 0.5d0) * bin_width(icomp)
   cauchy  = c / ( pi*(c**2 + bin_cen**2) )
   write(1, "(2(1x, e12.5))") bin_cen, cauchy
end do

end

!***************************************************************************************90
