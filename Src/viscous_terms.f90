!----------------------------------------------------------------------------------85

! This file contains routines for the treatment of viscosity, thermal conductivity, and
! viscous heating.

! Currently you can choose one of:
! (1) Laminar

! The sub-grid models due to:
! (2) Moin et. al., Phys. Fluids A, 3(11), p. 2746, (1991) ---  note we use the non-dynamic version.
! (3) Vreman,       Phys. Fluids   16(10), p. 3670, (2004).

! (4) Dilatation dependent shear viscosity to capture shear layers caused by shocks.

! For the internal energy eq. see Batchelor, p. 153.
! Our internal energy is rho Cv T.
! Heat conduction term on rhs of internal energy eq. = d/dx_i (k dT/dx_i).
! Viscous heating term is 2 mu S_ij^2 + (mu_b - 2/3 mu) Dilatation^2

!----------------------------------------------------------------------------------85

module viscosity_types
   integer, parameter :: molecular = 1, Moin_etal = 2, Vreman = 3, ddsv = 4, Moin_ddsv = 5, &
                         Vreman_ddsv = 6
end module viscosity_types

! Viscosity and conductivity arising from molecular terms or subgrid scale model or dilatation-dependent shear-viscosity.
module viscous
   logical :: apply_viscosity
   integer :: viscosity_type

   ! These are used so we don't separately calculate the dilatation if we can
   ! get it from the velocity gradient tensor or the strain tensor.
   logical :: have_velocity_gradient_tensor, have_strain_tensor

   real(8) :: C_DDSV ! coefficient for dilatation dependent shear viscosity.

   ! Constants of the LES models (Moin et al and Vreman):
   ! real(8), parameter :: C_Smag = 0.02d0, C_Yoshizawa = 0.05d0,  Pr_t = 1.0d0, C_Vreman = 2.5d0*C_Smag**2

   ! Standard value is 0.0484 = 0.22^2
   real(8) :: C_Smag

   real(8), parameter :: C_Yoshizawa = 0.05d0,  Pr_t = 1.0d0

   ! This is needed to go from viscosity to heat conductivity:
   real(8) :: gamma_over_Pr_molecular, gamma_over_Pr_t

   ! Molecular kinematic shear and bulk viscosities and Prandtl number.
   real(8) :: nu_molecular, nu_b_molecular, Pr_molecular
   
   ! Will also be used for the stress tensor to save space:
   real(8), allocatable, dimension(:, :, :) :: T_zz, T_zr, T_pz, T_rr, T_pr, T_pp

   ! Temperature gradient:
   real(8), allocatable, dimension(:, :, :) :: grad_cvT_z, grad_cvT_r, grad_cvT_p

   ! Viscosity for plotting purposes:
   real(8), allocatable, dimension(:, :, :) :: nu_plot   
   
   ! Viscous force = divergence of the stress
   real(8), allocatable, dimension(:, :, :) :: Fz, Fr, Fp
   
   real(8), allocatable, dimension(:, :, :, :) :: heat_flux
   ! This includes viscous heating + heat conduction.
   real(8), allocatable, dimension(:, :, :) :: heat_term

   ! Velocity gradient tensor for Vreman's subgrid model:
   real(8), allocatable, dimension(:, :, :, :, :) :: G
end module viscous

!----------------------------------------------------------------------------------85

subroutine activate_viscosity(apply_viscosity_arg, viscosity_type_arg, isothermal_arg, &
     nu_molecular_arg, nu_b_molecular_arg, Pr_molecular_arg, gamma_arg, C_DDSV_arg, C_Smag_arg)

! This routine activates one of four viscosity (and conductivity) treatments and should
! be called by the user's application subroutine.

! Set apply_viscosity = .true. if you want one of the four viscosities and conductivities
! or one of the two mass diffusivites.

! viscosity_type_arg : Choose one of the following integers (defined in module viscosity_types):
!                      Molecular = 1  : Molecular
!                      Moin_etal = 2  : Non-dynamic version of Moin etal sgs model
!                      Vreman    = 3  : Vreman sgs model + same compressibility additions as Moin et al
!                      ddsv      = 4  : Dilatation-dependent shear viscosity
!                      Moin_ddsv = 5  : 2 and 4.

! isothermal_arg : .true. or .false.  Tells us if we should allocate storage for thermal stuff.

! nu_molecular_arg   : Molecular kinematic shear viscosity
! nu_b_molecular_arg : Molecular kinematic bulk viscosity
! Pr_molecular_arg : Molecular Prandtl number
! gamma_arg      : cp/cv, ratio of specific heats.  This is needed to go from viscosity to
!                  heat conductivity for the molecular case.
! C_DDSV_arg : coefficient for the dilatation dependent shear viscosity used to capture shear
! layers caused by shocks in more than 1 dimension.  Currently I don't think this is a useful
! option.  Nevertheless I have left it in and may remove it in the future.

use grid
use viscous
use viscosity_types
use total_allocated_words, only: n_words_allocated
use partition_data
use dof_indices
use logical_units
implicit none
logical      :: apply_viscosity_arg
integer      :: viscosity_type_arg
logical      :: isothermal_arg
real(8)      :: nu_molecular_arg ! Molecular kinematic shear viscosity (constant).
real(8)      :: nu_b_molecular_arg ! Molecular kinematic bulk viscosity (constant).
real(8)      :: Pr_molecular_arg ! Molecular Prandtl number (constant).
real(8)      :: gamma_arg
real(8)      :: C_DDSV_arg
real(8)      :: C_Smag_arg

! Local:
integer :: ir, iz

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: activate_viscosity: first executable'
#endif

apply_viscosity = apply_viscosity_arg
viscosity_type  = viscosity_type_arg
nu_molecular    = nu_molecular_arg
nu_b_molecular  = nu_b_molecular_arg
Pr_molecular    = Pr_molecular_arg
C_DDSV          = C_DDSV_arg
C_Smag          = C_Smag_arg

if (my_node .eq. 0) then
   print *, ' node 0: subroutine activate_viscosity'
   print *, '    apply_viscosity = ', apply_viscosity
   print *, '    viscosity_type  = ', viscosity_type
   print *, '    nu_molecular    = ', nu_molecular
   print *, '    nu_b_molecular  = ', nu_b_molecular   
   print *, '    Pr_molecular    = ', Pr_molecular
   print *, '    C_DDSV          = ', C_DDSV
   print *, '    C_Smag          = ', C_Smag
end if

! See notes of Nov. 30, 2018.  We need these to calculate k/cv.
gamma_over_Pr_molecular = gamma_arg / Pr_molecular
gamma_over_Pr_t         = gamma_arg / Pr_molecular

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' node 0:'
   print *, ' viscosity_type = ', viscosity_type
   print *, ' C_DDSV         = ', C_DDSV
end if
#endif

! Check that viscosity_type is one of the allowable types:
if (my_node .eq. 0) then
   if (viscosity_type .eq. molecular) then
      print *, ' viscosity type = molecular'
   else if (viscosity_type .eq. Moin_etal) then
      print *, ' viscosity type = Moin et al LES model (non-dynamic)'
   else if (viscosity_type .eq. Vreman) then
      print *, ' viscosity type = Vreman LES model + compressibility additions'
   else if (viscosity_type .eq. ddsv) then
      print *, ' viscosity type = Dilatation-dependent shear viscosity'
   else if (viscosity_type .eq. Moin_ddsv) then
      print *, ' viscosity type = Moin et al + Dilatation-dependent shear viscosity'
   else if (viscosity_type .eq. Vreman_ddsv) then
      print *, ' viscosity type = Vreman + Dilatation-dependent shear viscosity'
   else
      print *, ' Unrecognized viscosity_type = ', viscosity_type
      call terminate_with_no_save(1)
   end if
end if

call allocate_viscous_quantities(isothermal_arg)

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: After activate_viscosity, Mb allocated = ', float(8*n_words_allocated)/1.d6
#endif

! Delta vector for Vreman model:
if ((viscosity_type .eq. Vreman) .or. (Viscosity_type .eq. Vreman_ddsv)) then
   allocate(Delta_vector(3, sr:er, nz)) ! This is in module grid
   if ((nz_actual .ne. 1) .and. (nr .ne. 1) .and. (nphi .ne. 1)) then
   ! 3D case:
      do iz = 1, nz
         do ir = sr, er
            Delta_vector(r_comp, ir, iz) = dr(ir)
            Delta_vector(z_comp, ir, iz) = dz(iz)
            Delta_vector(p_comp, ir, iz) = r_dphi(ir)
         end do
      end do
   else if ((nz_actual .eq. 1) .and. (nr .ne. 1) .and. (nphi .ne. 1)) then
      ! Planar case (vertically integrated, r, phi):   This needs to be redone.
      do iz = 1, nz
         do ir = sr, er
            Delta_vector(r_comp, ir, iz) = dr(ir)
            Delta_vector(z_comp, ir, iz) = 0.0d0
            Delta_vector(p_comp, ir, iz) = r_dphi(ir)         
         end do
      end do
   else if ((nz_actual .ne. 1) .and. (nr .ne. 1) .and. (nphi .eq. 1)) then
      ! Axisymmetric case (z, r):
      do iz = 1, nz
         do ir = sr, er
            Delta_vector(r_comp, ir, iz) = dr(ir)
            Delta_vector(z_comp, ir, iz) = dz(iz)
            Delta_vector(p_comp, ir, iz) = 0.0d0         
         end do
      end do
   else if ((nz_actual .ne. 1) .and. (nr .eq. 1) .and. (nphi .eq. 1)) then
      ! Vertical only case:
      do iz = 1, nz
         do ir = sr, er
            Delta_vector(r_comp, ir, iz) = 0.0d0
            Delta_vector(z_comp, ir, iz) = dz(iz)
            Delta_vector(p_comp, ir, iz) = 0.0d0         
         end do
      end do
   else
      print *, ' You are running a new case type for which you need to add coding for'
      print *, ' the minimum grid size in subroutine make_grid'
      call terminate_with_no_save(1)
   end if
end if

end subroutine activate_viscosity

!----------------------------------------------------------------------------------85

subroutine allocate_viscous_quantities(isothermal_arg)

use grid
use partition_data
use viscous
use viscosity_types
use total_allocated_words, only: n_words_allocated
implicit none
logical :: isothermal_arg

! Strain-rate then over-written by the shear-stress:
allocate(T_zz(sr:er, sphi:ephi, nz))
allocate(T_zr(sr:er, sphi:ephi, nz))
allocate(T_pz(sr:er, sphi:ephi, nz))
allocate(T_rr(sr:er, sphi:ephi, nz))
allocate(T_pr(sr:er, sphi:ephi, nz))
allocate(T_pp(sr:er, sphi:ephi, nz))
n_words_allocated = n_words_allocated + mr*mphi*nz*6

if (.not. isothermal_arg) then
   allocate(heat_flux(sr:er, sphi:ephi, nz, 3))
   allocate(heat_term(sr:er, sphi:ephi, nz))  ! viscous heating + divergence of heat flux   
   n_words_allocated = n_words_allocated + mr*mphi*nz*7
end if

! Viscous forces:
allocate(Fr(sr:er, sphi:ephi, nz))
allocate(Fz(sr:er, sphi:ephi, nz))
allocate(Fp(sr:er, sphi:ephi, nz))
n_words_allocated = n_words_allocated + mr*mphi*nz*3

! For plotting:
allocate(nu_plot(sr:er, sphi:ephi, nz))
n_words_allocated = n_words_allocated + mr*mphi*nz

! Need velocity gradient tensor for Vreman model:
if ((viscosity_type .eq. Vreman) .or. (viscosity_type .eq. Vreman_ddsv)) then
   allocate(G(sr:er, sphi:ephi, nz, 3, 3))
   n_words_allocated = n_words_allocated + 9*mr*mphi*nz      
end if

end subroutine allocate_viscous_quantities

!----------------------------------------------------------------------------------85

subroutine add_viscous_terms_to_qdot(t, q, qdot, lambda_max)

! This routine is called by rhs and adds the following contributions to qdot:
!    (1) Viscous force in momentum equation
!    (2) Heat conduction in energy equation (for the non-isothermal case)
!    (3) Viscous heating in energy equation (for the non-isothermal case)

! lambda_max is the eigenvalue for these terms (used for time-step selection).

use partition_data
use grid
use viscous
use viscosity_types
use dof_indices
use thermal_parameters
implicit none
real(8) :: t
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q, qdot
real(8) :: lambda_max ! eigenvalue for time step determination.

! Local:
integer :: iz, iphi, ir

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: add_viscous_terms_to_qdot has been called'
#endif

if (viscosity_type .eq. molecular) then
   call molecular_viscous_stress_and_heat_term(q, lambda_max)
else if (viscosity_type .eq. Moin_etal) then
   call Moin_etal_stress_and_heat_flux(q, lambda_max)
else if (viscosity_type .eq. Vreman) then            
   call Vreman_stress_and_heat_flux(q, lambda_max)      
else if (viscosity_type .eq. ddsv) then            
   call ddsv_stress_and_heat_flux(q, lambda_max)
else if (viscosity_type .eq. Moin_ddsv) then            
   call Moin_ddsv_stress_and_heat_flux(q, lambda_max)
else if (viscosity_type .eq. Vreman_ddsv) then            
   call Vreman_ddsv_stress_and_heat_flux(q, lambda_max)           
end if

! Take divergence of the stresses and heat flux:
call viscous_force_and_heat_term

! The mass equation does not receive a sub-grid contribution due to the use of Favre averages
! in Moin et al.

! Add terms to qdot:
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         qdot(ir,iphi,iz,zmom) = qdot(ir,iphi,iz,zmom) + Fz(ir,iphi,iz)
         qdot(ir,iphi,iz,rmom) = qdot(ir,iphi,iz,rmom) + Fr(ir,iphi,iz)
                                                         ! Viscous torque
         qdot(ir,iphi,iz,amom) = qdot(ir,iphi,iz,amom) + Fp(ir,iphi,iz)*rgrid(ir)                  
      end do
   end do
end do

! For the internal energy equation we use, see pg. 153 in Batchelor. eq. 3.4.4.
if (.not. isothermal) then
   ! heat_term = divergence of heat flux + viscous heating
   qdot(:,:,:,ener) = qdot(:,:,:,ener) + heat_term(:,:,:)                  
end if

end subroutine add_viscous_terms_to_qdot

!----------------------------------------------------------------------------------85

subroutine viscous_force_and_heat_term

! Computes the 3 viscous/subgrid force components in the momentum equation given the 6
! viscous/subgrid stresses.  Computes the viscous/subgrid heat conduction term.

! See notes of Sept. 6, 2018 for taking the divergence of the stress tensor.

use viscous
use partition_data
use grid
use fargo_or_plotting_shift
use dof_indices
use thermal_parameters
implicit none

! Local:
real(8), dimension(sr:er, sphi:ephi, nz) :: deriv_in_z_space

! These are dimensioned in z-space but will be used in all spaces (since pencils in
! all directions have the same size).
real(8), dimension(sr:er, sphi:ephi, nz) :: transposed, deriv

integer :: ir, iz, iphi

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: viscous_force_and_heat_term, first executable'
#endif

Fz = 0.0d0
Fr = 0.0d0
Fp = 0.0d0

! Certain stresses are multiplied by r:
do ir = sr, er
   T_zz(ir,:,:) = T_zz(ir,:,:) * rgrid(ir)
   T_zr(ir,:,:) = T_zr(ir,:,:) * rgrid(ir) 
   T_rr(ir,:,:) = T_rr(ir,:,:) * rgrid(ir)
end do

if (.not. isothermal) then
   do ir = sr, er
      heat_flux(ir,:,:,r_comp) = heat_flux(ir,:,:,r_comp) * rgrid(ir)
   end do
end if

! --------------
! z derivatives:
! --------------
if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
   ! Item 1: z-derivative of T_zz*r.  Paper item 1:
   call pade_diff_z(nbundle_z, T_zz, deriv_in_z_space)
   Fz = Fz + deriv_in_z_space

   ! Item 4: ddz of T_zr*r:
   call pade_diff_z(nbundle_z, T_zr, deriv_in_z_space)
   Fr = Fr + deriv_in_z_space
   
   ! Item 10: z-derivative of T_pz:
   call pade_diff_z(nbundle_z, T_pz, deriv_in_z_space)
   Fp = Fp + deriv_in_z_space

   ! Add z heat conduction term to the viscous heating term:
   if (.not. isothermal) then
      call pade_diff_z(nbundle_z, heat_flux(sr,sphi,1,z_comp), deriv_in_z_space)
      ! Recall that heat_term starts off as the viscous heating.
      heat_term = heat_term + deriv_in_z_space
   end if
end if

! ------------------------------
! r derivatives and Fargo terms:
! ------------------------------
if (nr .ne. 1) then
   ! Item 2: r-derivative of T_zr*r:   
   call transpose_z_to_r (1, T_zr, transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, deriv_in_z_space)
   Fz = Fz + deriv_in_z_space

   ! Item 5: r*T_rr:
   call transpose_z_to_r (1, T_rr, transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, deriv_in_z_space)
   Fr = Fr + deriv_in_z_space

   ! Item 11: T_pr
   call transpose_z_to_r (1, T_pr, transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, deriv_in_z_space)
   Fp = Fp + deriv_in_z_space

   if (.not. isothermal) then
      ! Heat term.  Item 14.
      call transpose_z_to_r (1, heat_flux(sr,sphi,1,r_comp), transposed)
      call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
      call transpose_r_to_z (1, deriv, deriv_in_z_space)
      do ir = sr, er
         heat_term(ir,:,:) = heat_term(ir,:,:) + deriv_in_z_space(ir,:,:)/rgrid(ir)
      end do
   end if

   ! Fargo terms:
   if ((add_fargo_extra_operator_now) .and. (nphi .ne. 1)) then
      ! fargo applied to item 2:
      call transpose_z_to_phi (1, T_zr, transposed)
      call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
      call transpose_phi_to_z (1, deriv, deriv_in_z_space)
      do ir = sr, er
         Fz(ir,:,:) = Fz(ir,:,:) + fargo_factor(ir)*deriv_in_z_space(ir,:,:)
      end do

      ! fargo applied to item 5:
      call transpose_z_to_phi (1, T_rr, transposed)
      call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
      call transpose_phi_to_z (1, deriv, deriv_in_z_space)
      do ir = sr, er
         Fr(ir,:,:) = Fr(ir,:,:) + fargo_factor(ir)*deriv_in_z_space(ir,:,:)
      end do

      ! fargo applied to item 11:
      call transpose_z_to_phi (1, T_pr, transposed)
      call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
      call transpose_phi_to_z (1, deriv, deriv_in_z_space)
      do ir = sr, er
         Fp(ir,:,:) = Fp(ir,:,:) + fargo_factor(ir)*deriv_in_z_space(ir,:,:)
      end do      

      if (.not. isothermal) then
         ! Item 15:
         call transpose_z_to_phi(1, heat_flux(sr,sphi,1,r_comp), transposed)
         call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
         call transpose_phi_to_z(1, deriv, deriv_in_z_space)   
         do ir = sr, er
            heat_term(ir,:,:) = heat_term(ir,:,:)+fargo_factor(ir)*deriv_in_z_space(ir,:,:)/rgrid(ir)
         end do
      end if
   end if ! Fargo
end if ! nr .ne. 1

! ----------------
! phi derivatives:
! ----------------
if (nphi .ne. 1) then
   ! Item 3: phi-derivative of T_zp:
   call transpose_z_to_phi(1, T_pz, transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z (1, deriv, deriv_in_z_space)
   Fz = Fz + deriv_in_z_space

   ! Item 6: phi-derivative of T_rp:
   call transpose_z_to_phi(1, T_pr, transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z (1, deriv, deriv_in_z_space)
   Fr = Fr + deriv_in_z_space

   ! Item 12: phi-derivative of T_pp and undifferentiated term:
   call transpose_z_to_phi(1, T_pp, transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z(1, deriv, deriv_in_z_space)
   do ir = sr, er
      Fp(ir,:,:) = Fp(ir,:,:) + deriv_in_z_space(ir,:,:)/rgrid(ir)
   end do

   if (.not. isothermal) then
      ! Heat term.  Item 16:
      call transpose_z_to_phi(1, heat_flux(sr,sphi,1,phi_comp), transposed)
      call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
      call transpose_phi_to_z(1, deriv, deriv_in_z_space)      
      do ir = sr, er
         heat_term(ir,:,:) = heat_term(ir,:,:) + deriv_in_z_space(ir,:,:)/rgrid(ir)
      end do
   end if
end if

! Item 3a: Divide by r and add undifferentiated terms:
do ir = sr, er
   Fz(ir,:,:) = Fz(ir,:,:)/rgrid(ir)
   Fr(ir,:,:) = (Fr(ir,:,:) - T_pp(ir,:,:)) / rgrid(ir)
   Fp(ir,:,:) = Fp(ir,:,:) + 2.d0/rgrid(ir)*T_pr(ir,:,:)   
end do

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: viscous_force_and_heat_term, returning'
#endif

end subroutine viscous_force_and_heat_term

!----------------------------------------------------------------------------------85

subroutine molecular_viscous_stress_and_heat_term(q, lambda_max)

use viscous
use partition_data
use grid
use dof_indices
use thermal_parameters
use math_constants
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: lambda_max

! Locals:
integer :: iz, iphi, ir, icomp
real(8) :: mu, mu_b, lambda_visc, S_mm, Sij_squared, k_over_cv

if (.not. have_strain_tensor) call strain_tensor(q)
if (.not. isothermal) call grad_cvT(q)
lambda_max = 0.0d0

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         mu   = q(ir,iphi,iz,irho) * nu_molecular   ! dynamic shear viscosity
         mu_b = q(ir,iphi,iz,irho) * nu_b_molecular ! dynamic bulk viscosity
         ! Another factor of pi comes later:
         lambda_max = MAX(lambda_max, pi * nu_molecular/min_grid_size(ir,iz)**2)
         S_mm = T_zz(ir,iphi,iz) + T_rr(ir,iphi,iz) + T_pp(ir,iphi,iz)
         ! This is a Lame parameter:
         lambda_visc = mu_b - 2.d0/3.d0*mu         
         
         ! Heat flux vector.  See notes of Nov. 9, 2018 in the "Governing Equations" tab.
         ! Note that we are overwriting the gradient of cvT which was sitting in heat_flux.
         if (.not. isothermal) then
            ! Conductivity divided by cv:
            k_over_cv = mu * gamma_over_Pr_molecular
            do icomp = 1, 3
               heat_flux(ir,iphi,iz,icomp) = k_over_cv * heat_flux(ir,iphi,iz,icomp)
            end do
            lambda_max = MAX(lambda_max, pi*nu_molecular*gamma_over_Pr_molecular/min_grid_size(ir,iz)**2)

            ! Viscous heating term:
            Sij_squared = 2.d0*T_zr(ir,iphi,iz)**2 + 2.d0*T_pz(ir,iphi,iz)**2 + 2.d0*T_pr(ir,iphi,iz)**2 + &
                               T_zz(ir,iphi,iz)**2 +      T_rr(ir,iphi,iz)**2 +      T_pp(ir,iphi,iz)**2
            
            heat_term(ir,iphi,iz) = 2.d0*mu*Sij_squared + lambda_visc*S_mm**2
         end if         
         
         ! This over-writes the strain tensor:
         T_zz(ir,iphi,iz) = 2.d0*mu *(T_zz(ir,iphi,iz) + lambda_visc*S_mm)
         T_rr(ir,iphi,iz) = 2.d0*mu *(T_rr(ir,iphi,iz) + lambda_visc*S_mm)
         T_pp(ir,iphi,iz) = 2.d0*mu *(T_pp(ir,iphi,iz) + lambda_visc*S_mm)
         T_zr(ir,iphi,iz) = 2.d0*mu * T_zr(ir,iphi,iz)
         T_pz(ir,iphi,iz) = 2.d0*mu * T_pz(ir,iphi,iz)
         T_pr(ir,iphi,iz) = 2.d0*mu * T_pr(ir,iphi,iz)
      end do
   end do
end do

! Since we have over-written it.
have_strain_tensor = .false.

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: molecular_stress_and_heat_flux: Returning'
#endif

end subroutine molecular_viscous_stress_and_heat_term

!----------------------------------------------------------------------------------85

subroutine ddsv_stress_and_heat_flux(q, lambda_max)

use viscous
use partition_data
use grid
use dof_indices
use thermal_parameters
use math_constants
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: lambda_max ! eigenvalue for time step determination.

! Locals:
integer :: iz, iphi, ir, icomp
real(8) :: nu, mu, S_mm, Sij_squared

if (.not. have_strain_tensor) call strain_tensor(q)
if (.not. isothermal        ) call grad_cvT(q)

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! See notes of Nov. 9, 2018:
         if (.not. isothermal) then
            Sij_squared = 2.d0*T_zr(ir,iphi,iz)**2 + 2.d0*T_pz(ir,iphi,iz)**2 + 2.d0*T_pr(ir,iphi,iz)**2 + &
                               T_zz(ir,iphi,iz)**2 +      T_rr(ir,iphi,iz)**2 +      T_pp(ir,iphi,iz)**2            
            heat_term(ir,iphi,iz) = 2.d0*mu*(Sij_squared - 1.d0/3.d0*S_mm**2)
            do icomp = 1, 3
               heat_flux(ir,iphi,iz,icomp) = 0.0d0
            end do
         end if
         
         S_mm = T_zz(ir,iphi,iz) + T_rr(ir,iphi,iz) + T_pp(ir,iphi,iz) ! = dilatation
         nu = C_DDSV*l_grid_squared(ir,iz) * abs(min(S_mm, 0.0d0))
         nu_plot(ir,iphi,iz) = nu
         mu = q(ir,iphi,iz,irho) * nu
         lambda_max = MAX(lambda_max, pi*nu/min_grid_size(ir,iz)**2)
         
         ! This over-writes the strain tensor:
         T_zz(ir,iphi,iz) = 2.d0*mu *(T_zz(ir,iphi,iz) - 1.d0/3.d0*S_mm)
         T_rr(ir,iphi,iz) = 2.d0*mu *(T_rr(ir,iphi,iz) - 1.d0/3.d0*S_mm)
         T_pp(ir,iphi,iz) = 2.d0*mu *(T_pp(ir,iphi,iz) - 1.d0/3.d0*S_mm)
         T_zr(ir,iphi,iz) = 2.d0*mu * T_zr(ir,iphi,iz)
         T_pz(ir,iphi,iz) = 2.d0*mu * T_pz(ir,iphi,iz)
         T_pr(ir,iphi,iz) = 2.d0*mu * T_pr(ir,iphi,iz)
      end do
   end do
end do

! Strain tensor has been over-written:
have_strain_tensor = .false.

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: add_dilatation_dependent_shear_viscosity: Returning'
#endif

end subroutine ddsv_stress_and_heat_flux

!----------------------------------------------------------------------------------85

subroutine Moin_etal_stress_and_heat_flux(q, lambda_max)

use viscous
use partition_data
use grid
use dof_indices
use gravity
use thermal_parameters
use math_constants
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: lambda_max ! eigenvalue for time step determination.

! Locals:
real(8) :: S_ij_squared, S_abs, S_mm, nu, mu, q2, k_over_cv
integer :: icomp, iz, iphi, ir

if (.not. have_strain_tensor) call strain_tensor(q)
if (.not. isothermal        ) call grad_cvT     (q)
lambda_max = 0.0d0

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         S_ij_squared = 2.d0*T_zr(ir,iphi,iz)**2 + 2.d0*T_pz(ir,iphi,iz)**2 + 2.d0*T_pr(ir,iphi,iz)**2 + &
                             T_zz(ir,iphi,iz)**2 +      T_rr(ir,iphi,iz)**2 +      T_pp(ir,iphi,iz)**2
         S_abs = SQRT(2.0d0 * S_ij_squared)
         S_mm = T_zz(ir,iphi,iz) + T_rr(ir,iphi,iz) + T_pp(ir,iphi,iz)

         ! Smagorinsky model:
         nu = C_Smag * l_grid_squared(ir,iz) * S_abs ! Smagorinsky model

         lambda_max = MAX(lambda_max, pi*nu/min_grid_size(ir,iz)**2)                  
         nu_plot(ir,iphi,iz) = nu
         mu = q(ir,iphi,iz,irho) * nu
         
         ! Yoshizawa model for turbulent kinetic energy.  Eq. 10 in Moin et al.
         q2 = 2.d0 * C_Yoshizawa * q(ir,iphi,iz,irho) * l_grid_squared(ir,iz) * S_abs**2

         ! Heat flux vector.  See notes of Nov. 9, 2018:
         if (.not. isothermal) then
            k_over_cv = mu * gamma_over_Pr_t
            do icomp = 1, 3
               heat_flux(ir,iphi,iz,icomp) = k_over_cv * heat_flux(ir,iphi,iz,icomp)
            end do
            lambda_max = MAX(lambda_max, pi*nu*gamma_over_Pr_t/min_grid_size(ir,iz)**2)                              

            ! Viscous heating:
            heat_term(ir,iphi,iz) = 2.d0*mu*(S_ij_squared - 1.d0/3.d0*S_mm**2)
         end if         

         ! Stresses:
         ! Note: We are overwriting the strain tensor.  The minus sign compared to Moin's formula
         ! accounts for the fact that we write +d/dx_l tau_kl on the RHS.
         T_zz(ir,iphi,iz) = -1.d0/3.d0*q2 + 2.d0*mu *(T_zz(ir, iphi, iz) - 1.d0/3.d0*S_mm)
         T_rr(ir,iphi,iz) = -1.d0/3.d0*q2 + 2.d0*mu *(T_rr(ir, iphi, iz) - 1.d0/3.d0*S_mm)
         T_pp(ir,iphi,iz) = -1.d0/3.d0*q2 + 2.d0*mu *(T_pp(ir, iphi, iz) - 1.d0/3.d0*S_mm)
         T_zr(ir,iphi,iz) =                 2.d0*mu *T_zr(ir, iphi, iz)
         T_pz(ir,iphi,iz) =                 2.d0*mu *T_pz(ir, iphi, iz)
         T_pr(ir,iphi,iz) =                 2.d0*mu *T_pr(ir, iphi, iz)
      end do
   end do
end do

! We need to do this since we have overwritten the strain tensor:
have_strain_tensor = .false.

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: Moin_etal_stress_and_heat_flux: Returning'
#endif

end subroutine Moin_etal_stress_and_heat_flux

!----------------------------------------------------------------------------------85

subroutine Moin_ddsv_stress_and_heat_flux(q, lambda_max)

use viscous
use partition_data
use grid
use dof_indices
use gravity
use thermal_parameters
use math_constants
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: lambda_max ! eigenvalue for time step determination.

! Locals:
real(8) :: S_ij_squared, S_abs, S_mm, nu, nu1, nu2, mu, q2, k_over_cv
integer :: icomp, iz, iphi, ir

if (.not. have_strain_tensor) call strain_tensor(q)
if (.not. isothermal        ) call grad_cvT     (q)
lambda_max = 0.0d0

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         S_ij_squared = 2.d0*T_zr(ir,iphi,iz)**2 + 2.d0*T_pz(ir,iphi,iz)**2 + 2.d0*T_pr(ir,iphi,iz)**2 + &
                             T_zz(ir,iphi,iz)**2 +      T_rr(ir,iphi,iz)**2 +      T_pp(ir,iphi,iz)**2
         S_abs = SQRT(2.0d0 * S_ij_squared)
         S_mm = T_zz(ir,iphi,iz) + T_rr(ir,iphi,iz) + T_pp(ir,iphi,iz) ! dilatation

         ! Smagorinsky model:
         nu1 = C_Smag * l_grid_squared(ir,iz) * S_abs ! Smagorinsky model

         ! DDSV:
         nu2 = C_DDSV*l_grid_squared(ir,iz) * abs(min(S_mm, 0.0d0))
         nu = nu1 + nu2
         
         lambda_max = MAX(lambda_max, pi*nu/min_grid_size(ir,iz)**2)                  
         nu_plot(ir,iphi,iz) = nu
         mu = q(ir,iphi,iz,irho) * nu
         
         ! Yoshizawa model for turbulent kinetic energy.  Eq. 10 in Moin et al.
         q2 = 2.d0 * C_Yoshizawa * q(ir,iphi,iz,irho) * l_grid_squared(ir,iz) * S_abs**2

         ! Heat flux vector.  See notes of Nov. 9, 2018:
         if (.not. isothermal) then
            k_over_cv = mu * gamma_over_Pr_t
            do icomp = 1, 3
               heat_flux(ir,iphi,iz,icomp) = k_over_cv * heat_flux(ir,iphi,iz,icomp)
            end do
            lambda_max = MAX(lambda_max, pi*nu*gamma_over_Pr_t/min_grid_size(ir,iz)**2)                              

            ! Viscous heating:
            heat_term(ir,iphi,iz) = 2.d0*mu*(S_ij_squared - 1.d0/3.d0*S_mm**2)
         end if         

         ! Stresses:
         ! Note: We are overwriting the strain tensor.  The minus sign compared to Moin's formula
         ! accounts for the fact that we write +d/dx_l tau_kl on the RHS.
         T_zz(ir,iphi,iz) = -1.d0/3.d0*q2 + 2.d0*mu *(T_zz(ir, iphi, iz) - 1.d0/3.d0*S_mm)
         T_rr(ir,iphi,iz) = -1.d0/3.d0*q2 + 2.d0*mu *(T_rr(ir, iphi, iz) - 1.d0/3.d0*S_mm)
         T_pp(ir,iphi,iz) = -1.d0/3.d0*q2 + 2.d0*mu *(T_pp(ir, iphi, iz) - 1.d0/3.d0*S_mm)
         T_zr(ir,iphi,iz) =                 2.d0*mu *T_zr(ir, iphi, iz)
         T_pz(ir,iphi,iz) =                 2.d0*mu *T_pz(ir, iphi, iz)
         T_pr(ir,iphi,iz) =                 2.d0*mu *T_pr(ir, iphi, iz)
      end do
   end do
end do

! We need to do this since we have overwritten the strain tensor:
have_strain_tensor = .false.

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: Moin_ddsv_stress_and_heat_flux: Returning'
#endif

end subroutine Moin_ddsv_stress_and_heat_flux

!----------------------------------------------------------------------------------85

subroutine Vreman_stress_and_heat_flux(q, lambda_max)

! Adds the stress, heat-flux, and viscous heating from Vreman's model. 

use viscous
use partition_data
use grid
use dof_indices
use gravity
use thermal_parameters
use math_constants
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: lambda_max

! Locals:
real(8) :: S_zz, S_rr, S_pp, S_zr, S_pz, S_pr
real(8) :: square, S_abs, S_mm, nu, mu, q2, alpha_squared, B_beta, beta(3,3), k_over_cv
integer :: iz, iphi, ir, i, j, m, icomp
real(8) :: C_Vreman

C_Vreman = 2.5d0 * C_Smag**2

if (.not. have_velocity_gradient_tensor) call velocity_gradient_tensor(q)
if (.not. isothermal                   ) call grad_cvT(q)

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         alpha_squared = 0.0d0
         do j = 1, 3
            do i = 1, 3
               beta(i, j) = 0.0d0
               do m = 1, 3
                  beta(i,j) = beta(i,j) + Delta_vector(m,ir,iz)*G(ir,iphi,iz,m,i)*G(ir,iphi,iz,m,j)
               end do
               alpha_squared = alpha_squared + G(ir,iphi,iz,i,j)
            end do
         end do

         ! Eq. (8) in Vreman:
         B_beta = beta(1,1)*beta(2,2) - beta(1,2)**2 + beta(1,1)*beta(3,3) - beta(1,3)**2 + &
              beta(2,2)*beta(3,3) - beta(2,3)**2

         nu = C_Vreman * (B_beta/alpha_squared)
         mu = q(ir,iphi,iz,irho) * nu
         nu_plot(ir,iphi,iz) = nu
         lambda_max = MAX(lambda_max, pi*nu/min_grid_size(ir,iz)**2)         

         ! S_ij
         S_zz = G(ir,iphi,iz,z_comp,z_comp)
         S_rr = G(ir,iphi,iz,r_comp,r_comp)
         S_pp = G(ir,iphi,iz,p_comp,p_comp)
         S_zr = 0.5d0*(G(ir,iphi,iz,z_comp,r_comp) + G(ir,iphi,iz,r_comp,z_comp))
         S_pz = 0.5d0*(G(ir,iphi,iz,p_comp,z_comp) + G(ir,iphi,iz,z_comp,p_comp))
         S_pr = 0.5d0*(G(ir,iphi,iz,p_comp,r_comp) + G(ir,iphi,iz,r_comp,p_comp))                  

         ! Sij^2:
         square = 2.d0*S_zr**2 + 2.d0*S_pz**2 + 2.d0*S_pr**2 + &
                       S_zz**2 +      S_rr**2 +      S_pp**2
         S_mm = S_zz + S_rr + S_pp
         S_abs = SQRT(2.0d0 * square)

         ! Yoshizawa model:
         q2 = 2.d0 * C_Yoshizawa * q(ir,iphi,iz,irho) * l_grid_squared(ir,iz) * S_abs**2 ! eq. 10

         T_zz(ir,iphi,iz) = - 1.d0/3.d0*q2 + 2.d0*mu *(S_zz - 1.d0/3.d0*S_mm)
         T_rr(ir,iphi,iz) = - 1.d0/3.d0*q2 + 2.d0*mu *(S_rr - 1.d0/3.d0*S_mm)
         T_pp(ir,iphi,iz) = - 1.d0/3.d0*q2 + 2.d0*mu *(S_pp - 1.d0/3.d0*S_mm)
         T_zr(ir,iphi,iz) =                  2.d0*mu *S_zr
         T_pz(ir,iphi,iz) =                  2.d0*mu *S_pz
         T_pr(ir,iphi,iz) =                  2.d0*mu *S_pr

         ! Heat flux vector.  See notes of Nov. 9, 2018:
         if (.not. isothermal) then
            k_over_cv = mu * gamma_over_Pr_t ! k is the conductivity
            do icomp = 1, 3
               heat_flux(ir,iphi,iz,icomp) = k_over_cv * heat_flux(ir,iphi,iz,icomp)
            end do
            lambda_max = MAX(lambda_max, pi*nu*gamma_over_Pr_t/min_grid_size(ir,iz)**2)

            ! Viscous heating term:
            heat_term(ir,iphi,iz) = heat_term(ir,iphi,iz) + 2.d0*mu*(square - 1.d0/3.d0*S_mm**2)
         end if
      end do
   end do
end do

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: Vreman_sgs_stress_and_heat_flux: Returning'
#endif

end subroutine Vreman_stress_and_heat_flux

!----------------------------------------------------------------------------------85

subroutine Vreman_ddsv_stress_and_heat_flux(q, lambda_max)

! Adds the stress, heat-flux, and viscous heating from Vreman's model. 

use viscous
use partition_data
use grid
use dof_indices
use gravity
use thermal_parameters
use math_constants
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: lambda_max

! Locals:
real(8) :: S_zz, S_rr, S_pp, S_zr, S_pz, S_pr
real(8) :: square, S_abs, S_mm, q2, alpha_squared, B_beta, beta(3,3), k_over_cv
integer :: iz, iphi, ir, i, j, m, icomp
real(8) :: C_Vreman
real(8) :: nu1, nu2, nu, mu

C_Vreman = 2.5d0 * C_Smag**2

if (.not. have_velocity_gradient_tensor) call velocity_gradient_tensor(q)
if (.not. isothermal                   ) call grad_cvT(q)

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         alpha_squared = 0.0d0
         do j = 1, 3
            do i = 1, 3
               beta(i, j) = 0.0d0
               do m = 1, 3
                  beta(i,j) = beta(i,j) + Delta_vector(m,ir,iz)*G(ir,iphi,iz,m,i)*G(ir,iphi,iz,m,j)
               end do
               alpha_squared = alpha_squared + G(ir,iphi,iz,i,j)
            end do
         end do

         ! Eq. (8) in Vreman:
         B_beta = beta(1,1)*beta(2,2) - beta(1,2)**2 + beta(1,1)*beta(3,3) - beta(1,3)**2 + &
              beta(2,2)*beta(3,3) - beta(2,3)**2

         nu1 = C_Vreman * (B_beta/alpha_squared)
         nu2 = C_DDSV*l_grid_squared(ir,iz) * abs(min(S_mm, 0.0d0)) ! DDSV
         nu = nu1 + nu2

         mu = q(ir,iphi,iz,irho) * nu
         nu_plot(ir,iphi,iz) = nu
         lambda_max = MAX(lambda_max, pi*nu/min_grid_size(ir,iz)**2)         

         ! S_ij
         S_zz = G(ir,iphi,iz,z_comp,z_comp)
         S_rr = G(ir,iphi,iz,r_comp,r_comp)
         S_pp = G(ir,iphi,iz,p_comp,p_comp)
         S_zr = 0.5d0*(G(ir,iphi,iz,z_comp,r_comp) + G(ir,iphi,iz,r_comp,z_comp))
         S_pz = 0.5d0*(G(ir,iphi,iz,p_comp,z_comp) + G(ir,iphi,iz,z_comp,p_comp))
         S_pr = 0.5d0*(G(ir,iphi,iz,p_comp,r_comp) + G(ir,iphi,iz,r_comp,p_comp))                  

         ! Sij^2:
         square = 2.d0*S_zr**2 + 2.d0*S_pz**2 + 2.d0*S_pr**2 + &
                       S_zz**2 +      S_rr**2 +      S_pp**2
         S_mm = S_zz + S_rr + S_pp
         S_abs = SQRT(2.0d0 * square)

         ! Yoshizawa model:
         q2 = 2.d0 * C_Yoshizawa * q(ir,iphi,iz,irho) * l_grid_squared(ir,iz) * S_abs**2 ! eq. 10

         T_zz(ir,iphi,iz) = - 1.d0/3.d0*q2 + 2.d0*mu *(S_zz - 1.d0/3.d0*S_mm)
         T_rr(ir,iphi,iz) = - 1.d0/3.d0*q2 + 2.d0*mu *(S_rr - 1.d0/3.d0*S_mm)
         T_pp(ir,iphi,iz) = - 1.d0/3.d0*q2 + 2.d0*mu *(S_pp - 1.d0/3.d0*S_mm)
         T_zr(ir,iphi,iz) =                  2.d0*mu *S_zr
         T_pz(ir,iphi,iz) =                  2.d0*mu *S_pz
         T_pr(ir,iphi,iz) =                  2.d0*mu *S_pr

         ! Heat flux vector.  See notes of Nov. 9, 2018:
         if (.not. isothermal) then
            k_over_cv = mu * gamma_over_Pr_t ! k is the conductivity
            do icomp = 1, 3
               heat_flux(ir,iphi,iz,icomp) = k_over_cv * heat_flux(ir,iphi,iz,icomp)
            end do
            lambda_max = MAX(lambda_max, pi*nu*gamma_over_Pr_t/min_grid_size(ir,iz)**2)

            ! Viscous heating term:
            heat_term(ir,iphi,iz) = heat_term(ir,iphi,iz) + 2.d0*mu*(square - 1.d0/3.d0*S_mm**2)
         end if
      end do
   end do
end do

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: Vreman_sgs_stress_and_heat_flux: Returning'
#endif

end subroutine Vreman_ddsv_stress_and_heat_flux

!----------------------------------------------------------------------------------85

subroutine zero_stress_bc(q)

! Zero stress conditions.  See notes of Dec. 11, 2018.  This is necessary for
! compatibility with the zero normal momentum condition. 

use partition_data
use boundary_condition_data
use boundary_condition_types
use dof_indices
use grid
use viscous
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: lambda_max ! eigenvalue for time step determination.

! Locals:
real(8) :: rho, ur, uz, uphi
real(8) :: square, S_abs, S_mm, mu, nu, q2
integer :: iz, iphi, ir
! For zero stress condition:
real(8) :: u_1, u_2, u_3, u_4, u_n, u_nm1, u_nm2, u_nm3, divisor, dur_dr, duz_dz

if ((sr .eq. 1) .and. (rmin_BC .eq. zero_normal_momentum)) then
   ir = 1
   do iz = 1, nz
      do iphi = sphi, ephi   
         ! (1) T_rr = 0 gives a condition on ur at ir = 1:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         square = 2.d0*T_zr(ir,iphi,iz)**2 + 2.d0*T_pz(ir,iphi,iz)**2 + 2.d0*T_pr(ir,iphi,iz)**2 + &
                       T_zz(ir,iphi,iz)**2 +      T_rr(ir,iphi,iz)**2 +      T_pp(ir,iphi,iz)**2
         S_mm = T_zz(ir,iphi,iz) + T_rr(ir,iphi,iz) + T_pp(ir,iphi,iz)
         S_abs = SQRT(2.0d0 * square)
         nu = C_Smag * l_grid_squared(ir,iz) * S_abs
         mu = nu * q(ir,iphi,iz,irho)
         q2 = 2.d0 * C_Yoshizawa * q(ir,iphi,iz,irho) * l_grid_squared(ir,iz) * S_abs**2 ! eq. 10

         dur_dr = 1.d0/6.d0*q2/mu + 1.d0/3.d0*S_mm
         
         u_2 = q(2, iphi, iz, zmom) / q(2, iphi, iz, irho) 
         u_3 = q(3, iphi, iz, zmom) / q(3, iphi, iz, irho) 
         u_4 = q(4, iphi, iz, zmom) / q(4, iphi, iz, irho)

         u_1 = 1.d0/11.d0 * (18.d0*u_2 - 9.d0*u_3 + 2.d0*u_4 - 6.d0*dur_dr/Ji_r(ir))
         q(1, iphi, iz, rmom) = q(1, iphi, iz, irho) * u_1         

         ! (2) T_rp = 0 gives a condition on uphi at ir = 1:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         u_2 = q(2, iphi, iz, amom) / q(2, iphi, iz, irho) / rgrid(2) 
         u_3 = q(3, iphi, iz, amom) / q(3, iphi, iz, irho) / rgrid(3)
         u_4 = q(4, iphi, iz, amom) / q(4, iphi, iz, irho) / rgrid(4)
         divisor = 11.d0/6.d0*Ji_r(ir) + 1.d0/rgrid(ir)
         u_1 = (18.d0*u_2 - 9.d0*u_3 + 2.d0*u_4) / 6.d0 * Ji_r(ir) / divisor 
         q(ir, iphi, iz, amom) = q(ir, iphi, iz, irho) * u_1 * rgrid(ir)         

         ! (3) T_rz = 0 gives a condition on uz at ir = 1:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         u_2 = q(2, iphi, iz, zmom) / q(2, iphi, iz, irho) 
         u_3 = q(3, iphi, iz, zmom) / q(3, iphi, iz, irho) 
         u_4 = q(4, iphi, iz, zmom) / q(4, iphi, iz, irho)
         u_1 = (18.d0*u_2 - 9.d0*u_3 + 2.d0*u_4) / 11.d0
         q(1, iphi, iz, zmom) = q(1, iphi, iz, irho) * u_1
      end do
   end do
end if

if ((er .eq. nr) .and. (rmax_BC .eq. zero_normal_momentum)) then
   ir = nr
   do iz = 1, nz
      do iphi = sphi, ephi   
         ! (4) T_rr = 0 gives a condition on ur at ir = nr:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Sij^2:
         square = 2.d0*T_zr(ir,iphi,iz)**2 + 2.d0*T_pz(ir,iphi,iz)**2 + 2.d0*T_pr(ir,iphi,iz)**2 + &
                       T_zz(ir,iphi,iz)**2 +      T_rr(ir,iphi,iz)**2 +      T_pp(ir,iphi,iz)**2
         S_mm = T_zz(ir,iphi,iz) + T_rr(ir,iphi,iz) + T_pp(ir,iphi,iz)
         S_abs = SQRT(2.0d0 * square)
               ! Subgrid viscosity:
         nu = C_Smag * l_grid_squared(ir,iz) * S_abs
         mu = nu * q(ir,iphi,iz,irho)
         q2 = 2.d0 * C_Yoshizawa * q(ir,iphi,iz,irho) * l_grid_squared(ir,iz) * S_abs**2 ! eq. 10

         dur_dr = 1.d0/6.d0*q2/mu + 1.d0/3.d0*S_mm
         
         u_nm1 = q(nr-1, iphi, iz, rmom) / q(nr-1, iphi, iz, irho) 
         u_nm2 = q(nr-2, iphi, iz, rmom) / q(nr-2, iphi, iz, irho) 
         u_nm3 = q(nr-3, iphi, iz, rmom) / q(nr-3, iphi, iz, irho)

         u_n = 1.d0/11.d0 * (18.d0*u_nm1 - 9.d0*u_nm1 + 2.d0*u_nm3 + 6.d0*dur_dr/Ji_r(nr))
         q(nr, iphi, iz, rmom) = q(nr, iphi, iz, irho) * u_n

         ! (5) T_rp = 0 gives a condition on uphi at ir = nr:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         u_nm1 = q(nr-1, iphi, iz, amom) / q(nr-1, iphi, iz, irho) / rgrid(nr-1) 
         u_nm2 = q(nr-2, iphi, iz, amom) / q(nr-2, iphi, iz, irho) / rgrid(nr-2)
         u_nm3 = q(nr-3, iphi, iz, amom) / q(nr-3, iphi, iz, irho) / rgrid(nr-3)
         divisor = 11.d0/6.d0*Ji_r(nr) - 1.d0/rgrid(nr)
         u_n = (18.d0*u_nm1 - 9.d0*u_nm2 + 2.d0*u_nm3) / 6.d0 * Ji_r(nr) / divisor 
         q(nr, iphi, iz, amom) = q(nr, iphi, iz, irho) * u_n * rgrid(nr)         

         ! (6) T_rz = 0 gives a condition on uz at ir = nr:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         
         u_nm1 = q(nr-1, iphi, iz, zmom) / q(nr-1, iphi, iz, irho) 
         u_nm2 = q(nr-2, iphi, iz, zmom) / q(nr-2, iphi, iz, irho) 
         u_nm3 = q(nr-3, iphi, iz, zmom) / q(nr-3, iphi, iz, irho)
         u_n = (18.d0*u_nm1 - 9.d0*u_nm2 + 2.d0*u_nm3) / 11.d0
         q(nr, iphi, iz, zmom) = q(nr, iphi, iz, irho) * u_n
      end do
   end do
end if

if (zmin_BC .eq. zero_normal_momentum) then
   iz = 1
   do iphi = sphi, ephi
      do ir = sr, er
         ! (7) T_zz = 0 at iz = 1 gives a condition on uz at iz = 1:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         square = 2.d0*T_zr(ir,iphi,iz)**2 + 2.d0*T_pz(ir,iphi,iz)**2 + 2.d0*T_pr(ir,iphi,iz)**2 + &
                       T_zz(ir,iphi,iz)**2 +      T_rr(ir,iphi,iz)**2 +      T_pp(ir,iphi,iz)**2
         S_mm = T_zz(ir,iphi,iz) + T_rr(ir,iphi,iz) + T_pp(ir,iphi,iz)
         S_abs = SQRT(2.0d0 * square)
         nu = C_Smag * l_grid_squared(ir,iz) * S_abs
         mu = nu * q(ir,iphi,iz,irho)
         q2 = 2.d0 * C_Yoshizawa * q(ir,iphi,iz,irho) * l_grid_squared(ir,iz) * S_abs**2 ! eq. 10

         duz_dz = 1.d0/6.d0*q2/mu + 1.d0/3.d0*S_mm
         
         u_2 = q(ir, iphi, 2, zmom) / q(ir, iphi, 2, irho) 
         u_3 = q(ir, iphi, 3, zmom) / q(ir, iphi, 3, irho) 
         u_4 = q(ir, iphi, 4, zmom) / q(ir, iphi, 4, irho)

         u_1 = 1.d0/11.d0 * (18.d0*u_2 - 9.d0*u_3 + 2.d0*u_4 - 6.d0*duz_dz/Ji_z(iz))
         q(ir, iphi, iz, zmom) = q(ir, iphi, iz, irho) * u_1

         ! (8) T_pz = 0 at iz = 1 gives a condition on uphi at iz = 1:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         u_2 = q(ir, iphi, 2, amom) / q(ir, iphi, 2, irho) / rgrid(ir) 
         u_3 = q(ir, iphi, 3, amom) / q(ir, iphi, 3, irho) / rgrid(ir)
         u_4 = q(ir, iphi, 4, amom) / q(ir, iphi, 4, irho) / rgrid(ir)

         u_1 = 1.d0/11.d0 * (18.d0*u_2 - 9.d0*u_3 + 2.d0*u_4)         
         q(ir, iphi, iz, amom) = q(ir, iphi, iz, irho) * u_1 * rgrid(ir)

         ! (9) T_zr = 0 at iz = 1 gives a condition on ur at iz = 1:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         u_2 = q(ir, iphi, 2, rmom) / q(ir, iphi, 2, irho) 
         u_3 = q(ir, iphi, 3, rmom) / q(ir, iphi, 3, irho) 
         u_4 = q(ir, iphi, 4, rmom) / q(ir, iphi, 4, irho)

         u_1 = 1.d0/11.d0 * (18.d0*u_2 - 9.d0*u_3 + 2.d0*u_4)
         q(ir, iphi, iz, rmom) = q(ir, iphi, iz, irho) * u_1         
      end do
   end do
end if

if (zmax_BC .eq. zero_normal_momentum) then
   iz = nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! (10) T_zz = 0 at iz = nz gives a condition on uz at iz = nz:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         square = 2.d0*T_zr(ir,iphi,iz)**2 + 2.d0*T_pz(ir,iphi,iz)**2 + 2.d0*T_pr(ir,iphi,iz)**2 + &
                       T_zz(ir,iphi,iz)**2 +      T_rr(ir,iphi,iz)**2 +      T_pp(ir,iphi,iz)**2
         S_mm = T_zz(ir,iphi,iz) + T_rr(ir,iphi,iz) + T_pp(ir,iphi,iz)
         S_abs = SQRT(2.0d0 * square)
         nu = C_Smag * l_grid_squared(ir,iz)**2 * S_abs
         mu = nu * q(ir,iphi,iz,irho)
         q2 = 2.d0 * C_Yoshizawa * q(ir,iphi,iz,irho) * l_grid_squared(ir,iz) * S_abs**2 ! eq. 10

         duz_dz = 1.d0/6.d0*q2/mu + 1.d0/3.d0*S_mm
         
         u_nm1 = q(ir, iphi, nz-1, zmom) / q(ir, iphi, nz-1, irho) 
         u_nm2 = q(ir, iphi, nz-2, zmom) / q(ir, iphi, nz-2, irho) 
         u_nm3 = q(ir, iphi, nz-3, zmom) / q(ir, iphi, nz-3, irho)

         u_n = 1.d0/11.d0 * (18.d0*u_nm1 - 9.d0*u_nm2 + 2.d0*u_nm3 + 6.d0*duz_dz/Ji_z(iz))
         q(ir, iphi, iz, zmom) = q(ir, iphi, iz, irho) * u_n

         ! (11) T_pz = 0 at iz = nz gives a condition on uphi at iz = nz:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         u_nm1 = q(ir, iphi, nz-1, amom) / q(ir, iphi, nz-1, irho) / rgrid(ir) 
         u_nm2 = q(ir, iphi, nz-2, amom) / q(ir, iphi, nz-2, irho) / rgrid(ir)
         u_nm3 = q(ir, iphi, nz-3, amom) / q(ir, iphi, nz-3, irho) / rgrid(ir)

         u_n = 1.d0/11.d0 * (18.d0*u_nm1 - 9.d0*u_nm2 + 2.d0*u_nm3)         
         q(ir, iphi, iz, amom) = q(ir, iphi, iz, irho) * u_n * rgrid(ir)

         ! (12) T_zr = 0 at iz = nz gives a condition on ur at iz = nz:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         u_nm1 = q(ir, iphi, nz-1, rmom) / q(ir, iphi, nz-1, irho)
         u_nm2 = q(ir, iphi, nz-2, rmom) / q(ir, iphi, nz-2, irho)
         u_nm3 = q(ir, iphi, nz-3, rmom) / q(ir, iphi, nz-3, irho)

         u_n = 1.d0/11.d0 * (18.d0*u_nm1 - 9.d0*u_nm2 + 2.d0*u_nm3)         
         q(ir, iphi, iz, rmom) = q(ir, iphi, iz, irho) * u_n
      end do
   end do
end if

! Debug: Recompute the strain-tensor:
!call strain_tensor(q)

end subroutine zero_stress_bc

!----------------------------------------------------------------------------------85

subroutine strain_tensor(q)

! Here I have 12 transposes.  Previously, I had ????? transposes but I computed
! the velocity again while in each space.  Hence I am going with the direct computation.
! The old coding will be stored as: VARIOUS_CODES/old_strain_rate_computation.f90

! This routine also allows me to insert the changes for Fargo.

use partition_data
use dof_indices
use grid
use viscous
use fargo_or_plotting_shift
use basic_state

implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
logical, parameter :: subtract_basic_state_uphi = .false.

! Local:
real(8), dimension(sr:er, sphi:ephi, nz) :: uz, ur, uphi
! These are dimensioned in z-space but will have data in other spaces since pencils in
! all the spaces have the same size:
real(8), dimension(sr:er, sphi:ephi, nz) :: transposed, deriv
real(8), dimension(sr:er, sphi:ephi, nz) :: deriv_in_z_space, T_cv

integer :: ir, iphi, iz

! Velocities:
if (subtract_basic_state_uphi) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            ur  (ir,iphi,iz) = q(ir,iphi,iz,rmom) / q(ir,iphi,iz,irho)
            uz  (ir,iphi,iz) = q(ir,iphi,iz,zmom) / q(ir,iphi,iz,irho)
            uphi(ir,iphi,iz) = q(ir,iphi,iz,amom) / q(ir,iphi,iz,irho) / rgrid(ir) - uphi_basic(ir, iz)         
         end do
      end do
   end do
else
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            ur  (ir,iphi,iz) = q(ir,iphi,iz,rmom) / q(ir,iphi,iz,irho)
            uz  (ir,iphi,iz) = q(ir,iphi,iz,zmom) / q(ir,iphi,iz,irho)
            uphi(ir,iphi,iz) = q(ir,iphi,iz,amom) / q(ir,iphi,iz,irho) / rgrid(ir)
         end do
      end do
   end do
end if
   

T_zz = 0.0d0; T_zr = 0.0d0; T_rr = 0.0d0; T_pr = 0.0d0; T_pp = 0.0d0; T_pz = 0.0d0

! The item numbers refer to your notes of 9/21/18: Strain tensor computation
! See also notes of April 23, 2021.

! --------------
! r derivatives:
! --------------
if (nr .ne. 1) then
   ! Item 4: duz/dr for T_zr.  Fargo is done in item 12 below and then
   ! T_zr completed in item 2:
   call transpose_z_to_r (1, uz, transposed) ! transpose 1
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, deriv_in_z_space) ! transpose 2
   T_zr = deriv_in_z_space ! T_zr partial.  Fargo still needs to be done.

   ! Item 5: dur/dr = T_rr partial.  Its Fargo is done in item 8 below
   ! which will complete T_rr:
   call transpose_z_to_r(1, ur, transposed) ! transpose 3
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, T_rr) ! transpose 4

   ! Item 6: duphi/dr = First term of T_pr.  Its Fargo is done in item 10 below.
   call transpose_z_to_r(1, uphi, transposed) ! transpose 5
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, deriv_in_z_space) ! transpose 6
   T_pr = deriv_in_z_space ! needs Fargo.  Fargo for this is done in item 10
end if

! ----------------
! phi derivatives:
! ----------------
if (nphi .ne. 1) then
   ! dur/dphi
   call transpose_z_to_phi(1, ur, transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z(1, deriv, deriv_in_z_space)

   ! Item 8: Complete Fargo for T_rr:
   if (add_fargo_extra_operator_now) then
      do ir = sr, er
         T_rr(ir,:,:) = T_rr(ir,:,:) + fargo_factor(ir)*deriv_in_z_space(ir,:,:) ! T_rr complete
      end do
   end if

   ! Item 7: Second term of T_pr
   do ir = sr, er                         ! dur/dphi
      T_pr(ir,:,:) = T_pr(ir,:,:) + deriv_in_z_space(ir,:,:)/rgrid(ir)
   end do

   ! Item 9: duphi/dphi
   call transpose_z_to_phi(1, uphi, transposed) ! 9
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z (1, deriv, deriv_in_z_space) ! 10th transpose

   do ir = sr, er
      T_pp(ir,:,:) = deriv_in_z_space(ir,:,:) ! division by r later
   end do

   ! Item 10: Fargo for item 6: Third term of T_pr (Fargo term).
   if (add_fargo_extra_operator_now) then
      do ir = sr, er
         ! This completes T_pr apart from an undifferentiated term and factor of 1/2:
         ! Note: deriv_in_z_space = duphi/dphi:
         T_pr(ir,:,:) = T_pr(ir,:,:) + fargo_factor(ir)*deriv_in_z_space(ir,:,:)
      end do
   end if

   ! Item 11: duz/dphi:
   call transpose_z_to_phi(1, uz, transposed) ! 11
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z (1, deriv, deriv_in_z_space) ! 12th transpose
   do ir = sr, er
      ! 1/r duz/dphi:
      T_pz(ir,:,:) = deriv_in_z_space(ir,:,:)/rgrid(ir)
   end do

   ! Item 12:  Adjust item 4 with Fargo since we have duz/dphi available in deriv_in_z_space:
   if (add_fargo_extra_operator_now) then
      do ir = sr, er
         T_zr(ir,:,:) = T_zr(ir,:,:) + fargo_factor(ir)*deriv_in_z_space(ir,:,:)
      end do
   end if   
end if
! -------------------------
! Done with phi derivatives
! -------------------------

! --------------
! z derivatives:
! --------------
if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
   ! Item 1: duz/dz = T_zz complete
   call pade_diff_z(nbundle_z, uz,   T_zz)
   ! Item 2: dur/dz gets added to T_zr and we get a factor of 0.5d0:
   call pade_diff_z(nbundle_z, ur, deriv_in_z_space) 
   T_zr = 0.5d0*(T_zr + deriv_in_z_space) ! Done
   ! Item 3: 0.5 duphi/dz get's added to T_pz:
   call pade_diff_z(nbundle_z, uphi, deriv_in_z_space)
   T_pz = 0.5d0*(T_pz + deriv_in_z_space)
end if

! -----------------------
! Undifferentiated terms:
! -----------------------
do ir = sr, er
   ! Item 13.  Final term of T_pr
   T_pr(ir,:,:) = 0.5d0*(T_pr(ir,:,:) - uphi(ir,:,:)/rgrid(ir))
   ! Item 14.  Coming here from item 9:
   T_pp(ir,:,:) = (T_pp(ir,:,:) + ur(ir,:,:)) / rgrid(ir)
end do

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: strain_tensor.  Returning'
#endif

end subroutine strain_tensor

!----------------------------------------------------------------------------------85

subroutine grad_cvT(q)

! Computes gradient of Cv*T and puts in into heat_flux

use partition_data
use dof_indices
use grid
use viscous
use fargo_or_plotting_shift
use thermal_parameters

implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
real(8), dimension(sr:er, sphi:ephi, nz) :: uz, ur, uphi
! These are dimensioned in z-space but will have data in other spaces since pencils in
! all the spaces have the same size:
real(8), dimension(sr:er, sphi:ephi, nz) :: transposed, deriv, cvT
integer :: iz, iphi, ir

cvT(:,:,:) = q(:,:,:,ener) / q(:,:,:,irho)
   
if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
   call pade_diff_z(nbundle_z, cvT(sr,sphi,1), heat_flux(sr,sphi,1,z_comp))
end if

if (nr .ne. 1) then
   call transpose_z_to_r(1, cvT(sr,sphi,1), transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z(1, deriv, heat_flux(sr,sphi,1,r_comp))
end if

if (nphi .ne. 1) then
   call transpose_z_to_phi(1, cvT(sr,sphi,1), transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z(1, deriv, heat_flux(sr,sphi,1,p_comp))

   ! Fargo for the r-derivative.  grad_rho_phi is d rho/dphi at the moment:
   if (add_fargo_extra_operator_now) then
      do iz = 1, nz   
         do iphi = sphi, ephi     
            do ir = sr, er
               heat_flux(ir,iphi,iz,r_comp) = heat_flux(ir,iphi,iz,r_comp) + fargo_factor(ir)*heat_flux(ir,iphi,iz,p_comp)
            end do
         end do
      end do
   end if
      
   ! Put in the r factor in the phi component:
   do iz = 1, nz   
      do iphi = sphi, ephi     
         do ir = sr, er
            heat_flux(ir,iphi,iz,p_comp) = heat_flux(ir,iphi,iz,p_comp)/rgrid(ir)
         end do
      end do
   end do      
end if

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: grad_CvT.  Returning'
#endif

end subroutine grad_cvT

!----------------------------------------------------------------------------------85

subroutine velocity_gradient_tensor(q)

! See notes dated Dec. 18, 2018 entitled "Velocity gradient tensor in cylindrical
! coordinates."

! This routine is used for the Vreman model (Phys. Fluids, v.16, p.3670, 2004).

use partition_data
use dof_indices
use grid
use viscous
use fargo_or_plotting_shift
use basic_state

implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
real(8), dimension(sr:er, sphi:ephi, nz) :: uz, ur, uphi
! These are dimensioned in z-space but will have data in other spaces since pencils in
! all the spaces have the same size:
real(8), dimension(sr:er, sphi:ephi, nz) :: transposed, deriv
real(8), dimension(sr:er, sphi:ephi, nz) :: deriv_in_z_space

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

! These is necessary because we could skip some sections below, for example,
! if the flow lacks a certain directional dependence.
G = 0.0d0

! --------------
! r derivatives:
! --------------
if (nr .ne. 1) then
   ! (1) G_zr.  Fargo is done below.
   call transpose_z_to_r (1, uz, transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, G(sr,sphi,1,z_comp,r_comp))

   ! (2) G_rr.  Fargo is done below.
   call transpose_z_to_r(1, ur, transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, G(sr,sphi,1,r_comp,r_comp))

   ! (3) G_pr.  Fargo is done below.
   call transpose_z_to_r(1, uphi, transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, G(sr,sphi,1,p_comp,r_comp))
end if

! ----------------
! phi derivatives:
! ----------------
if (nphi .ne. 1) then
   ! duz/dphi:
   call transpose_z_to_phi(1, uz, transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z (1, deriv, deriv_in_z_space)

   ! (4) Fargo for G_zr:
   if (add_fargo_extra_operator_now) then
      do ir = sr, er
         G(ir,:,:,z_comp,r_comp) = G(ir,:,:,z_comp,r_comp) + fargo_factor(ir)*deriv_in_z_space(ir,:,:)
      end do
   end if

   ! (5) G_zp:
   do ir = sr, er
      G(ir,:,:,z_comp,p_comp) = deriv_in_z_space(ir,:,:)/rgrid(ir)
   end do

   ! dur/dphi
   call transpose_z_to_phi(1, ur, transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z(1, deriv, deriv_in_z_space)

   ! (6) Fargo for G_rr:
   if (add_fargo_extra_operator_now) then
      do ir = sr, er
         G(ir,:,:,r_comp,r_comp) = G(ir,:,:,r_comp,r_comp) + fargo_factor(ir)*deriv_in_z_space(ir,:,:) ! G_rr complete
      end do
   end if

   ! (7) G_rp:
   do ir = sr, er     ! dur/dphi
      G(ir,:,:,r_comp,p_comp) = (deriv_in_z_space(ir,:,:) - uphi(ir,:,:))/rgrid(ir)
   end do

   ! duphi/dphi
   call transpose_z_to_phi(1, uphi, transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z (1, deriv, deriv_in_z_space)

   ! (8) Fargo for G_pr:   
   if (add_fargo_extra_operator_now) then
      do ir = sr, er
         G(ir,:,:,p_comp,r_comp) = G(ir,:,:,p_comp,r_comp) + fargo_factor(ir)*deriv_in_z_space(ir,:,:)
      end do
   end if

   ! (9) G_pp:
   do ir = sr, er
      G(ir,:,:,p_comp,p_comp) = (deriv_in_z_space(ir,:,:) + ur(ir,:,:))/rgrid(ir)
   end do
end if
! -------------------------
! Done with phi derivatives
! -------------------------

! --------------
! z derivatives:
! --------------
if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
   ! (10) G_zz:
   call pade_diff_z(nbundle_z, uz,   G(sr,sphi,1,z_comp,z_comp))
   ! (11) G_rz:
   call pade_diff_z(nbundle_z, ur,   G(sr,sphi,1,r_comp,z_comp)) 
   ! (12) G_pz:
   call pade_diff_z(nbundle_z, uphi, G(sr,sphi,1,p_comp,z_comp))
end if

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: velocity_gradient_tensor.  Returning'
#endif

have_velocity_gradient_tensor = .true.

end subroutine velocity_gradient_tensor

!----------------------------------------------------------------------------------85

subroutine tecplot_sgs_terms_in_meridional_plane (iphi, t)

! This subroutine works in either serial or parallel mode.  The direct access
! and specification of the record number for each data line allows each processor
! to write to the proper record.

use dof_indices
use logical_units
use grid
use thermal_parameters
use partition_data
use basic_state
use physical_constants
use viscous
implicit none
integer, intent(in) :: iphi
real(8), intent(in) :: t

! Local:
character(80) :: filename
integer :: ir, iz, irec

write (filename, 4) int(t), t - int(t)
4 format ('sgs_meridional_t_', i6.6, f0.4, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', access = 'direct', recl = 13*7 + 1)

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
   1  format ('TITLE = "t = ', e12.5, '"', a1, a1)
   write (lun_tecplot, 2, rec = 2) char(13), char(10)
   2  format('VARIABLES = "r", "z", "nu_t", "Fr", "Fz", "Fp", "grid min"', a1, a1)
   write (lun_tecplot, 3, rec = 3) nr, nz, char(13), char(10)
   3 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', a1, a1)
end if

do iz = 1, nz
   do ir = sr, er
      irec = (iz - 1)*nr + ir + 3

      write(lun_tecplot, "(7(1x, e12.5), a1)", rec = irec) rgrid(ir), zgrid(iz), &
           nu_plot(ir,iphi,iz), Fz(ir,iphi,iz), Fr(ir,iphi,iz), Fp(ir,iphi,iz), min_grid_size(ir,iz),char(10)
   end do
end do
close(lun_tecplot)

if (my_node .eq. 0) print *, ' node 0: Wrote tecplot file ', filename

end subroutine tecplot_sgs_terms_in_meridional_plane

!----------------------------------------------------------------------------------85

subroutine tecplot_sgs_force_and_heat_term(t)

use viscous
use grid
use xy_coordinates_of_grid
use logical_units, only: lun_tecplot
use partition_data
use dof_indices
implicit none
real(8) :: t 

character(80) :: filename
integer :: ir, iphi, iz, irec, icomp, iphi_rec, nz_plot

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: tecplot_Moin_etal_subgrid_force_and_heat_term has been called'
#endif

if (suppress_z_derivatives_when_nz_not_1) then
   nz_plot = 1
else
   nz_plot = nz
end if

write (filename, 7) int(t), t - int(t)
7 format ('sgs_force_and_heat_term_t_', i6.6, f0.4, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', &
     access = 'direct', recl = 15*7 + 1) ! 3 coordinates + 3 force components + mass term

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
1  format (' TITLE = "t = ', e12.5, '"', a1, a1)

   write (lun_tecplot, 2, rec = 2) char(13), char(10)
2  format(' VARIABLES = "x", "y", "z", "Fz", "Fr", "Fphi", "Q"', a1, a1)
   
   write (lun_tecplot, 6, rec = 3) nr, nphi+1, nz_plot, char(13), char(10)
6  format(' ZONE I=',i4,',',' J=',i4,',',' K=',i4,',',' DATAPACKING=POINT', a1, a1)
end if

do iz = 1, nz_plot
   do iphi = sphi, ephi
      do ir = sr, er
         ! nphi+1 leaves room for periodic completion in phi:
         irec = ir + (iphi-1)*nr + (iz-1)*nr*(nphi+1) + 3  
         write(lun_tecplot, "(7(1x, e14.7), a1)", rec = irec) &
              xgrid(ir,iphi), ygrid(ir,iphi), zgrid(iz), &
              Fz(ir, iphi, iz), Fr(ir, iphi, iz), Fp(ir, iphi, iz), heat_term(ir, iphi, iz), char(10)
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
         write(lun_tecplot, "(7(1x, e14.7), a1)", rec = irec) &
              xgrid(ir,0), ygrid(ir,0), zgrid(iz), &
              Fz(ir, iphi, iz), Fr(ir, iphi, iz), Fp(ir, iphi, iz), heat_term(ir, iphi, iz), char(10)
      end do
   end do
end if
close(lun_tecplot)

if (my_node .eq. 0) print *, ' node 0: Wrote sgs_force_and_heat_term tecplot file ', filename

end subroutine tecplot_sgs_force_and_heat_term

!----------------------------------------------------------------------------------85

