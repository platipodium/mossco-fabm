#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_mpb_cnp --- Fortran 2003 version of microphytobenthos (mpb) model
!
! !INTERFACE:
   module hzg_mpb_cnp
!
! !DESCRIPTION:
!
! The MPB model extracted from Hochard et al EcoMod 2010 adopted by Markus Kreus
! and Kai Wirtz.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_hzg_mpb_cnp
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day  = 1.0_rk / secs_pr_day
!
! !REVISION HISTORY:!
!  Original author(s): Markus Kreus, Richard Hofmeister & Kai Wirtz
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hzg_mpb_cnp
!     Variable identifiers
      type (type_state_variable_id)        :: id_no3, id_nh4, id_po4, id_oxy, id_ldet
      type (type_state_variable_id)        :: id_dic, id_zbC, id_zbN, id_zbP
      type (type_state_variable_id)        :: id_mpbCHL, id_mpbC, id_mpbN, id_mpbP, id_eps
      type (type_dependency_id)            :: id_temp, id_parz, id_porosity
      type (type_diagnostic_variable_id)   :: id_PrimProd, id_par, id_qNC, id_Q_chl, id_ThetaN
      type (type_diagnostic_variable_id)   :: id_expCProd, id_expNProd, id_expPProd
      type (type_diagnostic_variable_id)   :: id_NPP, id_SGR, id_TGR, id_SPR, id_SMR
! temporary for debugging
      type (type_diagnostic_variable_id)   :: id_MPB_din, id_MPB_no3, id_MPB_nh4, id_MPB_po4

!     Model parameters
      real(rk) :: Tref, Q10, alpha, mu_max, resp0, zeta_NO3, zeta_NH4, kPO4, kNO3, kNH4, kInNH4
      real(rk) :: theta_max, QNCmin, QNCmax, QNCupt, QPNmin, QPNmax, QPNupt, sigma_NC, sigma_PN
      real(rk) :: kO2resp, gamma, fracEPS, degrEPS, graz, aeff, exud, rzoo
      logical  :: use_no3, use_nh4, use_po4, use_oxy, use_ldet, use_dic, use_zbC, use_zbN, use_zbP

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do

   end type type_hzg_mpb_cnp
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the MPB model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the mpb namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hzg_mpb_cnp),intent(inout),target  :: self
   integer,             intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Markus Kreus, Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
!!------- Initial values of model mpb -------
   !Note: At Hochard et al (2010) units for mpbCHL was given as mgChla/l.
   !      MK assumed that to be a typo and defined its unit to mgChla m-3, which
   !      correspondents much better to other literature refering to the "Geider" model.
   real(rk)  :: mpbCHL_init   ! MicroPhytoBenthos chlorophyll               (mgChla m-3)
   real(rk)  :: mpbC_init     ! MicroPhytoBenthos carbon                    (mmolC m-3)
   real(rk)  :: mpbN_init     ! MicroPhytoBenthos nitrogen                  (mmolN m-3)
   real(rk)  :: mpbP_init     ! MicroPhytoBenthos phosphorus                (mmolP m-3)
   real(rk)  :: eps_init      ! Extracellular Polymeric Substances          (mmolC m-3)
!!------- Parameters for model mpb originating from omexdia_p -------
!    real(rk)  :: rLdet         ! decay rate labile detritus (fast decay)     (0.075  d-1)
!    real(rk)  :: rSdet         ! decay rate semilabile detritus (slow decay) (0.003  d-1)
!    real(rk)  :: rNCldet       ! N/C ratio labile detritus (fast decay)      (0.230  molN molC-1)
!!------- Parameters for model mpb -------  (reference Values refer to Kreus et al 2015)
   real(rk)  :: Tref          ! Reference Temperature for rate parameters   (15.0  degC)
   real(rk)  :: Q10           ! Temperature factor                          (1.5   -)
   real(rk)  :: alpha         ! Initial slope of the PI-curve               (1.2   m2 molC (gChla W d))
   real(rk)  :: mu_max        ! Maximum growth rate  (at Tref)              (2.00  d-1)
   real(rk)  :: resp0         ! Maintenance respiration/lysis rate          (0.01  d-1)
   real(rk)  :: zeta_NO3      ! Biosynthetic costs for NO3 uptake           (2.30  molC molN-1)
   real(rk)  :: zeta_NH4      ! Biosynthetic costs for NH4 uptake           (1.80  molC molN-1)
   real(rk)  :: kPO4          ! Half-saturation conc. for PO4 uptake        (0.30  mmolP m-3)
   real(rk)  :: kNO3          ! Half-saturation conc. for NO3 uptake        (3.00  mmolN m-3)
   real(rk)  :: kNH4          ! Half-saturation conc. for NH4 uptake        (3.00  mmolN m-3)
   real(rk)  :: kInNH4        ! Half-saturation conc. for NO3 uptake inhibition by NH4  (10.00  mmolN m-3)
   real(rk)  :: theta_max     ! Maximum Chla/N ratio                        (3.80  gChla molN-1 => 0.27 gChla gN-1)
   real(rk)  :: QNCmin        ! Minimum N/C ratio                           (0.05  molN molC-1)
   real(rk)  :: QNCmax        ! Maximum N/C ratio                           (0.20  molN molC-1)
   real(rk)  :: QNCupt        ! N uptake per carbon ratio                   (0.20  molN molC-1)
   real(rk)  :: QPNmin        ! Minimum P/N ratio                           (0.05  molP molN-1)
   real(rk)  :: QPNmax        ! Maximum P/N ratio                           (0.20  molP molN-1)
   real(rk)  :: QPNupt        ! P uptake per nitrogen ratio                 (0.20  molP molN-1)
   real(rk)  :: sigma_NC      ! Shape factor for N uptake regulation        (      -)
   real(rk)  :: sigma_PN      ! Shape factor for P uptake regulation        (      -)
   real(rk)  :: KO2resp       ! Half-saturation conc. O2 lim. for resp      (1.00  mmolO2 m-3)
   real(rk)  :: gamma         ! Mol O2 produced per mol C fixed by photosynthesis (1.0  molO2 molC-1)
   real(rk)  :: fracEPS       ! Fraction of primary prod. exudated as EPS   (0.20  -)
   real(rk)  :: degrEPS       ! Degradation rate for EPS                    (0.02  d-1)
   real(rk)  :: graz          ! Grazing rate                                (0.10  d-1)
   real(rk)  :: aeff          ! Assimilation efficiency (-> exportN)        (0.10  -)
   real(rk)  :: exud          ! Nitrogen exudation coeff. for zoobenthos    (0.15  -)
   real(rk)  :: rzoo          ! Respiration rate for zoobenthos             (0.10  d-1)
!
!!------- Optional external dependencies -------
   character(len=attribute_length) :: no3_variable  = ''  ! dissolved nitrate
   character(len=attribute_length) :: nh4_variable  = ''  ! dissolved ammonium
   character(len=attribute_length) :: po4_variable  = ''  ! dissolved phosphate
   character(len=attribute_length) :: oxy_variable  = ''  ! dissolved oxygen
   character(len=attribute_length) :: ldet_variable = ''  ! labile detritus carbon (fast decay)
   character(len=attribute_length) :: dic_variable  = ''  ! dissolved inorganic carbon
   character(len=attribute_length) :: zbC_variable  = ''  ! zoobenthos carbon
   character(len=attribute_length) :: zbN_variable  = ''  ! zoobenthos nitrogen
   character(len=attribute_length) :: zbP_variable  = ''  ! zoobenthos phosphorus

   namelist /hzg_mpb_cnp/ &
          Tref, Q10, alpha, mu_max, resp0, zeta_NO3, zeta_NH4, kPO4, kNO3, kNH4, kInNH4, &
          theta_max, QNCmin, QNCmax, QNCupt, QPNmin, QPNmax, QPNupt, sigma_NC, sigma_PN, &
          kO2resp, gamma, fracEPS, degrEPS, graz, aeff, exud, rzoo, &
          mpbCHL_init, mpbC_init, mpbN_init, mpbP_init, eps_init

   namelist /hzg_mpb_cnp_dependencies/  &
          no3_variable, nh4_variable, po4_variable, oxy_variable, ldet_variable, &
          dic_variable, zbC_variable, zbN_variable, zbP_variable

!EOP
!-----------------------------------------------------------------------
!BOC

   ! Read the namelist
   if (configunit>0) then
     read(configunit, nml=hzg_mpb_cnp, err=98, end=100)
     read(configunit, nml=hzg_mpb_cnp_dependencies, err=99, end=101)
   endif

   ! Store parameter values in our own derived type
   self%Tref       = Tref
   self%Q10        = Q10
   self%alpha      = alpha
   self%mu_max     = mu_max
   self%resp0      = resp0
   self%zeta_NO3   = zeta_NO3
   self%zeta_NH4   = zeta_NH4
   self%kPO4       = kPO4
   self%kNO3       = kNO3
   self%kNH4       = kNH4
   self%kInNH4     = kInNH4
   self%theta_max  = theta_max
   self%QNCmin     = QNCmin
   self%QNCmax     = QNCmax
   self%QNCupt     = QNCupt
   self%QPNmin     = QPNmin
   self%QPNmax     = QPNmax
   self%QPNupt     = QPNupt
   self%sigma_NC   = sigma_NC
   self%sigma_PN   = sigma_PN
   self%KO2resp    = KO2resp
   self%gamma      = gamma
   self%fracEPS    = fracEPS
   self%degrEPS    = degrEPS
   self%graz       = graz
   self%aeff       = aeff
   self%exud       = exud
   self%rzoo       = rzoo

   ! Register state variables
   call self%register_state_variable(self%id_mpbCHL, 'mpbCHL', 'mgChla m-3', &
         'MicroPhytoBenthos chlorophyll mpbCHL',                            &
         mpbCHL_init, minimum=0.0_rk, no_river_dilution=.true.)
   call self%set_variable_property(self%id_mpbCHL, 'particulate', .true.)
   !Note: At Hochard et al (2010) units for mpbCHL was given as mgChla/l.
   !      I (MK) assume that to be a typo and define its unit to mgChla m-3, which
   !      correspondents much better to other literature refering to the "Geider" model.

   call self%register_state_variable(self%id_mpbC,   'mpbC', 'mmolC m-3',   &
         'MicroPhytoBenthos carbon mpbC',                                   &
         mpbC_init, minimum=0.0_rk, no_river_dilution=.true.)
   call self%set_variable_property(self%id_mpbC, 'particulate', .true.)

   call self%register_state_variable(self%id_mpbN,   'mpbN', 'mmolN m-3',   &
         'MicroPhytoBenthos nitrogen mpbN',                                 &
         mpbN_init, minimum=0.0_rk, no_river_dilution=.true.)
   call self%set_variable_property(self%id_mpbN, 'particulate', .true.)

   call self%register_state_variable(self%id_mpbP,   'mpbP', 'mmolP m-3',   &
         'MicroPhytoBenthos phosphorus mpbP',                               &
         mpbP_init, minimum=0.0_rk, no_river_dilution=.true.)
   call self%set_variable_property(self%id_mpbP, 'particulate', .true.)

   call self%register_state_variable(self%id_eps,    'eps', 'mmolC m-3',    &
         'Extracellular Polymeric Substances eps',                          &
         eps_init, minimum=0.0_rk, no_river_dilution=.true.)
   call self%set_variable_property(self%id_eps, 'particulate', .false.)

   ! Register link to external dependencies, if variable names are provided in namelist.
   self%use_no3 = no3_variable/=''
   if (self%use_no3) then
      call self%register_state_dependency(self%id_no3, no3_variable, 'mmolN m-3',   &
            'dissolved nitrate', required=.true.)
      call self%request_coupling(self%id_no3, no3_variable)
   else
      call self%register_state_variable(self%id_no3, 'no3', 'mmolN m-3',    &
            'dissolved nitrate', 20._rk, minimum=0.0_rk,                    &
            standard_variable=standard_variables%mole_concentration_of_nitrate)
      call self%set_variable_property(self%id_no3, 'particulate', .false.)
   endif

   self%use_nh4 = nh4_variable/=''
   if (self%use_nh4) then
      call self%register_state_dependency(self%id_nh4, nh4_variable, 'mmolN m-3',   &
            'dissolved ammonium', required=.true.)
      call self%request_coupling(self%id_nh4, nh4_variable)
   else
      call self%register_state_variable(self%id_nh4, 'nh4', 'mmolN m-3',    &
            'dissolved ammonium', 40._rk, minimum=0.0_rk,                   &
            standard_variable=standard_variables%mole_concentration_of_ammonium)
      call self%set_variable_property(self%id_nh4, 'particulate', .false.)
   endif

   self%use_po4 = po4_variable/=''
   if (self%use_po4) then
      call self%register_state_dependency(self%id_po4, po4_variable, 'mmolP m-3',   &
            'dissolved phosphate', required=.true.)
      call self%request_coupling(self%id_po4, po4_variable)
   else
      call self%register_state_variable(self%id_po4, 'po4', 'mmolP m-3',    &
            'dissolved phosphate', 4._rk, minimum=0.0_rk,                   &
            standard_variable=standard_variables%mole_concentration_of_phosphate)
      call self%set_variable_property(self%id_po4, 'particulate', .false.)
   endif

   self%use_oxy  = oxy_variable/=''
   if (self%use_oxy) then
      call self%register_state_dependency(self%id_oxy, oxy_variable, 'mmolO2 m-3',  &
            'dissolved oxygen', required=.true.)
      call self%request_coupling(self%id_oxy, oxy_variable)
   else
      call self%register_state_variable(self%id_oxy, 'oxy', 'mmolO2 m-3',   &
            'dissolved oxygen', 100.0_rk, minimum=0.0_rk)
      call self%set_variable_property(self%id_oxy, 'particulate', .false.)
   endif

   self%use_ldet = ldet_variable/=''
   if (self%use_ldet) then
      call self%register_state_dependency(self%id_ldet, ldet_variable, 'mmolC m-3',  &
            'detritus labile carbon', required=.false.)
      call self%request_coupling(self%id_ldet, ldet_variable)
   else
      call self%register_state_variable(self%id_ldet, 'ldet', 'mmolC m-3',  &
            'detritus labile carbon', 4.e3_rk, minimum=0.0_rk)
      call self%set_variable_property(self%id_ldet, 'particulate', .true.)
   endif

   self%use_dic  = dic_variable/=''
   if (self%use_dic) then
      call self%register_state_dependency(self%id_dic, dic_variable, 'mmolC m-3',   &
            'dissolved inorganic carbon', required=.false.)
      call self%request_coupling(self%id_dic, dic_variable)
   endif

   self%use_zbC  = zbC_variable/=''
   if (self%use_zbC) then
      call self%register_state_dependency(self%id_zbC, zbC_variable, 'mmolC m-3',   &
            'ZooBenthos carbon', required=.false.)
      call self%request_coupling(self%id_zbC, zbC_variable)
   endif

   self%use_zbN  = zbN_variable/=''
   if (self%use_zbN) then
      call self%register_state_dependency(self%id_zbN, zbN_variable, 'mmolN m-3',   &
            'ZooBenthos nitrogen', required=.false.)
      call self%request_coupling(self%id_zbN, zbN_variable)
   endif

   self%use_zbP  = zbP_variable/=''
   if (self%use_zbP) then
      call self%register_state_dependency(self%id_zbP, zbP_variable, 'mmolP m-3',   &
            'ZooBenthos phosphorus', required=.false.)
      call self%request_coupling(self%id_zbP, zbP_variable)
   endif

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_PrimProd, 'MPB_PP', 'mmolC m-3 d-1',    &
         'MPB primary production rate PrimProd',     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_par,      'par', 'W m-2',               &
         'MPB photosynthetically active radiation',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_qNC,      'qNC', 'molN molC-1',         &
         'MPB nitrogen quota qNC',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Q_chl,    'Q_chl', 'gChla gC-1',        &
         'MPB CHL-to-C ratio Q_chl',                 output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_ThetaN,   'ThetaN', 'gChla molN-1',     &
         'MPB CHL-to-N ratio ThetaN',                output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_expCProd, 'expCProd', 'mmolC m-3 d-1',  &
         'MPB carbon export (zoobenthos grazing)',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_expNProd, 'expNProd', 'mmolN m-3 d-1',  &
         'MPB nitrogen export (zoobenthos grazing)', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_NPP, 'MPB_NPP',  'mmolC m-3 d-1',       &
         'MPB net primary production rate NPP',      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_SGR, 'MPB_SGR',  'd-1',                 &
         'MPB specific growth rate SGR',             output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_TGR, 'MPB_TGR',  'd-1',                 &
         'MPB total growth rate TGR',                output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_SPR, 'MPB_SPR',  'd-1',                 &
         'MPB specific photosynthesis rate SPR',     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_SMR, 'MPB_SMR',  'd-1',                 &
         'MPB specific mortality rate (grazing) SMR', output=output_instantaneous)
! temporary for debugging:
   call self%register_diagnostic_variable(self%id_MPB_DIN, 'MPB_DIN', 'mmolN m-3 d-1',    &
         'MPB DIN',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_no3, 'MPB_no3', 'mmolN m-3 d-1',    &
         'MPB no3',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_nh4, 'MPB_nh4', 'mmolN m-3 d-1',    &
         'MPB nh4',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_po4, 'MPB_po4', 'mmolP m-3 d-1',    &
         'MPB po4',  output=output_instantaneous)

   ! Register dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_parz, standard_variables%downwelling_photosynthetic_radiative_flux)

   return

98 call self%fatal_error('hzg_mpb_cnp_initialize','Error reading namelist hzg_mpb_cnp.')
99 call self%fatal_error('hzg_mpb_cnp_initialize','Error reading namelist hzg_mpb_cnp_dependencies')

100 call self%fatal_error('hzg_mpb_cnp_initialize','Namelist hzg_mpb_cnp was not found.')
101 call self%fatal_error('hzg_mpb_cnp_initialize','Namelist hzg_mpb_cnp_dependencies was not found.')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of MPB model
!
! !INTERFACE:
   subroutine do(self, _ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_hzg_mpb_cnp),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Markus Kreus, Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
   real(rk), parameter :: zero = 0.0_rk, one = 1.0_rk, TINY = 1.0e-3
   real(rk), parameter :: T0       = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk), parameter :: Q10b     = 1.5_rk    ! q10 temperature coefficient (source?)
   real(rk) :: mpbC, mpbN, mpbP, mpbCHL, eps, no3, nh4, po4, oxy, ldet
   real(rk) :: temp_celsius, temp_kelvin, f_T, E_a, tfac, parz, porosity
   real(rk) :: qNC, qPN, Q_chl, theta, mu_C, mu_N, mu_P, mu_chl, limN, limP, limO2
   real(rk) :: Pmax, f_IR, f_RCN, f_RPC, f_RPN, f_RNP, f_RNC, frac_NH4
   real(rk) :: prod, prodChl, prodO2, respphyto, uptP, uptN, uptNO3, uptNH4, prodEPS, degrEPS
   real(rk) :: grazingC, grazingN, grazingP, grazingChl, faecesC, faecesN, faecesP
   real(rk) :: respZoo, exudZoN, exudZoP, exportC, exportN, exportP
   real(rk) :: NPP, SGR, TGR, GRZ

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! initialize local external dependencies
   no3  = -9999.
   nh4  = -9999.
   po4  = -9999.
   oxy  = -9999.
   ldet = -9999.

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp, temp_celsius) ! sediment-water temperature in deg C
   _GET_(self%id_parz,   parz)    ! sediment light in W m-2
   _GET_(self%id_mpbCHL, mpbCHL)  ! MicroPhytoBenthos chlorophyll in mgChla m-3
   _GET_(self%id_mpbC,   mpbC)    ! MicroPhytoBenthos carbon in mmolC m-3
   _GET_(self%id_mpbN,   mpbN)    ! MicroPhytoBenthos nitrogen in mmolN m-3
   _GET_(self%id_mpbP,   mpbP)    ! MicroPhytoBenthos phosphorus in mmolP m-3
   _GET_(self%id_eps,    eps)     ! Extracellular Polymeric Substances in mmolC m-3
   ! Retrieve current external (local) dependencies if available
   _GET_(self%id_no3,    no3)     ! dissolved nitrate in mmolN m-3
   _GET_(self%id_nh4,    nh4)     ! dissolved ammonium in mmolN m-3
   _GET_(self%id_po4,    po4)     ! dissolved phosphate in mmolP m-3
   _GET_(self%id_oxy,    oxy)     ! dissolved oxygen in mmolO2 m-3
   _GET_(self%id_ldet,   ldet)    ! fast decaying detritus C in mmolC m-3

   ! global temperature dependency
   temp_kelvin = 273.15_rk + temp_celsius
   E_a = 0.1_rk*log(Q10b)*T0*(T0+10.0_rk)
   f_T = one*exp(-E_a*(one/temp_kelvin - one/T0)) ! temperature factor (-)


   !----- intracellular N:C:Chl stoichiometry
   qNC    = mpbN   / mpbC      ! molN  molC-1
   qPN    = mpbP   / mpbN      ! molP  molN-1
   Q_chl  = mpbCHL / mpbC      ! gChla molC-1
   theta  = mpbCHL / mpbN      ! gChla molN-1

   !----- microphytes temperature dependency
   tfac = exp( log(self%Q10)*((temp_celsius-10.)/10.) )                ! (-)

   !----- oxygen limitation factor
   limO2 = max( zero, oxy / (oxy+self%KO2resp) )                       ! (-)

   !----- photosynthesis and carbon assimilation
   f_RCN = max( zero, (qNC-self%QNCmin)/(self%QNCmax-self%QNCmin) )    ! (-)
   Pmax = self%mu_max *tfac *min(one,f_RCN)                            ! (d-1)
   print*,'#435 ', Pmax, tfac, f_RCN

   ! - for comparism: Pmax as formulated in i.e. Baumert, Pahlow etc.
   !Pmax = self%mu_max *tfac *(one- self%QNCmin/qNC)                   ! (d-1)
   ! light limitation
   f_IR = one -exp(-self%alpha *theta *parz /Pmax)                     ! (-)
   ! carbon assimilation rate
   mu_C = Pmax *f_IR *limO2                                            ! (d-1)
   prod = mu_C * mpbC                                                  ! (mmolC m-3 d-1)

   print*,'#443 ', prod, Pmax, f_IR, limO2

   !----- phosphorus assimilation and regulation
   limP  = max( zero, po4/(self%KPO4 + po4) )                          ! (-)
   ! down-regulation of P-uptake if cell is P-replete
   f_RPC = one                                                         ! (-)
   ! down-regulation of P-uptake due to higher P:N-quota caused by N-limitation
   f_RPN = one/( one + exp(-self%sigma_PN *(self%QPNmax-qPN)) )        ! (-)
   ! P-assimilation is not directly related to amount of proteins (as "expressed" by qNC)
   mu_P = self%QPNupt *self%mu_max *tfac *f_RPC *f_RPN *limP *limO2    ! (mmolP mmolN-1 d-1)
   ! net P uptake rate
   uptP = (mu_P * mpbN) - self%resp0*tfac*mpbP

   !----- nitrogen assimilation and regulation
   limN = no3/(no3+self%KNO3) *(one - nh4/(nh4+self%KinNH4) ) + nh4/(nh4+self%KNH4) ! (-)
   frac_NH4 = nh4/(nh4 + self%KNH4) /limN                              ! (-)
   f_RNP = max( zero, one - (self%QPNmin /qPN) )                       ! (-) (Droop)
   f_RNC = one /( one + exp(-self%sigma_NC *(self%QNCmax-qNC)) )       ! (-)
   ! carbon specific nitrogen uptake rate
   mu_N = self%QNCupt *self%mu_max *tfac *f_RNP *f_RNC *limN *limO2    ! (mmolN mmolC-1 d-1) (vgl. Geider)
   if (mu_N .lt. TINY .or. qNC .gt. self%QNCmax) mu_N = zero           ! (mmolN mmolC-1 d-1)
   ! net N uptake rates
   uptNH4  = ( (      frac_NH4)*mu_N - self%resp0*tfac ) *mpbC         ! (mmolN d-1)
   uptNO3  = ( (one - frac_NH4)*mu_N - self%resp0*tfac ) *mpbC         ! (mmolN d-1)
   uptN    = uptNO3 + uptNH4                                           ! (mmolN d-1)

   !----- chlorophyll synthesis
   if (mu_C .le. TINY) then
      mu_chl = mu_N *self%theta_max                                    ! (mgChla mmolC-1)
   else
      mu_chl = mu_N *self%theta_max *mu_C/(self%alpha *theta *parz)    ! (mgChla mmolC-1)
   endif
   ! net Chl synthesis
   prodchl = mu_chl *mpbC - self%resp0 *tfac *mpbCHL                   ! (mgChla m-3 d-1)

   !----- oxygen production and respiration loss:
   prodO2  = self%gamma * prod                                         ! (mmolO2 m-3 d-1)
   respphyto = self%resp0 *tfac *limO2 *mpbC +  &
               uptNO3 *self%zeta_NO3 + uptNH4 *self%zeta_NH4           ! (mmolC m-3 d-1)

   !----- Carbohydrate exudation:
   prodEPS = self%fracEPS *prod                                        ! (mmolC m-3 d-1)
   !CprodEPS = sqrt(self%rLdet*self%rSdet) * eps  ! Source of this formulation? Kai Wirtz ??
   degrEPS = max(zero, f_T *self%degrEPS *(eps-tiny) )                 ! (mmolC m-3 d-1)

   !----- Zoobenthos grazing and associated processes:
   grazingC   = self%graz *limO2 * mpbC                                ! (mmolC m-3 d-1)
   grazingN   = self%graz *limO2 * mpbN                                ! (mmolN m-3 d-1)
   grazingP   = self%graz *limO2 * mpbP                                ! (mmolP m-3 d-1)
   grazingChl = self%graz *limO2 * mpbCHL                              ! (mgChla m-3 d-1)
   faecesC    = (one -self%aeff) * grazingC                            ! (mmolC m-3 d-1)
   faecesN    = (one -self%aeff) * grazingN                            ! (mmolN m-3 d-1)
   faecesP    = (one -self%aeff) * grazingP                            ! (mmolN m-3 d-1)
   respZoo    = self%rzoo * grazingC                                   ! (mmolC m-3 d-1)
   exudZoN    = self%exud * grazingN                                   ! (mmolN m-3 d-1)
   exudZoP    = self%exud * grazingP                                   ! (mmolP m-3 d-1)
   exportC    = grazingC - faecesC - respZoo ! exported carbon   (open closure term)   ! (mmolC m-3 d-1)
   exportN    = grazingN - faecesN - exudZoN ! exported nitrogen (open closure term)   ! (mmolN m-3 d-1)
   exportP    = grazingP - faecesP - exudZoP ! exported nitrogen (open closure term)   ! (mmolN m-3 d-1)

   !----- Diagnostic parameters
   NPP =  prod - respphyto                                             ! (mmolC m-3 d-1)
   SGR = (prod - respphyto - prodeps ) /mpbC                           ! (d-1)
   TGR = (prod - grazingC - respphyto - prodeps) /mpbC                 ! (d-1)
   GRZ = (grazingC) /mpbC                                              ! (d-1)

#define _CONV_UNIT_ *one_pr_day
   ! reaction rates
   _SET_ODE_(self%id_mpbCHL, (prodchl                       - grazingChl)      _CONV_UNIT_)
   _SET_ODE_(self%id_mpbC,   (prod    - respphyto - prodeps - grazingC)        _CONV_UNIT_)
   _SET_ODE_(self%id_mpbN,   (uptN                          - grazingN)        _CONV_UNIT_)
   _SET_ODE_(self%id_eps,    (prodEPS - degrEPS)                               _CONV_UNIT_)
   ! external dependencies
   ! If externally maintained variables are present, change the pools accordingly
   _SET_ODE_(self%id_no3,    (        - uptNO3)                                _CONV_UNIT_)
   _SET_ODE_(self%id_nh4,    (exudZoN - uptNH4)                                _CONV_UNIT_)
   _SET_ODE_(self%id_po4,    (exudZoP - uptP)                                  _CONV_UNIT_)
   _SET_ODE_(self%id_oxy,    (prodO2  - respphyto - respZoo )                  _CONV_UNIT_)
   _SET_ODE_(self%id_ldet,   (faecesC)                                         _CONV_UNIT_)
   if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic , (respphyto + respZoo) _CONV_UNIT_)
   if (_AVAILABLE_(self%id_zbC)) _SET_ODE_(self%id_zbC , (exportC)             _CONV_UNIT_)
   if (_AVAILABLE_(self%id_zbN)) _SET_ODE_(self%id_zbN , (exportN)             _CONV_UNIT_)
   if (_AVAILABLE_(self%id_zbP)) _SET_ODE_(self%id_zbP , (exportP)             _CONV_UNIT_)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_SPR,      mu_C)           !instantaneous MPB specific photosynthesis rate (d-1)
   _SET_DIAGNOSTIC_(self%id_PrimProd, prod)           !instantaneous MPB primary production rate (molC m-3 d-1)
   _SET_DIAGNOSTIC_(self%id_par,      parz)           !instantaneous MPB photosynthetically active radiation (W m-2)
   _SET_DIAGNOSTIC_(self%id_qNC,      qNC)            !instantaneous MPB N:C quota (molN molC-1)
   _SET_DIAGNOSTIC_(self%id_Q_chl,    Q_chl/12.)      !instantaneous MPB CHL:C ratio (gChla gC-1)
   _SET_DIAGNOSTIC_(self%id_ThetaN,   theta)          !instantaneous MPB CHL:N ratio (gChla molN-1)
   _SET_DIAGNOSTIC_(self%id_expCProd, exportC)        !instantaneous MPB export production carbon (mmolC m-3 d-1)
   _SET_DIAGNOSTIC_(self%id_expNProd, exportN)        !instantaneous MPB export production nitrogen (mmolN m-3 d-1)
   _SET_DIAGNOSTIC_(self%id_expPProd, exportP)        !instantaneous MPB export production nitrogen (mmolN m-3 d-1)
   _SET_DIAGNOSTIC_(self%id_NPP,      NPP)            !instantaneous MPB net primary production rate (molC m-3 d-1)
   _SET_DIAGNOSTIC_(self%id_SGR,      SGR)            !instantaneous MPB specific growth rate (d-1)
   _SET_DIAGNOSTIC_(self%id_TGR,      TGR)            !instantaneous MPB total growth rate (d-1)
   _SET_DIAGNOSTIC_(self%id_SMR,      GRZ)            !instantaneous MPB specific mortality rate (d-1)
! temporary for debugging:
   _SET_DIAGNOSTIC_(self%id_mpb_din,  nh4+no3)        !instantaneous external DIN (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_no3,  no3)            !instantaneous external no3 (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_nh4,  nh4)            !instantaneous external nh4 (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_po4,  po4)            !instantaneous external nh4 (mmolN m-3)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
  end module hzg_mpb_cnp
!-----------------------------------------------------------------------

