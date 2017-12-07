#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_mpb --- Fortran 2003 version of microphytobenthos (mpb) model
!
! !INTERFACE:
   module hzg_mpb
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
   public type_hzg_mpb
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
   type,extends(type_base_model) :: type_hzg_mpb
!     Variable identifiers
      type (type_state_variable_id)        :: id_no3, id_nh4, id_po4, id_oxy, id_ldetC
      type (type_state_variable_id)        :: id_dic, id_zbC, id_zbN, id_zbP
      type (type_state_variable_id)        :: id_mpbCHL, id_mpbC, id_mpbN, id_eps
      type (type_dependency_id)            :: id_temp, id_parz, id_porosity
      type (type_diagnostic_variable_id)   :: id_PrimProd, id_par, id_Q_N, id_Q_chl, id_Theta_N
      type (type_diagnostic_variable_id)   :: id_expCProd, id_expNProd, id_mpbP
      type (type_diagnostic_variable_id)   :: id_NPP, id_SGR, id_TGR, id_SPR, id_SMR
! temporary for debugging
      type (type_diagnostic_variable_id)   :: id_MPB_din, id_MPB_no3, id_MPB_nh4, id_MPB_po4
      type (type_diagnostic_variable_id)   :: id_MPB_temp, id_MPB_totN, id_MPB_totP
      type (type_diagnostic_variable_id)   :: id_fIR, id_flimN, id_flimP, id_flimO2, id_uptN, id_uptP
      type (type_diagnostic_variable_id)   :: id_fRChl, id_fRCN, id_fRNC, id_fRNP, id_fRPN, id_ftfac

!     Model parameters
      real(rk) :: Tref, mumax, alpha, gamma, Qmin, Qmax, thetamax, uptmax, KNH4, KNO3
      real(rk) :: KinNH4, kEPS, rEPS, resp, Kresp, graz, kout, kexu, rzoo
      real(rk) :: k0, Achla, btemp
      logical  :: use_no3, use_nh4, use_po4, use_oxy, use_ldetC, use_dic, use_zbC, use_zbN

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do

   end type type_hzg_mpb
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
   class (type_hzg_mpb),intent(inout),target  :: self
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
   real(rk)  :: eps_init      ! Extracellular Polymeric Substances          (mmolC m-3)
!!------- Parameters for model mpb originating from omexdia_p -------
!    real(rk)  :: rLdet         ! decay rate labile detritus (fast decay)     (0.075  d-1)
!    real(rk)  :: rSdet         ! decay rate semilabile detritus (slow decay) (0.003  d-1)
!    real(rk)  :: rNCldet       ! N/C ratio labile detritus (fast decay)      (0.230  molN molC-1)
!!------- Parameters for model mpb -------  (reference Values refer to Hochard (2010))
   real(rk)  :: Tref          ! Reference Temperature for rate parameters   (0.0   degC)
   real(rk)  :: mumax         ! Maximum growth rate  (at Tref)              (0.70  d-1)
   real(rk)  :: btemp         ! Temperature factor                          (0.063  - )
   real(rk)  :: alpha         ! initial slope of the PI-curve               (1.296 m2 molC (gChla W d) := 1.5e-5  m2 molC (gChla J)-1)
   real(rk)  :: gamma         ! Mol O2 produced per mol C fixed by photosynthesis (1.0  molO2 molC-1)
   real(rk)  :: Qmin          ! Minimum N/C ratio                           (0.05  molN molC-1)
   real(rk)  :: Qmax          ! Maximum N/C ratio                           (0.20  molN molC-1)
   real(rk)  :: thetamax      ! Maximum Chla/N ratio                        (3.80  gChla molN-1 => 0.27 gChla gN-1)
   real(rk)  :: uptmax        ! Maximum N uptake rate per carbon unit (at Tref) (0.20  molN molC-1 d-1)
   real(rk)  :: KNH4          ! Half-saturation conc. for NH4 uptake        (3.00  mmolN m-3)
   real(rk)  :: KNO3          ! Half-saturation conc. for NO3 uptake        (3.00  mmolN m-3)
   real(rk)  :: KinNH4        ! Half-saturation conc. for NO3 uptake inhibition by NH4  (10.00  mmolN m-3)
   real(rk)  :: keps          ! fraction of primary production exudated as EPS  (0.20  -)
   real(rk)  :: rEPS          ! decay rate for EPS                          (0.02  d-1)
   real(rk)  :: resp          ! Respiration rate                            (0.10  d-1)
   real(rk)  :: Kresp         ! Half-saturation conc. O2 lim. for resp      (1.00  mmolO2 m-3)
   real(rk)  :: graz          ! Grazing rate                                (0.10  d-1)
   real(rk)  :: kout          ! Faeces production coeff. (-> exportN)       (0.10  -)
   real(rk)  :: kexu          ! Nitrogen exudation coeff.                   (0.15  -)
   real(rk)  :: rzoo          ! Respiration rate for zoobenthos             (0.10  d-1)
   real(rk)  :: k0            !*Extinction coefficient of sediment          (20.0  cm-1)         *obsolete/currently not used!!!
   real(rk)  :: Achla         !*Absorption factor of chlorophyll            (0.02  m2 mgChla-1)  *obsolete/currently not used!!!
!!------- Optional external dependencies -------
   character(len=attribute_length) :: no3_variable   = ''  ! dissolved nitrate
   character(len=attribute_length) :: nh4_variable   = ''  ! dissolved ammonium
   character(len=attribute_length) :: po4_variable   = ''  ! dissolved phosphorus
   character(len=attribute_length) :: oxy_variable   = ''  ! dissolved oxygen
   character(len=attribute_length) :: ldetC_variable = ''  ! labile detritus carbon (fast decay)
   character(len=attribute_length) :: dic_variable   = ''  ! dissolved inorganic carbon
   character(len=attribute_length) :: zbC_variable   = ''  ! zoobenthos carbon
   character(len=attribute_length) :: zbN_variable   = ''  ! zoobenthos nitrogen

   namelist /hzg_mpb/ &
          Tref, mumax, btemp, alpha, gamma, Qmin, Qmax, thetamax, uptmax, KNH4, KNO3,   &
          KinNH4, kEPS, rEPS, resp, Kresp, graz, kout, kexu, rzoo, k0, Achla,           &
          mpbCHL_init, mpbC_init, mpbN_init, eps_init

   namelist /hzg_mpb_dependencies/  &
          no3_variable, nh4_variable, po4_variable, oxy_variable, ldetC_variable, &
          dic_variable, zbC_variable, zbN_variable

!EOP
!-----------------------------------------------------------------------
!BOC

   ! Read the namelist
   if (configunit>0) then
     read(configunit, nml=hzg_mpb, err=98, end=100)
     read(configunit, nml=hzg_mpb_dependencies, err=99, end=101)
   endif

   ! Store parameter values in our own derived type
   self%Tref        = Tref
   self%mumax       = mumax
   self%btemp       = btemp
   self%alpha       = alpha
   self%gamma       = gamma
   self%Qmin        = Qmin
   self%Qmax        = Qmax
   self%thetamax    = thetamax
   self%uptmax      = uptmax
   self%KNH4        = KNH4
   self%KNO3        = KNO3
   self%KinNH4      = KinNH4
   self%keps        = keps
   self%resp        = resp
   self%Kresp       = Kresp
   self%graz        = graz
   self%kout        = kout
   self%kexu        = kexu
   self%rzoo        = rzoo
   self%k0          = k0
   self%Achla       = Achla

   ! Register state variables
   call self%register_state_variable(self%id_mpbCHL, 'mpbCHL', 'mgChla m-3', &
         'MicroPhytoBenthos chlorophyll mpbCHL',                             &
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

   call self%register_diagnostic_variable(self%id_mpbP,'mpbP','mmolP m-3',  &
         'MicroPhytoBenthos phosphorus mpbP',  output=output_instantaneous)

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
      call self%register_state_variable(self%id_no3,   'no3', 'mmolN m-3',   &
            'dissolved nitrate', 20.0_rk, minimum=0.0_rk,                    &
            standard_variable=standard_variables%mole_concentration_of_nitrate)
      call self%set_variable_property(self%id_no3, 'particulate', .false.)
   endif

   self%use_nh4 = nh4_variable/=''
   if (self%use_nh4) then
      call self%register_state_dependency(self%id_nh4, nh4_variable, 'mmolN m-3',   &
            'dissolved ammonium', required=.true.)
      call self%request_coupling(self%id_nh4, nh4_variable)
   else
      call self%register_state_variable(self%id_nh4,   'nh4', 'mmolN m-3',   &
            'dissolved ammonium', 40.0_rk, minimum=0.0_rk,                   &
            standard_variable=standard_variables%mole_concentration_of_ammonium)
      call self%set_variable_property(self%id_nh4, 'particulate', .false.)
   endif

   self%use_po4 = po4_variable/=''
   if (self%use_po4) then
      call self%register_state_dependency(self%id_po4, po4_variable, 'mmolP m-3',   &
            'dissolved phosphate', required=.true.)
      call self%request_coupling(self%id_po4, po4_variable)
   else
      call self%register_state_variable(self%id_po4,   'po4', 'mmolP m-3',   &
            'dissolved phosphate', 10.0_rk, minimum=0.0_rk,                  &
            standard_variable=standard_variables%mole_concentration_of_phosphate)
      call self%set_variable_property(self%id_po4, 'particulate', .false.)
   endif

   self%use_oxy  = oxy_variable/=''
   if (self%use_oxy) then
      call self%register_state_dependency(self%id_oxy, oxy_variable, 'mmolO2 m-3',  &
            'dissolved oxygen', required=.true.)
      call self%request_coupling(self%id_oxy, oxy_variable)
   else
      call self%register_state_variable(self%id_oxy,   'oxy', 'mmolO2 m-3',  &
            'dissolved oxygen', 100.0_rk, minimum=0.0_rk)
      call self%set_variable_property(self%id_oxy, 'particulate', .false.)
   endif

   self%use_ldetC = ldetC_variable/=''
   if (self%use_ldetC) then
      call self%register_state_dependency(self%id_ldetC, ldetC_variable, 'mmolC m-3',  &
            'detritus labile carbon', required=.false.)
      call self%request_coupling(self%id_ldetC, ldetC_variable)
   else
      call self%register_state_variable(self%id_ldetC,   'ldetC', 'mmolC m-3',  &
            'detritus labile carbon', 4.e3_rk, minimum=0.0_rk)
      call self%set_variable_property(self%id_ldetC, 'particulate', .true.)
   endif

   self%use_dic  = dic_variable/=''
   if (self%use_dic) then
      call self%register_state_dependency(self%id_dic, 'dic', 'mmolC m-3',   &
            'dissolved inorganic carbon', required=.false.)
      call self%request_coupling(self%id_dic, dic_variable)
   endif

   self%use_zbC  = zbC_variable/=''
   if (self%use_zbC) then
      call self%register_state_dependency(self%id_zbC, 'zbC', 'mmolC m-3',   &
            'ZooBenthos carbon', required=.false.)
      call self%request_coupling(self%id_zbC, zbC_variable)
   endif

   self%use_zbN  = zbN_variable/=''
   if (self%use_zbN) then
      call self%register_state_dependency(self%id_zbN, 'zbN', 'mmolN m-3',   &
            'ZooBenthos nitrogen', required=.false.)
      call self%request_coupling(self%id_zbN, zbN_variable)
   endif

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_par,      'MPB_par',    'W m-2',           &
         'MPB photosynthetically active radiation',         output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Q_N,      'MPB_Q_N',    'molN molC-1',     &
         'MPB nitrogen quota Q_N',                          output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Q_chl,    'MPB_Q_chl',  'gChla gC-1',      &
         'MPB CHL-to-C ratio Q_chl',                        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Theta_N,  'MPB_ThetaN', 'gChla molN-1',    &
         'MPB CHL-to-N ratio ThetaN',                       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_expCProd, 'MPB_expCProd', 'mmolC m-3 d-1', &
         'MPB carbon export (zoobenthos grazing)',          output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_expNProd, 'MPB_expNProd', 'mmolN m-3 d-1', &
         'MPB nitrogen export (zoobenthos grazing)',        output=output_instantaneous)
   ! production rates
   call self%register_diagnostic_variable(self%id_PrimProd, 'MPB_PP',     'mmolC m-3 d-1',   &
         'MPB primary production rate PrimProd',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_NPP,      'MPB_NPP',    'mmolC m-3 d-1',   &
         'MPB net primary production rate NPP',             output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_SGR,      'MPB_SGR',    'd-1',             &
         'MPB specific growth rate SGR',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_TGR,      'MPB_TGR',    'd-1',             &
         'MPB total growth rate TGR',                       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_SPR,      'MPB_SPR',    'd-1',             &
         'MPB specific photosynthesis rate SPR',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_SMR,      'MPB_SMR',    'd-1',             &
         'MPB specific mortality rate (grazing) SMR',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_uptN,     'MPB_uptN',   'd-1',             &
         'MPB specific nitrogen assimilation rate uptN',    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_uptP,     'MPB_uptP',   'd-1',             &
         'MPB specific phosphorus assimilation rate uptP',  output=output_instantaneous)
   ! limitation factors
   call self%register_diagnostic_variable(self%id_ftfac,    'MPB_ftfac',  '-',               &
         'MPB ftfac',                                       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_fIR,      'MPB_fIR',    '-',               &
         'MPB fIR',                                         output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_flimN,    'MPB_flimN',  '-',               &
         'MPB flimN',                                       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_flimP,    'MPB_flimP',  '-',               &
         'MPB flimP',                                       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_flimO2,   'MPB_flimO2', '-',               &
         'MPB flimO2',                                      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_fRChl,    'MPB_fRCHl',  '-',               &
         'MPB fRChl',                                       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_fRCN,     'MPB_fRCN',   '-',               &
         'MPB fRCN',                                        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_fRNC,     'MPB_fRNC',   '-',               &
         'MPB fRNC',                                        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_fRNP,     'MPB_fRNP',   '-',               &
         'MPB fRNP',                                        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_fRPN,     'MPB_fRPN',   '-',               &
         'MPB fRPN',                                        output=output_instantaneous)
  ! temporary for debugging:
   call self%register_diagnostic_variable(self%id_MPB_temp, 'MPB_temp',   'degC',            &
         'MPB temperature',                                 output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_DIN,  'MPB_DIN',    'mmolN m-3',       &
         'MPB DIN',                                         output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_no3,  'MPB_no3',    'mmolN m-3',       &
         'MPB no3',                                         output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_nh4,  'MPB_nh4',    'mmolN m-3',       &
         'MPB nh4',                                         output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_po4,  'MPB_po4',    'mmolP m-3',       &
         'MPB po4',                                         output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_totN, 'MPB_totN',   'mmolN m-3',       &
         'MPB totN',                                        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_totP, 'MPB_totP',   'mmolP m-3',       &
         'MPB totP',                                        output=output_instantaneous)

   ! Register dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_parz, standard_variables%downwelling_photosynthetic_radiative_flux)

   return

98 call self%fatal_error('hzg_mpb_initialize','Error reading namelist hzg_mpb.')
99 call self%fatal_error('hzg_mpb_initialize','Error reading namelist hzg_mpb_dependencies')

100 call self%fatal_error('hzg_mpb_initialize','Namelist hzg_mpb was not found.')
101 call self%fatal_error('hzg_mpb_initialize','Namelist hzg_mpb_dependencies was not found.')

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
   class (type_hzg_mpb),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Markus Kreus, Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
   real(rk), parameter :: zero = 0.0_rk, one = 1.0_rk, TINY = 1.0e-3, TINY_N = 0.1
   real(rk), parameter :: gammaO2  = 1.0_rk    ! molO2 used per molC respired
   real(rk), parameter :: T0       = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk), parameter :: Q10b     = 1.5_rk    ! q10 temperature coefficient (source?)
   real(rk), parameter :: NCrLdet  = 0.20_rk   ! N/C ratio labile detritus (fast decay)
   real(rk), parameter :: qNP      = 16._rk    ! N/P quota microphytobenthos
   real(rk) :: mpbC, mpbN, mpbP, mpbCHL, eps, no3, nh4, po4, oxy, ldetC, DIN
   real(rk) :: temp_celsius, temp_kelvin, f_T, E_a, parz, porosity, degrEPS
   real(rk) :: prodChl, theta, Q_N, Q_chl, Pmax, PP, prod, prodEPS, fac, tfac, totN, totP
   real(rk) :: prodO2, rhoChl, uptNH4, uptNO3, uptChl, uptN, respphyto, lossphyto
   real(rk) :: grazingC, grazingN, grazingChl, respzoo, exud, faecesC, faecesN
   real(rk) :: exportC, exportN, NPP, SGR, TGR, GRZ, ruptN, ruptP, x1, x2
   real(rk) :: f_IR, f_RCN, f_RPC, f_RPN, f_RNP, f_RNC, f_RChl
   real(rk) :: lim_NO3, lim_NH4, limN, limP, limO2, mu_N, frac_NH4

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! initialize local external dependencies
   temp_celsius = 15.0_rk
   parz         = 50.0_rk
   no3   = -9999._rk
   nh4   = -9999._rk
   po4   = -9999._rk
   oxy   = -9999._rk
   ldetC = -9999._rk

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp, temp_celsius) ! sediment-water temperature in deg C
   _GET_(self%id_parz,   parz)    ! sediment light in W m-2
   _GET_(self%id_mpbCHL, mpbCHL)  ! MicroPhytoBenthos chlorophyll in mgChla m-3
   _GET_(self%id_mpbC,   mpbC)    ! MicroPhytoBenthos carbon in mmolC m-3
   _GET_(self%id_mpbN,   mpbN)    ! MicroPhytoBenthos nitrogen in mmolN m-3
   mpbP = mpbN/qNP                ! MicroPhytoBenthos phosphorus in mmolP m-3
   _GET_(self%id_eps,    eps)     ! Extracellular Polymeric Substances  in mmolC m-3
   ! Retrieve current external (local) dependencies if available
   _GET_(self%id_no3,    no3)     ! dissolved nitrate in mmolN m-3
   _GET_(self%id_nh4,    nh4)     ! dissolved ammonium in mmolN m-3
   _GET_(self%id_po4,    po4)     ! dissolved phosphate in mmolP m-3
   _GET_(self%id_oxy,    oxy)     ! dissolved oxygen in mmolO2 m-3
   _GET_(self%id_ldetC,  ldetC)   ! fast decaying detritus C in mmolC m-3
   DIN = no3+nh4
   totN = DIN + mpbN
   totP = po4 + mpbP

   ! temperature dependency
   temp_kelvin = 273.15_rk + temp_celsius
   E_a = 0.1_rk *log(Q10b) *T0 *(T0+10.0_rk)
   f_T = one *exp( -E_a *(one/temp_kelvin - one/T0) ) ! temperature factor (-)

   ! light climate (obsolete)
   !par    = self%PAR0max * parz
   !print*,'par#386: ',par, parz
   !k      = self%k0 + self%Achla * mpbCHL

   !----- intracellular N:C:Chl stoichiometry
   Q_N    = mpbN   / mpbC      ! molN  molC-1
   Q_chl  = mpbCHL / mpbC      ! gChla molC-1
   theta  = mpbCHL / mpbN      ! gChla molN-1

   !----- microphytes temperature dependency
   tfac   = exp( self%btemp * (temp_celsius-self%Tref) )           ! (-) (@HochardEtAl Tref:=0degC) [eq. 19]

   !----- oxygen limitation factor
   limO2 = max( zero, oxy / (oxy+self%Kresp) )                     ! (-)               [eq. 28]

   !----- photosynthesis rate
   f_RCN  = max( zero, one - self%Qmin/Q_N)                        ! (-)               [eq. 19]
   !Pmax   = self%mumax *tfac *f_RCN                                ! (d-1)             [eq. 19]
   Pmax   = self%mumax *tfac *min(one,f_RCN)                       ! (d-1)             [eq. 19]
   f_IR   = zero
   PP     = zero
   if (Pmax > zero) then
     ! light limitation
     f_IR = (one - exp(-self%alpha *parz *Q_chl /Pmax))            ! (-)               [eq. 18]
     ! carbon assimilation rate
     PP   = Pmax * f_IR                                            ! (d-1)             [eq. 18]
   endif
   prod   = PP * mpbC                                              ! (mmolC m-3 d-1)   [eq. 17]

#ifdef HochardEtAl
   !----- chlorophyll synthesis (eq. 21,22,23)
   rhoChl = zero
   if (parz .gt. zero) then
      rhoChl = self%thetamax * PP/(self%alpha * parz * Q_chl)      ! (gChla molN-1)    [eq. 22]
   endif
   fac    = max(zero, one - theta/self%thetamax)                   ! (-)               [eq. 23]
   f_RChl = fac/(fac + 0.05_rk)                                    ! (-)               [eq. 23]
   uptNH4 = self%uptmax *nh4/(nh4 + self%KNH4)                     ! (molN molC-1 d-1) [eq. 23]
   uptNO3 = self%uptmax *no3/(no3 + self%KNO3)                     ! (molN molC-1 d-1) [eq. 23]
   uptChl = (uptNO3 + uptNH4) *f_RChl                              ! (molN molC-1 d-1) [eq. 23]
   prodChl= rhoChl *uptChl *mpbC                                   ! (mgChla d-1)      [eq. 21]

   !----- nutrient uptake (TODO add dependencies on P, Si)
   f_RNC  = max(zero, (self%Qmax - Q_N) / (self%Qmax - self%Qmin)) ! (-)               [eq. 25,26]
   uptNO3 = uptNO3 *f_RNC *mpbC *(one - nh4/(nh4 + self%KinNH4))   ! (mmolN d-1)       [eq. 26]
   uptNH4 = uptNH4 *f_RNC *mpbC                                    ! (mmolN d-1)       [eq. 25]
   uptN   = uptNO3 + uptNH4                                        ! (mmolN d-1)       [eq. 24]

   ! normalized limN
   limN = 0.5* ( no3/(no3 + self%KNO3) *(one - nh4/(nh4 + self%KinNH4)) &
               + nh4/(nh4 + self%KNH4) )                           ! (-)

#else
   !----- nutrient uptake (TODO add dependencies on P, Si)

!    f_RNC  = max(zero, (self%Qmax - Q_N) / (self%Qmax - self%Qmin)) ! (-)               [eq. 25,26]
!    lim_NO3 = no3/(no3 + self%KNO3) *(one - nh4/(nh4 + self%KinNH4))! (-)               [eq. 25]
!    lim_NH4 = nh4/(nh4 + self%KNH4)                                 ! (-)               [eq. 26]
!    !NOTE (mk): In Hochard et al (2010) no temperature dependency is considered for
!    !           eq. 25+26 even though it is in the model they refer to (Geider,1998) !!
!    uptNO3 = self%uptmax *tfac *f_RNC *lim_NO3 *mpbC                ! (mmolN d-1)       [eq. 25]
!    uptNH4 = self%uptmax *tfac *f_RNC *lim_NH4 *mpbC                ! (mmolN d-1)       [eq. 26]
!    uptN   = uptNO3 + uptNH4                                        ! (mmolN d-1)       [eq. 24]
!
!    ! normalized limN
!    limN = 0.5* (lim_NO3 + lim_NH4)
! (-)
   !----- nitrogen assimilation and regulation @ KreusEtAl
   frac_NH4 = max( zero, nh4 /(nh4 + self%KinNH4) )                    ! (-)
   x1 = max( zero, no3 /self%KNO3 )
   x2 = max( zero, nh4 /self%KNH4 )
   limN = (x1+x2) /(one+x1+x2)
   !limN = (no3+nh4) /(self%KNO3+no3+nh4)
   !limN = nh4/(nh4 + self%KNH4) + no3/(no3 + self%KNO3) *(one - frac_NH4) ! (-)
   if (limN > zero) then
      frac_NH4 = max( zero, one-( x1*(one-frac_NH4)/(one+x1+x2)) )     ! (-)
      !frac_NH4 = max( zero, one-( no3*(one-frac_NH4)/(self%KNO3+no3+nh4)) ) ! (-)
   endif
   if (nh4 <= zero) frac_NH4 = zero

   f_RNC = max(zero, (self%Qmax - Q_N) / (self%Qmax - self%Qmin)) ! (-)               [eq. 25,26]
   !f_RNC = one /( one + exp(-100. *(self%Qmax-Q_N)) )       ! (-) !for comparism @ KreusEtAl
   !f_RNC = max( zero, two/( one + exp(-self%sigma_NC *(self%QNCmax-qNC)) ) -one ) ! (-)
   ! carbon specific nitrogen uptake rate
   mu_N = self%Qmax *self%mumax *tfac *f_RNC *limN                     ! (mmolN mmolC-1 d-1) (vgl. Geider)
   if (mu_N .lt. TINY ) mu_N = zero                                    ! (mmolN mmolC-1 d-1)
   ! N uptake rates
   uptNH4  = (      frac_NH4)*mu_N*mpbC                                ! (mmolN d-1)
   uptNO3  = (one - frac_NH4)*mu_N*mpbC                                ! (mmolN d-1)
   uptN    = uptNO3 + uptNH4                                           ! (mmolN d-1)



   !----- Carbohydrate exudation:
   prodEPS = self%keps *prod                                       ! (mmolC m-3 d-1)   [eq. 27]
   !CprodEPS = sqrt(self%rLdet*self%rSdet) * eps  ! Source of this formulation? Kai Wirtz ??
   degrEPS = max( zero, self%rEPS *tfac *(eps-tiny) )              ! (mmolC m-3 d-1)

   !----- chlorophyll synthesis (eq. 21,22,23)
   !NOTE (mk): In Geider et al (1998) chlorophyll synthesis was related to nitrogen uptake !!
   rhoChl = zero
   if (parz .gt. tiny) then
      rhoChl = self%thetamax * PP/(self%alpha * parz * Q_chl)      ! (gChla molN-1)    [eq. 22]
   endif
   fac    = max(zero, one - theta/self%thetamax)                   ! (-)               [eq. 23]
   f_RChl = fac/(fac + 0.05_rk)                                    ! (-)               [eq. 23]
   uptChl = uptN *f_RChl                                           ! (mmolN  d-1)      [eq. 23]
   prodChl = rhoChl *uptChl                                        ! (mgChla d-1)      [eq. 21]
#endif

   !----- phosphorus uptake (TODO)
   limP   = one
   f_RNP  = one
   f_RPN  = one

   !----- oxygen production and respiration loss:
#ifdef HochardEtAl
   prodO2    = self%gamma * prod                                   ! (mmolO2 m-3 d-1)  [eq. 20]
   respphyto = self%resp *limO2 *mpbC                              ! (mmolC m-3 d-1)   [eq. 28]
   lossphyto = zero
#else
   prodO2    = self%gamma * prod                                   ! (mmolO2 m-3 d-1)  [eq. 20]
   respphyto = self%resp *tfac *limO2 *mpbC                        ! (mmolC m-3 d-1)   [eq. 28]
   !NOTE (mk): In Geider et al (1998) respiration accounts for costs of biosynthesis (Zeta)
   if (self%resp>0._rk .and. prod>TINY .and. uptN>TINY) then
     respphyto = respphyto + 2.3_rk*uptNO3 + 1.8_rk*uptNH4         ! (mmolC m-3 d-1)   [-]
   endif
   !NOTE (mk): In Geider et al (1998) loss due to cell lysis/damage was accounted for as well
   lossphyto = self%resp *tfac *mpbN                               ! (mmolN m-3 d-1)   [-]
#endif

   ! Zoobenthos grazing and associated processes:
   grazingC   = self%graz *mpbC                                    ! (mmolC m-3 d-1)   [eq. 29]
   grazingN   = self%graz *mpbN                                    ! (mmolN m-3 d-1)   [eq. 30]
   grazingChl = self%graz *mpbCHL                                  ! (mgChla m-3 d-1)  [eq. 31]
   faecesC    = self%kout *grazingC                                ! (mmolC m-3 d-1)   [eq. 32]
   faecesN    = faecesC * NCrLdet                                  ! (mmolN m-3 d-1)   [eq. 32]
   exud       = self%kexu *(grazingN-faecesN)                      ! (mmolN m-3 d-1)   [eq. 33]
   respzoo    = self%rzoo *(grazingC-faecesC) *limO2               ! (mmolC m-3 d-1)   [eq. 34]
   exportC    = grazingC - faecesC - respzoo ! exported carbon   (open closure term)   ! (mmolC m-3 d-1) [eq. 35]
   exportN    = grazingN - faecesN - exud    ! exported nitrogen (open closure term)   ! (mmolN m-3 d-1) [eq. 35]
   !NOTE (mk):
   !Die Formulierung für faecesN und exportN ist nicht ganz unproblematisch:
   !a) Dadurch, dass faeces stets ein Ratio "NCrLdet" (:= C/N=5) aufweisen welches i.d.R kleiner als das des MPB ist,
   !   wird u.U. der Export ein "überhöhtes" C/N aufweisen. Diesen Effekt könnte man ggfs. durch geeignete Wahl von
   !   rzoo und exud kompensieren. Aus meiner Sicht wäre es jedoch sinnvoller, eine "export stöchiometrie" festzulegen,
   !   die dem der Zoobenthos Biomasse entspricht. Die Größe "exud" ergibt sich dann diagnostisch aus der Summenbilanz.
   !b) Die gegenwärtige Formulierung kann dazu führen, dass exportN negativ wird. Das muß abgefangen werden!!!

   !----- Diagnostic parameters
   NPP =  prod - respphyto                                         ! (mmolC m-3 d-1)
   SGR = (prod - respphyto - prodEPS ) /mpbC                       ! (d-1)
   TGR = (prod - grazingC - respphyto - prodEPS) /mpbC             ! (d-1)
   GRZ = (grazingC) /mpbC                                          ! (d-1)

   ruptN = zero
   if (uptN>tiny) ruptN = uptN/mpbN                                ! (d-1)
   ruptP = ruptN/qNP                                               ! (d-1)

   if (mpbN<TINY) lossphyto = zero

   !if (uptN<tiny) print*,'mpb#560: (uptN,uptNO3,uptNH4)',uptN,uptNO3,uptNH4
   !if (lossphyto<zero) print*,'mpb#561: (lossphyto)',lossphyto
   !if (grazingN<zero) print*,'mpb#562: (grazingN)',grazingN

#define _CONV_UNIT_ *one_pr_day
   ! reaction rates
   _SET_ODE_(self%id_mpbCHL, (prodChl - lossphyto *theta - grazingChl)         _CONV_UNIT_)
   _SET_ODE_(self%id_mpbC,   (prod - respphyto - prodEPS - grazingC)           _CONV_UNIT_)
   _SET_ODE_(self%id_mpbN,   (uptN - lossphyto           - grazingN)           _CONV_UNIT_)
   _SET_ODE_(self%id_eps,    (prodEPS - degrEPS)                               _CONV_UNIT_)
   ! external dependencies
   ! If externally maintained variables are present, change the pools accordingly
   _SET_ODE_(self%id_no3 ,   (- uptNO3)                                        _CONV_UNIT_)
   _SET_ODE_(self%id_nh4 ,   (- uptNH4 + lossphyto + exud )                    _CONV_UNIT_)
   _SET_ODE_(self%id_po4,    ( -(uptN - lossphyto - grazingN) /qNP )           _CONV_UNIT_)
   _SET_ODE_(self%id_oxy ,   (prodO2   -(respphyto + respzoo) *gammaO2)        _CONV_UNIT_)
   _SET_ODE_(self%id_ldetC,  (faecesC)                                         _CONV_UNIT_)
   if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic , (respphyto + respzoo) _CONV_UNIT_)
   if (_AVAILABLE_(self%id_zbC)) _SET_ODE_(self%id_zbC , (exportC)             _CONV_UNIT_)
   if (_AVAILABLE_(self%id_zbN)) _SET_ODE_(self%id_zbN , (exportN)             _CONV_UNIT_)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_mpbP,     mpbP)          !instantaneous MPB phosphorus assuming redield N:P (mmolP m-3)
   _SET_DIAGNOSTIC_(self%id_par,      parz)          !instantaneous MPB photosynthetically active radiation (W m-2)
   _SET_DIAGNOSTIC_(self%id_Q_N,      Q_N)           !instantaneous MPB N:C quota (molN molC-1)
   _SET_DIAGNOSTIC_(self%id_Q_chl,    Q_chl/12._rk)  !instantaneous MPB CHL:C ratio (gChla gC-1)
   _SET_DIAGNOSTIC_(self%id_Theta_N,  theta)         !instantaneous MPB CHL:N ratio (gChla molN-1)
   _SET_DIAGNOSTIC_(self%id_expCProd, exportC)       !instantaneous MPB export production carbon (mmolC m-3 d-1)
   _SET_DIAGNOSTIC_(self%id_expNProd, exportN)       !instantaneous MPB export production nitrogen (mmolN m-3 d-1)
   ! production rates
   _SET_DIAGNOSTIC_(self%id_SPR,      PP)            !instantaneous MPB specific photosynthesis rate (d-1)
   _SET_DIAGNOSTIC_(self%id_PrimProd, prod)          !instantaneous MPB primary production rate (molC m-3 d-1)
   _SET_DIAGNOSTIC_(self%id_NPP,      NPP)           !instantaneous MPB net primary production rate (molC m-3 d-1)
   _SET_DIAGNOSTIC_(self%id_SGR,      SGR)           !instantaneous MPB specific growth rate (d-1)
   _SET_DIAGNOSTIC_(self%id_TGR,      TGR)           !instantaneous MPB total growth rate (d-1)
   _SET_DIAGNOSTIC_(self%id_SMR,      GRZ)           !instantaneous MPB specific mortality rate (d-1)
   _SET_DIAGNOSTIC_(self%id_uptN,     ruptN)         !instantaneous MPB specific N assimilation rate (d-1)
   _SET_DIAGNOSTIC_(self%id_uptP,     ruptP)         !instantaneous MPB specific P assimilation rate (d-1)
   ! limitation factors
   _SET_DIAGNOSTIC_(self%id_ftfac,    tfac)          !instantaneous temperature dependency
   _SET_DIAGNOSTIC_(self%id_fIR,      f_IR)          !instantaneous light limitation factor
   _SET_DIAGNOSTIC_(self%id_flimN,    limN)          !instantaneous DIN limitation factor
   _SET_DIAGNOSTIC_(self%id_flimP,    limP)          !instantaneous DIP limitation factor
   _SET_DIAGNOSTIC_(self%id_flimO2,   limO2)         !instantaneous oxygen limitation factor
   _SET_DIAGNOSTIC_(self%id_fRChl,    f_RChl)        !instantaneous Chla synthesis regulation
   _SET_DIAGNOSTIC_(self%id_fRCN,     f_RCN)         !instantaneous C uptake downregulation due to qNC
   _SET_DIAGNOSTIC_(self%id_fRNC,     f_RNC)         !instantaneous N uptake downregulation due to qNC
   _SET_DIAGNOSTIC_(self%id_fRNP,     f_RNP)         !instantaneous N uptake downregulation due to qNP
   _SET_DIAGNOSTIC_(self%id_fRPN,     f_RPN)         !instantaneous P uptake downregulation due to qNP
   ! temporary for debugging:
   _SET_DIAGNOSTIC_(self%id_mpb_temp, temp_celsius)  !instantaneous MPB temperature (degC)
   _SET_DIAGNOSTIC_(self%id_mpb_din,  DIN)           !instantaneous external DIN (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_no3,  no3)           !instantaneous external no3 (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_nh4,  nh4)           !instantaneous external nh4 (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_po4,  po4)           !instantaneous external po4 (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_totN, totN)          !instantaneous total N (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_totP, totP)          !instantaneous total P (mmolP m-3)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
  end module hzg_mpb
!-----------------------------------------------------------------------

