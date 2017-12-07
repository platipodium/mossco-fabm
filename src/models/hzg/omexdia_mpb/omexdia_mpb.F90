#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_omexdia_mpb --- Fortran 2003 version of OMEXDIA+CNP+MicroPhytoBenthos model
!
! !INTERFACE:
   module hzg_omexdia_mpb
!
! !DESCRIPTION:
!
! The OMEXDIA+CNP+MPB model is based on the model published by Hochard et al EcoMod 2010.
! This model is based on the OMEXDIA model (see Soetard et al. 1996a) and is intended to
! simulate early diagenesis in the sea sediments as well the benthic mircophytofauna.
! The major difference to the original OMEXDIA+MPB model is an added phosphorus cycle
! added by Kai Wirtz and Markus Kreus. Further modifications comprise efficient reaction and
! limitation terms to facilitate simple numerics (kai wirtz).
! OMEXDIA+CNP include an explizit formulation of detritus C:N:P content
! (no longer forced beeing constant) introduced by Markus Kreus.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_hzg_omexdia_mpb
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
   type,extends(type_base_model) :: type_hzg_omexdia_mpb
!     Variable identifiers
      ! for omexdia_p
      type (type_state_variable_id)        :: id_ldetC, id_sdetC, id_ldetN, id_sdetN, id_ldetP
      type (type_state_variable_id)        :: id_no3, id_nh3, id_oxy, id_po4, id_odu
      type (type_state_variable_id)        :: id_dic, id_zbC, id_zbN, id_zbP
      type (type_dependency_id)            :: id_temp
      type (type_diagnostic_variable_id)   :: id_denit, id_adsp
      ! for mpbC
      type (type_state_variable_id)        :: id_mpbCHL, id_mpbC, id_mpbN, id_mpbP, id_eps
      type (type_dependency_id)            :: id_parz, id_porosity
      type (type_diagnostic_variable_id)   :: id_PrimProd, id_par, id_qNC, id_Q_chl, id_ThetaN
      type (type_diagnostic_variable_id)   :: id_expCProd, id_expNProd, id_expPProd
      type (type_diagnostic_variable_id)   :: id_NPP, id_SGR, id_TGR, id_SPR, id_SMR
      ! temporary for debugging
      type (type_diagnostic_variable_id)   :: id_MPB_din, id_MPB_no3, id_MPB_nh4, id_MPB_po4
      type (type_diagnostic_variable_id)   :: id_MPB_temp, id_MPB_totN, id_MPB_totP
      type (type_diagnostic_variable_id)   :: id_fIR, id_flimN, id_flimP, id_flimO2, id_uptN, id_uptP
      type (type_diagnostic_variable_id)   :: id_fRChl, id_fRCN, id_fRNC, id_fRNP, id_fRPN, id_ftfac

!     Model parameters
      !for omexdia_p
      real(rk) :: rLabile, rSemilabile, NCrLdet, NCrSdet, PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox, CprodMax
      real(rk) :: ksO2oduox, ksO2oxic, ksNO3denit, kinO2denit, kinNO3anox, kinO2anox
      ! for mpbC
      logical  :: MPhytoBenOn   ! use MicroPhytoBenthos (MPB)
      real(rk) :: Tref, Q10, alpha, mu_max, resp0, zeta_NO3, zeta_NH4, kPO4, kNO3, kNH4, kInNH4
      real(rk) :: theta_max, QNCmin, QNCmax, QNCupt, QPNmin, QPNmax, QPNupt, sigma_NC, sigma_PN
      real(rk) :: kO2resp, gamma, fracEPS, degrEPS, graz, aeff, exud, rzoo
      ! for all
      logical  :: use_dic, use_zbC, use_zbN, use_zbP

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do

   end type type_hzg_omexdia_mpb
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the OMEXDIA+P+MPB model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the omexdia namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hzg_omexdia_mpb),intent(inout),target  :: self
   integer,                       intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Markus Kreus, Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
!!------- Initial values of model omexdia_p -------
   real(rk)  :: ldetC_init    ! labile detritus carbon (fast decay)         (        mmolC  m-3 solid)
   real(rk)  :: sdetC_init    ! semilabile detritus carbon (slow decay)     (        mmolC  m-3 solid)
   real(rk)  :: ldetN_init    ! labile detritus nitrogen (fast decay)       (        mmolN  m-3 solid)
   real(rk)  :: sdetN_init    ! semilabile detritus nitrogen (slow decay)   (        mmolN  m-3 solid)
   real(rk)  :: ldetP_init    ! labile detritus phosphorus                  (        mmolP  m-3 solid)
   real(rk)  :: no3_init      ! dissolved nitrate                           (        mmolN  m-3 liquid)
   real(rk)  :: nh3_init      ! dissolved ammonium                          (        mmolN  m-3 liquid)
   real(rk)  :: po4_init      ! dissolved phosphate                         (        mmolP  m-3 liquid)
   real(rk)  :: oxy_init      ! dissolved oxygen                            (        mmolO2 m-3 solid)
   real(rk)  :: odu_init      ! dissolved reduced substances                (        mmolP  m-3 liquid)
!!------- Parameters for model omexdia_p -------
   real(rk)  :: rLabile       ! decay rate labile detritus (fast decay)     (        d-1)
   real(rk)  :: rSemilabile   ! decay rate semilabile detritus (slow decay) (        d-1)
   real(rk)  :: NCrLdet       ! N/C ratio labile detritus (fast decay)      (0.1509  molN molC-1)
   real(rk)  :: NCrSdet       ! N/C ratio semilabile detritus (slow decay)  (0.1333  molN molC-1)
   real(rk)  :: PAds          ! Adsorption coeff phosphorus                 (        -)
   real(rk)  :: PAdsODU       ! PO4-Fe dissolution threshold in terms of [FeS]/ODU ( -)
   real(rk)  :: NH3Ads        ! Adsorption coeff ammonium                   (1.3     -)
   real(rk)  :: CprodMax=48.0 ! Max C-degrad rate for numeric stability     (        -)
   real(rk)  :: rnit          ! Max nitrification rate                      (20.0    d-1)
   real(rk)  :: ksO2nitri     ! half-sat for O2  limitation in nitrification                 (1.0  mmolO2 m-3)
   real(rk)  :: rODUox        ! Max rate oxidation rate of ODU              (20.0    d-1)
   real(rk)  :: ksO2oduox     ! half-sat for O2  limitation in oxidation of reduced nutriens (1.0  mmolO2 m-3)
   real(rk)  :: ksO2oxic      ! half-sat for O2  limitation in oxic mineralization           (3.0  mmolO2 m-3)
   real(rk)  :: ksNO3denit    ! half-sat for NO3 limitation in denitrification               (6.18 mmolN  m-3 := 30.0 mmolNO3 m-3)
   real(rk)  :: kinO2denit    ! half-sat for O2  inhibition in denitrification               (10.0 mmolO2 m-3)
   real(rk)  :: kinNO3anox    ! half-sat for NO3 inhibition in anoxic mineralization         (1.03 mmolN  m-3 := 5.0 mmolNO3 m-3)
   real(rk)  :: kinO2anox     ! half-sat for O2  inhibition in anoxic mineralization         (5.0  mmolO2 m-3)
!!------- Initial values of model mpb -------
   !Note: At Hochard et al (2010) units for mpbCHL was given as mgChla/l.
   !      MK assumed that to be a typo and defined its unit to mgChla m-3, which
   !      correspondents much better to other literature refering to the "Geider" model.
   real(rk)  :: mpbCHL_init   ! MicroPhytoBenthos chlorophyll               (mgChla m-3)
   real(rk)  :: mpbC_init     ! MicroPhytoBenthos carbon                    (mmolC m-3)
   real(rk)  :: mpbN_init     ! MicroPhytoBenthos nitrogen                  (mmolN m-3)
   real(rk)  :: mpbP_init     ! MicroPhytoBenthos phosphorus                (mmolP m-3)
   real(rk)  :: eps_init      ! Extracellular Polymeric Substances          (mmolC m-3)
!!------- Switch for MicroPhytoBenthos model -------
   logical   :: MPhytoBenOn   ! use MicroPhytoBenthos (MPB)
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
   character(len=attribute_length) :: dic_variable   = ''  ! dissolved inorganic carbon
   character(len=attribute_length) :: zbC_variable   = ''  ! zoobenthos carbon
   character(len=attribute_length) :: zbN_variable   = ''  ! zoobenthos nitrogen
   character(len=attribute_length) :: zbP_variable   = ''  ! zoobenthos phosphorus

   namelist /hzg_omexdia_cnp/  &
          rLabile, rSemilabile, NCrLdet, NCrSdet, PAds, PAdsODU, NH3Ads,       &
          CprodMax, rnit, ksO2nitri, rODUox, ksO2oduox, ksO2oxic, ksNO3denit,  &
          kinO2denit, kinNO3anox, kinO2anox,                                   &
          ldetC_init, sdetC_init, ldetN_init, sdetN_init, ldetP_init,          &
          no3_init, nh3_init, po4_init, oxy_init, odu_init

   namelist /hzg_omexdia_mpb/ MPhytoBenOn

   namelist /hzg_omexdia_mpb_dependencies/  &
          dic_variable, zbC_variable, zbN_variable, zbP_variable

   namelist /hzg_mpb_cnp/ &
          Tref, Q10, alpha, mu_max, resp0, zeta_NO3, zeta_NH4, kPO4, kNO3, kNH4, kInNH4, &
          theta_max, QNCmin, QNCmax, QNCupt, QPNmin, QPNmax, QPNupt, sigma_NC, sigma_PN, &
          kO2resp, gamma, fracEPS, degrEPS, graz, aeff, exud, rzoo, &
          mpbCHL_init, mpbC_init, mpbN_init, mpbP_init, eps_init

!EOP
!-----------------------------------------------------------------------
!BOC

   ! Read the namelist
   if (configunit>0) then
     read(configunit, nml=hzg_omexdia_cnp, err=96, end=100)
     read(configunit, nml=hzg_omexdia_mpb, err=97, end=101)
     read(configunit, nml=hzg_omexdia_mpb_dependencies, err=98, end=102)
     read(configunit, nml=hzg_mpb_cnp, err=99, end=103)
   endif

   ! Store parameter values in our own derived type
   self%rLabile    = rLabile
   self%rSemilabile= rSemilabile
   self%NCrLdet    = NCrLdet
   self%NCrSdet    = NCrSdet
   self%PAds       = PAds
   self%PAdsODU    = PAdsODU
   self%NH3Ads     = NH3Ads
   self%CprodMax   = CprodMax
   self%rnit       = rnit
   self%ksO2nitri  = ksO2nitri
   self%rODUox     = rODUox
   self%ksO2oduox  = ksO2oduox
   self%ksO2oxic   = ksO2oxic
   self%ksNO3denit = ksNO3denit
   self%kinO2denit = kinO2denit
   self%kinNO3anox = kinNO3anox
   self%kinO2anox  = kinO2anox

   ! Register state variables
   call self%register_state_variable(self%id_ldetC, 'ldetC', 'mmolC m-3',  &
                                    'detritus labile carbon', ldetC_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_ldetC, 'particulate', .true.)

   call self%register_state_variable(self%id_sdetC, 'sdetC', 'mmolC m-3',  &
                                    'detritus semilabile carbon', sdetC_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_sdetC, 'particulate', .true.)

   call self%register_state_variable(self%id_ldetN, 'ldetN', 'mmolN m-3',  &
                                    'detritus labile nitrogen', ldetN_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_ldetN, 'particulate', .true.)

   call self%register_state_variable(self%id_sdetN, 'sdetN', 'mmolN m-3',  &
                                    'detritus semilabile nitrogen', sdetN_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_sdetN, 'particulate', .true.)

   call self%register_state_variable(self%id_ldetP, 'ldetP', 'mmolP m-3',  &
                                    'detritus phosphorus', ldetP_init, minimum=0.0_rk)
!                                     'detritus labile phosphorus', ldetP_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_ldetP, 'particulate', .true.)

   call self%register_state_variable(self%id_no3, 'no3', 'mmolN m-3',  &
                                    'dissolved nitrate', no3_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_nitrate)
   call self%set_variable_property(self%id_no3, 'particulate', .false.)

   call self%register_state_variable(self%id_nh3, 'nh3', 'mmolN m-3',  &
                                    'dissolved ammonium', nh3_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_ammonium)
   call self%set_variable_property(self%id_nh3, 'particulate', .false.)

   call self%register_state_variable(self%id_po4, 'po4', 'mmolP m-3',  &
                                    'dissolved phosphate', po4_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_phosphate)
   call self%set_variable_property(self%id_po4, 'particulate', .false.)

   call self%register_state_variable(self%id_oxy, 'oxy', 'mmolO2 m-3',  &
                                    'dissolved oxygen', oxy_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_oxy, 'particulate', .false.)

   call self%register_state_variable(self%id_odu,'odu','mmolO2 m-3',  &
                                    'dissolved reduced substances', odu_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_odu,'particulate', .false.)

   if (MPhytoBenOn) then
      ! Store parameter values in our own derived type
      self%MPhytoBenOn = MPhytoBenOn
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
         call self%register_state_dependency(self%id_zbP, 'zbP_variable', 'mmolP m-3',   &
               'ZooBenthos phosphorus', required=.false.)
         call self%request_coupling(self%id_zbP, zbP_variable)
      endif

      ! Register diagnostic variables
!       call self%register_diagnostic_variable(self%id_par,      'MPB_CNP_par',    'W m-2',           &
!             'MPB-CNP photosynthetically active radiation',     output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_qNC,      'MPB_CNP_qNC',    'molN molC-1',     &
            'MPB-CNP nitrogen quota qNC',                      output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_Q_chl,    'MPB_CNP_Q_chl',  'gChla gC-1',      &
            'MPB-CNP CHL-to-C ratio Q_chl',                    output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_ThetaN,   'MPB_CNP_ThetaN', 'gChla molN-1',    &
            'MPB-CNP CHL-to-N ratio ThetaN',                   output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_expCProd, 'MPB_CNP_expCProd', 'mmolC m-3 d-1', &
            'MPB-CNP carbon export (zoobenthos grazing)',      output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_expNProd, 'MPB_CNP_expNProd', 'mmolN m-3 d-1', &
            'MPB-CNP nitrogen export (zoobenthos grazing)',    output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_expPProd, 'MPB_CNP_expPProd', 'mmolP m-3 d-1', &
            'MPB-CNP phosphorus export (zoobenthos grazing)',  output=output_instantaneous)
      ! production rates
      call self%register_diagnostic_variable(self%id_PrimProd, 'MPB_CNP_PP',     'mmolC m-3 d-1',   &
            'MPB-CNP primary production rate PrimProd',        output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_NPP,      'MPB_CNP_NPP',    'mmolC m-3 d-1',   &
            'MPB-CNP net primary production rate NPP',         output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_SGR,      'MPB_CNP_SGR',    'd-1',             &
            'MPB-CNP specific growth rate SGR',                output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_TGR,      'MPB_CNP_TGR',    'd-1',             &
            'MPB-CNP total growth rate TGR',                   output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_SPR,      'MPB_CNP_SPR',    'd-1',             &
            'MPB-CNP specific photosynthesis rate SPR',        output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_SMR,      'MPB_CNP_SMR',    'd-1',             &
            'MPB-CNP specific mortality rate (grazing) SMR',   output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_uptN,     'MPB_CNP_uptN',   'd-1',             &
            'MPB-CNP specific nitrogen assimilation rate uptN', output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_uptP,     'MPB_CNP_uptP',   'd-1',             &
            'MPB-CNP specific phosphorus assimilation rate uptP', output=output_instantaneous)
      ! limitation factors
      call self%register_diagnostic_variable(self%id_ftfac,    'MPB_CNP_ftfac',  '-',               &
            'MPB-CNP ftfac',                                   output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_fIR,      'MPB_CNP_fIR',    '-',               &
            'MPB-CNP fIR',                                     output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_flimN,    'MPB_CNP_flimN',  '-',               &
            'MPB-CNP flimN',                                   output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_flimP,    'MPB_CNP_flimP',  '-',               &
            'MPB-CNP flimP',                                   output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_flimO2,   'MPB_CNP_flimO2', '-',               &
            'MPB-CNP flimO2',                                  output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_fRChl,    'MPB_CNP_fRCHl',  '-',               &
            'MPB-CNP fRChl',                                   output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_fRCN,     'MPB_CNP_fRCN',   '-',               &
            'MPB-CNP fRCN',                                    output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_fRNC,     'MPB_CNP_fRNC',   '-',               &
            'MPB-CNP fRNC',                                    output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_fRNP,     'MPB_CNP_fRNP',   '-',               &
            'MPB-CNP fRNP',                                    output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_fRPN,     'MPB_CNP_fRPN',   '-',               &
            'MPB-CNP fRPN',                                    output=output_instantaneous)
   endif

   ! temporary for debugging:
   call self%register_diagnostic_variable(self%id_par,     'par',         'W m-2',      &
        'MPB-CNP photosynthetically active radiation',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_DIN, 'DIN',         'mmolN m-3',  &
        'MPB-CNP DIN',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_no3, 'MPB_CNP_no3', 'mmolN m-3',  &
        'MPB-CNP no3',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_nh4, 'MPB_CNP_nh4', 'mmolN m-3',  &
        'MPB-CNP nh4',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_po4, 'MPB_CNP_po4', 'mmolP m-3',  &
        'MPB-CNP po4',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_totN,'MPB_CNP_totN','mmolN m-3',  &
        'MPB-CNP totN', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_totP,'MPB_CNP_totP','mmolP m-3',  &
        'MPB-CNP totP', output=output_instantaneous)


   ! Register dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_parz, standard_variables%downwelling_photosynthetic_radiative_flux)

   return

96 call self%fatal_error('hzg_omexdia_mpb_initialize','Error reading namelist hzg_omexdia_p.')
97 call self%fatal_error('hzg_omexdia_mpb_initialize','Error reading namelist hzg_omexdia_mpb.')
98 call self%fatal_error('hzg_omexdia_mpb_initialize','Error reading namelist hzg_omexdia_mpb_dependencies.')
99 call self%fatal_error('hzg_omexdia_mpb_initialize','Error reading namelist hzg_mpb.')

100 call self%fatal_error('hzg_omexdia_mpb_initialize','Namelist hzg_omexdia_p was not found.')
101 call self%fatal_error('hzg_omexdia_mpb_initialize','Namelist hzg_omexdia_mpb was not found.')
102 call self%fatal_error('hzg_omexdia_mpb_initialize','Namelist hzg_omexdia_mpb_dependencies was not found.')
103 call self%fatal_error('hzg_omexdia_mpb_initialize','Namelist hzg_mpb was not found.')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of OMEXDIA+P model
!
! !INTERFACE:
   subroutine do(self, _ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_hzg_omexdia_mpb),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Markus Kreus, Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
!!------- for model omexdia_p -------
   real(rk), parameter :: gammaNO3 = 0.8_rk    ! molN  used per molC in denitrification     !RH 0.8-> ~104/106?
   real(rk), parameter :: gammaO2  = 1.0_rk    ! molO2 used per molC in oxic mineralization !RH 1.0->150/106*OxicMin (if [oxy]=mmolO2 m-3)
   real(rk), parameter :: gammaNH4 = 2.0_rk    ! molO2 needed to oxidize one molNH4 in nitrification
   real(rk), parameter :: T0       = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk), parameter :: Q10b     = 1.5_rk    ! q10 temperature coefficient (source?)
   real(rk), parameter :: relaxO2  = 0.04_rk   ! (source?)
   real(rk) :: ldetC, sdetC, ldetN, sdetN, ldetP, no3, nh4, po4, oxy, odu
   real(rk) :: temp_celsius, temp_kelvin, f_T, E_a
   real(rk) :: radsP, Oxicminlim, Denitrilim, Anoxiclim, Rescale, rP
   real(rk) :: CprodL, CprodS, Cprod, NprodL, NprodS, Nprod, Pprod
   real(rk) :: AnoxicMin, Denitrific, OxicMin, Nitri, OduDepo, OduOx, pDepo
!!------- for model mpb -------
   real(rk), parameter :: zero = 0.0_rk, one = 1.0_rk, two = 2.0_rk, TINY = 1.0e-3
   real(rk) :: mpbC, mpbN, mpbP, mpbCHL, eps, DIN
   real(rk) :: tfac, parz, porosity
   real(rk) :: qNC, qPN, Q_chl, theta, mu_C, mu_N, mu_P, mu_chl, limN, limP, limO2
   real(rk) :: Pmax, f_IR, f_RCN, f_RPC, f_RPN, f_RNP, f_RNC, frac_NH4
   real(rk) :: prod, prodChl, prodO2, respphyto, lossphytN, lossphytP
   real(rk) :: uptP, uptN, uptNO3, uptNH4, prodEPS, degrEPS, limGraz
   real(rk) :: grazingC, grazingN, grazingP, grazingChl, faecesC, faecesN, faecesP
   real(rk) :: respZoo, exudZoN, exudZoP, exportC, exportN, exportP
   real(rk) :: NPP, SGR, TGR, GRZ, ruptN, ruptP, x1, x2, fac, f_RChl, limNO3, limNH4, totN, totP
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp, temp_celsius) ! sediment-water temperature (degC)
   _GET_(self%id_ldetC, ldetC)       ! labile detritus C          (mmolC m-3)
   _GET_(self%id_sdetC, sdetC)       ! semilabile detritus C      (mmolC m-3)
   _GET_(self%id_ldetN, ldetN)       ! labile detritus N          (mmolN m-3)
   _GET_(self%id_sdetN, sdetN)       ! semilabile detritus N      (mmolN m-3)
   _GET_(self%id_ldetP, ldetP)       ! labile detritus P          (mmolP m-3)
   _GET_(self%id_oxy, oxy)           ! dissolved oxygen           (mmolO2 m-3)
   _GET_(self%id_odu, odu)           ! dissolved oxygen demand units (mmolO2 m-3)
   _GET_(self%id_no3, no3)           ! dissolved nitrate          (mmolN m-3)
   _GET_(self%id_nh3, nh4)           ! dissolved ammonium         (mmolN m-3)
   _GET_(self%id_po4, po4)           ! dissolved phosphate        (mmolP m-3)
   DIN = no3+nh4
   ldetN = ldetC*self%NCrLdet
   sdetN = sdetC*self%NCrSdet
   totN = DIN !+ ldetN + sdetN
   totP = po4 !+ detP
   if (self%MPhytoBenOn) then
      _GET_(self%id_parz,   parz)    ! sediment light in W m-2
      _GET_(self%id_mpbCHL, mpbCHL)  ! MicroPhytoBenthos chlorophyll in mgChla m-3
      _GET_(self%id_mpbC,   mpbC)    ! MicroPhytoBenthos carbon in mmolC m-3
      _GET_(self%id_mpbN,   mpbN)    ! MicroPhytoBenthos nitrogen in mmolN m-3
      _GET_(self%id_mpbP,   mpbP)    ! MicroPhytoBenthos phosphorus in mmolP m-3
      _GET_(self%id_eps,    eps)     ! Extracellular Polymeric Substances in mmolC m-3
      totN = totN + mpbN
      totP = totP + mpbP
   end if

   ! temperature dependency
   temp_kelvin = 273.15_rk + temp_celsius
   E_a = 0.1_rk*log(Q10b)*T0*(T0+10.0_rk)
   f_T = 1.0_rk*exp(-E_a*(1.0_rk/temp_kelvin - 1.0_rk/T0)) ! temperature factor (-)

   if (2*oxy < -self%kinO2anox) oxy = -self%kinO2anox/2    ! limit anoxic conditions

   ! limitation terms
   Oxicminlim = oxy/(oxy+self%ksO2oxic+relaxO2*(nh4+odu))                                ! (-)
   Denitrilim = (1.0_rk-oxy/(oxy+self%kinO2denit)) * NO3/(no3+self%ksNO3denit)           ! (-)
   Anoxiclim  = (1.0_rk-oxy/(oxy+self%kinO2anox)) * (1.0_rk-no3/(no3+self%kinNO3anox))   ! (-)
   if(Oxicminlim - 1E-3 < -Denitrilim-Anoxiclim) Oxicminlim = 1E-3-Denitrilim-Anoxiclim  ! (-)
   Rescale    = 1.0_rk/(Oxicminlim+Denitrilim+Anoxiclim)     !Soetaert eq 3.4

   CprodL = f_T * self%rLabile     * ldetC   ! (mmolC m-3 d-1)
   CprodS = f_T * self%rSemilabile * sdetC   ! (mmolC m-3 d-1)

! assume upper reactive surface area for POC hydrolysis (introduced by Kai Wirtz ?)
   if (CprodS > self%CprodMax) CprodS = self%CprodMax
   if (CprodL > self%CprodMax) CprodL = self%CprodMax
   Cprod  = CprodL + CprodS                               ! (mmolC m-3 d-1)

   NprodL = CprodL * (ldetN/ldetC)                        ! (mmolN m-3 d-1)
   NprodS = CprodS * (sdetN/sdetC)                        ! (mmolN m-3 d-1)
   Nprod  = NprodL + NprodS                               ! (mmolN m-3 d-1)

! PO4-adsorption ceases when critical capacity is reached
! [FeS] approximated by ODU
! TODO: temperature dependency
   radsP  = self%PAds  * po4 * 1.0_rk/(1.0_rk+exp(-5.0_rk+(odu-oxy)/self%PAdsODU)) ! (?)
   rP     = f_T * self%rLabile * 2.0_rk ! (1.0_rk - Oxicminlim)  ! (d-1)
   Pprod  = rP * ldetP                                           ! (mmolP m-3 d-1)

! Oxic mineralisation, denitrification, anoxic mineralisation
! then the mineralisation rates
   OxicMin    = Cprod*Oxicminlim*Rescale  ! carbon oxidized by O2              (mmolC m-3 d-1) !Soetaert eq 3.1
   Denitrific = Cprod*Denitrilim*Rescale  ! carbon oxidized by denitrification (mmolC m-3 d-1) !Soetaert eq 3.2
   AnoxicMin  = Cprod*Anoxiclim *Rescale  ! carbon oxidized by other oxidants  (mmolC m-3 d-1) !Soetaert eq 3.3

! Ammonium (NH4) nitrified (mmolN m-3 d-1) !Soetaert eq 3.7
   Nitri      = f_T * self%rnit   * nh4 * oxy/(oxy + self%ksO2nitri + relaxO2*(ldetC + odu)) !Soetaert eq 3.7, but where does the relax* term comes from?

! Oxygen consumed in the oxidation of reduced substances (mmolO2 m-3 d-1) !Soetaert eq 3.6
   OduOx      = f_T * self%rODUox * odu * oxy/(oxy + self%ksO2oduox + relaxO2*(nh4 + ldetC)) !Soetaert eq 3.6, but where does the relax* term comes from?
   if (OduOx > self%CprodMax) OduOx = self%CprodMax

! ODU deposition
   !pDepo      = min(1.0_rk,0.233_rk*(wDepo)**0.336_rk )
   pDepo      = 0.0_rk
   OduDepo    = AnoxicMin*pDepo           ! ODU deposited as solids (in O2-equivalents) (mmolO2 m-3 d-1) !Soetaert eq 3.5

   if (self%MPhytoBenOn) then
      ! light climate (obsolete)
      !par    = self%PAR0max * parz
      !print*,'par#386: ',par, parz
      !k      = self%k0 + self%Achla * mpbCHL

      !----- intracellular N:C:Chl stoichiometry
      qNC    = mpbN   / mpbC      ! molN  molC-1
      qPN    = mpbP   / mpbN      ! molP  molN-1
      Q_chl  = mpbCHL / mpbC      ! gChla molC-1
      theta  = mpbCHL / mpbN      ! gChla molN-1
      !print*,'mpb_cnp#425 ', mpbC, mpbN, mpbP, mpbChl
      !print*,'mpb_cnp#425 ', qNC, qPN, Q_chl, theta

      !----- microphytes temperature dependency
      tfac = exp( log(self%Q10)*((temp_celsius-10._rk)/10._rk) )          ! (-)

      !----- oxygen limitation factor
      limO2 = max( zero, oxy / (oxy+self%KO2resp) )                       ! (-)

      !----- photosynthesis and carbon assimilation
      !f_RCN = max( zero, (qNC-self%QNCmin)/(self%QNCmax-self%QNCmin) )    ! (-) [vgl. GeiderEtAl(1998), SchartauEtAl(2007), KreusEtAl(2015)
      f_RCN  = max( zero, one - self%QNCmin/qNC)      ! (-) [vgl. HochardEtAl(2010) eq. 19]
      Pmax = self%mu_max *tfac *min(one,f_RCN)                            ! (d-1)
      ! - for comparism: Pmax as formulated in i.e. Baumert, Pahlow etc.
      !Pmax = self%mu_max *tfac *(one- self%QNCmin/qNC)                   ! (d-1)

      !print*,'#435 ', Pmax, tfac, f_RCN

      f_IR = zero
      mu_C = zero
      if (Pmax > zero) then
        ! light limitation
        f_IR = one -exp(-self%alpha *Q_chl *parz /Pmax)                   ! (-)
        ! carbon assimilation rate
        mu_C = Pmax *f_IR                                                 ! (d-1)
      endif
      prod = mu_C * mpbC                                                  ! (mmolC m-3 d-1)

      !print*,'#443 ', prod, Pmax, f_IR, limO2

      !----- phosphorus assimilation and regulation
      limP  = max( zero, po4/(self%KPO4 + po4) )                          ! (-)
      ! down-regulation of P-uptake if cell is P-replete
      f_RPC = one                                                         ! (-)
      ! down-regulation of P-uptake due to higher P:N-quota caused by N-limitation
      f_RPN = one/( one + exp(-self%sigma_PN *(self%QPNmax-qPN)) )        ! (-)
      !f_RPN = max( zero, two/( one + exp(-self%sigma_PN *(self%QPNmax-qPN)) ) -one ) ! (-)
      ! P-assimilation is not directly related to amount of proteins (as "expressed" by qNC)
      mu_P = self%QPNupt *self%mu_max *tfac *f_RPC *f_RPN *limP           ! (mmolP mmolN-1 d-1)
      ! net P uptake rate
      uptP = mu_P * mpbN                                                  ! (mmolP d-1)

      !----- nitrogen assimilation and regulation
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

      f_RNP = max( zero, one-(self%QPNmin /qPN) )                         ! (-) (Droop)
      f_RNC = one /( one + exp(-self%sigma_NC *(self%QNCmax-qNC)) )       ! (-)
      !f_RNC = max( zero, two/( one + exp(-self%sigma_NC *(self%QNCmax-qNC)) ) -one ) ! (-)
      ! carbon specific nitrogen uptake rate
      mu_N = self%QNCupt *self%mu_max *tfac *f_RNP *f_RNC *limN           ! (mmolN mmolC-1 d-1) (vgl. Geider)
      if (mu_N .lt. TINY ) mu_N = zero                                    ! (mmolN mmolC-1 d-1)
      ! N uptake rates
      uptNH4  = (      frac_NH4)*mu_N*mpbC                                ! (mmolN d-1)
      uptNO3  = (one - frac_NH4)*mu_N*mpbC                                ! (mmolN d-1)
      uptN    = uptNO3 + uptNH4                                           ! (mmolN d-1)
      !uptP = uptN/16.

      !----- Carbohydrate exudation:
      prodEPS = self%fracEPS *prod                                        ! (mmolC m-3 d-1)
      !CprodEPS = sqrt(self%rLdet*self%rSdet) * eps  ! Source of this formulation? Kai Wirtz ??
      degrEPS = max( zero, f_T *self%degrEPS *(eps-TINY) )                ! (mmolC m-3 d-1)

      !----- chlorophyll synthesis
      if (mu_C .le. TINY) then
          mu_chl = mu_N *self%theta_max                                    ! (mgChla mmolC-1)
      else
          mu_chl = mu_N *self%theta_max *mu_C/(self%alpha *Q_chl *parz)    ! (mgChla mmolC-1)
      endif
      fac    = max(zero, one - theta/self%theta_max)                      ! (-)   [vgl. HochardEtAl(2010) eq. 23]
      f_RChl = fac/(fac + 0.05_rk)                                        ! (-)   [vgl. HochardEtAl(2010) eq. 23]
      prodchl = mu_chl *f_RChl *mpbC                                      ! (mgChla m-3 d-1)

      !----- oxygen production and respiration / maintenance loss:
      prodO2  = self%gamma * prod                                         ! (mmolO2 m-3 d-1)
      respphyto = self%resp0 *tfac *limO2 *mpbC                           ! (mmolC m-3 d-1)
      if (self%resp0>zero) then
          respphyto = respphyto + uptNO3 *self%zeta_NO3 + uptNH4 *self%zeta_NH4 ! (mmolC m-3 d-1)
      endif
      ! limit respiration according to N/P-losses if numbers get very low
      if (mu_C .lt. TINY) respphyto = self%resp0 *tfac *limO2 *mpbC       ! (mmolC m-3 d-1)
      lossphytN = self%resp0 *tfac *mpbN                                  ! (mmolN m-3 d-1)
      lossphytP = self%resp0 *tfac *mpbP                                  ! (mmolN m-3 d-1)

      !----- Zoobenthos grazing and associated processes:
      limGraz    = max( zero, (mpbC - TINY) /mpbC )
      grazingC   = self%graz *limO2 *limGraz *mpbC                        ! (mmolC m-3 d-1)
      grazingN   = self%graz *limO2 *limGraz *mpbN                        ! (mmolN m-3 d-1)
      grazingP   = self%graz *limO2 *limGraz *mpbP                        ! (mmolP m-3 d-1)
      grazingChl = self%graz *limO2 *limGraz *mpbCHL                      ! (mgChla m-3 d-1)
      faecesC    = (one -self%aeff) * grazingC                            ! (mmolC m-3 d-1)
      faecesN    = (one -self%aeff) * grazingN                            ! (mmolN m-3 d-1)
      faecesP    = (one -self%aeff) * grazingP                            ! (mmolN m-3 d-1)
      respZoo    = self%rzoo * (grazingC - faecesC)                       ! (mmolC m-3 d-1)
      exudZoN    = self%exud * (grazingN - faecesN)                       ! (mmolN m-3 d-1)
      exudZoP    = self%exud * (grazingP - faecesP)                       ! (mmolP m-3 d-1)
      exportC    = grazingC - faecesC - respZoo ! exported carbon   (open closure term)   ! (mmolC m-3 d-1)
      exportN    = grazingN - faecesN - exudZoN ! exported nitrogen (open closure term)   ! (mmolN m-3 d-1)
      exportP    = grazingP - faecesP - exudZoP ! exported nitrogen (open closure term)   ! (mmolN m-3 d-1)

      !----- Diagnostic parameters
      NPP =  prod - respphyto                                             ! (mmolC m-3 d-1)
      SGR = (prod - respphyto - prodeps ) /mpbC                           ! (d-1)
      TGR = (prod - grazingC - respphyto - prodeps) /mpbC                 ! (d-1)
      GRZ = (grazingC) /mpbC                                              ! (d-1)

      ruptN = zero
      if (uptN>TINY) ruptN = uptN/mpbN                                    ! (d-1)

      ruptP = zero
      if (uptP>TINY) ruptP = uptP/mpbP                                    ! (d-1)
      !if (uptN>TINY) ruptP = ruptN/16.                                    ! (d-1)

      if (mpbN<TINY) lossphytN = zero
      if (mpbP<TINY) lossphytP = zero

   else
!   end if
      prodChl    = 0.0_rk
      prod       = 0.0_rk
      uptN       = 0.0_rk
      uptP       = 0.0_rk
      uptNH4     = 0.0_rk
      uptNO3     = 0.0_rk
      respphyto  = 0.0_rk
      lossphytN  = 0.0_rk
      lossphytP  = 0.0_rk
      grazingChl = 0.0_rk
      grazingC   = 0.0_rk
      grazingN   = 0.0_rk
      grazingP   = 0.0_rk
      faecesC    = 0.0_rk
      faecesN    = 0.0_rk
      faecesP    = 0.0_rk
      respzoo    = 0.0_rk
      exudZoN    = 0.0_rk
      exudZoP    = 0.0_rk
      prodEPS    = 0.0_rk
      degrEPS    = 0.0_rk
      prodO2     = 0.0_rk
   end if

#define _CONV_UNIT_ *one_pr_day
! reaction rates
   _SET_ODE_(self%id_ldetC, (-CprodL)                                      _CONV_UNIT_)  ! (mmolC  m-3 d-1)
   _SET_ODE_(self%id_sdetC, (-CprodS)                                      _CONV_UNIT_)  ! (mmolC  m-3 d-1)
   _SET_ODE_(self%id_ldetN, (-NprodL)                                      _CONV_UNIT_)  ! (mmolN  m-3 d-1)
   _SET_ODE_(self%id_sdetN, (-NprodS)                                      _CONV_UNIT_)  ! (mmolN  m-3 d-1)
   _SET_ODE_(self%id_ldetP, (radsP - Pprod - self%NH3Ads*ldetP)            _CONV_UNIT_)  ! (mmolP  m-3 d-1)
   _SET_ODE_(self%id_no3,  (-Denitrific*gammaNO3 + Nitri)                  _CONV_UNIT_)  ! (mmolN  m-3 d-1)
!   _SET_ODE_(self%id_nh3,  (Nprod - Nitri) / (1.0_rk + self%NH3Ads)       _CONV_UNIT_)  ! (mmolN  m-3 d-1)
   _SET_ODE_(self%id_nh3,  (Nprod - Nitri)                                 _CONV_UNIT_)  ! (mmolN  m-3 d-1)
   _SET_ODE_(self%id_po4,  (Pprod - radsP)                                 _CONV_UNIT_)  ! (mmolP  m-3 d-1)
   _SET_ODE_(self%id_oxy,  (-OxicMin*gammaO2     - Nitri*gammaNH4 - OduOx) _CONV_UNIT_)  ! (mmolO2 m-3 d-1)
   _SET_ODE_(self%id_odu,  (AnoxicMin - OduOx - OduDepo)                   _CONV_UNIT_)  ! (mmolO2 m-3 d-1)
   if (self%MPhytoBenOn) then
      _SET_ODE_(self%id_mpbCHL, (prodChl - lossphytN *theta    - grazingChl)       _CONV_UNIT_)
      _SET_ODE_(self%id_mpbC,   (prod    - respphyto - prodeps - grazingC)         _CONV_UNIT_)
      _SET_ODE_(self%id_mpbN,   (uptN    - lossphytN           - grazingN)         _CONV_UNIT_)
      _SET_ODE_(self%id_mpbP,   (uptP    - lossphytP           - grazingP)         _CONV_UNIT_)
      _SET_ODE_(self%id_eps,    (prodEPS - degrEPS)                                _CONV_UNIT_)
      ! external dependencies
      _SET_ODE_(self%id_no3,    (-uptNO3)                                          _CONV_UNIT_)
      _SET_ODE_(self%id_nh3,    (-uptNH4 + lossphytN + exudZoN)                    _CONV_UNIT_)
      _SET_ODE_(self%id_po4,    (-uptP   + lossphytP + exudZoP)                    _CONV_UNIT_)
      _SET_ODE_(self%id_oxy,    (prodO2  -(respphyto + respZoo) *gammaO2)          _CONV_UNIT_)
      _SET_ODE_(self%id_ldetC,  (faecesC)                                          _CONV_UNIT_)
      _SET_ODE_(self%id_ldetN,  (faecesN)                                          _CONV_UNIT_)
      _SET_ODE_(self%id_ldetP,  (faecesP)                                          _CONV_UNIT_)
      ! If externally maintained variables are present, change the pools accordingly
      if (_AVAILABLE_(self%id_dic))   _SET_ODE_(self%id_dic, (respphyto + respZoo) _CONV_UNIT_)
      if (_AVAILABLE_(self%id_zbC))   _SET_ODE_(self%id_zbC, (exportC)             _CONV_UNIT_)
      if (_AVAILABLE_(self%id_zbN))   _SET_ODE_(self%id_zbN, (exportN)             _CONV_UNIT_)
      if (_AVAILABLE_(self%id_zbP))   _SET_ODE_(self%id_zbP, (exportP)             _CONV_UNIT_)
   end if

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit,       Denitrific*gammaNO3)  !last denitrification rate
   _SET_DIAGNOSTIC_(self%id_adsp,        radsP)                !instantaneous phosphate adsorption

   if (self%MPhytoBenOn) then
      !_SET_DIAGNOSTIC_(self%id_par,      parz)           !instantaneous MPB photosynthetically active radiation (W m-2)
      _SET_DIAGNOSTIC_(self%id_qNC,      qNC)            !instantaneous MPB N:C quota (molN molC-1)

      _SET_DIAGNOSTIC_(self%id_Q_chl,    Q_chl/12._rk)   !instantaneous MPB CHL:C ratio (gChla gC-1)
      _SET_DIAGNOSTIC_(self%id_ThetaN,   theta)          !instantaneous MPB CHL:N ratio (gChla molN-1)
      _SET_DIAGNOSTIC_(self%id_expCProd, exportC)        !instantaneous MPB export production carbon (mmolC m-3 d-1)
      _SET_DIAGNOSTIC_(self%id_expNProd, exportN)        !instantaneous MPB export production nitrogen (mmolN m-3 d-1)
      _SET_DIAGNOSTIC_(self%id_expPProd, exportP)        !instantaneous MPB export production nitrogen (mmolN m-3 d-1)
      ! production rates
      _SET_DIAGNOSTIC_(self%id_SPR,      mu_C)           !instantaneous MPB specific photosynthesis rate (d-1)
      _SET_DIAGNOSTIC_(self%id_PrimProd, prod)           !instantaneous MPB primary production rate (molC m-3 d-1)
      _SET_DIAGNOSTIC_(self%id_NPP,      NPP)            !instantaneous MPB net primary production rate (molC m-3 d-1)
      _SET_DIAGNOSTIC_(self%id_SGR,      SGR)            !instantaneous MPB specific growth rate (d-1)
      _SET_DIAGNOSTIC_(self%id_TGR,      TGR)            !instantaneous MPB total growth rate (d-1)
      _SET_DIAGNOSTIC_(self%id_SMR,      GRZ)            !instantaneous MPB specific mortality rate (d-1)
      _SET_DIAGNOSTIC_(self%id_uptN,     ruptN)          !instantaneous MPB specific N assimilation rate (d-1)
      _SET_DIAGNOSTIC_(self%id_uptP,     ruptP)          !instantaneous MPB specific P assimilation rate (d-1)
      ! limitation factors
      _SET_DIAGNOSTIC_(self%id_ftfac,    tfac)           !instantaneous temperature dependency
      _SET_DIAGNOSTIC_(self%id_fIR,      f_IR)           !instantaneous light limitation factor
      _SET_DIAGNOSTIC_(self%id_flimN,    limN)           !instantaneous DIN limitation factor
      _SET_DIAGNOSTIC_(self%id_flimP,    limP)           !instantaneous DIP limitation factor
      _SET_DIAGNOSTIC_(self%id_flimO2,   limO2)          !instantaneous oxygen limitation factor
      _SET_DIAGNOSTIC_(self%id_fRChl,    f_RChl)         !instantaneous Chla synthesis regulation
      _SET_DIAGNOSTIC_(self%id_fRCN,     f_RCN)          !instantaneous C uptake downregulation due to qNC
      _SET_DIAGNOSTIC_(self%id_fRNC,     f_RNC)          !instantaneous N uptake downregulation due to qNC
      _SET_DIAGNOSTIC_(self%id_fRNP,     f_RNP)          !instantaneous N uptake downregulation due to qNP
      _SET_DIAGNOSTIC_(self%id_fRPN,     f_RPN)          !instantaneous P uptake downregulation due to qNP
   endif

   ! temporary for debugging:
   _SET_DIAGNOSTIC_(self%id_mpb_temp, temp_celsius)   !instantaneous MPB temperature (degC)
   _SET_DIAGNOSTIC_(self%id_mpb_din,  DIN)            !instantaneous external DIN (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_no3,  no3)            !instantaneous external no3 (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_nh4,  nh4)            !instantaneous external nh4 (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_po4,  po4)            !instantaneous external po4 (mmolP m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_totN, totN)           !instantaneous total N (mmolN m-3)
   _SET_DIAGNOSTIC_(self%id_mpb_totP, totP)           !instantaneous total P (mmolP m-3)
   _SET_DIAGNOSTIC_(self%id_par,      parz)          !instantaneous MPB photosynthetically active radiation (W m-2)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
  end module hzg_omexdia_mpb
!-----------------------------------------------------------------------

