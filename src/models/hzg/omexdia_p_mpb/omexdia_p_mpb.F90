#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_omexdia_p_mpb --- Fortran 2003 version of OMEXDIA+P biogeochemical model
!
! !INTERFACE:
   module hzg_omexdia_p_mpb
!
! !DESCRIPTION:
!
! The OMEXDIA+P+MPB model is based on the OMEXDIA model (see Soetard et al. 1996a)
! and is intended to simulate early diagenesis in the sea sediments. The major
! difference to the original OMEXDIA is an added phosphorus cycle.
! P-cycle is added by kai wirtz
! efficient reaction and limitation terms to facilitate simple numerics added by kai wirtz
! MPB model from Hochard et al EcoMod 2010 added by kai wirtz
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_hzg_omexdia_p_mpb
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day  = 1.0_rk / secs_pr_day
!
! !REVISION HISTORY:!
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hzg_omexdia_p_mpb
!     Variable identifiers
      ! for omexdia_p
      type (type_state_variable_id)        :: id_ldetC, id_sdetC, id_detP
      type (type_state_variable_id)        :: id_no3, id_nh3, id_oxy, id_po4, id_odu
      type (type_state_variable_id)        :: id_dic, id_zbC, id_zbN
      type (type_dependency_id)            :: id_temp
      type (type_diagnostic_variable_id)   :: id_denit, id_adsp
      ! for mpbC
      type (type_state_variable_id)        :: id_mpbCHL, id_mpbC, id_mpbN, id_eps
      type (type_dependency_id)            :: id_parz, id_porosity
      type (type_diagnostic_variable_id)   :: id_PrimProd, id_par, id_Q_N, id_Q_chl
      type (type_diagnostic_variable_id)   :: id_expCProd, id_expNProd

!     Model parameters
      !for omexdia_p
      real(rk) :: rFast, rSlow, NCrLdet, NCrSdet, PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox, CprodMax
      real(rk) :: ksO2oduox, ksO2oxic, ksNO3denit, kinO2denit, kinNO3anox, kinO2anox
      ! for mpbC
      logical  :: MPhytoBenOn   ! use MicroPhytoBenthos (MPB)
      real(rk) :: rLdet, rSdet, rNCldet
      real(rk) :: mumax, alpha, gamma, Qmin, Qmax, thetamax, uptmax, KNH4, KNO3
      real(rk) :: KinNH4, keps, resp, Kresp, graz, kout, kexu, rzoo
      real(rk) :: PAR0max, k0, Achla, bTemp
      ! for all
      logical  :: use_dic, use_zbC, use_zbN

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do

   end type type_hzg_omexdia_p_mpb
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
   class (type_hzg_omexdia_p_mpb),intent(inout),target  :: self
   integer,                    intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
!!------- Initial values of model omexdia_p -------
   real(rk)  :: ldetC_init    ! labile detritus carbon (fast decay)
   real(rk)  :: sdetC_init    ! semilabile detritus carbon (slow decay)
   real(rk)  :: detP_init     ! detritus phosphorus
   real(rk)  :: po4_init      ! dissolved phosphate
   real(rk)  :: no3_init      ! dissolved nitrate
   real(rk)  :: nh3_init      ! dissolved ammonium
   real(rk)  :: oxy_init      ! dissolved oxygen
   real(rk)  :: odu_init      ! dissolved reduced substances
!!------- Parameters for model omexdia_p -------
   real(rk)  :: rFast         ! decay rate labile detritus (fast decay)
   real(rk)  :: rSlow         ! decay rate semilabile detritus (slow decay)
   real(rk)  :: NCrLdet       ! NC ratio labile detritus (fast decay)
   real(rk)  :: NCrSdet       ! NC ratio semilabile detritus (slow decay)
   real(rk)  :: PAds          ! Adsorption coeff phosphorus
   real(rk)  :: PAdsODU       ! PO4-Fe dissolution threshold in terms of [FeS]/ODU
   real(rk)  :: NH3Ads        ! Adsorption coeff ammonium
   real(rk)  :: CprodMax=48.0 ! Max C-degrad rate for numeric stability
   real(rk)  :: rnit          ! Max nitrification rate
   real(rk)  :: ksO2nitri     ! half-sat O2 in nitrification
   real(rk)  :: rODUox        ! Max rate oxidation of ODU
   real(rk)  :: ksO2oduox     ! half-sat O2 in oxidation of ODU
   real(rk)  :: ksO2oxic      ! half-sat O2 in oxic mineralis
   real(rk)  :: ksNO3denit    ! half-sat NO3 in denitrif
   real(rk)  :: kinO2denit    ! half-sat O2 inhib denitrif
   real(rk)  :: kinNO3anox    ! half-sat NO3 inhib anoxic min
   real(rk)  :: kinO2anox     ! half-sat O2 inhib anoxic min
!!------- Initial values of model mpb -------
   real(rk)  :: mpbCHL_init   ! MicroPhytoBenthos chlorophyll
   real(rk)  :: mpbC_init     ! MicroPhytoBenthos carbon
   real(rk)  :: mpbN_init     ! MicroPhytoBenthos nitrogen
   real(rk)  :: eps_init      ! Extracellular Polymeric Substances
!!------- Switch for MicroPhytoBenthos model -------
   logical   :: MPhytoBenOn   ! use MicroPhytoBenthos (MPB)
!!------- Parameters for model mpb originating from omexdia_p -------
   real(rk)  :: rLdet         ! decay rate labile detritus (fast decay)
   real(rk)  :: rSdet         ! decay rate semilabile detritus (slow decay)
   real(rk)  :: rNCldet       ! NC ratio labile detritus (fast decay)
!!------- Parameters for model mpb -------
   real(rk)  :: mumax         ! Maximum growth rate
   real(rk)  :: alpha         ! ific initial slope of the PI-curve
   real(rk)  :: gamma         ! Mol O2 produced per mol C fixed by photosynthesis
   real(rk)  :: Qmin          ! Minimum N/C ratio
   real(rk)  :: Qmax          ! Maximum N/C ratio
   real(rk)  :: thetamax      ! Maximum Chla/N ratio
   real(rk)  :: uptmax        ! Maximum uptake rate per carbon unit for NH4 & NO3
   real(rk)  :: KNH4          ! Half-saturation conc. for NH4 uptake
   real(rk)  :: KNO3          ! Half-saturation conc. for NO3 uptake
   real(rk)  :: KinNH4        ! Half-saturation conc. for NO3 uptake inhibition by ammonium
   real(rk)  :: keps          ! fraction of primary production exudated as EPS
   real(rk)  :: resp          ! Respiration rate
   real(rk)  :: Kresp         ! Half-saturation conc. O2 lim. for resp
   real(rk)  :: graz          ! Grazing rate
   real(rk)  :: kout          ! Faeces production coeff.
   real(rk)  :: kexu          ! Nitrogen exudation coeff.
   real(rk)  :: rzoo          ! Respiration rate for zoobenthos
   real(rk)  :: PAR0max       ! maximum light intensity
   real(rk)  :: k0            ! Extinction coefficient of sediment
   real(rk)  :: Achla         ! Absorption factor of chlorophyll
   real(rk)  :: bTemp         ! Temperature increase
!!------- Optional external dependencies -------
   character(len=attribute_length) :: dic_variable  = ''  ! dissolved inorganic carbon
   character(len=attribute_length) :: zbC_variable  = ''  ! zoobenthos carbon
   character(len=attribute_length) :: zbN_variable  = ''  ! zoobenthos nitrogen

   namelist /hzg_omexdia_p/  &
          rFast, rSlow, NCrLdet, NCrSdet, PAds, PAdsODU, NH3Ads, &
          CprodMax, rnit, ksO2nitri, rODUox, ksO2oduox, ksO2oxic, ksNO3denit,      &
          kinO2denit, kinNO3anox, kinO2anox,                                       &
          ldetC_init, sdetC_init, oxy_init, odu_init, no3_init, nh3_init,          &
          detP_init, po4_init

   namelist /hzg_omexdia_p_mpb/ MPhytoBenOn

   namelist /hzg_omexdia_p_mpb_dependencies/  &
          dic_variable, zbC_variable, zbN_variable

   namelist /hzg_mpb/ rLdet, rSdet, rNCldet, &
          mumax, alpha, gamma, Qmin, Qmax, thetamax, uptmax, KNH4, KNO3, KinNH4,   &
          keps, resp, Kresp, graz, kout, kexu, rzoo, PAR0max, k0, Achla, bTemp,    &
          mpbCHL_init, mpbC_init, mpbN_init, eps_init

!EOP
!-----------------------------------------------------------------------
!BOC

   ! Read the namelist
   if (configunit>0) then
     read(configunit, nml=hzg_omexdia_p, err=96, end=100)
     read(configunit, nml=hzg_omexdia_p_mpb, err=97, end=101)
     read(configunit, nml=hzg_omexdia_p_mpb_dependencies, err=98, end=102)
     read(configunit, nml=hzg_mpb, err=99, end=103)
   endif

   ! Store parameter values in our own derived type
   self%rFast      = rFast
   self%rSlow      = rSlow
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
   call self%register_state_variable(self%id_ldetC, 'ldetC', 'mmolC/m**3',  &
                                    'detritus labile carbon', ldetC_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_ldetC, 'particulate', .true.)

   call self%register_state_variable(self%id_sdetC, 'sdetC', 'mmolC/m**3',  &
                                    'detritus semilabile carbon', sdetC_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_sdetC, 'particulate', .true.)

   call self%register_state_variable(self%id_detP, 'detP', 'mmolP/m**3',  &
                                    'detritus phosphorus', detP_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_detP, 'particulate', .true.)

   call self%register_state_variable(self%id_po4, 'po4', 'mmolP/m**3',  &
                                    'dissolved phosphate', po4_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_phosphate)
   call self%set_variable_property(self%id_po4, 'particulate', .false.)

   call self%register_state_variable(self%id_no3, 'no3', 'mmolN/m**3',  &
                                    'dissolved nitrate', no3_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_nitrate)
   call self%set_variable_property(self%id_no3, 'particulate', .false.)

   call self%register_state_variable(self%id_nh3, 'nh3', 'mmolN/m**3',  &
                                    'dissolved ammonium', nh3_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_ammonium)
   call self%set_variable_property(self%id_nh3, 'particulate', .false.)

   call self%register_state_variable(self%id_oxy, 'oxy', 'mmolO2/m**3',  &
                                    'dissolved oxygen', oxy_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_oxy, 'particulate', .false.)

   call self%register_state_variable(self%id_odu,'odu','mmol/m**3',  &
                                    'dissolved reduced substances', odu_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_odu,'particulate', .false.)

   if (MPhytoBenOn) then
      ! Store parameter values in our own derived type
      self%MPhytoBenOn = MPhytoBenOn
      self%rLdet       = rLdet
      self%rSdet       = rSdet
      self%rNCldet     = rNCldet
      self%mumax       = mumax
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
      self%PAR0max     = PAR0max
      self%k0          = k0
      self%Achla       = Achla
      self%bTemp       = bTemp

      ! Register state variables
      call self%register_state_variable(self%id_mpbCHL, 'mpbCHL', 'mmolN/m**3', &
            'MicroPhytoBenthos chlorophyll mpbCHL',                             &
            mpbCHL_init, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_mpbCHL, 'particulate', .true.)

      call self%register_state_variable(self%id_mpbC,  'mpbC', 'mmolC/m**3',    &
            'MicroPhytoBenthos carbon mpbC',                                    &
            mpbC_init, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_mpbC, 'particulate', .true.)

      call self%register_state_variable(self%id_mpbN,  'mpbN', 'mmolN/m**3',    &
            'MicroPhytoBenthos nitrogen mpbN',                                  &
            mpbN_init, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_mpbN, 'particulate', .true.)

      call self%register_state_variable(self%id_eps,   'eps', 'mmolC/m**3',     &
            'Extracellular Polymeric Substances  eps',                          &
            eps_init, minimum=0.0_rk, no_river_dilution=.true.)
      call self%set_variable_property(self%id_eps, 'particulate', .false.)

      ! Register link to external dependencies, if variable names are provided in namelist.
      self%use_dic  = dic_variable/=''
      if (self%use_dic) then
          call self%register_state_dependency(self%id_dic, dic_variable, 'mmol/m**3',   &
                'dissolved inorganic carbon', required=.false.)
          call self%request_coupling(self%id_dic, dic_variable)
      endif
      self%use_zbC  = zbC_variable/=''
      if (self%use_zbC) then
          call self%register_state_dependency(self%id_zbC, zbC_variable, 'mmolC/m**3',   &
                'ZooBenthos carbon', required=.false.)
          call self%request_coupling(self%id_zbC, zbC_variable)
      endif
      self%use_zbN  = zbN_variable/=''
      if (self%use_zbN) then
          call self%register_state_dependency(self%id_zbN, zbN_variable, 'mmolN/m**3',   &
                'ZooBenthos nitrogen', required=.false.)
          call self%request_coupling(self%id_zbN, zbN_variable)
      endif

   endif

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_adsp, 'adsP', 'mmolP/m**3',  &
         'phosphate adsorption', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_denit, 'denit', 'mmol/m**3/d',  &
         'denitrification rate', output=output_instantaneous)

   self%MPhytoBenOn = MPhytoBenOn
   if (self%MPhytoBenOn) then
      call self%register_diagnostic_variable(self%id_PrimProd, 'PrimProd', 'mmolC/m**3/d',  &
            'MPB primary production rate PrimProd',    output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_par,      'par', 'w/m2',               &
            'photosynthetically active radiation par', output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_Q_N,      'Q_N', '-',                  &
            'MPB nitrogen quota Q_N',                  output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_Q_chl,    'Q_chl', '-',                &
            'MPB CHL:C ratio Q_chl',                   output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_expCProd, 'expCProd', 'mmolC/m**2/d',  &
            'MPB carbon export (zoobenthos grazing)',    output=output_instantaneous)
      call self%register_diagnostic_variable(self%id_expNProd, 'expNProd', 'mmolN/m**2/d',  &
            'MPB nitrogen export (zoobenthos grazing)',  output=output_instantaneous)
   endif

   ! Register dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_parz, standard_variables%downwelling_photosynthetic_radiative_flux)

   return

96 call self%fatal_error('hzg_omexdia_p_mpb_initialize','Error reading namelist hzg_omexdia_p.')
97 call self%fatal_error('hzg_omexdia_p_mpb_initialize','Error reading namelist hzg_omexdia_p_mpb.')
98 call self%fatal_error('hzg_omexdia_p_mpb_initialize','Error reading namelist hzg_omexdia_p_mpb_dependencies.')
99 call self%fatal_error('hzg_omexdia_p_mpb_initialize','Error reading namelist hzg_mpb.')

100 call self%fatal_error('hzg_omexdia_p_mpb_initialize','Namelist hzg_omexdia_p was not found.')
101 call self%fatal_error('hzg_omexdia_p_mpb_initialize','Namelist hzg_omexdia_p_mpb was not found.')
102 call self%fatal_error('hzg_omexdia_p_mpb_initialize','Namelist hzg_omexdia_p_mpb_dependencies was not found.')
103 call self%fatal_error('hzg_omexdia_p_mpb_initialize','Namelist hzg_mpb was not found.')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of OMEXDIA+P model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_hzg_omexdia_p_mpb),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
!!------- for model omexdia_p -------
   real(rk),parameter :: relaxO2=0.04_rk
   real(rk),parameter :: T0 = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk),parameter :: Q10b = 1.5_rk
   real(rk) :: ldetC, sdetC, oxy, odu, no3, nh3, detP, po4
   real(rk) :: temp_celsius, temp_kelvin, f_T, E_a
   real(rk) :: radsP, Oxicminlim, Denitrilim, Anoxiclim, Rescale, rP

   real(rk) :: CprodF, CprodS, Cprod, Nprod, Pprod
   real(rk) :: AnoxicMin, Denitrific, OxicMin, Nitri, OduDepo, OduOx, pDepo
!!------- for model mpb -------
   real(rk) :: mpbC, mpbN, mpbCHL, eps
   real(rk) :: par, parz, porosity, CprodEPS
   real(rk) :: prodChl, k, theta, Q_N, Q_chl, Pmax,  PP, prod, prodeps, fac
   real(rk) :: prodO2, rhochl, uptNH4, uptNO3, uptchl, uptN, respphyto, faecesC, faecesN
   real(rk) :: exud, grazingC, grazingN, grazingChl, respzoo, exportC, exportN

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp, temp_celsius) ! sediment-water temperature
   _GET_(self%id_ldetC, ldetC)       ! labile detritus C in mmolC/m**3
   _GET_(self%id_sdetC, sdetC)       ! semilabile detritus C in mmolC/m**3
   _GET_(self%id_detP, detP)         ! detritus P in mmolP/m**3
   _GET_(self%id_oxy, oxy)           ! dissolved oxygen in mmolO2/m**3
   _GET_(self%id_odu, odu)           ! dissolved reduced substances in mmolO2/m**3
   _GET_(self%id_no3, no3)           ! dissolved nitrate in mmolN/m**3
   _GET_(self%id_nh3, nh3)           ! dissolved ammonium in mmolN/m**3
   _GET_(self%id_po4, po4)           ! dissolved phosphate in mmolP/m**3
   if (self%MPhytoBenOn) then
      _GET_(self%id_parz,   parz)    ! sediment light normalized to one at top
      _GET_(self%id_mpbCHL, mpbCHL)  ! MicroPhytoBenthos chlorophyll in mmolN/m**3
      _GET_(self%id_mpbC,   mpbC)    ! MicroPhytoBenthos carbon in mmolC/m**3
      _GET_(self%id_mpbN,   mpbN)    ! MicroPhytoBenthos nitrogen in mmolN/m**3
      _GET_(self%id_eps,    eps)     ! Extracellular Polymeric Substances  in mmolC/m**3
   end if

   temp_kelvin = 273.15_rk + temp_celsius
   E_a = 0.1_rk*log(Q10b)*T0*(T0+10.0_rk)
   f_T = 1.0_rk*exp(-E_a*(1.0_rk/temp_kelvin - 1.0_rk/T0))

   if (2*oxy < -self%kinO2anox) oxy = -self%kinO2anox/2

   Oxicminlim = oxy/(oxy+self%ksO2oxic+relaxO2*(nh3+odu))                ! limitation terms

   Denitrilim = (1.0_rk-oxy/(oxy+self%kinO2denit)) * NO3/(no3+self%ksNO3denit)
   Anoxiclim  = (1.0_rk-oxy/(oxy+self%kinO2anox)) * (1.0_rk-no3/(no3+self%kinNO3anox))
   if(Oxicminlim - 1E-3 < -Denitrilim-Anoxiclim) Oxicminlim = 1E-3-Denitrilim-Anoxiclim
   Rescale    = 1.0_rk/(Oxicminlim+Denitrilim+Anoxiclim)

   CprodF = f_T * self%rFast * ldetC
   CprodS = f_T * self%rSlow * sdetC

   ! assume upper reactive surface area for POC hydrolysis
   if (CprodS > self%CprodMax) CprodS = self%CprodMax
   if (CprodF > self%CprodMax) CprodF = self%CprodMax

   Cprod  = CprodF + CprodS
   Nprod  = CprodF * self%NCrLdet + CprodS * self%NCrSdet


   ! PO4-adsorption ceases when critical capacity is reached
   ! [FeS] approximated by ODU
   ! TODO: temperature dependency
   radsP  = self%PAds  * po4 * 1.0_rk/(1.0_rk+exp(-5.0_rk+(odu-oxy)/self%PAdsODU))
   rP     = f_T * self%rFast * 2 ! (1.0_rk - Oxicminlim)
   Pprod  = rP * detP

   ! Oxic mineralisation, denitrification, anoxic mineralisation
   ! then the mineralisation rates
   OxicMin    = Cprod*Oxicminlim*Rescale        ! oxic mineralisation
   Denitrific = Cprod*Denitrilim*Rescale        ! Denitrification
   AnoxicMin  = Cprod*Anoxiclim *Rescale        ! anoxic mineralisation

   ! reoxidation and ODU deposition
   Nitri      = f_T * self%rnit   * nh3 * oxy/(oxy + self%ksO2nitri + relaxO2*(ldetC + odu))
   OduOx      = f_T * self%rODUox * odu * oxy/(oxy + self%ksO2oduox + relaxO2*(nh3 + ldetC))
   if (OduOx > self%CprodMax) OduOx = self%CprodMax

   !  pDepo      = min(1.0_rk,0.233_rk*(wDepo)**0.336_rk )
   pDepo      = 0.0_rk
   OduDepo    = AnoxicMin*pDepo

   if (self%MPhytoBenOn) then
      par    = self%PAR0max * parz
      !print*,'par#386: ',par, parz
      !k      = self%k0 + self%Achla * mpbCHL

      !----- intracellular N:C:Chl stoichiometry
      Q_N    = mpbN   / mpbC
      Q_chl  = mpbCHL / mpbC
      theta  = mpbCHL / mpbN

      !----- photosnythesis rate
      fac    = 1.0 - self%Qmin/Q_N
      if (fac .lt. 0.0) fac = 0.0
      Pmax   = 2*self%mumax * exp(self%bTemp * temp_celsius) * fac
      PP     = Pmax * (1.0 - exp(-self%alpha * parz*Q_chl/Pmax))
      prod   = PP * mpbC
      prodO2 = self%gamma * prod

      if (parz .gt. 0.0) then
        rhochl = self%thetamax * PP/(self%alpha * parz * Q_chl)
      else
        rhochl = 0.0
      endif
      !write(*,'(''#407'',6e20.10)') rhochl, self%thetamax, PP, self%alpha, parz, Q_chl

      !----- nutrient uptake (TODO add dependencies on P, Si)
      uptNO3 = self%uptmax * no3/(no3 + self%KNO3)
      uptNH4 = self%uptmax * nh3/(nh3 + self%KNH4)
      uptN   = uptNO3 + uptNH4
      !----- chlorophyll synthesis
      fac    = 1.0 - theta/self%thetamax
      if (fac .lt. 0.0) fac = 0.0
      uptchl = uptN * fac/(fac + 0.05)
      prodChl= rhochl * uptchl * mpbC
      ! write(*,'(''#525'',6e20.10)') prodChl, rhochl, uptchl, mpbC

      fac    = (self%Qmax - Q_N) / (self%Qmax - self%Qmin)
      if (fac .lt. 0.0) fac = 0.0
      uptNO3 = uptNO3 * fac * (1.0 - nh3/(nh3 + self%KinNH4)) * mpbC
      uptNH4 = uptNH4 * fac * mpbC
      uptN   = uptNO3 + uptNH4

      ! Carbohydrate exudation:
      prodeps = self%keps * prod
      CprodEPS = sqrt(self%rLdet*self%rSdet) * eps

      ! Respiration:
      respphyto = self%resp * mpbC * oxy/(oxy+self%Kresp)

      ! Zoobenthos grazing and associated processes:
      grazingC   = self%graz * mpbC
      grazingN   = self%graz * mpbN
      grazingChl = self%graz * mpbCHL
      faecesC    = self%kout * grazingC
      faecesN    = self%kout * grazingC * self%rNCldet ! Hochard et al 2010
      !faecesN    = self%kout * grazingN ! alternative proposed by Markus Kreus
      exud       = self%kexu * grazingN
      respzoo    = self%rzoo * grazingC * oxy/(oxy+self%Kresp)
      exportC    = grazingC - faecesC - respzoo ! exported carbon   (open closure term)
      exportN    = grazingN - faecesN - exud    ! exported nitrogen (open closure term)

   else
      faecesC   = 0.0_rk
      faecesN   = 0.0_rk
      respzoo   = 0.0_rk
      respphyto = 0.0_rk
      exud      = 0.0_rk
      uptNH4    = 0.0_rk
      uptNO3    = 0.0_rk
      CprodEPS  = 0.0_rk
   end if

#define _CONV_UNIT_ *one_pr_day
   ! reaction rates
   _SET_ODE_(self%id_ldetC, (-CprodF + faecesC) _CONV_UNIT_)
   _SET_ODE_(self%id_sdetC, (-CprodS) _CONV_UNIT_)
   _SET_ODE_(self%id_oxy , (-OxicMin - 2.0_rk* Nitri - OduOx - respzoo - respphyto) _CONV_UNIT_) !RH 1.0->150/106*OxicMin (if [oxy]=mmolO2/m**3)
   _SET_ODE_(self%id_no3 , (-0.8_rk*Denitrific + Nitri - uptNO3) _CONV_UNIT_)     !RH 0.8-> ~104/106?
!   _SET_ODE_(self%id_nh3 , (Nprod - Nitri + exud - uptNH4) / (1.0_rk + self%NH3Ads) _CONV_UNIT_)
   _SET_ODE_(self%id_nh3 , (Nprod - Nitri + exud - uptNH4) _CONV_UNIT_)
   _SET_ODE_(self%id_odu , (AnoxicMin - OduOx - OduDepo) _CONV_UNIT_)
   _SET_ODE_(self%id_po4 , (Pprod - radsP) _CONV_UNIT_)
   _SET_ODE_(self%id_detP, (radsP - Pprod - self%NH3Ads*detP) _CONV_UNIT_)
   if (self%MPhytoBenOn) then
      _SET_ODE_(self%id_mpbCHL, ( prodchl - grazingChl ) _CONV_UNIT_)
      _SET_ODE_(self%id_mpbC,   ( prod - grazingC - respphyto - prodeps ) _CONV_UNIT_)
      _SET_ODE_(self%id_mpbN,   ( uptN - grazingN ) _CONV_UNIT_)
      _SET_ODE_(self%id_eps,    ( prodeps - f_T * CprodEPS ) _CONV_UNIT_)
      if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic , (respzoo + respphyto) _CONV_UNIT_)
      if (_AVAILABLE_(self%id_zbC)) _SET_ODE_(self%id_zbC , (exportC) _CONV_UNIT_)
      if (_AVAILABLE_(self%id_zbN)) _SET_ODE_(self%id_zbN , (exportN) _CONV_UNIT_)
   end if

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit,0.8_rk*Denitrific)    !last denitrification rate
   _SET_DIAGNOSTIC_(self%id_adsp ,radsP)                !instantaneous phosphate adsorption
   if (self%MPhytoBenOn) then
      _SET_DIAGNOSTIC_(self%id_PrimProd, PP)            !instantaneous MPB primary production rate
      _SET_DIAGNOSTIC_(self%id_par, par)                !instantaneous photosynthetically active radiation
      _SET_DIAGNOSTIC_(self%id_Q_N, Q_N)                !instantaneous MPB N:C quota
      _SET_DIAGNOSTIC_(self%id_Q_chl, Q_chl)            !instantaneous MPB CHL:C ratio
      _SET_DIAGNOSTIC_(self%id_expCProd, exportC)       !instantaneous MPB export production carbon
      _SET_DIAGNOSTIC_(self%id_expNProd, exportN)       !instantaneous MPB export production nitrogen
   endif

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
  end module hzg_omexdia_p_mpb
!-----------------------------------------------------------------------

