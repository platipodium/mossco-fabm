#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_omexdia_p --- Fortran 2003 version of OMEXDIA+P biogeochemical model
!
! !INTERFACE:
   module hzg_omexdia_p
!
! !DESCRIPTION:
!
! The OMEXDIA+P model is based on the OMEXDIA model (see Soetard et al. 1996a)
! and is intended to simulate early diagenesis in the sea sediments. The major
! difference to the original OMEXDIA is an added phosphorus cycle.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_hzg_omexdia_p
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
   type,extends(type_base_model) :: type_hzg_omexdia_p
!     Variable identifiers
      type (type_state_variable_id)        :: id_ldetC, id_sdetC, id_detP
      type (type_state_variable_id)        :: id_no3, id_nh3, id_oxy, id_po4, id_odu
      type (type_dependency_id)            :: id_temp
      type (type_diagnostic_variable_id)   :: id_denit, id_adsp

!     Model parameters
      real(rk) :: rFast, rSlow, NCrLdet, NCrSdet, PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox, CprodMax
      real(rk) :: ksO2oduox, ksO2oxic, ksNO3denit, kinO2denit, kinNO3anox, kinO2anox

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do

   end type type_hzg_omexdia_p
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the OMEXDIA+P model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the omexdia namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hzg_omexdia_p), intent(inout),target  :: self
   integer,                    intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
!!------- Initial values -------
   real(rk)  :: ldetC_init    ! labile detritus carbon (fast decay)
   real(rk)  :: sdetC_init    ! semilabile detritus carbon (slow decay)
   real(rk)  :: detP_init     ! detritus phosphorus
   real(rk)  :: po4_init      ! dissolved phosphate
   real(rk)  :: no3_init      ! dissolved nitrate
   real(rk)  :: nh3_init      ! dissolved ammonium
   real(rk)  :: oxy_init      ! dissolved oxygen
   real(rk)  :: odu_init      ! dissolved reduced substances
!!------- Parameters -------
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

   namelist /hzg_omexdia_p/  &
          rFast, rSlow, NCrLdet, NCrSdet, PAds, PAdsODU, NH3Ads, &
          CprodMax, rnit, ksO2nitri, rODUox, ksO2oduox, ksO2oxic, ksNO3denit,      &
          kinO2denit, kinNO3anox, kinO2anox,                                       &
          ldetC_init, sdetC_init, oxy_init, odu_init, no3_init, nh3_init,          &
          detP_init, po4_init
!EOP
!-----------------------------------------------------------------------
!BOC

   ! Read the namelist
   if (configunit>0) read(configunit,nml=hzg_omexdia_p,err=99,end=100)

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

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_adsp, 'adsP', 'mmolP/m**3',  &
         'phosphate adsorption', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_denit, 'denit', 'mmol/m**3/d',  &
         'denitrification rate', output=output_instantaneous)

   ! Register dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)

   return

99 call self%fatal_error('hzg_omexdia_p_initialize', 'Error reading namelist hzg_omexdia_p.')

100 call self%fatal_error('hzg_omexdia_p_initialize', 'Namelist hzg_omexdia_p was not found.')

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
   class (type_hzg_omexdia_p),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
   real(rk), parameter :: relaxO2=0.04_rk
   real(rk), parameter :: T0 = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk), parameter :: Q10b = 1.5_rk
   real(rk) :: ldetC, sdetC, oxy, odu, no3, nh3, detP, po4
   real(rk) :: temp_celsius, temp_kelvin, f_T, E_a
   real(rk) :: radsP, Oxicminlim, Denitrilim, Anoxiclim, Rescale, rP
   real(rk) :: CprodF, CprodS, Cprod, Nprod, Pprod
   real(rk) :: AnoxicMin, Denitrific, OxicMin, Nitri, OduDepo, OduOx, pDepo
!   real(rk) :: CprodMax = 24.0 ! d^-1 maximal C-degrad rate for numeric stability
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

#define _CONV_UNIT_ *one_pr_day
! reaction rates
   _SET_ODE_(self%id_ldetC, -CprodF _CONV_UNIT_)
   _SET_ODE_(self%id_sdetC, -CprodS _CONV_UNIT_)
   _SET_ODE_(self%id_oxy,  (-OxicMin - 2.0_rk* Nitri - OduOx) _CONV_UNIT_) !RH 1.0->150/106*OxicMin (if [oxy]=mmolO2/m**3)
   _SET_ODE_(self%id_no3,  (-0.8_rk*Denitrific + Nitri) _CONV_UNIT_)       !RH 0.8-> ~104/106?
!   _SET_ODE_(self%id_nh3,  (Nprod - Nitri) / (1.0_rk + self%NH3Ads) _CONV_UNIT_)
   _SET_ODE_(self%id_nh3,  (Nprod - Nitri) _CONV_UNIT_)
   _SET_ODE_(self%id_odu,  (AnoxicMin - OduOx - OduDepo) _CONV_UNIT_)
   _SET_ODE_(self%id_po4,  (Pprod - radsP) _CONV_UNIT_)
   _SET_ODE_(self%id_detP, (radsP - Pprod - self%NH3Ads*detP) _CONV_UNIT_)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit, 0.8_rk*Denitrific)
   _SET_DIAGNOSTIC_(self%id_adsp, radsP)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
   end module hzg_omexdia_p
!-----------------------------------------------------------------------

