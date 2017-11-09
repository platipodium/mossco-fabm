#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_omexdia_cnp --- Fortran 2003 version of OMEXDIA+CNP biogeochemical model
!
! !INTERFACE:
   module hzg_omexdia_cnp
!
! !DESCRIPTION:
!
! The OMEXDIA+CNP model is based on the OMEXDIA model (see Soetard et al. 1996a)
! and is intended to simulate early diagenesis in the sea sediments.
! The major difference to the original OMEXDIA model is an added phosphorus cycle
! added by Kai Wirtz. Further modifications comprise efficient reaction and
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
   public type_hzg_omexdia_cnp
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
   type,extends(type_base_model) :: type_hzg_omexdia_cnp
!     Variable identifiers
      type (type_state_variable_id)        :: id_ldetC, id_sdetC, id_ldetN, id_sdetN, id_ldetP
      type (type_state_variable_id)        :: id_no3, id_nh3, id_oxy, id_po4, id_odu
      type (type_dependency_id)            :: id_temp
      type (type_diagnostic_variable_id)   :: id_denit, id_adsp

!     Model parameters
      real(rk) :: rLabile, rSemilabile, NCrLdet, NCrSdet, PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox, CprodMax
      real(rk) :: ksO2oduox, ksO2oxic, ksNO3denit, kinO2denit, kinNO3anox, kinO2anox

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do

   end type type_hzg_omexdia_cnp
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
   class (type_hzg_omexdia_cnp), intent(inout),target  :: self
   integer,                    intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
!!------- Initial values -------
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
!!------- Parameters -------  ! (reference Values refer to Soetaert (1996))
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

   namelist /hzg_omexdia_cnp/  &
          rLabile, rSemilabile, NCrLdet, NCrSdet, PAds, PAdsODU, NH3Ads,       &
          CprodMax, rnit, ksO2nitri, rODUox, ksO2oduox, ksO2oxic, ksNO3denit,  &
          kinO2denit, kinNO3anox, kinO2anox,                                   &
          ldetC_init, sdetC_init, ldetN_init, sdetN_init, ldetP_init,          &
          no3_init, nh3_init, po4_init, oxy_init, odu_init
!EOP
!-----------------------------------------------------------------------
!BOC

   ! Read the namelist
   if (configunit>0) read(configunit,nml=hzg_omexdia_cnp,err=99,end=100)

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
   call self%register_state_variable(self%id_ldetC, 'ldetC', 'mmolC/m**3',  &
                                    'detritus labile carbon', ldetC_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_ldetC, 'particulate', .true.)

   call self%register_state_variable(self%id_sdetC, 'sdetC', 'mmolC/m**3',  &
                                    'detritus semilabile carbon', sdetC_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_sdetC, 'particulate', .true.)

   call self%register_state_variable(self%id_ldetN, 'ldetN', 'mmolN/m**3',  &
                                    'detritus labile nitrogen', ldetN_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_ldetN, 'particulate', .true.)

   call self%register_state_variable(self%id_sdetN, 'sdetN', 'mmolN/m**3',  &
                                    'detritus semilabile nitrogen', sdetN_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_sdetN, 'particulate', .true.)

   call self%register_state_variable(self%id_ldetP, 'ldetP', 'mmolP/m**3',  &
                                    'detritus labile phosphorus', ldetP_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_ldetP, 'particulate', .true.)

   call self%register_state_variable(self%id_no3, 'no3', 'mmolN/m**3',  &
                                    'dissolved nitrate', no3_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_nitrate)
   call self%set_variable_property(self%id_no3, 'particulate', .false.)

   call self%register_state_variable(self%id_nh3, 'nh3', 'mmolN/m**3',  &
                                    'dissolved ammonium', nh3_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_ammonium)
   call self%set_variable_property(self%id_nh3, 'particulate', .false.)

   call self%register_state_variable(self%id_po4, 'po4', 'mmolP/m**3',  &
                                    'dissolved phosphate', po4_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_phosphate)
   call self%set_variable_property(self%id_po4, 'particulate', .false.)

   call self%register_state_variable(self%id_oxy, 'oxy', 'mmolO2/m**3',  &
                                    'dissolved oxygen', oxy_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_oxy, 'particulate', .false.)

   call self%register_state_variable(self%id_odu,'odu','mmol/m**3',  &
                                    'dissolved reduced substances', odu_init, minimum=0.0_rk)
   call self%set_variable_property(self%id_odu,'particulate', .false.)

   ! Register diagnostic variables
   !call self%register_diagnostic_variable(self%id_adsp, 'adsP', 'mmolP/m**3',  &
   !      'phosphate adsorption', output=output_instantaneous)
   !call self%register_diagnostic_variable(self%id_denit, 'denit', 'mmolN/m**3/d',  &
   !      'denitrification rate', output=output_instantaneous)

   ! Register dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)

   return

99 call self%fatal_error('hzg_omexdia_cnp_initialize', 'Error reading namelist hzg_omexdia_cnp.')

100 call self%fatal_error('hzg_omexdia_cnp_initialize', 'Namelist hzg_omexdia_cnp was not found.')

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
   class (type_hzg_omexdia_cnp),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
   real(rk), parameter :: gammaNO3 = 0.8_rk    ! molN  used per molC in denitrification     !RH 0.8-> ~104/106?
   real(rk), parameter :: gammaO2  = 1.0_rk    ! molO2 used per molC in oxic mineralization !RH 1.0->150/106*OxicMin (if [oxy]=mmolO2/m**3)
   real(rk), parameter :: gammaNH3 = 2.0_rk    ! molO2 needed to oxidize one molNH3 in nitrification
   real(rk), parameter :: T0       = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk), parameter :: Q10b     = 1.5_rk    ! q10 temperature coefficient (source?)
   real(rk), parameter :: relaxO2  = 0.04_rk   ! (source?)
   real(rk) :: ldetC, sdetC, ldetN, sdetN, ldetP, no3, nh3, po4, oxy, odu
   real(rk) :: temp_celsius, temp_kelvin, f_T, E_a
   real(rk) :: radsP, Oxicminlim, Denitrilim, Anoxiclim, Rescale, rP
   real(rk) :: CprodL, CprodS, Cprod, NprodL, NprodS, Nprod, Pprod
   real(rk) :: AnoxicMin, Denitrific, OxicMin, Nitri, OduDepo, OduOx, pDepo
!   real(rk) :: CprodMax = 24.0 ! d^-1 maximal C-degrad rate for numeric stability
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp, temp_celsius) ! sediment-water temperature (degC)
   _GET_(self%id_ldetC, ldetC)       ! labile detritus C          (mmolC/m**3)
   _GET_(self%id_sdetC, sdetC)       ! semilabile detritus C      (mmolC/m**3)
   _GET_(self%id_ldetN, ldetN)       ! labile detritus N          (mmolN/m**3)
   _GET_(self%id_sdetN, sdetN)       ! semilabile detritus N      (mmolN/m**3)
   _GET_(self%id_ldetP, ldetP)       ! labile detritus P          (mmolP/m**3)
   _GET_(self%id_oxy, oxy)           ! dissolved oxygen           (mmolO2/m**3)
   _GET_(self%id_odu, odu)           ! dissolved oxygen demand units (mmolO2/m**3)
   _GET_(self%id_no3, no3)           ! dissolved nitrate          (mmolN/m**3)
   _GET_(self%id_nh3, nh3)           ! dissolved ammonium         (mmolN/m**3)
   _GET_(self%id_po4, po4)           ! dissolved phosphate        (mmolP/m**3)

   ! temperature dependency
   temp_kelvin = 273.15_rk + temp_celsius
   E_a = 0.1_rk*log(Q10b)*T0*(T0+10.0_rk)
   f_T = 1.0_rk*exp(-E_a*(1.0_rk/temp_kelvin - 1.0_rk/T0)) ! temperature factor (-)

   if (2*oxy < -self%kinO2anox) oxy = -self%kinO2anox/2    ! limit anoxic conditions

   ! limitation terms
   Oxicminlim = oxy/(oxy+self%ksO2oxic+relaxO2*(nh3+odu))                                ! (-)
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
   rP     = f_T * self%rLabile * 2 ! (1.0_rk - Oxicminlim)  ! (d-1)
   Pprod  = rP * ldetP                                      ! (mmolP m-3 d-1)

! Oxic mineralisation, denitrification, anoxic mineralisation
! then the mineralisation rates
   OxicMin    = Cprod*Oxicminlim*Rescale  ! carbon oxidized by O2              (mmolC m-3 d-1) !Soetaert eq 3.1
   Denitrific = Cprod*Denitrilim*Rescale  ! carbon oxidized by denitrification (mmolC m-3 d-1) !Soetaert eq 3.2
   AnoxicMin  = Cprod*Anoxiclim *Rescale  ! carbon oxidized by other oxidants  (mmolC m-3 d-1) !Soetaert eq 3.3

! Ammonium (NH3) nitrified (mmolN m-3 d-1) !Soetaert eq 3.7
   Nitri      = f_T * self%rnit   * nh3 * oxy/(oxy + self%ksO2nitri + relaxO2*(ldetC + odu)) !Soetaert eq 3.7, but where does the relax* term comes from?

! Oxygen consumed in the oxidation of reduced substances (mmolO2 m-3 d-1) !Soetaert eq 3.6
   OduOx      = f_T * self%rODUox * odu * oxy/(oxy + self%ksO2oduox + relaxO2*(nh3 + ldetC)) !Soetaert eq 3.6, but where does the relax* term comes from?
   if (OduOx > self%CprodMax) OduOx = self%CprodMax

! ODU deposition
!  pDepo      = min(1.0_rk,0.233_rk*(wDepo)**0.336_rk )
   pDepo      = 0.0_rk
   OduDepo    = AnoxicMin*pDepo           ! ODU deposited as solids (in O2-equivalents) (mmolO2 m-3 d-1) !Soetaert eq 3.5

#define _CONV_UNIT_ *one_pr_day
! reaction rates
   _SET_ODE_(self%id_ldetC, -CprodL                                        _CONV_UNIT_)  ! (mmolC  m-3 d-1)
   _SET_ODE_(self%id_sdetC, -CprodS                                        _CONV_UNIT_)  ! (mmolC  m-3 d-1)
   _SET_ODE_(self%id_ldetN, -NprodL                                        _CONV_UNIT_)  ! (mmolN  m-3 d-1)
   _SET_ODE_(self%id_sdetN, -NprodS                                        _CONV_UNIT_)  ! (mmolN  m-3 d-1)
   _SET_ODE_(self%id_ldetP, (radsP - Pprod - self%NH3Ads*ldetP)            _CONV_UNIT_)  ! (mmolP  m-3 d-1)
   _SET_ODE_(self%id_no3,  (-Denitrific*gammaNO3 + Nitri)                  _CONV_UNIT_)  ! (mmolN  m-3 d-1)
!   _SET_ODE_(self%id_nh3,  (Nprod - Nitri) / (1.0_rk + self%NH3Ads)       _CONV_UNIT_)  ! (mmolN  m-3 d-1)
   _SET_ODE_(self%id_nh3,  (Nprod - Nitri)                                 _CONV_UNIT_)  ! (mmolN  m-3 d-1)
   _SET_ODE_(self%id_po4,  (Pprod - radsP)                                 _CONV_UNIT_)  ! (mmolP  m-3 d-1)
   _SET_ODE_(self%id_oxy,  (-OxicMin*gammaO2     - Nitri*gammaNH3 - OduOx) _CONV_UNIT_)  ! (mmolO2 m-3 d-1)
   _SET_ODE_(self%id_odu,  (AnoxicMin - OduOx - OduDepo)                   _CONV_UNIT_)  ! (mmolO2 m-3 d-1)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit, Denitrific*gammaNO3)  !last denitrification rate
   _SET_DIAGNOSTIC_(self%id_adsp, radsP)                 !instantaneous phosphate adsorption

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
   end module hzg_omexdia_cnp
!-----------------------------------------------------------------------

