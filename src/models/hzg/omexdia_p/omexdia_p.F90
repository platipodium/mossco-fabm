#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_omexdia_p --- Fortran 2003 version of OMEXDIA+P biogeochemical model
!
! !INTERFACE:
   module hzg_omexdia_p
! #define debugMK
!
! !DESCRIPTION:
!
! The OMEXDIA+P model is based on the OMEXDIA model (see Soetard et al. 1996a)
! and is intended to simulate early diagenesis in the sea sediments.
! The major difference to the original OMEXDIA model is an added phosphorus cycle
! added by Kai Wirtz. Further modifications comprise efficient reaction and
! limitation terms to facilitate simple numerics (kai wirtz).
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
!  Modified by author(s): Carsten Lemmen
!  Original author(s): Richard Hofmeister, Kai W. Wirtz
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
      real(rk) :: rLabile, rSemilabile, NCrLdet, NCrSdet, PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox, CprodMax
      real(rk) :: ksO2oduox, ksO2oxic, ksNO3denit, kinO2denit, kinNO3anox, kinO2anox

!     Variable identifiers for parameters (if requested), @todo rate constants
!     are not yet considered here (ksXXX and kinXXX)
      type (type_horizontal_dependency_id)  :: id_rLabile, id_rSemilabile, id_NCrLdet
      type (type_horizontal_dependency_id)  :: id_NCrSdet, id_PAds, id_PAdsODU, id_CprodMax
      type (type_horizontal_dependency_id)  :: id_NH3Ads, id_rnit, id_rODUox

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
!  Modified by author(s): Carsten Lemmen
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
!!------- Initial values -------
   real(rk)  :: ldetC_init    ! labile detritus carbon (fast decay)         (        mmolC  m-3 solid)
   real(rk)  :: sdetC_init    ! semilabile detritus carbon (slow decay)     (        mmolC  m-3 solid)
   real(rk)  :: detP_init     ! detritus phosphorus                         (        mmolP  m-3 liquid)
   real(rk)  :: po4_init      ! dissolved phosphate                         (        mmol   m-3 liquid)
   real(rk)  :: no3_init      ! dissolved nitrate                           (        mmolN  m-3 liquid)
   real(rk)  :: nh3_init      ! dissolved ammonium                          (        mmolN  m-3 liquid)
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

   namelist /hzg_omexdia_p/  &
          rLabile, rSemilabile, NCrLdet, NCrSdet, PAds, PAdsODU, NH3Ads, &
          CprodMax, rnit, ksO2nitri, rODUox, ksO2oduox, ksO2oxic, ksNO3denit,      &
          kinO2denit, kinNO3anox, kinO2anox,                                       &
          ldetC_init, sdetC_init, oxy_init, odu_init, no3_init, nh3_init,          &
          detP_init, po4_init
!EOP
!-----------------------------------------------------------------------
!BOC

   ! Read the namelist
   if (configunit>0) read(configunit, nml=hzg_omexdia_p, err=99, end=100)

   ! Store some parameter values in our own derived type
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

   ! Store other parameter values in our own derived type, if they are negative,
   ! then they are also registered as an external dependency (see below)
   self%rLabile    = rLabile
   self%rSemilabile= rSemilabile
   self%NCrLdet    = NCrLdet
   self%NCrSdet    = NCrSdet
   self%PAds       = PAds
   self%PAdsODU    = PAdsODU
   self%NH3Ads     = NH3Ads

   ! Register state variables
   call self%register_state_variable(self%id_ldetC, 'ldetC', 'mmolC/m**3',  &
                                    'detritus labile carbon', ldetC_init, &
                                    minimum=0.0_rk,missing_value=-1.e30_rk)
   call self%set_variable_property(self%id_ldetC, 'particulate', .true.)

   call self%register_state_variable(self%id_sdetC, 'sdetC', 'mmolC/m**3',  &
                                    'detritus semilabile carbon', sdetC_init, &
                                    minimum=0.0_rk,missing_value=-1.e30_rk)
   call self%set_variable_property(self%id_sdetC, 'particulate', .true.)

   call self%register_state_variable(self%id_detP, 'detP', 'mmolP/m**3',  &
                                    'detritus phosphorus', detP_init, &
                                    minimum=0.0_rk,missing_value=-1.e30_rk)
   call self%set_variable_property(self%id_detP, 'particulate', .true.)

   call self%register_state_variable(self%id_po4, 'po4', 'mmolP/m**3',  &
                                    'dissolved phosphate', po4_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_phosphate, &
                                    missing_value=-1.e30_rk)
   call self%set_variable_property(self%id_po4, 'particulate', .false.)

   call self%register_state_variable(self%id_no3, 'no3', 'mmolN/m**3',  &
                                    'dissolved nitrate', no3_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_nitrate, &
                                    missing_value=-1.e30_rk)
   call self%set_variable_property(self%id_no3, 'particulate', .false.)

   call self%register_state_variable(self%id_nh3, 'nh3', 'mmolN/m**3',  &
                                    'dissolved ammonium', nh3_init, minimum=0.0_rk,  &
                                    standard_variable=standard_variables%mole_concentration_of_ammonium, &
                                    missing_value=-1.e30_rk)
   call self%set_variable_property(self%id_nh3, 'particulate', .false.)

   call self%register_state_variable(self%id_oxy, 'oxy', 'mmolO2/m**3',  &
                                    'dissolved oxygen', oxy_init, &
                                    minimum=0.0_rk, missing_value=-1.e30_rk)
   call self%set_variable_property(self%id_oxy, 'particulate', .false.)

   call self%register_state_variable(self%id_odu,'odu','mmol/m**3',  &
                                    'dissolved reduced substances', odu_init, &
                                    minimum=0.0_rk, missing_value=-1.e30_rk)
   call self%set_variable_property(self%id_odu,'particulate', .false.)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_adsp, 'adsP', 'mmolP/m**3',  &
         'phosphate adsorption', output=output_instantaneous, missing_value=-1.e30_rk)
   call self%register_diagnostic_variable(self%id_denit, 'denit', 'mmolN/m**3/d',  &
         'denitrification rate', output=output_instantaneous, missing_value=-1.e30_rk)

   ! Register dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)

   ! Register dependencies for those parameters that are set to negative values
   if (rLabile < 0) call self%register_horizontal_dependency(self%id_rLabile, 'labile_detritus_decay_rate','d-1','labile_detritus_decay_rate')
   if (rSemiLabile < 0) call self%register_horizontal_dependency(self%id_rSemiLabile, 'semilabile_detritus_decay_rate','d-1','labile_detritus_decay_rate')
   if (NCrLdet < 0) call self%register_horizontal_dependency(self%id_NCrLdet, 'labile_detritus_n_to_c_ratio','d-1','labile_detritus_n_to_c_ratio')
   if (NCrSdet < 0) call self%register_horizontal_dependency(self%id_NCrSdet, 'semilabile_detritus_n_to_c_ratio','d-1','semilabile_detritus_n_to_c_ratio')
   if (PAds < 0) call self%register_horizontal_dependency(self%id_PAds, 'phosphorous_adsorption_coefficient','1','phosphorous_adsorption_coefficient')
   if (NH3Ads < 0) call self%register_horizontal_dependency(self%id_NH3Ads, 'ammonia_adsorption_coefficient','1','ammonia_adsorption_coefficient')
   if (PAdsODU < 0) call self%register_horizontal_dependency(self%id_PAdsODU, 'phosphate_iron_dissolution_threshold_ratio','mmolFe mmolODU-1','phosphate_iron_dissolution_threshold_ratio')

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
!  Modified by author(s): Carsten Lemmen
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
   real(rk), parameter :: gammaNO3 = 0.8_rk    ! molN  used per molC in denitrification     !RH 0.8-> ~104/106?
   real(rk), parameter :: gammaO2  = 1.0_rk    ! molO2 used per molC in oxic mineralization !RH 1.0->150/106*OxicMin (if [oxy]=mmolO2/m**3)
   real(rk), parameter :: gammaNH3 = 2.0_rk    ! molO2 needed to oxidize one molNH3 in nitrification
   real(rk), parameter :: T0       = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk), parameter :: Q10b     = 1.5_rk    ! q10 temperature coefficient (source?)
   real(rk), parameter :: relaxO2  = 0.04_rk   ! (source?)
   real(rk) :: ldetC, sdetC, oxy, odu, no3, nh3, detP, po4
   real(rk) :: temp_celsius, temp_kelvin, f_T, E_a
   real(rk) :: radsP, Oxicminlim, Denitrilim, Anoxiclim, Rescale, rP
   real(rk) :: CprodL, CprodS, Cprod, Nprod, Pprod
   real(rk) :: AnoxicMin, Denitrific, OxicMin, Nitri, OduDepo, OduOx, pDepo
   real(rk) :: rLabile, rSemiLabile, NCrLdet, NCrSdet, PAds, PAdsODU, NH3Ads

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
   _GET_(self%id_detP, detP)         ! detritus P                 (mmolP/m**3)
   _GET_(self%id_oxy, oxy)           ! dissolved oxygen           (mmolO2/m**3)
   _GET_(self%id_odu, odu)           ! dissolved oxygen demand units (mmolO2/m**3)
   _GET_(self%id_no3, no3)           ! dissolved nitrate          (mmolN/m**3)
   _GET_(self%id_nh3, nh3)           ! dissolved ammonium         (mmolN/m**3)
   _GET_(self%id_po4, po4)           ! dissolved phosphate        (mmolP/m**3)

   rLabile = self%rLabile
   if (rLabile < 0) _GET_HORIZONTAL_(self%id_rLabile, rLabile)
   rSemiLabile = self%rSemilabile
   if (rSemiLabile < 0) _GET_HORIZONTAL_(self%id_rSemiLabile, rSemiLabile)
   PAds = self%PAds
   if (PAds < 0) _GET_HORIZONTAL_(self%id_PAds, PAds)
   PAdsODU = self%PadsODU
   if (PAdsODU < 0) _GET_HORIZONTAL_(self%id_PAdsODU, PAdsODU)
   NH3Ads = self%NH3Ads
   if (NH3Ads < 0) _GET_HORIZONTAL_(self%id_NH3Ads, NH3Ads)
   NCrLdet = self%NCrLdet
   if (NCrLdet < 0) _GET_HORIZONTAL_(self%id_NCrLdet, NCrLdet)
   NCrSdet = self%NCrSdet
   if (NCrSdet < 0) _GET_HORIZONTAL_(self%id_NCrSdet, NCrSdet)

#ifdef debugMK
if ( temp_celsius /= temp_celsius ) then; write(0,*) 'omexdia#254 temp_celsius = ',temp_celsius; stop; endif
if ( ldetC        /= ldetC        ) then; write(0,*) 'omexdia#255 ldetC = ',ldetC     ; stop; endif
if ( sdetC        /= sdetC        ) then; write(0,*) 'omexdia#256 sdetC = ',sdetC     ; stop; endif
if ( detP         /= detP         ) then; write(0,*) 'omexdia#257 detP = ',detP       ; stop; endif
if ( oxy          /= oxy          ) then; write(0,*) 'omexdia#258 oxy = ',oxy         ; stop; endif
if ( odu          /= odu          ) then; write(0,*) 'omexdia#259 odu = ',odu         ; stop; endif
if ( no3          /= no3          ) then; write(0,*) 'omexdia#260 no3 = ',no3         ; stop; endif
if ( nh3          /= nh3          ) then; write(0,*) 'omexdia#261 nh3 = ',nh3         ; stop; endif
if ( po4          /= po4          ) then; write(0,*) 'omexdia#262 po4 = ',po4         ; stop; endif
#endif

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

   CprodL = f_T * rLabile     * ldetC   ! (mmolC m-3 d-1)
   CprodS = f_T * rSemilabile * sdetC   ! (mmolC m-3 d-1)

! assume upper reactive surface area for POC hydrolysis (introduced by Kai Wirtz ?)
   if (CprodS > self%CprodMax) CprodS = self%CprodMax
   if (CprodL > self%CprodMax) CprodL = self%CprodMax

   Cprod  = CprodL + CprodS                               ! (mmolC m-3 d-1)
   Nprod  = CprodL * NCrLdet + CprodS * NCrSdet ! (mmolN m-3 d-1) !MK: how to implement N/C-ratio from (coupled) external model?

! PO4-adsorption ceases when critical capacity is reached
! [FeS] approximated by ODU
! TODO: temperature dependency
   radsP  = PAds  * po4 * 1.0_rk/(1.0_rk+exp(-5.0_rk+(odu-oxy)/PAdsODU)) ! (?)
   rP     = f_T * rLabile * 2 ! (1.0_rk - Oxicminlim)  ! (d-1)
   Pprod  = rP * detP                                     ! (mmolP m-3 d-1)

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

#ifdef debugMK
if ( CprodL /= CprodL ) then; write(0,*) 'omexdia#302 CprodL = ',CprodL; stop; endif
if ( CprodS /= CprodS ) then; write(0,*) 'omexdia#302 CprodS = ',CprodS; stop; endif
if ( OxicMin /= OxicMin ) then; write(0,*) 'omexdia#302 OxicMin = ',OxicMin; stop; endif
if ( Nitri /= Nitri ) then; write(0,*) 'omexdia#302 Nitri = ',Nitri; stop; endif
if ( OduOx /= OduOx ) then; write(0,*) 'omexdia#302 OduOx = ',OduOx; stop; endif
if ( Denitrific /= Denitrific ) then; write(0,*) 'omexdia#302 Denitrific = ',Denitrific; stop; endif
if ( Nprod /= Nprod ) then; write(0,*) 'omexdia#302 Nprod = ',Nprod; stop; endif
if ( AnoxicMin /= AnoxicMin ) then; write(0,*) 'omexdia#302 AnoxicMin = ',AnoxicMin; stop; endif
if ( Nprod /= Nprod ) then; write(0,*) 'omexdia#302 Nprod = ',Nprod; stop; endif
if ( Pprod /= Pprod ) then; write(0,*) 'omexdia#302 Pprod = ',Pprod; stop; endif
if ( radsP /= radsP ) then; write(0,*) 'omexdia#302 radsP = ',radsP; stop; endif
#endif

#define _CONV_UNIT_ *one_pr_day
! reaction rates
   _SET_ODE_(self%id_ldetC, -CprodL                                        _CONV_UNIT_)  ! (mmolC  m-3 d-1)
   _SET_ODE_(self%id_sdetC, -CprodS                                        _CONV_UNIT_)  ! (mmolC  m-3 d-1)
   _SET_ODE_(self%id_oxy,  (-OxicMin*gammaO2     - Nitri*gammaNH3 - OduOx) _CONV_UNIT_)  ! (mmolO2 m-3 d-1)
   _SET_ODE_(self%id_no3,  (-Denitrific*gammaNO3 + Nitri)                  _CONV_UNIT_)  ! (mmolN  m-3 d-1)
!   _SET_ODE_(self%id_nh3,  (Nprod - Nitri) / (1.0_rk + self%NH3Ads)       _CONV_UNIT_)  ! (mmolN  m-3 d-1)
   _SET_ODE_(self%id_nh3,  (Nprod - Nitri)                                 _CONV_UNIT_)  ! (mmolN  m-3 d-1)
   _SET_ODE_(self%id_odu,  (AnoxicMin - OduOx - OduDepo)                   _CONV_UNIT_)  ! (mmolO2 m-3 d-1)
   _SET_ODE_(self%id_po4,  (Pprod - radsP)                                 _CONV_UNIT_)  ! (mmolP  m-3 d-1)
!   _SET_ODE_(self%id_detP, (radsP - Pprod - self%NH3Ads*detP)              _CONV_UNIT_)  ! (mmolP  m-3 d-1)
   _SET_ODE_(self%id_detP, (radsP - Pprod)                                 _CONV_UNIT_)  ! (mmolP  m-3 d-1)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit, Denitrific*gammaNO3)  !last denitrification rate
   _SET_DIAGNOSTIC_(self%id_adsp, radsP)                 !instantaneous phosphate adsorption

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
   end module hzg_omexdia_p
!-----------------------------------------------------------------------
