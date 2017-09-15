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
! The MPB model from Hochard et al EcoMod 2010 added (adopted) by kai wirtz
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
   real(rk), parameter :: one_pr_day = 1.0_rk / secs_pr_day
!
! !REVISION HISTORY:!
!  Original author(s): Markus Kreus, Richard Hofmeister & Kai Wirtz
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hzg_mpb
!     Variable identifiers
      type (type_state_variable_id)        :: id_no3, id_nh4, id_oxy, id_ldet
      type (type_state_variable_id)        :: id_dic, id_zbC, id_zbN
      type (type_state_variable_id)        :: id_mpbCHL, id_mpbC, id_mpbN, id_eps
      type (type_dependency_id)            :: id_temp, id_parz, id_porosity
      type (type_diagnostic_variable_id)   :: id_PrimProd, id_par, id_Q_N, id_Q_chl
      type (type_diagnostic_variable_id)   :: id_expCProd, id_expNProd
      type (type_diagnostic_variable_id)   :: id_NPP, id_SGR, id_TGR, id_SPR, id_SMR
! temporary for debugging
      type (type_diagnostic_variable_id)   :: id_MPB_din, id_MPB_no3, id_mpb_nh4

!     Model parameters
      real(rk) :: mumax, alpha, gamma, Qmin, Qmax, thetamax, uptmax, KNH4, KNO3
      real(rk) :: KinNH4, kEPS, rEPS, resp, Kresp, graz, kout, kexu, rzoo
      real(rk) :: PAR0max, k0, Achla, btemp
      logical  :: use_no3, use_nh4, use_oxy, use_ldet, use_dic, use_zbC, use_zbN

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
   real(rk)  :: mpbCHL_init   ! MicroPhytoBenthos chlorophyll               (g Chla m-3)
   real(rk)  :: mpbC_init     ! MicroPhytoBenthos carbon                    (mmolC m-3)
   real(rk)  :: mpbN_init     ! MicroPhytoBenthos nitrogen                  (mmolN m-3)
   real(rk)  :: eps_init      ! Extracellular Polymeric Substances          (mmolC m-3)
!!------- Parameters for model mpb originating from omexdia_p -------
!    real(rk)  :: rLdet         ! decay rate labile detritus (fast decay)     (0.075  d-1)
!    real(rk)  :: rSdet         ! decay rate semilabile detritus (slow decay) (0.003  d-1)
!    real(rk)  :: rNCldet       ! N/C ratio labile detritus (fast decay)      (0.230  molN molC-1)
!!------- Parameters for model mpb -------
   real(rk)  :: mumax         ! Maximum growth rate                         (0.70  d-1)
   real(rk)  :: btemp         ! Temperature factor                          (0.063  - )
   real(rk)  :: alpha         ! initial slope of the PI-curve               (1.5e-5  m2 molC (gChla J)-1)
   real(rk)  :: gamma         ! Mol O2 produced per mol C fixed by photosynthesis (1.0  molO2 molC-1)
   real(rk)  :: Qmin          ! Minimum N/C ratio                           (0.05  molN molC-1)
   real(rk)  :: Qmax          ! Maximum N/C ratio                           (0.20  molN molC-1)
   real(rk)  :: thetamax      ! Maximum Chla/N ratio                        (3.80  gChla molN-1)
   real(rk)  :: uptmax        ! Maximum N uptake rate per carbon unit       (0.20  molN molC-1)
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
   real(rk)  :: k0            ! Extinction coefficient of sediment          (20.0  cm-1)
   real(rk)  :: Achla         ! Absorption factor of chlorophyll            (0.02  m2 gChla-1)
!!------- Optional external dependencies -------
   character(len=attribute_length) :: no3_variable  = ''  ! dissolved nitrate
   character(len=attribute_length) :: nh4_variable  = ''  ! dissolved ammonium
   character(len=attribute_length) :: oxy_variable  = ''  ! dissolved oxygen
   character(len=attribute_length) :: ldet_variable = ''  ! labile detritus carbon (fast decay)
   character(len=attribute_length) :: dic_variable  = ''  ! dissolved inorganic carbon
   character(len=attribute_length) :: zbC_variable  = ''  ! zoobenthos carbon
   character(len=attribute_length) :: zbN_variable  = ''  ! zoobenthos nitrogen

   namelist /hzg_mpb/ &
          mumax, btemp, alpha, gamma, Qmin, Qmax, thetamax, uptmax, KNH4, KNO3,   &
          KinNH4, kEPS, rEPS, resp, Kresp, graz, kout, kexu, rzoo, k0, Achla,     &
          mpbCHL_init, mpbC_init, mpbN_init, eps_init

   namelist /hzg_mpb_dependencies/  &
          no3_variable, nh4_variable, oxy_variable, ldet_variable, &
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
   call self%register_state_variable(self%id_mpbCHL, 'mpbCHL', 'gChla m-3', &
         'MicroPhytoBenthos chlorophyll mpbCHL',                            &
         mpbCHL_init, minimum=0.0_rk, no_river_dilution=.true.)
   call self%set_variable_property(self%id_mpbCHL, 'particulate', .true.)

   call self%register_state_variable(self%id_mpbC,   'mpbC', 'mmolC m-3',   &
         'MicroPhytoBenthos carbon mpbC',                                   &
         mpbC_init, minimum=0.0_rk, no_river_dilution=.true.)
   call self%set_variable_property(self%id_mpbC, 'particulate', .true.)

   call self%register_state_variable(self%id_mpbN,   'mpbN', 'mmolN m-3',   &
         'MicroPhytoBenthos nitrogen mpbN',                                 &
         mpbN_init, minimum=0.0_rk, no_river_dilution=.true.)
   call self%set_variable_property(self%id_mpbN, 'particulate', .true.)

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

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_PrimProd, 'MPB_PP', 'mmolC m-3 d-1',  &
         'MPB primary production rate PrimProd',     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_par,      'par', 'W m-2',               &
         'MPB photosynthetically active radiation',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Q_N,      'Q_N', 'mmolN mmolC-1',       &
         'MPB nitrogen quota Q_N',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_Q_chl,    'Q_chl', 'gChla mmolC-1',     &
         'MPB CHL:C ratio Q_chl',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_expCProd, 'expCProd', 'mmolC m-2 d-1',  &
         'MPB carbon export (zoobenthos grazing)',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_expNProd, 'expNProd', 'mmolN m-2 d-1',  &
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
   call self%register_diagnostic_variable(self%id_MPB_DIN, 'DIN', 'mmolN m-2 d-1',        &
         'MPB DIN',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_no3, 'MPB_no3', 'mmolN m-2 d-1',    &
         'MPB no3',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_MPB_nh4, 'MPB_nh4', 'mmolN m-2 d-1',    &
         'MPB nh4',  output=output_instantaneous)

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
   real(rk), parameter :: T0       = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk), parameter :: Q10b     = 1.5_rk    ! q10 temperature coefficient (source?)
   real(rk) :: mpbC, mpbN, mpbCHL, eps, no3, nh4, oxy, ldet
   real(rk) :: temp_celsius, temp_kelvin, f_T, E_a, parz, porosity, CprodEPS
   real(rk) :: prodChl, k, theta, Q_N, Q_chl, Pmax,  PP, prod, prodeps, fac
   real(rk) :: prodO2, rhochl, uptNH4, uptNO3, uptchl, uptN, respphyto, faecesC, faecesN
   real(rk) :: exud, grazingC, grazingN, grazingChl, respzoo, exportC, exportN
   real(rk) :: NPP, SGR, TGR, GRZ

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! initialize local external dependencies
   no3  = -9999.
   nh4  = -9999.
   oxy  = -9999.
   ldet = -9999.

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp, temp_celsius) ! sediment-water temperature
   _GET_(self%id_parz,   parz)    ! sediment light
   _GET_(self%id_mpbCHL, mpbCHL)  ! MicroPhytoBenthos chlorophyll in mmolN m-3
   _GET_(self%id_mpbC,   mpbC)    ! MicroPhytoBenthos carbon in mmolC m-3
   _GET_(self%id_mpbN,   mpbN)    ! MicroPhytoBenthos nitrogen in mmolN m-3
   _GET_(self%id_eps,    eps)     ! Extracellular Polymeric Substances  in mmolC m-3
   ! Retrieve current external (local) dependencies if available
   _GET_(self%id_no3,    no3)     ! dissolved nitrate in mmolN m-3
   _GET_(self%id_nh4,    nh4)     ! dissolved ammonium in mmolN m-3
   _GET_(self%id_oxy,    oxy)     ! dissolved oxygen in mmolO2 m-3
   _GET_(self%id_ldet,   ldet)    ! fast decaying detritus C in mmolC m-3

   ! temperature dependency
   temp_kelvin = 273.15_rk + temp_celsius
   E_a = 0.1_rk*log(Q10b)*T0*(T0+10.0_rk)
   f_T = 1.0_rk*exp(-E_a*(1.0_rk/temp_kelvin - 1.0_rk/T0)) ! temperature factor (-)

   !----- intracellular N:C:Chl stoichiometry
   Q_N    = mpbN   / mpbC      ! mmolN mmolC-1
   Q_chl  = mpbCHL / mpbC      ! gChla mmolC-1
   theta  = mpbCHL / mpbN      ! gChla mmolN-1

   !----- photosnythesis rate
   fac    = max( 0.0, 1.0 - self%Qmin/Q_N)
   Pmax   = self%mumax * exp(self%btemp * temp_celsius) * fac
   PP     = Pmax * (1.0 - exp(-self%alpha * parz*Q_chl/Pmax))
   prod   = PP * mpbC
   prodO2 = self%gamma * prod

   rhochl = 0.0
   if (parz .gt. 0.0) rhochl = self%thetamax * PP/(self%alpha * parz * Q_chl)

   !----- nutrient uptake (TODO add dependencies on P, Si)
   uptNO3 = self%uptmax * no3/(no3 + self%KNO3)
   uptNH4 = self%uptmax * nh4/(nh4 + self%KNH4)
   uptN   = uptNO3 + uptNH4
   !----- chlorophyll synthesis
   fac    = max(0.0, 1.0 - theta/self%thetamax)
   uptchl = uptN * fac/(fac + 0.05)
   prodChl= rhochl * uptchl * mpbC

   fac    = max(0.0, (self%Qmax - Q_N) / (self%Qmax - self%Qmin))
   uptNO3 = uptNO3 * fac * (1.0 - nh4/(nh4 + self%KinNH4)) * mpbC
   uptNH4 = uptNH4 * fac * mpbC
   uptN   = uptNO3 + uptNH4

   ! Carbohydrate exudation:
   prodeps = self%keps * prod
   !CprodEPS = sqrt(self%rLdet*self%rSdet) * eps  ! Source of this formulation? Kai Wirtz ??
   CprodEPS = self%rEPS * eps

   ! Respiration:
   respphyto = self%resp * mpbC * oxy/(oxy+self%Kresp)

   ! Zoobenthos grazing and associated processes:
   grazingC   = self%graz * mpbC
   grazingN   = self%graz * mpbN
   grazingChl = self%graz * mpbCHL
   faecesC    = self%kout * grazingC
   !faecesN    = self%kout * grazingC * 0.23 !(:=rNCldet) ! Hochard et al 2010
   faecesN    = self%kout * grazingN ! alternative proposed by Markus Kreus
   exud       = self%kexu * grazingN
   respzoo    = self%rzoo * grazingC * oxy/(oxy+self%Kresp)
   exportC    = grazingC - faecesC - respzoo ! exported carbon   (open closure term)
   exportN    = grazingN - faecesN - exud    ! exported nitrogen (open closure term)

   NPP =  prod - respphyto
   SGR = (prod - respphyto - prodeps ) /mpbC
   TGR = (prod - grazingC - respphyto - prodeps) /mpbC
   GRZ = (grazingC) /mpbC

#define _CONV_UNIT_ *one_pr_day
   ! reaction rates
   _SET_ODE_(self%id_mpbCHL, (prodchl - grazingChl ) _CONV_UNIT_)
   _SET_ODE_(self%id_mpbC,   (prod - grazingC - respphyto - prodeps ) _CONV_UNIT_)
   _SET_ODE_(self%id_mpbN,   (uptN - grazingN ) _CONV_UNIT_)
   _SET_ODE_(self%id_eps,    (prodeps - f_T * CprodEPS ) _CONV_UNIT_)
   ! external dependencies
   ! If externally maintained variables are present, change the pools accordingly
   _SET_ODE_(self%id_no3 , (- uptNO3) _CONV_UNIT_)
   _SET_ODE_(self%id_nh4 , (exud - uptNH4) _CONV_UNIT_)
   _SET_ODE_(self%id_oxy , (- respzoo - respphyto) _CONV_UNIT_)
   _SET_ODE_(self%id_ldet, (faecesC) _CONV_UNIT_)
   if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic , (respzoo + respphyto) _CONV_UNIT_)
   if (_AVAILABLE_(self%id_zbC)) _SET_ODE_(self%id_zbC , (exportC) _CONV_UNIT_)
   if (_AVAILABLE_(self%id_zbN)) _SET_ODE_(self%id_zbN , (exportN) _CONV_UNIT_)

   ! Export diagnostic variables
   !_SET_DIAGNOSTIC_(self%id_PrimProd, PP)            !instantaneous MPB primary production rate (d-1)
   _SET_DIAGNOSTIC_(self%id_PrimProd, prod)          !instantaneous MPB primary production rate (molC m-3 d-1)
   _SET_DIAGNOSTIC_(self%id_par, parz)               !instantaneous MPB photosynthetically active radiation
   _SET_DIAGNOSTIC_(self%id_Q_N, Q_N)                !instantaneous MPB N:C quota
   _SET_DIAGNOSTIC_(self%id_Q_chl, Q_chl)            !instantaneous MPB CHL:C ratio
   _SET_DIAGNOSTIC_(self%id_expCProd, exportC)       !instantaneous MPB export production carbon
   _SET_DIAGNOSTIC_(self%id_expNProd, exportN)       !instantaneous MPB export production nitrogen
   _SET_DIAGNOSTIC_(self%id_NPP, NPP)                !instantaneous MPB net primary production rate (molC m-3 d-1)
   _SET_DIAGNOSTIC_(self%id_SGR, SGR)                !instantaneous MPB specific growth rate (d-1)
   _SET_DIAGNOSTIC_(self%id_TGR, TGR)                !instantaneous MPB total growth rate (d-1)
   _SET_DIAGNOSTIC_(self%id_SPR, PP)                 !instantaneous MPB specific photosynthesis rate (d-1)
   _SET_DIAGNOSTIC_(self%id_SMR, GRZ)                !instantaneous MPB specific mortality rate (d-1)
! temporary for debugging:
   _SET_DIAGNOSTIC_(self%id_mpb_din, nh4+no3)        !instantaneous external DIN
   _SET_DIAGNOSTIC_(self%id_mpb_no3, no3)            !instantaneous external no3
   _SET_DIAGNOSTIC_(self%id_mpb_nh4, nh4)            !instantaneous external nh4

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
  end module hzg_mpb
!-----------------------------------------------------------------------
