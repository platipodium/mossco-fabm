#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_npzd --- NPZD biogeochemical model based upon
! Fennel \& Neumann (1996, German Journal of Hydrography),
! with minor modifications by Burchard et al. (2005, Ocean Dynamics)
! taken from GOTM and adapted for FABM by Jorn Bruggeman
!
! !INTERFACE:
   module fabm_hzg_n2pzdq
!
! !DESCRIPTION:
! The NPZD (nutrient-phytoplankton-zooplankton-detritus) model described here
! consists of $I=4$ state variables.
! Nutrient uptake (phytoplankton growth) is limited by light and nutrient
! availability, the latter of which is modelled by means
! of Michaelis-Menten kinetics, see eq.\ (\ref{dnp}).
! The half-saturation nutrient concentration $\alpha$ used in this
! formulation has typically a value between 0.2 and 1.5 mmol N\, m$^{-3}$.
! Zooplankton grazing which is limited by the phytoplankton standing stock
! is modelled by means of an Ivlev formulation, see eq.\ (\ref{dpz}).
! All other processes are based on linear first-order kinematics,
! see eqs.\ (\ref{dpn}) - (\ref{dzd}).
! For all details of the NPZD model implemented here,
! see \cite{Burchardetal2005b}.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
!   public type_gotm_n2pzdq, gotm_n2pzdq_init, gotm_n2pzdq_do, & !gotm_n2pzdq_do_ppdd, &
!          gotm_n2pzdq_get_light_extinction, gotm_n2pzdq_get_conserved_quantities
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Lena Spruch, Kai Wirtz, Onur Kerimoglu
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public    ::  type_hzg_n2pzdq
!     Variable identifiers
      type (type_state_variable_id)        :: id_DIN,id_DIP,id_phyC,id_phyN,id_phyP,id_detN,id_detP,id_zooC
      type (type_state_variable_id)        :: id_dic
      type (type_dependency_id)            :: id_par
      type (type_dependency_id)            :: id_temp
      type (type_horizontal_dependency_id) :: id_I_0
      type (type_diagnostic_variable_id)   :: id_GPP,id_NCP,id_PPR,id_NPR,id_dPAR,id_dMort,id_dLlim,id_dNlim,id_dPlim,id_dqnc,id_dqpc,id_deN,id_deP,id_deC,id_dgraz,id_dmortz
      type (type_conserved_quantity_id)    :: id_totN,id_totP

!     Model parameters
      real(rk) :: p0,affin_par,upmax_N,upmax_P,grow_max,iv,halfsatN,halfsatP, & 
                  rem_N,rem_P,mort0_phy,mortpar_phy,qmax_N,qmax_P,qmin_N,qmin_P,kc,w_p,w_d,rpn,grazmax,mort_zoo,n,qzn,qzp,eff,e_C,k_detN,k_detP,mort_zoo2,n2,zexcdetfr 
      real(rk) :: dic_per_n
      logical  :: use_dic
      
            contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: get_light_extinction
      
   end type type_hzg_n2pzdq
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the NPZD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hzg_n2pzdq), intent(inout), target :: self
   integer,                 intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)          :: DIN_initial
   real(rk)          :: DIP_initial
   real(rk)          :: phyC_initial
   real(rk)          :: phyN_initial
   real(rk)          :: phyP_initial
   real(rk)          :: detN_initial
   real(rk)          :: detP_initial
   real(rk)          :: zooC_initial
   real(rk)          :: p0
   real(rk)          :: w_p
   real(rk)          :: w_d
   real(rk)          :: affin_par
   real(rk)          :: upmax_N
   real(rk)          :: upmax_P
   real(rk)          :: grow_max
   real(rk)          :: halfsatN
   real(rk)          :: halfsatP
   real(rk)          :: rpn
   real(rk)          :: rem_N
   real(rk)          :: rem_P
   real(rk)          :: mort0_phy
   real(rk)          :: mortpar_phy
   real(rk)          :: qmax_N
   real(rk)          :: qmax_P
   real(rk)          :: qmin_N
   real(rk)          :: qmin_P
   real(rk)          :: iv
   real(rk)          :: kc
   real(rk)          :: grazmax
   real(rk)          :: mort_zoo
   real(rk)          :: dic_per_n
   real(rk)          :: n
   real(rk)          :: mort_zoo2
   real(rk)          :: n2
   real(rk)          :: qzn
   real(rk)          :: qzp
   real(rk)          :: eff
   real(rk)          :: e_C
   real(rk)          :: zexcdetfr
   real(rk)          :: k_detN
   real(rk)          :: k_detP
   
   character(len=64) :: dic_variable

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   namelist /hzg_n2pzdq/ DIN_initial,DIP_initial,phyC_initial,phyN_initial,phyP_initial,detN_initial,detP_initial,zooC_initial,  &
                        p0,w_p,w_d,affin_par,upmax_N,upmax_P,grow_max,iv,halfsatN,halfsatP,rpn,  &
                        rem_N,rem_P,mort0_phy,mortpar_phy,qmax_N,qmax_P,qmin_N,qmin_P,kc,grazmax,mort_zoo,n,qzn,qzp,eff,e_C,k_detN,k_detP,mort_zoo2,n2,zexcdetfr
!EOP
!-----------------------------------------------------------------------
!BOC
   DIN_initial  = 4.5_rk        ! mmol-N/m³        , initial N concentration
   DIP_initial  = 4.5_rk        ! mmol-P/m³        , initial P concentration
   phyC_initial = 0.0_rk        ! mmol-C/m³        , initial phytoplankton C concentration
   phyN_initial = 0.0_rk        ! mmol-N/m³        , initial phytoplankton N concentration
   phyP_initial = 0.0_rk        ! mmol-N/m³        , initial phytoplankton N concentration
   detN_initial = 0.0_rk        ! mmol-N/m³        , initial phytoplankton N concentration
   detP_initial = 4.5_rk        ! mmol-N/m³        , initial detritus N concentration
   zooC_initial = 0.0_rk        ! mmol-C/m³        , initial zooplankton C concentration
   p0           = 0.0225_rk     ! mmol-C/m³        , minimum phytoplankton concentration
   w_p          = -1.0_rk       ! m/d              , settling velocity phytoplankton
   w_d          = -5.0_rk       ! m/d              , settling velocity detritus
   affin_par    = 25.0_rk       ! m²/(W*d)         , initial slope of light-prod-curve
   upmax_N      = 1.0_rk        ! mmol-N/(mmol-C*d), maximum nutrient(N) uptake rate
   upmax_P      = 1.0_rk        ! mmol-P/(mmol-C*d), maximum nutrient(P) uptake rate
   grow_max     = 0.5_rk        ! 1/d              , maximum phytoplankton growth rate
   iv           = 1.1_rk        !                  , ivlev constant
   halfsatN     = 0.3_rk        ! mmol-N/m³        , half saturation constant for N
   halfsatP     = 0.3_rk        ! mmol-N/m³        , half saturation constant for N
   rpn          = 0.01_rk       ! 1/d              , exudation rate of phytoplankton
   rem_N        = 0.003_rk      ! 1/d              , remineralisation rate, detN to N
   rem_P        = 0.02_rk       ! 1/d              , maximum phytoplankton mortality rate
   qmax_N       = 0.5_rk        ! mmol-N/mmol-C    , maximum N/C quota
   qmax_P       = 0.5_rk        ! mmol-N/mmol-C    , maximum N/C quota
   qmin_N       = 0.1_rk        ! mmol-N/mmol-C    , minimum N/C quota
   qmin_N       = 0.1_rk        ! mmol-N/mmol-C    , minimum N/C quota
   kc           = 0.03_rk       ! m²/mmol-C        , attenuation constant for self-shading
   mort0_phy    = 0.5_rk        ! 1/d              , maximum mortality phytoplankton mortality
   mortpar_phy  = 1.0_rk        !                  , mortality parameter phytoplankton
   grazmax      = 0.1_rk        ! 1/d              , maximum grazing rate of zooplankton
   mort_zoo     = 0.1_rk        ! 1/d              , maximum mortality rate for zooplankton 1
   n            = 1.0_rk        !                  , exponent for zooplankton mortality 1
   mort_zoo2    = 0.1_rk        ! 1/d              , maximum mortality rate for zooplankton 2
   n2           = 1.0_rk        !                  , exponent for zooplankton mortality 2
   qzn          = 0.5_rk        ! mmol-N/mmol-C    , zooplankton N/C quota
   qzp          = 0.5_rk        ! mmol-P/mmol-C    , zooplankton P/C quota
   e_C          = 0.4_rk        !                  , C grazing/uptake efficiency of zooplankton
   zexcdetfr    = 0.1_rk        !                  , detritus fraction of zooplankton excretion
   k_detN       = 1.0_rk        !
   k_detP       = 1.0_rk        !
   dic_per_n    = 106.0_rk/16.0_rk  ! Redfield C:N
   dic_variable = '' 
   
 
   ! Read the namelist
   read(configunit,nml=hzg_n2pzdq,err=99,end=100)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.

   call self%get_parameter(self%p0,         'p0',          default=p0)
   call self%get_parameter(self%affin_par,  'affinpar',    default=affin_par)
   call self%get_parameter(self%upmax_N,    'upmax_N',     default=upmax_N,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%upmax_P,    'upmax_P',     default=upmax_P,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%grow_max,   'grow_max',    default=grow_max,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%iv,         'iv',          default= iv)
   call self%get_parameter(self%halfsatN,   'halfsatN',    default=halfsatN)
   call self%get_parameter(self%halfsatP,   'halfsatP',    default=halfsatP)
   call self%get_parameter(self%rpn,        'rpn',         default=rpn,         scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%rem_N,      'rem_N',       default=rem_N,       scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter( self%rem_P,     'rem_P',       default=rem_P,       scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%mort0_phy,  'mort0_phy',   default=mort0_phy,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%mortpar_phy,'mortpar_phy', default=mortpar_phy)
   call self%get_parameter(self%qmax_N ,    'qmax_N',      default=qmax_N)
   call self%get_parameter(self%qmax_P,     'qmax_P',      default=qmax_P)
   call self%get_parameter(self%qmin_N,     'qmin_N',      default=qmin_N)
   call self%get_parameter(self%qmin_P,     'qmin_P',      default=qmin_P)
   call self%get_parameter(self%kc,         'kc',          default=kc)
   call self%get_parameter(self%grazmax,    'grazmax',     default=grazmax,     scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%mort_zoo ,  'mort_zoo',    default=mort_zoo,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%mort_zoo2,  'mort_zoo2',   default=mort_zoo2,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%n ,         'n',           default=n)
   call self%get_parameter(self%n2,         'n2',          default=n2)
   call self%get_parameter(self%dic_per_n,  'dic_per_n',   default=dic_per_n)
   call self%get_parameter(self%qzn,        'qzn',         default=qzn)
   call self%get_parameter(self%qzp,        'qzp',         default=qzp)
   call self%get_parameter(self%eff,        'eff',         default=eff)
   call self%get_parameter(self%e_C,        'e_C',         default=e_C)
   call self%get_parameter(self%k_detN,     'k_detN',      default=k_detN)
   call self%get_parameter(self%k_detP,     'k_detP',      default=k_detP)
   call self%get_parameter(self%zexcdetfr,  'zexcdetfr',   default=zexcdetfr)
   
   ! Register state variables
   call self%register_state_variable(self%id_DIN,'nut_N','mmol/m**3','DIN',   &
                                    DIN_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_DIP,'nut_P','mmol/m**3','DIP',   &
                                    DIP_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_phyC,'phy_C','mmol/m**3','phyC', &
                                    phyC_initial,minimum=0.0_rk,vertical_movement=w_p/secs_pr_day)
   call self%register_state_variable(self%id_phyN,'phy_N','mmol/m**3','phyN', &
                                    phyN_initial,minimum=0.0_rk,vertical_movement=w_p/secs_pr_day)
   call self%register_state_variable(self%id_phyP,'phy_P','mmol/m**3','phyP', &
                                    phyP_initial,minimum=0.0_rk,vertical_movement=w_p/secs_pr_day)
   call self%register_state_variable(self%id_detN,'det_N','mmol/m**3','detN', &
                                    detN_initial,minimum=0.0_rk,vertical_movement=w_d/secs_pr_day)
   call self%register_state_variable(self%id_detP,'det_P','mmol/m**3','detP', &
                                    detP_initial,minimum=0.0_rk,vertical_movement=w_d/secs_pr_day)
   call self%register_state_variable(self%id_zooC,'zoo_C','mmol/m**3','zooC', &
                                    zooC_initial,minimum=0.0_rk)
   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   !self%use_dic = dic_variable/=''
   !if (self%use_dic) call register_state_dependency(modelinfo,self%id_dic,dic_variable)
   if (self%use_dic) call self%register_state_dependency(self%id_dic,dic_variable)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_GPP,'GPP','mmol/m**3',  'gross primary production',            &
                     output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_dLlim,'Llim','',  'Llim',                          &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_dNlim,'Nlim','',  'Nlim',                       &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_dPlim,'Plim','',  'Plim',                     &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_dMort,'sp. phyto mortality rate','1/d',  'mort_phy',                         &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W/m**2',     'PAR',&
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_dqnc,'QNC','mmolN/mmolC',     'QNC',                     &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_dqpc,'QPC','mmolP/mmolC',     'QPC',                     &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_deN,'eN','',  'N_eff',           &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_deP,'eP','',  'P_eff',           &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_deC,'eC','',  'C_eff',           &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_dgraz,'grazing rate','mmolC/m**3/d',  'graz',           &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_dmortz,'zoo mortality rate','mmolC/m**3/d',  'mort_zoo',           &
                     output=output_time_step_averaged)
!   call register_diagnostic_variable(modelinfo,self%id_PPR,'PPR','mmol/m**3/d','gross primary production rate',      &
!                     time_treatment=time_treatment_averaged)
!   call register_diagnostic_variable(modelinfo,self%id_NPR,'NPR','mmol/m**3/d','net community production rate',      &
!                     time_treatment=time_treatment_averaged)
!   call register_diagnostic_variable(modelinfo,self%id_dPAR,'PAR','W/m**2',     'photosynthetically active radiation',&
!                     time_treatment=time_treatment_averaged)

   ! Add to conserved quantities
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_DIN)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_phyN)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_detN)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_zooC,scale_factor=self%qzn)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_DIP)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_phyP)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_detP)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_zooC,scale_factor=self%qzn)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)   
   call self%register_dependency(self%id_I_0,standard_variables%surface_downwelling_photosynthetic_radiative_flux)
    
   return

99 call self%fatal_error('hzg_n2pzdq_init','Error reading namelist hzg_n2pzdq.')

100 call self%fatal_error('hzg_n2pzdq_init','Namelist hzg_n2pzdq was not found.')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
! Seven processes expressed as sink terms are included in this
! conservative model, see eqs.\ (\ref{dnp}) - (\ref{dzd}). \\
!
! Nutrient uptake by phytoplankton:
! \begin{equation}\label{dnp}
! d_{np} = r_{\max}\frac{I_{PAR}}{I_{opt}}
! \exp\left(1-\frac{I_{PAR}}{I_{opt}}\right)
! \frac{c_n}{\alpha+c_n}c_p
! \end{equation}
!
! with
!
! \begin{equation}
! I_{opt}=\max\left(\frac14I_{PAR},I_{\min}\right).
! \end{equation}
!
! Grazing of zooplankton on phytoplankton:
! \begin{equation}\label{dpz}
! d_{pz}=g_{\max}\left(1-\exp\left(-I_v^2c_p^2\right)\right)c_z
! \end{equation}
!
! Phytoplankton excretion:
! \begin{equation}\label{dpn}
! d_{pn} = r_{pn} c_p
! \end{equation}
!
! Zooplankton excretion:
! \begin{equation}\label{dzn}
! d_{zn} = r_{zn} c_z
! \end{equation}
!
! Remineralisation of detritus into nutrients:
! \begin{equation}\label{ddn}
! d_{dn} = r_{dn} c_d
! \end{equation}
!
! Phytoplankton mortality:
! \begin{equation}\label{dpd}
! d_{pd} = r_{pd} c_p
! \end{equation}
!
! Zooplankton mortality:
! \begin{equation}\label{dzd}
! d_{zd} = r_{zd} c_z
! \end{equation}
!
! !INPUT PARAMETERS:
   class (type_hzg_n2pzdq),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
   
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)                   :: DIN,DIP,phyC,phyN,phyP,detN,detP,zooC,par,temp_fact,I_0
   real(rk)                   :: qnc,qpc,dn,reminN,reminP,mort_phy,mort_Z,uptakeN,uptakeP,grazing,primprod,Nlim,Plim,Llim,ee_C,e_N,e_P, &
   dp,dphyC,dphyN,dphyP,ddetN,ddetP,dzooC
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_DIN,DIN)   ! dissolved inorganic N
   _GET_(self%id_DIP,DIP)   ! dissolved inorganic P
   _GET_(self%id_phyC,phyC) ! phytoplankton C
   _GET_(self%id_phyN,phyN) ! phytoplankton N
   _GET_(self%id_phyP,phyP) ! phytoplankton P
   _GET_(self%id_detN,detN) ! detritus N
   _GET_(self%id_detP,detP) ! detritus P
   _GET_(self%id_zooC,zooC) ! zooplankton C

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)               ! local photosynthetically active radiation
   _GET_(self%id_temp,temp_fact)             ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,I_0)    ! surface short wave radiation

   ! Light acclimation formulation based on surface light intensity.
   !iopt = max(0.25*I_0,self%I_min)

   ! Loss rate of phytoplankton to detritus depends on local light intensity.
   !if (par .ge. self%I_min) then
   !  rpd = self%rpdu
   !else
   !  rpd = self%rpdl
   !end if

   ! N/C and P/C cell-quotas 
   qnc = phyN/phyC 
   qpc = phyP/phyC

   ! Define some intermediate quantities that will be reused multiple times.
   !primprod = fprod(self,par,qnc,qpc)

   ! phytoplankton processes
   call fprod(self, par,temp_fact,qnc,qpc, primprod,Nlim,Plim,Llim)
   !Nlim = 1-self%qmin_N/qnc
   !Plim = 1-self%qmin_P/qpc

   !primprod = self%grow_max**Llim
   !mort_phy = self%mort0_phy*exp(-3*(Nlim+Plim+Llim))
   mort_phy = self%mort0_phy*exp(-self%mortpar_phy*(Nlim*Plim))*detN!*Llim
   !mort_phy = self%mort0_phy*exp(-min(Nlim,Plim,Llim))
   !mort_phy = self%mort0_phy 
   uptakeN = fupN(self,DIN,qnc)*temp_fact
   uptakeP = fupP(self,DIP,qpc)*temp_fact
   
   ! pelagic remineralisation
   reminN = self%rem_N*temp_fact !*detN/(detN+self%k_detN)
   reminP = self%rem_P*temp_fact !*detP/(detP+self%k_detP)

   ! zooplankton processes
   !grazing = self%grazmax*(1.0_rk-exp(-self%iv*phyC))
   !grazing = self%grazmax*phyC**2.0/(phyC**2.0+self%iv**2.0)
   grazing = self%grazmax*(1.0_rk-exp(-self%iv*self%iv*phyC*phyC))*temp_fact
   mort_Z=(self%mort_zoo*(zooC**self%n)+self%mort_zoo2*(zooC**self%n2))*temp_fact

   
   ! define uptake efficiencies (e_N,e_P) for N and P, depending on C-efficiency and quotas
   e_N = self%e_C*self%qzn/qnc
   e_P = self%e_C*self%qzp/qpc
   ee_C = self%e_C


  if (max(e_N,e_P) > 1) then
     e_N = 1
     ee_C = qnc/self%qzn
     e_P = ee_C*self%qzp/qpc
    if (e_P > 1) then
       e_P = 1
       ee_C = qpc/self%qzp
       e_N = ee_C*self%qzn/qnc
      if (e_N > 1) then
        e_N = 1
        ee_C = qnc/self%qzn
        e_P = ee_C*self%qzp/qpc
      end if
    end if
  end if

   ! calculate rhs 
   dn    = reminN*detN - uptakeN*phyC+(1-self%zexcdetfr)*(1-e_N)*qnc*grazing*zooC
   dp    = reminP*detP - uptakeP*phyC+(1-self%zexcdetfr)*(1-e_P)*qpc*grazing*zooC
   dphyC = primprod*phyC - mort_phy*phyC - grazing*zooC
   dphyN = uptakeN*phyC - mort_phy*phyN - grazing*qnc*zooC
   dphyP = uptakeP*phyC - mort_phy*phyP - grazing*qpc*zooC
   ddetN = mort_phy*phyN - reminN*detN + mort_Z*self%qzn*zooC + self%zexcdetfr*(1-e_N)*qnc*grazing*zooC
   ddetP = mort_phy*phyP - reminP*detP + mort_Z*self%qzp*zooC + self%zexcdetfr*(1-e_P)*qpc*grazing*zooC
   dzooC = ee_C*grazing*zooC - mort_Z*zooC

   ! Set temporal derivatives
   _SET_ODE_(self%id_DIN,dn)
   _SET_ODE_(self%id_DIP,dp)
   _SET_ODE_(self%id_phyC,dphyC)
   _SET_ODE_(self%id_phyN,dphyN)
   _SET_ODE_(self%id_phyP,dphyP)
   _SET_ODE_(self%id_detN,ddetN)
   _SET_ODE_(self%id_detP,ddetP)
   _SET_ODE_(self%id_zooC,dzooC)

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   if (self%use_dic) then
     _SET_ODE_(self%id_dic,self%dic_per_n*dn)
   end if

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_GPP ,primprod)
   _SET_DIAGNOSTIC_(self%id_dLlim,Llim)
   _SET_DIAGNOSTIC_(self%id_dMort,mort_phy*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_dNlim,Nlim)
   _SET_DIAGNOSTIC_(self%id_dPlim,Plim)
   _SET_DIAGNOSTIC_(self%id_dqnc,qnc)
   _SET_DIAGNOSTIC_(self%id_dqpc,qpc)
   _SET_DIAGNOSTIC_(self%id_deN,e_N)
   _SET_DIAGNOSTIC_(self%id_deP,e_P)
   _SET_DIAGNOSTIC_(self%id_deC,ee_C)
   _SET_DIAGNOSTIC_(self%id_dgraz,grazing*secs_pr_day*ee_C*zooC)
   _SET_DIAGNOSTIC_(self%id_dmortz,mort_Z*secs_pr_day*zooC)
   !_SET_DIAGNOSTIC_(self%id_dexcr_N,ee_C)
   !_SET_DIAGNOSTIC_(self%id_NCP ,primprod - self%rpn*phyC)
   !_SET_DIAGNOSTIC_(self%id_PPR ,primprod*secs_pr_day)
   !_SET_DIAGNOSTIC_(self%id_NPR ,(primprod - self%rpn*phyC)*secs_pr_day)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_hzg_n2pzdq), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_

  !REVISION HISTORY:
  !Original author(s): Jorn Bruggeman

  !LOCAL VARIABLES:
  real(rk)                     :: phyC, detN
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_phyC,phyC) ! phytoplankton
   _GET_(self%id_detN,detN) ! detritus

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*(phyC+(detN/self%qmax_N)))

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Phytoplankton primary production
!
! !INTERFACE:
!  pure real(rk) function fprod(self,par,qnc,qpc)
!
! !DESCRIPTION:
! Here, the classical Michaelis-Menten formulation for nutrient uptake
! is formulated.
!
! !INPUT PARAMETERS:
!  type (type_gotm_n2pzdq), intent(in) :: self
!  real(rk), intent(in)         :: par,qnc,qpc
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  !fnp = self%gmax*(1.0_rk-self%q_min/qnc)*(par/self%I_opt*exp(1-par/self%I_opt))
!  fprod = max(self%grow_max*min((1-self%qmin_N/qnc),(1-self%qmin_P/qpc))*(1-exp(-self%affin_par*par)),0.0_rk)
!  end function fprod
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   subroutine fprod(self,par,temp_fact,qnc,qpc,primprod,Nlim,Plim,Llim)
!
! !IROUTINE: Phytoplankton primary production
! !DESCRIPTION:
! Here, the classical Michaelis-Menten formulation for nutrient uptake
! is formulated.
!
! !INPUT PARAMETERS:
   type (type_hzg_n2pzdq),    intent(in)   :: self
   real(rk), intent(in)         :: par,temp_fact,qnc,qpc
   real(rk), intent(out)        :: primprod,Nlim,Plim,Llim
   real(rk), parameter          :: secs_pr_day = 86400.0_rk

!EOP
!-----------------------------------------------------------------------
!BOC

   ! N- P- and light limitation 
   !Llim = 1-exp(-self%affin_par*par)
   !Llim = par/self%affin_par*exp(1-par/self%affin_par)
   Llim = (self%affin_par*par)/((((self%grow_max*secs_pr_day)**2.0)+(self%affin_par**2.0)*(par**2))**0.5)
   Nlim = 1-self%qmin_N/qnc
   Plim = 1-self%qmin_P/qpc

   primprod = self%grow_max*min(Nlim,Plim)*Llim*temp_fact
   !primprod = self%grow_max*min(Nlim,Plim)

   end subroutine fprod
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nutrient (N) uptake by phytoplankton
!
! !INTERFACE:
   pure real(rk) function fupN(self,DIN,qnc)
!
! !DESCRIPTION:
! 
!
! !INPUT PARAMETERS:
   type (type_hzg_n2pzdq), intent(in) :: self
   real(rk), intent(in)         :: qnc,DIN
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fupN = max(self%upmax_N-self%upmax_N/(self%qmax_N-self%qmin_N)*(qnc-self%qmin_N),0.0_rk)*DIN/(DIN+self%halfsatN)

   end function fupN
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nutrient (P) uptake by phytoplankton
!
! !INTERFACE:
   pure real(rk) function fupP(self,DIP,qpc)
!
! !DESCRIPTION:
! 
!
! !INPUT PARAMETERS:
   type (type_hzg_n2pzdq), intent(in) :: self
   real(rk), intent(in)         :: qpc,DIP
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fupP = max(self%upmax_P-self%upmax_P/(self%qmax_P-self%qmin_P)*(qpc-self%qmin_P),0.0_rk)*DIP/(DIP+self%halfsatP)

   end function fupP
!EOC

!-----------------------------------------------------------------------

!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
!  pure real(rk) function fgraz(self,phyC)
!
! !DESCRIPTION:
! Here, the classical Ivlev formulation for zooplankton grazing on
! phytoplankton is formulated.
!
! !INPUT PARAMETERS:
!  type (type_gotm_n2pzdq), intent(in) :: self
!  real(rk), intent(in)         :: phyC
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP 
!-----------------------------------------------------------------------
!BOC
!  fgraz = self%grazmax*(1.0_rk-exp(-self%iv*self%iv*phyC*phyC))

!  end function fgraz
!EOC
!-----------------------------------------------------------------------

   end module fabm_hzg_n2pzdq

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

