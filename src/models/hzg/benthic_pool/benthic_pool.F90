#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_benthic_predator --- test for benthic interfaces in FABM
!
! !INTERFACE:
   module fabm_hzg_benthic_pool
!
! !DESCRIPTION:
! This is a very simple model for a benthic predator, grazing according to a
! Monod/Michaelis-Menten functional response on a pelagic prey, and
! respiring/dying according to a linear loss term. Variables for the prey
! (e.g., phytoplankon or zooplankton) and the target for the losses (typically
! a detrital or mineral pool) must be provided by an external model, e.g.,
! gotm\_npzd.
!
! !USES:
   use fabm_types
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
!   public type_hzg_benthic_pool, hzg_benthic_pool_init, hzg_benthic_pool_do_benthos
!
! !PUBLIC TYPES:
   type,extends(type_base_model),public    ::   type_hzg_benthic_pool
!     Variable identifiers
      type (type_bottom_state_variable_id) :: id_det_ben,id_nut_ben
      type (type_state_variable_id)        :: id_det_pel,id_nut_pel
      type (type_diagnostic_variable_id)   :: id_dsink,id_ddiff,id_dremin,id_drhsdet,id_drhsnut,id_dloss

!     Model parameters: maximum grazing rate, half-saturation prey density, loss rate
      real(rk) :: g_max,K,h_const,remin_max,v_d,d_ben,diff,const_det,const_nut,k_remin,det_loss_max,k_loss
      logical  :: use_nut,use_det,sat_remin,det_loss
      
      !     Model procedures
      contains
      procedure :: initialize
      procedure :: do_bottom
      
   end type type_hzg_benthic_pool
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the benthic predator model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the examples\_benthic\_predator namelist is read and the variables
!  exported by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hzg_benthic_pool), intent(inout), target :: self
   integer,                          intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                  :: det_ben_initial=0.01, nut_ben_initial=0.01
   real(rk)                  :: g_max = 1., K=1., h_const=0.05,remin_max=0.5, v_d=0.5, d_ben=0.1, diff=1e-5,const_nut=10.,const_det=10.,k_remin=1.,det_loss_max=0.1,k_loss=1.
   character(len=64)         :: pelagic_nutrient_variable='',pelagic_detritus_variable=''
   logical                   :: sat_remin,det_loss
   real(rk), parameter :: secs_pr_day = 86400.
   namelist /hzg_benthic_pool/  pelagic_nutrient_variable,pelagic_detritus_variable, det_ben_initial,nut_ben_initial, &
g_max,K,h_const,remin_max,diff,v_d,d_ben,const_nut,const_det,k_remin,det_loss_max,k_loss,sat_remin,det_loss !pelagic_nutrient_variable,pelagic_detritus_variable,det_ben_initial,nut_ben_initial,g_max,K,h,diff 
                                    
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(configunit,nml=hzg_benthic_pool,err=99,end=100)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%g_max,         'g_max',         default=g_max,       scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%h_const,       'h_const',       default=h_const,     scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%remin_max,     'remin_max',     default=remin_max,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%K,             'K',             default=K)
   call self%get_parameter(self%v_d,           'v_d',           default= v_d,        scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%d_ben,         'd_ben',         default=d_ben)
   call self%get_parameter(self%diff,          'diff',          default=diff,        scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%const_nut,     'const_nut',     default=const_nut)
   call self%get_parameter(self%const_det,     'const_det',     default=const_det)
   call self%get_parameter(self%k_remin,       'k_remin',       default=k_remin)
   call self%get_parameter(self%det_loss_max,  'det_loss_max',  default=det_loss_max, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%k_loss,        'k_loss',        default= k_loss)
   call self%get_parameter(self%sat_remin,     'sat_remin',     default=sat_remin)
   call self%get_parameter(self%det_loss,      'det_loss',      default=det_loss)
   
   ! Register state variables
   ! NOTE the benthic=.true. argument, which specifies the variable is benthic.
   call self%register_state_variable(self%id_det_ben,'det_ben','mmol/m**2','det_ben', &
                                          det_ben_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_nut_ben,'nut_ben','mmol/m**2','nut_ben', &
                                          nut_ben_initial,minimum=_ZERO_)


   ! Register link to external pelagic detritus and mineral pools, if coupling to pelagic model
   self%use_nut = pelagic_nutrient_variable/=''
   self%use_det = pelagic_detritus_variable/=''
   if (self%use_nut) call self%register_state_dependency(self%id_det_pel,pelagic_detritus_variable)    
   if (self%use_det) call self%register_state_dependency(self%id_nut_pel,pelagic_nutrient_variable)
   
   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_dsink,'sinking','mmol/m**2/d',  'sinking rate',             &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_ddiff,'diffusion','',  'diffusion',                  &
                     output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_dremin,'remin_ben','',  'remin_ben',                 &
                     output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_dloss,'det_loss','',  'det_loss',                 &
                     output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_drhsdet,'ddet','',  'ddet',                          &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_drhsnut,'dnut','',  'dnut',                          &
                     output=output_time_step_averaged)
   return

99 call fatal_error('hzg_benthic_pool_init','Error reading namelist hzg_benthic_pool')
100 call fatal_error('hzg_benthic_pool_init','Namelist hzg_benthic_pool was not found')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of benthic_predator model
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
! This routine calculates the benthic sink and source terms, as well as
! (matching) bottom fluxes for pelagic variables. Both have units mmol/m**2/s.
!
! !INPUT PARAMETERS:
   !type (type_hzg_benthic_pool),       intent(in) :: self
   class (type_hzg_benthic_pool),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                   :: nut_ben,nut_pel,det_ben,det_pel,sink,diffusion,remin,ddet,dnut,det_loss_rate
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.

   if (self%use_nut) then
   _GET_(self%id_nut_pel,nut_pel)      ! nutrient concentration - pelagic
   else
   nut_pel = self%const_nut            ! no coupling - constant pelagic nutrient concentration 
   end if

   if (self%use_det) then
   _GET_(self%id_det_pel,det_pel)      ! detritus concentration - pelagic
   else
   det_pel = self%const_det            ! no coupling - constant pelagic detritus concentration
   end if

   _GET_HORIZONTAL_(self%id_det_ben,det_ben) ! detritus density - benthic
   _GET_HORIZONTAL_(self%id_nut_ben,nut_ben) ! nutrient density - benthic

   ! Calculate sinking of pelagic detritus
   sink = self%v_d*det_pel ! mmol/dm²

   ! Calculate remineralisation for benthic detritus, saturated remineralisation, if sat_remin=true
   if (self%sat_remin) then
   remin = self%remin_max*det_ben/(det_ben+self%k_remin)
   else
   remin = self%h_const  
   end if
   
   ! Calculate detritus loss, if det_loss=true
   if (self%det_loss) then
   det_loss_rate = self%det_loss_max*det_ben*det_ben/(det_ben+self%k_loss)
   else
   det_loss_rate = 0
   end if

   ! Calculate diffusive flux for particulate nutrients, has to be converted to get concentration
   diffusion = self%diff*(nut_ben/self%d_ben-nut_pel)/self%d_ben

   ! calculate rhs
   ddet = sink  - remin*det_ben - det_loss_rate
   dnut = remin*det_ben - diffusion

   ! Set local temporal derivatives of benthic variables
   _SET_ODE_BEN_(self%id_det_ben,ddet)
   _SET_ODE_BEN_(self%id_nut_ben,dnut)
   
   ! Set bottom fluxes of pelagic variables (these mirror local benthic derivatives), if coupled
   if (self%use_nut) then
   _SET_BOTTOM_EXCHANGE_(self%id_nut_pel,diffusion)
   end if 

   ! Set bottom fluxes of pelagic variables (these mirror local benthic derivatives), if coupled
   if (self%use_det) then
   _SET_BOTTOM_EXCHANGE_(self%id_det_pel,-sink)
   end if 

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dsink,sink)
   _SET_DIAGNOSTIC_(self%id_ddiff,diffusion*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_dremin,remin*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_dloss,det_loss_rate*secs_pr_day) 
   _SET_DIAGNOSTIC_(self%id_drhsdet ,ddet)
   _SET_DIAGNOSTIC_(self%id_drhsnut ,dnut)


   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

!-----------------------------------------------------------------------

   end module fabm_hzg_benthic_pool

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
