#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_benthic_predator --- test for benthic interfaces in FABM
!
! !INTERFACE:
   module hzg_benthic_pool
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
!   use fabm_driver

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
      type (type_horizontal_diagnostic_variable_id)   :: id_dsink,id_ddiff,id_dremin,id_drhsdet,id_drhsnut !,id_nloss
      !type (type_horizontal_dependency_id) :: id_depth
      !type (type_diagnostic_variable_id)   :: id_ninflux
!     Model parameters: maximum grazing rate, half-saturation prey density, loss rate
      real(rk) :: v_d,depth_ben,diff,const_det,const_nut,remin_const,remin_max,k_remin!,g_max,K,nut_loss_max,k_loss,nut_pel_influx
      character(len=16):: SVname
      logical  :: use_nut,use_det,do_sat_remin,do_nut_loss
      
      !     Model procedures
      contains
      procedure :: initialize
      !procedure :: do !to implement nut-reflux to counter nut-loss
      procedure :: do_bottom      
   end type type_hzg_benthic_pool
!
! !REVISION HISTORY:!
!  Original author(s): Lena Spruch, Kai Wirtz, Onur Kerimoglu
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
   real(rk)                  :: remin_const=0.05,remin_max=0.5, v_d=0.5, d_ben=0.1,  diff=1e-5,const_nut=10.,const_det=10.,k_remin=1.!,!g_max = 1., K=1., k_loss=1.,nut_loss_max=0.1,nut_pel_influx=0.10
   character(len=64)         :: SVname='',pelagic_nutrient_variable='',pelagic_detritus_variable=''
   logical                   :: do_sat_remin!,do_nut_loss
   real(rk), parameter :: secs_pr_day = 86400.
   namelist /hzg_benthic_pool/  SVname, pelagic_nutrient_variable,pelagic_detritus_variable, det_ben_initial,nut_ben_initial, &
diff,v_d,d_ben,const_nut,const_det,remin_const,remin_max,k_remin,do_sat_remin!,g_max,K,k_loss,do_nut_loss,nut_loss_max,nut_pel_influx 
                                    
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(configunit,nml=hzg_benthic_pool,err=99,end=100)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%SVname,        'SVname',             default=SVname)
   !call self%get_parameter(self%g_max,         'g_max',         default=g_max,       scale_factor=1.0_rk/secs_pr_day)
   !call self%get_parameter(self%K,             'K',             default=K)
   call self%get_parameter(self%v_d,           'v_d',           default= v_d,        scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%depth_ben,     'd_ben',         default=d_ben)
   call self%get_parameter(self%diff,          'diff',          default=diff,        scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%const_nut,     'const_nut',     default=const_nut)
   call self%get_parameter(self%const_det,     'const_det',     default=const_det)
   call self%get_parameter(self%remin_const,    'remin_const',       default=remin_const,     scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%remin_max,     'remin_max',     default=remin_max,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%k_remin,       'k_remin',       default=k_remin)
   call self%get_parameter(self%do_sat_remin,  'do_sat_remin',  default=do_sat_remin)
   !call self%get_parameter(self%k_loss,        'k_loss',        default=k_loss)
   !call self%get_parameter(self%nut_pel_influx,'nut_pel_influx',default=nut_pel_influx, scale_factor=1.0_rk/secs_pr_day)
   !call self%get_parameter(self%do_nut_loss,   'do_nut_loss',   default=do_nut_loss)
   !call self%get_parameter(self%nut_loss_max,  'nut_loss_max',  default=nut_loss_max, scale_factor=1.0_rk/secs_pr_day)
   
   ! Register state variables
   ! NOTE the benthic=.true. argument, which specifies the variable is benthic.
   !call self%register_state_variable(self%id_det_ben,'det_ben','mmol/m**2','det_ben', &
   call self%register_state_variable(self%id_det_ben,trim(SVname)//'_det','mmol/m**2',trim(SVname)//'-det', &
                                          det_ben_initial,minimum=_ZERO_)
   !call self%register_state_variable(self%id_nut_ben,'nut_ben','mmol/m**2','nut_ben', &
   call self%register_state_variable(self%id_nut_ben,trim(SVname)//'_nut','mmol/m**2',trim(SVname)//'-nut', &
                                          nut_ben_initial,minimum=_ZERO_)

   
   ! Register link to external pelagic detritus and mineral pools, if coupling to pelagic model
   self%use_nut = pelagic_nutrient_variable/=''
   self%use_det = pelagic_detritus_variable/=''
   if (self%use_det) call self%register_state_dependency(self%id_det_pel,pelagic_detritus_variable)    
   if (self%use_nut) call self%register_state_dependency(self%id_nut_pel,pelagic_nutrient_variable)
   !for implementing the nut_influx:
   ! call self%register_dependency(self%id_depth,standard_variables%bottom_depth) 
   
   ! Register diagnostic variables
   !id_dsink,id_ddiff,id_dremin,id_drhsdet,id_drhsnut,id_nloss
   call self%register_diagnostic_variable(self%id_dsink,'detsed','mmol/m**2/d',  'det_sed_flux',             &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_ddiff,'nutdif','mmol/m**2/d',  'nut_diff_flux',                  &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_dremin,'detrem','mmol/m**2/d',  'det_remin_rate',                 &
                     output=output_instantaneous) !output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_drhsdet,'ddet','mmol/m**2/d',  'det_RHS',                          &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_drhsnut,'dnut','mmol/m**2/d',  'nut_RHS',                          &
                     output=output_time_step_averaged)
   !call self%register_diagnostic_variable(self%id_nloss,'nutloss','mmol/m**2/d',  'nut_loss_rate',                 &
   !                  output=output_time_step_averaged)
   !call self%register_diagnostic_variable(self%id_ninflux,'ninflux','',  'ninflux',                          &
   !                  output=output_time_step_integrated)
   return

99 call self%fatal_error('hzg_benthic_pool_init','Error reading namelist hzg_benthic_pool')
100 call self%fatal_error('hzg_benthic_pool_init','Namelist hzg_benthic_pool was not found')

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
   real(rk)                   :: nut_ben,nut_pel,det_ben,det_pel,sink,diffusion,remin,ddet,dnut!,nut_loss_rate
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
   if (self%do_sat_remin) then
     remin = self%remin_max*det_ben/(det_ben+self%k_remin)
   else
     remin = self%remin_const  
   end if

   ! Calculate diffusive flux for particulate nutrients, the benthic variable (nut_ben. areal units) has to be converted to concentration using the depth of the benthic pool (d_ben) for being able to calculate the gradient. Then this gradient is assumed to be taking place within a distance equal to the depth of the benthic pool: 
   diffusion = self%diff*(nut_ben/self%depth_ben-nut_pel)/self%depth_ben
   ! m2/d * mmol/m3 *1/m = mmol/m2/d
   
   !    ! Calculate detritus loss, if do_nut_loss=true, this might be due to e.g., denitrification
!    if (self%do_nut_loss) then
!      !det_loss_rate = self%det_loss_max*det_ben*det_ben/(det_ben+self%k_loss)
!      nut_loss_rate = self%nut_loss_max*det_ben/(det_ben+self%k_loss)
!    else
!      nut_loss_rate = 0
!    end if

   ! calculate rhs
   ddet = sink  - remin*det_ben
   dnut = remin*det_ben - diffusion !- nut_loss_rate * nut_ben 

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
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dsink,sink*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_ddiff,diffusion*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dremin,remin*det_ben*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_drhsdet, ddet*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_drhsnut, dnut*secs_pr_day)
   !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_nloss,nut_loss_rate*nut_ben*secs_pr_day) 

   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC


! nut-reflux
! !-----------------------------------------------------------------------
! !BOP
! !
! ! !IROUTINE: Right hand sides of NPZD model
! !
! ! !INTERFACE:
!    subroutine do(self,_ARGUMENTS_DO_)
! !
! ! !INPUT PARAMETERS:
!    class (type_hzg_benthic_pool), intent(in) :: self
!    _DECLARE_ARGUMENTS_DO_
! !
! ! !LOCAL VARIABLES:
!    real(rk)                   :: depth_of_pelagic
! !EOP
! !-----------------------------------------------------------------------
! !BOC
!    ! Enter spatial loops (if any)
!    _LOOP_BEGIN_
! 
!    !in order to avoid drying out of nutrients from the system due to this loss rate,
!    !apply a constant rate of nutrient addition to the pelagic (eg., lateral inputs, nitrogen deposition, etc)
!    if (self%use_nut .AND. self%do_nut_loss) then
!    
!     _GET_HORIZONTAL_(self%id_depth,depth_of_pelagic)  
!     !print *, "depth = ", depth_of_pelagic 
!     _SET_ODE_(self%id_nut_pel, self%nut_pel_influx/depth_of_pelagic)
!     
!     !_SET_DIAGNOSTIC_(self%id_ninflux,self%nut_pel_influx/depth_of_pelagic)
!    end if
! 
!    ! Leave spatial loops (if any)
!    _LOOP_END_
! 
!    end subroutine do
! !EOC


!-----------------------------------------------------------------------

   end module hzg_benthic_pool

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
