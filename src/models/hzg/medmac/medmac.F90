#include "fabm_driver.h"

!todo: DO=f(det)
!todo: Pon, Non
!todo: OXb: refractory
   
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: hzg_medmac
!
! !INTERFACE:
   module hzg_medmac
!
! !DESCRIPTION:
! This is a simple model of early diagenesis for macronutrients.
!
! !USES:
   use fabm_types
!   use fabm_expressions
!   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
!   public type_hzg_medmac, hzg_medmac_init, hzg_medmac_do_benthos
!
! !PUBLIC TYPES:
   type,extends(type_base_model),public    ::   type_hzg_medmac
!     Variable identifiers
      type (type_bottom_state_variable_id) :: id_DINb,id_ONbl !nitrogen forms in benthos
      type (type_bottom_state_variable_id) :: id_DIPb,id_OPbl,id_sorpP !phosphorus forms in benthos

!     Fabm variables to couple
      type (type_state_variable_id)        :: id_DINp,id_PONp !nitrogen forms in pelagic
      type (type_state_variable_id)        :: id_DIPp,id_POPp !phosphorus forms in pelagic

!     External dependencies
      type (type_dependency_id)            :: id_temp
      
!     Diagnostic variables

      type (type_horizontal_diagnostic_variable_id)   :: id_f_do
      type (type_horizontal_diagnostic_variable_id)   :: id_dif_N,id_adv_N,id_rem_N,id_rhs_ONbl,id_rhs_DINb
      type (type_horizontal_diagnostic_variable_id)   :: id_dif_P,id_adv_P,id_rem_P,id_rhs_OPbl,id_rhs_DIPb
      type (type_horizontal_diagnostic_variable_id)   :: id_denit_lim,id_denit_rate                                                
      type (type_horizontal_diagnostic_variable_id)   :: id_sorp_rate,id_desorp_rate,id_sorpdesorp_rate
      
      !type (type_horizontal_dependency_id) :: id_depth
      !type (type_diagnostic_variable_id)   :: id_ninflux
!     Model parameters: maximum grazing rate, half-saturation prey density, loss rate
      real(rk) :: v_d,depth_ben,r_Q10,temp_ref
      integer  :: Rmeth_N,Rmeth_P
      logical  :: use_Temp,use_DINp,use_PONp,use_DIPp,use_POPp,do_denit,do_Psorp
      real(rk) :: DINp_presc,PONp_presc,DswN,rN,kN,K_ondo,K_denit
      real(rk) :: DIPp_presc,POPp_presc,DswP,rP,kP,Rsorp
      
      !     Model procedures
      contains
      procedure :: initialize
      !procedure :: do !to implement nut-reflux to counter nut-loss
      procedure :: do_bottom      
   end type type_hzg_medmac
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
! !IROUTINE: Initialise the medmac model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the examples\_benthic\_predator namelist is read and the variables
!  exported by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hzg_medmac), intent(inout), target :: self
   integer,                          intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Onur Kerimoglu
!
! !LOCAL VARIABLES:
   character(len=64) :: DINp_variable='',PONp_variable='',DIPp_variable='',POPp_variable=''
   real(rk) :: DINb0=16.0,ONbl0=16.0,DIPb0=1.0,OPbl0=1.0,sorpP0=1.0
   real(rk) :: v_d=0.5,depth_ben=0.1,r_Q10=2.0, temp_ref=10.0
   real(rk) :: DswN=1e-5,DINp_presc=16.0,PONp_presc=16.0,rN=0.05,kN=8.0,K_ondo=1000.0,K_denit=30.0 !kN 
   real(rk) :: DswP=1e-5,DIPp_presc=1.0,POPp_presc=1.0,rP=0.05,kP=0.5,Rsorp=0.5 
   integer  :: Rmeth_N=1,Rmeth_P=1          ! remin method
   logical  :: do_denit=.False., do_Psorp=.False.
 
   real(rk), parameter :: secs_pr_day = 86400.
   namelist /hzg_medmac/  &
   v_d,depth_ben,r_Q10,temp_ref, &
   DINp_variable, PONp_variable,DINp_presc,PONp_presc,&
   DINb0, ONbl0,DswN,rN,kN,Rmeth_N,do_denit,K_denit,K_ondo, &
   DIPp_variable, POPp_variable,DIPp_presc,POPp_presc,&
   DIPb0, OPbl0,DswP,rP,kP,Rmeth_P,do_Psorp,Rsorp,sorpP0
      
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(configunit,nml=hzg_medmac,err=99,end=100)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%v_d,           'v_d',           default= v_d,        scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%depth_ben,     'depth_ben',     default=depth_ben)
   call self%get_parameter(self%r_Q10,         'r_Q10',         default=r_Q10)
   call self%get_parameter(self%temp_ref,      'temp_ref',      default=temp_ref)
   !N-pars:
   call self%get_parameter(self%DswN,          'DswN',          default=DSwn,        scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%DINp_presc,    'DINp_presc',    default=DINp_presc)
   call self%get_parameter(self%PONp_presc,    'PONp_presc',    default=PONp_presc)
   call self%get_parameter(self%rN,            'rN',            default=rN,          scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kN,            'kN',            default=kN)
   call self%get_parameter(self%Rmeth_N,       'Rmeth_N',       default=Rmeth_N)
   call self%get_parameter(self%do_denit,      'do_denit',      default=do_denit)
   call self%get_parameter(self%K_denit,       'K_denit',       default=K_denit)
   call self%get_parameter(self%K_ondo,        'K_ondo',        default=K_ondo)
   !P-pars
   call self%get_parameter(self%DswP,          'DswP',          default=DSwP,        scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%DIPp_presc,    'DIPp_presc',    default=DIPp_presc)
   call self%get_parameter(self%POPp_presc,    'POPp_presc',    default=POPp_presc)
   call self%get_parameter(self%rP,            'rP',            default=rP,          scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kP,            'kP',            default=kP)
   call self%get_parameter(self%Rmeth_P,       'Rmeth_P',       default=Rmeth_P)
   call self%get_parameter(self%do_Psorp,      'do_Psorp',      default=do_Psorp)
   call self%get_parameter(self%Rsorp,         'Rsorp',         default=Rsorp,       scale_factor=1.0_rk/secs_pr_day)
   
   ! Register state variables
   !N-variables
   call self%register_state_variable(self%id_DINb, 'DIN','mmol/m**2',&
                                     'benthic DIN', DINb0, minimum=_ZERO_)
   call self%register_state_variable(self%id_ONbl, 'ON','mmol/m**2',&
                                     'benthic ON', ONbl0, minimum=_ZERO_)
   !P-variables
   call self%register_state_variable(self%id_DIPb, 'DIP','mmol/m**2',&
                                     'benthic DIP', DIPb0, minimum=_ZERO_)
   call self%register_state_variable(self%id_OPbl, 'OP','mmol/m**2',&
                                     'benthic OP', OPbl0, minimum=_ZERO_) 
   if (do_Psorp) call self%register_state_variable(self%id_sorpP, 'sorpP','mmol/m**2',&
                                     'benthic sorped P', sorpP0, minimum=_ZERO_) 
   
   ! Register links to external temperature field
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   
   ! Register links to external pelagic detritus and mineral pools, if coupling to pelagic model
   !N dependencies
   self%use_DINp = DINp_variable/=''
   if (self%use_DINp) call self%register_state_dependency(self%id_DINp,DINp_variable)
   self%use_PONp = PONp_variable/=''
   if (self%use_PONp) call self%register_state_dependency(self%id_PONp,PONp_variable)
   !P dependencies
   self%use_DIPp = DIPp_variable/=''
   if (self%use_DIPp) call self%register_state_dependency(self%id_DIPp,DIPp_variable)
   self%use_POPp = POPp_variable/=''
   if (self%use_POPp) call self%register_state_dependency(self%id_POPp,POPp_variable)
   
   !common diags
   call self%register_diagnostic_variable(self%id_f_do,'f_do','-', &
                                          'oxygen function', output=output_time_step_averaged)
   !N diags
   call self%register_diagnostic_variable(self%id_dif_N,'difN','mmol/m**2/d', &
                                          'diffusive flux rate of DIN at soil surface', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_adv_N,'advN','mmol/m**2/d', &
                                          'advective flux rate of PON', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_rem_N,'remN','mmol/m**2/d', &
                                          'remineralization of of PON in soil', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_rhs_ONbl,'rhsONbl','mmol/m**2/d', &
                                          'rhs of ON in soil', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_rhs_DINb,'rhsDINb','mmol/m**2/d', &
                                          'rhs of DIN in soil', output=output_time_step_averaged)
   if (do_denit) then                                        
     call self%register_diagnostic_variable(self%id_denit_rate,   'denit_rate','mmol/m**2/d', & 
				          'benthic denitrification rate', output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_denit_lim,   'denit_lim','-', & 
				          'denitrification limitation in soil', output=output_time_step_averaged)
   end if 
   
   !P diags
   call self%register_diagnostic_variable(self%id_dif_P,'difP','mmol/m**2/d', &
                                          'diffusive flux rate of DIP at soil surface', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_adv_P,'advP','mmol/m**2/d', &
                                          'advective flux rate of POP', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_rem_P,'remP','mmol/m**2/d', &
                                          'remineralization of of POP in soil', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_rhs_OPbl,'rhsOPbl','mmol/m**2/d', &
                                          'rhs of OP in soil', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_rhs_DIPb,'rhsDIPb','mmol/m**2/d', &
                                          'rhs of DIP in soil', output=output_time_step_averaged)
   
   if (do_Psorp) then
     call self%register_diagnostic_variable(self%id_sorp_rate,'sorp_rate','mmol/m**2/d', &
                                          'P-sorption rate', output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_desorp_rate,'desorp_rate','mmol/m**2/d', &
                                          'P-desorption rate', output=output_time_step_averaged)                                       
     call self%register_diagnostic_variable(self%id_sorpdesorp_rate,'delta_sorpdesorp','mmol/m**2/d', &
                                          'sorption-desorption', output=output_time_step_averaged)
   end if
   return

99 call self%fatal_error('hzg_medmac_init','Error reading namelist hzg_medmac')
100 call self%fatal_error('hzg_medmac_init','Namelist hzg_medmac was not found')

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
! (matching) bottom fluxes for pelagic variables. All have units mmol/m**2/s.
!
! !INPUT PARAMETERS:
   !type (type_hzg_medmac),       intent(in) :: self
   class (type_hzg_medmac),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !REVISION HISTORY:
!  Original author(s): Onur Kerimoglu
!
! !LOCAL VARIABLES:
   real(rk)                   ::temp,f_T,f_do
   real(rk)                   ::DINp,PONp,DINb,ONbl,denit_rate,denit_lim
   real(rk)                   ::DIPp,POPp,DIPb,OPbl,sorpP,sorp_rate,desorp_rate,sorpdesorp_rate
   real(rk)                   ::advN,difN,remN,rhs_DINb,rhs_ONbl
   real(rk)                   ::advP,difP,remP,rhs_DIPb,rhs_OPbl
   
   !,sink,diffusion,remin,ddet,dnut!,nut_loss_rate
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve external dependencies
   _GET_(self%id_temp,temp)
   
   ! Retrieve pelagic state variable values.
   ! N- variables
   if (self%use_DINp) then
    _GET_(self%id_DINp,DINp)      ! pelagic concentration
   else
    DINp = self%DINp_presc         ! no coupling - constant  concentration 
   end if
   if (self%use_PONp) then
    _GET_(self%id_PONp,PONp)      ! pelagic concentration
   else
    PONp = self%PONp_presc         ! no coupling - constant  concentration 
   end if
   ! P- variables

   if (self%use_DIPp) then
    _GET_(self%id_DIPp,DIPp)      ! pelagic concentration
   else
    DIPp = self%DIPp_presc         ! no coupling - constant  concentration 
   end if
   if (self%use_POPp) then
    _GET_(self%id_POPp,POPp)      ! pelagic concentration
   else
    POPp = self%POPp_presc         ! no coupling - constant  concentration 
   end if
   
   ! Retrieve local state variable values
   ! N-
   _GET_HORIZONTAL_(self%id_ONbl,ONbl) ! detritus density - benthic
   _GET_HORIZONTAL_(self%id_DINb,DINb) ! nutrient density - benthic
   ! P-
   _GET_HORIZONTAL_(self%id_OPbl,OPbl) ! detritus density - benthic
   _GET_HORIZONTAL_(self%id_DIPb,DIPb) ! nutrient density - benthic
   
   !write(*,*)'OPbl,DIPb:',OPbl,DIPb
   
   f_T=f_temp(self,temp)
    
   ! Calculate kinetic rates:
   
   !oxygen function:->1 in high DO, ->0 in low DO
   f_do=(1-exp(-ONbl/self%K_ondo))
   
   !N
   !todo: R as a function of OCbl
   if (self%Rmeth_N .eq. 1) then
    remN=ONbl * f_T * self%rN
   else if (self%Rmeth_N .eq. 2) then
    remN=ONbl * f_T * self%rN*ONbl/(ONbl+self%kN)
   end if
   
   !denitrification
   if (self%do_denit) then
     denit_lim = f_do * DINb/(DINb+self%K_denit) 
     denit_rate= denit_lim * remN    * 6.625       * 1         * 0.116
           !             mmolN/m2/s  * molC/molN  * molO/molC * molNdenitrified/molOconsumed 
                         !(Seitzinger & Giblin,1996)  
     !denit_rate=denitrate(self,ONbl,DINb,remN)
   else 
     denit_rate=0.0
   end if
   
   !P
   if (self%Rmeth_P .eq. 1) then
    remP=OPbl * f_T * self%rP
   else if (self%Rmeth_P .eq. 2) then
    remP=OPbl * f_T * self%rP*OPbl/(OPbl+self%kP)
   end if
   
   !phosphorus sorption
   if (self%do_Psorp) then
     _GET_HORIZONTAL_(self%id_sorpP,sorpP) ! nutrient density - benthic
     sorp_rate=self%Rsorp*f_do*DIPb
     desorp_rate=self%Rsorp*(1.0-f_do)*sorpP
   else
     sorp_rate=0.0
     desorp_rate=0.0
   end if
   
   ! sinking of pelagic POM
   advN = self%v_d * PONp ! mmol/d/m² !flux should be always towards soil
   advP = self%v_d * POPp ! mmol/d/m²
   
   ! Calculate diffusive flux for particulate nutrients, the benthic variable (nut_ben. areal units) has to be converted to concentration using the depth of the benthic pool (d_ben) for being able to calculate the gradient. Then this gradient is assumed to be taking place within a distance equal to the depth of the benthic pool: 
   difN = self%DswN*(DINp-DINb/self%depth_ben)/self%depth_ben !grad=+ means flux towards soil
   ! m2/d * mmol/m3 *1/m = mmol/m2/d
   difP = self%DswP*(DIPp-DIPb/self%depth_ben)/self%depth_ben
   
   
   ! calculate rhs
   rhs_DINb=difN+remN-denit_rate
   rhs_ONbl=advN-remN
   
   rhs_DIPb=difP+remP+desorp_rate-sorp_rate
   rhs_OPbl=advP-remP
   
   ! Set local temporal derivatives of benthic variables
   _SET_ODE_BEN_(self%id_DINb,rhs_DINb)
   _SET_ODE_BEN_(self%id_ONbl,rhs_ONbl)
   _SET_ODE_BEN_(self%id_DIPb,rhs_DIPb)
   _SET_ODE_BEN_(self%id_OPbl,rhs_OPbl)
   if (self%do_Psorp) _SET_ODE_BEN_(self%id_sorpP, sorp_rate-desorp_rate)
   
   ! Set bottom fluxes of pelagic variables (these mirror local benthic derivatives), if coupled
   if (self%use_DINp) then
    _SET_BOTTOM_EXCHANGE_(self%id_DINp,-difN)
   end if 
   if (self%use_PONp) then
    _SET_BOTTOM_EXCHANGE_(self%id_PONp,-advN)
   end if
   
   if (self%use_DIPp) then
    _SET_BOTTOM_EXCHANGE_(self%id_DIPp,-difP)
   end if 
   if (self%use_POPp) then
    _SET_BOTTOM_EXCHANGE_(self%id_POPp,-advP)
   end if 

   ! Export diagnostic variables
   !common
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f_do,f_do)
   !N
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_adv_N,advN*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dif_N,difN*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rem_N,remN*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rhs_DINb, rhs_DINb*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rhs_ONbl, rhs_ONbl*secs_pr_day)
   if (self%do_denit)then
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_denit_lim,denit_lim)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_denit_rate,denit_rate*secs_pr_day) 
   end if
   !P
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_adv_P,advP*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dif_P,difP*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rem_P,remP*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rhs_DIPb, rhs_DIPb*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rhs_OPbl, rhs_OPbl*secs_pr_day)
   if (self%do_Psorp) then
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sorp_rate,sorp_rate*secs_pr_day)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_desorp_rate,desorp_rate*secs_pr_day)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sorpdesorp_rate,(sorp_rate-desorp_rate)*secs_pr_day)
   end if

   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   !EOC
   !-----------------------------------------------------------------------
   end subroutine do_bottom
 

pure real(rk) function f_temp(self,temp_celc)
! Temperature response function
 
  !INPUT PARAMETERS:
  class (type_hzg_medmac), intent(in) :: self
  real(rk),intent(in)                 :: temp_celc
      
  f_temp  = self%r_Q10**((temp_celc-self%temp_ref)/10.0)
   
end function f_temp
 
 
 
pure real(rk) function denitrate(self,ONbl,DINb,remN)
! Denitrification function

  !INPUT PARAMETERS:
  class (type_hzg_medmac), intent(in) :: self
  real(rk),intent(in)                 :: ONbl,DINb,remN
  real(rk)                            :: denit_rate, denit_lim

  denit_lim = (1.0d0 - exp(-ONbl/self%K_ondo)) * DINb/(DINb+self%K_denit)
      
  denit_rate= denit_lim * remN    * 6.625       * 1         * 0.116
                     !mmolN/m2/s  * molC/molN  * molO/molC * molNdenitrified/molOconsumed (Seitzinger & Giblin,1996)  
                       
end function denitrate



end module hzg_medmac

!-----------------------------------------------------------------------
! Copyright by HZG under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------