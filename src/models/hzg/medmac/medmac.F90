#include "fabm_driver.h"

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
      type (type_state_variable_id)        :: id_POCp !id_DICp carbon forms in pelagic

!     External dependencies
      type (type_dependency_id)            :: id_temp
      
!     Diagnostic variables

      type (type_horizontal_diagnostic_variable_id)   :: id_tempd,id_do_est
      type (type_horizontal_diagnostic_variable_id)   :: id_dif_N,id_adv_N,id_rem_N,id_rhs_ONbl,id_rhs_DINb, id_f_ondo,id_f_doin_den
      type (type_horizontal_diagnostic_variable_id)   :: id_dif_P,id_adv_P,id_rem_P,id_rhs_OPbl,id_rhs_DIPb,id_f_do_sorp
      type (type_horizontal_diagnostic_variable_id)   :: id_adv_C
      type (type_horizontal_diagnostic_variable_id)   :: id_denit_lim,id_denit_rate                                                
      type (type_horizontal_diagnostic_variable_id)   :: id_sorp_rate,id_desorp_rate,id_net_sorp_rate,id_sorpPd,id_fsorbed
      
      !type (type_horizontal_dependency_id) :: id_depth
      !type (type_diagnostic_variable_id)   :: id_ninflux
!     Model parameters: maximum grazing rate, half-saturation prey density, loss rate
      real(rk) :: v_d,depth_ben,r_Q10,temp_ref,K_on2do,K_T2do
      integer  :: Rmeth_N,Rmeth_P,dometh,den_dometh,sorpmeth
      logical  :: use_Temp,couple_pelN,couple_pelP,couple_pelC,do_denit,do_Psorp
      real(rk) :: DINp_presc,PONp_presc,DswN,rN,kN,K_denit,K_doin_den,K_ondo
      real(rk) :: DIPp_presc,POPp_presc,DswP,rP,kP,Rsorp,K_sorp,K_do_sorp,do_sorpeq
      real(rk) :: POCp_presc !DICp_presc,DswC
      !     Model procedures
      contains
      procedure :: initialize
      !procedure :: do !to implement nut-reflux to counter nut-loss
      procedure :: do_bottom      
   end type type_hzg_medmac
!
! !REVISION HISTORY:!
!  Original author(s): Onur Kerimoglu
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
   character(len=64) :: DINp_variable='',PONp_variable='',DIPp_variable='',POPp_variable='',POCp_variable=''!,DICp_variable=''
   real(rk) :: DINb0=16.0,ONbl0=16.0,DIPb0=1.0,OPbl0=1.0,sorpP0=1.0
   real(rk) :: v_d=0.5,depth_ben=0.1,r_Q10=2.0, temp_ref=10.0,K_on2do=30.0,K_T2do=15.0
   real(rk) :: DswN=1e-5,DINp_presc=16.0,PONp_presc=16.0,rN=0.05,kN=8.0,K_denit=30.0,K_doin_den=10.0,K_ondo=1000.0
   real(rk) :: DswP=1e-5,DIPp_presc=1.0,POPp_presc=1.0,rP=0.05,kP=0.5,Rsorp=0.5,K_sorp=0.5,K_do_sorp=50.0,do_sorpeq=150.0
   real(rk) :: POCp_presc=1.0 !,DICp_presc=1.0,DswC=1e-5,
   integer  :: Rmeth_N=1,Rmeth_P=1,dometh=2,den_dometh=1,sorpmeth=2
   logical  :: do_denit=.False., do_Psorp=.False., couple_pelN=.False., couple_pelP=.False., couple_pelC=.False.
 
   real(rk), parameter :: secs_pr_day = 86400.
   !namelist /hzg_medmac/  &
   !v_d,depth_ben,r_Q10,temp_ref,K_on2do,K_T2do,dometh, &
   !DINp_variable, PONp_variable,DINp_presc,PONp_presc,&
   !DINb0, ONbl0,DswN,rN,kN,Rmeth_N,do_denit,K_denit,den_dometh,K_doin_den,K_ondo, &
   !DIPp_variable, POPp_variable,DIPp_presc,POPp_presc,&
   !DIPb0, OPbl0,DswP,rP,kP,Rmeth_P,do_Psorp,sorpmeth,Rsorp,sorpP0,K_sorp,K_do_sorp,do_sorpeq

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   !if (configunit .gt. 0) read(configunit,nml=hzg_medmac,err=99,end=100)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%v_d,           'v_d',           default= v_d,        scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%depth_ben,     'depth_ben',     default=depth_ben)
   call self%get_parameter(self%r_Q10,         'r_Q10',         default=r_Q10)
   call self%get_parameter(self%temp_ref,      'temp_ref',      default=temp_ref)
   call self%get_parameter(self%K_on2do,       'K_on2do',       default=K_on2do)
   call self%get_parameter(self%K_T2do,       'K_T2do',       default=K_T2do)
   call self%get_parameter(self%dometh,       'dometh',       default=dometh)
   !N-pars:
   call self%get_parameter(self%couple_pelN,   'couple_pelN',   default=couple_pelN)
   call self%get_parameter(self%DswN,          'DswN',          default=DSwn,        scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%DINp_presc,    'DINp_presc',    default=DINp_presc)
   call self%get_parameter(self%PONp_presc,    'PONp_presc',    default=PONp_presc)
   call self%get_parameter(self%rN,            'rN',            default=rN,          scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kN,            'kN',            default=kN)
   call self%get_parameter(self%Rmeth_N,       'Rmeth_N',       default=Rmeth_N)
   call self%get_parameter(self%do_denit,      'do_denit',      default=do_denit)
   call self%get_parameter(self%K_denit,       'K_denit',       default=K_denit)
   call self%get_parameter(self%den_dometh,    'den_dometh',    default=den_dometh)
   call self%get_parameter(self%K_doin_den,    'K_doin_den',    default=K_doin_den)
   call self%get_parameter(self%K_ondo,        'K_ondo',        default=K_ondo)
   !P-pars
   call self%get_parameter(self%couple_pelP,   'couple_pelP',   default=couple_pelP)
   call self%get_parameter(self%DswP,          'DswP',          default=DSwP,        scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%DIPp_presc,    'DIPp_presc',    default=DIPp_presc)
   call self%get_parameter(self%POPp_presc,    'POPp_presc',    default=POPp_presc)
   call self%get_parameter(self%rP,            'rP',            default=rP,          scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kP,            'kP',            default=kP)
   call self%get_parameter(self%Rmeth_P,       'Rmeth_P',       default=Rmeth_P)
   call self%get_parameter(self%do_Psorp,      'do_Psorp',      default=do_Psorp)
   call self%get_parameter(self%sorpmeth,      'sorpmeth',      default=sorpmeth)
   call self%get_parameter(self%Rsorp,         'Rsorp',         default=Rsorp,       scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%K_sorp,        'K_sorp',        default=K_sorp)
   call self%get_parameter(self%K_do_sorp,     'K_do_sorp',     default=K_do_sorp)
   call self%get_parameter(self%do_sorpeq,     'do_sorpeq',     default=do_sorpeq)
   !C-pars
   call self%get_parameter(self%couple_pelC,   'couple_pelC',   default=couple_pelC)
   call self%get_parameter(self%POCp_presc,    'POCp_presc',    default=POCp_presc)
   !call self%get_parameter(self%DswC,          'DswC',          default=DSwC,        scale_factor=1.0_rk/secs_pr_day)
   !call self%get_parameter(self%DIPp_presc,    'DICp_presc',    default=DIPp_presc)


   ! Register state variables
   !N-variables
   call self%register_state_variable(self%id_DINb, 'DIN','mmol/m**2',&
                                     'benthic DIN', 0.0_rk, minimum=_ZERO_)
   call self%register_state_variable(self%id_ONbl, 'ON','mmol/m**2',&
                                     'benthic ON', 0.0_rk, minimum=_ZERO_)
   !P-variables
   call self%register_state_variable(self%id_DIPb, 'DIP','mmol/m**2',&
                                     'benthic DIP', 0.0_rk, minimum=_ZERO_)
   call self%register_state_variable(self%id_OPbl, 'OP','mmol/m**2',&
                                     'benthic OP', 0.0_rk, minimum=_ZERO_) 
   if (self%do_Psorp) then
      if (self%sorpmeth .eq. 1) then
        call self%register_state_variable(self%id_sorpP, 'sorpP','mmol/m**2',&
                                     'benthic sorbed P', 0.0_rk, minimum=_ZERO_)
      else if (self%sorpmeth .eq. 2) then
        call self%register_diagnostic_variable(self%id_sorpPd,'sorpP','mmol/m**2/d', &
                                               'benthic sorbed P', output=output_time_step_averaged)
      end if
   end if
   
   ! Register links to external temperature field
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   
   ! Register links to external pelagic detritus and mineral pools, if coupling to pelagic model is requested
   !N dependencies
   if (self%couple_pelN) then
     call self%register_state_dependency(self%id_DINp,'DINp_variable','mmol m-3','pelagic DIN')
     call self%register_state_dependency(self%id_PONp,'PONp_variable','mmol m-3','pelagic PON')
   end if

   !P dependencies
   if (self%couple_pelP) then
     call self%register_state_dependency(self%id_DIPp,'DIPp_variable','mmol m-3','pelagic DIP')
     call self%register_state_dependency(self%id_POPp,'POPp_variable','mmol m-3','pelagic POP')
   end if 
 
   !C dependencies
   if (self%couple_pelC) then
     !call self%register_state_dependency(self%id_DICp,'DICp_variable','mmol m-3','pelagic DIC')
     call self%register_state_dependency(self%id_POCp,'POCp_variable','mmol m-3','pelagic POC')
   end if

!common diags
   call self%register_diagnostic_variable(self%id_tempd,'temp','-', &
                                          'temp at soil surface', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_do_est,'do_est','-', &
                                          'DO estimated from ON', output=output_time_step_averaged)

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
   if (self%do_denit) then
     call self%register_diagnostic_variable(self%id_f_doin_den,'f_doin_den','-', &
                                          'oxygen inhibition of denitrification', output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_f_ondo,'f_ondo','-', &
                                          'oxygen inhibition mimicked by ON', output=output_time_step_averaged)
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

   !C diags
   call self%register_diagnostic_variable(self%id_adv_C,'advC','mmol/m**2/d', &
                                          'advective flux rate of POC', output=output_time_step_averaged)

   if (self%do_Psorp) then
     if (self%sorpmeth .eq. 1) then
      call self%register_diagnostic_variable(self%id_f_do_sorp,'f_do_sorp','-', &
					  'oxygen limitation(inhibition) of sorption(desorption)', output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_sorp_rate,'sorp_rate','mmol/m**2/d', &
					  'P-sorption rate', output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_desorp_rate,'desorp_rate','mmol/m**2/d', &
					  'P-desorption rate', output=output_time_step_averaged)
      call self%register_diagnostic_variable(self%id_net_sorp_rate,'net_sorp_rate','mmol/m**2/d', &
					  'sorption-desorption', output=output_time_step_averaged)
     else if (self%sorpmeth .eq. 2) then
       call self%register_diagnostic_variable(self%id_fsorbed,'fsorbed','mmol/m**2/d', &
                                              'fraction of sorbed DIP', output=output_time_step_averaged)
     end if
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
   real(rk)                   ::temp,f_T,f_doin_den,f_do_sorp,do_est,den_do_in,f_ondo
   real(rk)                   ::DINp,PONp,DINb,ONbl,denit_rate,denit_lim
   real(rk)                   ::DIPp,POPp,DIPb,OPbl,sorpP,sorp_rate,desorp_rate,net_sorp_rate,fsorbed
   real(rk)                   ::POCp !,DICp,DICb,OCbl
   real(rk)                   ::advN,difN,remN,rhs_DINb,rhs_ONbl
   real(rk)                   ::advP,difP,remP,rhs_DIPb,rhs_OPbl
   real(rk)                   ::advC !,difP,remP,rhs_DIPb,rhs_OPbl
   
   !,sink,diffusion,remin,ddet,dnut!,nut_loss_rate
   real(rk), parameter        :: secs_pr_day = 86400.
   real(rk), parameter        :: do_max = 300.0 !max. DO concentration
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve external dependencies
   _GET_(self%id_temp,temp)
   
   ! Retrieve pelagic state variable values.
   ! N- variables
   if (self%couple_pelN) then
    _GET_(self%id_DINp,DINp)      ! pelagic concentration
    _GET_(self%id_PONp,PONp)      ! pelagic concentration
   else
    DINp = self%DINp_presc         ! no coupling - constant  concentration 
    PONp = self%PONp_presc         ! no coupling - constant  concentration 
   end if
   
   ! P- variables
   if (self%couple_pelP) then
    _GET_(self%id_DIPp,DIPp)      ! pelagic concentration
    _GET_(self%id_POPp,POPp)      ! pelagic concentration
   else
    DIPp = self%DIPp_presc         ! no coupling - constant  concentration 
    POPp = self%POPp_presc         ! no coupling - constant  concentration 
   end if
   
   ! C- variables
   if (self%couple_pelC) then
    !_GET_(self%id_DICp,DICp)      ! pelagic concentration
    _GET_(self%id_POCp,POCp)      ! pelagic concentration
   else
    !DICp = self%DICp_presc         ! no coupling - constant  concentration
    POCp = self%POCp_presc         ! no coupling - constant  concentration
   end if

   ! Retrieve local state variable values
   ! N-
   _GET_HORIZONTAL_(self%id_ONbl,ONbl) ! detritus density - benthic
   _GET_HORIZONTAL_(self%id_DINb,DINb) ! nutrient density - benthic
   ! P-
   _GET_HORIZONTAL_(self%id_OPbl,OPbl) ! detritus density - benthic
   _GET_HORIZONTAL_(self%id_DIPb,DIPb) ! nutrient density - benthic
   ! C-
   !_GET_HORIZONTAL_(self%id_OCbl,OCbl) ! detritus density - benthic
   !_GET_HORIZONTAL_(self%id_DICb,DICb) ! nutrient density - benthic
   
   !write(*,*)'OPbl,DIPb:',OPbl,DIPb
   
   ! Calculate common limitation functions/variables
   f_T=f_temp(self,temp)
   
   !estimate do from on (to be used for denitrification and/or sorption)
   if (self%dometh .eq. 1) then
     do_est=do_max*exp(-ONbl/self%K_on2do)
   else if (self%dometh .eq. 2) then
     do_est=max(0.0, do_max-self%K_T2do*temp)
   end if
   ! Calculate kinetic rates:
   !N
   !todo: R as a function of OCbl
   if (self%Rmeth_N .eq. 1) then
    remN=ONbl * f_T * self%rN
   else if (self%Rmeth_N .eq. 2) then
    remN=ONbl * f_T * self%rN*ONbl/(ONbl+self%kN)
   end if

   !denitrification
   if (self%do_denit) then
     if (self%den_dometh .eq. 1) then
       !ON (as a proxy of oxygen) inhibition: lowON -> high DO -> low f_do
       f_ondo=(1-exp(-ONbl/self%K_ondo))
       den_do_in= f_ondo
     else if (self%den_dometh .eq. 2) then
       !oxygen inhibition of denitrification
       f_doin_den=1.0-do_est/(do_est+self%K_doin_den)
       den_do_in = f_doin_den
     end if

     denit_lim=den_do_in*DINb/(DINb+self%K_denit)

     denit_rate= denit_lim * remN    * 6.625       * 1         * 0.116
           !             mmolN/m2/s  * molC/molN  * molO/molC * molNdenitrified/molOconsumed
                         !(Seitzinger & Giblin,1996)  
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
     if (self%sorpmeth .eq. 1) then
       _GET_HORIZONTAL_(self%id_sorpP,sorpP)
       f_do_sorp=do_est/(do_est+self%K_do_sorp)

       !i.e., higher the oxygen, higher the f_do_sorp, higher higher the sorption rate
       sorp_rate=self%Rsorp*f_do_sorp*DIPb/(self%K_sorp+DIPb)*DIPb
       !i.e., lower the oxygen, lower the f_do_sorp, lower the desorption rate
       desorp_rate=self%Rsorp*(1.0 -f_do_sorp)*sorpP/(self%K_sorp+sorpP)*sorpP
       net_sorp_rate=sorp_rate-desorp_rate
       !write(*,*)'f_do_sorp',f_do_sorp,'sorpR:',sorp_rate*secs_pr_day,'desorpR:',desorp_rate*secs_pr_day,'sorp-desorp',net_sorp_rate*secs_pr_day
       !write(*,*)'f_do_sorp',f_do_sorp,'DIP lim',DIPb/(self%K_sorp+DIPb)*DIPb,'relSR',sorp_rate/self%Rsorp,'(1-f_do_sorp)',(1.0 -f_do_sorp),'sorpP',sorpP/(self%K_sorp+sorpP)*sorpP, 'relDS', desorp_rate/self%Rsorp, 'net',(sorp_rate-desorp_rate)/self%Rsorp
     else if (self%sorpmeth .eq. 2) then
       fsorbed=1.0_rk/(1.0_rk+exp(0.05*(self%do_sorpeq-do_est)))
       sorpP=fsorbed*DIPb
       !update DIPb
       DIPb=(1.0-fsorbed)*DIPb
       net_sorp_rate=0.0
     end if
   else
     net_sorp_rate=0.0
   end if

   ! sinking of pelagic POM
   advN = self%v_d * PONp ! mmol/d/m² !flux should be always towards soil
   advP = self%v_d * POPp ! mmol/d/m²
   advC = self%v_d * POCp ! mmol/d/m²

   ! Calculate diffusive flux for particulate nutrients, the benthic variable (nut_ben. areal units) has to be converted to concentration using the depth of the benthic pool (d_ben) for being able to calculate the gradient. Then this gradient is assumed to be taking place within a distance equal to the depth of the benthic pool: 
   difN = self%DswN*(DINp-DINb/self%depth_ben)/self%depth_ben !grad=+ means flux towards soil
   ! m2/d * mmol/m3 *1/m = mmol/m2/d
   difP = self%DswP*(DIPp-DIPb/self%depth_ben)/self%depth_ben
   
   ! calculate rhs
   rhs_DINb=difN+remN-denit_rate
   rhs_ONbl=advN-remN

   rhs_DIPb=difP+remP-net_sorp_rate
   rhs_OPbl=advP-remP

   ! Set local temporal derivatives of benthic variables
   _SET_ODE_BEN_(self%id_DINb,rhs_DINb)
   _SET_ODE_BEN_(self%id_ONbl,rhs_ONbl)
   _SET_ODE_BEN_(self%id_DIPb,rhs_DIPb)
   _SET_ODE_BEN_(self%id_OPbl,rhs_OPbl)
   if (self%do_Psorp) then
     if (self%sorpmeth .eq. 1) then
       _SET_ODE_BEN_(self%id_sorpP, net_sorp_rate)
     else if (self%sorpmeth .eq. 2) then
       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sorpPd,sorpP)
     end if
   end if

   ! Set bottom fluxes of pelagic variables (these mirror local benthic derivatives), if coupled
   if (self%couple_pelN) then
    _SET_BOTTOM_EXCHANGE_(self%id_DINp,0.0-difN)
    _SET_BOTTOM_EXCHANGE_(self%id_PONp,0.0-advN)
   end if 

   if (self%couple_pelP) then
    _SET_BOTTOM_EXCHANGE_(self%id_DIPp,0.0-difP)
    _SET_BOTTOM_EXCHANGE_(self%id_POPp,0.0-advP)
   end if

   if (self%couple_pelC) then
    !_SET_BOTTOM_EXCHANGE_(self%id_DIPp,0.0-difP)
    _SET_BOTTOM_EXCHANGE_(self%id_POCp,0.0-advC)
   end if

   ! Export diagnostic variables
   !common
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tempd,temp)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_do_est,do_est)
   !N
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_adv_N,advN*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dif_N,difN*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rem_N,remN*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rhs_DINb, rhs_DINb*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rhs_ONbl, rhs_ONbl*secs_pr_day)
   if (self%do_denit)then
    if (self%den_dometh .eq. 1) then
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f_ondo,f_ondo)
    else if (self%den_dometh .eq. 2) then
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f_doin_den,f_doin_den)
    end if
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
     if (self%sorpmeth .eq. 1) then
       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f_do_sorp,f_do_sorp)
       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sorp_rate,sorp_rate*secs_pr_day)
       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_desorp_rate,desorp_rate*secs_pr_day)
       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_net_sorp_rate,(net_sorp_rate)*secs_pr_day)
     else if (self%sorpmeth .eq. 2) then
       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fsorbed,fsorbed)
     end if
   end if
   !C
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_adv_C,advC*secs_pr_day)

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

end module hzg_medmac

!-----------------------------------------------------------------------
! Copyright by HZG under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
