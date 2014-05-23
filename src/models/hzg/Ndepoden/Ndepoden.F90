#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_Ndepoden --- test for surface interfaces in FABM
!
! !INTERFACE:
   module fabm_hzg_Ndepoden
!
! !DESCRIPTION:
! This is a very simple model for pelagic consumption (eg., of det-N) & surface deposition (eg., of nut-N)
!
! !USES:
   use fabm_types
!   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
!   public type_hzg_Ndepoden, hzg_Ndepoden_init, hzg_Ndepoden_do_benthos
!
! !PUBLIC TYPES:
   type,extends(type_base_model),public    ::   type_hzg_Ndepoden
!     Variable identifiers
      type (type_bottom_state_variable_id) :: id_Ndef
      type (type_state_variable_id)        :: id_det_pel, id_nut_pel      
      type (type_diagnostic_variable_id)   :: id_peldenit
      type (type_horizontal_diagnostic_variable_id)   :: id_dNdef, id_deporate
      type (type_dependency_id)            :: id_temp !id_totNpel
      type (type_horizontal_dependency_id) :: id_depth, id_intpeldenit_dep, id_deporate_dep !id_totNpel
!     Model parameters:
      real(rk) :: Ntot0,Ndef0,const_det,rq10,T_ref,PON_denit,denit
      logical  :: use_det,NdepodenitON,Ndef_dyn
      
      !     Model procedures
      contains
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: do_surface
   end type type_hzg_Ndepoden
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
! !IROUTINE: Initialise the surface predator model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the examples\_surface\_predator namelist is read and the variables
!  exported by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hzg_Ndepoden), intent(inout), target    :: self
   integer,                          intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   character(len=64)         :: pelagic_detritus_variable='',pelagic_nutrient_variable=''
   logical                   :: NdepodenitON=.true., Ndef_dyn=.true.
   real(rk)                  :: Ntot0=100.,Ndef0=5.,const_det=10.,rq10=2.,T_ref=288.,PON_denit=5.,denit=0.024
   
   real(rk), parameter :: secs_pr_day = 86400.
   namelist /hzg_Ndepoden/  pelagic_detritus_variable,pelagic_nutrient_variable,Ntot0,Ndef0,NdepodenitON,const_det,rq10,T_ref,PON_denit,denit
                                    
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(configunit,nml=hzg_Ndepoden,err=99,end=100)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day and are converted here to values per second.
   call self%get_parameter(self%NdepodenitON,   'NdepodenitON',  default=NdepodenitON)
   call self%get_parameter(self%Ndef_dyn,        'Ndef_dyn',    default=Ndef_dyn)
   call self%get_parameter(self%Ntot0,          'Ntot0',        default=Ntot0)
   call self%get_parameter(self%Ndef0,          'Ndef0',        default=Ndef0)
   call self%get_parameter(self%const_det,      'const_det',     default=const_det)
   call self%get_parameter(self%rq10,           'rq10',          default=rq10)
   call self%get_parameter(self%T_ref,          'T_ref',          default=T_ref)
   call self%get_parameter(self%PON_denit,       'PON_denit',      default=PON_denit)
   call self%get_parameter(self%denit,          'denit',         default=denit, scale_factor=1.0_rk/secs_pr_day)
       
   ! Register link to external pelagic detritus  pool, if specified
   self%use_det = pelagic_detritus_variable/=''
   if (self%use_det) call self%register_state_dependency(self%id_det_pel,pelagic_detritus_variable)    
   if (self%NdepodenitON) call self%register_state_dependency(self%id_nut_pel,pelagic_nutrient_variable)
   
   
   ! Register state variables
   !? NOTE the benthic=.true. argument, which specifies the variable is benthic.
   call self%register_bottom_state_variable(self%id_Ndef,'Ndef','mmol/m**2','N deficit', &
                                          Ndef0,minimum=0.0_rk)
   
   !Register diagnostic variables:
   
   call self%register_horizontal_diagnostic_variable (self%id_dNdef, 'change in N deficit','mmol/m**3/d', 'change in N deficit', &
                                          output=output_time_step_averaged) 
   
   call self%register_horizontal_diagnostic_variable (self%id_deporate, 'deporate','mmol/m**3/d', 'surface deposition rate', &
                                          output=output_time_step_averaged) 
                                          
   call self%register_diagnostic_variable(self%id_peldenit,   'denitratepel','mmol/m**3/d', 'bulk pelagic denitrification rate', &
                                          output=output_time_step_averaged) 
   
                                          
   !dependencies
   ! call self%register_dependency(self%id_totNpel,standard_variables%total_nitrogen)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_depth,standard_variables%bottom_depth) !for implementing the nut_influx
   
   ! in order to be able to calculate the dNdef=total(denitrate)-deporate, we need to 
   call self%register_horizontal_dependency(self%id_deporate_dep,'hzg_Ndepoden_deporate') !,required=.false.
   !register the bulk pelagic denitrification as an integrated dependency
   !retreiving the integrated values of dependencies is not supported yet
   !call self%register_horizontal_dependency(self%id_intpeldenit_dep,vertical_integral(self%id_peldenit)) !,required=.false.vertical_integral(self%id_peldenit)
     
     
     !idNdef: horizontal state variable
     !_GET_HORIZONTAL_(self%id_Ndef,Ndef) 
     ! difference give the rate of change in Ndef
     !_SET_ODE_BEN_(self%id_Ndef, intdenitrate - deporate)
     
     !id_dNdef: horizontal diag
     !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_dNdef,(intdenitrate - deporate)*secs_pr_day)
     
     !id_det_pel: pelagic bulk state dependency
     !_GET_(self%id_det_pel,det_pel) 
     !_SET_ODE_(self%id_det_pel,-rhs_detpel)
     !_SET_SURFACE_EXCHANGE_(self%id_nut_pel,rhs_nutpel)
     
     !id_deporate: horizontal diag (+dependency?)
     !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_deporate,rhs_nutpel)
     !_GET_HORIZONTAL_(self%id_deporate, deporate)
     
     !id_peldenit: bulk diag 
     !_SET_DIAGNOSTIC_(self%id_peldenit,-rhs_detpel)
     ! id_intpeldenit: hor dependency (as water column integrated denitrification rate)
     !_GET_HORIZONTAL_(self%id_intpeldenit,intdenitrate)

   return

99 call self%fatal_error('hzg_Ndepoden_init','Error reading namelist hzg_Ndepoden')
100 call self%fatal_error('hzg_Ndepoden_init','Namelist hzg_Ndepoden was not found')

   end subroutine initialize
!EOC



!-----------------------------------------------------------------------
! Nitrogen defficiency, Ndef calculated as a 0-D dynamic variable 
! 
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   !type (type_hzg_Ndepoden),       intent(in) :: self
   class (type_hzg_Ndepoden),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: deporate, intdenitrate
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_
    

   ! turnover of long-term N-reservoir (denitrification - N-deposition)
   if (self%NdepodenitON) then
   
     ! deposition rate
     _GET_HORIZONTAL_(self%id_deporate_dep, deporate)
     
     ! water column integrated denitrification rate
     !retreiving the integrated values of dependencies is not supported yet
     !_GET_HORIZONTAL_(self%id_intpeldenit_dep,intdenitrate)
     intdenitrate=0.0_rk
     
     ! difference give the rate of change in Ndef
     _SET_ODE_BEN_(self%id_Ndef, intdenitrate - deporate)

   end if 


   ! Export diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dNdef,(intdenitrate - deporate)*secs_pr_day)

   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!-----------------------------------------------------------------------   

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: atmospheric deposition through pelagic-surface exchanges 
!
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !INPUT PARAMETERS:
   class (type_hzg_Ndepoden),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   real(rk)                   :: Ndef,totN,det_pel,rhs_nutpel,f_T,temp
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_
    
    ! if needed, calculate the N-deposition rate
    if (self%NdepodenitON) then
     
     !find the N deficit
     if (self%Ndef_dyn) then
       ! as  a dynamicaly calculated variable
       _GET_HORIZONTAL_(self%id_Ndef,Ndef) ! predator density - benthic
     else 
       ! as the difference between the initial state and the current state
       !_GET_(self%id_totNpel, totN)  ! water temperature
       !Ndef=self%totNini-totN
     end if
     
     !use the dynamic pelagic detritus concentrations or some constant
     if (self%use_det) then
       _GET_(self%id_det_pel,det_pel)      ! detritus concentration - pelagic
     else
       det_pel = self%const_det            ! no coupling - constant pelagic detritus concentration
     end if
     
     _GET_(self%id_temp, temp) !is this surface temperature? do we really need the surface temperature?
     !temp=15.0d0
     f_T=f_temp(self,temp + 273.d0)
     
     rhs_nutpel=deposit(self,det_pel,f_T,Ndef)     
     
   else
      rhs_nutpel = 0.0_rk
   endif 


   ! Set the surface flux of the pelagic nutrient variable
   if (self%use_det) then
     _SET_SURFACE_EXCHANGE_(self%id_nut_pel,rhs_nutpel)
   end if 

   ! save deposition rate as a diagnostic
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_deporate,rhs_nutpel)


   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC


!-----------------------------------------------------------------------
! !IROUTINE: pelagic denitrification through right hand sides of pelagic model

! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_hzg_Ndepoden), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: depth_of_pelagic
   real(rk)                   :: det_pel,rhs_detpel,f_T, temp
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

    !to check whether it loops over the spatial domain
    !_GET_HORIZONTAL_(self%id_depth,depth_of_pelagic)  
    !print *, "depth = ", depth_of_pelagic         
    
    if (self%NdepodenitON) then
      
      !use the dynamic pelagic detritus concentrations or some constant
      if (self%use_det) then
        _GET_(self%id_det_pel,det_pel)      ! detritus concentration - pelagic
      else
        det_pel = self%const_det            ! no coupling - constant pelagic detritus concentration
      end if
      
      _GET_(self%id_temp, temp)
      f_T=f_temp(self,temp + 273.d0)
      
      rhs_detpel=denitrify(self,det_pel,f_T)
    else
      rhs_detpel = 0.0_rk
    endif 
    
    !set the exchange rate 
    _SET_ODE_(self%id_det_pel,-rhs_detpel)
    
    !save denitrification as a pelagic variable
    _SET_DIAGNOSTIC_(self%id_peldenit,-rhs_detpel)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
pure real(rk) function denitrify(self,det_pel,f_T)
! pelagic N-loss by denitrification, emulating benthic pool and suboxic micro-environments

! !INPUT PARAMETERS:
   class (type_hzg_Ndepoden), intent(in) :: self
   real(rk), intent(in) :: det_pel,f_T
   
   !original code from Kai:
   !denitrate = self%denit * 4 * sens%f_T * (1.0d0 - exp(-det%N/self%PON_denit)) * det%N

   denitrify = self%denit * 4 * f_T * (1.0d0 - exp(-det_pel/self%PON_denit)) * det_pel
 end function denitrify
 !-----------------------------------------------------------------------
 
 
 !-----------------------------------------------------------------------
 pure real(rk) function deposit(self,det_pel,f_T,Ndef)
 ! pelagic deposition at the surface, slowly refueling N-losses 

! !INPUT PARAMETERS:
   class (type_hzg_Ndepoden), intent(in) :: self
   real(rk),intent(in) :: det_pel, f_T, Ndef
   
   !original code from Kai:
   !deporate  = self%denit * exp(-4*sens%f_T) * env%RNit

    deposit = self%denit * exp(-4*f_T)*Ndef
    
 end function deposit
 !-----------------------------------------------------------------------
 
 
 !-----------------------------------------------------------------------
 pure real(rk) function f_temp(self,T_Kelv)
 !Temperature response function
 
! !INPUT PARAMETERS:
   class (type_hzg_Ndepoden), intent(in) :: self
   real(rk),intent(in)                       :: T_Kelv
   real(rk) :: f_T
      
   !original code from calc_sensitivities in maecs_functions.f90:
   !sens%f_T  = exp( maecs%rq10*(T_Kelv-maecs%T_ref))
   f_T= exp( self%rq10*(T_Kelv-self%T_ref))

 end function f_temp
!-----------------------------------------------------------------------

 end module fabm_hzg_Ndepoden
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------