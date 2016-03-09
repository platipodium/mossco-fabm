#include "fabm_driver.h"
#define DEBUG 0
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_Ndepoden --- test for surface interfaces in FABM
!
! !INTERFACE:
   module hzg_Ndepoden
!
! !DESCRIPTION:
! This is a very simple model for pelagic consumption (eg., of det-N) & surface deposition (eg., of nut-N)
!
! !USES:
   use fabm_types
   use fabm_expressions
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
!      type (type_diagnostic_variable_id)   :: id_peldenit
      type (type_horizontal_diagnostic_variable_id)   :: id_dNdef, id_deporate,id_denitrate, id_denitrilim
      type (type_dependency_id)            :: id_depth, id_temp!,id_peldenit_dep !id_totNpel
      type (type_horizontal_dependency_id) :: id_deporate_dep!, id_denitrate_dep ! id_intpeldenit_dep, id_totNpel
      type (type_horizontal_dependency_id) :: id_bendetrem
!     Model parameters:
      real(rk) :: Ntot0,Ndef0,const_det,const_nut,rq10,T_ref,K_det,K_nut,denit,depocoef,depoconst
      logical  :: use_det,use_nut,NdepodenitON,Ndef_dyn
      integer  :: depometh,denitmeth
      
      !     Model procedures
      contains
      procedure :: initialize
      !procedure :: do
      procedure :: do_bottom
      procedure :: do_surface
   end type type_hzg_Ndepoden
!
! !REVISION HISTORY:!
!  Original author(s): Onur Kerimoglu
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
   character(len=64)         :: pelagic_detritus_variable='',pelagic_nutrient_variable='',benthic_detrem_variable=''
   logical                   :: NdepodenitON=.true., Ndef_dyn=.true.
   real(rk)                  :: Ntot0=100.,Ndef0=5.,const_det=10.,const_nut=10.,rq10=2.,T_ref=288.,K_det=3.,K_nut=5.,denit=0.024,depocoef=4.0,depoconst=0.1
   integer                   :: depometh=1, denitmeth=1
   real(rk), parameter :: secs_pr_day = 86400.
   namelist /hzg_Ndepoden/  pelagic_detritus_variable,pelagic_nutrient_variable,benthic_detrem_variable,Ntot0,Ndef0,NdepodenitON,const_det,const_nut,rq10,T_ref,K_det,K_nut,denit,denitmeth,depocoef,depoconst,depometh
                                    
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
   call self%get_parameter(self%K_det,       'K_det',      default=K_det)
   call self%get_parameter(self%K_nut,       'K_nut',      default=K_nut)
   call self%get_parameter(self%denit,          'denit',         default=denit, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%depocoef,       'depocoef',      default=depocoef)
   call self%get_parameter(self%depoconst,       'depoconst',      default=depoconst, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%depometh,       'depometh',      default=depometh)
   call self%get_parameter(self%denitmeth,       'denitmeth',      default=denitmeth)
   
   ! if (self%NdepodenitON) 
   ! Register link to external pelagic detritus  pool, if specified
   self%use_det = pelagic_detritus_variable/=''
   if (self%use_det) call self%register_state_dependency(self%id_det_pel,pelagic_detritus_variable)    
   
   self%use_nut = pelagic_nutrient_variable/=''
   if (self%use_nut) call self%register_state_dependency(self%id_nut_pel,pelagic_nutrient_variable)
   
   
   ! Register state variables
   !? NOTE the benthic=.true. argument, which specifies the variable is benthic.
   call self%register_bottom_state_variable(self%id_Ndef,'Ndef','mmol/m**2',& 
					  'N deficit', Ndef0,minimum=0.0_rk)
   
   
   !Register diagnostic variables:
   call self%register_horizontal_diagnostic_variable(self%id_dNdef,   'dNdef','mmol/m**2/d', & 
				  'change in N deficit', output=output_time_step_averaged)
				  
   call self%register_horizontal_diagnostic_variable (self%id_deporate, 'deporate','mmol/m**2/d', & 
				  'surface deposition rate', output=output_time_step_averaged) 
   
   call self%register_horizontal_diagnostic_variable(self%id_denitrate,   'denitrate','mmol/m**2/d', & 
				  'bottom denitrification rate', output=output_time_step_averaged) 
   
   call self%register_horizontal_diagnostic_variable(self%id_denitrilim,   'denitrilim','-', & 
				  'denitrification limitation', output=output_time_step_averaged) 
				  
   
   !call self%register_diagnostic_variable(self%id_peldenit,   'denitrate','mmol/m**3/d', 'bulk pelagic denitrification rate', &
   !                                       output=output_time_step_averaged) 
   
                                          
   !dependencies
   ! call self%register_dependency(self%id_totNpel,standard_variables%total_nitrogen)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_depth,standard_variables%depth) !for implementing the nut_influx
   ! in order to be able to calculate the dNdef=total(denitrate)-deporate, we need to 
   call self%register_horizontal_dependency(self%id_deporate_dep,'deporate') !,required=.false.
   
   if (self%denitmeth .eq. 3 .or. self%denitmeth .eq. 4) then
   call self%register_horizontal_dependency(self%id_bendetrem,benthic_detrem_variable) !,'hzg_benthic_pool01_detrem' required=.false.
   end if
   
   !call self%register_horizontal_dependency(self%id_denitrate_dep,'hzg_Ndepoden_denitrate') !,required=.false.
   
   !register the bulk pelagic denitrification as a bulk dependency
   !call self%register_dependency(self%id_peldenit_dep,'hzg_Ndepoden_denitrate') !for implementing the nut_influx
   
   !register the integrated bulk pelagic denitrification as a horizontal dependency: @TODO: doesn't work for some reason?
   !This is as instructed in the fabmlist:
   !call self%register_horizontal_dependency(self%id_intpeldenit_dep,vertical_integral(self%id_peldenit_dep)) !
   
   !This is how it feels right:
   !call !self%register_dependency(self%id_peldenit_dep,vertical_integral(self%id_peldenit_dep)) !,required=.false.  
     
   return

99 call self%fatal_error('hzg_Ndepoden_init','Error reading namelist hzg_Ndepoden')
100 call self%fatal_error('hzg_Ndepoden_init','Namelist hzg_Ndepoden was not found')

   end subroutine initialize
!EOC



!-----------------------------------------------------------------------
! Denitirifiation at the bottom layer,
! keep Nitrogen defficiency, Ndef calculated as a 0-D dynamic variable 
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   !type (type_hzg_Ndepoden),       intent(in) :: self
   class (type_hzg_Ndepoden),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: deporate, denitrate !_vertint
   real(rk)                   :: temp, f_T, det_pel,nut_pel,bendetrem,denitrilim
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_
    
   !if do Denitrification at the bottom
   if (self%NdepodenitON) then
     
     _GET_(self%id_temp, temp)
     f_T=f_temp(self,temp + 273.d0)
     
     !pelagic detritus concentration will be required for all methods
     if (self%use_det) then
       _GET_(self%id_det_pel,det_pel)      ! detritus concentration at the bottom (?)
     else
       det_pel = self%const_det            ! no coupling - constant pelagic detritus concentration
     end if
      
     
     
     if (self%denitmeth .eq. 1) then
       denitrate = self%denit * 4 * f_T * (1.0d0 - exp(-det_pel/self%K_det)) * det_pel
     else if (self%denitmeth .eq. 2) then
       !denitrate=denitrify(self,det_pel,f_T) !mmol/m2/s
       denitrate = self%denit * 4 * f_T * (1.0d0 - exp(-det_pel/self%K_det)) * det_pel*det_pel
     else if (self%denitmeth .eq. 3) then
       _GET_HORIZONTAL_(self%id_bendetrem, bendetrem) !mmolN/m2/d
       denitrate=bendetrem/secs_pr_day   *6.625       * 1         * 0.116
		 !mmolN/m2/d *d/s        * molC/molN  * molO/molC * molNdenitrified/molOconsumed (Seitzinger & Giblin,1996)
	!write(*,'(A,2(F16.14 ))')'bendetrem,denitrilim,denitrate:',bendetrem, denitrate	 
     else if (self%denitmeth .eq. 4) then
       _GET_HORIZONTAL_(self%id_bendetrem, bendetrem) !mmolN/m2/d
     if (self%use_nut) then
       _GET_(self%id_nut_pel,nut_pel)      ! detritus concentration at the bottom (?)
     else
       nut_pel = self%const_nut            ! no coupling - constant pelagic nutrient concentration
     end if

       !this is how omexdia calculates:
       !Denitrilim = (1.0_rk-oxy/(oxy+self%kinO2denit)) * NO3/(no3+self%ksNO3denit)
       !Denitrific = (self%rFast * fdet + self%rSlow * sdet)*Denitrilim*Rescale        ! Denitrification
       denitrilim = (1.0d0 - exp(-det_pel/self%K_det)) * nut_pel/(nut_pel+self%K_nut)
       denitrate= denitrilim*bendetrem/secs_pr_day   *6.625       * 1         * 0.116
       !write(*,'(A,2(F16.14))')'bendetrem,denitrilim,denitrate:',bendetrem,denitrilim, denitrate
     end if 
     
     ! Set the surface flux of the pelagic nutrient variable, if coupled
     if (self%use_det) then
        _SET_BOTTOM_EXCHANGE_(self%id_nut_pel,-denitrate)
     end if 
     
     !if (self%depofromdenit) then
     ! Calculate N-defficiency as a dynamic state variable, where Ndef/dt= denitrate-deporate
     !we know the denitrate, we need to retrieve the deporate:
     !!! THIS DOESN'T WORK PROPERLY IN GETM, AS THE DO_SURFACE IS CALLED AFTER DO_BOTTOM, 
     !!! SO THE DEPORATE IS NOT UPDATED FOR THE CURRENT LATERAL NODE.
     _GET_HORIZONTAL_(self%id_deporate_dep, deporate) !mmol/m2/d
     
     ! water column integrated denitrification rate
     !retreiving the integrated values of dependencies is not supported yet
     !_GET_HORIZONTAL_(self%id_intpeldenit_dep,peldenit_vertint) !mmol/m2/d
     !denitrate_verint=0.0_rk !mmol/m2/d
     !this should work too, but somehow it doesn't:
     !denitrate_vertint= vertical_dependency_integral(self%id_peldenit_dep)
     
     ! difference give the rate of change in Ndef
     _SET_ODE_BEN_(self%id_Ndef, denitrate - deporate/secs_pr_day)!mmol/m2/s
     !end if
   end if 

#if DEBUG   
   write (*,'(A,4(F10.5))') '(@do_bottom), det_pel,denitrate,deporate:', det_pel, denitrate*secs_pr_day, deporate 
#endif
   
   ! Export diagnostic variables
   
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_denitrilim,denitrilim) 
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_denitrate,denitrate*secs_pr_day) !mmol/m2/d
  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dNdef,denitrate*secs_pr_day - deporate) !mmol/m2/d

   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!-----------------------------------------------------------------------   

!-----------------------------------------------------------------------
!BOP
!
! atmospheric deposition through pelagic-surface exchanges 
! 
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !INPUT PARAMETERS:
   class (type_hzg_Ndepoden),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   real(rk)                   :: Ndef,totN,det_pel,deporate,f_T,temp
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_
    
    ! if needed, calculate the N-deposition rate
    if (self%NdepodenitON) then
     
     _GET_(self%id_temp, temp) !is this surface temperature? do we really need the surface temperature?
     !temp=15.0d0
     f_T=f_temp(self,temp + 273.d0)
     
     if (self%depometh .eq. 1) then
       deporate=self%depoconst
     else if (self%depometh .eq. 2) then
       !find the N deficit
       ! as  a dynamicaly calculated variable
       _GET_HORIZONTAL_(self%id_Ndef,Ndef) 
       deporate=deposit(self,f_T,Ndef)   
     else if (self%depometh .eq. 3) then
         ! as the difference between the initial state and the current state
         !_GET_(self%id_totNpel, totN)  ! water temperature
         !Ndef=self%totNini-totN  
     end if
     
       

#if DEBUG     
     write (*,'(A,4(F10.5))') '(@do_surface) temp,f_T,Ndef,depo.rate: ', temp,f_T,Ndef,deporate*secs_pr_day
#endif

     else
        deporate = 0.0_rk
     endif 

    
     ! Set the surface flux of the pelagic nutrient variable
     if (self%use_det) then
       _SET_SURFACE_EXCHANGE_(self%id_nut_pel,deporate)
     end if 

     ! save deposition rate as a diagnostic
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_deporate,deporate*secs_pr_day) !mmol/m2/d


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
   real(rk), parameter        :: secs_pr_day = 86400.
   real(rk)                   :: depth
   real(rk)                   :: det_pel,rhs_detpel,f_T, temp
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_
  
!  !applies pelagic denitrification
!     !to check whether it really loops over the spatial domain
!     _GET_(self%id_depth,depth)  !CHECK THIS        
!     
!     if (self%NdepodenitON) then
!       
!       !use the dynamic pelagic detritus concentrations or some constant
!       if (self%use_det) then
!         _GET_(self%id_det_pel,det_pel)      ! detritus concentration - pelagic
!       else
!         det_pel = self%const_det            ! no coupling - constant pelagic detritus concentration
!       end if
!       
!       _GET_(self%id_temp, temp)
!       f_T=f_temp(self,temp + 273.d0)
!       
!       rhs_detpel=denitrify(self,det_pel,f_T)

! !#if DEBUG       
!       write (*,'(A,5(F10.5))') '(@do) z,temp,f_T, det_pel,denit.rate:',depth, temp,f_T, det_pel, rhs_detpel*secs_pr_day
! !#endif

!     else
!       rhs_detpel = 0.0_rk
!     endif 
    
    !set the exchange rate 
    !_SET_ODE_(self%id_det_pel,-rhs_detpel)
    
    !save denitrification as a pelagic variable
    !_SET_DIAGNOSTIC_(self%id_denitrate,rhs_detpel*secs_pr_day) !mmol/m3/d

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
   !denitrate = self%denit * 4 * sens%f_T * (1.0d0 - exp(-det%N/self%K_det)) * det%N

   !denitrify = self%denit * 4 * f_T * (1.0d0 - exp(-det_pel/self%K_det)) * det_pel
   denitrify = self%denit * 4 * f_T * (1.0d0 - exp(-det_pel/self%K_det)) * det_pel*det_pel
 end function denitrify
 !-----------------------------------------------------------------------
 
 
 !-----------------------------------------------------------------------
 pure real(rk) function deposit(self,f_T,Ndef)
 ! pelagic deposition at the surface, slowly refueling N-losses 

! !INPUT PARAMETERS:
   class (type_hzg_Ndepoden), intent(in) :: self
   real(rk),intent(in) :: f_T, Ndef
   
   !original code from Kai:
   !deporate  = self%denit * exp(-4*sens%f_T) * env%RNit

    deposit = self%denit * exp(-self%depocoef*f_T)*Ndef
    
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
   !f_temp= exp( self%rq10*(T_Kelv-self%T_ref))
   f_temp  = self%rq10**((T_Kelv-self%T_ref)/10.0)
   
 end function f_temp
!-----------------------------------------------------------------------

 end module hzg_Ndepoden
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
