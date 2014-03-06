#include "fabm_driver.h"

!> @brief This is the main routine where right-hand-sides are calculated
!> @details
!> NOTE: Although this subroutine looks as if it's not a part of any module,
!! it is temporarily included in the fabm_hzg_maecs module (inside maecs.F90)
!! when compiling the documentation, such that the subroutine is documented 
!! under the 'Data Type Documentation' chapter, where the in-body docs are listed
!>
!> **Phytoplankton Equations** 
!> \n We distinguish between mass state variables 
!! (in units of carbon, nitrogen, & phosphorus) and property state variables. 
!! \latexonly For a textual narration and equations, see: Section \ref{sec:ModStr} \endlatexonly \n
!>  Current 'traits' are: 
!> - nitrogen allocated to rubisco (frac_Rub) 
!> - Chla content of chloroplasts  (theta) 
!>
!> **General code structure:**
!> 1. Calculation of quotas
!> 2. Calculation of fluxes and mass exchange rates & Specify rates of change of traits variables 
!> 3. Assign mass exchange rates ('rhs(j,i)')
!> 4. Assign rates of change of 'traits' property variables
!>
!> \n **Detailed Descriptions:**
! @todo: althougth HIDE_IN_BODY_DOCS=NO, the body-documentation is not included! HAS TO BE FIXED
! @todo: add equations
subroutine maecs_do(self,_ARGUMENTS_DO_)

use fabm_types
use maecs_types
use maecs_functions
use maecs_primprod 
use maecs_grazing

! !INPUT PARAMETERS:
 class (type_hzg_maecs),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
type (type_maecs_rhs)    :: rhsv
type (type_maecs_phy)    :: phy   ! phytoplankton type containing state and trait information
type (type_maecs_zoo)    :: zoo   ! zooplankton type 
type (type_maecs_om)     :: dom, det, nut, uptake, exud, lossZ, floppZ, nquot 
type (type_maecs_env)    :: env
type (type_maecs_switch) :: mswitch    
type (type_maecs_traitdyn)::acclim
type (type_maecs_sensitivities) :: sens

! --- LOCAL MODEL VARIABLES:
integer  :: i, j, iz, ihour, iloop
real(rk) :: reminT, degradT       ! Temp dependent remineralisation and hydrolysis rates
! --- QUOTA and FRACTIONS
real(rk) :: phys_status, dQN_dt, dRchl_phyC_dt=0.0_rk ! []

! --- ZOOPLANKTON GRAZING, INGESTION, MORTALITY, RESPIRATION... 
real(rk) :: graz_rate   ! carbon-specific grazing rate                          [d^{-1}]
real(rk) :: zoo_respC   ! temperature dependent carbon respiration rate  [mmolC m^{-3} d^{-1}]
real(rk) :: zoo_mort
real(rk) :: secs_pr_day = 86400.0_rk
! --- AGGREGATION 
real(rk) :: aggreg_rate ! aggregation among phytoplankton and between phytoplankton & detritus [d^{-1}]    
logical  :: out = .true.
!   if(36000.eq.secondsofday .and. mod(julianday,1).eq.0 .and. outn) out=.true.    
#define _KAI_ 0
#define _MARKUS_ 1
! #define UNIT / 86400
#define UNIT *1.1574074074E-5_rk

 _LOOP_BEGIN_
! First retrieve current (local) state  variable values
!#S_GET
!---------- GET for each state variable ----------
  _GET_(self%id_nutN, nut%N)  ! Dissolved Inorganic Nitrogen DIN in mmol-N/m**3
  _GET_(self%id_phyC, phy%C)  ! Phytplankton Carbon in mmol-C/m**3
  _GET_(self%id_phyN, phy%N)  ! Phytplankton Nitrogen in mmol-N/m**3
  _GET_(self%id_detC, det%C)  ! Detritus Carbon in mmol-C/m**3
  _GET_(self%id_detN, det%N)  ! Detritus Nitrogen in mmol-N/m**3
  _GET_(self%id_domC, dom%C)  ! Dissolved Organic Carbon in mmol-C/m**3
  _GET_(self%id_domN, dom%N)  ! Dissolved Organic Nitrogen in mmol-N/m**3
if (self%RubiscoOn) then
      _GET_(self%id_Rub, phy%Rub)  ! fraction of Rubisco in -
end if
if (self%PhotoacclimOn) then
      _GET_(self%id_chl, phy%chl)  ! Chl:C ratio in mg-Chla/mmol-C
end if
if (self%PhosphorusOn) then
      _GET_(self%id_nutP, nut%P)  ! Dissolved Inorganic Phosphorus DIP in mmol-P/m**3
      _GET_(self%id_phyP, phy%P)  ! Phytplankton Phosphorus in mmol-P/m**3
      _GET_(self%id_detP, det%P)  ! Detritus Phosphorus in mmol-P/m**3
      _GET_(self%id_domP, dom%P)  ! Dissolved Organic Phosphorus in mmol-P/m**3
end if
if (self%SiliconOn) then
      _GET_(self%id_nutS, nut%Si)  ! Dissolved Inorganic Silicon Si in mmol-Si/m**3
      _GET_(self%id_phyS, phy%Si)  ! Phytplankton Silicon in mmol-Si/m**3
      _GET_(self%id_detS, det%Si)  ! Detritus Silicon in mmol-Si/m**3
end if
if (self%GrazingOn) then
      _GET_(self%id_zooC, zoo%C)  ! Zooplankton Carbon in mmol-C/m**3
end if
!#E_GET
if (.not. self%GrazingOn) then
      zoo%C = self%zooC_initial
end if

!phy%Rub = Rub
!phy%chl = chl
! Retrieve current environmental conditions.
!#S_GED
  _GET_(self%id_temp, env%temp)  ! water temperature
  _GET_(self%id_par, env%par)  ! light photosynthetically active radiation
!#E_GED


! @ingroup main
!> @fn fabm_hzg_maecs::maecs_do () 
!> 1. Calculation of quotas
!>   - call min_mass with method=2, store phy\%C and \%N in phy\%reg
!>   - call calc_internal_states
!>   - if PhotoacclimOn=.false., calculate: 
!>     - phy\%chl=phy\%C * self\%frac_chl_ini 
!>     - phy\%frac\%theta = self\%frac_chl_ini * self\%itheta_max
!>     - phy\%theta= self\%frac_chl_ini / (self\%frac_Rub_ini * phy\%relQ\%N**self\%sigma)
!>   - call calc_sensitivities
!> @todo: min_mass correction of phy%\C and phy\%N at this stage requires specification of threshold values. What about back-calculating phy\%reg\%N from the smooth_small corrected phy\%Q\%N?

! --- checking and correcting extremely low state values  ------------  
call min_mass(self,phy,method=2) !_KAI_ minimal reasonable Phy-C and -Nitrogen

! --- stoichiometry of autotrophs (calculating QN_phy, frac_R, theta, and QP_phy)
call calc_internal_states(self,phy,det,dom,zoo)
!write (*,'(A,2(F10.3))') '2 P,chl=',phy%C,phy%chl

if (.not. self%PhotoacclimOn) then  
   phy%chl         = phy%C * self%frac_chl_ini   ! total Chl mg-CHL/m3
   phy%frac%theta  = self%frac_chl_ini * self%itheta_max
   phy%theta       = self%frac_chl_ini /(self%frac_Rub_ini*phy%relQ%N**self%sigma)
! g-CHL/mol-C*m3
! write (*,'(A,2(F10.3))') 'theta:',phy%relQ%N**self%sigma,phy%theta   
end if 

call calc_sensitivities(self,sens,phy,env,nut)


!> @fn fabm_hzg_maecs::maecs_do ()
!> 2. Calculation of fluxes and mass exchange rates 
!! & Specify rates of change of traits variables
!>   - call photosynthesis
!>   - if GrazingOn: 
!>     - call grazing
!>     - call grazing_losses
!>     - calculate graz_rate and zoo_mort
!>   - calc. aggreg_rate
!>   - calc. degradT and reminT

! --- ALGAL GROWTH and EXUDATION RATES, physiological trait dynamics ----------------
call photosynthesis(self,sens,phy,uptake,exud,acclim)


! ----------------       grazing        -------------------------------------------
if (self%GrazingOn) then
  call grazing(self%g_max * sens%f_T,self%k_grazC,phy%C,graz_rate)
  zoo%feeding = graz_rate
  zoo_respC   = self%basal_resp_zoo * sens%f_T  !  basal respiration of grazers
  nquot       = type_maecs_om(1.0_rk, phy%Q%N, phy%Q%P, phy%Q%Si )
  mswitch     = type_maecs_switch(self%PhosphorusOn,self%SiliconOn,.true. )
! --- calculates zooplankton loss rates (excretion->Nut, floppy+egestion->Det), specific to C
  call grazing_losses(zoo,zoo_respC,nquot,lossZ,floppZ, mswitch) 
!  --- transform from specific to bulk grazing rate
  graz_rate   = graz_rate * zoo%C 
!  --- quadratic closure term
  zoo_mort    = self%mort_zoo * sens%f_T  * zoo%C

else
  graz_rate   = 0.0_rk
  zoo_mort    = 0.0_rk
  lossZ       = type_maecs_om(0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk)
  floppZ      = type_maecs_om(0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk)
end if

! --- phytoplankton aggregation -------------------------------------------------
! If biovolume is primarily determined by the nitrogen content, also for detritus
!aggreg_rate = self%phi_agg * dom%C * (phy%N + det%N)                    ! [d^{-1}] 

!_GET_(self%id_fracR,phys_status )  
!write (*,'(A,1(F10.3))') 'phys=',phys_status

aggreg_rate = self%phi_agg * (1.0_rk - exp(-0.02*dom%C)) * (phy%N + det%N) 
!         vS * exp(-4*phys_status )                ! [d^{-1}] 
!aggreg_rate = aggreg_rate * exp(-4*phy%rel_phys ) 

! ------------------------------------------------------------------
!  ---  hydrolysis & remineralisation rate (temp dependent)
degradT     = self%hydrol * sens%f_T
reminT      = self%remin  * sens%f_T


!> @fn fabm_hzg_maecs::maecs_do ()
!> 3. Assign mass exchange rates ('rhs(j,i)')
!>   - phyC= uptake - dil - exud - aggreg\_rate - graz\_rate
!> @todo: add the rhs equations


! right hand side of ODE (rhs)    
!  #define UNIT /self%secs_pr_day
!__________________________________________________________________________
!
! PHYTOPLANKTON C
rhsv%phyC = uptake%C              * phy%C &
           - self%dil             * phy%C &
           - exud%C               * phy%C &  !TODO: move loss rates to mu, also checking for 
           - aggreg_rate          * phy%C &  !      trait dependencies
           - graz_rate                    

! write (*,'(A,3(F9.4))') 'flxc=',uptake%C, phy%C,nut%P 
! write (*,'(A,3(F9.4))') 'c=',phy%chl,phy%Rub, phy%C

!_____________________________________________________________________________
!
! PHYTOPLANKTON N
rhsv%phyN =  uptake%N             * phy%C &
           - exud%N               * phy%C & 
           - aggreg_rate          * phy%N &
           - self%dil             * phy%N &          
           - graz_rate * phy%Q%N       
   
!rhsv%phyN = 0.0_rk

!_____________________________________________________________________________

!> @fn fabm_hzg_maecs::maecs_do ()
!> 4. Assign rates of change of 'traits' property variables
!>    - if PhotoacclimOn: rhsv%chl
!>    - if RubiscoOn: rhsv%Rub

if (self%PhotoacclimOn) then
! PHYTOPLANKTON CHLa
     ! note that theta*rel_chloropl in units [mg Chla (mmol C)^{-1}] 
   dQN_dt        = (rhsv%phyN * phy%reg%C - rhsv%phyC * phy%reg%N) / (phy%reg%C*phy%reg%C)
! TODO: dangerous to work with RHS instead of net uptake rates (mortality has no physiological effect)

   dRchl_phyC_dt =  acclim%dRchl_dtheta * acclim%dtheta_dt   & 
                  + acclim%dRchl_dfracR * acclim%dfracR_dt   & 
                  + acclim%dRchl_dQN    * dQN_dt 

   rhsv%chl = phy%theta * phy%frac%Rub * phy%relQ%N**self%sigma * rhsv%phyC + dRchl_phyC_dt * phy%C

!write (*,'(A,4(F10.3))') 'rhs chl=', phy%theta * phy%frac%Rub * phy%relQ%N**self%sigma * rhsv%phyC,dRchl_phyC_dt * phy%reg%C*1E1,phy%relQ%N**self%sigma,phy%theta

!_____________________________________________ _________________________________
end if 

if (self%RubiscoOn) then 
   rhsv%Rub  = acclim%dfracR_dt * phy%C + phy%Rub/phy%reg%C * rhsv%phyC 
end if 

!________________________________________________________________________________
!
! ZOOPLANKTON zoo%feeding
if (self%GrazingOn) then  
   rhsv%zooC   =  zoo%yield * graz_rate       &
                - zoo_mort          * zoo%C   &
                - self%dil          * zoo%C   &         
                - lossZ%C           * zoo%C
else
   rhsv%zooC      = 0.0_rk
end if 
!________________________________________________________________________________
!
!  --- DETRITUS C
rhsv%detC   =  floppZ%C             * zoo%C   &
             + aggreg_rate          * phy%C   &
             + zoo_mort             * zoo%C   &
             - self%dil             * det%C   &             
             - degradT              * det%C                

!________________________________________________________________________________
!
!  --- DETRITUS N
rhsv%detN   = floppZ%N              * zoo%C   &
             + aggreg_rate          * phy%N   &
             - self%dil             * det%N   &
             + zoo_mort             * zoo%N   & 
             - degradT              * det%N 
!________________________________________________________________________________
!
!  --- DOC
rhsv%domC   = exud%C                * phy%C   & 
             + degradT              * det%C   &
             - self%dil             * dom%C   &        
             - reminT               * dom%C 
!________________________________________________________________________________
!
!  --- DON
rhsv%domN   = exud%N                * phy%C   &
             + degradT              * det%N   &
             - self%dil             * dom%N   &         
             - reminT               * dom%N
!________________________________________________________________________________
!
! DIC
!if (self%BioCarbochemOn) then
!  rhsv%dic     = -uptake%grossC      * phy%C   &
!                + uptake%lossC       * phy%C   &
!                + reminT             * dom%C   & 
!                + lossZ%C            * zoo%C
!
!_SET_ODE_(self%id_dic,rhsv%dic UNIT)
!end if
!________________________________________________________________________________
!
!  --- DIN
rhsv%nutN   = -uptake%N            * phy%C    &
             + reminT              * dom%N    &
             + lossZ%N             * zoo%C    &
             + self%dil * (self%nutN_initial - nut%N)
!________________________________________________________________________________
!
if (self%PhosphorusOn) then 
  ! ---  PHYTOPLANKTON P
   rhsv%phyP = uptake%P              * phy%C    & 
              - exud%P               * phy%C    & 
              - self%dil             * phy%P    &            
              - aggreg_rate          * phy%P    & 
              - graz_rate            * phy%Q%P    
  !  --- DETRITUS P 
   rhsv%detP = floppZ%P              * zoo%C    &
              + aggreg_rate          * phy%P    &
              - self%dil             * det%P    &         
              + zoo_mort             * zoo%P    & 
              - degradT              * det%P
  !  --- DOP
   rhsv%domP = exud%P                * phy%C    &
              + degradT              * det%P    &
              - self%dil             * dom%P    &              
              - reminT               * dom%P
  !  --- DIP
   rhsv%nutP = - uptake%P            * phy%C    & 
              + reminT               * dom%P    & 
              + lossZ%P              * zoo%C    &
              + self%dil * (self%nutP_initial - nut%P)
end if 
!________________________________________________________________________________
!
if (self%SiliconOn) then 
  ! ---  PHYTOPLANKTON Si
   rhsv%phyS = uptake%Si              * phy%C    & 
              - self%dil             * phy%Si    &            
              - aggreg_rate          * phy%Si    & 
              - graz_rate            * phy%Q%Si    
  !  --- DETRITUS Si
   rhsv%detS = floppZ%Si              * zoo%C    &
              + aggreg_rate          * phy%Si    &
              - self%dil             * det%Si    &         
              - degradT              * det%Si
  !  --- Dissolved Si
   rhsv%nutS = - uptake%Si            * phy%C    & 
              + degradT              * det%Si    & 
              + lossZ%Si              * zoo%C    &
              + self%dil * (self%nutS_initial - nut%Si)

end if 

!#S_ODE
!---------- ODE for each state variable ----------
  _SET_ODE_(self%id_nutN, rhsv%nutN UNIT)
  _SET_ODE_(self%id_phyC, rhsv%phyC UNIT)
  _SET_ODE_(self%id_phyN, rhsv%phyN UNIT)
  _SET_ODE_(self%id_detC, rhsv%detC UNIT)
  _SET_ODE_(self%id_detN, rhsv%detN UNIT)
  _SET_ODE_(self%id_domC, rhsv%domC UNIT)
  _SET_ODE_(self%id_domN, rhsv%domN UNIT)
if (self%RubiscoOn) then
      _SET_ODE_(self%id_Rub, rhsv%Rub UNIT)
end if
if (self%PhotoacclimOn) then
      _SET_ODE_(self%id_chl, rhsv%chl UNIT)
end if
if (self%PhosphorusOn) then
      _SET_ODE_(self%id_nutP, rhsv%nutP UNIT)
      _SET_ODE_(self%id_phyP, rhsv%phyP UNIT)
      _SET_ODE_(self%id_detP, rhsv%detP UNIT)
      _SET_ODE_(self%id_domP, rhsv%domP UNIT)
end if
if (self%SiliconOn) then
      _SET_ODE_(self%id_nutS, rhsv%nutS UNIT)
      _SET_ODE_(self%id_phyS, rhsv%phyS UNIT)
      _SET_ODE_(self%id_detS, rhsv%detS UNIT)
end if
if (self%GrazingOn) then
      _SET_ODE_(self%id_zooC, rhsv%zooC UNIT)
end if
!#E_ODE


!________________________________________________________________________________
! set diag variables, mostly from PrimProd module ______________
!#S_DIA
  _SET_DIAGNOSTIC_(self%id_chl2, phy%theta*phy%rel_chloropl) !last bulk chlorophyll concentration
  _SET_DIAGNOSTIC_(self%id_fracR, phy%frac%Rub)             !last 
  _SET_DIAGNOSTIC_(self%id_QN, phy%Q%N)                     !last 
  _SET_DIAGNOSTIC_(self%id_QP, phy%Q%P)                     !last 
  _SET_DIAGNOSTIC_(self%id_aVN, acclim%aV%N)                !last 
  _SET_DIAGNOSTIC_(self%id_aVP, acclim%aV%P)                !last 
  _SET_DIAGNOSTIC_(self%id_aVSi, acclim%aV%Si)              !last 
  _SET_DIAGNOSTIC_(self%id_rQSi, phy%relQ%Si)               !last 
  _SET_DIAGNOSTIC_(self%id_tmp, acclim%tmp)                 !last 
  _SET_DIAGNOSTIC_(self%id_fac1, acclim%fac1)               !last 
  _SET_DIAGNOSTIC_(self%id_fac2, acclim%fac2)               !last 
!#E_DIA

if (self%DebugDiagOn) then
!   _SET_DIAG_(self%id_chl_diag, phy%theta * phy%rel_chloropl ) !* phy%C
!   _SET_DIAG_(self%id_fracR, phy%frac%Rub) 
!   _SET_DIAG_(self%id_fracTheta, phy%frac%theta)
!   _SET_DIAG_(self%id_fracQN, phy%Q%N)
!   _SET_DIAG_(self%id_fracQP, phy%Q%P*1.0d3)
!   _SET_DIAG_(self%id_reg_VNC, uptake%N)
!   _SET_DIAG_(self%id_fac1,uptake%P  ) acclim%fac1
!   _SET_DIAG_(self%id_fac2,exud%P )
!   _SET_DIAG_(self%id_tmp,lossZ%P ) 

!   _SET_DIAG_(self%id_flxC, uptake%C*phy%C - reminT*dom%C- self%dil*totalC - lossZ%C*zoo%C) 
!   _SET_DIAG_(self%id_flxC,rhsv%phyP ) 

!   _SET_DIAG_(self%id_totC,totalC) 
!   _SET_DIAG_(self%id_totN, totalN) 
!   _SET_DIAG_(self%id_fac1, totalN - self%budget_N) 
!   if (self%PhosphorusOn) then
!     _SET_CONSERVED_QUANTITY_(self%id_totP,dip+phy%P+det%P+dom%P+zoo%P) 
!     _SET_DIAG_(self%id_totP,totalP) 
!     _SET_DIAG_(self%id_fac2, totalP - self%budget_P) 
!     _SET_DIAG_(self%id_fac1,floppZ%P*zoo%C) 
!     _SET_DIAG_(self%id_fac2,lossZ%P*zoo%C)   end if

end if

  _LOOP_END_

end subroutine maecs_do

!> @brief handles vertical movement for depth-varying movement rates
!> @details phyto sinking rate depends on the nutritional state, so for each node:
!! \n \f$ phy\%relQ \f$ obtained by calling calc_internal_states(self,phy,det,dom,zoo) 
!! \n then \f$ phyQstat=phy\%relQ\%N * phy\%relQ\%P \f$
!! \n finally, vsink = maecs_functions::sinking(self\%vS_phy, phyQstat, vsink)
subroutine maecs_get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)

use maecs_functions
use maecs_types

implicit none
!
! !INPUT PARAMETERS:
 class (type_hzg_maecs),intent(in) :: self
_DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_ 
 !   REALTYPE, intent(in)              ::vstokes 
type (type_maecs_phy):: phy !< maecs phytoplankton type
type (type_maecs_zoo) :: zoo
type (type_maecs_om):: det
type (type_maecs_om):: dom

!
! !LOCAL VARIABLES: 
REALTYPE    :: phyQstat,vsink
REALTYPE, parameter :: secs_pr_day = 86400.d0 
!EOP
!-----------------------------------------------------------------------
!BOC

_FABM_LOOP_BEGIN_
   
   ! Retrieve phtoplankton state
   
   !Retrieve the 'phyQstat' directly as a diagnostic variable: does not work yet.
   !fabm_get_bulk_diagnostic_data(self%id_phyqstat,phyQstatD) !where, phyQstat=relQ%N*relQ%P
   !_GET_(self%id_phyqstat,phyQstatD)
   
   !Calculate manually
   _GET_(self%id_phyC, phy%C)  ! Phytplankton Carbon in mmol-C/m**3
   _GET_(self%id_phyN, phy%N)  ! Phytplankton Nitrogen in mmol-N/m**3
   if (self%GrazingOn) then
     _GET_(self%id_zooC, zoo%C)  ! Zooplankton Carbon in mmol-C/m**3
   end if
   if (self%PhosphorusOn) then
     _GET_(self%id_phyP, phy%P)  ! Phytplankton Phosphorus in mmol-P/m**3
   end if

   !write (*,'(A,2(F10.3))') 'Before: phy%C, phy%N=', phy%C, phy%N
   call min_mass(self,phy,method=2) 
   !write (*,'(A,2(F10.3))') 'After: phy%C, phy%N=', phy%C, phy%N
   call calc_internal_states(self,phy,det,dom,zoo) 
   !write (*,'(A,2(F10.3))') 'phy%relQ%N, phy%relQ%P=', phy%relQ%N, phy%relQ%P
   
   !calculate Q state
   phyQstat = phy%relQ%N * phy%relQ%P 

  
   ! Calculate sinking

   call sinking(self%vS_phy, phyQstat, vsink)
   vsink = vsink / secs_pr_day
   !write (*,'(A,2(F10.3))') 'phyQstat, vsink=', phyQstat, vsink
   
   !set the rates
   _SET_VERTICAL_MOVEMENT_(self%id_detC,-1.0_rk*self%vS_det/secs_pr_day)
   _SET_VERTICAL_MOVEMENT_(self%id_detN,-1.0_rk*self%vs_det/secs_pr_day)

   _SET_VERTICAL_MOVEMENT_(self%id_phyN,vsink)
   _SET_VERTICAL_MOVEMENT_(self%id_phyC,vsink)
   if (self%PhosphorusOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_phyP,vsink)
      _SET_VERTICAL_MOVEMENT_(self%id_detP,-1.0_rk*self%vs_det/secs_pr_day)
   end if
   if (self%SiliconOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_phyS,vsink)
      _SET_VERTICAL_MOVEMENT_(self%id_detS,-1.0_rk*self%vs_det/secs_pr_day)
   end if
   if (self%PhotoacclimOn) then 
      _SET_VERTICAL_MOVEMENT_(self%id_chl,vsink)
      _SET_VERTICAL_MOVEMENT_(self%id_Rub,vsink)
   end if
  
_FABM_LOOP_END_
  
end subroutine maecs_get_vertical_movement