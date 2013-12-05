#include "fabm_driver.h"
!-----------------------------------------------------------------------------------
!          Model for Adaptive Ecosystems in Coastal Seas 
<<<<<<< HEAD
!  Original author(s): Richard Hofmeister, Markus Schartau & Kai Wirtz
!  HZG 2011-2013

=======
>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce
subroutine maecs_do(self,_ARGUMENTS_DO_)

use fabm_types
use fabm_driver
use maecs_types
use maecs_functions
use maecs_primprod 
use maecs_grazing

! !INPUT PARAMETERS:
 class (type_maecs_base_model),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
type (type_maecs_rhs) :: rhsv
type (type_maecs_phy) :: phy   ! phytoplankton type containing state and trait information
type (type_maecs_zoo) :: zoo   ! zooplankton type 
type (type_maecs_om)  :: dom, det, nut, uptake, exud, lossZ, floppZ 
type (type_maecs_env) :: env
type (type_maecs_traitdyn) :: acclim
type (type_maecs_sensitivities) :: sens

! --- LOCAL MODEL VARIABLES:
integer  :: i, j, iz, ihour, iloop
real(rk) :: reminT, degradT       ! Temp dependent remineralisation and hydrolysis rates
! --- QUOTA and FRACTIONS
real(rk) :: phys_status, dQN_dt, dRchl_phyC_dt=0.0_rk ! []

! --- ZOOPLANKTON GRAZING, INGESTION, MORTALITY, RESPIRATION... 
real(rk) :: graz_rate   ! carbon-specific grazing rate                          [d^{-1}]
<<<<<<< HEAD
real(rk) :: zoo_respC   ! temperature dependent carbon respiration rate  [mmolC m^{-3} d^{-1}]
real(rk) :: zoo_mort
real(rk) :: secs_pr_day = 86400.0_rk
=======
real(rk) :: yield_zoo   ! grazed fraction that is actually ingested            [dimensionless]
real(rk) :: zoo_respC   ! temperature dependent carbon respiration rate  [mmolC m^{-3} d^{-1}]
real(rk) :: zoo_mort,rub3

>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce
! --- AGGREGATION 
real(rk) :: aggreg_rate ! aggregation among phytoplankton and between phytoplankton & detritus [d^{-1}]    
logical  :: out = .true.
!   if(36000.eq.secondsofday .and. mod(julianday,1).eq.0 .and. outn) out=.true.    
#define _KAI_ 0
#define _MARKUS_ 1
<<<<<<< HEAD
! #define UNIT / 86400
#define UNIT *1.1574074074E-5_rk
=======
#define UNIT / 86400
!#define UNIT *1.1574074074E-5_rk
>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce

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
<<<<<<< HEAD
if (self%RubiscoOn) then
      _GET_(self%id_Rub, phy%Rub)  ! fraction of Rubisco in -
=======
if (self%GrazingOn) then
  _GET_(self%id_zooC, zoo%C)  ! Zooplankton Carbon in mmol-C/m**3
end if


if (self%RubiscoOn) then
      _GET_(self%id_Rub2, phy%Rub)  ! fraction of Rubisco in -

>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce
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
<<<<<<< HEAD
if (self%GrazingOn) then
      _GET_(self%id_zooC, zoo%C)  ! Zooplankton Carbon in mmol-C/m**3
end if
=======
>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce
!#E_GET
!phy%Rub = Rub
!phy%chl = chl
! Retrieve current environmental conditions.
!#S_GED
  _GET_(self%id_temp, env%temp)  ! water temperature
  _GET_(self%id_par, env%par)  ! light photosynthetically active radiation
!#E_GED

! *** PHYTOPLANKTON EQUATIONS 
! We distinguish between mass state variables (in units of carbon, nitrogen, & phosphorus) 
!   and property state variables. 
! Current 'traits' are: 1) nitrogen allocated to rubisco (frac_Rub) 
!                 2) Chla content of chloroplasts  (theta)
! General code structure:
!  A) Calculation of quotas, and mass exchange rates
!  B) Specify rates of change of property state variables ('traits') 
!  C) Assign mass exchange rates ('rhs(j,i)')
!  D) Assign rates of change of 'traits' property variables  

! *** BEGIN A) 

! --- checking and correcting extremely low state values  ------------  
<<<<<<< HEAD
call min_mass(self,phy,method=2) !_KAI_ minimal reasonable Phy-C and -Nitrogen
=======
call min_mass(self,phy,method=_KAI_) ! minimal reasonable Phy-C and -Nitrogen
>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce
!call min_mass(self,phy,method=3)
!write (*,'(A,3(F10.3))') '3 rub=',phy%Rub

! --- stoichiometry of autotrophs (calculating QN_phy, frac_R, theta, and QP_phy)
call calc_internal_states(self,phy,det,dom,zoo)
!write (*,'(A,2(F10.3))') '2 P,chl=',phy%C,phy%chl
!write (*,'(A,3(F10.3))') '2 Relchl=',phy%rel_chloropl,phy%C_reg,phy%QN
!write (*,'(A,4(F10.3))') '2 N,th=',phy%N,phy%theta,phy%frac%Rub,phy%rel_QN

if (.not. self%PhotoacclimOn) then  
   phy%chl = phy%C * self%frac_chl_ini   ! total Chl mg-CHL/m3
   phy%frac%theta  = self%frac_chl_ini * phy%QN  * self%itheta_max
   phy%theta       = self%frac_chl_ini * phy%QN / self%frac_Rub_ini! g-CHL/mol-N*m3
!   phy%frac%Rub   = phy%Rub / phy%C_reg
!   phy%frac%theta =phy%N*self%Chl_N_ini/phy%C_reg*self%itheta_max !
end if 
!write (*,'(A,2(F10.3))') 'chl theta:',phy%chl,phy%theta
!write (*,'(A,2(F10.3))') 'C rub=',phy%C,phy%frac%Rub

call calc_sensitivities(self,sens,phy,env,nut)

! --- ALGAL GROWTH and EXUDATION RATES, physiological trait dynamics ----------------
call photosynthesis(self,sens,phy,uptake,exud,acclim)
! write (*,'(A,4(F10.4))') 'upt C P=',uptake%C,uptake%P*1E3,phy%P*1E6,phy%QP*1E3
!write (*,'(A,4(F10.3))') '3 N,th=',phy%N,phy%theta,phy%frac%Rub,phy%rel_QN

! ----------------       grazing        -------------------------------------------
if (self%GrazingOn) then
  call grazing(self%g_max * sens%func_T,self%k_grazC,phy%C,graz_rate)
  zoo%feeding = graz_rate
  zoo_respC   = self%basal_resp_zoo * sens%func_T  !  basal respiration of grazers

! --- calculates zooplankton loss rates (excretion->Nut, floppy+egestion->Det), specific to C
  call grazing_losses(zoo,zoo_respC,phy%QN,phy%QP,lossZ,floppZ,self%PhosphorusOn, .true.) 
!  --- transform from specific to bulk grazing rate
  graz_rate   = graz_rate * zoo%C 
!  --- quadratic closure term
  zoo_mort    = self%mort_zoo * sens%func_T  * zoo%C

else
  graz_rate   = 0.0_rk
  lossZ       = type_maecs_om(0.0_rk, 0.0_rk, 0.0_rk)
  floppZ      = type_maecs_om(0.0_rk, 0.0_rk, 0.0_rk)
end if

! --- phytoplankton aggregation -------------------------------------------------
! If biovolume is primarily determined by the nitrogen content, also for detritus
!aggreg_rate = self%phi_agg * dom%C * (phy%N + det%N)                    ! [d^{-1}] 
<<<<<<< HEAD

!_GET_(self%id_fracR,phys_status )  
!write (*,'(A,1(F10.3))') 'phys=',phys_status
=======
! _SET_DIAGNOSTIC_(self%id_fracR,phy%frac%rel_phys)  

!_GET_(self%id_fracR,phys_status )  
!write (*,'(A,1(F10.3))') 'phys=',phys_status
! _SET_DIAGNOSTIC_(self%id_tmp,phy%frac%rel_phys ) !step_integrated bulk chlorophyll concentration
>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce

aggreg_rate = self%phi_agg * (1.0_rk - exp(-0.02*dom%C)) * (phy%N + det%N) 
!         vS * exp(-4*phys_status )                ! [d^{-1}] 
!aggreg_rate = aggreg_rate * exp(-4*phy%rel_phys ) 

! ------------------------------------------------------------------
!  ---  hydrolysis & remineralisation rate (temp dependent)
degradT     = self%hydrol * sens%func_T
reminT      = self%remin  * sens%func_T

! right hand side of ODE (rhs)    
<<<<<<< HEAD
!  #define UNIT /self%secs_pr_day
=======
!#define UNIT /self%secs_pr_day
>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce
!__________________________________________________________________________
!
! PHYTOPLANKTON C
rhsv%phyC = uptake%C              * phy%C &
           - self%dil             * phy%C &
           - exud%C               * phy%C &  !TODO: move loss rates to mu, also checking for 
           - aggreg_rate          * phy%C &  !      trait dependencies
           - graz_rate                    

!!! write (*,'(A,2(F9.4))') 'flxc=',uptake%C, phy%C

!_____________________________________________________________________________
!
! PHYTOPLANKTON N
rhsv%phyN =  uptake%N             * phy%C &
           - exud%N               * phy%C & 
           - aggreg_rate          * phy%N &
           - self%dil             * phy%N &          
           - graz_rate * phy%QN          

!_____________________________________________________________________________

if (self%PhotoacclimOn) then
! PHYTOPLANKTON CHLa
     ! note that theta*rel_chloropl in units [mg Chla (mmol C)^{-1}] 
     !             = phy%theta*(self%rel_chloropl_min+phy%frac%Rub*phy%rel_QN**self%sigma) 
   dQN_dt        = (rhsv%phyN * phy%C_reg - rhsv%phyC * phy%N_reg) / (phy%C_reg*phy%C_reg)
! TODO: dangerous to work with RHS instead of net uptake rates (mortality has no physiological effect)

   dRchl_phyC_dt =  acclim%dRchl_dtheta * acclim%dtheta_dt   & 
                  + acclim%dRchl_dfracR * acclim%dfracR_dt   & 
                  + acclim%dRchl_dQN    * dQN_dt 

<<<<<<< HEAD
   rhsv%chl = phy%theta * phy%frac%Rub * phy%rel_QN**self%sigma * rhsv%phyC + dRchl_phyC_dt * phy%C
=======
   rhsv%chl = phy%theta * phy%frac%Rub * phy%rel_QN**self%sigma * rhsv%phyC + dRchl_phyC_dt * phy%C_reg
>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce

!write (*,'(A,4(F10.3))') 'rhs chl=', phy%theta * phy%frac%Rub * phy%rel_QN**self%sigma * rhsv%phyC,dRchl_phyC_dt * phy%C_reg*1E1,phy%rel_QN**self%sigma,phy%theta

!_____________________________________________ _________________________________

if (self%RubiscoOn) then 
!        rhsv%rub  = acclim%dfracR_dt * phy%N_reg + phy%Rub / phy%N_reg * rhsv%phyN  !phy%Rub / phy%C_reg
<<<<<<< HEAD
    rhsv%Rub  = acclim%dfracR_dt * phy%C + phy%Rub/phy%C_reg * rhsv%phyC 
! add relaxation term to destroy unphysical Rub accumulation
!     if( phy%Rub .gt. 0.9d0*phy%C) rhsv%rub = rhsv%rub - phy%Rub/(1.d0+exp(-67.d0*(phy%Rub/ phy%C_reg  - 1.d0)))
!write (*,'(A,2(F10.3))') '2:',acc%dfracR_dt*1E3, grad_fracR *1E3
=======
    rhsv%Rub2  = acclim%dfracR_dt * phy%C_reg + phy%Rub/phy%C_reg * rhsv%phyC 
! add relaxation term to destroy unphysical Rub accumulation
!     if( phy%Rub .gt. 0.9d0*phy%C) rhsv%rub = rhsv%rub - phy%Rub/(1.d0+exp(-67.d0*(phy%Rub/ phy%C_reg  - 1.d0)))
!write (*,'(A,2(F10.3))') '2:',acc%dfracR_dt*1E3, grad_fracR *1E3

!write (*,'(A,5(F10.3))') 'rub=',phy%Rub,phy%C,phy%Rub/phy%C_reg,phy%frac%Rub,rhsv%Rub2-phy%Rub / phy%C_reg * rhsv%phyC
 _SET_DIAGNOSTIC_(self%id_tmp,phy%Rub ) !step_integrated bulk chlorophyll concentration
!  _SET_DIAGNOSTIC_(self%id_tmp,acclim%dfracR_dt ) !step_integrated bulk chlorophyll concentration

>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce
   end if 
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
!write (*,'(A,2(F10.3))') 'zoo=',zoo%C,rhsv%zooC
!if (zoo%C .gt. 10.d0) write (*,'(A,3(F10.3))') 'RHS zoo=',zoo%C,rhsv%zooC,zoo%yield * graz_rate
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
!zoo%flopp * phy%QN * graz_rate   &
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
              - graz_rate            * phy%QP    
!write (*,'(A,3(F10.3))') 'P=',phy%P,1E3*phy%P/phy%C_reg,1E3*phy%QP
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

!_SET_ODE_(self%id_phyC,rhsv%phyC UNIT)

!#S_ODE
!---------- ODE for each state variable ----------
  _SET_ODE_(self%id_nutN, rhsv%nutN UNIT)
  _SET_ODE_(self%id_phyC, rhsv%phyC UNIT)
  _SET_ODE_(self%id_phyN, rhsv%phyN UNIT)
  _SET_ODE_(self%id_detC, rhsv%detC UNIT)
  _SET_ODE_(self%id_detN, rhsv%detN UNIT)
  _SET_ODE_(self%id_domC, rhsv%domC UNIT)
  _SET_ODE_(self%id_domN, rhsv%domN UNIT)
<<<<<<< HEAD
if (self%RubiscoOn) then
      _SET_ODE_(self%id_Rub, rhsv%Rub UNIT)
=======

if (self%GrazingOn) then
  _SET_ODE_(self%id_zooC, rhsv%zooC UNIT)
end if 
if (self%RubiscoOn) then
      _SET_ODE_(self%id_Rub2, rhsv%Rub2 UNIT)
!      _SET_ODE_(self%id_Rub2, 0.0_rk  UNIT)
>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce
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
<<<<<<< HEAD
if (self%GrazingOn) then
      _SET_ODE_(self%id_zooC, rhsv%zooC UNIT)
end if
=======
>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce
!#E_ODE


!________________________________________________________________________________
! set diag variables, mostly from PrimProd module ______________
!#S_DIA
  _SET_DIAGNOSTIC_(self%id_chl2, phy%theta*phy%rel_chloropl) !step_integrated bulk chlorophyll concentration
  _SET_DIAGNOSTIC_(self%id_fracR, phy%frac%Rub)             !step_integrated 
<<<<<<< HEAD
  _SET_DIAGNOSTIC_(self%id_rhs_nutN, rhsv%nutN) !step_integrated RHS of Dissolved Inorganic Nitrogen DIN
  _SET_DIAGNOSTIC_(self%id_rhs_nutP, rhsv%nutP) !step_integrated RHS of Dissolved Inorganic Phosphorus DIP
  _SET_DIAGNOSTIC_(self%id_rhs_phyC, rhsv%phyC) !step_integrated RHS of Phytplankton Carbon
  _SET_DIAGNOSTIC_(self%id_rhs_phyN, rhsv%phyN) !step_integrated RHS of Phytplankton Nitrogen
  _SET_DIAGNOSTIC_(self%id_rhs_phyP, rhsv%phyP) !step_integrated RHS of Phytplankton Phosphorus
  _SET_DIAGNOSTIC_(self%id_rhs_zooC, rhsv%zooC) !step_integrated RHS of Zooplankton Carbon
  _SET_DIAGNOSTIC_(self%id_rhs_detC, rhsv%detC) !step_integrated RHS of Detritus Carbon
  _SET_DIAGNOSTIC_(self%id_rhs_detN, rhsv%detN) !step_integrated RHS of Detritus Nitrogen
  _SET_DIAGNOSTIC_(self%id_rhs_detP, rhsv%detP) !step_integrated RHS of Detritus Phosphorus
  _SET_DIAGNOSTIC_(self%id_rhs_domC, rhsv%domC) !step_integrated RHS of Dissolved Organic Carbon
  _SET_DIAGNOSTIC_(self%id_rhs_domN, rhsv%domN) !step_integrated RHS of Dissolved Organic Nitrogen
  _SET_DIAGNOSTIC_(self%id_rhs_domP, rhsv%domP) !step_integrated RHS of Dissolved Organic Phosphorus
  _SET_DIAGNOSTIC_(self%id_rhs_Rub, rhsv%Rub)   !step_integrated RHS of fraction of Rubisco
  _SET_DIAGNOSTIC_(self%id_rhs_chl, rhsv%chl)   !step_integrated RHS of Chl:C ratio
=======
>>>>>>> 560e3dbbe7b9557d60044e0e8f65e141e09252ce
!#E_DIA

if (self%DebugDiagOn) then
!   _SET_DIAG_(self%id_chl_diag, phy%theta * phy%rel_chloropl ) !* phy%C
!   _SET_DIAG_(self%id_fracR, phy%frac%Rub) 
!   _SET_DIAG_(self%id_fracTheta, phy%frac%theta)
!   _SET_DIAG_(self%id_fracQN, phy%QN)
!   _SET_DIAG_(self%id_fracQP, phy%QP*1.0d3)
!   _SET_DIAG_(self%id_reg_VNC, uptake%N)
!   _SET_DIAG_(self%id_fac1,uptake%P  ) 
!   _SET_DIAG_(self%id_fac2,exud%P )
!   _SET_DIAG_(self%id_tmp,lossZ%P ) 

   !_SET_DIAG_(self%id_flxC, uptake%C*phy%C - reminT*dom%C- self%dil*totalC - lossZ%C*zoo%C) 
!   _SET_DIAG_(self%id_flxC,rhsv%phyP ) 
!!!   write (*,'(A,4(F9.4))') 'flxc=',uptake%C, phy%C,uptake%C * phy%C,uptake%C*phy%C - reminT*dom%C- self%dil*totalC - lossZ%C*zoo%C

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
