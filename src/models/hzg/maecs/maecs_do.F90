! #include "fabm_driver.h"

!> @brief This is the main routine where right-hand-sides are calculated
!> @details
!> NOTE: Although this subroutine looks as if it's not a part of any module,
!! it is temporarily included in the fabm_hzg_maecs module (inside maecs.F90)
!! when compiling the documentation, such that the subroutine is documented 
!! under the 'Data Type Documentation' chapter, where the in-body docs are also listed
!>
!> **Phytoplankton Equations** 
!> \n We distinguish between mass state variables 
!! (in units of carbon, nitrogen, & phosphorus) and property state variables. 
!! \lref{For a textual narration and equations, see sec.,sec:ModStr,.}\n
!>  Current 'traits' are: 
!> - nitrogen allocated to rubisco [-] (frac_Rub) 
!> - Chla content of chloroplasts [chl-a/chl-C] (theta) 
!>
!> **General code structure:**
!> 1. Calculation of quotas, internal states, potential rates
!> 2. Calculation of fluxes, mass exchange rates &  rates of change of traits variables 
!> 3. Assign mass exchange rates ('rhs(j,i)')
!> 4. Assign rates of change of 'traits' property variables
!>
!> \n **Detailed Descriptions:**
! @todo: althougth HIDE_IN_BODY_DOCS=NO, the body-documentation is not included! HAS TO BE FIXED
! @todo: add equations
! @todo: why UNIT instead of secs_pr_day? 
! @todo: 'sensitivities' does not seem to be a proper name choice. Something more intuitive?
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
real(rk) :: denitrate   ! pelagic N-loss by denitrification, emulating benthic pool and suboxic micro-environments
real(rk) :: deporate    ! pelagic, "volumetric" deposition, slowly refueling N-losses
real(rk) :: qualPOM, ddegN, ddegP     !  POM quality -> degradation
real(rk) :: qualDOM, dremN, dremP     !  DOM quality -> degradation
real(rk) :: secs_pr_day = 86400.0_rk
! --- AGGREGATION 
real(rk) :: aggreg_rate ! aggregation among phytoplankton and between phytoplankton & detritus [d^{-1}]    
logical  :: out = .true.
!   if(36000.eq.secondsofday .and. mod(julianday,1).eq.0 .and. outn) out=.true.    
real(rk) :: pdet, no3
real(rk) :: det_prod, nh3f
real(rk) :: radsP,Oxicminlim,Denitrilim,Anoxiclim,Rescale,rP
real(rk),parameter :: relaxO2=0.04_rk
real(rk),parameter :: T0 = 288.15_rk ! reference Temperature fixed to 15 degC
real(rk),parameter :: Q10b = 1.5_rk
real(rk) :: Cprod, Nprod, Pprod
real(rk) :: AnoxicMin,Denitrific,OxicMin,Nitri,OduDepo,OduOx,pDepo, Anammox
real(rk) :: prodO2, rhochl, uptNH4, uptNO3, uptchl, uptN, respphyto,faeces

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
      _GET_(self%id_chl, phy%chl)  ! Chl in mg-Chla/mmol-C
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
if (self%BioOxyOn) then
      _GET_(self%id_nh3, env%nh3)  ! dissolved ammonium in mmolN/m**3
      _GET_(self%id_oxy, env%oxy)  ! dissolved oxygen in mmolO2/m**3
      _GET_(self%id_odu, env%odu)  ! dissolved reduced substances in mmolO2/m**3
end if
if (self%NResOn) then
      _GET_(self%id_RNit, env%RNit)  ! N-reservoir in mmol-N/m**3
end if
!#E_GET
if (.not. self%GrazingOn) then
      zoo%C = self%zooC_initial
end if

!phy%Rub = Rub
!phy%chl = chl
! Retrieve current environmental conditions.


!S_GED
  _GET_(self%id_temp, env%temp)  ! water temperature
  _GET_(self%id_par, env%par)  ! light photosynthetically active radiation
!E_GED  ! list outcommented due to different usage of zmax and doy (see light extinction)

! write (*,'(A,2(F10.3))') 'par/T:',env%par,env%temp


! @ingroup main
!> @fn fabm_hzg_maecs::maecs_do () 
!> 1. Calculation of quotas, internal states, potential rates
!>   - call min_mass with method=2, store phy\%C and \%N in phy\%reg 
!>   - call calc_internal_states: retrieve phy\%Q\%X, phy\%theta, phy\%frac\%X 
!>   - if PhotoacclimOn=.false., calculate: 
!>     - phy\%chl=phy\%C * self\%frac_chl_ini 
!>     - phy\%frac\%theta = self\%frac_chl_ini * self\%itheta_max
!>     - phy\%theta= self\%frac_chl_ini / (self\%frac_Rub_ini * phy\%relQ\%N**self\%sigma)
!>   - call calc_sensitivities: retrieve potential rates: @f$f_T@f$, sens\%upt\_pot\%C (=LH), sens\%upt\_pot\%X (= @f$ V_X @f$), sens\%P\_max
!> @todo: min_mass correction of phy%\C and phy\%N at this stage requires specification of threshold values. What about back-calculating phy\%reg\%N from the smooth_small corrected phy\%Q\%N?

! --- checking and correcting extremely low state values  ------------  
call min_mass(self,phy,method=2) !_KAI_ minimal reasonable Phy-C and -Nitrogen

! --- stoichiometry of autotrophs (calculating QN_phy, frac_R, theta, and QP_phy)
call calc_internal_states(self,phy,det,dom,zoo)

!write (*,'(A,2(F10.3))') 'PAR, chl=',env%par, phy%chl

if (.not. self%PhotoacclimOn) then  
   phy%chl         = phy%C * self%frac_chl_ini   ! total Chl mg-CHL/m3
   phy%frac%theta  = self%frac_chl_ini * self%itheta_max
   phy%theta       = self%frac_chl_ini /(self%frac_Rub_ini*phy%relQ%N**self%sigma)
! g-CHL/mol-C*m3
! write (*,'(A,2(F10.3))') 'theta:',phy%relQ%N**self%sigma,phy%theta   
end if 

call calc_sensitivities(self,sens,phy,env,nut,acclim)

!write (*,'(A,4(F10.3))') 'f=',phy%reg%C,phy%frac%Rub,self%small_finite + self%rel_chloropl_min,phy%frac%NutUpt
!write (*,'(A,4(F10.3))') 'PAR, T, th, P =',env%par,env%temp,phy%theta, sens%upt_pot%C 

!> @fn fabm_hzg_maecs::maecs_do ()
!> 2. Calculation of fluxes, mass exchange rates &  rates of change of traits variables 
!! & Specify rates of change of traits variables
!>   - call maecs_primprod::photosynthesis(): this is where everything happens!
!>   - if GrazingOn: 
!>     - graz\_rate=rate retrieved from call maecs_grazing::grazing()
!>     - lossZ\%X=lossZNut\%X, floppZ\%X=lossZDet\%X retrieved from call maecs_grazing::grazing_losses()
!>     - calculate graz_rate retr= graz_rate * zoo\%C and zoo_mort
!>   - calc. aggreg_rate @f$ = \mathrm{phi\_agg}*(1-e^{-0.02*\mathrm{dom\%C}}) * phy\%N *det\%N @f$
!>     - future work: aggreg_rate=f(size), \latexonly (see section \ref{sec:partagg}) \endlatexonly \n
!>   - calc. degradT=self\%hydrol * @f$ f_T @f$ and reminT=self\%remin * @f$ f_T @f$
!> @todo: aggregation equation: where does 0.02 come from? are the results sensitive to this par? 
!> @todo: specific graz_rate becomes pop. grazing rate. Do this at the rhs calculations
!> @todo: graz_rate: no temperature modification: forgotten?

! --- ALGAL GROWTH and EXUDATION RATES, physiological trait dynamics ----------------
call photosynthesis(self,sens,phy,uptake,exud,acclim)

! ----------------       grazing        -------------------------------------------
if (self%GrazingOn) then
  call grazing(self%g_max * sens%f_T,self%k_grazC,phy%C,graz_rate)
  zoo%feeding = graz_rate
  zoo_respC   = self%basal_resp_zoo * sens%f_T  !  basal respiration of grazers
  nquot       = type_maecs_om(1.0_rk, phy%Q%N, phy%Q%P, phy%Q%Si )
  mswitch     = type_maecs_switch(self%PhosphorusOn,self%SiliconOn,.true. )
                                  !isP, isSi, isTotIng
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

aggreg_rate = self%phi_agg * (1.0_rk - exp(-self%agg_doc*dom%C)) * (phy%N + det%N) 
!         vS * exp(-4*phys_status )                ! [d^{-1}] 
!aggreg_rate = aggreg_rate * exp(-4*phy%rel_phys ) TODO: DOM quality as proxy for TEP

!_____________________________________________________________________________
!
!      turnover of long-term N-reservoir (denitrification + wet N-deposition)
if (self%NResOn) then
!pelagic N-loss by denitrification, emulating benthic pool and suboxic micro-environments
 denitrate = self%denit * 4 * sens%f_T * (1.0d0 - exp(-det%N/self%PON_denit)) * det%N

!pelagic, "volumetric" deposition, slowly refueling N-losses 
 deporate  = self%denit * exp(-4*sens%f_T) * env%RNit

 rhsv%RNit = denitrate - deporate 
else
  deporate  = 0.0d0
  denitrate = 0.0d0
endif    


!> @fn fabm_hzg_maecs::maecs_do ()
!> 3. Assign mass exchange rates ('rhs(j,i)')
!>   - phyC= uptake - dil - exud - aggreg\_rate - graz\_rate
!> @todo: add the rhs equations


! right hand side of ODE (rhs)    
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
!>    - if PhotoacclimOn: rhsv%chl=A + B
!>      + A = rhsv\%phyC*phy\%theta*rel_chlorpl [= gC/m3/d * gchl/gchlorpl-C * gchloropl-C/gC]
!>      + B = dRchl/dtheta* dtheta/dt + dRchl/dfracR*dfracR/dt+dRchl/dQN * dQN/dt
!>      + all terms in B except dQN/dt are calculated in maecs_primprod::photosynthesis()
!>      + dQN/dt = (rhsv\%phyN* phyC - rhsv\*phyC*phyN )/(phyC^2)
!>    - if RubiscoOn: rhsv%Rub=A+B
!>      + A = rhsv\%phyC * phyRub/phyC
!>      + B = dfracR_dt is calculated in maecs_primprod::photosynthesis()

if (abs(phy%reg%C) .gt. 1d-4) then
 if (self%PhotoacclimOn ) then ! check for too small biomasses 

! PHYTOPLANKTON CHLa
     ! note that theta*rel_chloropl in units [mg Chla (mmol C)^{-1}] 
   dQN_dt        = (rhsv%phyN * phy%reg%C - rhsv%phyC * phy%reg%N) / (phy%reg%C*phy%reg%C)
! TODO: dangerous to work with RHS instead of net uptake rates (mortality has no physiological effect)

   dRchl_phyC_dt =  acclim%dRchl_dtheta * acclim%dtheta_dt   & 
                  + acclim%dRchl_dfracR * acclim%dfracR_dt   & 
                  + acclim%dRchl_dQN    * dQN_dt 

!   rhsv%chl = phy%theta * phy%frac%Rub * phy%relQ%N**self%sigma * rhsv%phyC + dRchl_phyC_dt * phy%C
   rhsv%chl = phy%theta * phy%frac%Rub * phy%relQ%N**self%sigma * rhsv%phyC + dRchl_phyC_dt * phy%reg%C

!write (*,'(A,4(F10.3))') 'rhs chl=', phy%theta * phy%frac%Rub * phy%relQ%N**self%sigma * rhsv%phyC,dRchl_phyC_dt * phy%reg%C*1E1,phy%relQ%N**self%sigma,phy%theta

!_____________________________________________ _________________________________
 end if 

 if (self%RubiscoOn) then 
   rhsv%Rub  = phy%Rub/phy%reg%C * rhsv%phyC + acclim%dfracR_dt * phy%C 
 end if 
else
  rhsv%Rub  = 0.0d0
  rhsv%chl  = 0.0d0
endif

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

! ------------------------------------------------------------------
!  ---  POM&DOM quality, relative to max N-quota of phytoplankton
 qualPOM     = (1.0d0-self%Nqual) + self%Nqual * det%N /(det%C + self%small_finite)  * self%iK_QN
 qualDOM     = (1.0d0-self%Nqual) + self%Nqual * dom%N /(dom%C + self%small_finite)  * self%iK_QN

!  ---  hydrolysis & remineralisation rate (temp dependent)
degradT     = self%hydrol * sens%f_T * qualPOM
reminT      = self%remin  * sens%f_T * qualDOM

!  ---  hydrolysis & remineralisation depend on quality, here propto N/C quota of OM
!  acceleration: rate difference for N-pool
ddegN       = self%hydrol * sens%f_T * smooth_small(1.0d0 - qualPOM, self%small_finite)
ddegP       = self%remNP * ddegN       
dremN       = self%remin * sens%f_T * smooth_small(1.0d0 - qualDOM, self%small_finite)
dremP       = self%remNP * dremN
!________________________________________________________________________________
!
!  --- DETRITUS C
det_prod    = floppZ%C              * zoo%C   &
             + aggreg_rate          * phy%C   &
             + zoo_mort             * zoo%C   

rhsv%detC   = det_prod                        &
             - self%dil             * det%C   &             
             - degradT              * det%C                

!________________________________________________________________________________
!
!  --- DETRITUS N
rhsv%detN   = floppZ%N              * zoo%C   &
             + aggreg_rate          * phy%N   &
             - self%dil             * det%N   &
             + zoo_mort             * zoo%N   & 
             - (degradT + ddegN)    * det%N   &
             - denitrate
  
!________________________________________________________________________________
!
!  --- DOC
 Cprod      = reminT                * dom%C
rhsv%domC   = exud%C                * phy%C   & 
             + degradT              * det%C   &
             - self%dil             * dom%C   &        
             - Cprod 
!________________________________________________________________________________
!
!  --- DON
Nprod       = (reminT + dremN)    * dom%N
rhsv%domN   = exud%N                * phy%C   &
             + (degradT + ddegN)    * det%N   &
             - self%dil             * dom%N   &         
             - Nprod
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
             + Nprod                          &
             + lossZ%N             * zoo%C    &
             + self%dil * (self%nutN_initial - nut%N) &
             + deporate
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
              - (degradT + ddegP)    * det%P    ! quality enhances P remin
  !  --- DOP
   Pprod     = (reminT + dremP)      * dom%P
   rhsv%domP = exud%P                * phy%C    &
              + (degradT + ddegP)    * det%P    &
              - self%dil             * dom%P    &              
              - Pprod
  !  --- DIP
   rhsv%nutP = - uptake%P            * phy%C    & 
              + Pprod                           & 
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

!---------- RHS for BGC/diagensis model part ----------
! Fortran 2003 version of OMEXDIA+P biogeochemical model
! The OMEXDIA+P+MPB model is based on the OMEXDIA model (see Soetard et al. 1996a)
! P-cycle is added by kai wirtz

if (self%BioOxyOn) then

! ---------- temperature    TODO: retrieve from existing temp variables 
   f_T    = sens%f_T

! ---------- manages overlapping state variables 
   no3    = smooth_small(nut%N - env%nh3,  self%small)
! ---------- remineralisation limitations 
   Oxicminlim = env%oxy/(env%oxy+self%ksO2oxic+relaxO2*(env%nh3+env%odu))              
   Denitrilim = (1.0_rk-env%oxy/(env%oxy+self%kinO2denit)) * no3/(no3+self%ksNO3denit)
   Anoxiclim  = (1.0_rk-env%oxy/(env%oxy+self%kinO2anox)) * (1.0_rk-no3/(no3+self%kinNO3anox))
   Rescale    = 1.0_rk/(Oxicminlim+Denitrilim+Anoxiclim)

! extra-omexdia P -dynamics  
  if (self%PhosphorusOn) then
! PO4-adsorption ceases when critical capacity is reached
! [FeS] approximated by ODU
!   po4    = nut%P
   radsP      = self%rPAds * degradT * nut%P * max(env%odu,self%PAdsODU)
   rhsv%nutP  = rhsv%nutP - radsP
   rhsv%detP  = rhsv%detP + radsP
!   rP     = self%rFast * (1.0_rk - Oxicminlim)
!   Pprod  = rP * pdet
endif

! Oxic mineralisation, denitrification, anoxic mineralisation
! then the mineralisation rates
   OxicMin    = Cprod*Oxicminlim*Rescale        ! oxic mineralisation
   Denitrific = Cprod*Denitrilim*Rescale        ! Denitrification
   AnoxicMin  = Cprod*Anoxiclim *Rescale        ! anoxic mineralisation

! reoxidation and ODU deposition
   Nitri      = f_T * self%rnit   * env%nh3 * env%oxy/(env%oxy + self%ksO2nitri &
                  + relaxO2*(dom%C + env%odu))
   OduOx      = f_T * self%rODUox * env%odu * env%oxy/(env%oxy + self%ksO2oduox &
                  + relaxO2*(env%nh3 + dom%C))

!  pDepo      = min(1.0_rk,0.233_rk*(wDepo)**0.336_rk )
   pDepo      = 0.0_rk
   OduDepo    = AnoxicMin * pDepo 

!  dynamics of env%oxy ~ dissolved oxygen
  rhsv%oxy    = - OxicMin - 2.0_rk* Nitri - OduOx &
                - lossZ%C * zoo%C + uptake%C * phy%C  &
                + self%dil * (self%oxy_initial - env%oxy)

!  dynamics of odu ~ dissolved reduced substances
  rhsv%odu    = (AnoxicMin - OduOx - OduDepo) 
!  dynamics of no3 ~ dissolved nitrate
!!  rhsv%no3 = (-0.8_rk*Denitrific + Nitri - uptNO3)

! Anammox: NH3 oxidation by nitrite, here related to NO3
  Anammox    = self%rAnammox * Anoxiclim *Rescale * env%nh3 * no3/(no3+self%ksNO3denit)

! preference for NH3 in DIN-uptake of autotrophs
  nh3f        = 1.0d0 - exp(-5*env%nh3/(nut%N+self%small))
!  dynamics of nh3 ~ dissolved ammonium
  rhsv%nh3    = Nprod - Nitri + lossZ%N * zoo%C  & !/ (1.0_rk + self%NH3Ads)
               + (exud%N - nh3f*uptake%N) * phy%C &!env%nh3/(nut%N+self%small) *
               + self%dil * (self%nh3_initial - env%nh3) &
               - Anammox

  rhsv%nutN   = rhsv%nutN - 0.8d0 * Denitrific - Anammox

!  dynamics of pdet ~ detritus-P
!    rhsv%pdet = (radsP - f_T * Pprod) 
!  dynamics of po4 ~ dissolved phosphate
!  rhsv%po4 = (f_T * Pprod - radsP) 
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
if (self%BioOxyOn) then
      _SET_ODE_(self%id_nh3, rhsv%nh3 UNIT)
      _SET_ODE_(self%id_oxy, rhsv%oxy UNIT)
      _SET_ODE_(self%id_odu, rhsv%odu UNIT)
end if
if (self%NResOn) then
      _SET_ODE_(self%id_RNit, rhsv%RNit UNIT)
end if
!#E_ODE

!________________________________________________________________________________
! set diag variables, mostly from PrimProd module

!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

if (self%DebugDiagOn) then
!#S_DIA
  _SET_DIAGNOSTIC_(self%id_GPPR, _REPLNAN_(phy%gpp*phy%C))   !average gross primary production
  _SET_DIAGNOSTIC_(self%id_Denitr, _REPLNAN_(0.8*Denitrific)) !average denitrification rate
  _SET_DIAGNOSTIC_(self%id_chl2C, _REPLNAN_(phy%theta*phy%rel_chloropl/12)) !average chlorophyll:carbon ratio 
  _SET_DIAGNOSTIC_(self%id_fracR, _REPLNAN_(phy%frac%Rub))   !average Rubisco fract
  _SET_DIAGNOSTIC_(self%id_fracT, _REPLNAN_(phy%frac%theta)) !average LHC fract
  _SET_DIAGNOSTIC_(self%id_Theta, _REPLNAN_(phy%theta))      !average Theta
  _SET_DIAGNOSTIC_(self%id_fracNU, _REPLNAN_(phy%frac%NutUpt)) !average Nut
  _SET_DIAGNOSTIC_(self%id_QN, _REPLNAN_(phy%Q%N))           !average N:C ratio
  _SET_DIAGNOSTIC_(self%id_QP, _REPLNAN_(phy%Q%P))           !average P:C ratio
  _SET_DIAGNOSTIC_(self%id_aVN, _REPLNAN_(acclim%aV%N))      !average N-uptake activity
  _SET_DIAGNOSTIC_(self%id_aVP, _REPLNAN_(acclim%aV%P))      !average P-uptake activity
  _SET_DIAGNOSTIC_(self%id_aVSi, _REPLNAN_(acclim%aV%Si))    !average Si-uptake activity
  _SET_DIAGNOSTIC_(self%id_faN, _REPLNAN_(acclim%fA%N))      !average N-uptake affinity allocation
  _SET_DIAGNOSTIC_(self%id_faP, _REPLNAN_(acclim%fA%P))      !average P-uptake affinity allocation
  _SET_DIAGNOSTIC_(self%id_faSi, _REPLNAN_(acclim%fA%Si))    !average Si-uptake affinity allocation
  _SET_DIAGNOSTIC_(self%id_rQSi, _REPLNAN_(phy%relQ%Si))     !average Relative Si-Quota
  _SET_DIAGNOSTIC_(self%id_tmp, _REPLNAN_(acclim%tmp))       !average Temporary diagnostic
  _SET_DIAGNOSTIC_(self%id_fac1, _REPLNAN_(dRchl_phyC_dt))   !average Auxiliary diagnostic 
  _SET_DIAGNOSTIC_(self%id_fac2, _REPLNAN_(acclim%dRchl_dfracR*acclim%dfracR_dt)) !average Auxiliary diagnostic
  _SET_DIAGNOSTIC_(self%id_fac3, _REPLNAN_(acclim%dRchl_dtheta*acclim%dtheta_dt)) !average Auxiliary diagnostic
  _SET_DIAGNOSTIC_(self%id_fac4, _REPLNAN_(acclim%fac1))     !average dtheta
  _SET_DIAGNOSTIC_(self%id_fac5, _REPLNAN_(acclim%fac2))     !average dtheta
  _SET_DIAGNOSTIC_(self%id_dPAR, _REPLNAN_(env%par))         !average Photosynthetically Active Radiation
  _SET_DIAGNOSTIC_(self%id_phyUR, _REPLNAN_(uptake%C))       !average Phytoplankton C Uptake Rate
  _SET_DIAGNOSTIC_(self%id_phyELR, _REPLNAN_(-exud%C))       !average Phytoplankton Exudation Loss Rate
  _SET_DIAGNOSTIC_(self%id_phyALR, _REPLNAN_(-aggreg_rate))  !average Phytoplankton Aggregation Loss Rate
  _SET_DIAGNOSTIC_(self%id_phyGLR, _REPLNAN_(-graz_rate/phy%reg%C)) !average Phytoplankton Grazing Loss Rate
  _SET_DIAGNOSTIC_(self%id_vsinkr, _REPLNAN_(exp(-self%sink_phys*phy%relQ%N*phy%relQ%P))) !average Relative Sinking Velocity
  _SET_DIAGNOSTIC_(self%id_qualPOM, _REPLNAN_(qualPOM))      !average Quality of POM 
  _SET_DIAGNOSTIC_(self%id_qualDOM, _REPLNAN_(qualDOM))      !average Quality of DOM 
  _SET_DIAGNOSTIC_(self%id_no3, _REPLNAN_(no3))              !average Nitrate
!#E_DIA
end if
!write (*,'(A,3(F11.5))') 'RN,depo,denit=',env%RNit,deporate,denitrate
                  
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
!  call sinking(self%vS_phy, phyQstat, vsink)

   !SINKING AS A FUNCTION OF INTERNAL STATES
   vsink = self%vS_phy * exp( -self%sink_phys * phyQstat)

   !CONSTANT SINKING
   !vsink = self%vS_phy
   
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

subroutine maecs_do_surface(self,_ARGUMENTS_DO_SURFACE_)
   use maecs_functions
   
   class (type_hzg_maecs), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: tot_vm_C,tot_vm_N,tot_vm_P,tot_vm_S, O2flux,O2airbl,oxy,tot_vm_GPPR,tot_vm_Denitr
   
!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

   _HORIZONTAL_LOOP_BEGIN_
      _GET_HORIZONTAL_(self%id_totN_vertint,tot_vm_N)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totN_vertint_diag,_REPLNAN_(tot_vm_N))
      _GET_HORIZONTAL_(self%id_totC_vertint,tot_vm_C)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totC_vertint_diag,_REPLNAN_(tot_vm_C))
      if (self%DiagOn) then
        _GET_HORIZONTAL_(self%id_GPPR_vertint,tot_vm_GPPR)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_GPPR_vertint_diag,_REPLNAN_(tot_vm_GPPR))
      end if
      if (self%BioOxyOn) then
        _GET_HORIZONTAL_(self%id_Denitr_vertint,tot_vm_Denitr)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Denitr_vertint_diag, _REPLNAN_(tot_vm_Denitr))
      end if
     
      if (self%PhosphorusOn) then
         _GET_HORIZONTAL_(self%id_totP_vertint,tot_vm_P)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totP_vertint_diag,_REPLNAN_(tot_vm_P))
! --- atmospheric deposition of PO4
         _SET_SURFACE_EXCHANGE_(self%id_nutP, self%P_depo UNIT)
      end if
      if (self%SiliconOn) then
         _GET_HORIZONTAL_(self%id_totS_vertint,tot_vm_S)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totS_vertint_diag,_REPLNAN_(tot_vm_S))
      end if

! --- wet and dry deposition of NO3 
      _SET_SURFACE_EXCHANGE_(self%id_nutN, self%N_depo UNIT)

! --- oxygen flux between sea water and air -----
      if (self%BioOxyOn) then
! O2 flux across the boundary layer
! O2airbl is the saturation concentration of O2
! airsea_ex is the average diffusivity coefficient (m2/sec) divided by the thickness of the boundary layer.
! for O2 in mmol m-3, the rate of exchange in mmol m-2 s-1).
! Positive values imply a flux into the water, negative: out of the water. 
         O2airbl = self%O2_sat
!        _GET_HORIZONTAL_(self%id_O2airbl, O2airbl)! boundary layer dissolved oxygen in mmolO2/m**3
        _GET_(self%id_oxy, oxy)   ! sea water dissolved oxygen in mmolO2/m**3

        O2flux  = self%ex_airsea * (O2airbl - oxy)!
        _SET_SURFACE_EXCHANGE_(self%id_oxy, O2flux )
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_O2flux_diag, _REPLNAN_(O2flux)) ! converts mmol/m2.s to mmol/m2.d
      endif

   _HORIZONTAL_LOOP_END_

end subroutine maecs_do_surface
! potential entries in maecs_deps.lst;but might work only using GOTM-input scheme
! O2airbl	mmol-C/m**3 horizontal_dependency O2airbl surface_molecular_oxygen_partial_pressure_difference_between_sea_water_and_air #BioOxyOn
!N2air	mmol-C/m**3 horizontal_dependency N2air mole_concentration_of_atomic_nitrogen_in_air  
!N2flux	mmol-N/m**2/d horizontal_diagnostic_variable  'N2flux','mmol-N/m**2/d','nitrogen_flux_between_sea_water_and_air' 
  