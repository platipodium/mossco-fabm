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
integer  :: i, j, iz, ihour, iloop, doy
real(rk) :: reminT, degradT       ! Temp dependent remineralisation and hydrolysis rates
! --- QUOTA and FRACTIONS
real(rk) :: phys_status, dQN_dt, dRchl_phyC_dt=0.0_rk ! []

! --- ZOOPLANKTON GRAZING, INGESTION, MORTALITY, RESPIRATION... 
real(rk) :: graz_rate   ! carbon-specific grazing rate                          [d^{-1}]
real(rk) :: zoo_respC   ! temperature dependent carbon respiration rate  [mmolC m^{-3} d^{-1}]
real(rk) :: zoo_mort
real(rk) :: decay       ! pigment-specific decay rate                          [d^{-1}]
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
real(rk) :: prodO2, rhochl, uptNH4, uptNO3, uptchl, uptN, respphyto,faeces, min_Cmass
logical  :: IsCritical = .false. ! phyC and phyN below reasonable range ?
#define _KAI_ 0
#define _MARKUS_ 1
#define _DEBUG_ 0
! #define UNIT / 86400
#define UNIT *1.1574074074E-5_rk

 _LOOP_BEGIN_

#if _DEBUG_
write(*,'(A)') 'begin DO'
#endif

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

!S_GED
  _GET_(self%id_temp, env%temp)  ! water temperature
  _GET_(self%id_par, env%par)  ! light photosynthetically active radiation
  if (self%ChemostatOn) then
    if (_AVAILABLE_(self%id_CO2)) then
      _GET_(self%id_CO2, env%CO2)  ! CO2
    else
      env%CO2 = 0.0_rk ! todo: throw an error, if necessary dependency cannot be found 
    end if
  end if

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
call min_mass(self,phy, min_Cmass, IsCritical, method=2) ! minimal reasonable Phy-C and -Nitrogen

if(self%maxVal .lt. 0.0d0) then 
  IsCritical=.false.
else
  if(phy%chl .gt. self%maxVal .or. phy%Rub .gt. self%maxVal) IsCritical=.true.
endif

if(IsCritical .and. .not. self%ChemostatOn) then
  rhsv%nutN=0.0d0
  rhsv%nutP=0.0d0
  rhsv%nutS=0.0d0
  rhsv%phyC=0.0d0
  rhsv%phyN=0.0d0
  rhsv%phyP=0.0d0
  rhsv%phyS=0.0d0
  rhsv%zooC=0.0d0
  rhsv%detC=0.0d0
  rhsv%detN=0.0d0
  rhsv%detP=0.0d0
  rhsv%detS=0.0d0
  rhsv%domC=0.0d0
  rhsv%domN=0.0d0
  rhsv%domP=0.0d0
  rhsv%RNit=0.0d0
  rhsv%Rub=0.0d0
  rhsv%chl=0.0d0
  rhsv%nh3=0.0d0
  rhsv%oxy=0.0d0
  rhsv%odu=0.0d0
else
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

!if (IsCritical .and. .false. ) then
!  write (*,'(A,4(F10.3))') 'fR=',phy%reg%C,phy%frac%Rub,self%small_finite + self%rel_chloropl_min,phy%frac%NutUpt
!end if
!if (phy%chl .lt. 0.01d0) then
!  phy%theta     = phy%chl / (phy%rel_chloropl * phy%reg%C)   ! trait variable
!  phy%frac%theta= phy%theta * phy%rel_chloropl * maecs%itheta_max ! []     no 
!  write (*,'(A,3(F10.3))') 'fT=',phy%chl,phy%reg%C,phy%Rub end if
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
call photosynthesis(self,sens,phy,nut,uptake,exud,acclim)

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
  zoo_mort    = self%mort_zoo * sens%f_T**self%fT_exp_mort  * zoo%C

else
  graz_rate   = 0.0_rk
!  if (self%ChemostatOn .and. .not. IsCritical) graz_rate = 0.2*phy%C
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

! additional mortality due to H2S stress (EC_50 :55 mmol-H2S/m3; Kuester2005  or for Skel cost 3 mmol-H2S/m3, Breteler1991)
if (self%BioOxyOn) then
   aggreg_rate = aggreg_rate + self%mort_ODU* env%odu
endif

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

!if (abs(phy%C) .gt. 1d-4) then
 if (self%PhotoacclimOn ) then ! check for too small biomasses %chl

! PHYTOPLANKTON CHLa
     ! note that theta*rel_chloropl in units [mg Chla (mmol C)^{-1}] 
   dQN_dt        = (rhsv%phyN * phy%reg%C - rhsv%phyC * phy%reg%N) / (phy%reg%C*phy%reg%C)
! TODO: dangerous to work with RHS instead of net uptake rates (mortality has no physiological effect)

   dRchl_phyC_dt =  acclim%dRchl_dtheta * acclim%dtheta_dt   & 
                  + acclim%dRchl_dfracR * acclim%dfracR_dt   & 
                  + acclim%dRchl_dQN    * dQN_dt 
! pigment decay to relieve from artificially high pigm:C ratios at very low phyC
   decay = self%decay_pigm * (exp(phy%frac%theta)-1.0d0)

!   rhsv%chl = phy%theta * phy%frac%Rub * phy%relQ%N**self%sigma * rhsv%phyC + dRchl_phyC_dt * phy%C
   rhsv%chl = phy%theta * phy%frac%Rub * phy%relQ%N**self%sigma * rhsv%phyC &
                   + dRchl_phyC_dt * phy%reg%C - decay * phy%chl

!if (rhsv%chl .lt. -200.d0) then
!  phy%theta     = phy%chl / (phy%rel_chloropl * phy%reg%C)   ! trait variable
!  phy%frac%theta= phy%theta * phy%rel_chloropl * maecs%itheta_max ! []     no 
!  write (*,'(A,4(F14.3))') 'dChl=',rhsv%chl,phy%chl,dRchl_phyC_dt,rhsv%phyC,dQN_dt
!  write (*,'(A,6(F14.3))') 'aa=',acclim%dRchl_dtheta,acclim%dtheta_dt,acclim%dRchl_dfracR,acclim%dfracR_dt,acclim%dRchl_dQN, dQN_dt
!  write (*,'(A,2(F14.3))') 'dmdx=',acclim%fac1,acclim%fac2
!end if

!write (*,'(A,4(F10.3))') 'rhs chl=', phy%theta * phy%frac%Rub * phy%relQ%N**self%sigma * rhsv%phyC,dRchl_phyC_dt * phy%reg%C*1E1,phy%relQ%N**self%sigma,phy%theta

!_____________________________________________ _________________________________
 end if ! PhotoacclimOn

 if (self%RubiscoOn) then 
   decay = self%decay_pigm * (exp(phy%frac%Rub)-1.0d0)
   rhsv%Rub  = phy%Rub/phy%reg%C * rhsv%phyC + acclim%dfracR_dt * phy%C - decay*phy%Rub
 end if 
!else  rhsv%Rub  = 0.0d0  rhsv%chl  = 0.0d0 
!endif !if (abs(phy%C) .gt. 1d-4)

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
!   PO4-adsorption ceases when critical capacity is reached
!   [FeS] approximated by ODU
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
  nh3f        = 1.0d0 - exp(-5*env%nh3/smooth_small(nut%N,self%small))
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
end if !BioOxyOn

end if !IsCritical

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

! artifical, serial nutrient input to illustrate co-limitation dynamics in 0D
!if (self%ChemostatOn) then
!  _GET_GLOBAL_ (self%id_doy,doy) !day of year
!  select case (doy)
!           case (89:92)
!            _SET_ODE_(self%id_nutP, 5*(1.0-nut%P) UNIT)
!           case (29:32)
!            _SET_ODE_(self%id_nutN, 5*(16.0-nut%N) UNIT)
!           case (59:62)
!            _SET_ODE_(self%id_nutS, 5*(16.0-nut%Si) UNIT)
!  end select
!endif

!________________________________________________________________________________
! set diag variables, mostly from PrimProd module

!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

!#S_DIA
if (self%DebugDiagOn) then
  _SET_DIAGNOSTIC_(self%id_tmp, _REPLNAN_(acclim%tmp))       !average Temporary_diagnostic_
  _SET_DIAGNOSTIC_(self%id_fac1, _REPLNAN_(dRchl_phyC_dt))   !average Auxiliary_diagnostic_
  _SET_DIAGNOSTIC_(self%id_fac2, _REPLNAN_(acclim%dRchl_dfracR*acclim%dfracR_dt)) !average Auxiliary_diagnostic_
  _SET_DIAGNOSTIC_(self%id_fac3, _REPLNAN_(acclim%dRchl_dtheta*acclim%dtheta_dt)) !average Auxiliary_diagnostic_
  _SET_DIAGNOSTIC_(self%id_fac4, _REPLNAN_(acclim%fac1))     !average dtheta_dt_due_to_flex_theta_
  _SET_DIAGNOSTIC_(self%id_fac5, _REPLNAN_(acclim%fac2))     !average dtheta_dt_due_to_grad_theta_
end if
if (self%BGC0DDiagOn) then
  _SET_DIAGNOSTIC_(self%id_GPPR, _REPLNAN_(phy%gpp*phy%C))   !average gross_primary_production_
  _SET_DIAGNOSTIC_(self%id_Denitr, _REPLNAN_(0.8*Denitrific)) !average denitrification_rate_
  _SET_DIAGNOSTIC_(self%id_dPAR, _REPLNAN_(env%par))         !average Photosynthetically_Active_Radiation_
  _SET_DIAGNOSTIC_(self%id_DNP, _REPLNAN_(nut%N/(nut%P+self%small)))   !average DIN:DIP_ratio_
  _SET_DIAGNOSTIC_(self%id_QNP, _REPLNAN_(phy%Q%N/phy%Q%P))  !average N:P_ratio_
  _SET_DIAGNOSTIC_(self%id_qualPOM, _REPLNAN_(qualPOM))      !average Quality_of_POM_
  _SET_DIAGNOSTIC_(self%id_qualDOM, _REPLNAN_(qualDOM))      !average Quality_of_DOM_
  _SET_DIAGNOSTIC_(self%id_no3, _REPLNAN_(no3))              !average Nitrate_
end if
if (self%PhysiolDiagOn) then
  _SET_DIAGNOSTIC_(self%id_chl2C, _REPLNAN_(phy%theta*phy%rel_chloropl/12)) !average chlorophyll:carbon_ratio_=_chl-a/chloroplast-C_*_chloroplast-C/phy-molC_*_1molC/12gC_
  _SET_DIAGNOSTIC_(self%id_Theta, _REPLNAN_(phy%theta))      !average Theta_
  _SET_DIAGNOSTIC_(self%id_fracR, _REPLNAN_(phy%frac%Rub))   !average Rubisco_fract._allocation_
  _SET_DIAGNOSTIC_(self%id_fracT, _REPLNAN_(phy%frac%theta)) !average LHC_fract._allocation_
  _SET_DIAGNOSTIC_(self%id_fracNU, _REPLNAN_(phy%frac%NutUpt)) !average Nut._Uptake_fract._allocation_
  _SET_DIAGNOSTIC_(self%id_QN, _REPLNAN_(phy%Q%N))           !average N:C_ratio_
  _SET_DIAGNOSTIC_(self%id_QP, _REPLNAN_(phy%Q%P))           !average P:C_ratio_
  _SET_DIAGNOSTIC_(self%id_QSi, _REPLNAN_(phy%Q%Si))         !average Si:C_ratio_
  _SET_DIAGNOSTIC_(self%id_aVN, _REPLNAN_(acclim%aV%N))      !average N-uptake_activity_
  _SET_DIAGNOSTIC_(self%id_aVP, _REPLNAN_(acclim%aV%P))      !average P-uptake_activity_
  _SET_DIAGNOSTIC_(self%id_aVSi, _REPLNAN_(acclim%aV%Si))    !average Si-uptake_activity_
  _SET_DIAGNOSTIC_(self%id_faN, _REPLNAN_(acclim%fA%N))      !average N-uptake_affinity_allocation_
  _SET_DIAGNOSTIC_(self%id_faP, _REPLNAN_(acclim%fA%P))      !average P-uptake_affinity_allocation_
  _SET_DIAGNOSTIC_(self%id_faSi, _REPLNAN_(acclim%fA%Si))    !average Si-uptake_affinity_allocation_
  _SET_DIAGNOSTIC_(self%id_rQN, _REPLNAN_(phy%relQ%N))       !average Relative_N-Quota_
  _SET_DIAGNOSTIC_(self%id_rQP, _REPLNAN_(phy%relQ%P))       !average Relative_P-Quota_
  _SET_DIAGNOSTIC_(self%id_rQSi, _REPLNAN_(phy%relQ%Si))     !average Relative_Si-Quota_
end if
if (self%RateDiagOn) then
  _SET_DIAGNOSTIC_(self%id_phyUR, _REPLNAN_(uptake%C))       !average Phytoplankton_C_Uptake_Rate_
  _SET_DIAGNOSTIC_(self%id_phyELR, _REPLNAN_(-exud%C))       !average Phytoplankton_Exudation_Loss_Rate_
  _SET_DIAGNOSTIC_(self%id_phyALR, _REPLNAN_(-aggreg_rate))  !average Phytoplankton_Aggregation_Loss_Rate_
  _SET_DIAGNOSTIC_(self%id_phyGLR, _REPLNAN_(-graz_rate/phy%reg%C)) !average Phytoplankton_Grazing_Loss_Rate_
!  _SET_DIAGNOSTIC_(self%id_vsinkr, _REPLNAN_(exp(-self%sink_phys*phy%relQ%N*phy%relQ%P))) !average Relative_Sinking_Rate_
end if
!#E_DIA


#if _DEBUG_
write(*,'(A)') 'end DO'
#endif

  _LOOP_END_

end subroutine maecs_do

!> @brief handles vertical movement for depth-varying movement rates
!> @details phyto sinking rate depends on the nutritional state, so for each node:
!! \n \f$ phy\%relQ \f$ obtained by calling calc_internal_states(self,phy,det,dom,zoo) 
!! \n then \f$ phyQstat=phy\%relQ\%N * phy\%relQ\%P \f$
!! \n finally, vs_phy = maecs_functions::sinking(self\%vS_phy, phyQstat, vs_phy)
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
logical  :: IsCritical = .false. ! phyC and phyN below reasonable range ?

!
! !LOCAL VARIABLES: 
REALTYPE    :: phyQstat,vs_phy,vs_det, phyEner, minPigm,min_Cmass, minc
!REALTYPE    :: aggf, agge=16.d0
REALTYPE, parameter :: secs_pr_day = 86400.d0 
!EOP
!-----------------------------------------------------------------------
!BOC

_FABM_LOOP_BEGIN_

#if _DEBUG_
write(*,'(A)') 'begin vert_move'
#endif

   ! Retrieve phtoplankton state
   
   !Retrieve the 'phyQstat' directly as a diagnostic variable: does not work yet.
   !fabm_get_bulk_diagnostic_data(self%id_phyqstat,phyQstatD) !where, phyQstat=relQ%N*relQ%P
   !_GET_(self%id_phyqstat,phyQstatD)
   
   !Calculate manually
   _GET_(self%id_phyC, phy%C)  ! Phytplankton Carbon in mmol-C/m**3
   _GET_(self%id_phyN, phy%N)  ! Phytplankton Nitrogen in mmol-N/m**3
!   _GET_(self%id_detC, det%C)  ! Detritus Nitrogen in mmol-N/m**3
!   _GET_(self%id_detN, det%N)  ! Detritus Nitrogen in mmol-N/m**3
   _GET_(self%id_domC, dom%C)  ! DONitrogen in mmol-N/m**3
!    aggf = det%C/106+det%N/16
!    if (self%PhosphorusOn) then
!      _GET_(self%id_detP, det%P)  ! Detritus Phosphorus in mmol-P/m**3
!      aggf = aggf + det%P
!    endif
!   aggf = 1.0_rk + 2*self%phi_agg * (phy%N + det%N) 
!   aggf = 0.1d0 + 1.0d0/(1.0d0+ exp(-3+agge*aggf))

   if (self%GrazingOn) then
     _GET_(self%id_zooC, zoo%C)  ! Zooplankton Carbon in mmol-C/m**3
   end if
   if (self%PhosphorusOn) then
     _GET_(self%id_phyP, phy%P)  ! Phytplankton Phosphorus in mmol-P/m**3
   end if

   !write (*,'(A,2(F10.3))') 'Before: phy%C, phy%N=', phy%C, phy%N
   call min_mass(self,phy, min_Cmass, IsCritical, method=2) 
   !write (*,'(A,2(F10.3))') 'After: phy%C, phy%N=', phy%C, phy%N
   call calc_internal_states(self,phy,det,dom,zoo) 
   !write (*,'(A,2(F10.3))') 'phy%relQ%N, phy%relQ%P=', phy%relQ%N, phy%relQ%P
   
   ! nutrient limitation ; TODO check product rule and add other elements such as Si
   phyQstat = phy%relQ%N * phy%relQ%P 

   ! energy limitation ; TODO check function and quantity
!   phyEner  = phy%gpp / (self%V_NC_max*self%zeta_CN)
! smoothed minimum of energy and nutrient limitation; 
!   phyQstat = phyQstat - smooth_small(phyQstat - phyEner, self%small)

   ! Calculate sinking
!  call sinking(self%vS_phy, phyQstat, vs_phy)

   !SINKING AS A FUNCTION OF INTERNAL STATES
   vs_phy = -self%vS_phy * exp( -self%sink_phys * phyQstat)
   if (self%RateDiagOn) then 
      _SET_DIAGNOSTIC_(self%id_vsinkr, _REPLNAN_(-vs_phy)) !average Relative Sinking Velocity
   end if

   !CONSTANT SINKING
   !vs_phy = self%vS_phy
   
   vs_phy = vs_phy / secs_pr_day
   !write (*,'(A,2(F10.3))') 'phyQstat, vs_phy=', phyQstat, vs_phy
!   vs_det = -self%vS_det*aggf/secs_pr_day
   vs_det = -1.0_rk*self%vS_det/secs_pr_day
   !set the rates
   _SET_VERTICAL_MOVEMENT_(self%id_detC,vs_det)
   _SET_VERTICAL_MOVEMENT_(self%id_detN,vs_det)
   _SET_VERTICAL_MOVEMENT_(self%id_phyN,vs_phy)
   _SET_VERTICAL_MOVEMENT_(self%id_phyC,vs_phy)
!   if (self%ZooSinkMeth .eq. 1) then
!    _SET_VERTICAL_MOVEMENT_(self%id_zooC,vs_phy)
!   endif
   if (self%PhosphorusOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_phyP,vs_phy)
      _SET_VERTICAL_MOVEMENT_(self%id_detP,vs_det)
   end if
   if (self%SiliconOn) then
      _SET_VERTICAL_MOVEMENT_(self%id_phyS,vs_phy)
      _SET_VERTICAL_MOVEMENT_(self%id_detS,vs_det)
   end if
   if (self%PhotoacclimOn) then 
      _SET_VERTICAL_MOVEMENT_(self%id_chl, vs_phy)
      _SET_VERTICAL_MOVEMENT_(self%id_Rub, vs_phy)
   end if

#if _DEBUG_
write(*,'(A)') 'end vert_move'
#endif

_FABM_LOOP_END_
  
end subroutine maecs_get_vertical_movement

subroutine maecs_do_surface(self,_ARGUMENTS_DO_SURFACE_)
   use maecs_functions
   
   class (type_hzg_maecs), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: tot_vi_C,tot_vi_N,tot_vi_P,tot_vi_S, O2flux,O2airbl,oxy,tot_vi_GPPR,tot_vi_Denitr
   
!define _REPLNAN_(X) X !changes back to original code
#define _REPLNAN_(X) nan_num(X)

   _HORIZONTAL_LOOP_BEGIN_

#if _DEBUG_
write(*,'(A)') 'begin surface_DO'
#endif

      if (self%BGC2DDiagOn) then
        _GET_HORIZONTAL_(self%id_GPPR_vertint,tot_vi_GPPR)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_GPPR_vertint_diag,_REPLNAN_(tot_vi_GPPR))
        if (self%BioOxyOn) then
          _GET_HORIZONTAL_(self%id_Denitr_vertint,tot_vi_Denitr)
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Denitr_vertint_diag, _REPLNAN_(tot_vi_Denitr))
        end if
      end if
      
      if (self%Budget2DDiagOn) then 
      _GET_HORIZONTAL_(self%id_totN_vertint,tot_vi_N)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totN_vertint_diag,_REPLNAN_(tot_vi_N))
      _GET_HORIZONTAL_(self%id_totC_vertint,tot_vi_C)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totC_vertint_diag,_REPLNAN_(tot_vi_C))
      if (self%PhosphorusOn) then
         _GET_HORIZONTAL_(self%id_totP_vertint,tot_vi_P)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totP_vertint_diag,_REPLNAN_(tot_vi_P))
      end if
      if (self%SiliconOn) then
         _GET_HORIZONTAL_(self%id_totS_vertint,tot_vi_S)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_totS_vertint_diag,_REPLNAN_(tot_vi_S))
      end if
      end if

! --- wet and dry deposition of NO3 
      _SET_SURFACE_EXCHANGE_(self%id_nutN, self%N_depo UNIT)
! --- atmospheric deposition of PO4
      if (self%PhosphorusOn) then
         _SET_SURFACE_EXCHANGE_(self%id_nutP, self%P_depo UNIT)
      end if

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
        if (self%BGC2DDiagOn) then
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_O2flux_diag, _REPLNAN_(O2flux)) ! converts mmol/m2.s to mmol/m2.d
        end if
      endif

#if _DEBUG_
write(*,'(A)') 'end surface_DO'
#endif

   _HORIZONTAL_LOOP_END_

end subroutine maecs_do_surface
! potential entries in maecs_deps.lst;but might work only using GOTM-input scheme
! O2airbl	mmol-C/m**3 horizontal_dependency O2airbl surface_molecular_oxygen_partial_pressure_difference_between_sea_water_and_air #BioOxyOn
!N2air	mmol-C/m**3 horizontal_dependency N2air mole_concentration_of_atomic_nitrogen_in_air  
!N2flux	mmol-N/m**2/d horizontal_diagnostic_variable  'N2flux','mmol-N/m**2/d','nitrogen_flux_between_sea_water_and_air' 
  