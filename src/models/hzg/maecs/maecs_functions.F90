#include "fabm_driver.h"
!---------------------------------------------------------
! !MODULE: MAECS_functions --- more to come
!          Model for Adaptive Ecosystems in Coastal Seas 
   module maecs_functions

! !USES:
   use fabm_types
   use maecs_types

   private
   public   uptflex       ,& 
            queuefunc, queuefunc0, queuederiv ,&
            smooth_small, sinking, min_mass, calc_rel_chloropl, &
            calc_sensitivities, calc_internal_states

 contains  

!---------------------------------------------------------
pure REALTYPE function smooth_small(x, eps)
!
!  continous smoothing function by kw Apr 2012
!  avoids step-like shifts at small lower boundaries
   implicit none
   REALTYPE, intent(in)          :: x, eps
   REALTYPE                      :: arg, larger
!--------------------------------------------------------------
   
   if (x .lt. 1.5d1*eps) then
     arg     = x/eps
     larger  = exp(arg)
     smooth_small  = eps + (x - eps) * larger / (larger + 1.0d0) 
   else
    smooth_small  = x
   endif
   end function smooth_small 

!-----------------------------------------------------------------------
pure REALTYPE function uptflex(Aff0, Vmax0, Nut, NutMin)

! !DESCRIPTION:
!
! !USES:
   implicit none
! !INPUT PARAMETERS:
   REALTYPE, intent(in)      :: Aff0, Vmax0, Nut, NutMin
   REALTYPE                  :: Aff,Vmax,fAv, Nutp

! optimal partitioning
   Nutp    = smooth_small(Nut, NutMin)
!   Nutp    = Nut
   fAv     = _ONE_/(sqrt(Aff0*Nutp/(Vmax0)) + _ONE_ );
! uptake regulation: sites vs. processing
   Aff     = fAv * Aff0    ! nutrient affinity 
   Vmax    = (_ONE_-fAv) * Vmax0   ! maximal N-uptake rate 

   uptflex = Vmax*Aff*Nut/(Aff*Nutp + Vmax)

   end function uptflex

!-----------------------------------------------------------------------
pure REALTYPE function queuefunc0(n,x)
!
!  response function from queing theory
!  synchrony of processing n->inf :liebig  n~1:product
!
! !USES:
   implicit none
! !INPUT PARAMETERS:
   REALTYPE, intent(in)          :: x, n
   REALTYPE                      :: px

   px = x**(n+_ONE_)
   if(abs(_ONE_-px) .lt. 1E-5) then
      queuefunc0 =  n/(n+_ONE_)
   else
      queuefunc0 =  (x-px)/(_ONE_-px)
   endif
   end function queuefunc0

!-----------------------------------------------
subroutine queuefunc(n,x,qfunc,qderiv)
!
!  response function from queing theory:  numerical approximation
!  synchrony of processing n->inf :liebig  n~1:product
!
   implicit none
   REALTYPE, intent(in)       :: x, n
   REALTYPE, intent(out)      :: qfunc, qderiv
   REALTYPE                   :: nn, hh,x0,en

   nn = n+1.
   hh = 1./nn
   x0 = (log(exp(nn)-1))*hh
   en = exp(-nn*(x/(1+hh*x) - x0))
   qfunc =  1. - hh*log(1.+ en)
   qderiv = 1. / ( (1. + hh * x)**2 * (1.+1./en))
   end subroutine queuefunc

!-------------------------------------------------------------
REALTYPE function queuederiv(n,x)
!
!  response function from queing theory
!  synchrony of processing n->inf :liebig  n~1:product
!
   implicit none
   REALTYPE, intent(in)          :: x, n
   REALTYPE                      :: nqueue_aa
   REALTYPE                      :: nqueue_bb = 1.3863d0
!-----------------------------------------------------------------------
   nqueue_aa = 6.9315d-1*(1.d0+1.0d0/(n))
   queuederiv = exp(-nqueue_aa*(x)**nqueue_bb) 
   end function queuederiv

!-------------------------------------------------------------
subroutine sinking(vS ,phys_status,sinkvel)
   implicit none
   REALTYPE, intent(in)     :: vS ,phys_status
   REALTYPE, intent(out)    :: sinkvel

   !   REALTYPE         :: 
!-----------------------------------------------!----------------------------------------------------------------
!------         sinking                --------------------------
!----------------------------------------------------------------
! ftempvisc = arrhenius(1.34,Temp); /* temperature dependence of viscosity */
! sv0  = ftempvisc*sink0 * STOKES0 ;
! assume proportionality between SPM and TKE */ 
! C=C0*exp(c*z) steady state: c=v_s/K_z
! critical shear stress tau TKE (Gordon & Dohne, 1973; Kim, Friedrichs, Maa, & Wright, 2000; Soulsby, 1983; Stapleton & Huntley, 1995)
! deposition rate given by Einstein and Krone [1962] ~ 1-tau;erosion rate by Partheniades [1965]  tau -1:*/
!-----------------------------------------------*/
!   density regulation: energy requirenment     */
!-----------------------------------------------*/
!!rest = smooth_small(resp,self%smallest);
!!px   = smooth_small(grossC_phy-rest)/rest,self%smallest);

!!relprod = 1.0d0 - exp(-px);
!-----------------------------------------------*/
!   density regulation: nutrient requirenment   */
!-----------------------------------------------*/
!qx=(nq-NQ0[i])/dpdQsMax;!dpdQsMax 0.003
!!relquot = 1.-exp(-rat_qN*rat_qP);
! Liebig approach for 'vitality' limitation ; TODO: change to queue*/
!!pv = relprod * relquot;

!rv0 = (2-exdens_exp)*lsz;ex = exp(1.0d0-lsz);
!iex = 1./(1.+ex);pv0 = 1.+3*exdens_exp*(lsz-l0)*iex;
!!rv = pv*pv0;
!--------------------------------------------------*/
!   effective sinking velocity:
!     Stokes - passive+active dens down-regulation   
!--------------------------------------------------*/
!!svel = sv0 * exp(rv0-rv);
!dsink_dQ = -svel*relprod*pv0*exp(-qx)/q0m;
!--------------------------------------------------------*/
! resuspension due higher TKE/K_z in shallow water reduces net sinking                                  
!   sres = resus0*tkef[b]*tkef[b];

! sink = fac*(svel-sres);  
!  dsink_ds = fac*svel*(2-exdens_exp-pv*3*exdens_exp*((1.+ex)-(lsz-l0)*ex)*iex*iex);
!------------------------------------------------------
!     stickiness -> aggregation -> enhanced sinking    
! epsc = stick0*(1.*Detritus[b]+0.*PConc[0][b]+1.*PConc[4][b])*ftemp[0];
! px =  ftempvisc*epsc*(PConc[i][b]+1.*Detritus[b]);
! stick = sink*px;
!dsink_ds *= 1 + 1.*px; dsink_dQ *= 1 + 0.5*px;
!sink += stick; ! Maerz&Wirtz, ECSS2009 Fig6 

   sinkvel   = -vS * exp(-4*phys_status ) ! /secs_pr_day
   end subroutine sinking

#define _KAI_ 0
#define _MARKUS_ 1
!------------------------------------------------------
subroutine min_mass(maecs,phy,method)

implicit none

type (type_maecs_base_model), intent(in)      :: maecs
type (type_maecs_phy), intent(inout)   :: phy
integer, intent(in), optional          :: method

REALTYPE                         :: min_Cmass,min_Nmass
integer                          :: ischanged, mm_method=_MARKUS_

if (present(method)) mm_method=method
ischanged = 0 

select case (mm_method)
 case (_MARKUS_)
   if (phy%N .le. 1.d-7) then
      phy%N = 1.d-7 
      phy%C = phy%N / maecs%aver_QN_phy
      if (maecs%PhosphorusOn) then
         phy%P = phy%C * maecs%aver_QP_phy
      end if
   end if
   phy%N_reg = phy%N
   phy%C_reg = phy%C
   phy%P_reg = phy%P

 case (_KAI_)
! -------------------------------------------------------------------------------
! --- Here, quota are smoothed as soon as phytoplankton biomass/biovolume 
!     approaches a lower threshold 'min_mass', in the order of o(1.d-4) / 
! set small boundary depending on numerical resolution
! TODO: insert h ~ level height, here 10cm 
   min_Cmass = maecs%small_finite * 1.0d-2 / maecs%a_spm  
   min_Nmass = min_Cmass * maecs%aver_QN_phy
   phy%C_reg = smooth_small( phy%C , min_Cmass)
   if (abs(phy%C-phy%C_reg) .gt. 1d-2*min_Cmass) then
      phy%N_reg = phy%C_reg * maecs%aver_QN_phy
      ischanged = 1 
   else 
      phy%N_reg = smooth_small( phy%N , min_Nmass)
      if (abs(phy%N-phy%N_reg) .gt. 1d-2*min_Nmass) then
         phy%C_reg = phy%N_reg / maecs%aver_QN_phy
         ischanged = 1 
      end if  
   end if
   if (ischanged==0) then ! retune P and Rub
      if (maecs%PhosphorusOn)  phy%P_reg =  phy%C_reg * maecs%aver_QP_phy
      if (maecs%RubiscoOn) then 
         phy%Rub =  phy%C_reg * maecs%frac_Rub_ini
!        phy%Rub = phy%N_reg * maecs%frac_Rub_ini
      end if
   else  ! additional check for Rub and P; TODO: omitt ??
       phy%Rub = smooth_small( phy%Rub , min_Cmass * maecs%frac_Rub_ini)
!       phy%Rub = smooth_small( phy%Rub , min_Cmass * maecs%aver_QN_phy* maecs%frac_Rub_ini)
      if (maecs%PhosphorusOn)  phy%P_reg =  smooth_small(phy%P,min_Cmass * maecs%aver_QP_phy)
   end if ! ischanged
 
 case (3)
! -------------------------------------------------------------------------------
! --- Here, quota are smoothed as soon as phytoplankton biomass/biovolume 
!     approaches a lower threshold 'min_mass', in the order of o(1.d-4) / 
! set small boundary depending on numerical resolution
! TODO: insert h ~ level height, here 10cm 
   min_Cmass = maecs%small_finite * 1.0d-3 / maecs%a_spm  
   phy%C_reg = smooth_small( phy%C , min_Cmass)
   phy%N_reg = smooth_small( phy%N , min_Cmass * maecs%aver_QN_phy)
   if (maecs%PhosphorusOn)  phy%P_reg = smooth_small( phy%P , min_Cmass * maecs%aver_QP_phy)

   if (maecs%RubiscoOn) then 
     phy%Rub = smooth_small( phy%Rub , min_Cmass * maecs%frac_Rub_ini)
!     phy%Rub = smooth_small( phy%Rub , min_Cmass * maecs%aver_QN_phy* maecs%frac_Rub_ini)
   end if    

end select
end subroutine min_mass

!------------------------------------------------------
subroutine calc_rel_chloropl(maecs,phy,method)

implicit none

type (type_maecs_base_model), intent(in)    :: maecs
type (type_maecs_phy), intent(inout) :: phy
integer, intent(in), optional        :: method
integer :: mm_method=_MARKUS_

if (present(method)) mm_method=method

select case (mm_method)
case (_MARKUS_)
   phy%rel_chloropl = maecs%rel_chloropl_min + phy%frac%Rub * phy%rel_QN**maecs%sigma
      
case (_KAI_)
!   phy%rel_chloropl = smooth_small(phy%frac%Rub * phy%rel_QN**maecs%sigma,maecs%rel_chloropl_min)
   phy%rel_chloropl = smooth_small(phy%frac%Rub * phy%rel_QN**maecs%sigma,maecs%rel_chloropl_min)

end select

end subroutine

!------------------------------------------------------
subroutine calc_sensitivities(maecs,sens,phy,env,nut)

implicit none
type (type_maecs_base_model), intent(in) :: maecs
type (type_maecs_sensitivities), intent(out) :: sens
type (type_maecs_phy),intent(in) :: phy
type (type_maecs_env),intent(in) :: env
type (type_maecs_om),intent(in) :: nut

REALTYPE :: par, T_Kelv

par          = maecs%frac_PAR * env%par ! use active  fraction frac_PAR of available radiation
T_Kelv       = env%Temp + 273.d0 ! temperature in Kelvin 
! ----------------------------------------------------------------------------
! +++ determine rates for photoautotophic growth ++++++++++++++++++++++++++++++++++++++++++++++++++
! --- temperature dependence of metabolic rates (with T_ref given in units [Kelvin]) --------------
sens%func_T  = exp(-maecs%AE_all *(1.0d0/T_Kelv - 1.0d0/maecs%T_ref ))

! --- (potential) maximum photosynthetic rate -----------------------------------------------------
sens%P_max_T = maecs%P_max * sens%func_T

! --- light response curve ------------------------------------------------------------------------
sens%a_light = maecs%alpha * par /(sens%P_max_T)    ! par NOCH BAUSTELLE, jetzt PAR in zellmitte 
sens%S_phot  = 1.0d0 - exp(- sens%a_light * phy%theta) ! [dimensionless]

! --- carbon specific N-uptake: sites vs. processing ----------------------------------
sens%up_NC   = uptflex(maecs%AffN,maecs%V_NC_max * sens%func_T,nut%N,maecs%small)

!  P-uptake coefficients
if (maecs%PhosphorusOn) then 
   sens%up_PC   = uptflex(maecs%AffP,maecs%V_PC_max * sens%func_T,nut%P,maecs%small)
end if
end subroutine

!------------------------------------------------------
subroutine calc_internal_states(self,phy,det,dom,zoo)

implicit none
type (type_maecs_base_model),intent(in)     :: self
type (type_maecs_phy), intent(inout) :: phy
type (type_maecs_om), intent(inout) :: det
type (type_maecs_om), intent(inout) :: dom
type (type_maecs_zoo), intent(inout) :: zoo
REALTYPE :: min_Cmass

min_Cmass = self%small_finite * 1.0d-3 / self%a_spm

! ------------------------------------------------------------------------------
!              calculate general quotas 
phy%QN    = phy%N_reg / phy%C_reg
! added for mixing effects in estuaries kw Jul, 15 2013
phy%QN  = smooth_small(phy%QN, self%QN_phy_0)

phy%frac%Rub=self%frac_Rub_ini

if (self%PhotoacclimOn) then
   if (self%RubiscoOn) then 
! trait + transporter needs division to become a trait again
     phy%frac%Rub = phy%Rub / phy%C_reg
!     phy%frac%Rub = phy%Rub / phy%N_reg
   end if    
end if

phy%frac%Rub = _ONE_ - smooth_small(_ONE_- phy%frac%Rub ,self%small_finite + self%rel_chloropl_min)

! --- stoichiometry of non-living organic matter  ---------------------------------
!dom%QN      = dom%N  /(dom%C + min_Cmass )  ! N:C ratio of dissolved organic matter (DOM)
!det%QN      = det%N /(det%C + min_Cmass)   ! N:C ratio of detritus
if (self%PhosphorusOn) then 
   phy%QP     = phy%P_reg / phy%C_reg
! added for mixing effects in estuaries kw Jul, 15 2013
   phy%QP  = smooth_small(phy%QP, self%QP_phy_0)

   phy%rel_QP = ( phy%QP - self%QP_phy_0 ) * self%iK_QP
! added for deep detritus traps with extreme quotas kw Jul, 16 2013
   phy%rel_QP  = _ONE_ - smooth_small(_ONE_- phy%rel_QP, self%small_finite)

   phy%QPN    = phy%P_reg / phy%N_reg
!   dom%QP     = dom%P  / (dom%C  + min_Cmass)  ! P:C ratio of DOM
!   det%QP     = det%P / (det%C + min_Cmass)  ! P:C ratio of detritus
end if   

! fraction of free (biochemically available) intracellular nitrogen
phy%rel_QN  = (phy%QN - self%QN_phy_0) * self%iK_QN
! added for deep detritus traps with extreme quotas kw Jul, 16 2013
phy%rel_QN  = _ONE_ - smooth_small(_ONE_- phy%rel_QN, self%small_finite)

if (self%PhotoacclimOn) then  
  ! calculate rel_chloropl
  call calc_rel_chloropl(self,phy,method=_KAI_)

  ! conversion of bulk chlorophyll concentration to chlorophyll content in chloroplasts  
  phy%theta     = phy%chl / (phy%rel_chloropl * phy%C_reg)   ! trait variable
! cell specific CHL:C ratio of chloroplasts / carbon bound to LHC per CHL-pigment
  phy%frac%theta= phy%chl/phy%C_reg * self%itheta_max ! []     no smaller than o(1.d-5)!
endif
! --- total pool-size of available/free proteins/enzymes and RNA -------------------   
phy%frac%TotFree= 1.0d0 
! -- remaining nitrogen fraction for uptake and nutrient processing --------------
phy%frac%NutUpt = smooth_small(phy%frac%TotFree - phy%frac%Rub - phy%frac%theta, self%small)

if (self%GrazingOn) then
  ! ---- herbivore stoichiometry ---------------------------
  zoo%QN    = self%const_NC_zoo
  zoo%N     = zoo%C * zoo%QN
  zoo%yield = self%yield_zoo
  zoo%flopp =  _ONE_ - self%yield_zoo
  if (self%PhosphorusOn) then 
    zoo%QP    = self%const_PC_zoo
    zoo%P     = zoo%C * zoo%QP
  endif
endif
end subroutine

end module maecs_functions
!------------------------------------------------------

