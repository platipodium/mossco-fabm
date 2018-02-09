!> @file maecs_functions.F90
!> @author Richard Hofmeister, Markus Schartau, Kai Wirtz, Onur Kerimoglu

#include "fabm_driver.h"
!---------------------------------------------------------
! !MODULE: MAECS_functions --- more to come
!> @brief  functions called by maecs_do, maecs_grazing and maecs_primprod
   module maecs_functions

! !USES:
   use fabm_types
   use maecs_types

   !private
   real(rk),private     ::zero   = 0.0_rk
   
   public   uptflex       ,& 
            queuefunc, queuefunc0, queuederiv ,&
            smooth_small, sinking, min_mass, &
            calc_sensitivities, calc_internal_states, nan_num
   
 contains  

 !------------------------------------------------------
!> @brief calculate the internal states 
!> @details 
!> @todo the theta-related calculations should obviously be related to \latexonly Eq. \ref{eq:ftheta} \endlatexonly but I get lost. See the Q's therein
subroutine calc_internal_states(maecs,phy,det,dom,zoo)

implicit none
class (type_maecs_base_model),intent(in)     :: maecs
type (type_maecs_phy), intent(inout) :: phy
type (type_maecs_om), intent(inout) :: det
type (type_maecs_om), intent(inout) :: dom
type (type_maecs_zoo), intent(inout) :: zoo

!> @fn maecs_functions::calc_internal_states()
!> 1. Calculate elemental absolute and relative quotas (Q and relQ):
!>   - phy\%Q\%X = phy\%X / phy\%C where x=N,P,Si
!>   - phy\%relQ\%X= (phy\%Q\%X - maecs\%qN_phy_0) / maecs\%iK_QN where x=N,P,Si
phy%Q%N    = phy%reg%N / phy%reg%C
! added for mixing effects in estuaries kw Jul, 15 2013
!write (*,'(A,2(E19.5))') 'QN-1:',phy%Q%N,maecs%small_finite* maecs%QN_phy_max

phy%Q%N  = smooth_small(phy%Q%N, maecs%QN_phy_0 + maecs%small_finite * maecs%QN_phy_max)

! fraction of free (biochemically available) intracellular nitrogen
phy%relQ%N  = (phy%Q%N - maecs%QN_phy_0) * maecs%iK_QN
phy%relQ%N  = smooth_small(phy%relQ%N,maecs%small)

! added for deep detritus traps with extreme quotas kw Jul, 16 2013
if(  phy%relQ%N .gt. 0.95d0*maecs%MaxRelQ ) then
   phy%relQ%N  = maecs%MaxRelQ - smooth_small(maecs%MaxRelQ- phy%relQ%N, maecs%small_finite)
endif

if (maecs%PhosphorusOn) then 
   phy%Q%P     = phy%reg%P / phy%reg%C
! added for mixing effects in estuaries kw Jul, 15 2013
   phy%Q%P     = smooth_small(phy%Q%P, maecs%QP_phy_0 + maecs%small_finite * maecs%QP_phy_max)

   phy%relQ%P = ( phy%Q%P - maecs%QP_phy_0 ) * maecs%iK_QP
   phy%relQ%P  = smooth_small(phy%relQ%P,maecs%small)

! added for deep detritus traps with extreme quotas kw Jul, 16 2013
   if(  phy%relQ%P .gt. 0.95d0*maecs%MaxRelQ ) then
     phy%relQ%P  = maecs%MaxRelQ - smooth_small(maecs%MaxRelQ- phy%relQ%P, maecs%small_finite)
   endif
else
   phy%Q%P     = maecs%QP_phy_0
   phy%relQ%P  = maecs%small_finite
end if 
  
if (maecs%SiliconOn) then 
   phy%Q%Si    = phy%Si / phy%reg%C
! added for mixing effects in estuaries kw Jul, 15 2013
   phy%Q%Si     = smooth_small(phy%Q%Si, maecs%QSi_phy_0 + maecs%small_finite *  maecs%QSi_phy_max)

   phy%relQ%Si = ( phy%Q%Si - maecs%QSi_phy_0 ) /(maecs%QSi_phy_max-maecs%QSi_phy_0) 

!TODO extend max-reQ to Si

! added for deep detritus traps with extreme quotas kw Jul, 16 2013
!   phy%relQ%Si = _ONE_ - smooth_small(_ONE_- phy%relQ%Si, maecs%small_finite)
!   phy%Q%SiN    = phy%Sii / phy%reg%N
end if   


!> @fn maecs_functions::calc_internal_states()
!> 2. Calculate Rubisco fraction (convert from the bulk variable) 
!>    - unpack phy\%frac\%Rub (=@f$ f_R @f$)= phy\%Rub / phy\%reg\%C
!>    - smooth 1-@f$ f_R @f$to (small\%finite + rel_chlropl_min) (both nml pars), such that @f$ f_R @f$ is always smaller than 1
!>    - @f$ f_R = \mathrm{phy\%Rub} / phy_C @f$

if (maecs%RubiscoOn) then 
! trait + transporter needs division to become a trait again
   phy%frac%Rub = phy%Rub / phy%reg%C
else
   phy%frac%Rub=maecs%frac_Rub_ini
end if   
 
!min-max correction: 
! note the  1-(1-x) structure
! which equals x for small x, and is slightly below x for x->1
! a counterpart of smooth-small, could be called smooth-at-one;
! ensures that f_R is always somewhat smaller than one, since it
! leaves a minimal fraction of resources to other compartments (f_V)
phy%frac%Rub = smooth_small(phy%frac%Rub , maecs%rel_chloropl_min)
phy%frac%Rub = _ONE_ - smooth_small(_ONE_- phy%frac%Rub ,maecs%small_finite + maecs%rel_chloropl_min)

!> @fn maecs_functions::calc_internal_states()
!> 3. Calculate @f$ \theta \mathrm{ and } f_{\theta} @f$ 
!\latexonly  See also: \ref{sec:uptsys} \endlatexonly
!>    - phy\%rel\_chloropl =  @f$  f_R*q_N^{\sigma} @f$
!>      - phy\%rel\_chloropl = factor that relates "chl-a per chloroplast" to "chl-a per cell-C"
!>        thus  chloroplast-C  over total intracellular C
!>    - @f$ \mathrm{phy\%theta} =  phy_{chl}/phy_C / \mathrm{phy\%rel\_chloropl} @f$
!>      - phy\%chl: bulk variable. biomass (phyc) times trait (theta) times factor)
!>      - thus, phy\%theta= (gchl/gC) * (gC/gchloroplastC) =gchl/ gchlroplastC
!>      - thus, with this transformation the transported trait variable (\%theta*phy_C) is translated into an observable, i.e. bulk CHL-a conc.
!>    - @f$ \mathrm{phy\%frac\%theta}= \mathrm{rel\_chloropl}*\theta* \mathrm{maecs\%itheta\_max} (=1/\theta_C) @f$
!>      - identically, \lref{eq. ,eq:ftheta,} says @f$ f_\theta = f_R*q_N^\sigma(=\mathrm{rel\_chloropl})*\theta / \theta_C @f$
!>    - @f$ \mathrm{phy\%frac\%NutUpt}=f_V = \mathrm{phy\%frac\%TotFree} - f_{\theta} - f_R  @f$, where phy\%frac\%TotFree=1.0

! calculate rel_chloropl
phy%rel_chloropl = smooth_small(phy%frac%Rub* phy%relQ%N**maecs%sigma,maecs%rel_chloropl_min)

if (maecs%PhotoacclimOn) then  

  ! conversion of bulk chlorophyll concentration to chlorophyll content in chloroplasts  
  phy%theta     = smooth_small(phy%chl / (phy%rel_chloropl * phy%reg%C),maecs%rel_chloropl_min*maecs%theta_LHC)   ! trait variable

! cell specific CHL:C ratio of chloroplasts / carbon bound to LHC per CHL-pigment
  phy%frac%theta= phy%theta * phy%rel_chloropl * maecs%itheta_max ! []     no smaller than o(1.d-5)!
!equally:
!  phy%frac%theta= (phy%chl / phy%reg%C) / theta_LHC
endif

!min-max correction: 
! note the  1-(1-x) structure
! which equals x for small x, and is slightly below x for x->1
! a counterpart of smooth-small, could be called smooth-at-one;
! ensures that f_R is always somewhat smaller than one, since it
! leaves a minimal fraction of resources to other compartments (f_V)
phy%frac%theta = _ONE_ - smooth_small(_ONE_- phy%frac%theta ,maecs%small_finite + maecs%rel_chloropl_min)


! --- total pool-size of available/free proteins/enzymes and RNA -------------------   
phy%frac%TotFree= 1.0d0 
! -- remaining nitrogen fraction for uptake and nutrient processing --------------
phy%frac%NutUpt = smooth_small(phy%frac%TotFree - phy%frac%Rub - phy%frac%theta, maecs%small)


!> @fn maecs_functions::calc_internal_states()
!> 4. Calculate zooplankton states:
!>    - @f$ zoo_{QX} = \mathrm{maecs\%const_NC_zoo} \mathrm{ , } zoo_{X} = zoo_C * zoo_{QN} \mathrm{ , } X=N,P @f$
!>    - @f$ zoo_{yield} = \mathrm{maecs\%yield_zoo} \mathrm{ , } zoo_{flopp} = 1-\mathrm{maecs\%yield_zoo} @f$
if (maecs%GrazingOn) then
  ! ---- herbivore stoichiometry ---------------------------
  zoo%Q%N   = maecs%const_NC_zoo
  zoo%N     = zoo%C * zoo%Q%N
  zoo%yield = maecs%yield_zoo
  zoo%flopp =  _ONE_ - maecs%yield_zoo
!  if (maecs%PhosphorusOn) then 
    zoo%Q%P   = maecs%const_PC_zoo
    zoo%P     = zoo%C * zoo%Q%P
! endif  
endif
end subroutine

!------------------------------------------------------
!> @brief calculate sensitivities
!> @details Details:
!> - sens\%f\_T \latexonly see eq. \ref{eq:arrhenius} \endlatexonly
!> - sens\%P\_max \latexonly see eq. \ref{eq:Pmax} \endlatexonly
!> - sens\%upt\_pot\%C \latexonly (=LH) see eq. \ref{eq:LH} \endlatexonly
!> - sens\%upt\_pot\%X, (@f$=V_X@f$), X=N,P,Si calculated by uptflex() \lref{see eq. ,eq:uptakecoeffcurr,.}
!\latexonly according to eq. \ref{eq:uptakecoeffcurr} \endlatexonly
!> @todo: a more intuitive name like calc_potentials?
!> @todo: Q: why maecs instead of self?
subroutine calc_sensitivities(maecs,sens,phy,env,nut,acc)

implicit none
class (type_maecs_base_model), intent(in) :: maecs
type (type_maecs_sensitivities), intent(out) :: sens
type (type_maecs_phy),intent(in) :: phy
type (type_maecs_env),intent(in) :: env
type (type_maecs_om),intent(in) :: nut
type (type_maecs_traitdyn), intent(out) :: acc

type (type_maecs_om) :: fA
real(rk) :: par, T_Kelv, NutF, affin, pmax, rqn
logical      ::IsAdap = .true.

if((maecs%adap_rub .lt. 0.001d0 .and. maecs%adap_theta.lt. 0.001d0) .or.  .not. maecs%PhotoacclimOn) IsAdap = .false. 
! switch off affinity-transport flexibility depending on overall flexibility; TODO: add new switch
if(IsAdap .and. maecs%tau_regV .gt. 100.0d0) IsAdap = .false. 
!> @fn maecs_functions::calc_sensitivities()
!> 1. calculate (sens\%) f\_T, P\_max\_T, a\_light, upt\_pot\%C 
par          = maecs%frac_PAR * env%par ! use active  fraction frac_PAR of available radiation
T_Kelv       = env%Temp + 273.d0 ! temperature in Kelvin 
! ----------------------------------------------------------------------------
! +++ determine rates for photoautotophic growth ++++++++++++++++++++++++++++++++++++++++++++++++++
! --- temperature dependence of metabolic rates (with T_ref given in units [Kelvin]) --------------
!standard Q10 RULE
sens%f_T     = maecs%rq10**((T_Kelv-maecs%T_ref)/10.0)
sens%f_T2    = (2*maecs%rq10)**((T_Kelv-maecs%T_ref)/10.0)

! --- (potential) maximum photosynthetic rate -----------------------------------------------------
sens%P_max_T = maecs%P_max * sens%f_T

! --- light response curve ------------------------------------------------------------------------
sens%a_light = maecs%alpha * par /(sens%P_max_T)  ! par NOCH BAUSTELLE, jetzt PAR in zellmitte 
sens%upt_pot%C  = 1.0d0 - exp(- sens%a_light * phy%theta) ! [dimensionless]
if (maecs%ChemostatOn) then
  if (maecs%rel_co2 .gt. 0.01d0) then 
! write (*,'(A,3(F10.4))') 'PAR CO2:',par,env%CO2 ,sens%upt_pot%C
    NutF    = max(env%CO2, 0.0d0)

! normalized affinity to DIC ([CO2]+[HCO3])
    pmax    = smooth_small(sens%upt_pot%C,maecs%small)
    affin   = pmax/maecs%rel_co2

! optimal partitioning between surface transporter and carboxylation/Rubisco
    fA%C    = fOptUpt(affin ,pmax, NutF, IsAdap)
    acc%fA%C= fA%C
! write (*,'(A,5(F10.4))') 'co2 A f Pm-> ',env%CO2,affin,fA%C, sens%upt_pot%C,uptflex(affin ,sens%upt_pot%C, NutF, fA%C)
    acc%fac3 = sens%upt_pot%C
    sens%upt_pot%C = uptflex(affin ,sens%upt_pot%C, NutF, fA%C)
  end if
end if

! --- carbon specific N-uptake: sites vs. processing ----------------------------------
!> @fn maecs_functions::calc_sensitivities()
!> 2. calculate fA\%X, (sens\%) upt\_pot\%X for each element, X
!>   - (acc\%)fA\%X= call: foptupt()
!>   - (sens\%)upt\_pot\%X= call: uptflex()
! non-zero nutrient concentration for regulation
!NutF          = max(nut%N,0.0d0)
NutF          = smooth_small(nut%N,maecs%small)
! optimal partitioning between
! surface uptake sites and internal enzymes (for assimilation)

! trait-hack: QmxaP ~ QN
if(maecs%mort_ODU .gt. 0.99) then
  rqn  = phy%relQ%N 
  if(rqn .gt. 1.0d0) rqn = 1.0d0
  if(rqn .lt. 0.1d0) rqn = 0.1d0
else
  rqn = 1.0d0
endif

fA%N          = fOptUpt(maecs%AffN/rqn,maecs%V_NC_max * sens%f_T2, NutF, IsAdap)
acc%fA%N      = fA%N
sens%upt_pot%N= uptflex(maecs%AffN/rqn,maecs%V_NC_max*sens%f_T2,NutF, fA%N)

!  P-uptake coefficients
if (maecs%PhosphorusOn) then 
   NutF          = smooth_small(nut%P,maecs%small)
! optimal partitioning 
   fA%P          = fOptUpt(maecs%AffP,maecs%V_PC_max * sens%f_T2, NutF, IsAdap)
   acc%fA%P      = fA%P
   sens%upt_pot%P= uptflex(maecs%AffP,maecs%V_PC_max*sens%f_T2,Nut%P,fA%P)
else
   acc%fA%P      = 0.5d0
   acc%Av%P      = 0.d0
end if

!  Si-uptake coefficients
if (maecs%SiliconOn) then 
   NutF          = smooth_small(nut%Si,maecs%small)
! optimal partitioning 
   fA%Si         = fOptUpt(maecs%AffSi,maecs%V_SiC_max * sens%f_T2, NutF, IsAdap)
   acc%fA%Si     = fA%Si
   sens%upt_pot%Si= uptflex(maecs%AffSi,maecs%V_SiC_max * sens%f_T2,nutF,fA%Si)
else
   acc%fA%Si     = 0.5d0
   acc%Av%Si     = 0.d0
end if
! TODO check temperature dependence of nutrient affinity


end subroutine

!-----------------------------------------------------------------------
!> @brief opt. partitioning between surf upt sites and intern. enzymes for nut assim.
!> @details 
!> calculates @f$ f_{A,X} @f$ \lref{see eq. ,eq:optutpalloc,.}
!> @f$ f_{A,X} = 1/ (1+ \sqrt(A_X^0*DIX/V_{max,X}^0)) @f$
pure real(rk) function fOptUpt(Aff0, Vmax0, Nut, IsOn)
   implicit none

   real(rk), intent(in)      :: Aff0, Vmax0, Nut
   logical, intent(in)       :: IsOn
   if (IsOn) then
! avoid fOptUpt=1 since Vmax=0 will induce a NaN in uptflex at Nut=0
     fOptUpt     = 1.d0/(sqrt(Aff0*Nut/(Vmax0)) + 1.001d0)
   else
     fOptUpt     = 0.5d0
   end if
   end function fOptUpt
   
!-----------------------------------------------------------------------
!> @brief calc's pot nut upt as f(external conc, allocations)
!> @details 
!> uptake regulation: sites vs. processing\n
!> \latexonly see eq. \ref{eq:uptake} - \ref{eq:optutpalloc} \endlatexonly
!> - sens%upt\_pot%X = @f$ (A_X*DIX*V_{max,X})/(A_X*DIX + V_{max,X}))^{-1} @f$
!> - @f$ A_X = f_{A,X}*A_X^0 @f$, where @f$ f_{A,X} @f$ =fOptUpt()
!> - @f$ V_{max}=(1-f_{A,X})*V_{max}^0*f_T @f$

pure real(rk) function uptflex(Aff0, Vmax0, Nut, fAv)
   implicit none
   real(rk), intent(in)      :: Aff0, Vmax0, Nut, fAv
   real(rk)                  :: Aff,Vmax

   Aff     = fAv * Aff0    ! nutrient affinity 
   Vmax    = (_ONE_-fAv) * Vmax0   ! maximal N-uptake rate 

   uptflex = Vmax*Aff*Nut/(Aff*Nut + Vmax + 1E-4)

   end function uptflex


!-----------------------------------------------------------------------
!> @brief the queue function 
!> @details 
!> provides both the queuing function and it's derivative 
!> with the parameter n->inf :liebig and n~1:product
!> \latexonly see: Section \ref{sec:colim} \endlatexonly \n
!> @todo: add equations
subroutine queuefunc(n,x,qfunc,dq_dx,dq_dn)

   implicit none
   real(rk), intent(in)          :: x, n
   real(rk), intent(out)         :: qfunc, dq_dx, dq_dn
   real(rk)                      :: px, dn

   if(abs(_ONE_-x) .lt. 1E-2) then
      qfunc = n/(n+_ONE_)
      dq_dx = qfunc/2 ! 1./(2*(1+hh)); 
      dq_dn = _ONE_/(n+_ONE_)**2
   else
      px    = x**(n+_ONE_)
      dn    = _ONE_ / (_ONE_-px)
      qfunc =  (x-px) * dn
      dq_dx = (_ONE_ -(n+_ONE_)*x**n+n*px)*dn*dn
      dq_dn = px*(x-_ONE_)*dn*dn * log( x + 1E-4)
   endif
   end subroutine queuefunc

!-----------------------------------------------
!> @brief numerical approximation of the queue function 
!> @details 
!> n->inf :liebig  n~1:product\n
!> Here is an example of adding a snip of code:
!> @snippet maecs_functions.F90 queuefunc1_snippet 
subroutine queuefunc1(n,x,qfunc,dq_dx)

   implicit none
   real(rk), intent(in)       :: x, n
   real(rk), intent(out)      :: qfunc, dq_dx
   real(rk)                   :: nn, hh,x0,en

   ! [queuefunc1_snippet]
   nn = n+1.
   hh = 1./nn
   x0 = (log(exp(nn)-1))*hh
   en = exp(-nn*(x/(1+hh*x) - x0))
   qfunc =  1. - hh*log(1.+ en)
   dq_dx = 1. / ( (1. + hh * x)**2 * (1.+1./en))
   ! [queuefunc1_snippet]
   
   end subroutine queuefunc1

!-------------------------------------------------------------
!> @brief derivative of the queue function ??
!> @details 
!> @return queuederiv
!> @todo: add equations
real(rk) function queuederiv(n,x)
!
   implicit none
   real(rk), intent(in)          :: x, n
   real(rk)                      :: nqueue_aa
   real(rk)                      :: nqueue_bb = 1.3863d0
!-----------------------------------------------------------------------
   nqueue_aa = 6.9315d-1*(1.d0+1.0d0/(n))
   queuederiv = exp(-nqueue_aa*(x)**nqueue_bb) 
   end function queuederiv

!-------------------------------------------------------------
!> @brief calculation of the sinking rate
!> @details 
!> - sinkvel   = @f$ -vS * e^{-4*\mathrm{phys\_status}} @f$
!>   - phys\_status calculated in fabm_hzg_maecs::maecs_get_vertical_movement()
!>   - Future work: sinkvel=f(size) \latexonly (see section \ref{sec:sink}) \endlatexonly \n
subroutine sinking(vS ,phys_status,sinkvel)
   implicit none
   real(rk), intent(in)     :: vS ,phys_status
   real(rk), intent(out)    :: sinkvel

   !   real(rk)         :: 
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
!> @brief minimum mass
!> @details  pushes the phyC,N and P to some lower boundary according to 4 different methods (controlled by the mm_method parameter):
!> phy\%N and phy\%C are stored in phy\%reg\%N and phy\%reg\%C, respectively
!> 1. if phy\%N <= 1e-7; phy\%N=1e-7, phy\%C=phy\%N/QN(aver), phy\%P=phy\%C*QP(aver)
!> 2. ..
!> 3. ..
!> 4. ..
!> @todo: assign some meaningful names to case numbers?
!> @todo: mm_method to be read from the nml?
!> @todo: Q: phy\%reg\%P either non existent or commented out for different cases. Why?
!> @todo: add equations
subroutine min_mass(maecs,phy,min_Cmass,iscritical,method)

implicit none

class (type_maecs_base_model), intent(in)      :: maecs
type (type_maecs_phy), intent(inout)   :: phy
real(rk), intent(out)                  :: min_Cmass
logical,intent(out)                    :: iscritical
integer, intent(in), optional          :: method

real(rk)     :: min_Nmass, delta_C, delta_N ! min_Cmass, 
integer      ::  mm_method=_KAI_
logical      ::  ischanged 

if (present(method)) mm_method=method
ischanged = .false.
iscritical= .false.

select case (mm_method)
 case (_MARKUS_)
   if (phy%N .le. 1.d-7) then
      phy%N = 1.d-7 
      phy%C = phy%N / maecs%aver_QN_phy
      if (maecs%PhosphorusOn) then
         phy%P = phy%C * maecs%aver_QP_phy
      end if
   end if
   phy%reg%N = phy%N
   phy%reg%C = phy%C
!   phy%reg%P = phy%P 

 case (_KAI_)
! -------------------------------------------------------------------------------
! --- Here, quota are smoothed as soon as phytoplankton biomass/biovolume 
!     approaches a lower threshold 'min_mass', in the order of o(1.d-4) / 
! set small boundary depending on numerical resolution
! TODO: insert h ~ level height, here 10cm 
   min_Cmass = maecs%small_finite * 1 !* 1.0d-2 / maecs%a_spm  
   min_Nmass = min_Cmass * maecs%aver_QN_phy
   phy%reg%C = smooth_small( phy%C , min_Cmass)
   if (abs(phy%C-phy%reg%C) .gt. 1d-2*min_Cmass) then
      phy%reg%N = phy%reg%C * maecs%aver_QN_phy
      ischanged = .true. 
   else 
      phy%reg%N = smooth_small( phy%N , min_Nmass)
      if (abs(phy%N-phy%reg%N) .gt. 1d-2*min_Nmass) then
         phy%reg%C = phy%reg%N / maecs%aver_QN_phy
         ischanged = .true.
      end if  
   end if
!write (*,'(A,4(F10.3))') 'P=',phy%P,phy%reg%C * maecs%aver_QP_phy,smooth_small(phy%P,min_Cmass * maecs%aver_QP_phy)
   if (ischanged) then ! retune P and Rub
!      if (maecs%PhosphorusOn)  phy%reg%P =  phy%reg%C * maecs%aver_QP_phy
      if (maecs%RubiscoOn) then 
         phy%Rub =  phy%reg%C * maecs%frac_Rub_ini
      end if
      if (maecs%PhotoacclimOn) then 
         phy%chl =  phy%reg%C * maecs%frac_chl_ini
      end if

   else  ! additional check for Rub and P; TODO: omitt ??
       phy%Rub = smooth_small( phy%Rub , min_Cmass * maecs%frac_Rub_ini)
!      if (maecs%PhosphorusOn)  phy%reg%P =  smooth_small(phy%P,min_Cmass * maecs%aver_QP_phy)
   end if ! ischanged
   if ((abs(phy%C-phy%reg%C) .gt. 1d-1*min_Cmass) .or. (abs(phy%N-phy%reg%N) .gt. 1d-1*min_Cmass* maecs%aver_QN_phy) ) iscritical= .true.

 case (2)
! -------------------------------------------------------------------------------
! --- Here, phyC and phyN are smoothed as soon as biomass approaches 'min_mass',  
! set small boundary depending on numerical resolution
! TODO: insert h ~ level height, here 10cm 
   min_Cmass = maecs%small_finite !* 1.0d-3 / maecs%a_spm  
   min_Nmass = min_Cmass * maecs%aver_QN_phy
   phy%reg%C = smooth_small( phy%C , min_Cmass)
   delta_C   = phy%reg%C - phy%C
   delta_N = 0.0_rk
  
   if (abs(delta_C) .gt. 1d-2*min_Cmass) then !.and. phy%N .lt. 1*min_Nmass
      phy%reg%N = smooth_small(phy%N + delta_C * maecs%aver_QN_phy, min_Nmass)
      if (maecs%PhosphorusOn) then
         phy%reg%P = phy%P + delta_C *  maecs%aver_QP_phy
      end if
      ischanged = .true. 
   else 
      phy%reg%N = smooth_small( phy%N , min_Nmass)
      if (maecs%PhosphorusOn) phy%reg%P = smooth_small( phy%P , min_Cmass * maecs%aver_QP_phy)
      delta_N   = phy%reg%N - phy%N
      if (abs(delta_N) .gt. 1d-2*min_Nmass) then
         phy%reg%C = smooth_small(phy%C + delta_N / maecs%aver_QN_phy, min_Cmass)
!         delta_C   = phy%reg%C - phy%C
         ischanged = .true.
      end if  
   end if
   if (ischanged) then
     if (abs(phy%C-phy%reg%C) .gt. 5d-1*min_Cmass .or. abs(phy%N-phy%reg%N) .gt. 5d-1*min_Cmass* maecs%aver_QN_phy ) iscritical= .true.
     phy%relax = max(0.0d0,delta_C)/phy%reg%C
   endif
!!      if (phy%reg%N .gt. 0.2 * phy%reg%C) phy%reg%N = 0.2 * phy%reg%C

 case (3)
! -------------------------------------------------------------------------------
! --- Here, quota are smoothed as soon as phytoplankton biomass/biovolume 
!     approaches a lower threshold 'min_mass', in the order of o(1.d-4) / 
! set small boundary depending on numerical resolution
! TODO: insert h ~ level height, here 10cm 
   min_Cmass = maecs%small_finite !* 1.0d-3 / maecs%a_spm  
   phy%reg%C = smooth_small( phy%C , min_Cmass)
   phy%reg%N = smooth_small( phy%N , min_Cmass * maecs%aver_QN_phy)
!   if (maecs%PhosphorusOn)  phy%reg%P = smooth_small( phy%P , min_Cmass * maecs%aver_QP_phy)

   if (maecs%RubiscoOn) then 
     phy%Rub = smooth_small( phy%Rub , min_Cmass * maecs%frac_Rub_ini)
   end if    
end select
end subroutine min_mass
   

   !---------------------------------------------------------
!> @brief  continous smoothing function by kw Apr 2012
!> @details 
!! smoothly converges x to **eps/2** for x<eps  
!! \f[ x=eps+(x-eps)*e^{x/eps}/(1+e^{x/eps}) \f]
pure real(rk) function smooth_small(x, eps)

   implicit none
   real(rk), intent(in)          :: x, eps
   real(rk)                      :: arg, larger, larger2
!   integer                       :: na = 2
   integer                       :: nb

!--------------------------------------------------------------
   nb      = 8
   if (x .lt. nb*eps) then
     arg     = x/(eps+1E-7)
     larger  = exp(2*arg)
     larger2 = exp(arg/2)   
     smooth_small  = ((nb+larger2)*eps + x*larger)/(nb+larger) 
   else
    smooth_small  = x
   endif
   end function smooth_small 

   !---------------------------------------------------------
!> @brief  converts nan-to be numbers to a real num (eg -9)
!> @details 
! \f[ x=eps+(x-eps)*e^{x/eps}/(1+e^{x/eps}) \f]
real(rk) function nan_num(x)

   implicit none
   real(rk), intent(in)          :: x
   !real(rk)            :: xnan,xinf,xinfneg
!--------------------------------------------------------------
   !xnan=zero/zero
   !xinf=1.0_rk/zero
   !xinfneg=-1.0_rk/zero
   !if ((x .eq. xinfneg) .or. (x .eq. xinf) .or. (x .eq. xnan)) then
   ! The above calculations trigger invalid in many compilers, thus the
   ! better formulation below
   if (      x /= x           & ! this is true for NaN 
      .or.   abs(x) > huge(x) & ! this is true for Inf 
      ) then  
     nan_num=-2.e20_rk !that's the default FABM missing value
   else 
     nan_num=x
   endif
   !write(0,*) xnan,x,nan_num
   
   end function nan_num 
   
end module maecs_functions
!------------------------------------------------------

