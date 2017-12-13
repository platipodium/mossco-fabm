!> @file maecs_primprod.F90
!> @author Kai Wirtz, Richard Hofmeister, Markus Schartau, Onur Kerimoglu

#include "fabm_driver.h"

!> @brief Primary production module
   module maecs_primprod

   use fabm_types
   use maecs_types
   use maecs_functions
   private
   public    photosynthesis     
 contains  

!> @brief  calculates photosynthesis
!> @details 
!> This is the subroutine, where the optimal regulation of phytoplankton traits
!> are described, which is central to the physiological-MAECS
subroutine photosynthesis(self,sens,phy,nut,uptake,exud,acc)
implicit none

class (type_maecs_base_model), intent(in)   :: self
type (type_maecs_sensitivities),intent(in), target :: sens
type (type_maecs_phy), intent(inout)       :: phy
type (type_maecs_om),intent(in) :: nut
type (type_maecs_om), intent(out), target  :: uptake
type (type_maecs_om), intent(out)          :: exud
type (type_maecs_traitdyn), intent(out), target    :: acc

integer  :: num_nut, i, j
real(rk) :: eps
! --- UPTAKE KINETICS
type (type_maecs_om), target :: upt_act ! free C-specific uptake rate [molX molC^{-1} d^{-1}]
! --- DERIVATIVES OF C-UPTAKE (FOR TRAIT-BASED ADAPTIVE DYNAMICS)
real(rk) :: dmu_dV, dmu_dV0   ! [...]
real(rk) :: dmuQ_dfracR,   dmuQ_dtheta
real(rk) :: dmu_dfracR, dmu_dtheta  ! [...]
real(rk) :: grad_theta , grad_fracR ! [...]
real(rk) :: flex_theta , flex_fracR ! [...]
real(rk) :: dfV_dfracR , dfV_dtheta
real(rk) :: grossC, Pmaxc
real(rk) :: dbal_dv, dmu_daV_tot
real(rk) :: syn_act, hh
real(rk) :: f_Lip, q_Lip, q_NoLip 
real(rk) :: c_h1, c_hq           ! correction terms of queuing function 
real(rk) :: fac_colim            ! colimitation factor (dimensionless)
real(rk) :: q_X, qp_Y, qp_X      ! placeholder for nutrient availability
real(rk) :: qfunc, d_qf_dx, d_qf_dn   ! queuing function and derivative
real(rk) :: prod_dq, prod_dn
real(rk) :: d_X, d_QX, dmu_daV, act_V, resp
real(rk) :: resC, darkf, ex, fmaxf, maxrq,safe
real(rk), dimension(8) :: sigmv ! relative quota-uptake feed-back 
real(rk), dimension(8) :: dqp_X_dq_X, dqp_X_dqp_Y, dqp_X_dn ! derivatives
real(rk), dimension(8) :: zeta_X  ! C-costs of assimilation of nutrient X
!real(rk), dimension(5) :: upt_act
type (stoich_pointer), dimension(8)::elem ! struct-pointer addressing elements wthin loops
!
!eps     =  self%small_finite ! just  a shorter namer for a small thing
eps     = 1E-2 ! just  a shorter namer for a small thing
maxrq   = 2.0d0 ! self%MaxRelQ
! TODO: energetic costs of P-assimilation \partial (\zeta_CN V_N) / \partial V_P not
!    resolved but assumed to be already included in protein synthesis
! prelim solution: stoichiometry in RNA (N:P ~ 4:1) and phospholipids (N:P~1:1)
! TODO: include proteins/mebranes (N:P >> 16:1) under low growth conditions 
!  q_NoLip  = 3.8  ! P-stoichiometry of active compounds (DNA, RNA)  
! q_NoLip  = self%zstoich_PN !   16 : Redfield-stoichiometry 
q_Lip    = 1.  ! 0.8 storage P-stoichiometry   
!f_Lip    = 1./(1.+exp(0*(1.-phy%relQ%P)))  f_Lip    = 0.5d0

!> @fn maecs_primprod::photosynthesis()
!> 1. Prepare loop over structure elements by assigning a pointer structure
!>    - Here, every possible nutrient is asked explicitely; thus first Si then P\n
!> @include 'maecs_stoichvars.F90p'
include './maecs_stoichvars.F90p' ! order of recursive colimitation scheme
num_nut  = self%nutind%nutnum  ! total number of nutrients

! --- relative amount of carbon invested into light harvesting complex (LHC) -------
! chlorophyll-to-carbon ratio of chloroplast * chloroplast concentration of cell
!> @fn maecs_primprod::photosynthesis()
!> 2. calculate the dchl/dtheta, dchl/dfracR, dchl/dQN, dbal/dv
!if (self%PhotoacclimOn) then  
! TODO complete derivatives for cases sigma .ne. 0 or 1
acc%dRchl_dtheta = phy%rel_chloropl 
acc%dRchl_dfracR = phy%theta * phy%relQ%N**self%sigma !* (1.0d0-self%sigma)
acc%dRchl_dQN    = phy%theta * phy%frac%Rub * self%sigma  * self%iK_QN !* phy%relQ%N**(self%sigma-1.0d0)
!endif

! --- gross carbon uptake by phytoplankton (gross-primary production) ---------------
!   phy%frac%Rub* sens%P_max_T* phy%relQ%N : carboxylation capacity = Rub * free proteins
!   fac_colim             : synchrony in  protein/RNA dynamics 
!   sens%upt_pot%C           : light harvesting (light limited growth)   

! --- synchrony in nutrient assimilation depends on growth cycle and N-quota
!> @fn maecs_primprod::photosynthesis()
!> 3. if \mathrm{self\%syn\_nut<0}: synergy is proportional to f_V 
!> @todo: explain & details         
if (self%syn_nut .lt. -0.001 ) then  !self%nutind%iP
!  synchrony increases with uptake machinery
!   syn_act = -self%syn_nut * phy%frac%NutUpt
!  if (self%MaxRelQ .lt. 1E3 ) then
!    syn_act = smooth_small(-self%syn_nut * phy%frac%NutUpt ,eps)
!  else
!   syn_act = smooth_small(-self%syn_nut * phy%frac%NutUpt* phy%frac%theta* phy%frac%Rub ,eps)
    syn_act = smooth_small(-self%syn_nut * (phy%relQ%N+0.5d0),eps)
!  endif
else
   syn_act  = self%syn_nut
endif 
acc%tmp   =   syn_act ! store

! $\climf = q_1\cdot g_h (q_2'/q_1) \cdot (1+h\cdot q_1 q_2' +c_h)$
! metabolic interdependence h is inverse of synergy syn_act
hh      = 1.0d0 / (syn_act + eps)

! one plus correction coefficient $c_h$ to ensure convergence to Liebig and product rule
 c_h1 = 1.0d0  + log( 1.0d0/(4**hh) + 0.5*hh );

! partial derivative of the balance eq. to uptake V
! equals ratio of gross to net primary production
dbal_dv = 1.0d0 + phy%Q%N * zeta_X(self%nutind%iN)  ! TODO P-uptake?
! e_N0    = 1.0d0 / (phy%Q%N * dbal_dv) 

!> @fn maecs_primprod::photosynthesis()
!> 4. retro-loop over structure elements with arbitrary number of nutrients
!> to calculate a process-based co-limitation factor and derivative terms needed for 
!> extended optimality functions
!> @todo: details
qp_Y    = elem(num_nut)%relQ

dqp_X_dq_X(num_nut)  = 1.0d0
dqp_X_dqp_Y(num_nut) = 0.0d0
dqp_X_dn(num_nut)    = 0.0d0

!if (num_nut .gt. 1) then
 ! loop over all nutrients starting from pre-final
do i = num_nut-1, 1, -1  

   q_X   = elem(i)%relQ
! correction of queuing function for symmetric/interdependent co-limitation 
   c_hq  = c_h1 + hh * qp_Y * q_X

! qfunc: relative nutrient limitation due to shortage in other nutrients 
   call queuefunc(syn_act, qp_Y / q_X, qfunc, d_qf_dx, d_qf_dn)

! $d_\mathrm{N} = \fracd{1}{\fQ{N}} - \pdiff{g_h}{x}\fracd{q_2'}{g_h\fQ{N}^2} +\fracd{h q_2'}{c_{hq}}$
   qp_X  = q_X * qfunc * c_hq

! +++ derivative of C-uptake rate with respect to quota ++++++++++++++++++++++++++++++
!   qfunc = smooth_small(qfunc,eps)
!   q_X   = smooth_small(q_X,eps)

   dqp_X_dq_X(i)  = qp_X* ( 1.0d0/q_X + qp_Y *( hh/c_hq - d_qf_dx/(qfunc*q_X*q_X+1E-4) ))

   dqp_X_dqp_Y(i) = qp_X *( hh*q_X/c_hq + d_qf_dx/(qfunc*q_X+1E-4) )

! +++ derivative of C-uptake rate with respect to exponent n ++++++++++++++++++++++++++++++
   dqp_X_dn(i) = q_X * d_qf_dn * c_hq  ! TODO make diff calculus complete (c_hq)

   qp_Y = qp_X

 end do
!else
!   qp_X = qp_Y   ! initial value for num_nut=1
!end if

!   final efficiency in interdependent multi-nutrient processing
fac_colim   = qp_X

!acc%fac1 = fac_colim
!acc%fac2 = elem(1)%relQ
!acc%fac3 = elem(2)%relQ

! phy%rel_phys= fac_colim * sens%upt_pot%C  ! auxiliary variable for sinking routine; not used?

! recursive product following from chain rule; first term : 1/LF
prod_dq     = 1.0d0/smooth_small(qp_X,eps)

! derivative of f_V (fraction of nutrient uptake machinery) on f_R (PS) and theta (LHC)

dfV_dfracR  = - 1*acc%dRchl_dfracR * self%itheta_max - 1.0d0 
dfV_dtheta  = - (acc%dRchl_dtheta + 0.*phy%frac%Rub+ 0*phy%frac%theta)  * self%itheta_max 
!dfV_dtheta  = - acc%dRchl_dtheta * self%itheta_max - phy%frac%Rub/(phy%theta + 1E-3)

! phy%frac%Rub  

! phy%frac%Rub* phy%relQ%N**maecs%sigma
!  auxiliary variable for the feed-back of N-quota and N-uptake change
! \todo stability   sigmv(1)    =   acc%dRchl_dQN  * self%itheta_max / phy%frac%NutUpt
sigmv(1:num_nut)      = 0.0d0
sigmv(self%nutind%iN) = acc%dRchl_dQN * self%itheta_max /(phy%frac%NutUpt+eps)

! --- gross carbon uptake by phytoplankton (gross-primary production) ---------------
! I DONT UNDERSTAND THE FOLLOWING:
! phy%frac%Rub * sens%P_max_T* phy%rel_QN : carboxylation capacity = Rub * free proteins

!> @fn maecs_primprod::photosynthesis()
!> 5. Calculate carbon uptake by phytoplankton:
!> - Pmaxc=fac\_colim-eps * sens\%Pmax\_T
!>   + fac\_colim : final synchrony in  protein/RNA dynamics (multi-nutrient processing) (obtained by the que function)
!> - grossC=Pmaxc*fR*sens\%upt\_pot\%C
!>   + sens\%upt\_pot\%C  : light harvesting (light limited growth)
!> - darkf= 1-exp(-grossC/self\%res0)

Pmaxc     = fac_colim * sens%P_max_T
grossC    = phy%frac%Rub * Pmaxc * sens%upt_pot%C  ! primary production

!acc%fac1 = phy%frac%Rub
!acc%fac2 = Pmaxc
!acc%fac3 = sens%upt_pot%C

! "darkness correction": marginal use should converge towards zero at very low light
!  offset (background respiration) derived from a_V(dmu_daV=0)*1/4*Vmax*zeta (f_A=f_V=1/2)

!darkf    = smooth_small(1.0d0 - exp(-2*fac_colim*sens%upt_pot%C),eps)
darkf = 1.0d0
!> @fn maecs_primprod::photosynthesis()
!> 6. Calculate (a first approximation of the ?!) activity @f$ a_{V,X} @f$, for each nutrient, X
!> - !! @f$ a_{V,X} @f$ are instantaneously optimized traits !!
!> @todo: details
! $\partial mu/\partial Q dQ/dV|_{tot}$ : initial zero berfore loop
dmuQ_dfracR = 0.0d0
dmuQ_dtheta = 0.0d0
dmu_daV_tot = 0.0d0
prod_dn     = 0.0d0

!acc%fac1 = prod_dq * dqp_X_dq_X(self%nutind%iN)
!acc%f = darkf

!do i = num_nut-1,1,-1 
do i = 1,num_nut-1 ! skip i=N:carbon
 ! ------------------ derivatives of co-limited growth function ---------------------
   d_X       = prod_dq * dqp_X_dq_X(i)

   prod_dn   = prod_dn + prod_dq * dqp_X_dn(i) ! integral diff wrt exponent for each quota
   prod_dq   = prod_dq * dqp_X_dqp_Y(i)
!   e_N       = e_N0 + d_X * elem(i)%iKQ / elem(i)%relQ
! contribution from synchrony
   if (i .eq. self%nutind%iN .and.  self%syn_nut .lt. -0.001.and. .true.) then
!     d_X     = d_X + prod_dn * (-self%syn_nut) * phy%frac%NutUpt 
     d_X     = d_X + prod_dn * (-self%syn_nut)  
   endif

   d_QX      = d_X* dbal_dv * elem(i)%iKQ + sigmv(i)* phy%Q%N * zeta_X(self%nutind%iN)  ! 
!TODO: include P-uptake resp  = zeta_X(self%nutind%iN) * (upt_act%N + self%zstoich_PN * upt_act%P)  
!if (phy%Q%P .gt. 0.01 .and. i .eq. self%nutind%iP) write (0,'(A,4(F9.4))') 'dQP=',d_X,dbal_dv,d_X* dbal_dv,sigmv(i)* phy%Q%N * zeta_X(self%nutind%iN)

! marginal use should converge towards zero at bad productivity conditions
!   (refine assumption Q\mu=V in derivation of d_QX)
   d_QX      = d_QX * darkf  
!if(i .eq. -self%nutind%iN) then
!  acc%fac1 = dqp_X_dq_X(i)
!  acc%fac2 = sens%upt_pot%C
!endif

!if(i .eq. self%nutind%iN) acc%fac2 = d_QX

!   steady-state down-regulation of uptake I: balance of respiration and indirect benefits 
   dmu_dV    = (1.0d0 + zeta_X(i) * elem(i)%Q) * d_QX/(1.0d0 + elem(i)%Q * (d_QX + sigmv(i)))

!   dmu_dV    = dmu_dV * e_N / (e_N + sigmv(i))

!   steady-state down-regulation of uptake I: balance of respiration and indirect benefits  
   safe      = 1.0d0 + max(0.0d0,(elem(i)%relQ-1.0d0)/maxrq)**3 * exp(-1*sens%P_max_T/self%res0)
   dmu_daV   = (-zeta_X(i)*safe + dmu_dV/safe) * phy%frac%NutUpt * elem(i)%upt_pot
!if (phy%Q%P .gt. 0.015 .and. i .eq. self%nutind%iP) write (0,'(A,5(F9.4))') 'aP=',dmu_daV,-zeta_X(i),dmu_dV,d_QX,d_QX/(1.0d0 + elem(i)%Q * (d_QX + sigmv(i)))
!if(i .eq. self%nutind%iN) then
!  acc%fac1 = -zeta_X(i)* phy%frac%NutUpt * elem(i)%upt_pot
!  acc%fac2 = dmu_dV* phy%frac%NutUpt * elem(i)%upt_pot
!endif
  
!   smoothed version of step function, uses marginal gain to emulate a continuous response 
!   act_V     = 1.0d0/(1.0d0 + exp( 3.1415d0 - self%tau_regV * dmu_daV));  ! 0.02
   act_V     = 1.0d0/(1.0d0 + exp( 1.0d0 - self%tau_regV * dmu_daV));  ! 0.02

   dmu_daV   = abs(dmu_daV)
   dmu_daV_tot = dmu_daV_tot + act_V * dmu_daV

   elem(i)%dmudV   = dmu_dV !* darkf!sqrt(darkf)
   elem(i)%dmudaV  = dmu_daV
   elem(i)%aV      = act_V 
!   write (*,'(A,4(F10.4))') 'aV:',dmu_dV,d_QX,elem(i)%aV,elem(i)%dmudaV
   
end do

!acc%fac1 = elem(1)%dmudV
!acc%fac2 = elem(self%nutind%iP)%relQ

!> @fn maecs_primprod::photosynthesis()
!> 7. (update !? - why not in the above loop?) the @f$ a_{V,X} @f$)
!> calculate nutrient uptake rates
!> calculate the derivative terms dmuQ/dfR, dmuQ/dtheta
dmuQ_dfracR = 0.0d0
dmuQ_dtheta = 0.0d0

do i = 1, num_nut-1 ! skip i=N:carbon
  if( self%adap_rub .gt. 0.001 .or. self%adap_theta .gt. 0.001) then
     act_V           = elem(i)%aV**2 * elem(i)%dmudaV/ (dmu_daV_tot + eps)
!    act_V           = elem(i)%aV
  else
     act_V           = max(0.0d0, 1.0d0 - elem(i)%relQ )
  endif
!  if (act_V .lt. 0.1 .and. i .eq. self%nutind%iN) act_V=0.1d0
! emulates passive Si diffusion through membrane (\todo not to be assimilated)
   if (self%SiliconOn) then
      if (i .eq. self%nutind%iSi .and. act_V .lt. 0.333d0 .and. elem(i)%upt_pot .gt. eps) act_V = 0.333d0  !num_nut  
   endif
   elem(i)%aV      = act_V 
   elem(i)%upt_act = act_V * elem(i)%upt_pot
   elem(i)%upt     = phy%frac%NutUpt * elem(i)%upt_act  ! [(molX) (molC)^{-1}
   ! d{-1}]

! emulates passive Fe diffusion through membrane (\todo not to be assimilated)
!TODO: check for similar process in Iron uptake
!   if (self%IronOn) then
!      if (i .eq. self%nutind%iFe .and. act_V .lt. 0.333d0 .and. elem(i)%upt_pot .gt. eps) act_V = 0.333d0  !num_nut  
!   endif
  
   elem(i)%aV      = act_V 
   elem(i)%upt_act = act_V * elem(i)%upt_pot
   elem(i)%upt     = phy%frac%NutUpt * elem(i)%upt_act  ! [(molX) (molC)^{-1} d{-1}]
!   write (*,'(A,2(F10.4))') 'upt:',elem(i)%aV,elem(i)%upt


! +++ derivative of C-uptake rate with respect to quota ++++++++++++++++++++++++++++++
   dmuQ_dfracR     = dmuQ_dfracR + elem(i)%dmudV * (elem(i)%upt_act+0*eps) * dfV_dfracR
   dmuQ_dtheta     = dmuQ_dtheta + elem(i)%dmudV * (elem(i)%upt_act+0*eps) * dfV_dtheta
! small *eps* correction at vanishing productivity since now aV=0 would entirely decouple regulation
! if (dmuQ_dfracR .lt. -20. .or. abs(dmuQ_dfracR+1.d0) .lt. 0.01) write (*,'(A,I3,10(F10.4))') 'Q',i,dmuQ_dfracR , elem(i)%dmudV , elem(i)%relQ,(elem(i)%upt_act+eps) , dfV_dfracR,act_V , elem(i)%upt_pot,Nut%P,Nut%N,grossC
end do
!acc%fac1 = elem(self%nutind%iP)%dmudV * elem(self%nutind%iP)%upt_pot* dfV_dfracR
!acc%fac1 = elem(self%nutind%iN)%dmudV * elem(self%nutind%iN)%upt_pot* dfV_dfracR

! account for differential effect of f_V on synchrony
if ( self%syn_nut .lt. -0.001 .and. .false.) then
!   acc%fac2 = prod_dn
   prod_dn = prod_dn * (-self%syn_nut)  * grossC *phy%relQ%N
   dmuQ_dfracR = dmuQ_dfracR + prod_dn * dfV_dfracR
   dmuQ_dtheta = dmuQ_dtheta + prod_dn * dfV_dtheta
!   acc%fac3 = prod_dn 
endif

!if (self%SiliconOn) uptake%Si=uptake%Si + 0.5*sens%upt_pot%Si
!if (self%IronOn) uptake%Fe=uptake%Fe + 0.5*sens%upt_pot%Fe

!> @fn maecs_primprod::photosynthesis()
!> 8. calculate phy%resp=f(uptake\%N), phy\%gpp=grossC-phy%resp, uptake\%C=grossC-phy%resp
! ---  respiration due to N & P assimilation --------------------------------------
phy%resp     = zeta_X(self%nutind%iN) * (uptake%N + self%zstoich_PN * uptake%P)     ! [d^{-1}]
!acc%fac4 = phy%resp

! --- relative growth rate RGR: gross production - exudation - uptake respiration --  
phy%gpp   = grossC
!phy%gpp   = fac_colim
uptake%C  = grossC - phy%resp     !* (1.0d0- self%exud_phy) ![d^{-1}]
!write (*,'(A,6(F10.4))') 'upC ',phy%frac%Rub , Pmaxc, sens%upt_pot%C,grossC,phy%resp,  uptake%C
!acc%fac4 = phy%resp

! --- photoacclimation and photosynthesis ------------------------------------            
!     differential coupling between pigment synthesis and costs due to N-uptake     
! ---  partitioning to chloroplast and rubisco  -------------------------------
!> @fn maecs_primprod::photosynthesis()
!> 9. calculate (acc\%)dfR/dt
! specific respiration, also accounting for P-uptake expenses
resp  = zeta_X(self%nutind%iN) * (upt_act%N + self%zstoich_PN * upt_act%P)  

if (self%RubiscoOn) then
! --- derivatives of C-uptake rate  ------------------------------------------         
   dmu_dfracR = Pmaxc * sens%upt_pot%C - 1*resp * dfV_dfracR 

!   dmuQ_dfracR = dmuQ_dfracR !* darkf**(0*exp(-phy%frac%Rub))

   grad_fracR = dmu_dfracR      &   ! marginal C gain of light independent processes 
              + dmuQ_dfracR         ! marginal loss due to reduced uptake 

 !  acc%fac3 = dmuQ_dfracR
!   acc%fac3 = dmu_dfracR
! upper boundary for fractional variables: smoothly integrates constrain from 2nd fractional variable
  fmaxf  = 1.0d0/(1.0d0+exp(0.0d0-2*grad_fracR/self%res0))
!  fmaxf  = 1.0d0
! upper boundary for fractional variables: resulting flexibility
  fmaxf  = max(1.0d0 - phy%frac%theta*fmaxf- phy%frac%Rub,0.0d0)

! --- regulation speed in Rubisco expression ------------------------------------ 
   flex_fracR = self%adap_Rub * fmaxf * (phy%frac%Rub-self%rel_chloropl_min)
! \todo check old version: (phy%frac%Rub - self%rel_chloropl_min-self%small_finite)

! *** ADAPTIVE EQUATION FOR 'frac_R'
   acc%dfracR_dt = flex_fracR * grad_fracR     
end if
!write (*,'(A,2(F10.3))') 'dfracR_dt: ',acc%dfracR_dt*1E3, grad_fracR *1E3


!> @fn maecs_primprod::photosynthesis()
!> 10. acc\%dtheta\_dt
!>   + reorganize \lref{eq. ,eq:ftheta,} as: @f$ \theta=\theta_C/chl_{r}* f_\theta @f$
!>     - where @f$ chl_{r}=f_R*Q_N^\sigma @f$ as set in subroutine calc_internal_states()
!>   + @f$ d\theta/dt= \theta_C/chl_{r}* df_\theta/dt  @f$
!>     - @f$ df_\theta/dt=\delta_\theta^* *f_\theta*(1-f_\theta)*(d\mu/df_\theta + dmuQ/df_\theta) @f$. \lref{see eq. ,eq:dmu_dtheta,.}
!>     - @f$ (d\mu/df_\theta + dmuQ/df_\theta) = (d\mu/d\theta+dmuQ/d\theta)*\theta_C/chl_{r}@f$
!>   + substituting, @f$ d\theta/dt= \delta_\theta^* * (\theta_C/chl_{r})^2* f_\theta*(1-f_\theta) * (d\mu/d\theta+dmuQ/d\theta)  @f$
!>     - @f$  d\mu/d\theta @f$ : \lref{see eq. ,eq:dmu_dtheta,.}
!>     - @f$  dmuQ/d\theta @f$ : see above

if (self%PhotoacclimOn) then
! --- derivatives of C-uptake rate  --------------------------------------------         
!     positive gradient term due to PAR adsorption by CHL 
   dmu_dtheta = Pmaxc* phy%frac%Rub * exp(- sens%a_light * phy%theta) *sens%a_light & 
                      -1*resp * dfV_dtheta 

!   acc%fac1 = Pmaxc* phy%frac%Rub
!   acc%fac2 = exp(- sens%a_light * phy%theta) *sens%a_light
!   acc%fac4 = -1*resp * dfV_dtheta
!   acc%fac4 = sens%a_light


   dmuQ_dtheta = dmuQ_dtheta !* darkf**(0*exp(-phy%frac%theta))

   grad_theta = dmu_dtheta + dmuQ_dtheta  ! marginal C gain and indirect costs of chloroplasts

  ! --- regulation speed in photoacclimation  -----------------------------------
  !flex_theta  = self%adap_theta * ( self%theta_LHC - phy%theta) * phy%theta
!  flex_theta  = self%adap_theta * (1- phy%frac%theta) * phy%frac%theta 

! upper boundary for fractional variables: smoothly integrates constrain from 2nd fractional variable
  fmaxf  = 1.0d0/(1.0d0+exp(0.0d0-2*grad_theta/self%res0))
!  fmaxf  = 1.0d0
! upper boundary for fractional variables: resulting flexibility
  fmaxf  = max(1.0d0 - phy%frac%theta- phy%frac%Rub*fmaxf,0.0d0)

  flex_theta  = self%adap_theta * (self%theta_LHC/phy%rel_chloropl)**2 * fmaxf * max((phy%frac%theta -self%rel_chloropl_min),0.0d0)

  ! *** ADAPTIVE EQUATION FOR 'theta'
  acc%dtheta_dt = flex_theta * grad_theta

   !for being able to save these intermediate quantities as diag vars:
!   acc%fac1 = dmu_dtheta+dmuQ_dtheta
!   acc%fac2 = phy%rel_chloropl
end if
! --- carbon exudation   -------------------------------------------------------
!  TODO: discuss and adjust; data?
exud%C      = self%exud_phy ! * grossC                        ![d^{-1}]

! --- carbon specific nitrogen & phosporus exudation   -------------------------
exud%N      = self%exud_phy * phy%Q%N ! * uptake%N   ! [(mmolN) (mmolC)^{-1} d^{-1}]
if (self%PhosphorusOn .and.  exud%N .gt. 0.0d0) then
  exud%P      =  self%exud_phy * phy%Q%P 
!  exud%P      = phy%P / phy%reg%N * exud%N  
else
  exud%P      = 0.0d0
end if

! ---- additional exudation to release unrealistic stoichiometry in depositional holes
ex = self%QN_phy_0*phy%reg%C - phy%N ! + self%small_finite*self%QN_phy_max
!  write (0,'(A,4(F9.4))') 'exC0=', self%QN_phy_0*phy%reg%C,phy%N,ex,(exp(ex * self%iK_QN/phy%reg%C)-1.0d0)
if( ex .gt. 0.0d0 ) then
  exud%C = exud%C + self%decay_nut * (exp(ex * self%iK_QN/phy%reg%C)-1.0d0) 
end if

if( phy%relQ%N .gt. maxrq ) then
  exud%N = exud%N + self%decay_nut * (exp(phy%relQ%N- maxrq)-1.0d0) * phy%Q%N
!  write (0,'(A,4(F9.4))') 'exN=',phy%relQ%N- maxrq,exp(phy%relQ%N- maxrq)-1.0d0,self%decay_nut * (exp(phy%relQ%N- maxrq)-1.0d0) * phy%Q%N,elem(self%nutind%iN)%aV 
end if

if (self%PhosphorusOn) then
  if( phy%relQ%P .gt. maxrq ) then
    exud%P = exud%P + self%decay_nut * (exp(phy%relQ%P - maxrq)-1.0d0) * phy%Q%P
  end if
end if

! set few volatile diag variables ___________________________________
!if (self%DebugDiagOn) then
!  acc%tmp    = sens%upt_pot%C
!  acc%fac1   = phy%theta * phy%rel_chloropl 
!  acc%fac2   = grad_fracR
! endif
end subroutine photosynthesis

end module maecs_primprod
