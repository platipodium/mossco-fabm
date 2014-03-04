!> @file maecs_primprod.F90
!> @author Richard Hofmeister, Markus Schartau, Kai Wirtz, Onur Kerimoglu

#include "fabm_driver.h"

!> @brief Primary production module
   module maecs_primprod

   use fabm_types
   use maecs_types
   use maecs_functions
   private
   public    photosynthesis     

 contains  

!> @brief  calculates grazing rate
!> @details 
!> This is the subroutine, where the optimal regulation of phytoplankton traits
!> are described, which is central to the physiological-MAECS
subroutine photosynthesis(self,sens,phy,uptake,exud,acc)
implicit none

type (type_maecs_base_model), intent(in)   :: self
type (type_maecs_sensitivities),intent(in), target :: sens
type (type_maecs_phy), intent(inout)       :: phy
type (type_maecs_om), intent(out), target  :: uptake
type (type_maecs_om), intent(out)          :: exud
type (type_maecs_traitdyn), intent(out), target    :: acc

integer  :: num_nut, i
integer  :: j, j2
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
real(rk) :: grossC, lossC, Pmaxc
real(rk) :: dbal_dv
real(rk) :: e_N0, e_N
real(rk) :: syn_act, hh
real(rk) :: f_Lip, q_Lip, q_NoLip 
real(rk) :: c_h1, c_hq           ! correction terms of queuing function 
real(rk) :: fac_colim            ! colimitation factor (dimensionless)
real(rk) :: q_X, qp_Y, qp_X      ! placeholder for nutrient availability
real(rk) :: qfunc, deriv_qfunc   ! queuing function and derivative
real(rk) :: prod_dq, d_X, d_QX, dmu_daV, act_V
real(rk), dimension(5) :: sigmv ! relative quota-uptake feed-back 
real(rk), dimension(5) :: dqp_X_dq_X, dqp_X_dqp_Y  ! derivatives
real(rk), dimension(5) :: zeta_X  ! C-costs of assimilation of nutrient X
!real(rk), dimension(5) :: upt_act
type (stoichiometry_pointer), dimension(5) :: elem ! pointer structure for addressing all elements wthin loops
! TODO: self%num_nut
! 
eps     =  self%small_finite ! just  a shorter namer


!> @fn maecs_primprod::photosynthesis()
!> *here some in-body-doc*
!> Prepare loop over structure elements by assigning a poinzer structure
!! here, every possible nutrient is asked explicitely; thus first Si then P
!>
!> First, the following is included from maecs_stoichvars.F90p:
!> @include 'maecs_stoichvars.F90p'
include './maecs_stoichvars.F90p'

num_nut  = i
! write (*,'(A,(I4))') 'num_nut:',num_nut

! prelim solution: elemiometry in RNA (N:P ~ 4:1) and phospholipids (N:P~1:1)
! TODO: include proteins/mebranes (N:P >> 16:1) under low growth conditions 
! TODO: energetic costs of P-assimilation not brought up as an extra term but 
!          assumed to be already included in protein synthesis
! \partial (\zeta_CN V_N) / \partial V_P
!    q_NoLip  = 3.8  ! P-stoichiometry of active compounds (DNA, RNA)  
    q_NoLip  = self%zstoich_PN
!    q_NoLip  = 16   ! Redfield-stoichiometry 
    q_Lip    = 0.8  ! storage P-stoichiometry   
    f_Lip    = 1./(1.+exp(10*(1.-phy%relQ%P)))

zeta_X(1)         = self%zeta_CN
if (self%PhosphorusOn) then
   if (num_nut .gt. 2) zeta_X(3:num_nut) = 0.0d0
   zeta_X(2)      = self%zeta_CN * ((1.-f_Lip)*q_NoLip + f_Lip*q_Lip)
!   zeta_X(2)      = self%zeta_CN * 0
end if
! --- relative amount of carbon invested into light harvesting complex (LHC) -------
! chlorophyll-to-carbon ratio of chloroplast * chloroplast concentration of cell
!if (self%PhotoacclimOn) then  
acc%dRchl_dtheta = phy%rel_chloropl 
acc%dRchl_dfracR = phy%theta * phy%relQ%N**self%sigma
acc%dRchl_dQN    = phy%theta * phy%frac%Rub * self%sigma  * self%iK_QN
                     ! * phy%relQ%N**(self%sigma-1)
!endif

dbal_dv = 1.0d0 + phy%Q%N * self%zeta_CN  ! partial derivative of the balance eq. to uptake V
                                          !   equals ratio of gross to net primary production

! --- gross carbon uptake by phytoplankton (gross-primary production) ---------------
!   phy%frac%Rub* sens%P_max_T* phy%relQ%N : carboxylation capacity = Rub * free proteins
!   fac_colim             : synchrony in  protein/RNA dynamics 
!   sens%upt_pot%C           : light harvesting (light limited growth)   

! --- synchrony in nutrient assimilation depends on growth cycle and N-quota
!         
if (self%syn_nut .le. _ZERO_ ) then  !
   syn_act = -self%syn_nut * phy%relQ%N
else
   syn_act  = self%syn_nut
endif 
acc%tmp   =   syn_act ! store

! $\climf = q_1\cdot g_h (q_2'/q_1) \cdot (1+h\cdot q_1 q_2' +c_h)$
! metabolic interdependence h is inverse of synergy syn_act
hh      = 1.0d0 / (syn_act + eps)

! one plus correction coefficient $c_h$ to ensure convergence to Liebig and product rule
c_h1 = 1.0d0  + log( 1.0d0/(4**hh) + 0.5*hh );

!  auxiliary variable expressing the ratio between gross and net production
dbal_dv = 1.0d0 + phy%Q%N * self%zeta_CN

e_N0    = 1.0d0 / (phy%Q%N * dbal_dv) 

qp_Y    = elem(num_nut)%relQ
dqp_X_dq_X(num_nut)  = 1.0d0
dqp_X_dqp_Y(num_nut) = 0.0d0

! retro-loop over structure elements with arbitrary number of nutrients
!   to calculate a process-based co-limitation factor and derivative terms needed for 
!   extended optimality functions

if (num_nut .gt. 1) then
 sigmv(2:num_nut)  = 0.0d0

 ! loop over all nutrients starting from pre-final
 do i = num_nut-1, 1, -1  

   q_X   = elem(i)%relQ

! correction of queuing function for symmetric/interdependent co-limitation 
   c_hq  = c_h1 + hh * qp_Y * q_X

! qfunc: relative nutrient limitation due to shortage in other nutrients 
   call queuefunc(syn_act, qp_Y / q_X, qfunc, deriv_qfunc)

! $d_\mathrm{N} = \fracd{1}{\fQ{N}} - \pdiff{g_h}{x}\fracd{q_2'}{g_h\fQ{N}^2} +\fracd{h q_2'}{c_{hq}}$

   qp_X  = q_X * qfunc * c_hq

! +++ derivative of C-uptake rate with respect to quota ++++++++++++++++++++++++++++++
   qfunc = smooth_small(qfunc,eps)
   q_X   = smooth_small(q_X,eps)

   dqp_X_dq_X(i)  = qp_X* ( 1.0d0/q_X + qp_Y *( hh/c_hq - deriv_qfunc/(qfunc*q_X*q_X) ))

   dqp_X_dqp_Y(i) = qp_X *( hh*q_X/c_hq + deriv_qfunc/(qfunc*q_X) )

   qp_Y  = qp_X

 end do
else
 qp_X = qp_Y   ! initial value for num_nut=1
! write (*,'(A,3(F11.5))') 'qp1:',phy%Q%N,elem(num_nut)%relQ,phy%relQ%N

end if

fac_colim   = qp_X

phy%rel_phys= fac_colim * sens%upt_pot%C  ! auxiliary variable for sinking routine

! recursive product following from chain rule; first term : 1/LF
prod_dq     = 1.0d0/smooth_small(qp_X,eps)

! derivative of f_V (fraction of nutrient uptake machinery) on f_R (PS) and theta
dfV_dfracR  = - acc%dRchl_dfracR * self%itheta_max - 1.0d0 
dfV_dtheta  = - acc%dRchl_dtheta * self%itheta_max

!  auxiliary variable for the feed-back of N-quota and N-uptake change
sigmv(1)    =   acc%dRchl_dQN  * self%itheta_max / phy%frac%NutUpt

! $\partial mu/\partial Q dQ/dV|_{tot}$ : initial zero berfore loop
dmuQ_dfracR = 0.0d0
dmuQ_dtheta = 0.0d0

do i = 1, num_nut 
 ! ------------------ derivatives of co-limited growth function ---------------------
   d_X       = prod_dq * dqp_X_dq_X(i)
   prod_dq   = prod_dq * dqp_X_dqp_Y(i)
   e_N       = e_N0 + d_X * elem(i)%iKQ / elem(i)%relQ

   d_QX      = d_X* dbal_dv * elem(i)%iKQ + sigmv(i)* phy%Q%N * self%zeta_CN  ! 

!   steady-state down-regulation of uptake I: balance of respiration and indirect benefits  
   dmu_dV    = (1.0d0 + zeta_X(i) * elem(i)%Q) * d_QX/(1.0d0 + elem(i)%Q * (d_QX + sigmv(i)))

   dmu_dV    = dmu_dV * e_N / (e_N + sigmv(i))

!   steady-state down-regulation of uptake I: balance of respiration and indirect benefits  
   dmu_daV   = (-zeta_X(i) + dmu_dV) * phy%frac%NutUpt * elem(i)%upt_pot 

!   smoothed version of step function, uses marginal gain to emulate a continuous response 
   act_V     = 1./(1.+exp(-self%tau_regV * dmu_daV));  ! 0.02

   elem(i)%aV      = act_V 
   elem(i)%upt_act = act_V * elem(i)%upt_pot
   elem(i)%upt     = phy%frac%NutUpt * elem(i)%upt_act  ! [(molX) (molC)^{-1} d{-1}]

! +++ derivative of C-uptake rate with respect to quota ++++++++++++++++++++++++++++++
   dmuQ_dfracR     = dmuQ_dfracR + dmu_dV * elem(i)%upt_act * dfV_dfracR
   dmuQ_dtheta     = dmuQ_dtheta + dmu_dV * elem(i)%upt_act * dfV_dtheta

!write (*,'(A,1(I2),5(F9.3))') 'df:',i,dmu_dV , phy%theta,phy%relQ%N**self%sigma, dfV_dfracR,dmuQ_dfracR
end do

! --- gross carbon uptake by phytoplankton (gross-primary production) ---------------
!   phy%frac%Rub* sens%P_max_T* phy%rel_QN : carboxylation capacity = Rub * free proteins
!   fac_colim       : final synchrony in  protein/RNA dynamics (multi-nutrient processing)
!   sens%upt_pot%C  : light harvesting (light limited growth)
Pmaxc     =  fac_colim * sens%P_max_T
grossC    =  phy%frac%Rub * Pmaxc * sens%upt_pot%C  ! primary production

! ---  respiration due to N assimilation --------------------------------------
lossC     = self%zeta_CN * uptake%N                           ! [d^{-1}]

! --- relative growth rate RGR: gross production - exudation - uptake respiration --  
uptake%C  = grossC - lossC     !* (1.0d0- self%exud_phy) ![d^{-1}]

! --- photoacclimation and photosynthesis ------------------------------------            
!     differential coupling between pigment synthesis and costs due to N-uptake     
! ---  partitioning to chloroplast and rubisco  -------------------------------
if (self%RubiscoOn) then
! --- derivatives of C-uptake rate  ------------------------------------------         
   dmu_dfracR = Pmaxc * sens%upt_pot%C - self%zeta_CN * upt_act%N * dfV_dfracR

   grad_fracR = dmu_dfracR      &   ! marginal C gain of light independent processes 
              + dmuQ_dfracR         ! marginal loss due to reduced uptake 

   acc%fac1 = dmu_dfracR
   acc%fac2 = dmuQ_dfracR

! --- regulation speed in Rubisco expression ------------------------------------ 
   flex_fracR = self%adap_Rub * (1.0d0 - phy%frac%Rub ) * phy%frac%Rub
! \todo check old version: (phy%frac%Rub - self%rel_chloropl_min-self%small_finite)

! *** ADAPTIVE EQUATION FOR 'frac_R'
   acc%dfracR_dt = flex_fracR * grad_fracR  
end if
!write (*,'(A,2(F10.3))') 'dfracR_dt: ',acc%dfracR_dt*1E3, grad_fracR *1E3

if (self%PhotoacclimOn) then

! --- derivatives of C-uptake rate  --------------------------------------------         
!     positive gradient term due to PAR adsorption by CHL 
   dmu_dtheta = Pmaxc* phy%frac%Rub * (1.0d0-sens%upt_pot%C)*sens%a_light & 
                      - self%zeta_CN * upt_act%N * dfV_dtheta

   grad_theta = dmu_dtheta + dmuQ_dtheta  ! marginal C gain and indirect costs of chloroplasts

  ! --- regulation speed in photoacclimation  -----------------------------------
   flex_theta  = self%adap_theta * ( self%theta_LHC - phy%theta ) * phy%theta

  ! *** ADAPTIVE EQUATION FOR 'theta'
   acc%dtheta_dt = flex_theta * grad_theta

end if
! --- carbon exudation   -------------------------------------------------------
!  TODO: discuss and adjust; data?
exud%C      = self%exud_phy * grossC                        ![d^{-1}]

! --- carbon specific nitrogen & phosporus exudation   -------------------------
exud%N      = self%exud_phy * uptake%N   ! [(mmolN) (mmolC)^{-1} d^{-1}]
exud%P      = phy%P / phy%reg%N * exud%N   ! [(mmolP) (mmolC)^{-1} d^{-1}]

! set few volatile diag variables ___________________________________
if (self%DebugDiagOn) then
!  acc%tmp    = sens%upt_pot%C
!  acc%fac1   = phy%theta * phy%rel_chloropl 
!  acc%fac2   = grad_fracR
endif

end subroutine photosynthesis

end module maecs_primprod