!
#include "fabm_driver.h"
!---------------------------------------------------------
!BOP
! !MODULE: MAECS_functions --- more to come
!  Model for Adaptive Ecosystems in Coastal Seas 
   module maecs_primprod

   use fabm_types
   use maecs_types
   use maecs_functions
   private
   public    photosynthesis     

 contains  

subroutine photosynthesis(self,sens,phy,uptake,exud,acc)
implicit none

type (type_maecs_base_model), intent(in)   :: self
type (type_maecs_sensitivities),intent(in) :: sens
type (type_maecs_phy), intent(inout)       :: phy
type (type_maecs_om), intent(out)          :: uptake
type (type_maecs_om), intent(out)          :: exud
type (type_maecs_traitdyn), intent(out)    :: acc

!integer, intent(in), optional          :: method
real(rk) :: ratio_rel               ! ratio of free P- & Si- over free N-reserves   (molP molN^{-1})
real(rk) :: fac_colim               ! colimitation due to 'ratio_free_PN' (dimensionless)
logical  :: IsVQP, IsVQSi           ! switching on/off additional gradient terms
integer  ::num_nut
! --- UPTAKE KINETICS
type (type_maecs_om)  :: upt_act ! actual C-specific nutrient uptake rate [molN molC^{-1} d^{-1}]
type (type_maecs_om)  :: dmu_dQ ! DERIVATIVE OF RGR WITH RESPECT TO QUOTAS  
real(rk) :: reg_V  ! 'fast trait': regulation of nutrient uptake
real(rk) :: gross_nue  ! gross carbon uptake rate per relative nutrient fraction ("Nitrogen Use Efficiency") [d^{-1}]
! ******************************************************************************  

! --- DERIVATIVES OF N-UPTAKE (FOR TRAIT-BASED ADAPTIVE DYNAMICS)
real(rk) :: dVNC_dfracR ! [...]
real(rk) :: dVNC_dQN    ! [...]
real(rk) :: dVNC_dtheta ! [...]

! --- DERIVATIVES OF P-UPTAKE (FOR TRAIT-BASED ADAPTIVE DYNAMICS)
real(rk) :: dVPC_dtheta, dV_dregV
real(rk) :: dVPC_dfracR ! [...]
real(rk) :: deriv_fac_colim   

! --- DERIVATIVES OF C-UPTAKE (FOR TRAIT-BASED ADAPTIVE DYNAMICS)
real(rk) :: dmu_dtheta  ! [...]
real(rk) :: dmu_dfracR  ! [...]

! --- DERIVATIVES FOR TRADE-OFF & ADAPTATION
real(rk) :: eps, rel_phys
real(rk) :: dmu, dmu_small       ! [...] 
real(rk) :: dQN_dtheta ! [...]
real(rk) :: dQN_dfracR ! [...]

real(rk) :: dmuP       ! [...]
real(rk) :: dQP_dtheta, dQP_dfracR ! [...]
real(rk) :: grossC, lossC
real(rk) :: dmudV
real(rk) :: grad_V
real(rk) :: dbal_dv
real(rk) :: lim_N, lim_P, lim_Si, lim_all
real(rk) :: grad_theta ! [...]
real(rk) :: grad_fracR ! [...]
real(rk) :: flex_theta ! [...]
real(rk) :: flex_fracR ! [...]
real(rk) :: a1, a2, a3
real(rk) :: feedb_vq
real(rk) :: syn_act, zeta_CP
real(rk) :: f_Lip, q_Lip, q_NoLip  
    
! 
eps     =  self%small_finite ! just  a shorter namer
! switching off phosphorus terms in the derivatives
IsVQP       = .true.  ! rethink dependence of P-uptake on N-regulation switching on phosphorus terms in the derivatives
!IsVQP       = .false. 
IsVQSi      = .false. ! TODO

! --- relative amount of carbon invested into light harvesting complex (LHC) -------
! chlorophyll-to-carbon ratio of chloroplast * chloroplast concentration of cell
if (self%PhotoacclimOn) then  
   acc%dRchl_dtheta = phy%rel_chloropl 
   acc%dRchl_dfracR = phy%theta * phy%rel_QN**self%sigma
   acc%dRchl_dQN    = phy%theta * phy%frac%Rub * self%sigma * phy%rel_QN**(self%sigma - 1.0d0) * self%iK_QN
endif

dbal_dv = 1.0d0 + phy%QN * self%zeta_CN

! --- gross carbon uptake by phytoplankton (gross-primary production) ---------------
!   phy%frac%Rub* sens%P_max_T* phy%rel_QN : carboxylation capacity = Rub * free proteins
!   fac_colim             : synchrony in  protein/RNA dynamics 
!   sens%S_phot           : light harvesting (light limited growth)   
 
! here, every possible nutrient is asked explicitely; thus first Si then P
! TODO: loop over structure elements with arbitrary number of nutrients
num_nut  = 0

! -- carbon specific Si-uptake with N, Si, and P co-limitation, according to queue model -
if (self%SiliconOn) then 
   ratio_rel  = phy%rel_QSi
   lim_Si     = 1.0d0
   num_nut    = num_nut + 1
else
   lim_Si     = 0
endif

! --- synchrony in nutrient assimilation depends on growth cycle and
!         P-quota
if (self%syn_nut .lt. _ZERO_ .and. self%PhosphorusOn) then
   syn_act = -self%syn_nut * phy%QP / self%QP_phy_max
else
   syn_act = self%syn_nut
endif 

if (self%PhosphorusOn) then 
  if (mod(num_nut,2) .eq. 1) then 
     ratio_rel = ratio_rel/phy%rel_QP

     call queuefunc(syn_act,ratio_rel,fac_colim,deriv_fac_colim)
 ! relative nutrient limitation due to P shortage 
     fac_colim = smooth_small(fac_colim,eps)
     ! relative nutrient limitation due to Si shortage
     lim_Si    = deriv_fac_colim*ratio_rel / fac_colim
     ! relative nutrient limitation due to P shortage
     lim_P     = smooth_small(1.0d0 - lim_Si, eps)
  else
     fac_colim = 1.0d0
     lim_P     = 1.0d0
  endif 
! contribution of P-status to gross assimilation
  ratio_rel = phy%rel_QP * fac_colim
  num_nut   = num_nut + 1
! uptake dependency on nutrient availability
else
  lim_P     = 0
endif

! -- always account for N-quota -
if (num_nut .gt. 0) then 
! takes co-limitation factor from all other nutrients
   ratio_rel  = ratio_rel/phy%rel_QN
   call queuefunc(syn_act,ratio_rel,fac_colim,deriv_fac_colim)
   fac_colim  = smooth_small(fac_colim,eps)
! deriv_fac_colim: numerical derivative of queue response  
   lim_all    = deriv_fac_colim*ratio_rel / fac_colim
! relative nutrient limitation due to N shortage
   lim_N      = smooth_small(1.0d0 - lim_all, eps)
   lim_P      = lim_P * lim_all
   lim_Si     = lim_Si* lim_all


else
   fac_colim   = 1.0d0
   deriv_fac_colim = 0.0d0
   lim_N       = 1.0d0

end if
! relative nutrient limitation due to P shortage

gross_nue   = fac_colim * phy%frac%Rub * sens%P_max_T * sens%S_phot  !Nitrogen use efficiency
grossC      = phy%rel_QN * gross_nue  ! [d^{-1}]

! ---------------------   derivatives of reg_V   -----------------------------------
dV_dregV    = phy%frac%NutUpt * sens%up_NC  !corrected Apr 17 kw

! +++ derivative of C-uptake rate with respect to N-quota ++++++++++++++++++++++++++++++
dmu_dQ%N    = lim_N * gross_nue * self%iK_QN 

a3          = phy%rel_QN 
a2          = lim_N + (dbal_dv -1.0d0)*(lim_N + a3)
a1          = self%zeta_CN*(phy%QN - self%QN_phy_0) + lim_N*phy%frac%NutUpt/(phy%frac%Rub + eps)* &
               self%K_QN_phy * (self%zeta_CN + 1.0d0/phy%QN)

! --- relative growth rate RGR: gross production - exudation - uptake respiration --  

feedb_vq    = 1.0d0 - dbal_dv * a3 / (1.0d0 - self%QN_phy_0/phy%QN + a2 + a3)
dmudV       = dbal_dv * feedb_vq/((phy%QN - self%QN_phy_0) * (1.0d0/a1 + 1.0d0/a2) + phy%QN )  

! approximated 1/(mu/(dmudQ*Q) +1 - dVdQ)
!   steady-state down-regulation of uptake I: balance of respiration and indirect benefits   
grad_V      = dV_dregV * (-self%zeta_CN + dmudV) 
!   smoothed version of step function, uses marginal gain to emulate a continuous response 
reg_V       = 1./(1.+exp(-self%tau_regV * grad_V));  ! 0.02

!------------------------------------------------------------------------
!     adjust all terms for down-regulated uptake      
!-------------------------------------------------------------------------
upt_act%N   = reg_V * sens%up_NC
uptake%N    = phy%frac%NutUpt * upt_act%N  ! [(mmolN) (mmolC)^{-1} d{-1}]
! ---  respiration due to N assimilation --------------------------------------
lossC       = self%zeta_CN * uptake%N                           ! [d^{-1}]

! --- relative growth rate RGR: gross production - exudation - uptake respiration --  
uptake%C    = grossC - lossC     !* (1.0d0- self%exud_phy) ![d^{-1}]
!  write (*,'("  fagN",5(F9.4))') par,,up_NC,mu,din

! exud propto uptake%grossC^2 ???
! +++ derivatives of N-uptake rate ++++++++++++++++++++++++++++++++++++++++++
dVNC_dfracR = -upt_act%N * (1.0d0 + phy%rel_QN**self%sigma * phy%theta * self%itheta_max)                        

dVNC_dQN    = -upt_act%N * acc%dRchl_dQN * self%itheta_max
dVNC_dtheta = -upt_act%N * phy%rel_chloropl * self%itheta_max 

! update dmu_dQ%N
dmu_dQ%N     = dmu_dQ%N - self%zeta_CN * dVNC_dQN 
!      if (.not. self%PhotoacclimOn) exit ! break the loop for constant theta and frac_R
! --- derivatives of P-uptake rate ----------------------------------------------------------------          
if (self%PhosphorusOn) then
    dV_dregV = phy%frac%NutUpt * sens%up_PC
    dmudV    = dbal_dv**2 / ((dbal_dv*lim_P+1.0d0) *phy%QP - self%QP_phy_0)  
 !TODO check costs !
! prelim solution: stoichiometry in RNA (N:P ~ 4:1) and phospholipids (N:P~1:1)
! TODO: include proteins/mebranes (N:P >> 16:1) under low growth conditions 
! TODO: energetic costs of P-assimilation not brought up as an extra term but 
!          assumed to be already included in protein synthesis
! \partial (\zeta_CN V_N) / \partial V_P
    q_NoLip  = 1./3.8  ! P-stoichiometry of active compounds (DNA, RNA)  
    q_Lip    = 1./0.8    ! storage P-stoichiometry 
    
    f_Lip    = 1./(1.+exp(10*(1.-phy%rel_QP)))
    zeta_CP  = ((1.-f_Lip)*q_NoLip + f_Lip*q_Lip) * self%zeta_CN

    grad_V   = dV_dregV * (-zeta_CP  + lim_P * dmudV) 
    reg_V    = 1./(1.+exp(-self%tau_regV * grad_V));  ! 0.02
    upt_act%P= reg_V * sens%up_PC 

! cell physiology regulating P - uptake is identical to N-uptake machinery
! TODO: explain, why (cf. Smith et al. 2009) 
! final, realized uptake = actual uptake rate per site times relative machinery
    uptake%P = phy%frac%NutUpt * upt_act%P    ! [(mmolP) (mmolC)^{-1} d{-1}]
    dmu_dQ%P = lim_P * gross_nue  * self%iK_QP
else
    dmu_dQ%P = 0.0d0
end if 

! --- derivatives of Si-uptake rate ----------------------------------------------------------------          
if (self%SiliconOn) then
!  uptake dependency on nutrient availability
    dV_dregV = phy%frac%NutUpt * sens%up_SiC

    dmudV    = dbal_dv**2 / ((dbal_dv*lim_Si+1.0d0) *phy%QSi - self%QSi_phy_0)  
    grad_V   = dV_dregV * (-0.*self%zeta_CN   + lim_Si * dmudV)
    reg_V    = 1./(1.+exp(-self%tau_regV * grad_V));  ! 0.02
    upt_act%S= reg_V * sens%up_SiC 

! final, realized uptake = actual uptake rate per site times relative machinery
    uptake%S = phy%frac%NutUpt * upt_act%S    ! [(mmolSi) (mmolC)^{-1} d{-1}]
    dmu_dQ%S = lim_Si * gross_nue /(self%QSi_phy_max-self%QSi_phy_0)
else
    dmu_dQ%S = 0.0d0
end if 

! --- derivatives of C-uptake rate  --------------------------------------------         
rel_phys    = phy%rel_QN * fac_colim 
dmu_dtheta  = phy%frac%Rub * rel_phys * sens%P_max_T * sens%a_light * (1.0d0 - sens%S_phot) & 
                      - self%zeta_CN * dVNC_dtheta               !OK

rel_phys= rel_phys * sens%S_phot 
!write (*,'(" phys ",4(F9.4))') rel_phys, phy%rel_QN, fac_colim, sens%S_phot 
phy%frac%rel_phys=rel_phys
dmu_dfracR  = rel_phys * sens%P_max_T - self%zeta_CN * dVNC_dfracR 

! +++ trade-off: effects on N uptake +++++++++++++++++++++++++++++++++++++++++++
! dmu is a central regulation parameter determining the coupling strength in ressource uptake; its position in the denominator makes it critical for frequently occuring low values
! TODO: replace with more robust formulation                                                
dmu         = dmu_dQ%N * phy%QN + uptake%C - dVNC_dQN
! rough estimate of maximal growth rate with frac_R=0.5 and q=1
dmu_small   = eps * sens%P_max_T 
dmu         = smooth_small(dmu , dmu_small); 

dQN_dtheta  = (dVNC_dtheta -1* dmu_dtheta * phy%QN) / dmu
dQN_dfracR  = (dVNC_dfracR -1* dmu_dfracR * phy%QN) / dmu
!dQN_dtheta  = (1.0d0  + 0*dQN_dtheta * dVNC_dtheta) * dQN_dtheta
!dQN_dfracR  = (1.0d0  + 0*dQN_dfracR * dVNC_dfracR) * dQN_dfracR

! --- photoacclimation and photosynthesis ------------------------------------------------            
!     differential coupling between pigment synthesis and costs due to N-uptake     
! positive gradient term due to PAR adsorption by CHL 
grad_theta  = dmu_dtheta + dmu_dQ%N * dQN_dtheta

! --- trade-off: effects on P uptake -----------------------------------------------------         
if (self%PhosphorusOn .and. IsVQP ) then
   dmuP       = dmu_dQ%P * phy%QP + uptake%C
   dmuP       = smooth_small(dmuP , dmu_small); 
! rough estimate of maximal growth rate with frac_R=0.5 and q=1

   dVPC_dfracR = -upt_act%P * (1.0d0 + phy%rel_QN**self%sigma * phy%theta * self%itheta_max)                        
   dVPC_dtheta = -upt_act%P * phy%rel_chloropl * self%itheta_max 

   dQP_dtheta = (dVPC_dtheta - 0.0d0* dmu_dtheta * phy%QP) /dmuP
   dQP_dfracR = (dVPC_dfracR - 0.0d0* dmu_dfracR * phy%QP) /dmuP
!   dQP_dtheta = (1.0d0  + 0*dQP_dtheta * dVPC_dtheta) * dQP_dtheta
!   dQP_dfracR = (1.0d0  + 0*dQP_dfracR * dVPC_dfracR) * dQP_dfracR
   grad_theta = grad_theta + dmu_dQ%P * dQP_dtheta

else
   dQP_dfracR = 0.0d0
end if ! PhosphorusOn

! --- nitrogen partitioning to chloroplast and rubisco --------------------------------------------
if (self%PhotoacclimOn .and. self%RubiscoOn) then
   grad_fracR = dmu_dfracR                   & ! marginal C gain of chloroplasts 
              + dmu_dQ%N * dQN_dfracR          ! marginal loss due to reduced uptake 
   if(self%PhosphorusOn .and. IsVQP) then
     grad_fracR = grad_fracR+ dmu_dQ%P * dQP_dfracR  ! marginal loss due to reduced uptake 
   end if !!!   acc%tmp,fac1,fac2 dmu_dQ%P

   acc%fac1 = dVPC_dfracR
   acc%fac2 = dQP_dfracR
   acc%tmp  = dmuP
! rubisco-N directly limits chl synthesisdmu_dQ%P * dQP_dfracR
! --- regulation speed in Rubisco expression -------------------------------------- 
   flex_fracR = self%adap_Rub * (1.0d0 - phy%frac%Rub ) * (phy%frac%Rub - self%rel_chloropl_min-self%small_finite)
! *** ADAPTIVE EQUATION FOR 'frac_R'
   acc%dfracR_dt = flex_fracR * grad_fracR  

!write (*,'(A,2(F10.3))') '2:',acc%dfracR_dt*1E3, grad_fracR *1E3
else
   grad_fracR = 0.0d0
end if

! --- regulation speed in photoacclimation   --------------------------------------
flex_theta  = self%adap_theta * ( self%theta_LHC - phy%theta ) * phy%theta

! *** ADAPTIVE EQUATION FOR 'theta'
acc%dtheta_dt = flex_theta * grad_theta
! *********************************

! --- carbon exudation   ----------------------------------------------------------
!  TODO: discuss and adjust; data?
exud%C      = self%exud_phy * grossC                        ![d^{-1}]

! --- carbon specific nitrogen & phosporus exudation   ----------------------------
exud%N      = self%exud_phy * uptake%N   ! [(mmolN) (mmolC)^{-1} d^{-1}]
exud%P      = phy%QPN * exud%N   ! [(mmolP) (mmolC)^{-1} d^{-1}]

! set few volatile diag variables ___________________________________
if (self%DebugDiagOn) then
  acc%tmp    = sens%S_phot
  acc%fac1   = phy%theta * phy%rel_chloropl 
  acc%fac2   = grad_fracR
endif

end subroutine photosynthesis

end module maecs_primprod
