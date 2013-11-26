!
#include "fabm_driver.h"
!---------------------------------------------------------
!BOP
!
! !MODULE: MAECS_functions --- more to come
!  Model for Adaptive Ecosystems in Coastal Seas 
   module maecs_primprod

! !USES:
   use fabm_types
!   use fabm_driver
   use fabm_hzg_maecs
   use maecs_types
   use maecs_functions
   private
   public    photosynthesis     

! --- local fixed parameter 
REALTYPE, parameter  :: n_queue   = 2.d0

 contains  

subroutine photosynthesis(self,sens,phy,uptake,exud,acc)
implicit none

type (type_hzg_maecs), intent(in)          :: self
type (type_maecs_sensitivities),intent(in) :: sens
type (type_maecs_phy), intent(inout)       :: phy
type (type_maecs_om), intent(out)          :: uptake
type (type_maecs_om), intent(out)          :: exud
type (type_maecs_traitdyn), intent(out)    :: acc

!integer, intent(in), optional          :: method
REALTYPE :: ratio_free_PN ! ratio of free P- over free N-reserves '(phy%rel_QP/phy%rel_QN)'     [molP molN^{-1}]
REALTYPE :: IsVQP         ! switching on/off additional gradient terms for phosphorus      [dimensionless]  
REALTYPE :: fac_PN_colim ! colimitation due to 'ratio_free_PN'                            [dimensionless]

! --- UPTAKE KINETICS

!REALTYPE :: V_NC_max_T ! temperature dependent, carbon-specific maximum N-uptake   [molN molC^{-1} d^{-1}]  
!REALTYPE :: V_PC_max_T ! temperature dependent, carbon-specific maximum P-uptake   [molP molC^{-1} d^{-1}]   
REALTYPE :: upt_N  ! kinematic, carbon-specific N-uptake rate                  [molN molC^{-1} d^{-1}]
REALTYPE :: upt_P      ! kinematic, carbon-specific P-uptake rate                  [molP molC^{-1} d^{-1}]
REALTYPE :: reg_V_NC   ! regulation of 'up_NC'                                             [dimensionless]
REALTYPE :: reg_V_PC   ! regulation of 'upt_P'                                             [dimensionless]
!REALTYPE :: V_NC       ! actual carbon-specific N-uptake rate                      [molN molC^{-1} d^{-1}]
!REALTYPE :: V_PC       ! actual carbon-specific P-uptake rate                      [molP molC^{-1} d^{-1}]
REALTYPE :: gross_nue  ! gross carbon uptake rate per relative nutrient fraction ("Nitrogen Use Efficiency") [d^{-1}]
!REALTYPE :: mu         ! relative growth rate                                                     [d^{-1}]
!REALTYPE :: resp,rest  ! respiration rate
! *************************************************************************************************  
! --- DERIVATIVES OF N-UPTAKE (FOR TRAIT-BASED ADAPTIVE DYNAMICS)

REALTYPE :: dVNC_dfracR ! [...]
REALTYPE :: dVNC_dfracP ! [...]
REALTYPE :: dVNC_dregV, dVPC_dregV ! [...]
REALTYPE :: dVNC_dQN    ! [...]
REALTYPE :: dVNC_dtheta ! [...]
! DERIVATIVE OF RELATIVE GROWTH RATE WITH RESPECT TO N:C QUOTA OF PHYTOPLANKTON 
REALTYPE :: dmu_dQN     ! [...]
REALTYPE :: dQ_dV
REALTYPE :: tmp
! --- DERIVATIVES OF P-UPTAKE (FOR TRAIT-BASED ADAPTIVE DYNAMICS)

REALTYPE :: dVPC_dtheta
REALTYPE :: dVPC_dfracR ! [...]
REALTYPE :: dVPC_dQN    ! [...] 
REALTYPE :: deriv_fac_PN_colim   
! DERIVATIVE OF RELATIVE GROWTH RATE WITH RESPECT TO P:C QUOTA OF PHYTOPLANKTON 
REALTYPE :: dmu_dQP     ! [...] 
! --- DERIVATIVES OF C-UPTAKE (FOR TRAIT-BASED ADAPTIVE DYNAMICS)
  
REALTYPE :: dmu_dtheta  ! [...]
REALTYPE :: dmu_dfracR  ! [...]
REALTYPE :: dmu_dfracPN_colim ! [...]
REALTYPE :: dmu_dregV   ! [...]

! --- DERIVATIVES FOR TRADE-OFF & ADAPTATION
REALTYPE :: eps
REALTYPE :: dmu, dmu_small       ! [...] 
REALTYPE :: dQN_dtheta ! [...]
REALTYPE :: dQN_dfracR ! [...]
REALTYPE :: dQN_dregV ! [...]
REALTYPE :: dQN_dfracP ! [...]

REALTYPE :: dmuP       ! [...]
REALTYPE :: dQP_dtheta ! [...]
REALTYPE :: dQP_dfracR ! [...]
REALTYPE :: dQP_dregV ! [...]
REALTYPE :: dQP_dfracP ! [...]
REALTYPE :: grossC,lossC
REALTYPE :: grad_VN, grad_VP     ! [...]
REALTYPE :: dbal_dv
REALTYPE :: lim_N, lim_P
REALTYPE :: dmudVN, dmudVP
REALTYPE :: grad_theta ! [...]
REALTYPE :: flex_theta ! [...]
REALTYPE :: a1, a2, a3
REALTYPE :: feedb_vq
REALTYPE :: phyCR, phyNR, phyPR
REALTYPE :: grad_fracR ! [...]
REALTYPE :: flex_fracR ! [...]
! 
eps     =  self%small_finite ! just  a shorter namer
! switching off phosphorus terms in the derivatives
IsVQP       = 1.0d0  ! rethink dependence of P-uptake on N-regulation switching on phosphorus terms in the derivatives
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
!   fac_PN_colim          : synchrony in  protein/RNA dynamics 
!   sens%S_phot           : light harvesting (light limited growth)    

! -- carbon specific P-uptake with N- and P co-limitation, according to queue model -
if (self%PhosphorusOn) then 
   ratio_free_PN= phy%rel_QP /phy%rel_QN       ! ratio of free P- over free N-reserves
        ! N_queue intermittency in protein/RNA dynamics creates weak co-limitation   
!          fac_PN_colim       = queuefunc(n_queue,ratio_free_PN,quefunc,queuederiv)
   call queuefunc(n_queue,ratio_free_PN,fac_PN_colim,deriv_fac_PN_colim)
          ! numerical derivative of ingestion response  
  ! deriv_fac_PN_colim = queuederiv(n_queue,ratio_free_PN) 
! relative nutrient limitation due to P shortage 
   fac_PN_colim = smooth_small(fac_PN_colim,eps)
   lim_P        = deriv_fac_PN_colim*ratio_free_PN / fac_PN_colim
   dVPC_dregV   = phy%frac%NutUpt * sens%up_PC
   dmudVP       = dbal_dv**2 / ((dbal_dv*lim_P+1.0d0) *phy%QP - self%QP_phy_0)  
  ! TODO smooth_small for Q_P -> Q_P0 ??
else
   fac_PN_colim = 1.0d0
   phy%rel_QP   = 1.0d0
   deriv_fac_PN_colim = 0.0d0
   dVPC_dregV   = 0.0d0
   lim_P        = 0.0d0
end if

! relative nutrient limitation due to N shortage
lim_N       = smooth_small(1.0d0 - lim_P, eps)

gross_nue   = phy%frac%Rub * sens%P_max_T * sens%S_phot  !Nitrogen use efficiency
phy%rel_QNP = phy%rel_QN * fac_PN_colim
grossC = phy%rel_QNP * gross_nue  ! [d^{-1}]


! ---------------------   derivatives of reg_V_NC   -----------------------------------
dVNC_dregV  = phy%frac%NutUpt * sens%up_NC  !corrected Apr 17 kw

! +++ derivative of C-uptake rate with respect to N-quota ++++++++++++++++++++++++++++++
! deriv_fac_PN_colim wahrscheinlich 'dmu_dfracPN_colim'
dmu_dQN     = fac_PN_colim * lim_N * gross_nue * self%iK_QN 

a3          = phy%rel_QN 
a2          = lim_N + (dbal_dv -1.0d0)*(lim_N + a3)
a1          = self%zeta_CN*(phy%QN - self%QN_phy_0) + lim_N*phy%frac%NutUpt/(phy%frac%Rub + eps)* &
               self%K_QN_phy * (self%zeta_CN + 1.0d0/phy%QN)

! --- relative growth rate RGR: gross production - exudation - uptake respiration --  

feedb_vq    = 1.0d0 - dbal_dv * a3 / (1.0d0 - self%QN_phy_0/phy%QN + a2 + a3)
dmudVN      = dbal_dv * feedb_vq/((phy%QN - self%QN_phy_0) * (1.0d0/a1 + 1.0d0/a2) + phy%QN )  

! approximated 1/(mu/(dmudQ*Q) +1 - dVdQ)
! --- down-regulation of uptake I: balance of respiration and indirect benefits -------------------  
grad_VN     = dVNC_dregV * (-self%zeta_CN + dmudVN) 
!   numerically smeared version of step function, uses basic
!              parameters to emulate a continuous response 
reg_V_NC    = 1./(1.+exp(-self%tau_regV * grad_VN));  ! 0.02

!------------------------------------------------------------------------
!     adajust all terms for down-regulated uptake      
!-------------------------------------------------------------------------
upt_N       = reg_V_NC * sens%up_NC
uptake%N    = phy%frac%NutUpt * upt_N  ! [(mmolN) (mmolC)^{-1} d{-1}]
! ---  respiration due to N assimilation --------------------------------------
lossC       = self%zeta_CN * uptake%N                           ! [d^{-1}]

! --- relative growth rate RGR: gross production - exudation - uptake respiration --  
uptake%C    = grossC - lossC     !* (1.0d0- self%exud_phy) ![d^{-1}]
!  write (*,'("  fagN",5(F9.4))') par,,up_NC,mu,din

! exud propto uptake%grossC^2 ???
! +++ derivatives of N-uptake rate ++++++++++++++++++++++++++++++++++++++++++
dVNC_dfracR = -upt_N * (1.0d0 + phy%rel_QN**self%sigma * phy%theta * self%itheta_max)                        

dVNC_dQN    = -upt_N * acc%dRchl_dQN * self%itheta_max
dVNC_dtheta = -upt_N * phy%rel_chloropl * self%itheta_max 

! update dmu_dQN
dmu_dQN     = dmu_dQN - self%zeta_CN * dVNC_dQN 
!      if (.not. self%PhotoacclimOn) exit ! break the loop for constant theta and frac_R
! --- derivatives of P-uptake rate ----------------------------------------------------------------          
if (self%PhosphorusOn) then
    grad_VP   = dVPC_dregV * (-2*self%zeta_CN   + lim_P * dmudVP)
    reg_V_PC  = 1./(1.+exp(-self%tau_regV * grad_VP));  ! 0.02
    upt_P     = reg_V_PC * sens%up_PC 
! cell physiology regulating P - uptake is identical to N-uptake machinery
! TODO: explain, why (cf. Smith et al. 2009)
!  
    uptake%P  = phy%frac%NutUpt * upt_P    ! [(mmolP) (mmolC)^{-1} d{-1}]
 ! dVPC_dQ  = V_PC * VPqN * self%iK_QN	with VPqn=1 ! TODO
    dVPC_dQN  = -upt_P * acc%dRchl_dQN * self%itheta_max
    dmu_dQP   = deriv_fac_PN_colim * gross_nue  * self%iK_QP
else
    dmu_dQP   = 0.0d0
end if 

! --- derivatives of C-uptake rate  --------------------------------------------         
dmu_dtheta  = phy%frac%Rub * phy%rel_QN * fac_PN_colim * sens%P_max_T * sens%a_light * (1.0d0 - sens%S_phot) & 
                      - self%zeta_CN * dVNC_dtheta                                                     !OK
dmu_dfracR  = phy%rel_QN * fac_PN_colim * sens%P_max_T * sens%S_phot - self%zeta_CN * dVNC_dfracR         

! +++ trade-off: effects on N uptake +++++++++++++++++++++++++++++++++++++++++++
! dmu is a central regulation parameter determining the coupling strength in ressource uptake; its position in the denominator makes it critical for frequently occuring low values
! TODO: replace with more robust formulation                                                
dmu         = dmu_dQN * phy%QN + uptake%C - dVNC_dQN
! rough estimate of maximal growth rate with frac_R=0.5 and q=1
dmu_small   = eps * sens%P_max_T 
dmu         = smooth_small(dmu , dmu_small); 

dQN_dtheta  = (dVNC_dtheta -1* dmu_dtheta * phy%QN) / dmu
dQN_dfracR  = (dVNC_dfracR -1* dmu_dfracR * phy%QN) / dmu
dQN_dtheta  = (1.0d0  + 0*dQN_dtheta * dVNC_dtheta) * dQN_dtheta
dQN_dfracR  = (1.0d0  + 0*dQN_dfracR * dVNC_dfracR) * dQN_dfracR

! --- photoacclimation and photosynthesis ------------------------------------------------            
!     differential coupling between pigment synthesis and costs due to N-uptake     
! positive gradient term due to PAR adsorption by CHL 
grad_theta  = dmu_dtheta + dmu_dQN * dQN_dtheta

! --- trade-off: effects on P uptake -----------------------------------------------------         
if (self%PhosphorusOn .and. IsVQP .gt. 0.0d0 ) then
   dmuP       = dmu_dQP * phy%QP + uptake%C
   dmuP       = smooth_small(dmuP , dmu_small); 
! rough estimate of maximal growth rate with frac_R=0.5 and q=1

   dVPC_dfracR = -upt_P * (1.0d0 + phy%rel_QN**self%sigma * phy%theta * self%itheta_max)                        
   dVPC_dtheta = -upt_P * phy%rel_chloropl * self%itheta_max 

   dQP_dtheta = (dVPC_dtheta - 0* dmu_dtheta * phy%QP) /dmuP
   dQP_dfracR = (dVPC_dfracR - 0* dmu_dfracR * phy%QP) /dmuP
!   dQP_dtheta = (1.0d0  + 0*dQP_dtheta * dVPC_dtheta) * dQP_dtheta
!   dQP_dfracR = (1.0d0  + 0*dQP_dfracR * dVPC_dfracR) * dQP_dfracR
   grad_theta = grad_theta + dmu_dQP * dQP_dtheta

else
   dQP_dfracR = 0.0d0
end if ! PhosphorusOn

! --- nitrogen partitioning to chloroplast and rubisco --------------------------------------------
if (self%PhotoacclimOn .and. self%RubiscoOn) then
   grad_fracR = dmu_dfracR                   & ! marginal C gain of chloroplasts 
              + dmu_dQN * dQN_dfracR         & ! marginal loss due to reduced uptake 
              + IsVQP * dmu_dQP * dQP_dfracR  ! marginal loss due to reduced uptake 
! rubisco-N directly limits chl synthesis
! --- regulation speed in Rubisco expression -------------------------------------- 
   flex_fracR = self%adap_Rub * (1.0d0 - phy%frac%Rub ) * phy%frac%Rub

          ! *** ADAPTIVE EQUATION FOR 'frac_R'
   acc%dfracR_dt = flex_fracR * grad_fracR  
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
