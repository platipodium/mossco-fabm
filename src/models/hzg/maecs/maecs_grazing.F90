!
#include "fabm_driver.h"
!---------------------------------------------------------
!BOP
! +++ zooplankton grazing +++++++++++++++++++++++++++++++++++++++++++++++++
!
! !MODULE: MAECS_functions --- more to come
!  Model for Adaptive Ecosystems in Coastal Seas 
   module maecs_grazing

! !USES:
   use fabm_types
   use maecs_types
   use maecs_functions
   private
   public   grazing, grazing_losses
 contains  

subroutine grazing(Imax,HalfSat,preyconc,rate)
!
  implicit none
  REALTYPE, intent(in)       :: Imax,HalfSat,preyconc
  REALTYPE, intent(out)      :: rate
!-------------------------------------------------------------
  ! Holling-III response function
  rate   = Imax * preyconc**2 /(HalfSat**2+preyconc**2)          ! [d^{-1}]
  end subroutine grazing


!---------------------------------------------------------
!subroutine grazing_exudates(rate,resC,NC_zoo,NC_prey,PC_zoo,PC_prey,lossC,lossN,lossP,isP)
! --------- loss rates of grazers to inorganic and   organic particulate
subroutine grazing_losses(zoo,resC,NC_prey,PC_prey,lossZNut,lossZDet,isP,isTotIng)
!
  implicit none
!  REALTYPE, intent(in)       :: rate,resC,NC_zoo,NC_prey,PC_prey,PC_zoo
  type (type_maecs_zoo), intent(in)		:: zoo ! zooplankton type containing quota information
  REALTYPE, intent(in)       			:: resC, NC_prey, PC_prey
  logical, intent(in)        			:: isP, isTotIng         ! switches
  type (type_maecs_om), intent(out)       	:: lossZNut,lossZDet 
! -------------- loss by floppy feeding and egestion to detritus
  lossZDet%C     = zoo%flopp  * zoo%feeding 
  lossZDet%N     = lossZDet%C * NC_prey

!----------------------------------------------------------------
! -------------- homeostasis for zooplankton
!     default rates, given that N:C and P:C of prey match 'const_NC_zoo' and 'const_PC_zoo'

! -------------- basal respiration  ; TODO: add activity respiration
  lossZNut%C     = resC  ! [d^{-1}]
! -------------- assuming homeostasis for zooplankton relaxation towards Redfield ratio
!                 nitrogen excretion (urea or ammonia)
   lossZNut%N     = lossZNut%C*zoo%QN + zoo%yield * zoo%feeding * ( NC_prey - zoo%QN ) ! [molN/molC d^{-1}]
! lossZNut%N     = zoo%yield * zoo%feeding * ( NC_prey - zoo%QN ) ! [molN/molC d^{-1}]
!  lossZNut%N  = lossZNut%C *  zoo%QN + excess_N_graz ! includes basal excretion

! --------------- respiration and/or excretion rates adjusted to maintain constant C:N:P ratio
  if (lossZNut%N .lt. 0.0d0) then
! ------ compensate with unused prey components  ?
     if (isTotIng) then  
        lossZDet%N  = lossZDet%N + lossZNut%N
        if (lossZDet%N .lt. 0.0d0 ) then ! compensate by respiration
          lossZNut%C  = lossZNut%C - lossZDet%N/zoo%QN
          lossZDet%N  = 0.0d0
        endif
     else
        lossZNut%C  = lossZNut%C - lossZNut%N/zoo%QN
     endif
     lossZNut%N  = 0.0d0
  end if 

  if (isP) then ! P-relaxation
!    lossZNut%P  = lossZNut%C * zoo%QP  + excess_P_graz ! includes basal excretion
    lossZNut%P   = lossZNut%C*zoo%QP + zoo%yield * zoo%feeding * ( PC_prey - zoo%QP ) ![mmolP m^{-3} d^{-1}] 
!    lossZNut%P    = zoo%yield * zoo%feeding * ( PC_prey - zoo%QP ) ![mmolP m^{-3} d^{-1}] 
! -------------- loss by floppy feeding to detritus
    lossZDet%P   = zoo%flopp * zoo%feeding* PC_prey 
!    lossZDet%P    = lossZDet%C * PC_prey
    if (lossZNut%P .lt. 0.0d0) then
      if (isTotIng) then  ! compensate with unused prey components
         lossZDet%P  =  lossZDet%P + lossZNut%P 
         if (lossZDet%P .lt. 0.0d0 ) then ! compensate by respiration & exudation
            lossZNut%C  = lossZNut%C - lossZDet%P/zoo%QP 
            lossZNut%N  = lossZNut%N - lossZDet%P/zoo%QP * zoo%QN
            lossZDet%P  = 0.0d0
         endif
      else
         lossZNut%C   = lossZNut%C - lossZNut%P/zoo%QP 
         lossZNut%N   = lossZNut%N - lossZNut%P/zoo%QP * zoo%QN
      endif ! isTotIng
      lossZNut%P   = 0.0d0
    end if ! lossZNut%P
  end if ! isP
  
  end subroutine grazing_losses
end module maecs_grazing
