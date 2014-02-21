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
  real(rk), intent(in)       :: Imax,HalfSat,preyconc
  real(rk), intent(out)      :: rate
!-------------------------------------------------------------
  ! Holling-III response function
  rate   = Imax * preyconc**2 /(HalfSat**2+preyconc**2)          ! [d^{-1}]
  end subroutine grazing


!---------------------------------------------------------
!subroutine grazing_exudates(rate,resC,NC_zoo,Q_prey%N,PC_zoo,Q_prey%P,lossC,lossN,lossP,isP)
! --------- loss rates of grazers to inorganic and   organic particulate
subroutine grazing_losses(zoo,resC,Q_prey,lossZNut,lossZDet,mswitch)
!
  implicit none
  type (type_maecs_zoo), intent(in)     :: zoo ! zooplankton type containing quota information
  real(rk), intent(in)                  :: resC
  type (type_maecs_switch), intent(in)  :: mswitch
  type (type_maecs_om), intent(in)      :: Q_prey
  type (type_maecs_om), intent(out)     :: lossZNut,lossZDet 
! -------------- loss by floppy feeding and egestion to detritus
  lossZDet%C     = zoo%flopp  * zoo%feeding 
  lossZDet%N     = lossZDet%C * Q_prey%N

!----------------------------------------------------------------
! -------------- homeostasis for zooplankton
!     default rates, given that N:C and P:C of prey match 'const_NC_zoo' and 'const_PC_zoo'

! -------------- basal respiration  ; TODO: add activity respiration
  lossZNut%C     = resC  ! [d^{-1}]
! -------------- assuming homeostasis for zooplankton relaxation towards Redfield ratio
!                 nitrogen excretion (urea or ammonia)
   lossZNut%N     = lossZNut%C*zoo%Q%N + zoo%yield * zoo%feeding * ( Q_prey%N - zoo%Q%N ) ! [molN/molC d^{-1}]
! lossZNut%N     = zoo%yield * zoo%feeding * ( Q_prey%N - zoo%Q%N ) ! [molN/molC d^{-1}]
!  lossZNut%N  = lossZNut%C *  zoo%Q%N + excess_N_graz ! includes basal excretion

! --------------- respiration and/or excretion rates adjusted to maintain constant C:N:P ratio
  if (lossZNut%N .lt. 0.0d0) then
! ------ compensate with unused prey components  ?
     if (mswitch%isTotIng) then  
        lossZDet%N  = lossZDet%N + lossZNut%N
        if (lossZDet%N .lt. 0.0d0 ) then ! compensate by respiration
          lossZNut%C  = lossZNut%C - lossZDet%N/zoo%Q%N
          lossZDet%N  = 0.0d0
        endif
     else
        lossZNut%C  = lossZNut%C - lossZNut%N/zoo%Q%N
     endif
     lossZNut%N  = 0.0d0
  end if 

  if (mswitch%isP) then ! P-relaxation
!    lossZNut%P  = lossZNut%C * zoo%Q%P  + excess_P_graz ! includes basal excretion
    lossZNut%P   = lossZNut%C*zoo%Q%P + zoo%yield * zoo%feeding * ( Q_prey%P - zoo%Q%P ) ![mmolP m^{-3} d^{-1}] 
!    lossZNut%P    = zoo%yield * zoo%feeding * ( Q_prey%P - zoo%Q%P ) ![mmolP m^{-3} d^{-1}] 
! -------------- loss by floppy feeding to detritus
    lossZDet%P   = zoo%flopp * zoo%feeding* Q_prey%P 
!    lossZDet%P    = lossZDet%C * Q_prey%P
    if (lossZNut%P .lt. 0.0d0) then
      if (mswitch%isTotIng) then  ! compensate with unused prey components
         lossZDet%P  =  lossZDet%P + lossZNut%P 
         if (lossZDet%P .lt. 0.0d0 ) then ! compensate by respiration & exudation
            lossZNut%C  = lossZNut%C - lossZDet%P/zoo%Q%P 
            lossZNut%N  = lossZNut%N - lossZDet%P/zoo%Q%P * zoo%Q%N
            lossZDet%P  = 0.0d0
         endif
      else
         lossZNut%C   = lossZNut%C - lossZNut%P/zoo%Q%P 
         lossZNut%N   = lossZNut%N - lossZNut%P/zoo%Q%P * zoo%Q%N
      endif ! mswitch%isTotIng
      lossZNut%P   = 0.0d0
    end if ! lossZNut%P
  end if ! mswitch%isP

  if (mswitch%isSi) then ! Si-release
! neglecting silicification in zooplankton such as in choanoflagellates, radiolarians
    lossZNut%Si   = zoo%yield * zoo%feeding * Q_prey%Si !(mmolSi m^{-3} d^{-1}) 
! -------------- loss by floppy feeding to detritus
    lossZDet%Si   = zoo%flopp * zoo%feeding* Q_prey%Si 
  end if ! mswitch%isSi
  
  end subroutine grazing_losses
end module maecs_grazing
