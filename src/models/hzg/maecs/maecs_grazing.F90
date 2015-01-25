!> @file maecs_grazing.F90
!> @author Richard Hofmeister, Markus Schartau, Kai Wirtz

#include "fabm_driver.h"

!> @brief Grazing module
module maecs_grazing

! !USES:
   use fabm_types
   use maecs_types
   use maecs_functions
   private
   public   grazing, grazing_losses
 contains  

!> @brief  calculates grazing rate
!> @details 
!> rate= @f$ I_{max} * F^2/(K^2+F^2) @f$
!> (Holling-III response function)
subroutine grazing(Imax,HalfSat,preyconc,rate)
!
  implicit none
  real(rk), intent(in)       :: Imax,HalfSat,preyconc
  real(rk), intent(out)      :: rate

  rate   = Imax * preyconc**2 /(HalfSat**2+preyconc**2)          ! [d^{-1}]
  end subroutine grazing


!---------------------------------------------------------
!grazing_exudates(rate,resC,NC_zoo,Q_prey%N,PC_zoo,Q_prey%P,lossC,lossN,lossP,isP)
!> @brief  loss rates of grazers to inorganic and  organic particulate
!> @details 
!> assumes a constant C:N:P unit feeding on variable C:N:P food 
subroutine grazing_losses(zoo,resC,Q_prey,lossZNut,lossZDet,mswitch)
!called as:grazing_losses(zoo,zoo_respC,nquot,lossZ,floppZ, mswitch) 
!
  implicit none
  type (type_maecs_zoo), intent(in)     :: zoo ! zooplankton type containing quota information
  real(rk), intent(in)                  :: resC
  type (type_maecs_switch), intent(in)  :: mswitch
  type (type_maecs_om), intent(in)      :: Q_prey
  type (type_maecs_om), intent(out)     :: lossZNut,lossZDet 
  
!>@fn maecs_grazing::grazing_losses()
!> 1. calculate lossZDet (floppZ in maecs_do())
!>    - lossZDet\%C=zoo\%floppI*zoo\%feeding
!>    - lossZDet\%X=zoo\%floppI*zoo\%feeding*Q_prey\%X
!> 2. calculate lossZNut (lossZ in maecs_do())
!>    - lossZNut\%C=resC (= self\%basal\_resp\_zoo * sens\%f\_T as calc in maecs_do() )
!>    - lossZNut\%X=resC*zoo\%Q\%X + zoo\%yield * zoo\%feeding * ( Q_prey\%X - zoo%Q%X )
!>      + idea is, zoo takes in Y*F*QX_{food}, keeps Y*F*QX_{zoo}, excretes the rest (+ bg)
!>      + zoo\%Q\%X calc. in maecs_functions::calc_internal_states()
!>      + zoo\%feeding=grazing() as set in maecs_do()

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

!>@fn maecs_grazing::grazing_losses()
!> 3. adjust lossZDet \&lossZNut to maintain constant C:N:P ratio
!>   + if lossZNut\%N<0 -> negative excretion! (would gain N)
!>     - lossZDet\%N = lossZDet\%N + lossZNut\%N;  lossZNut\%N=0 
!>     - if lossZDet\%N <0 -> negative loss. Then increase the C excr
!>       * lossZDet\%N=0, but still the N balance goes up -> compensate by increased lossZNut\%C
!>       * lossZNut\%C  = lossZNut\%C - (lossZNut\%N-lossZDet\%N(prev))
!>       * subst. lossZDet\%N(prev)=lossZDet\%N-lossZNut\%N
!>       * lossZNut\%C  = lossZNut\%C - (lossZNut\%N-(lossZDet\%N-lossZNut\%N))/zoo\%Q\%N
!>       * lossZNut\%C  = lossZNut\%C - lossZDet\%N/zoo\%Q\%N
!>   + if lossZNut\%P<0 -> negative excretion! (would gain P)
!>     - the same way, adjust lossZNut\%C(\&N) and lossZDet\%C(\&N) to maintain the const. C:N:P 
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
