!!-------------------------------------------------------------------
type type_maecs_om
   real(rk) :: C,N,P,Si  ! all elements resolved in MAECS are listed here; only here! 
end type

type type_maecs_basic_traits ! traits that apply to all -planktonic- life forms
   real(rk) :: lESD, vESD
   real(rk) :: trophy
end type

type type_maecs_allocation_fractions
   real(rk) :: Rub, theta
   real(rk) :: rel_phys ! physiological(=energetical&nutritional) status 0...1
end type

type,extends(type_maecs_om) type_maecs_life
      type (type_maecs_om)            :: quota       ! former QN,QP,QPN
      type (type_maecs_basic_traits)  :: btrait
end type

type,extends(type_maecs_life) type_maecs_phy
      real(rk)   :: Rub, theta       ! redundancy, cf. "Rub, theta" above
      real(rk)   :: chl, rel_chloropl
      type (type_maecs_om) :: rel_quota   ! former relQ%N, relQ%P, relQ%NP
      type (type_maecs_om) :: reg         ! former C_reg, N_reg, P_reg
      type (type_maecs_allocation_fractions) :: frac
end type

type,extends(type_maecs_life) type_maecs_zoo
      real(rk)   :: yield,flopp
      real(rk)   :: feeding
end type

! Time and Chl derivatives of trait variables
type type_maecs_traitdyn
      real(rk)   :: dtheta_dt
      real(rk)   :: dfracR_dt
      real(rk)   :: dRchl_dtheta
      real(rk)   :: dRchl_dfracR 
      real(rk)   :: dRchl_dQN
      real(rk)   :: tmp,fac1,fac2  ! for volatile diagnostics
end type


! TODO generalize matter structure C,N,P,QN,QP

type type_maecs_derivative
   real(rk) :: dregV
   real(rk) :: dtheta
   real(rk) :: dfracR
   real(rk) :: dfracP
   real(rk) :: dQN
end type

