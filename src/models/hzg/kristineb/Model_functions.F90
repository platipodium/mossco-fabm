subroutine ecophys_para(self,parsout)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
   real(rk),dimension(29),intent(out) :: parsout
   real(rk) :: mumax_incr2,b_Mumax_large2,b_Mumax_small2,a_Mumax_large2,a_Mumax_small2
   real(rk) :: b_Qmin_N2,a_Qmin_N2,b_Qmax_N2,a_Qmax_N2,b_Vmax_N2,a_Vmax_N2
   real(rk) :: a_affin_N2,a_affin_P2,a_carbon2
   real(rk) :: b_affin_N2,b_affin_P2,b_carbon2
   real(rk) :: b_Qmin_P2,a_Qmin_P2,b_Qmax_P2,a_Qmax_P2,b_Vmax_P2,a_Vmax_P2
!----------------
!--- Phytoplankton eco-physiological parameters ---
  b_qmin_N2=10.0**self%b_qmin_N
  b_qmax_N2=10.0**self%b_qmax_N
  b_vmax_N2=10.0**self%b_vmax_N
  !b_kn_N2=10.0**self%b_kn_N
  b_carbon2=10.0**self%b_carbon
  b_qmax_P2=10.0**self%b_qmax_P
  b_qmin_P2=10.0**self%b_qmin_P
  b_vmax_P2=10.0**self%b_vmax_P
 ! b_kn_P2=10.0**self%b_kn_P
  b_affin_N2=10.0**self%b_affin_N
  b_affin_P2=10.0**self%b_affin_P
  !      Nonlinear mumax
  b_mumax_large2=10.0**self%b_mumax_large
  b_Mumax_large2=b_mumax_large2*(acos(-1.0)/6.)**self%a_mumax_large
  a_Mumax_large2=3*self%a_mumax_large 
  b_Mumax_small2=10.0**self%b_mumax_small 
  b_Mumax_small2=b_mumax_small2*(acos(-1.0)/6.)**self%a_mumax_small
  a_Mumax_small2=3*self%a_mumax_small 
   !     
        ! Conversion of Phytoplankton eco-physiological parameters to ESD and mole-C base --- Nitorgen
        call convert_BGCparams(b_qmin_N2,self%a_qmin_N,self%a_carbon,b_carbon2,b_qmin_N2,a_qmin_N2)
        call convert_BGCparams(b_qmax_N2,self%a_qmax_N,self%a_carbon,b_carbon2,b_qmax_N2,a_qmax_N2)
        call convert_BGCparams(b_vmax_N2,self%a_vmax_N,self%a_carbon,b_carbon2,b_vmax_N2,a_vmax_N2)


      !  b_Kn_N2=b_kn_N2*(acos(-1.0)/6.)**self%a_kn_N
      !  a_Kn_N2=3*(self%a_kn_N)

        ! Nutrient affinity, m^3 mmol-C d^1
        !a_affin_N2= a_Vmax_N2/a_Kn_N2 !-1
        !b_affin_N2= b_Vmax_N2/b_Kn_N2 !0.4

        ! Conversion of Phytoplankton eco-physiological parameters to ESD and mole-C base --- Phosphorous
        call convert_BGCparams(b_qmin_P2,self%a_qmin_P,self%a_carbon,b_carbon2,b_qmin_P2,a_qmin_P2)
        call convert_BGCparams(b_qmax_P2,self%a_qmax_P,self%a_carbon,b_carbon2,b_qmax_P2,a_qmax_P2)
	call convert_BGCparams(b_vmax_P2,self%a_vmax_P,self%a_carbon,b_carbon2,b_vmax_P2,a_vmax_P2)
	a_Qmax_P2=0.0_rk
       ! b_Kn_P2=b_kn_P2*(acos(-1.0)/6.)**self%a_kn_P
       ! a_Kn_P2=3*(self%a_kn_P)

        call convert_BGCparams(b_affin_N2,self%a_affin_N,self%a_carbon,b_carbon2,b_affin_N2,a_affin_N2)
        call convert_BGCparams(b_affin_P2,self%a_affin_P,self%a_carbon,b_carbon2,b_affin_P2,a_affin_P2)
 parsout=0.0_rk
 parsout(1) = b_Mumax_small2
 parsout(2) = a_Mumax_small2
 parsout(3) = self%mumax_incr
 parsout(4) = b_Mumax_large2
 parsout(5) = a_Mumax_large2
 parsout(6) = b_Qmin_N2
 parsout(7) = a_Qmin_N2
 parsout(8) = b_Qmax_N2
 parsout(9) = a_Qmax_N2
 parsout(10) = b_Vmax_N2
 parsout(11) = a_Vmax_N2
 parsout(12) = b_affin_N2
 parsout(13) = a_affin_N2
 !parsout(14) = b_Kn_N2
 !parsout(15) = a_Kn_N2
 parsout(16) = b_Qmin_P2
 parsout(17) = a_Qmin_P2
 parsout(18) = b_Qmax_P2
 parsout(19) = a_Qmax_P2
 parsout(20) = b_Vmax_P2
 parsout(21) = a_Vmax_P2
 !parsout(22) = b_Kn_P2
 !parsout(23) = a_Kn_P2
 parsout(24) = b_affin_P2
 parsout(25) = a_affin_P2
 !parsout(26) = b_mumax2
 !parsout(27) = a_mumax2
 parsout(28) = b_carbon2
 parsout(29) = a_carbon2
 return
 end subroutine
!--------------------------------------------------------------------
subroutine f_T(self,deeptemp,f_T_out)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),dimension(2),intent(out):: f_T_out
 real(rk)::ft_Phy,ft_zoo
 real(rk),intent(in) :: deeptemp ! Interpolation functions for deep water temperature
! real(rk),intent(in) :: seatemp ! Interpolation functions for see surface temperature
 real(rk) :: meantemp
!        :return: Temperature dependency
        if  (self%T_forc .eqv. .true.) then
             meantemp=(deeptemp)!+seatemp)/2.0_rk
            fT_Phy= self%Tcons_phy**((meantemp-self%T_ref)/10._rk)
            fT_zoo= self%Tcons_zoo**((meantemp-self%T_ref)/10._rk)
            f_T_out=(/fT_Phy,fT_zoo/)
        else
            fT_Phy=1.0
            fT_zoo=1.0
            f_T_out=(/fT_Phy,fT_zoo/)
	end if
	return
end subroutine
!-----------------------------------------------------------------------------
subroutine  F_Co2sr(self,pCO2,f_co2)!f_LpCo2,f_HpCo2,log_ESD)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),intent(in) :: pCO2
integer::i
real(rk),dimension(self%phyto_num),intent(out)::f_co2
real(rk)::nom,dom
f_co2(:)=1.0_rk
       if (self%co2_forc .eqv. .true.) then
            do i=1,self%phyto_num
                nom=1.0_rk-exp(-self%a_co2*pco2)
                dom=1.0_rk+self%a_star*exp(self%log_ESD(i)-self%a_co2*pco2)
            	f_co2(i)=(nom/dom)
	    end do
        end if
        return
end subroutine

!-----------------------------------------------------------------------------
subroutine f_parsr(self, Phy, Q_N,f_co2,F_T,par,mixl,f_par)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 !real(rk),intent(in):: t! Time
 real(rk),dimension(self%phyto_num),intent(in):: Phy! Phytoplankton biomass concentration, mmol-C m^-3
 real(rk),dimension(self%phyto_num),intent(in):: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
 real(rk),dimension(self%phyto_num),intent(in):: f_co2! CO2 forcing
 real(rk),dimension(2),intent(in):: F_T! Temperature dependency for phytoplankton
 real(rk),dimension(self%phyto_num),intent(out)::f_par! PAR forcing
 real(rk),dimension(11) :: bgc_params
 real(rk),intent (in)::par
 integer::i
 real(rk)::k,par_w,phyNtot,par_tmp
 real(rk),dimension(self%phyto_num)::phyN
 real(rk),intent(in) :: mixl
 
 par_tmp=par!*3600*24
 f_par(:)=1._rk
        if (self%PAR_forc .eqv. .true.) then
!             phyN: phytoplankton concentration, mmol-N m^-3
	    do i=1,self%phyto_num
	            phyN(i) = Q_N(i)*Phy(i)
	    end do
!             phyNtot: Light attenuation due to phytoplankton biomass, m^2 mmol-N^-1
            phyNtot = sum(phyN(:)) * self%k_phyN
            k = self%kbg + phyNtot
!             par_w: Average light intensity within mixed layer depth,
            par_w = par / (mixl * k) * (1.0_rk - exp(-1.0_rk * k * mixl))
          !  par_w = par_tmp / (self%z * k) * (1._rk - exp(-1.0_rk * k * self%z))
            do i=1,self%phyto_num
		call bgc_parameters(self,self%log_ESD(i), bgc_params)
                f_par(i)=1.0_rk-exp(-(self%a_par*par_w*Q_N(i))/(bgc_params(1)*f_co2(i)*F_T(1)))
            !    f_par(i)=par_w/(par_w+((bgc_params(1)*f_co2(i)*F_T(1))/(self%alpha)))
	    end do
	end if
        return 
end subroutine
!-----------------------------------------------------------------------------
subroutine Phy_growth_rate(self, Q_N,Q_P,F_T,F_co2,F_par, P_growth_rate)
!use model_pars
implicit none
 class (type_hzg_kristineb),intent(in) :: self
   real(rk),dimension(self%phyto_num),intent(in) :: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
   real(rk),dimension(self%phyto_num),intent(in) :: Q_P! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
   real(rk),dimension(self%phyto_num),intent(in) :: F_T! Temperature dependency for phytoplankton
   real(rk),dimension(self%phyto_num),intent(in) :: F_co2! Co2 forcing
   real(rk),dimension(self%phyto_num),intent(in) :: F_par! PAR forcing
   real(rk),dimension(self%phyto_num),intent(out) :: P_growth_rate !Phytoplankton growth rate, d^-1
   integer::i,nutlim
   real(rk),dimension(11) :: bgc_params
   real(rk)::f_nut,r,n,g_N,mu_max,q_Nl,q_Pl
	nutlim=int(self%Nut_lim)
        Do i=1,self%phyto_num
	    call bgc_parameters(self,self%log_ESD(i),bgc_params)
            !if mu_inf is given instead of mu_max, we should correct:
            if (self%convert_mu .eqv. .true.) then
                mu_max=bgc_params(1)*(bgc_params(3)/(bgc_params(3)-bgc_params(2)))
                mu_max=bgc_params(1)*(bgc_params(8)/(bgc_params(8)-bgc_params(7)))
            else 
                mu_max=bgc_params(1)
	    end if
            !to use with mu_max:
            !q_Nl=(Q_N(i)-bgc_params(2))/(bgc_params(3)-bgc_params(2))
            !q_Pl=(Q_P(i)-bgc_params(7))/(bgc_params(8)-bgc_params(7))
            q_Nl=(Q_N(i)-bgc_params(2))/Q_N(i)
            q_Pl=(Q_P(i)-bgc_params(7))/Q_P(i)

            if (nutlim == 1) then
                r=q_Pl/q_Nl
                n=self%n_star*(1._rk+q_Nl)
                g_N=(r-r**(1._rk+n))/(1._rk-r**(1._rk+n))
                f_nut=q_Nl*g_N
            elseif (nutlim == 2) then
                f_nut= min(q_Nl,q_Pl)
            elseif (nutlim == 3) then
                f_nut= q_Nl*q_Pl/(q_Nl+q_Pl)
            elseif (nutlim == 4) then
                f_nut= q_Nl*q_Pl
    	    end if
            P_growth_rate(i) = ( mu_max * f_nut * F_T(1) * F_co2(i) * F_par(i))
	end do
        return
end subroutine
!------------------------------------------------------------------------------
real(rk) function aggr_rate(self, Phy,Q_N,D_N)! Aggregation rate, d^-1
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
 real(rk),dimension(self%phyto_num),intent(in) :: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
 real(rk),intent(in) :: D_N! Nitrogen content of detritus concentration,  mmol-N m^-3
 integer::i
 real(rk),dimension(self%phyto_num) :: phy_conc
  phy_conc=0.0_rk
	Do i=1,self%phyto_num
        	phy_conc(i)=Phy(i)*Q_N(i)
	end do
        !aggr_rate=self%A_star_opt*(product(phy_conc(:))+D_N)
        aggr_rate=self%A_star_opt*(sum(phy_conc(:))+D_N)
        return
end function
!------------------------------------------------------------------------------
subroutine N_uptake(self,N,Q_N,F_T,par,uptake_rate_N)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),intent(in) ::  N! Nitrogen concentration, mmol-N m^-3
 real(rk),dimension(self%phyto_num),intent(in) :: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
 real(rk),dimension(self%phyto_num),intent(in) :: F_T! Temperature dependency
 real(rk),dimension(self%phyto_num),intent(out) :: uptake_rate_N! Phytoplankton nitrogen uptake rate, mol-N mol-C^-1 d^-1
 real(rk),intent(in) :: par !PAR from data
 integer:: i
 real(rk),dimension(11) :: bgc_params
 real(rk)::nom_N,dom_N,q
  uptake_rate_N=0.0_rk
	if (par>=0.0) then  !No uptake during night
        Do i=1,self%phyto_num
	    call bgc_parameters(self,self%log_ESD(i), bgc_params)
            nom_N=bgc_params(4)*bgc_params(6)*N
            dom_N=bgc_params(4)+bgc_params(6)*N
            q=max(0.0_rk,(bgc_params(3)-Q_N(i))/(bgc_params(3)-bgc_params(2)))
            uptake_rate_N(i)=(nom_N/dom_N)*sqrt(F_T(1))*q
	end do
	end if
        return 
end subroutine
!------------------------------------------------------------------------------
subroutine P_uptake(self,P, Q_P,F_T,par,uptake_rate_P)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),intent(in) ::  P! Phsphorous concentration, mmol-P m^-3
 real(rk),dimension(self%phyto_num),intent(in) :: Q_P! Phytoplankton intracellular phosphorous cell quota, mol-P mol-C^-1
 real(rk),dimension(self%phyto_num),intent(in) :: F_T! Temperature dependency for phytoplankton
 real(rk),dimension(self%phyto_num),intent(out) :: uptake_rate_P! Phytoplankton nutrient uptake rate, mol-P mol-C^-1 d^-1
 real(rk),intent(in) :: par !PAR from data
 integer:: i
 real(rk) :: nom_P,dom_P,q
 real(rk),dimension(11) :: bgc_params
 	uptake_rate_P=0.0_rk
	if (par>=0.0) then !No uptake during night
	  Do i=1,self%phyto_num
	    call bgc_parameters(self,self%log_ESD(i), bgc_params)
            nom_P=bgc_params(9)*bgc_params(11)*P
            dom_P=bgc_params(9)+bgc_params(11)*P
        !    q=max(0.0_rk,(bgc_params(8)-Q_P(i))/(bgc_params(8)-bgc_params(7)))
            uptake_rate_P(i)=(nom_P/dom_P)*sqrt(F_T(1))
            !uptake_rate_P(i)=(nom_P/dom_P)*sqrt(F_T(1))*q
	end do
	end if
        return
end subroutine
!------------------------------------------------------------------------------
subroutine sink_rate(self,Q_N,Q_P,mixl,sinking)
!use model_pars
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),dimension(self%phyto_num),intent(in) :: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
 real(rk),dimension(self%phyto_num),intent(in) :: Q_P! Phytoplankton intracellular phosphorous cell quota, mol-P mol-C^-1
 real(rk),dimension(self%phyto_num),intent(out) :: sinking! sinking rate, d^-1
 integer:: i
 real(rk) :: physiol,qN,qP
 real(rk),dimension(11) :: bgc_params
 real(rk),intent(in)::mixl !mixing layer depth
 sinking=0.0_rk
        Do i=1,self%phyto_num
	    call bgc_parameters(self,self%log_ESD(i), bgc_params)
            qN=(Q_N(i)-bgc_params(2))/(bgc_params(3)-bgc_params(2))
            qP=(Q_P(i)-bgc_params(7))/(bgc_params(8)-bgc_params(7))
!         healthy, non-limited cells create buouyancy
 !           physiol = exp(-4.0_rk*qN*qP)
            physiol = exp(-0.5*((qP*qN)*16.0)**2) !Todo: Why?
!         size dependency: Stokes - vacuolation
            sinking(i)=physiol*exp(0.5_rk*self%log_ESD(i))* 0.06_rk/mixl
       !     sinking(i)=physiol*exp(0.25_rk*self%log_ESD(i))* 0.3_rk/mixl      
!z=mixed layer depth
	end do
        return
end subroutine
!------------------------------------------------------------------------------
subroutine Respiration(self, N_uptake,R_N)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),dimension(self%phyto_num),intent(in) :: N_uptake! Nitrogen uptake rate, mol-N mol-C^-1 d^-1
 real(rk),dimension(self%phyto_num),intent(out) :: R_N!Phytoplankton respiration rate, d^-1
 integer:: i
        do i=1,self%phyto_num
            R_N(i)=N_uptake(i)*self%mol_ratio
	end do
        return
end subroutine
!------------------------------------------------------------------------------
subroutine Grazing_forcing(self,Phy,F_T,Mean,zoo_pref,cop_pref,Zoo,grazing,Lz,I_max)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 integer:: i,j,c
 real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
 real(rk),dimension(self%zoo_num),intent(in) :: Zoo,Lz !Zooplankton biomass
 real(rk),dimension(2),intent(in) :: F_T!(2) Temperature dependency for zooplankton
 real(rk),intent(in) :: Mean! Community mean cell size, log_e ESD (mu m)
 real(rk),dimension(self%phyto_num),intent(out) :: grazing !Grazing forcing, mmol-C m^-3 d^-1
 real(rk)::eff_food_ciliat,nom,dom,a_zoo,I_max_star, x,g_x,G_min,mean_eff_food_ciliat,G,downreg
 real(rk),dimension(self%zoo_num,self%phyto_num):: eff_food,mean_eff_food,graz_j,graz_j_spec
 real(rk),dimension(self%zoo_num) :: eff_food_con,mean_eff_food_con,glob_graz
 real(rk),dimension(self%zoo_num), intent(out) :: I_max
 real(rk),dimension(self%zoo_num,self%phyto_num),intent(in)::zoo_pref
 real(rk),dimension(self%num_ciliat), intent(in)::cop_pref
 real(rk),parameter::eps=0.001_rk
 grazing=0.0_rk
	if (self%graz_forc .eqv. .true.) then
            ! Effective food concentration
	    eff_food_ciliat=0.0_rk
            do c=1,self%num_ciliat
	      eff_food_ciliat=eff_food_ciliat+(Zoo(c)*cop_pref(c))
	    end do
	    eff_food=0.0_rk
	    do j=1,self%zoo_num
	      do i=1,self%phyto_num
                eff_food(j,i)=Phy(i)*zoo_pref(j,i)
              end do
	      eff_food_con(j)=sum(eff_food(j,:))
	    end do
            !effective food con for copepoda-> Copepoda graze on phytoplankton + ciliates
            eff_food_con(self%zoo_num)=eff_food_con(self%zoo_num)+eff_food_ciliat

            ! Average effective food concentration
            mean_eff_food_ciliat=0.0_rk
            Do c=1,self%num_ciliat
	      mean_eff_food_ciliat=mean_eff_food_ciliat+Zoo(c)*cop_pref(c) * Lz(c)
	    end do
            mean_eff_food=0.0_rk
            do j=1,self%zoo_num
	      do i=1,self%phyto_num
                mean_eff_food(j,i)=Phy(i)*zoo_pref(j,i)*self%log_ESD(i)
              end do
            mean_eff_food_con(j)=sum(mean_eff_food(j,:))
            end do
            mean_eff_food_con(self%zoo_num)=mean_eff_food_con(self%zoo_num)+mean_eff_food_ciliat

            ! Maximum ingestion rate
            I_max=0.0_rk
            Do j=1,self%zoo_num
                a_zoo=self%a_Im0*(Lz(j)+self%Lz_star(j))
                I_max_star=self%I_max0*F_T(2)*exp(a_zoo+(2.0_rk-a_zoo)*self%Lz_star(j)+(a_zoo-3.0_rk)*Lz(j))
                I_max(j)=I_max_star*exp(-self%sel(j)*(self%Lz_star(j)-mean_eff_food_con(j)/eff_food_con(j))**2)
	    end do
            ! Calculation of grazing, x: food processing ratio, g_x: functional response
            glob_graz=0.0_rk
            do j=1,self%zoo_num
                x=self%a_gr * eff_food_con(j)/I_max(j)
                g_x=1.0_rk-exp(-x) !((1.0_rk-x**self%n_syn)/(1.0_rk-x**(self%n_syn+1.0_rk)))*x
                glob_graz(j)=self%graz_const*Zoo(j)*I_max(j)*g_x
	    end do
            !threshold
            G_min=self%a_gr*self%R_A/self%y

            graz_j=0.0_rk
            Do i=1,self%phyto_num
	      Do j=1,self%zoo_num
               ! graz_j_spec(j,i)=zoo_pref(j,i)*Phy(i)/(eff_food_con(j)+eps)
                graz_j(j,i)=glob_graz(j)*zoo_pref(j,i)*Phy(i)/(eff_food_con(j)+eps)
	      end do
	    end do
	    
	    Do j=1,self%zoo_num
	      G=sum(graz_j(j,:))
	      downreg=1.0_rk/(1.0_rk+exp(-(G-G_min)/0.05_rk))
	      graz_j(j,:)=downreg*graz_j(j,:)
	    end do
	    
	    do i=1,self%phyto_num
		grazing(i)=sum(graz_j(:,i))
	    end do
        else
            grazing(:)=0.0_rk
	end if
        return
end subroutine
!------------------------------------------------------------------------------
real(rk) function dD_N_dt(self, Phy, Q_N, D_N,aggregation,F_T,mixl,grazing_forc) 
implicit none
!Detritus concentration over time, mmol-N m^-3 d^-1
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
 real(rk),dimension(self%phyto_num),intent(in) :: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
 real(rk),intent(in) ::  D_N! Nitrogen content of detritus concentration,  mmol-N m^-3
 real(rk),intent(in) ::  aggregation! aggregation rate, d^-1
 real(rk),dimension(2),intent(in) :: F_T! Temperature dependency
 real(rk),dimension(self%phyto_num),intent(in) :: grazing_forc! Grazing forcing, mmol-C m^-3 d^-1
 integer::i
 real(rk) :: mtotal
 real(rk),dimension(self%phyto_num) :: source,gr,loss
  real(rk),intent(in) :: mixl
        mtotal = 0.0_rk
        do i=1,self%phyto_num
            !source(i)=(self%frac_md*self%m+aggregation) * Phy(i) * Q_N(i)+(1._rk-self%y)*grazing_forc(i)*Q_N(i)
            !mtotal = mtotal+source(i)
	    gr(i)=(1.0_rk-self%y)*grazing_forc(i)*Q_N(i)
	    loss(i)=(self%frac_md*self%m+aggregation)*Phy(i)*Q_N(i)
	end do
        !dD_N_dt=mtotal - self%r_dn * D_N*F_T(1)-(self%det_sink_r/self%z)*D_N
        dD_N_dt=sum(loss(:)+gr(:))-(self%r_dn*F_T(1)+(self%det_sink_r/mixl))*D_N
        return
end function
!------------------------------------------------------------------------------
real(rk) function dD_P_dt(self,Phy, Q_P, D_P,aggregation,F_T,mixl,grazing_forc) !Detritus concentration over time, mmol-P m^-3 d^-1
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
 real(rk),dimension(self%phyto_num),intent(in) :: Q_P! Phytoplankton intracellular phosphorous cell quota, mol-P mol-C^-1
 real(rk),intent(in) :: D_P! Phosphorous content of detritus concentration,  mmol-P m^-3
 real(rk),intent(in) :: aggregation! aggregation rate, d^-1
 real(rk),dimension(2),intent(in) :: F_T! Temperature dependency
 real(rk),dimension(self%phyto_num),intent(in) :: grazing_forc! Grazing forcing, mmol-C m^-3 d^-1
 integer::i
 real(rk) :: mtotal
 real(rk),dimension(self%phyto_num) :: source,gr,loss
  real(rk),intent(in) :: mixl !mixed layer depth
        mtotal = 0.0_rk
        do i=1,self%phyto_num
            !source(i)=(self%frac_md*self%m+aggregation) * Phy(i) * Q_P(i)+(1._rk-self%y)*grazing_forc(i)*Q_P(i)
            !mtotal = mtotal+source(i)
            gr(i)=(1.0_rk-self%y)*grazing_forc(i)*Q_P(i)
            loss(i)=(self%frac_md*self%m+aggregation)*Phy(i)*Q_P(i)
	end do
!        dD_P_dt=mtotal - self%r_dn * D_P*F_T(1) -(self%det_sink_r/self%z)*D_P
	dD_P_dt=sum(loss(:)+gr(:))-(self%r_dn*F_T(1)+(self%det_sink_r/mixl))*D_P
 return
end function
!------------------------------------------------------------------------------
real(rk) function dN_dt(self,N_uptake, Phy, D_N,grazing_forc,Q_N,F_T,N) ! Change of nitrogen concentration over time, mmol-N m^-3 d^-1
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),dimension(self%phyto_num),intent(in) :: N_uptake! Phytoplankton nitorgen uptake rate, mol-N mol-C^-1 d^-1
 real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
 real(rk),intent(in) ::  D_N! Nitrogen content of detritus concentration,  mmol-N m^-3
 real(rk),dimension(self%phyto_num),intent(in) :: grazing_forc! Grazing forcing, mmol-C m^-3 d^-1
 real(rk),dimension(self%phyto_num),intent(in) :: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^1
 real(rk),dimension(2),intent(in) :: F_T! Temperature dependency
 real(rk),intent(in) ::  N! Nitrogen concentration, mmol-N m^-3
! real(rk),intent(in) ::  r_mix! Mixing rate, d^-1
!        :param Nit_bot: Bottom layer nitrogen concentration, mmol-N m^-3
 integer::i
 real(rk),dimension(self%phyto_num) :: dyn_part,up,gr
 real(rk) :: common_part
        do i=1,self%phyto_num
            !dyn_part(i)=-N_uptake(i)*Phy(i)+self%y*grazing_forc(i)*Q_N(i)+self%frac_mn*self%m*Phy(i)*Q_N(i)
	    up(i)=N_uptake(i)*Phy(i)
	!    gr(i)=self%y*grazing_forc(i)*Q_N(i)
	end do
        !common_part =  self%r_dn*D_N*F_T(1)! #+ r_mix *(Nit_bot-N) 
	!dN_dt=common_part+sum(dyn_part(:))
        dN_dt=self%r_dn*D_N*F_T(1)-sum(up(:))!+sum(gr(:))
        return
end function
!------------------------------------------------------------------------------
real(rk) function dP_dt(self,P_uptake, Phy, D_P,grazing_forc,Q_P,F_T,P) ! Change of nutrient concentration over time
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),dimension(self%phyto_num),intent(in) :: P_uptake! Phytoplankton phosphorous uptake rate, mol-P mol-C^-1 d^-1
 real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
 real(rk),intent(in) :: D_P! Phosphorous content of detritus concentration,  mmol-P m^-3
 real(rk),dimension(self%phyto_num),intent(in) :: grazing_forc! Grazing forcing, mmol-C m^-3 d^-1
 real(rk),dimension(self%phyto_num),intent(in) :: Q_P! Phytoplankton intracellular phosphorous cell quota, mol-P mol-C^-1
 real(rk),dimension(2),intent(in) :: F_T! Temperature dependency
 real(rk),intent(in) :: P! Phosphorous concentration, mmol-P m^-3
 integer::i
 real(rk),dimension(self%phyto_num) :: dyn_part,up,gr
 real(rk) :: common_part
	dyn_part=0.0_rk
        do i=1,self%phyto_num
            !dyn_part(i)=-P_uptake(i)*Phy(i)+self%y*grazing_forc(i)*Q_P(i)+self%frac_mn*self%m*Phy(i)*Q_P(i)
!        #P_bot=Nit_bot*(self.pars['P0']/self.pars['N0'])
	    up(i)=P_uptake(i)*Phy(i)
!	    gr(i)=self%y*grazing_forc(i)*Q_P(i)
	end do
!        common_part = self%r_dn*D_P*F_T(1)!#+ r_mix*(P_bot-P)
!	dP_dt=common_part+sum(dyn_part(:))
	dP_dt=self%r_dn*D_P*F_T(1)-sum(up(:))!+sum(gr(:))
return
end function
!------------------------------------------------------------------------------
subroutine Rel_growth_rate_sr(self,Phy,aggregation,growth_rate,grazing_forc,respiration,sinking,rel_growth_rate)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
 real(rk),dimension(self%phyto_num),intent(out) :: rel_growth_rate
 real(rk),intent(in) ::  aggregation! Aggregation rate, d^-1
 real(rk),dimension(self%phyto_num),intent(in) :: growth_rate! Phytoplankton growth rate, d^-1
 real(rk),dimension(self%phyto_num),intent(in) :: grazing_forc! Grazing forcing, mmol-C m^-3 d^-1
 real(rk),dimension(self%phyto_num),intent(in) :: respiration! Phytoplankton respiration rate, d^-1
 real(rk),dimension(self%phyto_num),intent(in) :: sinking! Phytoplankton sinking rate, d^-1
 integer::i
 real(rk), parameter::eps=0.0001_rk !Small parameter to provide division by zero
 rel_growth_rate=0.0_rk
        do i=1,self%phyto_num
            rel_growth_rate(i)=growth_rate(i)-respiration(i)-sinking(i)-self%m-aggregation-grazing_forc(i)/(eps+Phy(i))
	end do
        return 
end subroutine
!------------------------------------------------------------------------------
real(rk) function chl_a(self,Phy,Q_N,par)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),intent(in) :: par
 real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
 real(rk),dimension(self%phyto_num),intent(in) :: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^1
!        :return: Chl_a concentration
 integer::i
 real(rk) :: phyto_conc_tot
 phyto_conc_tot=0.0_rk
	Do i=1,self%phyto_num
        phyto_conc_tot=phyto_conc_tot+Phy(i) * Q_N(i) !* (1._rk-par/self%I_opt)
	end do
        chl_a=phyto_conc_tot*self%chla_to_T_PhyN
        return
end function
!-----------------------------------------------------------------
real(rk) function mean_cell_size(self,Phy)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),dimension(self%phyto_num), intent(in) ::Phy
 real(rk),dimension(self%phyto_num) :: mean_nom
 real(rk) :: Phy_tot
 real(rk), dimension(self%phyto_num) :: Phy_tmp
 integer :: i
!        :param Phy: Phytoplankton biomass concentration, mmol-C m^-3
!        :return: community mean cell size, log_e ESD (mu m)
!        """
mean_nom=0.0_rk
mean_cell_size=0.0_rk
Phy_tmp=Phy
	do i=1,self%phyto_num
	!Larger Diatoms (L>3.5) were counted seperately
	        if(self%log_ESD(i)<self%log_ESD_crit) then! .and. self%log_ESD(i)/= 2.5) then
	        !Error in data for nan IV (2.6)
		  mean_nom(i)=Phy(i)*exp(self%log_ESD(i))
		!else if (self%log_ESD(i)==self%log_ESD_crit) then
		!  mean_nom(i)=0.5*Phy(i)*self%log_ESD(i)
		!  Phy_tmp(i)=0.5*Phy_tmp(i)
		else
		  Phy_tmp(i)=0.0_rk
		end if
	end do
        if (sum(Phy_tmp) /= 0.0) mean_cell_size=log(sum(mean_nom)/sum(Phy_tmp))
        return
end function
!------------------------------------------------------------------------------
real(rk) function size_diversity(self,Phy)
implicit none
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
!        :return: community size diversity, log_e ESD (mu m)^2
 integer::i
 real(rk),dimension(:),allocatable::div_nom
 real(rk)::mean_cellsize
 allocate(div_nom(self%phyto_num))
        do i=1,self%phyto_num
            div_nom(i)=(self%log_ESD(i)-mean_cell_size(self,Phy))**2*Phy(i)
	end do
        size_diversity=sum(div_nom)/sum(Phy)
        return
end function
!------------------------------------------------------------------------------
real(rk) function allometries_esd(beta,alpha,s)
 real(rk), intent(in) :: beta! Interception for trait
 real(rk), intent(in) :: alpha! Size scaling exponent for trait
 real(rk), intent(in) :: s! Equivalent spherical diamater, mu m
!    :return: trait=exp^(log_e(beta)+alpha*log_e(ESD))
    allometries_esd =beta*exp(alpha*s)
    return
end function
!------------------------------------------------------------------------------
subroutine bgc_parameters(self,s, bgc_params)
implicit none!    :return: Phytoplankton eco-physiological traits
 class (type_hzg_kristineb),intent(in) :: self
 real(rk),intent (in) :: s! Equivalent spherical diamater, mu m
 real(rk),dimension(11),intent (out)::bgc_params
 real(rk) :: mu_max, Qmin_N, Qmax_N, vmax_N, Kn_N, affinity_N
 real(rk) :: Qmin_P, Qmax_P, vmax_P, Kn_P, affinity_P,tmp
 
!     Nonlinear mumax
    mu_max=self%pars(3)*min(self%pars(1)*exp(s*self%pars(2)), &
     &  self%pars(4)*exp(s*self%pars(5)))

 !   tmp=self%pars(1)*exp(2.29_rk*self%pars(2))
 !   mu_max=self%pars(3)*min(self%pars(1)*exp(s*self%pars(2)),tmp*exp((s-tmp)*self%pars(5)))
    
  !linear mumax
    !mu_max=allometries_esd(self%pars(26),self%pars(27),s)
    
    Qmin_N=allometries_esd(self%pars(6),self%pars(7),s)
    Qmax_N=allometries_esd(self%pars(8),self%pars(9),s)
    vmax_N =allometries_esd(self%pars(10),self%pars(11),s)

    !if (log(s)<2.3_rk) then
    !    affinity_N=self%pars(12)*dexp(self%pars(13)*2.3_rk)
    !else

    !end if
    !Kn_N =allometries_esd(self%pars(14),self%pars(15),s)
    affinity_N=allometries_esd(self%pars(12),self%pars(13),s)
    Qmin_P=allometries_esd(self%pars(16),self%pars(17),s)
    Qmax_P=allometries_esd(self%pars(18),self%pars(19),s)
    vmax_P =allometries_esd(self%pars(20),self%pars(21),s)
    !Kn_P =allometries_esd(self%pars(22),self%pars(23),s)
    affinity_P=allometries_esd(self%pars(24),self%pars(25),s)

    bgc_params(1) = mu_max
    bgc_params(2) = Qmin_N
    bgc_params(3) = Qmax_N
    bgc_params(4) = vmax_N
    bgc_params(5) = 0.0!Kn_N
    bgc_params(6) = affinity_N

    bgc_params(7) = Qmin_P
    bgc_params(8) = Qmax_P
    bgc_params(9) = vmax_P
    bgc_params(10) = 0.0!Kn_P
    bgc_params(11) = affinity_P
    return
end subroutine
!------------------------------------------------------------------------------
subroutine convert_BGCparams(beta,alpha,a_carbon,b_carbon,ESD_beta,ESD_alpha)
implicit none
!    :param beta_params: Interception for trait
 real(rk),intent(in)::beta
!    :param alpha_params: Size scaling exponent for trait
 real(rk), intent(in)::alpha
 real(rk), intent(in):: a_carbon! size scaling exponent for cell carbon content
 real(rk), intent(in):: b_carbon! Interception for cell carbon content
 real(rk), intent(out)::ESD_beta,ESD_alpha
 real(rk) :: beta_cell_moleC
    beta_cell_moleC=(beta/b_carbon)*12._rk*10_rk**6!  #moleX/moleC
    ESD_beta=beta_cell_moleC*((acos(-1.0_rk)/6.0_rk)**(alpha-a_carbon))
    ESD_alpha=3.0_rk*(alpha-a_carbon)
    return
end subroutine
!------------------------------------------------------------------------------
!Umwandlsung von integer in Character Typ-----------------------------------
 Character*3 Function int2char(i)
implicit none
integer, intent(in) :: i
write(int2char,'(i3)') i
int2char = adjustl(int2char)
End Function int2char