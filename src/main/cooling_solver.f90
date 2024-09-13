!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_solver
!
! Generic module handling analytic cooling functions
!  with known derivatives. These can be solved either
!  explicitly, implicitly or with the Townsend (2009)
!  exact method. Implementation by Lionel Siess.
!
! :References:
!   Townsend (2009), ApJS 181, 391-397
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - T1_factor      : *factor by which T0 is increased (T1= T1_factor*T0)*
!   - bowen_Cprime   : *radiative cooling rate (g.s/cm³)*
!   - dust_collision : *dust collision (1=on/0=off)*
!   - excitation_HI  : *cooling via electron excitation of HI (1=on/0=off)*
!   - lambda_shock   : *Cooling rate parameter for analytic shock solution*
!   - relax_bowen    : *Bowen (diffusive) relaxation (1=on/0=off)*
!   - relax_stefan   : *radiative relaxation (1=on/0=off)*
!   - shock_problem  : *piecewise formulation for analytic shock solution (1=on/0=off)*
!
! :Dependencies: cooling_functions, infile_utils, io, physcon, timestep,
!   units
!

 use cooling_functions, only:bowen_Cprime,lambda_shock_cgs,T0_value,T1_factor
 USE physcon, ONLY:atomic_mass_unit
 implicit none
 character(len=*), parameter :: label = 'cooling_library'
 integer, public :: excitation_HI = 0, relax_Bowen = 0, dust_collision = 0, relax_Stefan = 0, shock_problem = 0
 integer, public :: cooltable_Chris = 0
 integer, public :: icool_method  = 0
 !integer, parameter :: nTg  = 64
 integer, parameter :: nTg  = 1000 ! cloudy
 !integer, parameter :: nTg  = 667 ! cloudy where first T bin spans T_floor=1e4K
 !integer, parameter :: nTg  = 201 ! GS07
 INTEGER :: nTg_Chris
 real,    parameter :: Tref = 1.d7 !higher value of the temperature grid (for exact cooling)
 real :: Tgrid(nTg)
 REAL :: LambdaTable(nTg),alphaTable(nTg),YkTable(nTg)
 REAL :: Tref_Chris
 REAL, DIMENSION(:), ALLOCATABLE :: Tgrid_Chris,LambdaTable_Chris,alphaTable_Chris,YkTable_Chris
 REAL, DIMENSION(:), ALLOCATABLE :: YTsegAlpha1_Chris,YTsegAlphaNot1_Chris
 REAL :: TNdivLN
 
 public :: init_cooling_solver,read_options_cooling_solver,write_options_cooling_solver
 public :: energ_cooling_solver,calc_cooling_rate, calc_Q
 public :: testfunc,print_cooling_rates
 public :: T0_value,lambda_shock_cgs ! expose to cooling module
 logical, public :: Townsend_test = .false. !for analysis_cooling

 private
 real,    parameter :: Tcap = 1.d3 !Townsend cap temperature

 REAL, PARAMETER :: habund=0.7
 REAL, PARAMETER :: mu_e=2.0*atomic_mass_unit/(1.0+habund)
 REAL, PARAMETER :: mu_h=atomic_mass_unit/habund
 !LOGICAL, PARAMETER :: methodLong=.TRUE. !Q-based method, native to Phantom
 LOGICAL, PARAMETER :: methodLong=.FALSE. !Lambda-based method, native to EIS algorithm


contains
!-----------------------------------------------------------------------
!+
!   Initialise cooling functions and check mix of options are sensible
!+
!-----------------------------------------------------------------------
subroutine init_cooling_solver(ierr)
 use io, only:error
 integer, intent(out) :: ierr

 ierr = 0
 !you can't have cool_relaxation_Stefan and cool_relaxation_Bowen at the same time
 if (relax_Bowen == 1 .and. relax_Stefan == 1) then
    call error(label,'you can"t have bowen and stefan cooling at the same time')
    ierr = 1
 endif
 !if no cooling flag activated, disable cooling
 if ( (excitation_HI+relax_Bowen+dust_collision+relax_Stefan+shock_problem) == 0) then
    print *,'ERROR: no cooling prescription activated'
    ierr = 2
    WRITE(*,*) 'NOTE: above is not actually an error -- need to fix this...'
    ierr=0
 endif
 !call set_Tgrid()
 cooltable_Chris=1
 IF(cooltable_Chris==1) THEN
    WRITE(*,*) 'methodLong = ',methodLong
    IF(methodLong) THEN
       WRITE(*,*) 
       WRITE(*,*) 'CALL set_Tgrid_cooltable_Chris()'
       CALL set_Tgrid_cooltable_Chris()
       WRITE(*,*) 'END CALL set_Tgrid_cooltable_Chris()'
    ELSE
       WRITE(*,*) 
       WRITE(*,*) 'CALL set_Tgrid_cooltable_Chris2()'
       CALL set_Tgrid_cooltable_Chris2()
       WRITE(*,*) 'END CALL set_Tgrid_cooltable_Chris2()'
    ENDIF
    WRITE(*,*) 
    WRITE(*,*) 'mu_e, mu_H =',mu_e,mu_H
    WRITE(*,*) 
 ELSE
    WRITE(*,*) 'call set_Tgrid()'
    call set_Tgrid()
 ENDIF

end subroutine init_cooling_solver

!-----------------------------------------------------------------------
!+
!   Get right hand side of energy equation for the desired
!   cooling prescription and choice of solver
!+
!-----------------------------------------------------------------------
subroutine energ_cooling_solver(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa)
 real, intent(in)  :: ui,rho,dt                ! in code units
 real, intent(in)  :: Tdust,mu,gamma,K2,kappa  ! in cgs
 real, intent(out) :: dudt                     ! in code units

 if (icool_method == 2) then
    !call exact_cooling   (ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa)
    IF(methodLong) THEN
       call exact_cooling_Chris(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa)
    ELSE
       call exact_cooling_Chris2(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa)
    ENDIF
 elseif (icool_method == 0) then
    call implicit_cooling(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa)
 else
    call explicit_cooling(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa)
 endif

end subroutine energ_cooling_solver

!-----------------------------------------------------------------------
!+
!   explicit cooling
!+
!-----------------------------------------------------------------------
subroutine explicit_cooling (ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa)

 use physcon, only:Rg
 use units,   only:unit_ergg

 real, intent(in)  :: ui, rho, dt, Tdust, mu, gamma !code units
 real, intent(in)  :: K2, kappa
 real, intent(out) :: dudt                         !code units

 real              :: u,Q,dlnQ_dlnT,T,T_on_u

 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T      = T_on_u*ui
 call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Tdust, mu, gamma, K2, kappa)
 if (ui + Q*dt < 0.) then   ! assume thermal equilibrium
    if (Townsend_test) then
       !special fix for Townsend benchmark
       u = Tcap/T_on_u
    else
       u = Tdust/T_on_u     ! set T=Tdust
    endif
    dudt = (u-ui)/dt
 else
    dudt = Q
 endif

end subroutine explicit_cooling

!-----------------------------------------------------------------------
!+
!   implicit cooling
!+
!-----------------------------------------------------------------------
subroutine implicit_cooling (ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa)

 use physcon, only:Rg
 use units,   only:unit_ergg

 real, intent(in)  :: ui, rho, dt, mu, gamma
 real, intent(in)  :: Tdust, K2, kappa
 real, intent(out) :: dudt

 real, parameter    :: tol = 1.d-6, Tmin = 1.
 integer, parameter :: iter_max = 40
 real               :: u,Q,dlnQ_dlnT,T_on_u,Qi,f0,fi,fmid,T,T0,dx,Tmid
 integer            :: iter

 u       = ui
 T_on_u  = (gamma-1.)*mu*unit_ergg/Rg
 T       = ui*T_on_u
 call calc_cooling_rate(Q,dlnQ_dlnT, rho, T, Tdust, mu, gamma, K2, kappa)
 !cooling negligible, return
 if (abs(Q) < tiny(0.)) then
    dudt = 0.
    return
 endif
 T0   = T
 f0   = -Q*dt*T_on_u
 fi   = f0
 iter = 0
 !define bisection interval for function f(T) = T^(n+1)-T^n-Q*dt*T_on_u
 do while (((f0 > 0. .and. fi > 0.) .or. (f0 < 0. .and. fi < 0.)) .and. iter < iter_max)
    Tmid = max(T+Q*dt*T_on_u,Tmin)
    call calc_cooling_rate(Qi,dlnQ_dlnT, rho, Tmid, Tdust, mu, gamma, K2, kappa)
    fi = Tmid-T0-Qi*dt*T_on_u
    T  = Tmid
    iter = iter+1
 enddo
 !Temperature is between T0 and Tmid
 if (iter > iter_max) stop '[implicit_cooling] cannot bracket cooling function'
 iter = 0
 if (Tmid > T0) then
    T = T0
 else
    if (Townsend_test) then
       !special fix for Townsend benchmark
       T = max(Tcap,Tmid)
    else
       T = Tmid
    endif
 endif
 dx = abs(Tmid-T0)
 do while (dx/T0 > tol .and. iter < iter_max)
    dx = dx*.5
    Tmid = T+dx
    call calc_cooling_rate(Qi,dlnQ_dlnT, rho, Tmid, Tdust, mu, gamma, K2, kappa)
    fmid = Tmid-T0-Qi*dt*T_on_u
    if (Townsend_test) then
       !special fix for Townsend benchmark
       if (fmid <= 0.) Tmid = max(Tcap,Tmid)
    else
       if (fmid <= 0.) T = Tmid
    endif
    iter = iter + 1
    !print *,iter,fmid,T,Tmid
 enddo
 u = Tmid/T_on_u
 dudt =(u-ui)/dt
 if (u < 0. .or. isnan(u)) then
    print *,u
    stop '[implicit_cooling] u<0'
 endif

end subroutine implicit_cooling


!-----------------------------------------------------------------------
!+
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!+
!-----------------------------------------------------------------------
subroutine exact_cooling(ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa)

 use physcon, only:Rg
 use units,   only:unit_ergg

 real, intent(in)  :: ui, rho, dt, Tdust, mu, gamma
 real, intent(in)  :: K2, kappa
 real, intent(out) :: dudt

 real, parameter :: tol = 1.d-12
 real            :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,T_on_u,T_floor,Qi
 integer         :: k

 if (Townsend_test) then
    T_floor = Tcap
 else
    T_floor = 10.
 endif
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T      = T_on_u*ui

 if (T < T_floor) then
    Temp = T_floor
 elseif (T > Tref) then
    call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Tdust, mu, gamma, K2, kappa)
    Temp = T+T_on_u*Q*dt
 else
    call calc_cooling_rate(Qref,dlnQref_dlnT, rho, Tref, Tdust, mu, gamma, K2, kappa)
    Qi = Qref
    Y         = 0.
    k         = nTg
    Q         = Qref          ! default value if Tgrid < T for all k
    dlnQ_dlnT = dlnQref_dlnT  ! default value if Tgrid < T for all k
    do while (Tgrid(k) > T)
       k = k-1
       call calc_cooling_rate(Q, dlnQ_dlnT, rho, Tgrid(k), Tdust, mu, gamma, K2, kappa)
       dlnQ_dlnT = log(Qi/Q)/log(Tgrid(k+1)/Tgrid(k))
       Qi = Q
       ! eqs A6 to get Yk
       if (abs(dlnQ_dlnT-1.) < tol) then
          y = y - Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/Tgrid(k+1))
       else
          y = y - Qref*Tgrid(k)/(Q*Tref*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
       endif
    enddo
    !eqs A5 for Y(T)
    yk = y
    if (abs(dlnQ_dlnT-1.) < tol) then
       y = yk + Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/T)
    else
       y = yk + Qref*Tgrid(k)/((Q*Tref)*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1.))
    endif
    !argument of Y^(-1) in eq 26
    dy = -Qref*dt*T_on_u/Tref
    y  = y + dy
    !find new k for eq A7 (not necessarily the same as k for eq A5)
    do while(y>yk .AND. k>1)
       k = k-1
       call calc_cooling_rate(Q, dlnQ_dlnT, rho, Tgrid(k), Tdust, mu, gamma, K2, kappa)
       dlnQ_dlnT = log(Qi/Q)/log(Tgrid(k+1)/Tgrid(k))
       Qi = Q
       ! eqs A6 to get Yk
       if (abs(dlnQ_dlnT-1.) < tol) then
          yk = yk - Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/Tgrid(k+1))
       else
          yk = yk - Qref*Tgrid(k)/(Q*Tref*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
       endif
    enddo
    !compute Yinv (eqs A7)
    if (abs(dlnQ_dlnT-1.) < tol) then
       Temp = max(Tgrid(k)*exp(-Q*Tref*(y-yk)/(Qref*Tgrid(k))),T_floor)
    else
       Yinv = 1.-(1.-dlnQ_dlnT)*Q*Tref/(Qref*Tgrid(k))*(y-yk)
       if (Yinv > 0.) then
          Temp = max(Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT))),T_floor)
       else
          Temp = T_floor
       endif
    endif
 endif

 dudt = (Temp-T)/T_on_u/dt
 !note that u = Temp/T_on_u

WRITE(*,*) ui,T,Temp,dudt,Q
end subroutine exact_cooling

!-----------------------------------------------------------------------
!+
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!+
!-----------------------------------------------------------------------
subroutine exact_cooling_Chris(ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa)

 !use physcon, only:Rg
 use physcon, only:Rg,kboltz
 !use units,   only:unit_ergg
 use units,   only:unit_ergg,unit_density,utime
 !use cooling, only:Tfloor

 real, intent(in)  :: ui, rho, dt, Tdust, mu, gamma
 real, intent(in)  :: K2, kappa
 real, intent(out) :: dudt

 real, parameter :: tol = 1.d-12
 real            :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,T_on_u,T_floor,Qi
 integer         :: k

 REAL :: rhocgsDivmueDivmuH
 REAL :: tcool,tcoolfac
 REAL :: Tref_Chris2
 LOGICAL :: opt0,opt1,opt2,withCorrection
 REAL, PARAMETER :: Tfloor=1.d4
 INTEGER :: ksave
 REAL :: yksave,ysave
 
 !Option 0: default calculation where ref=N
 opt0=.TRUE.
 opt1=.FALSE.
 !Option 1: new calcualiton where ref=k+1
 !opt0=.FALSE.
 !opt1=.TRUE.
 !Option 2: new calcualiton where ref=k
 !opt0=.FALSE.
 !opt1=.FALSE.
 !opt2=.TRUE.
 !Option 3: new calcualiton where ref=k,Y_k=0
 !opt0=.FALSE.
 !opt1=.FALSE.
 !opt2=.FALSE.
 
 !compute k' or not -- "withCorrection=.TRUE." means to copmpute k'
 !withCorrection=.FALSE.
 withCorrection=.TRUE.
 
 !WRITE(*,*) 'Rg, kB/amu, kB, amu =',Rg,kboltz/atomic_mass_unit,kboltz,atomic_mass_unit
 rhocgsDivmueDivmuH = rho*unit_density/mu_e/mu_H
 !tcoolfac = kboltz*mu_e*mu_H / ((gamma-1.)*rho*unit_density*mu*atomic_mass_unit)
 tcoolfac = Rg*mu_e*mu_H / ((gamma-1.)*rho*unit_density*mu)
 !WRITE(*,*) 'asdf: ',kboltz,mu_e,mu_H,gamma,rho,unit_density,mu,atomic_mass_unit,tcoolfac
 !WRITE(*,*) 'asdf2:',ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa

 if (Townsend_test) then
    T_floor = Tcap
 else
    !T_floor = 10.
    T_floor = Tfloor
 endif
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T      = T_on_u*ui
 !WRITE(*,*) 'asdf3:',T_on_u,gamma,mu,unit_ergg,Rg

 if (T < T_floor) then
    WRITE(*,*) 'T < T_floor'
    Temp = T_floor
 elseif (T > Tref_Chris) then
    WRITE(*,*) 'T > Tref_Chris'
    call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Tdust, mu, gamma, K2, kappa)
    Temp = T+T_on_u*Q*dt
 else
    !WRITE(*,*) 'else'
    !WRITE(*,*) 'else -- y=', y
    !WRITE(*,*) 'else -- mu=', mu
    !WRITE(*,*) 'else -- habund, mu_e, mu_h=', habund,mu_e, mu_h
    !WRITE(*,*) 'else -- rho, rho*unit_density, rhocgsDivmueDivmuH=',rho,rho*unit_density,rhocgsDivmueDivmuH,LambdaTable(nTg)
    WRITE(*,*) 'else -- opt0 =',opt0,', opt1 =',opt1,', opt2 =',opt2,', withCorrection =',withCorrection
    IF(opt0) THEN
       !call calc_cooling_rate(Qref,dlnQref_dlnT, rho, Tref_Chris, Tdust, mu, gamma, K2, kappa)
       !Qref=LambdaTable(nTg)*-2.e23
       Qref=-rhocgsDivmueDivmuH*LambdaTable(nTg) * utime/unit_ergg
       dlnQref_dlnT=alphaTable(nTg)
       !dlnQref_dlnT=0.
       Qi = Qref
       Y         = 0.
       k         = nTg
       Q         = Qref          ! default value if Tgrid < T for all k
       dlnQ_dlnT = dlnQref_dlnT  ! default value if Tgrid < T for all k
       do while (Tgrid(k) > T)
          k = k-1
          !call calc_cooling_rate(Q, dlnQ_dlnT, rho, Tgrid(k), Tdust, mu, gamma, K2, kappa)
          !Q=LambdaTable(k)*-2.e23
          Q=-rhocgsDivmueDivmuH*LambdaTable(k) * utime/unit_ergg
          dlnQ_dlnT=alphaTable(k)
          !dlnQ_dlnT=0.
          !dlnQ_dlnT = log(Qi/Q)/log(Tgrid(k+1)/Tgrid(k))
          Qi = Q
          ! eqs A6 to get Yk
          if (abs(dlnQ_dlnT-1.) < tol) then
             y = y - Qref*Tgrid(k)/(Q*Tref_Chris)*log(Tgrid(k)/Tgrid(k+1))
          else
             y = y - Qref*Tgrid(k)/(Q*Tref_Chris*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
          endif
       enddo
    ELSEIF(opt1) THEN !Y_{k+1}=0, T_N=T_{k+1}, Lambda_N=Lambda_{k+1}
       k         = nTg
       do while (Tgrid(k) > T)
          k = k-1
       enddo
       Qref=-rhocgsDivmueDivmuH*LambdaTable(k+1) * utime/unit_ergg
       dlnQref_dlnT=alphaTable(k+1)
       Tref_Chris2=Tgrid(k+1) !need to save since k might be changed later for Eq. A7
       !Y         = 0.
       Q=-rhocgsDivmueDivmuH*LambdaTable(k) * utime/unit_ergg
       dlnQ_dlnT=alphaTable(k)
       if (abs(dlnQ_dlnT-1.) < tol) then
          y = - Qref*Tgrid(k)/(Q*Tgrid(k+1))*log(Tgrid(k)/Tgrid(k+1))
       else
          y = - Qref*Tgrid(k)/(Q*Tgrid(k+1)*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
       endif
    ELSEIF(opt2) THEN !Y_{k+1}=0, T_N=T_k, Lambda_N=Lambda_k
       k         = nTg
       do while (Tgrid(k) > T)
          k = k-1
       enddo
       Qref=-rhocgsDivmueDivmuH*LambdaTable(k) * utime/unit_ergg
       dlnQref_dlnT=alphaTable(k)
       Tref_Chris2=Tgrid(k) !need to save since k might be changed later for Eq. A7
       !Y         = 0.
       Q=Qref
       dlnQ_dlnT=dlnQref_dlnT
       if (abs(dlnQ_dlnT-1.) < tol) then
          y = - log(Tgrid(k)/Tgrid(k+1))
       else
          y = - 1./(1.-dlnQ_dlnT)*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
       endif
    ELSE !Y_{k+1}=0, T_N=T_k, Lambda_N=Lambda_k
       k         = nTg
       do while (Tgrid(k) > T)
          k = k-1
       enddo
       Qref=-rhocgsDivmueDivmuH*LambdaTable(k) * utime/unit_ergg
       dlnQref_dlnT=alphaTable(k)
       Tref_Chris2=Tgrid(k) !need to save since k might be changed later for Eq. A7
       !Y         = 0.
       Q=Qref
       dlnQ_dlnT=dlnQref_dlnT
       !if (abs(dlnQ_dlnT-1.) < tol) then
       !   y = - log(Tgrid(k)/Tgrid(k+1))
       !else
       !   y = - 1./(1.-dlnQ_dlnT)*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
       !endif
       y=0.
    ENDIF
    tcool=tcoolfac*T/LambdaTable(k)
    !WRITE(*,*) 'tcool =',tcool,tcool/utime
    WRITE(*,*) 'tcool (s,yr,code) =',tcool,tcool/(365.25*24.*3600.),tcool/utime,k,LambdaTable(k),T,utime,dt,T_floor
    IF(opt0) THEN
       !eqs A5 for Y(T)
       yk = y
       if (abs(dlnQ_dlnT-1.) < tol) then
          y = yk + Qref*Tgrid(k)/(Q*Tref_Chris)*log(Tgrid(k)/T)
       else
          y = yk + Qref*Tgrid(k)/((Q*Tref_Chris)*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1.))
       endif
       !argument of Y^(-1) in eq 26
       dy = -Qref*dt*T_on_u/Tref_Chris
    ELSEIF(opt1) THEN !Y_{k+1}=0, T_N=T_{k+1}, Lambda_N=Lambda_{k+1}
       !eqs A5 for Y(T)
       yk = y
       if (abs(dlnQ_dlnT-1.) < tol) then
          !y = Qref*Tgrid(k)/(Q*Tref_Chris2)*(log(Tgrid(k)/T)-log(Tgrid(k)/Tref_Chris2))
          y = yk + Qref*Tgrid(k)/(Q*Tref_Chris2)*log(Tgrid(k)/T)
       else
          !y = Qref*Tgrid(k)/(Q*Tref_Chris2*(1.-dlnQ_dlnT))*((Tgrid(k)/Tref_Chris2)**(dlnQ_dlnT-1.)-(Tgrid(k)/T)**(dlnQ_dlnT-1.))
          y = yk + Qref*Tgrid(k)/(Q*Tref_Chris2*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1.))
       endif
       !argument of Y^(-1) in eq 26
       dy = -Qref*dt*T_on_u/Tref_Chris2
    ELSE !Y_{k+1}=0, T_N=T_k, Lambda_N=Lambda_k
       !eqs A5 for Y(T)
       yk = y
       if (abs(dlnQ_dlnT-1.) < tol) then
          y = yk + log(Tgrid(k)/T)
       else
          y = yk + 1./(1.-dlnQ_dlnT)*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1.))
       endif
       !argument of Y^(-1) in eq 26
       dy = -Qref*dt*T_on_u/Tref_Chris2
    ENDIF
    ksave = k
    yksave = yk
    ysave = y
    y  = y + dy

    !New Part -- Start
    IF(withCorrection) THEN
       do while(y>yk .AND. k>1)
          k = k-1
          !call calc_cooling_rate(Q, dlnQ_dlnT, rho, Tgrid(k), Tdust, mu, gamma, K2, kappa)
          !Q=LambdaTable(k)*-2.e23
          Q=-rhocgsDivmueDivmuH*LambdaTable(k) * utime/unit_ergg
          dlnQ_dlnT=alphaTable(k)
          !dlnQ_dlnT=0.
          !dlnQ_dlnT = log(Qi/Q)/log(Tgrid(k+1)/Tgrid(k))
          Qi = Q
          IF(opt0) THEN
             ! eqs A6 to get Yk
             if (abs(dlnQ_dlnT-1.) < tol) then
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris)*log(Tgrid(k)/Tgrid(k+1))
             else
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
             endif
          ELSEIF(opt1) THEN !Y_{k+1}=0, T_N=T_{k+1}, Lambda_N=Lambda_{k+1}
             ! eqs A6 to get Yk
             if (abs(dlnQ_dlnT-1.) < tol) then
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris2)*log(Tgrid(k)/Tgrid(k+1))
             else
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris2*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
             endif
          ELSE !Y_{k+1}=0, T_N=T_k, Lambda_N=Lambda_k
             ! eqs A6 to get Yk
             if (abs(dlnQ_dlnT-1.) < tol) then
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris2)*log(Tgrid(k)/Tgrid(k+1))
             else
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris2*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
             endif
          ENDIF
       enddo
    ENDIF
    !New Part -- End
    !WRITE(*,*) 'asdf',k,y,yk,Tgrid(k)
    
    IF(opt0) THEN
       !compute Yinv (eqs A7)
       if (abs(dlnQ_dlnT-1.) < tol) then
          Temp = max(Tgrid(k)*exp(-Q*Tref_Chris*(y-yk)/(Qref*Tgrid(k))),T_floor)
       else
          Yinv = 1.-(1.-dlnQ_dlnT)*Q*Tref_Chris/(Qref*Tgrid(k))*(y-yk)
          if (Yinv > 0.) then
             !Temp = Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT)))
             Temp = MAX(Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT))),T_floor)
             !WRITE(*,*) 'hmm, Temp =',Temp,Tgrid(k),k,Yinv
             !WRITE(*,*) 'confused ',Temp<T_floor
          else
             Temp = T_floor
             !WRITE(*,*) 'AtFloor',Temp,Yinv,k
          endif
       endif
    ELSEIF(opt1) THEN !Y_{k+1}=0, T_N=T_{k+1}, Lambda_N=Lambda_{k+1}
       !compute Yinv (eqs A7)
       if (abs(dlnQ_dlnT-1.) < tol) then
          Temp = max(Tgrid(k)*exp(-Q*Tref_Chris2*(y-yk)/(Qref*Tgrid(k))),T_floor)
       else
          Yinv = 1.-(1.-dlnQ_dlnT)*Q*Tref_Chris2/(Qref*Tgrid(k))*(y-yk)
          if (Yinv > 0.) then
             !Temp = Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT)))
             Temp = MAX(Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT))),T_floor)
             !WRITE(*,*) 'hmm, Temp =',Temp,Tgrid(k),k,Yinv
             !WRITE(*,*) 'confused ',Temp<T_floor
          else
             Temp = T_floor
             !WRITE(*,*) 'AtFloor',Temp,Yinv,k
          endif
       endif
    ELSE !Y_{k+1}=0, T_N=T_k, Lambda_N=Lambda_k
       !compute Yinv (eqs A7)
       if (abs(dlnQ_dlnT-1.) < tol) then
          Temp = max(Tgrid(k)*exp(-Q*Tref_Chris2*(y-yk)/(Qref*Tgrid(k))),T_floor)
       else
          Yinv = 1.-(1.-dlnQ_dlnT)*Q*Tref_Chris2/(Qref*Tgrid(k))*(y-yk)
          if (Yinv > 0.) then
             !Temp = Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT)))
             Temp = MAX(Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT))),T_floor)
             !WRITE(*,*) 'hmm, Temp =',Temp,Tgrid(k),k,Yinv
             !WRITE(*,*) 'confused ',Temp<T_floor
          else
             Temp = T_floor
             !WRITE(*,*) 'AtFloor',Temp,Yinv,k
          endif
       endif
    ENDIF
 endif

 IF(Temp<T_floor) WRITE(*,*) 'UH-HOH ',Temp,T_floor
 WRITE(*,*) 'T, Temp =',T,Temp
 dudt = (Temp-T)/T_on_u/dt
 !note that u = Temp/T_on_u

 !WRITE(*,'(A,6(1PE14.6))') 'C',ui,T,Temp,dudt,Q,Tref_Chris
 !WRITE(*,*) 'stuff: ',ksave,k,yksave,ysave,tcoolfac,dy,y,tcoolfac*T/(LambdaTable(ksave)*(T/Tgrid(ksave))**alphaTable(ksave)),tcoolfac*T/LambdaTable(ksave)
 !WRITE(*,*) 'ffuts: ',dy,dt,utime,tcoolfac,T_on_u,Tref_Chris,-Qref
end subroutine exact_cooling_Chris

!-----------------------------------------------------------------------
!+
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!+
!-----------------------------------------------------------------------
subroutine exact_cooling_Chris2(ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa)

 !use physcon, only:Rg
 use physcon, only:Rg,kboltz
 !use units,   only:unit_ergg
 use units,   only:unit_ergg,unit_density,utime
 !use cooling, only:Tfloor
!USE eos, only:gmw

 real, intent(in)  :: ui, rho, dt, Tdust, mu, gamma
 real, intent(in)  :: K2, kappa
 real, intent(out) :: dudt

 real, parameter :: tol = 1.d-12
 !real            :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,T_on_u,T_floor,Qi
 real            :: Yinv,Temp,T,T_on_u,T_floor
 integer         :: k

 !REAL :: rhocgsDivmueDivmuH
 REAL :: tcool,tcoolfac
 !REAL :: Tref_Chris2
 !LOGICAL :: opt0,opt1,opt2,withCorrection
 REAL, PARAMETER :: Tfloor=1.d4

 INTEGER :: kl,km,ku,kk
 REAL :: YT
 !REAL :: YTsave

 !The following 3 lines confirm that eos_vars(imu,i) is passed all the way to mu in this subroutine 
 !if (mu.ne.gmw) then
 !   WRITE(*,*) 'mu =',mu
 !endif

 !NOTE: When this cooling subroutine was written and tested,
 !      kboltz was not updated to enough sig figs, so computations using kboltz
 !      vs. computations using Rg*amu yielded different results.
 !At the time of testing, from src/main/physcon.f90:
 !   real(kind=8), parameter :: Rg = 8.31446261815324d7             !Gas constant              erg/K/g
 !   real(kind=8), parameter :: atomic_mass_unit = 1.660538921d-24  !Atomic mass unit          g
 !   real(kind=8), parameter :: kboltz = 1.38066d-16                !Boltzmann constant        erg/K
 !The following result should be 1.000000000, but it isn't:
 !   Rg*atomic_mass_unit/kboltz = 0.999991944768663
 !Though the difference is small, performing "diff" on files generated from kboltz 
 !   vs. files gerenated from Rg*amu will certainly generate differences.
 !
 !From the 2019 SI definitions, it should be
 !   kB = 1.380649d-16
 !   NA = 6.02214076d23
 !which therefore means that
 !   atomic_mass_unit = 1.66053906717385d-24 = 1/NA = amu
 !   Rg = 8.31446261815324d+07 = kB/amu
 !
 !Therefore, instead of using "kboltz" in calculations, use "Rg*atomic_mass_unit".
 !This is why tcoolfac below is currently written the way it is:
 !tcoolfac = kB mu_e mu_H / ((gamma-1) rho_cgs mu_cgs) = tcool*Lambda(T)/T, eq. 13
 !         = Rg*amu mu_e mu_H / ((gamma-1) rho_cgs mu_1*amu)
 !         = Rg mu_e mu_H / ((gamma-1) rho_cgs mu_1)

 !WRITE(*,*) 'Rg, kB/amu, kB, amu =',Rg,kboltz/atomic_mass_unit,kboltz,atomic_mass_unit
 !rhocgsDivmueDivmuH = rho*unit_density/mu_e/mu_H
 !tcoolfac = kboltz*mu_e*mu_H / ((gamma-1.)*rho*unit_density*mu*atomic_mass_unit) ! = tcool*Lambda(T)/T
 tcoolfac = Rg*mu_e*mu_H / ((gamma-1.)*rho*unit_density*mu) ! = tcool*Lambda(T)/T
 !WRITE(*,*) 'asdf: ',kboltz,mu_e,mu_H,gamma,rho,unit_density,mu,atomic_mass_unit,tcoolfac
 !WRITE(*,*) 'asdf2:',ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa

 if (Townsend_test) then
    T_floor = Tcap
 else
    !T_floor = 10.
    T_floor = Tfloor
 endif
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg ! = (gamma-1.)*mu*atomic_mass_unit/kB * unit_ergg = (gamma-1.)*mu_cgs/kB * unit_ergg
 T      = T_on_u*ui
 !WRITE(*,*) 'asdf3:',T_on_u,gamma,mu,unit_ergg,Rg

 if (T < T_floor) then
    !WRITE(*,*) 'T < T_floor'
    Temp = T_floor
 !elseif (T > Tref_Chris) then
 !   WRITE(*,*) 'T > Tref_Chris'
 !   call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Tdust, mu, gamma, K2, kappa)
 !   Temp = T+T_on_u*Q*dt
 !elseif (.true.) then !these two lines cool all particles to the floor temperature
 !   Temp = T_floor
 else
    !WRITE(*,*) 'else'

    !bisector to find k in Eq. A5
    !c----    Locate the interval in which temp is found.
    !c        This part is taken from locate.f in Numerical Recipes.
    !c
    !ignoring the initial conditionals, the results will be kl_init <= k <= ku_init-1
    !   therefore, choose kl_init = 1 and ku_init = nTg_Chris+1
    IF(T.LE.Tgrid_Chris(1)) THEN
       k = 1
    ELSEIF(T.GE.Tgrid_Chris(nTg_Chris)) THEN
       k = nTg_Chris
    ELSE
       kl = 1
       ku = nTg_Chris+1
       DO WHILE (ku-kl.GT.1)
          km = (ku+kl)/2
          IF(T.GE.Tgrid_Chris(km)) THEN
             kl = km
          ELSE
             ku = km
          ENDIF
       ENDDO
       k = kl
    ENDIF
    !WRITE(*,*) 'tcool = ',tcoolfac*T/(LambdaTable_Chris(k)*(T/Tgrid_Chris(k))**alphaTable_Chris(k)), &
    !                      tcoolfac*T/(LambdaTable_Chris(k)*(T/Tgrid_Chris(k))**alphaTable_Chris(k)) / utime, &
    !                      k,nTg_Chris, &
    !                      T,Tgrid_Chris(k), &
    !                      LambdaTable_Chris(k)*(T/Tgrid_Chris(k))**alphaTable_Chris(k),LambdaTable_Chris(k), &
    !                      alphaTable_Chris(k)

    !find Y(T), Eq. A5
    IF(ABS(alphaTable_Chris(k)-1.) < tol) THEN
       YT = YkTable_Chris(k) + YTsegAlpha1_Chris(k)*LOG(Tgrid_Chris(k)/T)
    ELSE
       YT = YkTable_Chris(k) + YTsegAlphaNot1_Chris(k)*(1.-(Tgrid_Chris(k)/T)**(alphaTable_Chris(k)-1.))
    ENDIF
    !YTsave = YT
    !WRITE(*,*) 'YT1 =',YT,YkTable_Chris(k),YTsegAlphaNot1_Chris(k),Tgrid_Chris(k),T,alphaTable_Chris(k), &
    !                   1.-(Tgrid_Chris(k)/T)**(alphaTable_Chris(k)-1.), &
    !                   YTsegAlphaNot1_Chris(k)*(1.-(Tgrid_Chris(k)/T)**(alphaTable_Chris(k)-1.))
    
    !argument of Eq. 26
    !cCool_v2b   delta_Y = (gamma-1)*rho*mu / (mu_e*mu_H*kB) * LambdaN/TempN * delta_t
    !         dyfunx_v2b = gamma1*DBLE(rho(ipart)*udens)
    !     &                *DBLE(amunit*gmw/(amue*amuh*boltz))
    !     &                /TNdivLN_v2b
    !     &                *DBLE(deltat*utime)
    !yfunx_v2b = yfunx_v2b + dyfunx_v2b
    !dYT = 1./(tcoolfac*TNdivLN) * dt*utime
    !yT = yT + dyT
    YT = YT + 1./(tcoolfac*TNdivLN) * dt*utime
    !WRITE(*,*) 'YT2 =',YT,YkTable_Chris(1),YkTable_Chris(nTg_Chris), &
    !                   tcoolfac,TNdivLN,dt,utime, &
    !                   1./(tcoolfac*TNdivLN),dt*utime, &
    !                   1./(tcoolfac*TNdivLN) * dt*utime

    !bisector to find kk in Eq. A7
    !cCool_v2b Bisection Method to find kk
    !cCool_v2b Note: 'k' is not used again, so 'k' could be used again for bisecting Y,
    !cCool_v2b        but for clarity, a new variable 'kk' is used here.
    !c----    Locate the interval in which temp is found.
    !c        This part is taken from locate.f in Numerical Recipes.
    !c
    !ignoring the initial conditionals, the results will be kl_init <= kk <= ku_init-1
    !   therefore, choose kl_init = 1 and ku_init = nTg_Chris+1
    !IF (yfunx_v2b.EQ.yfunc_v2b(1)) THEN
    IF(YT.GT.YkTable_Chris(1)) THEN
       kk = 0
    ELSEIF(YT.EQ.YkTable_Chris(1)) THEN
       kk = 1
    ELSEIF(YT.LE.YkTable_Chris(nTg_Chris)) THEN
       kk = nTg_Chris
    ELSE
       kl = 1
       ku = nTg_Chris+1
       DO WHILE(ku-kl.GT.1)
          km = (ku+kl)/2
          IF(YT.LE.YkTable_Chris(km)) THEN
             kl = km
          ELSE
             ku = km
          ENDIF
       ENDDO
       kk = kl
    ENDIF

    ! find Yinv, Eq. A7
    IF(kk==0) THEN
       Temp = T_floor
    ELSE
       IF(ABS(alphaTable_Chris(kk)-1.) < tol) THEN
          Temp = MAX(Tgrid_Chris(kk)*EXP(-(YT-YkTable_Chris(kk))/YTsegAlpha1_Chris(kk)) &
                     ,T_floor)
       ELSE
          Yinv = 1.-(YT-YkTable_Chris(kk))/YTsegAlphaNot1_Chris(kk)
          IF(Yinv.GT.0.0) THEN
             Temp = MAX(Tgrid_Chris(kk)*Yinv**(1./(1.-alphaTable_Chris(kk))) &
                        ,T_floor)
          ELSE
             Temp = T_floor
          ENDIF
       ENDIF
    ENDIF
 endif

 IF(Temp<T_floor) WRITE(*,*) 'UH-HOH ',Temp,T_floor
 !WRITE(*,*) 'T, Temp =',T,Temp
 dudt = (Temp-T)/T_on_u/dt
 !note that u = Temp/T_on_u

 !WRITE(*,'(A,6(1PE14.6))') 'C',ui,T,Temp,dudt,Q,Tref_Chris
 !WRITE(*,*) 'stuff: ',k,kk,YkTable_Chris(k),YTsave,tcoolfac,1./(tcoolfac*TNdivLN)*dt*utime,YT,tcoolfac*T/(LambdaTable_Chris(k)*(T/Tgrid_Chris(k))**alphaTable_Chris(k))
 !WRITE(*,*) 'ffuts: ',1./(tcoolfac*TNdivLN)*dt*utime,dt,utime,tcoolfac,T_on_u,Tref_Chris,TNdivLN
end subroutine exact_cooling_Chris2

!-----------------------------------------------------------------------
!+
!  calculate cooling rates
!+
!-----------------------------------------------------------------------
subroutine calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Teq, mu, gamma, K2, kappa)
 use units,   only:unit_ergg,unit_density,utime
 use physcon, only:mass_proton_cgs
 use cooling_functions, only:cooling_neutral_hydrogen,&
     cooling_Bowen_relaxation,cooling_dust_collision,&
     cooling_radiative_relaxation,piecewise_law,testing_cooling_functions
 !use cooling_molecular, only:do_molecular_cooling,calc_cool_molecular

 real, intent(in)  :: rho, T, Teq     !rho in code units
 real, intent(in)  :: mu, gamma
 real, intent(in)  :: K2, kappa       !cgs
 real, intent(out) :: Q, dlnQ_dlnT    !code units

 real :: Q_cgs,Q_H0, Q_relax_Bowen, Q_col_dust, Q_relax_Stefan, Q_molec, Q_shock
 real :: dlnQ_H0, dlnQ_relax_Bowen, dlnQ_col_dust, dlnQ_relax_Stefan, dlnQ_molec, dlnQ_shock
 real :: rho_cgs, ndens

 rho_cgs           = rho*unit_density
 ndens             = rho_cgs/mass_proton_cgs

 Q_H0              = 0.
 Q_relax_Bowen     = 0.
 Q_col_dust        = 0.
 Q_relax_Stefan    = 0.
 Q_shock           = 0.
 Q_molec           = 0.

 dlnQ_H0           = 0.
 dlnQ_relax_Bowen  = 0.
 dlnQ_col_dust     = 0.
 dlnQ_relax_Stefan = 0.
 dlnQ_shock        = 0.
 dlnQ_molec        = 0.

 if (excitation_HI  == 1) call cooling_neutral_hydrogen(T, rho_cgs, Q_H0, dlnQ_H0)
 if (relax_Bowen    == 1) call cooling_Bowen_relaxation(T, Teq, rho_cgs, mu, gamma, &
                                                        Q_relax_Bowen, dlnQ_relax_Bowen)
 if (dust_collision == 1 .and. K2 > 0.) call cooling_dust_collision(T, Teq, rho_cgs, K2,&
                                                        mu, Q_col_dust, dlnQ_col_dust)
 if (relax_Stefan   == 1) call cooling_radiative_relaxation(T, Teq, kappa, Q_relax_Stefan,&
                                                        dlnQ_relax_Stefan)
 if (shock_problem  == 1) call piecewise_law(T, T0_value, rho_cgs, ndens, Q_H0, dlnQ_H0)

 if (excitation_HI  == 99) call testing_cooling_functions(int(K2), T, Q_H0, dlnQ_H0)
 !if (do_molecular_cooling) call calc_cool_molecular(T, r, rho_cgs, Q_molec, dlnQ_molec)

 Q_cgs = Q_H0 + Q_relax_Bowen + Q_col_dust + Q_relax_Stefan + Q_molec + Q_shock
 if (Q_cgs == 0.) then
    dlnQ_dlnT = 0.
 else
    dlnQ_dlnT = (Q_H0*dlnQ_H0 + Q_relax_Bowen*dlnQ_relax_Bowen + Q_col_dust*dlnQ_col_dust&
   + Q_relax_Stefan*dlnQ_relax_Stefan + Q_molec*dlnQ_molec + Q_shock*dlnQ_shock)/Q_cgs
 endif
 !limit exponent to prevent overflow
 dlnQ_dlnT = sign(min(50.,abs(dlnQ_dlnT)),dlnQ_dlnT)
 Q         = Q_cgs/(unit_ergg/utime)

 !call testfunc()
 !call exit

end subroutine calc_cooling_rate

!-----------------------------------------------------------------------
!+
!  UTILITY: Total cooling function
!+
!-----------------------------------------------------------------------
real function calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 use cooling_functions
 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real, intent(in)  :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL

 calc_Q =  cool_dust_discrete_contact(T_gas, rho_gas, mu, T_dust, d2g, a, rho_grain, kappa_dust) &
!     + cool_dust_full_contact(T_gas, rho_gas, mu, T_dust, kappa_dust) &
!     + cool_dust_radiation(T_gas, kappa_gas, T_dust, kappa_dust) &
    + cool_coulomb(T_gas, rho_gas, mu, nH, nHe) &
    + cool_HI(T_gas, rho_gas, mu, nH, nHe) &
    + cool_H_ionisation(T_gas, rho_gas, mu, nH, nHe) &
    + cool_He_ionisation(T_gas, rho_gas, mu, nH, nHe) &
    + cool_H2_rovib(T_gas, nH, nH2) &
    + cool_H2_dissociation(T_gas, rho_gas, mu, nH, nH2) &
    + cool_CO_rovib(T_gas, rho_gas, mu, nH, nH2, nCO) &
    + cool_H2O_rovib(T_gas, rho_gas, mu, nH, nH2, nH2O) &
    + cool_OH_rot(T_gas, rho_gas, mu, nOH) &
    - heat_dust_friction(rho_gas, v_drift, d2g, a, rho_grain, kappa_dust) &
    - heat_dust_photovoltaic_soft(T_gas, rho_gas, mu, nH, nHe, kappa_dust) &
!     - heat_dust_photovoltaic_hard(T_gas, nH, d2g, kappa_dust, JL) &
    - heat_CosmicRays(nH, nH2)
!     - heat_H2_recombination(T_gas, rho_gas, mu, nH, nH2, T_dust)

!  call print_cooling_rates(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
!                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

end function calc_Q

!-----------------------------------------------------------------------
!+
!  UTILITY: numerical estimate of dlnQ/dlnT
!+
!-----------------------------------------------------------------------
real function calc_dlnQdlnT(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                            T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

 use timestep, only:bignumber

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real, intent(in)  :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust
 real, intent(in)  :: JL

 real, parameter    :: tolQ    = 1.d-4
 real               :: Qtot, dlnQ_dlnT, dT, Q1, Q2, dQdT
 integer, parameter :: itermax = 20
 integer            :: iter

 Qtot      = calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                    T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 dlnQ_dlnT = 0.

! centered finite order approximation for the numerical derivative
 dT = T_gas/100.
 Q1 = calc_Q(T_gas+dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
              T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 Q2 = calc_Q(T_gas-dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
              T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

 iter    = 0
 do while ( abs(Q1/Qtot-1.) > tolQ .and. abs(Q2/Qtot-1.) > tolQ .and. iter < itermax )
    dT   = dT/2.
    Q1   = calc_Q(T_gas+dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                  T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
    Q2   = calc_Q(T_gas-dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                  T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
    iter = iter + 1
 enddo

 dQdT      = (Q1-Q2)/(2.*dT)
 dlnQ_dlnT = (T_gas/Qtot)*dQdT

! gradient can become large at discontinuous physical temperature boundaries (see e.g. atomic and chemical cooling)
 if (dlnQ_dlnT > bignumber) then
    dlnQ_dlnT = 0.
 endif

 calc_dlnQdlnT = dlnQ_dlnT

end function calc_dlnQdlnT

!-----------------------------------------------------------------------
!+
!  Set Temperature grid for exact cooling: between 10K and Tref
!+
!-----------------------------------------------------------------------
subroutine set_Tgrid
 integer :: i
 real    :: dlnT

 dlnT = log(Tref/10.)/(nTg-1)
 do i = 1,nTg
    Tgrid(i) = 10.*exp((i-1)*dlnT)
    !print *,i,Tgrid(i)
 enddo

end subroutine set_Tgrid

SUBROUTINE set_Tgrid_cooltable_Chris
 INTEGER :: i,ierr,k
 REAL :: tol=1.d-12

 !WRITE(*,*) 'k, T(k), Lambda(k) in cgs'
 i = 1
 OPEN(UNIT=15,FILE='cooltable.dat',FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
 IF(ierr/=0) THEN
    WRITE(*,*)
    WRITE(*,*) 'ERROR ERROR ERROR'
    WRITE(*,*) 'cooltable.dat is missing'
    WRITE(*,*) 'Stopping...'
    STOP
 ENDIF
 READ(15,*,IOSTAT=ierr) Tgrid(i), LambdaTable(i)
 IF(ierr/=0) THEN
    WRITE(*,*)
    WRITE(*,*) 'ERROR ERROR ERROR'
    WRITE(*,*) 'cooltable.dat is not properly formatted -- error with the first entry'
    WRITE(*,*) 'Stopping...'
    STOP
 ENDIF
 WRITE(*,*) i, Tgrid(i), LambdaTable(i)
 DO WHILE(ierr==0 .AND. i<nTg)
    i = i+1
    READ(15,*,IOSTAT=ierr) Tgrid(i), LambdaTable(i)
    IF(ierr==0) WRITE(*,*) i, Tgrid(i), LambdaTable(i)
 ENDDO
 CLOSE(15)
 IF(ierr.NE.0) i = i-1
 IF(i/=nTg) THEN
    WRITE(*,*)
    WRITE(*,*) 'ERROR ERROR ERROR'
    IF(i<nTg) THEN
       WRITE(*,*) 'cooltable.dat is not properly formatted -- not enough entries'
    ELSE
       WRITE(*,*) 'cooltable.dat is not properly formatted -- too many entries'
    ENDIF
    WRITE(*,*) 'Stopping...'
    STOP
 ENDIF
 Tref_Chris = Tgrid(i)
 WRITE(*,*) 'set_Tgrid_cooltable_Chris: read in ',i,' temperatures and cooling values'
 WRITE(*,*) 'query: ',nTg,Tref_Chris,LambdaTable(nTg),Tref_Chris/LambdaTable(nTg),i

 DO k=1,i-1
    alphaTable(k) = LOG10(LambdaTable(k+1)/LambdaTable(k)) / LOG10(Tgrid(k+1)/Tgrid(k))
 ENDDO
 !Decision point: how to treat cooling for particles above the cooling curve
 !Option 1: continue the cooling curve at the same power law that connects the final two entries in the table
 alphaTable(i) = alphaTable(i-1)
 !Option 2: make the cooling constant, equal to the final value of the cooling table
 alphaTable(i) = 0.

 YkTable(i) = 0.
 DO k=i-1,1,-1
    IF(ABS(alphaTable(k)-1.) < tol) THEN
       YkTable(k) = YkTable(k+1) - LambdaTable(i)/Tgrid(i) &
                                 * Tgrid(k)/LambdaTable(k) * LOG(Tgrid(k)/Tgrid(k+1))
    ELSE
       YkTable(k) = YkTable(k+1) - LambdaTable(i)/Tgrid(i) &
                                 * Tgrid(k)/LambdaTable(k) / (1.-alphaTable(k)) &
                                 * (1. - (Tgrid(k)/Tgrid(k+1))**(alphaTable(k)-1.))
    ENDIF
 ENDDO

 WRITE(*,*) 'k, alpha(k), Y_k(k)'
 DO k=1,i
    WRITE(*,*) k,alphaTable(k),YkTable(k)
 ENDDO
END SUBROUTINE set_Tgrid_cooltable_Chris

SUBROUTINE set_Tgrid_cooltable_Chris2
 INTEGER :: k,ierr
 REAL :: tol=1.d-12,dummy1(10000),dummy2(10000)

 WRITE(*,*) 'k, T(k), Lambda(k) in cgs'
 k = 1
 OPEN(UNIT=15,FILE='cooltable.dat',FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
 IF(ierr/=0) THEN
    WRITE(*,*)
    WRITE(*,*) 'ERROR ERROR ERROR'
    WRITE(*,*) 'cooltable.dat is missing'
    WRITE(*,*) 'Stopping...'
    STOP
 ENDIF
 READ(15,*,IOSTAT=ierr) dummy1(k), dummy2(k)
 IF(ierr/=0) THEN
    WRITE(*,*)
    WRITE(*,*) 'ERROR ERROR ERROR'
    WRITE(*,*) 'cooltable.dat is not properly formatted -- error with the first entry'
    WRITE(*,*) 'Stopping...'
    STOP
 ENDIF
 DO WHILE(ierr==0) ! .AND. k<nTg)
    !WRITE(*,*) k, dummy1(k), dummy2(k) !now commented out -- cooling table is outputted later
    k = k+1
    IF(k>10000) THEN
       WRITE(*,*) 'ERROR -- increase size of dummy1 and dummy2 in set_Tgrid_cooltable_Chris2'
       STOP
    ENDIF
    READ(15,*,IOSTAT=ierr) dummy1(k), dummy2(k)
 ENDDO
 CLOSE(15)
 IF(ierr.NE.0) k = k-1
 IF(k<2) THEN
    WRITE(*,*)
    WRITE(*,*) 'ERROR ERROR ERROR'
    WRITE(*,*) 'cooltable.dat is not properly formatted -- not enough entries'
    WRITE(*,*) 'Stopping...'
    STOP
 ENDIF
 nTg_Chris = k
 
 IF(ALLOCATED(Tgrid_Chris         )) DEALLOCATE(Tgrid_Chris         )
 ALLOCATE(    Tgrid_Chris(        nTg_Chris))
 IF(ALLOCATED(LambdaTable_Chris   )) DEALLOCATE(LambdaTable_Chris   )
 ALLOCATE(    LambdaTable_Chris(  nTg_Chris))
 IF(ALLOCATED(alphaTable_Chris    )) DEALLOCATE(alphaTable_Chris    )
 ALLOCATE(    alphaTable_Chris(   nTg_Chris))
 IF(ALLOCATED(YkTable_Chris       )) DEALLOCATE(YkTable_Chris       )
 ALLOCATE(    YkTable_Chris(      nTg_Chris))
 IF(ALLOCATED(YTsegAlpha1_Chris   )) DEALLOCATE(YTsegAlpha1_Chris   )
 ALLOCATE(    YTsegAlpha1_Chris(  nTg_Chris))
 IF(ALLOCATED(YTsegAlphaNot1_Chris)) DEALLOCATE(YTsegAlphaNot1_Chris)
 ALLOCATE(    YTsegAlphaNot1_Chris(nTg_Chris))
 
 !cooling table -- put dummy variables into allocated arrays
 DO k=1,nTg_Chris
    Tgrid_Chris(k) = dummy1(k)
    LambdaTable_Chris(k) = dummy2(k)
 ENDDO
 
 !reference temperature T_ref=T_N, which is the final entry in the cooling table
 Tref_Chris = Tgrid_Chris(nTg_Chris)
 WRITE(*,*) 'set_Tgrid_cooltable_Chris2: read in ',nTg_Chris,' temperatures and cooling values'
 
 !frequently needed quantity -- precompute for optimization
 TNdivLN = Tgrid_Chris(nTg_Chris) / LambdaTable_Chris(nTg_Chris)
 !WRITE(*,*) 'query: ',nTg_Chris,Tref_Chris,LambdaTable_Chris(nTg_Chris),TNdivLN
 WRITE(*,*) 'nTg_Chris                             : ',nTg_Chris
 WRITE(*,*) 'Tref_Chris [= Tgrid_Chris(nTg_Chris)] : ',Tref_Chris 
 WRITE(*,*) 'LambdaTableref_Chris                  : ',LambdaTable_Chris(nTg_Chris)
 WRITE(*,*) 'TNdivLN [= Tref / Lambdaref]          : ',TNdivLN

 !Eq. A4, piecewise power law
 DO k=1,nTg_Chris-1
    alphaTable_Chris(k) = LOG10(LambdaTable_Chris(k+1)/LambdaTable_Chris(k)) / LOG10(Tgrid_Chris(k+1)/Tgrid_Chris(k))
 ENDDO
 !Decision point: how to treat cooling for particles above the cooling curve
 !   Option 1: continue the cooling curve at the same power law that connects the final two entries in the table
 alphaTable_Chris(nTg_Chris) = alphaTable_Chris(nTg_Chris-1)
 !   Option 2: make the cooling constant, equal to the final value of the cooling table
 !alphaTable_Chris(nTg_Chris) = 0.

 !Eq. A6
 YkTable_Chris(nTg_Chris) = 0.
 DO k=nTg_Chris-1,1,-1
    IF(ABS(alphaTable_Chris(k)-1.) < tol) THEN
       !YkTable_Chris(k) = YkTable_Chris(k+1) - LambdaTable_Chris(nTg_Chris)/Tgrid_Chris(nTg_Chris) &
       !                          * Tgrid_Chris(k)/LambdaTable_Chris(k) * LOG(Tgrid_Chris(k)/Tgrid_Chris(k+1))
       YkTable_Chris(k) = YkTable_Chris(k+1) - 1./TNdivLN &
                                 * Tgrid_Chris(k)/LambdaTable_Chris(k) * LOG(Tgrid_Chris(k)/Tgrid_Chris(k+1))
    ELSE
       !YkTable_Chris(k) = YkTable_Chris(k+1) - LambdaTable_Chris(nTg_Chris)/Tgrid_Chris(nTg_Chris) &
       !                          * Tgrid_Chris(k)/LambdaTable_Chris(k) / (1.-alphaTable_Chris(k)) &
       !                          * (1. - (Tgrid_Chris(k)/Tgrid_Chris(k+1))**(alphaTable_Chris(k)-1.))
       YkTable_Chris(k) = YkTable_Chris(k+1) - 1./TNdivLN &
                                 * Tgrid_Chris(k)/LambdaTable_Chris(k) / (1.-alphaTable_Chris(k)) &
                                 * (1. - (Tgrid_Chris(k)/Tgrid_Chris(k+1))**(alphaTable_Chris(k)-1.))
    ENDIF
 ENDDO

 !Precomputed combination of Y-related constants for optimization (from Eq. A5 in Townsend+09)        
 !YTsegAlpha1_Chris: for alpha=1
 !YTsegAlphaNot1_Chris: for alpha/=1
 DO k=1,nTg_Chris
    YTsegAlpha1_Chris(k) = Tgrid_Chris(k)/LambdaTable_Chris(k) / TNdivLN
    YTsegAlphaNot1_Chris(k) = YTsegAlpha1_Chris(k)/(1.-alphaTable_Chris(k))
 ENDDO

 WRITE(*,*) 'full cooling table and computed quantities:'
 WRITE(*,*) '          k   T_k (K)                Lambda_k (erg*cm^3/s)    alpha_k                 Y_k'
 DO k=1,nTg_Chris
    WRITE(*,*) k,Tgrid_Chris(k),LambdaTable_Chris(k),alphaTable_Chris(k),YkTable_Chris(k)
 ENDDO
END SUBROUTINE set_Tgrid_cooltable_Chris2

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling_solver(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 !call write_options_molecularcooling(iunit)
 call write_inopt(icool_method,'icool_method',&
                 'integration method (0=implicit, 1=explicit, 2=exact solution)',iunit)
 call write_inopt(excitation_HI,'excitation_HI','cooling via electron excitation of HI (1=on/0=off)',iunit)
 call write_inopt(relax_bowen,'relax_bowen','Bowen (diffusive) relaxation (1=on/0=off)',iunit)
 call write_inopt(relax_stefan,'relax_stefan','radiative relaxation (1=on/0=off)',iunit)
 call write_inopt(dust_collision,'dust_collision','dust collision (1=on/0=off)',iunit)
 call write_inopt(shock_problem,'shock_problem','piecewise formulation for analytic shock solution (1=on/0=off)',iunit)
 if (shock_problem == 1) then
    call write_inopt(lambda_shock_cgs,'lambda_shock','Cooling rate parameter for analytic shock solution',iunit)
    call write_inopt(T1_factor,'T1_factor','factor by which T0 is increased (T1= T1_factor*T0)',iunit)
    call write_inopt(T0_value,'T0','temperature to cool towards (do not modify! set by setup)',iunit)
 endif
 call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)

end subroutine write_options_cooling_solver

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling_solver(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 integer :: nn

 imatch        = .true.
 igotall       = .false.  ! cooling options are compulsory
 select case(trim(name))
 case('icool_method')
    read(valstring,*,iostat=ierr) icool_method
    ngot = ngot + 1
 case('excitation_HI')
    read(valstring,*,iostat=ierr) excitation_HI
    ngot = ngot + 1
 case('relax_bowen')
    read(valstring,*,iostat=ierr) relax_bowen
    ngot = ngot + 1
 case('relax_stefan')
    read(valstring,*,iostat=ierr) relax_stefan
    ngot = ngot + 1
 case('dust_collision')
    read(valstring,*,iostat=ierr) dust_collision
    ngot = ngot + 1
 case('shock_problem')
    read(valstring,*,iostat=ierr) shock_problem
    ngot = ngot + 1
 case('lambda_shock')
    read(valstring,*,iostat=ierr) lambda_shock_cgs
    ngot = ngot + 1
 case('T1_factor')
    read(valstring,*,iostat=ierr) T1_factor
    ngot = ngot + 1
 case('T0')
    read(valstring,*,iostat=ierr) T0_value
    ngot = ngot + 1
 case('bowen_Cprime')
    read(valstring,*,iostat=ierr) bowen_Cprime
    ngot = ngot + 1
 case default
    imatch = .false.
    ierr = 0
 end select
 if (shock_problem == 1) then
    nn = 10
 else
    nn = 7
 endif
 if (ngot >= nn) igotall = .true.

end subroutine read_options_cooling_solver

!=======================================================================
!=======================================================================
!=======================================================================
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!
!  Test routine for cooling functions
!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!=======================================================================
!=======================================================================
!=======================================================================

subroutine testfunc()

 use physcon, only: mass_proton_cgs

 real :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL
 real :: n_gas

 ! evaluate parameters using plausible values to test the cooling functions

 T_gas      = 4000.
 rho_gas    = 2.d-16
 mu         = 1.78

 n_gas      = rho_gas/(mu*mass_proton_cgs)

 nH         = 0.5 *n_gas
 nH2        = 0.5 *n_gas
 nHe        = 0.1 *n_gas
 nCO        = 1.d-4*n_gas
 nH2O       = 5.d-5*n_gas
 nOH        = 1.d-7*n_gas
 kappa_gas  = 2.d-4

 T_dust     = 1500.
 v_drift    = 1.d6
 d2g        = 1./200.
 a          = 1.d-5
 rho_grain  = 2.
 kappa_dust = 1.d-4

 JL = 2.5d-12     ! Value taken from Barstow et al. 1997
 call print_cooling_rates(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

end subroutine testfunc

subroutine print_cooling_rates(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 use cooling_functions
 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real, intent(in)  :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust
 real, intent(in)  :: JL
 real :: Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15, Q16, Q17, Qtot, dlnQ_dlnT, nH_tot

 !nH_tot = nH+2.*nH2
 nH_tot = 1.
 print*, ' '
 print*, ' '
 print*, '-----------------------------------------------------------------------------'
 print*, ' '
 print*, ' '
 Q1  = cool_dust_discrete_contact(T_gas, rho_gas, mu, T_dust, d2g, a, rho_grain, kappa_dust)
 print*, 'Q1   = ', Q1/nH_tot
 Q2  = cool_dust_full_contact(T_gas, rho_gas, mu, T_dust, kappa_dust)
 print*, 'Q2   = ', Q2/nH_tot
 Q3  = cool_dust_radiation(T_gas, kappa_gas, T_dust, kappa_dust)
 print*, 'Q3   = ', Q3/nH_tot
 Q4  = cool_coulomb(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q4   = ', Q4/nH_tot
 Q5  = cool_HI(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q5   = ', Q5/nH_tot
 Q6  = cool_H_ionisation(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q6   = ', Q6/nH_tot
 Q7  = cool_He_ionisation(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q7   = ', Q7/nH_tot
 Q8  = cool_H2_rovib(T_gas, nH, nH2)
 print*, 'Q8   = ', Q8/nH_tot
 Q9  = cool_H2_dissociation(T_gas, rho_gas, mu, nH, nH2)
 print*, 'Q9   = ', Q9/nH_tot
 Q10 = cool_CO_rovib(T_gas, rho_gas, mu, nH, nH2, nCO)
 print*, 'Q10  = ', Q10/nH_tot
 Q11 = cool_H2O_rovib(T_gas, rho_gas, mu, nH, nH2, nH2O)
 print*, 'Q11  = ', Q11/nH_tot
 Q12 = cool_OH_rot(T_gas, rho_gas, mu, nOH)
 print*, 'Q12  = ', Q12/nH_tot
 Q13 = heat_dust_friction(rho_gas, v_drift, d2g, a, rho_grain, kappa_dust)
 print*, 'Q13  = ', Q13/nH_tot
 Q14 = heat_dust_photovoltaic_soft(T_gas, rho_gas, mu, nH, nHe, kappa_dust)
 print*, 'Q14  = ', Q14/nH_tot
 Q15 = heat_dust_photovoltaic_hard(T_gas, nH, d2g, kappa_dust, JL)
 print*, 'Q15  = ', Q15/nH_tot
 Q16 = heat_CosmicRays(nH, nH2)
 print*, 'Q16  = ', Q16/nH_tot
 Q17 = heat_H2_recombination(T_gas, rho_gas, mu, nH, nH2, T_dust)
 print*, 'Q17  = ', Q17/nH_tot

 Qtot = calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 print*, 'Qtot = ', Qtot/nH_tot

 dlnQ_dlnT = calc_dlnQdlnT(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                            T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 print*, 'dlnQdlnT = ', dlnQ_dlnT

 print*, ' '
 print*, ' '
 print*, '------------------- exit in calc_cooling_rate --------------------------------'
 print*, ' '
 print*, ' '

end subroutine print_cooling_rates

end module cooling_solver
