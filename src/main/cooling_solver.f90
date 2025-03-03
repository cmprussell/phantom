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
! :Owner: Christopher Russell
!
! :Runtime parameters:
!   - T1_factor      : *factor by which T0 is increased (T1= T1_factor*T0)*
!   - bowen_Cprime   : *radiative cooling rate (g.s/cm³)*
!   - dust_collision : *dust collision (1=on/0=off)*
!   - excitation_HI  : *cooling via electron excitation of HI (1=on/0=off)*
!   - lambda_shock   : *Cooling rate parameter for analytic shock solution*
!   - lambda_table   : *Lambda(T) cooling table (1=on/0=off)*
!   - relax_bowen    : *Bowen (diffusive) relaxation (1=on/0=off)*
!   - relax_stefan   : *radiative relaxation (1=on/0=off)*
!   - shock_problem  : *piecewise formulation for analytic shock solution (1=on/0=off)*
!
! :Dependencies: cooling_functions, eos, infile_utils, io, part, physcon,
!   timestep, units
!

 use cooling_functions, only:bowen_Cprime,lambda_shock_cgs,T0_value,T1_factor
 use physcon, only:atomic_mass_unit
 implicit none
 character(len=*), parameter :: label = 'cooling_library'
 integer, public :: excitation_HI = 0, relax_Bowen = 0, dust_collision = 0, relax_Stefan = 0, lambda_table = 0, shock_problem = 0
 integer, public :: icool_method  = 0
 !integer, parameter :: nTg  = 64
 integer, parameter :: nTg  = 1000 ! cloudy
 !integer, parameter :: nTg  = 667 ! cloudy where first T bin spans T_floor=1e4K
 !integer, parameter :: nTg  = 201 ! GS07
 !integer :: nTg_Chris
 integer, allocatable :: nTg_Chris(:)
 real,    parameter :: Tref = 1.d7 !higher value of the temperature grid (for exact cooling)
 real :: Tgrid(nTg)
 real :: LambdaTable(nTg),alphaTable(nTg),YkTable(nTg)
 !real :: Tref_Chris,TNdivLN
 real, dimension(:), allocatable :: Tref_Chris,TNdivLN
 !real, dimension(:), allocatable :: Tgrid_Chris,LambdaTable_Chris,alphaTable_Chris,YkTable_Chris
 !real, dimension(:), allocatable :: YTsegAlphaEQ1_Chris,YTsegAlphaNE1_Chris
 real, dimension(:,:), allocatable :: Tgrid_Chris,LambdaTable_Chris,alphaTable_Chris,YkTable_Chris
 real, dimension(:,:), allocatable :: YTsegAlphaEQ1_Chris,YTsegAlphaNE1_Chris

 public :: init_cooling_solver,read_options_cooling_solver,write_options_cooling_solver
 public :: energ_cooling_solver,calc_cooling_rate,calc_Q
 public :: testfunc,print_cooling_rates
 public :: T0_value,lambda_shock_cgs ! expose to cooling module
 logical, public :: Townsend_test = .false. !for analysis_cooling

 private
 real,    parameter :: Tcap = 1.d3 !Townsend cap temperature

 real, parameter :: habund_solar = 0.7
 real, parameter :: mu_e_solar = 2.0*atomic_mass_unit/(1.0+habund_solar)
 real, parameter :: mu_H_solar = atomic_mass_unit/habund_solar
 real :: mu_e_Arr(60),mu_H_Arr(60)
 !logical, parameter :: methodLong=.true. !Q-based method, native to Phantom
 logical, parameter :: methodLong=.false. !Lambda-based method, native to EIS algorithm


contains
!-----------------------------------------------------------------------
!+
!   Initialise cooling functions and check mix of options are sensible
!+
!-----------------------------------------------------------------------
subroutine init_cooling_solver(ierr)
 use io, only:error,fatal
 integer, intent(out) :: ierr

 print*
 ierr = 0
 !you can't have cool_relaxation_Stefan and cool_relaxation_Bowen at the same time
 if (relax_Bowen == 1 .and. relax_Stefan == 1) then
    call error(label,'you can"t have bowen and stefan cooling at the same time')
    ierr = 1
 endif
 !if no cooling flag activated, disable cooling
 if ( (excitation_HI+relax_Bowen+dust_collision+relax_Stefan+lambda_table+shock_problem) == 0) then
    print *,'ERROR: no cooling prescription activated'
    ierr = 2
 endif
 !call set_Tgrid()
 if (lambda_table == 1) then
    if (methodLong) then
       call set_Tgrid_cooltable_Chris()
    else
       call set_Tgrid_cooltable_Chris2()
    endif
    print "(/,a)", ' Important! Lambda(T) cooling curves are parameterized based on the solar value for mu_H=amu/X, so make sure that the'
    print "(a,f11.9)", '   solar hydrogen abundance mass fraction (X) that was used for creating the cooling tables is the following value:'
    print "(a,f11.9)", '   X = habund_solar = ',habund_solar
    print "(a)", ' If this is not the hydrogen mass fraction used to create the cooling tables, update "habund_solar" in cooling_solver.f90.'
    print "(a)", ' Using habund_solar, the solar value for mu_e and mu_H are:'
    print "(a,es21.14,a)", ' mu_e_solar =',mu_e_solar,' g'
    print "(a,es21.14,a)", ' mu_H_solar =',mu_H_solar,' g'
    print*
    !check for an invalid habund_solar -- if invalid, kill the sim after the above text appears in hopes that it will help solve the issue
    if (habund_solar<tiny(habund_solar)) call fatal('cooling_solver','error with habund_solar value -- must be positive',var='habund_solar',val=habund_solar)
 else
    call set_Tgrid()
 endif

end subroutine init_cooling_solver

!-----------------------------------------------------------------------
!+
!   Get right hand side of energy equation for the desired
!   cooling prescription and choice of solver
!+
!-----------------------------------------------------------------------
subroutine energ_cooling_solver(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa,ict)
 real, intent(in)    :: ui,rho,dt                ! in code units
 real, intent(in)    :: Tdust,mu,gamma,K2,kappa  ! in cgs
 integer, intent(in) :: ict                      ! cooling table's integer value
 real, intent(out)   :: dudt                     ! in code units

 if (icool_method == 2) then
    !call exact_cooling   (ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa)
    if (methodLong) then
       call exact_cooling_Chris(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa)
    else
       call exact_cooling_Chris2(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa,ict)
    endif
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

 !print*, ui,T,Temp,dudt,Q
end subroutine exact_cooling

!-----------------------------------------------------------------------
!+
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!+
!-----------------------------------------------------------------------
subroutine exact_cooling_Chris(ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa)

 use physcon, only:Rg,kboltz
 use units,   only:unit_ergg,unit_density,utime
 !use cooling, only:Tfloor

 real, intent(in)  :: ui, rho, dt, Tdust, mu, gamma
 real, intent(in)  :: K2, kappa
 real, intent(out) :: dudt

 real, parameter :: tol = 1.d-12
 real            :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,T_on_u,T_floor,Qi
 integer         :: k

 real :: rhocgsDivmueDivmuH
 real :: tcool,tcoolfac
 real :: Tref_Chris2
 logical :: opt0,opt1,opt2,withCorrection
 real, parameter :: Tfloor=1.d4
 integer :: ksave
 real :: yksave,ysave
 integer :: ict ! cooling table's integer value

 !Option 0: default calculation where ref=N
 opt0=.true.
 opt1=.false.
 !Option 1: new calcualiton where ref=k+1
 !opt0=.false.
 !opt1=.true.
 !Option 2: new calcualiton where ref=k
 !opt0=.false.
 !opt1=.false.
 !opt2=.true.
 !Option 3: new calcualiton where ref=k,Y_k=0
 !opt0=.false.
 !opt1=.false.
 !opt2=.false.

 !compute k' or not -- "withCorrection=.true." means to copmpute k'
 !withCorrection=.false.
 withCorrection=.true.

 ict = 1

 !print*, 'Rg, kB/amu, kB, amu =',Rg,kboltz/atomic_mass_unit,kboltz,atomic_mass_unit
 rhocgsDivmueDivmuH = rho*unit_density/mu_e_solar/mu_H_solar
 !tcoolfac = kboltz*mu_e_solar*mu_H_solar / ((gamma-1.)*rho*unit_density*mu*atomic_mass_unit)
 tcoolfac = Rg*mu_e_solar*mu_H_solar / ((gamma-1.)*rho*unit_density*mu)
 !print*, 'asdf: ',kboltz,mu_e,mu_H,gamma,rho,unit_density,mu,atomic_mass_unit,tcoolfac
 !print*, 'asdf2:',ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa

 if (Townsend_test) then
    T_floor = Tcap
 else
    !T_floor = 10.
    T_floor = Tfloor
 endif
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T      = T_on_u*ui
 !print*, 'asdf3:',T_on_u,gamma,mu,unit_ergg,Rg

 if (T < T_floor) then
    print*, 'T < T_floor'
    Temp = T_floor
 elseif (T > Tref_Chris(ict)) then
    print*, 'T > Tref_Chris'
    call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Tdust, mu, gamma, K2, kappa)
    Temp = T+T_on_u*Q*dt
 else
    !print*, 'else'
    !print*, 'else -- y=', y
    !print*, 'else -- mu=', mu
    !print*, 'else -- habund, mu_e_solar, mu_H_solar=', habund,mu_e_solar, mu_H_solar
    !print*, 'else -- rho, rho*unit_density, rhocgsDivmueDivmuH=',rho,rho*unit_density,rhocgsDivmueDivmuH,LambdaTable(nTg)
    print*, 'else -- opt0 =',opt0,', opt1 =',opt1,', opt2 =',opt2,', withCorrection =',withCorrection
    if (opt0) then
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
             y = y - Qref*Tgrid(k)/(Q*Tref_Chris(ict))*log(Tgrid(k)/Tgrid(k+1))
          else
             y = y - Qref*Tgrid(k)/(Q*Tref_Chris(ict)*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
          endif
       enddo
    elseif (opt1) then !Y_{k+1}=0, T_N=T_{k+1}, Lambda_N=Lambda_{k+1}
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
    elseif (opt2) then !Y_{k+1}=0, T_N=T_k, Lambda_N=Lambda_k
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
    else !Y_{k+1}=0, T_N=T_k, Lambda_N=Lambda_k
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
    endif
    tcool=tcoolfac*T/LambdaTable(k)
    !print*, 'tcool =',tcool,tcool/utime
    print*, 'tcool (s,yr,code) =',tcool,tcool/(365.25*24.*3600.),tcool/utime,k,LambdaTable(k),T,utime,dt,T_floor
    if (opt0) then
       !eqs A5 for Y(T)
       yk = y
       if (abs(dlnQ_dlnT-1.) < tol) then
          y = yk + Qref*Tgrid(k)/(Q*Tref_Chris(ict))*log(Tgrid(k)/T)
       else
          y = yk + Qref*Tgrid(k)/((Q*Tref_Chris(ict))*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1.))
       endif
       !argument of Y^(-1) in eq 26
       dy = -Qref*dt*T_on_u/Tref_Chris(ict)
    elseif (opt1) then !Y_{k+1}=0, T_N=T_{k+1}, Lambda_N=Lambda_{k+1}
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
    else !Y_{k+1}=0, T_N=T_k, Lambda_N=Lambda_k
       !eqs A5 for Y(T)
       yk = y
       if (abs(dlnQ_dlnT-1.) < tol) then
          y = yk + log(Tgrid(k)/T)
       else
          y = yk + 1./(1.-dlnQ_dlnT)*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1.))
       endif
       !argument of Y^(-1) in eq 26
       dy = -Qref*dt*T_on_u/Tref_Chris2
    endif
    ksave = k
    yksave = yk
    ysave = y
    y  = y + dy

    !New Part -- Start
    if (withCorrection) then
       do while(y>yk .AND. k>1)
          k = k-1
          !call calc_cooling_rate(Q, dlnQ_dlnT, rho, Tgrid(k), Tdust, mu, gamma, K2, kappa)
          !Q=LambdaTable(k)*-2.e23
          Q=-rhocgsDivmueDivmuH*LambdaTable(k) * utime/unit_ergg
          dlnQ_dlnT=alphaTable(k)
          !dlnQ_dlnT=0.
          !dlnQ_dlnT = log(Qi/Q)/log(Tgrid(k+1)/Tgrid(k))
          Qi = Q
          if (opt0) then
             ! eqs A6 to get Yk
             if (abs(dlnQ_dlnT-1.) < tol) then
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris(ict))*log(Tgrid(k)/Tgrid(k+1))
             else
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris(ict)*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
             endif
          elseif (opt1) then !Y_{k+1}=0, T_N=T_{k+1}, Lambda_N=Lambda_{k+1}
             ! eqs A6 to get Yk
             if (abs(dlnQ_dlnT-1.) < tol) then
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris2)*log(Tgrid(k)/Tgrid(k+1))
             else
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris2*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
             endif
          else !Y_{k+1}=0, T_N=T_k, Lambda_N=Lambda_k
             ! eqs A6 to get Yk
             if (abs(dlnQ_dlnT-1.) < tol) then
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris2)*log(Tgrid(k)/Tgrid(k+1))
             else
                yk = yk - Qref*Tgrid(k)/(Q*Tref_Chris2*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
             endif
          endif
       enddo
    endif
    !New Part -- End
    !print*, 'asdf',k,y,yk,Tgrid(k)

    if (opt0) then
       !compute Yinv (eqs A7)
       if (abs(dlnQ_dlnT-1.) < tol) then
          Temp = max(Tgrid(k)*exp(-Q*Tref_Chris(ict)*(y-yk)/(Qref*Tgrid(k))),T_floor)
       else
          Yinv = 1.-(1.-dlnQ_dlnT)*Q*Tref_Chris(ict)/(Qref*Tgrid(k))*(y-yk)
          if (Yinv > 0.) then
             !Temp = Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT)))
             Temp = max(Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT))),T_floor)
             !print*, 'hmm, Temp =',Temp,Tgrid(k),k,Yinv
             !print*, 'confused ',Temp<T_floor
          else
             Temp = T_floor
             !print*, 'AtFloor',Temp,Yinv,k
          endif
       endif
    elseif (opt1) then !Y_{k+1}=0, T_N=T_{k+1}, Lambda_N=Lambda_{k+1}
       !compute Yinv (eqs A7)
       if (abs(dlnQ_dlnT-1.) < tol) then
          Temp = max(Tgrid(k)*exp(-Q*Tref_Chris2*(y-yk)/(Qref*Tgrid(k))),T_floor)
       else
          Yinv = 1.-(1.-dlnQ_dlnT)*Q*Tref_Chris2/(Qref*Tgrid(k))*(y-yk)
          if (Yinv > 0.) then
             !Temp = Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT)))
             Temp = max(Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT))),T_floor)
             !print*, 'hmm, Temp =',Temp,Tgrid(k),k,Yinv
             !print*, 'confused ',Temp<T_floor
          else
             Temp = T_floor
             !print*, 'AtFloor',Temp,Yinv,k
          endif
       endif
    else !Y_{k+1}=0, T_N=T_k, Lambda_N=Lambda_k
       !compute Yinv (eqs A7)
       if (abs(dlnQ_dlnT-1.) < tol) then
          Temp = max(Tgrid(k)*exp(-Q*Tref_Chris2*(y-yk)/(Qref*Tgrid(k))),T_floor)
       else
          Yinv = 1.-(1.-dlnQ_dlnT)*Q*Tref_Chris2/(Qref*Tgrid(k))*(y-yk)
          if (Yinv > 0.) then
             !Temp = Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT)))
             Temp = max(Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT))),T_floor)
             !print*, 'hmm, Temp =',Temp,Tgrid(k),k,Yinv
             !print*, 'confused ',Temp<T_floor
          else
             Temp = T_floor
             !print*, 'AtFloor',Temp,Yinv,k
          endif
       endif
    endif
 endif

 if (Temp<T_floor) print*, 'UH-HOH ',Temp,T_floor
 print*, 'T, Temp =',T,Temp
 dudt = (Temp-T)/T_on_u/dt
 !note that u = Temp/T_on_u

 !print "(A,6(1PE14.6))", 'C',ui,T,Temp,dudt,Q,Tref_Chris
 !print*, 'stuff: ',ksave,k,yksave,ysave,tcoolfac,dy,y,tcoolfac*T/(LambdaTable(ksave)*(T/Tgrid(ksave))**alphaTable(ksave)),tcoolfac*T/LambdaTable(ksave)
 !print*, 'ffuts: ',dy,dt,utime,tcoolfac,T_on_u,Tref_Chris,-Qref
end subroutine exact_cooling_Chris

!-----------------------------------------------------------------------
!+
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!+
!-----------------------------------------------------------------------
subroutine exact_cooling_Chris2(ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa, ict)

 use physcon, only:Rg,kboltz
 use units,   only:unit_ergg,unit_density,utime
 !use cooling, only:Tfloor
 use eos,  only:use_var_comp
 use part, only:n_startypes,mu_startypes
 use io,   only:fatal

 real, intent(in)    :: ui, rho, dt, Tdust, mu, gamma
 real, intent(in)    :: K2, kappa
 integer, intent(in) :: ict ! cooling table's integer value
 real, intent(out)   :: dudt

 real, parameter :: tol = 1.d-12
 real            :: Yinv,Temp,T,T_on_u,T_floor
 integer         :: k

 !real :: rhocgsDivmueDivmuH
 real :: tcoolfac
 !real :: Tref_Chris2
 !logical :: opt0,opt1,opt2,withCorrection
 real, parameter :: Tfloor=1.d4

 integer :: kl,km,ku,kk
 real :: YT
 !real :: YTsave
 integer :: ict2
 real :: tolmu=1.d-6

 !determine which cooling table to use for this particular particle
 !for now, this is just a consistency check. Eventually, this mu-dependent ict2 method will no longer be used
 if (use_var_comp) then
    !find the correct cooling table based on mu -- if the particle's mu is not found, default to the first cooling table, i.e. ict=1
    !   Note: It would be faster if iwindorig(i) was available to determine the cooling table, but this variable would need to
    !   be added to the calls of energ_cooling() to make it here, which involves making many modification to the code.
    !   So for now, we will use the following mu-to-mu_startypes comparison to determine which cooling table to use.
    ict2 = n_startypes
    do while (abs(mu-mu_startypes(ict2))>tolmu .and. ict2>1)
       ict2 = ict2-1
    enddo
    !print*, ict2,mu
    if (ict<=0) call fatal ('exact_cooling_Chris2','particle has undefined cooling table -- ict is too low',var='ict',ival=ict)
    if (ict>n_startypes) call fatal ('exact_cooling_Chris2','particle has undefined cooling table -- ict is too high',var='ict',ival=ict)
    if (ict/=ict2) print "(2(a,i0))", 'WEIRD: ict/=ict2 -- ict = ',ict,', ict2 = ',ict2
    !else  !ict=1 is now the defualt option, so this isn't needed
    !   single composition, so use the only cooling table available
    !   ict = 1
 endif

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

 !print*, 'Rg, kB/amu, kB, amu =',Rg,kboltz/atomic_mass_unit,kboltz,atomic_mass_unit
 !rhocgsDivmueDivmuH = rho*unit_density/mu_e_solar/mu_H_solar
 !tcoolfac = kboltz*mu_e_solar*mu_H_solar / ((gamma-1.)*rho*unit_density*mu*atomic_mass_unit) ! = tcool*Lambda(T)/T
 !tcoolfac = Rg*mu_e_solar*mu_H_solar / ((gamma-1.)*rho*unit_density*mu) ! = tcool*Lambda(T)/T
 !if (use_var_comp) then
 tcoolfac = Rg*mu_e_Arr(ict)*mu_H_solar / ((gamma-1.)*rho*unit_density*mu) ! = tcool*Lambda(T)/T
 !else
 !   tcoolfac = Rg*mu_e_solar*mu_H_solar / ((gamma-1.)*rho*unit_density*mu) ! = tcool*Lambda(T)/T
 !endif

 if (Townsend_test) then
    T_floor = Tcap
 else
    !T_floor = 10.
    T_floor = Tfloor
 endif
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg ! = (gamma-1.)*mu*atomic_mass_unit/kB * unit_ergg = (gamma-1.)*mu_cgs/kB * unit_ergg
 T      = T_on_u*ui

 if (T < T_floor) then
    Temp = T_floor
    !elseif (T > Tref_Chris) then
    !   print*, 'T > Tref_Chris'
    !   call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Tdust, mu, gamma, K2, kappa)
    !   Temp = T+T_on_u*Q*dt
    !elseif (.true.) then !these two lines cool all particles to the floor temperature
    !   Temp = T_floor
 else
    !bisector to find k in Eq. A5
    !c----    Locate the interval in which temp is found.
    !c        This part is taken from locate.f in Numerical Recipes.
    !c
    !ignoring the initial conditionals, the results will be kl_init <= k <= ku_init-1
    !   therefore, choose kl_init = 1 and ku_init = nTg_Chris+1
    if (T<=Tgrid_Chris(1,ict)) then
       k = 1
    elseif (T>=Tgrid_Chris(nTg_Chris(ict),ict)) then
       k = nTg_Chris(ict)
    else
       kl = 1
       ku = nTg_Chris(ict)+1
       do while (ku-kl>1)
          km = (ku+kl)/2
          if (T>=Tgrid_Chris(km,ict)) then
             kl = km
          else
             ku = km
          endif
       enddo
       k = kl
    endif
    !print*, 'tcool = ',tcoolfac*T/(LambdaTable_Chris(k)*(T/Tgrid_Chris(k))**alphaTable_Chris(k)), &
    !                      tcoolfac*T/(LambdaTable_Chris(k)*(T/Tgrid_Chris(k))**alphaTable_Chris(k)) / utime, &
    !                      k,nTg_Chris, &
    !                      T,Tgrid_Chris(k), &
    !                      LambdaTable_Chris(k)*(T/Tgrid_Chris(k))**alphaTable_Chris(k),LambdaTable_Chris(k), &
    !                      alphaTable_Chris(k)

    !find Y(T), Eq. A5
    if (abs(alphaTable_Chris(k,ict)-1.) < tol) then
       YT = YkTable_Chris(k,ict) + YTsegAlphaEQ1_Chris(k,ict)*log(Tgrid_Chris(k,ict)/T)
    else
       YT = YkTable_Chris(k,ict) + YTsegAlphaNE1_Chris(k,ict)*(1.-(Tgrid_Chris(k,ict)/T)**(alphaTable_Chris(k,ict)-1.))
    endif
    !YTsave = YT
    !print*, 'YT1 =',YT,YkTable_Chris(k),YTsegAlphaNE1_Chris(k),Tgrid_Chris(k),T,alphaTable_Chris(k), &
    !                   1.-(Tgrid_Chris(k)/T)**(alphaTable_Chris(k)-1.), &
    !                   YTsegAlphaNE1_Chris(k)*(1.-(Tgrid_Chris(k)/T)**(alphaTable_Chris(k)-1.))

    !argument of Eq. 26
    !cCool_v2b   delta_Y = (gamma-1)*rho*mu / (mu_e*mu_H*kB) * LambdaN/TempN * delta_t
    !         dyfunx_v2b = gamma1*dble(rho(ipart)*udens)
    !     &                *dble(amunit*gmw/(amue*amuh*boltz))
    !     &                /TNdivLN_v2b
    !     &                *dble(deltat*utime)
    !yfunx_v2b = yfunx_v2b + dyfunx_v2b
    !dYT = 1./(tcoolfac*TNdivLN) * dt*utime
    !yT = yT + dyT
    YT = YT + 1./(tcoolfac*TNdivLN(ict)) * dt*utime
    !print*, 'YT2 =',YT,YkTable_Chris(1),YkTable_Chris(nTg_Chris), &
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
    if (YT>YkTable_Chris(1,ict)) then
       kk = 0
    elseif (YT==YkTable_Chris(1,ict)) then
       kk = 1
    elseif (YT<=YkTable_Chris(nTg_Chris(ict),ict)) then
       kk = nTg_Chris(ict)
    else
       kl = 1
       ku = nTg_Chris(ict)+1
       do while(ku-kl>1)
          km = (ku+kl)/2
          if (YT<=YkTable_Chris(km,ict)) then
             kl = km
          else
             ku = km
          endif
       enddo
       kk = kl
    endif

    ! find Yinv, Eq. A7
    if (kk==0) then
       Temp = T_floor
    else
       if (abs(alphaTable_Chris(kk,ict)-1.) < tol) then
          Temp = max(Tgrid_Chris(kk,ict)*exp(-(YT-YkTable_Chris(kk,ict))/YTsegAlphaEQ1_Chris(kk,ict)) &
                     ,T_floor)
       else
          Yinv = 1.-(YT-YkTable_Chris(kk,ict))/YTsegAlphaNE1_Chris(kk,ict)
          if (Yinv>0.0) then
             Temp = max(Tgrid_Chris(kk,ict)*Yinv**(1./(1.-alphaTable_Chris(kk,ict))) &
                        ,T_floor)
          else
             Temp = T_floor
          endif
       endif
    endif
 endif

 if (Temp<T_floor) print*, 'UH-HOH ',Temp,T_floor
 dudt = (Temp-T)/T_on_u/dt
 !note that u = Temp/T_on_u

 !print "(a,6(es14.6))", 'C',ui,T,Temp,dudt,Q,Tref_Chris
 !print* 'stuff: ',k,kk,YkTable_Chris(k),YTsave,tcoolfac,1./(tcoolfac*TNdivLN)*dt*utime,YT,tcoolfac*T/(LambdaTable_Chris(k)*(T/Tgrid_Chris(k))**alphaTable_Chris(k))
 !print*, 'ffuts: ',1./(tcoolfac*TNdivLN)*dt*utime,dt,utime,tcoolfac,T_on_u,Tref_Chris,TNdivLN
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

!-----------------------------------------------------------------------
!+
!  Read in cooling table for exact cooling
!+
!-----------------------------------------------------------------------
subroutine set_Tgrid_cooltable_Chris
 use io, only:fatal
 integer :: i,ierr,k
 real    :: tol=1.d-12

 print "(a)", ' k, T(k), Lambda(k) in cgs'
 i = 1
 open(UNIT=15,file='cooltable.dat',form='formatted',status='old',iostat=ierr)
 if (ierr/=0) then
    print "(/,a)", ' ERROR ERROR ERROR'
    print "(a)", ' cooltable.dat is missing'
    call fatal('cooling','cooltable.dat is missing')
 endif
 read(15,*,iostat=ierr) Tgrid(i), LambdaTable(i)
 if (ierr/=0) then
    print "(/,a)", ' ERROR ERROR ERROR'
    print "(a)", ' cooltable.dat is not properly formatted -- error with the first entry'
    call fatal('cooling','cooltable.dat is not properly formatted -- error with the first entry')
 endif
 print*, i, Tgrid(i), LambdaTable(i)
 do while(ierr==0 .AND. i<nTg)
    i = i+1
    read(15,*,iostat=ierr) Tgrid(i), LambdaTable(i)
    if (ierr==0) print "(i6,2(es16.8))", i, Tgrid(i), LambdaTable(i)
 enddo
 close(15)
 if (ierr/=0) i = i-1
 if (i/=nTg) then
    print "(/,a)", ' ERROR ERROR ERROR'
    if (i<nTg) then
       print "(a,i0)", ' cooltable.dat is not properly formatted -- not enough entries for nTg = ',nTg
       call fatal('cooling','cooltable.dat is not properly formatted -- not enough entries',var='i',ival=i)
    else
       print "(a,i0)", ' cooltable.dat is not properly formatted -- too many entries for nTg = ',nTg
       call fatal('cooling','cooltable.dat is not properly formatted -- too many entries',var='i',ival=i)
    endif
 endif
 Tref_Chris = Tgrid(i)
 print "(a,i0,a)", ' set_Tgrid_cooltable_Chris: read in ',i,' temperatures and cooling values'
 print "(a,i0,a,2(es16.8,a),i0)", 'query: ',nTg,' ',Tref_Chris,' ',LambdaTable(nTg),' ',Tref_Chris/LambdaTable(nTg),' ',i

 do k=1,i-1
    alphaTable(k) = log10(LambdaTable(k+1)/LambdaTable(k)) / log10(Tgrid(k+1)/Tgrid(k))
 enddo
 !Decision point: how to treat cooling for particles above the cooling curve
 !Option 1: continue the cooling curve at the same power law that connects the final two entries in the table
 alphaTable(i) = alphaTable(i-1)
 !Option 2: make the cooling constant, equal to the final value of the cooling table
 alphaTable(i) = 0.

 YkTable(i) = 0.
 do k=i-1,1,-1
    if (abs(alphaTable(k)-1.) < tol) then
       YkTable(k) = YkTable(k+1) - LambdaTable(i)/Tgrid(i) &
                                 * Tgrid(k)/LambdaTable(k) * log(Tgrid(k)/Tgrid(k+1))
    else
       YkTable(k) = YkTable(k+1) - LambdaTable(i)/Tgrid(i) &
                                 * Tgrid(k)/LambdaTable(k) / (1.-alphaTable(k)) &
                                 * (1. - (Tgrid(k)/Tgrid(k+1))**(alphaTable(k)-1.))
    endif
 enddo

 print "(a)", ' k, alpha(k), Y_k(k)'
 do k=1,i
    print "(i6,2(es16.8))", k,alphaTable(k),YkTable(k)
 enddo
end subroutine set_Tgrid_cooltable_Chris

!-----------------------------------------------------------------------
!+
!  Read in cooling table(s) for exact cooling
!+
!-----------------------------------------------------------------------
subroutine set_Tgrid_cooltable_Chris2
 use io,   only:fatal,error
 use eos,  only:use_var_comp,num_var_comp
 use part, only:name_startypes,habund_startypes
 integer            :: k,ierr,ict
 integer, parameter :: max_ct_entries=10000
 real               :: tol=1.d-12
 real, allocatable  :: dummy1(:,:),dummy2(:,:)
 integer            :: nCoolTables=0,nTg_Chris_max=0

 nCoolTables = max(num_var_comp,1)
 print "(2(a,i0))", ' set_Tgrid_cooltable_Chris2: num_var_comp = ',num_var_comp,', nCoolTables = ',nCoolTables
 allocate(dummy1(max_ct_entries,nCoolTables))
 allocate(dummy2(max_ct_entries,nCoolTables))
 dummy1 = 0.
 dummy2 = 0.
 allocate(nTg_Chris(nCoolTables)) !currently this variable in not deallocated -- should I just define it to be the same size as "mu_startypes"?
 nTg_Chris = 0 !initialize to be used in max statement for multiple cooling tables

 !read in cooling tables -- first using dummy variables, then in correct-size allocated arrays
 do ict=1,nCoolTables
    k = 1
    if (nCoolTables==1) then
       open(unit=15,file='cooltable.dat',form='formatted',status='old',iostat=ierr)
       if (ierr/=0) then
          print "(/,a)", ' ERROR ERROR ERROR'
          print "(a)", ' cooltable.dat is missing'
          call fatal('cooling','cooltable.dat is missing',var='ict',ival=ict)
       endif
    else
       print "(a,i0,a)", ' opening cooling table ',ict,': cooltable_'//trim(name_startypes(ict))//'.dat'
       open(unit=15,file='cooltable_'//trim(name_startypes(ict))//'.dat',form='formatted',status='old',iostat=ierr)
       if (ierr/=0) then
          print "(/,a)", ' ERROR ERROR ERROR'
          print "(a)", ' cooltable_'//trim(name_startypes(ict))//'.dat is missing'
          call fatal('cooling','cooltable_'//trim(name_startypes(ict))//'.dat is missing',var='ict',ival=ict)
       endif
    endif
    read(15,*,iostat=ierr) dummy1(k,ict), dummy2(k,ict)
    if (ierr/=0) then
       print "(/,a)", ' ERROR ERROR ERROR'
       print "(a)", ' cooltable.dat is not properly formatted -- error with the first entry'
       call fatal('cooling','cooltable_'//trim(name_startypes(ict))//'.dat is not properly formatted -- error with the first entry',var='ict',ival=ict)
       stop
    endif
    do while(ierr==0) ! .AND. k<nTg)
       k = k+1
       if (k>max_ct_entries) then
          print "(/,a)", ' ERROR -- increase size of dummy1 and dummy2 in set_Tgrid_cooltable_Chris2'
          print "(/,a,i0)", ' ERROR -- max_ct_entries = ',max_ct_entries
          call fatal('cooling','increase size of dummy1 and dummy2 in set_Tgrid_cooltable_Chris2',var='ict',ival=ict)
       endif
       read(15,*,iostat=ierr) dummy1(k,ict), dummy2(k,ict)
       if (ierr==0 .and. dummy1(k,ict)<=dummy1(k-1,ict)) then
          !cooling table is not correct
          print "(/,a,i0,a)", ' ERROR -- temperatures for cooling table ',ict,' need to be increasing'
          print "(2(a,es22.14))", ' ERROR -- instead, T(k) = ',dummy1(k,ict),' <= T(k-1) = ',dummy1(k-1,ict)
          print "(a,i0)", ' ERROR -- offending entry is k = ',k
          call fatal('cooling','temperatures for cooling tables need to be increasing',var='ict',ival=ict)
       endif
    enddo
    close(15)
    if (ierr/=0) k = k-1
    if (k<2) then
       print "(/,a)", ' ERROR ERROR ERROR'
       print "(a)", ' cooltable.dat is not properly formatted -- not enough entries'
       call fatal('cooling','cooltable_'//trim(name_startypes(ict))//'.dat is not properly formatted -- not enough entries',var='ict',ival=ict)
       stop
    endif
    nTg_Chris(ict) = k
 enddo

 do ict=1,nCoolTables
    print "(2(a,i0),a)", ' set_Tgrid_cooltable_Chris2: cooling table ',ict,' -- read in ',nTg_Chris(ict),' temperatures and cooling values'
 enddo
 nTg_Chris_max = maxval(nTg_Chris)
 !print "(a,1000(1x,i0))", ' nTg_Chris =',nTg_Chris
 print "(2(a,i0),a)", ' allocated array size for various cooling variables: (nTg_Chris_max, nCoolTables) = (',nTg_Chris_max,', ',nCoolTables,')'

 if (allocated(Tgrid_Chris        )) deallocate(Tgrid_Chris        )
 allocate(     Tgrid_Chris(        nTg_Chris_max,nCoolTables))
 if (allocated(LambdaTable_Chris  )) deallocate(LambdaTable_Chris  )
 allocate(     LambdaTable_Chris(  nTg_Chris_max,nCoolTables))
 if (allocated(alphaTable_Chris   )) deallocate(alphaTable_Chris   )
 allocate(     alphaTable_Chris(   nTg_Chris_max,nCoolTables))
 if (allocated(YkTable_Chris      )) deallocate(YkTable_Chris      )
 allocate(     YkTable_Chris(      nTg_Chris_max,nCoolTables))
 if (allocated(YTsegAlphaEQ1_Chris)) deallocate(YTsegAlphaEQ1_Chris)
 allocate(     YTsegAlphaEQ1_Chris(nTg_Chris_max,nCoolTables))
 if (allocated(YTsegAlphaNE1_Chris)) deallocate(YTsegAlphaNE1_Chris)
 allocate(     YTsegAlphaNE1_Chris(nTg_Chris_max,nCoolTables))

 if (allocated(Tref_Chris)) deallocate(Tref_Chris)
 allocate(     Tref_Chris(nCoolTables))
 if (allocated(TNdivLN   )) deallocate(TNdivLN   )
 allocate(     TNdivLN(   nCoolTables))

 !cooling table -- put dummy variables into allocated arrays
 do ict=1,nCoolTables
    do k=1,nTg_Chris(ict)
       Tgrid_Chris(k,ict) = dummy1(k,ict)
       LambdaTable_Chris(k,ict) = dummy2(k,ict)
    enddo
 enddo

 !reference temperature T_ref=T_N, which is the final entry in the cooling table
 do ict=1,nCoolTables
    Tref_Chris(ict) = Tgrid_Chris(nTg_Chris(ict),ict)
 enddo

 !mean molecular weights for electrons and hydrogen atoms
 if (use_var_comp) then
    do ict=1,nCoolTables
       mu_e_Arr(ict) = 2.0*atomic_mass_unit/(1.0+habund_startypes(ict))
       if (habund_startypes(ict)>tiny(habund_startypes)) then
          mu_H_Arr(ict) = atomic_mass_unit/habund_startypes(ict)
       else
          mu_H_Arr(ict) = 0.
       endif
    enddo
 else
    mu_e_Arr(1) = 2.0*atomic_mass_unit/(1.0+habund_solar)
    if (habund_solar>tiny(habund_startypes)) then
       mu_H_Arr(1) = atomic_mass_unit/habund_solar
    else
       mu_H_Arr(1) = 0
       call error('cooling_solver','incorrect value for habund_solar',var='habund_solar',val=habund_solar)
    endif
 endif

 !frequently needed quantity -- precompute for optimization
 print "(/,a)", ' select quantities for each cooling table:'
 print "(1x,67('-'))"
 do ict=1,nCoolTables
    TNdivLN(ict) = Tgrid_Chris(nTg_Chris(ict),ict) / LambdaTable_Chris(nTg_Chris(ict),ict)
    !print*, 'query: ',nTg_Chris,Tref_Chris,LambdaTable_Chris(nTg_Chris),TNdivLN
    if (nCoolTables==1) then
       print "(a,i0,a)", ' cooling table ',ict,': cooltable.dat'
    else
       print "(a,i0,a)", ' cooling table ',ict,': cooltable_'//trim(name_startypes(ict))//'.dat'
    endif
    print "(a,i0,a,i0)",      ' nTg_Chris(',ict,')                             : ',nTg_Chris(ict)
    print "(a,i0,a,es22.14)", ' Tref_Chris(',ict,') [= Tgrid_Chris(nTg_Chris)] : ',Tref_Chris(ict)
    print "(a,i0,a,es22.14)", ' LambdaTableref_Chris(',ict,')                  : ',LambdaTable_Chris(nTg_Chris(ict),ict)
    print "(a,i0,a,es22.14)", ' TNdivLN(',ict,') [= Tref / Lambdaref]          : ',TNdivLN(ict)
    print "(a,i0,a,es22.14)", ' mu_e_Arr(',ict,') (grams)                      : ',mu_e_Arr(ict)
    print "(a,i0,a,es22.14)", ' mu_H_Arr(',ict,') (grams)                      : ',mu_H_Arr(ict)
    print "(a,i0,a,es22.14)", ' habund_startypes(',ict,') [=X_H] (mass frac)   : ',habund_startypes(ict)
    print "(1x,67('-'))"
 enddo

 !Eq. A4, piecewise power law
 do ict=1,nCoolTables
    do k=1,nTg_Chris(ict)-1
       alphaTable_Chris(k,ict) = log10(LambdaTable_Chris(k+1,ict)/LambdaTable_Chris(k,ict)) / log10(Tgrid_Chris(k+1,ict)/Tgrid_Chris(k,ict))
    enddo
    !Decision point: how to treat cooling for particles above the cooling curve
    !   Option 1: continue the cooling curve at the same power law that connects the final two entries in the table
    alphaTable_Chris(nTg_Chris(ict),ict) = alphaTable_Chris(nTg_Chris(ict)-1,ict)
    !   Option 2: make the cooling constant, equal to the final value of the cooling table
    !alphaTable_Chris(nTg_Chris(ict),ict) = 0.
 enddo

 !Eq. A6
 do ict=1,nCoolTables
    YkTable_Chris(nTg_Chris(ict),ict) = 0.
    do k=nTg_Chris(ict)-1,1,-1
       if (abs(alphaTable_Chris(k,ict)-1.) < tol) then
          !YkTable_Chris(k) = YkTable_Chris(k+1) - LambdaTable_Chris(nTg_Chris)/Tgrid_Chris(nTg_Chris) &
          !                          * Tgrid_Chris(k)/LambdaTable_Chris(k) * log(Tgrid_Chris(k)/Tgrid_Chris(k+1))
          YkTable_Chris(k,ict) = YkTable_Chris(k+1,ict) - 1./TNdivLN(ict) &
                                    * Tgrid_Chris(k,ict)/LambdaTable_Chris(k,ict) * log(Tgrid_Chris(k,ict)/Tgrid_Chris(k+1,ict))
       else
          !YkTable_Chris(k) = YkTable_Chris(k+1) - LambdaTable_Chris(nTg_Chris)/Tgrid_Chris(nTg_Chris) &
          !                          * Tgrid_Chris(k)/LambdaTable_Chris(k) / (1.-alphaTable_Chris(k)) &
          !                          * (1. - (Tgrid_Chris(k)/Tgrid_Chris(k+1))**(alphaTable_Chris(k)-1.))
          YkTable_Chris(k,ict) = YkTable_Chris(k+1,ict) - 1./TNdivLN(ict) &
                                    * Tgrid_Chris(k,ict)/LambdaTable_Chris(k,ict) / (1.-alphaTable_Chris(k,ict)) &
                                    * (1. - (Tgrid_Chris(k,ict)/Tgrid_Chris(k+1,ict))**(alphaTable_Chris(k,ict)-1.))
       endif
    enddo
 enddo

 !Precomputed combination of Y-related constants for optimization (from Eq. A5 in Townsend+09)
 !YTsegAlphaEQ1_Chris: for alpha=1
 !YTsegAlphaNE1_Chris: for alpha/=1
 do ict=1,nCoolTables
    do k=1,nTg_Chris(ict)
       YTsegAlphaEQ1_Chris(k,ict) = Tgrid_Chris(k,ict)/LambdaTable_Chris(k,ict) / TNdivLN(ict)
       YTsegAlphaNE1_Chris(k,ict) = YTsegAlphaEQ1_Chris(k,ict)/(1.-alphaTable_Chris(k,ict))
    enddo
 enddo

 print "(/,a)", ' cooling table(s) and computed quantities:'
 print "(1x,99('-'))"
 do ict=1,nCoolTables
    print "(a,i0)", ' cooling table ',ict
    print "(a)", '     k   T_k (K)                Lambda_k (erg*cm^3/s)  alpha_k                Y_k'
    if (nTg_Chris(ict)>30) then
       do k=1,10
          print "(i6,4(es23.14))", k,Tgrid_Chris(k,ict),LambdaTable_Chris(k,ict),alphaTable_Chris(k,ict),YkTable_Chris(k,ict)
       enddo
       print "(a)", '   ...'
       do k=100,nTg_Chris(ict)-10,100
          print "(i6,4(es23.14))", k,Tgrid_Chris(k,ict),LambdaTable_Chris(k,ict),alphaTable_Chris(k,ict),YkTable_Chris(k,ict)
       enddo
       print "(a)", '   ...'
       do k=nTg_Chris(ict)-9,nTg_Chris(ict)
          print "(i6,4(es23.14))", k,Tgrid_Chris(k,ict),LambdaTable_Chris(k,ict),alphaTable_Chris(k,ict),YkTable_Chris(k,ict)
       enddo
    else
       do k=1,nTg_Chris(ict)
          print "(i6,4(es23.14))", k,Tgrid_Chris(k,ict),LambdaTable_Chris(k,ict),alphaTable_Chris(k,ict),YkTable_Chris(k,ict)
       enddo
    endif
    print "(1x,99('-'))"
 enddo

 if (allocated(dummy1)) deallocate(dummy1)
 if (allocated(dummy2)) deallocate(dummy2)
end subroutine set_Tgrid_cooltable_Chris2

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
 call write_inopt(lambda_table,'lambda_table','Lambda(T) cooling table (1=on/0=off)',iunit)
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
 case('lambda_table')
    read(valstring,*,iostat=ierr) lambda_table
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
