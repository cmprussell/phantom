!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module extern_windaccel
!
! Implementation of external forces from wind-acceleration model
!
! :References: Madura et al. (2013), 436, 4, 3820, Appendix A1
!
! :Owner: Christopher Russell
!
! :Runtime parameters:
!   - kappa_follows_gas    : *kappa follows gas or radiation field*
!   - radiative_inhibition : *radiative inhibition on or off*
!
! :Dependencies: infile_utils, io, kernel, part
!
 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !
 logical, public :: kappa_follows_gas=.true.,radiative_inhibition=.true.

 public :: get_windaccel_force
 public :: write_options_windaccel, read_options_windaccel
 private

contains

!------------------------------------------------
!+
!  Compute wind acceleration where the wind follows a
!  beta velocity law in the absence of companion stars
!+
!------------------------------------------------
subroutine get_windaccel_force(xi,yi,zi,hi,fxi,fyi,fzi,phi,iwindorigj)!,jj)
 use part,   only:xyzmh_ptmass,nptmass,ihsoft,ihacc
 use part,   only:windaccel,ixantgrav,ikappac,ikappar,ivbeta
 use kernel, only:kernel_softening,radkern

 real, intent(in)    :: xi,yi,zi,hi
 real, intent(out)   :: fxi,fyi,fzi,phi
 integer, intent(in) :: iwindorigj!,jj

 real :: ddinv,xantigr,kappa_windaccel,kappa_windaccel_phi
 integer :: iorder(2),k

 real    :: dx,dy,dz,rr2,ddr,dr3,f1,pmassj
 integer :: j,nptmass_max

 !Determine the order that ptmass computations are performed, which is the order that kappa computations are performed.
 !   This is needed when kappa is tied to gas since the star iwindorigj needs to got first when computing kappa_windaccel.
 !   This is not needed when kappa is tied to the radiation field, in which both kappa(1) and kappa(2) need to be computed.
 if (iwindorigj==1) then
    iorder(1) = 1
    iorder(2) = 2
 elseif (iwindorigj==2) then
    iorder(1) = 2
    iorder(2) = 1
 elseif (iwindorigj==0) then
    !choose the first point mass to be the default case for initialization particles
    iorder(1) = 1
    iorder(2) = 2
 else
    print*, 'Fatal Error --> extern_windaccel.f90 --> get_windaccel_force() --> iwindorigj =',iwindorigj, &
            iwindorigj==1, iwindorigj==2, iwindorigj==0 !,jj
    print*, kind(iwindorigj),kind(2),kind(iorder)
    stop
 endif
 !this case statement was made when the above nested if statements were behaving unpredictably
 !select case(iwindorigj)
 !   case(1)
 !      iorder(1) = 1
 !      iorder(2) = 2
 !   case(2)
 !      iorder(1) = 2
 !      iorder(2) = 1
 !   case(0)
 !      !choose the first point mass to be the default case for initialization particles
 !      iorder(1) = 1
 !      iorder(2) = 2
 !   case default
 !      print*, 'Fatal Error --> extern_windaccel.f90 --> get_windaccel_force() --> iwindorigj =',iwindorigj, &
 !              iwindorigj==1, iwindorigj==2, iwindorigj==0 !,jj
 !      print*, kind(iwindorigj),kind(2),kind(iorder)
 !      stop
 !end select
 !print*,'iwindorigj =',iwindorigj,', iorder = ',iorder

 !only loop over one star if radiative inhibition is turned off; loop over all stars if radiative inhibition is turned on
 if (radiative_inhibition) then
    nptmass_max = nptmass
 else
    nptmass_max = 1
 endif

 fxi = 0.
 fyi = 0.
 fzi = 0.
 phi = 0.

 do k=1,nptmass_max
    j = iorder(k)
    !if (j == iwindorigj) cycle !used for gravity/gravity-canceling tests
    pmassj = xyzmh_ptmass(4,j)
    if (pmassj < 0.0) cycle

    dx = xi - xyzmh_ptmass(1,j)
    dy = yi - xyzmh_ptmass(2,j)
    dz = zi - xyzmh_ptmass(3,j)

    rr2 = dx*dx + dy*dy + dz*dz + epsilon(rr2)
    ddr = 1./sqrt(rr2)
    dr3 = ddr*ddr*ddr

    ! acceleration of gas due to point mass particle
    !This is the gravity component, which is how the wind acceleration is scaled.
    f1 = pmassj*dr3

    !wind acceleration
    !!!!!if (k==1) then !use this if statement for no radiative inhibition; comment out this if statement for using radiative inhibition
    !if (k==1) then !use this if statement for radiative inhibition where kappa follows the gas
    !   !kappa follows gas
    !   ddinv = xyzmh_ptmass(ihacc,iwindorigj) * ddr ! Radius/distance = R/r
    !   !kappa follows radiation field
    !   !ddinv = xyzmh_ptmass(ihacc,j) * ddr ! Radius/distance = R/r
    !   if (ddinv<1.d0) then
    !      !kappa follows gas
    !      kappa_windaccel = windaccel(ikappac,iwindorigj) + &
    !                        windaccel(ikappar,iwindorigj)*(1.0d0-ddinv)**(2.d0*windaccel(ivbeta,iwindorigj)-1.d0)
    !      !kappa follows radiation field
    !      !kappa_windaccel = windaccel(ikappac,j) + &
    !      !                  windaccel(ikappar,j)*(1.0d0-ddinv)**(2.d0*windaccel(ivbeta,j)-1.d0)
    !      !elseif (ddinv>1.0d0) then !use for debugging of exact location of particle; not needed for regular use
    !      !kappa_windaccel = 0.d0
    !      !print*,'particle originating from star ',iwindorigj,' is within the radius of star ',j,': ',xyzmh_ptmass(ihacc,iwindorigj),rr2
    !   else
    !      kappa_windaccel = 0.d0
    !      !print*,'particle originating from star ',iwindorigj,' is exactly at the radius of star ',j,': ',xyzmh_ptmass(ihacc,iwindorigj),rr2
    !   endif
    !endif
    !Compute kappa_windaccel for the first star or if kappa follows the radiation field.
    !If kappa follows the gas, then kappa_windaccel from k=1 is used for k>1 as well,
    !   so kappa_windaccel does not need to be recomputed.
    if (k==1 .or. (.not.kappa_follows_gas)) then
       ddinv = xyzmh_ptmass(ihacc,j) * ddr ! radius/distance = R/r
       if (ddinv<1.d0) then
          !force component, used for all coordinate directions
          kappa_windaccel = windaccel(ikappac,j) + &
                            windaccel(ikappar,j)*(1.0d0-ddinv)**(2.d0*windaccel(ivbeta,j)-1.d0)

          !potential component, equivalent method to force component
          !The expression for general beta involves a hypergeometric function, which is not included here.
          !   Instead, the following specific beta values have their specific formula entered below.
          !   For now, uncomment the formula corresponding to the appropriate beta value.  If multiple beta
          !   values are in use, add an if statement to distinguish between the different potential formulas.
          !for beta = 1
          kappa_windaccel_phi = windaccel(ikappac,j)+windaccel(ikappar,j)*(1.0d0-0.5d0*ddinv) !beta=1
          !for beta = 0.5
          !kappa_windaccel_phi = windaccel(ikappac,j)+windaccel(ikappar,j) !beta=0.5
          !for beta = 0.8
          !kappa_windaccel_phi = (1.0d0-1.0d0/ddinv)*(windaccel(ikappac,j) + &
          !                      windaccel(ikappar,j)*0.625d0*(1.0d0-ddinv)**0.6d0) !beta=0.8
          !for beta = 1.5
          !kappa_windaccel_phi = windaccel(ikappac,j) - &
          !                      windaccel(ikappar,j)/(3.0d0*ddinv)*(1.0d0-ddinv)**3 !beta=1.5
          !for beta = 2
          !kappa_windaccel_phi = windaccel(ikappac,j) + &
          !                      windaccel(ikappar,j)*(1.0d0-1.5d0*ddinv+ddinv**2-0.25d0*ddinv**3) !beta=2
          !for beta = 2.5
          !kappa_windaccel_phi = windaccel(ikappac,j) - &
          !                      windaccel(ikappar,j)/(5.0d0*ddinv)*(1.0d0-ddinv)**5 !beta=2.5
          !for beta = 3
          !kappa_windaccel_phi = windaccel(ikappac,j) + &
          !                      windaccel(ikappar,j)*(1.0d0-2.5d0*ddinv+10.0d0/3.0d0*ddinv**2-&
          !                                            2.5d0*ddinv**3+ddinv**4-1.0d0/6.0d0*ddinv**5) !beta=3
          !elseif (ddinv>1.0d0) then !use for debugging of exact location of particle; not needed for regular use
          !kappa_windaccel = 0.d0
          !kappa_windaccel_phi = 0.d0
          !print*,'particle originating from star ',iwindorigj,' is within the radius of star ',j,': ',xyzmh_ptmass(ihacc,j),rr2
       else
          kappa_windaccel = 0.d0
          kappa_windaccel_phi = 0.d0
          !print*,'particle originating from star ',iwindorigj,' is exactly at the radius of star ',j,': ',xyzmh_ptmass(ihacc,j),rr2
       endif
    endif
    xantigr = windaccel(ixantgrav,j)*kappa_windaccel

    fxi = fxi + xantigr*dx*f1
    fyi = fyi + xantigr*dy*f1
    fzi = fzi + xantigr*dz*f1
    phi = phi + windaccel(ixantgrav,j)*kappa_windaccel_phi*pmassj*ddr

    !if (jj>33552) then
    !   print*, jj,kappa_windaccel,windaccel(ixantgrav,j),xantigr,xantigr*dx*f1,dx*f1
    !endif
    !if (mod(jj,1000)==0) then
    !   print*, jj,k,j,kappa_windaccel,xantigr*dx*f1,ftmpxi
    !endif
    !if (mod(jj,100)==0 .and. -0.01<yi .and. yi<0.01 .and. -0.01<zi .and. zi<0.01) then
    !   !print*, jj,k,j,xi,xyzmh_ptmass(1,j),dx,f1,kappa_windaccel,windaccel(ixantgrav,j),xantigr*dx*f1,ftmpxi
    !   print*, jj,k,j,xi,xyzmh_ptmass(1,j),dx,f1,kappa_windaccel,windaccel(ixantgrav,j),xantigr*dx*f1,fxi,phi
    !endif

    !!keep this old method around to see how to mix opacities for the kappa_bar implementation via Eq. A6
    !!--- if kappa_bar is implemented, this would need to be pulled out of the current loop over point masses
    !!       since the kappa for each particle interacting with each point mass would need to be computed
    !!--- additionally, the density contribution from each star's wind would need to be computed, so the
    !!       density loop would need to be modified as well
    !!--- this code snippet was in the density loop
    !do j=1,nptmass
    !   if (akappar(iptm) /= 0.0) then
    !      ddinv = rptmas(j) / sqrt((xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2)
    !      akapp(j) = akappac(j) + akappar(j)*(1.0d0-ddinv)**(2.0d0*vbeta(j)-1.0d0)
    !   else
    !      akapp(j) = akappac(j)
    !   endif
    !enddo
    !!if (iexf==7) then
    !akappa(ipart) = akapp(ibelong(ipart))
    !!else
    !!   if (rhoi /= 0.0) then
    !!      akappa(ipart) = akapp(1)+(akapp(2)-akapp(1))*rhoip2/rhoi
    !!   else
    !!      akappa(ipart) = 0.0
    !!   endif
    !!endif
 enddo

end subroutine get_windaccel_force

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_windaccel(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(kappa_follows_gas,'kappa_follows_gas','kappa follows gas or radiation field',iunit)
 call write_inopt(radiative_inhibition,'radiative_inhibition','radiative inhibition on or off',iunit)

end subroutine write_options_windaccel

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_windaccel(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 !character(len=30), parameter :: label = 'read_options_windaccel'

 igotall = .false.
 imatch = .true.
 select case(trim(name))
 case('kappa_follows_gas')
    read(valstring,*,iostat=ierr) kappa_follows_gas
    ngot = ngot + 1
 case('radiative_inhibition')
    read(valstring,*,iostat=ierr) radiative_inhibition
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 2)

end subroutine read_options_windaccel

end module extern_windaccel
