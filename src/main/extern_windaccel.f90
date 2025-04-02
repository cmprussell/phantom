!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
!
! :Dependencies: infile_utils, io
!
 implicit none
 !
 !--code input parameters: these are the default values
 !  and can be changed in the input file
 !

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
 real, intent(inout) :: fxi,fyi,fzi
 real, intent(inout) :: phi
 integer, intent(in) :: iwindorigj!,jj

 real :: ddinv,xantigr,kappa_windaccel 
 integer :: iorder(2),k

 real                             :: ftmpxi,ftmpyi,ftmpzi
 real                             :: dx,dy,dz,rr2,ddr,dr3,f1,f2,pmassj
 real                             :: hsoft,hsoft1,hsoft21,q2i,qi,psoft,fsoft
 real                             :: fxj,fyj,fzj,dsx,dsy,dsz
 integer                          :: j
 logical                          :: tofrom,extrap,kappa

 !
 ! Determine if acceleration is from/to gas, or to gas
 !
 !if (present(pmassi) .and. present(fxyz_ptmass) .and. present(fonrmax)) then
 !   tofrom  = .true.
 !   fonrmax = 0.
 !else
 tofrom  = .false.
 !endif

 !! check if it is a force computed using Omelyan extrapolation method for FSI
 !if (present(extrapfac)) then
 !   extrap = .true.
 !else
 extrap = .false.
 !endif

 !if (present(bin_info)) then
 !   kappa = .true.
 !else
 kappa = .false.
 !endif

 !order in which ptmass computations are performed, which is the order in which kappa computations are performed
 !   needed when kappa is tied to gas, so only kappa(iwindorigi) needs to be computed
 !   not needed when kappa is tied to radiation field, in which kappa(1) and kappa(2) need to be computed
 !if(.TRUE.) then
 !if(.FALSE.) then
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
    print*, KIND(iwindorigj),KIND(2),KIND(iorder)
    stop
 endif
 !else !this case statement was made when the above nested if statements were behaving unpredictably
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
 !      print*, KIND(iwindorigj),KIND(2),KIND(iorder)
 !      stop
 !end select
 !endif
 !print*,'iwindorigj =',iwindorigj,', iorder = ',iorder

 ftmpxi = 0.  ! use temporary summation variable
 ftmpyi = 0.  ! (better for round-off, plus we need this bit of
 ftmpzi = 0.  ! the force to calculate the dtphi timestep)
 phi    = 0.
 f2     = 0.

 !do j=1,nptmass
 do k=1,nptmass
    !if (k>1) cycle !turn off wind accel on gas particles from non-originating (i.e. companion) star
    j = iorder(k)
    !if(j == iwindorigj) cycle !used for gravity/gravity-canceling tests
    !if (extrap) then
    !   dx     = xi - (xyzmh_ptmass(1,j) + extrapfac*fsink_old(1,j))
    !   dy     = yi - (xyzmh_ptmass(2,j) + extrapfac*fsink_old(2,j))
    !   dz     = zi - (xyzmh_ptmass(3,j) + extrapfac*fsink_old(3,j))
    !else
    dx     = xi - xyzmh_ptmass(1,j)
    dy     = yi - xyzmh_ptmass(2,j)
    dz     = zi - xyzmh_ptmass(3,j)
    !endif
    pmassj = xyzmh_ptmass(4,j)
    hsoft  = xyzmh_ptmass(ihsoft,j)
    if (hsoft > 0.0) hsoft = max(hsoft,hi)
    if (pmassj < 0.0) cycle

    rr2    = dx*dx + dy*dy + dz*dz + epsilon(rr2)
!#ifdef FINVSQRT
!    ddr    = finvsqrt(rr2)
!#else
    ddr    = 1./sqrt(rr2)
!#endif
    dsx = 0.
    dsy = 0.
    dsz = 0.
    fxj = 0.
    fyj = 0.
    fzj = 0.
    if (rr2 < (radkern*hsoft)**2) then
       !
       ! if the sink particle is given a softening length, soften the
       ! force and potential if r < radkern*hsoft
       !
       hsoft1 = 1.0/hsoft
       hsoft21= hsoft1**2
       q2i    = rr2*hsoft21
       qi     = sqrt(q2i)
       call kernel_softening(q2i,qi,psoft,fsoft)  ! Note: psoft < 0

       ! acceleration of gas due to point mass particle
       f1     = pmassj*fsoft*hsoft21*ddr
       !ftmpxi = ftmpxi - dx*f1
       !ftmpyi = ftmpyi - dy*f1
       !ftmpzi = ftmpzi - dz*f1
       !phi    = phi + pmassj*psoft*hsoft1  ! potential (spline-softened)

       !! acceleration of sink from gas
       !if (tofrom) f2 = pmassi*fsoft*hsoft21*ddr
    else
       ! no softening on the sink-gas interaction
       dr3  = ddr*ddr*ddr

       ! acceleration of gas due to point mass particle
       f1     = pmassj*dr3

       !! acceleration of sink from gas
       !if (tofrom) f2 = pmassi*dr3

       !ftmpxi = ftmpxi - dx*f1
       !ftmpyi = ftmpyi - dy*f1
       !ftmpzi = ftmpzi - dz*f1
       !phi    = phi    - pmassj*ddr      ! potential (GM/r)
    endif

    !wind acceleration
    !if (k==1) then !use this if statement for no radiative inhibition; comment out this if statement for using radiative inhibition
       !kappa follows gas
       !ddinv = xyzmh_ptmass(ihacc,iwindorigj) * ddr ! Radius/distance = R/r
       !!kappa follows radiation field
       ddinv = xyzmh_ptmass(ihacc,j) * ddr ! Radius/distance = R/r
       if (ddinv<1.d0) then
          !kappa follows gas
          !kappa_windaccel = windaccel(ikappac,iwindorigj) + &
          !                  windaccel(ikappar,iwindorigj)*(1.0d0-ddinv)**(2.d0*windaccel(ivbeta,iwindorigj)-1.d0)
          !!kappa follows radiation field
          kappa_windaccel = windaccel(ikappac,j) + &
                            windaccel(ikappar,j)*(1.0d0-ddinv)**(2.d0*windaccel(ivbeta,j)-1.d0)
       elseif (ddinv>1.0d0) then
          kappa_windaccel = 0.d0
          !print*,'particle originating from star ',iwindorigj,' is within the radius of star ',j,': ',xyzmh_ptmass(ihacc,iwindorigj),rr2
       else
          kappa_windaccel = 0.d0
          !print*,'particle originating from star ',iwindorigj,' is exactly at the radius of star ',j,': ',xyzmh_ptmass(ihacc,iwindorigj),rr2
       endif 
    !endif
    xantigr = windaccel(ixantgrav,j)*kappa_windaccel
    ftmpxi = ftmpxi + xantigr*dx*f1
    ftmpyi = ftmpyi + xantigr*dy*f1
    ftmpzi = ftmpzi + xantigr*dz*f1
    
    !if (jj>33552) then
    !   print*, jj,kappa_windaccel,windaccel(ixantgrav,j),xantigr,xantigr*dx*f1,dx*f1
    !endif

    !!keep this old method around to see how to mix opacities for the kappa_bar implementation via Eq. A6
    !!--- if kappa_bar is implemented, this would need to be pulled out of the current loop over point masses since the kappa for each particle interacting with each point mass would need to be computed
    !!--- additionally, the density contribution from each star's wind would need to be computed, so the density loop would need to be modified as well
    !!--- this code snippet was in the density loop
    !DO j=1,nptmass
    !   IF (akappar(iptm).NE.0.0) THEN
    !      ddinv = rptmas(j) / SQRT((xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2)
    !      akapp(j) = akappac(j) + akappar(j)*(1.0d0-ddinv)**(2.0d0*vbeta(j)-1.0d0)
    !   ELSE
    !      akapp(j) = akappac(j)
    !   ENDIF
    !ENDDO
    !!IF (iexf.EQ.7) THEN
    !akappa(ipart) = akapp(ibelong(ipart))
    !!ELSE
    !!   IF (rhoi.NE.0.0) THEN
    !!      akappa(ipart) = akapp(1)+(akapp(2)-akapp(1))*rhoip2/rhoi
    !!   ELSE
    !!      akappa(ipart) = 0.0
    !!   ENDIF
    !!ENDIF

    !if (tofrom) then
    !   ! backreaction of gas onto sink
    !   fxyz_ptmass(1,j) = fxyz_ptmass(1,j) + dx*f2 + fxj*pmassi/pmassj
    !   fxyz_ptmass(2,j) = fxyz_ptmass(2,j) + dy*f2 + fyj*pmassi/pmassj
    !   fxyz_ptmass(3,j) = fxyz_ptmass(3,j) + dz*f2 + fzj*pmassi/pmassj
    !
    !   ! backreaction torque of gas onto oblate sink
    !   dsdt_ptmass(1,j) = dsdt_ptmass(1,j) + pmassi*dsx
    !   dsdt_ptmass(2,j) = dsdt_ptmass(2,j) + pmassi*dsy
    !   dsdt_ptmass(3,j) = dsdt_ptmass(3,j) + pmassi*dsz
    !
    !   ! timestep is sqrt(separation/force)
    !   fonrmax = max(f1,f2,fonrmax)
    !   if (kappa) then
    !      if (abs(bin_info(isemi,j))>tiny(f2)) then
    !         bin_info(ipert,j) = bin_info(ipert,j) + f2
    !      endif
    !   endif
    !endif
 enddo
 !
 ! external force timestep based on sqrt(phi)/accel
 !
 !if (present(dtphi2)) then
 !   if (abs(phi) > epsilon(phi)) then
 !      f2     = ftmpxi*ftmpxi + ftmpyi*ftmpyi + ftmpzi*ftmpzi
 !      !dtphi is sqrt of this, but for optimisation we take the sqrt outside of
 !      !the loop
 !      dtphi2 = dtfacphi2*abs(phi)/f2
 !   else
 !      dtphi2 = huge(dtphi2)
 !   endif
 !endif
 !
 ! add temporary sums to existing force on gas particle
 !
 fxi = fxi + ftmpxi
 fyi = fyi + ftmpyi
 fzi = fzi + ftmpzi

end subroutine get_windaccel_force

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_windaccel(iunit)
 !use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

end subroutine write_options_windaccel

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_windaccel(name,valstring,imatch,igotall,ierr)
 !use io,      only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 !integer, save :: ngot = 0
 !character(len=30), parameter :: label = 'read_options_windaccel'

 igotall=.true.
 imatch = .false.
 ierr=0

end subroutine read_options_windaccel

end module extern_windaccel
