!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module metric_tools
!
! This module contains wrapper subroutines to get:
!      - The metric (covariant and contravariant)
!      - Derivatives of the covariant metric
!  As well as some general tools that are not specfic to each metric:
!      - Numerical metric derivatives
!      - Tensor transformations
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: inverse4x4, metric
!
 use metric, only:imetric
 implicit none

 !--- List of coordinates
 integer, public, parameter :: &
    icoord_cartesian  = 1,     &    ! Cartesian coordinates
    icoord_spherical  = 2           ! Spherical coordinates

 !--- List of metrics
 integer, public, parameter :: &
    imet_minkowski      = 1,   &    ! Minkowski metric
    imet_schwarzschild  = 2,   &    ! Schwarzschild metric
    imet_kerr           = 3,   &    ! Kerr metric
    imet_et             = 6         ! Tabulated metric from Einstein toolkit

 !--- Choice of coordinate system
 !    (When using this with PHANTOM, it should always be set to cartesian)
 integer, public, parameter :: icoordinate = icoord_cartesian

 !--- Choice for contravariant metric
 !    false  ->  use analytic contravariant metric
 !    true   ->  invert the covariant metric
 logical, private, parameter :: useinv4x4 = .true.

 public :: get_metric, get_metric_derivs, print_metricinfo, init_metric, pack_metric, unpack_metric
 public :: pack_metricderivs
 public :: imetric
 public :: numerical_metric_derivs

 private

contains

!-------------------------------------------------------------------------------
!+
!  This is a wrapper subroutine to get the metric tensor in both
!  covariant (gcov) and contravariant (gcon) form.
!+
!-------------------------------------------------------------------------------
pure subroutine get_metric(position,gcov,gcon,sqrtg)
 use metric,     only:get_metric_cartesian,get_metric_spherical,cartesian2spherical
 use inverse4x4, only:inv4x4
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
 real :: det

 select case(icoordinate)
 case(icoord_cartesian)
    if (useinv4x4) then
       call get_metric_cartesian(position,gcov,sqrtg=sqrtg)
       call inv4x4(gcov,gcon,det)
    else
       call get_metric_cartesian(position,gcov,gcon=gcon,sqrtg=sqrtg)
    endif
 case(icoord_spherical)
    if (useinv4x4) then
       call get_metric_spherical(position,gcov,sqrtg=sqrtg)
       call inv4x4(gcov,gcon,det)
    else
       call get_metric_spherical(position,gcov,gcon=gcon,sqrtg=sqrtg)
    endif
 end select

end subroutine get_metric

!-------------------------------------------------------------------------------
!+
!  This is a wrapper subroutine to get the derivatives of the covariant metric tensor.
!  The actual analytic metric derivaties are in the metric module, which are different for each type
!  of metric.
!+
!-------------------------------------------------------------------------------
subroutine get_metric_derivs(position,dgcovdx1,dgcovdx2,dgcovdx3)
 use metric, only:metric_cartesian_derivatives,metric_spherical_derivatives,imetric
 real, intent(in)  :: position(3)
 real, intent(out) :: dgcovdx1(0:3,0:3), dgcovdx2(0:3,0:3), dgcovdx3(0:3,0:3)

 select case(icoordinate)
 case(icoord_cartesian)
    call metric_cartesian_derivatives(position,dgcovdx1, dgcovdx2, dgcovdx3)
 case(icoord_spherical)
    call metric_spherical_derivatives(position,dgcovdx1, dgcovdx2, dgcovdx3)
 end select

end subroutine get_metric_derivs

!-------------------------------------------------------------------------------
!+
!  Numerical derivatives of the covariant metric tensor
!+
!-------------------------------------------------------------------------------
pure subroutine numerical_metric_derivs(position,dgcovdx,dgcovdy,dgcovdz)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdx,dgcovdy,dgcovdz
 real :: gblah(0:3,0:3), temp(3), gplus(0:3,0:3),gminus(0:3,0:3),dx,dy,dz,di,sqrtgblag
 di = 1.e-8
 dx = di
 dy = di
 dz = di

 temp      = position
 temp(1)   = temp(1)+dx
 call get_metric(temp,gplus,gblah,sqrtgblag)
 temp      = position
 temp(1)   = temp(1)-dx
 call get_metric(temp,gminus,gblah,sqrtgblag)
 dgcovdx = 0.5*(gplus-gminus)/dx

 temp      = position
 temp(2)   = temp(2)+dy
 call get_metric(temp,gplus,gblah,sqrtgblag)
 temp      = position
 temp(2)   = temp(2)-dy
 call get_metric(temp,gminus,gblah,sqrtgblag)
 dgcovdy = 0.5*(gplus-gminus)/dy

 temp      = position
 temp(3)   = temp(3)+dz
 call get_metric(temp,gplus,gblah,sqrtgblag)
 temp      = position
 temp(3)   = temp(3)-dz
 call get_metric(temp,gminus,gblah,sqrtgblag)
 dgcovdz = 0.5*(gplus-gminus)/dz

end subroutine numerical_metric_derivs

!-------------------------------------------------------------------------------
!+
!  print the metric type
!+
!-------------------------------------------------------------------------------
subroutine print_metricinfo(iprint)
 use metric, only:metric_type
 integer, intent(in) :: iprint

 write(iprint,*) 'Metric = ',trim(metric_type)

end subroutine print_metricinfo

!-------------------------------------------------------------------------------
!+
!  initialise arrays for the metric and metric derivatives
!+
!-------------------------------------------------------------------------------
subroutine init_metric(npart,xyzh,metrics,metricderivs)
 integer,         intent(in)  :: npart
 real,            intent(in)  :: xyzh(:,:)
 real,            intent(out) :: metrics(:,:,:,:)
 real, optional,  intent(out) :: metricderivs(:,:,:,:)
 integer :: i

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,metrics) &
 !$omp private(i)
 do i=1,npart
    call pack_metric(xyzh(1:3,i),metrics(:,:,:,i))
 enddo
 !omp end parallel do

 if (present(metricderivs)) then
    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,metricderivs) &
    !$omp private(i)
    do i=1,npart
       call pack_metricderivs(xyzh(1:3,i),metricderivs(:,:,:,i))
    enddo
    !omp end parallel do
 endif

end subroutine init_metric

!-------------------------------------------------------------------------------
!+
!  subroutine to pack the metric into a 4x4x2 array
!+
!-------------------------------------------------------------------------------
pure subroutine pack_metric(xyz,metrici)
 real, intent(in)  :: xyz(3)
 real, intent(out) :: metrici(:,:,:)
 real :: sqrtg

 call get_metric(xyz,gcov=metrici(:,:,1),gcon=metrici(:,:,2),sqrtg=sqrtg)

end subroutine pack_metric

!-------------------------------------------------------------------------------
!+
!  subroutine to pack the metric derivatives into a 4x4x3 array
!+
!-------------------------------------------------------------------------------
subroutine pack_metricderivs(xyzi,metricderivsi)
 real, intent(in)  :: xyzi(3)
 real, intent(out) :: metricderivsi(0:3,0:3,3)

 call get_metric_derivs(xyzi,metricderivsi(:,:,1),metricderivsi(:,:,2),metricderivsi(:,:,3))

end subroutine pack_metricderivs

!-------------------------------------------------------------------------------
!+
!  Subroutine to return metric/components from metrici array
!+
!-------------------------------------------------------------------------------
pure subroutine unpack_metric(metrici,gcov,gcon,gammaijdown,gammaijUP,alpha,betadown,betaUP)
 real, intent(in), dimension(0:3,0:3,2) :: metrici
 real, intent(out), dimension(0:3,0:3), optional :: gcov,gcon
 real, intent(out), dimension(1:3,1:3), optional :: gammaijdown,gammaijUP
 real, intent(out),                     optional :: alpha,betadown(3),betaUP(3)
 integer :: i

 if (present(alpha)) alpha  = sqrt(-1./metrici(0,0,2))

 if (present(betaUP)) betaUP = metrici(0,1:3,2) * (-1./metrici(0,0,2)) ! = gcon(0,1:3)*alpha**2

 if (present(gammaijUP)) then
    gammaijUP = 0.
    do i=1,3
       gammaijUP(:,i) = metrici(1:3,i,2) + metrici(1:3,0,2)*metrici(i,0,2)*(-1./metrici(0,0,2)) ! = gcon(i,j) + betaUP(i)*betaUP(j)/alpha**2
    enddo
 endif

 if (present(gcov))        gcov        = metrici(0:3,0:3,1)
 if (present(gcon))        gcon        = metrici(0:3,0:3,2)
 if (present(gammaijdown)) gammaijdown = metrici(1:3,1:3,1)
 if (present(betadown))    betadown    = metrici(0,1:3,1)

end subroutine unpack_metric

end module metric_tools
