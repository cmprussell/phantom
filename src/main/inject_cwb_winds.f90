!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module inject
!
!!!!! Wind injection from galactic centre stars
!!!!!   Written by Daniel Price, Jorge Cuadra, and Christopher Russell
! Wind injection from colliding wind binary stars
!   Written by Christopher Russell, which is based on the Galactic Center code that was
!   Written by Daniel Price, Jorge Cuadra, and Christopher Russell
!
! :References: Cuadra et al. (2008), MNRAS 383, 458
!
! :Owner: Daniel Price and Christopher Russell
!
! :Runtime parameters:
!   - datafile       : *name of data file for wind injection*
!   - outer_boundary : *kill gas particles outside this radius*
!
! :Dependencies: dim, eos, infile_utils, io, part, partinject, physcon,
!   random, units
!
 use dim,  only:maxptmass
 use part, only:nptmass
 implicit none
 character(len=*), parameter, public :: inject_type = 'cwb_winds'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,update_injected_par

 real :: outer_boundary = 20.
 character(len=120) :: datafile = 'winddata.txt'

 ! enumerated type for wind properties
 integer, parameter :: n_wind_prop = 2
 integer, parameter :: i_vel   = 1, &
                       i_Mdot  = 2

 REAL :: dtinject_cwb = MIN(1.e99,0.1*HUGE(dtinject_cwb)) !prevents 3-digit exponents, which fortran writes weirdly
 INTEGER :: ninjectmax_cwb = 100

 ! array containing properties of the wind from each star
 real,    private :: wind(n_wind_prop,maxptmass)
 integer, private :: total_particles_injected(maxptmass) = 0
 logical, private :: first_iteration = .true.
 integer, private :: iseed = -666

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
use part, only:massoftype,igas
use units, only:umass,utime
use physcon, only:solarm,years
 integer, intent(out) :: ierr
 !
 ! return without error
 !
 ierr = 0

 !limit the max number of particles injected per timestep per star to ninjectmax_cwb by specifiying idtinject_cwb, 
 !   the time interval it takes to inject ninjectmax_cwb particles for the wind with the highest Mdot
 !dtinject_cwb is written to cwb.in
 !dtinject_cwb becomes dtinject in the evol(...) subroutine
 !time interval = max mass to inject per timestep / max mass-loss rate [in code units]
 dtinject_cwb = ninjectmax_cwb*massoftype(igas) / (MAXVAL(wind(i_Mdot,1:nptmass)) * (solarm/umass)*(utime/years))
WRITE(*,*) 'CWB-specific code: ninjectmax_cwb =',ninjectmax_cwb
WRITE(*,*) 'CWB-specific code: dtinject_cwb =',dtinject_cwb

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling injection at the L1 point.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use io,        only:fatal,iverbose
 use part,      only:massoftype,igas,ihacc,i_tlast
 use partinject,only:add_or_update_particle
 use physcon,   only:pi,solarm,seconds,years,km,kb_on_mH
 use units,     only:umass,udist,utime,unit_velocity
 use random,    only:ran2
 use eos,       only:gmw,gamma
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 !real,    intent(out)   :: dtinject !for CWB, dtinject is not set here
 real,    intent(in)    :: dtinject
 real :: r2,Mcut,Mdot_fac,vel_fac,Minject,Mdot_code,tlast
 real :: xyzi(3),vxyz(3),xyz_star(3),vxyz_star(3),dir(3)
 real :: rr,phi,theta,cosphi,sinphi,costheta,sintheta
 real :: deltat,h,u,vinject,temp_inject,uu_inject,gam1
 integer :: i,j,k,nskip,i_part,ninject
!
! kill particles outside some outer radius
!
 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,outer_boundary) &
 !$omp private(i,r2)
 do i=1,npart
    r2 = xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2
    if (r2 > outer_boundary**2) xyzh(4,i) = -abs(xyzh(4,i))
 enddo
 !$omp end parallel do

 Mcut = 1000.*(solarm/umass)
 !nskip = 1 was set to 1 to account for SMBH; should now stay nskip=0
 nskip = 0
 do while(xyzmh_ptmass(4,nskip) > Mcut)
    nskip = nskip + 1
 enddo
 if (iverbose >= 2) print*,' skipping ',nskip,' point masses'
!WRITE(*,*) "nskip = ",nskip
!
! convert mass loss rate from Msun/yr to code units
!
 Mdot_fac = (solarm/umass)*(utime/years)
 vel_fac  = (km/udist)*(utime/seconds)

!
! If restarting, compute the number of particles already injected from each star.
!    This overestimates the total number injected by 1 timestep since 'time'
!    should be the last value before the restart dump was written, which is less
!    than its current value.  Therefore, the first time through no particles will
!    be injected.  This error is small though. A better idea is to add
!    'total_particles_injected' to the dump file, which will eliminate this
!    error altogether.
! Note: I imagine there's a better place to put this.  This if statement will
!    evaluate to false an awfully large number of times.  Needs to be after
!    dumpfile ('time') and wind data (Mdots) have been read in.
!
 if (first_iteration) then
    if (time /= 0) then   ! only if restarting
       do i=nskip+1,nptmass
          j = i - nskip ! position in wind table
          total_particles_injected(i) = int(wind(i_Mdot,j)*Mdot_fac * time / massoftype(igas))
       enddo
       print*
       print*, 'cwb initialization: wind particles already injected (total_particles_injected) =',&
               total_particles_injected(1:nptmass)
       print*
    endif
    first_iteration = .false.
 endif

 temp_inject = 1.e4
 gam1 = gamma - 1.
 if (gam1 <= 0) call fatal('inject','require gamma > 1 for wind injection')
!
! convert from temperature to thermal energy
! P/rho = kT/(mu m_H) = (gam-1)*u
!
 uu_inject = temp_inject * (((kb_on_mh) / unit_velocity)/unit_velocity) / (gmw*gam1)
 !print*,' uu_inject = ',uu_inject,kb_on_mh,unit_velocity,gmw,gam1
!
! loop over all wind particles
!
 !!$omp parallel do default(none) &
 !!$omp shared(nptmass)
 do i=nskip+1,nptmass
    !
    ! extract current position, velocity and injection radius of star
    !
    xyz_star  = xyzmh_ptmass(1:3,i)
    rr        = 1.0001*xyzmh_ptmass(ihacc,i)
    tlast     = xyzmh_ptmass(i_tlast,i)
    vxyz_star = vxyz_ptmass(1:3,i)

    !
    ! calculate how much mass to inject based on
    ! time interval since last injection
    !
    j = i - nskip ! position in wind table
    Mdot_code = wind(i_Mdot,j)*Mdot_fac
    vinject   = wind(i_vel,j)*vel_fac
    deltat    = time - tlast
    Minject   = Mdot_code*time
    !
    ! divide by mass of gas particles
    !
    ninject = int(Minject/massoftype(igas))-total_particles_injected(i)
    if (iverbose >= 2) print*,' point mass ',i,j,' injecting ',&
                       ninject,Minject-total_particles_injected(i)*massoftype(igas),massoftype(igas),time,tlast
!WRITE(*,*) 'Mdot, Mdot_code =',wind(i_Mdot,j), Mdot_code
!WRITE(*,*) 'v, vinject =',wind(i_vel,j), vinject
!WRITE(*,*) 'time, tlast =',time,tlast
!WRITE(*,*) 'point mass ',i,j,' injecting ',ninject,' -- total_particles_injected(',i,') = ',total_particles_injected(i)
!WRITE(*,'(4(A,I0))') 'point mass ',j,' injecting ',ninject,' -- total_particles_injected(',i,') = ',total_particles_injected(i)

    !
    ! this if statement is no longer essential for more accurate mass-loss rates,
    !    but it should help with setting h since tlast --> deltat is more accurate
    !
    ! don't update tlast for a particular star unless that star injected
    !    particles this timestep; this way, fractional particles/timestep can
    !    accumulate and eventually inject a particle, making Mdot more accurate
    !
    if (ninject > 0) then
!WRITE(*,'(4(A,I0),A,1PE15.8)') 'point mass ',j,' injecting ',ninject,' -- total_particles_injected(',i,')+ninject = ',total_particles_injected(i)+ninject,', time = ',time
!i_part=1
       do k=1,ninject
          !
          ! get random position on sphere
          !
          phi = 2.*pi*(ran2(iseed) - 0.5)
          theta = acos(2.*ran2(iseed) - 1.)
          sintheta = sin(theta)
          costheta = cos(theta)
          sinphi   = sin(phi)
          cosphi   = cos(phi)
          dir  = (/sintheta*cosphi,sintheta*sinphi,costheta/)

          xyzi = rr*dir + xyz_star
          vxyz = vinject*dir + vxyz_star
          !print*,' v = ',vinject,vxyz_star
          !print*,rr,vinject*deltat,100.*rr
          h = max(rr,10.*vinject*deltat) !/ninject**(1./3.)

          u = uu_inject

          i_part = npart + 1 ! all particles are new
!DO WHILE(xyzh(4,i_part)>0 .AND. i_part<npart+1)
!i_part=i_part+1
!ENDDO
!WRITE(*,*) 'i_part = ',i_part,', npart = ',npart
          call add_or_update_particle(igas, xyzi, vxyz, h, u, i_part, npart, npartoftype, xyzh, vxyzu)
       enddo
       !
       ! update tlast to current time
       !
       xyzmh_ptmass(i_tlast,i) = time
       !
       ! update total particles injected for this star
       !
       total_particles_injected(i) = total_particles_injected(i) + ninject
    endif
 enddo
 if (iverbose >= 2) then
    print*,'npart = ',npart
    print*,'tpi = ',total_particles_injected(1:nptmass)
 endif
 !
 !-- no constraint on timestep
 !
 !dtinject = huge(dtinject)
!WRITE(*,*) 'CWB-specific code: Did not write dtinject here; dtinject = ',dtinject

end subroutine inject_particles

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

!WRITE(*,*) 'Started write_options_inject'
 write(iunit,"(/,a)") '# options controlling particle injection'
 call write_inopt(trim(datafile),'datafile','name of data file for wind injection',iunit)
 call write_inopt(outer_boundary,'outer_boundary','kill gas particles outside this radius',iunit)
 call write_inopt(dtinject_cwb,'dtinject_cwb','timestep limit based on number of particles injected per timestep',iunit)
 call write_inopt(ninjectmax_cwb,'ninjectmax_cwb','max particles to inject per timestep per star',iunit)
!WRITE(*,*) 'Ended   write_options_inject'

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal,error,warning
 use physcon, only:solarm,years
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'
 integer :: nstars

!WRITE(*,*) 'Started read_options_inject'
 imatch  = .true.
 select case(trim(name))
 case('outer_boundary')
    read(valstring,*,iostat=ierr) outer_boundary
 case('datafile')
    read(valstring,*,iostat=ierr) datafile
    call read_wind_data(datafile,nstars)
    !nstars is determined here from startrun --> read_infile
    !nptmass is read in from startrun --> read_dump, which has yet to be executed
    !so since nptmass has not been defined yet, comparing it to nstars is pointless
    !if (nstars /= nptmass) then
    !   call warning('read_options_inject','number of stars /= number of wind sources')
    !    WRITE(*,*) 'nstars  =',nstars,', nptmass =',nptmass
    !endif
    ngot = ngot + 1
 case('dtinject_cwb')
    WRITE(*,*) 'CWB-specific code: pre-read  dtinject_cwb = ',dtinject_cwb
    read(valstring,*,iostat=ierr) dtinject_cwb
    WRITE(*,*) 'CWB-specific code: post-read dtinject_cwb = ',dtinject_cwb
 case('ninjectmax_cwb')
    WRITE(*,*) 'CWB-specific code: pre-read  ninjectmax_cwb = ',ninjectmax_cwb
    read(valstring,*,iostat=ierr) ninjectmax_cwb
    WRITE(*,*) 'CWB-specific code: post-read ninjectmax_cwb = ',ninjectmax_cwb
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

!WRITE(*,*) 'Ended   read_options_inject'
end subroutine read_options_inject

!-----------------------------------------------------------------------
!+
!  Reads wind input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_wind_data(filename,nstars)
 use io, only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: nstars
 integer :: iunit,ierr,idum,i

 nstars = 0
 wind = 0.
 open(newunit=iunit,file=trim(filename),status='old',action='read',iostat=ierr)
 do while(ierr == 0)
    nstars = nstars + 1
    if (nstars <= maxptmass) then
       read(iunit,*,iostat=ierr) idum,wind(i_vel,nstars),wind(i_Mdot,nstars)
       if (ierr /= 0) nstars = nstars - 1
    else
       call error('read_wind_data','array bounds exceeded')
    endif
 enddo

 if (nstars > 0) print "(1x,37('-'),/,1x,a,'|',2(a15,1x,'|'),/,1x,37('-'))",&
                        'ID',' Wind Vel(km/s)',' Mdot(Msun/yr)'
 do i=1,nstars
    print "(i3,'|',2(1pg15.4,1x,'|'))",i,wind(i_vel,i),wind(i_Mdot,i)
 enddo
 if (nstars > 0) print "(1x,37('-'))"

end subroutine read_wind_data

end module inject
