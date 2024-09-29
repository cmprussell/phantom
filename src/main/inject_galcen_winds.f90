!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Wind injection from galactic centre stars
!   Written by Daniel Price, Jorge Cuadra, and Christopher Russell
!
! :References: Cuadra et al. (2008), MNRAS 383, 458
!
! :Owner: Daniel Price
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
 character(len=*), parameter, public :: inject_type = 'galcen_winds'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      set_default_options_inject,update_injected_par

 real :: outer_boundary = 20.
 character(len=120) :: datafile = 'winddata.txt'

 ! enumerated type for wind properties
 integer, parameter :: n_wind_prop = 2
 integer, parameter :: i_vel   = 1, &
                       i_Mdot  = 2

 ! array containing properties of the wind from each star
 real,    private :: wind(n_wind_prop,maxptmass)
 integer, private :: total_particles_injected(maxptmass) = 0
 logical, private :: first_iteration = .true.
 integer, private :: iseed = -666
 integer, private :: nskip_ptmass
 real,    private :: temp_inject=0. !set to value that should cause an error if value in init_inject() is not run properly
 real,    private :: uu_inject=0.   !set to value that should cause an error if value in init_inject() is not run properly
 real,    private :: Mdot_fac=0.    !set to value that should cause an error if value in init_inject() is not run properly
 real,    private :: vel_fac=0.     !set to value that should cause an error if value in init_inject() is not run properly


contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use io,        only:fatal
 use part,      only:xyzmh_ptmass
 use physcon,   only:solarm,seconds,years,km,kb_on_mH
 use units,     only:umass,udist,utime,unit_velocity
 use eos,       only:gmw,gamma
 integer, intent(out) :: ierr
 integer :: i
 real    :: Mcut

 !
 ! return without error
 !
 ierr = 0

 !
 ! correlate pointmass postion/velocity/mass table and pointmass wind table
 ! for wind injection, skip SMBH and IMBH -- which are above 200Msun -- since
 ! black holes don't appear in the wind table
 !
 Mcut = 200.*(solarm/umass)
 nskip_ptmass = 0
 do while(xyzmh_ptmass(4,nskip_ptmass+1) > Mcut .and. nskip_ptmass<nptmass)
    nskip_ptmass = nskip_ptmass + 1
 enddo
 if (nskip_ptmass==nptmass) then
    ierr = ierr + 1
    print*,' ERROR: no winds since all read-in point masses are skipped'
 endif
 print "(a,i0,a,f0.1,a)", ' Skipping ',nskip_ptmass,' point masses for wind injection since their point masses are >200Msun,'
 print "(a)", '   which is presently interpreted to mean that these point masses are black holes.'
 !
 ! convert mass loss rate from Msun/yr to code units
 !
 Mdot_fac = (solarm/umass)*(utime/years)
 vel_fac  = (km/udist)*(utime/seconds)
 ! 
 ! verification of Mdots and vinfs
 !
 print "(/a)", ' Pointmass table relating the following quantities, which is assembled from two different input tables:'
 print "(a)",  '    1. index, which is used in xyzmh_ptmass -- from position/velocity/mass table'
 print "(a)",  '    2. mass in a variety of units -- from position/velocity/mass table'
 print "(a)",  '    3. wind properties -- from wind table (which excludes entries for black holes)'
 print "(1x,90('-'))"
 print "(a)",  '  i |  M(code units)  |  M(g)           |  M(Msun)        | Mdot(Msun/yr)  |  vinf(km/s) |'
 print "(1x,90('-'))"
 do i=1,nptmass
    if (i<=nskip_ptmass) then
       print "(i3,1x,'|',3(es16.8,1x,'|'),a)", i,xyzmh_ptmass(4,i),xyzmh_ptmass(4,i)*umass, &
          xyzmh_ptmass(4,i)*umass/solarm,'  skipped -- no wind injection for this pointmass'
    else
       print "(i3,1x,'|',4(es16.8,1x,'|'),3x,f9.1,1x,'|')", i,xyzmh_ptmass(4,i),xyzmh_ptmass(4,i)*umass, &
          xyzmh_ptmass(4,i)*umass/solarm,wind(i_Mdot,i-nskip_ptmass),wind(i_vel,i-nskip_ptmass)
    endif
 enddo
 print "(1x,90('-'))"

 !
 ! set more global wind-injection properties
 !
 temp_inject = 1.e4
 print "(/,a,es14.6,a)", ' injected wind particles will have a temperature of temp_inject =',temp_inject,'K'
 if (gamma<=1.) call fatal('inject','require gamma > 1 for wind injection')
 !
 ! convert from temperature to thermal energy
 ! P/rho = kT/(mu m_H) = (gam-1)*u
 !
 uu_inject = temp_inject * kb_on_mh / unit_velocity**2 / (gmw*(gamma-1.))

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection from the stars.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use io,        only:fatal,iverbose
 use part,      only:massoftype,igas,ihacc,i_tlast,iwindorig,isdead_or_accreted,iphase,iunknown
 use partinject,only:add_or_update_particle,updated_particle
 use physcon,   only:pi,solarm,seconds,years,km,kb_on_mH
 use random,    only:ran2
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real :: r2,Mdot_fac,vel_fac,Minject,Mdot_code,tlast
 real :: xyzi(3),vxyz(3),xyz_star(3),vxyz_star(3),dir(3)
 real :: rr,phi,theta,cosphi,sinphi,costheta,sintheta
 real :: deltat,h,u,vinject,uu_inject
 integer :: i,j,k,i_part,ninject

 !print*,'init: tpi = ',total_particles_injected(1:nptmass)
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

 if (iverbose >= 2) print*,' skipping ',nskip_ptmass,' point masses'

 !!
 !! Verification of Mdots and vinfs
 !!
 !do i=1,nptmass
 !   if (i<=nskip_ptmass) then
 !      write(*,*) 'm(',i,') = ',xyzmh_ptmass(4,i),xyzmh_ptmass(4,i)*umass,xyzmh_ptmass(4,i)*umass/solarm
 !   else
 !      write(*,*) 'm(',i,') = ',xyzmh_ptmass(4,i),xyzmh_ptmass(4,i)*umass,xyzmh_ptmass(4,i)*umass/solarm,wind(i_Mdot,i-nskip_ptmass),wind(i_vel,i-nskip_ptmass)
 !   endif
 !enddo
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

!
! loop over all wind particles
!
 i_part=1 !for tracking resuse particles
 !!$omp parallel do default(none) &
 !!$omp shared(nptmass)
 do i=nskip_ptmass+1,nptmass
    !
    ! calculate how much mass to inject based on
    ! time interval since last injection
    !
    j = i - nskip_ptmass ! position in wind table
    Mdot_code = wind(i_Mdot,j)*Mdot_fac
    vinject   = wind(i_vel,j)*vel_fac
    tlast     = xyzmh_ptmass(i_tlast,i)
    deltat    = time - tlast
    Minject   = Mdot_code*time
    !
    ! divide by mass of gas particles
    !
    ninject = int(Minject/massoftype(igas))-total_particles_injected(i)
    if (iverbose >= 2) print*,' point mass ',i,j,' injecting ',&
                       ninject,Minject-total_particles_injected(i)*massoftype(igas),massoftype(igas),time,tlast
    !
    ! this if statement is no longer essential for more accurate mass-loss rates,
    !    but it should help with setting h since tlast --> deltat is more accurate
    !
    ! don't update tlast for a particular star unless that star injected
    !    particles this timestep; this way, fractional particles/timestep can
    !    accumulate and eventually inject a particle, making Mdot more accurate
    !
    if (ninject > 0) then
       !
       ! extract current position, velocity and injection radius of star
       !
       xyz_star  = xyzmh_ptmass(1:3,i)
       rr        = 1.0001*xyzmh_ptmass(ihacc,i)
       vxyz_star = vxyz_ptmass(1:3,i)

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

          !!
          !! original method -- all particles are new
          !!
          !i_part = npart + 1 ! all particles are new

          !
          ! reuse particles
          ! find first dead or accreted particle
          !
          do while((.not.isdead_or_accreted(xyzh(4,i_part))) .and. i_part<npart+1)
             i_part=i_part+1
          enddo
          !
          ! only inject particles that won't be immediately killed due to being beyond the outer_boundary
          ! this can be used with the original all-particles-are-new method as well
          !
          if (sqrt(xyzmh_ptmass(1,i)**2+xyzmh_ptmass(2,i)**2+xyzmh_ptmass(3,i)**2) < outer_boundary) then
             !!
             !! track particles that are actually injected
             !!
             !ninject_actual = ninject_actual+1
             !
             ! add or update the particle via built-in method within partinject
             !
             call add_or_update_particle(igas, xyzi, vxyz, h, u, i_part, npart, npartoftype, xyzh, vxyzu)
             !
             ! star from where this wind particle originated
             !
             iwindorig(i_part) = i
             !note: add_or_update_particle increased npart by 1 if a new particle was added, so now the
             !         comparison is with npart (whereas above the comparison was with npart+1)
             if (i_part<npart) then
                !
                ! particle was updated, not added
                !
                updated_particle=.true.
                !
                ! flag this particle to update its timestep -- this overrides
                !    "call set_particle_type(particle_number,itype)" in partinject-->add_or_update_particle
                !
                iphase(i_part) = iunknown
                !call set_particle_type(i_part,iunknown) !alternative/equivalent to above line
                !
                ! begin the search for the next accreted or dead particle to reuse with the next particle
                !
                i_part = i_part + 1
             endif
             !else (not needed, but kept here for the following comment)
             !
             ! do not inject the particle, but keep track of it via total_particles_injected
             !    in case this star comes back within outer_boundary and then starts injecting particles
             !
          endif
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
 dtinject = huge(dtinject)

end subroutine inject_particles

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(trim(datafile),'datafile','name of data file for wind injection',iunit)
 call write_inopt(outer_boundary,'outer_boundary','kill gas particles outside this radius',iunit)

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

 imatch  = .true.
 select case(trim(name))
 case('outer_boundary')
    read(valstring,*,iostat=ierr) outer_boundary
 case('datafile')
    read(valstring,*,iostat=ierr) datafile
    call read_wind_data(datafile,nstars)
    !note: nptmass=0 here, so can't compare nstars and nptmass
    !if (nstars /= nptmass) then
    !   call warning('read_options_inject','number of stars /= number of wind sources')
    !endif
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 1)

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

subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
end subroutine set_default_options_inject

end module inject
