!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
! :Owner: Christopher Russell
!
! :Runtime parameters:
!   - datafile_mhn   : *filename for wind composition (mu,habund,name)*
!   - datafile_wind  : *filename for wind injection (m,x,y,z,vx,vy,vz)*
!   - outer_boundary : *kill gas particles outside this radius*
!
! :Dependencies: dim, eos, infile_utils, io, options, part, partinject,
!   physcon, random, timestep, units
!
 use dim,  only:maxptmass
 use part, only:nptmass
 implicit none
 character(len=*), parameter, public :: inject_type = 'galcen_winds'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      set_default_options_inject,update_injected_par,read_use_var_comp_data

 real :: outer_boundary = 20.
 character(len=120) :: datafile_wind = 'winddata.txt'
 character(len=120) :: datafile_mhn  = 'mu_habund_name.txt'

 ! enumerated type for wind properties that are reals
 integer, parameter :: n_wind_prop = 3
 integer, parameter :: i_vel  = 1, &
                       i_Mdot = 2, &
                       i_gmw  = 3
 ! enumerated type for wind properties that are integers
 integer, parameter :: n_wind_prop_int = 2
 integer, parameter :: i_astroID      = 1, &
                       i_comp         = 2
 ! enumerated type for wind properties that are characters
 integer, parameter :: n_wind_prop_char = 2
 integer, parameter :: i_comp_char = 1, &
                       i_fulltype_char = 2

 ! array containing properties of the wind from each star
 real,              private :: wind(     n_wind_prop     ,maxptmass) = 0.
 integer,           private :: wind_int( n_wind_prop_int ,maxptmass) = 0
 character(len=20), private :: wind_char(n_wind_prop_char,maxptmass) = ''
 integer, private :: total_particles_injected(maxptmass) = 0
 !integer, private :: total_particles_injected_actual(maxptmass) = 0
 integer, private :: iseed = -666
 integer, private :: nskip_ptmass
 real,    private :: temp_inject = 0. !set to value that should cause an error if value in init_inject() is not run properly
 real,    private :: uu_inject = 0.   !set to value that should cause an error if value in init_inject() is not run properly
 real,    private :: Mdot_fac = 0.    !set to value that should cause an error if value in init_inject() is not run properly
 real,    private :: vel_fac = 0.     !set to value that should cause an error if value in init_inject() is not run properly
 integer, private :: i_unknowninit_startypes !index in the mu/name table of the startype for unknown and initialization particles

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use io,        only:fatal
 use timestep,  only:time
 use part,      only:massoftype,igas,xyzmh_ptmass,iwindorig_to_ict,mu_startypes
 use physcon,   only:solarm,seconds,years,km,kb_on_mH
 use units,     only:umass,udist,utime,unit_velocity
 use eos,       only:gmw,gamma
 use timestep,  only:dtmax
 integer, intent(out) :: ierr
 integer :: i,j,j_corrupt
 real    :: Mcut
 logical :: iexist
 integer :: ierr_tpi=0
 real    :: time_tpi(10000) = 0.
 integer :: nptmass_tpi(10000) = 0
 integer :: total_particles_injected_tpi(maxptmass,10000) = 0
 integer :: i_curr,i_curr2,i_first
 logical :: tpi_read_from_file=.false.
 real    :: time_test
 integer :: nbinmax_test
 integer :: total_particles_injected_test(maxptmass) = 0
 integer :: total_diff
 real    :: tol_init_inject !tolerance for error when comparing time values
 logical :: tpi_file_rewrite_add = .false.
 logical :: tpi_file_rewrite_delete = .false.
 logical :: tpi_file_rewrite_reorder = .false.
 integer :: iunit_tpi

 !
 ! initialize to return without error
 !
 ierr = 0

 !
 ! display time to immediately verify whether this is a new sim or a restarting sim
 !
 if (time<1) then
    print "(/,a,g0/)", ' init_inject : time = 0',time
 else
    print "(/,a,g0/)", ' init_inject : time = ',time
 endif

 !
 ! correlate pointmass postion/velocity/mass table and pointmass wind table
 ! for wind injection, skip SMBH and IMBH -- which are above 200Msun -- since black holes don't appear in the wind table
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
 ! convert mass loss rate from Msun/yr to code units
 Mdot_fac = (solarm/umass)*(utime/years)
 vel_fac  = (km/udist)*(utime/seconds)
 ! make table that correlates iwindorig to ict (i.e. originating star to cooling tables)
 if (allocated(iwindorig_to_ict)) deallocate(iwindorig_to_ict)
 allocate(iwindorig_to_ict(0:nptmass)) !indez zero is for unknown/init particles
 ! verification of Mdots and vinfs
 print "(/a)", ' Pointmass table relating the following quantities, which is assembled from two different input tables:'
 print "(a)",  '    1. index, which is used in xyzmh_ptmass -- from position/velocity/mass table'
 print "(a)",  '    2. mass in a variety of units -- from position/velocity/mass table'
 print "(a)",  '    3. wind properties -- from wind table (which excludes entries for black holes)'
 print "(1x,108('-'))"
 print "(a)",  ' |   i |  M(code units)  |  M(g)           |  M(Msun)        |  Mdot(Msun/yr)  |  vinf(km/s) |  composition |'
 print "(a)",  ' |-----|-----------------|-----------------|-----------------|-----------------|-------------|--------------|'
 !print "(1x,108('-'))"
 do i=1,nptmass
    if (i<=nskip_ptmass) then
       print "(1x'|',1x,i3,1x,'|',3(es16.8,1x,'|'),a)", &
          i,xyzmh_ptmass(4,i),xyzmh_ptmass(4,i)*umass,xyzmh_ptmass(4,i)*umass/solarm,' --- skipped: no wind from this pointmass --- |'
       iwindorig_to_ict(i) = -1 !this is not a wind, so it does not have a cooling table -- this value should cause an error if ever used
    else
       print "(1x,'|',1x,i3,1x,'|',4(es16.8,1x,'|'),3x,f9.1,1x,'|',1x,i12,1x,'|')", i,xyzmh_ptmass(4,i),xyzmh_ptmass(4,i)*umass, &
          xyzmh_ptmass(4,i)*umass/solarm,wind(i_Mdot,i-nskip_ptmass),wind(i_vel,i-nskip_ptmass),wind_int(i_comp,i-nskip_ptmass)
       wind(i_gmw,i-nskip_ptmass) = mu_startypes(wind_int(i_comp,i-nskip_ptmass))
       !print*, 'wind(i_gmw,i-nskip_ptmass) = ',wind(i_gmw,i-nskip_ptmass)
       iwindorig_to_ict(i) = wind_int(i_comp,i-nskip_ptmass) !this is a wind, so it does have a cooling table
    endif
 enddo
 print "(1x,108('-'))"
 !unknown and initialization particles are stored in 0th index
 i = 0
 iwindorig_to_ict(i) = wind_int(i_comp,nptmass-nskip_ptmass+1)
 print "(1x,'|',1x,i3,1x,'|',1x,a,46(1x),1x,'|',1x,i12,1x,'|')", i,' unknown and initialization particles',iwindorig_to_ict(0)
 print "(1x,108('-'))"

 !
 ! if this is a new simulation, create and populate total_particles_injected.dat
 !
 ! if this is a restarting simulation, either
 !    A. read in total_particles_injected from total_particles_injected.dat
 !    B. compute total_particles_injected, with the slight discrepancy
 !       that the last time this would have been computed is one fraction of dtmax
 !       before "time" and therefore might have a slight error assiciated with it
 !
 if (time < tiny(time)) then
    ! new sim --> create total_particles_injected.dat
    print "(/,a)", ' New simulation: Creating total_particles_injected.dat to track total_particles_injected,'
    print "(a)", '     which will negate any injected-particle errors upon restarting.'
    total_particles_injected(1:nptmass) = 0
    open(newunit=iunit_tpi,file='total_particles_injected.dat',form='formatted')
    write(iunit_tpi,*) time,nptmass,total_particles_injected(1:nptmass)
    close(iunit_tpi)
 else
    ! compute particles that should have been injected already
    ! hopefully, this is just for comparison with total_particles_injected.dat values
    do i=nskip_ptmass+1,nptmass
       j = i - nskip_ptmass ! position in wind table
       total_particles_injected(i) = int(wind(i_Mdot,j)*Mdot_fac * time / massoftype(igas))
    enddo
    ! see if total_particles_injected.dat exists
    ! it should always exist for future galcen sims, but this verification will
    !    allow backwards compatibility
    inquire(file='total_particles_injected.dat',exist=iexist)
    if (iexist) then
       ! read-in total_particles_injected
       open(newunit=iunit_tpi,file='total_particles_injected.dat',form='formatted',status='old')
       ierr_tpi=0
       j_corrupt = 0
       j = 0
       do while(ierr_tpi==0)
          j = j+1
          read(iunit_tpi,*,iostat=ierr_tpi) time_tpi(j),nptmass_tpi(j),total_particles_injected_tpi(1:nptmass,j)
          if (ierr_tpi>0) then
             ! Corrupt entry: Skip it but keep reading from the file in hopes of one or
             !    more well-formatted entries later on in the file.
             ! Not sure if all compilers have unreadable error codes > 0 and EOF error codes < 0, upon
             !    which this distinction relies -- take out this if statement if things seem to be going
             !    haywire when reading in total_particles_injected.dat when not using ifort.
             j = j-1
             ierr_tpi = 0
             j_corrupt = j_corrupt+1
          endif
       enddo
       j = j-1
       close(iunit_tpi)
       print "(/,a,i0,a)", ' Read in ',j,' entries from total_particles_injected.dat.'
       if (j_corrupt>0) print "(a,i0,a)", ' Read in ',j_corrupt,' corrupt entries from total_particles_injected.dat, which will be removed from the file.'
       ! at least one entry was read in
       if (j>0) then
          ! find entry in total_particles_injected.dat that corresponds to this specific restart time
          tol_init_inject = dtmax*1.d-3 !tolerance for error when comparing time values
          !                             !make it a small number relative to dtmax since the successive entries in
          !                             !   total_particles_injected.dat should be separated by dtmax
          i_curr = 1
          do while(abs(time_tpi(i_curr)-time)>tol_init_inject .and. i_curr<j)
             i_curr = i_curr+1
          enddo
          ! verify that the correct entry was read -- if not, use the already
          !    computed values for total_particle_injected (or kill the sim)
          if (abs(time_tpi(i_curr)-time)>tol_init_inject) then
             ! correct entry in total_particles_injected.dat for "time" was not found
             print "(a)", ' Warning: the correct entry in total_particles_injected.dat does not seem to exist!'
             print "(a,g0)", '    time =',time
             do i=1,i_curr
                print "(a,i0,a,g0)", '    time_tpi(',i,') =',time_tpi(i)
             enddo
             print "(a)", ' total_particles_injected will be computed using "time".'

             !! kill the sim via reporting an error if the correct entry in total_particles_injected appears to be missing
             !ierr = ierr+1
             !print "(a)", ' ERROR: the correct entry in total_particles_injected.dat does not seem to exist!'

             ! attempt to make total_particles_injected.dat as useful as possible
             !    for future restarts by removing all entries past the current time
             !    and ensuring that the entries are in sequential order
             tpi_file_rewrite_add = .false.
             tpi_file_rewrite_delete = .false.
             tpi_file_rewrite_reorder = .false.
             do i=1,j
                if (time_tpi(i)<time) tpi_file_rewrite_add = .true.
                if (time_tpi(i)>time) tpi_file_rewrite_delete = .true.
             enddo
             do i=2,j
                if (time_tpi(i)<time_tpi(i-1)) tpi_file_rewrite_reorder = .true.
             enddo
             if (tpi_file_rewrite_add) then
                ! there are entries that should be in total_particles_injected.dat
                if (tpi_file_rewrite_delete .or. tpi_file_rewrite_reorder) then
                   ! rewriting total_particles_injected.dat is needed due to either deletion or reordering
                   ! Note: The reordering simply ensures that sequential entries are increasing.
                   !    At present, there is no fancy sort algorithm to bring a drastically
                   !    out-of-order list into sequential order while keeping all unique entries.
                   !
                   ! find the first entry to write
                   i_first = 1
                   do while (time_tpi(i_first)>time)
                      i_first = i_first+1
                   enddo
                   ! rewrite total_particles_injected.dat
                   print "(/,a)", ' total_particles_injected.dat seems to have irregular entries given this simulation''s restart variable "time".'
                   print "(a)", ' Rewrite this file with only the correct entries, which are values less than "time" and in sequential order.'
                   open(newunit=iunit_tpi,file='total_particles_injected.dat',form='formatted')
                   write(iunit_tpi,*,iostat=ierr_tpi) time_tpi(i_first),nptmass_tpi(i_first),total_particles_injected_tpi(1:nptmass,i_first)
                   i_curr2 = i_first
                   do i = i_first+1,j
                      if (time_tpi(i) < time .and. time_tpi(i)>time_tpi(i_curr2)) then
                         write(iunit_tpi,*,iostat=ierr_tpi) time_tpi(i),nptmass_tpi(i),total_particles_injected_tpi(1:nptmass,i)
                         i_curr2 = i
                      endif
                   enddo
                   close(iunit_tpi)
                   !else (not needed, but kept here for the following comment)
                   ! total_particles_injected.dat only has good entries -- no deleting or reordering needed -- so leave the file as is
                endif
             else
                ! there are no current entries in total_particles_injected.dat that should still be there
                ! replace total_particles_injected.dat with a file that has an entry only for the start
                !    of the simulation
                print "(/,a)", ' total_particles_injected.dat seems to have no relevant entries given this simulation''s restart variable "time".'
                print "(a)", ' Rewrite this file with only the correct entry for the start of the simulation.'
                total_particles_injected_tpi(1:nptmass,1) = 0
                open(newunit=iunit_tpi,file='total_particles_injected.dat',form='formatted')
                write(iunit_tpi,*,iostat=ierr_tpi) 0.0d0,nptmass,total_particles_injected_tpi(1:nptmass,1)
                close(iunit_tpi)
             endif
          else
             ! correct entry in total_particles_injected.dat for "time" was found
             tpi_read_from_file = .true. !correct entry for current restart time was found in total_particles_injected.dat
             !                           !   --> don't use total_particles_injected values computed above
             print "(a)", ' correct entry in total_particles_injected.dat has been found'
             print "(a,g0,a,g0,2(a,i0),a)", ' time = ',time,', time_tpi(correct) = ',time_tpi(i_curr),', correct index i_curr = ',i_curr,' (out of ',j,')'
             print "(a,i0)", ' nptmass_tpi(correct) = ',nptmass_tpi(i_curr)
             if (nptmass/=nptmass_tpi(i_curr)) then
                ! kill the sim via reporting an error if the number of pointmasses in total_particles_injected is incorrect
                ierr = ierr+1
                print "(a)", ' ERROR: nptmass from total_particles_injected.dat does not equal nptmass!'
                print "(a)", ' ERROR: nptmass = ',nptmass,', nptmass_tpi(i_curr) = ',nptmass_tpi(i_curr)
             endif
             print "(/,a)", ' total_particles_injected values for current restart time (correct) compared to the'
             print "(a)", '    calculation with the "time" variable (error prone) and their associated errors.'
             print "(a)", '    These errors are not in the computation due to using total_particles_injected.dat.'
             print "(a)", ' Note: These might not be the actual numbers of particles injected if stars beyond'
             print "(a)", '    outer_boundary do not have any wind particles injected [check inject_particles()].'
             print "(a)", '    These non-injected particles are still tracked so that if the star goes within'
             print "(a)", '    outer_boundary, then that star should start injecting wind particles normally.'
             print "(1x,55('-'))"
             print "(a)", '  i |    tpi_time(i) |    tpi_file(i) |          error |'
             print "(a)", '    |  (error prone) |      (correct) |    (time-file) |'
             print "(1x,55('-'))"
             do i=1,nptmass
                print "(i3,1x,'|',3(i15,1x,'|'))", i,total_particles_injected(i),total_particles_injected_tpi(i,i_curr),total_particles_injected(i)-total_particles_injected_tpi(i,i_curr)
             enddo
             print "(1x,55('-'))"
             total_particles_injected(1:nptmass) = total_particles_injected_tpi(1:nptmass,i_curr)
             ! see if total_particles_injected.dat has any out-of-order entries
             tpi_file_rewrite_reorder = .false.
             do i=2,j
                if (time_tpi(i)<time_tpi(i-1)) tpi_file_rewrite_reorder = .true.
             enddo
             ! if not restarting from the end of the last simulation run, remove any
             !    extra entries in total_particles_injected.dat by rewriting the file
             ! if entries are out of order, reorder them
             ! Note: The reordering simply ensures that sequential entries are increasing.
             !    At present, there is no fancy sort algorithm to bring a drastically
             !    out-of-order list into sequential order while keeping all unique entries.
             if (i_curr<j .or. tpi_file_rewrite_reorder) then
                print "(/,a)", ' Rewriting total_particles_injected.dat to remove unnecessary entries'
                print "(a)", '    from times that are after this current sim''s restart time,'
                print "(a)", '    to reorder the entries so they are in sequential order, or both.'
                open(newunit=iunit_tpi,file='total_particles_injected.dat',form='formatted')
                ! find the first entry to write
                i_first = 1
                do while (time_tpi(i_first)>time)
                   i_first = i_first+1
                enddo
                ! rewrite the file in sequential order ignoring entries beyond the sim's restart variable "time"
                write(iunit_tpi,*,iostat=ierr_tpi) time_tpi(i_first),nptmass_tpi(i_first),total_particles_injected_tpi(1:nptmass,i_first)
                i_curr2 = i_first
                do i = i_first+1,i_curr
                   if (time_tpi(i)<time+tol_init_inject .and. time_tpi(i)>time_tpi(i_curr2)) then
                      write(iunit_tpi,*,iostat=ierr_tpi) time_tpi(i),nptmass_tpi(i),total_particles_injected_tpi(1:nptmass,i)
                      i_curr2 = i
                   endif
                enddo
                close(iunit_tpi)
             endif
          endif
       else
          ! total_particles_injected.dat is either a blank or corrupt file since
          !    no entries were read in correctly.
          ! To be sure that it is not corrupt, make total_particles_injected.dat
          !    have a single entry for the start of the simulation.
          print "(/,a)", ' Rewriting total_particles_injected.dat to include the'
          print "(a)", '    entry for the start of the simulation.'
          total_particles_injected_tpi(1:nptmass,1) = 0
          open(newunit=iunit_tpi,file='total_particles_injected.dat',form='formatted')
          write(iunit_tpi,*,iostat=ierr_tpi) 0.0d0,nptmass,total_particles_injected_tpi(1:nptmass,1)
          close(iunit_tpi)
       endif
    else
       print "(/,a)", ' Warning: file total_particles_injected.dat does not exist!'
    endif
    if (.not.tpi_read_from_file) then
       ! compute particles that should have been injected already, noting that
       !    there is an error associated with "time" vs. "time-1/(2**nbinmax)*dtmax";
       !    the former is the current calculation and the latter is the actual
       !    last time value used for injecting particles prior to the restart
       print "(/,a)", ' ****************************************************************************'
       print "(a)", ' * Restarting sim without the correct entry in total_particles_injected.dat *'
       print "(a)", ' ****************************************************************************'
       print "(a)", ' Here are the wind particles that  should have been injected already for each star'
       print "(a)", '    (ignoring any errors in the injected-particle numbers based on the last'
       print "(a)", '    fractional timestep of dtmax).'
       print "(a)", ' Note: These might not be the actual numbers of particles injected if stars beyond'
       print "(a)", '    outer_boundary do not have any wind particles injected [check inject_particles()].'
       print "(a)", '    These non-injected particles are still tracked so that if the star goes within'
       print "(a)", '    outer_boundary, then that star should start injecting wind particles normally.'
       print "(1x,36('-'))"
       print "(a)", '  i |   total_particles_injected(i) |'
       print "(1x,36('-'))"
       do i=nskip_ptmass+1,nptmass
          !j = i - nskip_ptmass ! position in wind table
          !total_particles_injected(i) = int(wind(i_Mdot,j)*Mdot_fac * time / massoftype(igas))
          if (total_particles_injected(i)<0) then
             ierr = ierr + 1 !flag an error in the calculation of particles already injected
             print*,'ERROR: star ',i,' has negative particles already injected, which is not possible.'
          endif
          print "(i3,1x,'|',i30,1x,'|')",i,total_particles_injected(i)
       enddo
       print "(1x,36('-'))"
       ! investigate the errors that could be introduced since the last
       !    fractional timestep is not known (or at least I don't think this info
       !    is stored in the full dump files)
       print "(/,a)", ' Timestep discrepancy investigation: Once the simulation restarts,'
       print "(a)", '    find the initial nbinmax and compare it to the appropriate table in'
       print "(a)", '    order to estimate the errors in the number of wind particles injected.'
       print "(a)", '    If the last nbinmax before the restart (not stored) equals the first new'
       print "(a)", '    nbinmax upon restarting, then the injected-particle exact errors are known.'
       print "(1x,48('-'))"
       do nbinmax_test=3,20
          time_test = time - 1./(2**nbinmax_test) * dtmax
          print "(2(a,i0),a,f14.10)", ' nbinmax_test = ',nbinmax_test,', 2^nbinmax_test = ',2**nbinmax_test,', time_test = ',time_test
          total_diff = 0
          do i=nskip_ptmass+1,nptmass
             j = i - nskip_ptmass ! position in wind table
             total_particles_injected_test(i) = int(wind(i_Mdot,j)*Mdot_fac * time_test / massoftype(igas))
             total_diff = total_diff + total_particles_injected(i)-total_particles_injected_test(i)
             print "(i3,1x,'|',i30,1x,'|',i10,1x,'|')",i,total_particles_injected_test(i),total_particles_injected(i)-total_particles_injected_test(i)
          enddo
          print "(a,i0)", ' total missing particles from timestep discrepancy is: ',total_diff
          print "(1x,48('-'))"
          if (total_diff == 0) exit !finer timesteps will also show no errors in total_particles_injected, so cease this calculation
       enddo
    endif
 endif

 !
 ! set more global wind-injection properties
 !
 temp_inject = 1.e4
 print "(/,a,es14.6,a)", ' injected wind particles will have a temperature of temp_inject =',temp_inject,'K'
 if (gamma<=1.) call fatal('inject','require gamma > 1 for wind injection')
 ! convert from temperature to thermal energy
 !    P/rho = kT/(mu m_H) = (gam-1)*u
 uu_inject = temp_inject * kb_on_mh / unit_velocity**2 / (gmw*(gamma-1.))
 if (gmw<1) then
    print "(a,es21.14,a,g0)", ' init_inject: uu_inject =',uu_inject,' based on gmw = 0',gmw
 else
    print "(a,es21.14,a,g0)", ' init_inject: uu_inject =',uu_inject,' based on gmw = ',gmw
 endif

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection from the stars.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use io,        only:fatal,iverbose
 use part,      only:massoftype,igas,ihacc,i_tlast,iwindorig,isdead_or_accreted,iphase,iunknown,eos_vars,imu,itemp
 use partinject,only:add_or_update_particle,updated_particle
 use physcon,   only:pi,solarm,seconds,years,km,kb_on_mH
 use random,    only:ran2
 use eos,       only:gmw
 use options,   only:use_var_comp
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real :: r2,Minject,Mdot_code,tlast
 real :: xyzi(3),vxyz(3),xyz_star(3),vxyz_star(3),dir(3)
 real :: rr,phi,theta,cosphi,sinphi,costheta,sintheta
 real :: deltat,h,u,vinject
 real :: uu_inject_withgmw,mui
 integer :: i,j,k,i_part,ninject
 !integer :: ninject_actual

 !print*,'init: tpi = ',total_particles_injected(1:nptmass)

 ! flag to signal to partinject-->update_injected_particles(...) that a particle
 !    was updated (rather than added), which necessitates twas being computed
 updated_particle=.false.

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

 ! for reusing particles
 i_part = 1

 if (use_var_comp) uu_inject_withgmw = uu_inject

!
! loop over all wind particles
!
 !!$omp parallel do default(none) &
 !!$omp shared(nptmass)
 do i=nskip_ptmass+1,nptmass
    ! calculate how much mass to inject based on time interval since last injection
    j = i - nskip_ptmass ! position in wind table
    Mdot_code = wind(i_Mdot,j)*Mdot_fac
    vinject   = wind(i_vel,j)*vel_fac
    tlast     = xyzmh_ptmass(i_tlast,i)
    deltat    = time - tlast
    Minject   = Mdot_code*time
    ! divide by mass of gas particles
    ninject = int(Minject/massoftype(igas))-total_particles_injected(i)
    if (iverbose >= 2) print*,' point mass ',i,j,' injecting ',&
                       ninject,Minject-total_particles_injected(i)*massoftype(igas),massoftype(igas),time,tlast

    !! track when wind particles are actually injected, excluding particles that
    !!   are not actually injected due to the star being beyond outer_boundary
    !ninject_actual = 0

    if (ninject > 0) then
       ! extract current position of star
       xyz_star  = xyzmh_ptmass(1:3,i)
       ! only inject particles that won't be immediately killed due to being beyond outer_boundary
       ! this can be used with the original all-particles-are-new method as well
       if (sqrt(xyz_star(1)**2+xyz_star(2)**2+xyz_star(3)**2) < outer_boundary) then
          ! extract current velocity and injection radius of star
          rr        = 1.0001*xyzmh_ptmass(ihacc,i)
          vxyz_star = vxyz_ptmass(1:3,i)

          !set abundance of new particle based on origin-star's abundances, which is used in u and mu->eos_vars
          if (use_var_comp) then
             mui = wind(i_gmw,j) !mu is based on the injecting star's spectral type
             !print*, 'j, mui =',j,mui
          endif

          do k=1,ninject
             ! get random position on sphere
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

             !set abundance of new particle based on origin-star's abundances, which is used in u and mu->eos_vars
             if (use_var_comp) then
                u = uu_inject_withgmw * gmw/mui
             else
                u = uu_inject
             endif

             !! original method -- all particles are new
             !i_part = npart + 1 ! all particles are new

             ! reuse particles
             ! find first dead or accreted particle
             do while((.not.isdead_or_accreted(xyzh(4,i_part))) .and. i_part<npart+1)
                i_part=i_part+1
             enddo

             !! track particles that are actually injected
             !ninject_actual = ninject_actual+1

             ! add or update the particle via built-in method within partinject
             call add_or_update_particle(igas, xyzi, vxyz, h, u, i_part, npart, npartoftype, xyzh, vxyzu)
             ! point mass from where this wind particle originated
             !    not the star indexing -- the full point-mass indexing, which includes non-stars as well
             iwindorig(i_part) = i
             !composition for this particular particle based on its origin-star's abundances
             if (use_var_comp) then
                eos_vars(imu,i_part) = mui
                eos_vars(itemp,i_part) = temp_inject
                !if (mod(i_part,101)==0) then
                !   print "(2(a,i0),a,f0.2,3(a,es14.6))", 'inject_particles: eos_vars(imu,i_part) = eos_vars(',imu,',',i_part,') = ',eos_vars(imu,i_part), &
                !                                         ', eos_vars(itemp,i_part) = ',eos_vars(itemp,i_part),', u = ',u,', u*mu = ',u*mui
                !endif
             endif

             !note: add_or_update_particle increased npart by 1 if a new particle was added, so now the
             !         comparison is with npart (whereas above the comparison was with npart+1)
             if (i_part<npart) then
                ! particle was updated, not added; needed for partinject-->update_injected_particles(...)
                updated_particle=.true.
                ! flag this particle to update its timestep -- this overrides
                !    "call set_particle_type(particle_number,itype)" in partinject-->add_or_update_particle
                iphase(i_part) = iunknown
                !call set_particle_type(i_part,iunknown) !alternative/equivalent to above line
                ! begin the search for the next accreted or dead particle to reuse with the next particle
                i_part = i_part + 1
             endif
             !else (not needed, but kept here for the following comment)
             ! do not inject the particle, but keep track of it via total_particles_injected
             !    in case this star comes back within outer_boundary and then starts injecting particles
             !endif
          enddo
       endif
       ! update tlast to the current time
       ! tlast is only updated for a particular star that should have injected particles this timestep;
       !    this way, fractional particles-per-timestep can accumulate and eventually
       !    inject a particle, making Mdot more accurate
       xyzmh_ptmass(i_tlast,i) = time
       ! update total particles injected for this star
       ! this tally includes particles that would have been injected had their star
       !    not been beyond outer_boundary; this is because these not-actually-injected particles
       !    need to be tracked in case the star orbits back within outer_boundary and starts
       !    actually injecting particles
       total_particles_injected(i) = total_particles_injected(i) + ninject
       !! track particles that are actually injected
       !total_particles_injected_actual(i) = total_particles_injected_actual(i) + ninject_actual
    endif
 enddo
 if (iverbose >= 2) then
    print*,'npart = ',npart
    print*,'tpi = ',total_particles_injected(1:nptmass)
    !print*,'tpi_actual = ',total_particles_injected_actual(1:nptmass)
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
 use timestep,     only:time
 use options,      only:use_var_comp
 integer, intent(in) :: iunit
 integer             :: iunit_tpi

 if (use_var_comp) then
    call write_inopt(trim(datafile_mhn),'datafile_mhn','filename for wind composition (mu,habund,name)',iunit)
 endif
 call write_inopt(trim(datafile_wind),'datafile_wind','filename for wind injection (m,x,y,z,vx,vy,vz)',iunit)
 call write_inopt(outer_boundary,'outer_boundary','kill gas particles outside this radius',iunit)

 if (.not.isnan(time)) then !prevents issues when a simulation starts while debugging flags (which are stricter than normal flags) are turned on
    if (time>tiny(time)) then
       !
       ! write new entry in total_particles_injected.dat when each full dump is written,
       !    which is needed for eliminating discrepancies in injected-particle numbers when restarting
       !
       open(newunit=iunit_tpi,file='total_particles_injected.dat',form='formatted',position='append')
       write(iunit_tpi,*) time,nptmass,total_particles_injected(1:nptmass)
       close(iunit_tpi)
    endif
 endif

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 !use io,      only:fatal,error,warning
 use physcon, only:solarm,years
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'
 integer :: nstars

 imatch  = .true.
 select case(trim(name))
 case('datafile_mhn')
    read(valstring,*,iostat=ierr) datafile_mhn
    call read_use_var_comp_data(datafile_mhn)
 case('datafile_wind')
    read(valstring,*,iostat=ierr) datafile_wind
    call read_wind_data(datafile_wind,nstars)
    !note: nptmass=0 here, so can't compare nstars and nptmass
    !if (nstars /= nptmass) then
    !   call warning('read_options_inject','number of stars /= number of wind sources')
    !endif
    ngot = ngot + 1
 case('outer_boundary')
    read(valstring,*,iostat=ierr) outer_boundary
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
subroutine read_wind_data(filename_wind,nstars)
 use io,      only:error
 use options, only:use_var_comp
 use part,    only:n_startypes,mu_startypes,habund_startypes,name_startypes
 character(len=*), intent(in)  :: filename_wind
 integer,          intent(out) :: nstars
 integer :: iunit,ierr,i,j

 nstars = 0
 wind = 0.
 open(newunit=iunit,file=trim(filename_wind),status='old',action='read',iostat=ierr)

 if (use_var_comp) then
    !read in star's astronomy-specific ID, wind velocity, wind Mdot, and composition

    do while(ierr == 0)
       nstars = nstars + 1
       if (nstars <= maxptmass) then
          read(iunit,*,iostat=ierr) wind_int(i_astroID,nstars),wind(i_vel,nstars),wind(i_Mdot,nstars),wind_char(i_comp_char,nstars),wind_char(i_fulltype_char,nstars)
          if (ierr /= 0) then
             !entry for this wind was not read in properly, signaling the end of reading in wind entries
             !reset entries in the wind arrays that were unsuccessfully read in
             wind_int(i_astroID,nstars) = 0
             wind(i_vel,nstars) = 0.
             wind(i_Mdot,nstars) = 0.
             wind_char(i_comp_char,nstars) = ''
             wind_char(i_fulltype_char,nstars) = ''
             !decrement nstars so it corresponds to the number of stars with successfully read-in winds
             nstars = nstars - 1
          endif
       else
          call error('read_wind_data','array bounds exceeded')
       endif
    enddo

    if (nstars > 0) then
       print "(a,i0,a)", ' read wind data for ',nstars,' stars from '//trim(filename_wind)

       !correlation between wind_char and name_startypes
       do i=1,nstars
          j = 1
          !print "(a)", trim(name_startypes(j))//' -- '//trim(wind_char(i_comp_char,i))
          do while (trim(wind_char(i_comp_char,i))/=trim(name_startypes(j)) .and. j<=n_startypes)
             j = j+1
             !print "(a)", trim(name_startypes(j))//' -- '//trim(wind_char(i_comp_char,i))
          enddo
          if (j<=n_startypes) then
             !print "(a,i0)", 'star type found: wind_char(i_comp_char,i) = '//trim(wind_char(i_comp_char,i))//', name_startypes(j) = '//trim(name_startypes(j))
             wind_int(i_comp,i) = j
          else
             !print "(a)", 'star type not found: wind_char(i_comp_char,i) = '//trim(wind_char(i_comp_char,i))
             wind_int(i_comp,i) = 1
          endif
       enddo

       !entries for the unknown/init startype
       i = nstars+1
       wind_int(i_comp,i) = i_unknowninit_startypes
       wind_char(i_comp_char,i) = trim(name_startypes(wind_int(i_comp,i)))

       print "(1x,124('-'))"
       print "(a)", ' | ID           | ID       | Wind Vel | Mdot          | Composition | Composition | mu          | habund      | Spec type   |'
       print "(a)", ' | (sequential) | (astro)  | (km/s)   | (Msun/yr)     | (code ID)   | (name)      | (amu)       | (mass frac) | (full name) |'
       print "(a)", ' |--------------|----------|----------|---------------|-------------|-------------|-------------|-------------|-------------|'
       do i=1,nstars
          write(*,'(1x,"|",5x,i8,1x,"|",1x,i8,1x,"|",1x,f8.2,1x,"|",1x,es13.6,1x,"|",1x,3x,i8,1x,"|",1x,a)',advance='no') &
             i,wind_int(i_astroID,i),wind(i_vel,i),wind(i_Mdot,i),wind_int(i_comp,i),trim(wind_char(i_comp_char,i))
          do j=10,len(trim(wind_char(i_comp_char,i))),-1
             write(*,'(a)',advance='no') ' '
          enddo
          write(*,'(1x,2("|",1x,f11.9,1x),"|",1x,a)',advance='no') mu_startypes(wind_int(i_comp,i)),habund_startypes(wind_int(i_comp,i)),trim(wind_char(i_fulltype_char,i))
          do j=10,len(trim(wind_char(i_fulltype_char,i))),-1
             write(*,'(a)',advance='no') ' '
          enddo
          write(*,'(1x,"|")')
       enddo
       print "(a)", ' |--------------|-------------------------------------|-------------|-------------|-------------|-------------|-------------|'
       i = nstars+1
       write(*,'(1x,"|",5x,i8,1x,"|",1x,a,1x,"|",1x,3x,i8,1x,"|",1x,a)',advance='no') &
          i,'unknown & initialization particles ',wind_int(i_comp,i),trim(name_startypes(wind_int(i_comp,i)))
       do j=10,len(trim(wind_char(i_comp_char,i))),-1
          write(*,'(a)',advance='no') ' '
       enddo
       write(*,'(1x,2("|",1x,f11.9,1x),"| N/A         |")') mu_startypes(wind_int(i_comp,i)),habund_startypes(wind_int(i_comp,i))
       print "(1x,124('-'))"
    else
       print "(a)", ' *******************************************'
       print "(a)", ' * No wind data read in from '//trim(filename_wind)//'  *'
       print "(a)", ' * --> setting Mdot''s and vel''s to zero.   *'
       print "(a)", ' * --> setting composition indices to one. *'
       print "(a)", ' *******************************************'
       wind(i_vel ,:) = 0.
       wind(i_Mdot,:) = 0.
       wind_int(i_comp,:) = 1
    endif

 else
    !read in star's astronomy-specific ID, wind velocity, and wind Mdot --> no composition

    print "(a)", ' single composition for stars, so composition data is not being read in'
    print "(a,i0)", ' --> all stars have default composition value of wind_int(i_comp,:) = ',wind_int(i_comp,1)
    do while(ierr == 0)
       nstars = nstars + 1
       if (nstars <= maxptmass) then
          read(iunit,*,iostat=ierr) wind_int(i_astroID,nstars),wind(i_vel,nstars),wind(i_Mdot,nstars)
          if (ierr /= 0) then
             !entry for this wind was not read in properly, signaling the end of reading in wind entries
             !reset entries in the wind arrays that were unsuccessfully read in
             wind_int(i_astroID,nstars) = 0
             wind(i_vel,nstars) = 0.
             wind(i_Mdot,nstars) = 0.
             !decrement nstars so it corresponds to the number of stars with successfully read-in winds
             nstars = nstars - 1
          endif
       else
          call error('read_wind_data','array bounds exceeded')
       endif
    enddo

    if (nstars > 0) then
       print "(1x,54('-'))"
       print "(a)", ' | ID           | ID       | Wind Vel | Mdot          |'
       print "(a)", ' | (sequential) | (astro)  | (km/s)   | (Msun/yr)     |'
       print "(a)", ' |--------------|----------|----------|---------------|'
       do i=1,nstars
          print "(1x,'|',5x,i8,1x,'|',1x,i8,1x,'|',1x,f8.2,1x,'|',1x,es13.6,1x,'|')", &
             i,wind_int(i_astroID,i),wind(i_vel,i),wind(i_Mdot,i)
       enddo
       print "(1x,54('-'))"
    else
       print "(a)", ' ******************************************'
       print "(a)", ' * No wind data read in from '//trim(filename_wind)//' *'
       print "(a)", ' * --> setting Mdot''s and vel''s to zero.  *'
       print "(a)", ' ******************************************'
       wind(i_vel ,:) = 0.
       wind(i_Mdot,:) = 0.
    endif

 endif
 close(iunit)

 print*
end subroutine read_wind_data

!----------------------------------------------------------------
!+
!  read various composition data for stars from file
!+
!----------------------------------------------------------------
subroutine read_use_var_comp_data(filename_mhn)
 use io,   only:error,warning
 use part, only:n_startypes,mu_startypes,habund_startypes,name_startypes
 use eos,  only:gmw,num_var_comp
 character(len=*), intent(in) :: filename_mhn
 integer :: iunit,ierr
 integer :: i,j
 character(len=12) :: unknown_init_char
 character(len=20) :: unknown_init_name_startypes
 logical :: read_unknown_startypes=.false.,found_unknown_startypes=.false.
 real :: tol_mu = 1.d-6

 n_startypes = 0
 open(newunit=iunit,file=filename_mhn,status='old',action='read',iostat=ierr)
 if (ierr /= 0) then
    print "(2(/,a))",' ERROR opening "'//trim(filename_mhn)//'" for read of various composition data', &
                       ' -> this file should contain mu,name for each point mass, one per line'
    print "(a)", ' no startypes were read in from '//trim(filename_mhn)
    print "(a,g0)", ' --> mu for all particles will come from module eos-->gmw, which is ',gmw
    name_startypes(1) = ''
    mu_startypes(1) = gmw
    habund_startypes(1) = 0.7
    num_var_comp = 1
    i_unknowninit_startypes = 1
    close(iunit)
    call error('read_use_var_comp_data','no startypes were read in -- create '//trim(filename_mhn))
    return
 endif

 !read startype for how to treat unknown stellar-wind particles and initialization particles
 !to make the user aware of this, the first line should be "unknown_init <startype>" so it is explicitly stated
 !   that the first entry is for unknown startypes and initialization particles
 read(iunit,*,iostat=ierr) unknown_init_char,unknown_init_name_startypes
 !print "(a,i0)", '"'//trim(unknown_init_char)//'" "'//trim(unknown_init_name_startypes)//'" ',ierr
 if (ierr==0 .and. trim(unknown_init_char)=='unknown_init') then
    print "(a)", ' startype for unknown and initialization particles is startype = '//trim(unknown_init_name_startypes)
    read_unknown_startypes =.true.
 else
    print "(a)", ' error reading startype for unknown and initialization particles -- first line of '//trim(filename_mhn)//' should be "unknown_init <startype>"'
    print "(a)", ' --> setting startype for unknown and initialization particles to be the first correct entry of '//trim(filename_mhn)
    rewind(iunit)
    read(iunit,*,iostat=ierr) mu_startypes(1),habund_startypes(1),name_startypes(1)
    if (ierr==0) then
       rewind(iunit) !first entry is likely a mu/startype combination, so it seems like the "unknown_init" line was skipped
    else
       ierr = 0 !reset ierr in hopes that the rest of the file is correct even though the "unknown_init" line was not correct
    endif
    call error('read_use_var_comp_data','startype for unknown and initialization particles not found -- fix '//trim(filename_mhn))
 endif

 !read mu and startype entries for each type of stellar wind; check for uniqueness
 do while(ierr==0)
    n_startypes = n_startypes + 1
    if (n_startypes > size(mu_startypes)) then
       ierr = 66
    else
       read(iunit,*,iostat=ierr) mu_startypes(n_startypes),habund_startypes(n_startypes),name_startypes(n_startypes)
       if (ierr==0) then
          !new startype entry was read in
          do j=1,n_startypes-1
             !check for uniqueness among name_startypes
             if (trim(name_startypes(n_startypes))==trim(name_startypes(j))) then
                print "(a)", ' *****************************************************************************'
                print "(2(a,i0),a)", ' * duplicate name entry found in '//trim(filename_mhn)//': entry ',n_startypes,' is the same as entry ',j
                print "(a,i0,a,g0,a)", ' * entry ',n_startypes,': mu = ',mu_startypes(n_startypes),' and name = '//name_startypes(n_startypes)
                print "(a,i0,a,g0,a)", ' * entry ',j,': mu = ',mu_startypes(j),' and name = '//name_startypes(j)
                print "(a)", ' *****************************************************************************'
                call error('read_use_var_comp_data','duplication of name for startypes -- fix '//trim(filename_mhn))
             endif
             !check for uniqueness among mu_startypes -- not always an error, but could be, so make a warning
             if (abs(mu_startypes(n_startypes)-mu_startypes(j))<tol_mu) then
                print "(a)", ' ***************************************************************************'
                print "(2(a,i0),a)", ' * duplicate mu entry found in '//trim(filename_mhn)//': entry ',n_startypes,' is the same as entry ',j
                print "(a,i0,a,g0,a)", ' * entry ',n_startypes,': mu = ',mu_startypes(n_startypes),' and name = '//name_startypes(n_startypes)
                print "(a,i0,a,g0,a)", ' * entry ',j,': mu = ',mu_startypes(j),' and name = '//name_startypes(j)
                print "(a)", ' ***************************************************************************'
                call warning('read_use_var_comp_data','duplication of mu for startypes -- might need to  fix '//trim(filename_mhn))
             endif
          enddo
       else
          !entry for this startype was not read in properly, signaling the end of reading in startype entries
          !reset entries in the startype arrays that were unsuccessfully read in
          mu_startypes(n_startypes) = 0.
          habund_startypes(n_startypes) = 0.
          name_startypes(n_startypes) = ''
          !decrement startypes so it corresponds to the number of startypes successfully read-in
          n_startypes = n_startypes - 1
       endif
    endif
 enddo
 print "(a,i0,a)",' read ',n_startypes,' various compositions for point masses from '//trim(filename_mhn)
 if (ierr==66) then
    call error('read_use_var_comp_data','array size exceeded in read_use_var_comp_data, increase size of mu_startypes',var='n_startypes',ival=n_startypes+1)
 endif
 close(iunit)

 if (n_startypes>0) then
    !at least one startype was read in
    if (n_startypes==1) then
       !one startype was read in, so that will become the unknown type
       mu_startypes(n_startypes+1) = mu_startypes(1)
       habund_startypes(n_startypes+1) = habund_startypes(1)
       name_startypes(n_startypes+1) = name_startypes(1)
       i_unknowninit_startypes = 1
    else
       !multiple startypes were read-in, so determining the unknown type is more involved
       if (read_unknown_startypes) then
          !correlate the unknown startype with the read-in startypes
          do i=1,n_startypes
             if (index(trim(unknown_init_name_startypes),trim(name_startypes(i)))==1 .and. &
                 index(trim(name_startypes(i)),trim(unknown_init_name_startypes))==1) then
                if (.not.found_unknown_startypes) then
                   mu_startypes(n_startypes+1) = mu_startypes(i)
                   habund_startypes(n_startypes+1) = habund_startypes(i)
                   name_startypes(n_startypes+1) = name_startypes(i)
                   i_unknowninit_startypes = i
                   found_unknown_startypes = .true.
                else
                   print "(a)", ' ERROR: multiple startypes match the startype for unknown and initialization particles -- rewrite '//trim(filename_mhn)//' to make sure this isn''t happening'
                endif
             endif
          enddo
       endif
       if (.not.found_unknown_startypes) then
          !unknown startype was not found --> set to the first startype
          if (read_unknown_startypes) then
             print "(a)", ' startype for unknown and initialization particles was not found among the read-in startype entries'
          else
             print "(a)", ' startype for unknown and initialization particles was not read in'
          endif
          print "(a)", ' --> setting startype for unknown and initialization particles to be the first read-in startype'
          mu_startypes(n_startypes+1) = mu_startypes(1)
          habund_startypes(n_startypes+1) = habund_startypes(1)
          name_startypes(n_startypes+1) = name_startypes(1)
          i_unknowninit_startypes = 1
       endif
    endif

    print "(1x,56('-'))"
    print "(a)", ' | n_startype | mu          | habund      | name        |'
    print "(a)", ' |            | (amu)       | (mass frac) |             |'
    print "(a)", ' |------------|-------------|-------------|-------------|'
    do i=1,n_startypes
       write(*,'(1x,"|",1x,i10,1x,"|",2(1x,f11.9,1x,"|"),1x,a)',advance='no') i,mu_startypes(i),habund_startypes(i),'"'//trim(name_startypes(i))//'"'
       do j=8,len(trim(name_startypes(i))),-1
          write(*,'(a)',advance='no') ' '
       enddo
       if (abs(mu_startypes(n_startypes+1)-mu_startypes(i))<tol_mu .and. &
           abs(habund_startypes(n_startypes+1)-habund_startypes(i))<tol_mu .and. &
           trim(name_startypes(n_startypes+1))==trim(name_startypes(i))) then
          write(*,'(1x,"|",a)')  '  <-- unknown/init startype'
       else
          write(*,'(1x,"|")')
       endif
    enddo
    print "(1x,56('-'))"
    num_var_comp = n_startypes
    print "(a,i0)", ' num_var_comp has been set to n_startypes -- num_var_comp = ',num_var_comp
 else
    print "(a)", ' no startypes were read in from '//trim(filename_mhn)
    print "(a,g0)", ' --> mu for all particles will come from module eos-->gmw, which is ',gmw
    print "(a,g0)", ' --> habund for all particles will be 0.7'
    mu_startypes(1) = gmw
    habund_startypes(1) = 0.7
    name_startypes(1) = ''
    num_var_comp = 1
    i_unknowninit_startypes = 1
    call error('read_use_var_comp_data','no startypes were read in -- fix '//trim(filename_mhn))
 endif

 ! end of file error is OK
 if (ierr < 0) ierr = 0

 print*
end subroutine read_use_var_comp_data

subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
end subroutine set_default_options_inject

end module inject
