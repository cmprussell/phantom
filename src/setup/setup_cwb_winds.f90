!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for simulations of Colliding Wind Binaries (CWBs)
!    Done by Christopher Russell via adapting the Galactic Center setup, which was
!    Adapted by Daniel Price in collaboration with Jorge Cuadra
!
! :References: Madura et al. (2013), MNRAS, 436, 4, 3820
!              Russell et al. (2016), MNRAS, 458, 3, 2275
!
! :Owner: Christopher Russell
!
! :Runtime parameters:
!   - datafile_mhn       : *filename for various-composition data (mu,habund,name)*
!   - datafile_mpv       : *filename for star data (m,x,y,z,vx,vy,vz)*
!   - m_gas              : *gas mass resolution in solar masses*
!   - ninjectmax         : *max number of particles to inject in a single timestep per star (0=no limit)*
!   - use_var_comp_local : *whether or not to use various compositions*
!
! :Dependencies: cooling, cooling_solver, datafiles, dim, eos,
!   externalforces, infile_utils, inject, io, options, part, physcon,
!   prompting, spherical, timestep, units
!
 implicit none
 public :: setpart

 !
 ! setup options and default values for these
 !
 character(len=120) :: datafile_mhn = 'mu_habund_name.txt'
 character(len=120) :: datafile_mpv = 'stars_mass_pos_vel.txt'
 real :: m_gas = 1.e-10
 integer :: ninjectmax = 0
 logical :: use_var_comp_local=.false.

contains

!----------------------------------------------------------------
!+
!  setup for colliding wind binaries
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas,iwindorig,eos_vars,imu,mu_startypes,n_startypes
 use units,     only:set_units,umass,udist,utime
 use physcon,   only:solarm,kpc,pi,au,solarr
 use io,        only:fatal,iprint,master
 use eos,       only:gmw,use_var_comp
 use timestep,  only:dtmax,tmax
 use spherical, only:set_sphere
 use datafiles, only:find_phantom_datafile
 use inject,    only:ninjectmax_cwb
 use options,        only:icooling,nfulldump,iexternalforce
 use cooling_solver, only:icool_method,lambda_table
 use cooling,        only:Tfloor
 use inject,         only:read_use_var_comp_data
 use externalforces, only:iext_windaccel
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=len(fileprefix)+6) :: setupfile
 character(len=len(datafile_mpv)) :: filename_mpv
 character(len=len(datafile_mhn)) :: filename_mhn
 integer :: ierr,i,iexist
 real    :: scale,psep
 integer :: j,nTooClose
 real    :: dis1,dis2,dis,vinit
 integer :: iclo

 print "(/,a)", ' start setpart(...)'
!
!set some CWB-specific parameters in case cwb.in does not exist
!
 inquire(file=trim(fileprefix)//'.in',exist=iexist)
 if (.not.iexist) then
    print "(a)", ' the input file '//trim(fileprefix)//'.in does not exist, so some convenient CWB-specific values for '//trim(fileprefix)//'.in are being specified at the beginning of setup/setup_cwb_winds.f90 --> setpart(...)'
    tmax      = 0.2
    dtmax     = 0.1
    nfulldump = 1
    lambda_table = 1 !use Lambda(T) cooling table for radiative cooling
    iexternalforce = iext_windaccel
 endif
!
! units (mass = solar mass, length = AU, time set so that G=1)
!
 scale = au
 call set_units(mass=solarm,dist=scale,G=1.d0)
!
! general parameters
!
 time = 0.
 hfact = 1.2
 polyk = 0.
 gamma = 5./3.
 gmw = 0.6  ! completely ionized, solar abu; eventually needs to be WR abu
 !
 ! read setup parameters from the .setup file
 ! if file does not exist, then ask for user input
 !
 setupfile = trim(fileprefix)//'.setup'
 call read_setupfile(setupfile,iprint,ierr)
 if (ierr /= 0 .and. id==master) then
    call interactive_setup()               ! read setup options from user
    call write_setupfile(setupfile,iprint) ! write .setup file with defaults
 endif
!
! space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 massoftype = m_gas*(solarm/umass)  ! mass resolution

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
!
! Read positions, masses and velocities of stars from file
!
 !print*, 'filename_mpv =',trim(filename_mpv)
 filename_mpv=find_phantom_datafile(datafile_mpv,'cwb')
 print "(a)", ' filename_mpv = '//trim(filename_mpv)
 !print "(a,i0)", ' nptmass E = ',nptmass !note: nptmass=0 here
 call read_ptmass_data(filename_mpv,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr) !nptmass is determined in this subroutine
 do i=1,nptmass
    xyzmh_ptmass(1:3,i) = xyzmh_ptmass(1:3,i)
    xyzmh_ptmass(4,i)   = xyzmh_ptmass(4,i)
    !xyzmh_ptmass(ihacc,i)  = h_sink*solarr/au
    !xyzmh_ptmass(ihsoft,i) = h_sink*solarr/au
    xyzmh_ptmass(ihacc,i)  = 0.1**solarr/au !for now, only for removing gas particles too close to stars; set later in init_inject()
    !xyzmh_ptmass(ihsoft,i) = 0.1*solarr/au
    xyzmh_ptmass(5,i) = 0.1*solarr/au !for now, only for removing gas particles too close to stars; set in init_inject()
    vxyz_ptmass(1:3,i)  = vxyz_ptmass(1:3,i)*1.e5*utime/udist
 enddo
 print "(a,i0)", ' nptmass = ',nptmass
 do i=1,nptmass
    print "(a,i0)", ' ptmass i = ',i
    print*, 'xyzmh_ptmass(:,i) =',xyzmh_ptmass(:,i)
    print*, 'vxyz_ptmass(:,i) =',vxyz_ptmass(:,i)
 enddo
 !print*, 'h_sink (Rsun) =',h_sink
 !print*, 'h_sink (AU) =',h_sink*solarr/au
 print*, 'umass, udist, utime =',umass,udist,utime
!
! setup initial sphere of particles to prevent initialisation problems
!
 psep = 1.0
 call set_sphere('cubic',id,master,0.,20.,psep,hfact,npart,xyzh)
!
! initialize mean molecular weight if needed
!
 print "(/,a,l)", ' setpart: use_var_comp =',use_var_comp
 if (use_var_comp) then
    !if cwb.in does not already exist, then the following call to read_use_var_comp_data() will be run
    if (n_startypes<=0) then !if already called from inject_cwb_winds()-->read_options_inject()-->read_use_var_comp_data(), then don't call this again
       filename_mhn = find_phantom_datafile(datafile_mhn,'cwb')
       print '(a)', ' filename_mhn = '//trim(filename_mhn)
       print "(a)", ' calling read_use_var_comp_data(filename_mhn) from setpart()'
       call read_use_var_comp_data(filename_mhn)
    endif
    if (n_startypes>=1) then
       gmw = mu_startypes(n_startypes+1)
    else
       gmw = 0.6 !completely ionized solar-abundance gas
    endif
    eos_vars(imu,:) = gmw
 else
    gmw = 0.6 !completely ionized solar-abundance gas
 endif
 if (gmw<1.) then
    print "(a,g0)", ' setpart: gmw = 0',gmw
 else
    print "(a,g0)", ' setpart: gmw = ',gmw
 endif
!
! set temperature of initial particles
! initialize internal energy based on 1. temperature -- either read from a file or 1e4K
!                                     2. mean molecular weight -- gmw
!
 !this commented-out code is useful for radiative cooling tests since it is easy to change the temperature of the initialization particles
 !this code is not needed for full GC sims, but keep this for future coolig tests
 !open(unit=61,file='Tinit.dat',form='formatted',status='old',iostat=ierr_Tinit)
 !if (ierr_Tinit==0) then
 !   read(61,*,iostat=ierr_Tinit) vxyzu(4,1)
 !   if (ierr_Tinit==0) then
 !      !vxyzu(4,1)=5.356136065348470d-4 * vxyzu(4,1)/1.d4
 !      vxyzu(4,1) = 5.356136065348470d-4 * 0.6d0/gmw * vxyzu(4,1)/1.d4
 !      vxyzu(4,:) = vxyzu(4,1)
 !      print "(a,g0,a,/)", ' Tinit.dat entry successfully found -- initial particle energy set to T = ',vxyzu(4,1)/5.356136065348470d-4*gmw/0.6d0*1.d4,' K'
 !   endif
 !endif
 !if (ierr_Tinit/=0) then
 !   !vxyzu(4,:) = 5.317e-4 ! T_init=1e4
 !   !vxyzu(4,:) = 5.317e-4 * 1.e2 ! T_init=1e6
 !   !vxyzu(4,:) = 5.356136065348470d-4 ! T_init=1e4K to more accuracy
 !   vxyzu(4,:) = 5.356136065348470d-4 * 0.6d0/gmw ! T_init=1e4K to more accuracy
 !   !vxyzu(4,:) = 5.356136065348470d-4 * 1.d1 ! T_init=1e5K to more accuracy
 !   !vxyzu(4,:) = 5.356136065348470d-4 * 1.d2 ! T_init=1e6K to more accuracy
 !   !vxyzu(4,:) = 5.356136065348470d-4 * 1.d3 ! T_init=1e7K to more accuracy
 !   !vxyzu(4,:) = 5.356136065348470d-4 * 1.d4 ! T_init=1e8K to more accuracy
 !   !print "(a,g0,a)", ' Tinit.dat not found or not configured correctly -- initial particle energy set to T = ',vxyzu(4,1)/5.356136065348470d-4*1.d4,' K'
 !   print "(a,g0,a,/)", ' Tinit.dat not found or not configured correctly -- initial particle energy set to T = ',vxyzu(4,1)/5.356136065348470d-4*gmw/0.6d0*1.d4,' K'
 !endif
 !close(61)
 vxyzu(4,:) = 5.356136065348470d-4 * 0.6d0/gmw ! T_init=1e4K to more accuracy
 print "(a,g0,a)", ' initial particle energy set so that T_init = ',vxyzu(4,1)/5.356136065348470d-4*gmw/0.6d0*1.d4,' K'

 npartoftype(igas) = npart

 ! Radial flow away from nearest star
 do i=1,npart
    dis1=(xyzh(1,i)-xyzmh_ptmass(1,1))**2+(xyzh(2,i)-xyzmh_ptmass(2,1))**2+(xyzh(3,i)-xyzmh_ptmass(3,1))**2
    dis2=(xyzh(1,i)-xyzmh_ptmass(1,2))**2+(xyzh(2,i)-xyzmh_ptmass(2,2))**2+(xyzh(3,i)-xyzmh_ptmass(3,2))**2
    iclo=1
    if (dis2<dis1) iclo=2
    dis=sqrt(min(dis1,dis2))
    vinit=1000*1.e5/(udist/utime)
    !vxyzu(1,i)=(xyzh(1,i)-xyzmh_ptmass(1,iclo))/dis*vinit
    vxyzu(1:3,i)=(xyzh(1:3,i)-xyzmh_ptmass(1:3,iclo))/dis*vinit
 enddo

 ! Remove particles near stars
 nTooClose=0
 do i=1,npart
    do j=1,nptmass
       if (sqrt((xyzh(1,i)-xyzmh_ptmass(1,j))**2+(xyzh(2,i)-xyzmh_ptmass(2,j))**2+(xyzh(3,i)-xyzmh_ptmass(3,j))**2) < 3.*xyzmh_ptmass(5,j)) then
          !xyzh(4,i) = -xyzh(4,i) !causes an error since the setup procedure shouldn't yield any particle with h<0
          xyzh(3,i) = 100.+nTooClose !instead, move the particle way out of the domain -- eventually the simulation boundary should be incorporated into this computation
          nTooClose=nTooClose+1
          print*, i,j,xyzh(1,i),xyzmh_ptmass(1,j),xyzh(2,i),xyzmh_ptmass(2,j),xyzh(3,i),xyzmh_ptmass(3,j),sqrt((xyzh(1,i)-xyzmh_ptmass(1,j))**2+(xyzh(2,i)-xyzmh_ptmass(2,j))**2+(xyzh(3,i)-xyzmh_ptmass(3,j))**2),xyzmh_ptmass(5,j)
       endif
    enddo
 enddo
 print '(a,i0)', ' Particle Removal near stars: nTooClose = ',nTooClose
 !npartoftype(igas) = npart-nTooClose

 ninjectmax_cwb = ninjectmax
 print '(2(a,i0))', ' ninjectmax, ninjectmax_cwb = ',ninjectmax,', ',ninjectmax_cwb

 !denote these particles (and all not-yet-injected particles) as initial particles
 !later this variable will be the star from which each wind particle originated
 iwindorig = 0

 !ensure EIS is used for cooling
 icooling = 1
 icool_method = 2
 Tfloor = 1.d4

 if (nptmass == 0) call fatal('setup','no particles setup')
 if (ierr /= 0) call fatal('setup','ERROR during setup')

 print*
end subroutine setpart

!----------------------------------------------------------------
!+
!  read sink particle masses, positions and velocities from file
!+
!----------------------------------------------------------------
subroutine read_ptmass_data(filename,xyzmh_ptmass,vxyz_ptmass,n,ierr)
 use io, only:error
 character(len=*), intent(in) :: filename
 real,    intent(out)   :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: n
 integer, intent(out)   :: ierr
 integer :: iunit,n_input,i
 logical :: fulloutput=.false.

 !print "(a,i0)", ' n A = ',n !note: n=0 here for CWBs (since all point masses are stars)
 n_input = n
 open(newunit=iunit,file=filename,status='old',action='read',iostat=ierr)
 if (ierr /= 0) then
    print "(/,2(a,/))",' ERROR opening "'//trim(filename)//'" for read of point mass data', &
                       ' -> this file should contain m,x,y,z,vx,vy,vz for each point mass, one per line'
 endif
 if (fulloutput) then
    print "(1x,175('-'))"
    print "(a)", ' |  ID |                  mass |' // &
                         '                     x |                     y |                     z |' // &
                         '                    vx |                    vy |                    vz |'
    print "(1x,'|',5('-'),'|',7(23('-'),'|'))"
 else
    print "(1x,133('-'))"
    print "(a)", ' |  ID |            mass |' // &
                         '               x |               y |               z |' // &
                         '              vx |              vy |              vz |'
    print "(1x,'|',5('-'),'|',7(17('-'),'|'))"
 endif
 do while(ierr==0)
    n = n + 1
    if (n > size(xyzmh_ptmass(1,:))) then
       ierr = 66
    else
       read(iunit,*,iostat=ierr) xyzmh_ptmass(4,n),xyzmh_ptmass(1:3,n),vxyz_ptmass(1:3,n)
    endif
    if (ierr /= 0) n = n - 1
 enddo
 !do i=n_input,n
 do i=n_input+1,n
    if (fulloutput) then
       print "(1x,'|',1x,i3,1x,'|',7(es22.14,1x,'|'))", i,xyzmh_ptmass(4,i),xyzmh_ptmass(1:3,i),vxyz_ptmass(1:3,i)
    else
       print "(1x,'|',1x,i3,1x,'|',7(es16.8 ,1x,'|'))", i,xyzmh_ptmass(4,i),xyzmh_ptmass(1:3,i),vxyz_ptmass(1:3,i)
    endif
 enddo
 if (fulloutput) then
    print "(1x,175('-'))"
 else
    print "(1x,133('-'))"
 endif
 if (n-n_input==1) then
    print "(a,i0,a)",' read mass, position, and velocity for ',n - n_input,' point mass from '//trim(filename)
 else
    print "(a,i0,a)",' read mass, position, and velocity for ',n - n_input,' point masses from '//trim(filename)
 endif
 print "(a)",' Note: The above mass, position, and velocity values are in code units.'
 if (ierr==66) then
    call error('read_ptmass_data','array size exceeded in read_ptmass_data, recompile with MAXPTMASS=n',var='n',ival=n+1)
 endif

 ! end of file error is OK
 if (ierr < 0) ierr = 0
 !print "(a,i0)", ' n B = ',n

 print*
end subroutine read_ptmass_data

!------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!------------------------------------------
subroutine write_setupfile(filename,iprint)
 use infile_utils, only:write_inopt
 use dim,          only:tagline
 character(len=*), intent(in) :: filename
 integer,          intent(in) :: iprint
 integer                      :: lu,ierr1,ierr2

 print "(a)", ' Writing '//trim(filename)//' with setup options'
 open(newunit=lu,file=filename,status='replace',form='formatted')
 write(lu,"(a)") '# '//trim(tagline)
 write(lu,"(a)") '# input file for Phantom colliding wind binary setup'

 write(lu,"(/,a)") '# datafile for stars'
 call write_inopt(datafile_mpv,'datafile_mpv','filename for star data (m,x,y,z,vx,vy,vz)',lu,ierr1)

 write(lu,"(/,a)") '# resolution'
 call write_inopt(m_gas, 'm_gas','gas mass resolution in solar masses',lu,ierr2)
 call write_inopt(ninjectmax, 'ninjectmax','max number of particles to inject in a single timestep per star (0=no limit)',lu,ierr2)

 write(lu,"(/,a)") '# use various compositions'
 call write_inopt(use_var_comp_local, 'use_var_comp_local','whether or not to use various compositions',lu,ierr2)
 if (use_var_comp_local) then
    write(lu,"(/,a)") '# datafile for mu and stellar classification names'
    call write_inopt(datafile_mhn,'datafile_mhn','filename for various-composition data (mu,habund,name)',lu,ierr1)
 endif

 close(lu)

end subroutine write_setupfile

!------------------------------------------
!+
!  Read setup parameters from input file
!+
!------------------------------------------
subroutine read_setupfile(filename,iprint,ierr)
 use infile_utils, only:open_db_from_file,inopts,close_db,read_inopt
 use dim,          only:maxvxyzu
 use eos,          only:use_var_comp
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(in)  :: iprint
 integer,          intent(out) :: ierr
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 print "(/,a)", ' start read_setupfile(filename,...) = read_setupfile('//trim(filename)//')'
 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) then
    print "(/,1x,a,/)", 'END read_setupfile -- '//trim(filename)//' not read properly, so will try interactive mode'
    return
 endif

 print "(a)", ' These are the values read-in from '//trim(filename)//' that will be used for the ensuing calculation:'
 nerr = 0
 call read_inopt(datafile_mpv,'datafile_mpv',db,errcount=nerr)
 print "(a)", ' datafile_mpv = '//trim(datafile_mpv)
 call read_inopt(m_gas,'m_gas',db,errcount=nerr)
 print "(a,es21.14)", ' m_gas =',m_gas
 call read_inopt(ninjectmax,'ninjectmax',db,errcount=nerr)
 print "(a,i0)", ' ninjectmax = ',ninjectmax

 call read_inopt(use_var_comp_local,'use_var_comp_local',db,errcount=nerr)
 print "(a,l)", ' use_var_comp_local =',use_var_comp_local
 use_var_comp = use_var_comp_local
 if (use_var_comp) then
    call read_inopt(datafile_mhn,'datafile_mhn',db,errcount=nerr)
    print "(a)", ' datafile_mhn = '//trim(datafile_mhn)
 endif

 if (nerr > 0) then
    print "(1x,a,i2,a)",'Setup_cwb: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif
 call close_db(db)

 print *
end subroutine read_setupfile

!------------------------------------------
!+
!  Prompt user for setup options
!+
!------------------------------------------
subroutine interactive_setup()
 use prompting, only:prompt
 use eos,       only:use_var_comp

 print "(2(/,a),/)",'*** Welcome to your friendly colliding-wind binary setup,',&
                    '    where the stars are so bright that their winds sail on starlight. ***'
 call prompt('Enter filename for mass, position, and velocity of star data',datafile_mpv,noblank=.true.)
 call prompt('Enter mass resolution of injected gas particles in Msun',m_gas,1.e-15,1.)
 call prompt('Enter max number of particles to inject in a single timestep per star (0=no limit)',ninjectmax,0,1000000000)
 call prompt('Enter logical value for use_var_comp',use_var_comp_local)
 use_var_comp = use_var_comp_local
 if (use_var_comp_local) then
    call prompt('Enter filename for various compositions data',datafile_mhn,noblank=.true.)
 endif

 print*
end subroutine interactive_setup

end module setup
