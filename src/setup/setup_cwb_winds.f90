!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
!!!! Setup for simulations of the Galactic Centre
!!!!    Adapted by Daniel Price in collaboration with Jorge Cuadra
! Setup for simulations of Collidinw Wind Binaries (CWBs)
!    Done by Christopher Russell via adaptinge the Galactic Center setup, which was
!    Adapted by Daniel Price in collaboration with Jorge Cuadra
!
! :References: Paumard et al. (2006)
!
! :Owner: Christopher Russell
!
! :Runtime parameters:
!   - datafile : *filename for star data (m,x,y,z,vx,vy,vz)*
!   - h_sink   : *sink particle radii in Rsun*
!   - m_gas    : *gas mass resolution in solar masses*
!   - ninjectmax : *max number of particles to inject in a single timestep per star (0=no limit)*
!
! :Dependencies: datafiles, dim, eos, infile_utils, io, part, physcon,
!   prompting, spherical, timestep, units
!
 implicit none
 public :: setpart

 !
 ! setup options and default values for these
 !
 character(len=120) :: datafile = 'stars_mass_pos_vel.txt'
 real :: m_gas = 1.e-6 ! gas mass resolution in Msun
 real :: h_sink = 10. ! sink particle radii in Rsun
 INTEGER :: ninjectmax = 0 ! max number of particles to inject in a single timestep per star (0=no limit)

 private

contains

!----------------------------------------------------------------
!+
!  setup for colliding wind binaries
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas
 use units,     only:set_units,umass,udist,utime
 use physcon,   only:solarm,kpc,pi,au,solarr
 use io,        only:fatal,iprint,master
 use eos,       only:gmw
 use timestep,  only:dtmax
 use spherical, only:set_sphere
 use datafiles, only:find_phantom_datafile
 USE inject,    ONLY:ninjectmax_cwb
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
 character(len=len(datafile)) :: filename
 integer :: ierr,i
 real    :: scale,psep
 INTEGER :: j,nTooClose
 REAL    :: dis1,dis2,dis,vinit
 INTEGER :: iclo
!
! units (mass = mass of black hole, length = 1 arcsec at 8kpc)
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
 dtmax = 0.01
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
WRITE(*,*) 'filename =',TRIM(filename)
 filename=find_phantom_datafile(datafile,'cwb') !will this work? does search in curent directory
WRITE(*,*) 'filename =',TRIM(filename)
WRITE(*,'(A,I0)') 'nptmass E = ',nptmass
 call read_ptmass_data(filename,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr) !nptmass is determined in this subroutine
 do i=1,nptmass
    xyzmh_ptmass(1:3,i)  = xyzmh_ptmass(1:3,i)
    xyzmh_ptmass(4,i)    = xyzmh_ptmass(4,i)
    xyzmh_ptmass(ihacc,i)  = h_sink*solarr/au
    xyzmh_ptmass(ihsoft,i) = h_sink*solarr/au
    vxyz_ptmass(1:3,i)  = vxyz_ptmass(1:3,i)*1.e5*utime/udist
 enddo
WRITE(*,'(A,I0)') 'nptmass = ',nptmass
DO i=1,nptmass
WRITE(*,*) i
WRITE(*,*) xyzmh_ptmass(:,i)
WRITE(*,*) vxyz_ptmass(:,i)
ENDDO
WRITE(*,*) 'h_sink (Rsun) =',h_sink
WRITE(*,*) 'h_sink (AU) =',h_sink*solarr/au
WRITE(*,*) 'umass, udist, utime =',umass,udist,utime
!
! setup initial sphere of particles to prevent initialisation problems
!
 psep = 1.0
 call set_sphere('cubic',id,master,0.,20.,psep,hfact,npart,xyzh)
 vxyzu(4,:) = 5.317e-4
 npartoftype(igas) = npart
!!! Radial flow away from nearest star
 DO i=1,npart
    dis1=(xyzh(1,i)-xyzmh_ptmass(1,1))**2+(xyzh(2,i)-xyzmh_ptmass(2,1))**2+(xyzh(3,i)-xyzmh_ptmass(3,1))**2
    dis2=(xyzh(1,i)-xyzmh_ptmass(1,2))**2+(xyzh(2,i)-xyzmh_ptmass(2,2))**2+(xyzh(3,i)-xyzmh_ptmass(3,2))**2
    iclo=1
    IF(dis2<dis1) iclo=2
    dis=SQRT(MIN(dis1,dis2))
    vinit=1000*1.e5/(udist/utime)
    !vxyzu(1,i)=(xyzh(1,i)-xyzmh_ptmass(1,iclo))/dis*vinit
    vxyzu(1:3,i)=(xyzh(1:3,i)-xyzmh_ptmass(1:3,iclo))/dis*vinit
 ENDDO
!!! Radial flow away from nearest star
!!! Remove particles near star
 nTooClose=0
 DO i=1,npart
    DO j=1,nptmass
       IF(SQRT((xyzh(1,i)-xyzmh_ptmass(1,j))**2+(xyzh(2,i)-xyzmh_ptmass(2,j))**2+(xyzh(3,i)-xyzmh_ptmass(3,j))**2) < 3.*xyzmh_ptmass(5,j)) THEN
          !xyzh(4,i) = -xyzh(4,i) !causes an error since the setup procedure shouldn't yield any particle with h<0
          xyzh(3,i) = 100.+nTooClose !instead, move the particle way out of the domain -- eventually the simulation boundary should be incorporated into this computation
          nTooClose=nTooClose+1
WRITE(*,*) i,j,xyzh(1,i),xyzmh_ptmass(1,j),xyzh(2,i),xyzmh_ptmass(2,j),xyzh(3,i),xyzmh_ptmass(3,j),SQRT((xyzh(1,i)-xyzmh_ptmass(1,j))**2+(xyzh(2,i)-xyzmh_ptmass(2,j))**2+(xyzh(3,i)-xyzmh_ptmass(3,j))**2),xyzmh_ptmass(5,j)
       ENDIF
    ENDDO
 ENDDO
WRITE(*,*) 'Particle Removal near stars: nTooClose =',nTooClose
 !npartoftype(igas) = npart-nTooClose
!!! Remove particles near star

 ninjectmax_cwb = ninjectmax
WRITE(*,*) 'ninjectmax, ninjectmax_cwb =',ninjectmax,ninjectmax_cwb

 if (nptmass == 0) call fatal('setup','no particles setup')
 if (ierr /= 0) call fatal('setup','ERROR during setup')

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
 integer :: iunit,n_input

WRITE(*,*) 'n_input =',n_input
WRITE(*,*) 'n A =',n
 n_input = n
 open(newunit=iunit,file=filename,status='old',action='read',iostat=ierr)
 if (ierr /= 0) then
    print "(/,2(a,/))",' ERROR opening "'//trim(filename)//'" for read of point mass data', &
                       ' -> this file should contain m,x,y,z,vx,vy,vz for each point mass, one per line'
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
 print "(a,i4,a)",' READ',n - n_input,' point masses from '//trim(filename)
 if (ierr==66) then
    call error('read_ptmass_data','array size exceeded in read_ptmass_data, recompile with MAXPTMASS=n',var='n',ival=n+1)
 endif

 ! end of file error is OK
 if (ierr < 0) ierr = 0
WRITE(*,*) 'n B =',n

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

 write(iprint,"(a)") ' Writing '//trim(filename)//' with setup options'
 open(newunit=lu,file=filename,status='replace',form='formatted')
 write(lu,"(a)") '# '//trim(tagline)
 write(lu,"(a)") '# input file for Phantom colliding wind binary setup'

 write(lu,"(/,a)") '# datafile'
 call write_inopt(datafile,'datafile','filename for star data (m,x,y,z,vx,vy,vz)',lu,ierr1)

 write(lu,"(/,a)") '# resolution'
 call write_inopt(m_gas, 'm_gas','gas mass resolution in solar masses',lu,ierr2)
 call write_inopt(h_sink, 'h_sink','sink particle radii in Rsun',lu,ierr2)
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
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(in)  :: iprint
 integer,          intent(out) :: ierr
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 write(iprint, '(1x,2a)') 'Setup_cwb: Reading setup options from ',trim(filename)

 nerr = 0
 call read_inopt(datafile,'datafile',db,errcount=nerr)
 call read_inopt(m_gas,'m_gas',db,errcount=nerr)
 call read_inopt(h_sink,'h_sink',db,errcount=nerr)
 call read_inopt(ninjectmax,'ninjectmax',db,errcount=nerr)

 if (nerr > 0) then
    print "(1x,a,i2,a)",'Setup_cwb: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif
 call close_db(db)

end subroutine read_setupfile

!------------------------------------------
!+
!  Prompt user for setup options
!+
!------------------------------------------
subroutine interactive_setup()
 use prompting, only:prompt

 print "(2(/,a),/)",'*** Welcome to your friendly colliding-wind binary setup',&
                 '    ...where the stars are so bright that their winds sail on starlight ***'
 call prompt('Enter filename for star data',datafile,noblank=.true.)
 call prompt('Enter mass resolution of injected gas particles in Msun',m_gas,1.e-15,1.)
 call prompt('Enter sink particle radii in Rsun',h_sink,1.e-5,100.)
 call prompt('Enter max number of particles to inject in a single timestep per star (0=no limit)',ninjectmax,0,1000000000)
 print "(a)"

end subroutine interactive_setup

end module setup
