      subroutine Inverse(buffst1ns,buffst1ew,buffst2ns,buffst2ew,
     &edist,faz)
c
c********1*********2*********3*********4*********5*********6*********7**
c
c name:      inverse
c version:   200208.19
c author:    stephen j. frakes
c purpose:   to compute a geodetic inverse  
c            and then display output information
c
c input parameters:
c -----------------
c
c output parameters:
c ------------------
c
c local variables and constants:
c ------------------------------
c answer           user prompt response
c b                semiminor axis polar (in meters)
c baz              azimuth back (in radians)
c buff             input char buffer array
c dd,dm,ds         temporary values for degrees, minutes, seconds
c dlon             temporary value for difference in longitude (radians)   
c dmt,d_mt         char constants for meter units         
c edist            ellipsoid distance (in meters)
c elips            ellipsoid choice
c esq              eccentricity squared for reference ellipsoid
c faz              azimuth forward (in radians)
c filout           output file name
c finv             reciprocal flattening
c hem              hemisphere flag for lat & lon entry  
c ierror           error condition flag with d,m,s conversion
c lgh              length of buff() array
c option           user prompt response             
c r1,r2            temporary variables    
c ss               temporary variable     
c tol              tolerance for conversion of seconds
c
c name1            name of station one
c ld1,lm1,sl1      latitude  sta one - degrees,minutes,seconds
c ald1,alm1,sl1    latitude  sta one - degrees,minutes,seconds
c lat1sn           latitude  sta one - sign (+/- 1)
c d_ns1            latitude  sta one - char ('N','S')
c lod1,lom1,slo1   longitude sta one - degrees,minutes,seconds
c alod1,alom1,slo1 longitude sta one - degrees,minutes,seconds
c lon1sn           longitude sta one - sign (+/- 1)
c d_ew1            longitude sta one - char ('E','W')
c iaz1,maz1,saz1   forward azimuth   - degrees,minutes,seconds
c isign1           forward azimuth   - flag  (+/- 1)
c glat1,glon1      station one       - (lat & lon in radians )
c p1,e1            standpoint one    - (lat & lon in radians )
c
c name2            name of station two
c ld2,lm2,sl2      latitude  sta two - degrees,minutes,seconds
c ald2,alm2,sl2    latitude  sta two - degrees,minutes,seconds
c lat2sn           latitude  sta two - sign (+/- 1)
c d_ns2            latitude  sta two - char ('N','S')
c lod2,lom2,slo2   longitude sta two - degrees,minutes,seconds
c alod2,alom2,slo2 longitude sta one - degrees,minutes,seconds
c lon2sn           longitude sta two - sign (+/- 1)
c d_ew2            longitude sta two - char ('E','W')
c iaz2,maz2,saz2   back azimuth      - degrees,minutes,seconds
c isign2           back azimuth      - flag  (+/- 1)
c glat2,glon2      station two       - (lat & lon in radians )
c p2,e2            forepoint two     - (lat & lon in radians )
c
c global variables and constants:
c -------------------------------
c a                semimajor axis equatorial (in meters)
c f                flattening
c pi               constant 3.14159....
c rad              constant 180.0/pi  
c
c    this module called by:  n/a
c
c    this module calls:      elipss, getrad, inver1, todmsp
c    gethem, trim,   bufdms, gvalr8, gvali4, fixdms, gpnhri
c    datan,  write,  read,   dabs,   open,   stop
c
c    include files used:     n/a
c
c    common blocks used:     const, elipsoid
c
c    references:             see comments within subroutines
c
c    comments:
c
c********1*********2*********3*********4*********5*********6*********7**
c::modification history
c::1990mm.dd, sjf, creation of program           
c::199412.15, bmt, creation of program on viper
c::200203.08, crs, modified by c.schwarz to correct spelling of Clarke
c::200207.15, rws, modified i/o & standardized program documentation
c::                added subs trim, bufdms, gethem, gvali4, gvalr8      
c::200207.23, rws, replaced sub inver1 with gpnarc, gpnloa, gpnhri
c::200208.15, rws, fixed an error in bufdms
c::              - renamed ellips to elipss "common error" with dirct2
c::              - added FAZ & BAZ to printed output      
c::200208.19, rws, added more error flags for web interface code
c::              - added logical nowebb                             
c::200208.xx, sjf, program version number 2.0          
c********1*********2*********3*********4*********5*********6*********7**
ce::inverse
c
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
c
      logical  nowebb
c
      character*1  buffst1ns(50),buffst1ew(50)
      character*1  buffst2ns(50),buffst2ew(50)
      character*1  dmt,hem
      character*6  d_ns1, d_ew1, d_ns2, d_ew2, d_mt
      character*30 elips
c
      integer*4    ierror
      integer*4    lgh
c
      common/const/pi,rad
      common/elipsoid/a,f
c
c     ms_unix      0 = web version
c                  1 = ms_dos or unix version
c
      pi   = 4.d0*datan(1.d0)
      rad  = 180.d0/pi
      dmt  = 'm'
      d_mt = 'Meters'
      nowebb = .false.
c
      a=6378137.d0
      f=1.d0/298.25722210088d0
      elips='GRS80 / WGS84  (NAD83)'
c
      esq = f*(2.0d0-f)
c
      hem='N'
      call trim(buffst1ns,lgh,hem)
      call bufdms (buffst1ns,lgh,hem,dd,dm,ds,ierror)
c
      if( ierror.eq.0 )then
        irlat1 = 0
      else
        irlat1 = 1
        write(0,*) ' Invalid Latitude ... Program stopped '
        stop
      endif
c
      ald1 = dd
      alm1 = dm
      sl1  = ds
c
      if( hem.eq.'N' ) then
        lat1sn = +1
      else
        lat1sn = -1
      endif
c
      hem='E'
      call trim (buffst1ew,lgh,hem)
      call bufdms (buffst1ew,lgh,hem,dd,dm,ds,ierror)
c
      if( ierror.eq.0 )then
        irlon1 = 0
      else
        irlon1 = 1
        write(0,*) ' Invalid Longitude ... Program stopped '
        stop
      endif
c
      alod1 = dd
      alom1 = dm
      slo1  = ds
c
      if( hem.eq.'E' ) then
        lon1sn = +1
      else
        lon1sn = -1
      endif
c
      call getrad(ald1, alm1, sl1, lat1sn, glat1)
      call getrad(alod1,alom1,slo1,lon1sn, glon1)
c
      hem='N'
      call trim (buffst2ns,lgh,hem)
      call bufdms (buffst2ns,lgh,hem,dd,dm,ds,ierror)
c
      if( ierror.eq.0 )then
        irlat2 = 0
      else
        irlat2 = 1
        write(0,*) ' Invalid Latitude ... Program stopped '
        stop
      endif
c
      ald2 = dd
      alm2 = dm
      sl2  = ds
c
      if( hem.eq.'N' ) then
        lat2sn = +1
      else
        lat2sn = -1
      endif
c
c
      hem='E'
      call trim (buffst2ew,lgh,hem)
      call bufdms (buffst2ew,lgh,hem,dd,dm,ds,ierror)
c
      if( ierror.eq.0 )then
        irlon2 = 0
      else
        irlon2 = 1
        write(0,*) ' Invalid Longitude ... Program stopped '
        stop
      endif
c
      alod2 = dd
      alom2 = dm
      slo2  = ds
c
      if( hem.eq.'E' )then
        lon2sn = +1
      else
        lon2sn = -1
      endif
c
      call getrad(ald2, alm2, sl2, lat2sn, glat2)
      call getrad(alod2,alom2,slo2,lon2sn, glon2)
c
      p1 = glat1
      e1 = glon1
      p2 = glat2
      e2 = glon2
c
      if( e1.lt.0.0d0 )then
        e1 = e1+2.0d0*pi
      endif
      if( e2.lt.0.0d0 )then
        e2 = e2+2.0d0*pi
      endif
c
c     compute the geodetic inverse
c
c ************************************************************
c *   replaced subroutine inver1 with gpnhri
c *  
c *   call inver1 (glat1,glon1,glat2,glon2,faz,baz,edist)
c *
c ************************************************************
c
      call gpnhri (a,f,esq,pi,p1,e1,p2,e2,faz,baz,edist)
c
c     check for a non distance ... p1,e1 & p2,e2 equal zero ?
c
      if( edist.lt.0.00005d0 )then
        faz = 0.0d0
        baz = 0.0d0
      endif
c
c     set the tolerance (in seconds) for the azimuth conversion
c
      tol = 0.00005d0
c
      call todmsp(faz,iaz1,maz1,saz1,isign1)
      if(isign1.lt.0) then
        iaz1=359-iaz1
        maz1=59-maz1
        saz1=60.d0-saz1
      endif
      call fixdms ( iaz1, maz1, saz1, tol )
c
      call todmsp(baz,iaz2,maz2,saz2,isign2)
      if(isign2.lt.0) then
        iaz2=359-iaz2
        maz2=59-maz2
        saz2=60.d0-saz2
      endif
      call fixdms ( iaz2, maz2, saz2, tol ) 
c
      call todmsp(glat1, ld1,  lm1,  sl1,  lat1sn)
      call todmsp(glon1, lod1, lom1, slo1, lon1sn)
      call todmsp(glat2, ld2,  lm2,  sl2,  lat2sn)
      call todmsp(glon2, lod2, lom2, slo2, lon2sn)
c
      call hem_ns ( lat1sn, d_ns1 )
      call hem_ew ( lon1sn, d_ew1 )
      call hem_ns ( lat2sn, d_ns2 )
      call hem_ew ( lon2sn, d_ew2 )
c 
      finv=1.d0/f
      b=a*(1.d0-f)

      if(faz.lt.0)then
        faz=(faz+(2*pi))*rad
      elseif(faz.ge.(2*pi))then
        faz=(faz-(2*pi))*rad
      else
        faz=faz*rad
      endif 
      return 
      end subroutine Inverse



      subroutine gpnhri (a,f,esq,pi,p1,e1,p2,e2,az1,az2,s)      
c
c********1*********2*********3*********4*********5*********6*********7*
c
c name:        gpnhri
c version:     200208.09
c written by:  robert (sid) safford
c purpose:     subroutine to compute helmert rainsford inverse problem 
c 
c     solution of the geodetic inverse problem after t. vincenty
c     modified rainsford's method with helmert's elliptical terms
c     effective in any azimuth and at any distance short of antipocal
c     from/to stations must not be the geographic pole.
c     parameter a is the semi-major axis of the reference ellipsoid
c     finv=1/f is the inverse flattening of the reference ellipsoid
c     latitudes and longitudes in radians positive north and west
c     forward and back azimuths returned in radians clockwise from south
c     geodesic distance s returned in units of semi-major axis a
c     programmed for ibm 360-195   09/23/75
c
c     note - note - note -
c     1. do not use for meridional arcs and be careful on the equator.
c     2. azimuths are from north(+) clockwise and 
c     3. longitudes are positive east(+) 
c
c input parameters:
c -----------------
c a            semi-major axis of reference ellipsoid      meters
c f            flattening (0.0033528...)
c esq          eccentricity squared 
c pi           3.14159...
c p1           lat station 1                               radians
c e1           lon station 1                               radians
c p2           lat station 2                               radians
c e2           lon station 2                               radians
c
c output parameters:
c ------------------
c az1          azi at sta 1 -> sta 2                       radians
c az2          azi at sta 2 -> sta 1                       radians
c s            geodetic dist between sta(s) 1 & 2          meters
c
c local variables and constants:
c ------------------------------
c aa               constant from subroutine gpnloa                    
c alimit           equatorial arc distance along the equator   (radians)
c arc              meridional arc distance latitude p1 to p2 (in meters)      
c az1              azimuth forward                          (in radians)
c az2              azimuth back                             (in radians)
c bb               constant from subroutine gpnloa                    
c dlon             temporary value for difference in longitude (radians)   
c equ              equatorial distance                       (in meters)
c r1,r2            temporary variables    
c s                ellipsoid distance                        (in meters)
c sms              equatorial - geodesic distance (S - s) "Sms"       
c ss               temporary variable     
c tol0             tolerance for checking computation value         
c tol1             tolerance for checking a real zero value         
c tol2             tolerance for close to zero value  
c twopi            two times constant pi               
c
c global variables and constants:
c -------------------------------
c
c    module called by:    general 
c
c    this module calls:   gpnarc, gpnloa
c       llibfore/ dsin,   dcos,   dsqrt,  dabs,  datan2, write
c
c    include files used:
c    common blocks used:  
c
c    references: microsoft fortran 4.10 optimizing compiler, 1988
c                ms-dos operating system
c    comments:
c********1*********2*********3*********4*********5*********6*********7*
c::modification history
c::197507.05, rws, ver 00 tencol released for field use
c::198311.20, rws, ver 01 mten   released to field
c::198411.26, rws, ver 07 mten2  released to field
c::198506.10, rws, wrk    enhancements released to field
c::198507.22, rws, code   modified for mten3
c::198509.01, rws, ver 11 mten3  released to field
c::198708.10, rws, code   modified to use new mten4 gpn record format
c::199112.31, rws, ver 20 mten4 released to field
c::200001.13, rws, ver 21 mten4 released to field
c::200005.26, rws, code   restructured & documentation added             
c::200012.31, rws, ver 23 mten5 released                                 
c::200104.09, rws, code   added to calblin program                       
c::200208.09, rws, code   added subroutines gpnarc & gpnloa              
c********1*********2*********3*********4*********5*********6*********7*
ce::gpnhri
c  -------------------------------
c     m t e n  (version 3)
c              (version 4.22)
c              (version 5.23)
c  -------------------------------
c
      implicit real*8 (a-h,o-z)
c
      data tol0 /5.0d-15/
      data tol1 /5.0d-14/
      data tol2 /7.0d-03/
c
      twopi = 2.0d0*pi
c
c     test the longitude difference with tol1
c     tol1 is approximately 0.000000001 arc seconds
c
      ss = e2-e1
      if( dabs(ss).lt.tol1 )then
        e2 = e2+tol1
        write(*,*) ' longitudal difference is near zero '
c                 
        r2 = p2
        r1 = p1
        call gpnarc ( a, f, esq, pi, r1, r2, arc )
        s  = dabs( arc )
c
        if( p2.gt.p1 )then
          az1 = 0.0d0
          az2 = pi
        else
          az1 = pi   
          az2 = 0.0d0
        endif
        return 
      endif
c
c     test for longitude over 180 degrees
c
      dlon = e2-e1
c
      if( dlon.ge.0.0d0 )then
        if( pi.le.dlon .and. dlon.lt.twopi )then
          dlon = dlon-twopi
        endif
      else
        ss = dabs(dlon)
        if( pi.le.ss .and. ss.lt.twopi )then
          dlon = dlon+twopi
        endif
      endif
c
      ss = dabs( dlon )
      if( ss.gt.pi )then
c::     write(*,*) '  '
c::     write(*,*) ' Longitude difference over 180 degrees  '  
c::     write(*,*) ' Turn it around '
        ss = twopi-ss
      endif
c
c     compute the limit in longitude (alimit), it is equal 
c     to twice the distance from the equator to the pole,
c     as measured along the equator (east/ewst)
c
      alimit = pi*(1.0d0-f)
c
c     test for anti-nodal difference      
c
      if( ss.ge.alimit )then
        r1 = dabs(p1)
        r2 = dabs(p2)
c
c       latitudes r1 & r2 are not near the equator
c
        if( r1.gt.tol2 .and. r2.gt.tol2 )then
          goto 60
        endif
c
c       longitude difference is greater than lift-off point
c       now check to see if  "both"  r1 & r2 are on equator
c
        if( r1.lt.tol1 .and. r2.gt.tol2 )then
          goto 60
        endif
        if( r2.lt.tol1 .and. r1.gt.tol2 )then
          goto 60
        endif
c
c       check for either r1 or r2 just off the equator but < tol2
c
        if( r1.gt.tol1. or. r2.gt.tol1 )then
          az1 = 0.0d0
          az2 = 0.0d0
          s   = 0.0d0
          return 
        endif
c
c       compute the azimuth to anti-nodal point
c
c::     write(*,*) '  '
c::     write(*,*) ' Longitude difference beyond lift-off point '  
c::     write(*,*) '  '
c
        call gpnloa (a,f,esq,pi,dlon,az1,az2,aa,bb,sms)
c
c       compute the equatorial distance & geodetic
c
        equ = a*dabs(dlon)
        s   = equ-sms
        return 
      endif
c
   60 continue
c
      f0   = (1.0d0-f)
      b    = a*f0
      epsq = esq/(1.0d0-esq)
      f2   = f*f     
      f3   = f*f2    
      f4   = f*f3    
c
c     the longitude difference 
c
      dlon  = e2-e1   
      ab    = dlon      
      kount = 0    
c
c     the reduced latitudes    
c
      u1    = f0*dsin(p1)/dcos(p1)     
      u2    = f0*dsin(p2)/dcos(p2)
c
      u1    = datan(u1)
      u2    = datan(u2)
c
      su1   = dsin(u1)    
      cu1   = dcos(u1)    
c
      su2   = dsin(u2)
      cu2   = dcos(u2)
c
c     counter for the iteration operation
c
    1 kount = kount+1     
c
      clon  = dcos(ab)   
      slon  = dsin(ab)   
c
      csig  = su1*su2+cu1*cu2*clon  
      ssig  = dsqrt((slon*cu2)**2+(su2*cu1-su1*cu2*clon)**2)  
c
      sig   = datan2(ssig,csig)
      sinalf=cu1*cu2*slon/ssig
c
      w   = (1.0d0-sinalf*sinalf)
      t4  = w*w   
      t6  = w*t4   
c
c     the coefficients of type a      
c
      ao  = f-f2*(1.0d0+f+f2)*w/4.0d0+3.0d0*f3*(1.0d0+
     1        9.0d0*f/4.0d0)*t4/16.0d0-25.0d0*f4*t6/128.0d0
      a2  = f2*(1.0d0+f+f2)*w/4.0d0-f3*(1.0d0+9.0d0*f/4.0d0)*t4/4.0d0+
     1        75.0d0*f4*t6/256.0d0
      a4  = f3*(1.0d0+9.0d0*f/4.0d0)*t4/32.0d0-15.0d0*f4*t6/256.0d0
      a6  = 5.0d0*f4*t6/768.0d0
c
c     the multiple angle functions    
c
      qo  = 0.0d0
      if( w.gt.tol0 )then
        qo = -2.0d0*su1*su2/w
      endif     
c
      q2  = csig+qo
      q4  = 2.0d0*q2*q2-1.0d0    
      q6  = q2*(4.0d0*q2*q2-3.0d0)      
      r2  = 2.0d0*ssig*csig      
      r3  = ssig*(3.0d0-4.0d0*ssig*ssig) 
c
c     the longitude difference 
c
      s   = sinalf*(ao*sig+a2*ssig*q2+a4*r2*q4+a6*r3*q6)    
      xz  = dlon+s   
c
      xy  = dabs(xz-ab)    
      ab  = dlon+s   
c
      if( xy.lt.0.5d-13 )then
        goto 4
      endif
c
      if( kount.le.7 )then
        goto 1
      endif
c
c     the coefficients of type b      
c
    4 z   = epsq*w
c
      bo  = 1.0d0+z*(1.0d0/4.0d0+z*(-3.0d0/64.0d0+z*(5.0d0/256.0d0-
     1         z*175.0d0/16384.0d0)))      
      b2  = z*(-1.0d0/4.0d0+z*(1.0d0/16.0d0+z*(-15.0d0/512.0d0+
     1         z*35.0d0/2048.0d0)))  
      b4  = z*z*(-1.0d0/128.0d0+z*(3.0d0/512.0d0-z*35.0d0/8192.0d0))
      b6  = z*z*z*(-1.0d0/1536.0d0+z*5.0d0/6144.0d0)    
c
c     the distance in meters   
c
      s   = b*(bo*sig+b2*ssig*q2+b4*r2*q4+b6*r3*q6) 
c
c     first compute the az1 & az2 for along the equator
c
      if( dlon.gt.pi )then
        dlon = (dlon-2.0d0*pi)
      endif
c
      if( dabs(dlon).gt.pi )then
        dlon = (dlon+2.0d0*pi)
      endif
c
      az1 = pi/2.0d0
      if( dlon.lt.0.0d0 )then
        az1 = 3.0d0*az1
      endif
c
      az2 = az1+pi
      if( az2.gt.2.0d0*pi )then
        az2 = az2-2.0d0*pi
      endif
c
c     now compute the az1 & az2 for latitudes not on the equator
c
      if( .not.(dabs(su1).lt.tol0 .and. dabs(su2).lt.tol0) )then
        tana1 =  slon*cu2/(su2*cu1-clon*su1*cu2)  
        tana2 =  slon*cu1/(su1*cu2-clon*su2*cu1)  
        sina1 =  sinalf/cu1
        sina2 = -sinalf/cu2      
c
c       azimuths from north,longitudes positive east  
c
        az1   = datan2(sina1,sina1/tana1)   
        az2   = pi-datan2(sina2,sina2/tana2)
      endif
c
      if( az1.lt.0.0d0 )then
        az1 = az1+2.0d0*pi   
      endif
c
      if( az2.lt.0.0d0 )then
        az2 = az2+2.0d0*pi
      endif
c
      return     
      end 

CB::GPNLOA
C
      SUBROUTINE GPNLOA (AMAX,FLAT,ESQ,PI,DL,AZ1,AZ2,AO,BO,SMS)
C
C********1*********2*********3*********4*********5*********6*********7*
C
C NAME:        GPNLOA
C VERSION:     200005.26
C WRITTEN BY:  ROBERT (Sid) SAFFORD
C PURPOSE:     SUBROUTINE TO COMPUTE THE LIFF-OFF-AZIMUTH CONSTANTS
C
C INPUT PARAMETERS:
C -----------------
C AMAX         SEMI-MAJOR AXIS OF REFERENCE ELLIPSOID
C FLAT         FLATTENING (0.0033528 ... )
C ESQ          ECCENTRICITY SQUARED FOR REFERENCE ELLIPSOID
C PI           3.14159...
C DL           LON DIFFERENCE
C AZ1          AZI AT STA 1 -> STA 2
C
C OUTPUT PARAMETERS:
C ------------------
C AZ2          AZ2 AT STA 2 -> STA 1
C AO           CONST
C BO           CONST
C SMS          DISTANCE ... EQUATORIAL - GEODESIC  (S - s)   "SMS"
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C GLOBAL VARIABLES AND CONSTANTS:
C -------------------------------
C
C    MODULE CALLED BY:    GENERAL 
C
C    THIS MODULE CALLS:   
C       LLIBFORE/ DSIN,   DCOS,   DABS,   DASIN 
C
C    INCLUDE FILES USED:
C    COMMON BLOCKS USED:  
C
C    REFERENCES: Microsoft FORTRAN 4.10 Optimizing Compiler, 1988
C                MS-DOS Operating System
C    COMMENTS:
C********1*********2*********3*********4*********5*********6*********7*
C::MODIFICATION HISTORY
C::1985xx.xx, RWS, CODE   CREATED               
C::198506.10, RWS, WRK    ENHANCEMENTS RELEASED TO FIELD
C::198509.01, RWS, VER 11 MTEN3  RELEASED TO FIELD
C::198512.18, RWS, CODE   MODIFIED FOR MTEN3
C::198708.10, RWS, CODE   MODIFIED TO USE NEW MTEN4 GPN RECORD FORMAT
C::199112.31, RWS, VER 20 MTEN4 RELEASED TO FIELD
C::200001.13, RWS, VER 21 MTEN4 RELEASED TO FIELD
C::200005.26, RWS, CODE   RESTRUCTURED & DOCUMENTATION ADDED             
C::200012.31, RWS, VER 23 MTEN5 RELEASED                                 
C********1*********2*********3*********4*********5*********6*********7*
CE::GPNLOA
C ---------------------------
C     M T E N  (VERSION 3)
C              (VERSION 4.22)
C              (VERSION 5.23)
C ---------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DATA TT/5.0D-13/
C
      DLON = DABS(DL)
      CONS = (PI-DLON)/(PI*FLAT)
      F    = FLAT
C
C     COMPUTE AN APPROXIMATE AZ
C
      AZ   = DASIN(CONS)
C
      T1   =    1.0D0
      T2   =  (-1.0D0/4.0D0)*F*(1.0D0+F+F*F)
      T4   =    3.0D0/16.0D0*F*F*(1.0D0+(9.0D0/4.0D0)*F)
      T6   = (-25.0D0/128.0D0)*F*F*F
C
      ITER = 0
    1 ITER = ITER+1
      S    = DCOS(AZ)
      C2   = S*S
C
C     COMPUTE NEW AO
C
      AO   = T1 + T2*C2 + T4*C2*C2 + T6*C2*C2*C2
      CS   = CONS/AO
      S    = DASIN(CS)
      IF( DABS(S-AZ).LT.TT )THEN
        GOTO 2
      ENDIF
C
      AZ   = S
      IF( ITER.LE.6 )THEN
        GOTO 1
      ENDIF
C
    2 AZ1  = S
      IF( DL.LT.0.0D0 )THEN
        AZ1 = 2.0D0*PI-AZ1
      ENDIF
C
      AZ2  = 2.0D0*PI-AZ1
C
C     EQUATORIAL - GEODESIC  (S - s)   "SMS"
C
      ESQP = ESQ/(1.0D0-ESQ)
      S    = DCOS(AZ1)
C
      U2   = ESQP*S*S
      U4   = U2*U2
      U6   = U4*U2
      U8   = U6*U2
C
      T1   =     1.0D0
      T2   =    (1.0D0/4.0D0)*U2
      T4   =   (-3.0D0/64.0D0)*U4
      T6   =    (5.0D0/256.0D0)*U6
      T8   = (-175.0D0/16384.0D0)*U8
C
      BO   = T1 + T2 + T4 + T6 + T8
      S    = DSIN(AZ1)
      SMS  = AMAX*PI*(1.0D0 - FLAT*DABS(S)*AO - BO*(1.0D0-FLAT))
C
      RETURN
      END
