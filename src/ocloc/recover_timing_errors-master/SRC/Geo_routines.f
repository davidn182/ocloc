      subroutine bufdms (buff,lgh,hem,dd,dm,ds,ierror)
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
c
      logical     done,flag
c
      character*1 buff(*),abuf(21)
      character*1 ch
      character*1 hem
      integer*4   ll,lgh
      integer*4   i4,id,im,is,icond,ierror
      real*8      x(5)
c
c     set the "error flag" 
c
      ierror = 0
      icond  = 0
c
c     set defaults for dd,dm,ds
c
      dd = 0.0d0
      dm = 0.0d0
      ds = 0.0d0
c
c     set default limits for "hem" flag
c
      if(     hem.eq.'N' .or. hem.eq.'S' )then
        ddmax = 90.0d0
      elseif( hem.eq.'E' .or. hem.eq.'W' )then
        ddmax = 360.0d0
      elseif( hem.eq.'A' )then
        ddmax = 360.0d0
      elseif( hem.eq.'Z' )then
        ddmax = 180.0d0
      elseif( hem.eq.'*' )then
        ddmax  = 0.0d0
        ierror = 1
      else
        ddmax = 360.0d0
      endif
c
      do 1 i=1,5
        x(i) = 0.0d0
    1 continue
c
      icolon = 0
      ipoint = 0
      icount = 0
      flag   = .true.
      jlgh   = lgh
c
      do 2 i=1,jlgh
        if( buff(i).eq.':' )then
          icolon = icolon+1
        endif
        if( buff(i).eq.'.' )then
          ipoint = ipoint+1
          flag   = .false.
        endif
        if( flag )then
          icount = icount+1
        endif
    2 continue
c
      if( ipoint.eq.1 .and. icolon.eq.0 )then
c
c       load temp buffer
c
        do 3 i=1,jlgh
          abuf(i) = buff(i)
    3   continue
        abuf(jlgh+1) = '$'
        ll = jlgh
c
        call gvalr8 (abuf,ll,r8,icond)
c
        if( icount.ge.5 )then
c
c         value is a packed decimal of ==>  DDMMSS.sssss       
c
          ss = r8/10000.0d0
          id = idint( ss )
c
          r8 = r8-10000.0d0*dble(float(id))
          ss = r8/100.0d0
          im = idint( ss )
c
          r8 = r8-100.0d0*dble(float(im))
        else
c
c         value is a decimal of ==>  .xx   X.xxx   X.  
c
          id = idint( r8 )
          r8 = (r8-id)*60.0d0
          im = idint( r8 )
          r8 = (r8-im)*60.0d0
        endif
c
c       account for rounding error
c
        is = idnint( r8*1.0d5 )
        if( is.ge.6000000 )then
           r8 = 0.0d0
           im = im+1
        endif
c
        if( im.ge.60 )then
          im = 0
          id = id+1
        endif
c
        dd = dble( float( id ) )
        dm = dble( float( im ) )
        ds = r8
      else
c
c       buff() value is a d,m,s of ==>  NN:NN:XX.xxx    
c
        k    = 0
        next = 1
        done = .false.
        ie   = jlgh
c
        do 100 j=1,5
          ib = next
          do 90 i=ib,ie
            ch   = buff(i)
            last = i
            if( i.eq.jlgh .or. ch.eq.':' )then
              if( i.eq.jlgh )then
                done = .true.
              endif
              if( ch.eq.':' )then
                last = i-1
              endif
              goto 91
            endif
   90     continue
          goto 98
c
   91     ipoint = 0
          ik     = 0
          do 92 i=next,last
            ik = ik+1
            ch = buff(i)
            if( ch.eq.'.' )then
              ipoint = ipoint+1
            endif
            abuf(ik) = buff(i) 
   92     continue
          abuf(ik+1) = '$' 
c
          ll = ik
          if( ipoint.eq.0 )then
            call gvali4 (abuf,ll,i4,icond)
            r8 = dble(float( i4 )) 
          else
            call gvalr8 (abuf,ll,r8,icond)
          endif
c
          k    = k+1
          x(k) = r8
c
   98     if( done )then
            goto 101
          endif
c
          next = last
   99     next = next+1     
          if( buff(next).eq.':' )then
            goto 99
          endif
  100   continue
c
c       load dd,dm,ds
c
  101   if( k.ge.1 )then
          dd = x(1)
        endif
c
        if( k.ge.2 )then
          dm = x(2)
        endif
c
        if( k.ge.3 )then
          ds = x(3)
        endif
      endif
c
      if( dd.gt.ddmax  .or.
     1    dm.ge.60.0d0 .or.
     1    ds.ge.60.0d0 )then
        ierror = 1
        dd = 0.0d0
        dm = 0.0d0
        ds = 0.0d0
      endif
c
      if( icond.ne.0 )then
        ierror = 1
      endif
c
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine elipss (elips)
      implicit double precision(a-h,o-z)
      character*1  answer
      character*30 elips
      common/elipsoid/a,f
      write(*,*) '  Other Ellipsoids.'
      write(*,*) '  -----------------'
      write(*,*) '  '
      write(*,*) '  A) Airy 1858'
      write(*,*) '  B) Airy Modified'
      write(*,*) '  C) Australian National'
      write(*,*) '  D) Bessel 1841'
      write(*,*) '  E) Clarke 1880'
      write(*,*) '  F) Everest 1830'
      write(*,*) '  G) Everest Modified'
      write(*,*) '  H) Fisher 1960'
      write(*,*) '  I) Fisher 1968'
      write(*,*) '  J) Hough 1956'
      write(*,*) '  K) International (Hayford)'
      write(*,*) '  L) Krassovsky 1938'
      write(*,*) '  M) NWL-9D (WGS 66)'
      write(*,*) '  N) South American 1969'
      write(*,*) '  O) Soviet Geod. System 1985'
      write(*,*) '  P) WGS 72'
      write(*,*) '  Q-Z) User defined.'
      write(*,*) '  '
      write(*,*) '  Enter choice : '
      read(*,10) answer
   10 format(a1)
c
      if(answer.eq.'A'.or.answer.eq.'a') then
        a=6377563.396d0
        f=1.d0/299.3249646d0
        elips='Airy 1858'
      elseif(answer.eq.'B'.or.answer.eq.'b') then
        a=6377340.189d0
        f=1.d0/299.3249646d0
        elips='Airy Modified'
      elseif(answer.eq.'C'.or.answer.eq.'c') then
        a=6378160.d0
        f=1.d0/298.25d0
        elips='Australian National'
      elseif(answer.eq.'D'.or.answer.eq.'d') then
        a=6377397.155d0
        f=1.d0/299.1528128d0
        elips='Bessel 1841'
      elseif(answer.eq.'E'.or.answer.eq.'e') then
        a=6378249.145d0
        f=1.d0/293.465d0
        elips='Clarke 1880'
      elseif(answer.eq.'F'.or.answer.eq.'f') then
        a=6377276.345d0
        f=1.d0/300.8017d0
        elips='Everest 1830'
      elseif(answer.eq.'G'.or.answer.eq.'g') then
        a=6377304.063d0
        f=1.d0/300.8017d0
        elips='Everest Modified'
      elseif(answer.eq.'H'.or.answer.eq.'h') then
        a=6378166.d0
        f=1.d0/298.3d0
        elips='Fisher 1960'
      elseif(answer.eq.'I'.or.answer.eq.'i') then
        a=6378150.d0
        f=1.d0/298.3d0
        elips='Fisher 1968'
      elseif(answer.eq.'J'.or.answer.eq.'j') then
        a=6378270.d0
        f=1.d0/297.d0
        elips='Hough 1956'
      elseif(answer.eq.'K'.or.answer.eq.'k') then
        a=6378388.d0
        f=1.d0/297.d0
        elips='International (Hayford)'
      elseif(answer.eq.'L'.or.answer.eq.'l') then
        a=6378245.d0
        f=1.d0/298.3d0
        elips='Krassovsky 1938'
      elseif(answer.eq.'M'.or.answer.eq.'m') then
        a=6378145.d0
        f=1.d0/298.25d0
        elips='NWL-9D  (WGS 66)'
      elseif(answer.eq.'N'.or.answer.eq.'n') then
        a=6378160.d0
        f=1.d0/298.25d0
        elips='South American 1969'
      elseif(answer.eq.'O'.or.answer.eq.'o') then
        a=6378136.d0
        f=1.d0/298.257d0
        elips='Soviet Geod. System 1985'
      elseif(answer.eq.'P'.or.answer.eq.'p') then
        a=6378135.d0
        f=1.d0/298.26d0
        elips='WGS 72'
      else
        elips = 'User defined.'
c
        write(*,*) '  Enter Equatorial axis,   a : '
        read(*,*) a
        a  = dabs(a)
c
        write(*,*) '  Enter either Polar axis, b or '
        write(*,*) '  Reciprocal flattening,   1/f : '
        read(*,*) ss
        ss = dabs(ss)
c
        f = 0.0d0
        if( 200.0d0.le.ss .and. ss.le.310.0d0 )then
          f = 1.d0/ss  
        elseif( 6000000.0d0.lt.ss .and. ss.lt.a )then
          f = (a-ss)/a
        else
          elips = 'Error: default GRS80 used.'
          a     = 6378137.0d0
          f     = 1.0d0/298.25722210088d0
        endif
      endif
c
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine fixdms (ideg, min, sec, tol )
c
      implicit double precision (a-h, o-z)
      implicit integer (i-n)
c
c     test for seconds near 60.0-tol
c
      if( sec.ge.( 60.0d0-tol ) )then
        sec  = 0.0d0
        min  = min+1
      endif
c
c     test for minutes near 60
c
      if( min.ge.60 )then
        min  = 0
        ideg = ideg+1
      endif 
c
c     test for degrees near 360
c
      if( ideg.ge.360 )then
        ideg = 0
      endif 
c
      return
      end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hem_ns ( lat_sn, hem )
      implicit integer (i-n)
      character*6  hem
c
      if( lat_sn.eq.1 ) then
        hem = 'North '
      else
        hem = 'South '
      endif
c
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hem_ew ( lon_sn, hem )
      implicit integer (i-n)
      character*6  hem
c
      if( lon_sn.eq.1 ) then
        hem = 'East  '
      else
        hem = 'West  '
      endif
c
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine getrad(d,m,sec,isign,val)
 
*** comvert deg, min, sec to radians
 
      implicit double precision(a-h,j-z)
      common/const/pi,rad
 
      val=(d+m/60.d0+sec/3600.d0)/rad
      val=dble(isign)*val
 
      return
      end

      SUBROUTINE GPNARC (AMAX,FLAT,ESQ,PI,P1,P2,ARC)
C
C********1*********2*********3*********4*********5*********6*********7*
C
C NAME:        GPNARC
C VERSION:     200005.26
C WRITTEN BY:  ROBERT (Sid) SAFFORD
C PURPOSE:     SUBROUTINE TO COMPUTE THE LENGTH OF A MERIDIONAL ARC 
C              BETWEEN TWO LATITUDES
C
C INPUT PARAMETERS:
C -----------------
C AMAX         SEMI-MAJOR AXIS OF REFERENCE ELLIPSOID
C FLAT         FLATTENING (0.0033528 ... )
C ESQ          ECCENTRICITY SQUARED FOR REFERENCE ELLIPSOID
C PI           3.14159...
C P1           LAT STATION 1
C P2           LAT STATION 2
C
C OUTPUT PARAMETERS:
C ------------------
C ARC          GEODETIC DISTANCE 
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C GLOBAL VARIABLES AND CONSTANTS:
C -------------------------------
C
C    MODULE CALLED BY:    GENERAL 
C
C    THIS MODULE CALLS:   
C       LLIBFORE/ OPEN,   CLOSE,  READ,   WRITE,  INQUIRE
C                 DABS,   DBLE,   FLOAT,  IABS,   CHAR,   ICHAR
C
C    INCLUDE FILES USED:
C    COMMON BLOCKS USED:  
C
C    REFERENCES: Microsoft FORTRAN 4.10 Optimizing Compiler, 1988
C                MS-DOS Operating System
C    COMMENTS:
C********1*********2*********3*********4*********5*********6*********7*
C::MODIFICATION HISTORY
C::197507.05, RWS, VER 00 TENCOL RELEASED FOR FIELD USE
C::198311.20, RWS, VER 01 MTEN   RELEASED TO FIELD
C::198411.26, RWS, VER 07 MTEN2  RELEASED TO FIELD
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
CE::GPNARC
C ---------------------------
C     M T E N  (VERSION 3)
C     M T E N  (VERSION 5.23)
C ---------------------------
C 
      IMPLICIT REAL*8 (A-H,O-Z)
C
      LOGICAL  FLAG
C
      DATA TT/5.0D-15/
C
C     CHECK FOR A 90 DEGREE LOOKUP
C
      FLAG = .FALSE.
C
      S1 = DABS(P1)
      S2 = DABS(P2)
C
      IF( (PI/2.0D0-TT).LT.S2 .AND. S2.LT.(PI/2.0D0+TT) )THEN
        FLAG = .TRUE.
      ENDIF
C
      IF( S1.GT.TT )THEN
        FLAG = .FALSE.
      ENDIF
C
      DA = (P2-P1)
      S1 = 0.0D0
      S2 = 0.0D0
C
C     COMPUTE THE LENGTH OF A MERIDIONAL ARC BETWEEN TWO LATITUDES
C
      E2 = ESQ
      E4 = E2*E2
      E6 = E4*E2
      E8 = E6*E2
      EX = E8*E2
C
      T1 = E2*(003.0D0/4.0D0)
      T2 = E4*(015.0D0/64.0D0)
      T3 = E6*(035.0D0/512.0D0)
      T4 = E8*(315.0D0/16384.0D0)
      T5 = EX*(693.0D0/131072.0D0)
C
      A  = 1.0D0+T1+3.0D0*T2+10.0D0*T3+35.0D0*T4+126.0D0*T5
C
      IF( FLAG )THEN
        GOTO 1
      ENDIF
C
      B  = T1+4.0D0*T2+15.0D0*T3+56.0D0*T4+210.0D0*T5
      C  = T2+06.0D0*T3+28.0D0*T4+120.0D0*T5
      D  = T3+08.0D0*T4+045.0D0*T5
      E  = T4+010.0D0*T5
      F  = T5
C
      DB = DSIN(P2*2.0D0)-DSIN(P1*2.0D0)
      DC = DSIN(P2*4.0D0)-DSIN(P1*4.0D0)
      DD = DSIN(P2*6.0D0)-DSIN(P1*6.0D0)
      DE = DSIN(P2*8.0D0)-DSIN(P1*8.0D0)
      DF = DSIN(P2*10.0D0)-DSIN(P1*10.0D0)
C
C     COMPUTE THE S2 PART OF THE SERIES EXPANSION
C
      S2 = -DB*B/2.0D0+DC*C/4.0D0-DD*D/6.0D0+DE*E/8.0D0-DF*F/10.0D0
C
C     COMPUTE THE S1 PART OF THE SERIES EXPANSION
C
    1 S1 = DA*A
C
C     COMPUTE THE ARC LENGTH
C
      ARC = AMAX*(1.0D0-ESQ)*(S1+S2)
C
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine gvali4 (buff,ll,vali4,icond)
      implicit     integer (i-n)
c
      logical      plus,sign,done,error
      character*1  buff(*)
      character*1  ch
c
c     integer*2    i
c     integer*2    l1
c
      integer*4    ich,icond
      integer*4    ll    
      integer*4    vali4
c
      l1    = ll
      vali4 = 0
      icond = 0
      plus  = .true.
      sign  = .false.
      done  = .false.
      error = .false.
c
      i = 0
   10 i = i+1
      if( i.gt.l1 .or. done )then
        go to 1000
      else
        ch  = buff(i)
        ich = ichar( buff(i) )
      endif
c
      if(     ch.eq.'+' )then
c
c       enter on plus sign
c
        if( sign )then
          goto 150
        else 
          sign = .true.
          goto 10
        endif
      elseif( ch.eq.'-' )then
c
c       enter on minus sign
c
        if( sign )then
          goto 150
        else
          sign = .true.
          plus = .false.
          goto 10
        endif
      elseif( ch.ge.'0' .and. ch.le.'9' )then
        goto 100
      elseif( ch.eq.' ' )then
c
c       enter on space -- ignore leading spaces
c
        if( .not.sign )then
          goto 10
        else
          buff(i) = '0'
          ich = 48
          goto 100
        endif
      elseif( ch.eq.':' )then
c
c       enter on colon -- ignore 
c
        if( .not.sign )then
          goto 10
        else
          goto 1000
        endif
      elseif( ch.eq.'$' )then
c
c       enter on dollar "$"      
c
        done = .true.
        goto 10
      else
c
c       something wrong
c
        goto 150
      endif
c
c     enter on numeric
c
  100 vali4 = 10*vali4+(ich-48)
      sign  = .true.
      goto 10
c
c     treat illegal character
c
  150 buff(i) = '0'
      vali4 = 0
      icond = 1
c
 1000 if( .not.plus )then
        vali4 = -vali4
      endif
c
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine gvalr8 (buff,ll,valr8,icond)
      implicit     integer (i-n)
c
      logical      plus,sign,dpoint,done
c
      character*1  buff(*)
      character*1  ch
c
c     integer*2    i, ip
c     integer*2    l1
c     integer*2    nn, num, n48
c
      integer*4    ich,icond
      integer*4    ll
c
      real*8       ten
      real*8       valr8
      real*8       zero
c
      data zero,ten/0.0d0,10.0d0/
c
      n48     =  48
      l1      =  ll
      icond   =   0
      valr8   =  zero  
      plus    = .true.
      sign    = .false.
      dpoint  = .false.
      done    = .false.
c
c     start loop thru buffer
c
      i = 0
   10 i = i+1
      if( i.gt.l1 .or. done )then
        go to 1000
      else 
        ch  = buff(i)
        nn  = ichar( ch )
        ich = nn
      endif 
c
      if(     ch.eq.'+' )then
c
c       enter on plus sign
c
        if( sign )then
          goto 150
        else
          sign = .true.
          goto 10
        endif
      elseif( ch.eq.'-' )then
c
c       enter on minus sign
c
        if( sign )then
          goto 150
        else
          sign = .true.
          plus = .false.
          goto 10
        endif
      elseif( ch.eq.'.' )then
c
c       enter on decimal point
c
        ip     = 0
        sign   = .true.
        dpoint = .true.
        goto 10
      elseif( ch.ge.'0' .and. ch.le.'9' )then
        goto 100
      elseif( ch.eq.' ' )then
c
c       enter on space
c
        if( .not.sign )then
          goto 10
        else
          buff(i) = '0'
          ich = 48
          goto 100
        endif
      elseif( ch.eq.':' .or. ch.eq.'$' )then
c
c       enter on colon or "$" sign
c
        done = .true.
        goto 10
      else
c
c       something wrong
c
        goto 150
      endif
c
c     enter on numeric
c
  100 sign = .true.
      if( dpoint )then
        ip = ip+1
      endif
c
      num   = ich
      valr8 = ten*valr8+dble(float( num-n48 ))
      goto 10
c
c     treat illegal character
c
  150 buff(i) = '0'
      valr8   =  0.0d0
      icond   =  1
c
 1000 if( dpoint )then
        valr8 =  valr8/(ten**ip)
      endif
c
      if( .not.plus )then
        valr8 = -valr8
      endif
c
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine todmsp(val,id,im,s,isign)
 
*** convert position radians to deg,min,sec
*** range is [-pi to +pi]
 
      implicit double precision(a-h,o-z)
      common/const/pi,rad
 
    1 if(val.gt.pi) then
        val=val-pi-pi
        go to 1
      endif
 
    2 if(val.lt.-pi) then
        val=val+pi+pi
        go to 2
      endif
 
      if(val.lt.0.d0) then
        isign=-1
      else
        isign=+1
      endif
 
      s=dabs(val*rad)
      id=idint(s)
      s=(s-id)*60.d0
      im=idint(s)
      s=(s-im)*60.d0
 
*** account for rounding error
 
      is=idnint(s*1.d5)
      if(is.ge.6000000) then
        s=0.d0
        im=im+1
      endif
      if(im.ge.60) then
        im=0
        id=id+1
      endif
 
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine trim (buff,lgh,hem)
c
      implicit integer (i-n)
      character*1 ch,hem
      character*1 buff(*)
      integer*4   lgh
c
      ibeg = 1
      do 10 i=1,50
        if( buff(i).ne.' ' )then
          goto 11
        endif
        ibeg = ibeg+1
   10 continue
   11 continue
      if( ibeg.ge.50 )then
        ibeg = 1
        buff(ibeg) = '0'
      endif
c
      iend = 50
      do 20 i=1,50
        j = 51-i
        if( buff(j).eq.' ' )then
          iend = iend-1
        else
          goto 21
        endif
   20 continue
   21 continue
c
      ch = buff(ibeg)
      if( hem.eq.'N' )then
        if( ch.eq.'N' .or. ch.eq.'n' .or. ch.eq.'+' )then
          hem = 'N'
          ibeg = ibeg+1
        endif
        if( ch.eq.'S' .or. ch.eq.'s' .or. ch.eq.'-' )then
          hem = 'S'
          ibeg = ibeg+1
        endif
c
c       check for wrong hemisphere entry
c
        if( ch.eq.'E' .or. ch.eq.'e' )then
          hem = '*'
          ibeg = ibeg+1
        endif
        if( ch.eq.'W' .or. ch.eq.'w' )then
          hem = '*'
          ibeg = ibeg+1
        endif
      elseif( hem.eq.'W' )then
        if( ch.eq.'E' .or. ch.eq.'e' .or. ch.eq.'+' )then
          hem = 'E'
          ibeg = ibeg+1
        endif
        if( ch.eq.'W' .or. ch.eq.'w' .or. ch.eq.'-' )then
          hem = 'W'
          ibeg = ibeg+1
        endif
c
c       check for wrong hemisphere entry
c
        if( ch.eq.'N' .or. ch.eq.'n' )then
          hem = '*'
          ibeg = ibeg+1
        endif
        if( ch.eq.'S' .or. ch.eq.'s' )then
          hem = '*'
          ibeg = ibeg+1
        endif
      elseif( hem.eq.'A' )then
        if( .not.('0'.le.ch .and. ch.le.'9') )then
          hem = '*'
          ibeg = ibeg+1
        endif
      else
c        do nothing
      endif
c
c
      do 30 i=ibeg,iend
        ch = buff(i)
c
        if(     ch.eq.':' .or. ch.eq.'.' )then
          goto 30
        elseif( ch.eq.' ' .or. ch.eq.',' )then
          buff(i) = ':'
        elseif( '0'.le.ch .and. ch.le.'9' )then
          goto 30      
        else
          buff(i) = ':'
        endif
c
   30 continue
c
c     left justify buff() array to its first character position
c     also check for a ":" char in the starting position,
c     if found!!  skip it
c
      j  = 0
      ib = ibeg
      ie = iend
c
      do 40 i=ib,ie
        if( i.eq.ibeg .and. buff(i).eq.':' )then
c
c         move the 1st position pointer to the next char &
c         do not put ":" char in buff(j) array where j=1    
c
          ibeg = ibeg+1
          goto 40
        endif
        j = j+1
        buff(j) = buff(i)
   40 continue
c
c
      lgh = iend-ibeg+1
      j   = lgh+1
      buff(j) = '$'
c
c     clean-up the rest of the buff() array
c
      do 50 i=j+1,50   
        buff(i) = ' '    
   50 continue
c
c     save a maximum of 20 characters
c
      if( lgh.gt.20 )then
        lgh = 20
        j   = lgh+1
        buff(j) = '$'
      endif
c
      return
      end

