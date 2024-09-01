      SUBROUTINE upc(string)
      IMPLICIT NONE
      CHARACTER(*) :: string
      CHARACTER          :: ch
      INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
      INTEGER            :: i
      DO i = 1,LEN(string)
        ch = string(i:i)
        IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
        string(i:i) = ch
      END DO
      END

      SUBROUTINE openwr (lu,filename,st,ch1,ch2,i1,i2,ier)
      IMPLICIT NONE
      INTEGER*4 lu,i1,i2,ier
      CHARACTER(*) filename,st,ch1,ch2
      OPEN(lu,file=filename,status=st,iostat=ier)
      RETURN
      END

      FUNCTION lenact(string)
      implicit none
      INTEGER*4 lenact
      CHARACTER(*) string
      lenact=len_trim(string)
      RETURN
      END

      SUBROUTINE rmvlbk(string)
      implicit none 
      CHARACTER(*) string
      string = adjustl(string)
      RETURN
      END

      SUBROUTINE rdforn(string,length)
      implicit none 
      INTEGER*4 length,in
      CHARACTER(*) string
      CHARACTER*800 str2
      CALL GET_COMMAND(str2,length)
      in = index(str2,' ')
      length = len_trim(str2)
      string=str2(in+1:length)
      length = len_trim(string)
      RETURN
      END
c
      SUBROUTINE getlun(newunit)
c This is a simple function to search for an available unit.
c LUN_MIN and LUN_MAX define the range of possible LUNs to check.
c The UNIT value is returned by the function, and also by the optional
c argument. This allows the function to be used directly in an OPEN
c statement, and optionally save the result in a local variable.
c If no units are available, -1 is returned.
      IMPLICIT none
      integer :: newunit
      integer, parameter :: LUN_MIN=10, LUN_MAX=1000
      logical :: opened
      integer :: lun
! begin
      newunit=-1
      DO lun=LUN_MIN,LUN_MAX
        inquire(unit=lun,opened=opened)
        IF (.not. opened) THEN
          newunit=lun
          RETURN
        ENDIF 
      ENDDO 
      RETURN
      END

      SUBROUTINE frelun(unit)
      IMPLICIT NONE
      INTEGER*4 unit
      RETURN
      END

      SUBROUTINE sla_CALDJ (IY, IM, ID, DJM, J)
*+
*     - - - - - -
*      C A L D J
*     - - - - - -
*
*  Gregorian Calendar to Modified Julian Date
*
*  (Includes century default feature:  use sla_CLDJ for years
*   before 100AD.)
*
*  Given:
*     IY,IM,ID     int    year, month, day in Gregorian calendar
*
*  Returned:
*     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs
*     J            int    status:
*                           0 = OK
*                           1 = bad year   (MJD not computed)
*                           2 = bad month  (MJD not computed)
*                           3 = bad day    (MJD computed)
*
*  Acceptable years are 00-49, interpreted as 2000-2049,
*                       50-99,     "       "  1950-1999,
*                       100 upwards, interpreted literally.
*
*  Called:  sla_CLDJ
*
*  P.T.Wallace   Starlink   November 1985
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*
*  License:
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program (see SLA_CONDITIONS); if not, write to the 
*    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
*    Boston, MA  02111-1307  USA
*
*-

      IMPLICIT NONE

      INTEGER IY,IM,ID
      DOUBLE PRECISION DJM
      INTEGER J

      INTEGER NY




*  Default century if appropriate
      IF (IY.GE.0.AND.IY.LE.49) THEN
         NY=IY+2000
      ELSE IF (IY.GE.50.AND.IY.LE.99) THEN
         NY=IY+1900
      ELSE
         NY=IY
      END IF

*  Modified Julian Date
      CALL sla_CLDJ(NY,IM,ID,DJM,J)

      END


      SUBROUTINE sla_CLDJ (IY, IM, ID, DJM, J)
*+
*     - - - - -
*      C L D J
*     - - - - -
*
*  Gregorian Calendar to Modified Julian Date
*
*  Given:
*     IY,IM,ID     int    year, month, day in Gregorian calendar
*
*  Returned:
*     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs
*     J            int    status:
*                           0 = OK
*                           1 = bad year   (MJD not computed)
*                           2 = bad month  (MJD not computed)
*                           3 = bad day    (MJD computed)
*
*  The year must be -4699 (i.e. 4700BC) or later.
*
*  The algorithm is derived from that of Hatcher 1984
*  (QJRAS 25, 53-55).
*
*  P.T.Wallace   Starlink   11 March 1998
*
*  Copyright (C) 1998 Rutherford Appleton Laboratory
*
*  License:
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program (see SLA_CONDITIONS); if not, write to the 
*    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
*    Boston, MA  02111-1307  USA
*
*-

      IMPLICIT NONE

      INTEGER IY,IM,ID
      DOUBLE PRECISION DJM
      INTEGER J

*  Month lengths in days
      INTEGER MTAB(12)
      DATA MTAB / 31,28,31,30,31,30,31,31,30,31,30,31 /



*  Preset status
      J=0

*  Validate year
      IF (IY.LT.-4699) THEN
         J=1
      ELSE

*     Validate month
         IF (IM.GE.1.AND.IM.LE.12) THEN

*        Allow for leap year
            IF (MOD(IY,4).EQ.0) THEN
               MTAB(2)=29
            ELSE
               MTAB(2)=28
            END IF
            IF (MOD(IY,100).EQ.0.AND.MOD(IY,400).NE.0)
     :         MTAB(2)=28

*        Validate day
            IF (ID.LT.1.OR.ID.GT.MTAB(IM)) J=3

*        Modified Julian Date
            DJM=DBLE((1461*(IY-(12-IM)/10+4712))/4
     :               +(306*MOD(IM+9,12)+5)/10
     :               -(3*((IY-(12-IM)/10+4900)/100))/4
     :               +ID-2399904)

*        Bad month
         ELSE
            J=2
         END IF

      END IF

      END
c
c
c
      SUBROUTINE sort2(n,arr,brr)  
      INTEGER n,M,NSTACK  
      REAL*8 arr(n),brr(n)  
      PARAMETER (M=7,NSTACK=50)  
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)  
      REAL a,b,temp  
      jstack=0  
      l=1  
      ir=n  
1     if(ir-l.lt.M)then  
        do 12 j=l+1,ir  
          a=arr(j)  
          b=brr(j)  
          do 11 i=j-1,l,-1  
            if(arr(i).le.a)goto 2  
            arr(i+1)=arr(i)  
            brr(i+1)=brr(i)  
11        continue  
          i=l-1  
2         arr(i+1)=a  
          brr(i+1)=b  
12      continue  
        if(jstack.eq.0)return  
        ir=istack(jstack)  
        l=istack(jstack-1)  
        jstack=jstack-2  
      else  
        k=(l+ir)/2  
        temp=arr(k)  
        arr(k)=arr(l+1)  
        arr(l+1)=temp  
        temp=brr(k)  
        brr(k)=brr(l+1)  
        brr(l+1)=temp  
        if(arr(l).gt.arr(ir))then  
          temp=arr(l)  
          arr(l)=arr(ir)  
          arr(ir)=temp  
          temp=brr(l)  
          brr(l)=brr(ir)  
          brr(ir)=temp  
        endif  
        if(arr(l+1).gt.arr(ir))then  
          temp=arr(l+1)  
          arr(l+1)=arr(ir)  
          arr(ir)=temp  
          temp=brr(l+1)  
          brr(l+1)=brr(ir)  
          brr(ir)=temp  
        endif  
        if(arr(l).gt.arr(l+1))then  
          temp=arr(l)  
          arr(l)=arr(l+1)  
          arr(l+1)=temp  
          temp=brr(l)  
          brr(l)=brr(l+1)  
          brr(l+1)=temp  
        endif  
        i=l+1  
        j=ir  
        a=arr(l+1)  
        b=brr(l+1)  
3       continue  
          i=i+1  
        if(arr(i).lt.a)goto 3  
4       continue  
          j=j-1  
        if(arr(j).gt.a)goto 4  
        if(j.lt.i)goto 5  
        temp=arr(i)  
        arr(i)=arr(j)  
        arr(j)=temp  
        temp=brr(i)  
        brr(i)=brr(j)  
        brr(j)=temp  
        goto 3  
5       arr(l+1)=arr(j)  
        arr(j)=a  
        brr(l+1)=brr(j)  
        brr(j)=b  
        jstack=jstack+2  
        if(jstack.gt.NSTACK) print *, 'NSTACK too small in sort2'  
        if(ir-i+1.ge.j-l)then  
          istack(jstack)=ir  
          istack(jstack-1)=i  
          ir=j-1  
        else  
          istack(jstack)=j-1  
          istack(jstack-1)=l  
          l=i  
        endif  
      endif  
      goto 1  
      END  
c
      SUBROUTINE indexx(n,arr,indx)  
      INTEGER*4 n,indx(n),M,NSTACK  
      REAL*8 arr(n)  
      PARAMETER (M=7,NSTACK=50)  
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)  
      REAL a  
      do 11 j=1,n  
        indx(j)=j  
11    continue  
      jstack=0  
      l=1  
      ir=n  
1     if(ir-l.lt.M)then  
        do 13 j=l+1,ir  
          indxt=indx(j)  
          a=arr(indxt)  
          do 12 i=j-1,l,-1  
            if(arr(indx(i)).le.a)goto 2  
            indx(i+1)=indx(i)  
12        continue  
          i=l-1  
2         indx(i+1)=indxt  
13      continue  
        if(jstack.eq.0)return  
        ir=istack(jstack)  
        l=istack(jstack-1)  
        jstack=jstack-2  
      else  
        k=(l+ir)/2  
        itemp=indx(k)  
        indx(k)=indx(l+1)  
        indx(l+1)=itemp  
        if(arr(indx(l)).gt.arr(indx(ir)))then  
          itemp=indx(l)  
          indx(l)=indx(ir)  
          indx(ir)=itemp  
        endif  
        if(arr(indx(l+1)).gt.arr(indx(ir)))then  
          itemp=indx(l+1)  
          indx(l+1)=indx(ir)  
          indx(ir)=itemp  
        endif  
        if(arr(indx(l)).gt.arr(indx(l+1)))then  
          itemp=indx(l)  
          indx(l)=indx(l+1)  
          indx(l+1)=itemp  
        endif  
        i=l+1  
        j=ir  
        indxt=indx(l+1)  
        a=arr(indxt)  
3       continue  
          i=i+1  
        if(arr(indx(i)).lt.a)goto 3  
4       continue  
          j=j-1  
        if(arr(indx(j)).gt.a)goto 4  
        if(j.lt.i)goto 5  
        itemp=indx(i)  
        indx(i)=indx(j)  
        indx(j)=itemp  
        goto 3  
5       indx(l+1)=indx(j)  
        indx(j)=indxt  
        jstack=jstack+2  
        if(jstack.gt.NSTACK) print *, 'NSTACK too small in indexx'  
        if(ir-i+1.ge.j-l)then  
          istack(jstack)=ir  
          istack(jstack-1)=i  
          ir=j-1  
        else  
          istack(jstack)=j-1  
          istack(jstack-1)=l  
          l=i  
        endif  
      endif  
      goto 1  
      END  

      SUBROUTINE SCHRA(Ra,string)
C
C CHANGE ra in degs to a string containing hr, min and sec
C
c
      character*(*) string
      character(80) xalpha
      REAL*8 Ra
      INTEGER*4 status
c
      string = xalpha (Ra,status)
      if(status.ne.0)then
      string='** ** ******'
      endif
      RETURN
      END

      SUBROUTINE SCHDEC(Dec,String)
c
c change decimal degrees to a string containing deg mn sec
c
c
      character*(*) String
      character(80) xdelta
      integer*4 status
      REAL*8 Dec
      string = xdelta(Dec,status)
      status = 0
      if(status.ne.0)then
      string='*** ** *****'
      endif
c
      RETURN
      END

*- xalpha - alpha string from ra in degrees
      character(80) function xalpha(ra,status)

      include 'status.codes'

* Import :
      real*8 ra
* Status :
      integer status
* Local variables :
      integer hour,minute
      real*8 second
      character(80) a
*-
      if(status.ne.ok__)return

      hour=int(ra/15d0)
      minute=int(4d0*(ra-15d0*dble(hour)))
      second=6d1*4d0*(ra-15d0*(dble(hour)+dble(minute)/6d1))

      a=' '
      write(a(1:2),'(i2.2)',iostat=status)hour
      write(a(4:5),'(i2.2)',iostat=status)minute
      write(a(7:12),'(f6.3)',iostat=status)second
      if(a(7:7).eq.' ')a(7:7)='0'

      xalpha=a

      return

      end

*- xdelta - delta string from declination in degrees
      character(80) function xdelta(dec,status)

      include 'status.codes'

* Import :
      real*8 dec
* Status :
      integer status
* Local variables :
      real*8 d
      integer degree,arcmin
      real*8 arcsec
      character(80) a
*-
      if(status.ne.ok__)return

      if(dec.le.0d0)then
      a='-'
      d=-dec
      else
      a=' '
      d=dec
      endif
      degree=int(d)
      arcmin=int(6d1*(d-dble(degree)))
      arcsec=6d1*6d1*(d-dble(degree)-dble(arcmin)/6d1)

      write(a(2:3),'(i2.2)',iostat=status)degree
      write(a(5:6),'(i2.2)',iostat=status)arcmin
      write(a(8:12),'(f5.2)',iostat=status)arcsec
      if(a(8:8).eq.' ')a(8:8)='0'

      xdelta=a

      return

      end

      SUBROUTINE CHRA(Ra,Irh,Irm,Rsec,Iflag)
C
C CHANGE ra hr, min and sec to decimal degrees (iflag .ne. 1)
C          or decimal degrees to hr, mn, sec (iflag = 1 )
C
c
      REAL*8 Ra , rra , ram , rasec
      REAL*4 Rsec
      INTEGER*4 Irh , Irm , Iflag
c
      IF ( Iflag.EQ.1 ) THEN
      rra = Ra/15.
      Irh = rra
      ram = rra - Irh
      ram = ram*60.
      Irm = ram
      rasec = ram - Irm
      Rsec = rasec*60.
      IF ( Rsec.EQ.60.0 ) THEN
      Rsec = 0.0
      Irm = Irm + 1
      ENDIF
      IF ( Irm.EQ.60.0 ) THEN
      Irm = 0.0
      Irh = Irh + 1
      ENDIF
      IF ( Irh.GE.24.0 ) Irh = Irh - 24.0
      ELSE
      Ra = (DFLOAT(Irh)+(DFLOAT(Irm)+DBLE(Rsec)/60.)/60.)*15.D0
      ENDIF
      RETURN
      END

      SUBROUTINE CHDEC(Dec,Idd,Idm,Dsec,Iflag)
c
c change decimal degrees to deg mn sec (iflag = 1)
c or change decimal degrees to deg mn sec (iflag .ne. 1)
c
c
      character(1) sign
      REAL*8 Dec , rdm
      REAL*4 Dsec
      INTEGER*4 Iflag , Idd , iddd , Idm
c
      IF ( Iflag.EQ.1 ) THEN
      if(dec.lt.0.0)then
      sign = '-'
      else
      sign = '+'
      endif
c
      Idd = Dec
      iddd = ABS(Idd)
      rdm = ABS(Dec) - iddd
      rdm = rdm*60
      Idm = rdm
      Dsec = rdm - Idm
      Dsec = Dsec*60.
      IF ( Dsec.EQ.60.0 ) THEN
      Dsec = 0.0
      Idm = Idm + 1
      ENDIF
      IF ( Idm.EQ.60.0 ) THEN
      Idm = 0.0
      IF ( Idd.GT.0 ) THEN
      Idd = Idd + 1
      ELSE
      Idd = Idd - 1
      ENDIF
      ENDIF
c
      if(sign.eq.'-')then
      if(idd.gt.0)then
      idd = -idd
      elseif(idm.gt.0.and.idd.eq.0)then
      idm = -idm
      elseif(idm.eq.0.and.idd.eq.0)then
      dsec = -dsec
      endif
      endif
c
      ELSEIF ( Dsec.LT.0. ) THEN
      Dec = DBLE(Dsec)/3600.D0
      ELSEIF ( Idm.LT.0 ) THEN
      Dec = DFLOAT(Idm)/60.D0 - DBLE(Dsec)/3600.D0
      ELSEIF ( Idd.LT.0 ) THEN
      Dec = DFLOAT(Idd) - DFLOAT(Idm)/60.D0 - DBLE(Dsec)/3600.D0
      ELSE
      Dec = DFLOAT(Idd) + DFLOAT(Idm)/60.D0 + DBLE(Dsec)/3600.D0
      ENDIF
      RETURN
      END

      subroutine nhdeabsorb2(flag_nh,emin,emax,alpha,anh,reduce,iflag)
c
c    Calculates count-rate to flux conversion and nh absorption correction factors
c    for various satellites, instruments and filters
c
      IMPLICIT none
      INTEGER*4  max_sat , length
      PARAMETER (max_sat=200)
      CHARACTER*200 file_eff_area , string
      CHARACTER*200 f_eff_area(max_sat)
      CHARACTER*40 satellite(max_sat), instrument(max_sat)
      CHARACTER*1 yesno
      INTEGER*4 lu_input , model , ifear , ifl , itype
      INTEGER*4 i , j
      INTEGER*4 lu_pf, in, lenact
      INTEGER*4 iflag, ii, flag_nh
      REAL*4 nh , emin_area , emax_area , emin , emax
      REAL*4 emin_areas, emax_areas, emin_aream, emax_aream,nufnu
      REAL*4 emin_areah, emax_areah, fekev, conv
      REAL*4 energy
      REAL*4 eff_area , anh , alpha , alpha1 , bbreak_energy
      REAL*4 tt , ak , aa7 , alpha2 , gamm , bk , t , gamma1 , gamma2
      REAL*4 break_energy, umalpha
      REAL*4 t1 , ekev,reduce
      REAL*4 eemin(max_sat),eemax(max_sat), aa8
      REAL*4 eemin_area(max_sat), eemax_area(max_sat)
      REAL*4 eemin_areas(max_sat), eemax_areas(max_sat)
      REAL*4 eemin_aream(max_sat), eemax_aream(max_sat)
      REAL*4 eemin_areah(max_sat), eemax_areah(max_sat)
      CHARACTER*200 webprograms,filename, stringin
      LOGICAL ok
      EXTERNAL espec
      EXTERNAL eacma
      COMMON webprograms
      COMMON /filt  / energy(5000) , eff_area(5000)
      COMMON /esp   / nh , gamm , bk , ifl , t , itype , gamma1 ,
     &                gamma2 , break_energy
      ok = .TRUE.
         ekev = 1.
         model = 1
      filename = webprograms(1:lenact(webprograms)) //
     &'/count_rate/countrates2.cf'
      OPEN (15,FILE=filename,STATUS='old')
      i=0
c read configuration file
      DO WHILE(ok)
         string=' '
         READ (15,'(a)',end=100) string
         CALL rmvlbk(string)
         in =index(string,' ')
         CALL upc(string(1:in-1))
         IF (string(1:in-1).EQ.'SATELLITE') THEN
            i=i+1
            satellite(i)=string(in+1:lenact(string))
            call rmvlbk(satellite(i))
         ENDIF
         IF (string(1:in-1).EQ.'INSTRUMENT') THEN
            instrument(i)=string(in+1:lenact(string))
            call rmvlbk(instrument(i))
         ENDIF
         IF (string(1:in-1).EQ.'FILE_EFF_AREA') THEN
            f_eff_area(i)=webprograms(1:lenact(webprograms)) //
     &                    string(in+1:lenact(string))
         ENDIF
         IF (string(1:in-1).EQ.'FLUX_ENERGY_RANGE') THEN
            read(string(in+1:lenact(string)),*) eemin(i), eemax(i)
         ENDIF
         IF (string(1:in-1).EQ.'COUNTS_ENERGY_RANGE') THEN
            read(string(in+1:lenact(string)),*)
     &                       eemin_area(i), eemax_area(i)
         ENDIF
         IF (string(1:in-1).EQ.'SOFT_ENERGY_BAND') THEN
            read(string(in+1:lenact(string)),*)
     &                       eemin_areas(i), eemax_areas(i)
         ENDIF
         IF (string(1:in-1).EQ.'MEDIUM_ENERGY_BAND') THEN
            read(string(in+1:lenact(string)),*)
     &                       eemin_aream(i), eemax_aream(i)
         ENDIF
         IF (string(1:in-1).EQ.'HARD_ENERGY_BAND') THEN
            read(string(in+1:lenact(string)),*)
     &                       eemin_areah(i), eemax_areah(i)
         ENDIF
      ENDDO
100   CONTINUE
      IF (flag_nh .NE. 1) THEN
         file_eff_area=f_eff_area(iflag)
         emin_area=emin!eemin_area(iflag)
         emax_area=emax!eemax_area(iflag)
         emin_areas=eemin_areas(iflag)
         emax_areas=eemax_areas(iflag)
         emin_aream=eemin_aream(iflag)
         emax_aream=eemax_aream(iflag)
         emin_areah=eemin_areah(iflag)
         emax_areah=eemax_areah(iflag)
         OPEN (16,FILE=file_eff_area,STATUS='old')
         ok = .TRUE.
         i = 1
         DO WHILE (ok)
            READ (16,*,END=200) energy(i),eff_area(i)
            i = i + 1
         ENDDO
 200     CONTINUE
         DO j = i , 5000
           energy(j) = energy(i-1)
           eff_area(j) = eff_area(i-1)
         ENDDO
         CLOSE (16)
         CALL frelun(16)
      ENDIF
      IF ( alpha.EQ.1 ) alpha = alpha + 0.0001
      umalpha=1.-alpha
      IF (umalpha.eq.0.) THEN
         conv=ekev**(-alpha)*1./log(emax/emin)/(ekev*2.418E-12)
      ELSE
         conv=umalpha*ekev**(-alpha)/(emax**umalpha-emin**umalpha)/
     &                                             (ekev*2.418E-12)
      ENDIF
      IF (flag_nh == 1) THEN
c
c NH correction factor here
c
         ifear = 1
c ifear = 1 => flux NOT corrected for Galactic absorption
         CALL energycr(-1.,alpha,t1,anh,0.,model,5,ak,emin,emax,
     &      emin_area,emax_area,aa7,ifear,alpha1,
     &      alpha2,bbreak_energy)
            ifear = 0
c ifear = 0 => flux corrected for Galactic absorption
         CALL energycr(-1.,alpha,t1,anh,0.,model,5,aa8,emin,emax,
     &      emin_area,emax_area,aa7,ifear,alpha1,
     &      alpha2,bbreak_energy)
c         WRITE(*,'(''Integrated flux correction factor ='',f8.2)') aa8/ak
         reduce=aa8/ak
         IF (lenact(stringin) == 0 ) THEN
            WRITE(*,'("Conversion factor for nuFnu flux at ",f5.2,"keV :",
     &             1pe10.3)') ekev,conv*ekev*ekev*2.418e-12
         ENDIF
      ELSE
cccccccccccccccc
         ifear = 1
c ifear = 1 => flux NOT corrected for Galactic absorption
            ifear = 0
         CALL energycr(1.,alpha,tt,anh,0.,model,5,ak,emin,emax,
     &             emin_area,emax_area,aa7,ifear,alpha1,alpha2,
     &             bbreak_energy)
         fekev=ak*conv
         nufnu = fekev*ekev*ekev*2.418e-12
         reduce=ak
      ENDIF
99001 FORMAT (/,' Model type',/,'   1 -> power law',/,
     &        '   2 -> black body',/,'   3 -> thermal bremsstrahlung',/,
     &        '   7 -> double power law     (d/f=1) : ',$)
99999 CONTINUE
      close(15)
      return
      END


      SUBROUTINE energycr(crate,ggamm,tt,nhydr,hys,iifl,iitype,
     &                 flxog,emin_flux,emax_flux,emin_area,emax_area,
     &                 aa8,ifear,ggamma1,ggamma2,bbreak_energy)
C
C  this subroutine converts count rates into fluxes
C
C  crate - source count rate -
C
C  ggamm - power law energy slope -
C
C  tt    - temperature in kelvin degrees -
C
C  nhydr - hydrogen column density in our galaxy  (real*4) -
C
C  hys   - intrinsic absorption  (real*4) -
C
C  iifl  - flag for spectral type     ifl = 1  :  power law
C                                     ifl = 2  :  black body
C                                     ifl = 3  :  thermal brems
C                                     ifl = 4  :  exponential
C                                     ifl = 7  :  double power law
C
C  iitype- flag for absorption model  itype = 0 : brown & gould
C                                     itype = 1 : gas model
C                                     itype = 2 : 0.15 micron model
C                                     itype = 3 : 0.60 micron model
C                                     itype = 4 : brown & gould (xuv ext
C                                     itype = 5 : morrison & mccammon
C
C  flxog - source flux outside the galaxy ( on source side of galaxy)
C
C  emin_flux, emax_flux lower and upper boundary of energy band where flux is
C             calculated
C
C  ifear = 1 for fluxes at the earth
C  ifear = 0 for fluxes corrected for galactic absorption
C
c
c  emin_area minumum value of energy in effective area table
c  emax_area maximum value of energy in effective area table
c
      implicit none
      EXTERNAL espec
      EXTERNAL eacma
c      integer*4 ifl, iifl, itype, nmax, ifear, n4, n7
      integer*4 ifl, iifl, itype, nmax, ifear
      integer*4 iitype
      REAL*4 nhydr, nh, t, tt, gamm, ggamm, gamma1, ggamma1, gamma2
      real*4 ggamma2, break_energy, bk, hys , aa4, aa8, ra48
      real*4 emin_flux, emax_flux, emin_area, emax_area, flxog, crate
      real*4 bbreak_energy
      real*4 espec,eacma
      COMMON /esp   / nh, gamm, bk, ifl, t, itype, gamma1,
     &                gamma2, break_energy
      ifl = iifl
      t = tt
      itype = iitype
      gamm = ggamm
      gamma1 = ggamma1
      gamma2 = ggamma2
      break_energy = bbreak_energy
      nmax = 200
      bk = 8.6171E-8
      IF ( ifear.EQ.1 ) THEN
          nh = (nhydr+hys)*1.E-22
      ELSE
          nh = hys*1.E-22
      END IF
      CALL stup1(espec,emin_flux,emax_flux,aa4,nmax)
      nh = (nhydr+hys)*1.E-22
      gamm = gamm + 1.
      gamma1 = gamma1 + 1.
      gamma2 = gamma2 + 1.
      IF (crate > 0.) THEN
c crate < 0. forGalactice NH correction only, it avoids integral with effective area that
c is not necessary
         CALL stup1(eacma,emin_area,emax_area,aa8,nmax)
      ELSE
         aa8=1.
      ENDIF
      gamm = gamm - 1.
      gamma1 = gamma1 - 1.
      gamma2 = gamma2 - 1.
      ra48 = aa4/aa8*1.602192E-9
      flxog = ra48*crate
      RETURN
      END

**==eacma.f
      FUNCTION eacma(ener)
      integer*4 itype, ifl
      REAL*4 nh, absor , ener, atten, effa, gamm, eacma, es, bk, t
      real*4 emean, sigma, break_energy, gamma1, gamma2, effective_area
      COMMON /esp   / nh, gamm, bk, ifl, t, itype,
     &                gamma1, gamma2, break_energy
      absor = atten(ener,nh,itype,1.,1.)
      effa = effective_area(ener)
      IF ( ifl.EQ.1 ) THEN
         IF ( gamm.LE.-20. ) THEN
            eacma = absor*ener**20*effa
            RETURN
         END IF
         IF ( gamm.GT.20 ) THEN
            eacma = 0.
            RETURN
         END IF
         eacma = absor*ener**(-gamm)*effa
         RETURN
      ELSE IF ( ifl.EQ.2 ) THEN
         es = ener/(bk*t)
         IF ( es.GT.60. ) THEN
            eacma = absor*effa*exp(2.*alog(ener)-es)
            RETURN
         END IF
         eacma = absor*effa*ener*ener/(exp(es)-1.)
         RETURN
      ELSE IF ( ifl.EQ.3 ) THEN
         es = ener/(bk*t)
         eacma = absor*effa*exp(-.4*alog(es)-alog(ener)-es)
         RETURN
      ELSE IF ( ifl.EQ.4 ) THEN
         es = ener/(bk*t)
         eacma = absor*effa*exp(-es)
         RETURN
      ELSE IF ( ifl.EQ.5 ) THEN
         emean = gamm
         sigma = 3.E-4*emean
         es = (ener-emean)*(ener-emean)/(2.*sigma*sigma)
         IF ( es.GE.88 ) THEN
            eacma = 0.
            RETURN
         ELSE
            eacma = absor*effa/sigma*exp(-es)
            RETURN
         END IF
      ELSE IF ( ifl.EQ.6 ) THEN
         es = ener/bk/t
         eacma = absor*effa*ener**2.12*exp(-es)
         RETURN
      ELSE IF ( ifl.EQ.7 ) THEN
         IF ( ener.LT.break_energy ) THEN
            eacma = absor*ener**(-gamma1)*effa
         ELSE
            eacma = absor*ener**(-gamma2)
     &              *effa*break_energy**(gamma2-gamma1)
         END IF
         RETURN
      END IF
      END

**==STUP1.FOR
      SUBROUTINE stup1(func,xmin,xmax,result,nloops)
c
C ** this rutine calculates definite integrals of the function "func"
C ** in the interval xmin-xmax using the stupid trapezoidal rule.
C ** the only difference between this and other "smarter" routines
C ** is that this one always works while others too often don't!
c
      INTEGER*4 nloops , i , last
      REAL*4 aloo , h , xmax , xmin , first , func , result , xw
      REAL*4 fi , xl
      aloo = nloops - 1
      h = (xmax-xmin)/aloo
      first = func(xmin)/2.
      result = first
      xw = xmin
c
      DO 100 i = 1 , nloops - 2
         xw = xw + h
         fi = func(xw)
         result = result + fi
 100  CONTINUE
c
      xl = xw + h
      last = func(xl)/2.
      result = (result+last)*h
      RETURN
      END


**==ATTEN.FOR
      FUNCTION atten(e,enx,itype,raox,rafe)
c
c  interstellar absorption
c
      INTEGER*4 itype
      REAL*4 ax , enx , e , raox , rafe , axs , atten
      ax = -enx*axs(e,itype,raox,rafe)
      IF ( ax.GE.0 ) THEN
         atten = 1.
         RETURN
      ELSEIF ( ax.GT.-40. ) THEN
         atten = exp(ax)
         RETURN
      ELSE
         atten = 0.
         RETURN
      ENDIF
      END


**==EFFECTIVE_AREA.FOR
      FUNCTION effective_area(e)
c
C   this function returns the effective area at a given energy e
c
      INTEGER*4 j , k , kk
      REAL*4 r , effective_area, e,energy, area
      COMMON /filt  / energy(5000) , area(5000)
      j=1
      DO WHILE ( e.GE.energy(j) )
         k = j
         j=j+1
      ENDDO
      kk = k - 1
      IF ( kk.LE.0 ) THEN
         k = k + 1
         kk = kk + 1
      ENDIF
      r = alog(area(kk)/area(k))/alog(energy(kk)/energy(k))
c      effective_area = area(k)*((e/energy(k))**r)
      effective_area = area(k)
      RETURN
      END

**==AXS.FOR
C=============================================================================
C modified mssl interstellar absorption package  (rjb+pg)
C with addition adapted from new
C
      FUNCTION axs(e,itype,raox,rafe)
C ism cross section
C this function calculates interstellar attenuation coefficients
C itype=0     brown & gould model
C itype=1     fireman gas model
C itype=2     fireman 0.15 micr model
C itype=3     fireman 0.60 micr model
C itype=4     brown & gould (extended to xuv) (likely cruddace)
C itype=5     morrison and mccammon (adapted from new)
C raox is the abundance of oxygen required rel to 'cosmic'
C rafe is the abundance of iron required rel to 'cosmic'
C mods to be made ... correct edge energies .....
C  ni  8.33165
C  fe  7.1112
C  cr  5.9892
C  ca  4.0381
C  a   3.2029
C  s   2.47048
C  si  1.8400
C  mg  1.30339
C  ne  0.866889
C  o   0.5317
C  n   0.4000
C  c   0.28384
C  he  0.0243
C  h   0.0136
      INTEGER*4 j1 , itype , j , k , kk , i
      REAL*4 e , ex , r , xs , sig , axs , raox , oxy , eng , absr2 ,
     &       coef , e1 , xd
      REAL*4 sg , dex , rafe , fe , fct
      DIMENSION ex(40) , sig(40) , e1(22) , sg(22,3)
      DIMENSION coef(15,3) , eng(16)
      DATA ex/.0136 , .0177 , .02 , .0243 , .0243 , .03 , .035 , .045 ,
     &     .07 , .1 , .2 , .2838 , .2838 , .4016 , .4016 , .5320 ,
     &     .5320 , .8669 , .8669 , 1.3050 , 1.3050 , 1.8389 , 1.8389 ,
     &     2.4720 , 2.4720 , 3.2029 , 3.2029 , 4. , 4. , 5. , 5. , 6. ,
     &     6. , 7. , 7. , 8. , 8. , 9. , 10. , 12./
                                   !mmcc
      DATA sig/.15 , .172 , .175 , .17 , .25 , .32 , .35 , .41 , .51 ,
     &     .58 , .63 , .60 , .70 , .67 , .75 , .74 , 1.63 , 1.72 ,
     &     2.01 , 2.16 , 2.28 , 2.4 , 2.74 , 2.89 , 3.3 , 3.45 , 3.71 ,
     &     3.86 , 3.86 , 3.92 , 3.92 , 3.95 , 3.95 , 3.95 , 3.95 ,
     &     3.94 , 3.94 , 3.92 , 3.90 , 3.90/
      DATA e1/.3 , .5320 , .5320 , .8669 , .8669 , 1.3050 , 1.3050 ,
     &     1.8389 , 1.8389 , 2.4720 , 2.4720 , 3.2029 , 3.2029 ,
     &     4.0381 , 4.0381 , 5.9892 , 5.9892 , 7.1120 , 7.1120 ,
     &     8.3328 , 8.3328 , 10./
      DATA sg/0.8 , 0.9 , 1.9 , 2.0 , 2.4 , 2.6 , 2.8 , 2.9 , 3.2 ,
     &     3.5 , 3.9 , 4.1 , 4.3 , 4.6 , 4.7 , 5.0 , 5.1 , 5.3 , 10.6 ,
     &     11.2 , 11.8 , 12.3 , 0.7 , 0.9 , 1.4 , 1.7 , 2.0 , 2.4 ,
     &     2.6 , 2.8 , 3.0 , 3.4 , 3.9 , 4.1 , 4.3 , 4.6 , 4.7 , 5.0 ,
     &     5.1 , 5.3 , 10.6 , 11.2 , 11.8 , 12.3 , 0.6 , 0.7 , 0.8 ,
     &     1.1 , 1.5 , 1.9 , 2.1 , 2.4 , 2.7 , 3.3 , 3.9 , 4.1 , 4.3 ,
     &     4.6 , 4.7 , 5.0 , 5.1 , 5.3 , 10.6 , 11.2 , 11.8 , 12.3/
      DATA coef/17.3 , 34.6 , 78.1 , 71.4 , 95.5 , 308.9 , 120.6 ,
     &     141.3 , 202.7 , 342.7 , 352.2 , 433.9 , 629.0 , 701.2 ,
     &     953.0 , 608.1 , 267.9 , 18.8 , 66.8 , 145.8 , -380.6 ,
     &     169.3 , 146.8 , 104.7 , 18.7 , 18.7 , -2.4 , 30.9 , 25.2 ,
     &     0.0 , -2150. , -476.1 , 4.3 , -51.4 , -61.1 , 294.0 , -47.7 ,
     &     -31.5 , -17.0 , 0.0 , 0.0 , 0.75 , 0.0 , 0.0 , 0.0/
      DATA eng/0.03 , 0.1 , 0.284 , 0.4 , 0.532 , 0.707 , 0.867 ,
     &     1.303 , 1.84 , 2.471 , 3.21 , 4.038 , 7.111 , 8.331 , 10.0 ,
     &     100.0/
      j1 = 9
      IF ( itype.EQ.0 .OR. itype.EQ.4 ) THEN
C
C.....brown and gould ( 0.0136 or 0.07 to 12 kev) ... only o supported........
C
         IF ( itype.EQ.4 ) j1 = 1
         DO 50 j = j1 , 40
            k = j
            IF ( e.LE.ex(j) ) GOTO 100
 50      CONTINUE
 100     kk = k - 1
         IF ( kk.LE.0 ) THEN
            k = k + 1
            kk = kk + 1
         ENDIF
         r = (e-ex(kk))/(ex(k)-ex(kk))
         xs = r*(sig(k)-sig(kk)) + sig(kk)
         xs = xs/e/e/e
         xs = xs
         axs = xs + 1.082*(raox-1.)*oxy(e)
         RETURN
      ELSEIF ( itype.EQ.5 ) THEN
C
C.....morrison and mc cammon ( 0.03-100 kev) ... no o and fe supported........
C     in units of 10**22
C     also added in is thomson xsection
C
         DO 150 i = 1 , 16
            IF ( e.LT.eng(i) ) GOTO 200
 150     CONTINUE
 200     i = i - 1
         IF ( i.EQ.0 ) i = 1
         absr2 = (coef(i,1)+coef(i,2)*e+coef(i,3)*e*e)/e**3
                      !protection added by lc (le range)
         axs = absr2*0.01 + 0.0067
         RETURN
      ELSE
C
C.....fireman ( 0.3-10 kev ) ... o and fe supported ..........................
C
         DO 250 j = 1 , 22
            k = j
            IF ( e.LE.e1(j) ) GOTO 300
 250     CONTINUE
 300     kk = k - 1
         IF ( kk.LE.0 ) THEN
            k = k + 1
            kk = kk + 1
         ENDIF
         r = (e-e1(kk))/(e1(k)-e1(kk))
         xd = r*(sg(k,itype)-sg(kk,itype)) + sg(kk,itype)
         xd = xd/e/e/e
         xd = xd
         dex = 1.022*(rafe-1.)*fe(e)
         axs = xd + 1.216*(raox-1.)*oxy(e)*fct(e,itype) + dex
         RETURN
      ENDIF
      END

**==OXY.FOR
      FUNCTION oxy(e)
c
c  oxygen
c
      integer*4 j,k,kk
      real*4 eox(8), sox(8), r, e , oxy
      DATA sox, eox/1.79, 2.74, 3.30, 85.53, 104.36, 123.76, 144.71,
     &     155.37, 0.10, 0.30, 0.5320, 0.5320, 1., 2., 5., 10./
      DO j = 1, 8
         k = j
         IF ( e.LE.eox(j) ) GO TO 100
      END DO
 100  CONTINUE
      kk = k - 1
      IF ( kk.LE.0 ) THEN
         k = k + 1
         kk = kk + 1
      END IF
      r = (e-eox(kk))/(eox(k)-eox(kk))
      oxy = r*(sox(k)-sox(kk)) + sox(kk)
      oxy = oxy/e/e/e
      oxy = oxy*1.E-2
      RETURN
      END


**==FE.FOR
      FUNCTION fe(e)
c
c  iron
c
      integer*4 j,k,kk
      real*4 efe(12), sfe(12), e, fe, r
      DATA sfe, efe/2.15, 7.74, 23.9, 34.9, 50.6, 59.7, 69.7, 74.5,
     &     593., 626.9, 686.1, 772.5, .3, .6, 1., 1.5, 2.4, 3.5, 5.0,
     &     7.112, 7.112, 8., 10., 15./
      IF ( e.LE..2 ) THEN
         fe = 0.
         RETURN
      ELSE
         DO j = 1, 12
            k = j
            IF ( e.LE.efe(j) ) GO TO 50
         END DO
 50      CONTINUE
         kk = k - 1
         IF ( kk.LE.0 ) THEN
            k = k + 1
            kk = kk + 1
         END IF
         r = (e-efe(kk))/(efe(k)-efe(kk))
         fe = r*(sfe(k)-sfe(kk)) + sfe(kk)
         fe = fe/e/e/e
         fe = fe*1.E-2
         RETURN
      END IF
      END

**==FCT.FOR
      FUNCTION fct(e,itype)
c
c  auxiliary ism routine
c
      integer*4 j,k,kk,itype
      real*4 eq(7), f(7,3) , e , r , fct
      DATA f/7*1., .155, .70, .86, .93, .98, 1., 1., .187, .317, .587,
     &     .77, .93, 1., 1./
      DATA eq/.3, .6, 1., 1.5, 2.4, 3., 4./
      DO j = 1, 7
         k = j
         IF ( e.LE.eq(j) ) GO TO 100
      END DO
 100  CONTINUE
      kk = k - 1
      IF ( kk.LE.0 ) THEN
         k = k + 1
         kk = kk + 1
      END IF
      r = (e-eq(kk))/(eq(k)-eq(kk))
      fct = r*(f(k,itype)-f(kk,itype)) + f(kk,itype)
      IF ( fct.LT.0. ) fct = 0.
      RETURN
      END

**==ESPEC.FOR
      FUNCTION espec(ener)
      INTEGER*4 itype , ifl
      REAL*4 nh , absor , ener , atten , espec , gamm , es , bk , t ,
     &       emean
      REAL*4 sigma , break_energy , gamma1 , gamma2
      COMMON /esp   / nh , gamm , bk , ifl , t , itype , gamma1 ,
     &                gamma2 , break_energy
      absor = atten(ener,nh,itype,1.,1.)
      IF ( ifl.EQ.1 ) THEN
         espec = absor*ener**(-gamm)
         RETURN
      ELSEIF ( ifl.EQ.2 ) THEN
         es = ener/(bk*t)
         IF ( es.GT.60. ) THEN
            espec = absor*exp(3*alog(ener)-es)
            RETURN
         ENDIF
         espec = absor*ener*ener*ener/(exp(es)-1.)
         RETURN
      ELSEIF ( ifl.EQ.3 ) THEN
         es = ener/(bk*t)
         espec = absor*exp(-.4*alog(es)-es)
         RETURN
      ELSEIF ( ifl.EQ.4 ) THEN
         es = ener/(bk*t)
         espec = absor*ener*exp(-es)
         RETURN
      ELSEIF ( ifl.EQ.5 ) THEN
         emean = gamm
         sigma = 3.E-4*emean
         es = (ener-emean)*(ener-emean)/(2.*sigma*sigma)
         IF ( es.GE.88 ) THEN
            espec = 0.
            RETURN
         ELSE
            espec = absor*ener/sigma*exp(-es)
            RETURN
         ENDIF
      ELSEIF ( ifl.EQ.6 ) THEN
         es = ener/(bk*t)
         IF ( es.GT.60 ) THEN
            espec = ener**3.12
         ELSE
            espec = absor*ener**3.12*exp(-es)
         ENDIF
         RETURN
      ELSEIF ( ifl.EQ.7 ) THEN
         IF ( ener.LT.break_energy ) THEN
            espec = absor*ener**(-gamma1)
         ELSE
            espec = absor*ener**(-gamma2)*break_energy**(gamma2-gamma1)
         ENDIF
         RETURN
      ELSE
         WRITE (*,'('' sub ESPEC called with a wrong flag '')')
         RETURN
      ENDIF
      END

c=========================EBL

      SUBROUTINE ebl_flux(emin,emax,ggg,zzz,reduction_factor)
c
c
c Calculates ebl attenuation factor for VHE flux from spectra integrated between emin (TeV) and emax (TeV)
c
      IMPLICIT NONE
      INTEGER*4 i , lu_input, n
      INTEGER*4 j, length,lenact
      INTEGER*4 nmax
      REAL*4 redshift,binz(39),e(50),tau(50,39)
      REAL*4 emin, emax,ggg,zzz
      REAL*4 gamma1,gamma2,break_energy,out1,out2,reduction_factor
      REAL*4 accur, r
      REAL*4 delta_slope,tau_ebl
      CHARACTER*200 input_file,webprograms,string
c      DATA binz/ 0.01, 0.02526316, 0.04052632, 0.05078947, 0.07105263,
c     & 0.08631579, 0.10157895, 0.11684211, 0.13210526, 0.14736842,
c     & 0.16263158, 0.17789474, 0.19315789, 0.20842105, 0.22368421,
c     & 0.23894737, 0.25421053, 0.26947368, 0.28473684, 0.3, 0.35, 0.4,
c     & 0.45, 0.5, 0.50, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
c     &  1.,1.2,1.4,1.6,1.8,2./
      EXTERNAL especdebl
      EXTERNAL noabs
      real*4 especdebl,noabs
      COMMON /espebl   / gamma1,gamma2,break_energy,redshift
      COMMON /bins  / binz,e,tau
      common webprograms
      DATA binz/ 0.01, 0.02526316, 0.04052632, 0.05078947, 0.07105263,
     & 0.08631579, 0.10157895, 0.11684211, 0.13210526, 0.14736842,
     & 0.16263158, 0.17789474, 0.19315789, 0.20842105, 0.22368421,
     & 0.23894737, 0.25421053, 0.26947368, 0.28473684, 0.3, 0.35, 0.4,
     & 0.45, 0.5, 0.50, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
     &  1.,1.2,1.4,1.6,1.8,2./
      input_file=webprograms(1:lenact(webprograms)) //
     &'/count_rate/tau_dominguez11.txt'

c      CALL getlun(lu_input)
      OPEN (17,file=input_file,status='old')
c      DO i = 1,5
c         read (10,'(a)') string
c      ENDDO
      DO i = 1,50
      read (17,*,end=999) e(i),tau(i,1:39)
c      write(*,*) i,e(i)
      ENDDO
999   CONTINUE
      CLOSE(17)
      redshift=zzz
c      write(*,*) redshift,emin,emax,ggg
      gamma1 = ggg - 1.
      delta_slope=0.5
      gamma2 = gamma1 + delta_slope
      nmax = 2000
      accur = 1.e-2
      break_energy = 0.2
c      break_energy = 0.1
      IF (redshift .eq. 0.e0) redshift = 0.3
c      write(*,*) redshift
      IF (redshift .LT. 0.01 ) redshift = 0.01
      IF (redshift .GT. 2.0 )  redshift = 2.0
      if (emin .eq. emax) then
         reduction_factor=1./exp(-tau_ebl(emin))
      else
         CALL  simp1(especdebl,emin,emax,out1,r,n,accur,nmax)
         CALL  simp1(noabs,emin,emax,out2,r,n,accur,nmax)
         reduction_factor=out2/out1 !no abs / abs
      endif
c      write (*,*) ' Reduction factor ' , reduction_factor
      return
      END



      FUNCTION tau_ebl(energy)
c Calculates the amount of obsorption due to EBL given energy and redshift
      IMPLICIT none
      INTEGER*4 i, j
      REAL*4 energy, redshift, binz(39),e(50),tau(50,39),tau_ebl
      REAL*4 interpolated_tau,x0,x1,y0,y1,y0ie,y1ie
      REAL*4 gamma1,gamma2,break_energy
      LOGICAL ok
      COMMON /espebl   / gamma1,gamma2,break_energy,redshift
      COMMON /bins  / binz,e,tau
      IF (energy.LT.e(1) ) THEN
      tau_ebl = 0.0
      RETURN
      ENDIF
      IF (energy.LT.e(1) .OR. energy.GT.e(50) ) THEN
      write (*,*) 'Energy ',energy,' outside allowed range (',e(1),e(50),')'
      stop
      endif
      ok = .TRUE.
      i=1
      DO WHILE (ok)
      i = i +1
      IF (energy.LE.e(i)) ok = .FALSE.
      ENDDO
      IF (redshift.LT.binz(1)) THEN
      j = 1
      ELSE IF (redshift.GT.binz(39)) THEN
      j = 39
      ELSE
      j=1
      ok = .TRUE.
      DO WHILE (ok)
      j = j +1
      IF (redshift.LE.binz(j)) ok = .FALSE.
      ENDDO
      ENDIF
      y0 = tau(i-1,j-1)
      y1 = tau(i,j-1)
      x0 = e(i-1)
      x1 = e(i)
      y0ie = y0 + (y1-y0)*(energy-x0)/(x1-x0)
      y0 = tau(i-1,j)
      y1 = tau(i,j)
      y1ie = y0 + (y1-y0)*(energy-x0)/(x1-x0)
      x0 = binz(j-1)
      x1 = binz(j)
      interpolated_tau = y0ie + (y1ie-y0ie)*(redshift-x0)/(x1-x0)
c       write (*,*) 'tau(',e(i-1),binz(j-1),'),tau(',e(i-1),binz(j),')',
c     &              tau(i-1,j-1),tau(i-1,j)
c       write (*,*) 'tau(',e(i),binz(j-1),'),  tau(',e(i),binz(j),')  ',
c     &              tau(i,j-1),tau(i,j)
c       write (*,*) ' interpolated_tau ', interpolated_tau
      tau_ebl = interpolated_tau
      RETURN
      END
c
c
c
      FUNCTION especdebl(energy)
      REAL*4 absor , energy , especdebl , break_energy , gamma1 , gamma2
      REAL*4 tau_ebl
      REAL*4 redshift, binz(39),e(50),tau(50,39)
      COMMON /espebl   / gamma1,gamma2,break_energy,redshift
      COMMON /bins  / binz,e,tau
      absor = exp(-tau_ebl(energy))
c      write (*,*) ' absor energy ', absor, energy
      IF ( energy.LT.break_energy ) THEN
      especdebl = absor*energy**(-gamma1)
      ELSE
      especdebl = absor*energy**(-gamma2)*break_energy**(gamma2-gamma1)
      ENDIF
      RETURN
      END
c
c
c
      FUNCTION noabs(energy)
      REAL*4 energy , noabs , break_energy , gamma1 , gamma2, redshift
      COMMON /espebl   / gamma1,gamma2,break_energy,redshift
      IF ( energy.LT.break_energy ) THEN
      noabs = energy**(-gamma1)
      ELSE
      noabs = energy**(-gamma2)*break_energy**(gamma2-gamma1)
      ENDIF
      RETURN
      END
c
**==SIMP1.FOR
c
      SUBROUTINE simp1(akern,a,b,aint,r,n,accur,nmax)
C
C      Calculate definite integrals of a function
C      using Simpson's rule
C      Inputs :
C      akern  : function to be integrated, must be a function subroutine
C               and must be declared 'external' in main program.
C      a,b    : limits of integration
C      accur  : relative accuracy required
C      nmax   : max number of iterations allowed
C      Outputs:
C      aint   : integral of function
C      r      : actual relative accuracy achieved
C      n      : actual no of iterations performed
C      ERROR messages
C      exit 14 : no convergence to accuracy requird achieved
C      exit 15 : a > b
      INTEGER*4 n , i , nmax
      REAL*4 dh , a , b , aint , akern , aone , atwo , afour , t ,
     &       oldint , r , accur
      dh = (b-a)/4.D0
      IF ( dh.EQ.0.D0 ) THEN
      aint = 0.D0
      RETURN
      ENDIF
      IF ( dh.LT.0.D0 ) THEN
      WRITE (*,*) ' Simp1 exit 15 , a > b ! ' , dh
      STOP
      ENDIF
      n = 2
      aone = akern(a) + akern(b)
      atwo = akern(a+dh+dh)
      afour = akern(a+dh) + akern(b-dh)
      aint = dh/3.D0*(aone+2.*atwo+4.D0*afour)
100   dh = dh/2.D0
      n = n + n
      atwo = atwo + afour
      afour = 0.D0
      t = a - dh
      DO 200 i = 1 , n
      t = t + 2.D0*dh
      afour = afour + akern(t)
200   CONTINUE
      oldint = aint
      aint = dh/3.D0*(aone+2.D0*atwo+4.D0*afour)
      r = abs((aint-oldint)/aint)
      IF ( r.LE.accur ) RETURN
      IF ( n.LT.nmax ) GOTO 100
      WRITE (*,*) ' Simp1 exit 14, accuracy required not achieved ' ,
     &            n , aint , r
      STOP
      END
