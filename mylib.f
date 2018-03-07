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
      CHARACTER*200 str2
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
      string = xdelta(Dec,String)
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
