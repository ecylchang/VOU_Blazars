      PROGRAM plot_lc

c This program plot the SED for candidate

      implicit none
      integer*4 ier,pgbeg,length,ns,j,rah, ram, id, dm,in,im,pttest
      integer*4 i,sfound,npt(5000),rtype,stype(5000),spectype(5000,1000)
      real*4 frequency(5000,1000),flux(5000,1000),uflux(5000,1000),lflux(5000,1000)
      real*4 rasec,decsec,testflux,mjdstart(5000,1000),mjdend(5000,1000),mjdavg(5000,1000)
      real*4 mjdlow(2),mjdup(2),lcup(2),lclow(2)
      real*8 rra,rdec,ra(1000),dec(1000)
      character*160 string
      character*100 title
      character*80 input_file,output_file,refs(5000,1000)
      character*14 stringin
      character*6 number
      character*1 sign
      logical ok,there
      ok = .true.

      call rdforn(string,length)
      call rmvlbk(string)
      in=index(string(1:length),' ')
      input_file=string(1:in-1)
      im=index(string(in+1:length),' ')+in
      output_file=string(in+1:im-1)
      number=string(im+1:length)
      !write(*,*) number
      if (number == 'finish') then
         Stop '!!!!Exit the source exploring routine!!!!'
      else if (number == 'sed') then
         ns=99
      else
         read(number,*) ns
      endif

c      write(*,*) output_file
      INQUIRE (FILE=input_file,EXIST=there)
      IF (.NOT.there) THEN
         write (*,'('' file Sed.txt not found '')')
         STOP
      ENDIF

      npt(1:1000)=0
      open(10,file=input_file,status='old')
      !if (ns .ne. 1) .and. read(10,'(a)') stringin
101   continue
      read(10,'(i4,2x,a,2(2x,f9.5),2x,i2)',end=100,err=100) sfound,stringin,rra,rdec,rtype
      ra(sfound)=rra
      dec(sfound)=rdec
      stype(sfound)=rtype
      do while (ok)
         npt(sfound)=npt(sfound)+1
         read(10,*,end=99,err=99) frequency(npt(sfound),sfound),flux(npt(sfound),sfound),
     &      uflux(npt(sfound),sfound),lflux(npt(sfound),sfound),mjdstart(npt(sfound),sfound),
     &      mjdend(npt(sfound),sfound),spectype(npt(sfound),sfound)
c     &      refs(npt(sfound),sfound)
      enddo
99    continue
      npt(sfound)=npt(sfound)-1
      write(*,*) 'number of pts read for light curve:',npt(sfound)!,spectype(npt(sfound),sfound)
      goto 101
100   continue

      if (npt(sfound) .eq. 0) STOP 'No data for SED'

      if (ns .eq. 99) then
         i=1
      else
         i=ns
      endif
      call chra(ra(i),rah,ram,rasec,1)
      call chdec(abs(dec(i)),id,dm,decsec,1)
      sign(1:1) =' '
      if (dec(i) < 0.0) sign(1:1)='-'

      !write(*,*) i
c      IER = PGBEG(0,"/xwindow",1,1)
c      IER = PGBEG(0,"/xs",1,1)
      IER = PGBEG(0,output_file,1,2)
      call pgslw(4)
c      if ( ns .eq. 99 ) THEN
c         write(title,'(a,i2.2,1x,i2.2,1x,f4.1,a,a,i2.2,1x,i2.2,1x,f4.1)')
c     &     'Source position:  ',rah,ram,rasec,' , ',sign,id,dm,decsec
c      else
c         write(title,'(a,i2,a,i2.2,1x,i2.2,1x,f4.1,a,a,i2.2,1x,i2.2,1x,f4.1)')
c     &     'Source',i,' position:  ',rah,ram,rasec,' , ',sign,id,dm,decsec
c      endif
      CALL PGSCRN(0, 'White', IER)
      CALL PGSCRN(1, 'Black', IER)
      call pgsch(1.3)
      mjdavg=(mjdstart+mjdend)/2.
c      write(*,*) mjdavg(1:8,1)

c define the upper limit and lower limit for LC
      lcup(1)=max(uflux(1,i),abs(flux(1,i)))
      lclow(1)=min(lflux(1,i),abs(flux(1,i)))
      mjdlow(1)=mjdstart(1,i)
      mjdup(1)=mjdend(1,i)
      do j=1,npt(i)
         if (spectype(j,i) .eq. 11) then
            testflux=max(uflux(j,i),abs(flux(j,i)))
            if ((testflux .gt. lcup(1)) .and. (testflux .gt. 0.d0)) lcup(1)=testflux
            testflux=min(lflux(j,i),abs(flux(j,i)))
            if ((testflux .lt. lclow(1)) .and. (testflux .gt. 0.d0)) lclow(1)=testflux
            if (mjdstart(j,i) .lt. mjdlow(1)) mjdlow(1)=mjdstart(j,i)
            if (mjdend(j,i) .gt. mjdup(1)) mjdup(1)=mjdend(j,i)
         else if (spectype(j,i) .eq. 00) then
            if (spectype(j,i) .ne. spectype(j-1,i)) then
               lcup(2)=max(uflux(j,i),abs(flux(j,i)))
               lclow(2)=min(lflux(j,i),abs(flux(j,i)))
               mjdup(2)=mjdend(j,i)
               mjdlow(2)=mjdstart(j,i)
            endif
            testflux=max(uflux(j,i),abs(flux(j,i)))
            if ((testflux .gt. lcup(2)) .and. (testflux .gt. 0.d0)) lcup(2)=testflux
            testflux=min(lflux(j,i),abs(flux(j,i)))
            if ((testflux .lt. lclow(2)) .and. (testflux .gt. 0.d0)) lclow(2)=testflux
            if (mjdstart(j,i) .lt. mjdlow(2)) mjdlow(2)=mjdstart(j,i)
            if (mjdend(j,i) .gt. mjdup(2)) mjdup(2)=mjdend(j,i)
         endif
      enddo
      lcup=alog10(lcup)+0.5
      lclow=alog10(lclow)-0.5
      mjdlow=mjdlow-20.
      mjdup=mjdup+20.

      CALL PGENV(mjdlow(1),mjdup(1),lclow(1),lcup(1),0,1)
      CALL PGLAB('', 'Log \gnf\d\gn\u (erg/s/cm\u2\d)','1 keV')
      call pgsci(4)
      do j=1,npt(i)
         if (spectype(j,i) .eq. 11.) then
            if ((flux(j,i) .eq. lflux(j,i)) .and. (flux(j,i) .eq. uflux(j,i))) then
               call pgsch(1.5)
               CALL PGPT(1,mjdavg(j,i),log10(uflux(j,i)),45)
               call PGPT(1,mjdavg(j,i),log10(uflux(j,i))-0.07,31)
               CALL PGERRX(1,mjdstart(j,i),mjdend(j,i),log10(uflux(j,i)),1.0)
            else if ((lflux(j,i) .eq. 0.) .and. (uflux(j,i) .ne. 0.) ) then
               call pgsch(1.5)
               CALL PGPT(1,mjdavg(j,i),log10(uflux(j,i)),45)
               call PGPT(1,mjdavg(j,i),log10(uflux(j,i))-0.07,31)
               CALL PGERRX(1,mjdstart(j,i),mjdend(j,i),log10(uflux(j,i)),1.0)
            else if (flux(j,i) .lt. 0.) then
               call pgsch(1.2)
               CALL PGPT(1,mjdavg(j,i),log10(-flux(j,i)),13)
               CALL PGERRY(1,mjdavg(j,i),log10(uflux(j,i)),log10(lflux(j,i)),1.0)
               CALL PGERRX(1,mjdstart(j,i),mjdend(j,i),log10(-flux(j,i)),1.0)
            else
               call pgsch(1.2)
               CALL PGPT(1,mjdavg(j,i),log10(flux(j,i)),-17)
               CALL PGERRY(1,mjdavg(j,i),log10(uflux(j,i)),log10(lflux(j,i)),1.0)
               CALL PGERRX(1,mjdstart(j,i),mjdend(j,i),log10(flux(j,i)),1.0)
            endif
         endif
      enddo

c second panel for LC
      call pgsci(1)
      call pgsch(1.3)
      CALL PGENV(mjdlow(2),mjdup(2),lclow(2),lcup(2),0,1)
      CALL PGLAB('MJD', 'Log \gnf\d\gn\u (erg/s/cm\u2\d)','1 TeV')
      call pgsci(8)
      do j=1,npt(i)
         if (spectype(j,i) .eq. 00.) then
            if ((flux(j,i) .eq. lflux(j,i)) .and. (flux(j,i) .eq. uflux(j,i))) then
               call pgsch(1.5)
               CALL PGPT(1,mjdavg(j,i),log10(uflux(j,i)),45)
               call PGPT(1,mjdavg(j,i),log10(uflux(j,i))-0.07,31)
               CALL PGERRX(1,mjdstart(j,i),mjdend(j,i),log10(uflux(j,i)),1.0)
            else if ((lflux(j,i) .eq. 0.) .and. (uflux(j,i) .ne. 0.) ) then
               call pgsch(1.5)
               CALL PGPT(1,mjdavg(j,i),log10(uflux(j,i)),45)
               call PGPT(1,mjdavg(j,i),log10(uflux(j,i))-0.07,31)
               CALL PGERRX(1,mjdstart(j,i),mjdend(j,i),log10(uflux(j,i)),1.0)
            else if (flux(j,i) .lt. 0.) then
               call pgsch(1.2)
               CALL PGPT(1,mjdavg(j,i),log10(-flux(j,i)),13)
               CALL PGERRY(1,mjdavg(j,i),log10(uflux(j,i)),log10(lflux(j,i)),1.0)
               CALL PGERRX(1,mjdstart(j,i),mjdend(j,i),log10(-flux(j,i)),1.0)
            else
               call pgsch(1.2)
               CALL PGPT(1,mjdavg(j,i),log10(flux(j,i)),-17)
               CALL PGERRY(1,mjdavg(j,i),log10(uflux(j,i)),log10(lflux(j,i)),1.0)
               CALL PGERRX(1,mjdstart(j,i),mjdend(j,i),log10(flux(j,i)),1.0)
            endif
         endif
      enddo
      call pgsci(1)
      CALL PGEND

      close (10)
      END

