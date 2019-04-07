      PROGRAM plot_lc

c This program plot the SED for candidate

      implicit none
      integer*4 ier,pgbeg,length,ns,j,rah, ram, id, dm,in,im,pttest
      integer*4 i,sfound,npt(5000),rtype,stype(5000),ivhe,ixray,iilc
      real*4 frequency(5000,1000),flux(5000,1000),uflux(5000,1000),lflux(5000,1000)
      real*4 rasec,decsec,testflux,mjdstart(5000,1000),mjdend(5000,1000),mjdavg(5000)
      real*4 mjdlow(2),mjdup(2),lcup(2),lclow(2),flux_lc(5000),uflux_lc(5000),lflux_lc(5000)
      real*4 mjdst_lc(5000),mjded_lc(5000),freq_lc(5000),fq1tev,sloperat
      real*8 rra,rdec,ra(1000),dec(1000)
      character*160 string
      character*100 title
      character*80 input_file,output_file,refs(5000,1000)
      character*14 stringin
      character*10 spectype(5000,1000)
      character*6 number,xtitle
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

c define the upper limit and lower limit for LC
c need to define the range first then plot...
      ixray=0
      ivhe=0
      iilc=0
      do j=1,npt(i)
         if ((spectype(j,i) == 'OUSXB') .and. (frequency(j,i) .eq. 2.418E17)) then
            ixray=ixray+1
            if (ixray .eq. 1) then
               lcup(1)=max(uflux(j,i),abs(flux(j,i)))
               lclow(1)=min(lflux(j,i),abs(flux(j,i)))
               mjdlow(1)=mjdstart(j,i)
               mjdup(1)=mjdend(j,i)
            endif
            iilc=iilc+1
            testflux=max(uflux(j,i),abs(flux(j,i))) ! range for 1kev
            if ((testflux .gt. lcup(1)) .and. (testflux .gt. 0.d0)) lcup(1)=testflux
            testflux=min(lflux(j,i),abs(flux(j,i)))
            if ((testflux .lt. lclow(1)) .and. (testflux .gt. 0.d0)) lclow(1)=testflux
            if (mjdstart(j,i) .lt. mjdlow(1)) mjdlow(1)=mjdstart(j,i)
            if (mjdend(j,i) .gt. mjdup(1)) mjdup(1)=mjdend(j,i)
            freq_lc(iilc)=frequency(j,i)
            flux_lc(iilc)=flux(j,i)
            lflux_lc(iilc)=lflux(j,i)
            uflux_lc(iilc)=uflux(j,i)
            mjdst_lc(iilc)=mjdstart(j,i)
            mjded_lc(iilc)=mjdend(j,i)
         else if (spectype(j,i) == 'VERITAS') then
c begin to find 1TeV point for every era
c   frequency_vhe(ivhe)=(1.602E-19)*(frequency_vhe(ivhe)*1.e12)/(6.626e-34)
            fq1tev=(1.602E-19)*(1.e12)/(6.626e-34)
            ivhe=ivhe+1
            if (ivhe .eq. 1) then
               lcup(2)=max(uflux(j,i),abs(flux(j,i)))
               lclow(2)=min(lflux(j,i),abs(flux(j,i)))
               mjdup(2)=mjdend(j,i)
               mjdlow(2)=mjdstart(j,i)
            endif
            if (frequency(j,i) .gt. fq1tev) then
               if (ivhe .eq. 1) then
                  iilc=iilc+1
                  sloperat=((flux(j+1,i)/frequency(j+1,i))-(flux(j,i)/frequency(j,i)))/
     &                       (frequency(j+1,i)/frequency(j,i))
                  flux_lc(iilc)=(flux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                  flux_lc(iilc)=flux_lc(iilc)*fq1tev
                  uflux_lc(iilc)=(uflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                  uflux_lc(iilc)=uflux_lc(iilc)*fq1tev
                  lflux_lc(iilc)=(lflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                  lflux_lc(iilc)=lflux_lc(iilc)*fq1tev
                  mjdst_lc(iilc)=mjdstart(j,i)
                  mjded_lc(iilc)=mjdend(j,i)
               else
                  if (mjdstart(j,i) .ne. mjdstart(j-1,i)) then
                     if (iilc-ixray .ne. 0) then
                        if (mjdstart(j,i) .eq. mjdst_lc(iilc)) then
                           goto 600
                        endif
                     endif
                     iilc=iilc+1
                     sloperat=((flux(j+1,i)/frequency(j+1,i))-(flux(j,i)/frequency(j,i)))/
     &                       (frequency(j+1,i)/frequency(j,i))
                     flux_lc(iilc)=(flux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     flux_lc(iilc)=flux_lc(iilc)*fq1tev
                     uflux_lc(iilc)=(uflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     uflux_lc(iilc)=uflux_lc(iilc)*fq1tev
                     lflux_lc(iilc)=(lflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     lflux_lc(iilc)=lflux_lc(iilc)*fq1tev
                     mjdst_lc(iilc)=mjdstart(j,i)
                     mjded_lc(iilc)=mjdend(j,i)
                  else
                     if (iilc-ixray .ne. 0) then
                        if (mjdstart(j,i) .eq. mjdst_lc(iilc)) then
                           goto 600
                        endif
                     endif
                     iilc=iilc+1
                     sloperat=((flux(j,i)/frequency(j,i))-(flux(j-1,i)/frequency(j-1,i)))
     &                              /(frequency(j,i)/frequency(j-1,i))
                     flux_lc(iilc)=(flux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     flux_lc(iilc)=flux_lc(iilc)*fq1tev
                     uflux_lc(iilc)=(uflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     uflux_lc(iilc)=uflux_lc(iilc)*fq1tev
                     lflux_lc(iilc)=(lflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     lflux_lc(iilc)=lflux_lc(iilc)*fq1tev
                     mjdst_lc(iilc)=mjdstart(j,i)
                     mjded_lc(iilc)=mjdend(j,i)
                  endif
               endif
c               write(*,*) ixray,iilc,flux_lc(iilc),uflux_lc(iilc),mjdst_lc(iilc),mjded_lc(iilc)
               testflux=max(uflux_lc(iilc),flux_lc(iilc))
               if ((testflux .gt. lcup(2)) .and. (testflux .gt. 0.d0)) lcup(2)=testflux
               testflux=min(lflux_lc(iilc),flux_lc(iilc))
               if ((testflux .lt. lclow(2)) .and. (testflux .gt. 0.d0)) lclow(2)=testflux
               if (mjdst_lc(iilc) .lt. mjdlow(2)) mjdlow(2)=mjdst_lc(iilc)
               if (mjded_lc(iilc) .gt. mjdup(2)) mjdup(2)=mjdst_lc(iilc)
            endif
  600       continue
         endif
      enddo

      if (iilc .eq. 0) then
         write(*,*) 'NO Light Curve Plot'
         stop
      endif
      write(*,*) 'number of pts read for light curve:',iilc!,spectype(npt(sfound),sfound)
      if (ixray .eq. 0) then
         mjdlow(1)=0
         mjdup(1)=0
         lcup(1)=0
         lclow(1)=0
      endif
      if (iilc-ixray .eq. 0) then
         mjdlow(2)=0
         mjdup(2)=0
         lcup(2)=0
         lclow(2)=0
      endif
      lcup=alog10(lcup)+0.5
      lclow=alog10(lclow)-0.5
      mjdlow=mjdlow-20.
      mjdup=mjdup+20.
      mjdavg=(mjdst_lc+mjded_lc)/2.
c      write(*,*) mjdlow(1),mjdup(1),mjdlow(2),mjdup(2)

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
c      write(*,*) mjdavg(1:8,1)
      if (ixray .ne. 0) then
      CALL PGENV(mjdlow(1),mjdup(1),lclow(1),lcup(1),0,1)
      if (iilc-ixray .eq. 0) then
         xtitle='MJD'
      else
         xtitle=''
      endif
      CALL PGLAB(xtitle, 'Log \gnf\d\gn\u (erg/s/cm\u2\d)','1 keV')
      call pgsci(4)
      do j=1,iilc
         if (j .le. ixray) then
            if ((flux_lc(j) .eq. lflux_lc(j)) .and. (flux_lc(j) .eq. uflux_lc(j))) then
               call pgsch(1.5)
               CALL PGPT(1,mjdavg(j),log10(uflux_lc(j)),45)
               call PGPT(1,mjdavg(j),log10(uflux_lc(j))-0.07,31)
               CALL PGERRX(1,mjdst_lc(j),mjded_lc(j),log10(uflux_lc(j)),1.0)
            else if ((lflux_lc(j) .eq. 0.) .and. (uflux_lc(j) .ne. 0.) ) then
               call pgsch(1.5)
               CALL PGPT(1,mjdavg(j),log10(uflux_lc(j)),45)
               call PGPT(1,mjdavg(j),log10(uflux_lc(j))-0.07,31)
               CALL PGERRX(1,mjdst_lc(j),mjded_lc(j),log10(uflux_lc(j)),1.0)
            else if (flux_lc(j) .lt. 0.) then
               call pgsch(1.2)
               CALL PGPT(1,mjdavg(j),log10(-flux_lc(j)),13)
               CALL PGERRY(1,mjdavg(j),log10(uflux_lc(j)),log10(lflux_lc(j)),1.0)
               CALL PGERRX(1,mjdst_lc(j),mjded_lc(j),log10(-flux_lc(j)),1.0)
            else
               call pgsch(1.2)
               CALL PGPT(1,mjdavg(j),log10(flux_lc(j)),-17)
               CALL PGERRY(1,mjdavg(j),log10(uflux_lc(j)),log10(lflux_lc(j)),1.0)
               CALL PGERRX(1,mjdst_lc(j),mjded_lc(j),log10(flux_lc(j)),1.0)
            endif
         endif
      enddo
      endif

c second panel for LC
      call pgsci(1)
      call pgsch(1.3)
      if (iilc-ixray .ne. 0) then
      CALL PGENV(mjdlow(2),mjdup(2),lclow(2),lcup(2),0,1)
      CALL PGLAB('MJD', 'Log \gnf\d\gn\u (erg/s/cm\u2\d)','1 TeV')
      call pgsci(8)
      do j=1,iilc
         if (j .gt. ixray) then
            if ((flux_lc(j) .eq. lflux_lc(j)) .and. (flux_lc(j) .eq. uflux_lc(j))) then
               call pgsch(1.5)
               CALL PGPT(1,mjdavg(j),log10(uflux_lc(j)),45)
               call PGPT(1,mjdavg(j),log10(uflux_lc(j))-0.07,31)
               CALL PGERRX(1,mjdst_lc(j),mjded_lc(j),log10(uflux_lc(j)),1.0)
            else if ((lflux_lc(j) .eq. 0.) .and. (uflux_lc(j) .ne. 0.) ) then
               call pgsch(1.5)
               CALL PGPT(1,mjdavg(j),log10(uflux_lc(j)),45)
               call PGPT(1,mjdavg(j),log10(uflux_lc(j))-0.07,31)
               CALL PGERRX(1,mjdst_lc(j),mjded_lc(j),log10(uflux_lc(j)),1.0)
            else if (flux_lc(j) .lt. 0.) then
               call pgsch(1.2)
               CALL PGPT(1,mjdavg(j),log10(-flux_lc(j)),13)
               CALL PGERRY(1,mjdavg(j),log10(uflux_lc(j)),log10(lflux_lc(j)),1.0)
               CALL PGERRX(1,mjdst_lc(j),mjded_lc(j),log10(-flux_lc(j)),1.0)
            else
               call pgsch(1.2)
               CALL PGPT(1,mjdavg(j),log10(flux_lc(j)),-17)
               CALL PGERRY(1,mjdavg(j),log10(uflux_lc(j)),log10(lflux_lc(j)),1.0)
               CALL PGERRX(1,mjdst_lc(j),mjded_lc(j),log10(flux_lc(j)),1.0)
            endif
         endif
      enddo
      endif
      call pgsci(1)
      CALL PGEND

      close (10)
      END

