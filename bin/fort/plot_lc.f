      PROGRAM plot_lc

c This program plot the SED for candidate

      implicit none
      integer*4 ier,pgbeg,length,ns,j,rah, ram, id, dm,in,im,pttest,iflcuv,iir,IERGAM
      integer*4 i,sfound,npt(15000),rtype,stype(1000),ivhe,ixray,iilc,ialma,lctype(15000)
      real*4 frequency(15000,1000),flux(15000,1000),uflux(15000,1000),lflux(15000,1000)
      real*4 rasec,decsec,testflux,mjdstart(15000,1000),mjdend(15000,1000),mjdavg(15000)
      real*4 mjdlow,mjdup,lcup(5),lclow(5),flux_lc(15000),uflux_lc(15000),lflux_lc(15000)
      real*4 mjdst_lc(15000),mjded_lc(15000),freq_lc(15000),fq1tev,sloperat,mjdstgam,mjdedgam
      real*8 rra,rdec,ra(1000),dec(1000)
      character*160 string
      character*100 title
      character*200 input_file,output_file,output_file2,refs(15000,1000),output_file3
      character*14 stringin
      character*10 spectype(15000,1000)
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
      in=im
      im=index(string(in+1:length),' ')+in
      output_file2=string(in+1:im-1)
      in=im
      im=index(string(in+1:length),' ')+in
      output_file3=string(in+1:im-1)
      number=string(im+1:length)
      !write(*,*) number
      if (number == 'finish') then
         Stop '!!!!Exit the source exploring routine!!!!'
      else if (number == 'sed') then
         ns=99
      else
         read(number,*) ns
      endif

      !write(*,*) output_file,output_file2,output_file3
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
      read(10,*) string
      read(10,*) string
      read(10,*) string
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
      ialma=0
      iflcuv=0
      iir=0
      iilc=0
      do j=1,npt(i)
         if (((spectype(j,i) == 'OUSXB') .or. (spectype(j,i) == 'OUSPEC') .or. (spectype(j,i) == 'OUSXG'))
     &              .and. (frequency(j,i) .eq. 2.418E17)) then
            ixray=ixray+1
            iilc=iilc+1
            if (ixray .eq. 1) then
               lcup(3)=max(uflux(j,i),abs(flux(j,i)))
               lclow(3)=min(lflux(j,i),abs(flux(j,i)))
               if (iilc .eq. 1) mjdlow=mjdstart(j,i)
               if (iilc .eq. 1) mjdup=mjdend(j,i)
            endif
            if (spectype(j,i) == 'OUSPEC') then
               lctype(iilc)=10
            else
               lctype(iilc)=20
            endif
            testflux=max(uflux(j,i),abs(flux(j,i))) ! range for 1kev
            if ((testflux .gt. lcup(3)) .and. (testflux .gt. 0.d0)) lcup(3)=testflux
            testflux=min(lflux(j,i),abs(flux(j,i)))
            if ((testflux .lt. lclow(3)) .and. (testflux .gt. 0.d0)) lclow(3)=testflux
            if (mjdstart(j,i) .lt. mjdlow) mjdlow=mjdstart(j,i)
            if (mjdend(j,i) .gt. mjdup) mjdup=mjdend(j,i)
            freq_lc(iilc)=frequency(j,i)
            flux_lc(iilc)=flux(j,i)
            lflux_lc(iilc)=lflux(j,i)
            uflux_lc(iilc)=uflux(j,i)
            mjdst_lc(iilc)=mjdstart(j,i)
            mjded_lc(iilc)=mjdend(j,i)
         else if (spectype(j,i) == 'ALMA') then
            ialma=ialma+1
            iilc=iilc+1
            if (ialma .eq. 1) then
               lcup(1)=max(uflux(j,i),abs(flux(j,i)))
               lclow(1)=min(lflux(j,i),abs(flux(j,i)))
               if (iilc .eq. 1) mjdlow=mjdstart(j,i)
               if (iilc .eq. 1) mjdup=mjdend(j,i)
            endif
            lctype(iilc)=30
            testflux=max(uflux(j,i),abs(flux(j,i))) ! range for 1kev
            if ((testflux .gt. lcup(1)) .and. (testflux .gt. 0.d0)) lcup(1)=testflux
            testflux=min(lflux(j,i),abs(flux(j,i)))
            if ((testflux .lt. lclow(1)) .and. (testflux .gt. 0.d0)) lclow(1)=testflux
            if (mjdstart(j,i) .lt. mjdlow) mjdlow=mjdstart(j,i)
            if (mjdend(j,i) .gt. mjdup) mjdup=mjdend(j,i)
            freq_lc(iilc)=frequency(j,i)
            flux_lc(iilc)=flux(j,i)
            lflux_lc(iilc)=lflux(j,i)
            uflux_lc(iilc)=uflux(j,i)
            mjdst_lc(iilc)=mjdstart(j,i)
            mjded_lc(iilc)=mjdend(j,i)
          else if (((spectype(j,i) == 'WISELC') .or. (spectype(j,i) == 'NEOWISE'))
     &      .and. (frequency(j,i) .eq. 6.517e13)) then
            iir=iir+1
            iilc=iilc+1
            if (iir .eq. 1) then
               lcup(2)=max(uflux(j,i),abs(flux(j,i)))
               lclow(2)=min(lflux(j,i),abs(flux(j,i)))
               if (iilc .eq. 1) mjdlow=mjdstart(j,i)
               if (iilc .eq. 1) mjdup=mjdend(j,i)
               if (lcup(2) .eq. 0.) lcup(2)=5.e-10
               if (lclow(2) .eq. 0.) lclow(2)=5.e-12
            endif
            lctype(iilc)=40
            testflux=max(uflux(j,i),abs(flux(j,i))) ! range for 1kev
            if ((testflux .gt. lcup(2)) .and. (testflux .gt. 0.d0)) lcup(2)=testflux
            testflux=min(lflux(j,i),abs(flux(j,i)))
            if ((testflux .lt. lclow(2)) .and. (testflux .gt. 0.d0)) lclow(2)=testflux
            if (mjdstart(j,i) .lt. mjdlow) mjdlow=mjdstart(j,i)
            if (mjdend(j,i) .gt. mjdup) mjdup=mjdend(j,i)
            freq_lc(iilc)=frequency(j,i)
            flux_lc(iilc)=flux(j,i)
            lflux_lc(iilc)=lflux(j,i)
            uflux_lc(iilc)=uflux(j,i)
            mjdst_lc(iilc)=mjdstart(j,i)
            mjded_lc(iilc)=mjdend(j,i)
          else if (spectype(j,i) == 'FMonLC') then
            iflcuv=iflcuv+1
            iilc=iilc+1
            if (iflcuv .eq. 1) then
               lcup(4)=max(uflux(j,i),abs(flux(j,i)))
               lclow(4)=min(lflux(j,i),abs(flux(j,i)))
               if (iilc .eq. 1) mjdlow=mjdstart(j,i)
               if (iilc .eq. 1) mjdup=mjdend(j,i)
               if (lcup(4) .eq. 0.) lcup(4)=5.e-10
               if (lclow(4) .eq. 0.) lclow(4)=5.e-12
               mjdstgam=mjdstart(j,i)
               mjdedgam=mjdend(j,i)
            endif
c            write(*,*) lcup(3),lclow(3)
            if (frequency(j,i) .eq. 2.418e23) then
               lctype(iilc)=50
            else if (frequency(j,i) .eq. 2.418e24) then
               lctype(iilc)=52
            else
               lctype(iilc)=51
            endif
            testflux=max(uflux(j,i),abs(flux(j,i))) ! range for 1kev
            if ((testflux .gt. lcup(4)) .and. (testflux .gt. 0.d0)) lcup(4)=testflux
            testflux=min(lflux(j,i),abs(flux(j,i)))
            if ((testflux .lt. lclow(4)) .and. (testflux .gt. 0.d0)) lclow(4)=testflux
            if (mjdstart(j,i) .lt. mjdlow) mjdlow=mjdstart(j,i)
            if (mjdend(j,i) .gt. mjdup) mjdup=mjdend(j,i)
            if (mjdstart(j,i) .lt. mjdstgam) mjdstgam=mjdstart(j,i)
            if (mjdend(j,i) .gt. mjdedgam) mjdedgam=mjdend(j,i)
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
               lcup(5)=max(uflux(j,i),abs(flux(j,i)))
               lclow(5)=min(lflux(j,i),abs(flux(j,i)))
               if (iilc .eq. 0) mjdup=mjdend(j,i)
               if (iilc .eq. 0) mjdlow=mjdstart(j,i)
            endif
            if (frequency(j,i) .gt. fq1tev) then
               if (ivhe .eq. 1) then
                  iilc=iilc+1
                  sloperat=(log10(flux(j+1,i)/frequency(j+1,i))-(flux(j,i)/frequency(j,i)))/
     &                       (log10(frequency(j+1,i))-log10(frequency(j,i)))
                  flux_lc(iilc)=(flux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                  flux_lc(iilc)=flux_lc(iilc)*fq1tev
                  uflux_lc(iilc)=(uflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                  uflux_lc(iilc)=uflux_lc(iilc)*fq1tev
                  lflux_lc(iilc)=(lflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                  lflux_lc(iilc)=lflux_lc(iilc)*fq1tev
                  mjdst_lc(iilc)=mjdstart(j,i)
                  mjded_lc(iilc)=mjdend(j,i)
                  freq_lc(iilc)=fq1tev
                  lctype(iilc)=60
               else
                  if (mjdstart(j,i) .ne. mjdstart(j-1,i)) then
                     if (iilc-ixray-ialma-iflcuv .ne. 0) then !exclude the same one
                        if (mjdstart(j,i) .eq. mjdst_lc(iilc)) then
                           goto 600
                        endif
                     endif
                     iilc=iilc+1
                     sloperat=(log10(flux(j+1,i)/frequency(j+1,i))-log10(flux(j,i)/frequency(j,i)))/
     &                       (log10(frequency(j+1,i))-log10(frequency(j,i)))
                     flux_lc(iilc)=(flux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     flux_lc(iilc)=flux_lc(iilc)*fq1tev
                     uflux_lc(iilc)=(uflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     uflux_lc(iilc)=uflux_lc(iilc)*fq1tev
                     lflux_lc(iilc)=(lflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     lflux_lc(iilc)=lflux_lc(iilc)*fq1tev
                     mjdst_lc(iilc)=mjdstart(j,i)
                     mjded_lc(iilc)=mjdend(j,i)
                     freq_lc(iilc)=fq1tev
                     lctype(iilc)=60
                  else
                     if (iilc-ixray-ialma-iflcuv .ne. 0) then
                        if (mjdstart(j,i) .eq. mjdst_lc(iilc)) then
                           goto 600
                        endif
                     endif
                     iilc=iilc+1
                     sloperat=(log10(flux(j,i)/frequency(j,i))-log10(flux(j-1,i)/frequency(j-1,i)))
     &                              /(log10(frequency(j,i))-log10(frequency(j-1,i)))
                     flux_lc(iilc)=(flux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     flux_lc(iilc)=flux_lc(iilc)*fq1tev
                     uflux_lc(iilc)=(uflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     uflux_lc(iilc)=uflux_lc(iilc)*fq1tev
                     lflux_lc(iilc)=(lflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     lflux_lc(iilc)=lflux_lc(iilc)*fq1tev
                     mjdst_lc(iilc)=mjdstart(j,i)
                     mjded_lc(iilc)=mjdend(j,i)
                     freq_lc(iilc)=fq1tev
                     lctype(iilc)=60
                  endif
               endif
c               write(*,*) ixray,iilc,flux_lc(iilc),uflux_lc(iilc),mjdst_lc(iilc),mjded_lc(iilc)
               testflux=max(uflux_lc(iilc),flux_lc(iilc))
               if ((testflux .gt. lcup(5)) .and. (testflux .gt. 0.d0)) lcup(5)=testflux
               testflux=min(lflux_lc(iilc),flux_lc(iilc))
               if ((testflux .lt. lclow(5)) .and. (testflux .gt. 0.d0)) lclow(5)=testflux
               if (mjdst_lc(iilc) .lt. mjdlow) mjdlow=mjdst_lc(iilc)
               if (mjded_lc(iilc) .gt. mjdup) mjdup=mjdst_lc(iilc)
            else
               if (j .ne. npt(i)) then
                  if ( mjdstart(j,i) .ne. mjdstart(j+1,i) ) then
                     iilc=iilc+1
                     sloperat=(log10(flux(j,i)/frequency(j,i))-log10(flux(j-1,i)/frequency(j-1,i)))
     &                              /(log10(frequency(j,i))-log10(frequency(j-1,i)))
                     flux_lc(iilc)=(flux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     flux_lc(iilc)=flux_lc(iilc)*fq1tev
                     uflux_lc(iilc)=(uflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     uflux_lc(iilc)=uflux_lc(iilc)*fq1tev
                     lflux_lc(iilc)=(lflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                     lflux_lc(iilc)=lflux_lc(iilc)*fq1tev
                     mjdst_lc(iilc)=mjdstart(j,i)
                     mjded_lc(iilc)=mjdend(j,i)
                     freq_lc(iilc)=fq1tev
                     lctype(iilc)=60
                  endif
               else
                  iilc=iilc+1
                  sloperat=(log10(flux(j,i)/frequency(j,i))-log10(flux(j-1,i)/frequency(j-1,i)))
     &                              /(log10(frequency(j,i))-log10(frequency(j-1,i)))
                  flux_lc(iilc)=(flux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                  flux_lc(iilc)=flux_lc(iilc)*fq1tev
                  uflux_lc(iilc)=(uflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                  uflux_lc(iilc)=uflux_lc(iilc)*fq1tev
                  lflux_lc(iilc)=(lflux(j,i)/frequency(j,i))*(fq1tev/frequency(j,i))**sloperat
                  lflux_lc(iilc)=lflux_lc(iilc)*fq1tev
                  mjdst_lc(iilc)=mjdstart(j,i)
                  mjded_lc(iilc)=mjdend(j,i)
                  freq_lc(iilc)=fq1tev
                  lctype(iilc)=60
               endif
               testflux=max(uflux_lc(iilc),flux_lc(iilc))
               if ((testflux .gt. lcup(5)) .and. (testflux .gt. 0.d0)) lcup(5)=testflux
               testflux=min(lflux_lc(iilc),flux_lc(iilc))
               if ((testflux .lt. lclow(5)) .and. (testflux .gt. 0.d0)) lclow(5)=testflux
               if (mjdst_lc(iilc) .lt. mjdlow) mjdlow=mjdst_lc(iilc)
               if (mjded_lc(iilc) .gt. mjdup) mjdup=mjdst_lc(iilc)
            endif
  600       continue
         endif
      enddo

      write(*,*) iilc,ialma,iir,ixray,iflcuv
      !write(*,*) freq_lc(1:iilc),flux_lc(1:iilc)
      open(11,file=output_file2,status='unknown',iostat=ier)

      if (iilc .eq. 0) then
         write(*,*) 'NO Light Curve Plot'
         stop
      endif
      write(*,*) 'number of pts read for light curve:',iilc!,spectype(npt(sfound),sfound)
      if (iilc .eq. 0) then
         mjdlow=0
         mjdup=0
      endif
      if (ialma .eq. 0) then
         lcup(1)=0
         lclow(1)=0
      endif
      if (iir .eq. 0) then
         lcup(2)=0
         lclow(2)=0
      endif
      if (ixray .eq. 0) then
         lcup(3)=0
         lclow(3)=0
      endif
      if (iflcuv .eq. 0) then
         lcup(4)=0
         lclow(4)=0
      endif
      if (iilc-ixray-ialma-iflcuv .eq. 0) then
         lcup(5)=0
         lclow(5)=0
      endif
c      write(*,*) 'Fermi range',lcup(4),lclow(4)

      lcup=alog10(lcup)+0.5
      lclow=alog10(lclow)-0.5
      mjdlow=mjdlow-20.
      mjdup=mjdup+20.
      mjdstgam=mjdstgam-20.
      mjdedgam=mjdedgam+20.
      mjdavg=(mjdst_lc+mjded_lc)/2.

      write(11,'(i4,2x,a,2(2x,f9.5),2x,i2)') i,"matched source",ra(i),dec(i),stype(i)

!write(*,*) i
c      IER = PGBEG(0,"/xwindow",1,1)
c      IER = PGBEG(0,"/xs",1,1)
      IER = PGBEG(0,output_file,1,4)
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
      call pgsch(2.2)
c      write(*,*) mjdavg(1:8,1)

      if (ialma .ne.0) then
      CALL PGENV(mjdlow,mjdup,lclow(1),lcup(1),0,1)
      if ((iilc-ixray-ialma-iflcuv .eq. 0) .and. (ixray .eq. 0) .and. (iflcuv .eq. 0) .and. (iir .eq. 0)) then
         xtitle='MJD'
      else
         xtitle=''
      endif
      CALL PGLAB(xtitle, 'Log \gnf\d\gn\u (erg/s/cm\u2\d)','100 GHz')
      call pgsci(3)
      do j=1,iilc
         if (freq_lc(j) .lt. 5.e12 ) then
            write(11,'(4(es10.3,2x),2(f10.4,2x),i2)') freq_lc(j),flux_lc(j),
     *        uflux_lc(j),lflux_lc(j),mjdst_lc(j),mjded_lc(j),lctype(j)
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

c panel for IR
      call pgsci(1)
      call pgsch(2.2)
      if (iir .ne. 0) then
      CALL PGENV(mjdlow,mjdup,lclow(2),lcup(2),0,1)
      if ((iilc-ixray-ialma-iflcuv .eq. 0) .and. (iflcuv .eq. 0)
     &    .and.  (ixray .eq. 0)  ) then
         xtitle='MJD'
      else
         xtitle=''
      endif
      CALL PGLAB(xtitle, 'Log \gnf\d\gn\u (erg/s/cm\u2\d)','W2')
      do j=1,iilc
         if ((freq_lc(j) .gt. 5.e12) .and. (freq_lc(j) .lt. 1.e15)) then
            write(11,'(4(es10.3,2x),2(f10.4,2x),i2)') freq_lc(j),flux_lc(j),
     *        uflux_lc(j),lflux_lc(j),mjdst_lc(j),mjded_lc(j),lctype(j)
            call pgsci(8)
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

c panel for X-ray
      call pgsci(1)
      call pgsch(2.2)
      if (ixray .ne. 0) then
      CALL PGENV(mjdlow,mjdup,lclow(3),lcup(3),0,1)
      if ((iilc-ixray-ialma-iflcuv .eq. 0) .and. (iflcuv .eq. 0)) then
         xtitle='MJD'
      else
         xtitle=''
      endif
      CALL PGLAB(xtitle, 'Log \gnf\d\gn\u (erg/s/cm\u2\d)','1 keV')
      do j=1,iilc
         if ((freq_lc(j) .gt. 1.e14) .and. (freq_lc(j) .lt. 1.e20)) then
            write(11,'(4(es10.3,2x),2(f10.4,2x),i2)') freq_lc(j),flux_lc(j),
     *        uflux_lc(j),lflux_lc(j),mjdst_lc(j),mjded_lc(j),lctype(j)
            if (lctype(j) .eq. 10) then
               call pgsci(5) !light blue, for OUSPEC
            else
               call pgsci(4) !Dark blue for OUSXB
            endif
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

c   Fermi light curve
      call pgsci(1)
      call pgsch(2.2)
      if (iflcuv .ne. 0) then
      CALL PGENV(mjdlow,mjdup,lclow(4),lcup(4),0,1)
      if (iilc-ixray-ialma-iflcuv .eq. 0) then
        xtitle='MJD'
      else
        xtitle=''
      endif
      CALL PGLAB(xtitle, 'Log \gnf\d\gn\u (erg/s/cm\u2\d)','1 GeV')
      do j=1,iilc
         if (freq_lc(j) .eq. 2.418e23) then
            write(11,'(4(es10.3,2x),2(f10.4,2x),i2)') freq_lc(j),flux_lc(j),
     *        uflux_lc(j),lflux_lc(j),mjdst_lc(j),mjded_lc(j),lctype(j)
            call pgsci(12) !Purple for Fermi
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

c TeV panel for LC
      call pgsci(1)
      call pgsch(2.2)
      if (ivhe .ne. 0) then
      CALL PGENV(mjdlow,mjdup,lclow(5),lcup(5),0,1)
      CALL PGLAB('MJD', 'Log \gnf\d\gn\u (erg/s/cm\u2\d)','1 TeV')
      call pgsci(2)
      do j=1,iilc
         if (freq_lc(j) .gt. 1.e24) then
            write(11,'(4(es10.3,2x),2(f10.4,2x),i2)') freq_lc(j),flux_lc(j),
     *        uflux_lc(j),lflux_lc(j),mjdst_lc(j),mjded_lc(j),lctype(j)
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

      if (iflcuv .ne. 0) then
      IERGAM = PGBEG(1,output_file3,1,3)
      call pgslw(4)
      CALL PGSCRN(0, 'White', IERGAM)
      CALL PGSCRN(1, 'Black', IERGAM)
      call pgsch(1.7)
      call pgsci(1)
      CALL PGENV(mjdstgam,mjdedgam,lclow(4),lcup(4),0,1)
      CALL PGLAB('MJD','Log \gnf\d\gn\u (erg/s/cm\u2\d)','500 MeV')
      do j=1,iilc
         if (lctype(j) .eq. 51) then
            write(11,'(4(es10.3,2x),2(f10.4,2x),i2)') freq_lc(j),flux_lc(j),
     *        uflux_lc(j),lflux_lc(j),mjdst_lc(j),mjded_lc(j),lctype(j)
            call pgsci(12) !Purple for Fermi
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

      call pgsch(1.7)
      call pgsci(1)
      CALL PGENV(mjdstgam,mjdedgam,lclow(4),lcup(4),0,1)
      CALL PGLAB('','Log \gnf\d\gn\u (erg/s/cm\u2\d)','1 GeV')
      do j=1,iilc
         if (lctype(j) .eq. 50) then
c            write(11,'(4(es10.3,2x),2(f10.4,2x),i2)') freq_lc(j),flux_lc(j),
c     *        uflux_lc(j),lflux_lc(j),mjdst_lc(j),mjded_lc(j),lctype(j)
            call pgsci(12) !Purple for Fermi
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

      call pgsch(1.7)
      call pgsci(1)
      CALL PGENV(mjdstgam,mjdedgam,lclow(4),lcup(4),0,1)
      CALL PGLAB('','Log \gnf\d\gn\u (erg/s/cm\u2\d)','10 GeV')
      do j=1,iilc
         if (lctype(j) .eq. 52) then
            write(11,'(4(es10.3,2x),2(f10.4,2x),i2)') freq_lc(j),flux_lc(j),
     *        uflux_lc(j),lflux_lc(j),mjdst_lc(j),mjded_lc(j),lctype(j)
            call pgsci(12) !Purple for Fermi
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
      call pgsci(1)
      CALL PGEND
      endif
 
      close (10)
      END

