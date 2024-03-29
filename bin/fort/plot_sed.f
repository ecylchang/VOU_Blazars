      PROGRAM plot_sed

c This program plot the SED for candidate

      implicit none
      integer*4 ier,pgbeg,length,ns,j,rah, ram, id, dm,in,im,is
      integer*4 i,sfound,rtype,iconfig,iarr,arrsize(30),lenact
      real*4 sedup,sedlow,rasec,decsec,testflux
      real*8 rra,rdec
      character*400 string
      character*100 title
      character*200 input_file,output_file,webprograms,array_size
      character*14 stringin
      character*6 number
      character*1 sign
      logical ok,there

      integer*4,dimension(:),allocatable :: npt,stype
      real*4,dimension(:,:),allocatable :: frequency,flux,uflux,lflux,mjdend,mjdstart
      character*200,dimension(:,:),allocatable :: refs
      character*15,dimension(:,:),allocatable :: spectype
      character*2,dimension(:,:),allocatable :: flag
      real*8,dimension(:),allocatable :: ra,dec

      ok = .true.

      call rdforn(string,length)
      call rmvlbk(string)
      in=index(string(1:length),' ')
      input_file=string(1:in-1)
      im=index(string(in+1:length),' ')+in
      output_file=string(in+1:im-1)
      in=im
      im=index(string(in+1:length),' ')+in
      webprograms=string(in+1:im-1)
      number=string(im+1:length)
      !write(*,*) number
      if (number == 'finish') then
         Stop '!!!!Exit the source exploring routine!!!!'
      else if ((number == 'sed') .or. (number == 'lcurve')) then
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
      
      iarr=0
      iconfig=0
      arrsize(1:30)=0
      array_size=webprograms(1:lenact(webprograms))//'/array_size.cf'
      open(18,file=array_size,status='old',iostat=ier)
      do while(ok)
      read(18,'(a)',end=700) string
      if (string(1:3) == '---' ) then
         iconfig=iconfig+1
      else
         if (iconfig .eq. 1) then
            iarr=iarr+1
            is=index(string(1:len(string)),':')
            !write(*,*) iconfig,iarr,is,lenact(string)
            read(string(is+1:lenact(string)),*) arrsize(iarr)
         endif
      endif
      enddo
700   continue
      close(18)
      !write(*,*) arrsize

      allocate(stype(arrsize(5)),npt(arrsize(5)))
      allocate(ra(arrsize(1)),dec(arrsize(1)))
      allocate(frequency(arrsize(5),arrsize(1)),flux(arrsize(5),arrsize(1)),uflux(arrsize(5),arrsize(1)),lflux(arrsize(5),arrsize(1)),mjdend(arrsize(5),arrsize(1)),mjdstart(arrsize(5),arrsize(1)))
      allocate(refs(arrsize(5),arrsize(1)),spectype(arrsize(5),arrsize(1)),flag(arrsize(5),arrsize(1)))

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
     &     uflux(npt(sfound),sfound),lflux(npt(sfound),sfound),mjdstart(npt(sfound),sfound),
     &     mjdend(npt(sfound),sfound),flag(npt(sfound),sfound),spectype(npt(sfound),sfound)!,
c     &      refs(npt(sfound),sfound)
c         write(*,*) spectype(npt(sfound),sfound)
      enddo
99    continue
      npt(sfound)=npt(sfound)-1
      write(*,*) 'number of pts read',npt(sfound)!,spectype(npt(sfound),sfound)
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
      IER = PGBEG(0,output_file,1,1)
      call pgslw(4)
      if ( ns .eq. 99 ) THEN
         write(title,'(a,i2.2,1x,i2.2,1x,f4.1,a,a,i2.2,1x,i2.2,1x,f4.1)')
     &     'Source position:  ',rah,ram,rasec,' , ',sign,id,dm,decsec
      else
         write(title,'(a,i2,a,i2.2,1x,i2.2,1x,f4.1,a,a,i2.2,1x,i2.2,1x,f4.1)')
     &     'Source',i,' position:  ',rah,ram,rasec,' , ',sign,id,dm,decsec
      endif
      CALL PGSCRN(0, 'White', IER)
      CALL PGSCRN(1, 'Black', IER)
      call pgsch(1.3)
      sedup=max(uflux(1,i),abs(flux(1,i)))
      sedlow=5.e-16

      do j=1,npt(i)
         testflux=max(uflux(j,i),abs(flux(j,i)))
         if ((frequency(j,i) .lt. 1.e10) .and. (testflux .gt. 0.d0)) sedlow=testflux
      enddo
      do j=1,npt(i)
         testflux=max(uflux(j,i),abs(flux(j,i)))
         if ((testflux .gt. sedup) .and. (testflux .gt. 0.d0)) sedup=testflux
         if ((testflux .lt. sedlow) .and. (testflux .gt. 0.d0)) sedlow=testflux
      enddo
      write(*,*) 'SED range upper limit',sedup,'SED range lower limit',sedlow
         sedup=alog10(sedup)+0.5
         sedlow=alog10(sedlow)-0.5
      CALL PGENV(7.5,27.,sedlow,sedup,0,1)
      CALL PGLAB('Log \gn (Hz)', 'Log \gnf\d\gn\u (erg/s/cm\u2\d)',title)
c      if (stype(i) .eq. 1) call pgsci(8)
c      if (stype(i) .eq. 2) call pgsci(5)
c      if (stype(i) .eq. 3) call pgsci(4)
c      if (stype(i) .eq. 4) call pgsci(3)
c      if (stype(i) .eq. 5) call pgsci(1)
c      if (stype(i) .eq. 0) call pgsci(12)
      do j=1,npt(i)
         call pgsci(1)
c         if ((spectype(j,i) /= 'FMonLC') .and. (spectype(j,i) /= 'WISELC')
c     &     .and. (spectype(j,i) /= 'FTAptLC') .and. (spectype(j,i) /= 'NEOWISE')) then
         if (spectype(j,i) == 'DEBL') then
            call pgsch(1.2)
            call pgsci(8)
            CALL PGPT(1,log10(frequency(j,i)),log10(flux(j,i)),11)
            CALL PGERRY(1,log10(frequency(j,i)),log10(uflux(j,i)),log10(lflux(j,i)),1.0)
            goto 200
         endif
         if (flag(j,i) == 'UL' ) then
            call pgsch(1.5)
            CALL PGPT(1,log10(frequency(j,i)),log10(uflux(j,i)),45)
            call PGPT(1,log10(frequency(j,i)),log10(uflux(j,i))-0.07,31)
c         else if ((lflux(j,i) .eq. 0.) .and. (uflux(j,i) .ne. 0.) ) then
c            call pgsch(1.5)
c            CALL PGPT(1,log10(frequency(j,i)),log10(uflux(j,i)),45)
c            call PGPT(1,log10(frequency(j,i)),log10(uflux(j,i))-0.07,31)
         else if (flux(j,i) .lt. 0.) then
            call pgsch(1.2)
            CALL PGPT(1,log10(frequency(j,i)),log10(-flux(j,i)),13)
            CALL PGERRY(1,log10(frequency(j,i)),log10(uflux(j,i)),log10(lflux(j,i)),1.0)
         else
            if ((spectype(j,i) == 'XRTSPEC') .or. (spectype(j,i) == 'OUSPEC') .or. (spectype(j,i) == 'ALMA')
     &         .or. (spectype(j,i) == 'FMonLC') .or. (spectype(j,i) == 'FTAptLC')
     &         .or. (spectype(j,i) == 'WISELC') .or. (spectype(j,i) == 'NEOWISE')) then
               call pgsch(0.5)
               CALL PGPT(1,log10(frequency(j,i)),log10(flux(j,i)),3)
               CALL PGERRY(1,log10(frequency(j,i)),log10(uflux(j,i)),log10(lflux(j,i)),1.0)
            else
               call pgsch(0.8)
               CALL PGPT(1,log10(frequency(j,i)),log10(flux(j,i)),-17)
               CALL PGERRY(1,log10(frequency(j,i)),log10(uflux(j,i)),log10(lflux(j,i)),1.0)
            endif
         endif
c         endif
200      continue
      enddo
      call pgsci(1)
      CALL PGEND

      deallocate(npt,stype)
      deallocate(frequency,flux,uflux,lflux,mjdend,mjdstart)
      deallocate(refs,spectype,flag)
      deallocate(ra,dec)

      close (10)
      END

