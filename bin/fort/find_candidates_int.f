      PROGRAM find_candidates_int

      IMPLICIT none
      integer i,j,k,m,s,iuv,igam,i4p8,irr,ixx,in,im,lenact,length,lu_in,ier,icand,filen,is,ie,xpts(2000)
      integer ii1,ii2,ipass,isource,pass2(4000),rr_type(2000),xx_type(2000),xxot_type(10,2000)
      integer nrepxx(2000),backxx(2000,2000),ixxss,ixxrep,trackxx(2000),repnumberxx(2000),posindxx(2000)
      integer nreprr(2000),backrr(2000,2000),irrss,irrrep,trackrr(2000),repnumberrr(2000),posindrr(2000)
      integer nnuvx,nnruv,nnralpha,xxss_type(2000),rrss_type(2000),code
      real*8 ra,dec,dist,ra_uv(10000),dec_uv(10000),ra_gam(50),dec_gam(50),ra_4p8(500),dec_4p8(500)
      real*8 ra_rr(2000),dec_rr(2000),ra_xx(2000),dec_xx(2000),ra_center,dec_center
      real*8 ra_xxss(2000),dec_xxss(2000),ra_rrss(2000),dec_rrss(200)
      real*4 nh,flux,uflux,lflux,epos,freq,radius,flux2nufnu_4p8,fdens,nudens
      real*4 frequency_rr(2000),flux_rr(2000),FluxL_rr(2000),FluxU_rr(2000),poserr_rr(2000)
      real*4 frequency_xx(2000),flux_xx(2000),FluxL_xx(2000),FluxU_xx(2000),poserr_xx(2000)
      real*4 frequency_xxot(10,2000),flux_xxot(10,2000),FluxL_xxot(10,2000),FluxU_xxot(10,2000)
      real*4 mjdst_rr(2000),mjded_rr(2000),mjdst_xx(2000),mjded_xx(2000)
      real*4 frequency_4p8(500),flux_4p8(500),FluxL_4p8(500),FluxU_4p8(500)
      real*4 poserr_4p8(500),Ferr_4p8(500),poserr_uv(10000),poserr_xxss(2000),flux_rrss(2000)
      real*4 frequency_uv(10000,2),flux_uv(10000,2),FluxL_uv(10000,2),FluxU_uv(10000,2)
      real*4 uvmag(10000,2),uvmagerr(10000,2),slope_gam(50),specerr_gam(50)
      real*4 frequency_gam(50),flux_gam(50),FluxL_gam(50),FluxU_gam(50),poserr_gam(50),Ferr_gam(50)
      real*4 aruv,auvx,min_dist,alphar,mjdstart,mjdend
      character*80 input_file,input_file2,input_file3,output_file
      character*10 catalog,type_4p8(500)
      character*800 string
      LOGICAL there,ok,found
      ok = .TRUE.
      found = .FALSE.
c      min_dist_uv=8./3600.
c      min_dist_gam=5./60.
c      min_dist_4p8=30./3600.
      flux2nufnu_4p8=4.85E9*1.E-26

      CALL rdforn(string,length)
      CALL rmvlbk(string)
      in=index(string(1:length),' ')
      input_file=string(1:in-1)
      im=index(string(in+1:length),' ')+in
      input_file2=string(in+1:im-1)
      in=im
      im=index(string(in+1:length),' ')+in
      input_file3=string(in+1:im-1)
      output_file=string(im+1:length)
      in = index(input_file(1:lenact(input_file)),'.')
      IF (in == 0) input_file(lenact(input_file)+1:lenact(input_file)+4) = '.csv'
      INQUIRE (FILE=input_file,EXIST=there)
      IF (.NOT.there) THEN
         write (*,'('' file '',a,'' not found '')')
     &     input_file(1:lenact(input_file))
         STOP
      ENDIF

      lu_in = 10
      open(lu_in,file=input_file,status='old',iostat=ier)
      open(11,file=input_file2,status='old',iostat=ier)
      open(13,file=input_file3,status='old',iostat=ier)
      IF (ier.NE.0) THEN
         write (*,*) ' Error ',ier,' opening file ', input_file
      ENDIF

      icand=0
      ixx=0
      irr=0
      Do WHILE(ok)
         read(11,*,end=100) freq,flux,uflux,lflux,ra,dec,epos,mjdstart,mjdend,code
         if (code .lt. 0) THEN
            icand=icand+1
            irr=irr+1
            ra_rr(irr)=ra
            dec_rr(irr)=dec
            poserr_rr(irr)=epos
            frequency_rr(irr)=freq
            flux_rr(irr)=flux
            FluxU_rr(irr)=uflux
            FluxL_rr(irr)=lflux
            mjdst_rr(irr)=mjdstart
            mjded_rr(irr)=mjdend
            rr_type(irr)=abs(code)
         else if (code .gt. 50) THEN
            icand=icand+1
            ixx=ixx+1
            xpts(ixx)=0
            ra_xx(ixx)=ra
            dec_xx(ixx)=dec
            poserr_xx(ixx)=epos
            frequency_xx(ixx)=freq
            flux_xx(ixx)=flux
            FluxU_xx(ixx)=uflux
            FluxL_xx(ixx)=lflux
            mjdst_xx(ixx)=mjdstart
            mjded_xx(ixx)=mjdend
            xx_type(ixx)=code-50
         else
            xpts(ixx)=xpts(ixx)+1
            frequency_xxot(xpts(ixx),ixx)=freq
            flux_xxot(xpts(ixx),ixx)=flux
            FluxU_xxot(xpts(ixx),ixx)=uflux
            FluxL_xxot(xpts(ixx),ixx)=lflux
            xxot_type(xpts(ixx),ixx)=code
         endif
      ENDDO
100   continue
c      write(*,*) icand,irr,ixx

      isource=0
      do while(ok)
         read(13,*,end=98) ra,dec,code
         if ((code .gt. 0.) .or. (code .le. -50000.) .or. (code .eq. -9999)) then
            isource=isource+1
         endif
      enddo
98    continue
c      write(*,*) isource

      igam=0
      iuv=0
      i4p8=0
      READ(lu_in,'(a)') string
      is = index(string(1:len(string)),'=')
      ie = index(string(is+5:len(string)),' ') +is+4
      read(string(is+1:ie-1),*) ra_center
      is = index(string(ie+1:len(string)),'=') +ie
      ie = index(string(is+5:len(string)),' ') +is+4
      read(string(is+1:ie-1),*) dec_center
      is = index(string(ie+1:len(string)),'=') +ie
      ie = index(string(is+5:len(string)),' ') +is+4
      read(string(is+1:ie-1),*) radius
      READ(lu_in,'(a)') string !begin reading nh
      is = index(string(1:len(string)),'=')
      ie = index(string(is+5:len(string)),' ') +is+4
      read(string(is+1:ie-1),*) nh
      DO WHILE(ok)
         READ(lu_in,'(a)',end=99) string
         ie=index(string(1:len(string)),',')
         read(string(1:ie-1),*) filen
         is=ie
         ie=index(string(is+1:len(string)),',')+is
         read(string(is+1:ie-1),'(a)') catalog
         is=ie
         ie=index(string(is+1:len(string)),',')+is
         read(string(is+1:ie-1),*) ra
         is=ie
         ie=index(string(is+1:len(string)),',')+is
         read(string(is+1:ie-1),*) dec
         IF ( (catalog(1:3) == 'pmn') .OR. (catalog(1:3) == 'gb6') ) then
            i4p8=i4p8+1
            ra_4p8(i4p8)=ra
            dec_4p8(i4p8)=dec
            if (catalog(1:3) == 'gb6') then
               type_4p8(i4p8)='GB6'
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_4p8(i4p8)
            else
               type_4p8(i4p8)='PMN'
               poserr_4p8(i4p8)=sqrt((15.*15.)+100.)
            endif
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_4p8(i4p8)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_4p8(i4p8)
            FluxU_4p8(i4p8)=flux_4p8(i4p8)+Ferr_4p8(i4p8)
            FluxL_4p8(i4p8)=flux_4p8(i4p8)-Ferr_4p8(i4p8)
            flux_4p8(i4p8)=flux_4p8(i4p8)*flux2nufnu_4p8
            FluxU_4p8(i4p8)=FluxU_4p8(i4p8)*flux2nufnu_4p8
            FluxL_4p8(i4p8)=FluxL_4p8(i4p8)*flux2nufnu_4p8
            frequency_4p8(i4p8)=4.85e9
         else if (catalog(1:5) == 'galex') then
            iuv=iuv+1
            ra_uv(iuv)=ra
            dec_uv(iuv)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            read(string(is+1:ie-1),*) uvmag(iuv,1)
            if (uvmag(iuv,1) .le. -999.) uvmag(iuv,1)=0.
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            read(string(is+1:ie-1),*) uvmag(iuv,2)
            if (uvmag(iuv,2) .le. -999.) uvmag(iuv,2)=0.
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            read(string(is+1:ie-1),*) uvmagerr(iuv,1)
            if (uvmagerr(iuv,1) .le. -999.) uvmagerr(iuv,1)=0.
            is=ie
            ie=index(string(is+1:len(string)),' ')+is
            read(string(is+1:ie-1),*) uvmagerr(iuv,2)
            if (uvmagerr(iuv,2) .le. -999.) uvmagerr(iuv,2)=0.
            CALL  mag2flux (nh,uvmag(iuv,2),'fuv',flux_uv(iuv,2),frequency_uv(iuv,2))
         CALL  mag2flux (nh,uvmag(iuv,2)-uvmagerr(iuv,2),'fuv',FluxU_uv(iuv,2),frequency_uv(iuv,2))
         CALL  mag2flux (nh,uvmag(iuv,2)+uvmagerr(iuv,2),'fuv',FluxL_uv(iuv,2),frequency_uv(iuv,2))
            CALL  mag2flux (nh,uvmag(iuv,1),'nuv',flux_uv(iuv,1),frequency_uv(iuv,1))
         CALL  mag2flux (nh,uvmag(iuv,1)-uvmagerr(iuv,1),'nuv',FluxU_uv(iuv,1),frequency_uv(iuv,1))
         CALL  mag2flux (nh,uvmag(iuv,1)+uvmagerr(iuv,1),'nuv',FluxL_uv(iuv,1),frequency_uv(iuv,1))
            poserr_uv(iuv)=1.
         else if (catalog(1:8) == 'fermi8yr') then
            igam=igam+1
            ra_gam(igam)=ra
            dec_gam(igam)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_gam(igam)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) slope_gam(igam)
            is=ie
            ie=index(string(is+1:len(string)),' ')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) specerr_gam(igam)
            FluxU_gam(igam)=flux_gam(igam)+Ferr_gam(igam)
            FluxL_gam(igam)=flux_gam(igam)-Ferr_gam(igam)
            call fluxtofdens(slope_gam(igam),1.,100.,flux_gam(igam),1.,fdens,nudens)
            flux_gam(igam)=fdens
            frequency_gam(igam)=nudens
            call fluxtofdens(slope_gam(igam),1.,100.,FluxU_gam(igam),1.,fdens,nudens)
            FluxU_gam(igam)=fdens
            call fluxtofdens(slope_gam(igam),1.,100.,FluxL_gam(igam),1.,fdens,nudens)
            FluxL_gam(igam)=fdens
         endif
      ENDDO
99    continue
      close(lu_in)
      close(11)
c      write(*,*) i4p8,iuv,igam

      nrepxx(1:ixx)=0
      ixxss=0
      ixxrep=0
      do i=1,ixx
         if (i .ne. 1) then
            do j=1,i-1
               call DIST_SKY(abs(ra_xx(i)),dec_xx(i),abs(ra_xx(j)),dec_xx(j),dist)
               if (dist*3600. .lt. 18.) then
                  !write(*,*) 'nearby X-ray',ra_xx(i),dec_xx(i),xx_type(ixx)
                  ixxrep=ixxrep+1
c                  xrayrepeat(nrep(i),i)=j
                  trackxx(i)=trackxx(j)
                  nrepxx(trackxx(i))=nrepxx(trackxx(i))+1
                  repnumberxx(i)=nrepxx(trackxx(i))
c                  write(*,*) i,track(i),xx_type(i),poserr_xx(i),nrep(track(i))
                  goto 97
               endif
            enddo
         endif
         ixxss=ixxss+1
         trackxx(i)=ixxss
         nrepxx(trackxx(i))=1
         repnumberxx(i)=nrepxx(trackxx(i))
c         write(*,*) i,track(i),xx_type(i),poserr_xx(i)!,back(i,repnumber(i))
97    continue
      enddo
      if (ixxss+ixxrep .ne. ixx ) write(*,*) 'Warning! may have the wrong number.'

c xrt 1sxps chandra
c bmw xmm xmmsl
c rass wga ipc
c      write(*,*) '========================'
      do i=1,ixxss
         do j=1,ixx
            if (trackxx(j) .eq. i) then
               backxx(i,repnumberxx(j))=j
               if (repnumberxx(j) .eq. 1) then
                  poserr_xxss(trackxx(j))=poserr_xx(j)
                  ra_xxss(trackxx(j))=ra_xx(j)
                  dec_xxss(trackxx(j))=dec_xx(j)
                  posindxx(trackxx(j))=j
                  xxss_type(trackxx(j))=xx_type(j)
               else
c                  write(*,*) j,track(j),xx_type(j),poserr_xx(j),repnumber(j)
                  if ((xx_type(j) .eq. 55) .or. (xx_type(j) .eq. 59) .or. (xx_type(j) .eq. 58) .or.
     &                (xx_type(j) .eq. 61) .or. (xx_type(j) .eq. 64)) then
                     if ((xxss_type(trackxx(j)) .ne. 55) .or. (xxss_type(trackxx(j)) .ne. 59)
     &                .or. (xxss_type(trackxx(j)) .ne. 58) .or. (xxss_type(trackxx(j)) .eq. 61)
     &                .or. (xxss_type(trackxx(j)) .eq. 64)) then
                        ra_xxss(trackxx(j))=ra_xx(j)
                        dec_xxss(trackxx(j))=dec_xx(j)
                        posindxx(trackxx(j))=j
                     else
                        if (poserr_xx(j) .lt. poserr_xxss(trackxx(j)) ) then
                           ra_xxss(trackxx(j))=ra_xx(j)
                           dec_xxss(trackxx(j))=dec_xx(j)
                           posindxx(trackxx(j))=j
                        endif
                     endif
                  else if ((xx_type(j) .eq. 51) .or. (xx_type(j) .eq. 52) .or. (xx_type(j) .eq. 57)) then
                     if ((xxss_type(trackxx(j)) .eq. 53) .or. (xxss_type(trackxx(j)) .eq. 54)
     &                 .or. (xxss_type(trackxx(j)) .eq. 56) .or. (xxss_type(trackxx(j)) .eq. 60)
     &                 .or. (xxss_type(trackxx(j)) .eq. 62) .or. (xxss_type(trackxx(j)) .eq. 63)) then
                        ra_xxss(trackxx(j))=ra_xx(j)
                        dec_xxss(trackxx(j))=dec_xx(j)
                        posindxx(trackxx(j))=j
                     else if ((xxss_type(trackxx(j)) .eq. 51) .or. (xxss_type(trackxx(j)) .eq. 52)
     &                 .or. (xxss_type(trackxx(j)) .eq. 57)) then
                        if (poserr_xx(j) .lt. poserr_xxss(trackxx(j)) ) then
                           ra_xxss(trackxx(j))=ra_xx(j)
                           dec_xxss(trackxx(j))=dec_xx(j)
                           posindxx(trackxx(j))=j
                        endif
                     endif
                  else if ((xx_type(j) .eq. 53) .or. (xx_type(j) .eq. 54) .or. (xx_type(j) .eq. 62)
     &              .or. (xx_type(j) .eq. 56) .or. (xx_type(j) .eq. 60) .or. (xx_type(j) .eq. 63)) then
                     if ((xxss_type(trackxx(j)) .eq. 53) .or. (xxss_type(trackxx(j)) .eq. 54)
     &                 .or. (xxss_type(trackxx(j)) .eq. 56) .or. (xxss_type(trackxx(j)) .eq. 60)
     &                 .or. (xxss_type(trackxx(j)) .eq. 62) .or. (xxss_type(trackxx(j)) .eq. 63)) then
                        if (poserr_xx(j) .lt. poserr_xxss(trackxx(j)) ) then
                           ra_xxss(trackxx(j))=ra_xx(j)
                           dec_xxss(trackxx(j))=dec_xx(j)
                           posindxx(trackxx(j))=j
                        endif
                     endif
                  endif
               endif
            endif
         enddo
c         write(*,*) i,posind(i),ra_xxss(i),dec_xxss(i),nrep(i),back(i,1:nrep(i))
      enddo

      nreprr(1:irr)=0
      irrss=0
      irrrep=0
      do i=1,irr
         if (i .ne. 1) then
            do j=1,i-1
               call DIST_SKY(ra_rr(i),dec_rr(i),ra_rr(j),dec_rr(j),dist)
               if (dist*3600. .lt. 18.) then
                  !write(*,*) 'nearby X-ray',ra_xx(i),dec_xx(i),xx_type(ixx)
                  irrrep=irrrep+1
c                  xrayrepeat(nrep(i),i)=j
                  trackrr(i)=trackrr(j) !for repeated sources track back to the number
                  nreprr(trackrr(i))=nreprr(trackrr(i))+1 !repeated number
                  repnumberrr(i)=nreprr(trackrr(i))
c                  write(*,*) i,trackrr(i),rr_type(i),poserr_rr(i),nreprr(trackrr(i))
                  goto 96
               endif
            enddo
         endif
         irrss=irrss+1
         trackrr(i)=irrss
         nreprr(trackrr(i))=1
         repnumberrr(i)=nreprr(trackrr(i))
c         write(*,*) i,trackrr(i),rr_type(i),poserr_rr(i),nreprr(trackrr(i))!,back(i,repnumber(i))
96    continue
      enddo
      if (irrss+irrrep .ne. irr ) write(*,*) 'Warning! may have the wrong number.'

      do i=1,irrss
         do j=1,irr
            if (trackrr(j) .eq. i) then
               backrr(i,repnumberrr(j))=j !matrix for repeated source id
               if (repnumberrr(j) .eq. 1) then
                  flux_rrss(trackrr(j))=flux_rr(j)
                  ra_rrss(trackrr(j))=ra_rr(j)
                  dec_rrss(trackrr(j))=dec_rr(j)
                  posindrr(trackrr(j))=j
                  rrss_type(trackrr(j))=rr_type(j)
               else
c                  write(*,*) j,trackrr(j),rr_type(j),poserr_rr(j),repnumberrr(j)
                  if (rr_type(j) .eq. 1) then
                     if (rrss_type(trackrr(j)) .ne. 1) then
                        ra_rrss(trackrr(j))=ra_rr(j)
                        dec_rrss(trackrr(j))=dec_rr(j)
                        posindrr(trackrr(j))=j
                     else
                        if (flux_rr(j) .gt. flux_rrss(trackrr(j)) ) then
                           ra_rrss(trackrr(j))=ra_rr(j)
                           dec_rrss(trackrr(j))=dec_rr(j)
                           posindrr(trackrr(j))=j
                        endif
                     endif
                  else if (rr_type(j) .eq. 2) then
                     if (rrss_type(trackrr(j)) .eq. 3) then
                        ra_rrss(trackrr(j))=ra_rr(j)
                        dec_rrss(trackrr(j))=dec_rr(j)
                        posindrr(trackrr(j))=j
                     else if (rrss_type(trackrr(j)) .eq. 2) then
                        if (flux_rr(j) .gt. flux_rrss(trackrr(j)) ) then
                           ra_rrss(trackrr(j))=ra_rr(j)
                           dec_rrss(trackrr(j))=dec_rr(j)
                           posindrr(trackrr(j))=j
                        endif
                     endif
                  else if (rr_type(j) .eq. 3) then
                     if (rrss_type(trackrr(j)) .eq. 3) then
                        if (flux_rr(j) .gt. flux_rrss(trackrr(j)) ) then
                           ra_rrss(trackrr(j))=ra_rr(j)
                           dec_rrss(trackrr(j))=dec_rr(j)
                           posindrr(trackrr(j))=j
                        endif
                     endif
                  endif
               endif
            endif
         enddo
c         write(*,*) i,posindrr(i),ra_rrss(i),dec_rrss(i),nreprr(i),backrr(i,1:nreprr(i))
      enddo

      open(12,file=output_file,status='unknown',iostat=ier)
      ipass=0
      do i=1,irrss
         ii1=0
         ii2=0
         do m=1,nreprr(i)
            k=backrr(i,m)
            nnruv=0
            nnralpha=0
            do j=1,iuv
               call DIST_SKY(ra_rr(k),dec_rr(k),ra_uv(j),dec_uv(j),dist)
               min_dist=sqrt(poserr_uv(j)**2+poserr_rr(k)**2)
               if ((dist*3600. .lt. max(min_dist,2.)) .and. (flux_rr(k)/frequency_rr(k) .ge. 10*1.e-26)) then
                  if (flux_uv(j,1) .ne. 0.) then
                     aruv = 1.-log10(flux_rr(k)/flux_uv(j,1))/log10(frequency_rr(k)/frequency_uv(j,1))
                  else
                     aruv = 1.-log10(flux_rr(k)/flux_uv(j,2))/log10(frequency_rr(k)/frequency_uv(j,2))
                  endif
                  if (aruv .le. 0.75 ) then
                     ii1=ii1+1 !!!remove UV-r slope strange source 0.85
                     nnruv=nnruv+1
                     if (ii1 .eq. 1) then
                         write(*,'(a,i4,a)') "----------------------------------------"
                         write (*,'(a,i4)') "source nr.", isource+ii1
                     endif
                     if ((ii1 .eq. 1) .or. ((m .ne. 1) .and. (nnruv .eq. 1)))
     &                 write(*,'(f9.5,2x,f9.5,a,f9.3,a,2x,i2)') ra_rr(k),dec_rr(k)," radio source ",
     &                 flux_rr(k)/frequency_rr(k)/1.E-26," mJy",rr_type(k)
c                     write(*,'("GALEX : ",2(f6.3,2x),10x,f7.3," arcsec away")') uvmag(j,1),uvmag(j,2),dist*3600.
c                     write(*,'(6x,"radio-UV slope: ",f6.3)') aruv
                  endif
               endif
            enddo
            do j=1,i4p8
               call DIST_SKY(ra_rr(k),dec_rr(k),ra_4p8(j),dec_4p8(j),dist)
               min_dist=sqrt(poserr_4p8(j)**2+poserr_rr(k)**2)
               if ((dist*3600. .lt. min_dist ) .and. (flux_rr(k)/frequency_rr(k) .ge. 10*1.e-26)) then
                  alphar = 1.-log10(flux_rr(k)/flux_4p8(j))/log10(frequency_rr(k)/frequency_4p8(j))
                  if (alphar .le. 0.7 ) then
                     ii2=ii2+1 !!!remove radio extended sources
                     nnralpha=nnralpha+1
                     if ((ii2 .eq. 1) .and. (ii1 .eq. 0)) then
                         write(*,'(a,i4,a)') "----------------------------------------"
                         write (*,'(a,i4)') "source nr.", isource+ii1+ii2
                     endif
                     if (((ii2 .eq. 1) .and. (ii1 .eq. 0)) .or. ((m .ne. 1) .and. (nnruv .eq. 0)
     &                .and. (nnralpha .eq. 1))) write(*,'(f9.5,2x,f9.5,a,f9.3,a,2x,i2)')  ra_rr(k),
     &                dec_rr(k)," radio source ",flux_rr(k)/frequency_rr(k)/1.E-26," mJy",rr_type(k)
c                     write(*,'(a,f9.3," mJy",10x,f7.3," arcsec away")')
c     &                type_4p8(j),flux_4p8(j)/flux2nufnu_4p8,dist*3600.
c                     write(*,'("radio slope: ",f6.3)') alphar
                  endif
               endif
            enddo
         enddo
         if (ii1+ii2 .gt. 0.) then
            ipass=ipass+1
            CALL RXgraphic_code(flux_rr(posindrr(i))/frequency_rr(posindrr(i))/1.E-26,'R',code)
            write (12,'(f9.5,2x,f9.5,2x,i6)') ra_rr(posindrr(i)),dec_rr(posindrr(i)),int(code)
            pass2(ipass)=i
         endif
      enddo

      do i=1,ixxss
         ii1=0
         do m=1,nrepxx(i)
            k=backxx(i,m)
            nnuvx=0
            do j=1,iuv
               call DIST_SKY(abs(ra_xx(k)),dec_xx(k),ra_uv(j),dec_uv(j),dist)
               min_dist=sqrt(poserr_uv(j)**2+poserr_xx(k)**2)
               if ((dist*3600. .le. max(min_dist,2.) ) .and. (poserr_xx(k) .le. 20.)
     &              .and. (flux_xx(k) .ne. 0. ) ) then
                  if (flux_uv(j,1) .ne. 0.) then
                     auvx = 1.-log10(flux_uv(j,1)/flux_xx(k))/log10(frequency_uv(j,1)/frequency_xx(k))
                  else
                     auvx = 1.-log10(flux_uv(j,2)/flux_xx(k))/log10(frequency_uv(j,2)/frequency_xx(k))
                  endif
                  if (auvx .le. 1.4) then
                     ii1=ii1+1 !!!remove UV-x slope strange source 0.85
                     nnuvx=nnuvx+1
                     if (ii1 .eq. 1) then
                        write(*,'(a,i4,a)') "----------------------------------------"
                        write(*,'(a,i4)') "source nr.", ipass+ii1+isource
                     endif
                     if ((ii1 .eq. 1) .or. ((m .ne. 1) .and. (nnuvx .eq. 1)))
     $                  write(*,'(f9.5,2x,f9.5,a,es10.3,2x,i2)')  abs(ra_xx(k)),
     $                  dec_xx(k)," X-ray source with 1 keV flux ",flux_xx(k),xx_type(k)
c                     write(*,'("GALEX : ",2(f6.3,2x),10x,f7.3," arcsec away")') uvmag(j,1),uvmag(j,2),dist*3600.
c                     write(*,'(6x,"UV-X-ray slope: ",f6.3)') auvx
                  endif
               endif
            enddo
         enddo
         if (ii1 .gt. 0.) then
             ipass=ipass+1
             CALL RXgraphic_code(flux_xx(posindxx(i)),'X',code)
             write (12,'(f9.5,2x,f9.5,2x,i6)') abs(ra_xx(posindxx(i))),dec_xx(posindxx(i)),int(code)
             pass2(ipass)=i+irr
         endif
      enddo

      if ((ipass .eq. 0) .and. (isource .ne. 0)) then
         print *,achar(27),'[35;1m No Candidates Found in Intermediate Phase',achar(27),'[0m'
         stop
      endif
      if (ipass+isource .eq. 0) then
         print *,achar(27),'[31;1m No Candidates Found！',achar(27),'[0m'
         stop
      endif
c      write(*,*) pass2(1:ipass)

      do i=1,ipass
         k=pass2(i)
         if (k .le. irr) then
            do s=1,nreprr(k)
               m=backrr(k,s)
               if ((i+isource .ne. 1) .or. (s .ne. 1)) write(12,*) "===================="
               write(12,'(i4,2x,a,2(2x,f9.5),2x,a,2x,i2)') i+isource,"matched source",
     &            ra_rr(m),dec_rr(m),'source type',int(code/10000)
               write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,2(f10.4,2x),i2)') frequency_rr(m),flux_rr(m),
     &         FluxU_rr(m),FluxL_rr(m),ra_rr(m),dec_rr(m),poserr_rr(m),mjdst_rr(m),mjded_rr(m),rr_type(m)
            enddo
         else
            do s=1,nrepxx(k-irr)
               m=backxx(k-irr,s)
               if ((i+isource .ne. 1) .or. (s .ne. 1)) write(12,*) "===================="
               write(12,'(i4,2x,a,2(2x,f10.5),2x,a,2x,i2)') i+isource,"matched source",
     &           abs(ra_xx(m)),dec_xx(m),'source type',int(code/10000)
               write(12,'(4(es10.3,2x),2(f10.5,2x),f7.3,2x,2(f10.4,2x),i2)') frequency_xx(m),flux_xx(m),
     &         FluxU_xx(m),FluxL_xx(m),ra_xx(m),dec_xx(m),poserr_xx(m),mjdst_xx(m),mjded_xx(m),xx_type(m)
               do j=1,xpts(m)
                  write(12,'(4(es10.3,2x),i2)') frequency_xxot(j,m),flux_xxot(j,m),
     &             FluxU_xxot(j,m),FluxL_xxot(j,m),xxot_type(j,m)
               enddo
            enddo
         endif
      enddo
      write(12,*) ipass

      close(12)
      end

      SUBROUTINE DIST_SKY(alpha1,delta1,alpha2,delta2,dist)
      IMPLICIT NONE
      REAL*8 dist,alpha1,alpha2,delta1,delta2,costheta
      REAL*8 radian
      radian=57.2957795
      costheta=sin(delta1/radian)*sin(delta2/radian)+
     &         cos(delta1/radian)*cos(delta2/radian)*
     &         cos((alpha1-alpha2)/radian)
      dist=acos(costheta)*radian
      RETURN
      END

      SUBROUTINE fluxtofdens(gamma,bandl,bandu,flux,gev,fdens,nudens)
      real*4 bandu,bandl,flux,nudens,fdens,conval,gev,gamma
      !write(*,*) alpha,flux,kev,bandu,bandl
      if (gamma .ne. 2. ) then
      conval=(1./(-gamma+1.))*((bandu)**(-gamma+1.)-(bandl)**(-gamma+1.))
      else
      conval=log(bandu/bandl)
      endif
      fdens=gev*(flux/conval)*((gev)**(-gamma))!!!!To photon flux at gev
      fdens=fdens*1.602E-19*1.E7*(gev*1.E9)
      nudens=(1.602E-19)*(gev*1.E9)/(6.626e-34)
      RETURN
      end

      SUBROUTINE RXgraphic_code(flux,RX,code)
      IMPLICIT none
      REAL*4 flux,rfl_max,rfl_min
      REAL*4 xfl_min,xfl_max
      INTEGER*4 radio_component,x_ray_component,code
      CHARACTER*1 RX
      code = 0
      rfl_min=0.8 ! 0.8 mJy
      rfl_max=8000. ! 8 Jy
      xfl_min = 1.e-16 ! 1.e-16 erg/cm2/s, nufnu
      xfl_max = 5.e-11
      radio_component = 0.
      x_ray_component = 0.
      IF (RX == 'R') THEN
      radio_component=int(alog10(flux/rfl_min)/alog10(rfl_max/rfl_min)*99.)
      IF (radio_component .GE. 99) THEN
      radio_component = 99
      ELSE IF (radio_component .LE. 1) THEN
      radio_component = 1
      ENDIF
      code = -90000
      ELSE IF (RX == 'X') THEN
      IF (flux > xfl_min) THEN
      x_ray_component=int(alog10(flux/xfl_min)/alog10(xfl_max/xfl_min)*99.)
      ELSE
      x_ray_component = 1
      ENDIF
      IF (x_ray_component .GE. 99) THEN
      x_ray_component = 99
      ELSE IF (x_ray_component .Lt. 1) THEN
      x_ray_component = 1
      ENDIF
      code = -80000
      ENDIF
      code = code -radio_component-100*x_ray_component
      RETURN
      END

      SUBROUTINE mag2flux (nh,m_band,filter,flux,frequency)
c
c  converts u,v,i,h,b,r,j,k magnitudes into monochromatic fluxes
c  in units of erg/cm2/s for nufnu vs nu plots
c
      IMPLICIT none
      REAL*4 nh, flux , av , m_band, a_band, Rv
      REAL*4 c, lambda, const, frequency, a
      REAL*4 x,aa,bb,c1,c2,dx,px,ebv
      CHARACTER*3 filter
c        print *,' nh, m_band, filter ', nh, m_band,filter
c        call upc(filter)
      IF ((filter(1:3).NE.'fuv') .and. (filter(1:3).NE.'nuv')) THEN
         write (*,*) ' mag2flux: Filter not supported  '
         stop
      ENDIF
c extintion law taken from Cardelli et al. 1989 ApJ 345, 245
c in UV apply the UV relation from Fitzpatrick 1999
      Rv=3.1
      av = Rv*(-0.055+nh*1.987e-22) !!dust map from BH1978, assumed constant gas-to-dust ratio
      ebv=av/Rv
      if (av < 0.) av=0.
      if (filter(1:3) == 'fuv') then
         lambda=1528.
         const=log10(3631.)-23.
      else if (filter(1:3) == 'nuv') then
         lambda=2271.
         const=log10(3631.)-23.
      endif
c lambda from Amstrongs to microns
      x=10000./lambda
      if ((x .le. 1.1) .and. (x .ge. 0.3)) then
      a_band=(0.574*(x**1.61)-0.527*(x**1.61)/Rv)*av ! the a_lambda
      else if ((x .le. 3.3) .and. (x .ge. 1.1)) then
      x=x-1.82
      aa=1+(0.17699*x)-(0.50447*x**2)-(0.02427*x**3)+(0.73085*x**4)
     &     +(0.01979*x**5)-(0.77530*x**6)+(0.32999*x**7)
      bb=1.41338*x+(2.28305*x**2)+(1.07233*x**3)-(5.38434*x**4)
     &    -(0.662251*x**5)+(5.30260*x**6)-(2.09002*x**7)
      a_band=(aa+(bb/Rv))*av
      else if ((x .le. 10.) .and. (x .ge. 3.3)) then
      c2=-0.824+4.717/Rv
      c1=2.03-3.007*c2
      dx=(x*x)/((x**2-4.596**2)+(x*0.99)**2)
      px=0.5392*(x-5.9)**2+0.05644*(x-5.9)**2
      if (x .le. 5.9) px=0.
      a_band=(c1+c2*x+3.23*dx+0.41*px)*ebv+av
      endif
      c=2.9979e10
      a=1.0
      frequency=c/(lambda*1.e-8)
      flux = 10.**(-0.4*(m_band -a_band)+const)*frequency
      if (m_band .le. 0.) flux=0.
!m_band=0.
!lambda=0.
      RETURN
      END
