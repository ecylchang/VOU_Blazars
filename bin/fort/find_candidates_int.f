      PROGRAM find_candidates_int

      IMPLICIT none
      integer i,j,k,m,s,iuv,ihst,igaia,iunwise,igam,i4p8,irr,ixx,in,im,lenact,length,lu_in,ier,is,ie
      integer ii1,ii2,ipass,isource,icand,filen,iconfig,iarr,arrsize(30),offset,l,ipanstarrs,iradio
      integer nnuvx,nnox,nnruv,nnralpha,code,irrrep,ixxrep,irrss,ixxss,iaox,file_size
      real*8 ra,dec,dist,ra_center,dec_center,ragood,decgood
      real*4 nh,flux,uflux,lflux,epos,freq,radius,flux2nufnu_4p8,airo,p_err,arx,arg,aro
      real*4 aruv,auvx,aox,max_dist,alphar,mjdstart,mjdend,airx,hstbmag,aw1w2
      character*200 input_file,input_file2,input_file3,output_file,output_file2
      character*200 webprograms,array_size
      character*10 catalog
      character*800 string

      real*8,dimension(:),allocatable :: ra_uv,dec_uv,ra_hst,dec_hst,ra_unwise,dec_unwise,ra_gaia,dec_gaia
      real*8,dimension(:),allocatable :: ra_panstarrs,dec_panstarrs,ra_radio,dec_radio
      real*4,dimension(:),allocatable :: poserr_uv,poserr_hst,poserr_unwise,slopegaia 
      real*4,dimension(:),allocatable :: frequency_gaiag,flux_gaiag,frequency_gaiab,flux_gaiab
      real*4,dimension(:),allocatable :: frequency_gaiar,flux_gaiar
      real*4,dimension(:),allocatable :: frequency_panstarrsr,flux_panstarrsr,frequency_radio,flux_radio
      real*4,dimension(:),allocatable :: frequency_panstarrsg,flux_panstarrsg
      real*4,dimension(:),allocatable :: gaiagmag,gaiabmag,gaiarmag,panstarrsgmag,panstarrsrmag
      real*4,dimension(:,:),allocatable :: frequency_uv,flux_uv,FluxL_uv,FluxU_uv,uvmag,uvmagerr,flux_ir,irmag,frequency_ir
      real*4,dimension(:,:),allocatable :: frequency_hst,flux_hst,hstvmag

      real*8,dimension(:),allocatable :: ra_4p8,dec_4p8
      real*4,dimension(:),allocatable :: frequency_4p8,flux_4p8,FluxL_4p8,FluxU_4p8,poserr_4p8,Ferr_4p8
      character*10,dimension(:),allocatable :: type_4p8

      real*8,dimension(:),allocatable :: ra_gam,dec_gam
      real*4,dimension(:),allocatable :: slope_gam,specerr_gam,poserr_gam
      real*4,dimension(:),allocatable :: frequency_gam,flux_gam,FluxL_gam,FluxU_gam,Ferr_gam

      integer*4,dimension(:),allocatable :: xpts,pass2,rr_type,xx_type,xxss_type,rrss_type
      integer*4,dimension(:),allocatable :: nreprr,trackrr,repnumberrr,posindrr
      integer*4,dimension(:),allocatable :: nrepxx,trackxx,repnumberxx,posindxx
      real*8,dimension(:),allocatable :: ra_rr,dec_rr,ra_rrss,dec_rrss
      real*8,dimension(:),allocatable :: ra_xx,dec_xx,ra_xxss,dec_xxss

      real*4,dimension(:),allocatable :: frequency_rr,flux_rr,FluxL_rr,FluxU_rr,poserr_rr
      real*4,dimension(:),allocatable :: frequency_xx,flux_xx,FluxL_xx,FluxU_xx,poserr_xx
      real*4,dimension(:),allocatable :: mjdst_rr,mjded_rr,mjdst_xx,mjded_xx,poserr_xxss,flux_rrss

      integer*4,dimension(:,:),allocatable :: backxx,backrr
      integer*4,dimension(:,:),allocatable :: xxot_type
      real*4,dimension(:,:),allocatable :: frequency_xxot,flux_xxot,FluxL_xxot,FluxU_xxot

      LOGICAL there,ok,found,written
      common webprograms
      ok = .TRUE.
      found = .FALSE.
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
      in=im
      im=index(string(in+1:length),' ')+in
      output_file=string(in+1:im-1)
      webprograms=string(im+1:length)
      in = index(input_file(1:lenact(input_file)),'.')
      IF (in == 0) input_file(lenact(input_file)+1:lenact(input_file)+4) = '.csv'
      open(13,file=input_file3,status='old',iostat=ier)
      inquire(13, size=file_size)
      INQUIRE (FILE=input_file,EXIST=there)
      IF (.NOT.there .OR. file_size .EQ. 0) THEN
         !print *,'Here !'
         if (.NOT.there) write (*,'(''No data found in Intermediate phase '')')
         INQUIRE (FILE=input_file2,EXIST=there)
         if (there) then
            output_file='tmp/Sed_temp.txt'
            open(12,file=output_file,status='old',iostat=ier)
            open(11,file=input_file2,status='old',iostat=ier)
            i = 0
            do while(ok)
                i = i +1
                read(11,*,end=801) freq,flux,uflux,lflux,ra,dec,epos,mjdstart,mjdend,code
                if (i .EQ. 1) write(12,'("   1  matched source  ",f10.5,2x,f10.5,"  source type   2")')ra,dec
                write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') freq,flux,uflux,lflux,ra,dec,
     &                     epos,mjdstart,mjdend,abs(code)
            enddo
         endif
801      continue
         STOP
      ENDIF
      output_file2='tmp/candidates_int.txt'
      array_size=webprograms(1:lenact(webprograms))//'/array_size.cf'
      lu_in = 10
      open(lu_in,file=input_file,status='old',iostat=ier)
      open(11,file=input_file2,status='old',iostat=ier)
      open(18,file=array_size,status='old',iostat=ier)
      IF (ier.NE.0) THEN
         write (*,*) ' Error ',ier,' opening file ', input_file
      ENDIF

      iarr=0
      iconfig=0
      arrsize(1:30)=0
      do while(ok)
      read(18,'(a)',end=700) string
      if (string(1:3) == '---' ) then
         iconfig=iconfig+1
      else
         if (iconfig .eq. 4) then
            iarr=iarr+1
            is=index(string(1:lenact(string)),':')
            read(string(is+1:lenact(string)),*) arrsize(iarr)
         endif
      endif
      enddo
700   continue
      close(18)

      allocate(xxot_type(arrsize(3),arrsize(1)))
      allocate(frequency_xxot(arrsize(3),arrsize(1)),flux_xxot(arrsize(3),arrsize(1)),FluxL_xxot(arrsize(3),arrsize(1)),FluxU_xxot(arrsize(3),arrsize(1)))
      allocate(ra_rr(arrsize(2)),dec_rr(arrsize(2)),ra_xx(arrsize(1)),dec_xx(arrsize(1)))
      allocate(rr_type(arrsize(2)),xx_type(arrsize(1)),xpts(arrsize(1)))
      allocate(frequency_rr(arrsize(2)),flux_rr(arrsize(2)),FluxL_rr(arrsize(2)),FluxU_rr(arrsize(2)),poserr_rr(arrsize(2)))
      allocate(frequency_xx(arrsize(1)),flux_xx(arrsize(1)),FluxL_xx(arrsize(1)),FluxU_xx(arrsize(1)),poserr_xx(arrsize(1)))
      allocate(mjdst_rr(arrsize(2)),mjded_rr(arrsize(2)),mjdst_xx(arrsize(1)),mjded_xx(arrsize(1)))


      icand=0
      ixx=0
      irr=0
      Do WHILE(ok)
         read(11,*,end=100) freq,flux,uflux,lflux,ra,dec,epos,mjdstart,mjdend,code
         !print *,'-->',freq,flux,uflux,lflux,ra,dec,epos,mjdstart,mjdend,code
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

      isource=0
      do while(ok)
         read(13,'(a)',end=98) string 
         !print *,'string ',string(1:lenact(string))
         read(string,*,end=98) ra,dec,code
         !read(13,*,end=98) ra,dec,code
         if ((code .gt. 0.) .or. (code .le. -50000.) .or. (code .eq. -9999)) then
            isource=isource+1
         endif
      enddo
98    continue
      offset = isource
      !write(*,*) isource

      allocate(ra_hst(arrsize(6)),dec_hst(arrsize(6)),poserr_hst(arrsize(6)))
      allocate(ra_gaia(arrsize(6)),dec_gaia(arrsize(6)))
      allocate(ra_panstarrs(arrsize(6)),dec_panstarrs(arrsize(6)))
      allocate(ra_radio(arrsize(6)),dec_radio(arrsize(6)))
      allocate(ra_unwise(arrsize(6)),dec_unwise(arrsize(6)),poserr_unwise(arrsize(6)))
      allocate(frequency_hst(arrsize(6),2),flux_hst(arrsize(6),2),hstvmag(arrsize(6),2))
      allocate(frequency_gaiag(arrsize(6)),flux_gaiag(arrsize(6)),gaiagmag(arrsize(6)))
      allocate(frequency_gaiab(arrsize(6)),flux_gaiab(arrsize(6)),gaiabmag(arrsize(6)))
      allocate(frequency_gaiar(arrsize(6)),flux_gaiar(arrsize(6)),gaiarmag(arrsize(6)))
      allocate(frequency_panstarrsg(arrsize(6)),flux_panstarrsg(arrsize(6)),panstarrsgmag(arrsize(6)))
      allocate(frequency_panstarrsr(arrsize(6)),flux_panstarrsr(arrsize(6)),panstarrsrmag(arrsize(6)))
      allocate(frequency_radio(arrsize(6)),flux_radio(arrsize(6)))
      allocate(flux_ir(arrsize(6),2),irmag(arrsize(6),2),frequency_ir(arrsize(6),2),slopegaia(6))   

      allocate(ra_uv(arrsize(6)),dec_uv(arrsize(6)),poserr_uv(arrsize(6)))
      allocate(frequency_uv(arrsize(6),2),flux_uv(arrsize(6),2),FluxL_uv(arrsize(6),2),FluxU_uv(arrsize(6),2),uvmag(arrsize(6),2),uvmagerr(arrsize(6),2))

      allocate(ra_4p8(arrsize(5)),dec_4p8(arrsize(5)))
      allocate(frequency_4p8(arrsize(5)),flux_4p8(arrsize(5)),FluxL_4p8(arrsize(5)),FluxU_4p8(arrsize(5)),poserr_4p8(arrsize(5)),Ferr_4p8(arrsize(5)))
      allocate(type_4p8(arrsize(5)))

      allocate(ra_gam(arrsize(4)),dec_gam(arrsize(4)))
      allocate(slope_gam(arrsize(4)),specerr_gam(arrsize(4)),poserr_gam(arrsize(4)))
      allocate(frequency_gam(arrsize(4)),flux_gam(arrsize(4)),FluxL_gam(arrsize(4)),FluxU_gam(arrsize(4)),Ferr_gam(arrsize(4)))

      igam=0
      iuv=0
      ihst=0
      igaia=0
      ipanstarrs=0
      iunwise=0
      iradio=0
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
         !print *,'string ',string(1:lenact(string))
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
         else if (catalog(1:4) == 'nvss') then
            iradio = iradio + 1
            ra_radio(iradio)=ra
            dec_radio(iradio)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            !is=ie
            !ie=index(string(is+1:len(string)),',')+is
            read(string(is+1:ie-1),*) flux_radio(iradio)
            !print *,'flux_radio mjy NVSS iradio ',flux_radio(iradio),iradio
            flux_radio(iradio) = flux_radio(iradio)*1.4e9*1.e-26
            frequency_radio(iradio) = 1.4e9
         else if (catalog(1:5) == 'first') then
            iradio = iradio + 1
            ra_radio(iradio)=ra
            dec_radio(iradio)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            read(string(is+1:ie-1),*) flux_radio(iradio)
            !print *,'FIRST flux_radio mJy iradio ',flux_radio(iradio),iradio
            flux_radio(iradio) = flux_radio(iradio)*1.4e9*1.e-26
            frequency_radio(iradio) = 1.4e9
            !print *,'flux_radio(iradio) ',flux_radio(iradio)
         else if (catalog(1:4) == 'racs') then
            iradio = iradio + 1
            ra_radio(iradio)=ra
            dec_radio(iradio)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            read(string(is+1:ie-1),*) flux_radio(iradio)
            flux_radio(iradio) = flux_radio(iradio)*8.9e8*1.e-26
            frequency_radio(iradio) = 8.9e8
         else if (catalog(1:4) == 'sums') then
            iradio = iradio + 1
            ra_radio(iradio)=ra
            dec_radio(iradio)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            read(string(is+1:ie-1),*) flux_radio(iradio)
            flux_radio(iradio) = flux_radio(iradio)*8.9e8*1.e-26
            frequency_radio(iradio) = 8.9e8
         else if (catalog(1:7) == 'vlassql') then
            iradio = iradio + 1
            ra_radio(iradio)=ra
            dec_radio(iradio)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            read(string(is+1:ie-1),*) flux_radio(iradio)
            !print *,'VLASS flux_radio mJy iradio ',flux_radio(iradio), iradio
            flux_radio(iradio) = flux_radio(iradio)*3.0e9*1.e-26
            frequency_radio(iradio) = 8.9e8
         else if (catalog(1:9) == 'panstarrs') then
            ipanstarrs=ipanstarrs+1
            ra_panstarrs(ipanstarrs)=ra
            dec_panstarrs(ipanstarrs)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            is=ie
            read(string(is+1:lenact(string)),*) panstarrsgmag(ipanstarrs)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            is=ie
            read(string(is+1:lenact(string)),*) panstarrsrmag(ipanstarrs)
            if ( panstarrsgmag(ipanstarrs) < -99. ) then 
                if (panstarrsrmag(ipanstarrs) > -10.) then
                    panstarrsgmag(ipanstarrs) = panstarrsrmag(ipanstarrs)+1.0
                else
                   panstarrsgmag(ipanstarrs) = 25.
                endif
            endif
            CALL  mag2flux (nh,panstarrsgmag(ipanstarrs),'psg',flux_panstarrsg(ipanstarrs),frequency_panstarrsg(ipanstarrs))
            CALL  mag2flux (nh,panstarrsrmag(ipanstarrs),'psr',flux_panstarrsr(ipanstarrs),frequency_panstarrsr(ipanstarrs))
         else if (catalog(1:4) == 'gaia') then
            igaia=igaia+1
            ra_gaia(igaia)=ra
            dec_gaia(igaia)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            is=ie
            read(string(is+1:lenact(string)),*) gaiagmag(igaia)
            ie=index(string(is+1:len(string)),',')+is
            is=ie
            read(string(is+1:lenact(string)),*) gaiabmag(igaia)
            ie=index(string(is+1:len(string)),',')+is
            is=ie
            read(string(is+1:lenact(string)),*) gaiarmag(igaia)
            !print *,'gaiagmag gaiabmag gaiarmag ',gaiagmag(igaia),gaiabmag(igaia),gaiarmag(igaia)
            CALL  mag2flux (nh,gaiagmag(igaia),'gG ',flux_gaiag(igaia),frequency_gaiag(igaia))
            CALL  mag2flux (nh,gaiabmag(igaia),'bpG',flux_gaiab(igaia),frequency_gaiab(igaia))
            CALL  mag2flux (nh,gaiarmag(igaia),'rpG',flux_gaiar(igaia),frequency_gaiar(igaia))
            !print *,'flux_gaiar flux_gaiab ',flux_gaiar(igaia),flux_gaiab(igaia)
            !print *,'frequency_gaiar frequency_gaiab ',frequency_gaiar(igaia),frequency_gaiab(igaia)
         else if (catalog(1:3) == 'hst') then
            ihst=ihst+1
            ra_hst(ihst)=ra
            dec_hst(ihst)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            !print *,'is,is ',is,ie
            if (ie > is+1) then 
                read(string(is+1:ie-1),*) poserr_hst(ihst) 
            else
                poserr_hst(ihst) = 2.
            endif
            do i = 1,2
               is=ie
               ie=index(string(is+1:len(string)),',')+is
            enddo
            read(string(is+1:ie-1),*) hstbmag 
            do i = 1,2
               is=ie
               ie=index(string(is+1:len(string)),',')+is
            enddo
            read(string(is+1:ie-1),*) hstvmag(ihst,1)
            if (hstvmag(ihst,1) > 25.) then 
                if (hstbmag < 25.) then
                   hstvmag(ihst,1) = hstbmag
                else
                   hstvmag(ihst,1) = 25.
                endif
            endif
            CALL  mag2flux (nh,hstvmag(ihst,1),'V  ',flux_hst(ihst,1),frequency_hst(ihst,1))
         else if (catalog(1:6) == 'unwise') then
            iunwise = iunwise +1
            ra_unwise(iunwise)=ra
            dec_unwise(iunwise)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            read(string(is+1:ie-1),*) irmag(iunwise,1)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            read(string(is+1:ie-1),*) irmag(iunwise,2)
            if ( (irmag(iunwise,1) > 0.0) .and. (irmag(iunwise,2) > 0.0) ) then
               irmag(iunwise,1) = 22.5-2.5*log10(irmag(iunwise,1)) ! this conversion from vizier flux to mag comes from a note in vizier catalog
               irmag(iunwise,2) = 22.5-2.5*log10(irmag(iunwise,2)) ! this conversion from vizier flux to mag comes from a note in vizier catalog
               call mag2flux (nh,irmag(iunwise,1),'ww1',flux_ir(iunwise,1),frequency_ir(iunwise,1))
               call mag2flux (nh,irmag(iunwise,2),'ww2',flux_ir(iunwise,2),frequency_ir(iunwise,2))
               !print *,'irmag(iunwise,1),irmag(iunwise,2)',irmag(iunwise,1),irmag(iunwise,2)
               !print *,'flux_ir(iunwise,1) flux_ir(iunwise,2)',flux_ir(iunwise,1),flux_ir(iunwise,2)
            else
               flux_ir(iunwise,1) = 0.
               flux_ir(iunwise,2) = 0.
            endif
         endif
      ENDDO
99    continue
      close(lu_in)
      !close(11)

      allocate(xxss_type(arrsize(1)),rrss_type(arrsize(2)),pass2(arrsize(1)+arrsize(2)))
      allocate(nreprr(arrsize(2)),trackrr(arrsize(2)),repnumberrr(arrsize(2)),posindrr(arrsize(2)))
      allocate(nrepxx(arrsize(1)),trackxx(arrsize(1)),repnumberxx(arrsize(1)),posindxx(arrsize(1)))
      allocate(ra_rrss(arrsize(2)),dec_rrss(arrsize(2)),flux_rrss(arrsize(2)))
      allocate(ra_xxss(arrsize(1)),dec_xxss(arrsize(1)),poserr_xxss(arrsize(1)))
      allocate(backxx(arrsize(1),arrsize(1)),backrr(arrsize(2),arrsize(2)))

      nrepxx(1:ixx)=0
      ixxss=0
      ixxss=0
      ixxrep=0
      do i=1,ixx
         if (i .ne. 1) then
            do j=1,i-1
               call DIST_SKY(abs(ra_xx(i)),dec_xx(i),abs(ra_xx(j)),dec_xx(j),dist)
               if (dist*3600. .lt. 18.) then
                  ixxrep=ixxrep+1
                  trackxx(i)=trackxx(j)
                  nrepxx(trackxx(i))=nrepxx(trackxx(i))+1
                  repnumberxx(i)=nrepxx(trackxx(i))
                  goto 97
               endif
            enddo
         endif
         ixxss=ixxss+1
         trackxx(i)=ixxss
         nrepxx(trackxx(i))=1
         repnumberxx(i)=nrepxx(trackxx(i))
97       continue
      enddo
      if (ixxss+ixxrep .ne. ixx ) write(*,*) 'Warning! may have the wrong number.'

c xrt 1sxps chandra
c bmw xmm xmmsl
c rass wga ipc

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
      enddo

      nreprr(1:irr)=0
      irrss=0
      irrrep=0
      do i=1,irr
         if (i .ne. 1) then
            do j=1,i-1
               call DIST_SKY(ra_rr(i),dec_rr(i),ra_rr(j),dec_rr(j),dist)
               if (dist*3600. .lt. 18.) then
                  irrrep=irrrep+1
                  trackrr(i)=trackrr(j) !for repeated sources track back to the number
                  nreprr(trackrr(i))=nreprr(trackrr(i))+1 !repeated number
                  repnumberrr(i)=nreprr(trackrr(i))
                  goto 96
               endif
            enddo
         endif
         irrss=irrss+1
         trackrr(i)=irrss
         nreprr(trackrr(i))=1
         repnumberrr(i)=nreprr(trackrr(i))
96    continue
      enddo
      !print *,'irrss - ',irrss
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
      enddo

      open(12,file=output_file,status='unknown',iostat=ier)
      open(15,file=output_file2,status='unknown',iostat=ier)
      ipass=0
      !print *,'irrss ',irrss
      do i=1,irrss
         ii1=0
         ii2=0
         do m=1,nreprr(i)
            k=backrr(i,m)
            nnruv=0
            nnralpha=0
            do j=1,iuv
               call DIST_SKY(ra_rr(k),dec_rr(k),ra_uv(j),dec_uv(j),dist)
               max_dist=sqrt(poserr_uv(j)**2+poserr_rr(k)**2)
               if ((dist*3600. .lt. max(max_dist,2.)) .and. (flux_rr(k)/frequency_rr(k) .ge. 10*1.e-26)) then
                  if (flux_uv(j,1) .ne. 0.) then
                     aruv = 1.-log10(flux_rr(k)/flux_uv(j,1))/log10(frequency_rr(k)/frequency_uv(j,1))
                  else
                     aruv = 1.-log10(flux_rr(k)/flux_uv(j,2))/log10(frequency_rr(k)/frequency_uv(j,2))
                  endif
                  if (aruv .le. 0.75 ) then
                     ii1=ii1+1 !!!remove UV-r slope strange source 0.85
                     nnruv=nnruv+1
                  endif
               endif
            enddo
            do j=1,i4p8
               call DIST_SKY(ra_rr(k),dec_rr(k),ra_4p8(j),dec_4p8(j),dist)
               max_dist=sqrt(poserr_4p8(j)**2+poserr_rr(k)**2)
               if ((dist*3600. .lt. max_dist ) .and. (flux_rr(k)/frequency_rr(k) .ge. 10*1.e-26)) then
                  alphar = 1.-log10(flux_rr(k)/flux_4p8(j))/log10(frequency_rr(k)/frequency_4p8(j))
                  if (alphar .le. 0.7 ) then
                     ii2=ii2+1 !!!remove radio extended sources
                     nnralpha=nnralpha+1
                     if ((ii2 .eq. 1) .and. (ii1 .eq. 0)) then
                         !write(*,'(a,i4,a)') "----------------------------------------"
                         !write (*,'(a,i4)') "source nr.", isource+ii1+ii2
                         offset = offset +ii1+ii2
                     endif
c                     !if (((ii2 .eq. 1) .and. (ii1 .eq. 0)) .or. ((m .ne. 1) .and. (nnruv .eq. 0)
c     !&                .and. (nnralpha .eq. 1))) write(*,'(f9.5,2x,f9.5,a,f9.3,a,2x,i2)')  ra_rr(k),
c     !&                dec_rr(k)," radio source ",flux_rr(k)/frequency_rr(k)/1.E-26," mJy",rr_type(k)
                  endif
               endif
            enddo
         enddo
         if (ii1+ii2 .gt. 0.) then
            ipass=ipass+1
            CALL RXgraphic_code(flux_rr(posindrr(i))/frequency_rr(posindrr(i))/1.E-26,'R',code)
            write (12,'(f9.5,2x,f9.5,2x,i6," N.A.")') ra_rr(posindrr(i)),dec_rr(posindrr(i)),int(code)
            pass2(ipass)=i
         endif
      enddo
      !print *,'ixxss ',ixxss
      do i=1,ixxss
         ii1=0
         do m=1,nrepxx(i)
            k=backxx(i,m)
            nnuvx=0
            do j=1,iuv
               call DIST_SKY(abs(ra_xx(k)),dec_xx(k),ra_uv(j),dec_uv(j),dist)
               max_dist=sqrt(poserr_uv(j)**2+poserr_xx(k)**2)
               if ((dist*3600. .le. max(max_dist,2.) ) .and. (poserr_xx(k) .le. 20.)
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
                        !write(*,'(a,i4,a)') "----------------------------------------"
                        !write(*,'(a,i4)') "source nr.", ipass+ii1+isource
                        offset = offset + ipass+ii1
                     endif
                  endif
               endif
            enddo
         enddo
         if (ii1 .gt. 0.) then
             ipass=ipass+1
             CALL RXgraphic_code(flux_xx(posindxx(i)),'X',code)
             write (12,'(f9.5,2x,f9.5,2x,i6," N.A.")') abs(ra_xx(posindxx(i))),dec_xx(posindxx(i)),int(code)
             pass2(ipass)=i+irr
         endif
      enddo
      write(15,'(''alpha_ox,alpha_irx,alpha_iro,alpha_ro,alpha_rx,alpha_w1w2,alpha_rmag-gmag,ra,dec,fx_1kev'')')
      do i=1,ixxss
         ii1=0
         iaox = 0
         airx = 0
         arx = 0
         do m=1,nrepxx(i)
           k=backxx(i,m)
           nnox=0
           if (flux_xx(k) > 0. ) then 
              if ( (ipanstarrs > 0) .and. (igaia < 1) ) then
                 do j=1,ipanstarrs
                    aox = -99.
                    arg = -99.
                    call DIST_SKY(abs(ra_xx(k)),dec_xx(k),ra_panstarrs(j),dec_panstarrs(j),dist)
                    if (poserr_xx(k) < 3.0) poserr_xx(k) = 3.0
                    p_err=1.0
                    max_dist=sqrt(p_err**2+poserr_xx(k)**2)*1.25 ! increase matching radius by 25% 
                    if ( (dist*3600. .le. max_dist) .and. (flux_xx(k) .ne. 0. ) ) then
                        if ( (flux_panstarrsg(j) > 0.) .and. (panstarrsgmag(j) < 24.9) ) then
                           aox  = 1.-log10(flux_panstarrsg(j)/flux_xx(k))/log10(frequency_panstarrsg(j)/frequency_xx(k))
                        endif
                        airx = -99.
                        airo = -99.
                        aw1w2 = -99.
                        if (iunwise > 0) then 
                           do l=1,iunwise
                              call DIST_SKY(ra_panstarrs(j),dec_panstarrs(j),ra_unwise(l),dec_unwise(l),dist)
                              if ( (dist*3600. < 3.) .and. (flux_ir(l,2) > 8.e-14) .and. (flux_xx(k) .ne. 0. )) then
                                  airx = 1.-log10(flux_ir(l,2)/flux_xx(k))/log10(frequency_ir(l,2)/frequency_xx(k)) 
                                  aw1w2 = log10(flux_ir(l,1)/flux_ir(l,2))/log10(frequency_ir(l,1)/frequency_ir(l,2))
                                  if ( flux_panstarrsg(j) > 0.) then
                                     airo  = 1.-log10(flux_ir(l,2)/flux_panstarrsg(j))/log10(frequency_ir(l,2)/frequency_panstarrsg(j))
                                  endif
                              endif
                           enddo
                        endif
                        written = .FALSE.
                        if (iradio > 0) then
                           do l=1,iradio
                              arx = -99.
                              aro = -99.
                              call DIST_SKY(ra_panstarrs(j),dec_panstarrs(j),ra_radio(l),dec_radio(l),dist)
                              if (flux_radio(l) > 1.e-14) then
                                 max_dist = 8
                              else
                                 max_dist = 6
                              endif
                              if (dist*3600. .le. max_dist) then
                                 if ( ( flux_radio(l) > 0. ) .and. (flux_xx(k) .ne. 0. ) ) then
                                    arx  = 1.-log10(flux_radio(l)/flux_xx(k))/log10(frequency_radio(l)/frequency_xx(k))
                                    aro  = 1.-log10(flux_radio(l)/flux_panstarrsg(j))/log10(frequency_radio(l)/frequency_panstarrsg(j))
                                 endif
                                 if ( (aox > -90) .OR. (airx > -90) .OR. (airo > -90) .OR. (arx > -90.) .or. (arg > -90.) ) then
                                    !print *,'aox airx airo arx ',aox,airx,airo,arx
                                    ragood = ra_panstarrs(j)
                                    decgood = dec_panstarrs(j)
                                    write (15,'(f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",
     &                                          f6.2,", ",f6.2,", ",f9.5,", ",f9.5,", ",e9.3)') 
     &                                          aox,airx,airo,aro,arx,aw1w2,arg,ragood,decgood,flux_xx(k)
                                    written = .TRUE.
                                 endif
                              endif
                           enddo
                        endif
                        if (.NOT. written) then 
                           if ( (aox > -90) .OR. (airx > -90) .OR. (airo > -90) .OR. (arx > -90.) .OR. (aro > -90.)) then
                              ragood = ra_panstarrs(j)
                              decgood = dec_panstarrs(j)
                              write (15,'(f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f9.5,", ",f9.5,", ",e9.3)')
     &                                    aox,airx,airo,aro,arx,aw1w2,arg,ragood,decgood,flux_xx(k)
                           endif
                        endif
                    endif
                 enddo
              else if (igaia > 0) then
                 do j=1,igaia
                    aox = -99.
                    call DIST_SKY(abs(ra_xx(k)),dec_xx(k),ra_gaia(j),dec_gaia(j),dist)
                    if (poserr_xx(k) < 3.) poserr_xx(k) = 3.0
                    p_err=1.0
                    max_dist=sqrt(p_err**2+poserr_xx(k)**2)*1.25 ! increase matching radius by 25% 
                    !print *,'gaiagmag(j) dist*3600. max_dist ',gaiagmag(j),dist*3600.,max_dist
                    if ( (dist*3600. .le. max_dist) .and. (flux_xx(k) .ne. 0. ) ) then
                        if ( gaiagmag(j) < 24. ) then
                           aox  = 1.-log10(flux_gaiag(j)/flux_xx(k))/log10(frequency_gaiag(j)/frequency_xx(k))
                           if ( flux_gaiab(j) > 0. .AND. flux_gaiar(j) > 0.) then
                              arg  = 1.-log10(flux_gaiar(j)/flux_gaiab(j))/
     &                                  log10(frequency_gaiar(j)/frequency_gaiab(j))
                              !print *,'arg from gaia ',arg
                           endif
                        endif
                        airx = -99.
                        airo = -99.
                        aw1w2 = -99.
                        if (iunwise > 0) then 
                           do l=1,iunwise
                              call DIST_SKY(ra_gaia(j),dec_gaia(j),ra_unwise(l),dec_unwise(l),dist)
                              if ( (dist*3600. < 3.) .and. (flux_ir(l,2) > 1.e-13) .and. (flux_xx(k) .ne. 0. )) then
                                  airx = 1.-log10(flux_ir(l,2)/flux_xx(k))/log10(frequency_ir(l,2)/frequency_xx(k)) 
                                  aw1w2 = log10(flux_ir(l,1)/flux_ir(l,2))/log10(frequency_ir(l,1)/frequency_ir(l,2))
                                  if ( gaiagmag(j) < 24. ) then
                                     airo  = 1.-log10(flux_ir(l,2)/flux_gaiag(j))/log10(frequency_ir(l,2)/frequency_gaiag(j))
                                     !print *,'aox airx airo',aox,airx,airo
                                  endif
                              endif
                           enddo
                        endif
                        written = .FALSE.
                        if (iradio > 0) then
                           do l=1,iradio
                              arx = -99.
                              aro = -99.
                              call DIST_SKY(ra_gaia(j),dec_gaia(j),ra_radio(l),dec_radio(l),dist)
                              if (flux_radio(l) > 1.e-14) then
                                 max_dist = 8
                              else
                                 max_dist = 5
                              endif
                              dist = dist * 3600.
                              if (dist .le. max_dist) then
                                 if ( ( flux_radio(l) > 0. ) .and. (flux_xx(k) > 0. ) ) then
                                    arx  = 1.-log10(flux_radio(l)/flux_xx(k))/log10(frequency_radio(l)/frequency_xx(k))
                                    aro  = 1.-log10(flux_radio(l)/flux_gaiag(j))/log10(frequency_radio(l)/frequency_gaiag(j))
                                 endif
                                 if ( (aox > -90) .OR. (airx > -90) .OR. (airo > -90) .OR. (arx > -90.) ) then
                                    ragood = ra_gaia(j)
                                    decgood = dec_gaia(j)
                                    write (15,'(f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",
     &                                          f6.2,", ",f6.2,", ",f9.5,", ",f9.5,", ",e9.3)') 
     &                                          aox,airx,airo,aro,arx,aw1w2,arg,ragood,decgood,flux_xx(k)
                                    written = .TRUE.
                                 endif
                              endif
                           enddo
                        endif
                        if (.NOT. written) then 
                           if ( (aox > -90) .OR. (airx > -90) .OR. (airo > -90) .OR. (arx > -90.) .OR. (aro > -90.)) then
                              ragood = ra_gaia(j)
                              decgood = dec_gaia(j)
                              write (15,'(f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f9.5,", ",f9.5,", ",e9.3)')
     &                                    aox,airx,airo,aro,arx,aw1w2,arg,ragood,decgood,flux_xx(k)
                           endif
                        endif
                    endif

                 enddo
              else if (ihst > 0) then
                do j=1,ihst
                  aox = -99.
                  call DIST_SKY(abs(ra_xx(k)),dec_xx(k),ra_hst(j),dec_hst(j),dist)
                  if (poserr_xx(k) < 1.2) poserr_xx(k) = 1.2
                  max_dist=sqrt(poserr_hst(j)**2+poserr_xx(k)**2)*1.25 ! increase matching radius by 25% 
                  !print *,'dist max_dist poserr_xx(k) ',dist*3600.,max_dist,poserr_xx(k)
                  if ( (dist*3600. .le. max_dist) .and. (flux_xx(k) .ne. 0. ) ) then
                      if (aox .EQ. -99.) then 
                         if (hstvmag(j,1) > 24.) then
                            aox = -99.
                         else if (flux_xx(k) > 0.) then
                            aox  = 1.-log10(flux_hst(j,1)/flux_xx(k))/log10(frequency_hst(j,1)/frequency_xx(k))
                         endif
                      endif
                      airx = -99.
                      airo = -99.
                      aw1w2 = -99.
                      do l=1,iunwise
                         call DIST_SKY(ra_hst(j),dec_hst(j),ra_unwise(l),dec_unwise(l),dist)
                         if ( (dist*3600. < 2.) .and. (flux_ir(l,2) > 1.e-13) .and. (flux_xx(k) > 0.) ) then
                            airx = 1.-log10(flux_ir(l,2)/flux_xx(k))/log10(frequency_ir(l,2)/frequency_xx(k)) 
                            aw1w2 = log10(flux_ir(l,1)/flux_ir(l,2))/log10(frequency_ir(l,1)/frequency_ir(l,2))
                            if (hstvmag(j,1) > 24.) then
                               airo = -99.
                            else
                               airo = 1.-log10(flux_ir(l,2)/flux_hst(j,1))/log10(frequency_ir(l,2)/frequency_hst(j,1)) 
                            endif
                         endif
                      enddo
                      if (iradio > 0) then
                         do l=1,iradio
                            arx = -99.
                            aro = -99.
                            call DIST_SKY(ra_hst(j),dec_hst(j),ra_radio(l),dec_radio(l),dist)
                            if (flux_radio(l) > 1.e-14) then
                               max_dist = 8
                            else
                               max_dist = 4
                            endif
                            if (dist*3600. .le. max_dist) then
                               if ( ( flux_radio(l) > 0. ) .and. (flux_xx(k) .ne. 0. ) ) then
                                  arx  = 1.-log10(flux_radio(l)/flux_xx(k))/log10(frequency_radio(l)/frequency_xx(k))
                                  aro  = 1.-log10(flux_radio(l)/flux_hst(j,1))/log10(frequency_radio(l)/frequency_hst(j,1))
                               endif
                               if ( (aox > -90) .OR. (airx > -90) .OR. (airo > -90) .OR. (arx > -90.) .OR. (aro > -90.)) then
                               !print *,'aox airx airo arx ',aox,airx,airo,arx
                                  ragood = ra_gaia(j)
                                  decgood = dec_gaia(j)
                                  write (15,'(f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f9.5,", ",f9.5,", ",e9.3)') 
     &                                        aox,airx,airo,aro,arx,aw1w2,arg,ragood,decgood,flux_xx(k)
                               endif
                            endif
                         enddo
                      endif
                      if (aox .le. 1.1) then
                         ii1=ii1+1 !!!remove o-x slope strange source 0.85
                         nnox=nnox+1
                         if (ii1 .eq. 1) then
                            !write(*,'("----------------------------------------")')
                            !write(*,'(a,i4)') "source nr.", ipass+ii1+isource
                            offset = offset + ipass+ii1
                         endif
                      endif
                      ragood = ra_hst(j)
                      decgood = dec_hst(j)
                      write (15,'(f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f6.2,", ",f9.5,", ",f9.5,", ",e9.3)') 
     &                            aox,airx,airo,aro,arx,aw1w2,arg,ragood,decgood,flux_xx(k)
                  endif
               enddo
            endif
           endif
         enddo
         if (ii1 .gt. 0.) then
             ipass=ipass+1
             CALL RXgraphic_code(flux_xx(posindxx(i)),'X',code)
             write (12,'(f9.5,2x,f9.5,2x,i6," N.A.")') abs(ra_xx(posindxx(i))),dec_xx(posindxx(i)),int(code)
             pass2(ipass)=i+irr
         endif
      enddo

      deallocate(ra_uv,dec_uv,poserr_uv)
      deallocate(frequency_uv,flux_uv,FluxL_uv,FluxU_uv,uvmag,uvmagerr)

      deallocate(ra_4p8,dec_4p8)
      deallocate(frequency_4p8,flux_4p8,FluxL_4p8,FluxU_4p8,poserr_4p8,Ferr_4p8)
      deallocate(type_4p8)

      deallocate(ra_gam,dec_gam)
      deallocate(slope_gam,specerr_gam,poserr_gam)
      deallocate(frequency_gam,flux_gam,FluxL_gam,FluxU_gam,Ferr_gam)

      if ((ipass .eq. 0) .and. (isource .ne. 0)) then
         print *,achar(27),'[35;1m No Candidates Found in the Intermediate Phase',achar(27),'[0m'
         stop
      endif
      if (ipass+isource .eq. 0) then
         print *,achar(27),'[31;1m No Candidates Found！',achar(27),'[0m'
         write(12,'("No Blazar Candidates Found")') 
         stop
      endif

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
     &         FluxU_xx(m),FluxL_xx(m),ra_xx(m),dec_xx(m),poserr_xx(m),mjdst_xx(m),mjded_xx(m),xx_type(m)+50
               do j=1,xpts(m)
                  write(12,'(4(es10.3,2x),i2)') frequency_xxot(j,m),flux_xxot(j,m),
     &             FluxU_xxot(j,m),FluxL_xxot(j,m),xxot_type(j,m)
               enddo
            enddo
         endif
      enddo
      write(12,*) ipass
      deallocate(ra_rr,dec_rr,ra_xx,dec_xx)
      deallocate(rr_type,xx_type,xpts)
      deallocate(frequency_rr,flux_rr,FluxL_rr,FluxU_rr,poserr_rr)
      deallocate(frequency_xx,flux_xx,FluxL_xx,FluxU_xx,poserr_xx)
      deallocate(mjdst_rr,mjded_rr,mjdst_xx,mjded_xx)

      deallocate(backxx,backrr)
      deallocate(xxot_type)
      deallocate(frequency_xxot,flux_xxot,FluxL_xxot,FluxU_xxot)

      deallocate(xxss_type,rrss_type)
      deallocate(nreprr,trackrr,repnumberrr,posindrr)
      deallocate(nrepxx,trackxx,repnumberxx,posindxx)
      deallocate(ra_rrss,dec_rrss,flux_rrss)
      deallocate(ra_xxss,dec_xxss,poserr_xxss)

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
      REAL*4 nh, flux, av , m_band, a_band, Rv
      REAL*4 c, lambda, frequency, a
      REAL*4 x,aa,bb,c1,c2,dx,px,ebv
      REAL*8 const
      CHARACTER*3 filter
      IF ( (filter(1:3).NE.'fuv') .and. (filter(1:3).NE.'nuv') .and. 
     &     (filter(1:1).NE.'V')   .and. (filter(1:3).NE.'ww1') .and. 
     &     (filter(1:3).NE.'psg') .and. (filter(1:3).NE.'psr')   .and. 
     &     (filter(1:3).NE.'psi') .and. (filter(1:3).NE.'psz')   .and. 
     &     (filter(1:3).NE.'psy') .and. 
     &     (filter(1:1).NE.'g')   .and. (filter(1:1).NE.'u')   .and. 
     &     (filter(1:2).NE.'gG') .and. (filter(1:3).NE.'bpG') .and.
     &     (filter(1:2).NE.'rpG') .and. (filter(1:3).NE.'rvs') .and.
     &     (filter(1:1).NE.'r')   .and. (filter(1:3).NE.'ww2') ) THEN
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
      else if (filter(1:1) == 'V') then
         lambda=5500.d0
         const=log10(3640.d0)-23.d0
      else if (filter(1:3) == 'ww1') then
         lambda=34000.
         const=log10(309.540)-23.
      else if (filter(1:3) == 'ww2') then 
         lambda=46000.
         const=log10(171.787)-23.
      else if (filter(1:3) == 'u  ') then !effective wavelength from SDSS, Doi et al. 2010 ApJ 139, 1628
         !m_band=m_band-0.04 !calibrate of the SDSS u band to AB mag
         lambda=3568.d0
         const=log10(3631.)-23.d0 !3631 is the 0 mag flux of AB mag system
      else if (filter(1:3) == 'g  ') then
         lambda=4653.d0
         const=log10(3631.d0)-23.d0
      else if (filter(1:3) == 'psg') then !tonry et al. 2012 Panstarrs mags are in the AB magnitude system
         lambda=4814.d0
         const=log10(3631.d0)-23.d0
      else if (filter(1:3) == 'psr') then
         lambda=6176.d0
         const=log10(3631.d0)-23.d0
      else if (filter(1:3) == 'psi') then
         lambda=7522.d0
         const=log10(3631.d0)-23.d0
      else if (filter(1:3) == 'psz') then
         lambda=8665.d0
         const=log10(3631.d0)-23.d0
      else if (filter(1:3) == 'psy') then
         lambda=9620.d0
         const=log10(3631.d0)-23.d0
      else if (filter(1:3) == 'rpG') then
         lambda=7770.d0 ! GAIA DR3
         const=log10(2555.d0)-23.d0 !!!!the zero mag. flux are estimated from Vega flux(22.3061)!!!
      else if (filter(1:3) == 'bpG') then
         lambda=5111.d0 ! GAIA DR3
         const=log10(3552.d0)-23.d0 !!!!the zero mag. flux are estimated from Vega flux(23.4318)!!!
      else if (filter(1:2) == 'gG') then
         lambda=6730.d0 !!!Jordi et al. 2010
         !lambda=6405.d0 ! Gaia DR2 in flight calibration
         const=log10(3296.d0)-23.d0 !!!!the zero mag. flux are estimated from Vega flux(25.6885)!!!
      else if (filter(1:3) == 'r  ') then
         lambda=6400.d0
         const=log10(3080.d0)-23.d0
      endif
c lambda from Amstrongs to microns
      x=10000./lambda
      a_band = 0.
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
      if (m_band .le. 0.) then
         flux=0. 
      else
         flux = 10.**(-0.4*(m_band -a_band)+const)*frequency
      endif
      RETURN
      END
