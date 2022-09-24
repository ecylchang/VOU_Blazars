       PROGRAM find_candidates1
c
c This program reads a query results file from the ASDC error circle explorer tool
c and plots potential blazars of different types based on radio/X-ray ratios 
c
c The output file should be used as input to the program gnomo_plot_types 
c that generates a postscript file with a plot in gnomonic coordinates. 
c In this plot each radio-X/ray match appears as a filled and an open circle with size that is 
c proportional to the radio flux density (filled circles) and to X-ray flux (open circles)
c the color coding is as follows
c orange : possible blazar of the HBL type
c Cyan   : possible blazar of the IBL type 
c Blue   : possible blazar of the LBL type 
c Green  : possible non-jetted AGN  
c black  : unknown type
c
c Sources included in blazar catalogs are plotted as:
c Red filled triangle if included in the BZCAT5 catalog
c Red open star symbol if included in the 2WHSP catalog
c
c A red "question mark" symbol is plotted if the Radio/X-ray match is 
c within a knwon cluster of galaxy. This is to warn that the X-ray could be 
c extended and due to the cluster rather than from the radio source.
c
      IMPLICIT none
      INTEGER*4 ier, lu_in, in,k, length,im,is,ie, i, j,l
      INTEGER*4 lenact,source_type,type_average,types(0:5)
      INTEGER*4 no_found,sfound,rfound,s,aim,xray_type,ncat,ifound
      INTEGER*4 iradio,ixmm,irosat,iswift,iipc,iother,ichandra,ibmw
      INTEGER*4 rah, ram, id, dm ,filen,ix,xpts,iarr,iconfig
      INTEGER*4 igrb,ierosita,imaxi,ibigb,igam,arrsize(30)
      REAL*8 ra, dec,dist,ra_center,dec_center,radius
      real*8 ra_bary,dec_bary,ra_cattemp,dec_cattemp
      REAL*4 min_dist_rosat,min_dist_xmm,min_dist_ipc,min_dist_cluster
      REAL*4 min_dist_other,min_dist_swift,min_dist_bmw,min_dist_chandra
      REAL*4 flux2nufnu_vlass,flux2nufnu_nvss,flux2nufnu_rosat
      REAL*4 flux2nufnu_swift,flux2nufnu_ipc,flux2nufnu_bmw,flux2nufnu_sumss
      REAL*4 radian,rasec,decsec,code,fdens,nudens,ratio,reduce
      REAL*4 flux_x,nh,errfrx,totweight,erraxis,totxerr,min_dist
      real*4 major,minor,posang,posxerr,posyerr
      real*4 errrad,errmaj,errmin,errang,mjdavg
      CHARACTER*1 sign
      CHARACTER*200 input_file,output_file,output_file2,output_file3
      character*200 output_file4,webprograms,array_size!,output_file5
      CHARACTER*15 catalog
      CHARACTER*800 string

      integer*4,dimension(:,:),allocatable :: spec_type,spec_xpts
      real*8,dimension(:,:),allocatable :: ra_1kev,dec_1kev,distrx
      real*4,dimension(:,:),allocatable :: flux_1kev,uflux_1kev,lflux_1kev,uflux_xpts,lflux_xpts,flux_xpts,frequency_xpts
      real*4,dimension(:,:),allocatable :: poserr_1kev,mjdstart,mjdend
      real*8,dimension(:),allocatable :: ra_radio,dec_radio
      integer*4,dimension(:),allocatable :: radio_type,ra_index
      real*4,dimension(:),allocatable :: ppss,const,Ferr_radio,FluxU_radio,FluxL_radio,poserr_radio,flux_radio,frequency_radio

      integer*4,dimension(:),allocatable :: xmm_type,xrt_type,rosat_type
      real*8,dimension(:),allocatable :: ra_xmm,dec_xmm,ra_swift,dec_swift
      real*4,dimension(:),allocatable :: poserr_xmm,poserr_swift,mjdst_swift,mjded_swift
      real*4,dimension(:,:),allocatable :: flux_xmm,Ferr_xmm,FluxU_xmm,FluxL_xmm,frequency_xmm
      real*4,dimension(:,:),allocatable :: flux_swift,Ferr_swift,FluxU_swift,FluxL_swift,frequency_swift

      real*8,dimension(:),allocatable :: ra_rosat,dec_rosat,ra_chandra,dec_chandra
      real*4,dimension(:),allocatable :: poserr_rosat,poserr_chandra
      real*4,dimension(:),allocatable :: flux_rosat,Ferr_rosat,FluxU_rosat,FluxL_rosat,frequency_rosat
      real*4,dimension(:,:),allocatable :: flux_chandra,FluxU_chandra,FluxL_chandra,frequency_chandra

      real*8,dimension(:),allocatable :: ra_erosita,dec_erosita,ra_bmw,dec_bmw
      real*4,dimension(:),allocatable :: poserr_erosita,poserr_bmw
      real*4,dimension(:,:),allocatable :: flux_erosita,fluxL_erosita,fluxU_erosita,Ferr_erosita,frequency_erosita
      real*4,dimension(:),allocatable :: flux_bmw,fluxL_bmw,fluxU_bmw,Ferr_bmw,frequency_bmw

      integer*4,dimension(:),allocatable :: ipc_type,maxi_type
      real*8,dimension(:),allocatable :: ra_ipc,dec_ipc,ra_maxi,dec_maxi
      real*4,dimension(:),allocatable :: poserr_ipc,poserr_maxi
      real*4,dimension(:),allocatable :: flux_ipc,Ferr_ipc,fluxL_ipc,fluxU_ipc,frequency_ipc
      real*4,dimension(:,:),allocatable :: flux_maxi,Ferr_maxi,FluxL_maxi,FluxU_maxi,frequency_maxi

      integer*4,dimension(:),allocatable :: rtype_source,nrep,t,track,ttsource,bary,rank,priority
      real*8,dimension(:),allocatable :: ra_source,dec_source,ra_xx,dec_xx
      real*4,dimension(:),allocatable :: xxerr,poserr_source,flux_source,xflux,rflux,rrconst,zsource
      character*30,dimension(:),allocatable :: nnsource

      integer*4,dimension(:),allocatable :: track2,bigbind
      real*8,dimension(:),allocatable :: ra_gam,dec_gam,ra_cat,dec_cat
      character*30,dimension(:),allocatable :: name_cat,namegam
      real*8,dimension(:),allocatable :: ra_other,dec_other
      real*4,dimension(:),allocatable :: savemjy,zz
      CHARACTER*30,dimension(:),allocatable :: name_other,vlasssrnm

      LOGICAL there,ok,found,catsrc
      common webprograms
      ok = .TRUE. 
      found = .FALSE.
      catsrc=.false.

      iconfig=0
      iarr=0
      arrsize(1:30)=0
      ifound = 0
      sfound = 0
      rfound = 0
      ncat=0
      iradio=0
      ixmm=0
      irosat=0
      iswift=0
      iipc=0
      ibmw = 0
      iother=0
      ichandra=0
      imaxi=0
      igam=0
      ibigb=0
      igrb=0
      ierosita=0
      radian = 45.0/atan(1.0)
      flux2nufnu_nvss=1.4e9*1.e-26
      flux2nufnu_sumss=8.43e8*1.e-26 !assumed radio alpha=0.2 !f_0.8 to f_1.4
      flux2nufnu_vlass=3.e9*1.e-26
c approximate flux conversions from cts/s to erg/cm2/s at 1 kev (NH=5.e20)
      flux2nufnu_rosat=7.e-12!!!!
      !flux2nufnu_xmm=7.e-13
      flux2nufnu_swift=9.e-12 !!!!
      flux2nufnu_ipc=1.2e-11
      flux2nufnu_bmw=1.8e-11
      !flux2nufnu_rxs=1.0
      sign=' '
      min_dist_ipc=50./3600.
      min_dist_rosat=40./3600.
c 40 arcsecs
      min_dist_xmm=15./3600.
      !min_dist_xmmsl=15./3600.
      min_dist_swift=7./3600.
      min_dist_bmw=10./3600.
      min_dist_chandra=5./3600.
c 10 arcsecs
      min_dist_other=15./3600. !check 30 extended CRATES no radio, check 30 instand.
      min_dist_cluster=60./3600.
c 15 arcsecs
      mjdavg=55000.

      CALL rdforn(string,length)
      IF ( length.NE.0 ) THEN
         CALL rmvlbk(string)
c         write(*,*) string
c         write(*,*) length
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
         in=im
         im=index(string(in+1:length),' ')+in
         output_file4=string(in+1:im-1)
c         in=im
c         im=index(string(in+1:length),' ')+in
c         output_file5=string(in+1:im-1)
         in=im
         im=index(string(in+1:length),' ')+in
         webprograms=string(in+1:im-1)
c         write(*,*) webprograms
         read(string(im+1:length),*) aim
c         write(*,*) 'the aim',aim
      ENDIF
      lu_in = 10
      !lu_output=11
      in = index(input_file(1:lenact(input_file)),'.')
      IF (in == 0) input_file(lenact(input_file)+1:lenact(input_file)+4) = '.csv' 
      INQUIRE (FILE=input_file,EXIST=there)
      IF (.NOT.there) THEN
         write (*,'('' file '',a,'' not found. No data found in phase 1 '')')
     &     input_file(1:lenact(input_file))
         STOP
      ENDIF

      array_size=webprograms(1:lenact(webprograms))//'/array_size.cf'
      open(lu_in,file=input_file,status='old',iostat=ier)
      open(11,file=output_file,status='unknown',iostat=ier)
      open(13,file=output_file2,status='unknown',iostat=ier)
      open(14,file=output_file4,status='unknown',iostat=ier)
      open(12,file=output_file3,status='unknown',iostat=ier)
      open(18,file=array_size,status='old',iostat=ier)
      IF (ier.NE.0) THEN
        write (*,*) ' Error ',ier,' opening file ', input_file
      ENDIF

      do while(ok)
      read(18,'(a)',end=700) string
      if (string(1:3) == '---' ) then
         iconfig=iconfig+1
      else
         if (iconfig .eq. 1) then
            iarr=iarr+1
            if (iarr .le. 3) then
               is=index(string(1:len(string)),':')
               read(string(is+1:lenact(string)),*) arrsize(iarr)
            else
               iarr=iarr-1
            endif
         else if (iconfig .eq. 2) then
            iarr=iarr+1
            is=index(string(1:len(string)),':')
            read(string(is+1:lenact(string)),*) arrsize(iarr)
         endif
      endif
      enddo
700   continue
      close(18)
c      write(*,*) arrsize

      allocate(nrep(arrsize(1)),zsource(arrsize(1)),nnsource(arrsize(1)))
      allocate(ra_xmm(arrsize(6)),dec_xmm(arrsize(6)),xmm_type(arrsize(6)),poserr_xmm(arrsize(6)))
      allocate(ra_swift(arrsize(7)),dec_swift(arrsize(7)),xrt_type(arrsize(7)),poserr_swift(arrsize(7)),mjdst_swift(arrsize(7)),mjded_swift(arrsize(7)))
      allocate(Ferr_swift(arrsize(7),5),FluxU_swift(arrsize(7),5),FluxL_swift(arrsize(7),5),flux_swift(arrsize(7),5),frequency_swift(arrsize(7),5))
      allocate(Ferr_xmm(arrsize(6),6),FluxU_xmm(arrsize(6),6),FluxL_xmm(arrsize(6),6),flux_xmm(arrsize(6),6),frequency_xmm(arrsize(6),6))
      allocate(ra_rosat(arrsize(8)),dec_rosat(arrsize(8)),poserr_rosat(arrsize(8)),rosat_type(arrsize(8)))
      allocate(poserr_chandra(arrsize(9)),ra_chandra(arrsize(9)),dec_chandra(arrsize(9)))
      allocate(flux_rosat(arrsize(8)),Ferr_rosat(arrsize(8)),FluxU_rosat(arrsize(8)),FluxL_rosat(arrsize(8)),frequency_rosat(arrsize(8)))
      allocate(flux_chandra(arrsize(9),5),FluxU_chandra(arrsize(9),5),FluxL_chandra(arrsize(9),5),frequency_chandra(arrsize(9),5))
      allocate(poserr_erosita(arrsize(10)),ra_erosita(arrsize(10)),dec_erosita(arrsize(10)))
      allocate(poserr_bmw(arrsize(11)),ra_bmw(arrsize(11)),dec_bmw(arrsize(11)))
      allocate(flux_erosita(arrsize(10),8),fluxL_erosita(arrsize(10),8),fluxU_erosita(arrsize(10),8),Ferr_erosita(arrsize(10),8),frequency_erosita(arrsize(10),8))
      allocate(flux_bmw(arrsize(11)),fluxL_bmw(arrsize(11)),fluxU_bmw(arrsize(11)),Ferr_bmw(arrsize(11)),frequency_bmw(arrsize(11)))
      allocate(ipc_type(arrsize(12)),poserr_ipc(arrsize(12)),ra_ipc(arrsize(12)),dec_ipc(arrsize(12)))
      allocate(poserr_maxi(arrsize(13)),ra_maxi(arrsize(13)),dec_maxi(arrsize(13)),maxi_type(arrsize(13)))
      allocate(flux_ipc(arrsize(12)),Ferr_ipc(arrsize(12)),fluxL_ipc(arrsize(12)),fluxU_ipc(arrsize(12)),frequency_ipc(arrsize(12)))
      allocate(flux_maxi(arrsize(13),4),Ferr_maxi(arrsize(13),4),FluxL_maxi(arrsize(13),4),FluxU_maxi(arrsize(13),4),frequency_maxi(arrsize(13),4))
      allocate(ra_gam(arrsize(5)),dec_gam(arrsize(5)),ra_cat(arrsize(5)),dec_cat(arrsize(5)))
      allocate(name_cat(arrsize(5)),namegam(arrsize(5)),track2(arrsize(5)),bigbind(arrsize(5)))
      allocate(ra_other(arrsize(4)),dec_other(arrsize(4)),zz(arrsize(4)))
      allocate(name_other(arrsize(4)),vlasssrnm(arrsize(3)))
      allocate(ra_radio(arrsize(3)),dec_radio(arrsize(3)),radio_type(arrsize(3)),ra_index(arrsize(3)))
      allocate(ppss(arrsize(3)),const(arrsize(3)),Ferr_radio(arrsize(3)),FluxU_radio(arrsize(3)),FluxL_radio(arrsize(3)),poserr_radio(arrsize(3)),flux_radio(arrsize(3)),frequency_radio(arrsize(3)))

      nrep(1:arrsize(1))=1
      zsource(1:arrsize(1))=0.
      nnsource(1:arrsize(1))='NONAME'

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
      is = index(string(ie+1:len(string)),'=') +ie
      read(string(is+1:len(string)),*) errrad,errmaj,errmin,errang
      !write(*,*) nh,errrad,errmaj,errmin,errang

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
         IF ( (catalog(1:4) == 'nvss') .OR. (catalog(1:5) == 'first') .or.
     &        (catalog(1:7) == 'vlassql') .OR. (catalog(1:5) == 'sumss') ) THEN
            iradio=iradio+1
            !write(*,*) "Nr. radio",iradio
            IF (iradio > arrsize(3)) Stop 'Too many NVSS/SUMSS points'
            ra_radio(iradio)=ra
            dec_radio(iradio)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            IF (catalog(1:4) == 'nvss') radio_type(iradio) = 3
            IF (catalog(1:5) == 'first') radio_type(iradio) = 2
            IF (catalog(1:5) == 'sumss') radio_type(iradio) = 4
            IF (catalog(1:7) == 'vlassql') radio_type(iradio) = 1
            IF ((catalog(1:5) == 'sumss') .or. (catalog(1:4) == 'nvss')
     &               .or. (catalog(1:7) == 'vlassql')) THEN
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_radio(iradio)
               poserr_radio(iradio)=poserr_radio(iradio)*2.*(1./0.95) !95% error
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_radio(iradio)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               !if (catalog(1:7) == 'vlassql') ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_radio(iradio)
               if (catalog(1:4) == 'nvss') then
                  erraxis=0.
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) major
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) erraxis
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) minor
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if ((is .ne. ie-1) .and. (erraxis .ne. 0.)) then
                     read(string(is+1:ie-1),*) erraxis
                  else
                     erraxis=0.
                  endif
                  is=ie
                  ie=index(string(is+1:len(string)),' ')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) posang
                  posxerr=sqrt(((sin(posang)*major)**2)+((cos(posang)*minor)**2))
                  posyerr=sqrt(((cos(posang)*major)**2)+((sin(posang)*minor)**2))
c                  write(*,*) ra_radio(iradio),dec_radio(iradio),major,minor,erraxis
                  if (erraxis .ne. 0.) poserr_radio(iradio)=max(posxerr,posyerr)
               else if (catalog(1:5) == 'sumss') then
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) major
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) minor
                  is=ie
                  ie=index(string(is+1:len(string)),' ')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) posang
                  posxerr=sqrt(((sin(posang)*major)**2)+((cos(posang)*minor)**2))
                  posyerr=sqrt(((cos(posang)*major)**2)+((sin(posang)*minor)**2))
c                  poserr_radio(iradio)=max(posxerr,posyerr)
               else if (catalog(1:7) == 'vlassql') then
                  is=ie
                  ie=index(string(is+1:len(string)),' ')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) vlasssrnm(iradio)
                  !write(*,*) iradio,catalog,vlasssrnm(iradio)(1:4)
               endif
            ELSE
               if (is .ne. ie-1) read(string(is+1:ie-1),*) ppss(iradio)
               !ppss(iradio)=0.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_radio(iradio)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_radio(iradio)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) major
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) minor
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) posang
               posxerr=sqrt(((sin(posang)*major)**2)+((cos(posang)*minor)**2))
               posyerr=sqrt(((cos(posang)*major)**2)+((sin(posang)*minor)**2))
               poserr_radio(iradio)=max(posxerr,posyerr)
               if (poserr_radio(iradio) .lt. 0.1 ) poserr_radio(iradio)=0.1
c               write(*,*) 'large',poserr_radio(iradio)
c               poserr_radio(iradio)=major*((Ferr_radio(iradio)/(flux_radio(iradio)-0.25))+0.05 )
c               if (poserr_radio(iradio) .lt. 0.1 ) poserr_radio(iradio)=0.1
c               write(*,*) 'act',poserr_radio(iradio)
            ENDIF
            IF (catalog(1:5) == 'sumss') then
               const(iradio) = flux2nufnu_sumss
               frequency_radio(iradio)=8.43E8
            else if (catalog(1:7) == 'vlassql') then
               const(iradio) = flux2nufnu_vlass
               frequency_radio(iradio)=3.e9
            else
               const(iradio) = flux2nufnu_nvss
               frequency_radio(iradio)=1.4e9
            endif
            !write(*,*) flux_radio(iradio),Ferr_radio(iradio),catalog
            FluxU_radio(iradio)=(flux_radio(iradio)+Ferr_radio(iradio))*const(iradio)
            FluxL_radio(iradio)=(flux_radio(iradio)-Ferr_radio(iradio))*const(iradio)
c PG
            CALL RXgraphic_code(flux_radio(iradio),'R',code)
            write (13,'(f9.5,2x,f9.5,2x,i6)') ra_radio(iradio),dec_radio(iradio),int(code)
c PG end
            flux_radio(iradio)=flux_radio(iradio)*const(iradio)
            !write(*,*) iradio,radio_type(iradio)
            if (radio_type(iradio) .eq. 2) then
               if (ppss(iradio) .gt. 0.15) iradio=iradio-1
            else if (radio_type(iradio) .eq. 1) then
               if (vlasssrnm(iradio)(1:1) .ne. 'J') iradio=iradio-1
            endif
c            if (iradio .ge. 1) then
c               write(*,*) iradio,catalog,radio_type(iradio),flux_radio(iradio)/const(iradio)
c            endif
            !if (iradio .eq. 14) iradio=iradio-1
         ELSE IF ( (catalog(1:5) == 'xmmsl') .OR.
     &             (catalog(1:4) == '4xmm') )  THEN
            ixmm=ixmm+1
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            ra_xmm(ixmm)=ra
            dec_xmm(ixmm)=dec
            IF (ixmm > arrsize(6)) Stop 'Too many XMM points'
            !write(*,*) FluxU_xmm(ixmm,1),flux_xmm(ixmm,1),FluxL_xmm(ixmm,1)
            IF (catalog(1:5) == 'xmmsl') THEN
               xmm_type(ixmm)=1
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xmm(ixmm,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xmm(ixmm,2)
               FluxU_xmm(ixmm,2)=flux_xmm(ixmm,2)+Ferr_xmm(ixmm,2)
               FluxL_xmm(ixmm,2)=flux_xmm(ixmm,2)-Ferr_xmm(ixmm,2)
               if (Ferr_xmm(ixmm,2) .gt. flux_xmm(ixmm,2)) then
                  FluxU_xmm(ixmm,2)=0.!Ferr_xmm(ixmm,1)*3.
                  FluxL_xmm(ixmm,2)=0.
                  flux_xmm(ixmm,2)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xmm(ixmm,1) !6
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xmm(ixmm,1)
               FluxU_xmm(ixmm,1)=flux_xmm(ixmm,1)+Ferr_xmm(ixmm,1)
               FluxL_xmm(ixmm,1)=flux_xmm(ixmm,1)-Ferr_xmm(ixmm,1)
               if (Ferr_xmm(ixmm,1) .gt. flux_xmm(ixmm,1)) then
                  FluxU_xmm(ixmm,1)=0.!Ferr_xmm(ixmm,2)*3.
                  FluxL_xmm(ixmm,1)=0.
                  Flux_xmm(ixmm,1)=0.
               endif
               !write(*,*) FluxU_xmm(ixmm,2),flux_xmm(ixmm,2),FluxL_xmm(ixmm,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xmm(ixmm,3) !7
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xmm(ixmm,3)
               FluxU_xmm(ixmm,3)=flux_xmm(ixmm,3)+Ferr_xmm(ixmm,3)
               FluxL_xmm(ixmm,3)=flux_xmm(ixmm,3)-Ferr_xmm(ixmm,3)
               if (Ferr_xmm(ixmm,3) .gt. flux_xmm(ixmm,3)) then
                  FluxU_xmm(ixmm,3)=0.!Ferr_xmm(ixmm,3)*3.
                  FluxL_xmm(ixmm,3)=0.
                  Flux_xmm(ixmm,3)=0.
               endif
               !write(*,*) FluxU_xmm(ixmm,3),flux_xmm(ixmm,3),FluxL_xmm(ixmm,3)
               call nhdeabsorb2 (1,0.2,2.,0.9,nh,reduce,100)
               flux_xmm(ixmm,1)=flux_xmm(ixmm,1)*reduce
               FluxU_xmm(ixmm,1)=FluxU_xmm(ixmm,1)*reduce
               FluxL_xmm(ixmm,1)=FluxL_xmm(ixmm,1)*reduce
               call fluxtofdens(0.9,0.2,2.,flux_xmm(ixmm,1),1.,fdens,nudens)
               flux_xmm(ixmm,1)=fdens
               frequency_xmm(ixmm,1)=nudens
               call fluxtofdens(0.9,0.2,2.,FluxU_xmm(ixmm,1),1.,fdens,nudens)
               FluxU_xmm(ixmm,1)=fdens
               call fluxtofdens(0.9,0.2,2.,FluxL_xmm(ixmm,1),1.,fdens,nudens)
               FluxL_xmm(ixmm,1)=fdens
               if (flux_xmm(ixmm,1) .eq. 0.) then
                  if (flux_xmm(ixmm,2) .ne. 0) then
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,2)
                     Ferr_xmm(ixmm,1)=Ferr_xmm(ixmm,2)
                     FluxU_xmm(ixmm,1)=flux_xmm(ixmm,1)+Ferr_xmm(ixmm,1)
                     FluxL_xmm(ixmm,1)=flux_xmm(ixmm,1)-Ferr_xmm(ixmm,1)
                     call nhdeabsorb2 (1,0.2,12.,0.9,nh,reduce,100)
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,1)*reduce
                     FluxU_xmm(ixmm,1)=FluxU_xmm(ixmm,1)*reduce
                     FluxL_xmm(ixmm,1)=FluxL_xmm(ixmm,1)*reduce
                     call fluxtofdens(0.9,0.2,12.,flux_xmm(ixmm,1),1.0,fdens,nudens)
                     flux_xmm(ixmm,1)=fdens
                     frequency_xmm(ixmm,1)=nudens
                     call fluxtofdens(0.9,0.2,12.,FluxU_xmm(ixmm,1),1.0,fdens,nudens)
                     FluxU_xmm(ixmm,1)=fdens
                     call fluxtofdens(0.9,0.2,12.,FluxL_xmm(ixmm,1),1.0,fdens,nudens)
                     FluxL_xmm(ixmm,1)=fdens
                  ELSE if (flux_xmm(ixmm,3) .ne. 0) then
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,3)
                     Ferr_xmm(ixmm,1)=Ferr_xmm(ixmm,3)
                     FluxU_xmm(ixmm,1)=flux_xmm(ixmm,1)+Ferr_xmm(ixmm,1)
                     FluxL_xmm(ixmm,1)=flux_xmm(ixmm,1)-Ferr_xmm(ixmm,1)
                     call nhdeabsorb2 (1,2.,12.,0.9,nh,reduce,100)
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,1)*reduce
                     FluxU_xmm(ixmm,1)=FluxU_xmm(ixmm,1)*reduce
                     FluxL_xmm(ixmm,1)=FluxL_xmm(ixmm,1)*reduce
                     call fluxtofdens(0.9,2.,12.,flux_xmm(ixmm,1),1.,fdens,nudens)
                     flux_xmm(ixmm,1)=fdens
                     frequency_xmm(ixmm,1)=nudens
                     call fluxtofdens(0.9,2.,12.,FluxU_xmm(ixmm,1),1.,fdens,nudens)
                     FluxU_xmm(ixmm,1)=fdens
                     call fluxtofdens(0.9,2.,12.,FluxL_xmm(ixmm,1),1.,fdens,nudens)
                     FluxL_xmm(ixmm,1)=fdens
                  endif
               endif
               call nhdeabsorb2 (1,0.2,12.,0.9,nh,reduce,100)
               !write(*,*) reduce
               flux_xmm(ixmm,2)=flux_xmm(ixmm,2)*reduce
               FluxU_xmm(ixmm,2)=FluxU_xmm(ixmm,2)*reduce
               FluxL_xmm(ixmm,2)=FluxL_xmm(ixmm,2)*reduce
               call fluxtofdens(0.9,0.2,12.,flux_xmm(ixmm,2),3.,fdens,nudens)
               flux_xmm(ixmm,2)=fdens
               frequency_xmm(ixmm,2)=nudens
               call fluxtofdens(0.9,0.2,12.,FluxU_xmm(ixmm,2),3.,fdens,nudens)
               FluxU_xmm(ixmm,2)=fdens
               call fluxtofdens(0.9,0.2,12.,FluxL_xmm(ixmm,2),3.,fdens,nudens)
               FluxL_xmm(ixmm,2)=fdens
               call nhdeabsorb2 (1,2.,12.,0.9,nh,reduce,100)
               flux_xmm(ixmm,3)=flux_xmm(ixmm,3)*reduce
               FluxU_xmm(ixmm,3)=FluxU_xmm(ixmm,3)*reduce
               FluxL_xmm(ixmm,3)=FluxL_xmm(ixmm,3)*reduce
               call fluxtofdens(0.9,2.,12.,flux_xmm(ixmm,3),5.,fdens,nudens)
               flux_xmm(ixmm,3)=fdens
               frequency_xmm(ixmm,3)=nudens
               call fluxtofdens(0.9,2.,12.,FluxU_xmm(ixmm,3),5.,fdens,nudens)
               FluxU_xmm(ixmm,3)=fdens
               call fluxtofdens(0.9,2.,12.,FluxL_xmm(ixmm,3),5.,fdens,nudens)
               FluxL_xmm(ixmm,3)=fdens
            ELSE IF (catalog(1:4) == '4xmm') THEN
               xmm_type(ixmm)=2
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xmm(ixmm,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xmm(ixmm,1)
               FluxU_xmm(ixmm,1)=flux_xmm(ixmm,1)+Ferr_xmm(ixmm,1)
               FluxL_xmm(ixmm,1)=flux_xmm(ixmm,1)-Ferr_xmm(ixmm,1)
               if (Ferr_xmm(ixmm,1) .gt. flux_xmm(ixmm,1)) then
                  FluxU_xmm(ixmm,1)=0.!Ferr_xmm(ixmm,1)*3.
                  FluxL_xmm(ixmm,1)=0.
                  flux_xmm(ixmm,1)=0.
               endif
               call nhdeabsorb2 (1,0.2,12.,0.9,nh,reduce,100)
               flux_xmm(ixmm,1)=flux_xmm(ixmm,1)*reduce
               FluxU_xmm(ixmm,1)=FluxU_xmm(ixmm,1)*reduce
               FluxL_xmm(ixmm,1)=FluxL_xmm(ixmm,1)*reduce
               call fluxtofdens(0.9,0.2,12.,flux_xmm(ixmm,1),1.,fdens,nudens)
               flux_xmm(ixmm,1)=fdens
               frequency_xmm(ixmm,1)=nudens
               call fluxtofdens(0.9,0.2,12.,FluxU_xmm(ixmm,1),1.,fdens,nudens)
               FluxU_xmm(ixmm,1)=fdens
               call fluxtofdens(0.9,0.2,12.,FluxL_xmm(ixmm,1),1.,fdens,nudens)
               FluxL_xmm(ixmm,1)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xmm(ixmm,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xmm(ixmm,2)
               FluxU_xmm(ixmm,2)=flux_xmm(ixmm,2)+Ferr_xmm(ixmm,2)
               FluxL_xmm(ixmm,2)=flux_xmm(ixmm,2)-Ferr_xmm(ixmm,2)
               if (Ferr_xmm(ixmm,2) .gt. flux_xmm(ixmm,2)) then
                  FluxU_xmm(ixmm,2)=0.!Ferr_xmm(ixmm,2)*3.
                  FluxL_xmm(ixmm,2)=0.
                  flux_xmm(ixmm,2)=0.
               endif
               !write(*,*) FluxU_xmm(ixmm,2),flux_xmm(ixmm,2),FluxL_xmm(ixmm,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xmm(ixmm,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xmm(ixmm,3)
               FluxU_xmm(ixmm,3)=flux_xmm(ixmm,3)+Ferr_xmm(ixmm,3)
               FluxL_xmm(ixmm,3)=flux_xmm(ixmm,3)-Ferr_xmm(ixmm,3)
               if (Ferr_xmm(ixmm,3) .gt. flux_xmm(ixmm,3)) then
                  FluxU_xmm(ixmm,3)=0.!Ferr_xmm(ixmm,3)*3.
                  FluxL_xmm(ixmm,3)=0.
                  flux_xmm(ixmm,3)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xmm(ixmm,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xmm(ixmm,4)
               FluxU_xmm(ixmm,4)=flux_xmm(ixmm,4)+Ferr_xmm(ixmm,4)
               FluxL_xmm(ixmm,4)=flux_xmm(ixmm,4)-Ferr_xmm(ixmm,4)
               if (Ferr_xmm(ixmm,4) .gt. flux_xmm(ixmm,4)) then
                  FluxU_xmm(ixmm,4)=0.!Ferr_xmm(ixmm,4)*3.
                  FluxL_xmm(ixmm,4)=0.
                  flux_xmm(ixmm,4)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xmm(ixmm,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xmm(ixmm,5)
               FluxU_xmm(ixmm,5)=flux_xmm(ixmm,5)+Ferr_xmm(ixmm,5)
               FluxL_xmm(ixmm,5)=flux_xmm(ixmm,5)-Ferr_xmm(ixmm,5)
               if (Ferr_xmm(ixmm,5) .gt. flux_xmm(ixmm,5)) then
                  FluxU_xmm(ixmm,5)=0.!Ferr_xmm(ixmm,5)*3.
                  FluxL_xmm(ixmm,5)=0.
                  flux_xmm(ixmm,5)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xmm(ixmm,6)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xmm(ixmm,6)
               FluxU_xmm(ixmm,6)=flux_xmm(ixmm,6)+Ferr_xmm(ixmm,6)
               FluxL_xmm(ixmm,6)=flux_xmm(ixmm,6)-Ferr_xmm(ixmm,6)
               if (Ferr_xmm(ixmm,6) .gt. flux_xmm(ixmm,6)) then
                  FluxU_xmm(ixmm,6)=0.!Ferr_xmm(ixmm,6)*3.
                  FluxL_xmm(ixmm,6)=0.
                  flux_xmm(ixmm,6)=0.
               endif
               !write(*,*) FluxU_xmm(ixmm,6),flux_xmm(ixmm,6),FluxL_xmm(ixmm,6)
               if (flux_xmm(ixmm,1) .eq. 0.) then
                  if (flux_xmm(ixmm,4) .ne. 0.) then
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,4)
                     Ferr_xmm(ixmm,1)=Ferr_xmm(ixmm,4)
                     FluxU_xmm(ixmm,1)=flux_xmm(ixmm,1)+Ferr_xmm(ixmm,1)
                     FluxL_xmm(ixmm,1)=flux_xmm(ixmm,1)-Ferr_xmm(ixmm,1)
                     call nhdeabsorb2 (1,1.,2.,0.9,nh,reduce,100)
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,1)*reduce
                     FluxU_xmm(ixmm,1)=FluxU_xmm(ixmm,1)*reduce
                     FluxL_xmm(ixmm,1)=FluxL_xmm(ixmm,1)*reduce
                     call fluxtofdens(0.9,1.,2.,flux_xmm(ixmm,1),1.,fdens,nudens)
                     flux_xmm(ixmm,1)=fdens
                     call fluxtofdens(0.9,1.,2.,FluxU_xmm(ixmm,1),1.,fdens,nudens)
                     FluxU_xmm(ixmm,1)=fdens
                     call fluxtofdens(0.9,1.,2.,FluxL_xmm(ixmm,1),1.,fdens,nudens)
                     FluxL_xmm(ixmm,1)=fdens
                  else if (flux_xmm(ixmm,3) .ne. 0) then
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,3)
                     Ferr_xmm(ixmm,1)=Ferr_xmm(ixmm,3)
                     FluxU_xmm(ixmm,1)=flux_xmm(ixmm,1)+Ferr_xmm(ixmm,1)
                     FluxL_xmm(ixmm,1)=flux_xmm(ixmm,1)-Ferr_xmm(ixmm,1)
                     call nhdeabsorb2 (1,0.5,1.,0.9,nh,reduce,100)
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,1)*reduce
                     FluxU_xmm(ixmm,1)=FluxU_xmm(ixmm,1)*reduce
                     FluxL_xmm(ixmm,1)=FluxL_xmm(ixmm,1)*reduce
                     call fluxtofdens(0.9,0.5,1.,flux_xmm(ixmm,1),1.,fdens,nudens)
                     flux_xmm(ixmm,1)=fdens
                     frequency_xmm(ixmm,1)=nudens
                     call fluxtofdens(0.9,0.5,1.,FluxU_xmm(ixmm,1),1.,fdens,nudens)
                     FluxU_xmm(ixmm,1)=fdens
                     call fluxtofdens(0.9,0.5,1.,FluxL_xmm(ixmm,1),1.,fdens,nudens)
                     FluxL_xmm(ixmm,1)=fdens
                  else if (flux_xmm(ixmm,2) .ne. 0) then
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,2)
                     Ferr_xmm(ixmm,1)=Ferr_xmm(ixmm,2)
                     FluxU_xmm(ixmm,1)=flux_xmm(ixmm,1)+Ferr_xmm(ixmm,1)
                     FluxL_xmm(ixmm,1)=flux_xmm(ixmm,1)-Ferr_xmm(ixmm,1)
                     call nhdeabsorb2 (1,0.2,0.5,0.9,nh,reduce,100)
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,1)*reduce
                     FluxU_xmm(ixmm,1)=FluxU_xmm(ixmm,1)*reduce
                     FluxL_xmm(ixmm,1)=FluxL_xmm(ixmm,1)*reduce
                     call fluxtofdens(0.9,0.2,0.5,flux_xmm(ixmm,2),1.,fdens,nudens)
                     frequency_xmm(ixmm,1)=nudens
                     call fluxtofdens(0.9,0.2,0.5,FluxU_xmm(ixmm,2),1.,fdens,nudens)
                     FluxU_xmm(ixmm,1)=fdens
                     call fluxtofdens(0.9,0.2,0.5,FluxL_xmm(ixmm,2),1.,fdens,nudens)
                     FluxL_xmm(ixmm,1)=fdens
                  else if (flux_xmm(ixmm,5) .ne. 0) then
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,5)
                     Ferr_xmm(ixmm,1)=Ferr_xmm(ixmm,5)
                     FluxU_xmm(ixmm,1)=flux_xmm(ixmm,1)+Ferr_xmm(ixmm,1)
                     FluxL_xmm(ixmm,1)=flux_xmm(ixmm,1)-Ferr_xmm(ixmm,1)
                     call nhdeabsorb2 (1,2.,4.5,0.9,nh,reduce,100)
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,1)*reduce
                     FluxU_xmm(ixmm,1)=FluxU_xmm(ixmm,1)*reduce
                     FluxL_xmm(ixmm,1)=FluxL_xmm(ixmm,1)*reduce
                     call fluxtofdens(0.9,2.,4.5,flux_xmm(ixmm,1),1.,fdens,nudens)
                     frequency_xmm(ixmm,1)=nudens
                     call fluxtofdens(0.9,2.,4.5,FluxU_xmm(ixmm,1),1.,fdens,nudens)
                     FluxU_xmm(ixmm,1)=fdens
                     call fluxtofdens(0.9,2.,4.5,FluxL_xmm(ixmm,1),1.,fdens,nudens)
                     FluxL_xmm(ixmm,1)=fdens
                  else if (flux_xmm(ixmm,6) .ne. 0) then
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,6)
                     Ferr_xmm(ixmm,1)=Ferr_xmm(ixmm,6)
                     FluxU_xmm(ixmm,1)=flux_xmm(ixmm,1)+Ferr_xmm(ixmm,1)
                     FluxL_xmm(ixmm,1)=flux_xmm(ixmm,1)-Ferr_xmm(ixmm,1)
                     call nhdeabsorb2 (1,4.5,12.,0.9,nh,reduce,100)
                     !write(*,*) reduce
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,1)*reduce
                     FluxU_xmm(ixmm,1)=FluxU_xmm(ixmm,1)*reduce
                     FluxL_xmm(ixmm,1)=FluxL_xmm(ixmm,1)*reduce
                     call fluxtofdens(0.9,4.5,12.,flux_xmm(ixmm,1),1.,fdens,nudens)
                     flux_xmm(ixmm,1)=fdens
                     call fluxtofdens(0.9,4.5,12.,FluxU_xmm(ixmm,1),1.,fdens,nudens)
                     FluxU_xmm(ixmm,1)=fdens
                     call fluxtofdens(0.9,4.5,12.,FluxL_xmm(ixmm,1),1.,fdens,nudens)
                     FluxL_xmm(ixmm,1)=fdens
                  endif
               endif
               call nhdeabsorb2 (1,0.2,0.5,0.9,nh,reduce,100)
               flux_xmm(ixmm,2)=flux_xmm(ixmm,2)*reduce
               FluxU_xmm(ixmm,2)=FluxU_xmm(ixmm,2)*reduce
               FluxL_xmm(ixmm,2)=FluxL_xmm(ixmm,2)*reduce
               call fluxtofdens(0.9,0.2,0.5,flux_xmm(ixmm,2),0.3,fdens,nudens)
               flux_xmm(ixmm,2)=fdens
               frequency_xmm(ixmm,2)=nudens
               call fluxtofdens(0.9,0.2,0.5,FluxU_xmm(ixmm,2),0.3,fdens,nudens)
               FluxU_xmm(ixmm,2)=fdens
               call fluxtofdens(0.9,0.2,0.5,FluxL_xmm(ixmm,2),0.3,fdens,nudens)
               FluxL_xmm(ixmm,2)=fdens
               call nhdeabsorb2 (1,0.5,1.,0.9,nh,reduce,100)
               flux_xmm(ixmm,3)=flux_xmm(ixmm,3)*reduce
               FluxU_xmm(ixmm,3)=FluxU_xmm(ixmm,3)*reduce
               FluxL_xmm(ixmm,3)=FluxL_xmm(ixmm,3)*reduce
               call fluxtofdens(0.9,0.5,1.,flux_xmm(ixmm,3),0.75,fdens,nudens)
               flux_xmm(ixmm,3)=fdens
               frequency_xmm(ixmm,3)=nudens
               call fluxtofdens(0.9,0.5,1.,FluxU_xmm(ixmm,3),0.75,fdens,nudens)
               FluxU_xmm(ixmm,3)=fdens
               call fluxtofdens(0.9,0.5,1.,FluxL_xmm(ixmm,3),0.75,fdens,nudens)
               FluxL_xmm(ixmm,3)=fdens
               call nhdeabsorb2 (1,1.,2.,0.9,nh,reduce,100)
               flux_xmm(ixmm,4)=flux_xmm(ixmm,4)*reduce
               FluxU_xmm(ixmm,4)=FluxU_xmm(ixmm,4)*reduce
               FluxL_xmm(ixmm,4)=FluxL_xmm(ixmm,4)*reduce
               call fluxtofdens(0.9,1.,2.,flux_xmm(ixmm,4),1.5,fdens,nudens)
               flux_xmm(ixmm,4)=fdens
               frequency_xmm(ixmm,4)=nudens
               call fluxtofdens(0.9,1.,2.,FluxU_xmm(ixmm,4),1.5,fdens,nudens)
               FluxU_xmm(ixmm,4)=fdens
               call fluxtofdens(0.9,1.,2.,FluxL_xmm(ixmm,4),1.5,fdens,nudens)
               FluxL_xmm(ixmm,4)=fdens
               call nhdeabsorb2 (1,2.,4.5,0.9,nh,reduce,100)
               flux_xmm(ixmm,5)=flux_xmm(ixmm,5)*reduce
               FluxU_xmm(ixmm,5)=FluxU_xmm(ixmm,5)*reduce
               FluxL_xmm(ixmm,5)=FluxL_xmm(ixmm,5)*reduce
               call fluxtofdens(0.9,2.,4.5,flux_xmm(ixmm,5),3.,fdens,nudens)
               flux_xmm(ixmm,5)=fdens
               frequency_xmm(ixmm,5)=nudens
               call fluxtofdens(0.9,2.,4.5,FluxU_xmm(ixmm,5),3.,fdens,nudens)
               FluxU_xmm(ixmm,5)=fdens
               call fluxtofdens(0.9,2.,4.5,FluxL_xmm(ixmm,5),3.,fdens,nudens)
               FluxL_xmm(ixmm,5)=fdens
               call nhdeabsorb2 (1,4.5,12.,0.9,nh,reduce,100)
               !write(*,*) reduce
               flux_xmm(ixmm,6)=flux_xmm(ixmm,6)*reduce
               FluxU_xmm(ixmm,6)=FluxU_xmm(ixmm,6)*reduce
               FluxL_xmm(ixmm,6)=FluxL_xmm(ixmm,6)*reduce
               call fluxtofdens(0.9,4.5,12.,flux_xmm(ixmm,6),7.,fdens,nudens)
               flux_xmm(ixmm,6)=fdens
               call fluxtofdens(0.9,4.5,12.,FluxU_xmm(ixmm,6),7.,fdens,nudens)
               FluxU_xmm(ixmm,6)=fdens
               call fluxtofdens(0.9,4.5,12.,FluxL_xmm(ixmm,6),7.,fdens,nudens)
               FluxL_xmm(ixmm,6)=fdens
               frequency_xmm(ixmm,6)=nudens
            ENDIF
            is=ie
            ie=index(string(is+1:len(string)),' ')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_xmm(ixmm)
            if ((xmm_type(ixmm) == 1 )) poserr_xmm(ixmm)=sqrt((8**2)+(poserr_xmm(ixmm)**2))
            poserr_xmm(ixmm)=poserr_xmm(ixmm)*2
c PG
            CALL RXgraphic_code(flux_xmm(ixmm,1),'X',code)
            write (13,'(f9.5,2x,f9.5,2x,i6)') ra_xmm(ixmm),dec_xmm(ixmm),int(code)
c end PG
         ELSE IF ( (catalog(1:4) == 'rass') .OR.
     &             (catalog(1:3) == 'wga') ) THEN
            irosat=irosat+1
            ra_rosat(irosat)=ra
            dec_rosat(irosat)=dec
            IF (irosat > arrsize(8)) Stop 'Too many RASS points'
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_rosat(irosat)
            IF (catalog(1:4) == 'rass') THEN
               rosat_type(irosat)=1
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_rosat(irosat)
               if (flux_rosat(irosat) .gt. 0.1) then
                  poserr_rosat(irosat)=20.
               else
                  poserr_rosat(irosat)=40.
               endif
               FluxU_rosat(irosat)=flux_rosat(irosat)+Ferr_rosat(irosat)
               FluxL_rosat(irosat)=flux_rosat(irosat)-Ferr_rosat(irosat)
               call nhdeabsorb2 (0,0.1,2.4,0.9,nh,reduce,1)
               flux_rosat(irosat)=flux_rosat(irosat)*reduce
               FluxU_rosat(irosat)=FluxU_rosat(irosat)*reduce
               FluxL_rosat(irosat)=FluxL_rosat(irosat)*reduce
               call fluxtofdens(0.9,0.1,2.4,flux_rosat(irosat),1.,fdens,nudens)
               flux_rosat(irosat)=fdens
               frequency_rosat(irosat)=nudens
               call fluxtofdens(0.9,0.1,2.4,FluxU_rosat(irosat),1.,fdens,nudens)
               FluxU_rosat(irosat)=fdens
               call fluxtofdens(0.9,0.1,2.4,FluxL_rosat(irosat),1.,fdens,nudens)
               FluxL_rosat(irosat)=fdens
            ELSE IF (catalog(1:3) == 'wga') THEN
               rosat_type(irosat)=2
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_rosat(irosat)
               FluxU_rosat(irosat)=flux_rosat(irosat)+Ferr_rosat(irosat)
               FluxL_rosat(irosat)=flux_rosat(irosat)-Ferr_rosat(irosat)
               call nhdeabsorb2 (0,0.24,2.,0.9,nh,reduce,1)
               flux_rosat(irosat)=flux_rosat(irosat)*reduce
               FluxU_rosat(irosat)=FluxU_rosat(irosat)*reduce
               FluxL_rosat(irosat)=FluxL_rosat(irosat)*reduce
               call fluxtofdens(0.9,0.24,2.,flux_rosat(irosat),1.,fdens,nudens)
               flux_rosat(irosat)=fdens
               frequency_rosat(irosat)=nudens
               call fluxtofdens(0.9,0.24,2.,FluxU_rosat(irosat),1.,fdens,nudens)
               FluxU_rosat(irosat)=fdens
               call fluxtofdens(0.9,0.24,2.,FluxL_rosat(irosat),1.,fdens,nudens)
               FluxL_rosat(irosat)=fdens
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_rosat(irosat)
            ENDIF
            !write(*,*) catalog,FluxU_rosat(irosat),flux_rosat(irosat),FluxL_rosat(irosat),poserr_rosat(irosat)
c PG
            CALL RXgraphic_code(flux_rosat(irosat),'X',code)
            write (13,'(f9.5,2x,f9.5,2x,i6)') ra_rosat(irosat),dec_rosat(irosat),int(code)
c end PG
         ELSE IF ((catalog(1:5) == '2sxps') .or. (catalog(1:7) == 'xrtdeep')
     &           .or. (catalog(1:5) == 'sds82') .or. (catalog(1:5) == 'ousxb')
     &           .or.  (catalog(1:5) == 'ousxg') .or. (catalog(1:5) == '1ousx'))THEN
            iswift=iswift+1
            IF (iswift > arrsize(7)) Stop 'Too many swift points'
            ra_swift(iswift)=ra
            dec_swift(iswift)=dec
            frequency_swift(iswift,5)=999
            if (catalog(1:5) == '2sxps') then
               xrt_type(iswift)=1
               ra_swift(iswift)=-ra_swift(iswift)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_swift(iswift)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_swift(iswift,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,1)
               FluxU_swift(iswift,1)=flux_swift(iswift,1)+Ferr_swift(iswift,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,1)
               FluxL_swift(iswift,1)=flux_swift(iswift,1)+Ferr_swift(iswift,1)
               if (flux_swift(iswift,1) .ne. 0.) then
               if ((FluxL_swift(iswift,1) .le. 0.) .or.
     &             ((FluxU_swift(iswift,1)-flux_swift(iswift,1)) .gt. flux_swift(iswift,1))) then
                  flux_swift(iswift,1)=0.
                  FluxU_swift(iswift,1)=0.
                  FluxL_swift(iswift,1)=0.
               endif
               endif
               call nhdeabsorb2 (0,0.3,10.,0.9,nh,reduce,4)
               flux_swift(iswift,1)=flux_swift(iswift,1)*reduce
               FluxU_swift(iswift,1)=FluxU_swift(iswift,1)*reduce
               FluxL_swift(iswift,1)=FluxL_swift(iswift,1)*reduce
               call fluxtofdens(0.9,0.3,10.,flux_swift(iswift,1),3.,fdens,nudens)
               flux_swift(iswift,1)=fdens
               frequency_swift(iswift,1)=2.418e17 !nudens 3 keV
               call fluxtofdens(0.9,0.3,10.,FluxU_swift(iswift,1),3.,fdens,nudens)
               FluxU_swift(iswift,1)=fdens
               call fluxtofdens(0.9,0.3,10.,FluxL_swift(iswift,1),3.,fdens,nudens)
               FluxL_swift(iswift,1)=fdens
               if ((FluxU_swift(iswift,1) .ne. 0.) .and. (flux_swift(iswift,1) .eq. 0. )) then
                   FluxU_swift(iswift,1)=FluxU_swift(iswift,1)*3.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_swift(iswift,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_swift(iswift,2)
               FluxU_swift(iswift,2)=flux_swift(iswift,2)+Ferr_swift(iswift,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,2)
               FluxL_swift(iswift,2)=flux_swift(iswift,2)+Ferr_swift(iswift,2)
               if (flux_swift(iswift,2) .ne. 0.) then
               if ((FluxL_swift(iswift,2) .le. 0.) .or.
     &             ((FluxU_swift(iswift,2)-flux_swift(iswift,2)) .gt. flux_swift(iswift,2))) then
                  flux_swift(iswift,2)=0.
                  FluxU_swift(iswift,2)=0.
                  FluxL_swift(iswift,2)=0.
               endif
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_swift(iswift,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_swift(iswift,3)
               FluxU_swift(iswift,3)=flux_swift(iswift,3)+Ferr_swift(iswift,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,3)
               FluxL_swift(iswift,3)=flux_swift(iswift,3)+Ferr_swift(iswift,3)
               if (flux_swift(iswift,3) .ne. 0.) then
               if ((FluxL_swift(iswift,3) .le. 0.) .or.
     &             ((FluxU_swift(iswift,3)-flux_swift(iswift,3)) .gt. flux_swift(iswift,3))) then
                  flux_swift(iswift,3)=0.
                  FluxU_swift(iswift,3)=0.
                  FluxL_swift(iswift,3)=0.
               endif
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_swift(iswift,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_swift(iswift,4)
               FluxU_swift(iswift,4)=flux_swift(iswift,4)+Ferr_swift(iswift,4)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,4)
               FluxL_swift(iswift,4)=flux_swift(iswift,4)+Ferr_swift(iswift,4)
               if (flux_swift(iswift,4) .ne. 0.) then
               if ((FluxL_swift(iswift,4) .le. 0.) .or.
     &             ((FluxU_swift(iswift,4)-flux_swift(iswift,4)) .gt. flux_swift(iswift,4))) then
                  flux_swift(iswift,4)=0.
                  FluxU_swift(iswift,4)=0.
                  FluxL_swift(iswift,4)=0.
               endif
               endif
               if (flux_swift(iswift,1) .eq. 0.) then
                  if (flux_swift(iswift,3) .ne. 0.) then
                     flux_swift(iswift,1)=flux_swift(iswift,3)
                     Ferr_swift(iswift,1)=Ferr_swift(iswift,3)
                     FluxU_swift(iswift,1)=flux_swift(iswift,1)+Ferr_swift(iswift,1)
                     FluxL_swift(iswift,1)=flux_swift(iswift,1)-Ferr_swift(iswift,1)
                     call nhdeabsorb2 (0,1.,2.,0.9,nh,reduce,4)
                     flux_swift(iswift,1)=flux_swift(iswift,1)*reduce
                     FluxU_swift(iswift,1)=FluxU_swift(iswift,1)*reduce
                     FluxL_swift(iswift,1)=FluxL_swift(iswift,1)*reduce
                     call fluxtofdens(0.9,1.,2.,flux_swift(iswift,1),3.,fdens,nudens)
                     flux_swift(iswift,1)=fdens
                     frequency_swift(iswift,1)=2.418e17 !nudens 3 keV
                     call fluxtofdens(0.9,1.,2.,FluxU_swift(iswift,1),3.,fdens,nudens)
                     FluxU_swift(iswift,1)=fdens
                     call fluxtofdens(0.9,1.,2.,FluxL_swift(iswift,1),3.,fdens,nudens)
                     FluxL_swift(iswift,1)=fdens
                  else if (flux_swift(iswift,2) .ne. 0.) then
                     flux_swift(iswift,1)=flux_swift(iswift,2)
                     Ferr_swift(iswift,1)=Ferr_swift(iswift,2)
                     FluxU_swift(iswift,1)=flux_swift(iswift,1)+Ferr_swift(iswift,1)
                     FluxL_swift(iswift,1)=flux_swift(iswift,1)-Ferr_swift(iswift,1)
                     call nhdeabsorb2 (0,0.3,1.,0.9,nh,reduce,4)
                     flux_swift(iswift,1)=flux_swift(iswift,1)*reduce
                     FluxU_swift(iswift,1)=FluxU_swift(iswift,1)*reduce
                     FluxL_swift(iswift,1)=FluxL_swift(iswift,1)*reduce
                     call fluxtofdens(0.9,0.3,1.,flux_swift(iswift,1),3.,fdens,nudens)
                     flux_swift(iswift,1)=fdens
                     frequency_swift(iswift,1)=2.418e17 !nudens 3 keV
                     call fluxtofdens(0.9,0.3,1.,FluxU_swift(iswift,1),3.,fdens,nudens)
                     FluxU_swift(iswift,1)=fdens
                     call fluxtofdens(0.9,0.3,1.,FluxL_swift(iswift,1),3.,fdens,nudens)
                     FluxL_swift(iswift,1)=fdens
                  else if (flux_swift(iswift,4) .ne. 0.) then
                     flux_swift(iswift,1)=flux_swift(iswift,4)
                     Ferr_swift(iswift,1)=Ferr_swift(iswift,4)
                     FluxU_swift(iswift,1)=flux_swift(iswift,1)+Ferr_swift(iswift,1)
                     FluxL_swift(iswift,1)=flux_swift(iswift,1)-Ferr_swift(iswift,1)
                     call nhdeabsorb2 (0,2.,10.,0.9,nh,reduce,4)
                     flux_swift(iswift,1)=flux_swift(iswift,1)*reduce
                     FluxU_swift(iswift,1)=FluxU_swift(iswift,1)*reduce
                     FluxL_swift(iswift,1)=FluxL_swift(iswift,1)*reduce
                     call fluxtofdens(0.9,2.,10.,flux_swift(iswift,1),3.,fdens,nudens)
                     flux_swift(iswift,1)=fdens
                     frequency_swift(iswift,1)=2.418e17 !nudens 3 keV
                     call fluxtofdens(0.9,2.,10.,FluxU_swift(iswift,1),3.,fdens,nudens)
                     FluxU_swift(iswift,1)=fdens
                     call fluxtofdens(0.9,2.,10.,FluxL_swift(iswift,1),3.,fdens,nudens)
                     FluxL_swift(iswift,1)=fdens
                  endif
               endif
               call nhdeabsorb2 (0,0.3,1.,0.9,nh,reduce,4)
               flux_swift(iswift,2)=flux_swift(iswift,2)*reduce
               FluxU_swift(iswift,2)=FluxU_swift(iswift,2)*reduce
               FluxL_swift(iswift,2)=FluxL_swift(iswift,2)*reduce
               call fluxtofdens(0.9,0.3,1.,flux_swift(iswift,2),0.5,fdens,nudens)
               flux_swift(iswift,2)=fdens
               frequency_swift(iswift,2)=nudens
               call fluxtofdens(0.9,0.3,1.,FluxU_swift(iswift,2),0.5,fdens,nudens)
               FluxU_swift(iswift,2)=fdens
               call fluxtofdens(0.9,0.3,1.,FluxL_swift(iswift,2),0.5,fdens,nudens)
               FluxL_swift(iswift,2)=fdens
               call nhdeabsorb2 (0,1.,2.,0.9,nh,reduce,4)
               flux_swift(iswift,3)=flux_swift(iswift,3)*reduce
               FluxU_swift(iswift,3)=FluxU_swift(iswift,3)*reduce
               FluxL_swift(iswift,3)=FluxL_swift(iswift,3)*reduce
               call fluxtofdens(0.9,1.,2.,flux_swift(iswift,3),1.5,fdens,nudens)
               flux_swift(iswift,3)=fdens
               frequency_swift(iswift,3)=nudens
               call fluxtofdens(0.9,1.,2.,FluxU_swift(iswift,3),1.5,fdens,nudens)
               FluxU_swift(iswift,3)=fdens
               call fluxtofdens(0.9,1.,2.,FluxL_swift(iswift,3),1.5,fdens,nudens)
               FluxL_swift(iswift,3)=fdens
               call nhdeabsorb2 (0,2.,10.,0.9,nh,reduce,4)
               flux_swift(iswift,4)=flux_swift(iswift,4)*reduce
               FluxU_swift(iswift,4)=FluxU_swift(iswift,4)*reduce
               FluxL_swift(iswift,4)=FluxL_swift(iswift,4)*reduce
               call fluxtofdens(0.9,2.,10.,flux_swift(iswift,4),4.5,fdens,nudens)
               flux_swift(iswift,4)=fdens
               frequency_swift(iswift,4)=nudens
               call fluxtofdens(0.9,2.,10.,FluxU_swift(iswift,4),4.5,fdens,nudens)
               FluxU_swift(iswift,4)=fdens
               call fluxtofdens(0.9,2.,10.,FluxL_swift(iswift,4),4.5,fdens,nudens)
               FluxL_swift(iswift,4)=fdens
               if ((FluxU_swift(iswift,2) .ne. 0.) .and. (flux_swift(iswift,2) .eq. 0.)) then
                  FluxU_swift(iswift,2)=FluxU_swift(iswift,2)*3.
               endif
               if ((FluxU_swift(iswift,3) .ne. 0.) .and. (flux_swift(iswift,3) .eq. 0.)) then
                  FluxU_swift(iswift,3)=FluxU_swift(iswift,3)*3.
               endif
               if ((FluxU_swift(iswift,4) .ne. 0.) .and. (flux_swift(iswift,4) .eq. 0.)) then
                  FluxU_swift(iswift,4)=FluxU_swift(iswift,4)*3.
               endif
            else if ((catalog(1:7) == 'xrtdeep') .or. (catalog(1:5) == 'sds82') .or. (catalog(1:5) == '1ousx')
     &        .or.  (catalog(1:5) == 'ousxb') .or. (catalog(1:5) == 'ousxg')) then
               if ((catalog(1:5) == '1ousx')) then
                  xrt_type(iswift)=3
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdst_swift(iswift)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) mjded_swift(iswift)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_swift(iswift,2)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,2)
                  FluxU_swift(iswift,2)=flux_swift(iswift,2)+Ferr_swift(iswift,2)
                  FluxL_swift(iswift,2)=flux_swift(iswift,2)-Ferr_swift(iswift,2)
                  frequency_swift(iswift,2)=(1.602E-19)*(5.e2)/(6.626e-34)
                  if ((Ferr_swift(iswift,2) .lt. 0) .or. (FluxL_swift(iswift,2) .lt. 0)) then
c                     if (flux_swift(iswift,2) .gt. 0.) then
                        FluxU_swift(iswift,2)=flux_swift(iswift,2)
                        flux_swift(iswift,2)=0.
                        FluxL_swift(iswift,2)=0.
c                     else
c                        FluxL_swift(iswift,2)=0.
c                        FluxU_swift(iswift,2)=0.
c                        flux_swift(iswift,2)=0.
c                     endif
                  endif
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_swift(iswift,1)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,1)
                  FluxU_swift(iswift,1)=flux_swift(iswift,1)+Ferr_swift(iswift,1)
                  FluxL_swift(iswift,1)=flux_swift(iswift,1)-Ferr_swift(iswift,1)
                  frequency_swift(iswift,1)=(1.602E-19)*(1.e3)/(6.626e-34)
                  if ((Ferr_swift(iswift,1) .lt. 0) .or. (FluxL_swift(iswift,1) .lt. 0)) then
c                     if (flux_swift(iswift,1) .gt. 0.) then
                        FluxU_swift(iswift,1)=Ferr_swift(iswift,1)
                        flux_swift(iswift,1)=0.
                        FluxL_swift(iswift,1)=0.
c                     else
c                        FluxL_swift(iswift,1)=0.
c                        FluxU_swift(iswift,1)=0.
c                        flux_swift(iswift,1)=0.
c                     endif
                  endif
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_swift(iswift,3)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,3)
                  FluxU_swift(iswift,3)=flux_swift(iswift,3)+Ferr_swift(iswift,3)
                  FluxL_swift(iswift,3)=flux_swift(iswift,3)-Ferr_swift(iswift,3)
                  frequency_swift(iswift,3)=(1.602E-19)*(1.5e3)/(6.626e-34)
                  if ((Ferr_swift(iswift,3) .lt. 0) .or. (FluxL_swift(iswift,3) .lt. 0)) then
c                     if (flux_swift(iswift,3) .gt. 0.) then
                        FluxU_swift(iswift,3)=Ferr_swift(iswift,3)
                        flux_swift(iswift,3)=0.
                        FluxL_swift(iswift,3)=0.
c                     else
c                        FluxL_swift(iswift,3)=0.
c                        FluxU_swift(iswift,3)=0.
c                        flux_swift(iswift,3)=0.
c                     endif
                  endif
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_swift(iswift,4)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,4)
                  FluxU_swift(iswift,4)=flux_swift(iswift,4)+Ferr_swift(iswift,4)
                  FluxL_swift(iswift,4)=flux_swift(iswift,4)-Ferr_swift(iswift,4)
                  frequency_swift(iswift,4)=(1.602E-19)*(3.e3)/(6.626e-34)
                  if ((Ferr_swift(iswift,4) .lt. 0) .or. (FluxL_swift(iswift,4) .lt. 0)) then
c                     if (flux_swift(iswift,4) .gt. 0.) then
                        FluxU_swift(iswift,4)=Ferr_swift(iswift,4)
                        flux_swift(iswift,4)=0.
                        FluxL_swift(iswift,4)=0.
c                     else
c                        FluxL_swift(iswift,4)=0.
c                        FluxU_swift(iswift,4)=0.
c                        flux_swift(iswift,4)=0.
c                     endif
                  endif
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_swift(iswift,5)
                  is=ie
                  ie=index(string(is+1:len(string)),' ')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,5)
                  FluxU_swift(iswift,5)=flux_swift(iswift,5)+Ferr_swift(iswift,5)
                  FluxL_swift(iswift,5)=flux_swift(iswift,5)-Ferr_swift(iswift,5)
                  frequency_swift(iswift,5)=(1.602E-19)*(4.5e3)/(6.626e-34)
                  if ((Ferr_swift(iswift,5) .lt. 0) .or. (FluxL_swift(iswift,5) .lt. 0)) then
c                     if (flux_swift(iswift,5) .gt. 0.) then
                        FluxU_swift(iswift,5)=Ferr_swift(iswift,5)
                        flux_swift(iswift,5)=0.
                        FluxL_swift(iswift,5)=0.
c                     else
c                        FluxL_swift(iswift,5)=0.
c                        FluxU_swift(iswift,5)=0.
c                        flux_swift(iswift,5)=0.
c                     endif
                  endif
                  if (flux_swift(iswift,1) .eq. 0.) then
                     if (flux_swift(iswift,3) .ne. 0.) then
                        flux_swift(iswift,1)=(flux_swift(iswift,3)/frequency_swift(iswift,3))*(1./1.5)**(-0.9)
                        FluxU_swift(iswift,1)=(FluxU_swift(iswift,3)/frequency_swift(iswift,3))*(1./1.5)**(-0.9)
                        FluxL_swift(iswift,1)=(FluxL_swift(iswift,3)/frequency_swift(iswift,3))*(1./1.5)**(-0.9)
                     else
                        flux_swift(iswift,1)=(flux_swift(iswift,4)/frequency_swift(iswift,4))*(1./3.)**(-0.9)
                        FluxU_swift(iswift,1)=(FluxU_swift(iswift,4)/frequency_swift(iswift,4))*(1./3.)**(-0.9)
                        FluxL_swift(iswift,1)=(FluxL_swift(iswift,4)/frequency_swift(iswift,4))*(1./3.)**(-0.9)
                     endif
                     flux_swift(iswift,1)=flux_swift(iswift,1)*frequency_swift(iswift,1)
                     FluxU_swift(iswift,1)=FluxU_swift(iswift,1)*frequency_swift(iswift,1)
                     FluxL_swift(iswift,1)=FluxL_swift(iswift,1)*frequency_swift(iswift,1)
                  endif
               else
                  xrt_type(iswift)=2
                  if (catalog(1:7) == 'xrtdeep') then
                     is=ie
                     ie=index(string(is+1:len(string)),',')+is !nh
                     is=ie
                     ie=index(string(is+1:len(string)),',')+is !slope
                     is=ie
                     ie=index(string(is+1:len(string)),',')+is !slope err
                     is=ie
                     ie=index(string(is+1:len(string)),',')+is !exp
                  endif
                  ra_swift(iswift)=-ra_swift(iswift)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is !3kev
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_swift(iswift,1)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,1)
                  FluxU_swift(iswift,1)=flux_swift(iswift,1)+Ferr_swift(iswift,1)
                  FluxL_swift(iswift,1)=flux_swift(iswift,1)-Ferr_swift(iswift,1)
                  frequency_swift(iswift,1)=(1.602E-19)*(1.e3)/(6.626e-34)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is !0.5kev
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_swift(iswift,2)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,2)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  FluxU_swift(iswift,2)=flux_swift(iswift,2)+Ferr_swift(iswift,2)
                  FluxL_swift(iswift,2)=flux_swift(iswift,2)-Ferr_swift(iswift,2)
                  if ((Ferr_swift(iswift,2) .lt. 0) .or. (FluxL_swift(iswift,2) .lt. 0)) then
                     if (is .ne. ie-1) read(string(is+1:ie-1),*) FluxU_swift(iswift,2)
                     if (FluxU_swift(iswift,2) .gt. 0.) then
                        FluxL_swift(iswift,2)=0.
                     else
                        FluxL_swift(iswift,2)=0.
                        FluxU_swift(iswift,2)=0.!3.*Ferr_swift(iswift,2)
                        flux_swift(iswift,2)=0.
                     endif
                  endif
                  frequency_swift(iswift,2)=(1.602E-19)*(5.e2)/(6.626e-34)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is !1.5kev
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_swift(iswift,3)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is !1.5kev
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,3)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  FluxU_swift(iswift,3)=flux_swift(iswift,3)+Ferr_swift(iswift,3)
                  FluxL_swift(iswift,3)=flux_swift(iswift,3)-Ferr_swift(iswift,3)
                  if ((Ferr_swift(iswift,3) .lt. 0) .or. (FluxL_swift(iswift,3) .lt. 0)) then
                     if (is .ne. ie-1) read(string(is+1:ie-1),*) FluxU_swift(iswift,3)
                     if (FluxU_swift(iswift,3) .gt. 0.) then
                        FluxL_swift(iswift,3)=0.
                     else
                        FluxL_swift(iswift,3)=0.
                        FluxU_swift(iswift,3)=0.!3.*Ferr_swift(iswift,3)
                        flux_swift(iswift,3)=0.
                     endif
                  endif
                  frequency_swift(iswift,3)=(1.602E-19)*(1.5E3)/(6.626e-34)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is !4.5kev
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_swift(iswift,4)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is !4.5kev
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_swift(iswift,4)
                  is=ie
                  ie=index(string(is+1:len(string)),' ')+is
                  FluxU_swift(iswift,4)=flux_swift(iswift,4)+Ferr_swift(iswift,4)
                  FluxL_swift(iswift,4)=flux_swift(iswift,4)-Ferr_swift(iswift,4)
                  if ((Ferr_swift(iswift,4) .lt. 0) .or. (FluxL_swift(iswift,4) .lt. 0)) then
                     if (is .ne. ie-1) read(string(is+1:ie-1),*) FluxU_swift(iswift,4)
                     if (FluxU_swift(iswift,4) .gt. 0.) then
                        FluxL_swift(iswift,4)=0.
                     else
                        FluxL_swift(iswift,4)=0.
                        FluxU_swift(iswift,4)=0.!3.*Ferr_swift(iswift,4)
                        flux_swift(iswift,4)=0.
                     endif
                  endif
                  frequency_swift(iswift,4)=(1.602E-19)*(4.5e3)/(6.626e-34)
               endif
               poserr_swift(iswift)=8.!!!!!!!!
            endif
c PG
               CALL RXgraphic_code(flux_swift(iswift,1),'X',code)
               write (13,'(f9.5,2x,f9.5,2x,i6)') abs(ra_swift(iswift)),dec_swift(iswift),int(code)
c            write(*,*) frequency_swift(iswift,5)
c end PG
         ELSE IF (catalog(1:3) == 'ipc') THEN
            iipc=iipc+1
            IF (iipc > arrsize(12)) Stop 'Too many Einstein IPC points'
            ra_ipc(iipc)=ra
            dec_ipc(iipc)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_ipc(iipc)
            if (catalog(1:5) == 'ipcsl') then
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_ipc(iipc)
               poserr_ipc(iipc)=72. !!!!90% error
               ipc_type(iipc)=1
            else
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_ipc(iipc)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_ipc(iipc)
               ipc_type(iipc)=2
            endif
            FluxU_ipc(iipc)=flux_ipc(iipc)+Ferr_ipc(iipc)
            FluxL_ipc(iipc)=flux_ipc(iipc)-Ferr_ipc(iipc)
            call nhdeabsorb2 (0,0.2,3.5,0.9,nh,reduce,3)
            flux_ipc(iipc)=flux_ipc(iipc)*reduce
            FluxU_ipc(iipc)=FluxU_ipc(iipc)*reduce
            FluxL_ipc(iipc)=FluxL_ipc(iipc)*reduce
            call fluxtofdens(0.9,0.2,3.5,flux_ipc(iipc),1.,fdens,nudens)
            flux_ipc(iipc)=fdens
            frequency_ipc(iipc)=nudens
            call fluxtofdens(0.9,0.2,3.5,FluxU_ipc(iipc),1.,fdens,nudens)
            FluxU_ipc(iipc)=fdens
            call fluxtofdens(0.9,0.2,3.5,FluxL_ipc(iipc),1.,fdens,nudens)
            FluxL_ipc(iipc)=fdens
            !write(*,*) catalog,FluxU_ipc(iipc),flux_ipc(iipc),FluxL_ipc(iipc),poserr_ipc(iipc)
c PG
            CALL RXgraphic_code(flux_ipc(iipc),'X',code)
            write (13,'(f9.5,2x,f9.5,2x,i6)') ra_ipc(iipc),dec_ipc(iipc),int(code)
c end PG
         ELSE IF (catalog(1:3) == 'bmw') THEN
            ibmw=ibmw+1
            IF (ibmw > arrsize(11)) Stop 'Too many ROSAT-BMW points'
            ra_bmw(ibmw)=ra
            dec_bmw(ibmw)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_bmw(ibmw)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_bmw(ibmw)
            is=ie
            ie=index(string(is+1:len(string)),' ')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_bmw(ibmw)
            !call fluxtofdens(0.9,0.5,7.,flux_bmw(ibmw,1),1.,fdens,nudens)
            FluxU_bmw(ibmw)=flux_bmw(ibmw)+Ferr_bmw(ibmw)
            FluxL_bmw(ibmw)=flux_bmw(ibmw)-Ferr_bmw(ibmw)
            call nhdeabsorb2 (0,0.1,2.4,0.9,nh,reduce,2)
            flux_bmw(ibmw)=flux_bmw(ibmw)*reduce
            FluxU_bmw(ibmw)=FluxU_bmw(ibmw)*reduce
            FluxL_bmw(ibmw)=FluxL_bmw(ibmw)*reduce
            call fluxtofdens(0.9,0.1,2.4,flux_bmw(ibmw),1.,fdens,nudens)
            flux_bmw(ibmw)=fdens
            frequency_bmw(ibmw)=nudens
            call fluxtofdens(0.9,0.1,2.4,FluxU_bmw(ibmw),1.,fdens,nudens)
            FluxU_bmw(ibmw)=fdens
            call fluxtofdens(0.9,0.1,2.4,FluxL_bmw(ibmw),1.,fdens,nudens)
            FluxL_bmw(ibmw)=fdens
c PG
            CALL RXgraphic_code(flux_bmw(ibmw),'X',code)
            write (13,'(f9.5,2x,f9.5,2x,i6)') ra_bmw(ibmw),dec_bmw(ibmw),int(code)
c end PG
         ELSE IF (catalog(1:11) == 'chandracsc2') THEN
            ichandra=ichandra+1
            IF (ichandra > arrsize(9)) Stop 'Too many Chandra points'
            ra_chandra(ichandra)=ra
            dec_chandra(ichandra)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_chandra(ichandra)
            poserr_chandra(ichandra)=sqrt((poserr_chandra(ichandra)**2)+(0.8**2))
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_chandra(ichandra,4)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_chandra(ichandra,2) !h
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_chandra(ichandra,3) !m
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_chandra(ichandra,1) !s
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_chandra(ichandra,5) !us
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxU_chandra(ichandra,4)
            if ((flux_chandra(ichandra,4) .eq. 0.) .and. (FluxU_chandra(ichandra,4) .ne. 0.))
     &             FluxU_chandra(ichandra,4)=FluxU_chandra(ichandra,4)*3.
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxU_chandra(ichandra,2) !h
            if ((flux_chandra(ichandra,2) .eq. 0.) .and. (FluxU_chandra(ichandra,2) .ne. 0.))
     &             FluxU_chandra(ichandra,2)=FluxU_chandra(ichandra,2)*3.
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxU_chandra(ichandra,3) !m
            if ((flux_chandra(ichandra,3) .eq. 0.) .and. (FluxU_chandra(ichandra,3) .ne. 0.))
     &             FluxU_chandra(ichandra,3)=FluxU_chandra(ichandra,3)*3.
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxU_chandra(ichandra,1) !s
            if ((flux_chandra(ichandra,1) .eq. 0.) .and. (FluxU_chandra(ichandra,1) .ne. 0.))
     &             FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,1)*3.
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxU_chandra(ichandra,5) !us
            if ((flux_chandra(ichandra,5) .eq. 0.) .and. (FluxU_chandra(ichandra,5) .ne. 0.))
     &             FluxU_chandra(ichandra,5)=FluxU_chandra(ichandra,5)*3.
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxL_chandra(ichandra,4)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxL_chandra(ichandra,2) !h
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxL_chandra(ichandra,3) !m
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxL_chandra(ichandra,1) !s
            is=ie
            ie=index(string(is+1:len(string)),' ')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) FluxL_chandra(ichandra,5) !us
            call nhdeabsorb2 (1,0.5,1.2,0.9,nh,reduce,100)
            flux_chandra(ichandra,1)=flux_chandra(ichandra,1)*reduce
            call fluxtofdens(0.9,0.5,1.2,flux_chandra(ichandra,1),1.,fdens,nudens)
            flux_chandra(ichandra,1)=fdens
            frequency_chandra(ichandra,1)=nudens
            FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,1)*reduce
            call fluxtofdens(0.9,0.5,1.2,FluxU_chandra(ichandra,1),1.,fdens,nudens)
            FluxU_chandra(ichandra,1)=fdens
            FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,1)*reduce
            call fluxtofdens(0.9,0.5,1.2,FluxL_chandra(ichandra,1),1.,fdens,nudens)
            FluxL_chandra(ichandra,1)=fdens
            if (flux_chandra(ichandra,1) .eq. 0.) then
               if (flux_chandra(ichandra,4) .ne. 0.) then
                  flux_chandra(ichandra,1)=flux_chandra(ichandra,4)
                  FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,4)
                  FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,4)
                  call nhdeabsorb2 (1,0.5,7.,0.9,nh,reduce,100)
                  flux_chandra(ichandra,1)=flux_chandra(ichandra,1)*reduce
                  FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,1)*reduce
                  FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,1)*reduce
                  call fluxtofdens(0.9,0.5,7.,flux_chandra(ichandra,1),1.,fdens,nudens)
                  flux_chandra(ichandra,1)=fdens
                  frequency_chandra(ichandra,1)=nudens
                  call fluxtofdens(0.9,0.5,7.,fluxU_chandra(ichandra,1),1.,fdens,nudens)
                  FluxU_chandra(ichandra,1)=fdens
                  call fluxtofdens(0.9,0.5,7.,fluxL_chandra(ichandra,1),1.,fdens,nudens)
                  FluxL_chandra(ichandra,1)=fdens
               ELSE if (flux_chandra(ichandra,3) .ne. 0) then
                  flux_chandra(ichandra,1)=flux_chandra(ichandra,3)
                  FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,3)
                  FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,3)
                  call nhdeabsorb2 (1,1.2,2.,0.9,nh,reduce,100)
                  flux_chandra(ichandra,1)=flux_chandra(ichandra,1)*reduce
                  FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,1)*reduce
                  FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,1)*reduce
                  call fluxtofdens(0.9,1.2,2.,flux_chandra(ichandra,1),1.,fdens,nudens)
                  flux_chandra(ichandra,1)=fdens
                  frequency_chandra(ichandra,1)=nudens
                  call fluxtofdens(0.9,1.2,2.,fluxU_chandra(ichandra,1),1.,fdens,nudens)
                  FluxU_chandra(ichandra,1)=fdens
                  call fluxtofdens(0.9,1.2,2.,fluxL_chandra(ichandra,1),1.,fdens,nudens)
                  FluxL_chandra(ichandra,1)=fdens
               ELSE if (flux_chandra(ichandra,5) .ne. 0) then
                  flux_chandra(ichandra,1)=flux_chandra(ichandra,5)
                  FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,5)
                  FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,5)
                  call nhdeabsorb2 (1,0.2,0.5,0.9,nh,reduce,100)
                  flux_chandra(ichandra,1)=flux_chandra(ichandra,1)*reduce
                  FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,1)*reduce
                  FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,1)*reduce
                  call fluxtofdens(0.9,0.2,0.5,flux_chandra(ichandra,1),1.,fdens,nudens)
                  flux_chandra(ichandra,1)=fdens
                  frequency_chandra(ichandra,1)=nudens
                  call fluxtofdens(0.9,0.2,0.5,fluxU_chandra(ichandra,1),1.,fdens,nudens)
                  FluxU_chandra(ichandra,1)=fdens
                  call fluxtofdens(0.9,0.2,0.5,fluxL_chandra(ichandra,1),1.,fdens,nudens)
                  FluxL_chandra(ichandra,1)=fdens
               ELSE
                  flux_chandra(ichandra,1)=flux_chandra(ichandra,2)
                  FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,2)
                  FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,2)
                  call nhdeabsorb2 (1,2.,7.,0.9,nh,reduce,100)
                  flux_chandra(ichandra,1)=flux_chandra(ichandra,1)*reduce
                  FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,1)*reduce
                  FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,1)*reduce
                  call fluxtofdens(0.9,2.,7.,flux_chandra(ichandra,1),1.,fdens,nudens)
                  flux_chandra(ichandra,1)=fdens
                  frequency_chandra(ichandra,1)=nudens
                  call fluxtofdens(0.9,2.,7.,fluxU_chandra(ichandra,1),1.,fdens,nudens)
                  FluxU_chandra(ichandra,1)=fdens
                  call fluxtofdens(0.9,2.,7.,fluxL_chandra(ichandra,1),1.,fdens,nudens)
                  FluxL_chandra(ichandra,1)=fdens
               endif
            endif
            call nhdeabsorb2 (1,2.,7.,0.9,nh,reduce,100)
            flux_chandra(ichandra,2)=flux_chandra(ichandra,2)*reduce
            call fluxtofdens(0.9,2.,7.,flux_chandra(ichandra,2),4.5,fdens,nudens)
            flux_chandra(ichandra,2)=fdens
            frequency_chandra(ichandra,2)=nudens
            FluxU_chandra(ichandra,2)=FluxU_chandra(ichandra,2)*reduce
            call fluxtofdens(0.9,2.,7.,fluxU_chandra(ichandra,2),4.5,fdens,nudens)
            FluxU_chandra(ichandra,2)=fdens
            FluxL_chandra(ichandra,2)=FluxL_chandra(ichandra,2)*reduce
            call fluxtofdens(0.9,2.,7.,fluxL_chandra(ichandra,2),4.5,fdens,nudens)
            FluxL_chandra(ichandra,2)=fdens
            call nhdeabsorb2 (1,1.2,2.,0.9,nh,reduce,100)
            flux_chandra(ichandra,3)=flux_chandra(ichandra,3)*reduce
            call fluxtofdens(0.9,1.2,2.,flux_chandra(ichandra,3),1.5,fdens,nudens)
            flux_chandra(ichandra,3)=fdens
            frequency_chandra(ichandra,3)=nudens
            FluxU_chandra(ichandra,3)=FluxU_chandra(ichandra,3)*reduce
            call fluxtofdens(0.9,1.2,2.,fluxU_chandra(ichandra,3),1.5,fdens,nudens)
            FluxU_chandra(ichandra,3)=fdens
            FluxL_chandra(ichandra,3)=FluxL_chandra(ichandra,3)*reduce
            call fluxtofdens(0.9,1.2,2.,fluxL_chandra(ichandra,3),1.5,fdens,nudens)
            FluxL_chandra(ichandra,3)=fdens
            call nhdeabsorb2 (1,0.5,7.,0.9,nh,reduce,100)
            flux_chandra(ichandra,4)=flux_chandra(ichandra,4)*reduce
            call fluxtofdens(0.9,0.5,7.,flux_chandra(ichandra,4),3.,fdens,nudens)
            flux_chandra(ichandra,4)=fdens
            frequency_chandra(ichandra,4)=nudens
            FluxU_chandra(ichandra,4)=FluxU_chandra(ichandra,4)*reduce
            call fluxtofdens(0.9,0.5,7.,fluxU_chandra(ichandra,4),3.,fdens,nudens)
            FluxU_chandra(ichandra,4)=fdens
            FluxL_chandra(ichandra,4)=FluxL_chandra(ichandra,4)*reduce
            call fluxtofdens(0.9,0.5,7.,fluxL_chandra(ichandra,4),3.,fdens,nudens)
            FluxL_chandra(ichandra,4)=fdens
            call nhdeabsorb2 (1,0.2,0.5,0.9,nh,reduce,100)
            flux_chandra(ichandra,5)=flux_chandra(ichandra,5)*reduce
            call fluxtofdens(0.9,0.2,0.5,flux_chandra(ichandra,5),0.3,fdens,nudens)
            flux_chandra(ichandra,5)=fdens
            frequency_chandra(ichandra,5)=nudens
            FluxU_chandra(ichandra,5)=FluxU_chandra(ichandra,5)*reduce
            call fluxtofdens(0.9,0.2,0.5,fluxU_chandra(ichandra,5),0.3,fdens,nudens)
            FluxU_chandra(ichandra,5)=fdens
            FluxL_chandra(ichandra,5)=FluxL_chandra(ichandra,5)*reduce
            call fluxtofdens(0.9,0.2,0.5,fluxL_chandra(ichandra,5),0.3,fdens,nudens)
            FluxL_chandra(ichandra,5)=fdens
            !write(*,*) fluxU_chandra(ichandra,1),flux_chandra(ichandra,1),fluxL_chandra(ichandra,1)
            !write(*,*) fluxU_chandra(ichandra,5),flux_chandra(ichandra,5),fluxL_chandra(ichandra,5)
c PG
            CALL RXgraphic_code(flux_chandra(ichandra,1),'X',code)
            write (13,'(f9.5,2x,f9.5,2x,i6)') ra_chandra(ichandra),dec_chandra(ichandra),int(code)
c end PG
         ELSE IF ((catalog(1:7) == 'maxissc') .or. (catalog(1:7) == 'maxigsc')) THEN
            imaxi=imaxi+1
            ra_maxi(imaxi)=ra
            dec_maxi(imaxi)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (catalog(1:7) == 'maxigsc') then
               ra_maxi(imaxi)=-ra_maxi(imaxi)
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_maxi(imaxi)
               poserr_maxi(imaxi)=2.*poserr_maxi(imaxi)*3600.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_maxi(imaxi,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_maxi(imaxi,3)
               FluxU_maxi(imaxi,3)=flux_maxi(imaxi,3)+Ferr_maxi(imaxi,3)
               FluxL_maxi(imaxi,3)=flux_maxi(imaxi,3)-Ferr_maxi(imaxi,3)
               call nhdeabsorb2 (1,4.,10.,0.9,nh,reduce,100)
               flux_maxi(imaxi,3)=flux_maxi(imaxi,3)*reduce*1.e-12
               FluxU_maxi(imaxi,3)=FluxU_maxi(imaxi,3)*reduce*1.e-12
               FluxL_maxi(imaxi,3)=FluxL_maxi(imaxi,3)*reduce*1.e-12
               call fluxtofdens(0.9,4.,10.,flux_maxi(imaxi,3),7.,fdens,nudens)
               flux_maxi(imaxi,3)=fdens
               frequency_maxi(imaxi,3)=nudens
               call fluxtofdens(0.9,4.,10.,FluxU_maxi(imaxi,3),7.,fdens,nudens)
               FluxU_maxi(imaxi,3)=fdens
               call fluxtofdens(0.9,4.,10.,FluxL_maxi(imaxi,3),7.,fdens,nudens)
               FluxL_maxi(imaxi,3)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_maxi(imaxi,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_maxi(imaxi,2)
               FluxU_maxi(imaxi,2)=flux_maxi(imaxi,2)+Ferr_maxi(imaxi,2)
               FluxL_maxi(imaxi,2)=flux_maxi(imaxi,2)-Ferr_maxi(imaxi,2)
               if (FluxL_maxi(imaxi,2) .lt. 0.) then
                  FluxU_maxi(imaxi,2)=Ferr_maxi(imaxi,2)*3.
                  FluxL_maxi(imaxi,2)=0.
               endif
               call nhdeabsorb2 (1,3.,4.,0.9,nh,reduce,100)
               flux_maxi(imaxi,2)=flux_maxi(imaxi,2)*reduce*1.e-12
               FluxU_maxi(imaxi,2)=FluxU_maxi(imaxi,2)*reduce*1.e-12
               FluxL_maxi(imaxi,2)=FluxL_maxi(imaxi,2)*reduce*1.e-12
               call fluxtofdens(0.9,3.,4.,flux_maxi(imaxi,2),3.5,fdens,nudens)
               flux_maxi(imaxi,2)=fdens
               frequency_maxi(imaxi,2)=nudens
               call fluxtofdens(0.9,3.,4.,FluxU_maxi(imaxi,2),3.5,fdens,nudens)
               FluxU_maxi(imaxi,2)=fdens
               call fluxtofdens(0.9,3.,4.,FluxL_maxi(imaxi,2),3.5,fdens,nudens)
               FluxL_maxi(imaxi,2)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_maxi(imaxi,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_maxi(imaxi,4)
               FluxU_maxi(imaxi,4)=flux_maxi(imaxi,4)+Ferr_maxi(imaxi,4)
               FluxL_maxi(imaxi,4)=flux_maxi(imaxi,4)-Ferr_maxi(imaxi,4)
               if (FluxL_maxi(imaxi,4) .lt. 0.) then
                  FluxU_maxi(imaxi,4)=Ferr_maxi(imaxi,4)*3.
                  FluxL_maxi(imaxi,4)=0.
               endif
               call nhdeabsorb2 (1,10.,20.,0.9,nh,reduce,100)
               flux_maxi(imaxi,4)=flux_maxi(imaxi,4)*reduce*1.e-12
               FluxU_maxi(imaxi,4)=FluxU_maxi(imaxi,4)*reduce*1.e-12
               FluxL_maxi(imaxi,4)=FluxL_maxi(imaxi,4)*reduce*1.e-12
               call fluxtofdens(0.9,10.,20.,flux_maxi(imaxi,4),15.,fdens,nudens)
               flux_maxi(imaxi,4)=fdens
               frequency_maxi(imaxi,4)=nudens
               call fluxtofdens(0.9,10.,20.,FluxU_maxi(imaxi,4),15.,fdens,nudens)
               FluxU_maxi(imaxi,4)=fdens
               call fluxtofdens(0.9,10.,20.,FluxL_maxi(imaxi,4),15.,fdens,nudens)
               FluxL_maxi(imaxi,4)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_maxi(imaxi,1)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_maxi(imaxi,1)
               FluxU_maxi(imaxi,1)=flux_maxi(imaxi,1)+Ferr_maxi(imaxi,1)
               FluxL_maxi(imaxi,1)=flux_maxi(imaxi,1)-Ferr_maxi(imaxi,1)
               call nhdeabsorb2 (1,3.,10.,0.9,nh,reduce,100)
               flux_maxi(imaxi,1)=flux_maxi(imaxi,1)*reduce*1.e-12
               FluxU_maxi(imaxi,1)=FluxU_maxi(imaxi,1)*reduce*1.e-12
               FluxL_maxi(imaxi,1)=FluxL_maxi(imaxi,1)*reduce*1.e-12
               call fluxtofdens(0.9,3.,10.,flux_maxi(imaxi,1),5.,fdens,nudens)
               flux_maxi(imaxi,1)=fdens
               frequency_maxi(imaxi,1)=2.418E17!nudens check the flux value later
               call fluxtofdens(0.9,3.,10.,FluxU_maxi(imaxi,1),5.,fdens,nudens)
               FluxU_maxi(imaxi,1)=fdens
               call fluxtofdens(0.9,3.,10.,FluxL_maxi(imaxi,1),5.,fdens,nudens)
               FluxL_maxi(imaxi,1)=fdens
               maxi_type(imaxi)=1
            else
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_maxi(imaxi,1) !s
               call nhdeabsorb2 (1,0.7,1.85,0.9,nh,reduce,100)
               flux_maxi(imaxi,1)=flux_maxi(imaxi,1)*reduce
               call fluxtofdens(0.9,0.7,1.85,flux_maxi(imaxi,1),1.,fdens,nudens)
               flux_maxi(imaxi,1)=fdens
               frequency_maxi(imaxi,1)=nudens
               FluxU_maxi(imaxi,1)=0.
               FluxL_maxi(imaxi,1)=0.
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_maxi(imaxi,2) !h
               call nhdeabsorb2 (1,1.85,7.,0.9,nh,reduce,100)
               flux_maxi(imaxi,2)=flux_maxi(imaxi,2)*reduce
               call fluxtofdens(0.9,1.85,7.,flux_maxi(imaxi,2),4.4,fdens,nudens)
               flux_maxi(imaxi,2)=fdens
               frequency_maxi(imaxi,2)=nudens
               FluxU_maxi(imaxi,2)=0.
               FluxL_maxi(imaxi,2)=0.
               poserr_maxi(imaxi)=sqrt(0.2**2+0.2**2)*3600.
               maxi_type(imaxi)=2
            endif
            CALL RXgraphic_code(flux_maxi(imaxi,1),'X',code)
            write (13,'(f9.5,2x,f9.5,2x,i6)') abs(ra_maxi(imaxi)),dec_maxi(imaxi),int(code)
         ELSE IF (catalog(1:7) == 'erosita') THEN
            ierosita=ierosita+1
            ra_erosita(ierosita)=ra
            dec_erosita(ierosita)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_erosita(ierosita)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_erosita(ierosita,1)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_erosita(ierosita,1)
            fluxU_erosita(ierosita,1)=flux_erosita(ierosita,1)+Ferr_erosita(ierosita,1)
            fluxL_erosita(ierosita,1)=flux_erosita(ierosita,1)-Ferr_erosita(ierosita,1)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
ccc            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_erosita(ierosita,2)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
ccc            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,2)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_erosita(ierosita,2)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_erosita(ierosita,3)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_erosita(ierosita,4)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_erosita(ierosita,5)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_erosita(ierosita,6)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_erosita(ierosita,7)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_erosita(ierosita,8)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,2)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,3)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,4)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,5)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,6)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,7)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,8)
            fluxL_erosita(ierosita,2)=flux_erosita(ierosita,2)-Ferr_erosita(ierosita,2)
            fluxL_erosita(ierosita,3)=flux_erosita(ierosita,3)-Ferr_erosita(ierosita,3)
            fluxL_erosita(ierosita,4)=flux_erosita(ierosita,4)-Ferr_erosita(ierosita,4)
            fluxL_erosita(ierosita,5)=flux_erosita(ierosita,5)-Ferr_erosita(ierosita,5)
            fluxL_erosita(ierosita,6)=flux_erosita(ierosita,6)-Ferr_erosita(ierosita,6)
            fluxL_erosita(ierosita,7)=flux_erosita(ierosita,7)-Ferr_erosita(ierosita,7)
            fluxL_erosita(ierosita,8)=flux_erosita(ierosita,8)-Ferr_erosita(ierosita,8)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,2)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,3)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,4)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,5)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,6)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,7)
            is=ie
            ie=index(string(is+1:len(string)),' ')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_erosita(ierosita,8)
            fluxU_erosita(ierosita,2)=flux_erosita(ierosita,2)+Ferr_erosita(ierosita,2)
            fluxU_erosita(ierosita,3)=flux_erosita(ierosita,3)+Ferr_erosita(ierosita,3)
            fluxU_erosita(ierosita,4)=flux_erosita(ierosita,4)+Ferr_erosita(ierosita,4)
            fluxU_erosita(ierosita,5)=flux_erosita(ierosita,5)+Ferr_erosita(ierosita,5)
            fluxU_erosita(ierosita,6)=flux_erosita(ierosita,6)+Ferr_erosita(ierosita,6)
            fluxU_erosita(ierosita,7)=flux_erosita(ierosita,7)+Ferr_erosita(ierosita,7)
            fluxU_erosita(ierosita,8)=flux_erosita(ierosita,8)+Ferr_erosita(ierosita,8)

            call nhdeabsorb2 (0,0.2,2.3,0.9,nh,reduce,5)
            flux_erosita(ierosita,1)=flux_erosita(ierosita,1)*reduce
            fluxU_erosita(ierosita,1)=fluxU_erosita(ierosita,1)*reduce
            fluxL_erosita(ierosita,1)=fluxL_erosita(ierosita,1)*reduce
            call nhdeabsorb2 (0,0.2,0.5,0.9,nh,reduce,5)
            flux_erosita(ierosita,2)=flux_erosita(ierosita,2)*reduce
            fluxU_erosita(ierosita,2)=fluxU_erosita(ierosita,2)*reduce
            fluxL_erosita(ierosita,2)=fluxL_erosita(ierosita,2)*reduce
            call nhdeabsorb2 (0,0.5,1.,0.9,nh,reduce,5)
            flux_erosita(ierosita,3)=flux_erosita(ierosita,3)*reduce
            fluxU_erosita(ierosita,3)=fluxU_erosita(ierosita,3)*reduce
            fluxL_erosita(ierosita,3)=fluxL_erosita(ierosita,3)*reduce
            call nhdeabsorb2 (0,1.,2.,0.9,nh,reduce,5)
            flux_erosita(ierosita,4)=flux_erosita(ierosita,4)*reduce
            fluxU_erosita(ierosita,4)=fluxU_erosita(ierosita,4)*reduce
            fluxL_erosita(ierosita,4)=fluxL_erosita(ierosita,4)*reduce
            call nhdeabsorb2 (0,2.,4.5,0.9,nh,reduce,5)
            flux_erosita(ierosita,5)=flux_erosita(ierosita,5)*reduce
            fluxU_erosita(ierosita,5)=fluxU_erosita(ierosita,5)*reduce
            fluxL_erosita(ierosita,5)=fluxL_erosita(ierosita,5)*reduce
            call nhdeabsorb2 (0,0.5,2.,0.9,nh,reduce,5)
            flux_erosita(ierosita,6)=flux_erosita(ierosita,6)*reduce
            fluxU_erosita(ierosita,6)=fluxU_erosita(ierosita,6)*reduce
            fluxL_erosita(ierosita,6)=fluxL_erosita(ierosita,6)*reduce
            call nhdeabsorb2 (0,2.3,5.,0.9,nh,reduce,5)
            flux_erosita(ierosita,7)=flux_erosita(ierosita,7)*reduce
            fluxU_erosita(ierosita,7)=fluxU_erosita(ierosita,7)*reduce
            fluxL_erosita(ierosita,7)=fluxL_erosita(ierosita,7)*reduce
            call nhdeabsorb2 (0,5.,8.,0.9,nh,reduce,5)
            flux_erosita(ierosita,8)=flux_erosita(ierosita,8)*reduce
            fluxU_erosita(ierosita,8)=fluxU_erosita(ierosita,8)*reduce
            fluxL_erosita(ierosita,8)=fluxL_erosita(ierosita,8)*reduce

            call fluxtofdens(0.9,0.2,2.3,flux_erosita(ierosita,1),1.,fdens,nudens)
            flux_erosita(ierosita,1)=fdens
            frequency_erosita(ierosita,1)=nudens
            call fluxtofdens(0.9,0.2,2.3,fluxU_erosita(ierosita,1),1.,fdens,nudens)
            fluxU_erosita(ierosita,1)=fdens
            call fluxtofdens(0.9,0.2,2.3,fluxL_erosita(ierosita,1),1.,fdens,nudens)
            fluxL_erosita(ierosita,1)=fdens
            call fluxtofdens(0.9,0.2,0.5,flux_erosita(ierosita,2),0.3,fdens,nudens)
            flux_erosita(ierosita,2)=fdens
            frequency_erosita(ierosita,2)=nudens
            call fluxtofdens(0.9,0.2,0.5,fluxU_erosita(ierosita,2),0.3,fdens,nudens)
            fluxU_erosita(ierosita,2)=fdens
            call fluxtofdens(0.9,0.2,0.5 ,fluxL_erosita(ierosita,2),0.3,fdens,nudens)
            fluxL_erosita(ierosita,2)=fdens
            call fluxtofdens(0.9,0.5,1.,flux_erosita(ierosita,3),0.7,fdens,nudens)
            flux_erosita(ierosita,3)=fdens
            frequency_erosita(ierosita,3)=nudens
            call fluxtofdens(0.9,0.5,1.,fluxU_erosita(ierosita,3),0.7,fdens,nudens)
            fluxU_erosita(ierosita,3)=fdens
            call fluxtofdens(0.9,0.5,1.,fluxL_erosita(ierosita,3),0.7,fdens,nudens)
            fluxL_erosita(ierosita,3)=fdens
            call fluxtofdens(0.9,1.,2.,flux_erosita(ierosita,4),1.5,fdens,nudens)
            flux_erosita(ierosita,4)=fdens
            frequency_erosita(ierosita,4)=nudens
            call fluxtofdens(0.9,1.,2.,fluxU_erosita(ierosita,4),1.5,fdens,nudens)
            fluxU_erosita(ierosita,4)=fdens
            call fluxtofdens(0.9,1.,2.,fluxL_erosita(ierosita,4),1.5,fdens,nudens)
            fluxL_erosita(ierosita,4)=fdens
            call fluxtofdens(0.9,2.,4.5,flux_erosita(ierosita,5),3.,fdens,nudens)
            flux_erosita(ierosita,5)=fdens
            frequency_erosita(ierosita,5)=nudens
            call fluxtofdens(0.9,2.,4.5,fluxU_erosita(ierosita,5),3.,fdens,nudens)
            fluxU_erosita(ierosita,5)=fdens
            call fluxtofdens(0.9,2.,4.5 ,fluxL_erosita(ierosita,5),3.,fdens,nudens)
            fluxL_erosita(ierosita,5)=fdens
            call fluxtofdens(0.9,0.5,2.,flux_erosita(ierosita,6),1.2,fdens,nudens)
            flux_erosita(ierosita,6)=fdens
            frequency_erosita(ierosita,6)=nudens
            call fluxtofdens(0.9,0.5,2.,fluxU_erosita(ierosita,6),1.2,fdens,nudens)
            fluxU_erosita(ierosita,6)=fdens
            call fluxtofdens(0.9,0.5,2.,fluxL_erosita(ierosita,6),1.2,fdens,nudens)
            fluxL_erosita(ierosita,6)=fdens
            call fluxtofdens(0.9,2.3,5.,flux_erosita(ierosita,7),4.,fdens,nudens)
            flux_erosita(ierosita,7)=fdens
            frequency_erosita(ierosita,7)=nudens
            call fluxtofdens(0.9,2.3,5.,fluxU_erosita(ierosita,7),4.,fdens,nudens)
            fluxU_erosita(ierosita,7)=fdens
            call fluxtofdens(0.9,2.3,5.,fluxL_erosita(ierosita,7),4.,fdens,nudens)
            fluxL_erosita(ierosita,7)=fdens
            call fluxtofdens(0.9,5.,8.,flux_erosita(ierosita,8),7.,fdens,nudens)
            flux_erosita(ierosita,8)=fdens
            frequency_erosita(ierosita,8)=nudens
            call fluxtofdens(0.9,5.,8.,fluxU_erosita(ierosita,8),7.,fdens,nudens)
            fluxU_erosita(ierosita,8)=fdens
            call fluxtofdens(0.9,5.,8.,fluxL_erosita(ierosita,8),7.,fdens,nudens)
            fluxL_erosita(ierosita,8)=fdens

            if (fluxL_erosita(ierosita,1) .le. 0) then
                flux_erosita(ierosita,1)=0
                fluxL_erosita(ierosita,1)=0
                fluxU_erosita(ierosita,1)=fluxU_erosita(ierosita,1)*3.
            endif
            if (fluxL_erosita(ierosita,2) .le. 0) then
                flux_erosita(ierosita,2)=0
                fluxL_erosita(ierosita,2)=0
                fluxU_erosita(ierosita,2)=fluxU_erosita(ierosita,2)*3.
            endif
            if (fluxL_erosita(ierosita,3) .le. 0) then
                flux_erosita(ierosita,3)=0
                fluxL_erosita(ierosita,3)=0
                fluxU_erosita(ierosita,3)=fluxU_erosita(ierosita,3)*3.
            endif
            if (fluxL_erosita(ierosita,4) .le. 0) then
                flux_erosita(ierosita,4)=0
                fluxL_erosita(ierosita,4)=0
                fluxU_erosita(ierosita,4)=fluxU_erosita(ierosita,4)*3.
            endif
            if (fluxL_erosita(ierosita,5) .le. 0) then
                flux_erosita(ierosita,5)=0
                fluxL_erosita(ierosita,5)=0
                fluxU_erosita(ierosita,5)=fluxU_erosita(ierosita,5)*3.
            endif
            if (fluxL_erosita(ierosita,6) .le. 0) then
                flux_erosita(ierosita,6)=0
                fluxL_erosita(ierosita,6)=0
                fluxU_erosita(ierosita,6)=fluxU_erosita(ierosita,6)*3.
            endif
            if (fluxL_erosita(ierosita,7) .le. 0) then
                flux_erosita(ierosita,7)=0
                fluxL_erosita(ierosita,7)=0
                fluxU_erosita(ierosita,7)=fluxU_erosita(ierosita,7)*3.
            endif
            if (fluxL_erosita(ierosita,8) .le. 0) then
                flux_erosita(ierosita,8)=0
                fluxL_erosita(ierosita,8)=0
                fluxU_erosita(ierosita,8)=fluxU_erosita(ierosita,8)*3.
            endif
c         write(*,*) flux_erosita(ierosita,1),fluxU_erosita(ierosita,1),fluxL_erosita(ierosita,1)
c         write(*,*) flux_erosita(ierosita,2),fluxU_erosita(ierosita,2),fluxL_erosita(ierosita,2)
c         write(*,*) flux_erosita(ierosita,8),fluxU_erosita(ierosita,8),fluxL_erosita(ierosita,8)

            CALL RXgraphic_code(flux_erosita(ierosita,1),'X',code)
            write (13,'(f9.5,2x,f9.5,2x,i6)') ra_erosita(ierosita),dec_erosita(ierosita),int(code)

         ELSE IF ((catalog(1:4) == '3fhl') .or. (catalog(1:7) == '4fgldr3') .or.
     &       (catalog(1:4) == '3fgl') .or. (catalog(1:5) == '2bigb') .or. (catalog(1:7) == 'f357cat')
     &           .or.  (catalog(1:5) == 'mst9y') .or. (catalog(1:5) == '2agile')
     &           .or.  (catalog(1:5) == 'fmev') .or. (catalog(1:4) == 'fgrb')) then
            igam=igam+1
            ra_gam(igam)=ra
            dec_gam(igam)=dec
            is=ie
            if ((catalog(1:5) == 'mst9y') .or. (catalog(1:4) == 'fgrb') .or. (catalog(1:7) == 'f357cat')) then
               ie=index(string(is+10:len(string)),' ')+is+9
            else
               ie=index(string(is+1:len(string)),',')+is
            endif
            if (catalog(1:5) == '2bigb') then
               ibigb=ibigb+1
cc               bigbind(igam)=MOD(ibigb,10)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') namegam(igam)
               namegam(igam)(6:lenact(namegam(igam))+6)=namegam(igam)(1:lenact(namegam(igam)))
               namegam(igam)(1:5)='2BIGB'
               if (ibigb .eq. 1) then
                  bigbind(igam)=1
               else
                  !write(*,*) namegam(igam),namegam(igam-1)
                  if (namegam(igam) == namegam(igam-1)) then
                     bigbind(igam)=bigbind(igam-1)+1
                  else
                     bigbind(igam)=1
                  endif
               endif
               !write(*,*) namegam(igam),ibigb,bigbind(igam)
            else if (catalog(1:4) == 'fgrb') then
               igrb=igrb+1
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') namegam(igam)
c               namegam(igam)(5:lenact(namegam(igam))+5)=namegam(igam)(1:lenact(namegam(igam)))
c               namegam(igam)(1:4)='GRB '
            else
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') namegam(igam)
            endif
            if ((bigbind(igam) .eq. 1) .or. (namegam(igam)(1:5) /= '2BIGB')) then
               if (namegam(igam)(1:3) =='GRB') then
                  write (11,'(f9.5,2x,f9.5,2x,i6,2x,a,a,a)') ra_gam(igam),dec_gam(igam),int(-2222),'"',trim(namegam(igam)),'"'
               else if (catalog(1:7) == 'f357cat') then
                  namegam(igam)(7:lenact(namegam(igam))+7)=namegam(igam)(1:lenact(namegam(igam)))
                  namegam(igam)(1:6)='Radio-'
                  write (11,'(f9.5,2x,f9.5,2x,i6,2x,a,a,a)') ra_gam(igam),dec_gam(igam),int(-3333),'"',trim(namegam(igam)),'"'
                  write (13,'(f9.5,2x,f9.5,2x,i6)') ra_gam(igam),dec_gam(igam),int(-3333)
               else
                  write (13,'(f9.5,2x,f9.5,2x,i6)') ra_gam(igam),dec_gam(igam),int(-1111)
                  write (11,'(f9.5,2x,f9.5,2x,i6,2x,a,a,a)') ra_gam(igam),dec_gam(igam),int(-1111),'"',trim(namegam(igam)),'"'
               endif
            endif
c            write(*,*) namegam(igam)
         ELSE
            iother=iother+1
            IF (iother > arrsize(4)) Stop 'Too many catalogued sources'
            ra_other(iother)=ra
            dec_other(iother)=dec
            is=ie
            if (catalog(1:5) == 'mquas') then
               ie=index(string(is+1:len(string)),',')+is
               read(string(is+1:ie-1),'(a)') name_other(iother)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               !read(string(is+1:ie-1),'(a)') classmq(iother)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) zz(iother)
            else
               ie=len(string)
               read(string(is+1:ie-1),'(a)') name_other(iother)
            endif
            if ((catalog(1:4) == 'mcxc') .or. (catalog(1:2) == 'zw') .or.
     &           (catalog(1:3) == 'whl')) then
               name_other(iother)(5:lenact(name_other(iother))+5)=name_other(iother)(1:lenact(name_other(iother)))
               name_other(iother)(1:4)=catalog(1:4)
            endif
            if (catalog(1:5) == 'mquas') then
               name_other(iother)(4:lenact(name_other(iother))+4)=name_other(iother)(1:lenact(name_other(iother)))
               name_other(iother)(1:3)='MQ '
            endif
            if (iother .ne. 1) then
               do i=1,iother-1
                  if (name_other(i) == name_other(i+1)) iother=iother-1
               enddo
            endif
            if (name_other(iother)(1:2) == 'MQ') then
c               write(*,*) name_other(iother),'   ',classmq(iother)
               write (13,'(f9.5,2x,f9.5,2x,i6)') ra_other(iother),dec_other(iother),int(-7777)
c               ra_other(iother) = -ra_other(iother)
            endif
         ENDIF
      ENDDO
 99   CONTINUE
      CLOSE (lu_in)
      if (aim .eq. 0) goto 501

      deallocate(vlasssrnm)
      allocate(rtype_source(arrsize(1)),t(arrsize(1)),track(arrsize(1)),ttsource(arrsize(1)),bary(arrsize(1)),rank(arrsize(1)),priority(arrsize(1)))
      allocate(ra_source(arrsize(1)),dec_source(arrsize(1)),ra_xx(arrsize(1)),dec_xx(arrsize(1)))
      allocate(xxerr(arrsize(1)),poserr_source(arrsize(1)),flux_source(arrsize(1)),xflux(arrsize(1)),rflux(arrsize(1)),rrconst(arrsize(1)))
      allocate(savemjy(arrsize(4)))
      allocate(spec_type(arrsize(2),arrsize(3)),spec_xpts(arrsize(2),arrsize(3)))
      allocate(ra_1kev(arrsize(2),arrsize(3)),dec_1kev(arrsize(2),arrsize(3)),distrx(arrsize(2),arrsize(3)))
      allocate(flux_1kev(arrsize(2),arrsize(3)),uflux_1kev(arrsize(2),arrsize(3)),lflux_1kev(arrsize(2),arrsize(3)),uflux_xpts(arrsize(2),arrsize(3)),lflux_xpts(arrsize(2),arrsize(3)),flux_xpts(arrsize(2),arrsize(3)),frequency_xpts(arrsize(2),arrsize(3)))
      allocate(poserr_1kev(arrsize(2),arrsize(3)),mjdstart(arrsize(2),arrsize(3)),mjdend(arrsize(2),arrsize(3)))

      CALL indexx (iradio,ra_radio,ra_index)
      !write(*,*) ra_radio(ra_index(1:iradio))
c      if (iradio .eq. 14) iradio=iradio-1
c      write(*,*) "Nr. of radio", iradio+1
      DO j=1,iradio
         DO i =0,5
           types(i) = 0
         ENDDO
         found = .FALSE.
         flux_x = 0.
         type_average = 0.
         ix = 0
         xpts=0
         k = ra_index(j)
         call chra(ra_radio(k),rah,ram,rasec,1)
         call chdec(dec_radio(k),id,dm,decsec,1)
         DO i=1,ixmm
c            if (xmm_type(i) == 2) min_dist_xmm=4./3600.
c            if (xmm_type(i) == 1) min_dist_xmm=15./3600.
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_xmm(i),dec_xmm(i),dist)
c            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
c               min_dist = sqrt(min_dist_xmm**2+(5./3600.)**2)
c            ELSE
            min_dist = sqrt(poserr_xmm(i)**2+poserr_radio(k)**2)/3600.
            !write(*,*) dist*3600.,min_dist*3600.,radio_type(k)
c            ENDIF
            IF (dist < max(min_dist,2./3600.)) THEN
               found = .TRUE.
               IF (xmm_type(i) == 1) THEN 
                 xray_type = 1
               ELSE IF (xmm_type(i) == 2) THEN
                 xray_type = 2
               ENDIF 
               !IF (flux_xmm(i,1) > 0. ) THEN
                  flux_x = flux_x + flux_xmm(i,1)
                  ix = ix +1
                  flux_1kev(ix,k)=flux_xmm(i,1)
                  uflux_1kev(ix,k)=FluxU_xmm(i,1)
                  lflux_1kev(ix,k)=FluxL_xmm(i,1)
                  ra_1kev(ix,k)=ra_xmm(i)
                  dec_1kev(ix,k)=dec_xmm(i)
                  poserr_1kev(ix,k)=poserr_xmm(i)
                  distrx(ix,k)=dist*3600.
                  spec_type(ix,k)=xray_type+50
                  mjdend(ix,k)=55000.
                  mjdstart(ix,k)=55000.
               !ENDIF
               if (xray_type == 2) then
                  Do l=2,6
                     xpts=xpts+1
                     flux_xpts(xpts,k)=flux_xmm(i,l)
                     uflux_xpts(xpts,k)=FluxU_xmm(i,l)
                     lflux_xpts(xpts,k)=FluxL_xmm(i,l)
                     frequency_xpts(xpts,k)=frequency_xmm(i,l)
                     spec_xpts(xpts,k)=xray_type
                  enddo
               else
                  do l=2,3
                     xpts=xpts+1
                     flux_xpts(xpts,k)=flux_xmm(i,l)
                     uflux_xpts(xpts,k)=FluxU_xmm(i,l)
                     lflux_xpts(xpts,k)=FluxL_xmm(i,l)
                     frequency_xpts(xpts,k)=frequency_xmm(i,l)
                     spec_xpts(xpts,k)=xray_type
                  ENDDO
               endif
               !if (flux_xmm(i,1) .le. 1.e-13) write(*,*) i,flux_xmm(i,1:6),dist*3600.
               if (ix .eq. 1) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),xray_type,
     &                             flux_xmm(i,1),const(k),ra_center,dec_center,source_type)
               else if ((i .eq. 1) .or. (spec_type(ix,k) .ne. spec_type(ix-1,k))) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),xray_type,
     &                             flux_xmm(i,1),const(k),ra_center,dec_center,source_type)
               endif
               IF (source_type .GE. 0) THEN 
                  types(source_type) = types(source_type) + 1
               ENDIF
            ENDIF
         ENDDO
         DO i=1,irosat
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_rosat(i),dec_rosat(i),dist)
c            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
            min_dist = sqrt(poserr_rosat(i)**2+poserr_radio(k)**2)/3600.
c            ELSE
c               min_dist = min_dist_rosat
c            ENDIF
            IF (dist < min_dist) THEN 
               found = .TRUE.
               IF (rosat_type(i) == 1) THEN 
                 xray_type = 3
               ELSE IF (rosat_type(i) == 2) THEN
                 xray_type = 4
               ENDIF 
               flux_x = flux_x + flux_rosat(i) 
               ix = ix +1
               flux_1kev(ix,k)=flux_rosat(i)
               uflux_1kev(ix,k)=FluxU_rosat(i)
               lflux_1kev(ix,k)=FluxL_rosat(i)
               ra_1kev(ix,k)=ra_rosat(i)
               dec_1kev(ix,k)=dec_rosat(i)
               poserr_1kev(ix,k)=poserr_rosat(i)
               distrx(ix,k)=dist*3600.
               spec_type(ix,k)=xray_type+50
               mjdend(ix,k)=55000.
               mjdstart(ix,k)=55000.
               if (ix .eq. 1) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_rosat(i),const(k),ra_center,dec_center,source_type)
               else if ((i .eq. 1) .or. (spec_type(ix,k) .ne. spec_type(ix-1,k))) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_rosat(i),const(k),ra_center,dec_center,source_type)
               endif
               IF (source_type .GE. 0) THEN 
                  types(source_type) = types(source_type) + 1
               ENDIF
            ENDIF
         ENDDO
         DO i=1,iswift
            CALL DIST_SKY(ra_radio(k),dec_radio(k),abs(ra_swift(i)),dec_swift(i),dist)
c            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
            min_dist = sqrt(poserr_swift(i)**2+poserr_radio(k)**2)/3600.
            !write(*,*) dist*3600.,min_dist*3600.,radio_type(k)
c            ELSE
c               min_dist = min_dist_swift
c            ENDIF
c            write(*,*) 'XRT',dist*3600.,min_dist*3600.
            IF (dist < max(min_dist,2./3600.)) THEN
               IF (xrt_type(i) == 1) THEN
                  xray_type = 5 !2sxps
               ELSE IF (xrt_type(i) == 2) THEN
                  xray_type = 9 !deep
               ELSE If (xrt_type(i) == 3) THEN
                  xray_type = 11 !ousx
               ENDIF
               found = .TRUE.
               flux_x = flux_x + flux_swift(i,1)
               ix = ix +1
               flux_1kev(ix,k)=flux_swift(i,1)
               uflux_1kev(ix,k)=FluxU_swift(i,1)
               lflux_1kev(ix,k)=FluxL_swift(i,1)
               ra_1kev(ix,k)=ra_swift(i)
               dec_1kev(ix,k)=dec_swift(i)
               poserr_1kev(ix,k)=poserr_swift(i)
               distrx(ix,k)=dist*3600.
               spec_type(ix,k)=xray_type+50
               if ((spec_type(ix,k) .eq. 61)) then
                  mjdstart(ix,k)=mjdst_swift(i)
                  mjdend(ix,k)=mjded_swift(i)
               else
                  mjdend(ix,k)=55000.
                  mjdstart(ix,k)=55000.
               endif
               do l=2,4
                  xpts=xpts+1
                  flux_xpts(xpts,k)=flux_swift(i,l)
                  uflux_xpts(xpts,k)=FluxU_swift(i,l)
                  lflux_xpts(xpts,k)=FluxL_swift(i,l)
                  frequency_xpts(xpts,k)=frequency_swift(i,l)
                  spec_xpts(xpts,k)=xray_type
               enddo
               if (frequency_swift(i,5) .ne. 999.) THEN
                  xpts=xpts+1
                  flux_xpts(xpts,k)=flux_swift(i,5)
                  uflux_xpts(xpts,k)=FluxU_swift(i,5)
                  lflux_xpts(xpts,k)=FluxL_swift(i,5)
                  frequency_xpts(xpts,k)=frequency_swift(i,5)
                  spec_xpts(xpts,k)=xray_type
               endif
               if (ix .eq. 1) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_swift(i,1),const(k),ra_center,dec_center,source_type)
               else if ((i .eq. 1) .or. (spec_type(ix,k) .ne. spec_type(ix-1,k))) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_swift(i,1),const(k),ra_center,dec_center,source_type)
               endif
               IF (source_type .GE. 0) THEN
                  types(source_type) = types(source_type) + 1
               ENDIF
            ENDIF
         ENDDO
         DO i=1,iipc
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_ipc(i),dec_ipc(i),dist)
c            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
            min_dist = sqrt(poserr_ipc(i)**2+poserr_radio(k)**2)/3600.
c            ELSE
c               min_dist = min_dist_ipc
c            ENDIF
            IF (dist < min_dist) THEN 
               found = .TRUE.
               if (ipc_type(i) == 1) then
                 xray_type = 12
               else
                 xray_type = 6
               endif
               flux_x = flux_x + flux_ipc(i) 
               ix = ix +1
               flux_1kev(ix,k)=flux_ipc(i)
               uflux_1kev(ix,k)=FluxU_ipc(i)
               lflux_1kev(ix,k)=FluxL_ipc(i)
               ra_1kev(ix,k)=ra_ipc(i)
               dec_1kev(ix,k)=dec_ipc(i)
               poserr_1kev(ix,k)=poserr_ipc(i)
               distrx(ix,k)=dist*3600.
               spec_type(ix,k)=xray_type+50
               mjdend(ix,k)=55000.
               mjdstart(ix,k)=55000.
               if (ix .eq. 1) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_ipc(i),const(k),ra_center,dec_center,source_type)
               else if ((i .eq. 1) .or. (spec_type(ix,k) .ne. spec_type(ix-1,k))) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_ipc(i),const(k),ra_center,dec_center,source_type)
               endif
               IF (source_type .GE. 0) THEN 
                  types(source_type) = types(source_type) + 1
               ENDIF
            ENDIF
         ENDDO
         DO i=1,ibmw
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_bmw(i),dec_bmw(i),dist)
c            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
            min_dist = sqrt(poserr_bmw(i)**2+poserr_radio(k)**2)/3600.
c            ELSE
c               min_dist = min_dist_bmw
c            ENDIF
            !write(*,*) dist,max(min_dist,2./3600.)
            IF (dist < max(min_dist,2./3600.)) THEN
               found = .TRUE.
               xray_type = 7
               flux_x = flux_x + flux_bmw(i)
               ix = ix +1
               flux_1kev(ix,k)=flux_bmw(i)
               uflux_1kev(ix,k)=FluxU_bmw(i)
               lflux_1kev(ix,k)=FluxL_bmw(i)
               ra_1kev(ix,k)=ra_bmw(i)
               dec_1kev(ix,k)=dec_bmw(i)
               poserr_1kev(ix,k)=poserr_bmw(i)
               distrx(ix,k)=dist*3600.
               spec_type(ix,k)=xray_type+50
               mjdend(ix,k)=55000.
               mjdstart(ix,k)=55000.
               if (ix .eq. 1) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_bmw(i),const(k),ra_center,dec_center,source_type)
               else if ((i .eq. 1) .or. (spec_type(ix,k) .ne. spec_type(ix-1,k))) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_bmw(i),const(k),ra_center,dec_center,source_type)
               endif
               IF (source_type .GE. 0) THEN 
                  types(source_type) = types(source_type) + 1
               ENDIF
            ENDIF
         ENDDO
         DO i=1,ichandra
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_chandra(i),dec_chandra(i),dist)
c            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
            min_dist = sqrt(poserr_chandra(i)**2+poserr_radio(k)**2)/3600.
c            ELSE
c            min_dist = min_dist_chandra
c            ENDIF
            IF (dist < max(min_dist,2./3600.)) THEN
               found = .TRUE.
               xray_type = 8
               flux_x = flux_x + flux_chandra(i,1)
               ix = ix +1
               flux_1kev(ix,k)=flux_chandra(i,1)
               uflux_1kev(ix,k)=FluxU_chandra(i,1)
               lflux_1kev(ix,k)=FluxL_chandra(i,1)
               ra_1kev(ix,k)=ra_chandra(i)
               dec_1kev(ix,k)=dec_chandra(i)
               poserr_1kev(ix,k)=poserr_chandra(i)
               distrx(ix,k)=dist*3600.
               spec_type(ix,k)=xray_type+50
               mjdend(ix,k)=55000.
               mjdstart(ix,k)=55000.
               if (ix .eq. 1) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                xray_type,flux_chandra(i,1),const(k),ra_center,dec_center,source_type)
               else if ((i .eq. 1) .or. (spec_type(ix,k) .ne. spec_type(ix-1,k))) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                xray_type,flux_chandra(i,1),const(k),ra_center,dec_center,source_type)
               endif
               IF (source_type .GE. 0) THEN
                  types(source_type) = types(source_type) + 1
               ENDIF
               do l=2,5
                  xpts=xpts+1
                  flux_xpts(xpts,k)=flux_chandra(i,l)
                  frequency_xpts(xpts,k)=frequency_chandra(i,l)
                  uflux_xpts(xpts,k)=FluxU_chandra(i,l)
                  lflux_xpts(xpts,k)=FluxL_chandra(i,l)
                  spec_xpts(xpts,k)=xray_type
               enddo
            ENDIF
         ENDDO
         Do i=1,imaxi
            CALL DIST_SKY(ra_radio(k),dec_radio(k),abs(ra_maxi(i)),dec_maxi(i),dist)
            min_dist = sqrt(poserr_maxi(i)**2+poserr_radio(k)**2)/3600.
            IF (dist < min_dist) THEN
               found = .TRUE.
               IF (maxi_type(i) == 1) THEN
                  xray_type = 10
               ELSE IF (maxi_type(i) == 2) THEN
                  xray_type = 13
               ENDIF
               flux_x = flux_x + flux_maxi(i,1)
               ix = ix +1
               flux_1kev(ix,k)=flux_maxi(i,1)
               uflux_1kev(ix,k)=FluxU_maxi(i,1)
               lflux_1kev(ix,k)=FluxL_maxi(i,1)
               ra_1kev(ix,k)=ra_maxi(i)
               dec_1kev(ix,k)=dec_maxi(i)
               poserr_1kev(ix,k)=poserr_maxi(i)
               distrx(ix,k)=dist*3600.
               spec_type(ix,k)=xray_type+50
               mjdend(ix,k)=55000.
               mjdstart(ix,k)=55000.
               if (ix .eq. 1) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                xray_type,flux_maxi(i,1),const(k),ra_center,dec_center,source_type)
               else if ((i .eq. 1) .or. (spec_type(ix,k) .ne. spec_type(ix-1,k))) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                xray_type,flux_maxi(i,1),const(k),ra_center,dec_center,source_type)
               endif
               IF (source_type .GE. 0) THEN
                  types(source_type) = types(source_type) + 1
               ENDIF
               if ( xray_type == 13 ) then
                  do l=2,2
                     xpts=xpts+1
                     flux_xpts(xpts,k)=flux_maxi(i,l)
                     frequency_xpts(xpts,k)=frequency_maxi(i,l)
                     uflux_xpts(xpts,k)=FluxU_maxi(i,l)
                     lflux_xpts(xpts,k)=FluxL_maxi(i,l)
                     spec_xpts(xpts,k)=xray_type
                  enddo
               ELSE
                  do l=2,4
                     xpts=xpts+1
                     flux_xpts(xpts,k)=flux_maxi(i,l)
                     frequency_xpts(xpts,k)=frequency_maxi(i,l)
                     uflux_xpts(xpts,k)=FluxU_maxi(i,l)
                     lflux_xpts(xpts,k)=FluxL_maxi(i,l)
                     spec_xpts(xpts,k)=xray_type
                  enddo
               endif
            ENDIF
         ENDDO
         !write(*,*) const

         DO i=1,ierosita
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_erosita(i),dec_erosita(i),dist)
            min_dist = sqrt(poserr_erosita(i)**2+poserr_radio(k)**2)/3600.
            IF (dist < max(min_dist,2./3600.)) THEN
               found = .TRUE.
               xray_type = 14
               flux_x = flux_x + flux_erosita(i,1)
               ix = ix +1
               flux_1kev(ix,k)=flux_erosita(i,1)
               uflux_1kev(ix,k)=fluxU_erosita(i,1)
               lflux_1kev(ix,k)=fluxL_erosita(i,1)
               ra_1kev(ix,k)=ra_erosita(i)
               dec_1kev(ix,k)=dec_erosita(i)
               poserr_1kev(ix,k)=poserr_erosita(i)
               distrx(ix,k)=dist*3600.
               spec_type(ix,k)=xray_type+50
               mjdend(ix,k)=55000.
               mjdstart(ix,k)=55000.
               if (ix .eq. 1) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                xray_type,flux_erosita(i,1),const(k),ra_center,dec_center,source_type)
               else if ((i .eq. 1) .or. (spec_type(ix,k) .ne. spec_type(ix-1,k))) then
                  CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                xray_type,flux_erosita(i,1),const(k),ra_center,dec_center,source_type)
               endif
               IF (source_type .GE. 0) THEN
                  types(source_type) = types(source_type) + 1
               ENDIF
               do l=2,8
                  xpts=xpts+1
                  flux_xpts(xpts,k)=flux_erosita(i,l)
                  frequency_xpts(xpts,k)=frequency_erosita(i,l)
                  uflux_xpts(xpts,k)=fluxU_erosita(i,l)
                  lflux_xpts(xpts,k)=fluxL_erosita(i,l)
                  spec_xpts(xpts,k)=xray_type
               enddo
            ENDIF
         ENDDO

cq      write(*,*) iradio
         IF (found) THEN 
            ifound = ifound +1
c            write(*,*) 'number of matched',ifound
            no_found = 0
            DO i = 0,5
               IF (types(i) > no_found) THEN
                  no_found=types(i)
                  type_average = i
               ENDIF
            ENDDO
            if (ifound .ne. 1) then
               do i = 1,ifound-1
                  call DIST_SKY(ra_radio(k),dec_radio(k),ra_radio(t(i)),dec_radio(t(i)),dist)
                  !if (dist*60 .lt. 0.8) then
                  if (dist*3600. .lt. 12.) then
                     rfound=rfound+1
                     IF ( ix.NE.0 ) THEN
                     flux_x = flux_x/float(ix)
                     ELSE
                     flux_x = 0.
                     ENDIF
                     if (rfound .eq. 1) then
                        track(ifound)=i !the source matched number
                     else
                        track(ifound)=track(i)
                     endif
                     if (flux_x .ne. xflux(track(ifound)) )  write(*,*) '!!!Warning, check X-ray counterpart!'
                     !rflux(track(ifound))=rflux(track(ifound))+flux_radio(k)
                     nrep(track(ifound))=nrep(track(ifound))+1
                     !rflux(track(ifound))=rflux(track(ifound))/nrep(track(ifound))
                     write(*,'(15x,"Repeated radio counterpart, ",f6.3,2x,"arcsec away from the matched nr.",2x,i2)')
     &                   dist*3600, track(ifound)
                     write(*,*) '     '
                     write(12,*) "===================="
                     write(12,'(i4,2x,a,2(2x,f10.5),2x,a,2x,i2)') track(ifound),"matched source",
     &                  ra_radio(k),dec_radio(k),"source type",type_average
                     write(12,'(4(es10.3,2x),2(f10.5,2x),f7.3,2(2x,f10.4),2x,i2)') frequency_radio(k),
     &                flux_radio(k),FluxU_radio(k),FluxL_radio(k),ra_radio(k),dec_radio(k),
     &                poserr_radio(k),mjdavg,mjdavg,radio_type(k)
                     t(ifound)=k
                     goto 97 !repeated radio counterpart, so end the loop
                  endif
               enddo
            endif
            sfound=sfound+1
            track(ifound)=sfound
            rflux(sfound)=flux_radio(k)
            print '(a,a,i4,6x,a,a,f9.5,a,f9.5)',  achar(27),'[31;1m Match nr.',sfound,achar(27),
     &           '[0m ra dec: ',ra_radio(k),',',dec_radio(k)
            IF ( ix.NE.0 ) THEN 
               flux_x = flux_x/float(ix)
            ELSE
               flux_x = 0.
            ENDIF
            xflux(sfound)=flux_x
            if (sfound .ne. 1 ) write(12,*) "===================="
c            if (sfound .ne. 1 ) write(17,*) "===================="
            write(12,'(i4,2x,a,2(2x,f10.5),2x,a,2x,i2)') sfound,"matched source",
     &         ra_radio(k),dec_radio(k),'source type',type_average
c            write(17,'(i4,2x,a,2(2x,f10.5),2x,a,2x,i2)') sfound,"matched source",
c     &         ra_radio(k),dec_radio(k),'source type',type_average
            write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_radio(k),flux_radio(k),
     &      FluxU_radio(k),FluxL_radio(k),ra_radio(k),dec_radio(k),poserr_radio(k),mjdavg,mjdavg,radio_type(k)
            do i=1,ix
               write(12,'(" 2.418E+17",3(2x,es10.3),2(2x,f10.5),2x,f8.3,2(2x,f10.4),2x,i2)')
     &            flux_1kev(i,k),uflux_1kev(i,k),lflux_1kev(i,k),ra_1kev(i,k),dec_1kev(i,k),
     &            poserr_1kev(i,k),mjdstart(i,k),mjdstart(i,k),spec_type(i,k)
c               if (spec_type(i,k) .eq. 61) THEN
c                  write(17,'(" 2.418E+17",3(2x,es10.3),2(2x,f10.4),2x,i2)')
c     &            flux_1kev(i,k),uflux_1kev(i,k),lflux_1kev(i,k),mjdstart(i,k),mjdend(i,k),spec_type(i,k)-50
c               endif
            enddo
            !write(*,*) 'how many other x-ray pts',xpts
            do i=1,xpts
               write(12,'(4(es10.3,2x),i2)') frequency_xpts(i,k),flux_xpts(i,k),
     &             uflux_xpts(i,k),lflux_xpts(i,k),spec_xpts(i,k)
            enddo
cccccccccccc check radio repeted!!!!!!!!!!!!!!!!!!!!
            !write(*,*) frequency_radio(k),flux_radio(k),radio_type(k)
            ra_source(sfound)=ra_radio(k)
            dec_source(sfound)=dec_radio(k)
            rtype_source(sfound)=radio_type(k)
            flux_source(sfound)=flux_radio(k)
            poserr_source(sfound)=poserr_radio(k)
            rrconst(sfound)=const(k)
            ttsource(sfound)=type_average
            bary(sfound)=0
 97         continue
            !write(*,*) "repeated sources",nrep(sfound),frequency_radio(k),flux_radio(k),radio_type(k)
            if (nrep(sfound) .gt. 1.) then
               !write(*,*) "repeated sources",rtype_source(sfound),radio_type(k)
               !write(*,*) flux_source(sfound),flux_radio(k)
               if (radio_type(k) .lt. rtype_source(sfound) ) then
                  !write(*,*) "replace type"
                  ra_source(sfound)=ra_radio(k)
                  dec_source(sfound)=dec_radio(k)
                  flux_source(sfound)=flux_radio(k)
                  poserr_source(sfound)=poserr_radio(k)
                  rrconst(sfound)=const(k)
               else if ((radio_type(k) .eq. rtype_source(sfound)) .and.
     &                      (flux_radio(k) .gt. flux_source(sfound)))then
                  !write(*,*) "replace flux"
                  ra_source(sfound)=ra_radio(k)
                  dec_source(sfound)=dec_radio(k)
                  flux_source(sfound)=flux_radio(k)
                  poserr_source(sfound)=poserr_radio(k)
                  rrconst(sfound)=const(k)
               endif
               goto 98  !skip
               !write(*,*) "Final radio",ra_source(sfound),dec_source(sfound),flux_source(sfound)
            endif
            write(*,*) '................Cataloged sources.................'
            catsrc=.false.
            DO i=1,iother
               CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_other(i),
     &                       dec_other(i),dist)
c               write(*,*) 'CHECK CAT.',ra_radio(k),dec_radio(k),ra_other(i),dec_other(i),dist*3600.
               IF (dist*3600. < max(poserr_radio(k),10.)) THEN !!!!!
                  IF (name_other(i)(1:4) == '3HSP') THEN
                     type_average = -1
                     write(*,'(2x,a,1x,a)') name_other(i)
                     ra_other(i) = -ra_other(i)
                     ra_cattemp=-ra_other(i)
                     dec_cattemp=dec_other(i)
                     if (nnsource(sfound)(1:6) == 'NONAME') nnsource(sfound)=name_other(i)
                     catsrc=.true.
                  ELSE IF (name_other(i)(1:3) == '5BZ') THEN
                     type_average = -2
                     write(*,'(2x,a,1x,a)') name_other(i)
                     ra_other(i) = -ra_other(i)
                     ra_cattemp=-ra_other(i)
                     dec_cattemp=dec_other(i)
                     if (nnsource(sfound)(1:6) == 'NONAME') nnsource(sfound)=name_other(i)
                     catsrc=.true.
                  ELSE IF ((name_other(i)(1:6) == 'CRATES') .or. (name_other(i)(1:4) == 'BROS')) THEN
                     type_average = -3
                     write(*,'(2x,a,1x,a)') name_other(i)
                     ra_other(i) = -ra_other(i)
                     if (nnsource(sfound)(1:6) == 'NONAME') nnsource(sfound)=name_other(i)
                  ELSE IF (name_other(i)(1:3) == 'PSR') THEN
                     type_average = -88
                     code=-8888
                     write(*,'(2x,a,1x,a)') name_other(i)
                     ra_other(i) = -ra_other(i)
                     if (nnsource(sfound)(1:6) == 'NONAME') nnsource(sfound)=name_other(i)
                  ElSE IF (name_other(i)(1:2) == 'MQ') THEN
                     type_average = -70
                     code=-7000
                     write(*,'(2x,a,1x,a)') name_other(i)
                     ra_other(i) = -ra_other(i)
                     zsource(sfound)=zz(i)
                     if (nnsource(sfound)(1:6) == 'NONAME') nnsource(sfound)=name_other(i)(4:lenact(name_other(i)))
                  ENDIF
               ENDIF
               IF (dist < min_dist_cluster) THEN
                  IF ( (name_other(i)(1:5) == 'ABELL') .OR.
     &                      (name_other(i)(1:4) == 'PSZ2') .OR.
     &                      (name_other(i)(1:4) == 'mcxc') .OR. !!!
     &                      (name_other(i)(1:5) == 'SWXCS') .OR.
     &                      (name_other(i)(1:2) == 'zw') .OR. !!!
     &                      (name_other(i)(1:3) == 'whl') ) THEN !!!
                    type_average = -4
                    write(*,'(2x,a,1x,a)') name_other(i)
                    if (nnsource(sfound)(1:6) == 'NONAME') nnsource(sfound)=name_other(i)
                  ENDIF
               ENDIF
               IF ((type_average .gt. -20) .and. (type_average .lt. 0)) THEN
                  CALL graphic_code (flux_x,flux_radio(k),type_average,code)
                  write(11,'(f9.5,2x,f9.5,2x,i6,2x,a,a,a)') abs(ra_other(i)),dec_other(i),int(code),'"',trim(name_other(i)),'"'
               ELSE if (code .eq. -8888) then
                  write(11,'(f9.5,2x,f9.5,2x,i6,2x,a,a,a)') abs(ra_other(i)),dec_other(i),int(code),'"',trim(name_other(i)),'"'
               ELSE if (code .eq. -7000) then
                  write(11,'(f9.5,2x,f9.5,2x,i6,2x,a,a,a)') abs(ra_other(i)),dec_other(i),int(code),'"',trim(name_other(i)(4:lenact(name_other(i)))),'"'
               ENDIF
               type_average=0
               code=0
            ENDDO
            if (nnsource(sfound)(1:6) == 'NONAME') then
               call chra(ra_radio(k),rah,ram,rasec,1)
               call chdec(dec_radio(k),id,dm,decsec,1)
               sign(1:1) ='+'
               if (dec < 0.0) sign(1:1)='-'
               write(nnsource(sfound),"(a,i2.2,i2.2,i2.2,a,i2.2,i2.2,i2.2)")
     &              "J",rah,ram,int(rasec),sign(1:1),abs(id),abs(dm),int(abs(decsec))
            endif

            t(ifound)=k !!!recourd the former index
cccccccccccccccccccccccc  check X-ray points cccccccccccccccccc
c            ra_xx(sfound)=0
c            dec_xx(sfound)=0
            rank(sfound)=0
c            xxerr(sfound)=0.
            totweight=0
            totxerr=0
            if (ix .gt. 1) rank(sfound)=rank(sfound)+1
            totweight=sum(1./poserr_1kev(1:ix,k))  !!!!temperately weighting...
c            write(*,*) 'WEIGHT TOTAL',sfound,totweight
            do i=1,ix
               ra_xx(sfound)=abs(ra_xx(sfound))+abs(ra_1kev(i,k))*((1./poserr_1kev(i,k))/totweight)
               dec_xx(sfound)=dec_xx(sfound)+dec_1kev(i,k)*((1./poserr_1kev(i,k))/totweight)
               totxerr=totxerr+poserr_1kev(i,k)
            enddo
            xxerr(sfound)=totxerr/float(ix)
c            write(*,*) 'CHECK X-ray',ra_xx(sfound),dec_xx(sfound)
            call DIST_SKY (ra_source(sfound),dec_source(sfound),ra_xx(sfound),dec_xx(sfound),dist)
            if ((xxerr(sfound) .lt. 15.) .and. (poserr_source(sfound) .lt. 15.)) then
               rank(sfound)=rank(sfound)+1
c               write(*,*) 'RANK ERROR',xxerr(sfound),poserr_source(sfound)
            endif
            if (dist*3600. .lt. 10.) then
               rank(sfound)=rank(sfound)+1
c               write(*,*) 'RANK DIST',dist*3600.
            endif
            if (((xxerr(sfound) .gt. 15.) .or. (poserr_source(sfound) .gt. 15.)) .and. (dist*3600. .lt. 15.)) then
               ra_xx(sfound)=ra_radio(k)
               dec_xx(sfound)=dec_radio(k)
               bary(sfound)=1
            endif
            if (catsrc) then
               rank(sfound)=rank(sfound)+1
               ra_xx(sfound)=ra_cattemp
               dec_xx(sfound)=dec_cattemp
               bary(sfound)=1
c               write(*,*) 'RANK CATS'
            endif
c            write(*,*) 'TEST X-ray position',ra_xx(sfound),dec_xx(sfound)
            write(*,*) '        '
         ELSE !!check radio without matched
            do i=1,iother
               if ( ( (name_other(i)(1:3) == '5BZ') .OR. (name_other(i)(1:4) == '3HSP') .or.
     &           (name_other(i)(1:6) == 'CRATES') .or. (name_other(i)(1:3) == 'PSR') .or.
     &           (name_other(i)(1:5) == 'mquas') .or. (name_other(i)(1:4) == 'BROS'))
     &            .AND. (ra_other(i) .gt. 0.) ) THEN
                  CALL DIST_SKY(ra_other(i),dec_other(i),ra_radio(k),dec_radio(k),dist)
                  if (dist*3600. .lt. max(poserr_radio(k),2.) ) found=.true.
               endif
            enddo
            if (.not. found) then
               CALL DIST_SKY(ra_center,dec_center,ra_radio(k),dec_radio(k),dist)
               if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
                  if (dist .le. errrad/60.) then
                     write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_radio(k),
     &                flux_radio(k),FluxU_radio(k),FluxL_radio(k),ra_radio(k),dec_radio(k),
     &                poserr_radio(k),mjdavg,mjdavg,-radio_type(k)
                  endif
               else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
                  if (dist .le. errmaj/60.) then
                     write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_radio(k),
     &               flux_radio(k),FluxU_radio(k),FluxL_radio(k),ra_radio(k),dec_radio(k),
     &               poserr_radio(k),mjdavg,mjdavg,-radio_type(k)
                  endif
               endif
            endif
         ENDIF
  98     continue
      ENDDO

c      if (iradio .eq. 14) iradio=iradio-1
c      write(*,*) "Nr. of radio", iradio+1
      deallocate(spec_type,spec_xpts)
      deallocate(ra_1kev,dec_1kev,distrx)
      deallocate(flux_1kev,uflux_1kev,lflux_1kev,uflux_xpts,lflux_xpts,flux_xpts,frequency_xpts)
      deallocate(poserr_1kev,mjdstart,mjdend)


      Do i=1,ixmm
         found=.false.
c         if (xmm_type(i) == 2) min_dist_xmm=4./3600.
c         if (xmm_type(i) == 1) min_dist_xmm=15./3600.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),ra_xmm(i),dec_xmm(i),dist)
            min_dist=sqrt(poserr_radio(j)**2+poserr_xmm(i)**2)/3600.
            if (dist .lt. max(min_dist,2./3600.))  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .OR. (name_other(j)(1:3) == 'PSR') .or.
     &      (name_other(j)(1:5) == 'mquas') .or. (name_other(j)(1:4) == 'BROS')) .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_xmm(i),dec_xmm(i),dist)
               if (dist*3600. .lt. max(poserr_xmm(i),10.)) found=.true.
             endif
         enddo
         if (.not. found) THEN
            if (xmm_type(i) == 1 ) then
               xray_type=1
            else
               xray_type=2
            endif
            CALL DIST_SKY(ra_center,dec_center,ra_xmm(i),dec_xmm(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_xmm(i,1),
     &            flux_xmm(i,1),FluxU_xmm(i,1),FluxL_xmm(i,1),ra_xmm(i),dec_xmm(i),
     &            poserr_xmm(i),mjdavg,mjdavg,xray_type+50
                  if (xray_type .eq. 1) then
                     do s=2,3
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_xmm(i,s),
     $                   flux_xmm(i,s),FluxU_xmm(i,s),FluxL_xmm(i,s),ra_xmm(i),dec_xmm(i),poserr_xmm(i),
     $                   mjdavg,mjdavg,xray_type
                     enddo
                  else
                     do s=2,6
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_xmm(i,s),
     $                   flux_xmm(i,s),FluxU_xmm(i,s),FluxL_xmm(i,s),ra_xmm(i),dec_xmm(i),
     &                   poserr_xmm(i),mjdavg,mjdavg,xray_type
                     enddo
                  endif
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_xmm(i,1),
     &            flux_xmm(i,1),FluxU_xmm(i,1),FluxL_xmm(i,1),ra_xmm(i),dec_xmm(i),poserr_xmm(i),
     &            mjdavg,mjdavg,xray_type+50
                  if (xray_type .eq. 1) then
                     do s=2,3
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_xmm(i,s),
     &                  flux_xmm(i,s),FluxU_xmm(i,s),FluxL_xmm(i,s),ra_xmm(i),dec_xmm(i),
     &                  poserr_xmm(i),mjdavg,mjdavg,xray_type
                     enddo
                  else
                     do s=2,6
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_xmm(i,s),
     &                     flux_xmm(i,s),FluxU_xmm(i,s),FluxL_xmm(i,s),ra_xmm(i),dec_xmm(i),
     &                     poserr_xmm(i),mjdavg,mjdavg,xray_type
                     enddo
                  endif
               endif
            endif
         endif
      enddo

      Do i=1,irosat
         found=.false.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),ra_rosat(i),dec_rosat(i),dist)
            min_dist=sqrt(poserr_radio(j)**2+poserr_rosat(i)**2)/3600.
            if (dist .lt. min_dist)  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(j)(1:3) == 'PSR') .or.
     &      (name_other(j)(1:5) == 'mquas') .or. (name_other(j)(1:4) == 'BROS') ) .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_rosat(i),dec_rosat(i),dist)
               if (dist*3600. .lt. max(poserr_rosat(i),10.)) found=.true.
            endif
         enddo
         if (.not. found) THEN
            if (rosat_type(i) == 1 ) then
               xray_type=3
            else
               xray_type=4
            endif
            CALL DIST_SKY(ra_center,dec_center,ra_rosat(i),dec_rosat(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_rosat(i),
     &             flux_rosat(i),FluxU_rosat(i),FluxL_rosat(i),ra_rosat(i),dec_rosat(i),
     &             poserr_rosat(i),mjdavg,mjdavg,xray_type+50
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_rosat(i),
     &             flux_rosat(i),FluxU_rosat(i),FluxL_rosat(i),ra_rosat(i),dec_rosat(i),
     &             poserr_rosat(i),mjdavg,mjdavg,xray_type+50
               endif
            endif
         endif
      enddo

      Do i=1,iswift
         found=.false.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),abs(ra_swift(i)),dec_swift(i),dist)
            min_dist=sqrt(poserr_radio(j)**2+poserr_swift(i)**2)/3600.
            if (dist .lt. max(min_dist,2./3600.))  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(j)(1:3) == 'PSR') .or.
     &      (name_other(j)(1:5) == 'mquas') .or. (name_other(j)(1:4) == 'BROS')) .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),abs(ra_swift(i)),dec_swift(i),dist)
               if (dist*3600. .lt. max(poserr_swift(i),10.)) found=.true.
            endif
         enddo
         if (.not. found) THEN
            if (xrt_type(i) == 1 ) then
               xray_type=5
            else if (xrt_type(i) == 3 ) then
               xray_type=11
            else
               xray_type=9
            endif
            CALL DIST_SKY(ra_center,dec_center,abs(ra_swift(i)),dec_swift(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  if ((xray_type .eq. 11) )then
                     write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(i,1),
     &                flux_swift(i,1),FluxU_swift(i,1),FluxL_swift(i,1),ra_swift(i),dec_swift(i),
     &                poserr_swift(i),mjdst_swift(i),mjded_swift(i),xray_type+50
                     do s=2,5
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(i,s),
     &                  flux_swift(i,s),FluxU_swift(i,s),FluxL_swift(i,s),abs(ra_swift(i)),dec_swift(i),
     &                  poserr_swift(i),mjdst_swift(i),mjded_swift(i),xray_type
                     enddo
                  else
                     write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(i,1),
     &                flux_swift(i,1),FluxU_swift(i,1),FluxL_swift(i,1),ra_swift(i),dec_swift(i),
     &                poserr_swift(i),mjdavg,mjdavg,xray_type+50
                     do s=2,4
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(i,s),
     &                  flux_swift(i,s),FluxU_swift(i,s),FluxL_swift(i,s),abs(ra_swift(i)),dec_swift(i),
     &                  poserr_swift(i),mjdavg,mjdavg,xray_type
                     enddo
                  endif
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  if ((xray_type .eq. 11)) then
                     write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(i,1),
     &               flux_swift(i,1),FluxU_swift(i,1),FluxL_swift(i,1),ra_swift(i),dec_swift(i),
     &               poserr_swift(i),mjdst_swift(i),mjded_swift(i),xray_type+50
                     do s=2,5
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(i,s),
     &                  flux_swift(i,s),FluxU_swift(i,s),FluxL_swift(i,s),abs(ra_swift(i)),dec_swift(i),
     &                  poserr_swift(i),mjdst_swift(i),mjded_swift(i),xray_type
                     enddo
                  else
                     write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(i,1),
     &               flux_swift(i,1),FluxU_swift(i,1),FluxL_swift(i,1),ra_swift(i),dec_swift(i),
     &                poserr_swift(i),mjdavg,mjdavg,xray_type+50
                     do s=2,4
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(i,s),
     &                  flux_swift(i,s),FluxU_swift(i,s),FluxL_swift(i,s),abs(ra_swift(i)),dec_swift(i),
     &                  poserr_swift(i),mjdavg,mjdavg,xray_type
                     enddo
                  endif
               endif
            endif
         endif
      enddo

      Do i=1,iipc
         found=.false.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),ra_ipc(i),dec_ipc(i),dist)
            min_dist=sqrt(poserr_radio(j)**2+poserr_ipc(i)**2)/3600.
            if (dist .lt. min_dist)  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(j)(1:3) == 'PSR') .or.
     &      (name_other(j)(1:5) == 'mquas') .or. (name_other(j)(1:4) == 'BROS')) .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_ipc(i),dec_ipc(i),dist)
               if (dist*3600. .lt. max(poserr_ipc(i),10.)) found=.true.
            endif
         enddo
         if (.not. found) THEN
            if (ipc_type(i) == 1) THEN
               xray_type=12
            else
               xray_type=6
            endif
            CALL DIST_SKY(ra_center,dec_center,ra_ipc(i),dec_ipc(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_ipc(i),
     &             flux_ipc(i),FluxU_ipc(i),FluxL_ipc(i),ra_ipc(i),dec_ipc(i),
     &             poserr_ipc(i),mjdavg,mjdavg,xray_type+50
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_ipc(i),
     &              flux_ipc(i),FluxU_ipc(i),FluxL_ipc(i),ra_ipc(i),dec_ipc(i),
     &              poserr_ipc(i),mjdavg,mjdavg,xray_type+50
               endif
            endif
         endif
      enddo

      Do i=1,ibmw
         found=.false.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),ra_bmw(i),dec_bmw(i),dist)
            min_dist=sqrt(poserr_radio(j)**2+poserr_bmw(i)**2)/3600.
            if (dist .lt. max(min_dist,2./3600.))  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(j)(1:3) == 'PSR') .or.
     &       (name_other(j)(1:5) == 'mquas') .or. (name_other(j)(1:4) == 'BROS')) .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_bmw(i),dec_bmw(i),dist)
               if (dist*3600. .lt. max(poserr_bmw(i),10.)) found=.true.
            endif
         enddo
         if (.not. found) THEN
            xray_type=7
            CALL DIST_SKY(ra_center,dec_center,ra_bmw(i),dec_bmw(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_bmw(i),
     &              flux_bmw(i),FluxU_bmw(i),FluxL_bmw(i),ra_bmw(i),dec_bmw(i),
     &              poserr_bmw(i),mjdavg,mjdavg,xray_type+50
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_bmw(i),
     &              flux_bmw(i),FluxU_bmw(i),FluxL_bmw(i),ra_bmw(i),dec_bmw(i),
     &              poserr_bmw(i),mjdavg,mjdavg,xray_type+50
               endif
            endif
         endif
      enddo

      Do i=1,ichandra
         found=.false.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),ra_chandra(i),dec_chandra(i),dist)
            min_dist=sqrt(poserr_radio(j)**2+poserr_chandra(i)**2)/3600.
            if (dist .lt. max(min_dist,2./3600.))  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(j)(1:3) == 'PSR') .or.
     &      (name_other(j)(1:5) == 'mquas') .or. (name_other(j)(1:4) == 'BROS')) .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_chandra(i),dec_chandra(i),dist)
               if (dist*3600. .lt. max(poserr_chandra(i),10.)) found=.true.
            endif
         enddo
         if (.not. found) THEN
            xray_type=8
            CALL DIST_SKY(ra_center,dec_center,ra_chandra(i),dec_chandra(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_chandra(i,1),
     &            flux_chandra(i,1),FluxU_chandra(i,1),FluxL_chandra(i,1),ra_chandra(i),dec_chandra(i),
     &            poserr_chandra(i),mjdavg,mjdavg,xray_type+50
                  do s=2,5
                     write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_chandra(i,s),
     &               flux_chandra(i,s),FluxU_chandra(i,s),FluxL_chandra(i,s),ra_chandra(i),dec_chandra(i),
     &               poserr_chandra(i),mjdavg,mjdavg,xray_type
                  enddo
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_chandra(i,1),
     &            flux_chandra(i,1),FluxU_chandra(i,1),FluxL_chandra(i,1),ra_chandra(i),dec_chandra(i),
     &            poserr_chandra(i),mjdavg,mjdavg,xray_type+50
                  do s=2,5
                     write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_chandra(i,s),
     &               flux_chandra(i,s),FluxU_chandra(i,s),FluxL_chandra(i,s),ra_chandra(i),dec_chandra(i),
     &               poserr_chandra(i),mjdavg,mjdavg,xray_type
                  enddo
               endif
            endif
         endif
      enddo

      Do i=1,imaxi
         found=.false.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),abs(ra_maxi(i)),dec_maxi(i),dist)
            min_dist=sqrt(poserr_radio(j)**2+poserr_maxi(i)**2)/3600.
            if (dist .lt. min_dist)  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(j)(1:3) == 'PSR') .or.
     &       (name_other(j)(1:5) == 'mquas') .or. (name_other(j)(1:4) == 'BROS')) .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),abs(ra_maxi(i)),dec_maxi(i),dist)
               if (dist*3600. .lt. max(poserr_maxi(i),10.)) found=.true.
            endif
         enddo
         if (.not. found) THEN
            IF (maxi_type(i) == 1) THEN
               xray_type = 10
            ELSE IF (maxi_type(i) == 2) THEN
               xray_type = 13
            ENDIF
            CALL DIST_SKY(ra_center,dec_center,abs(ra_maxi(i)),dec_maxi(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_maxi(i,1),
     &            flux_maxi(i,1),FluxU_maxi(i,1),FluxL_maxi(i,1),ra_maxi(i),dec_maxi(i),
     &            poserr_maxi(i),mjdavg,mjdavg,xray_type+50
                  if (xray_type .eq. 13) then
                     do s=2,2
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_maxi(i,s),
     &                  flux_maxi(i,s),FluxU_maxi(i,s),FluxL_maxi(i,s),abs(ra_maxi(i)),dec_maxi(i),
     &                  poserr_maxi(i),mjdavg,mjdavg,xray_type
                     enddo
                  else
                     do s=2,4
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_maxi(i,s),
     &                  flux_maxi(i,s),FluxU_maxi(i,s),FluxL_maxi(i,s),abs(ra_maxi(i)),dec_maxi(i),
     &                  poserr_maxi(i),mjdavg,mjdavg,xray_type
                     enddo
                  endif
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_maxi(i,1),
     &             flux_maxi(i,1),FluxU_maxi(i,1),FluxL_maxi(i,1),ra_maxi(i),dec_maxi(i),
     &             poserr_maxi(i),mjdavg,mjdavg,xray_type+50
                  if (xray_type .eq. 13) then
                     do s=2,2
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_maxi(i,s),
     &                   flux_maxi(i,s),FluxU_maxi(i,s),FluxL_maxi(i,s),abs(ra_maxi(i)),dec_maxi(i),
     &                   poserr_maxi(i),mjdavg,mjdavg,xray_type
                     enddo
                  else
                     do s=2,4
                        write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_maxi(i,s),
     &                   flux_maxi(i,s),FluxU_maxi(i,s),FluxL_maxi(i,s),abs(ra_maxi(i)),dec_maxi(i),
     &                   poserr_maxi(i),mjdavg,mjdavg,xray_type
                     enddo
                  endif
               endif
            endif
         endif
      enddo

      Do i=1,ierosita
         found=.false.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),ra_erosita(i),dec_erosita(i),dist)
            min_dist=sqrt(poserr_radio(j)**2+poserr_erosita(i)**2)/3600.
            if (dist .lt. max(min_dist,2./3600.))  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(j)(1:3) == 'PSR') .or.
     &      (name_other(j)(1:5) == 'mquas') .or. (name_other(j)(1:4) == 'BROS')) .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_erosita(i),dec_erosita(i),dist)
               if (dist*3600. .lt. max(poserr_erosita(i),10.)) found=.true.
            endif
         enddo
         if (.not. found) THEN
            xray_type=14
            CALL DIST_SKY(ra_center,dec_center,ra_erosita(i),dec_erosita(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_erosita(i,1),
     &            flux_erosita(i,1),FluxU_erosita(i,1),FluxL_erosita(i,1),ra_erosita(i),dec_erosita(i),
     &            poserr_erosita(i),mjdavg,mjdavg,xray_type+50
                  do s=2,8
                     write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_erosita(i,s),
     &               flux_erosita(i,s),FluxU_erosita(i,s),FluxL_erosita(i,s),ra_erosita(i),dec_erosita(i),
     &               poserr_erosita(i),mjdavg,mjdavg,xray_type
                  enddo
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_erosita(i,1),
     &            flux_erosita(i,1),FluxU_erosita(i,1),FluxL_erosita(i,1),ra_erosita(i),dec_erosita(i),
     &            poserr_erosita(i),mjdavg,mjdavg,xray_type+50
                  do s=2,8
                     write(14,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_erosita(i,s),
     &               flux_erosita(i,s),FluxU_erosita(i,s),FluxL_erosita(i,s),ra_erosita(i),dec_erosita(i),
     &               poserr_erosita(i),mjdavg,mjdavg,xray_type
                  enddo
               endif
            endif
         endif
      enddo

      !write(*,*) sfound,rfound,ifound
      if (ifound .ne. sfound+rfound ) stop 'Warning, might have wrong matched number'
c      write(*,*) 'RANK=',rank(1:sfound)
      Do i=1,sfound
         CALL graphic_code (xflux(i),flux_source(i)/rrconst(i),ttsource(i),code)
         if (bary(i) .eq. 0) then
            errfrx=(poserr_source(i))/(poserr_source(i)+xxerr(i))
            call DIST_SKY(ra_source(i),dec_source(i),ra_xx(i),dec_xx(i),dist)
c            write(*,*) ra_source(i),dec_source(i),ra_xx(i),dec_xx(i)
c            write(*,*) poserr_source(i),xxerr(i),dist,errfrx
            call int_great_circle(ra_source(i),dec_source(i),ra_xx(i),dec_xx(i),errfrx,dist,ra_bary,dec_bary)
         else if (bary(i) .eq. 1) then
            ra_bary=ra_xx(i)
            dec_bary=dec_xx(i)
         endif
c            write(*,*) ra_bary,dec_bary
         write(11,'(f9.5,2x,f9.5,2x,i6,2x,a,a,a)') ra_bary,dec_bary,int(code),'"',trim(nnsource(i)),'"'
         do j=1,sfound
            if (rank(j) .eq. 4) then
c               write(*,*) i,j
               rank(j)=-1
               priority(i)=j
               goto 300
            endif
         enddo
         do j=1,sfound
            if (rank(j) .eq. 3) then
c               write(*,*) i,j
               rank(j)=-1
               priority(i)=j
               goto 300
            endif
         enddo
         do j=1,sfound
            if (rank(j) .eq. 2) then
c              write(*,*) i,j
               rank(j)=-1
               priority(i)=j
               goto 300
            endif
         enddo
         do j=1,sfound
            if (rank(j) .eq. 1) then
c               write(*,*) i,j
               rank(j)=-1
               priority(i)=j
               goto 300
            endif
         enddo
         do j=1,sfound
            if (rank(j) .ge. 0) then
c               write(*,*) rank(j),i,j
               rank(j)=-1
               priority(i)=j
               goto 300
            endif
         enddo
300      continue
      enddo
      IF (ifound  ==  0) print *,achar(27),'[31;1m No radio/X-ray matches were found.',achar(27),'[0m'

      savemjy(1:iother)=0.
      DO l=1,iother
         IF ( ( (name_other(l)(1:3) == '5BZ') .OR. (name_other(l)(1:4) == '3HSP') .or.
     &  (name_other(l)(1:6) == 'CRATES') .or. (name_other(l)(1:3) == 'PSR') .or.
     &  (name_other(l)(1:5) == 'mquas') .or. (name_other(l)(1:4) == 'BROS')) .AND. (ra_other(l) .gt. 0.)) THEN
            ncat=ncat+1
            ra_cat(ncat)=ra_other(l)
            dec_cat(ncat)=dec_other(l)
            name_cat(ncat)=name_other(l)
            if (ncat .ne. 1) then ! check repeat cataloged sources
               do j=1,ncat-1
                  call DIST_SKY (ra_other(l),dec_other(l),ra_cat(j),dec_cat(j),dist)
                  if (dist*3600. .lt. 6. ) then
                     track2(ncat)=track2(j)
                     IF (name_other(l)(1:3) == '5BZ') type_average = -2
                     IF (name_other(l)(1:4) == '3HSP') type_average = -1
                     IF ((name_other(l)(1:6) == 'CRATES').or. (name_other(l)(1:4) == 'BROS')) type_average = -3
                     IF (name_other(l)(1:3) == 'PSR') then
                        type_average = -88
                        code=-8888
                     endif
c                     IF (name_other(l)(1:5) == 'mquas') then
c                        type_average = -70
c                        code=-7000
c                     endif
                     write(*,'(a,a,i4,2x,a)') name_other(l)(1:lenact(name_other(l))),
     &                   ", repeated with candidate nr.",track2(ncat),name_cat(j)(1:lenact(name_cat(j)))
                     goto 100
                  endif
               enddo
            endif
            IF (name_other(l)(1:3) == '5BZ') type_average = -6
            IF (name_other(l)(1:4) == '3HSP') type_average = -5
            IF ((name_other(l)(1:6) == 'CRATES') .or. (name_other(l)(1:4) == 'BROS')) type_average = -7
            IF (name_other(l)(1:3) == 'PSR') type_average = -99
c            IF (name_other(l)(1:5) == 'mquas') type_average = -77
            sfound=sfound+1
            track2(ncat)=sfound
            if (sfound .ne. 1 ) write(12,*) "===================="
            write(12,'(i4,2x,a,2(2x,f10.5),2x,a,2x,i2)') sfound,"matched source",
     &         ra_other(l),dec_other(l),'source type',type_average
c            write(17,'(i4,2x,a,2(2x,f10.5),2x,a,2x,i2)') sfound,"matched source",
c     &         ra_other(l),dec_other(l),'source type',type_average
            CALL DIST_SKY(ra_other(l),dec_other(l),ra_center,dec_center,dist)
            dist = dist*60
            if (type_average .eq. -7) then
               print '(a,a,i4,a,a,a,a,a,f7.3,a)', achar(27),'[35;1m Candidate nr.',sfound,
     &          ', Known flat spectrum radio source with no radio/X-ray match: ',
     &          achar(27),'[0m',name_other(l)(1:lenact(name_other(l))),' found at a distance of ',
     &                          real(dist),' arcmin '
            else if (type_average .eq. -99) then
               code=-9999
               !write(*,*) code,type_average
               write(*,'("Pulsar",2x,a,2x,f7.3,2x,"arcmin away")') name_other(l)(1:lenact(name_other(l))),dist
c            else if (type_average .eq. -77) then
c               code=-7777
c               write(*,'("Million Quasar",2x,a,2x,f7.3,2x,"arcmin away")') name_other(l)(1:lenact(name_other(l))),dist
            else
               print '(a,a,i4,a,a,a,a,a,f7.3,a)', achar(27),'[35;1m Candidate nr.',sfound,
     &          ', Known blazar with no radio/X-ray match: ',
     &          achar(27),'[0m',name_other(l)(1:lenact(name_other(l))),' found at a distance of ',
     &                          real(dist),' arcmin '
            endif
            do j=1,iradio
               call DIST_SKY(ra_other(l),dec_other(l),ra_radio(j),dec_radio(j),dist)
               if (dist*3600. < max(poserr_radio(j),10.)) THEN
                  write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_radio(j),
     &            flux_radio(j),FluxU_radio(j),FluxL_radio(j),ra_radio(j),dec_radio(j),
     &            poserr_radio(j),mjdavg,mjdavg,radio_type(j)
                  savemjy(ncat)=flux_radio(j)/const(j)
               endif
            enddo
            do j=1,ixmm
c               if (xmm_type(i) == 1) min_dist_xmm=15./3600.
c               if (xmm_type(i) == 2) min_dist_xmm=4./3600.
               call DIST_SKY(ra_other(l),dec_other(l),ra_xmm(j),dec_xmm(j),dist)
               if (dist*3600. < max(poserr_xmm(j),10.)) THEN
                  if (xmm_type(j) == 1) then
                     xray_type=1
                     write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_xmm(j,1),
     &                flux_xmm(j,1),FluxU_xmm(j,1),FluxL_xmm(j,1),ra_xmm(j),dec_xmm(j),
     &                poserr_xmm(j),mjdavg,mjdavg,xray_type+50
                     do s=2,3
                        write(12,'(4(es10.3,2x),i2)') frequency_xmm(j,s),flux_xmm(j,s),FluxU_xmm(j,s),
     &                     FluxL_xmm(j,s),xray_type
                     enddo
                  else
                     xray_type=2
                     write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_xmm(j,1),
     &               flux_xmm(j,1),FluxU_xmm(j,1),FluxL_xmm(j,1),ra_xmm(j),dec_xmm(j),
     &               poserr_xmm(j),mjdavg,mjdavg,xray_type+50
                     do s=2,6
                        write(12,'(4(es10.3,2x),i2)') frequency_xmm(j,s),flux_xmm(j,s),FluxU_xmm(j,s),
     &                     FluxL_xmm(j,s),xray_type
                     enddo
                  endif
               endif
            enddo
            do j=1,irosat
               call DIST_SKY(ra_other(l),dec_other(l),ra_rosat(j),dec_rosat(j),dist)
               if (dist*3600. < max(poserr_rosat(j),10.)) THEN
                  if (rosat_type(j) == 1) THEN
                     xray_type=3
                     write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_rosat(j),
     &                   flux_rosat(j),FluxU_rosat(j),FluxL_rosat(j),ra_rosat(j),dec_rosat(j),
     &                   poserr_rosat(j),mjdavg,mjdavg,xray_type+50
                  else
                     xray_type=4
                     write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_rosat(j),
     &                  flux_rosat(j),FluxU_rosat(j),FluxL_rosat(j),ra_rosat(j),dec_rosat(j),
     &                   poserr_rosat(j),mjdavg,mjdavg,xray_type+50
                  endif
               endif
            enddo
            do j=1,iswift
               call DIST_SKY(ra_other(l),dec_other(l),abs(ra_swift(j)),dec_swift(j),dist)
               if (dist*3600. < max(poserr_swift(j),10.)) THEN
                  if (xrt_type(j) == 1) then
                     xray_type=5
                  else if (xrt_type(j) == 3) then
                     xray_type=11
                  else
                     xray_type=9
                  endif
                  if ((xray_type .eq. 11))then
c                     write(17,'(" 2.418E+17",3(2x,es10.3),2(2x,f10.4),2x,i2)')
c     &            flux_swift(j,1),FluxU_swift(j,1),FluxL_swift(j,1),mjdst_swift(j),mjded_swift(j),xray_type
                     write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(j,1),
     &                     flux_swift(j,1),FluxU_swift(j,1),FluxL_swift(j,1),ra_swift(j),dec_swift(j),
     &                     poserr_swift(j),mjdst_swift(j),mjded_swift(j),xray_type+50
                  else
                     write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(j,1),
     &                     flux_swift(j,1),FluxU_swift(j,1),FluxL_swift(j,1),ra_swift(j),dec_swift(j),
     &                     poserr_swift(j),mjdavg,mjdavg,xray_type+50
                  endif
                  do s=2,4
                     write(12,'(4(es10.3,2x),i2)') frequency_swift(j,s),flux_swift(j,s),FluxU_swift(j,s),
     &                     FluxL_swift(j,s),xray_type
                  enddo
                  if (frequency_swift(j,5) .ne. 999.) then
                     write(12,'(4(es10.3,2x),i2)') frequency_swift(j,5),flux_swift(j,5),
     &              FluxU_swift(j,5),FluxL_swift(j,5),xray_type
                  endif
               endif
            enddo
            do j=1,iipc
               call DIST_SKY(ra_other(l),dec_other(l),ra_ipc(j),dec_ipc(j),dist)
               if (dist*3600. < max(poserr_ipc(j),10.)) THEN
                  if (ipc_type(j) == 1) then
                     xray_type=12
                  else
                     xray_type=6
                  endif
                  write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_ipc(j),flux_ipc(j),
     &             FluxU_ipc(j),FluxL_ipc(j),ra_ipc(j),dec_ipc(j),poserr_ipc(j),mjdavg,mjdavg,xray_type+50
               endif
            enddo
            do j=1,ibmw
               call DIST_SKY(ra_other(l),dec_other(l),ra_bmw(j),dec_bmw(j),dist)
               if (dist*3600. < max(poserr_bmw(j),10.)) THEN
                  xray_type=7
                  write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_bmw(j),flux_bmw(j),
     &              FluxU_bmw(j),FluxL_bmw(j),ra_bmw(j),dec_bmw(j),poserr_bmw(j),mjdavg,mjdavg,xray_type+50
               endif
            enddo
            do j=1,ichandra
               call DIST_SKY(ra_other(l),dec_other(l),ra_chandra(j),dec_chandra(j),dist)
               if (dist*3600. < max(poserr_chandra(j),10.)) THEN
                  xray_type=8
                  write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_chandra(j,1),
     &                 flux_chandra(j,1),FluxU_chandra(j,1),FluxL_chandra(j,1),ra_chandra(j),dec_chandra(j),
     &                   poserr_chandra(j),mjdavg,mjdavg,xray_type+50
                  do s=2,5
                     write(12,'(4(es10.3,2x),i2)') frequency_chandra(j,s),flux_chandra(j,s),FluxU_chandra(j,s),
     &                     FluxL_chandra(j,s),xray_type
                  enddo
               endif
            enddo
            do j=1,imaxi
               call DIST_SKY(ra_other(l),dec_other(l),abs(ra_maxi(j)),dec_maxi(j),dist)
               if (dist*3600. < max(poserr_maxi(j),10.)) THEN
                  IF (maxi_type(i) == 1) THEN
                     xray_type = 10
                  ELSE IF (maxi_type(i) == 2) THEN
                     xray_type = 13
                  ENDIF
                  write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_maxi(j,1),
     &                 flux_maxi(j,1),FluxU_maxi(j,1),FluxL_maxi(j,1),ra_maxi(j),dec_maxi(j),
     &                   poserr_maxi(j),mjdavg,mjdavg,xray_type+50
                  if (xray_type .eq. 13) then
                     do s=2,2
                        write(12,'(4(es10.3,2x),i2)') frequency_maxi(j,s),flux_maxi(j,s),FluxU_maxi(j,s),
     &                     FluxL_maxi(j,s),xray_type
                     enddo
                  ELSE
                     do s=2,4
                        write(12,'(4(es10.3,2x),i2)') frequency_maxi(j,s),flux_maxi(j,s),FluxU_maxi(j,s),
     &                     FluxL_maxi(j,s),xray_type
                     enddo
                  endif
               endif
            enddo
            do j=1,ierosita
               call DIST_SKY(ra_other(l),dec_other(l),ra_erosita(j),dec_erosita(j),dist)
               if (dist*3600. < max(poserr_erosita(j),10.)) THEN
                  xray_type=14
                  write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_erosita(j,1),
     &                 flux_erosita(j,1),FluxU_erosita(j,1),FluxL_erosita(j,1),ra_erosita(j),dec_erosita(j),
     &                   poserr_erosita(j),mjdavg,mjdavg,xray_type+50
                  do s=2,8
                     write(12,'(4(es10.3,2x),i2)') frequency_erosita(j,s),flux_erosita(j,s),FluxU_erosita(j,s),
     &                     FluxL_erosita(j,s),xray_type
                  enddo
               endif
            enddo
100   continue
            if ((type_average .lt. -3) .and. (type_average .gt. -20)) then
!for no X-ray blazars and CRATES sources
                if (savemjy(ncat) .gt. 0.) then
                   if (savemjy(ncat) .lt. 20.) savemjy(ncat)=20.
                   call RXgraphic_code(savemjy(ncat),'R',code)
                   write(11,'(f9.5,2x,f9.5,2x,i6,2x,a,a,a)') ra_other(l),dec_other(l),int(code+60000),'"',trim(name_other(l)),'"' !radio source to CRATES
                endif
            endif
            if (type_average .gt. -20) CALL graphic_code (1.,1.,type_average,code) !produce code
            !if (type_average .eq. -7) CALL graphic_code (1.,savemjy(l),type_average,code)
            write(11,'(f9.5,2x,f9.5,2x,i6,2x,a,a,a)') ra_other(l),dec_other(l),int(code),'"',trim(name_other(l)),'"'
         ENDIF
      ENDDO
      WRITE (*,*) '      '
      if (sfound .gt. 0) write(12,*) "===================="
      if (sfound .gt. 0) write(12,'(i4,2x,a)') 1000,"Redshift"
      if (sfound .gt. 0) write(12,*) zsource(1:sfound)

      do i=1,igam
         if (namegam(i)(1:3) == 'GRB') then
            CALL DIST_SKY(ra_gam(i),dec_gam(i),ra_center,dec_center,dist)
            dist = dist*60
            write(*,'(a,2x,f7.3,2x,"arcmin away")') namegam(i)(1:lenact(namegam(i))),dist
         endif
      enddo
      WRITE (*,*) '      '
      if (igam-igrb .gt. 0) then
         write(*,*) 'Gamma-ray Counterparts'
         do i=1,igam
            if ((bigbind(i) .eq. 1) .or. (namegam(i)(1:5) /= '2BIGB')) then
               if ((namegam(i)(1:3) /= 'GRB') .and. (namegam(i)(1:5) /= 'Radio')) write(*,*) namegam(i),ra_gam(i),dec_gam(i)
            endif
         enddo
      endif
      WRITE (*,*) '      '
      !write(*,*) 'Suggest Priority Order of Source nr.:',priority(1:ifound-rfound)
      !WRITE (*,*) '      '

      deallocate(ra_other,dec_other)
      deallocate(name_other)
      deallocate(savemjy,zz)
      deallocate(ra_gam,dec_gam,ra_cat,dec_cat)
      deallocate(name_cat,namegam,track2,bigbind)
      deallocate(rtype_source,nrep,t,track,ttsource,bary,rank,priority)
      deallocate(ra_source,dec_source,ra_xx,dec_xx)
      deallocate(xxerr,poserr_source,flux_source,xflux,rflux,rrconst,zsource)
      deallocate(nnsource)

      goto 502

cccccc for skip the phase 1
501   continue
      sfound=1
      type_average=99
      write(11,'(f9.5,2x,f9.5,2x,a)') ra_center,dec_center,"99"
      write(12,'(i4,2x,a,2(2x,f10.5),2x,a,2x,i2)') sfound,"matched source",
     &         ra_center,dec_center,'source type',type_average
c      write(17,'(i4,2x,a,2(2x,f10.5),2x,a,2x,i2)') sfound,"matched source",
c     &         ra_center,dec_center,'source type',type_average
      do j=1,iradio
         call DIST_SKY(ra_center,dec_center,ra_radio(j),dec_radio(j),dist)
         !write(*,*) 'RADIO',ra_center,dec_center,ra_radio(j),dec_radio(j),dist*3600.,poserr_radio(j)
         if (dist*3600. .lt. max(poserr_radio(j),2.) ) then !18 arcsec for radio sources
            write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_radio(j),flux_radio(j),
     &      FluxU_radio(j),FluxL_radio(j),ra_radio(j),dec_radio(j),poserr_radio(j),mjdavg,mjdavg,radio_type(j)
         endif
      enddo
      do j=1,ixmm
c         if (xmm_type(j) == 1) min_dist_xmm=15./3600.
c         if (xmm_type(j) == 2) min_dist_xmm=4./3600.
         call DIST_SKY(ra_center,dec_center,ra_xmm(j),dec_xmm(j),dist)
c         write(*,*) 'XMM',ra_center,dec_center,ra_xmm(j),dec_xmm(j),dist*3600.,poserr_xmm(j)
         if ( dist*3600. .lt. max(poserr_xmm(j),2.) ) then
            if (xmm_type(j) == 1) then
               xray_type=1
               write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_xmm(j,1),flux_xmm(j,1),
     &         FluxU_xmm(j,1),FluxL_xmm(j,1),ra_xmm(j),dec_xmm(j),poserr_xmm(j),mjdavg,mjdavg,xray_type+50
               do s=2,3
                  write(12,'(4(es10.3,2x),i2)') frequency_xmm(j,s),flux_xmm(j,s),FluxU_xmm(j,s),
     &                     FluxL_xmm(j,s),xray_type
               enddo
            else
               xray_type=2
               write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_xmm(j,1),flux_xmm(j,1),
     &         FluxU_xmm(j,1),FluxL_xmm(j,1),ra_xmm(j),dec_xmm(j),poserr_xmm(j),mjdavg,mjdavg,xray_type+50
               do s=2,6
                  write(12,'(4(es10.3,2x),i2)') frequency_xmm(j,s),flux_xmm(j,s),FluxU_xmm(j,s),
     &                     FluxL_xmm(j,s),xray_type
               enddo
            endif
         endif
      enddo
      do j=1,irosat
         call DIST_SKY(ra_center,dec_center,ra_rosat(j),dec_rosat(j),dist)
         if ( dist*3600. .lt. poserr_rosat(j) ) then
            if (rosat_type(j) == 1) THEN
               xray_type=3
               write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_rosat(j),
     &           flux_rosat(j),FluxU_rosat(j),FluxL_rosat(j),ra_rosat(j),dec_rosat(j),
     &           poserr_rosat(j),mjdavg,mjdavg,xray_type+50
            else
               xray_type=4
               write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_rosat(j),
     &            flux_rosat(j),FluxU_rosat(j),FluxL_rosat(j),ra_rosat(j),dec_rosat(j),
     &            poserr_rosat(j),mjdavg,mjdavg,xray_type+50
            endif
         endif
      enddo
      do j=1,iswift
         call DIST_SKY(ra_center,dec_center,abs(ra_swift(j)),dec_swift(j),dist)
c         write(*,*) 'SXPS',ra_center,dec_center,ra_swift(j),dec_swift(j),dist*3600.,poserr_swift(j)
         if ( dist*3600. .lt. max(poserr_swift(j),2.)) then
            if (xrt_type(j) == 1) THEN
               xray_type=5
            else if (xrt_type(j) == 3) THEN
               xray_type=11
            else
               xray_type=9
            endif
            if ((xray_type .eq. 11)) then
               write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(j,1),
     &                     flux_swift(j,1),FluxU_swift(j,1),FluxL_swift(j,1),ra_swift(j),dec_swift(j),
     &                     poserr_swift(j),mjdst_swift(j),mjded_swift(j),xray_type+50
c               write(17,'(" 2.418E+17",3(2x,es10.3),2(2x,f10.4),2x,i2)')
c     &            flux_swift(j,1),FluxU_swift(j,1),FluxL_swift(j,1),mjdst_swift(j),mjded_swift(j),xray_type
            else
               write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_swift(j,1),
     &                     flux_swift(j,1),FluxU_swift(j,1),FluxL_swift(j,1),ra_swift(j),dec_swift(j),
     &                     poserr_swift(j),mjdavg,mjdavg,xray_type+50
            endif
            do s=2,4
               write(12,'(4(es10.3,2x),i2)') frequency_swift(j,s),flux_swift(j,s),FluxU_swift(j,s),
     &                     FluxL_swift(j,s),xray_type
            enddo
            if (frequency_swift(j,5) .ne. 999.) then
               write(12,'(4(es10.3,2x),i2)') frequency_swift(j,5),flux_swift(j,5),
     &              FluxU_swift(j,5),FluxL_swift(j,5),xray_type
            endif
         endif
      enddo
      do j=1,iipc
         call DIST_SKY(ra_center,dec_center,ra_ipc(j),dec_ipc(j),dist)
         if ( dist*3600. .lt. poserr_ipc(j) ) then
            if (ipc_type(j) == 1) then
               xray_type=12
            else
               xray_type=6
            endif
            write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_ipc(j),flux_ipc(j),
     &      FluxU_ipc(j),FluxL_ipc(j),ra_ipc(j),dec_ipc(j),poserr_ipc(j),mjdavg,mjdavg,xray_type+50
         endif
      enddo
      do j=1,ibmw
         call DIST_SKY(ra_center,dec_center,ra_bmw(j),dec_bmw(j),dist)
         if ( dist*3600. .lt. max(poserr_bmw(j),2.) ) then
            xray_type=7
            write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_bmw(j),flux_bmw(j),
     &      FluxU_bmw(j),FluxL_bmw(j),ra_bmw(j),dec_bmw(j),poserr_bmw(j),mjdavg,mjdavg,xray_type+50
         endif
      enddo
      do j=1,ichandra
         call DIST_SKY(ra_center,dec_center,ra_chandra(j),dec_chandra(j),dist)
         if ( dist*3600. .lt. max(poserr_chandra(j),2.) ) then
            xray_type=8
            write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_chandra(j,1),
     &                 flux_chandra(j,1),FluxU_chandra(j,1),FluxL_chandra(j,1),ra_chandra(j),dec_chandra(j),
     &                   poserr_chandra(j),mjdavg,mjdavg,xray_type+50
            do s=2,5
               write(12,'(4(es10.3,2x),i2)') frequency_chandra(j,s),flux_chandra(j,s),FluxU_chandra(j,s),
     &                     FluxL_chandra(j,s),xray_type
            enddo
         endif
      enddo
      do j=1,imaxi
         call DIST_SKY(ra_center,dec_center,abs(ra_maxi(j)),dec_maxi(j),dist)
         if ( dist*3600. .lt. poserr_maxi(j) ) then
            IF (maxi_type(i) == 1) THEN
               xray_type = 10
            ELSE IF (maxi_type(i) == 2) THEN
               xray_type = 13
            ENDIF
            write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_maxi(j,1),
     &                 flux_maxi(j,1),FluxU_maxi(j,1),FluxL_maxi(j,1),ra_maxi(j),dec_maxi(j),
     &                   poserr_maxi(j),mjdavg,mjdavg,xray_type+50
            if (xray_type .eq. 13) then
               do s=2,2
                  write(12,'(4(es10.3,2x),i2)') frequency_maxi(j,s),flux_maxi(j,s),FluxU_maxi(j,s),
     &                     FluxL_maxi(j,s),xray_type
               enddo
            ELSE
               do s=2,4
                  write(12,'(4(es10.3,2x),i2)') frequency_maxi(j,s),flux_maxi(j,s),FluxU_maxi(j,s),
     &                     FluxL_maxi(j,s),xray_type
               enddo
            endif
         endif
      enddo
      do j=1,ierosita
         call DIST_SKY(ra_center,dec_center,ra_erosita(j),dec_erosita(j),dist)
         if ( dist*3600. .lt. max(poserr_erosita(j),2.) ) then
            xray_type=14
            write(12,'(4(es10.3,2x),2(f10.5,2x),f8.3,2(2x,f10.4),2x,i2)') frequency_erosita(j,1),
     &                 flux_erosita(j,1),FluxU_erosita(j,1),FluxL_erosita(j,1),ra_erosita(j),dec_erosita(j),
     &                   poserr_erosita(j),mjdavg,mjdavg,xray_type+50
            do s=2,8
               write(12,'(4(es10.3,2x),i2)') frequency_erosita(j,s),flux_erosita(j,s),FluxU_erosita(j,s),
     &                     FluxL_erosita(j,s),xray_type
            enddo
         endif
      enddo
502   close(11)
      close(12)
      close(13)
      close(14)

      deallocate(ra_xmm,dec_xmm,xmm_type,poserr_xmm)
      deallocate(ra_swift,dec_swift,xrt_type,poserr_swift,mjdst_swift,mjded_swift)
      deallocate(Ferr_swift,FluxU_swift,FluxL_swift,flux_swift,frequency_swift)
      deallocate(flux_xmm,Ferr_xmm,FluxU_xmm,FluxL_xmm,frequency_xmm)
      deallocate(ra_rosat,dec_rosat,poserr_rosat,rosat_type)
      deallocate(poserr_chandra,ra_chandra,dec_chandra)
      deallocate(flux_rosat,Ferr_rosat,FluxU_rosat,FluxL_rosat,frequency_rosat)
      deallocate(flux_chandra,FluxU_chandra,FluxL_chandra,frequency_chandra)
      deallocate(poserr_erosita,ra_erosita,dec_erosita)
      deallocate(poserr_bmw,ra_bmw,dec_bmw)
      deallocate(flux_erosita,fluxL_erosita,fluxU_erosita,Ferr_erosita,frequency_erosita)
      deallocate(flux_bmw,fluxL_bmw,fluxU_bmw,Ferr_bmw,frequency_bmw)
      deallocate(ipc_type,poserr_ipc,ra_ipc,dec_ipc)
      deallocate(poserr_maxi,ra_maxi,dec_maxi,maxi_type)
      deallocate(flux_ipc,Ferr_ipc,fluxL_ipc,fluxU_ipc,frequency_ipc)
      deallocate(flux_maxi,Ferr_maxi,FluxL_maxi,FluxU_maxi,frequency_maxi)
      deallocate(ra_radio,dec_radio,radio_type,ra_index)
      deallocate(ppss,const,Ferr_radio,FluxU_radio,FluxL_radio,poserr_radio,flux_radio,frequency_radio)

      END
c
      SUBROUTINE print_results (ratio,ra,dec,flux_radio,radio_type,xray_type,
     &                          flux,constsub,ra_center,dec_center,source_type)
      IMPLICIT none
      INTEGER*4 rah,ram,id,dm,radio_type,xray_type,source_type
      REAL*4 ratio,rasec,decsec,constsub,arx,lognupeak,flux_radio
      REAL*4 flux
      REAL*8 ra,dec,ra_center,dec_center,dist
      CHARACTER*1 sign
      CHARACTER*15 xmission,radio_survey
      CHARACTER*80 type
      ratio=flux/flux_radio
      call chra(ra,rah,ram,rasec,1)
      call chdec(dec,id,dm,decsec,1)
      sign(1:1) =' '
      if (dec < 0.0) sign(1:1)='-'
      type = ' (type unknown) '
c      IF (ratio < 0.) RETURN 
      IF (radio_type == 3) THEN
         radio_survey='NVSS'
      ELSE IF (radio_type == 2) THEN
         radio_survey='FIRST'
      ELSE IF (radio_type == 4) THEN
         radio_survey='SUMSS'
      Else IF (radio_type == 1) then
         radio_survey='VLASSQL'
      ELSE
         radio_survey='UNKNOWN'
      ENDIF
      IF (xray_type == 1) THEN
         xmission='XMMSLEW2'
      ELSE IF (xray_type == 2) THEN
         xmission='4XMM-DR11'
      ELSE IF (xray_type == 3) THEN
         xmission='RASS'
      ELSE IF (xray_type == 4) THEN
         xmission='WGA'
      ELSE IF (xray_type == 5) THEN
         xmission='2SXPS'
      ELSE IF (xray_type == 6) THEN
         xmission='IPC'
      ELSE IF (xray_type == 7) THEN
         xmission='BMW'
      ELSE IF (xray_type == 8) THEN
         xmission='Chandra-CSC2'
      ELSE IF (xray_type == 9) THEN
         xmission='XRTDEEP'
      ELSE IF (xray_type == 10) THEN
         xmission='MAXIGSC'
      ELSE IF (xray_type == 11) THEN
         xmission='1OUSX'
      ELSE IF (xray_type == 12) THEN
         xmission='IPCSLEW'
      ELSE IF (xray_type == 13) THEN
         xmission='MAXISSC'
      ELSE IF (xray_type == 14) THEN
         xmission='eROSITA-EDR'
      ELSE
         xmission='UNKNOWN'
      ENDIF
      source_type = 0
c      print *,' xmission, flux-x ',xmission,flux
      CALL DIST_SKY(ra,dec,ra_center,dec_center,dist)
      dist=dist*60.
         !write(*,*) 'test flux',flux
      IF (flux > 0.) THEN 
         arx = 1.-log10(ratio)/log10(2.41e17/1.4e9)
         !write(*,*) arx,'test arx'
         IF ( (abs(arx) > 0.42).AND.(abs(arx).LE.0.78) ) THEN
            lognupeak=(1.44-arx)/0.05
            IF (lognupeak > 15.5) THEN 
               type(1:1)=achar(27) 
               type(21:21)=achar(27) 
               type(2:20) = '[31;1m possible HBL' 
               type(22:24)  ='[0m'
               source_type = 1
            ELSE
               type(1:1)=achar(27) 
               type(21:21)=achar(27) 
               type(2:20) = '[36;1m possible IBL' 
               type(22:24)  ='[0m'
               source_type = 2
            ENDIF
            WRITE(*,'(1x,a,''/'',a,'' ra dec '',i2.2,1x,i2.2,1x,f4.1,
     &      '','',a,i2.2,1x,i2.2,1x,f4.1,'' radio flux d.'',
     &      f9.3,'' X-ray/radio flux-ratio'',f8.0,'' arx '',f6.3,
     &      '' Log(nu_p) '',f4.1,''+/-~1 '',a, 
     &      '' Dist. '',f7.3,'' arcmin'')')
     &      xmission(1:len_trim(xmission)),radio_survey(1:len_trim(radio_survey)),
     &      rah,ram,rasec,sign(1:1),abs(id),abs(dm),abs(decsec),
     &      flux_radio/constsub,ratio,arx,lognupeak,type(1:len_trim(type)),dist
         ELSE  
          IF ((abs(arx) > 0.78).AND.(abs(arx) < 0.95)) THEN
            type(1:30) = '           [34;1m possible LBL' 
            type(11:11)=achar(27) 
            type(32:32)=achar(27) 
            type(33:35)  ='[0m'
            source_type = 3
          ELSE IF (abs(arx) < 0.43) THEN
            type(1:1)=achar(27) 
            type(2:31) = '[32;1m possible non-jetted AGN'
            type(32:32)=achar(27) 
            type(33:35)  ='[0m'
            source_type = 4
          ELSE
            type(1:1)=achar(27) 
            type(21:21)=achar(27) 
            type(2:20) = '[38;1m type unknown'
            type(22:24)  ='[0m'
            source_type = 5 
          ENDIF
          WRITE(*,'(1x,a,''/'',a,'' ra dec '',i2.2,1x,i2.2,1x,f4.1,
     &       '','',a,i2.2,1x,i2.2,1x,f4.1,'' radio flux d. '',
     &       f9.3,'' flux-ratio'',f8.0,'' arx '',f6.3,13x,a,
     &       '' Dist. '',f7.3,'' arcmin'')')
     &       xmission(1:len_trim(xmission)),radio_survey(1:len_trim(radio_survey)),
     &       rah,ram,rasec,sign(1:1),abs(id),abs(dm),abs(decsec),
     &       flux_radio/constsub,ratio,arx,type(1:len_trim(type)),dist
         ENDIF
      ELSE
           source_type = 5
           type(1:1)=achar(27) 
           type(33:33)=achar(27) 
           type(2:32) = '[38;1m No X-ray flux available!'
           type(34:36) = '[0m'
         WRITE(*,'(1x,a,''/'',a,'' ra dec '',i2.2,1x,i2.2,1x,f4.1,
     &       1x,'','',a,i2.2,1x,i2.2,1x,f4.1,'' radio flux d. '',
     &       f9.3,a,'' Dist. '',f7.3,'' arcmin'')')
     &       xmission(1:len_trim(xmission)),radio_survey(1:len_trim(radio_survey)),
     &       rah,ram,rasec,sign(1:1),abs(id),abs(dm),abs(decsec),
     &       flux_radio/constsub,type(1:len_trim(type)),dist

      ENDIF
      RETURN
      END
C
      SUBROUTINE DIST_SKY(alpha1,delta1,alpha2,delta2,dist)
      IMPLICIT NONE
      REAL*8 dist,alpha1,alpha2,delta1,delta2,costheta
      REAL*8 radian
      radian=57.2957795 !(180/pi), degree to radian
      costheta=sin(delta1/radian)*sin(delta2/radian)+
     &         cos(delta1/radian)*cos(delta2/radian)*
     &         cos((alpha1-alpha2)/radian)
      dist=acos(costheta)*radian !radian to degree.
      RETURN
      END
c
ccc
      subroutine int_great_circle(lon1,lat1,lon2,lat2,f,d,lon,lat)
      implicit none
      real*8 lat1,lon1,lat2,lon2,lat,lon,d,radian
      real*8 A,B,x,y,z
      real*4 f

      radian=57.2957795
      lat1=(lat1)/radian
      lat2=lat2/radian
      lon1=lon1/radian
      lon2=lon2/radian
      d=d/radian
c      write(*,*) 'TEST value',lat1,lat2,lon1,lon2,d
      A=sin((1-f)*d)/sin(d)
      B=sin(f*d)/sin(d)
      x = A*cos(lat1)*cos(lon1) +  B*cos(lat2)*cos(lon2)
      y = A*cos(lat1)*sin(lon1) +  B*cos(lat2)*sin(lon2)
      z = A*sin(lat1)           +  B*sin(lat2)
c      write(*,*) 'TEST value',A,B,x,y,z
      lat=atan2(z,sqrt(x**2+y**2))*radian
      lon=atan2(y,x)*radian
      if (lon .lt. 0.) lon=lon+360.
      RETURN
      end
c

      SUBROUTINE graphic_code(flux_x,flux_r,source_type,code)
      IMPLICIT none 
      REAL*4 flux_x,flux_r,code,rfl_max,rfl_min
      REAL*4 xfl_min,xfl_max
      INTEGER*4 radio_component,x_ray_component,source_type,temp
      IF (source_type < 0) THEN
         code = source_type*10000
         RETURN
      ENDIF 
      temp = source_type
      rfl_min=0.8 ! 0.8 mJy
      rfl_max=8000. ! 8 Jy 
      xfl_min = 1.e-16 ! 1.e-16 erg/cm2/s, nufnu
      xfl_max = 5.e-11
      radio_component=int(alog10(flux_r/rfl_min)/alog10(rfl_max/rfl_min)*99.)
      IF (radio_component .GE. 99) THEN 
         radio_component = 99
      ELSE IF (radio_component .LE. 1) THEN 
         radio_component = 1
      ENDIF
      IF (flux_x > xfl_min) THEN  
         x_ray_component=int(alog10(flux_x/xfl_min)/alog10(xfl_max/xfl_min)*99.)
      ELSE
         x_ray_component = 1
      ENDIF
      IF (x_ray_component .GE. 99) THEN 
         x_ray_component = 99
      ELSE IF (x_ray_component .Lt. 1) THEN 
         x_ray_component = 1
      ENDIF
      code = source_type*10000.+radio_component+100.*x_ray_component
      !IF (code  <  10000 ) THEN
      !  print *,' source type radio_component, x_ray_component ', source_type,radio_component,x_ray_component
      !  stop
      !ENDIF
      source_type = temp
      RETURN  
      END
c PG
      SUBROUTINE RXgraphic_code(flux,RX,code)
      IMPLICIT none
      REAL*4 flux,code,rfl_max,rfl_min
      REAL*4 xfl_min,xfl_max
      INTEGER*4 radio_component,x_ray_component
      CHARACTER*1 RX 
      code = 0.
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
         code = -90000.
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
         code = -80000.
      ENDIF
      code = code -radio_component-100.*x_ray_component
      RETURN 
      END

c
      SUBROUTINE fluxtofdens(alpha,bandl,bandu,flux,kev,fdens,nudens)
      real*4 alpha,bandu,bandl,flux,nudens,fdens,conval,kev!,nuu,nul
      !write(*,*) alpha,flux,kev,bandu,bandl
      if (alpha .ne. 1. ) then
         conval=(1./(-alpha+1.))*((bandu)**(-alpha+1.)-(bandl)**(-alpha+1.))
      else
         conval=log(bandu/bandl)
      endif
      !write(*,*) kev,(1./conval)*kev*kev**(-alpha)!/(kev*2.418E-12)
      fdens=kev*(flux/conval)*((kev)**(-alpha))!!!!
      nudens=(1.602E-19)*(kev*1.e3)/(6.626e-34)
      RETURN
      end
