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
      INTEGER*4 ier, lu_in, ia,xray_type,lu_output, in,k, length,spec_type(5000,5000)
      INTEGER*4 radio_type(5000),xmm_type(5000),rosat_type(1000),rtype_source(5000),utype,iifound
      INTEGER*4 lenact,source_type,type_average,ix,types(0:5),xpts,spec_xpts(5000,5000)
      INTEGER*4 no_found,sfound,nrep(1000),rfound,s,track(1000),t(1000),aim,xrt_type(5000),ncat
      INTEGER*4 iradio,ixmm,irosat,iswift,iipc,iother,ichandra,ibmw,ifound,exits,iuv,isuv,xpts2,iuvx
      INTEGER*4 rah, ram, id, dm ,is,ie, i, j,ra_index(5000),l,filen,ttsource(5000),ihighpeak,track2(100)
      REAL*8 ra_other(10000),dec_other(10000),ra, dec,dist,ra_center, dec_center,radius,dec_1kev(5000,5000)
      REAL*8 ra_radio(5000),dec_radio(5000),ra_xmm(5000),dec_xmm(5000),ra_rosat(1000),dec_rosat(1000)
      REAL*8 ra_swift(5000),dec_swift(5000),ra_bmw(500),dec_bmw(500),ra_ipc(200),dec_ipc(200)
      REAL*8 ra_chandra(500),dec_chandra(500),ra_source(5000),dec_source(5000),ra_1kev(5000,5000)
      real*8 ra_cat(100),dec_cat(100)
      REAL*4 flux_radio(5000),flux_xmm(5000,6),flux_rosat(1000),flux_chandra(500,5),radian
      REAL*4 flux_swift(5000,4),flux_ipc(200),flux_bmw(500),flux_x,nh,ppss(5000)
      REAL*4 frequency_xmm(5000,6),frequency_bmw(500),frequency_rosat(1000)
      REAL*4 frequency_chandra(500,5),frequency_swift(5000,4),frequency_ipc(200)
      REAL*4 min_dist_rosat,min_dist_xmm,rasec,decsec,min_dist_ipc,min_dist_cluster
      REAL*4 min_dist_other,min_dist_swift,min_dist_bmw,min_dist_chandra
      REAL*4 flux2nufnu_nvss,flux2nufnu_rosat,flux2nufnu_xmm,min_dist,reduce
      REAL*4 flux2nufnu_swift,flux2nufnu_ipc,code,fdens,nudens,flux_source(5000),rrconst(5000)
      REAL*4 flux2nufnu_bmw,flux2nufnu_rxs,ratio,const(5000),flux2nufnu_sumss
      REAL*4 flux_1kev(5000,5000),frequency_radio(5000),uflux_1kev(5000,5000),lflux_1kev(5000,5000)
      REAL*4 xflux(5000),rflux(5000),flux_xpts(5000,5000),frequency_xpts(500,5000),poserr_1kev(5000,5000)
      real*4 major,minor,posang,posxerr,posyerr,uflux_xpts(5000,5000),lflux_xpts(5000,5000),distrx(5000,5000)
      real*4 Ferr_radio(5000),FluxU_radio(5000),FluxL_radio(5000),poserr_radio(5000),distrx2(5000,5000)
      real*4 Ferr_xmm(5000,6),FluxU_xmm(5000,6),FluxL_xmm(5000,6),poserr_xmm(5000)
      real*4 Ferr_rosat(1000),FluxU_rosat(1000),FluxL_rosat(1000),poserr_rosat(1000)
      real*4 Ferr_swift(5000,4),FluxU_swift(5000,4),FluxL_swift(5000,4),poserr_swift(5000)
      real*4 Ferr_ipc(200),FluxU_ipc(200),FluxL_ipc(200),poserr_ipc(200)
      real*4 Ferr_bmw(500),FluxU_bmw(500),FluxL_bmw(500),poserr_bmw(500)
      real*4 Ferr_chandra(500,5),FluxU_chandra(500,5),FluxL_chandra(500,5),poserr_chandra(500)
      real*4 errrad,errmaj,errmin,errang,savemjy(10000)
      CHARACTER*1 sign
      CHARACTER*30 name_other(10000),name_cat(10000)
      CHARACTER*80 input_file,output_file,output_file2
      CHARACTER*8 catalog,uv_type(20000)
      CHARACTER*800 string,repflux
      LOGICAL there,ok,found 
      ok = .TRUE. 
      found = .FALSE.
      nrep(1:500)=1
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
      radian = 45.0/atan(1.0)
      flux2nufnu_nvss=1.4e9*1.e-26
      flux2nufnu_sumss=8.43e8*1.e-26 !assumed radio alpha=0.2 !f_0.8 to f_1.4
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

      CALL rdforn(string,length)
      IF ( length.NE.0 ) THEN
         CALL rmvlbk(string)
         in=index(string(1:length),' ')
         input_file=string(1:in-1)
         read(string(in+1:length),*) aim
         !write(*,*) 'the aim',aim
      ELSE 
         WRITE (*,'('' Enter query results file '',$)')
         READ (*,'(a)') input_file
      ENDIF
      output_file='find_out_temp.txt'
      output_file2='RX_temp.txt'
      lu_in = 10
      !lu_output=11
      in = index(input_file(1:lenact(input_file)),'.')
      IF (in == 0) input_file(lenact(input_file)+1:lenact(input_file)+4) = '.csv' 
      INQUIRE (FILE=input_file,EXIST=there)
      IF (.NOT.there) THEN
         write (*,'('' file '',a,'' not found '')')
     &     input_file(1:lenact(input_file))
         STOP
      ENDIF

      open(lu_in,file=input_file,status='old',iostat=ier)
      open(11,file=output_file,status='unknown',iostat=ier)
      open(13,file=output_file2,status='unknown',iostat=ier)
      open(14,file='no_matched_temp.txt',status='unknown',iostat=ier)
      !write(*,*) ier
      open(12,file='Sed_temp.txt',status='unknown')
      IF (ier.NE.0) THEN
        write (*,*) ' Error ',ier,' opening file ', input_file
      ENDIF

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
         IF ( (catalog(1:4) == 'nvss') .OR.
     &        (catalog(1:5) == 'first') .OR.
     &        (catalog(1:5) == 'sumss') ) THEN
            iradio=iradio+1
            IF (iradio > 5000) Stop 'Too many NVSS/SUMSS points'
            ra_radio(iradio)=ra
            dec_radio(iradio)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            IF (catalog(1:4) == 'nvss') radio_type(iradio) = 2
            IF (catalog(1:5) == 'first') radio_type(iradio) = 1
            IF (catalog(1:5) == 'sumss') radio_type(iradio) = 3
            IF ((catalog(1:5) == 'sumss') .or. (catalog(1:4) == 'nvss')) THEN
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_radio(iradio)
               poserr_radio(iradio)=poserr_radio(iradio)*2. !95% error
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_radio(iradio)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_radio(iradio)
            ELSE
               if (is .ne. ie-1) read(string(is+1:ie-1),*) ppss(iradio)
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
               poserr_radio(iradio)=major*((Ferr_radio(iradio)/(flux_radio(iradio)-0.25))+0.05 )
               if (poserr_radio(iradio) .lt. 0.1 ) poserr_radio(iradio)=0.1
               !write(*,*) flux_radio(iradio)
            ENDIF
            IF (catalog(1:5) == 'sumss') then
               const(iradio) = flux2nufnu_sumss
               frequency_radio(iradio)=8.43E8
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
            if ((ppss(iradio) .gt. 0.15) .and. (radio_type(iradio) .eq. 1)) iradio=iradio-1
            !write(*,*) catalog,flux_radio(iradio)/const(iradio),iradio
         ELSE IF ( (catalog(1:5) == 'xmmsl') .OR.
     &             (catalog(1:4) == '3xmm') )  THEN
            ixmm=ixmm+1
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            ra_xmm(ixmm)=ra
            dec_xmm(ixmm)=dec
            IF (ixmm > 5000) Stop 'Too many XMM points'
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
            !write(*,*) FluxU_xmm(ixmm,1),flux_xmm(ixmm,1),FluxL_xmm(ixmm,1)
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
            IF (catalog(1:5) == 'xmmsl') THEN
               xmm_type(ixmm)=1
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xmm(ixmm,2) !6
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xmm(ixmm,2)
               FluxU_xmm(ixmm,2)=flux_xmm(ixmm,2)+Ferr_xmm(ixmm,2)
               FluxL_xmm(ixmm,2)=flux_xmm(ixmm,2)-Ferr_xmm(ixmm,2)
               if (Ferr_xmm(ixmm,2) .gt. flux_xmm(ixmm,2)) then
                  FluxU_xmm(ixmm,2)=0.!Ferr_xmm(ixmm,2)*3.
                  FluxL_xmm(ixmm,2)=0.
                  Flux_xmm(ixmm,2)=0.
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
               if (flux_xmm(ixmm,1) .eq. 0.) then
                  if (flux_xmm(ixmm,2) .ne. 0) then
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,2)
                     Ferr_xmm(ixmm,1)=Ferr_xmm(ixmm,2)
                     FluxU_xmm(ixmm,1)=flux_xmm(ixmm,1)+Ferr_xmm(ixmm,1)
                     FluxL_xmm(ixmm,1)=flux_xmm(ixmm,1)-Ferr_xmm(ixmm,1)
                     call nhdeabsorb2 (1,0.2,2.,0.9,nh,reduce,100)
                     flux_xmm(ixmm,1)=flux_xmm(ixmm,1)*reduce
                     FluxU_xmm(ixmm,1)=FluxU_xmm(ixmm,1)*reduce
                     FluxL_xmm(ixmm,1)=FluxL_xmm(ixmm,1)*reduce
                     call fluxtofdens(0.9,0.2,2.,flux_xmm(ixmm,1),1.0,fdens,nudens)
                     flux_xmm(ixmm,1)=fdens
                     frequency_xmm(ixmm,1)=nudens
                     call fluxtofdens(0.9,0.2,2.,FluxU_xmm(ixmm,1),1.0,fdens,nudens)
                     FluxU_xmm(ixmm,1)=fdens
                     call fluxtofdens(0.9,0.2,2.,FluxL_xmm(ixmm,1),1.0,fdens,nudens)
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
               call nhdeabsorb2 (1,0.2,2.,0.9,nh,reduce,100)
               !write(*,*) reduce
               flux_xmm(ixmm,2)=flux_xmm(ixmm,2)*reduce
               FluxU_xmm(ixmm,2)=FluxU_xmm(ixmm,2)*reduce
               FluxL_xmm(ixmm,2)=FluxL_xmm(ixmm,2)*reduce
               call fluxtofdens(0.9,0.2,2.,flux_xmm(ixmm,2),1.1,fdens,nudens)
               flux_xmm(ixmm,2)=fdens
               frequency_xmm(ixmm,2)=nudens
               call fluxtofdens(0.9,0.2,2.,FluxU_xmm(ixmm,2),1.1,fdens,nudens)
               FluxU_xmm(ixmm,2)=fdens
               call fluxtofdens(0.9,0.2,2.,FluxL_xmm(ixmm,2),1.1,fdens,nudens)
               FluxL_xmm(ixmm,2)=fdens
               call nhdeabsorb2 (1,2.,12.,0.9,nh,reduce,100)
               flux_xmm(ixmm,3)=flux_xmm(ixmm,3)*reduce
               FluxU_xmm(ixmm,3)=FluxU_xmm(ixmm,3)*reduce
               FluxL_xmm(ixmm,3)=FluxL_xmm(ixmm,3)*reduce
               call fluxtofdens(0.9,2.,12.,flux_xmm(ixmm,3),7.,fdens,nudens)
               flux_xmm(ixmm,3)=fdens
               frequency_xmm(ixmm,3)=nudens
               call fluxtofdens(0.9,2.,12.,FluxU_xmm(ixmm,3),7.,fdens,nudens)
               FluxU_xmm(ixmm,3)=fdens
               call fluxtofdens(0.9,2.,12.,FluxL_xmm(ixmm,3),7.,fdens,nudens)
               FluxL_xmm(ixmm,3)=fdens
            ELSE IF (catalog(1:4) == '3xmm') THEN
               xmm_type(ixmm)=2
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
               call fluxtofdens(0.9,0.2,0.5,flux_xmm(ixmm,2),0.35,fdens,nudens)
               flux_xmm(ixmm,2)=fdens
               frequency_xmm(ixmm,2)=nudens
               call fluxtofdens(0.9,0.2,0.5,FluxU_xmm(ixmm,2),0.35,fdens,nudens)
               FluxU_xmm(ixmm,2)=fdens
               call fluxtofdens(0.9,0.2,0.5,FluxL_xmm(ixmm,2),0.35,fdens,nudens)
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
               call fluxtofdens(0.9,2.,4.5,flux_xmm(ixmm,5),3.25,fdens,nudens)
               flux_xmm(ixmm,5)=fdens
               frequency_xmm(ixmm,5)=nudens
               call fluxtofdens(0.9,2.,4.5,FluxU_xmm(ixmm,5),3.25,fdens,nudens)
               FluxU_xmm(ixmm,5)=fdens
               call fluxtofdens(0.9,2.,4.5,FluxL_xmm(ixmm,5),3.25,fdens,nudens)
               FluxL_xmm(ixmm,5)=fdens
               call nhdeabsorb2 (1,4.5,12.,0.9,nh,reduce,100)
               !write(*,*) reduce
               flux_xmm(ixmm,6)=flux_xmm(ixmm,6)*reduce
               FluxU_xmm(ixmm,6)=FluxU_xmm(ixmm,6)*reduce
               FluxL_xmm(ixmm,6)=FluxL_xmm(ixmm,6)*reduce
               call fluxtofdens(0.9,4.5,12.,flux_xmm(ixmm,6),8.25,fdens,nudens)
               flux_xmm(ixmm,6)=fdens
               call fluxtofdens(0.9,4.5,12.,FluxU_xmm(ixmm,6),8.25,fdens,nudens)
               FluxU_xmm(ixmm,6)=fdens
               call fluxtofdens(0.9,4.5,12.,FluxL_xmm(ixmm,6),8.25,fdens,nudens)
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
            IF (irosat > 1000) Stop 'Too many RASS points'
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
         ELSE IF ((catalog(1:4) == 'sxps') .or. (catalog(1:7) == 'xrtdeep'))THEN
            iswift=iswift+1
            IF (iswift > 5000) Stop 'Too many swift points'
            ra_swift(iswift)=ra
            dec_swift(iswift)=dec
            if (catalog(1:4) == 'sxps') then
               xrt_type(iswift)=1
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
               call fluxtofdens(0.9,0.3,10.,flux_swift(iswift,1),1.,fdens,nudens)
               flux_swift(iswift,1)=fdens
               frequency_swift(iswift,1)=nudens
               call fluxtofdens(0.9,0.3,10.,FluxU_swift(iswift,1),1.,fdens,nudens)
               FluxU_swift(iswift,1)=fdens
               call fluxtofdens(0.9,0.3,10.,FluxL_swift(iswift,1),1.,fdens,nudens)
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
                     call fluxtofdens(0.9,1.,2.,flux_swift(iswift,1),1.,fdens,nudens)
                     flux_swift(iswift,1)=fdens
                     frequency_swift(iswift,1)=nudens
                     call fluxtofdens(0.9,1.,2.,FluxU_swift(iswift,1),1.,fdens,nudens)
                     FluxU_swift(iswift,1)=fdens
                     call fluxtofdens(0.9,1.,2.,FluxL_swift(iswift,1),1.,fdens,nudens)
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
                     call fluxtofdens(0.9,0.3,1.,flux_swift(iswift,1),1.,fdens,nudens)
                     flux_swift(iswift,1)=fdens
                     frequency_swift(iswift,1)=nudens
                     call fluxtofdens(0.9,0.3,1.,FluxU_swift(iswift,1),1.,fdens,nudens)
                     FluxU_swift(iswift,1)=fdens
                     call fluxtofdens(0.9,0.3,1.,FluxL_swift(iswift,1),1.,fdens,nudens)
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
                     call fluxtofdens(0.9,2.,10.,flux_swift(iswift,1),1.,fdens,nudens)
                     flux_swift(iswift,1)=fdens
                     frequency_swift(iswift,1)=nudens
                     call fluxtofdens(0.9,2.,10.,FluxU_swift(iswift,1),1.,fdens,nudens)
                     FluxU_swift(iswift,1)=fdens
                     call fluxtofdens(0.9,2.,10.,FluxL_swift(iswift,1),1.,fdens,nudens)
                     FluxL_swift(iswift,1)=fdens
                  endif
               endif
               call nhdeabsorb2 (0,0.3,1.,0.9,nh,reduce,4)
               flux_swift(iswift,2)=flux_swift(iswift,2)*reduce
               FluxU_swift(iswift,2)=FluxU_swift(iswift,2)*reduce
               FluxL_swift(iswift,2)=FluxL_swift(iswift,2)*reduce
               call fluxtofdens(0.9,0.3,1.,flux_swift(iswift,2),0.65,fdens,nudens)
               flux_swift(iswift,2)=fdens
               frequency_swift(iswift,2)=nudens
               call fluxtofdens(0.9,0.3,1.,FluxU_swift(iswift,2),0.65,fdens,nudens)
               FluxU_swift(iswift,2)=fdens
               call fluxtofdens(0.9,0.3,1.,FluxL_swift(iswift,2),0.65,fdens,nudens)
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
               call fluxtofdens(0.9,2.,10.,flux_swift(iswift,4),6.,fdens,nudens)
               flux_swift(iswift,4)=fdens
               frequency_swift(iswift,4)=nudens
               call fluxtofdens(0.9,2.,10.,FluxU_swift(iswift,4),6.,fdens,nudens)
               FluxU_swift(iswift,4)=fdens
               call fluxtofdens(0.9,2.,10.,FluxL_swift(iswift,4),6.,fdens,nudens)
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
            else if (catalog(1:7) == 'xrtdeep') then
               xrt_type(iswift)=2
               is=ie
               ie=index(string(is+1:len(string)),',')+is !nh
               is=ie
               ie=index(string(is+1:len(string)),',')+is !slope
               is=ie
               ie=index(string(is+1:len(string)),',')+is !slope err
               is=ie
               ie=index(string(is+1:len(string)),',')+is !exp
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
               poserr_swift(iswift)=7.!!!!!!!!
            endif
c PG
            CALL RXgraphic_code(flux_swift(iswift,1),'X',code)
            write (13,'(f9.5,2x,f9.5,2x,i6)') ra_swift(iswift),dec_swift(iswift),int(code)
c end PG
         ELSE IF (catalog(1:3) == 'ipc') THEN
            iipc=iipc+1
            IF (iipc > 200) Stop 'Too many Einstein IPC points'
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
            else
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_ipc(iipc)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_ipc(iipc)
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
            IF (ibmw > 500) Stop 'Too many ROSAT-BMW points'
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
         ELSE IF (catalog(1:7) == 'chandra') THEN
            ichandra=ichandra+1
            IF (ichandra > 500) Stop 'Too many Chandra points'
            ra_chandra(ichandra)=ra
            dec_chandra(ichandra)=dec
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_chandra(ichandra)
            poserr_chandra(ichandra)=sqrt((poserr_chandra(ichandra)**2)+(0.8**2))
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_chandra(ichandra,1)
            call nhdeabsorb2 (1,0.5,7.,0.9,nh,reduce,100)
            flux_chandra(ichandra,1)=flux_chandra(ichandra,1)*reduce
            call fluxtofdens(0.9,0.5,7.,flux_chandra(ichandra,1),1.,fdens,nudens)
            flux_chandra(ichandra,1)=fdens
            frequency_chandra(ichandra,1)=nudens
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_chandra(ichandra,2) !h
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_chandra(ichandra,3) !m
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_chandra(ichandra,4) !s
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_chandra(ichandra,5) !us
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxU_chandra(ichandra,1)
            if ((flux_chandra(ichandra,1) .eq. 0.) .and. (FluxU_chandra(ichandra,1) .ne. 0.))
     &             FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,1)*3.
            call nhdeabsorb2 (1,0.5,7.,0.9,nh,reduce,100)
            FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,1)*reduce
            call fluxtofdens(0.9,0.5,7.,FluxU_chandra(ichandra,1),1.,fdens,nudens)
            FluxU_chandra(ichandra,1)=fdens
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
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxU_chandra(ichandra,4) !s
            if ((flux_chandra(ichandra,4) .eq. 0.) .and. (FluxU_chandra(ichandra,4) .ne. 0.))
     &             FluxU_chandra(ichandra,4)=FluxU_chandra(ichandra,4)*3.
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxU_chandra(ichandra,5) !us
            if ((flux_chandra(ichandra,5) .eq. 0.) .and. (FluxU_chandra(ichandra,5) .ne. 0.))
     &             FluxU_chandra(ichandra,5)=FluxU_chandra(ichandra,5)*3.
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxL_chandra(ichandra,1)
            call nhdeabsorb2 (1,0.5,7.,0.9,nh,reduce,100)
            FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,1)*reduce
            call fluxtofdens(0.9,0.5,7.,FluxL_chandra(ichandra,1),1.,fdens,nudens)
            FluxL_chandra(ichandra,1)=fdens
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxL_chandra(ichandra,2) !h
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxL_chandra(ichandra,3) !m
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*)FluxL_chandra(ichandra,4) !s
            is=ie
            ie=index(string(is+1:len(string)),' ')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) FluxL_chandra(ichandra,5) !us
            if (flux_chandra(ichandra,1) .eq. 0.) then
               if (flux_chandra(ichandra,4) .ne. 0.) then
                  flux_chandra(ichandra,1)=flux_chandra(ichandra,4)
                  FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,4)
                  FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,4)
                  call nhdeabsorb2 (1,0.5,1.2,0.9,nh,reduce,100)
                  flux_chandra(ichandra,1)=flux_chandra(ichandra,1)*reduce
                  FluxU_chandra(ichandra,1)=FluxU_chandra(ichandra,1)*reduce
                  FluxL_chandra(ichandra,1)=FluxL_chandra(ichandra,1)*reduce
                  call fluxtofdens(0.9,0.5,1.2,flux_chandra(ichandra,1),1.,fdens,nudens)
                  flux_chandra(ichandra,1)=fdens
                  frequency_chandra(ichandra,1)=nudens
                  call fluxtofdens(0.9,0.5,1.2,fluxU_chandra(ichandra,1),1.,fdens,nudens)
                  FluxU_chandra(ichandra,1)=fdens
                  call fluxtofdens(0.9,0.5,1.2,fluxL_chandra(ichandra,1),1.,fdens,nudens)
                  FluxL_chandra(ichandra,1)=fdens
               ELSE if (flux_xmm(ixmm,3) .ne. 0) then
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
               ELSE if (flux_xmm(ixmm,5) .ne. 0) then
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
            call nhdeabsorb2 (1,1.2,2.,0.9,nh,reduce,100)
            flux_chandra(ichandra,3)=flux_chandra(ichandra,3)*reduce
            call fluxtofdens(0.9,1.2,2.,flux_chandra(ichandra,3),1.6,fdens,nudens)
            flux_chandra(ichandra,3)=fdens
            frequency_chandra(ichandra,3)=nudens
            call nhdeabsorb2 (1,0.5,1.2,0.9,nh,reduce,100)
            flux_chandra(ichandra,4)=flux_chandra(ichandra,4)*reduce
            call fluxtofdens(0.9,0.5,1.2,flux_chandra(ichandra,4),0.85,fdens,nudens)
            flux_chandra(ichandra,4)=fdens
            frequency_chandra(ichandra,4)=nudens
            call nhdeabsorb2 (1,0.2,0.5,0.9,nh,reduce,100)
            flux_chandra(ichandra,5)=flux_chandra(ichandra,5)*reduce
            call fluxtofdens(0.9,0.2,0.5,flux_chandra(ichandra,5),0.35,fdens,nudens)
            flux_chandra(ichandra,5)=fdens
            frequency_chandra(ichandra,5)=nudens
            call nhdeabsorb2 (1,2.,7.,0.9,nh,reduce,100)
            FluxU_chandra(ichandra,2)=FluxU_chandra(ichandra,2)*reduce
            call fluxtofdens(0.9,2.,7.,fluxU_chandra(ichandra,2),4.5,fdens,nudens)
            FluxU_chandra(ichandra,2)=fdens
            call nhdeabsorb2 (1,1.2,2.,0.9,nh,reduce,100)
            FluxU_chandra(ichandra,3)=FluxU_chandra(ichandra,3)*reduce
            call fluxtofdens(0.9,1.2,2.,fluxU_chandra(ichandra,3),1.6,fdens,nudens)
            FluxU_chandra(ichandra,3)=fdens
            call nhdeabsorb2 (1,0.5,1.2,0.9,nh,reduce,100)
            FluxU_chandra(ichandra,4)=FluxU_chandra(ichandra,4)*reduce
            call fluxtofdens(0.9,0.5,1.2,fluxU_chandra(ichandra,4),0.85,fdens,nudens)
            FluxU_chandra(ichandra,4)=fdens
            call nhdeabsorb2 (1,0.2,0.5,0.9,nh,reduce,100)
            FluxU_chandra(ichandra,5)=FluxU_chandra(ichandra,5)*reduce
            call fluxtofdens(0.9,0.2,0.5,fluxU_chandra(ichandra,5),0.35,fdens,nudens)
            FluxU_chandra(ichandra,5)=fdens
            call nhdeabsorb2 (1,2.,7.,0.9,nh,reduce,100)
            FluxL_chandra(ichandra,2)=FluxL_chandra(ichandra,2)*reduce
            call fluxtofdens(0.9,2.,7.,fluxL_chandra(ichandra,2),4.5,fdens,nudens)
            FluxL_chandra(ichandra,2)=fdens
            call nhdeabsorb2 (1,1.2,2.,0.9,nh,reduce,100)
            FluxL_chandra(ichandra,3)=FluxL_chandra(ichandra,3)*reduce
            call fluxtofdens(0.9,1.2,2.,fluxL_chandra(ichandra,3),1.6,fdens,nudens)
            FluxL_chandra(ichandra,3)=fdens
            call nhdeabsorb2 (1,0.5,1.2,0.9,nh,reduce,100)
            FluxL_chandra(ichandra,4)=FluxL_chandra(ichandra,4)*reduce
            call fluxtofdens(0.9,0.5,1.2,fluxL_chandra(ichandra,4),0.85,fdens,nudens)
            FluxL_chandra(ichandra,4)=fdens
            call nhdeabsorb2 (1,0.2,0.5,0.9,nh,reduce,100)
            FluxL_chandra(ichandra,5)=FluxL_chandra(ichandra,5)*reduce
            call fluxtofdens(0.9,0.2,0.5,fluxL_chandra(ichandra,5),0.35,fdens,nudens)
            FluxL_chandra(ichandra,5)=fdens
            !write(*,*) fluxU_chandra(ichandra,1),flux_chandra(ichandra,1),fluxL_chandra(ichandra,1)
            !write(*,*) fluxU_chandra(ichandra,5),flux_chandra(ichandra,5),fluxL_chandra(ichandra,5)
c PG
            CALL RXgraphic_code(flux_chandra(ichandra,1),'X',code)
            write (13,'(f9.5,2x,f9.5,2x,i6)') ra_chandra(ichandra),dec_chandra(ichandra),int(code)
c end PG
         ELSE
            iother=iother+1
            IF (iother > 10000) Stop 'Too many catalogued sources'
            ra_other(iother)=ra
            dec_other(iother)=dec
            is=ie
            ie=len(string)
            read(string(is+1:ie),'(a)') name_other(iother)
            if ((catalog(1:4) == 'mcxc') .or. (catalog(1:2) == 'zw') .or.
     &           (catalog(1:3) == 'whl')) then
               name_other(iother)(5:lenact(name_other(iother))+5)=name_other(iother)(1:lenact(name_other(iother)))
               name_other(iother)(1:4)=catalog(1:4)
            else if (catalog(1:6) == 'pulsar') then
               name_other(iother)(5:lenact(name_other(iother))+5)=name_other(iother)(1:lenact(name_other(iother)))
               name_other(iother)(1:4)='PSR '
            endif
            !write(*,*) name_other(iother),iother
            if (iother .ne. 1) then
               do i=1,iother-1
                  if (name_other(i) == name_other(i+1)) iother=iother-1
               enddo
            endif
            !write(*,*) iother
         ENDIF
      ENDDO
 99   CONTINUE
      CLOSE (lu_in)
      if (aim .eq. 0) goto 501

      CALL indexx (iradio,ra_radio,ra_index)
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
            if (xmm_type(i) == 2) min_dist_xmm=4./3600.
            if (xmm_type(i) == 1) min_dist_xmm=15./3600.
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_xmm(i),dec_xmm(i),dist)
            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
               min_dist = sqrt(min_dist_xmm**2+(5./3600.)**2)
            ELSE
               min_dist = min_dist_xmm
            ENDIF
            IF (dist < min_dist) THEN 
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
                  spec_type(ix,k)=xray_type+10
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
               CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),xray_type,
     &                             flux_xmm(i,1),const(k),ra_center,dec_center,source_type)
               IF (source_type .GE. 0) THEN 
                  types(source_type) = types(source_type) + 1
               ENDIF
            ENDIF
         ENDDO
         DO i=1,irosat
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_rosat(i),dec_rosat(i),dist)
            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
               min_dist = sqrt(min_dist_rosat**2+(5./3600.)**2)
            ELSE
               min_dist = min_dist_rosat
            ENDIF
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
               spec_type(ix,k)=xray_type+10
               CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_rosat(i),const(k),ra_center,dec_center,source_type)
               IF (source_type .GE. 0) THEN 
                  types(source_type) = types(source_type) + 1
               ENDIF
            ENDIF
         ENDDO
         DO i=1,iswift
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_swift(i),dec_swift(i),dist)
            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
               min_dist = sqrt(min_dist_swift**2+(5./3600.)**2)
            ELSE
               min_dist = min_dist_swift
            ENDIF
            IF (dist < min_dist) THEN
               IF (xrt_type(i) == 1) THEN
                  xray_type = 5
               ELSE IF (xrt_type(i) == 2) THEN
                  xray_type = 9
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
               spec_type(ix,k)=xray_type+10
               do l=2,4
                  xpts=xpts+1
                  flux_xpts(xpts,k)=flux_swift(i,l)
                  uflux_xpts(xpts,k)=FluxU_swift(i,l)
                  lflux_xpts(xpts,k)=FluxL_swift(i,l)
                  frequency_xpts(xpts,k)=frequency_swift(i,l)
                  spec_xpts(xpts,k)=xray_type
               enddo
               CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_swift(i,1),const(k),ra_center,dec_center,source_type)
               IF (source_type .GE. 0) THEN
                  types(source_type) = types(source_type) + 1
               ENDIF
            ENDIF
         ENDDO
         DO i=1,iipc
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_ipc(i),dec_ipc(i),dist)
            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
               min_dist = sqrt(min_dist_ipc**2+(5./3600.)**2)
            ELSE
               min_dist = min_dist_ipc
            ENDIF
            IF (dist < min_dist) THEN 
               found = .TRUE.
               xray_type = 6
               flux_x = flux_x + flux_ipc(i) 
               ix = ix +1
               flux_1kev(ix,k)=flux_ipc(i)
               uflux_1kev(ix,k)=FluxU_ipc(i)
               lflux_1kev(ix,k)=FluxL_ipc(i)
               ra_1kev(ix,k)=ra_ipc(i)
               dec_1kev(ix,k)=dec_ipc(i)
               poserr_1kev(ix,k)=poserr_ipc(i)
               distrx(ix,k)=dist*3600.
               spec_type(ix,k)=xray_type+10
               CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_ipc(i),const(k),ra_center,dec_center,source_type)
               IF (source_type .GE. 0) THEN 
                  types(source_type) = types(source_type) + 1
               ENDIF
            ENDIF
         ENDDO
         DO i=1,ibmw
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_bmw(i),dec_bmw(i),dist)
            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
               min_dist = sqrt(min_dist_bmw**2+(5./3600.)**2)
            ELSE
               min_dist = min_dist_bmw
            ENDIF
            IF (dist < min_dist) THEN
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
               spec_type(ix,k)=xray_type+10
               CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_bmw(i),const(k),ra_center,dec_center,source_type)
               IF (source_type .GE. 0) THEN 
                  types(source_type) = types(source_type) + 1
               ENDIF
            ENDIF
         ENDDO
         DO i=1,ichandra
            CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_chandra(i),dec_chandra(i),dist)
            IF (radio_type(k) == 3) THEN ! 5 arcsec increase of min_dist for the case of SUMSS
            min_dist = sqrt(min_dist_chandra**2+(5./3600.)**2)
            ELSE
            min_dist = min_dist_chandra
            ENDIF
            IF (dist < min_dist) THEN
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
               spec_type(ix,k)=xray_type+10
               CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                xray_type,flux_chandra(i,1),const(k),ra_center,dec_center,source_type)
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
         !write(*,*) const
         IF (found) THEN 
            ifound = ifound +1
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
                  if (dist*3600 .lt. 6.) then
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
                     write(12,'(i4,2x,a,2(2x,f9.5),2x,a,2x,i2)') track(ifound),"matched source",
     &                  ra_radio(k),dec_radio(k),"source type",type_average
                     write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_radio(k),flux_radio(k),
     &                FluxU_radio(k),FluxL_radio(k),ra_radio(k),dec_radio(k),poserr_radio(k),radio_type(k)
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
            write(12,'(i4,2x,a,2(2x,f9.5),2x,a,2x,i2)') sfound,"matched source",
     &         ra_radio(k),dec_radio(k),'source type',type_average
            write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_radio(k),flux_radio(k),FluxU_radio(k),
     &          FluxL_radio(k),ra_radio(k),dec_radio(k),poserr_radio(k),radio_type(k)
            do i=1,ix
               write(12,'(" 2.418E+17",3(2x,es10.3),2(2x,f9.5),2x,f7.3,2x,i2)') flux_1kev(i,k),uflux_1kev(i,k),
     &             lflux_1kev(i,k),ra_1kev(i,k),dec_1kev(i,k),poserr_1kev(i,k),spec_type(i,k)
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
            rrconst(sfound)=const(k)
            ttsource(sfound)=type_average
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
                  rrconst(sfound)=const(k)
               else if ((radio_type(k) .eq. rtype_source(sfound)) .and.
     &                      (flux_radio(k) .gt. flux_source(sfound)))then
                  !write(*,*) "replace flux"
                  ra_source(sfound)=ra_radio(k)
                  dec_source(sfound)=dec_radio(k)
                  flux_source(sfound)=flux_radio(k)
                  rrconst(sfound)=const(k)
               endif
               goto 98
               !write(*,*) "Final radio",ra_source(sfound),dec_source(sfound),flux_source(sfound)
            endif
            write(*,*) '................Cataloged sources.................'
            DO i=1,iother
               CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_other(i),
     &                       dec_other(i),dist)
               IF (dist < min_dist_other) THEN 
                  write(*,'(2x,a,1x,a)') name_other(i)
                  IF (name_other(i)(1:4) == '3HSP') THEN
                    type_average = -1
                  ELSE IF (name_other(i)(1:3) == '5BZ') THEN
                    type_average = -2
                  ELSE IF (name_other(i)(1:6) == 'CRATES') THEN
                    type_average = -3
                  ELSE IF (name_other(i)(1:3) == 'PSR') THEN
                    type_average = 0
                    code=-8888
                  ENDIF
                  ra_other(i) = -ra_other(i)
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
                  ENDIF
               ENDIF
               IF (type_average < 0) THEN
                  CALL graphic_code (flux_x,flux_radio(k),type_average,code)
                  write(11,'(f9.5,2x,f9.5,2x,i6)') ra_radio(k),dec_radio(k),int(code)
               ELSE if (code .eq. -8888) then
                  write(11,'(f9.5,2x,f9.5,2x,i6)') ra_radio(k),dec_radio(k),int(code)
               ENDIF
               type_average=0
            ENDDO
            t(ifound)=k !!!recourd the former index
            write(*,*) '        '
         ELSE !!check radio without matched
            do i=1,iother
               if ( ( (name_other(i)(1:3) == '5BZ') .OR. (name_other(i)(1:4) == '3HSP') .or.
     &             (name_other(i)(1:6) == 'CRATES') .or. (name_other(i)(1:3) == 'PSR')) .AND.
     &                  (ra_other(i) .gt. 0.) ) THEN
                  CALL DIST_SKY(ra_other(i),dec_other(i),ra_radio(k),dec_radio(k),dist)
                  if (dist .le. min_dist_other ) found=.true.
               endif
            enddo
            if (.not. found) then
               CALL DIST_SKY(ra_center,dec_center,ra_radio(k),dec_radio(k),dist)
               if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
                  if (dist .le. errrad/60.) then
                     write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_radio(k),flux_radio(k),
     &               FluxU_radio(k),FluxL_radio(k),ra_radio(k),dec_radio(k),poserr_radio(k),radio_type(k)
                  endif
               else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
                  if (dist .le. errmaj/60.) then
                     write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_radio(k),flux_radio(k),
     &               FluxU_radio(k),FluxL_radio(k),ra_radio(k),dec_radio(k),poserr_radio(k),radio_type(k)
                  endif
               endif
            endif
         ENDIF
  98     continue
      ENDDO

      Do i=1,ixmm
         found=.false.
         if (xmm_type(i) == 2) min_dist_xmm=4./3600.
         if (xmm_type(i) == 1) min_dist_xmm=15./3600.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),ra_xmm(i),dec_xmm(i),dist)
            if (dist .le. min_dist_xmm)  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .OR. (name_other(i)(1:3) == 'PSR'))
     &                  .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_xmm(i),dec_xmm(i),dist)
               if (dist .le. min_dist_other) found=.true.
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
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_xmm(i,1),flux_xmm(i,1),
     &            FluxU_xmm(i,1),FluxL_xmm(i,1),ra_xmm(i),dec_xmm(i),poserr_xmm(i),xray_type+10
                  if (xray_type .eq. 1) then
                     do s=2,3
                        write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_xmm(i,s),flux_xmm(i,s),
     $                   FluxU_xmm(i,s),FluxL_xmm(i,s),ra_xmm(i),dec_xmm(i),poserr_xmm(i),xray_type+50
                     enddo
                  else
                     do s=2,6
                        write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_xmm(i,s),flux_xmm(i,s),
     $                   FluxU_xmm(i,s),FluxL_xmm(i,s),ra_xmm(i),dec_xmm(i),poserr_xmm(i),xray_type+50
                     enddo
                  endif
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_xmm(i,1),flux_xmm(i,1),
     &            FluxU_xmm(i,1),FluxL_xmm(i,1),ra_xmm(i),dec_xmm(i),poserr_xmm(i),xray_type+10
                  if (xray_type .eq. 1) then
                     do s=2,3
                        write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_xmm(i,s),flux_xmm(i,s),
     &                  FluxU_xmm(i,s),FluxL_xmm(i,s),ra_xmm(i),dec_xmm(i),poserr_xmm(i),xray_type+50
                     enddo
                  else
                     do s=2,6
                        write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_xmm(i,s),flux_xmm(i,s),
     &                     FluxU_xmm(i,s),FluxL_xmm(i,s),ra_xmm(i),dec_xmm(i),poserr_xmm(i),xray_type+50
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
            if (dist .le. min_dist_rosat)  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(i)(1:3) == 'PSR'))
     &               .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_rosat(i),dec_rosat(i),dist)
               if (dist .le. min_dist_other) found=.true.
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
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_rosat(i),flux_rosat(i),
     &             FluxU_rosat(i),FluxL_rosat(i),ra_rosat(i),dec_rosat(i),poserr_rosat(i),xray_type+10
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_rosat(i),flux_rosat(i),
     &             FluxU_rosat(i),FluxL_rosat(i),ra_rosat(i),dec_rosat(i),poserr_rosat(i),xray_type+10
               endif
            endif
         endif
      enddo

      Do i=1,iswift
         found=.false.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),ra_swift(i),dec_swift(i),dist)
            if (dist .le. min_dist_swift)  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(i)(1:3) == 'PSR'))
     &               .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_swift(i),dec_swift(i),dist)
               if (dist .le. min_dist_other) found=.true.
            endif
         enddo
         if (.not. found) THEN
            if (xrt_type(i) == 1 ) then
               xray_type=5
            else
               xray_type=9
            endif
            CALL DIST_SKY(ra_center,dec_center,ra_swift(i),dec_swift(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_swift(i,1),flux_swift(i,1),
     &             FluxU_swift(i,1),FluxL_swift(i,1),ra_swift(i),dec_swift(i),poserr_swift(i),xray_type+10
                  do s=2,4
                     write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_swift(i,s),flux_swift(i,s),
     &              FluxU_swift(i,s),FluxL_swift(i,s),ra_swift(i),dec_swift(i),poserr_swift(i),xray_type+50
                  enddo
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_swift(i,1),flux_swift(i,1),
     &              FluxU_swift(i,1),FluxL_swift(i,1),ra_swift(i),dec_swift(i),poserr_swift(i),xray_type+10
                  do s=2,4
                     write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_swift(i,s),flux_swift(i,s),
     &               FluxU_swift(i,s),FluxL_swift(i,s),ra_swift(i),dec_swift(i),poserr_swift(i),xray_type+50
                  enddo
               endif
            endif
         endif
      enddo

      Do i=1,iipc
         found=.false.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),ra_ipc(i),dec_ipc(i),dist)
            if (dist .le. min_dist_ipc)  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(i)(1:3) == 'PSR') )
     &            .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_ipc(i),dec_ipc(i),dist)
               if (dist .le. min_dist_other) found=.true.
            endif
         enddo
         if (.not. found) THEN
            xray_type=6
            CALL DIST_SKY(ra_center,dec_center,ra_ipc(i),dec_ipc(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_ipc(i),flux_ipc(i),
     &              FluxU_ipc(i),FluxL_ipc(i),ra_ipc(i),dec_ipc(i),poserr_ipc(i),xray_type+10
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_ipc(i),flux_ipc(i),
     &              FluxU_ipc(i),FluxL_ipc(i),ra_ipc(i),dec_ipc(i),poserr_ipc(i),xray_type+10
               endif
            endif
         endif
      enddo

      Do i=1,ibmw
         found=.false.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),ra_bmw(i),dec_bmw(i),dist)
            if (dist .le. min_dist_bmw)  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(i)(1:3) == 'PSR'))
     &           .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_bmw(i),dec_bmw(i),dist)
               if (dist .le. min_dist_other) found=.true.
            endif
         enddo
         if (.not. found) THEN
            xray_type=7
            CALL DIST_SKY(ra_center,dec_center,ra_bmw(i),dec_bmw(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_bmw(i),flux_bmw(i),
     &              FluxU_bmw(i),FluxL_bmw(i),ra_bmw(i),dec_bmw(i),poserr_bmw(i),xray_type+10
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_bmw(i),flux_bmw(i),
     &              FluxU_bmw(i),FluxL_bmw(i),ra_bmw(i),dec_bmw(i),poserr_bmw(i),xray_type+10
               endif
            endif
         endif
      enddo

      Do i=1,ichandra
         found=.false.
         do j=1,iradio
            call dist_sky(ra_radio(j),dec_radio(j),ra_chandra(i),dec_chandra(i),dist)
            if (dist .le. min_dist_chandra)  found=.true.
         enddo
         do j=1,iother
            if ( ( (name_other(j)(1:3) == '5BZ') .OR. (name_other(j)(1:4) == '3HSP') .or.
     &             (name_other(j)(1:6) == 'CRATES') .or. (name_other(i)(1:3) == 'PSR') )
     &           .AND. (ra_other(j) .gt. 0.) ) THEN
               CALL DIST_SKY(ra_other(j),dec_other(j),ra_chandra(i),dec_chandra(i),dist)
               if (dist .le. min_dist_other) found=.true.
            endif
         enddo
         if (.not. found) THEN
            xray_type=8
            CALL DIST_SKY(ra_center,dec_center,ra_chandra(i),dec_chandra(i),dist)
            if ((errrad .ne. 0.) .and. (errmaj .eq. 0.)) then
               if (dist .le. errrad/60.) then
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_chandra(i,1),flux_chandra(i,1),
     &     FluxU_chandra(i,1),FluxL_chandra(i,1),ra_chandra(i),dec_chandra(i),poserr_chandra(i),xray_type+10
                  do s=2,5
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_chandra(i,s),flux_chandra(i,s),
     &     FluxU_chandra(i,s),FluxL_chandra(i,s),ra_chandra(i),dec_chandra(i),poserr_chandra(i),xray_type+50
                  enddo
               endif
            else if ((errrad .eq. 0.) .and. (errmaj .ne. 0.)) then
               if (dist .le. errmaj/60.) then
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_chandra(i,1),flux_chandra(i,1),
     &     FluxU_chandra(i,1),FluxL_chandra(i,1),ra_chandra(i),dec_chandra(i),poserr_chandra(i),xray_type+10
                  do s=2,5
                  write(14,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_chandra(i,s),flux_chandra(i,s),
     &     FluxU_chandra(i,s),FluxL_chandra(i,s),ra_chandra(i),dec_chandra(i),poserr_chandra(i),xray_type+50
                  enddo
               endif
            endif
         endif
      enddo

      !write(*,*) sfound,rfound,ifound
      if (ifound .ne. sfound+rfound ) stop 'Warning, might have wrong matched number'
      Do i=1,sfound
         CALL graphic_code (xflux(i),flux_source(i)/rrconst(i),ttsource(i),code)
         write(11,'(f9.5,2x,f9.5,2x,i6)') ra_source(i),dec_source(i),int(code)
      enddo
      IF (ifound  ==  0) print *,achar(27),'[31;1m No radio/X-ray matches were found.',achar(27),'[0m'

      savemjy(1:iother)=25.
      DO l=1,iother
         IF ( ( (name_other(l)(1:3) == '5BZ') .OR. (name_other(l)(1:4) == '3HSP') .or.
     &             (name_other(l)(1:6) == 'CRATES') .or. (name_other(l)(1:3) == 'PSR'))
     &             .AND. (ra_other(l) .gt. 0.)) THEN
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
                     IF (name_other(l)(1:6) == 'CRATES') type_average = -3
                     IF (name_other(l)(1:3) == 'PSR') then
                        type_average = 0
                        code=-8888
                     endif
                     write(*,'(a,a,i4,2x,a)') name_other(l)(1:lenact(name_other(l))),
     &                   ", repeated with candidate nr.",track2(ncat),name_cat(j)(1:lenact(name_cat(j)))
                     goto 100
                  endif
               enddo
            endif
            IF (name_other(l)(1:3) == '5BZ') type_average = -6
            IF (name_other(l)(1:4) == '3HSP') type_average = -5
            IF (name_other(l)(1:6) == 'CRATES') type_average = -7
            IF (name_other(l)(1:3) == 'PSR') type_average = 0
            sfound=sfound+1
            track2(ncat)=sfound
            if (sfound .ne. 1 ) write(12,*) "===================="
            write(12,'(i4,2x,a,2(2x,f9.5),2x,a,2x,i2)') sfound,"matched source",
     &         ra_other(l),dec_other(l),'source type',type_average
            CALL DIST_SKY(ra_other(l),dec_other(l),ra_center,dec_center,dist)
            dist = dist*60
            if (type_average .eq. -7) then
               print '(a,a,i4,a,a,a,a,a,f7.3,a)', achar(27),'[35;1m Candidate nr.',sfound,
     &          ', Known flat spectrum radio source with no radio/X-ray match: ',
     &          achar(27),'[0m',name_other(l)(1:lenact(name_other(l))),' found at a distance of ',
     &                          real(dist),' arcmin '
            else if (type_average .eq. 0) then
               code=-9999
               !write(*,*) code,type_average
               write(*,'("Pulsar",2x,a,2x,f7.3,2x,"arcmin away")') name_other(l)(1:lenact(name_other(l))),dist
            else
               print '(a,a,i4,a,a,a,a,a,f7.3,a)', achar(27),'[35;1m Candidate nr.',sfound,
     &          ', Known blazar with no radio/X-ray match: ',
     &          achar(27),'[0m',name_other(l)(1:lenact(name_other(l))),' found at a distance of ',
     &                          real(dist),' arcmin '
            endif
            do j=1,iradio
               call DIST_SKY(ra_other(l),dec_other(l),ra_radio(j),dec_radio(j),dist)
               if (dist < min_dist_other) THEN
                  write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_radio(j),flux_radio(j),FluxU_radio(j),
     &          FluxL_radio(j),ra_radio(j),dec_radio(j),poserr_radio(j),radio_type(j)
                savemjy(l)=flux_radio(j)/const(j)
               endif
            enddo
            do j=1,ixmm
               if (xmm_type(i) == 1) min_dist_xmm=15./3600.
               if (xmm_type(i) == 2) min_dist_xmm=4./3600.
               call DIST_SKY(ra_other(l),dec_other(l),ra_xmm(j),dec_xmm(j),dist)
               if (dist < min_dist_xmm) THEN
                  if (xmm_type(j) == 1) then
                     xray_type=1
                     write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_xmm(j,1),
     &                flux_xmm(j,1),FluxU_xmm(j,1),FluxL_xmm(j,1),ra_xmm(j),dec_xmm(j),poserr_xmm(j),xray_type+10
                     do s=2,3
                        write(12,'(4(es10.3,2x),i2)') frequency_xmm(j,s),flux_xmm(j,s),FluxU_xmm(j,s),
     &                     FluxL_xmm(j,s),xray_type
                     enddo
                  else
                     xray_type=2
                     write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_xmm(j,1),
     &               flux_xmm(j,1),FluxU_xmm(j,1),FluxL_xmm(j,1),ra_xmm(j),dec_xmm(j),poserr_xmm(j),xray_type+10
                     do s=2,6
                        write(12,'(4(es10.3,2x),i2)') frequency_xmm(j,s),flux_xmm(j,s),FluxU_xmm(j,s),
     &                     FluxL_xmm(j,s),xray_type
                     enddo
                  endif
               endif
            enddo
            do j=1,irosat
               call DIST_SKY(ra_other(l),dec_other(l),ra_rosat(j),dec_rosat(j),dist)
               if (dist < min_dist_rosat) THEN
                  if (rosat_type(j) == 1) THEN
                     xray_type=3
                     write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_rosat(j),
     &                    flux_rosat(j),FluxU_rosat(j),FluxL_rosat(j),ra_rosat(j),dec_rosat(j),
     &                   poserr_rosat(j),xray_type+10
                  else
                     xray_type=4
                     write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_rosat(j),
     &                  flux_rosat(j),FluxU_rosat(j),FluxL_rosat(j),ra_rosat(j),dec_rosat(j),
     &                   poserr_rosat(j),xray_type+10
                  endif
               endif
            enddo
            do j=1,iswift
               call DIST_SKY(ra_other(l),dec_other(l),ra_swift(j),dec_swift(j),dist)
               if (dist < min_dist_swift) THEN
                  if (xrt_type(j) == 1) then
                     xray_type=5
                  else
                     xray_type=9
                  endif
                  write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_swift(j,1),
     &                     flux_swift(j,1),FluxU_swift(j,1),FluxL_swift(j,1),ra_swift(j),dec_swift(j),
     &                     poserr_swift(j),xray_type+10
                  do s=2,4
                     write(12,'(4(es10.3,2x),i2)') frequency_swift(j,s),flux_swift(j,s),FluxU_swift(j,s),
     &                     FluxL_swift(j,s),xray_type
                  enddo
               endif
            enddo
            do j=1,iipc
               call DIST_SKY(ra_other(l),dec_other(l),ra_ipc(j),dec_ipc(j),dist)
               if (dist < min_dist_ipc) THEN
                  xray_type=6
                  write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_ipc(j),flux_ipc(j),
     &               FluxU_ipc(j),FluxL_ipc(j),ra_ipc(j),dec_ipc(j),poserr_ipc(j),xray_type+10
               endif
            enddo
            do j=1,ibmw
               call DIST_SKY(ra_other(l),dec_other(l),ra_bmw(j),dec_bmw(j),dist)
               if (dist < min_dist_bmw) THEN
                  xray_type=7
                  write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_bmw(j),flux_bmw(j),
     &                 FluxU_bmw(j),FluxL_bmw(j),ra_bmw(j),dec_bmw(j),poserr_bmw(j),xray_type+10
               endif
            enddo
            do j=1,ichandra
               call DIST_SKY(ra_other(l),dec_other(l),ra_chandra(j),dec_chandra(j),dist)
               if (dist < min_dist_chandra) THEN
                  xray_type=8
                  write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_chandra(j,1),
     &                 flux_chandra(j,1),FluxU_chandra(j,1),FluxL_chandra(j,1),ra_chandra(j),dec_chandra(j),
     &                   poserr_chandra(j),xray_type+10
                  do s=2,5
                     write(12,'(4(es10.3,2x),i2)') frequency_chandra(j,s),flux_chandra(j,s),FluxU_chandra(j,s),
     &                     FluxL_chandra(j,s),xray_type
                  enddo
               endif
            enddo
100   continue
            if ((type_average .lt. -3)) then
                call RXgraphic_code(savemjy(l),'R',code)
                write(11,'(f9.5,2x,f9.5,2x,i6)') ra_other(l),dec_other(l),int(code+60000)
            endif
            if (type_average .ne. 0) CALL graphic_code (1.,1.,type_average,code)
            !if (type_average .eq. -7) CALL graphic_code (1.,savemjy(l),type_average,code)
            write(11,'(f9.5,2x,f9.5,2x,i6)') ra_other(l),dec_other(l),int(code)
         ENDIF
      ENDDO
      WRITE (*,*) '      '
      goto 502

cccccc for skip the phase 1
501   continue
      sfound=1
      type_average=99
      write(11,'(f9.5,2x,f9.5,2x,a)') ra_center,dec_center,"99"
      write(12,'(i4,2x,a,2(2x,f9.5),2x,a,2x,i2)') sfound,"matched source",
     &         ra_center,dec_center,'source type',type_average
      do j=1,iradio
         call DIST_SKY(ra_center,dec_center,ra_radio(j),dec_radio(j),dist)
         if (dist .lt. 18./3600. ) then !18 arcsec for radio sources
            write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_radio(j),flux_radio(j),FluxU_radio(j),
     &          FluxL_radio(j),ra_radio(j),dec_radio(j),poserr_radio(j),radio_type(j)
         endif
      enddo
      do j=1,ixmm
         if (xmm_type(j) == 1) min_dist_xmm=15./3600.
         if (xmm_type(j) == 2) min_dist_xmm=4./3600.
         call DIST_SKY(ra_center,dec_center,ra_xmm(j),dec_xmm(j),dist)
         if ( dist .lt. min_dist_xmm ) then
            if (xmm_type(j) == 1) then
               xray_type=1
               write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_xmm(j,1),
     &                flux_xmm(j,1),FluxU_xmm(j,1),FluxL_xmm(j,1),ra_xmm(j),dec_xmm(j),poserr_xmm(j),xray_type+10
               do s=2,3
                  write(12,'(4(es10.3,2x),i2)') frequency_xmm(j,s),flux_xmm(j,s),FluxU_xmm(j,s),
     &                     FluxL_xmm(j,s),xray_type
               enddo
            else
               xray_type=2
               write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_xmm(j,1),
     &               flux_xmm(j,1),FluxU_xmm(j,1),FluxL_xmm(j,1),ra_xmm(j),dec_xmm(j),poserr_xmm(j),xray_type+10
               do s=2,6
                  write(12,'(4(es10.3,2x),i2)') frequency_xmm(j,s),flux_xmm(j,s),FluxU_xmm(j,s),
     &                     FluxL_xmm(j,s),xray_type
               enddo
             endif
          endif
       enddo
       do j=1,irosat
         call DIST_SKY(ra_center,dec_center,ra_rosat(j),dec_rosat(j),dist)
         if ( dist .lt. min_dist_rosat ) then
            if (rosat_type(j) == 1) THEN
               xray_type=3
               write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_rosat(j),
     &                    flux_rosat(j),FluxU_rosat(j),FluxL_rosat(j),ra_rosat(j),dec_rosat(j),
     &                   poserr_rosat(j),xray_type+10
            else
               xray_type=4
               write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_rosat(j),
     &                  flux_rosat(j),FluxU_rosat(j),FluxL_rosat(j),ra_rosat(j),dec_rosat(j),
     &                   poserr_rosat(j),xray_type+10
            endif
         endif
      enddo
      do j=1,iswift
         call DIST_SKY(ra_center,dec_center,ra_swift(j),dec_swift(j),dist)
         if ( dist .lt. min_dist_swift ) then
            if (xrt_type(j) == 1) THEN
               xray_type=5
            else
               xray_type=9
            endif
            write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_swift(j,1),
     &                     flux_swift(j,1),FluxU_swift(j,1),FluxL_swift(j,1),ra_swift(j),dec_swift(j),
     &                     poserr_swift(j),xray_type+10
            do s=2,4
               write(12,'(4(es10.3,2x),i2)') frequency_swift(j,s),flux_swift(j,s),FluxU_swift(j,s),
     &                     FluxL_swift(j,s),xray_type
            enddo
         endif
      enddo
      do j=1,iipc
         call DIST_SKY(ra_center,dec_center,ra_ipc(j),dec_ipc(j),dist)
         if ( dist .lt. min_dist_ipc ) then
            xray_type=6
            write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_ipc(j),flux_ipc(j),
     &               FluxU_ipc(j),FluxL_ipc(j),ra_ipc(j),dec_ipc(j),poserr_ipc(j),xray_type+10
         endif
      enddo
      do j=1,ibmw
         call DIST_SKY(ra_center,dec_center,ra_bmw(j),dec_bmw(j),dist)
         if ( dist .lt. min_dist_bmw ) then
            xray_type=7
            write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_bmw(j),flux_bmw(j),
     &                 FluxU_bmw(j),FluxL_bmw(j),ra_bmw(j),dec_bmw(j),poserr_bmw(j),xray_type+10
         endif
      enddo
      do j=1,ichandra
         call DIST_SKY(ra_center,dec_center,ra_chandra(j),dec_chandra(j),dist)
         if ( dist .lt. min_dist_chandra ) then
            xray_type=8
            write(12,'(4(es10.3,2x),2(f9.5,2x),f7.3,2x,i2)') frequency_chandra(j,1),
     &                 flux_chandra(j,1),FluxU_chandra(j,1),FluxL_chandra(j,1),ra_chandra(j),dec_chandra(j),
     &                   poserr_chandra(j),xray_type+10
            do s=2,5
               write(12,'(4(es10.3,2x),i2)') frequency_chandra(j,s),flux_chandra(j,s),FluxU_chandra(j,s),
     &                     FluxL_chandra(j,s),xray_type
            enddo
         endif
      enddo
502   close(11)
      close(12)
      close(13)
      close(14)
      END
c
      SUBROUTINE print_results (ratio,ra,dec,flux_radio,radio_type,xray_type,
     &                          flux,const,ra_center,dec_center,source_type)
      IMPLICIT none
      INTEGER*4 rah,ram,id,dm,radio_type,xray_type,source_type
      REAL*4 ratio,rasec,decsec,const,arx,lognupeak,flux_radio
      REAL*4 flux
      REAL*8 ra,dec,ra_center,dec_center,dist
      CHARACTER*1 sign
      CHARACTER*10 xmission,radio_survey
      CHARACTER*80 type
      ratio=flux/flux_radio
      call chra(ra,rah,ram,rasec,1)
      call chdec(dec,id,dm,decsec,1)
      sign(1:1) =' '
      if (dec < 0.0) sign(1:1)='-'
      type = ' (type unknown) '
c      IF (ratio < 0.) RETURN 
      IF (radio_type == 2) THEN
         radio_survey='NVSS'
      ELSE IF (radio_type == 1) THEN
         radio_survey='FIRST'
      ELSE IF (radio_type == 3) THEN
         radio_survey='SUMSS'
      ELSE
         radio_survey='UNKNOWN'
      ENDIF
      IF (xray_type == 1) THEN
         xmission='XMMSLEW'
      ELSE IF (xray_type == 2) THEN
         xmission='3XMM'
      ELSE IF (xray_type == 3) THEN
         xmission='RASS'
      ELSE IF (xray_type == 4) THEN
         xmission='WGA'
      ELSE IF (xray_type == 5) THEN
         xmission='XRT'
      ELSE IF (xray_type == 6) THEN
         xmission='IPC'
      ELSE IF (xray_type == 7) THEN
         xmission='BMW'
      ELSE IF (xray_type == 8) THEN
         xmission='CHANDRA'
      ELSE IF (xray_type == 9) THEN
         xmission='XRTDEEP'
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
     &      flux_radio/const,ratio,arx,lognupeak,type(1:len_trim(type)),dist
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
     &       flux_radio/const,ratio,arx,type(1:len_trim(type)),dist
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
     &       flux_radio/const,type(1:len_trim(type)),dist

      ENDIF
      RETURN
      END
C
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
      INTEGER*4 radio_component,x_ray_component,source_type,temp
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


