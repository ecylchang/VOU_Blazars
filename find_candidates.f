      PROGRAM find_candidates
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
      INTEGER*4 ier, lu_in, ia, radio_type(5000),xray_type,lu_output, in,rfound,ir100found,s
      INTEGER*4 xmm_type(5000),rosat_type(1000),swift_type(5000),no_found,sfound,nrep(5000)
      INTEGER*4 bmw_type(5000),ipc_type(200),lenact,source_type,type_average,track(5000)
      INTEGER*4 iradio,ixmm,irosat,iswift,iipc,iother,k,ix,types(0:5),i4p8,pccs100_type(200)
      INTEGER*4 is2,ie2,iir,iuv,ixrt,igam,iuvfound,iirfound,ixrtfound,igamfound
      INTEGER*4 rah, ram, id, dm ,is,ie, i, j,ibmw,ifound,ra_index(5000),l,t(5000)
      INTEGER*4 iusno, iofound, length,ialphar,iofound_index(100),ipccs100,ioptmatch
      REAL*8 ra_other(10000),dec_other(10000),ra_usno(1000),dec_usno(1000)
      REAL*8 ra_radio(5000),dec_radio(5000),ra_xmm(5000),dec_xmm(5000)
      REAL*8 ra_rosat(1000),dec_rosat(1000),ra, dec,min_dist_gam
      REAL*8 ra_swift(5000),dec_swift(5000),ra_bmw(500),dec_bmw(500)
      REAL*8 ra_ipc(200),dec_ipc(200),dist,ra_center, dec_center
      REAL*8 ra_opt_cand(100),dec_opt_cand(100)
      REAL*8 ra_pccs100(200),dec_pccs100(200),ra_gam(100),dec_gam(100)
      REAL*8 ra_4p8(1000),dec_4p8(1000),dist_opt_cand(100)
      REAL*8 ra_ir(1000),dec_ir(1000),ra_uv(1000),dec_uv(1000),ra_xrt(5000),dec_xrt(5000)
      REAL*4 flux_radio(5000),flux_xmm(5000),flux_rosat(1000),radian,aox,a100x
      REAL*4 flux_swift(5000),flux_ipc(200),flux_bmw(500),flux_x,nh,aro,arx,alpho
      REAL*4 flux_4p8(1000),alphar,flux_usno(1000,5),frequency_usno(1000,5)
      REAL*4 min_dist_rosat,min_dist_xmm,rasec,decsec,min_dist_ipc,min_dist2opt
      REAL*4 min_dist_other,min_dist_swift,min_dist_bmw,min_dist_4p8,min_dist_uv,min_dist_ir
      REAL*4 flux2nufnu_nvss,flux2nufnu_rosat,flux2nufnu_xmm,min_dist
      REAL*4 flux2nufnu_swift,flux2nufnu_ipc,code,flux2nufnu_4p8,aalphar
      REAL*4 flux2nufnu_bmw,flux2nufnu_rxs,ratio,const(5000),flux2nufnu_sumss
      REAL*4 aro_opt_cand(100),aox_opt_cand(100),mag_opt_cand(100),mag_usno(1000,5)
      REAL*4 min_dist_pccs100,flux2nufnu_pccs100,flux_pccs100(200),min_dist_cluster
      REAL*4 Umag,Bmag,Vmag,Rmag,Imag,uumag,ggmag,rrmag,iimag,zzmag,mag_uv(1000,6)
      REAL*4 flux_ir(1000,4),w1mag,w2mag,w3mag,w4mag,Jmag,Hmag,Kmag,frequency_ir(1000,4)
      REAL*4 fuvmag,nuvmag,uvvmag,uvbmag,uvumag,uvw1mag,uvm2mag,uvw2mag,flux_uv(1000,6)
      REAL*4 flux_xrt(5000,4),flux_gam(100,7),slope_gam(100,2),frequency_uv(1000,6),frequency_pccs100(200)
      REAL*4 frequency_gam(100,7),flux_1kev(5000,5000),frequency_xrt(4),frequency_radio(5000)
      REAL*4 auvx,aruv,airx,arir,aswift,alphauv,xflux(5000),rflux(5000),frequency_4p8(1000)
      CHARACTER*1 sign
      CHARACTER*30 name_other(10000),name2_other(10000)
      CHARACTER*8 opt_type(1000),opt_type_cand(100),uv_type(1000),ir_type(1000),gam_type(100)
      CHARACTER*8 uv_type_cand(100),ir_type_cand(100)
      CHARACTER*300 input_file,string,output_file
      LOGICAL there,ok,found 
      ok = .TRUE. 
      found = .FALSE.
      nrep(1:5000)=1
      ifound = 0
      sfound = 0
      rfound = 0
      iradio=0
      ixmm=0
      i4p8=0
      iusno=0
      irosat=0
      iswift=0
      iipc=0
      ibmw = 0
      iother=0
      ipccs100=0
      iir=0
      iuv=0
      ixrt=0
      igam=0
      radian = 45.0/atan(1.0)
      flux2nufnu_nvss=1.4e9*1.e-26
      flux2nufnu_sumss=8.43e8*1.e-26 !assumed radio alpha=0.2 !f_0.8 to f_1.4
      flux2nufnu_4p8=4.8e9*1.e-26
c approximate flux conversions from cts/s to erg/cm2/s at 1 kev (NH=5.e20)
      flux2nufnu_rosat=7.e-12
      flux2nufnu_xmm=7.e-13
      flux2nufnu_swift=9.e-12
      flux2nufnu_ipc=1.2e-11
      flux2nufnu_bmw=1.8e-11
      flux2nufnu_rxs=1.0
      flux2nufnu_pccs100=1.e-15
      sign=' '
      min_dist_ipc=50./3600.
      min_dist_rosat=40./3600.
c 40 arcsecs
      min_dist_xmm=15./3600.
      min_dist_swift=10./3600.
      min_dist_bmw=10./3600.
      min_dist_4p8=30./3600.
      min_dist_pccs100=180./3600.
c 10 arcsecs
      min_dist_other=15./3600.
      min_dist_cluster=60./3600.
c 15 arcsecs
      min_dist2opt = 5./3600.
      min_dist_uv=10./3600.
      min_dist_ir=10./3600.
      min_dist_gam=10./60. !20 arcmin
      CALL rdforn(string,length)
      IF ( length.NE.0 ) THEN
         CALL rmvlbk(string)
         input_file=string(1:len_trim(string))
      ELSE 
         WRITE (*,'('' Enter query results file '',$)')
         READ (*,'(a)') input_file
      ENDIF
c      WRITE (*,'('' Enter output file '',$)')
c      READ (*,'(a)') output_file
      output_file='find_out.txt'
c      CALL getlun(lu_in)
      lu_in = 10
      lu_output = 11
      in = index(input_file(1:lenact(input_file)),'.')
      IF (in == 0) input_file(lenact(input_file)+1:lenact(input_file)+4) = '.csv' 
      INQUIRE (FILE=input_file,EXIST=there)
      IF (.NOT.there) THEN
         write (*,'('' file '',a,'' not found '')')
     &     input_file(1:lenact(input_file))
         STOP
      ENDIF
      open(lu_in,file=input_file,status='old',iostat=ier)
c      CALL getlun(lu_output)
      open(lu_output,file=output_file(1:lenact(output_file)),status='unknown',iostat=ier)
      open(12,file='Sed.txt')
      IF (ier.NE.0) THEN
        write (*,*) ' Error ',ier,' opening file ', input_file
      ENDIF
      ok = . TRUE.
      READ(lu_in,'(a)') string 
      READ(lu_in,'(a)') string 
      is = index(string(1:len(string)),'=') +1
      ie = index(string(1:len(string)),'(') -1
      read(string(is:ie),*) ra_center
      is = index(string(is:len(string)),'=') +is  
      ie = index(string(ie+2:len(string)),'(') +ie 
      read(string(is:ie),*) dec_center
c      READ(lu_in,'(a)') string !begin reading nh
c      is = index(string(1:len(string)),'=') +1
c      ie = index(string(1:len(string)),'(') -1
c      read(string(is:ie),*) nh
c      write(*,*) nh ! print out nh, end of reading nh
      nh=5.e20
      READ(lu_in,'(a)') string
      DO WHILE(ok)
         READ(lu_in,'(a)',end=99) string 
         is=index(string(1:len(string)),';') 
         IF (is == 0) GOTO 97
         ie=index(string(is+1:len(string)),';')+is
         read(string(is+1:is+2),'(i2)')rah
         read(string(is+4:is+5),'(i2)')ram
         read(string(is+7:ie-1),*) rasec
         is=ie
         ie=index(string(is+1:len(string)),';')+is
         sign=string(is+1:is+1)
         IF (sign(1:1) == '-') is=is+1 
         IF (sign(1:1) == '+') is=is+1 
         read(string(is+1:is+2),'(i2)') id
         read(string(is+4:is+5),'(i2)') dm
         read(string(is+7:ie-1),*) decsec
         call chra(ra,rah,ram,rasec,0)
         call chdec(dec,id,dm,decsec,0)
         IF (sign(1:1) == '-') dec=-abs(dec) 
         is=index(string(1:len(string)),'=')
         ie=index(string(is+2:len(string)),' ')+is+1
         IF ( (string(1:4) == 'NVSS') .OR.
     &        (string(1:5) == 'FIRST') .OR.
     &        (string(1:5) == 'SUMSS') ) THEN
            iradio=iradio+1
            IF (string(1:4) == 'NVSS') radio_type(iradio) = 1
            IF (string(1:5) == 'FIRST') radio_type(iradio) = 2
            IF (string(1:5) == 'SUMSS') radio_type(iradio) = 3
            IF (iradio > 5000) Stop 'Too many NVSS/SUMSS points'
            ra_radio(iradio)=ra
            dec_radio(iradio)=dec
            IF (is == 0) THEN 
               flux_radio(iradio)=-1. 
            ELSE
               read(string(is+1:ie-1),*)flux_radio(iradio) 
            ENDIF
            IF (string(1:5) == 'SUMSS') THEN
               const(iradio) = flux2nufnu_sumss
               flux_radio(iradio)=flux_radio(iradio)*flux2nufnu_sumss
               frequency_radio(iradio)=8.43E8
            ELSE
               const(iradio) = flux2nufnu_nvss
               flux_radio(iradio)=flux_radio(iradio)*flux2nufnu_nvss
               frequency_radio(iradio)=1.4e9
            ENDIF
         ELSE IF ( (string(1:3) == 'PMN') .OR. 
     &             (string(1:6) == 'ATPMN') .OR.
     &             (string(1:4) == 'GB87') .OR.
     &             (string(1:3) == 'GB6') )  THEN
            i4p8=i4p8+1
            IF (i4p8 > 1000) Stop 'Too many PMN points'
            ra_4p8(i4p8)=ra
            dec_4p8(i4p8)=dec
            IF (is == 0) THEN 
               flux_4p8(i4p8)=-1. 
            ELSE
               read(string(is+1:ie-1),*)flux_4p8(i4p8) 
            ENDIF
            flux_4p8(i4p8)=flux_4p8(i4p8)*flux2nufnu_4p8
         ELSE IF ( (string(1:9) == 'PCCS1F100') .OR.
     &             (string(1:9) == 'PCCS2F100') .OR.
     &             (string(1:9) == 'PCCS1F143') .OR.
     &             (string(1:9) == 'PCCS2F143') )  THEN
            ipccs100=ipccs100+1
            IF (ipccs100 > 200) Stop 'Too many PCCS points'
            ra_pccs100(ipccs100)=ra
            dec_pccs100(ipccs100)=dec
            IF (is == 0) THEN 
               flux_pccs100(ipccs100)=-1. 
            ELSE if ( (string(1:9) == 'PCCS1F100') .OR. (string(1:9) == 'PCCS2F100')) then
               read(string(is+1:ie-1),*)flux_pccs100(ipccs100)
               frequency_pccs100(ipccs100)=1.0e11
               pccs100_type(ipccs100)=100
            ELSE if ( (string(1:9) == 'PCCS1F143') .OR. (string(1:9) == 'PCCS2F143')) then
               read(string(is+1:ie-1),*)flux_pccs100(ipccs100)
               flux_pccs100(ipccs100)=flux_pccs100(ipccs100)*1.43 !100 to 143
               frequency_pccs100(ipccs100)=1.43e11
               pccs100_type(ipccs100)=143
            ENDIF
            flux_pccs100(ipccs100)=flux_pccs100(ipccs100)*flux2nufnu_pccs100
         ELSE IF ((string(1:4) == 'WISE') .OR.
     &             (string(1:5) == '2MASS') ) THEN
            iir=iir+1
            IF (iir > 1000) Stop 'Too many Infrared points'
            ra_ir(iir)=ra
            dec_ir(iir)=dec
            IF (is == 0) THEN
            flux_ir(iir,1:4)=-1.
            ELSE IF (string(1:4) == 'WISE') then
               read(string(is+1:ie-1),*) w1mag
               is2=index(string(ie+1:len(string)),'=')+ie
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) w2mag
               is2=index(string(ie2:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) w3mag
               is2=index(string(ie2:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) w4mag
               call mag2flux (nh,w1mag,'ww1',flux_ir(iir,1),frequency_ir(iir,1))
               call mag2flux (nh,w2mag,'ww2',flux_ir(iir,2),frequency_ir(iir,2))
               call mag2flux (nh,w3mag,'ww3',flux_ir(iir,3),frequency_ir(iir,3))
               call mag2flux (nh,w4mag,'ww4',flux_ir(iir,4),frequency_ir(iir,4))
               ir_type(iir)='WISE'
            ELSE IF (string(1:5) == '2MASS') then
               read(string(is+1:ie-1),*) Jmag
               is2=index(string(ie+1:len(string)),'=')+ie
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) Hmag
               is2=index(string(ie2:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) Kmag
               ir_type(iir)='2MASS'
               CALL  mag2flux (nh,Jmag,'J  ',flux_ir(iir,4),frequency_ir(iir,4))
               call  mag2flux (nh,Hmag,'H  ',flux_ir(iir,3),frequency_ir(iir,3))
               call  mag2flux (nh,Kmag,'K  ',flux_ir(iir,2),frequency_ir(iir,2))
            ENDIF
         ELSE IF ( (string(1:4) == 'USNO') .OR.
     &             (string(1:4) == 'SDSS') ) THEN
            iusno = iusno + 1 
            IF (iusno > 1000) Stop 'Too many USNO points'
            ra_usno(iusno)=ra
            dec_usno(iusno)=dec
            is2=index(string(ie+1:len(string)),'=')+ie
            ie2=index(string(is2+2:len(string)),' ')+is2+1
            IF (is == 0) THEN
               flux_usno(iusno,1:2)=-1.
            ELSE IF (string(1:4) == 'USNO') THEN !tell if the data is from SDSS or USNO
               read(string(is+1:ie-1),*) Rmag
               mag_usno(iusno,1) = Rmag
               read(string(is2+1:ie2-1),*) Bmag
               mag_usno(iusno,2) = Bmag
               CALL  mag2flux (nh,Rmag,'R  ',flux_usno(iusno,1),frequency_usno(iusno,1))
               CALL  mag2flux (nh,Bmag,'B  ',flux_usno(iusno,2),frequency_usno(iusno,2))
               opt_type(iusno)='USNO'
            ELSE
               read(string(is+1:ie-1),*) rrmag
               mag_usno(iusno,1) = rrmag
               read(string(is2+1:ie2-1),*) uumag
               mag_usno(iusno,2) = uumag
               CALL  mag2flux (nh,rrmag,'r  ',flux_usno(iusno,1),frequency_usno(iusno,1))
               CALL  mag2flux (nh,uumag,'u  ',flux_usno(iusno,2),frequency_usno(iusno,2))
               opt_type(iusno)='SDSS'
            ENDIF
         ELSE IF  ( (string(1:5) == 'GALEX') .OR.
     &              (string(1:4) == 'UVOT') )  THEN
            iuv=iuv+1
            if (iuv > 1000) Stop 'Too many UV points'
            ra_uv(iuv)=ra
            dec_uv(iuv)=dec
            IF (is == 0) THEN
               flux_uv(iuv,1:6)=-1.
            ELSE IF (string(1:5) == 'GALEX') then
              read(string(is+1:ie-1),*) fuvmag
               is2=index(string(ie+1:len(string)),'=')+ie
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) nuvmag
               CALL  mag2flux (nh,17.,'fuv',flux_uv(iuv,4),frequency_uv(iuv,4))
               flux_uv(iuv,4)=fuvmag
               uv_type(iuv)='GALEX'
            ELSE IF (string(1:4) == 'UVOT') then
               read(string(is+1:ie-1),*) uvvmag
               is2=index(string(ie+1:len(string)),'=')+ie
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) uvbmag
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) uvumag
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) uvw1mag
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) uvm2mag
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) uvw2mag
               CALL  mag2flux (nh,18.,'su ',flux_uv(iuv,1),frequency_uv(iuv,1))
               CALL  mag2flux (nh,18.,'sb ',flux_uv(iuv,2),frequency_uv(iuv,2))
               CALL  mag2flux (nh,18.,'sv ',flux_uv(iuv,3),frequency_uv(iuv,3))
               CALL  mag2flux (nh,18.,'sw1',flux_uv(iuv,4),frequency_uv(iuv,4))
               CALL  mag2flux (nh,18.,'sm2',flux_uv(iuv,5),frequency_uv(iuv,5))
               CALL  mag2flux (nh,18.,'sw2',flux_uv(iuv,6),frequency_uv(iuv,6))
               flux_uv(iuv,1)=uvumag
               flux_uv(iuv,2)=uvbmag
               flux_uv(iuv,3)=uvvmag
               flux_uv(iuv,4)=uvw1mag
               flux_uv(iuv,5)=uvm2mag
               flux_uv(iuv,6)=uvw2mag
               uv_type(iuv)='UVOT'
            ENDIF
         ELSE IF ( (string(1:5) == 'XMMSL') .OR.
     &             (string(1:6) == 'TWOXMM') )  THEN
            ixmm=ixmm+1
            IF (string(1:5) == 'XMMSL') THEN
              xmm_type(ixmm)=1
            ELSE IF (string(1:6) == 'TWOXMM') THEN
              xmm_type(ixmm)=2
            ENDIF
            IF (ixmm > 5000) Stop 'Too many XMM points'
            ra_xmm(ixmm)=ra
            dec_xmm(ixmm)=dec
            IF (is == 0) THEN 
               flux_xmm(ixmm)=-1. 
            ELSE
               read(string(is+1:ie-1),*)flux_xmm(ixmm) 
            ENDIF
            flux_xmm(ixmm)=flux_xmm(ixmm)*flux2nufnu_xmm
         ELSE IF ( (string(1:4) == 'RASS') .OR.
     &             (string(1:3) == 'RXS')  .OR. 
     &             (string(1:3) == 'WGA') ) THEN
            irosat=irosat+1
            IF (irosat > 1000) Stop 'Too many RASS points'
            IF (string(1:4) == 'RASS') THEN 
               rosat_type(irosat)=1
            ELSE IF (string(1:3) == 'RXS') THEN
               rosat_type(irosat)=2
            ELSE IF (string(1:3) == 'WGA') THEN
               rosat_type(irosat)=3
            ENDIF
            ra_rosat(irosat)=ra
            dec_rosat(irosat)=dec
            IF (is == 0) THEN 
               flux_rosat(irosat)=-1. 
            ELSE
               read(string(is+1:ie-1),*)flux_rosat(irosat) 
            ENDIF
            IF (string(1:3) == 'RXS') THEN
               flux_rosat(irosat)=flux_rosat(irosat)*flux2nufnu_rxs
            ELSE
               flux_rosat(irosat)=flux_rosat(irosat)*flux2nufnu_rosat
            ENDIF
         ELSE IF ( (string(1:5) == 'SWXRT') .OR. 
     &             (string(1:3) == 'XRT') ) THEN
            iswift=iswift+1
            IF (string(1:5) == 'SWXRT') THEN 
               swift_type(iswift)=1
            ELSE IF (string(1:3) == 'XRT') THEN
               swift_type(iswift)=2
            ENDIF 
            IF (iswift > 5000) THEN 
               write(*,*) 'Too many swift points', iswift
               Stop 
            ENDIF
            ra_swift(iswift)=ra
            dec_swift(iswift)=dec
            IF (is == 0) THEN 
               flux_swift(iswift)=-1. 
            ELSE
               read(string(is+1:ie-1),*)flux_swift(iswift) 
            ENDIF
            flux_swift(iswift)=flux_swift(iswift)*flux2nufnu_swift
         ELSE IF (string(1:3) == 'IPC') THEN
            iipc=iipc+1
            ipc_type(iipc)=1
            IF (iipc > 200) Stop 'Too many Einstein IPC points'
            ra_ipc(iipc)=ra
            dec_ipc(iipc)=dec
            IF (is == 0) THEN 
               flux_ipc(iipc)=-1. 
            ELSE
               read(string(is+1:ie-1),*)flux_ipc(iipc) 
            ENDIF
            flux_ipc(iipc)=flux_ipc(iipc)*flux2nufnu_ipc
         ELSE IF (string(1:3) == 'BMW') THEN
            ibmw=ibmw+1
            bmw_type(ibmw)=1
            IF (ibmw > 500) Stop 'Too many ROSAT-BMW points'
            ra_bmw(ibmw)=ra
            dec_bmw(ibmw)=dec
            IF (is == 0) THEN 
               flux_bmw(ibmw)=-1. 
            ELSE
               read(string(is+1:ie-1),*)flux_bmw(ibmw) 
            ENDIF
            flux_bmw(ibmw)=flux_bmw(ibmw)*flux2nufnu_bmw
         ELSE IF (string(1:4) == 'SXPS') then
            ixrt=ixrt+1
            If (ixrt > 5000) stop 'Too many 1SXPS points'
            ra_xrt(ixrt)=ra
            dec_xrt(ixrt)=dec
            frequency_xrt(1:4)=(/4.188E+17,1.325E+17,3.419E+17,1.081E+18/)
            if (is == 0) then
               flux_xrt(ixrt,1:4)=-1.
            ELSE
               read(string(is+1:ie-1),*) flux_xrt(ixrt,1)
               is2=index(string(ie+1:len(string)),'=')+ie
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) flux_xrt(ixrt,2)
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) flux_xrt(ixrt,3)
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) flux_xrt(ixrt,4)
            endif
         ELSE IF ((string(1:4) == '2FHL') .or.
     &      (string(1:4) == '3FGL'))then
            igam=igam+1
            If (igam > 100) stop 'Too many Gamma-ray points'
            ra_gam(igam)=ra
            dec_gam(igam)=dec
            if (is == 0) then
               flux_gam(igam,1:4)=-1
               slope_gam(igam,1:2)=-1
            Else if (string(1:4) == '2FHL') then
               read(string(is+1:ie-1),*) flux_gam(igam,1)
               is2=index(string(ie+1:len(string)),'=')+ie
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) flux_gam(igam,2)
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) flux_gam(igam,3)
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) flux_gam(igam,4)
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) slope_gam(igam,1)
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) slope_gam(igam,2)
               frequency_gam(igam,1:7)=(/50.,92.,316.,1081.,-1.,-1.,-1./)!GeV
               frequency_gam(igam,1:7)=frequency_gam(igam,1:7)*1.e9*(1.4e-19/6.6e-34)
               gam_type(igam)='2FHL'
            ELSE IF (string(1:4) == '3FGL') then
               read(string(is+1:ie-1),*) flux_gam(igam,1)
               is2=index(string(ie+1:len(string)),'=')+ie
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) flux_gam(igam,2)
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) flux_gam(igam,3)
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) flux_gam(igam,4)
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) flux_gam(igam,5)
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) flux_gam(igam,6)
               is2=index(string(ie2+1:len(string)),'=')+ie2
               ie2=index(string(is2+2:len(string)),' ')+is2+1
               read(string(is2+1:ie2-1),*) slope_gam(igam,1)
               frequency_gam(igam,1:7)=(/1.,0.2,0.6,2.,6.,60.,-1./) !GeV
               frequency_gam(igam,1:7)=frequency_gam(igam,1:7)*1.e9*(1.4e-19/6.6e-34)
               gam_type(igam)='3FGL'
            ENDIF
         ELSE
            iother=iother+1
            IF (iother > 10000) Stop 'Too many catalogued sources'
            ra_other(iother)=ra
            dec_other(iother)=dec
            is = index(string(1:len(string)),';')
            ie = 1 
            ia = 1 
            DO WHILE(ia.NE.0) 
               ia = index(string(ie:len(string)),';')
               ie = ia+ie
            ENDDO
            IF (is.LE.30) THEN 
              name_other(iother)(1:is) = string(1:is-1)
            ELSE
              name_other(iother)(1:30) = string(1:30)
            ENDIF
            ia = min((len(string)-ie),30) 
            name2_other(iother)(1:ia+1) = string(ie:ie+ia+1)
         ENDIF
 97      CONTINUE
      ENDDO 
 99   CONTINUE
      CLOSE (lu_in)
      CALL indexx (iradio,ra_radio,ra_index)
      DO j=1,iradio
         DO i =0,5
           types(i) = 0
         ENDDO
         found = .FALSE.
         flux_x = 0.
         type_average = 0.
         ix = 0
         k = ra_index(j)
         call chra(ra_radio(k),rah,ram,rasec,1)
         call chdec(dec_radio(k),id,dm,decsec,1)
         DO i=1,ixmm
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
               IF (flux_xmm(i) > 0. ) THEN 
                  flux_x = flux_x + flux_xmm(i)
                  ix = ix +1
                  flux_1kev(ix,k)=flux_xmm(i)
               ENDIF
               CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),xray_type,
     &                             flux_xmm(i),const(k),ra_center,dec_center,source_type)
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
               ELSE IF (rosat_type(i) == 3) THEN
                 xray_type = 5
               ENDIF 
               flux_x = flux_x + flux_rosat(i) 
               ix = ix +1
               flux_1kev(ix,k)=flux_rosat(i)
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
               found = .TRUE.
               IF (swift_type(i) == 1) THEN 
                 xray_type = 6
               ELSE IF (swift_type(i) == 2) THEN 
                 xray_type = 7
               ENDIF
               flux_x = flux_x + flux_swift(i) 
               ix = ix +1
               flux_1kev(ix,k)=flux_swift(i)
               CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_swift(i),const(k),ra_center,dec_center,source_type)
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
               IF (ipc_type(i) == 1) THEN 
                 xray_type = 8
               ENDIF
               flux_x = flux_x + flux_ipc(i) 
               ix = ix +1
               flux_1kev(ix,k)=flux_ipc(i)
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
               IF (bmw_type(i) == 1) THEN 
                 xray_type = 9
               ENDIF
               flux_x = flux_x + flux_bmw(i) 
               ix = ix +1
               flux_1kev(ix,k)=flux_bmw(i)
               CALL print_results (ratio,ra_radio(k),dec_radio(k),flux_radio(k),radio_type(k),
     &                             xray_type,flux_bmw(i),const(k),ra_center,dec_center,source_type)
               IF (source_type .GE. 0) THEN 
                  types(source_type) = types(source_type) + 1
               ENDIF
            ENDIF
         ENDDO
         !write(*,*) const
         IF (found) THEN 
            ifound = ifound +1
            write(*,*) 'radio counterpart',j,k,t(ifound-1),radio_type(k),flux_radio(k)/const(k)
            if (ifound .ne. 1) then
               do i = 1,ifound-1
                  call DIST_SKY(ra_radio(k),dec_radio(k),ra_radio(t(i)),dec_radio(t(i)),dist)
                  if (dist*3600 .lt. 6.) then
                     rfound=rfound+1
                     IF ( ix.NE.0 ) THEN
                     flux_x = flux_x/float(ix)
                     ELSE
                     flux_x = 0.
                     ENDIF
                     if (rfound .eq. 1) then
                        track(ifound)=i
                     else
                        track(ifound)=track(i)
                     endif
                     if (flux_x .ne. xflux(track(ifound)) )  write(*,*) '!!!Warning, check X-ray counterpart!'
                     !rflux(track(ifound))=rflux(track(ifound))+flux_radio(k)
                     !nrep(track(ifound))=nrep(track(ifound))+1
                     !rflux(track(ifound))=rflux(track(ifound))/nrep(track(ifound))
                     write(*,'(15x,"Repeated radio counterpart, ",f5.3,2x,"arcsec away from the matched nr.",2x,i2)')
     &                   dist*3600, track(ifound)
                     write(*,*) '     '
                     write(12,*) "===================="
                     write(12,'(i4,2x,a,2(2x,f8.4))') track(ifound),"matched source",ra_radio(k),dec_radio(k)!,type_average
                     write(12,'(es9.3,2x,es9.3)') frequency_radio(k),flux_radio(k)
                     t(ifound)=k
                     goto 98 !repeated radio counterpart, so end the loop
                  endif
               enddo
            endif
            sfound=sfound+1
            track(ifound)=sfound
            rflux(sfound)=flux_radio(k)
            print *,  achar(27),'[31;1m Match nr. ',sfound,achar(27),'[0m ra dec',ra_radio(k),',',dec_radio(k)
            IF ( ix.NE.0 ) THEN 
               flux_x = flux_x/float(ix)
            ELSE
               flux_x = 0.
            ENDIF
            xflux(sfound)=flux_x
            no_found = 0
            DO i = 0,5
              IF (types(i) > no_found) THEN
                 no_found=types(i) 
                 type_average = i 
              ENDIF
            ENDDO
            alphar = -99.0
            a100x = -99.9
            ialphar = 0
            aalphar = 0.
            ir100found=0
            if (sfound .ne. 1 ) write(12,*) "===================="
            write(12,'(i4,2x,a,2(2x,f8.4),2x,i2)') sfound,"matched source",ra_radio(k),dec_radio(k),type_average
            write(12,'(es9.3,2x,es9.3)') frequency_radio(k),flux_radio(k)
            do i=1,ix
               write(12,'("2.418E+17",2x,es9.3)') flux_1kev(i,k)
            enddo
            write(*,*) '.......................100 GHz Radio........................'
            DO i=1,ipccs100
               CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_pccs100(i),
     &                       dec_pccs100(i),dist)
               IF (dist < min_dist_pccs100) THEN 
                 a100x = 1. - log10(flux_pccs100(i)/flux_x)/log10(1.e11/2.418e17)
                 write(*,'(6x,''100 GHz - X-ray slope: '',f5.2,",",2x,f6.3," arcmin away")') a100x,dist*60
                 ir100found=ir100found+1
                    write(12,'(es9.3,2x,es9.3)') frequency_pccs100(i),flux_pccs100(i)
               ENDIF
            ENDDO
            if (ir100found .eq. 0) write(*,'(" No 100 GHz detection within",f5.0,2x,"arcmin")') min_dist_pccs100*60
            write(*,*) '........................5 GHz Radio......................'
            DO i=1,i4p8
               CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_4p8(i),dec_4p8(i),dist)
               IF (dist < min_dist_4p8) THEN 
                  ialphar = ialphar +1
                  alphar = 1. - log10(flux_4p8(i)/flux_radio(k))/log10(4.8/1.4)
                  write(*,'(6x,''Radio alpha: '',f5.2,",",2x,f6.3," arcsec away")') alphar,dist*3600
                  write(12,'("4.800E+09",2x,es9.3)') flux_4p8(i)
                  aalphar = aalphar + alphar
               ENDIF
            ENDDO
            IF (ialphar .ne. 0) then
               alphar = aalphar/float(ialphar)
            else
               write(*,'(" No 5 GHz detection within",f5.0,2x,"arcsec")') min_dist_4p8*3600
            endif
            iofound = 0
            iuvfound = 0
            iirfound = 0
            write(*,*) '.....................IR.............................'
            do i=1,iir
               call DIST_SKY(ra_radio(k),dec_radio(k),ra_ir(i),dec_ir(i),dist)
               IF ( dist < min_dist_ir ) THEN
                  iirfound=iirfound+1
                  If (iirfound > 100) Stop 'Too many IR candidate'
                  airx = 1.-log10(flux_ir(i,2)/flux_x)/log10(frequency_ir(i,2)/2.418e17) !k band, w2 band
                  arir = 1.-log10(flux_radio(k)/flux_ir(i,2))/log10(1.4e9/frequency_ir(i,2))
                  write(*,'(a,"IR-X-ray slope: ",f6.3,",",2x,"radio-IR slope: ",f6.3,",",2x,f6.3," arcsec away")')
     &               ir_type(i),airx,arir,dist*3600
                  if (ir_type(i) == 'WISE') THEN
                     write(*,*) ir_type(i)
                     do s=1,4
                        write(12,'(es9.3,2x,es9.3)') frequency_ir(i,s),flux_ir(i,s)
                     enddo
                  else if (ir_type(i) == '2MASS') THEN
                     do s=2,4
                        write(12,'(es9.3,2x,es9.3)') frequency_ir(i,s),flux_ir(i,s)
                     enddo
                  endif
               endif
            enddo
            IF (iirfound == 0) then
            write(*,'('' NO IR object within '',f5.0,'' arcsec'')') min_dist_ir*3600.
            endif
            write(*,*) '.......................Optical........................'
            DO i=1,iusno
               CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_usno(i),dec_usno(i),dist)
               IF ( dist < min_dist2opt ) THEN 
                    iofound = iofound +1
                    IF (iofound > 100 ) Stop 'Too many optical candidates'
                    aox = 1.-log10(flux_usno(i,1)/flux_x)/log10(frequency_usno(i,1)/2.418e17)
c calculate the optical slope
                    alpho= 1.-log10(flux_usno(i,2)/flux_usno(i,1))/
     &                    log10(frequency_usno(i,2)/frequency_usno(i,1))
                    write(*,'(a,"Optical slope: ",f6.3,",",2x,f6.3," arcsec away")')
     &                 opt_type(i),alpho,dist*3600
                    IF (flux_x > 0. ) THEN 
                       arx = 1.-log10(flux_x/flux_radio(k))/log10(2.41e17/1.4e9)
                       aro = 1.-log10(flux_radio(k)/flux_usno(i,1))/log10(1.4e9/frequency_usno(i,1))
                    ELSE
                       arx = 0.
                       aro = 0.
                    ENDIF
                    do s=1,2
                       write(12,'(es9.3,2x,es9.3)') frequency_usno(i,s),flux_usno(i,s)
                    enddo
                    ra_opt_cand(iofound)   = ra_usno(i)
                    dec_opt_cand(iofound)  = dec_usno(i)
                    aro_opt_cand(iofound)  = aro
                    aox_opt_cand(iofound)  = aox
                    dist_opt_cand(iofound) = dist*3600.
                    mag_opt_cand(iofound)  = mag_usno(i,1)
                    opt_type_cand(iofound) = opt_type(i)
               ENDIF
            ENDDO
            IF (iofound == 0) THEN
            write(*,'('' NO optical object within '', f5.0,
     &                   '' arcsec'')') min_dist2opt*3600.
            ELSE
            CALL indexx (iofound,dist_opt_cand,iofound_index)
            DO i = 1, iofound
            l = iofound_index(i)
            write(*,'(a,''Optical object with mag = '',f5.2,'' found at'',f4.1,
     &                      '' arcsec, aro ='',
     &                      f5.2,'' aox ='',f5.2,'' <arx> = '',f5.2,'' a100x ='',f6.2,
     &                      '' <alpha_r> ='',f6.2,'' RA, DEC ='',f12.7,1x,'','',1x,f12.7)')
     &                      opt_type_cand(l),mag_opt_cand(l),dist_opt_cand(l),aro_opt_cand(l),
     &                      aox_opt_cand(l),arx,a100x,alphar,ra_opt_cand(l),dec_opt_cand(l)
            ENDDO
            ENDIF
            write(*,*) '.......................UV...........................'
            Do i=1,iuv
               CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_uv(i),dec_uv(i),dist)
               IF ( dist < min_dist_uv ) THEN
                  iuvfound = iuvfound +1
                  If (iuvfound > 100) Stop 'Too many UV candidates'
                  auvx = 1.-log10(flux_uv(i,4)/flux_x)/log10(frequency_uv(i,4)/2.418e17) !w1 Fuv
                  aruv = 1.-log10(flux_radio(k)/flux_uv(i,4))/log10(1.4e9/frequency_uv(i,4))
                  write(*,'(a,"UV-X-ray slope: ",f6.3,",",2x,"radio-UV slope: ",f6.3,",",2x,f6.3," arcsec away")')
     &              uv_type(i),auvx,aruv,dist*3600
                  if (uv_type(i) == 'UVOT') then
                     alphauv = 1.-log10(flux_uv(i,1)/flux_uv(i,6))/log10(frequency_uv(i,1)/frequency_uv(i,6)) !u to w2
                     write(*,*) 'UVOT slope',alphauv
                     do s=1,6
                        write(12,'(es9.3,2x,es9.3)') frequency_uv(i,s),flux_uv(i,s)
                     enddo
                  else
                     alphauv = 0.
                     write(12,'(es9.3,2x,es9.3)') frequency_uv(i,4),flux_uv(i,4)
                     !write(12,'(es9.3,2x,es9.3)') frequency_uv(i,s),flux_uv(i,s)
                  endif
               endif
            enddo
            IF (iuvfound == 0) then
            write(*,'('' NO UV object within '',f5.0,'' arcsec'')') min_dist_uv*3600.
            endif
            write(*,*) '........................XRT........................'
            ixrtfound=0
            do i=1,ixrt
               call DIST_SKY(ra_radio(k),dec_radio(k),ra_xrt(i),dec_xrt(i),dist)
               if (dist < min_dist_swift) then
                  ixrtfound=ixrtfound+1
                  if (ixrtfound > 100) stop 'Too many XRT candidate'
                  aswift = 1.-log10(flux_xrt(i,2)/flux_xrt(i,4))/log10(1.325e17/1.091e18)
                  write(*,*) 'xrt slope',aswift
                  do s=2,4
                     write(12,'(es9.3,2x,es9.3)') frequency_xrt(s),flux_xrt(i,s)
                  enddo
               endif
            enddo
            IF (ixrtfound == 0) write(*,'('' NO XRT object within '',f5.0,'' arcsec'')') min_dist_swift*3600.
            igamfound=0
            write(*,*) '.......................Gamma-ray..................'
            do i=1,igam
               call Dist_sky(ra_radio(k),dec_radio(k),ra_gam(i),dec_gam(i),dist)
               if (dist < min_dist_gam) then
                  igamfound=igamfound+1
                  if (igamfound > 20) stop 'Too many Gamma-ray candidate'
                  write(*,'(a,"photon index: ",f5.3,",",2x,f6.3," arcmin away")') gam_type(i),slope_gam(i,1),dist*60
                  if (gam_type(i) == '2FHL' ) then
                     do s=2,4
                        write(12,'(es9.3,2x,es9.3)') frequency_gam(i,s),flux_gam(i,s)
                     enddo
                  else if (gam_type(i) == '3FGL') then
                     do s=2,6
                        write(12,'(es9.3,2x,es9.3)') frequency_gam(i,s),flux_gam(i,s)
                     enddo
                  endif
               endif
            enddo
            IF (igamfound == 0) write(*,'('' NO Gamma-ray detection within '',f5.0,'' arcmin'')') min_dist_gam*60.
            CALL graphic_code (flux_x,flux_radio(k)/const(k),type_average,code)
            write(lu_output,*) ra_radio(k),dec_radio(k),code
            write(*,*) '................Cataloged sources.................'
            DO i=1,iother
               CALL DIST_SKY(ra_radio(k),dec_radio(k),ra_other(i),
     &                       dec_other(i),dist)
               IF (dist < min_dist_other) THEN 
                  write(*,'(2x,a,1x,a)') name_other(i)(1:len(name_other(i))),
     &                           name2_other(i)(1:len(name2_other(i)))
                  IF (name_other(i)(4:7) == 'WHSP') THEN 
                    type_average = -1
                  ELSE IF (name_other(i)(1:5) == 'BZCAT') THEN
                    type_average = -3
                  ENDIF
                  ra_other(i) = -ra_other(i)
               ENDIF
               IF (dist < min_dist_cluster) THEN 
                  IF ( (name_other(i)(1:5) == 'ABELL') .OR. 
     &                      (name_other(i)(1:4) == 'PSZ2') .OR.
     &                      (name_other(i)(1:4) == 'MCXC') .OR.
     &                      (name_other(i)(1:5) == 'SWXCS') .OR.
     &                      (name_other(i)(1:5) == 'ZWCLU') .OR.
     &                      (name_other(i)(1:7) == 'SDSSWHL') ) THEN
                    type_average = -2
                  ENDIF
               ENDIF
               IF (type_average < 0) THEN
                  CALL graphic_code (flux_x,flux_radio(k),type_average,code)
                  write(lu_output,*) ra_radio(k),dec_radio(k),code
               ENDIF
               type_average=0
            ENDDO
            t(ifound)=k
            write(*,*) '        '
         ENDIF
  98     continue
      ENDDO
      if (ifound .ne. sfound+rfound ) stop 'Warning, might have wrong matched number'
      DO i=1,iother
         IF ( ( (name_other(i)(1:5) == 'BZCAT') .OR. (name_other(i)(4:7) == 'WHSP') ) .AND.
     &        (ra_other(i) > 0.) ) THEN
            type_average = -3
            IF (name_other(i)(4:7) == 'WHSP') type_average = -1  
            CALL DIST_SKY(ra_other(i),dec_other(i),ra_center,dec_center,dist)
            dist = dist*60
            print *, achar(27),'[35;1m Known blazar with no radio/X-ray match: ',achar(27),'[0m',
     &                          name2_other(i)(1:lenact(name2_other(i))),' found at a distance of ',
     &                          real(dist),' arcmin '
            CALL graphic_code (1.,1.,-3,code)
            write(lu_output,*) ra_other(i),dec_other(i),code
         ENDIF
      ENDDO
      WRITE (*,*) '      '
      WRITE(*,'(''gnomo_plot_types '',a,'',test.ps/vcps, '',f12.5,2x,f11.4,'' 0. 90.'',
     &           f12.5,2x,f11.4,'' 20.0 10.0 0.3 0. 0.4 0.3 2.1 0. 0.'')') 
     &           output_file(1:lenact(output_file)),ra_center,dec_center,
     &           ra_center,dec_center
      IF (ifound  ==  0) print *,achar(27),'[38;1m Sorry, no radio/X-ray matches were found ',achar(27),'[0m'
      close(lu_output)
      close(12)
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
      IF (radio_type == 1) THEN
         radio_survey='NVSS'
      ELSE IF (radio_type == 2) THEN
         radio_survey='FIRST'
      ELSE IF (radio_type == 3) THEN
         radio_survey='SUMSS'
      ELSE
         radio_survey='UNKNOWN'
      ENDIF
      IF (xray_type == 1) THEN
         xmission='XMMSLEW'
      ELSE IF (xray_type == 2) THEN
         xmission='TWOXMM'
      ELSE IF (xray_type == 3) THEN
         xmission='RASS'
      ELSE IF (xray_type == 4) THEN
         xmission='RXS'
      ELSE IF (xray_type == 5) THEN
         xmission='WGA'
      ELSE IF (xray_type == 6) THEN
         xmission='XRT'
      ELSE IF (xray_type == 7) THEN
         xmission='XRT'
      ELSE IF (xray_type == 8) THEN
         xmission='IPC'
      ELSE IF (xray_type == 9) THEN
         xmission='BMW'
      ELSE
         xmission='UNKNOWN'
      ENDIF
      source_type = 0
c      print *,' xmission, flux-x ',xmission,flux
      CALL DIST_SKY(ra,dec,ra_center,dec_center,dist)
      dist=dist*60.
      IF (flux > 0.) THEN 
         arx = 1.-log10(ratio)/log10(2.41e17/1.4e9)
         IF ( (abs(arx) > 0.42).AND.(abs(arx).LE.0.76) ) THEN 
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
     &      f8.1,'' X-ray/radio flux-ratio'',f8.0,'' arx '',f4.2, 
     &      '' Log(nu_p) '',f4.1,''+/-~1 '',a, 
     &      '' Dist. '',f5.1,'' arcmin'')')
     &      xmission(1:len_trim(xmission)),radio_survey(1:len_trim(radio_survey)),
     &      rah,ram,rasec,sign(1:1),abs(id),abs(dm),abs(decsec),
     &      flux_radio/const,ratio,arx,lognupeak,type(1:len_trim(type)),dist
         ELSE  
          IF ((abs(arx) > 0.76).AND.(abs(arx) < 0.95)) THEN
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
     &       f8.1,'' flux-ratio'',f8.0,'' arx '',f4.2,13x,a,
     &       '' Dist. '',f5.1,'' arcmin'')')
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
     &       f7.1,a,'' Dist. '',f5.1,'' arcmin'')')
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
      IF (source_type < 0 ) THEN
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
      IF (code  <  10000 ) THEN
        print *,' source type radio_component, x_ray_component ', source_type,radio_component,x_ray_component
        stop
      ENDIF
      source_type = temp
      RETURN  
      END
c
c
c

        SUBROUTINE mag2flux (nh,m_band,filter,flux,frequency)
c
c  converts u,v,i,h,b,r,j,k magnitudes into monochromatic fluxes
c  in units of erg/cm2/s for nufnu vs nu plots 
c
        IMPLICIT none
        REAL*4 nh, flux , av , m_band, a_band, Rv 
        REAL*4 c, lambda, const, frequency, a
        REAL*8 x,aa,bb,c1,c2,dx,px,ebv
        CHARACTER*3 filter
c        print *,' nh, m_band, filter ', nh, m_band,filter
c        call upc(filter)
        IF ( (filter(1:1).NE.'U  ') .AND. (filter(1:1).NE.'B  ') .AND.
     &       (filter(1:1).NE.'V  ') .AND. (filter(1:1).NE.'R  ') .AND.
     &       (filter(1:1).NE.'I  ') .AND. (filter(1:1).NE.'J  ') .AND.
     &       (filter(1:1).NE.'H  ') .AND. (filter(1:1).NE.'K  ') .AND.
     &       (filter(1:1).NE.'u  ') .AND. (filter(1:1).NE.'g  ') .AND.
     &       (filter(1:1).NE.'r  ') .AND. (filter(1:1).NE.'i  ') .AND.
     &       (filter(1:1).NE.'z  ') .AND.
     &       (filter(1:3).NE.'fuv') .and. (filter(1:3).NE.'nuv') .AND.
     &       (filter(1:3).NE.'su ') .and. (filter(1:3).NE.'sv ') .AND.
     &       (filter(1:3).NE.'sb ') .and. (filter(1:3).NE.'sw1') .AND.
     &       (filter(1:3).NE.'sm2') .and. (filter(1:3).NE.'sw2') .AND.
     &       (filter(1:3).NE.'ww1') .and. (filter(1:3).NE.'ww2') .AND.
     &       (filter(1:3).NE.'ww3') .and. (filter(1:3).NE.'ww4')) THEN
         write (*,*) ' mag2flux: Filter not supported  '
         stop
        ENDIF
c extintion law taken from Cardelli et al. 1989 ApJ 345, 245
c in UV apply the UV relation from Fitzpatrick 1999
        Rv=3.1
        av = Rv*(-0.055+nh*1.987e-22) !!dust map from BH1978, assumed constant gas-to-dust ratio
        ebv=av/Rv
        if (av < 0.) av=0.
        if (filter(1:1) == 'U') then
           lambda=3550 
           const=log10(1810.)-23.
        else if (filter(1:1) == 'B') then
           lambda=4400
           const=log10(4260.)-23.
        else if (filter(1:1) == 'V') then
           lambda=5500
           const=log10(3640.)-23.
        else if (filter(1:1) == 'R') then 
           lambda=6400 
           const=log10(3080.)-23.
        else if (filter(1:1) == 'I') then 
           lambda=7900 
           const=log10(2550.)-23.
        else if (filter(1:1) == 'J') then !2MASS
           lambda=12350.
           const=log10(1594.)-23.
        else if (filter(1:1) == 'H') then 
           lambda=16620.
           const=log10(1024.)-23.
        else if (filter(1:1) == 'K') then 
           lambda=21590.
           const=log10(666.7)-23.
        else if (filter(1:1) == 'u') then !effective wavelength from SDSS, Doi et al. 2010 ApJ 139, 1628
           lambda=3568.
           const=log10(3631.)-23. !3631 is the 0 mag flux of AB mag system
        else if (filter(1:1) == 'g') then
           m_band=m_band-0.04 !calibrate of the SDSS u band to AB mag
           lambda=4653.
           const=log10(3631.)-23.
        else if (filter(1:1) == 'r') then !effective wavelength from SDSS, Doi et al. 2010 ApJ 139, 1628
           lambda=6148.
           const=log10(3631.)-23. !3631 is the 0 mag flux of AB mag system
        else if (filter(1:1) == 'i') then
           m_band=m_band-0.04 !calibrate of the SDSS u band to AB mag
           lambda=7468.
           const=log10(3631.)-23.
        else if (filter(1:1) == 'z') then
           m_band=m_band-0.04 !calibrate of the SDSS u band to AB mag
           lambda=8863.
           const=log10(3631.)-23.
        else if (filter(1:3) == 'fuv') then
           lambda=1528.
           const=log10(3631.)-23.
        else if (filter(1:3) == 'nuv') then
           lambda=2271.
           const=log10(3631.)-23.
        else if (filter(1:3) == 'su ') then
           lambda=3501. !from Poole et al. (2008) effective wavelength
           const=log10(3631.)-23
        else if (filter(1:3) == 'sb ') then
           lambda=4329. !from Poole et al. (2008) effective wavelength
           const=log10(3631.)-23
        else if (filter(1:3) == 'sv ') then
           lambda=5402. !from Poole et al. (2008) effective wavelength
           const=log10(3631.)-23
        else if (filter(1:3) == 'sw1') then
           lambda=2634. !from Poole et al. (2008) effective wavelength
           const=log10(3631.)-23
        else if (filter(1:3) == 'sm2') then
           lambda=2231. !from Poole et al. (2008) effective wavelength
           const=log10(3631.)-23
        else if (filter(1:3) == 'sw2') then
           lambda=2030. !from Poole et al. (2008) effective wavelength
           const=log10(3631.)-23
        else if (filter(1:3) == 'ww1') then
           lambda=34000
           const=log10(309.540)-23
        else if (filter(1:3) == 'ww2') then
           lambda=46000
           const=log10(171.787)-23
        else if (filter(1:3) == 'ww3') then
           lambda=120000
           const=log10(31.674)-23
        else if (filter(1:3) == 'ww4') then
           lambda=220000
           const=log10(8.363)-23
        endif
c lambda from Amstrongs to microns
        x=10000./lambda
        if ((x .le. 1.1) .and. (x .ge. 0.3)) then
        a_band=(0.574*(x**1.61)-0.527*(x**1.61)/Rv)*av ! the a_lambda
        else if ((x .le. 3.3) .and. (x .ge. 1.1)) then
        x=x-1.82
        aa=1+(0.17699*x)-(0.50447*x**2)-(0.02427*x**3)+(0.73085*x**4)
     &  +(0.01979*x**5)-(0.77530*x**6)+(0.32999*x**7)
        bb=1.41338*x+(2.28305*x**2)+(1.07233*x**3)-(5.38434*x**4)
     &  -(0.662251*x**5)+(5.30260*x**6)-(2.09002*x**7)
        a_band=(aa+(bb/Rv))*av
        else if ((x .le. 10.) .and. (x .ge. 3.3)) then
        c2=-0.824+4.717/Rv
        c1=2.03-3.007*c2
        dx=(x*x)/((x**2-4.596**2)+(x*0.99)**2)
        px=0.5392*(x-5.9)**2+0.05644*(x-5.9)**2
        if (x .le. 5.9) px=0.
        a_band=(c1+c2*x+3.23*dx+0.41*px)*ebv+av
        endif
        c=3.e10
        a=1.0
        frequency=c/(lambda*1.e-8)
        flux = 10.**(-0.4*(m_band -a_band)+const)*frequency
        RETURN 
        END

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
