      PROGRAM find_candidates2
c
c This program reads a query results file from the ASDC error circle explorer tool
c and plots potential blazars of different types based on radio/X-ray ratios 
c
c The output file should be used as input to the program gnomo_plot_types 
c that generates a postscript file with a plot in gnomonic coordinates. 
c In this plot each radio-X/ray match appears as a filled and an open circle with size that ispart
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
      INTEGER*4 ier, lu_in, lu_output,in,im,s,is,ie, i, j,iarr,iconfig
      INTEGER*4 sfound,lenact,type_average,ns,iradio,icat,ix,ir
      INTEGER*4 i4p8,ifar,ilowr,ixrtsp,ipccs100,iusno,iousxb,iswort
      INTEGER*4 iir,iuv,ixray,igam,iiswort,i4fgl,idebl,iref
      INTEGER*4 ialma,ir100found,iverit,ivhe,imagic,ibigb
      INTEGER*4 iofound, length,ialphar,ifarfound,iflcuvfound
      integer*4 isource,filen,sourceu,sourcel,ikuehr,ioptlc
      integer*4 ii1,ii2,ii3,ii4,ii5,iircheck,iirlc,iflcuv,igamfound
      integer*4 r,drop,ixxfound,ilowrfound,iuvfound,iirfound
      integer*4 year,month,date,hour,minute,second,arrsize(30)
      REAL*8 ra, dec,ra_center, dec_center,radius
      REAL*8 dist,min_dist_gam,mjdtest,irdistval
      REAL*4 min_dist2opt,min_dist_at,min_dist_4p8,min_dist_uv,min_dist_ir
      REAL*4 min_dist,min_dist_other,min_dist_pccs100,min_dist_far
      REAL*4 radian,aox,a100x,flux_x,nh,aro,arx,alpho,flux_r
      REAL*4 alphar,sigma,flux2nufnu_pccs100,farirx,code,flux2nufnu_4p8
      REAL*4 auvx,aruv,airx,arir,alphauv,fdens,nudens,aalphar,pccconv
      real*4 posxerr,posyerr,posang,major,minor,intensity,like
      real*4 mjdavg,redshift,reduction_factor,zzfermi,gammatev,rms_lowr
      real*4 engmax,engmin,fluxind,Ufluxind,Lfluxind,zzinput
      CHARACTER*1 sign
      character*6 aim
      CHARACTER*30 irlcid(2)
      character*200 input_file,output_file,input_file2,input_file3,output_file2
      character*200 webprograms,refsfile,array_size
      CHARACTER*15 catalog
      CHARACTER*800 string,repflux
      LOGICAL ok,found,debl

      integer*4,dimension(:),allocatable :: typer
      real*8,dimension(:),allocatable :: ra_source,dec_source
      real*8,dimension(:),allocatable :: ra_cat,dec_cat
      real*4,dimension(:),allocatable :: type_cat
      character*15,dimension(:),allocatable :: name_cat
      character*200,dimension(:),allocatable :: refs

      integer*4,dimension(:),allocatable :: npt
      real*4,dimension(:),allocatable :: zsource
      integer*4,dimension(:,:),allocatable :: spec_type,rrxx_ref,recordmjd
      real*4,dimension(:,:),allocatable :: mjdst_rrxx,mjded_rrxx,frequency,flux,uflux,lflux,epos
      real*8,dimension(:,:),allocatable :: ra_rrxx,dec_rrxx
      character*15,dimension(:,:),allocatable :: rrxx_type
      character*4,dimension(:,:),allocatable :: rrxx_flag

      integer*4,dimension(:),allocatable :: filen_v,vhe_ref
      real*8,dimension(:),allocatable :: ra_vhe,dec_vhe
      real*4,dimension(:),allocatable :: frequency_vhe,flux_vhe,FluxU_vhe,FluxL_vhe
      real*4,dimension(:),allocatable :: poserr_vhe,Ferr_vhe,mjdstart,mjdend
      character*4,dimension(:),allocatable :: vhe_flag
      character*15,dimension(:),allocatable :: vhe_type

      integer*4,dimension(:),allocatable :: veritind,magicind,eblnn
      real*4,dimension(:,:),allocatable :: flux_debl,FluxU_debl,FluxL_debl,frequency_debl
      character*4,dimension(:,:),allocatable :: debl_flag

      integer*4,dimension(:),allocatable :: filen_g,gampart,gam_ref,bigbind
      real*8,dimension(:),allocatable :: ra_gam,dec_gam
      real*4,dimension(:),allocatable :: gamlike,poserr_gam
      real*4,dimension(:,:),allocatable :: frequency_gam,slope_gam,specerr_gam
      real*4,dimension(:,:),allocatable :: flux_gam,FluxL_gam,FluxU_gam,Ferr_gam
      character*4,dimension(:,:),allocatable :: gam_flag
      character*30,dimension(:),allocatable :: namegam
      character*15,dimension(:),allocatable :: name_g,gam_type

      integer*4,dimension(:),allocatable :: xrtspind,xraypart,xray_ref,filen_x
      real*8,dimension(:),allocatable :: ra_xray,dec_xray
      real*4,dimension(:),allocatable :: slope_xray,poserr_xray,mjdst_xrt,mjded_xrt
      real*4,dimension(:,:),allocatable :: frequency_xray,flux_xray,FluxU_xray,FluxL_xray,Ferr_xray
      character*4,dimension(:,:),allocatable :: xray_flag
      character*15,dimension(:),allocatable :: name_x,xray_type

      integer*4,dimension(:),allocatable :: filen_u
      real*8,dimension(:),allocatable :: ra_uv,dec_uv
      real*4,dimension(:),allocatable :: poserr_uv
      real*4,dimension(:,:),allocatable :: uvmag,flux_uv,frequency_uv,FluxU_uv,FluxL_uv,uvmagerr
      character*15,dimension(:),allocatable :: uv_type,name_u

      integer*4,dimension(:),allocatable :: uv_ref
      real*8,dimension(:),allocatable :: ra_uvcand,dec_uvcand
      real*4,dimension(:),allocatable :: uvdist,epos_uvcand
      real*4,dimension(:,:),allocatable :: flux_uvcand,uvmag_cand,freq_uvcand,uflux_uvcand,lflux_uvcand
      character*4,dimension(:,:),allocatable :: uv_flag
      character*15,dimension(:),allocatable :: uvcand_type

      integer*4,dimension(:),allocatable :: flcuvpart,filen_a
      real*8,dimension(:),allocatable :: ra_flcuv,dec_flcuv
      real*4,dimension(:),allocatable :: mjdst_flcuv,mjded_flcuv,ts,duration,slope_flcuv,poserr_flcuv
      real*4,dimension(:,:),allocatable :: frequency_flcuv,flux_flcuv,FluxU_flcuv,FluxL_flcuv,Ferr_flcuv
      character*4,dimension(:,:),allocatable :: flcuv_flag
      character*15,dimension(:),allocatable :: name_a,flcuv_type

      integer*4,dimension(:),allocatable :: optpart,filen_o,indoptlc,optlc_ref
      real*8,dimension(:),allocatable :: ra_usno,dec_usno
      real*4,dimension(:),allocatable :: poserr_usno,mjdst_optlc
      real*4,dimension(:,:),allocatable :: flux_usno,frequency_usno,FluxU_usno,FluxL_usno,usnomagerr,usnomag
      character*4,dimension(:),allocatable :: optlc_flag
      character*15,dimension(:),allocatable :: opt_type,name_o

      integer*4,dimension(:),allocatable :: filen_i,ir_ref,indirlc
      real*8,dimension(:),allocatable :: ra_ir,dec_ir,ra_ircand,dec_ircand
      real*4,dimension(:),allocatable :: irdist,poserr_ir,epos_ircand
      real*4,dimension(:),allocatable :: mjdst_irlc,mjded_irlc,mjdst_irlccand,mjded_irlccand
      real*4,dimension(:,:),allocatable :: flux_ir,irmag,frequency_ir,FluxU_ir,FluxL_ir,irmagerr
      real*4,dimension(:,:),allocatable :: flux_ircand,irmag_cand,freq_ircand,uflux_ircand,lflux_ircand
      character*4,dimension(:,:),allocatable :: flag_ir,ir_flag
      character*15,dimension(:),allocatable :: ir_type,ircand_type,name_i
      character*30,dimension(:),allocatable :: irlc_name

      integer*4,dimension(:),allocatable :: filen_f,far_ref,farpart
      real*8,dimension(:),allocatable :: ra_far,dec_far
      real*4,dimension(:),allocatable :: farlike,poserr_far
      real*4,dimension(:,:),allocatable :: frequency_far,flux_far,Ferr_far,FluxU_far,FluxL_far
      character*4,dimension(:,:),allocatable :: far_flag
      character*15,dimension(:),allocatable :: name_f,far_type

      integer*4,dimension(:),allocatable :: filen_p,pccspart,pccs100_ref,almaind
      real*8,dimension(:),allocatable :: ra_pccs100,dec_pccs100
      real*4,dimension(:),allocatable :: mjdst_alma,mjded_alma,poserr_pccs100,pccslike
      real*4,dimension(:,:),allocatable :: flux_pccs100,frequency_pccs100,flux2_pccs100,snr_pccs100
      real*4,dimension(:,:),allocatable :: Ferr_pccs100,FluxU_pccs100,FluxL_pccs100,Ferr2_pccs100
      character*4,dimension(:,:),allocatable :: pccs100_flag
      character*15,dimension(:),allocatable :: name_p,pccs100_type
      character*30,dimension(:),allocatable :: date_alma

      integer*4,dimension(:),allocatable :: filen_r,f4p8part,f4p8_ref,kuehrind
      real*8,dimension(:),allocatable :: ra_4p8,dec_4p8
      real*4,dimension(:),allocatable :: f4p8like,poserr_4p8
      real*4,dimension(:,:),allocatable :: flux_4p8,frequency_4p8,Ferr_4p8,FluxU_4p8,FluxL_4p8
      character*1,dimension(:,:),allocatable :: flag_4p8
      character*4,dimension(:,:),allocatable :: f4p8_flag
      character*15,dimension(:),allocatable :: f4p8_type,name_r

      integer*4,dimension(:),allocatable :: filen_l,lowr_ref
      real*8,dimension(:),allocatable :: ra_lowr,dec_lowr
      real*4,dimension(:),allocatable :: frequency_lowr,flux_lowr,Ferr_lowr,FluxU_lowr,FluxL_lowr,poserr_lowr
      character*4,dimension(:),allocatable :: lowr_flag
      character*15,dimension(:),allocatable :: name_l,lowr_type
      
      integer*4,dimension(:),allocatable :: opt_ref
      real*8,dimension(:),allocatable :: ra_lowrcand,ra_usnocand,ra_gamslp
      real*8,dimension(:),allocatable :: dec_lowrcand,dec_usnocand,dec_gamslp
      real*4,dimension(:),allocatable :: optdist,epos_usnocand,lowrdist,gamslp,poserr_gamslp
      real*4,dimension(:),allocatable :: freq_lowrcand,flux_lowrcand,uflux_lowrcand,lflux_lowrcand,epos_lowrcand
      real*4,dimension(:,:),allocatable :: freq_usnocand,flux_usnocand,uflux_usnocand,lflux_usnocand,usnomag_cand
      character*4,dimension(:,:),allocatable :: opt_flag
      character*15,dimension(:),allocatable :: optcand_type,lowrcand_type

      common webprograms
      ok = .TRUE.
      found = .FALSE.

      arrsize(1:30)=0
      iconfig=0
      iarr=0
      isource=0
      icat=0
      iswort=0
      iousxb=0
      sfound = 0
      iradio=0
      ilowr=0
      i4p8=0
      iusno=0
      ipccs100=0
      iir=0
      iuv=0
      igam=0
      ixray=0
      ixrtsp=0
      ikuehr=0
      ibigb=0
      imagic=0
      iverit=0
      ivhe=0
      ialma=0
      idebl=0
      iflcuv=0
      iirlc=0
      i4fgl=0
      ioptlc=0
      radian = 45.0/atan(1.0)
c approximate flux conversions from cts/s to erg/cm2/s at 1 kev (NH=5.e20)
      sign=' '
      flux2nufnu_4p8=4.8E9*1.E-26
      flux2nufnu_pccs100=1.e11*1.E-26
      min_dist_4p8=30./3600.
      min_dist_pccs100=180./3600.
      min_dist_at=10./3600.
      min_dist_far=30./3600.
      min_dist_other=15./3600.
      !min_dist_cluster=60./3600.
c 15 arcsecs
      min_dist2opt = 6./3600.
      min_dist_uv=8./3600.
      min_dist_ir=8./3600.
      min_dist_gam=20./60. !10 arcmin
      mjdavg=55000.
      irlcid(1)='First'
      irlcid(2)='First'

      CALL rdforn(string,length)
      IF ( length.NE.0 ) THEN
         CALL rmvlbk(string)
c         write(*,*) string,length
         in=index(string(1:length),' ')
         input_file=string(1:in-1)
         im=index(string(in+1:length),' ')+in
c         write(*,*) in,im
         input_file2=string(in+1:im-1)
         in=im
         im=index(string(in+1:length),' ')+in
         input_file3=string(in+1:im-1)
c         in=im
c         im=index(string(in+1:length),' ')+in
c         input_file4=string(in+1:im-1)
         in=im
         im=index(string(in+1:length),' ')+in
         output_file=string(in+1:im-1)
         in=im
         im=index(string(in+1:length),' ')+in
         output_file2=string(in+1:im-1)
         in=im
         im=index(string(in+1:length),' ')+in
         refsfile=string(in+1:im-1)
         in=im
         im=index(string(in+1:length),' ')+in
         webprograms=string(in+1:im-1)
         in=im
         im=index(string(in+1:length),' ')+in
         read(string(in+1:im-1),*) zzinput
         read(string(im+1:length),'(a)') aim
c         write(*,*) zzinput,aim
         if (aim == 'finish') then
            Stop '!!!!Exit the source exploring routine!!!!'
         else if ((aim == 'sed') .or. (aim == 'lcurve')) then
            ns=1
         else
            read(aim,*) ns
         endif
      ELSE 
         WRITE (*,'('' Enter query results file '',$)')
         READ (*,'(a)') input_file
      ENDIF

      lu_in = 10
      lu_output = 11
      in = index(input_file(1:lenact(input_file)),'.')
      IF (in == 0) input_file(lenact(input_file)+1:lenact(input_file)+4) = '.csv'
      array_size=webprograms(1:lenact(webprograms))//'/array_size.cf'
      !write(*,*) input_file
c      INQUIRE (FILE=input_file,EXIST=there)
c      IF (.NOT.there) THEN
c         write (*,'('' file '',a,'' not found. No data found in phase 2.'')')
c     &     input_file(1:lenact(input_file))
c         open(18,file=input_file,status='unknown',iostat=ier)
c         write(18,'("RA= ",f9.5,2x,"Dec= ",f9.5,2x,"Searching radius= ",f6.2)')
c     &      0.,0.,0.
c         write(18,'("nH= ",es9.3,2x,"Error circle/elliptical= ",4(f6.2,2x))') 0.,0.,0.,0.,0.
c         write(18,'("0,nvss,41.21762,11.53044,1.097,33,1.4,26.7,2.7,19.3,,-69.5")')
c         close(18)
c      ENDIF
      open(lu_in,file=input_file,status='old',iostat=ier)
      open(13,file=input_file2,status='old',iostat=ier)
      open(12,file=input_file3,status='old',iostat=ier)
      open(15,file=refsfile,status='old',iostat=ier)
      open(18,file=array_size,status='old',iostat=ier)
c      open(16,file=input_file4,status='old',iostat=ier)
      IF (ier.NE.0) THEN
        write (*,*) ' Error ',ier,' opening file '
      ENDIF

      do while(ok)
      read(18,'(a)',end=700) string
      if (string(1:3) == '---' ) then
         iconfig=iconfig+1
      else
         if (iconfig .eq. 1) then
            iarr=iarr+1
            if (iarr .le. 2) then
               is=index(string(1:len(string)),':')
               read(string(is+1:lenact(string)),*) arrsize(iarr)
            else
               iarr=iarr-1
            endif
         else if (iconfig .eq. 3) then
            iarr=iarr+1
            is=index(string(1:len(string)),':')
            read(string(is+1:lenact(string)),*) arrsize(iarr)
         endif
      endif
      enddo
700   continue
      close(18)
      !write(*,*) arrsize

      allocate(typer(arrsize(3)))
      allocate(ra_source(arrsize(3)),dec_source(arrsize(3)))
      allocate(ra_cat(arrsize(4)),dec_cat(arrsize(4)))
      allocate(type_cat(arrsize(4)),name_cat(arrsize(4)),refs(arrsize(4)))

c read the find_out.txt first
      Do WHILE(ok)
         read(13,*,end=100) ra,dec,code
         if ((code .ge. 0.) .or. (code .lt. -40000.) .or. (code .eq. -9999)) then
            isource=isource+1
            ra_source(isource)=ra
            dec_source(isource)=dec
            typer(isource)=int(code/10000.)
            if (code .eq. 99) typer(isource)=99
         else
            icat=icat+1
            ra_cat(icat)=ra
            dec_cat(icat)=dec
            type_cat(icat)=int(code/10000.)
         ENDIF
      enddo
100   continue
      !write(*,*) icat,isource  !!!!!!the number of source
      close(13)
c      do i=1,icat
c         write(lu_output,*) ra_cat(i),dec_cat(i),type_cat(i)*10000.
c      enddo

      allocate(mjdst_rrxx(arrsize(2),arrsize(1)),mjded_rrxx(arrsize(2),arrsize(1)),frequency(arrsize(2),arrsize(1)),flux(arrsize(2),arrsize(1)),uflux(arrsize(2),arrsize(1)),lflux(arrsize(2),arrsize(1)),epos(arrsize(2),arrsize(1)))
      allocate(spec_type(arrsize(2),arrsize(1)),ra_rrxx(arrsize(2),arrsize(1)),dec_rrxx(arrsize(2),arrsize(1)),rrxx_type(arrsize(2),arrsize(1)))
      allocate(recordmjd(3,arrsize(2)),npt(arrsize(1)),zsource(arrsize(1)))

c read the sed.txt first
      npt(1:arrsize(1))=0
202   continue
      read(12,'(i4,a)',end=200,err=200) sfound,string
c      write(*,*) sfound,string
      do while(ok)
      if (sfound .ge. 1000) then
         read(12,*,end=201,err=201) zsource(1:isource)
c         write(*,*) "redshift",zsource(1:isource)
      else
         npt(sfound)=npt(sfound)+1
         read(12,'(es10.3,a)',end=201,err=201) frequency(npt(sfound),sfound),string
         !write(*,*) frequency(npt(sfound),sfound)
         if ((frequency(npt(sfound),sfound) .lt. 1.e10) .or. (frequency(npt(sfound),sfound) .eq. 2.418e17)) then
            read(string(1:lenact(string)),*) flux(npt(sfound),sfound),uflux(npt(sfound),sfound),lflux(npt(
     &      sfound),sfound),ra_rrxx(npt(sfound),sfound),dec_rrxx(npt(sfound),sfound),epos(npt(sfound),sfound),
     &      mjdst_rrxx(npt(sfound),sfound),mjded_rrxx(npt(sfound),sfound),spec_type(npt(sfound),sfound)
            if ((spec_type(npt(sfound),sfound) .eq. 61) .or. (spec_type(npt(sfound),sfound) .eq. 66)) then
               iousxb=iousxb+1
               recordmjd(1:3,iousxb)=[iousxb,npt(sfound),sfound]
            endif
         else
            read(string(1:lenact(string)),*) flux(npt(sfound),sfound),uflux(npt(sfound),sfound),
     &          lflux(npt(sfound),sfound),spec_type(npt(sfound),sfound)
            if ((spec_type(npt(sfound),sfound) .eq. 11) .or. (spec_type(npt(sfound),sfound) .eq. 16)) then
               iswort=iswort+1
               iiswort=iswort/4
               if (MOD(iswort,4) .ne. 0) iiswort=iiswort+1
               mjdst_rrxx(npt(sfound),sfound)=mjdst_rrxx(recordmjd(2,iiswort),recordmjd(3,iiswort))
               mjded_rrxx(npt(sfound),sfound)=mjded_rrxx(recordmjd(2,iiswort),recordmjd(3,iiswort))
            else
               mjdst_rrxx(npt(sfound),sfound)=55000.
               mjded_rrxx(npt(sfound),sfound)=55000.
            endif
         endif
         if (((spec_type(npt(sfound),sfound) .eq. 59) .or. (spec_type(npt(sfound),sfound) .eq. 55)) .and. (ra_rrxx(npt(sfound),sfound) .lt. 0.)) then
            frequency(npt(sfound),sfound)=(1.602E-19)*(3.e3)/(6.626e-34)
            ra_rrxx(npt(sfound),sfound)=abs(ra_rrxx(npt(sfound),sfound))
         endif
      !write(*,*) sfound,frequency(npt(sfound),sfound),spec_type(npt(sfound),sfound),epos(npt(sfound),sfound)
         if ((spec_type(npt(sfound),sfound) .eq. 60) .and. (ra_rrxx(npt(sfound),sfound) .lt. 0.)) then
            frequency(npt(sfound),sfound)=(1.602E-19)*(5.e3)/(6.626e-34)
            ra_rrxx(npt(sfound),sfound)=abs(ra_rrxx(npt(sfound),sfound))
         endif
         if (frequency(npt(sfound),sfound) .gt. 1.E11 ) then
         if ((spec_type(npt(sfound),sfound) .eq. 51) .or. (spec_type(npt(sfound),sfound) .eq. 1))
     &      rrxx_type(npt(sfound),sfound)='XMMSL2'
         if ((spec_type(npt(sfound),sfound) .eq. 52) .or. (spec_type(npt(sfound),sfound) .eq. 2))
     &      rrxx_type(npt(sfound),sfound)='4XMM-DR11'
         if ((spec_type(npt(sfound),sfound) .eq. 53) .or. (spec_type(npt(sfound),sfound) .eq. 3))
     &      rrxx_type(npt(sfound),sfound)='RASS'
         if ((spec_type(npt(sfound),sfound) .eq. 54) .or. (spec_type(npt(sfound),sfound) .eq. 4))
     &      rrxx_type(npt(sfound),sfound)='WGACAT'
         if ((spec_type(npt(sfound),sfound) .eq. 55) .or. (spec_type(npt(sfound),sfound) .eq. 5))
     &      rrxx_type(npt(sfound),sfound)='2SXPS'
         if ((spec_type(npt(sfound),sfound) .eq. 56) .or. (spec_type(npt(sfound),sfound) .eq. 6))
     &      rrxx_type(npt(sfound),sfound)='IPC2E'
         if ((spec_type(npt(sfound),sfound) .eq. 57) .or. (spec_type(npt(sfound),sfound) .eq. 7))
     &      rrxx_type(npt(sfound),sfound)='BMW'
         if ((spec_type(npt(sfound),sfound) .eq. 58) .or. (spec_type(npt(sfound),sfound) .eq. 8))
     &      rrxx_type(npt(sfound),sfound)='Chandra-CSC2'
         if ((spec_type(npt(sfound),sfound) .eq. 59) .or. (spec_type(npt(sfound),sfound) .eq. 9))
     &      rrxx_type(npt(sfound),sfound)='XRTDEEP'
         if ((spec_type(npt(sfound),sfound) .eq. 60) .or. (spec_type(npt(sfound),sfound) .eq. 10))
     &      rrxx_type(npt(sfound),sfound)='MAXIGSC'
         if ((spec_type(npt(sfound),sfound) .eq. 61) .or. (spec_type(npt(sfound),sfound) .eq. 11))
     &      rrxx_type(npt(sfound),sfound)='1OUSX'
         if ((spec_type(npt(sfound),sfound) .eq. 62) .or. (spec_type(npt(sfound),sfound) .eq. 12))
     &      rrxx_type(npt(sfound),sfound)='IPCSL'
         if ((spec_type(npt(sfound),sfound) .eq. 63) .or. (spec_type(npt(sfound),sfound) .eq. 13))
     &      rrxx_type(npt(sfound),sfound)='MAXISSC'
         if ((spec_type(npt(sfound),sfound) .eq. 64) .or. (spec_type(npt(sfound),sfound) .eq. 14))
     &      rrxx_type(npt(sfound),sfound)='eROSITA-EDR'

         else
            if (spec_type(npt(sfound),sfound) .eq. 2) rrxx_type(npt(sfound),sfound)='FIRST'
            if (spec_type(npt(sfound),sfound) .eq. 3) rrxx_type(npt(sfound),sfound)='NVSS'
            if (spec_type(npt(sfound),sfound) .eq. 4) rrxx_type(npt(sfound),sfound)='SUMSS'
            if (spec_type(npt(sfound),sfound) .eq. 1) rrxx_type(npt(sfound),sfound)='VLASSQL'
         endif
      endif
      enddo
201   continue
      npt(sfound)=npt(sfound)-1
      goto 202
200   continue
      close(12)

c         write(*,*) "redshift",zsource(1:isource)
c read the reference file
      iref=0
      Do WHILE(ok)
         read(15,'(a)',end=400) string
         iref=iref+1
         ie = index(string(1:len(string)),' ')
         read(string(1:ie-1),'(a)') name_cat(iref)
         is = ie
         ie = index(string(1:len(string)),'!')
         read(string(is+1:ie-1),'(a)') refs(iref)
c         write(*,*) name_cat(iref),refs(iref)(1:lenact(refs(iref)))
      enddo
400   continue
      close(15)

      allocate(filen_v(arrsize(5)),vhe_ref(arrsize(5)),eblnn(arrsize(5)+arrsize(7)))
      allocate(ra_vhe(arrsize(5)),dec_vhe(arrsize(5)))
      allocate(frequency_vhe(arrsize(5)),flux_vhe(arrsize(5)),FluxU_vhe(arrsize(5)),FluxL_vhe(arrsize(5)))
      allocate(poserr_vhe(arrsize(5)),Ferr_vhe(arrsize(5)),mjdstart(arrsize(5)),mjdend(arrsize(5)))
      allocate(vhe_flag(arrsize(5)),vhe_type(arrsize(5)))

      allocate(veritind(arrsize(6)),magicind(arrsize(6)))
      allocate(flux_debl(arrsize(6),5),FluxU_debl(arrsize(6),5),FluxL_debl(arrsize(6),5),frequency_debl(arrsize(6),5))
      allocate(debl_flag(arrsize(6),5))

      allocate(filen_g(arrsize(7)),gampart(arrsize(7)),gam_ref(arrsize(7)),bigbind(arrsize(7)))
      allocate(ra_gam(arrsize(7)),dec_gam(arrsize(7)),gamlike(arrsize(7)),poserr_gam(arrsize(7)))
      allocate(frequency_gam(arrsize(7),9),slope_gam(arrsize(7),2),specerr_gam(arrsize(7),2))
      allocate(flux_gam(arrsize(7),9),FluxL_gam(arrsize(7),9),FluxU_gam(arrsize(7),9),Ferr_gam(arrsize(7),9))
      allocate(gam_flag(arrsize(7),9))
      allocate(namegam(arrsize(7)),name_g(arrsize(7)),gam_type(arrsize(7)))

      allocate(xray_ref(arrsize(9)),xray_flag(arrsize(9),5))
      allocate(ra_xray(arrsize(9)),dec_xray(arrsize(9)),xrtspind(arrsize(9)),xraypart(arrsize(9)))
      allocate(slope_xray(arrsize(9)),poserr_xray(arrsize(9)),mjdst_xrt(arrsize(9)),mjded_xrt(arrsize(9)))
      allocate(frequency_xray(arrsize(9),5),flux_xray(arrsize(9),5),FluxU_xray(arrsize(9),5),FluxL_xray(arrsize(9),5),Ferr_xray(arrsize(9),5))
      allocate(name_x(arrsize(9)),xray_type(arrsize(9)),filen_x(arrsize(9)))

      allocate(filen_u(arrsize(10)),uv_type(arrsize(10)),name_u(arrsize(10)))
      allocate(ra_uv(arrsize(10)),dec_uv(arrsize(10)),poserr_uv(arrsize(10)))
      allocate(uvmag(arrsize(10),6),flux_uv(arrsize(10),6),frequency_uv(arrsize(10),6),FluxU_uv(arrsize(10),6),FluxL_uv(arrsize(10),6),uvmagerr(arrsize(10),6))

      allocate(uv_ref(arrsize(11)))
      allocate(ra_uvcand(arrsize(11)),dec_uvcand(arrsize(11)),uvdist(arrsize(11)),epos_uvcand(arrsize(11)))
      allocate(flux_uvcand(arrsize(11),6),uvmag_cand(arrsize(11),6),freq_uvcand(arrsize(11),6),uflux_uvcand(arrsize(11),6),lflux_uvcand(arrsize(11),6))
      allocate(uv_flag(arrsize(11),6),uvcand_type(arrsize(11)))

      allocate(ra_flcuv(arrsize(8)),dec_flcuv(arrsize(8)))
      allocate(mjdst_flcuv(arrsize(8)),mjded_flcuv(arrsize(8)),ts(arrsize(8)),duration(arrsize(8)),slope_flcuv(arrsize(8)),poserr_flcuv(arrsize(8)))
      allocate(frequency_flcuv(arrsize(8),4),flux_flcuv(arrsize(8),4),FluxU_flcuv(arrsize(8),4),FluxL_flcuv(arrsize(8),4),Ferr_flcuv(arrsize(8),4))
      allocate(flcuv_flag(arrsize(8),4))
      allocate(flcuv_type(arrsize(8)),name_a(arrsize(8)),flcuvpart(arrsize(8)),filen_a(arrsize(8)))

      allocate(optpart(arrsize(12)),filen_o(arrsize(12)),indoptlc(arrsize(12)),optlc_ref(arrsize(12)))
      allocate(optlc_flag(arrsize(12)),opt_type(arrsize(12)),name_o(arrsize(12)))
      allocate(ra_usno(arrsize(12)),dec_usno(arrsize(12)),poserr_usno(arrsize(12)),mjdst_optlc(arrsize(12)))
      allocate(flux_usno(arrsize(12),5),frequency_usno(arrsize(12),5),FluxU_usno(arrsize(12),5),FluxL_usno(arrsize(12),5),usnomagerr(arrsize(12),5),usnomag(arrsize(12),5))

      allocate(filen_i(arrsize(13)),ir_ref(arrsize(13)),indirlc(arrsize(13)))
      allocate(ra_ir(arrsize(13)),dec_ir(arrsize(13)),ra_ircand(arrsize(13)),dec_ircand(arrsize(13)))
      allocate(irdist(arrsize(13)),poserr_ir(arrsize(13)),epos_ircand(arrsize(13)))
      allocate(mjdst_irlc(arrsize(13)),mjded_irlc(arrsize(13)),mjdst_irlccand(arrsize(13)),mjded_irlccand(arrsize(13)))
      allocate(flux_ir(arrsize(13),4),irmag(arrsize(13),4),frequency_ir(arrsize(13),4),FluxU_ir(arrsize(13),4),FluxL_ir(arrsize(13),4),irmagerr(arrsize(13),4))
      allocate(flux_ircand(arrsize(13),4),irmag_cand(arrsize(13),4),freq_ircand(arrsize(13),4),uflux_ircand(arrsize(13),4),lflux_ircand(arrsize(13),4))
      allocate(flag_ir(arrsize(13),4),ir_flag(arrsize(13),4))
      allocate(ir_type(arrsize(13)),ircand_type(arrsize(13)),name_i(arrsize(13)),irlc_name(arrsize(13)))

      allocate(filen_f(arrsize(14)),far_ref(arrsize(14)),farpart(arrsize(14)))
      allocate(ra_far(arrsize(14)),dec_far(arrsize(14)),farlike(arrsize(14)),poserr_far(arrsize(14)))
      allocate(frequency_far(arrsize(14),5),flux_far(arrsize(14),5),Ferr_far(arrsize(14),5),FluxU_far(arrsize(14),5),FluxL_far(arrsize(14),5))
      allocate(far_flag(arrsize(14),5),name_f(arrsize(14)),far_type(arrsize(14)))

      allocate(filen_p(arrsize(15)),pccspart(arrsize(15)),pccs100_ref(arrsize(15)),almaind(arrsize(15)))
      allocate(ra_pccs100(arrsize(15)),dec_pccs100(arrsize(15)))
      allocate(mjdst_alma(arrsize(15)),mjded_alma(arrsize(15)),poserr_pccs100(arrsize(15)),pccslike(arrsize(15)))
      allocate(flux_pccs100(arrsize(15),9),frequency_pccs100(arrsize(15),9),flux2_pccs100(arrsize(15),9),snr_pccs100(arrsize(15),9))
      allocate(Ferr_pccs100(arrsize(15),9),FluxU_pccs100(arrsize(15),9),FluxL_pccs100(arrsize(15),9),Ferr2_pccs100(arrsize(15),9))
      allocate(pccs100_flag(arrsize(15),9))
      allocate(name_p(arrsize(15)),pccs100_type(arrsize(15)),date_alma(arrsize(15)))

      allocate(filen_r(arrsize(16)),f4p8part(arrsize(16)),f4p8_ref(arrsize(16)),kuehrind(arrsize(16)))
      allocate(ra_4p8(arrsize(16)),dec_4p8(arrsize(16)),f4p8like(arrsize(16)),poserr_4p8(arrsize(16)))
      allocate(flux_4p8(arrsize(16),3),frequency_4p8(arrsize(16),3),Ferr_4p8(arrsize(16),3),FluxU_4p8(arrsize(16),3),FluxL_4p8(arrsize(16),3))
      allocate(flag_4p8(arrsize(16),4),f4p8_flag(arrsize(16),3))
      allocate(f4p8_type(arrsize(16)),name_r(arrsize(16)))

      allocate(filen_l(arrsize(17)),lowr_ref(arrsize(17)))
      allocate(ra_lowr(arrsize(17)),dec_lowr(arrsize(17)))
      allocate(frequency_lowr(arrsize(17)),flux_lowr(arrsize(17)),Ferr_lowr(arrsize(17)),FluxU_lowr(arrsize(17)),FluxL_lowr(arrsize(17)),poserr_lowr(arrsize(17)))
      allocate(lowr_flag(arrsize(17)),name_l(arrsize(17)),lowr_type(arrsize(17)))

      allocate(opt_ref(arrsize(18)))
      allocate(ra_lowrcand(arrsize(18)),ra_usnocand(arrsize(18)),ra_gamslp(arrsize(18)))
      allocate(dec_lowrcand(arrsize(18)),dec_usnocand(arrsize(18)),dec_gamslp(arrsize(18)))
      allocate(optdist(arrsize(18)),epos_usnocand(arrsize(18)),lowrdist(arrsize(18)),gamslp(arrsize(18)),poserr_gamslp(arrsize(18)))
      allocate(freq_lowrcand(arrsize(18)),flux_lowrcand(arrsize(18)),uflux_lowrcand(arrsize(18)),lflux_lowrcand(arrsize(18)),epos_lowrcand(arrsize(18)))
      allocate(freq_usnocand(arrsize(18),arrsize(18)),flux_usnocand(arrsize(18),arrsize(18)),uflux_usnocand(arrsize(18),5),lflux_usnocand(arrsize(18),5),usnomag_cand(arrsize(18),5))
      allocate(opt_flag(arrsize(18),5),optcand_type(arrsize(18)),lowrcand_type(arrsize(18)) )

c read the data file
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
      if (nh .eq. 0.d0) THEN
         write (*,'('' file '',a,'' not found. No data found in phase 2.'')')
     &     input_file(1:lenact(input_file))
      endif
c      write(*,*) ra_center,dec_center,radius
c      write(*,*) nh
      !nh=5.e20

c        write(*,*) zsource(ns),zzinput,redshift
      if (zsource(ns) .gt. 0.) redshift=zsource(ns)
      if (zzinput .gt. 0.) redshift=zzinput

      DO WHILE(ok)
300    continue
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
         IF ( (catalog(1:3) == 'pmn') .OR. (catalog(1:6) == 'crates') .OR.
     &       (catalog(1:2) == 'at') .OR. (catalog(1:4) == 'gb87') .or.
     &       (catalog(1:3) == 'gb6') .or.(catalog(1:7) == 'north20') .or.
     &       (catalog(1:7) == 'f357det') .or.(catalog(1:5) == 'kuehr'))  THEN
            i4p8=i4p8+1
            IF (i4p8 > arrsize(16)) Stop 'Too many PMN points'
            if (i4p8 .ne. 1) THEN
               do j=1,i4p8-1
                  if ((ra_4p8(j) .eq. ra) .and. (dec_4p8(j) .eq. dec)) THEN
                     if ((name_r(j) == catalog ) .and. (filen_r(j) .ne. filen)) then
c                     write(*,'(4x,"The counterpart",i4,3x,2(f9.5,2x),3(f9.3,2x))')
c     &                 j,ra_4p8(j),dec_4p8(j),flux_4p8(j,1:3)/(frequency_4p8(j,1:3)*1.E-26)
                     is=ie
                     ie=index(string(is+1:len(string)),' ')+is
                     read(string(is+1:ie-1),'(a)') repflux
c                     write(*,'(a,i4,"_",a,2x,2(f9.5,2x),a)') ' Repeated one:',
c     &                   filen,catalog,ra,dec,repflux(1:lenact(repflux))
                     i4p8=i4p8-1
                     goto 300
                     endif
                  endif
               enddo
            endif
            ra_4p8(i4p8)=ra
            dec_4p8(i4p8)=dec
            name_r(i4p8)=catalog
            filen_r(i4p8)=filen
            if  (catalog(1:4) == 'gb87')  THEN
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_4p8(i4p8)
               poserr_4p8(i4p8)=poserr_4p8(i4p8)*2.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_4p8(i4p8,1)
               FluxU_4p8(i4p8,1)=flux_4p8(i4p8,1)+Ferr_4p8(i4p8,1)
               FluxL_4p8(i4p8,1)=flux_4p8(i4p8,1)-Ferr_4p8(i4p8,1)
               flux_4p8(i4p8,1)=flux_4p8(i4p8,1)*flux2nufnu_4p8
               FluxU_4p8(i4p8,1)=FluxU_4p8(i4p8,1)*flux2nufnu_4p8
               FluxL_4p8(i4p8,1)=FluxL_4p8(i4p8,1)*flux2nufnu_4p8
               frequency_4p8(i4p8,1)=4.8e9
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,2)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,3)
               f4p8_type(i4p8)='GB87'
               if (flag_4p8(i4p8,2) == 'W') flux_4p8(i4p8,1)=-flux_4p8(i4p8,1)
            else if  (catalog(1:3) == 'gb6')  THEN
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_4p8(i4p8)
               poserr_4p8(i4p8)=2.*poserr_4p8(i4p8)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_4p8(i4p8,1)
               FluxU_4p8(i4p8,1)=flux_4p8(i4p8,1)+Ferr_4p8(i4p8,1)
               FluxL_4p8(i4p8,1)=flux_4p8(i4p8,1)-Ferr_4p8(i4p8,1)
               flux_4p8(i4p8,1)=flux_4p8(i4p8,1)*flux2nufnu_4p8
               FluxU_4p8(i4p8,1)=FluxU_4p8(i4p8,1)*flux2nufnu_4p8
               FluxL_4p8(i4p8,1)=FluxL_4p8(i4p8,1)*flux2nufnu_4p8
               frequency_4p8(i4p8,1)=4.8e9
               f4p8_type(i4p8)='GB6'
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,2)
               if (flag_4p8(i4p8,2) == 'W') flux_4p8(i4p8,1)=-flux_4p8(i4p8,1)
            else if (catalog(1:5) == 'kuehr') then
               ikuehr=ikuehr+1
               kuehrind(i4p8)=ikuehr
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) frequency_4p8(i4p8,1)
               frequency_4p8(i4p8,1)=frequency_4p8(i4p8,1)*1.e6
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_4p8(i4p8,1)
               FluxU_4p8(i4p8,1)=flux_4p8(i4p8,1)+Ferr_4p8(i4p8,1)
               FluxL_4p8(i4p8,1)=flux_4p8(i4p8,1)-Ferr_4p8(i4p8,1)
               flux_4p8(i4p8,1)=flux_4p8(i4p8,1)*1.e-26*frequency_4p8(i4p8,1)
               FluxU_4p8(i4p8,1)=FluxU_4p8(i4p8,1)*1.e-26*frequency_4p8(i4p8,1)
               FluxL_4p8(i4p8,1)=FluxL_4p8(i4p8,1)*1.e-26*frequency_4p8(i4p8,1)
               poserr_4p8(i4p8)=30.
               f4p8_type(i4p8)='KUEHR'
            else if (catalog(1:3) == 'pmn') THEN
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_4p8(i4p8,1)
               posxerr=(1300./flux_4p8(i4p8,1))**2+(6**2)
               posyerr=(1100./flux_4p8(i4p8,1))**2+(4**2)
               if (Ferr_4p8(i4p8,1) .eq. 99.) then
                  if ((dec_4p8(i4p8) .ge. -87.5) .and. (dec_4p8(i4p8) .lt. -37.)) then
                     Ferr_4p8(i4p8,1)=sqrt((12.3+0.085*real(dec_4p8(i4p8)))**2+(0.052*flux_4p8(i4p8,1))**2)
                  else if ((dec_4p8(i4p8) .ge. 37.) .and. (dec_4p8(i4p8) .lt. -29.)) then
                     Ferr_4p8(i4p8,1)=sqrt((11.**2)+(0.052*flux_4p8(i4p8,1))**2)
                  else if ((dec_4p8(i4p8) .ge. -29.) .and. (dec_4p8(i4p8) .lt. -9.5)) then
                     Ferr_4p8(i4p8,1)=sqrt((10.5**2)+(0.052*flux_4p8(i4p8,1))**2)
                  else if ((dec_4p8(i4p8) .ge. -9.5) .and. (dec_4p8(i4p8) .lt. -10.)) then
                     Ferr_4p8(i4p8,1)=sqrt((9.1**2)+(0.052*flux_4p8(i4p8,1))**2)
                  endif
               endif
               FluxU_4p8(i4p8,1)=flux_4p8(i4p8,1)+Ferr_4p8(i4p8,1)
               FluxL_4p8(i4p8,1)=flux_4p8(i4p8,1)-Ferr_4p8(i4p8,1)
               flux_4p8(i4p8,1)=flux_4p8(i4p8,1)*flux2nufnu_4p8
               FluxU_4p8(i4p8,1)=FluxU_4p8(i4p8,1)*flux2nufnu_4p8
               FluxL_4p8(i4p8,1)=FluxL_4p8(i4p8,1)*flux2nufnu_4p8
               frequency_4p8(i4p8,1)=4.8e9
               poserr_4p8(i4p8)=2.*sqrt(posxerr+posyerr)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,1)
               f4p8_type(i4p8)='PMN'
               if (flag_4p8(i4p8,2) == 'X') flux_4p8(i4p8,1)=-flux_4p8(i4p8,1)
            else if (catalog(1:5) == 'atpmn') THEN
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_4p8(i4p8,2)
               FluxU_4p8(i4p8,2)=flux_4p8(i4p8,2)+Ferr_4p8(i4p8,2)
               FluxL_4p8(i4p8,2)=flux_4p8(i4p8,2)-Ferr_4p8(i4p8,2)
               flux_4p8(i4p8,2)=flux_4p8(i4p8,2)*flux2nufnu_4p8*(8.6/4.8)
               FluxU_4p8(i4p8,2)=FluxU_4p8(i4p8,2)*flux2nufnu_4p8*(8.6/4.8)
               FluxL_4p8(i4p8,2)=FluxL_4p8(i4p8,2)*flux2nufnu_4p8*(8.6/4.8)
               frequency_4p8(i4p8,2)=8.6e9
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_4p8(i4p8,1)
               FluxU_4p8(i4p8,1)=flux_4p8(i4p8,1)+Ferr_4p8(i4p8,1)
               FluxL_4p8(i4p8,1)=flux_4p8(i4p8,1)-Ferr_4p8(i4p8,1)
               flux_4p8(i4p8,1)=flux_4p8(i4p8,1)*flux2nufnu_4p8
               FluxU_4p8(i4p8,1)=FluxU_4p8(i4p8,1)*flux2nufnu_4p8
               FluxL_4p8(i4p8,1)=FluxL_4p8(i4p8,1)*flux2nufnu_4p8
               frequency_4p8(i4p8,1)=4.8e9
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) posxerr
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,2)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) posyerr
               posxerr=2*150.*10**(-posxerr)
               posyerr=2*10.*10**(-posyerr)
               poserr_4p8(i4p8)=sqrt((posxerr*posxerr)+(posyerr*posyerr))
               f4p8_type(i4p8)='ATPMN'
               if (flag_4p8(i4p8,2) == '1') flux_4p8(i4p8,1:2)=-flux_4p8(i4p8,1:2)
            else if (catalog(1:4) == 'at20') THEN
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_4p8(i4p8,3)
               FluxU_4p8(i4p8,3)=flux_4p8(i4p8,3)+Ferr_4p8(i4p8,3)
               FluxL_4p8(i4p8,3)=flux_4p8(i4p8,3)-Ferr_4p8(i4p8,3)
               flux_4p8(i4p8,3)=flux_4p8(i4p8,3)*flux2nufnu_4p8*(20/4.8)
               FluxU_4p8(i4p8,3)=FluxU_4p8(i4p8,3)*flux2nufnu_4p8*(20/4.8)
               FluxL_4p8(i4p8,3)=FluxL_4p8(i4p8,3)*flux2nufnu_4p8*(20/4.8)
               frequency_4p8(i4p8,3)=2.e10
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_4p8(i4p8,2)
               FluxU_4p8(i4p8,2)=flux_4p8(i4p8,2)+Ferr_4p8(i4p8,2)
               FluxL_4p8(i4p8,2)=flux_4p8(i4p8,2)-Ferr_4p8(i4p8,2)
               flux_4p8(i4p8,2)=flux_4p8(i4p8,2)*flux2nufnu_4p8*(8/4.8)
               FluxU_4p8(i4p8,2)=FluxU_4p8(i4p8,2)*flux2nufnu_4p8*(8/4.8)
               FluxL_4p8(i4p8,2)=FluxL_4p8(i4p8,2)*flux2nufnu_4p8*(8/4.8)
               frequency_4p8(i4p8,2)=8.e9
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_4p8(i4p8,1)
               FluxU_4p8(i4p8,1)=flux_4p8(i4p8,1)+Ferr_4p8(i4p8,1)
               FluxL_4p8(i4p8,1)=flux_4p8(i4p8,1)-Ferr_4p8(i4p8,1)
               flux_4p8(i4p8,1)=flux_4p8(i4p8,1)*flux2nufnu_4p8*(5/4.8)
               FluxU_4p8(i4p8,1)=FluxU_4p8(i4p8,1)*flux2nufnu_4p8*(5/4.8)
               FluxL_4p8(i4p8,1)=FluxL_4p8(i4p8,1)*flux2nufnu_4p8*(5/4.8)
               frequency_4p8(i4p8,1)=5.e9
               poserr_4p8(i4p8)=sqrt(1.81)*2. !average 0.9 1
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,2)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,1)
               f4p8_type(i4p8)='AT20G'
               if (flag_4p8(i4p8,2) == 'p') flux_4p8(i4p8,1:3)=-flux_4p8(i4p8,1:3)
            else if (catalog(1:6) == 'crates') THEN
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,2)
               flux_4p8(i4p8,1)=flux_4p8(i4p8,1)*flux2nufnu_4p8
               flux_4p8(i4p8,2)=flux_4p8(i4p8,2)*flux2nufnu_4p8*(8.4/4.8)
               frequency_4p8(i4p8,1)=4.8E9
               frequency_4p8(i4p8,2)=8.4E9
               poserr_4p8(i4p8)=5. !!!!!!!
               FluxU_4p8(i4p8,1)=0.
               FluxL_4p8(i4p8,1)=0.
               FluxU_4p8(i4p8,2)=0.
               FluxL_4p8(i4p8,2)=0.
               f4p8_type(i4p8)='CRATES'
            else if (catalog(1:7) == 'north20') THEN
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_4p8(i4p8,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,2)
               flux_4p8(i4p8,2)=flux_4p8(i4p8,2)*flux2nufnu_4p8*(1.4/4.8)
               frequency_4p8(i4p8,2)=1.4E9
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flag_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,1)
               flux_4p8(i4p8,1)=flux_4p8(i4p8,1)*flux2nufnu_4p8*(4.85/4.8)
               frequency_4p8(i4p8,1)=4.85E9
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if ((flag_4p8(i4p8,1) /= '*') .and. (is .ne. ie-1)) read(string(is+1:ie-1),*)flag_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,3)
               flux_4p8(i4p8,3)=flux_4p8(i4p8,3)*flux2nufnu_4p8*(0.365/4.8)
               frequency_4p8(i4p8,3)=3.65e8
               poserr_4p8(i4p8)=160. !90 accuracy
               FluxU_4p8(i4p8,1)=0.
               FluxL_4p8(i4p8,1)=0.
               FluxU_4p8(i4p8,2)=0.
               FluxL_4p8(i4p8,2)=0.
               FluxU_4p8(i4p8,3)=0.
               FluxL_4p8(i4p8,3)=0.
               f4p8_type(i4p8)='NORTH20'
          else if (catalog(1:7) == 'f357det') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_4p8(i4p8)
               poserr_4p8(i4p8)=2.*poserr_4p8(i4p8)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_4p8(i4p8,2)
               FluxU_4p8(i4p8,2)=flux_4p8(i4p8,2)+Ferr_4p8(i4p8,2)
               FluxL_4p8(i4p8,2)=flux_4p8(i4p8,2)-Ferr_4p8(i4p8,2)
               flux_4p8(i4p8,2)=flux_4p8(i4p8,2)*flux2nufnu_4p8*(3./4.8)
               FluxU_4p8(i4p8,2)=FluxU_4p8(i4p8,2)*flux2nufnu_4p8*(3./4.8)
               FluxL_4p8(i4p8,2)=FluxL_4p8(i4p8,2)*flux2nufnu_4p8*(3./4.8)
               frequency_4p8(i4p8,2)=3.E9
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_4p8(i4p8,1)
               FluxU_4p8(i4p8,1)=flux_4p8(i4p8,1)+Ferr_4p8(i4p8,1)
               FluxL_4p8(i4p8,1)=flux_4p8(i4p8,1)-Ferr_4p8(i4p8,1)
               flux_4p8(i4p8,1)=flux_4p8(i4p8,1)*flux2nufnu_4p8*(5./4.8)
               FluxU_4p8(i4p8,1)=FluxU_4p8(i4p8,1)*flux2nufnu_4p8*(5./4.8)
               FluxL_4p8(i4p8,1)=FluxL_4p8(i4p8,1)*flux2nufnu_4p8*(5./4.8)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flag_4p8(i4p8,1)
               if (flag_4p8(i4p8,1) == '1' ) then
                  flux_4p8(i4p8,1)=0.
                  FluxL_4p8(i4p8,1)=0.
               endif
               frequency_4p8(i4p8,1)=5.E9
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_4p8(i4p8,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_4p8(i4p8,3)
               FluxU_4p8(i4p8,3)=flux_4p8(i4p8,3)+Ferr_4p8(i4p8,3)
               FluxL_4p8(i4p8,3)=flux_4p8(i4p8,3)-Ferr_4p8(i4p8,3)
               flux_4p8(i4p8,3)=flux_4p8(i4p8,3)*flux2nufnu_4p8*(7./4.8)
               FluxU_4p8(i4p8,3)=FluxU_4p8(i4p8,3)*flux2nufnu_4p8*(7./4.8)
               FluxL_4p8(i4p8,3)=FluxL_4p8(i4p8,3)*flux2nufnu_4p8*(7./4.8)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flag_4p8(i4p8,3)
               if (flag_4p8(i4p8,3) == '1' ) then
                  flux_4p8(i4p8,3)=0.
                  FluxL_4p8(i4p8,3)=0.
               endif
               frequency_4p8(i4p8,3)=7.E9
               f4p8_type(i4p8)='F357det'
            endif
               !write(*,*) catalog,FluxU_4p8(i4p8,1),flux_4p8(i4p8,1),FluxL_4p8(i4p8,1)
               !write(*,*) catalog,"Flag: ",flag_4p8(i4p8,1),flag_4p8(i4p8,2),flag_4p8(i4p8,3),flag_4p8(i4p8,4)
          ELSE IF ((catalog(1:7) == 'wish352') .or. (catalog(1:7) == 'tgss150') .or. (catalog(1:5) == 'gleam')
     &           .or. (catalog(1:5) == 'lotss') .or. (catalog(1:5) == 'vlssr')) then
            ilowr=ilowr+1
            ra_lowr(ilowr)=ra
            dec_lowr(ilowr)=dec
            name_l(ilowr)=catalog
            filen_l(ilowr)=filen
            if (ilowr .ne. 1) THEN
               do j=1,ilowr-1
                  if ((ra_lowr(j) .eq. ra) .and. (dec_lowr(j) .eq. dec)) THEN
                     if ((name_l(j) == catalog ) .and. (filen_l(j) .ne. filen)) then
c                        write(*,'(4x,a,i4,3x,2(f9.5,2x),f9.3)') 'The counterpart',
c     &                  j,ra_lowr(j),dec_lowr(j),flux_lowr(j)/(frequency_lowr(j)*1.E-26)
                        is=ie
                        ie=index(string(is+1:len(string)),' ')+is
                        read(string(is+1:ie-1),'(a)') repflux
c                        write(*,'(a,i4,"_",a,2x,2(f9.5,2x),a)') ' Repeated one:',
c     &                   filen,catalog,ra,dec,repflux(1:lenact(repflux))
                        ilowr=ilowr-1
                        goto 300
                     endif
                  endif
               enddo
            endif
            if (catalog(1:7) == 'wish352') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_lowr(ilowr)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_lowr(ilowr)
               FluxU_lowr(ilowr)=flux_lowr(ilowr)+Ferr_lowr(ilowr)
               FluxL_lowr(ilowr)=flux_lowr(ilowr)-Ferr_lowr(ilowr)
               frequency_lowr(ilowr)=3.52e8
               flux_lowr(ilowr)=flux_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
               FluxU_lowr(ilowr)=FluxU_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
               FluxL_lowr(ilowr)=FluxL_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
               posxerr=sqrt(1.5**2+0.12**2)
               posyerr=sqrt(1+0.09**2)
               poserr_lowr(ilowr)=2*sqrt(posxerr**2+posyerr**2)
               lowr_type(ilowr)='WISH352'
           else IF (catalog(1:7) == 'tgss150') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)poserr_lowr(ilowr)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_lowr(ilowr)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_lowr(ilowr)
               FluxU_lowr(ilowr)=flux_lowr(ilowr)+Ferr_lowr(ilowr)
               FluxL_lowr(ilowr)=flux_lowr(ilowr)-Ferr_lowr(ilowr)
               frequency_lowr(ilowr)=1.5e8
               flux_lowr(ilowr)=flux_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
               FluxU_lowr(ilowr)=FluxU_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
               FluxL_lowr(ilowr)=FluxL_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)rms_lowr
               poserr_lowr(ilowr)=poserr_lowr(ilowr)+rms_lowr
               lowr_type(ilowr)='TGSS150'
           else IF (catalog(1:5) == 'gleam') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)poserr_lowr(ilowr)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_lowr(ilowr)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_lowr(ilowr)
               FluxU_lowr(ilowr)=flux_lowr(ilowr)+Ferr_lowr(ilowr)
               FluxL_lowr(ilowr)=flux_lowr(ilowr)-Ferr_lowr(ilowr)
               frequency_lowr(ilowr)=2.e8
               flux_lowr(ilowr)=flux_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
               FluxU_lowr(ilowr)=FluxU_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
               FluxL_lowr(ilowr)=FluxL_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)rms_lowr
               !poserr_lowr(ilowr)=poserr_lowr(ilowr)+rms_lowr
               lowr_type(ilowr)='GLEAM'
           else IF (catalog(1:5) == 'vlssr') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_lowr(ilowr)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_lowr(ilowr)
               FluxU_lowr(ilowr)=flux_lowr(ilowr)+Ferr_lowr(ilowr)
               FluxL_lowr(ilowr)=flux_lowr(ilowr)-Ferr_lowr(ilowr)
               frequency_lowr(ilowr)=7.4e7
               flux_lowr(ilowr)=flux_lowr(ilowr)*frequency_lowr(ilowr)*1.E-23
               FluxU_lowr(ilowr)=FluxU_lowr(ilowr)*frequency_lowr(ilowr)*1.E-23
               FluxL_lowr(ilowr)=FluxL_lowr(ilowr)*frequency_lowr(ilowr)*1.E-23
               poserr_lowr(ilowr)=sqrt((3**2)+(3.4**2))
               lowr_type(ilowr)='VLSSR'
           else IF (catalog(1:5) == 'lotss') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)poserr_lowr(ilowr)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_lowr(ilowr)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_lowr(ilowr)
               FluxU_lowr(ilowr)=flux_lowr(ilowr)+Ferr_lowr(ilowr)
               FluxL_lowr(ilowr)=flux_lowr(ilowr)-Ferr_lowr(ilowr)
               frequency_lowr(ilowr)=1.44e8
               flux_lowr(ilowr)=flux_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
               FluxU_lowr(ilowr)=FluxU_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
               FluxL_lowr(ilowr)=FluxL_lowr(ilowr)*frequency_lowr(ilowr)*1.E-26
           !poserr_lowr(ilowr)=poserr_lowr(ilowr)+rms_lowr
               lowr_type(ilowr)='LoTSS'
           endif
         ELSE IF ( (catalog(1:6) == 'pccs44') .OR. (catalog(1:6) == 'pccs70') .or.
     &             (catalog(1:7) == 'pccs143') .or. (catalog(1:7) == 'pccs100') .or.
     &             (catalog(1:7) == 'pccs217') .or. (catalog(1:7) == 'pccs353') .or.
     &             (catalog(1:4) == 'pcnt') .or. (catalog(1:4) == 'alma'))  THEN
            ipccs100=ipccs100+1
            IF (ipccs100 > arrsize(15)) Stop 'Too many PCCS points'
            ra_pccs100(ipccs100)=ra
            dec_pccs100(ipccs100)=dec
            name_p(ipccs100)=catalog
            filen_p(ipccs100)=filen
            if (ipccs100 .ne. 1) THEN
               do j=1,ipccs100-1
                  if ((ra_pccs100(j) .eq. ra) .and. (dec_pccs100(j) .eq. dec)) THEN
                     if ((name_p(j) == catalog ) .and. (filen_p(j) .ne. filen)) then
c                     write(*,'(4x,a,i4,3x,2(f9.5,2x),f9.3)') 'The counterpart',
c     &                  j,ra_pccs100(j),dec_pccs100(j),flux_pccs100(j,1)/(frequency_pccs100(j,1)*1.E-26)
                     is=ie
                     ie=index(string(is+1:len(string)),' ')+is
                     read(string(is+1:ie-1),'(a)') repflux
c                     write(*,'(a,i4,"_",a,2x,2(f9.5,2x),a)') ' Repeated one:',
c     &                   filen,catalog,ra,dec,repflux(1:lenact(repflux))
                     ipccs100=ipccs100-1
                     goto 300
                     endif
                  endif
               enddo
            endif
            if ((catalog(1:4) == 'pcnt')) then
               frequency_pccs100(ipccs100,1)=3.e10
               frequency_pccs100(ipccs100,2)=4.4e10
               frequency_pccs100(ipccs100,3)=7.e10
               frequency_pccs100(ipccs100,4)=1.e11
               frequency_pccs100(ipccs100,5)=1.43e11
               frequency_pccs100(ipccs100,6)=2.17e11
               frequency_pccs100(ipccs100,7)=3.53e11
               frequency_pccs100(ipccs100,8)=5.45e11
               frequency_pccs100(ipccs100,9)=8.57e11
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_pccs100(ipccs100,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_pccs100(ipccs100,1)
               FluxU_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)+Ferr_pccs100(ipccs100,1)
               FluxL_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)-Ferr_pccs100(ipccs100,1)
               if (Ferr_pccs100(ipccs100,1) .gt. flux_pccs100(ipccs100,1)) then
                  FluxU_pccs100(ipccs100,1)=0.!Ferr_xmm(ixmm,3)*3.
                  FluxL_pccs100(ipccs100,1)=0.
                  flux_pccs100(ipccs100,1)=0.
               endif
               flux_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)*flux2nufnu_pccs100*0.3
               FluxU_pccs100(ipccs100,1)=FluxU_pccs100(ipccs100,1)*flux2nufnu_pccs100*0.3
               FluxL_pccs100(ipccs100,1)=FluxL_pccs100(ipccs100,1)*flux2nufnu_pccs100*0.3
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_pccs100(ipccs100,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_pccs100(ipccs100,2)
               FluxU_pccs100(ipccs100,2)=flux_pccs100(ipccs100,2)+Ferr_pccs100(ipccs100,2)
               FluxL_pccs100(ipccs100,2)=flux_pccs100(ipccs100,2)-Ferr_pccs100(ipccs100,2)
               if (Ferr_pccs100(ipccs100,2) .gt. flux_pccs100(ipccs100,2)) then
                  FluxU_pccs100(ipccs100,2)=0.!Ferr_xmm(ixmm,3)*3.
                  FluxL_pccs100(ipccs100,2)=0.
                  flux_pccs100(ipccs100,2)=0.
               endif
               flux_pccs100(ipccs100,2)=flux_pccs100(ipccs100,2)*flux2nufnu_pccs100*0.44
               FluxU_pccs100(ipccs100,2)=FluxU_pccs100(ipccs100,2)*flux2nufnu_pccs100*0.44
               FluxL_pccs100(ipccs100,2)=FluxL_pccs100(ipccs100,2)*flux2nufnu_pccs100*0.44
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_pccs100(ipccs100,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_pccs100(ipccs100,3)
               FluxU_pccs100(ipccs100,3)=flux_pccs100(ipccs100,3)+Ferr_pccs100(ipccs100,3)
               FluxL_pccs100(ipccs100,3)=flux_pccs100(ipccs100,3)-Ferr_pccs100(ipccs100,3)
               if (Ferr_pccs100(ipccs100,3) .gt. flux_pccs100(ipccs100,3)) then
                  FluxU_pccs100(ipccs100,3)=0.!Ferr_xmm(ixmm,3)*3.
                  FluxL_pccs100(ipccs100,3)=0.
                  flux_pccs100(ipccs100,3)=0.
               endif
               flux_pccs100(ipccs100,3)=flux_pccs100(ipccs100,3)*flux2nufnu_pccs100*0.7
               FluxU_pccs100(ipccs100,3)=FluxU_pccs100(ipccs100,3)*flux2nufnu_pccs100*0.7
               FluxL_pccs100(ipccs100,3)=FluxL_pccs100(ipccs100,3)*flux2nufnu_pccs100*0.7
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_pccs100(ipccs100,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_pccs100(ipccs100,4)
               FluxU_pccs100(ipccs100,4)=flux_pccs100(ipccs100,4)+Ferr_pccs100(ipccs100,4)
               FluxL_pccs100(ipccs100,4)=flux_pccs100(ipccs100,4)-Ferr_pccs100(ipccs100,4)
               if (Ferr_pccs100(ipccs100,4) .gt. flux_pccs100(ipccs100,4)) then
                  FluxU_pccs100(ipccs100,4)=0.!Ferr_xmm(ixmm,3)*3.
                  FluxL_pccs100(ipccs100,4)=0.
                  flux_pccs100(ipccs100,4)=0.
               endif
               flux_pccs100(ipccs100,4)=flux_pccs100(ipccs100,4)*flux2nufnu_pccs100
               FluxU_pccs100(ipccs100,4)=FluxU_pccs100(ipccs100,4)*flux2nufnu_pccs100
               FluxL_pccs100(ipccs100,4)=FluxL_pccs100(ipccs100,4)*flux2nufnu_pccs100
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_pccs100(ipccs100,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_pccs100(ipccs100,5)
               FluxU_pccs100(ipccs100,5)=flux_pccs100(ipccs100,5)+Ferr_pccs100(ipccs100,5)
               FluxL_pccs100(ipccs100,5)=flux_pccs100(ipccs100,5)-Ferr_pccs100(ipccs100,5)
               if (Ferr_pccs100(ipccs100,5) .gt. flux_pccs100(ipccs100,5)) then
                  FluxU_pccs100(ipccs100,5)=0.!Ferr_xmm(ixmm,3)*3.
                  FluxL_pccs100(ipccs100,5)=0.
                  flux_pccs100(ipccs100,5)=0.
               endif
               flux_pccs100(ipccs100,5)=flux_pccs100(ipccs100,5)*flux2nufnu_pccs100*1.43
               FluxU_pccs100(ipccs100,5)=FluxU_pccs100(ipccs100,5)*flux2nufnu_pccs100*1.43
               FluxL_pccs100(ipccs100,5)=FluxL_pccs100(ipccs100,5)*flux2nufnu_pccs100*1.43
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_pccs100(ipccs100,6)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_pccs100(ipccs100,6)
               FluxU_pccs100(ipccs100,6)=flux_pccs100(ipccs100,6)+Ferr_pccs100(ipccs100,6)
               FluxL_pccs100(ipccs100,6)=flux_pccs100(ipccs100,6)-Ferr_pccs100(ipccs100,6)
               if (Ferr_pccs100(ipccs100,6) .gt. flux_pccs100(ipccs100,6)) then
                  FluxU_pccs100(ipccs100,6)=0.!Ferr_xmm(ixmm,3)*3.
                  FluxL_pccs100(ipccs100,6)=0.
                  flux_pccs100(ipccs100,6)=0.
               endif
               flux_pccs100(ipccs100,6)=flux_pccs100(ipccs100,6)*flux2nufnu_pccs100*2.17
               FluxU_pccs100(ipccs100,6)=FluxU_pccs100(ipccs100,6)*flux2nufnu_pccs100*2.17
               FluxL_pccs100(ipccs100,6)=FluxL_pccs100(ipccs100,6)*flux2nufnu_pccs100*2.17
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_pccs100(ipccs100,7)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_pccs100(ipccs100,7)
               FluxU_pccs100(ipccs100,7)=flux_pccs100(ipccs100,7)+Ferr_pccs100(ipccs100,7)
               FluxL_pccs100(ipccs100,7)=flux_pccs100(ipccs100,7)-Ferr_pccs100(ipccs100,7)
               if (Ferr_pccs100(ipccs100,7) .gt. flux_pccs100(ipccs100,7)) then
                  FluxU_pccs100(ipccs100,7)=0.!Ferr_xmm(ixmm,3)*3.
                  FluxL_pccs100(ipccs100,7)=0.
                  flux_pccs100(ipccs100,7)=0.
               endif
               flux_pccs100(ipccs100,7)=flux_pccs100(ipccs100,7)*flux2nufnu_pccs100*3.53
               FluxU_pccs100(ipccs100,7)=FluxU_pccs100(ipccs100,7)*flux2nufnu_pccs100*3.53
               FluxL_pccs100(ipccs100,7)=FluxL_pccs100(ipccs100,7)*flux2nufnu_pccs100*3.53
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_pccs100(ipccs100,8)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_pccs100(ipccs100,8)
               FluxU_pccs100(ipccs100,8)=flux_pccs100(ipccs100,8)+Ferr_pccs100(ipccs100,8)
               FluxL_pccs100(ipccs100,8)=flux_pccs100(ipccs100,8)-Ferr_pccs100(ipccs100,8)
               if (Ferr_pccs100(ipccs100,8) .gt. flux_pccs100(ipccs100,8)) then
                  FluxU_pccs100(ipccs100,8)=0.!Ferr_xmm(ixmm,3)*3.
                  FluxL_pccs100(ipccs100,8)=0.
                  flux_pccs100(ipccs100,8)=0.
               endif
               flux_pccs100(ipccs100,8)=flux_pccs100(ipccs100,8)*flux2nufnu_pccs100*5.45
               FluxU_pccs100(ipccs100,8)=FluxU_pccs100(ipccs100,8)*flux2nufnu_pccs100*5.45
               FluxL_pccs100(ipccs100,8)=FluxL_pccs100(ipccs100,8)*flux2nufnu_pccs100*5.45
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_pccs100(ipccs100,9)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_pccs100(ipccs100,9)
               FluxU_pccs100(ipccs100,9)=flux_pccs100(ipccs100,9)+Ferr_pccs100(ipccs100,9)
               FluxL_pccs100(ipccs100,9)=flux_pccs100(ipccs100,9)-Ferr_pccs100(ipccs100,9)
               if (Ferr_pccs100(ipccs100,9) .gt. flux_pccs100(ipccs100,9)) then
                  FluxU_pccs100(ipccs100,9)=0.!Ferr_xmm(ixmm,3)*3.
                  FluxL_pccs100(ipccs100,9)=0.
                 flux_pccs100(ipccs100,9)=0.
               endif
               flux_pccs100(ipccs100,9)=flux_pccs100(ipccs100,9)*flux2nufnu_pccs100*8.57
               FluxU_pccs100(ipccs100,9)=FluxU_pccs100(ipccs100,9)*flux2nufnu_pccs100*8.57
               FluxL_pccs100(ipccs100,9)=FluxL_pccs100(ipccs100,9)*flux2nufnu_pccs100*8.57
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)snr_pccs100(ipccs100,1)
               if ((snr_pccs100(ipccs100,1) .lt. 3.) .and. (flux_pccs100(ipccs100,1) .ne. 0.)) then
                  flux_pccs100(ipccs100,1)=0.
                  FluxU_pccs100(ipccs100,1)=Ferr_pccs100(ipccs100,1)*flux2nufnu_pccs100*0.3*3.
                  FluxL_pccs100(ipccs100,1)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)snr_pccs100(ipccs100,2)
               if ((snr_pccs100(ipccs100,2) .lt. 3.) .and. (flux_pccs100(ipccs100,2) .ne. 0.)) then
                  flux_pccs100(ipccs100,2)=0.
                  FluxU_pccs100(ipccs100,2)=Ferr_pccs100(ipccs100,2)*flux2nufnu_pccs100*0.44*3.
                  FluxL_pccs100(ipccs100,2)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)snr_pccs100(ipccs100,3)
               if ((snr_pccs100(ipccs100,3) .lt. 3.) .and. (flux_pccs100(ipccs100,3) .ne. 0.)) then
                  flux_pccs100(ipccs100,3)=0.
                  FluxU_pccs100(ipccs100,3)=Ferr_pccs100(ipccs100,3)*flux2nufnu_pccs100*0.7*3.
                  FluxL_pccs100(ipccs100,3)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)snr_pccs100(ipccs100,4)
               if ((snr_pccs100(ipccs100,4) .lt. 3.) .and. (flux_pccs100(ipccs100,4) .ne. 0.)) then
                  flux_pccs100(ipccs100,4)=0.
                  FluxU_pccs100(ipccs100,4)=Ferr_pccs100(ipccs100,4)*flux2nufnu_pccs100*3.
                  FluxL_pccs100(ipccs100,4)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)snr_pccs100(ipccs100,5)
               if ((snr_pccs100(ipccs100,5) .lt. 3.) .and. (flux_pccs100(ipccs100,5) .ne. 0.)) then
                  flux_pccs100(ipccs100,5)=0.
                  FluxU_pccs100(ipccs100,5)=Ferr_pccs100(ipccs100,5)*flux2nufnu_pccs100*1.43*3.
                  FluxL_pccs100(ipccs100,5)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)snr_pccs100(ipccs100,6)
               if ((snr_pccs100(ipccs100,6) .lt. 3.) .and. (flux_pccs100(ipccs100,6) .ne. 0.)) then
                  flux_pccs100(ipccs100,6)=0.
                  FluxU_pccs100(ipccs100,6)=Ferr_pccs100(ipccs100,6)*flux2nufnu_pccs100*2.17*3.
                  FluxL_pccs100(ipccs100,6)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)snr_pccs100(ipccs100,7)
               if ((snr_pccs100(ipccs100,7) .lt. 3.) .and. (flux_pccs100(ipccs100,7) .ne. 0.)) then
                  flux_pccs100(ipccs100,7)=0.
                  FluxU_pccs100(ipccs100,7)=Ferr_pccs100(ipccs100,7)*flux2nufnu_pccs100*3.53*3.
                  FluxL_pccs100(ipccs100,7)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)snr_pccs100(ipccs100,8)
               if ((snr_pccs100(ipccs100,8) .lt. 3.) .and. (flux_pccs100(ipccs100,8) .ne. 0.)) then
                  flux_pccs100(ipccs100,8)=0.
                  FluxU_pccs100(ipccs100,8)=Ferr_pccs100(ipccs100,8)*flux2nufnu_pccs100*5.45*3.
                  FluxL_pccs100(ipccs100,8)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)snr_pccs100(ipccs100,9)
               if ((snr_pccs100(ipccs100,9) .lt. 3.) .and. (flux_pccs100(ipccs100,9) .ne. 0.)) then
                  flux_pccs100(ipccs100,9)=0.
                  FluxU_pccs100(ipccs100,9)=Ferr_pccs100(ipccs100,9)*flux2nufnu_pccs100*8.57*3.
                  FluxL_pccs100(ipccs100,9)=0.
               endif
               poserr_pccs100(ipccs100)=50.
               pccs100_type(ipccs100)='PCNT'
            else if (catalog(1:4) == 'alma') then
               ialma=ialma+1
               almaind(ipccs100)=ialma
c               if ((ra_pccs100(ipccs100) .eq. ra_alma ) .and. (dec_pccs100(ipccs100) .eq. dec_alma)) then
c                  almaind(ipccs100)=0
c               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_pccs100(ipccs100,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_pccs100(ipccs100,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)frequency_pccs100(ipccs100,1)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)')date_alma(ipccs100)
               read(date_alma(ipccs100)(1:4),*) year
               read(date_alma(ipccs100)(6:7),*) month
               read(date_alma(ipccs100)(9:10),*) date
               read(date_alma(ipccs100)(12:13),*) hour
               read(date_alma(ipccs100)(15:16),*) minute
               read(date_alma(ipccs100)(18:19),*) second
               call date_to_mjd(year,month,date,hour,minute,second,mjdtest)
c               is=ie
c               ie=index(string(is+1:len(string)),' ')+is
c               if (is .ne. ie-1) read(string(is+1:ie-1),*)mjdst_alma(ipccs100)
c               mjdst_alma(ipccs100)=mjdst_alma(ipccs100)-2400000.5
c               mjded_alma(ipccs100)=mjdst_alma(ipccs100)
               mjded_alma(ipccs100)=real(mjdtest)
               mjdst_alma(ipccs100)=real(mjdtest)
c               call mjd_to_date(mjdtest,year,month,date,hour)
               FluxU_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)+Ferr_pccs100(ipccs100,1)
               FluxL_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)-Ferr_pccs100(ipccs100,1)
               frequency_pccs100(ipccs100,1)=frequency_pccs100(ipccs100,1)*1.e9
               flux_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)*frequency_pccs100(ipccs100,1)*1.e-23
               FluxU_pccs100(ipccs100,1)=FluxU_pccs100(ipccs100,1)*frequency_pccs100(ipccs100,1)*1.e-23
               FluxL_pccs100(ipccs100,1)=FluxL_pccs100(ipccs100,1)*frequency_pccs100(ipccs100,1)*1.e-23
               poserr_pccs100(ipccs100)=10.
               pccs100_type(ipccs100)='ALMA'
c               ra_alma=ra_pccs100(ipccs100)
c               dec_alma=dec_pccs100(ipccs100)
            else
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux_pccs100(ipccs100,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr_pccs100(ipccs100,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)flux2_pccs100(ipccs100,1)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*)Ferr2_pccs100(ipccs100,1)
               if (Ferr2_pccs100(ipccs100,1) .lt. Ferr_pccs100(ipccs100,1)) then
                  flux_pccs100(ipccs100,1)=flux2_pccs100(ipccs100,1)
                  Ferr_pccs100(ipccs100,1)=Ferr2_pccs100(ipccs100,1)
               endif
               FluxU_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)+Ferr_pccs100(ipccs100,1)
               FluxL_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)-Ferr_pccs100(ipccs100,1)
               flux_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)*flux2nufnu_pccs100
               FluxU_pccs100(ipccs100,1)=FluxU_pccs100(ipccs100,1)*flux2nufnu_pccs100
               FluxL_pccs100(ipccs100,1)=FluxL_pccs100(ipccs100,1)*flux2nufnu_pccs100
               if (FluxL_pccs100(ipccs100,1) .le. 0) then
                  FluxU_pccs100(ipccs100,1)=0.
                  FluxL_pccs100(ipccs100,1)=0.
                  flux_pccs100(ipccs100,1)=0.
               endif
               If (catalog(1:7) == 'pccs100') then
                  frequency_pccs100(ipccs100,1)=1.0e11
               ELSE if (catalog(1:7) == 'pccs44') then
                  flux_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)*0.44 !100 to 143
                  FluxU_pccs100(ipccs100,1)=FluxU_pccs100(ipccs100,1)*0.44
                  FluxL_pccs100(ipccs100,1)=FluxL_pccs100(ipccs100,1)*0.44
                  frequency_pccs100(ipccs100,1)=4.4e10
               ELSE if (catalog(1:7) == 'pccs70') then
                  flux_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)*0.7 !100 to 143
                  FluxU_pccs100(ipccs100,1)=FluxU_pccs100(ipccs100,1)*0.7
                  FluxL_pccs100(ipccs100,1)=FluxL_pccs100(ipccs100,1)*0.7
                  frequency_pccs100(ipccs100,1)=7.e10
               ELSE if (catalog(1:7) == 'pccs143') then
                  flux_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)*1.43 !100 to 143
                  FluxU_pccs100(ipccs100,1)=FluxU_pccs100(ipccs100,1)*1.43
                  FluxL_pccs100(ipccs100,1)=FluxL_pccs100(ipccs100,1)*1.43
                  frequency_pccs100(ipccs100,1)=1.43e11
               ELSE if (catalog(1:7) == 'pccs217') then
                  flux_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)*2.17 !100 to 143
                  FluxU_pccs100(ipccs100,1)=FluxU_pccs100(ipccs100,1)*2.17
                  FluxL_pccs100(ipccs100,1)=FluxL_pccs100(ipccs100,1)*2.17
                  frequency_pccs100(ipccs100,1)=2.17e11
               ELSE if (catalog(1:7) == 'pccs353') then
                  flux_pccs100(ipccs100,1)=flux_pccs100(ipccs100,1)*3.53 !100 to 143
                  FluxU_pccs100(ipccs100,1)=FluxU_pccs100(ipccs100,1)*3.53
                  FluxL_pccs100(ipccs100,1)=FluxL_pccs100(ipccs100,1)*3.53
                  frequency_pccs100(ipccs100,1)=3.53e11
               ENDIF
               if ((flux_pccs100(ipccs100,1) .ne. 0 ) .and.
     &             (flux_pccs100(ipccs100,1)/Ferr_pccs100(ipccs100,1) .gt. 20.)) then
                  if (frequency_pccs100(ipccs100,1) .eq. 4.4E10) poserr_pccs100(ipccs100)=44.35*2.
                  if (frequency_pccs100(ipccs100,1) .eq. 7.E10) poserr_pccs100(ipccs100)=39.69*2.
                  if (frequency_pccs100(ipccs100,1) .eq. 1.E11) poserr_pccs100(ipccs100)=45.8*2.
                  if (frequency_pccs100(ipccs100,1) .eq. 1.43E11) poserr_pccs100(ipccs100)=39.53*2.
                  if (frequency_pccs100(ipccs100,1) .eq. 2.17E11) poserr_pccs100(ipccs100)=38.33*2.
                  if (frequency_pccs100(ipccs100,1) .eq. 3.53E11) poserr_pccs100(ipccs100)=38.57*2.
               ELSE
                  if (frequency_pccs100(ipccs100,1) .eq. 4.4E10) poserr_pccs100(ipccs100)=59.57*2.
                  if (frequency_pccs100(ipccs100,1) .eq. 7.E10) poserr_pccs100(ipccs100)=44.07*2.
                  if (frequency_pccs100(ipccs100,1) .eq. 1.E11) poserr_pccs100(ipccs100)=51.96*2.
                  if (frequency_pccs100(ipccs100,1) .eq. 1.43E11) poserr_pccs100(ipccs100)=43.68*2.
                  if (frequency_pccs100(ipccs100,1) .eq. 2.17E11) poserr_pccs100(ipccs100)=39.94*2.
                  if (frequency_pccs100(ipccs100,1) .eq. 3.53E11) poserr_pccs100(ipccs100)=39.59*2.
               endif
               pccs100_type(ipccs100)='PCCS2'
            !write(*,*) catalog,poserr_pccs100(ipccs100)
            endif
         ELSE IF ((catalog(1:5) == 'spire') .OR. (catalog(1:6) == 'hatlas') .OR.
     &             (catalog(1:8) == 'akaribsc') .OR. (catalog(1:7) == 'iraspsc')) THEN
            ifar=ifar+1
            IF (ifar > arrsize(14)) Stop 'Too many Hershel SPIRE points'
            ra_far(ifar)=ra
            dec_far(ifar)=dec
            name_f(ifar)=catalog
            filen_f(ifar)=filen
            if (ifar .ne. 1) THEN
               do j=1,ifar-1
                  if ((ra_far(j) .eq. ra) .and. (dec_far(j) .eq. dec)) THEN
                     if ((name_f(j) == catalog ) .and. (filen_f(j) .ne. filen)) then
c                     write(*,'(4x,a,i4,3x,2(f9.5,2x),f9.3)') 'The counterpart',
c     &                  j,ra_far(j),dec_far(j),flux_far(j)/(frequency_far(j)*1.E-26)
                     is=ie
                     ie=index(string(is+1:len(string)),' ')+is
                     read(string(is+1:ie-1),'(a)') repflux
c                     write(*,'(a,i4,"_",a,2x,2(f9.5,2x),a)') ' Repeated one:',
c     &                   filen,catalog,ra,dec,repflux(1:lenact(repflux))
                     ifar=ifar-1
                     goto 300
                     endif
                  endif
               enddo
            endif
            if (catalog(1:5) == 'spire') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_far(ifar)
               poserr_far(ifar)=poserr_far(ifar)*2.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,1)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,1)
               if (catalog(1:8) == 'spire250') frequency_far(ifar,1)=(3.E8/2.5e-4)
               if (catalog(1:8) == 'spire350') frequency_far(ifar,1)=(3.E8/3.5e-4)
               if (catalog(1:8) == 'spire500') frequency_far(ifar,1)=(3.E8/5.e-4)
               FluxU_far(ifar,1)=flux_far(ifar,1)+Ferr_far(ifar,1)
               FluxL_far(ifar,1)=flux_far(ifar,1)-Ferr_far(ifar,1)
               flux_far(ifar,1)=flux_far(ifar,1)*frequency_far(ifar,1)*1.E-26
               FluxU_far(ifar,1)=FluxU_far(ifar,1)*frequency_far(ifar,1)*1.E-26
               FluxL_far(ifar,1)=FluxL_far(ifar,1)*frequency_far(ifar,1)*1.E-26
               far_type(ifar)='SPIRE'
            !write(*,*) frequency_far(ifar),FluxU_far(ifar),flux_far(ifar),FluxL_far(ifar)
            Else if (catalog(1:7) == 'iraspsc') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_far(ifar)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,4)
               Ferr_far(ifar,1)=flux_far(ifar,1)*Ferr_far(ifar,1)*1.e-2
               FluxU_far(ifar,1)=flux_far(ifar,1)+Ferr_far(ifar,1)
               FluxL_far(ifar,1)=flux_far(ifar,1)-Ferr_far(ifar,1)
               frequency_far(ifar,1)=(3.E8/1.2e-5)
               flux_far(ifar,1)=flux_far(ifar,1)*frequency_far(ifar,1)*1.E-23
               FluxU_far(ifar,1)=FluxU_far(ifar,1)*frequency_far(ifar,1)*1.E-23
               FluxL_far(ifar,1)=FluxL_far(ifar,1)*frequency_far(ifar,1)*1.E-23
               Ferr_far(ifar,2)=flux_far(ifar,2)*Ferr_far(ifar,2)*1.e-2
               FluxU_far(ifar,2)=flux_far(ifar,2)+Ferr_far(ifar,2)
               FluxL_far(ifar,2)=flux_far(ifar,2)-Ferr_far(ifar,2)
               frequency_far(ifar,2)=(3.E8/2.5e-5)
               flux_far(ifar,2)=flux_far(ifar,2)*frequency_far(ifar,2)*1.E-23
               FluxU_far(ifar,2)=FluxU_far(ifar,2)*frequency_far(ifar,2)*1.E-23
               FluxL_far(ifar,2)=FluxL_far(ifar,2)*frequency_far(ifar,2)*1.E-23
               Ferr_far(ifar,3)=flux_far(ifar,3)*Ferr_far(ifar,3)*1.e-2
               FluxU_far(ifar,3)=flux_far(ifar,3)+Ferr_far(ifar,3)
               FluxL_far(ifar,3)=flux_far(ifar,3)-Ferr_far(ifar,3)
               frequency_far(ifar,3)=(3.E8/6.e-5)
               flux_far(ifar,3)=flux_far(ifar,3)*frequency_far(ifar,3)*1.E-23
               FluxU_far(ifar,3)=FluxU_far(ifar,3)*frequency_far(ifar,3)*1.E-23
               FluxL_far(ifar,3)=FluxL_far(ifar,3)*frequency_far(ifar,3)*1.E-23
               Ferr_far(ifar,4)=flux_far(ifar,4)*Ferr_far(ifar,4)*1.e-2
               FluxU_far(ifar,4)=flux_far(ifar,4)+Ferr_far(ifar,4)
               FluxL_far(ifar,4)=flux_far(ifar,4)-Ferr_far(ifar,4)
               frequency_far(ifar,4)=(3.E8/1.e-4)
               flux_far(ifar,4)=flux_far(ifar,4)*frequency_far(ifar,4)*1.E-23
               FluxU_far(ifar,4)=FluxU_far(ifar,4)*frequency_far(ifar,4)*1.E-23
               FluxL_far(ifar,4)=FluxL_far(ifar,4)*frequency_far(ifar,4)*1.E-23
               far_type(ifar)='IRAS-PSC'
            Else if (catalog(1:8) == 'akaribsc') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_far(ifar)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,3)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,4)
               FluxU_far(ifar,1)=flux_far(ifar,1)+Ferr_far(ifar,1)
               FluxL_far(ifar,1)=flux_far(ifar,1)-Ferr_far(ifar,1)
               frequency_far(ifar,1)=(3.E8/6.e-5)
               flux_far(ifar,1)=flux_far(ifar,1)*frequency_far(ifar,1)*1.E-23
               FluxU_far(ifar,1)=FluxU_far(ifar,1)*frequency_far(ifar,1)*1.E-23
               FluxL_far(ifar,1)=FluxL_far(ifar,1)*frequency_far(ifar,1)*1.E-23
               FluxU_far(ifar,2)=flux_far(ifar,2)+Ferr_far(ifar,2)
               FluxL_far(ifar,2)=flux_far(ifar,2)-Ferr_far(ifar,2)
               frequency_far(ifar,2)=(3.E8/9.e-5)
               flux_far(ifar,2)=flux_far(ifar,2)*frequency_far(ifar,2)*1.E-23
               FluxU_far(ifar,2)=FluxU_far(ifar,2)*frequency_far(ifar,2)*1.E-23
               FluxL_far(ifar,2)=FluxL_far(ifar,2)*frequency_far(ifar,2)*1.E-23
               FluxU_far(ifar,3)=flux_far(ifar,3)+Ferr_far(ifar,3)
               FluxL_far(ifar,3)=flux_far(ifar,3)-Ferr_far(ifar,3)
               frequency_far(ifar,3)=(3.E8/1.4e-4)
               flux_far(ifar,3)=flux_far(ifar,3)*frequency_far(ifar,3)*1.E-23
               FluxU_far(ifar,3)=FluxU_far(ifar,3)*frequency_far(ifar,3)*1.E-23
               FluxL_far(ifar,3)=FluxL_far(ifar,3)*frequency_far(ifar,3)*1.E-23
               FluxU_far(ifar,4)=flux_far(ifar,4)+Ferr_far(ifar,4)
               FluxL_far(ifar,4)=flux_far(ifar,4)-Ferr_far(ifar,4)
               frequency_far(ifar,4)=(3.E8/1.6e-4)
               flux_far(ifar,4)=flux_far(ifar,4)*frequency_far(ifar,4)*1.E-23
               FluxU_far(ifar,4)=FluxU_far(ifar,4)*frequency_far(ifar,4)*1.E-23
               FluxL_far(ifar,4)=FluxL_far(ifar,4)*frequency_far(ifar,4)*1.E-23
               if ((FluxL_far(ifar,1) .lt. 0) .or. (flux_far(ifar,1) .le. 0)) then
                  flux_far(ifar,1)=0.
                  FluxL_far(ifar,1)=0.
               endif
               if ((FluxL_far(ifar,2) .lt. 0) .or. (flux_far(ifar,2) .le. 0)) then
                  flux_far(ifar,2)=0.
                  FluxL_far(ifar,2)=0.
               endif
               if ((FluxL_far(ifar,3) .lt. 0) .or. (flux_far(ifar,3) .le. 0)) then
                  flux_far(ifar,3)=0.
                  FluxL_far(ifar,3)=0.
               endif
               if ((FluxL_far(ifar,4) .lt. 0) .or. (flux_far(ifar,4) .le. 0)) then
                  flux_far(ifar,4)=0.
                  FluxL_far(ifar,4)=0.
               endif
               far_type(ifar)='AKARIBSC'
            ELSE
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,1)
               FluxU_far(ifar,1)=flux_far(ifar,1)+Ferr_far(ifar,1)
               FluxL_far(ifar,1)=flux_far(ifar,1)-Ferr_far(ifar,1)
               frequency_far(ifar,1)=(3.E8/2.5e-4)
               flux_far(ifar,1)=flux_far(ifar,1)*frequency_far(ifar,1)*1.E-23
               FluxU_far(ifar,1)=FluxU_far(ifar,1)*frequency_far(ifar,1)*1.E-23
               FluxL_far(ifar,1)=FluxL_far(ifar,1)*frequency_far(ifar,1)*1.E-23
               if ((FluxL_far(ifar,1) .le. 0.) .or. (FluxL_far(ifar,1) .eq. flux_far(ifar,1))) then
                  flux_far(ifar,1)=0.
                  fluxL_far(ifar,1)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,2)
               FluxU_far(ifar,2)=flux_far(ifar,2)+Ferr_far(ifar,2)
               FluxL_far(ifar,2)=flux_far(ifar,2)-Ferr_far(ifar,2)
               frequency_far(ifar,2)=(3.E8/3.5e-4)
               flux_far(ifar,2)=flux_far(ifar,2)*frequency_far(ifar,2)*1.E-23
               FluxU_far(ifar,2)=FluxU_far(ifar,2)*frequency_far(ifar,2)*1.E-23
               FluxL_far(ifar,2)=FluxL_far(ifar,2)*frequency_far(ifar,2)*1.E-23
               if ((FluxL_far(ifar,2) .le. 0.) .or. (FluxL_far(ifar,2) .eq. flux_far(ifar,2))) then
                  flux_far(ifar,2)=0.
                  fluxL_far(ifar,2)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,3)
               FluxU_far(ifar,3)=flux_far(ifar,3)+Ferr_far(ifar,3)
               FluxL_far(ifar,3)=flux_far(ifar,3)-Ferr_far(ifar,3)
               frequency_far(ifar,3)=(3.E8/5.e-4)
               flux_far(ifar,3)=flux_far(ifar,3)*frequency_far(ifar,3)*1.E-23
               FluxU_far(ifar,3)=FluxU_far(ifar,3)*frequency_far(ifar,3)*1.E-23
               FluxL_far(ifar,3)=FluxL_far(ifar,3)*frequency_far(ifar,3)*1.E-23
               if ((FluxL_far(ifar,3) .le. 0.) .or. (FluxL_far(ifar,3) .eq. flux_far(ifar,3))) then
                  flux_far(ifar,3)=0.
                  fluxL_far(ifar,3)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,4)
               FluxU_far(ifar,4)=flux_far(ifar,4)+Ferr_far(ifar,4)
               FluxL_far(ifar,4)=flux_far(ifar,4)-Ferr_far(ifar,4)
               frequency_far(ifar,4)=(3.E8/1.e-4)
               flux_far(ifar,4)=flux_far(ifar,4)*frequency_far(ifar,4)*1.E-23
               FluxU_far(ifar,4)=FluxU_far(ifar,4)*frequency_far(ifar,4)*1.E-23
               FluxL_far(ifar,4)=FluxL_far(ifar,4)*frequency_far(ifar,4)*1.E-23
               if ((FluxL_far(ifar,4) .le. 0.) .or. (FluxL_far(ifar,4) .eq. flux_far(ifar,4))) then
                  flux_far(ifar,4)=0.
                  fluxL_far(ifar,4)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_far(ifar,5)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_far(ifar,5)
               FluxU_far(ifar,5)=flux_far(ifar,5)+Ferr_far(ifar,5)
               FluxL_far(ifar,5)=flux_far(ifar,5)-Ferr_far(ifar,5)
               frequency_far(ifar,5)=(3.E8/1.6e-4)
               flux_far(ifar,5)=flux_far(ifar,5)*frequency_far(ifar,5)*1.E-23
               FluxU_far(ifar,5)=FluxU_far(ifar,5)*frequency_far(ifar,5)*1.E-23
               FluxL_far(ifar,5)=FluxL_far(ifar,5)*frequency_far(ifar,5)*1.E-23
               if ((FluxL_far(ifar,5) .le. 0.) .or. (FluxL_far(ifar,5) .eq. flux_far(ifar,5))) then
                  flux_far(ifar,5)=0.
                  fluxL_far(ifar,5)=0.
               endif
               poserr_far(ifar)=3.
               if (catalog(1:7) == 'hatlas1') far_type(ifar)='H-ATLAS-DR1'
               if (catalog(1:7) == 'hatlas2') far_type(ifar)='H-ATLAS-DR2'
            endif
         ELSE IF ((catalog(1:4) == 'wise') .OR. (catalog(1:7) == 'neowise')
     &            .or.  (catalog(1:5) == '2mass') ) THEN
            iir=iir+1
            IF (iir > arrsize(13)) Stop 'Too many Infrared points'
            if (iir .ne. 1) THEN
               do j=1,iir-1
                  if ((ra_ir(j) .eq. ra) .and. (dec_ir(j) .eq. dec)) THEN
                     if ((name_i(j) == catalog ) .and. (filen_i(j) .ne. filen)) then
c                     write(*,'(4x,"The counterpart",i4,3x,2(f9.5,2x),4(f6.3,2x))')
c     &                  j,ra_ir(j),dec_ir(j),irmag(j,1:4)
                     is=ie
                     ie=index(string(is+1:len(string)),' ')+is
                     read(string(is+1:ie-1),'(a)') repflux
c                     write(*,'(a,i4,"_",a,2x,2(f9.5,2x),a)') ' Repeated one:',
c     &                   filen,catalog,ra,dec,repflux(1:lenact(repflux))
                     iir=iir-1
                     goto 300
                     endif
                  endif
               enddo
            endif
            ra_ir(iir)=ra
            dec_ir(iir)=dec
            name_i(iir)=catalog
            filen_i(iir)=filen
            IF (catalog(1:4) == 'wise') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (catalog(1:5) == 'wise ') then
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_ir(iir)
                  poserr_ir(iir)=poserr_ir(iir)*2.
               else
                  irlc_name(iir)='nan'
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) irlc_name(iir)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdst_irlc(iir)
                  mjded_irlc(iir)=mjdst_irlc(iir)
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmag(iir,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmagerr(iir,1)
               call mag2flux (nh,irmag(iir,1),'ww1',flux_ir(iir,1),frequency_ir(iir,1))
               call mag2flux (nh,irmag(iir,1)-irmagerr(iir,1),'ww1',FluxU_ir(iir,1),frequency_ir(iir,1))
               call mag2flux (nh,irmag(iir,1)+irmagerr(iir,1),'ww1',FluxL_ir(iir,1),frequency_ir(iir,1))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmag(iir,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmagerr(iir,2)
               call mag2flux (nh,irmag(iir,2),'ww2',flux_ir(iir,2),frequency_ir(iir,2))
               call mag2flux (nh,irmag(iir,2)-irmagerr(iir,2),'ww2',FluxU_ir(iir,2),frequency_ir(iir,2))
               call mag2flux (nh,irmag(iir,2)+irmagerr(iir,2),'ww2',FluxL_ir(iir,2),frequency_ir(iir,2))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmag(iir,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmagerr(iir,3)
               call mag2flux (nh,irmag(iir,3),'ww3',flux_ir(iir,3),frequency_ir(iir,3))
               call mag2flux (nh,irmag(iir,3)-irmagerr(iir,3),'ww3',FluxU_ir(iir,3),frequency_ir(iir,3))
               call mag2flux (nh,irmag(iir,3)+irmagerr(iir,3),'ww3',FluxL_ir(iir,3),frequency_ir(iir,3))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmag(iir,4)
               if (catalog(1:5) == 'wise ') then
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
               else
                  is=ie
                  ie=index(string(is+1:len(string)),' ')+is
               ENDIF
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmagerr(iir,4)
               call mag2flux (nh,irmag(iir,4),'ww4',flux_ir(iir,4),frequency_ir(iir,4))
               call mag2flux (nh,irmag(iir,4)-irmagerr(iir,4),'ww4',FluxU_ir(iir,4),frequency_ir(iir,4))
               call mag2flux (nh,irmag(iir,4)+irmagerr(iir,4),'ww4',FluxL_ir(iir,4),frequency_ir(iir,4))
               if (catalog(1:5) == 'wise ') then
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_ir(iir,2)
                  is=ie
                  ie=index(string(is+1:len(string)),' ')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') flag_ir(iir,1)
                  ir_type(iir)='WISE'
                  if (flag_ir(iir,1)(1:1) == 'C') flux_ir(iir,1)=-flux_ir(iir,1)
                  if (flag_ir(iir,1)(2:2) == 'C') flux_ir(iir,2)=-flux_ir(iir,2)
                  if (flag_ir(iir,1)(3:3) == 'C') flux_ir(iir,3)=-flux_ir(iir,3)
                  if (flag_ir(iir,1)(4:4) == 'C') flux_ir(iir,4)=-flux_ir(iir,4)
                else
                  ir_type(iir)='WISEME'
                  if ((irlcid(1) == 'First')  .and. (irlc_name(iir) /= 'nan') ) irlcid(1)=irlc_name(iir)
                  iirlc=iirlc+1
                  indirlc(iirlc)=iir
                endif
               !write(*,*) catalog(1:4),flux_ir(iir,1),flux_ir(iir,2),flux_ir(iir,3),flux_ir(iir,4)
            ELSE if (catalog(1:7) == 'neowise') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmag(iir,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmagerr(iir,1)
               call mag2flux (nh,irmag(iir,1),'ww1',flux_ir(iir,1),frequency_ir(iir,1))
               call mag2flux (nh,irmag(iir,1)-irmagerr(iir,1),'ww1',FluxU_ir(iir,1),frequency_ir(iir,1))
               call mag2flux (nh,irmag(iir,1)+irmagerr(iir,1),'ww1',FluxL_ir(iir,1),frequency_ir(iir,1))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmag(iir,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmagerr(iir,2)
               call mag2flux (nh,irmag(iir,2),'ww2',flux_ir(iir,2),frequency_ir(iir,2))
               call mag2flux (nh,irmag(iir,2)-irmagerr(iir,2),'ww2',FluxU_ir(iir,2),frequency_ir(iir,2))
               call mag2flux (nh,irmag(iir,2)+irmagerr(iir,2),'ww2',FluxL_ir(iir,2),frequency_ir(iir,2))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdst_irlc(iir)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               irlc_name(iir)='nan'
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irlc_name(iir)
               mjded_irlc(iir)=mjdst_irlc(iir)
               if ((irlcid(2) == 'First')  .and. (irlc_name(iir) /= 'nan') ) irlcid(2)=irlc_name(iir)
               ir_type(iir)='NEOWISE'
               iirlc=iirlc+1
               indirlc(iirlc)=iir
            ELSE IF (catalog(1:5) == '2mass') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_ir(iir)
               poserr_ir(iir)=poserr_ir(iir)*2.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmag(iir,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmagerr(iir,4)
               if (irmagerr(iir,4) .le. -99.) irmagerr(iir,4)=0.
               CALL  mag2flux (nh,irmag(iir,4),'J  ',flux_ir(iir,4),frequency_ir(iir,4))
               CALL  mag2flux (nh,irmag(iir,4)-irmagerr(iir,4),'J  ',FluxU_ir(iir,4),frequency_ir(iir,4))
               CALL  mag2flux (nh,irmag(iir,4)+irmagerr(iir,4),'J  ',FluxL_ir(iir,4),frequency_ir(iir,4))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmag(iir,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmagerr(iir,3)
               if (irmagerr(iir,3) .le. -99.) irmagerr(iir,3)=0.
               call  mag2flux (nh,irmag(iir,3),'H  ',flux_ir(iir,3),frequency_ir(iir,3))
               call  mag2flux (nh,irmag(iir,3)-irmagerr(iir,3),'H  ',FluxU_ir(iir,3),frequency_ir(iir,3))
               call  mag2flux (nh,irmag(iir,3)+irmagerr(iir,3),'H  ',FluxL_ir(iir,3),frequency_ir(iir,3))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmag(iir,2)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) irmagerr(iir,2)
               if (irmagerr(iir,2) .le. -99.) irmagerr(iir,2)=0.
               call  mag2flux (nh,irmag(iir,2),'K  ',flux_ir(iir,2),frequency_ir(iir,2))
               call  mag2flux (nh,irmag(iir,2)-irmagerr(iir,2),'K  ',FluxU_ir(iir,2),frequency_ir(iir,2))
               call  mag2flux (nh,irmag(iir,2)+irmagerr(iir,2),'K  ',FluxL_ir(iir,2),frequency_ir(iir,2))
               ir_type(iir)='2MASS'
            ENDIF
            !write(*,*) catalog,FluxU_ir(iir,3),flux_ir(iir,3),FluxL_ir(iir,3)
         ELSE IF ( (catalog(1:4) == 'usno') .OR. (catalog(1:3) == 'hst') .or.
     &           (catalog(1:4) == 'sdss') .or. (catalog(1:9) == 'panstarrs') .or.
     &           (catalog(1:4) == 'gaia') .or. (catalog(1:6) == 'smarts')) THEN
            iusno = iusno + 1
            IF (iusno > arrsize(12)) Stop 'Too many USNO points'
            if (iusno .ne. 1) THEN
               do j=1,iusno-1
                  if ((ra_usno(j) .eq. ra) .and. (dec_usno(j) .eq. dec)) THEN
                     if ((name_o(j) == catalog ) .and. (filen_o(j) .ne. filen)) then
c                     write(*,'(4x,"The counterpart",i4,3x,2(f9.5,2x),4(f6.3,2x))')
c     &                  j,ra_usno(j),dec_usno(j),usnomag(j,1:5)
                     is=ie
                     ie=index(string(is+1:len(string)),' ')+is
                     read(string(is+1:ie-1),'(a)') repflux
c                     write(*,'(a,i4,"_",a,2x,2(f9.5,2x),a)') ' Repeated one:',
c     &                   filen,catalog,ra,dec,repflux(1:lenact(repflux))
                     iusno=iusno-1
                     goto 300
                     endif
                  endif
               enddo
            endif
            ra_usno(iusno)=ra
            dec_usno(iusno)=dec
            name_o(iusno)=catalog
            filen_o(iusno)=filen
            IF (catalog(1:4) == 'usno') THEN !tell if the data is from SDSS or USNO
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,2)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,4)
               CALL  mag2flux (nh,usnomag(iusno,4),'R  ',flux_usno(iusno,4),frequency_usno(iusno,4))
               CALL  mag2flux (nh,usnomag(iusno,2),'B  ',flux_usno(iusno,2),frequency_usno(iusno,2))
               FluxU_usno(iusno,1:5)=0.
               FluxL_usno(iusno,1:5)=0.
               opt_type(iusno)='USNO'
               poserr_usno(iusno)=0.5
            ELSE if (catalog(1:4) == 'sdss') THEN
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,4)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,5)
              CALL  mag2flux (nh,usnomag(iusno,1),'u  ',flux_usno(iusno,1),frequency_usno(iusno,1))
         CALL  mag2flux (nh,usnomag(iusno,1)-usnomagerr(iusno,1),'u  ',FluxU_usno(iusno,1),frequency_usno(iusno,1))
         CALL  mag2flux (nh,usnomag(iusno,1)+usnomagerr(iusno,1),'u  ',FluxL_usno(iusno,1),frequency_usno(iusno,1))
              CALL  mag2flux (nh,usnomag(iusno,2),'g  ',flux_usno(iusno,2),frequency_usno(iusno,2))
         CALL  mag2flux (nh,usnomag(iusno,2)-usnomagerr(iusno,2),'g  ',FluxU_usno(iusno,2),frequency_usno(iusno,2))
         CALL  mag2flux (nh,usnomag(iusno,2)+usnomagerr(iusno,2),'g  ',FluxL_usno(iusno,2),frequency_usno(iusno,2))
              CALL  mag2flux (nh,usnomag(iusno,3),'r  ',flux_usno(iusno,3),frequency_usno(iusno,3))
         CALL  mag2flux (nh,usnomag(iusno,3)-usnomagerr(iusno,3),'r  ',FluxU_usno(iusno,3),frequency_usno(iusno,3))
         CALL  mag2flux (nh,usnomag(iusno,3)+usnomagerr(iusno,3),'r  ',FluxL_usno(iusno,3),frequency_usno(iusno,3))
              CALL  mag2flux (nh,usnomag(iusno,4),'i  ',flux_usno(iusno,4),frequency_usno(iusno,4))
         CALL  mag2flux (nh,usnomag(iusno,4)-usnomagerr(iusno,4),'i  ',FluxU_usno(iusno,4),frequency_usno(iusno,4))
         CALL  mag2flux (nh,usnomag(iusno,4)+usnomagerr(iusno,4),'i  ',FluxL_usno(iusno,4),frequency_usno(iusno,4))
               CALL  mag2flux (nh,usnomag(iusno,5),'z  ',flux_usno(iusno,5),frequency_usno(iusno,5))
         CALL  mag2flux (nh,usnomag(iusno,5)-usnomagerr(iusno,5),'z  ',FluxU_usno(iusno,5),frequency_usno(iusno,5))
         CALL  mag2flux (nh,usnomag(iusno,5)+usnomagerr(iusno,5),'z  ',FluxL_usno(iusno,5),frequency_usno(iusno,5))
c checked photometric quality for SDSS ! no upper limit for SDSS
               if (usnomag(iusno,1) .gt. 24.63) THEN
                  flux_usno(iusno,1)=0.
                  fluxU_usno(iusno,1)=0.
                  fluxL_usno(iusno,1)=0.
               else if ((usnomag(iusno,1) .gt. 22.12) .and. (usnomag(iusno,1) .lt. 24.63)) then
                  flux_usno(iusno,1)=-flux_usno(iusno,1)
               endif
               if (usnomag(iusno,2) .gt. 25.11) THEN
                  flux_usno(iusno,2)=0.
                  fluxU_usno(iusno,2)=0.
                  fluxL_usno(iusno,2)=0.
               else if ((usnomag(iusno,2) .gt. 22.6) .and. (usnomag(iusno,2) .lt. 25.11)) then
                  flux_usno(iusno,2)=-flux_usno(iusno,2)
               endif
               if (usnomag(iusno,3) .gt. 24.8) THEN
                  flux_usno(iusno,3)=0.
                  fluxU_usno(iusno,3)=0.
                  fluxL_usno(iusno,3)=0.
               else if ((usnomag(iusno,3) .gt. 22.29) .and. (usnomag(iusno,3) .lt. 24.8)) then
                  flux_usno(iusno,3)=-flux_usno(iusno,3)
               endif
               if (usnomag(iusno,4) .gt. 24.36) THEN
                  flux_usno(iusno,4)=0.
                  fluxU_usno(iusno,4)=0.
                  fluxL_usno(iusno,4)=0.
               else if ((usnomag(iusno,4) .gt. 21.85) .and. (usnomag(iusno,4) .lt. 24.36)) then
                  flux_usno(iusno,4)=-flux_usno(iusno,4)
               endif
               if (usnomag(iusno,5) .gt. 22.83) THEN
                  flux_usno(iusno,5)=0.
                  fluxU_usno(iusno,5)=0.
                  fluxL_usno(iusno,5)=0.
               else if ((usnomag(iusno,5) .gt. 20.32) .and. (usnomag(iusno,5) .lt. 22.83)) then
                  flux_usno(iusno,5)=-flux_usno(iusno,5)
               endif
               opt_type(iusno)='SDSS'
               poserr_usno(iusno)=0.5
            else if (catalog(1:3) == 'hst') THEN
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_usno(iusno)
               poserr_usno(iusno)=poserr_usno(iusno)*2.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,1)
               if (usnomag(iusno,1) .gt. 99.) usnomag(iusno,1)=0.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,1)
               if (usnomagerr(iusno,1) .gt. 99.) usnomagerr(iusno,1)=0.
               CALL  mag2flux (nh,usnomag(iusno,1),'U  ',flux_usno(iusno,1),frequency_usno(iusno,1))
         CALL  mag2flux (nh,usnomag(iusno,1)-usnomagerr(iusno,1),'U  ',FluxU_usno(iusno,1),frequency_usno(iusno,1))
         CALL  mag2flux (nh,usnomag(iusno,1)+usnomagerr(iusno,1),'U  ',FluxL_usno(iusno,1),frequency_usno(iusno,1))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,2)
               if (usnomag(iusno,2) .gt. 99.) usnomag(iusno,2)=0.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,2)
               if (usnomagerr(iusno,2) .gt. 99.) usnomagerr(iusno,2)=0.
               CALL  mag2flux (nh,usnomag(iusno,2),'B  ',flux_usno(iusno,2),frequency_usno(iusno,2))
         CALL  mag2flux(nh,usnomag(iusno,2)-usnomagerr(iusno,2),'B  ',FluxU_usno(iusno,2),frequency_usno(iusno,2))
         CALL  mag2flux(nh,usnomag(iusno,2)+usnomagerr(iusno,2),'B  ',FluxL_usno(iusno,2),frequency_usno(iusno,2))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,3)
               if (usnomag(iusno,3) .gt. 99.) usnomag(iusno,3)=0.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,3)
               if (usnomagerr(iusno,3) .gt. 99.) usnomagerr(iusno,3)=0.
               CALL  mag2flux (nh,usnomag(iusno,3),'V  ',flux_usno(iusno,3),frequency_usno(iusno,3))
         CALL  mag2flux(nh,usnomag(iusno,3)-usnomagerr(iusno,3),'V  ',FluxU_usno(iusno,3),frequency_usno(iusno,3))
         CALL  mag2flux(nh,usnomag(iusno,3)+usnomagerr(iusno,3),'V  ',FluxL_usno(iusno,3),frequency_usno(iusno,3))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,4)
               if (usnomag(iusno,4) .gt. 99.) usnomag(iusno,4)=0.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,4)
               if (usnomagerr(iusno,4) .gt. 99.) usnomagerr(iusno,4)=0.
               CALL  mag2flux (nh,usnomag(iusno,4),'R  ',flux_usno(iusno,4),frequency_usno(iusno,4))
         CALL  mag2flux (nh,usnomag(iusno,4)-usnomagerr(iusno,4),'R  ',FluxU_usno(iusno,4),frequency_usno(iusno,4))
         CALL  mag2flux (nh,usnomag(iusno,4)+usnomagerr(iusno,4),'R  ',FluxL_usno(iusno,4),frequency_usno(iusno,4))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,5)
               if (usnomag(iusno,5) .gt. 99.) usnomag(iusno,5)=0.
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,5)
               if (usnomagerr(iusno,5) .gt. 99.) usnomagerr(iusno,5)=0.
               CALL  mag2flux (nh,usnomag(iusno,5),'I  ',flux_usno(iusno,5),frequency_usno(iusno,5))
         CALL  mag2flux (nh,usnomag(iusno,5)-usnomagerr(iusno,5),'I  ',FluxU_usno(iusno,5),frequency_usno(iusno,5))
         CALL  mag2flux (nh,usnomag(iusno,5)+usnomagerr(iusno,5),'I  ',FluxL_usno(iusno,5),frequency_usno(iusno,5))
               opt_type(iusno)='HSTGSC'
            else if (catalog(1:9) == 'panstarrs') THEN
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_usno(iusno)
               poserr_usno(iusno)=2.*poserr_usno(iusno)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,1)
               CALL  mag2flux (nh,usnomag(iusno,1),'psg',flux_usno(iusno,1),frequency_usno(iusno,1))
        CALL  mag2flux (nh,usnomag(iusno,1)-usnomagerr(iusno,1),'psg',FluxU_usno(iusno,1),frequency_usno(iusno,1))
        CALL  mag2flux (nh,usnomag(iusno,1)+usnomagerr(iusno,1),'psg',FluxL_usno(iusno,1),frequency_usno(iusno,1))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,2)
               CALL  mag2flux (nh,usnomag(iusno,2),'psr',flux_usno(iusno,2),frequency_usno(iusno,2))
        CALL  mag2flux(nh,usnomag(iusno,2)-usnomagerr(iusno,2),'psr',FluxU_usno(iusno,2),frequency_usno(iusno,2))
        CALL  mag2flux(nh,usnomag(iusno,2)+usnomagerr(iusno,2),'psr',FluxL_usno(iusno,2),frequency_usno(iusno,2))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,3)
               CALL  mag2flux (nh,usnomag(iusno,3),'psi',flux_usno(iusno,3),frequency_usno(iusno,3))
        CALL  mag2flux(nh,usnomag(iusno,3)-usnomagerr(iusno,3),'psi',FluxU_usno(iusno,3),frequency_usno(iusno,3))
        CALL  mag2flux(nh,usnomag(iusno,3)+usnomagerr(iusno,3),'psi',FluxL_usno(iusno,3),frequency_usno(iusno,3))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,4)
               CALL  mag2flux (nh,usnomag(iusno,4),'psz',flux_usno(iusno,4),frequency_usno(iusno,4))
        CALL  mag2flux (nh,usnomag(iusno,4)-usnomagerr(iusno,4),'psz',FluxU_usno(iusno,4),frequency_usno(iusno,4))
        CALL  mag2flux (nh,usnomag(iusno,4)+usnomagerr(iusno,4),'psz',FluxL_usno(iusno,4),frequency_usno(iusno,4))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,5)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,5)
               CALL  mag2flux (nh,usnomag(iusno,5),'psy',flux_usno(iusno,5),frequency_usno(iusno,5))
        CALL  mag2flux (nh,usnomag(iusno,5)-usnomagerr(iusno,5),'psy',FluxU_usno(iusno,5),frequency_usno(iusno,5))
        CALL  mag2flux (nh,usnomag(iusno,5)+usnomagerr(iusno,5),'psy',FluxL_usno(iusno,5),frequency_usno(iusno,5))
               opt_type(iusno)='PanSTARRS'
            else if (catalog(1:4) == 'gaia') THEN
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_usno(iusno)
               poserr_usno(iusno)=2.*sqrt(poserr_usno(iusno)**2+(0.0003**2))
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomag(iusno,1)
               !write(*,*) 'TEST GAIA',poserr_usno(iusno),usnomag(iusno,1)
               CALL  mag2flux (nh,usnomag(iusno,1),'bbG',flux_usno(iusno,1),frequency_usno(iusno,1))
               FluxU_usno(iusno,1)=0.
               FluxL_usno(iusno,1)=0.
               opt_type(iusno)='GAIA'
            else if (catalog(1:6) == 'smarts') THEN
               ioptlc=ioptlc+1
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) frequency_usno(iusno,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_usno(iusno,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) usnomagerr(iusno,1)
               FluxU_usno(iusno,1)=flux_usno(iusno,1)+usnomagerr(iusno,1)
               FluxL_usno(iusno,1)=flux_usno(iusno,1)-usnomagerr(iusno,1)
               if (FluxL_usno(iusno,1) .lt. 0) then
                  flux_usno(iusno,1)=0.
                  FluxL_usno(iusno,1)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdst_optlc(iusno)
               poserr_usno(iusno)=10.
               indoptlc(ioptlc)=iusno
               !write(*,*) 'TEST GAIA',poserr_usno(iusno),usnomag(iusno,1)
               opt_type(iusno)='SMARTS'
            ENDIF
            !write(*,*) usnomag(iusno,1)-usnomagerr(iusno,1),usnomag(iusno,1),usnomag(iusno,1)+usnomagerr(iusno,1)
            !write(*,*) catalog,FluxU_usno(iusno,1),flux_usno(iusno,1),FluxL_usno(iusno,1),poserr_usno(iusno)
            !testmag1=15.688
            !testmag2=15.691
            !CALL  mag2flux (nh,testmag1,'g  ',testmag1,ftest1)
            !CALL  mag2flux (nh,testmag2,'g  ',testmag2,ftest1)
            !write(*,*) 'test SDSS',testmag1,testmag2
         ELSE IF  ( (catalog(1:5) == 'galex') .OR. (catalog(1:5) == 'xmmom') .or.
     &              (catalog(1:4) == 'uvot') )  THEN
            iuv=iuv+1
            if (iuv > arrsize(10)) Stop 'Too many UV points'
            if (iuv .ne. 1) THEN
               do j=1,iuv-1
                  if ((ra_uv(j) .eq. ra) .and. (dec_uv(j) .eq. dec)) THEN
                     if ((name_u(j) == catalog ) .and. (filen_u(j) .ne. filen)) then
c                     write(*,'(4x,"The counterpart",i4,3x,2(f9.5,2x),6(f6.3,2x))')
c     &                  j,ra_uv(j),dec_uv(j),uvmag(j,1:6)
                     is=ie
                     ie=index(string(is+1:len(string)),' ')+is
                     read(string(is+1:ie-1),'(a)') repflux
c                     write(*,'(a,i4,"_",a,2x,2(f9.5,2x),a)') ' Repeated one:',
c     &                   filen,catalog,ra,dec,repflux(1:lenact(repflux))
                     iuv=iuv-1
                     goto 300
                     endif
                  endif
               enddo
            endif
            ra_uv(iuv)=ra
            dec_uv(iuv)=dec
            name_u(iuv)=catalog
            filen_u(iuv)=filen
            IF (catalog(1:5) == 'galex') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               read(string(is+1:ie-1),*) uvmag(iuv,4)
               if (uvmag(iuv,4) .le. -999.) uvmag(iuv,4)=0.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               read(string(is+1:ie-1),*) uvmag(iuv,3)
               if (uvmag(iuv,3) .le. -999.) uvmag(iuv,3)=0.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               read(string(is+1:ie-1),*) uvmagerr(iuv,4)
               if (uvmagerr(iuv,4) .le. -999.) uvmagerr(iuv,4)=0.
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               read(string(is+1:ie-1),*) uvmagerr(iuv,3)
               if (uvmagerr(iuv,3) .le. -999.) uvmagerr(iuv,3)=0.
               CALL  mag2flux (nh,uvmag(iuv,3),'fuv',flux_uv(iuv,3),frequency_uv(iuv,3))
               CALL  mag2flux (nh,uvmag(iuv,3)-uvmagerr(iuv,3),'fuv',FluxU_uv(iuv,3),frequency_uv(iuv,3))
               CALL  mag2flux (nh,uvmag(iuv,3)+uvmagerr(iuv,3),'fuv',FluxL_uv(iuv,3),frequency_uv(iuv,3))
               CALL  mag2flux (nh,uvmag(iuv,4),'nuv',flux_uv(iuv,4),frequency_uv(iuv,4))
               CALL  mag2flux (nh,uvmag(iuv,4)-uvmagerr(iuv,4),'nuv',FluxU_uv(iuv,4),frequency_uv(iuv,4))
               CALL  mag2flux (nh,uvmag(iuv,4)+uvmagerr(iuv,4),'nuv',FluxL_uv(iuv,4),frequency_uv(iuv,4))
               uv_type(iuv)='GALEX'
               poserr_uv(iuv)=1.
            ELSE IF (catalog(1:4) == 'uvot') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_uv(iuv)
               poserr_uv(iuv)=2*sqrt(poserr_uv(iuv)**2+(0.26**2)) !systematic error 0.26
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,6)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,6)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,3)
               !write(*,*) uvmag(iuv,1),uvmag(iuv,2),uvmag(iuv,3),uvmag(iuv,4),uvmag(iuv,5),uvmag(iuv,6)
               CALL  mag2flux (nh,uvmag(iuv,1),'su ',flux_uv(iuv,1),frequency_uv(iuv,1))
               CALL  mag2flux (nh,uvmag(iuv,2),'sb ',flux_uv(iuv,2),frequency_uv(iuv,2))
               CALL  mag2flux (nh,uvmag(iuv,3),'sv ',flux_uv(iuv,3),frequency_uv(iuv,3))
               CALL  mag2flux (nh,uvmag(iuv,4),'sw1',flux_uv(iuv,4),frequency_uv(iuv,4))
               CALL  mag2flux (nh,uvmag(iuv,5),'sm2',flux_uv(iuv,5),frequency_uv(iuv,5))
               CALL  mag2flux (nh,uvmag(iuv,6),'sw2',flux_uv(iuv,6),frequency_uv(iuv,6))
               CALL  mag2flux (nh,uvmag(iuv,1)-uvmagerr(iuv,1),'su ',FluxU_uv(iuv,1),frequency_uv(iuv,1))
               CALL  mag2flux (nh,uvmag(iuv,2)-uvmagerr(iuv,2),'sb ',FluxU_uv(iuv,2),frequency_uv(iuv,2))
               CALL  mag2flux (nh,uvmag(iuv,3)-uvmagerr(iuv,3),'sv ',FluxU_uv(iuv,3),frequency_uv(iuv,3))
               CALL  mag2flux (nh,uvmag(iuv,4)-uvmagerr(iuv,4),'sw1',FluxU_uv(iuv,4),frequency_uv(iuv,4))
               CALL  mag2flux (nh,uvmag(iuv,5)-uvmagerr(iuv,5),'sm2',FluxU_uv(iuv,5),frequency_uv(iuv,5))
               CALL  mag2flux (nh,uvmag(iuv,6)-uvmagerr(iuv,6),'sw2',FluxU_uv(iuv,6),frequency_uv(iuv,6))
               CALL  mag2flux (nh,uvmag(iuv,1)+uvmagerr(iuv,1),'su ',FluxL_uv(iuv,1),frequency_uv(iuv,1))
               CALL  mag2flux (nh,uvmag(iuv,2)+uvmagerr(iuv,2),'sb ',FluxL_uv(iuv,2),frequency_uv(iuv,2))
               CALL  mag2flux (nh,uvmag(iuv,3)+uvmagerr(iuv,3),'sv ',FluxL_uv(iuv,3),frequency_uv(iuv,3))
               CALL  mag2flux (nh,uvmag(iuv,4)+uvmagerr(iuv,4),'sw1',FluxL_uv(iuv,4),frequency_uv(iuv,4))
               CALL  mag2flux (nh,uvmag(iuv,5)+uvmagerr(iuv,5),'sm2',FluxL_uv(iuv,5),frequency_uv(iuv,5))
               CALL  mag2flux (nh,uvmag(iuv,6)+uvmagerr(iuv,6),'sw2',FluxL_uv(iuv,6),frequency_uv(iuv,6))
               uv_type(iuv)='UVOT'
            ELSE IF (catalog(1:5) == 'xmmom') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_uv(iuv)
               poserr_uv(iuv)=2.*sqrt(poserr_uv(iuv)*poserr_uv(iuv)+(0.7**2))
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,6)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,6)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmag(iuv,3)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) uvmagerr(iuv,3)
               CALL  mag2flux (nh,uvmag(iuv,1),'xu ',flux_uv(iuv,1),frequency_uv(iuv,1))
               CALL  mag2flux (nh,uvmag(iuv,2),'xb ',flux_uv(iuv,2),frequency_uv(iuv,2))
               CALL  mag2flux (nh,uvmag(iuv,3),'xv ',flux_uv(iuv,3),frequency_uv(iuv,3))
               CALL  mag2flux (nh,uvmag(iuv,4),'xw1',flux_uv(iuv,4),frequency_uv(iuv,4))
               CALL  mag2flux (nh,uvmag(iuv,5),'xm2',flux_uv(iuv,5),frequency_uv(iuv,5))
               CALL  mag2flux (nh,uvmag(iuv,6),'xw2',flux_uv(iuv,6),frequency_uv(iuv,6))
               CALL  mag2flux (nh,uvmag(iuv,1)-uvmagerr(iuv,1),'xu ',FluxU_uv(iuv,1),frequency_uv(iuv,1))
               !write(*,*) 'test1',FluxU_uv(iuv,1)
               CALL  mag2flux (nh,uvmag(iuv,2)-uvmagerr(iuv,2),'xb ',FluxU_uv(iuv,2),frequency_uv(iuv,2))
               CALL  mag2flux (nh,uvmag(iuv,3)-uvmagerr(iuv,3),'xv ',FluxU_uv(iuv,3),frequency_uv(iuv,3))
               !write(*,*) 'test2',FluxU_uv(iuv,1)
               CALL  mag2flux (nh,uvmag(iuv,4)-uvmagerr(iuv,4),'xw1',FluxU_uv(iuv,4),frequency_uv(iuv,4))
               CALL  mag2flux (nh,uvmag(iuv,5)-uvmagerr(iuv,5),'xm2',FluxU_uv(iuv,5),frequency_uv(iuv,5))
               CALL  mag2flux (nh,uvmag(iuv,6)-uvmagerr(iuv,6),'xw2',FluxU_uv(iuv,6),frequency_uv(iuv,6))
               CALL  mag2flux (nh,uvmag(iuv,1)+uvmagerr(iuv,1),'xu ',FluxL_uv(iuv,1),frequency_uv(iuv,1))
               CALL  mag2flux (nh,uvmag(iuv,2)+uvmagerr(iuv,2),'xb ',FluxL_uv(iuv,2),frequency_uv(iuv,2))
               CALL  mag2flux (nh,uvmag(iuv,3)+uvmagerr(iuv,3),'xv ',FluxL_uv(iuv,3),frequency_uv(iuv,3))
               CALL  mag2flux (nh,uvmag(iuv,4)+uvmagerr(iuv,4),'xw1',FluxL_uv(iuv,4),frequency_uv(iuv,4))
               CALL  mag2flux (nh,uvmag(iuv,5)+uvmagerr(iuv,5),'xm2',FluxL_uv(iuv,5),frequency_uv(iuv,5))
               CALL  mag2flux (nh,uvmag(iuv,6)+uvmagerr(iuv,6),'xw2',FluxL_uv(iuv,6),frequency_uv(iuv,6))
               uv_type(iuv)='XMMOM'
               !write(*,*) flux_uv(iuv,1),FluxL_uv(iuv,1)
            ENDIF
            !write(*,*) catalog,FluxU_uv(iuv,6),flux_uv(iuv,6),FluxL_uv(iuv,6),frequency_uv(iuv,6)
         else if ((catalog(1:7) == 'xrtspec') .or. (catalog(1:6) == 'bat105') .or.
     &    (catalog(1:4) == 'oulc') .or. (catalog(1:8) == 'bepposax') .or. (catalog(1:8) == 'nublazar')
     &    .or. (catalog(1:5) == 'ousxg') .or. (catalog(1:5) == 'ousxb')) then
            ixray=ixray+1
            if (ixray .ne. 1) THEN
               do j=1,ixray-1
                  if ((ra_xray(j) .eq. ra) .and. (dec_xray(j) .eq. dec)) THEN
                     if ((name_x(j) == catalog ) .and. (filen_x(j) .ne. filen)) then
c                        write(*,'(4x,"The counterpart",i4,3x,2(f9.5,2x),6(es10.3,2x))')
c     &                        j,ra_xray(j),dec_xray(j),flux_xray(j,1:2)
                        is=ie
                        ie=index(string(is+1:len(string)),' ')+is
                        read(string(is+1:ie-1),'(a)') repflux
c                        write(*,'(a,i4,"_",a,2x,2(f9.5,2x),a)') ' Repeated one:',
c     &                         filen,catalog,ra,dec,repflux(1:lenact(repflux))
                        ixray=ixray-1
                        goto 300
                     endif
                  endif
               enddo
            endif
            ra_xray(ixray)=ra
            dec_xray(ixray)=dec
            name_x(ixray)=catalog
            filen_x(ixray)=filen
            if ((catalog(1:7) == 'xrtspec') .or. (catalog(1:8) == 'nublazar')) then
               if (catalog(1:7) == 'xrtspec') then
                  ixrtsp=ixrtsp+1
                  xrtspind(ixray)=ixrtsp
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) frequency_xray(ixray,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is !!!frequency error
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xray(ixray,1)
               !if (flux_xray(ixray) .le. 0.) write(*,*) ixray,frequency_xray(ixray),flux_xray(ixray)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xray(ixray,1)
               if (catalog(1:7) == 'xrtspec') then
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdst_xrt(ixray)
               endif
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) mjded_xrt(ixray)
               if (catalog(1:8) == 'nublazar') mjdst_xrt(ixray)=mjded_xrt(ixray)
               FluxU_xray(ixray,1)=flux_xray(ixray,1)+Ferr_xray(ixray,1)
               FluxL_xray(ixray,1)=flux_xray(ixray,1)-Ferr_xray(ixray,1)
               if ((flux_xray(ixray,1) .lt. 0.) .or. (FluxL_xray(ixray,1) .lt. 0.)) then
                  flux_xray(ixray,1)=0.
                  FluxU_xray(ixray,1)=0.
                  FluxL_xray(ixray,1)=0.
               endif
               if (catalog(1:7) == 'xrtspec') then
                  poserr_xray(ixray)=5.
                  xray_type(ixray)='XRTSPEC'
               else
                  poserr_xray(ixray)=15.
                  xray_type(ixray)='NuBlazar'
               endif
            else if (catalog(1:8) == 'bepposax') then
c               ixrtsp=ixrtsp+1
c               xrtspind(ixray)=ixrtsp
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) frequency_xray(ixray,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xray(ixray,1)
               !if (flux_xray(ixray) .le. 0.) write(*,*) ixray,frequency_xray(ixray),flux_xray(ixray)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xray(ixray,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdst_xrt(ixray)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) mjded_xrt(ixray)
               FluxU_xray(ixray,1)=flux_xray(ixray,1)+Ferr_xray(ixray,1)
               FluxL_xray(ixray,1)=flux_xray(ixray,1)-Ferr_xray(ixray,1)
               if ((flux_xray(ixray,1) .lt. 0.) .or. (FluxL_xray(ixray,1) .lt. 0.)) then
                  flux_xray(ixray,1)=0.
                  FluxU_xray(ixray,1)=0.
                  FluxL_xray(ixray,1)=0.
               endif
               poserr_xray(ixray)=10.
               xray_type(ixray)='BEPPOSAX'
            else if (catalog(1:4) == 'oulc') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xray(ixray,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xray(ixray,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xray(ixray,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xray(ixray,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdst_xrt(ixray)
               mjded_xrt(ixray)=mjdst_xrt(ixray)
               FluxU_xray(ixray,1)=flux_xray(ixray,1)+Ferr_xray(ixray,1)
               FluxL_xray(ixray,1)=flux_xray(ixray,1)-Ferr_xray(ixray,1)
               if ((flux_xray(ixray,1) .lt. 0.) .or. (FluxL_xray(ixray,1) .lt. 0.) .or. (Ferr_xray(ixray,1) .lt. 0. )) then
                  flux_xray(ixray,1)=0.
                  FluxU_xray(ixray,1)=0.
                  FluxL_xray(ixray,1)=0.
               endif
               call fluxtofdens2(slope_xray(ixray),0.5,2.,flux_xray(ixray,1),1.,fdens,nudens)
               frequency_xray(ixray,1)=nudens
               flux_xray(ixray,1)=fdens
               call fluxtofdens2(slope_xray(ixray),0.5,2.,FluxU_xray(ixray,1),1.,fdens,nudens)
               FluxU_xray(ixray,1)=fdens
               call fluxtofdens2(slope_xray(ixray),0.5,2.,FluxL_xray(ixray,1),1.,fdens,nudens)
               FluxL_xray(ixray,1)=fdens
               FluxU_xray(ixray,2)=flux_xray(ixray,2)+Ferr_xray(ixray,2)
               FluxL_xray(ixray,2)=flux_xray(ixray,2)-Ferr_xray(ixray,2)
               if ((flux_xray(ixray,2) .lt. 0.) .or. (FluxL_xray(ixray,2) .lt. 0.) .or. (Ferr_xray(ixray,2) .lt. 0. )) then
                  flux_xray(ixray,2)=0.
                  FluxU_xray(ixray,2)=0.
                  FluxL_xray(ixray,2)=0.
               endif
               call fluxtofdens2(slope_xray(ixray),2.,10.,flux_xray(ixray,2),4.5,fdens,nudens)
               frequency_xray(ixray,2)=nudens
               flux_xray(ixray,2)=fdens
               call fluxtofdens2(slope_xray(ixray),2.,10.,FluxU_xray(ixray,2),4.5,fdens,nudens)
               FluxU_xray(ixray,2)=fdens
               call fluxtofdens2(slope_xray(ixray),2.,10.,FluxL_xray(ixray,2),4.5,fdens,nudens)
               FluxL_xray(ixray,2)=fdens
               poserr_xray(ixray)=5.
               xray_type(ixray)='OULC'
            else if ((catalog(1:5) == 'ousxb') .or. (catalog(1:5) == 'ousxg')) then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdst_xrt(ixray)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) mjded_xrt(ixray)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xray(ixray,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xray(ixray,2)
               FluxU_xray(ixray,2)=flux_xray(ixray,2)+Ferr_xray(ixray,2)
               FluxL_xray(ixray,2)=flux_xray(ixray,2)-Ferr_xray(ixray,2)
               frequency_xray(ixray,2)=(1.602E-19)*(5.e2)/(6.626e-34)
               if ((Ferr_xray(ixray,2) .lt. 0) .or. (FluxL_xray(ixray,2) .lt. 0)) then
                  FluxU_xray(ixray,2)=Ferr_xray(ixray,2)
                  flux_xray(ixray,2)=0.
                  FluxL_xray(ixray,2)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xray(ixray,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xray(ixray,1)
               FluxU_xray(ixray,1)=flux_xray(ixray,1)+Ferr_xray(ixray,1)
               FluxL_xray(ixray,1)=flux_xray(ixray,1)-Ferr_xray(ixray,1)
               frequency_xray(ixray,1)=(1.602E-19)*(1.e3)/(6.626e-34)
               if ((Ferr_xray(ixray,1) .lt. 0) .or. (FluxL_xray(ixray,1) .lt. 0)) then
                  FluxU_xray(ixray,1)=Ferr_xray(ixray,1)
                  flux_xray(ixray,1)=0.
                  FluxL_xray(ixray,1)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xray(ixray,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xray(ixray,3)
               FluxU_xray(ixray,3)=flux_xray(ixray,3)+Ferr_xray(ixray,3)
               FluxL_xray(ixray,3)=flux_xray(ixray,3)-Ferr_xray(ixray,3)
               frequency_xray(ixray,3)=(1.602E-19)*(1.5e3)/(6.626e-34)
               if ((Ferr_xray(ixray,3) .lt. 0) .or. (FluxL_xray(ixray,3) .lt. 0)) then
                  FluxU_xray(ixray,3)=Ferr_xray(ixray,3)
                  flux_xray(ixray,3)=0.
                  FluxL_xray(ixray,3)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xray(ixray,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xray(ixray,4)
               FluxU_xray(ixray,4)=flux_xray(ixray,4)+Ferr_xray(ixray,4)
               FluxL_xray(ixray,4)=flux_xray(ixray,4)-Ferr_xray(ixray,4)
               frequency_xray(ixray,4)=(1.602E-19)*(3.e3)/(6.626e-34)
               if ((Ferr_xray(ixray,4) .lt. 0) .or. (FluxL_xray(ixray,4) .lt. 0)) then
                  FluxU_xray(ixray,4)=Ferr_xray(ixray,4)
                  flux_xray(ixray,4)=0.
                  FluxL_xray(ixray,4)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xray(ixray,5)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xray(ixray,5)
               FluxU_xray(ixray,5)=flux_xray(ixray,5)+Ferr_xray(ixray,5)
               FluxL_xray(ixray,5)=flux_xray(ixray,5)-Ferr_xray(ixray,5)
               frequency_xray(ixray,5)=(1.602E-19)*(4.5e3)/(6.626e-34)
               if ((Ferr_xray(ixray,5) .lt. 0) .or. (FluxL_xray(ixray,5) .lt. 0)) then
                  FluxU_xray(ixray,5)=flux_xray(ixray,5)
                  flux_xray(ixray,5)=0.
                  FluxL_xray(ixray,5)=0.
               endif
               if (flux_xray(ixray,1) .eq. 0.) then
                   if (flux_xray(ixray,3) .ne. 0.) then
                      flux_xray(ixray,1)=(flux_xray(ixray,3)/frequency_xray(ixray,3))*(1./1.5)**(-0.9)
                      FluxU_xray(ixray,1)=(FluxU_xray(ixray,3)/frequency_xray(ixray,3))*(1./1.5)**(-0.9)
                      FluxL_xray(ixray,1)=(FluxL_xray(ixray,3)/frequency_xray(ixray,3))*(1./1.5)**(-0.9)
                   else
                      flux_xray(ixray,1)=(flux_xray(ixray,4)/frequency_xray(ixray,4))*(1./3.)**(-0.9)
                      FluxU_xray(ixray,1)=(FluxU_xray(ixray,4)/frequency_xray(ixray,4))*(1./3.)**(-0.9)
                      FluxL_xray(ixray,1)=(FluxL_xray(ixray,4)/frequency_xray(ixray,4))*(1./3.)**(-0.9)
                   endif
                  flux_xray(ixray,1)=flux_xray(ixray,1)*frequency_xray(ixray,1)
                  FluxU_xray(ixray,1)=FluxU_xray(ixray,1)*frequency_xray(ixray,1)
                  FluxL_xray(ixray,1)=FluxL_xray(ixray,1)*frequency_xray(ixray,1)
               endif
               if (catalog(1:5) == 'ousxb') then
                  xray_type(ixray)='OUSXB'
               else
                  xray_type(ixray)='OUSXG'
               endif
               poserr_xray(ixray)=7.
               !write(*,*) 'xrt'
            else if (catalog(1:6) == 'bat105') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_xray(ixray,1) !SNR
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_xray(ixray,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) FluxL_xray(ixray,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) FluxU_xray(ixray,1)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) slope_xray(ixray)
               slope_xray(ixray)=slope_xray(ixray)-1
               poserr_xray(ixray)=2*60.*sqrt(0.01+(30.5/Ferr_xray(ixray,1))**2)
               call fluxtofdens2(slope_xray(ixray),14.,195.,flux_xray(ixray,1)*1.e-12,100.,fdens,nudens)
               frequency_xray(ixray,2)=nudens
               flux_xray(ixray,2)=fdens
               call fluxtofdens2(slope_xray(ixray),14.,195.,FluxU_xray(ixray,1)*1.e-12,100.,fdens,nudens)
               FluxU_xray(ixray,2)=fdens
               call fluxtofdens2(slope_xray(ixray),14.,195.,FluxL_xray(ixray,1)*1.e-12,100.,fdens,nudens)
               FluxL_xray(ixray,2)=fdens
               call fluxtofdens2(slope_xray(ixray),14.,195.,flux_xray(ixray,1)*1.e-12,50.,fdens,nudens)
               frequency_xray(ixray,1)=nudens
               flux_xray(ixray,1)=fdens
               call fluxtofdens2(slope_xray(ixray),14.,195.,FluxU_xray(ixray,1)*1.e-12,50.,fdens,nudens)
               FluxU_xray(ixray,1)=fdens
               call fluxtofdens2(slope_xray(ixray),14.,195.,FluxL_xray(ixray,1)*1.e-12,50.,fdens,nudens)
               FluxL_xray(ixray,1)=fdens
               xray_type(ixray)='BAT105m'
               !write(*,*) 'BAT',flux_xray(ixray,1),FluxU_xray(ixray,1),FluxL_xray(ixray,1),poserr_xray(ixray)
            endif
         ELSE IF ((catalog(1:4) == '2fhl') .or. (catalog(1:7) == '4fgldr3') .or.
     &      (catalog(1:4) == '3fgl') .or. (catalog(1:4) == '3fhl') .or. (catalog(1:5) == '2bigb')
     &      .or. (catalog(1:4) == 'fmev') .or. (catalog(1:5) == '2agile')) then
            igam=igam+1
            If (igam > arrsize(7)) stop 'Too many Gamma-ray points'
            if (igam .ne. 1) THEN
               do j=1,igam-1
                  if ((ra_gam(j) .eq. ra) .and. (dec_gam(j) .eq. dec)) THEN
                     if ((name_g(j) == catalog ) .and. (filen_g(j) .ne. filen)) then
                     write(*,'(4x,"The counterpart",i4,3x,2(f9.5,2x),6(es10.3,2x))')
     &                  j,ra_gam(j),dec_gam(j),flux_gam(j,1:6)
                     is=ie
                     ie=index(string(is+1:len(string)),' ')+is
                     read(string(is+1:ie-1),'(a)') repflux
                     write(*,'(a,i4,"_",a,2x,2(f9.5,2x),a)') ' Repeated one:',
     &                   filen,catalog,ra,dec,repflux(1:lenact(repflux))
                     igam=igam-1
                     goto 300
                     endif
                  endif
               enddo
            endif
            ra_gam(igam)=ra
            dec_gam(igam)=dec
            name_g(igam)=catalog
            filen_g(igam)=filen
            
            If (catalog(1:4) == '2fhl') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_gam(igam)
               poserr_gam(igam)=poserr_gam(igam)*3600.*2.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) slope_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) specerr_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) slope_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) specerr_gam(igam,1)
               if (slope_gam(igam,1) .eq. 0.) slope_gam(igam,1)=slope_gam(igam,2)
               if (specerr_gam(igam,1) .eq. 0.) specerr_gam(igam,1)=specerr_gam(igam,2)
               if (slope_gam(igam,2) .gt. 5.) slope_gam(igam,2)=2.
               !write(*,*) 'Gamma slope',slope_gam(igam,1),slope_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,1) !!!!EBL deabsorption
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,1)
               FluxU_gam(igam,1)=flux_gam(igam,1)+Ferr_gam(igam,1)
               FluxL_gam(igam,1)=flux_gam(igam,1)-Ferr_gam(igam,1)
               if ((FluxL_gam(igam,1) .lt. 0.) .or. (FluxL_gam(igam,1) .eq. flux_gam(igam,1))) then
                  FluxU_gam(igam,1)=Ferr_gam(igam,1)*3.
                  FluxL_gam(igam,1)=0.
               endif
               if  (redshift .gt. 0) then
                  idebl=idebl+1
                  eblnn(igam)=idebl
                  call ebl_flux(0.05,2.,slope_gam(igam,2),redshift,reduction_factor)
                  flux_debl(idebl,1)=flux_gam(igam,1)*reduction_factor
                  FluxU_debl(idebl,1)=FluxU_gam(igam,1)*reduction_factor
                  FluxL_debl(idebl,1)=FluxL_gam(igam,1)*reduction_factor
               endif
c               write(*,*) FluxU_gam(igam,1),Flux_gam(igam,1),FluxL_gam(igam,1),slope_gam(igam,2)
               call fluxtofdens(slope_gam(igam,2),50.,2000.,flux_gam(igam,1),50.,fdens,nudens)
               flux_gam(igam,1)=fdens
               frequency_gam(igam,1)=nudens
               call fluxtofdens(slope_gam(igam,2),50.,2000.,FluxU_gam(igam,1),50.,fdens,nudens)
               FluxU_gam(igam,1)=fdens
               call fluxtofdens(slope_gam(igam,2),50.,2000.,FluxL_gam(igam,1),50.,fdens,nudens)
               FluxL_gam(igam,1)=fdens
               if (redshift .gt. 0) then
               if ((FluxL_gam(igam,1) .eq. 0.) .and. (FluxU_gam(igam,1) .gt. 0.)) then
                  flux_debl(idebl,1)=0
                  FluxU_debl(idebl,1)=0
                  FluxL_debl(idebl,1)=0
               else
                  call fluxtofdens(slope_gam(igam,2),50.,2000.,flux_debl(idebl,1),50.,fdens,nudens)
                  flux_debl(idebl,1)=fdens
                  frequency_debl(idebl,1)=nudens
                  call fluxtofdens(slope_gam(igam,2),50.,2000.,FluxU_debl(idebl,1),50.,fdens,nudens)
                  FluxU_debl(idebl,1)=fdens
                  call fluxtofdens(slope_gam(igam,2),50.,2000.,FluxL_debl(idebl,1),50.,fdens,nudens)
                  FluxL_debl(idebl,1)=fdens
               endif
               endif
c               write(*,*) FluxU_gam(igam,1),Flux_gam(igam,1),FluxL_gam(igam,1),slope_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,2)
               FluxL_gam(igam,2)=flux_gam(igam,2)+Ferr_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,2)
               FluxU_gam(igam,2)=flux_gam(igam,2)+Ferr_gam(igam,2)
               if ((FluxL_gam(igam,2) .lt. 0.) .or. (FluxL_gam(igam,2) .eq. flux_gam(igam,2))) then
                  FluxU_gam(igam,2)=Ferr_gam(igam,2)*3.
                  FluxL_gam(igam,2)=0.
               endif
               if  (redshift .gt. 0) then
                  call ebl_flux(0.05,0.171,slope_gam(igam,2),redshift,reduction_factor)
                  flux_debl(idebl,2)=flux_gam(igam,2)*reduction_factor
                  FluxU_debl(idebl,2)=FluxU_gam(igam,2)*reduction_factor
                  FluxL_debl(idebl,2)=FluxL_gam(igam,2)*reduction_factor
               endif
               call fluxtofdens(slope_gam(igam,2),50.,171.,flux_gam(igam,2),100.,fdens,nudens)
               flux_gam(igam,2)=fdens
               frequency_gam(igam,2)=nudens
               call fluxtofdens(slope_gam(igam,2),50.,171.,FluxU_gam(igam,2),100.,fdens,nudens)
               FluxU_gam(igam,2)=fdens
               call fluxtofdens(slope_gam(igam,2),50.,171.,FluxL_gam(igam,2),100.,fdens,nudens)
               FluxL_gam(igam,2)=fdens
               if (redshift .gt. 0) then
               if ((FluxL_gam(igam,2) .eq. 0.) .and. (FluxU_gam(igam,2) .gt. 0.)) then
                  flux_debl(idebl,2)=0
                  FluxU_debl(idebl,2)=0
                  FluxL_debl(idebl,2)=0
               else
                  call fluxtofdens(slope_gam(igam,2),50.,171.,flux_debl(idebl,2),100.,fdens,nudens)
                  flux_debl(idebl,2)=fdens
                  frequency_debl(idebl,2)=nudens
                  call fluxtofdens(slope_gam(igam,2),50.,171.,FluxU_debl(idebl,2),100.,fdens,nudens)
                  FluxU_debl(idebl,2)=fdens
                  call fluxtofdens(slope_gam(igam,2),50.,171.,FluxL_debl(idebl,2),100.,fdens,nudens)
                  FluxL_debl(idebl,2)=fdens
               endif
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,3)
               FluxL_gam(igam,3)=flux_gam(igam,3)+Ferr_gam(igam,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,3)
               FluxU_gam(igam,3)=flux_gam(igam,3)+Ferr_gam(igam,3)
               if ((FluxL_gam(igam,3) .lt. 0.) .or. (FluxL_gam(igam,3) .eq. flux_gam(igam,3))) then
                  FluxU_gam(igam,3)=Ferr_gam(igam,3)*3.
                  FluxL_gam(igam,3)=0.
               endif
               if  (redshift .gt. 0) then
                  call ebl_flux(0.171,0.585,slope_gam(igam,2),redshift,reduction_factor)
                  flux_debl(idebl,3)=flux_gam(igam,3)*reduction_factor
                  FluxU_debl(idebl,3)=FluxU_gam(igam,3)*reduction_factor
                  FluxL_debl(idebl,3)=FluxL_gam(igam,3)*reduction_factor
               endif
               call fluxtofdens(slope_gam(igam,2),171.,585.,flux_gam(igam,3),300.,fdens,nudens)
               flux_gam(igam,3)=fdens
               frequency_gam(igam,3)=nudens
               call fluxtofdens(slope_gam(igam,2),171.,585.,FluxU_gam(igam,3),300.,fdens,nudens)
               FluxU_gam(igam,3)=fdens
               call fluxtofdens(slope_gam(igam,2),171.,585.,FluxL_gam(igam,3),300.,fdens,nudens)
               FluxL_gam(igam,3)=fdens
               if (redshift .gt. 0) then
               if ((FluxL_gam(igam,3) .eq. 0.) .and. (FluxU_gam(igam,3) .gt. 0.)) then
                  flux_debl(idebl,3)=0
                  FluxU_debl(idebl,3)=0
                  FluxL_debl(idebl,3)=0
               else
                  call fluxtofdens(slope_gam(igam,2),171.,585.,flux_debl(idebl,3),300.,fdens,nudens)
                  flux_debl(idebl,3)=fdens
                  frequency_debl(idebl,3)=nudens
                  call fluxtofdens(slope_gam(igam,2),171.,585.,FluxU_debl(idebl,3),300.,fdens,nudens)
                  FluxU_debl(idebl,3)=fdens
                  call fluxtofdens(slope_gam(igam,2),171.,585.,FluxL_debl(idebl,3),300.,fdens,nudens)
                  FluxL_debl(idebl,3)=fdens
               endif
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,4)
               FluxL_gam(igam,4)=flux_gam(igam,4)+Ferr_gam(igam,4)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,4)
               FluxU_gam(igam,4)=flux_gam(igam,4)+Ferr_gam(igam,4)
               if ((FluxL_gam(igam,4) .lt. 0.) .or. (FluxL_gam(igam,4) .eq. flux_gam(igam,4))) then
                  FluxU_gam(igam,4)=Ferr_gam(igam,4)*3.
                  FluxL_gam(igam,4)=0.
               endif
               if  (redshift .gt. 0) then
                  call ebl_flux(0.585,2.,slope_gam(igam,2),redshift,reduction_factor)
                  flux_debl(idebl,4)=flux_gam(igam,4)*reduction_factor
                  FluxU_debl(idebl,4)=FluxU_gam(igam,4)*reduction_factor
                  FluxL_debl(idebl,4)=FluxL_gam(igam,4)*reduction_factor
               endif
               call fluxtofdens(slope_gam(igam,2),585.,2000.,flux_gam(igam,4),1000.,fdens,nudens)
               flux_gam(igam,4)=fdens
               call fluxtofdens(slope_gam(igam,2),585.,2000.,FluxU_gam(igam,4),1000.,fdens,nudens)
               FluxU_gam(igam,4)=fdens
               call fluxtofdens(slope_gam(igam,2),585.,2000.,FluxL_gam(igam,4),1000.,fdens,nudens)
               FluxL_gam(igam,4)=fdens
               frequency_gam(igam,4)=nudens
               if (redshift .gt. 0) then
               if ((FluxL_gam(igam,4) .eq. 0.) .and. (FluxU_gam(igam,4) .gt. 0.)) then
                  flux_debl(idebl,4)=0
                  FluxU_debl(idebl,4)=0
                  FluxL_debl(idebl,4)=0
               else
                  call fluxtofdens(slope_gam(igam,2),585.,2000.,flux_debl(idebl,4),1000.,fdens,nudens)
                  flux_debl(idebl,4)=fdens
                  frequency_debl(idebl,4)=nudens
                  call fluxtofdens(slope_gam(igam,2),585.,2000.,FluxU_debl(idebl,4),1000.,fdens,nudens)
                  FluxU_debl(idebl,4)=fdens
                  call fluxtofdens(slope_gam(igam,2),585.,2000.,FluxL_debl(idebl,4),1000.,fdens,nudens)
                  FluxL_debl(idebl,4)=fdens
               endif
               endif
               gam_type(igam)='2FHL'
               !write(*,*) flux_gam(igam,1),flux_gam(igam,2),flux_gam(igam,3),flux_gam(igam,4)
            ELSE IF (catalog(1:4) == '3fgl') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) slope_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) specerr_gam(igam,1)
               if (slope_gam(igam,1) .gt. 5.) slope_gam(igam,1)=2.
               FluxL_gam(igam,1)=flux_gam(igam,1)-Ferr_gam(igam,1)
               FluxU_gam(igam,1)=flux_gam(igam,1)+Ferr_gam(igam,1)
               if (Ferr_gam(igam,1) .gt. flux_gam(igam,1)) then
                  FluxU_gam(igam,1)=Ferr_gam(igam,1)*3.
                  FluxL_gam(igam,1)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),1.,100.,flux_gam(igam,1),1.,fdens,nudens)
               flux_gam(igam,1)=fdens
               frequency_gam(igam,1)=nudens
               call fluxtofdens(slope_gam(igam,1),1.,100.,FluxU_gam(igam,1),1.,fdens,nudens)
               FluxU_gam(igam,1)=fdens
               call fluxtofdens(slope_gam(igam,1),1.,100.,FluxL_gam(igam,1),1.,fdens,nudens)
               FluxL_gam(igam,1)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,4)
               FluxL_gam(igam,4)=flux_gam(igam,4)+Ferr_gam(igam,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,4)
               FluxU_gam(igam,4)=flux_gam(igam,4)+Ferr_gam(igam,4)
               if ((FluxL_gam(igam,4) .lt. 0.) .or. (FluxL_gam(igam,4) .eq. flux_gam(igam,4))) then
                  FluxU_gam(igam,4)=Ferr_gam(igam,4)*3.
                  FluxL_gam(igam,4)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),0.3,1.,flux_gam(igam,4),0.6,fdens,nudens) !flux at 0.6 GeV
               flux_gam(igam,4)=fdens
               frequency_gam(igam,4)=nudens
               call fluxtofdens(slope_gam(igam,1),0.3,1.,FluxU_gam(igam,4),0.6,fdens,nudens) !flux at 0.6 GeV
               FluxU_gam(igam,4)=fdens
               call fluxtofdens(slope_gam(igam,1),0.3,1.,FluxL_gam(igam,4),0.6,fdens,nudens) !flux at 0.6 GeV
               FluxL_gam(igam,4)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,3)
               FluxL_gam(igam,3)=flux_gam(igam,3)+Ferr_gam(igam,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,3)
               FluxU_gam(igam,3)=flux_gam(igam,3)+Ferr_gam(igam,3)
               if ((FluxL_gam(igam,3) .lt. 0.) .or. (FluxL_gam(igam,3) .eq. flux_gam(igam,3))) then
                  FluxU_gam(igam,3)=Ferr_gam(igam,3)*3.
                  FluxL_gam(igam,3)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,flux_gam(igam,3),0.2,fdens,nudens) !flux at 0.2 GeV
               flux_gam(igam,3)=fdens
               frequency_gam(igam,3)=nudens
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,FluxU_gam(igam,3),0.2,fdens,nudens) !flux at 0.2 GeV
               FluxU_gam(igam,3)=fdens
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,FluxL_gam(igam,3),0.2,fdens,nudens) !flux at 0.2 GeV
               FluxL_gam(igam,3)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,7)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,7)
               FluxL_gam(igam,7)=flux_gam(igam,7)+Ferr_gam(igam,7)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,7)
               FluxU_gam(igam,7)=flux_gam(igam,7)+Ferr_gam(igam,7)
               if ((FluxL_gam(igam,7) .lt. 0.) .or. (FluxL_gam(igam,7) .eq. flux_gam(igam,7))) then
                  FluxU_gam(igam,7)=Ferr_gam(igam,7)*3.
                  FluxL_gam(igam,7)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),10.,100.,flux_gam(igam,7),60.,fdens,nudens)
               flux_gam(igam,7)=fdens
               frequency_gam(igam,7)=nudens
               call fluxtofdens(slope_gam(igam,1),10.,100.,FluxU_gam(igam,7),60.,fdens,nudens)
               FluxU_gam(igam,7)=fdens
               call fluxtofdens(slope_gam(igam,1),10.,100.,FluxL_gam(igam,7),60.,fdens,nudens)
               FluxL_gam(igam,7)=fdens
               if (redshift .gt. 0.) then ! convert the debl flux
                  idebl=idebl+1
                  eblnn(igam)=idebl
                  call ebl_flux(0.06,0.06,slope_gam(igam,1),redshift,reduction_factor)
                  if ((FluxL_gam(igam,7) .eq. 0.) .and. (FluxU_gam(igam,7) .gt. 0.)) then
                     flux_debl(idebl,1)=0. !correction of EBL
                     FluxU_debl(idebl,1)=0.
                     FluxL_debl(idebl,1)=0.
                  else
                     frequency_debl(idebl,1)=frequency_gam(igam,7)
                     flux_debl(idebl,1)=flux_gam(igam,7)*reduction_factor !correction of EBL
                     FluxU_debl(idebl,1)=FluxU_gam(igam,7)*reduction_factor
                     FluxL_debl(idebl,1)=FluxL_gam(igam,7)*reduction_factor
                  endif
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,5)
               !write(*,*) 'test',Ferr_gam(igam,5)
               FluxL_gam(igam,5)=flux_gam(igam,5)+Ferr_gam(igam,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,5)
               FluxU_gam(igam,5)=flux_gam(igam,5)+Ferr_gam(igam,5)
               if ((FluxL_gam(igam,5) .lt. 0.) .or. (FluxL_gam(igam,5) .eq. flux_gam(igam,5))) then
                  FluxU_gam(igam,5)=Ferr_gam(igam,5)*3.
                  FluxL_gam(igam,5)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),1.,3.,flux_gam(igam,5),2.,fdens,nudens)
               flux_gam(igam,5)=fdens
               frequency_gam(igam,5)=nudens
               call fluxtofdens(slope_gam(igam,1),1.,3.,FluxU_gam(igam,5),2.,fdens,nudens)
               FluxU_gam(igam,5)=fdens
               call fluxtofdens(slope_gam(igam,1),1.,3.,FluxL_gam(igam,5),2.,fdens,nudens)
               FluxL_gam(igam,5)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,2)
               FluxL_gam(igam,2)=flux_gam(igam,2)+Ferr_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,2)
               FluxU_gam(igam,2)=flux_gam(igam,2)+Ferr_gam(igam,2)
               if ((FluxL_gam(igam,2) .lt. 0.) .or. (FluxL_gam(igam,2) .eq. flux_gam(igam,2))) then
                  FluxU_gam(igam,2)=Ferr_gam(igam,2)*3.
                  FluxL_gam(igam,2)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),0.03,0.1,flux_gam(igam,2),0.06,fdens,nudens)
               flux_gam(igam,2)=fdens
               frequency_gam(igam,2)=nudens
               call fluxtofdens(slope_gam(igam,1),0.03,0.1,FluxU_gam(igam,2),0.06,fdens,nudens)
               FluxU_gam(igam,2)=fdens
               call fluxtofdens(slope_gam(igam,1),0.03,0.1,FluxL_gam(igam,2),0.06,fdens,nudens)
               FluxL_gam(igam,2)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,6)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,6)
               FluxL_gam(igam,6)=flux_gam(igam,6)+Ferr_gam(igam,6)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,6)
               FluxU_gam(igam,6)=flux_gam(igam,6)+Ferr_gam(igam,6)
               if ((FluxL_gam(igam,6) .lt. 0.) .or. (FluxL_gam(igam,6) .eq. flux_gam(igam,6))) then
                  FluxU_gam(igam,6)=Ferr_gam(igam,6)*3.
                  FluxL_gam(igam,6)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),3.,10.,flux_gam(igam,6),6.,fdens,nudens)
               flux_gam(igam,6)=fdens
               frequency_gam(igam,6)=nudens
               call fluxtofdens(slope_gam(igam,1),3.,10.,FluxU_gam(igam,6),6.,fdens,nudens)
               FluxU_gam(igam,6)=fdens
               call fluxtofdens(slope_gam(igam,1),3.,10.,FluxL_gam(igam,6),6.,fdens,nudens)
               FluxL_gam(igam,6)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) posang
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) major
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) minor
               posxerr=sqrt(((sin(posang)*major)**2)+((cos(posang)*minor)**2))
               posyerr=sqrt(((cos(posang)*major)**2)+((sin(posang)*minor)**2))
               poserr_gam(igam)=max(posxerr,posyerr)*3600.
               gam_type(igam)='3FGL'
            ELSE IF (catalog(1:7) == '4fgldr3') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) major
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) minor
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) posang
               posxerr=sqrt(((sin(posang)*major)**2)+((cos(posang)*minor)**2))
               posyerr=sqrt(((cos(posang)*major)**2)+((sin(posang)*minor)**2))
               poserr_gam(igam)=max(posxerr,posyerr)*3600.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) slope_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) specerr_gam(igam,1)
               if (slope_gam(igam,1) .gt. 5.) slope_gam(igam,1)=2.
               FluxU_gam(igam,1)=flux_gam(igam,1)+Ferr_gam(igam,1)
               FluxL_gam(igam,1)=flux_gam(igam,1)-Ferr_gam(igam,1)
               call fluxtofdens(slope_gam(igam,1),1.,100.,flux_gam(igam,1),1.,fdens,nudens)
               flux_gam(igam,1)=fdens
               frequency_gam(igam,1)=nudens
               call fluxtofdens(slope_gam(igam,1),1.,100.,FluxU_gam(igam,1),1.,fdens,nudens)
               FluxU_gam(igam,1)=fdens
               call fluxtofdens(slope_gam(igam,1),1.,100.,FluxL_gam(igam,1),1.,fdens,nudens)
               FluxL_gam(igam,1)=fdens
               !write(*,*) 'Fermi',poserr_gam(igam),FluxU_gam(igam,1),flux_gam(igam,1),FluxL_gam(igam,1)
               !write(*,*) 'Fermi',poserr_gam(igam),FluxU_gam(igam,2),flux_gam(igam,2),FluxL_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,6)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,7)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,8)
c               write(*,*) flux_gam(igam,7)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,9)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,2)
               FluxL_gam(igam,2)=flux_gam(igam,2)+Ferr_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,2)
               FluxU_gam(igam,2)=flux_gam(igam,2)+Ferr_gam(igam,2)
               if ((FluxL_gam(igam,2) .lt. 0.) .or. (FluxL_gam(igam,2) .eq. flux_gam(igam,2))) then
                  FluxU_gam(igam,2)=Ferr_gam(igam,2)*3.
                  FluxL_gam(igam,2)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),0.05,0.1,flux_gam(igam,2),0.075,fdens,nudens)
               flux_gam(igam,2)=fdens
               frequency_gam(igam,2)=nudens
               call fluxtofdens(slope_gam(igam,1),0.05,0.1,FluxU_gam(igam,2),0.075,fdens,nudens)
               FluxU_gam(igam,2)=fdens
               call fluxtofdens(slope_gam(igam,1),0.05,0.1,FluxL_gam(igam,2),0.075,fdens,nudens)
               FluxL_gam(igam,2)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,3)
               FluxL_gam(igam,3)=flux_gam(igam,3)+Ferr_gam(igam,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,3)
               FluxU_gam(igam,3)=flux_gam(igam,3)+Ferr_gam(igam,3)
               if ((FluxL_gam(igam,3) .lt. 0.) .or. (FluxL_gam(igam,3) .eq. flux_gam(igam,3))) then
                  FluxU_gam(igam,3)=Ferr_gam(igam,3)*3.
                  FluxL_gam(igam,3)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,flux_gam(igam,3),0.2,fdens,nudens)
               flux_gam(igam,3)=fdens
               frequency_gam(igam,3)=nudens
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,FluxU_gam(igam,3),0.2,fdens,nudens)
               FluxU_gam(igam,3)=fdens
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,FluxL_gam(igam,3),0.2,fdens,nudens)
               FluxL_gam(igam,3)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,4)
               FluxL_gam(igam,4)=flux_gam(igam,4)+Ferr_gam(igam,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,4)
               FluxU_gam(igam,4)=flux_gam(igam,4)+Ferr_gam(igam,4)
               if ((FluxL_gam(igam,4) .lt. 0.) .or. (FluxL_gam(igam,4) .eq. flux_gam(igam,4))) then
                  FluxU_gam(igam,4)=Ferr_gam(igam,4)*3.
                  FluxL_gam(igam,4)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),0.3,1.,flux_gam(igam,4),0.6,fdens,nudens)
               flux_gam(igam,4)=fdens
               frequency_gam(igam,4)=nudens
               call fluxtofdens(slope_gam(igam,1),0.3,1.,FluxU_gam(igam,4),0.6,fdens,nudens)
               FluxU_gam(igam,4)=fdens
               call fluxtofdens(slope_gam(igam,1),0.3,1.,FluxL_gam(igam,4),0.6,fdens,nudens)
               FluxL_gam(igam,4)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,5)
               FluxL_gam(igam,5)=flux_gam(igam,5)+Ferr_gam(igam,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,5)
               FluxU_gam(igam,5)=flux_gam(igam,5)+Ferr_gam(igam,5)
               if ((FluxL_gam(igam,5) .lt. 0.) .or. (FluxL_gam(igam,5) .eq. flux_gam(igam,5))) then
               FluxU_gam(igam,5)=Ferr_gam(igam,5)*3.
               FluxL_gam(igam,5)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),1.,3.,flux_gam(igam,5),2.,fdens,nudens)
               flux_gam(igam,5)=fdens
               frequency_gam(igam,5)=nudens
               call fluxtofdens(slope_gam(igam,1),1.,3.,FluxU_gam(igam,5),2.,fdens,nudens)
               FluxU_gam(igam,5)=fdens
               call fluxtofdens(slope_gam(igam,1),1.,3.,FluxL_gam(igam,5),2.,fdens,nudens)
               FluxL_gam(igam,5)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,6)
               FluxL_gam(igam,6)=flux_gam(igam,6)+Ferr_gam(igam,6)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,6)
               FluxU_gam(igam,6)=flux_gam(igam,6)+Ferr_gam(igam,6)
               if ((FluxL_gam(igam,6) .lt. 0.) .or. (FluxL_gam(igam,6) .eq. flux_gam(igam,6))) then
                  FluxU_gam(igam,6)=Ferr_gam(igam,6)*3.
                  FluxL_gam(igam,6)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),3.,10.,flux_gam(igam,6),6.,fdens,nudens)
               flux_gam(igam,6)=fdens
               frequency_gam(igam,6)=nudens
               call fluxtofdens(slope_gam(igam,1),3.,10.,FluxU_gam(igam,6),6.,fdens,nudens)
               FluxU_gam(igam,6)=fdens
               call fluxtofdens(slope_gam(igam,1),3.,10.,FluxL_gam(igam,6),6.,fdens,nudens)
               FluxL_gam(igam,6)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,7)
               FluxL_gam(igam,7)=flux_gam(igam,7)+Ferr_gam(igam,7)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,7)
               FluxU_gam(igam,7)=flux_gam(igam,7)+Ferr_gam(igam,7)
               if ((FluxL_gam(igam,7) .lt. 0.) .or. (FluxL_gam(igam,7) .eq. flux_gam(igam,7))) then
                  FluxU_gam(igam,7)=Ferr_gam(igam,7)*3.
                  FluxL_gam(igam,7)=0.
               endif
               call fluxtofdens(slope_gam(igam,1),10.,30.,flux_gam(igam,7),20.,fdens,nudens)
               flux_gam(igam,7)=fdens
               frequency_gam(igam,7)=nudens
               call fluxtofdens(slope_gam(igam,1),10.,30.,FluxU_gam(igam,7),20.,fdens,nudens)
               FluxU_gam(igam,7)=fdens
               call fluxtofdens(slope_gam(igam,1),10.,30.,FluxL_gam(igam,7),20.,fdens,nudens)
               FluxL_gam(igam,7)=fdens
               if (redshift .gt. 0.) then ! convert the debl flux
                  idebl=idebl+1
                  eblnn(igam)=idebl
                  call ebl_flux(0.02,0.02,slope_gam(igam,1),redshift,reduction_factor)
                  if ((FluxL_gam(igam,7) .eq. 0.) .and. (FluxU_gam(igam,7) .gt. 0.)) then
                     frequency_debl(idebl,1)=frequency_gam(igam,7)
                     flux_debl(idebl,1)=0. !correction of EBL
                     FluxU_debl(idebl,1)=0.
                     FluxL_debl(idebl,1)=0.
                  else
                     frequency_debl(idebl,1)=frequency_gam(igam,7)
                     flux_debl(idebl,1)=flux_gam(igam,7)*reduction_factor !correction of EBL
                     FluxU_debl(idebl,1)=FluxU_gam(igam,7)*reduction_factor
                     FluxL_debl(idebl,1)=FluxL_gam(igam,7)*reduction_factor
                  endif
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,8)
               FluxL_gam(igam,8)=flux_gam(igam,8)+Ferr_gam(igam,8)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,8)
               FluxU_gam(igam,8)=flux_gam(igam,8)+Ferr_gam(igam,8)
               if ((FluxL_gam(igam,8) .lt. 0.) .or. (FluxL_gam(igam,8) .eq. flux_gam(igam,8))) then
                 FluxU_gam(igam,8)=Ferr_gam(igam,8)*3.
                 FluxL_gam(igam,8)=0.
               endif
               if  (redshift .gt. 0) then
                  call ebl_flux(0.03,0.1,slope_gam(igam,1),redshift,reduction_factor)
                  flux_debl(idebl,2)=flux_gam(igam,8)*reduction_factor
                  FluxU_debl(idebl,2)=FluxU_gam(igam,8)*reduction_factor
                  FluxL_debl(idebl,2)=FluxL_gam(igam,8)*reduction_factor
               endif
               call fluxtofdens(slope_gam(igam,1),30.,100.,flux_gam(igam,8),60.,fdens,nudens)
               flux_gam(igam,8)=fdens
               frequency_gam(igam,8)=nudens
               call fluxtofdens(slope_gam(igam,1),30.,100.,FluxU_gam(igam,8),60.,fdens,nudens)
               FluxU_gam(igam,8)=fdens
               call fluxtofdens(slope_gam(igam,1),30.,100.,FluxL_gam(igam,8),60.,fdens,nudens)
               FluxL_gam(igam,8)=fdens
               if (redshift .gt. 0) then
               if ((FluxL_gam(igam,8) .eq. 0.) .and. (FluxU_gam(igam,8) .gt. 0.)) then
                  flux_debl(idebl,2)=0
                  FluxU_debl(idebl,2)=0
                  FluxL_debl(idebl,2)=0
               else
                  call fluxtofdens(slope_gam(igam,1),30.,100.,flux_debl(idebl,2),60.,fdens,nudens)
                  flux_debl(idebl,2)=fdens
                  frequency_debl(idebl,2)=nudens
                  call fluxtofdens(slope_gam(igam,1),30.,100.,FluxU_debl(idebl,2),60.,fdens,nudens)
                  FluxU_debl(idebl,2)=fdens
                  call fluxtofdens(slope_gam(igam,1),30.,100.,FluxL_debl(idebl,2),60.,fdens,nudens)
                  FluxL_debl(idebl,2)=fdens
               endif
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,9)
               FluxL_gam(igam,9)=flux_gam(igam,9)+Ferr_gam(igam,9)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,9)
               FluxU_gam(igam,9)=flux_gam(igam,9)+Ferr_gam(igam,9)
               if ((FluxL_gam(igam,9) .lt. 0.) .or. (FluxL_gam(igam,9) .eq. flux_gam(igam,9))) then
                 FluxU_gam(igam,9)=Ferr_gam(igam,9)*3.
                 FluxL_gam(igam,9)=0.
               endif
               if  (redshift .gt. 0) then
                  call ebl_flux(0.1,1.,slope_gam(igam,1),redshift,reduction_factor)
                  flux_debl(idebl,3)=flux_gam(igam,9)*reduction_factor
                  FluxU_debl(idebl,3)=FluxU_gam(igam,9)*reduction_factor
                  FluxL_debl(idebl,3)=FluxL_gam(igam,9)*reduction_factor
               endif
               call fluxtofdens(slope_gam(igam,1),100.,1000.,flux_gam(igam,9),600.,fdens,nudens)
               flux_gam(igam,9)=fdens
               frequency_gam(igam,9)=nudens
               call fluxtofdens(slope_gam(igam,1),100.,1000.,FluxU_gam(igam,9),600.,fdens,nudens)
               FluxU_gam(igam,9)=fdens
               call fluxtofdens(slope_gam(igam,1),100.,1000.,FluxL_gam(igam,9),600.,fdens,nudens)
               FluxL_gam(igam,9)=fdens
               if (redshift .gt. 0) then
               if ((FluxL_gam(igam,9) .eq. 0.) .and. (FluxU_gam(igam,9) .gt. 0.)) then
                  flux_debl(idebl,3)=0
                  FluxU_debl(idebl,3)=0
                  FluxL_debl(idebl,3)=0
               else
                  call fluxtofdens(slope_gam(igam,1),100.,1000.,flux_debl(idebl,3),600.,fdens,nudens)
                  flux_debl(idebl,3)=fdens
                  frequency_debl(idebl,3)=nudens
                  call fluxtofdens(slope_gam(igam,1),100.,1000.,FluxU_debl(idebl,3),600.,fdens,nudens)
                  FluxU_debl(idebl,3)=fdens
                  call fluxtofdens(slope_gam(igam,1),100.,1000.,FluxL_debl(idebl,3),600.,fdens,nudens)
                  FluxL_debl(idebl,3)=fdens
               endif
               endif

               i4fgl=i4fgl+1
               ra_gamslp(i4fgl)=ra_gam(igam)
               dec_gamslp(i4fgl)=dec_gam(igam)
               poserr_gamslp(i4fgl)=poserr_gam(igam)
               gamslp(i4fgl)=slope_gam(igam,1)
               gam_type(igam)='4FGL-DR3'
            ELSE IF (catalog(1:4) == '3fhl') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_gam(igam)
               poserr_gam(igam)=poserr_gam(igam)*3600.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) slope_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) specerr_gam(igam,1)
               if (slope_gam(igam,1) .gt. 5.) slope_gam(igam,1)=2.
               FluxU_gam(igam,1)=flux_gam(igam,1)+Ferr_gam(igam,1)
               FluxL_gam(igam,1)=flux_gam(igam,1)-Ferr_gam(igam,1)
               if ((FluxL_gam(igam,1) .lt. 0.) .or. (FluxL_gam(igam,1) .eq. flux_gam(igam,1))) then
                  FluxU_gam(igam,1)=Ferr_gam(igam,1)*3.
                  FluxL_gam(igam,1)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,2)
               FluxL_gam(igam,2)=flux_gam(igam,2)+Ferr_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,2)
               FluxU_gam(igam,2)=flux_gam(igam,2)+Ferr_gam(igam,2)
               if ((FluxL_gam(igam,2) .lt. 0.) .or. (FluxL_gam(igam,2) .eq. flux_gam(igam,2))) then
                  FluxU_gam(igam,2)=Ferr_gam(igam,2)*3.
                  FluxL_gam(igam,2)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,3)
               FluxL_gam(igam,3)=flux_gam(igam,3)+Ferr_gam(igam,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,3)
               FluxU_gam(igam,3)=flux_gam(igam,3)+Ferr_gam(igam,3)
               if ((FluxL_gam(igam,3) .lt. 0.) .or. (FluxL_gam(igam,3) .eq. flux_gam(igam,3))) then
                  FluxU_gam(igam,3)=Ferr_gam(igam,3)*3.
                  FluxL_gam(igam,3)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,4)
               FluxL_gam(igam,4)=flux_gam(igam,4)+Ferr_gam(igam,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,4)
               FluxU_gam(igam,4)=flux_gam(igam,4)+Ferr_gam(igam,4)
               if ((FluxL_gam(igam,4) .lt. 0.) .or. (FluxL_gam(igam,4) .eq. flux_gam(igam,4))) then
                  FluxU_gam(igam,4)=Ferr_gam(igam,4)*3.
                  FluxL_gam(igam,4)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,5)
               FluxL_gam(igam,5)=flux_gam(igam,5)+Ferr_gam(igam,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,5)
               FluxU_gam(igam,5)=flux_gam(igam,5)+Ferr_gam(igam,5)
               if ((FluxL_gam(igam,5) .lt. 0.) .or. (FluxL_gam(igam,5) .eq. flux_gam(igam,5))) then
                  FluxU_gam(igam,5)=Ferr_gam(igam,5)*3.
                  FluxL_gam(igam,5)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,6)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,6)
               FluxL_gam(igam,6)=flux_gam(igam,6)+Ferr_gam(igam,6)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,6)
               FluxU_gam(igam,6)=flux_gam(igam,6)+Ferr_gam(igam,6)
               if ((FluxL_gam(igam,6) .lt. 0.) .or. (FluxL_gam(igam,6) .eq. flux_gam(igam,6))) then
                  FluxU_gam(igam,6)=Ferr_gam(igam,6)*3.
                  FluxL_gam(igam,6)=0.
               endif
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) zzfermi
               !write(*,*) '4fgl=',redshift,'input=',zzinput,'million=',zsource(2)
               if (redshift .eq. 0) then
                  if (zzfermi .gt. 0) redshift=zzfermi
               endif
c               if (zsource(ns) .gt. 0.) redshift=zsource(ns)
c               if (zzinput .gt. 0.) redshift=zzinput
               !write(*,*) redshift,idebl
c               if (((FluxL_gam(igam,1) .eq. 0.) .and. (FluxU_gam(igam,1) .gt. 0.)) .or. (redshift .eq. 0.)) then
c                  idebl=idebl+1
c                  call ebl_flux(0.01,1.,slope_gam(igam,1),redshift,reduction_factor)
c                  flux_debl(idebl,1)=flux_gam(igam,1)*reduction_factor !correction of EBL
c                  FluxU_debl(idebl,1)=FluxU_gam(igam,1)*reduction_factor
c                  FluxL_debl(idebl,1)=FluxL_gam(igam,1)*reduction_factor
c               endif
               call fluxtofdens(slope_gam(igam,1),10.,1000.,flux_gam(igam,1),50.,fdens,nudens) !convert the flux
               flux_gam(igam,1)=fdens
               frequency_gam(igam,1)=nudens
               call fluxtofdens(slope_gam(igam,1),10.,1000.,FluxU_gam(igam,1),50.,fdens,nudens)
               FluxU_gam(igam,1)=fdens
               call fluxtofdens(slope_gam(igam,1),10.,1000.,FluxL_gam(igam,1),50.,fdens,nudens)
               FluxL_gam(igam,1)=fdens
               if (redshift .gt. 0.) then ! convert the debl flux
                  idebl=idebl+1
                  eblnn(igam)=idebl
                  call ebl_flux(0.05,0.05,slope_gam(igam,1),redshift,reduction_factor)
                  if ((FluxL_gam(igam,1) .eq. 0.) .and. (FluxU_gam(igam,1) .gt. 0.)) then
                     flux_debl(idebl,1)=0. !correction of EBL
                     FluxU_debl(idebl,1)=0.
                     FluxL_debl(idebl,1)=0.
                  else
                     frequency_debl(idebl,1)=frequency_gam(igam,1)
                     flux_debl(idebl,1)=flux_gam(igam,1)*reduction_factor !correction of EBL
                     FluxU_debl(idebl,1)=FluxU_gam(igam,1)*reduction_factor
                     FluxL_debl(idebl,1)=FluxL_gam(igam,1)*reduction_factor
                  endif
c                  call fluxtofdens(slope_gam(igam,1),10.,1000.,flux_debl(idebl,1),50.,fdens,nudens)
c                  flux_debl(idebl,1)=fdens
c                  frequency_debl(idebl,1)=nudens
c                  call fluxtofdens(slope_gam(igam,1),10.,1000.,FluxU_debl(idebl,1),50.,fdens,nudens)
c                  FluxU_debl(idebl,1)=fdens
c                  call fluxtofdens(slope_gam(igam,1),10.,1000.,FluxL_debl(idebl,1),50.,fdens,nudens)
c                  FluxL_debl(idebl,1)=fdens
               endif
               call fluxtofdens(slope_gam(igam,1),10.,20.,flux_gam(igam,2),15.,fdens,nudens)
               flux_gam(igam,2)=fdens
               frequency_gam(igam,2)=nudens
               call fluxtofdens(slope_gam(igam,1),10.,20.,FluxU_gam(igam,2),15.,fdens,nudens)
               FluxU_gam(igam,2)=fdens
               call fluxtofdens(slope_gam(igam,1),10.,20.,FluxL_gam(igam,2),15.,fdens,nudens)
               FluxL_gam(igam,2)=fdens
               if (redshift .gt. 0) then
                  call ebl_flux(0.02,0.05,slope_gam(igam,1),redshift,reduction_factor)
                  flux_debl(idebl,2)=flux_gam(igam,3)*reduction_factor
                  FluxU_debl(idebl,2)=FluxU_gam(igam,3)*reduction_factor
                  FluxL_debl(idebl,2)=FluxL_gam(igam,3)*reduction_factor
               endif
               call fluxtofdens(slope_gam(igam,1),20.,50.,flux_gam(igam,3),35.,fdens,nudens)
               flux_gam(igam,3)=fdens
               frequency_gam(igam,3)=nudens
               call fluxtofdens(slope_gam(igam,1),20.,50.,FluxU_gam(igam,3),35.,fdens,nudens)
               FluxU_gam(igam,3)=fdens
               call fluxtofdens(slope_gam(igam,1),20.,50.,FluxL_gam(igam,3),35.,fdens,nudens)
               FluxL_gam(igam,3)=fdens
               if (redshift .gt. 0) then
               if ((FluxL_gam(igam,3) .eq. 0.) .and. (FluxU_gam(igam,3) .gt. 0.)) then
                  flux_debl(idebl,2)=0
                  FluxU_debl(idebl,2)=0
                  FluxL_debl(idebl,2)=0
               else
                  call fluxtofdens(slope_gam(igam,1),20.,50.,flux_debl(idebl,2),35.,fdens,nudens)
                  flux_debl(idebl,2)=fdens
                  frequency_debl(idebl,2)=nudens
                  call fluxtofdens(slope_gam(igam,1),20.,50.,FluxU_debl(idebl,2),35.,fdens,nudens)
                  FluxU_debl(idebl,2)=fdens
                  call fluxtofdens(slope_gam(igam,1),20.,50.,FluxL_debl(idebl,2),35.,fdens,nudens)
                  FluxL_debl(idebl,2)=fdens
               endif
               endif
               if (redshift .gt. 0) then
                  call ebl_flux(0.05,0.15,slope_gam(igam,1),redshift,reduction_factor)
                  flux_debl(idebl,3)=flux_gam(igam,4)*reduction_factor
                  FluxU_debl(idebl,3)=FluxU_gam(igam,4)*reduction_factor
                  FluxL_debl(idebl,3)=FluxL_gam(igam,4)*reduction_factor
               endif
               call fluxtofdens(slope_gam(igam,1),50.,150.,flux_gam(igam,4),100.,fdens,nudens)
               flux_gam(igam,4)=fdens
               frequency_gam(igam,4)=nudens
               call fluxtofdens(slope_gam(igam,1),50.,150.,FluxU_gam(igam,4),100.,fdens,nudens)
               FluxU_gam(igam,4)=fdens
               call fluxtofdens(slope_gam(igam,1),50.,150.,FluxL_gam(igam,4),100.,fdens,nudens)
               FluxL_gam(igam,4)=fdens
               if (redshift .gt. 0) then
               if ((FluxL_gam(igam,4) .eq. 0.) .and. (FluxU_gam(igam,4) .gt. 0.)) then
                  flux_debl(idebl,3)=0
                  FluxU_debl(idebl,3)=0
                  FluxL_debl(idebl,3)=0
               else
                  call fluxtofdens(slope_gam(igam,1),50.,150.,flux_debl(idebl,3),100.,fdens,nudens)
                  flux_debl(idebl,3)=fdens
                  frequency_debl(idebl,3)=nudens
                  call fluxtofdens(slope_gam(igam,1),50.,150.,FluxU_debl(idebl,3),100.,fdens,nudens)
                  FluxU_debl(idebl,3)=fdens
                  call fluxtofdens(slope_gam(igam,1),50.,150.,FluxL_debl(idebl,3),100.,fdens,nudens)
                  FluxL_debl(idebl,3)=fdens
               endif
               endif
               if (redshift .gt. 0) then
                  call ebl_flux(0.15,0.5,slope_gam(igam,1),redshift,reduction_factor)
                  flux_debl(idebl,4)=flux_gam(igam,5)*reduction_factor
                  FluxU_debl(idebl,4)=FluxU_gam(igam,5)*reduction_factor
                  FluxL_debl(idebl,4)=FluxL_gam(igam,5)*reduction_factor
               endif
               call fluxtofdens(slope_gam(igam,1),150.,500.,flux_gam(igam,5),300.,fdens,nudens)
               flux_gam(igam,5)=fdens
               frequency_gam(igam,5)=nudens
               call fluxtofdens(slope_gam(igam,1),150.,500.,FluxU_gam(igam,5),300.,fdens,nudens)
               FluxU_gam(igam,5)=fdens
               call fluxtofdens(slope_gam(igam,1),150.,500.,FluxL_gam(igam,5),300.,fdens,nudens)
               FluxL_gam(igam,5)=fdens
               if (redshift .gt. 0) then
               if ((FluxL_gam(igam,5) .eq. 0.) .and. (FluxU_gam(igam,5) .gt. 0.)) then
                  flux_debl(idebl,4)=0
                  FluxU_debl(idebl,4)=0
                  FluxL_debl(idebl,4)=0
               else
                  call fluxtofdens(slope_gam(igam,1),150.,500.,flux_debl(idebl,4),300.,fdens,nudens)
                  flux_debl(idebl,4)=fdens
                  frequency_debl(idebl,4)=nudens
                  call fluxtofdens(slope_gam(igam,1),150.,500.,FluxU_debl(idebl,4),300.,fdens,nudens)
                  FluxU_debl(idebl,4)=fdens
                  call fluxtofdens(slope_gam(igam,1),150.,500.,FluxL_debl(idebl,4),300.,fdens,nudens)
                  FluxL_debl(idebl,4)=fdens
               endif
               endif
               if (redshift .gt. 0) then
                  call ebl_flux(0.5,2.,slope_gam(igam,1),redshift,reduction_factor)
                  flux_debl(idebl,5)=flux_gam(igam,6)*reduction_factor
                  FluxU_debl(idebl,5)=FluxU_gam(igam,6)*reduction_factor
                  FluxL_debl(idebl,5)=FluxL_gam(igam,6)*reduction_factor
               endif
               call fluxtofdens(slope_gam(igam,1),500.,2000.,flux_gam(igam,6),1000.,fdens,nudens)
               flux_gam(igam,6)=fdens
               frequency_gam(igam,6)=nudens
               call fluxtofdens(slope_gam(igam,1),500.,2000.,FluxU_gam(igam,6),1000.,fdens,nudens)
               FluxU_gam(igam,6)=fdens
               call fluxtofdens(slope_gam(igam,1),500.,2000.,FluxL_gam(igam,6),1000.,fdens,nudens)
               FluxL_gam(igam,6)=fdens
               if (redshift .gt. 0) then
               if ((FluxL_gam(igam,6) .eq. 0.) .and. (FluxU_gam(igam,6) .gt. 0.)) then
                  flux_debl(idebl,5)=0
                  FluxU_debl(idebl,5)=0
                  FluxL_debl(idebl,5)=0
               else
                  call fluxtofdens(slope_gam(igam,1),500.,2000.,flux_debl(idebl,5),1000.,fdens,nudens)
                  flux_debl(idebl,5)=fdens
                  frequency_debl(idebl,5)=nudens
                  call fluxtofdens(slope_gam(igam,1),500.,2000.,FluxU_debl(idebl,5),1000.,fdens,nudens)
                  FluxU_debl(idebl,5)=fdens
                  call fluxtofdens(slope_gam(igam,1),500.,2000.,FluxL_debl(idebl,5),1000.,fdens,nudens)
                  FluxL_debl(idebl,5)=fdens
               endif
               endif
               gam_type(igam)='3FHL'
               if (zsource(ns)+zzinput .eq. 0) redshift=0
            ELSE IF (catalog(1:5) =='2bigb') then
               ibigb=ibigb+1
c               bigbind(igam)=MOD(ibigb,10)
               poserr_gam(igam)=10. !!!set to 10 arcsec
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) frequency_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,1)
               FluxU_gam(igam,1)=flux_gam(igam,1)+Ferr_gam(igam,1)
               FluxL_gam(igam,1)=flux_gam(igam,1)-Ferr_gam(igam,1)
               if ((FluxL_gam(igam,1) .lt. 0.) .or. (FluxL_gam(igam,1) .eq .flux_gam(igam,1))) then
                  FluxU_gam(igam,1)=flux_gam(igam,1)
                  flux_gam(igam,1)=0.
                  FluxL_gam(igam,1)=0.
               endif
               gammatev=frequency_gam(igam,1)/(2.418e26)
               slope_gam(igam,1)=1.8 !!!!!!!!!!!change later
               if ((redshift .gt. 0.) .and. (gammatev .ge. 0.02 )) then ! convert the debl flux
                  idebl=idebl+1
                  eblnn(igam)=idebl
                  call ebl_flux(gammatev,gammatev,slope_gam(igam,1),redshift,reduction_factor)
                  if ((FluxL_gam(igam,1) .eq. flux_gam(igam,1)) .and. (FluxU_gam(igam,1) .ge. flux_gam(igam,1))) then
                     flux_debl(idebl,1)=0. !correction of EBL
                     FluxU_debl(idebl,1)=0.
                     FluxL_debl(idebl,1)=0.
                  else
                     frequency_debl(idebl,1)=frequency_gam(igam,1)
                     flux_debl(idebl,1)=flux_gam(igam,1)*reduction_factor !correction of EBL
                     FluxU_debl(idebl,1)=FluxU_gam(igam,1)*reduction_factor
                     FluxL_debl(idebl,1)=FluxL_gam(igam,1)*reduction_factor
                  endif
               endif
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),'(a)') namegam(igam)
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
               gam_type(igam)='2BIGB'
               !write(*,*) namegam(igam),ibigb,bigbind(igam),poserr_gam(igam)
            ELSE IF (catalog(1:5) =='2agile') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_gam(igam)
               poserr_gam(igam)=sqrt((0.1+poserr_gam(igam)**2))*3600.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) slope_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) specerr_gam(igam,1)
               slope_gam(igam,1)=slope_gam(igam,1)!+1
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,1) !0.1-10 GeV
               FluxU_gam(igam,1)=flux_gam(igam,1)+Ferr_gam(igam,1)
               FluxL_gam(igam,1)=flux_gam(igam,1)-Ferr_gam(igam,1)
               flux_gam(igam,1)=flux_gam(igam,1)*1.e-8
               FluxU_gam(igam,1)=FluxU_gam(igam,1)*1.e-8
               FluxL_gam(igam,1)=FluxL_gam(igam,1)*1.e-8
               call fluxtofdens(slope_gam(igam,1),0.1,10.,flux_gam(igam,1),1.,fdens,nudens)
               flux_gam(igam,1)=fdens
               frequency_gam(igam,1)=nudens
               call fluxtofdens(slope_gam(igam,1),0.1,10.,FluxU_gam(igam,1),1.,fdens,nudens)
               FluxU_gam(igam,1)=fdens
               call fluxtofdens(slope_gam(igam,1),0.1,10.,FluxL_gam(igam,1),1.,fdens,nudens)
               FluxL_gam(igam,1)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,2) !0.1-0.3 GeV
               FluxL_gam(igam,2)=flux_gam(igam,2)-Ferr_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,2) !0.1-0.3 GeV
               FluxU_gam(igam,2)=flux_gam(igam,2)+Ferr_gam(igam,2)
               flux_gam(igam,2)=flux_gam(igam,2)*1.e-8
               FluxU_gam(igam,2)=FluxU_gam(igam,2)*1.e-8
               FluxL_gam(igam,2)=FluxL_gam(igam,2)*1.e-8
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,flux_gam(igam,2),0.2,fdens,nudens)
               flux_gam(igam,2)=fdens
               frequency_gam(igam,2)=nudens
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,FluxU_gam(igam,2),0.2,fdens,nudens)
               FluxU_gam(igam,2)=fdens
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,FluxL_gam(igam,2),0.2,fdens,nudens)
               FluxL_gam(igam,2)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,3) !0.1-0.3 GeV
               FluxL_gam(igam,3)=flux_gam(igam,3)-Ferr_gam(igam,3)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,3) !0.1-0.3 GeV
               FluxU_gam(igam,3)=flux_gam(igam,3)+Ferr_gam(igam,3)
               flux_gam(igam,3)=flux_gam(igam,3)*1.e-8
               FluxU_gam(igam,3)=FluxU_gam(igam,3)*1.e-8
               FluxL_gam(igam,3)=FluxL_gam(igam,3)*1.e-8
               call fluxtofdens(slope_gam(igam,1),0.3,1.,flux_gam(igam,3),0.6,fdens,nudens)
               flux_gam(igam,3)=fdens
               frequency_gam(igam,3)=nudens
               call fluxtofdens(slope_gam(igam,1),0.3,1.,FluxU_gam(igam,3),0.6,fdens,nudens)
               FluxU_gam(igam,3)=fdens
               call fluxtofdens(slope_gam(igam,1),0.3,1.,FluxL_gam(igam,3),0.6,fdens,nudens)
               FluxL_gam(igam,3)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
                if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,4) !0.1-0.3 GeV
               FluxL_gam(igam,4)=flux_gam(igam,4)-Ferr_gam(igam,4)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,4) !0.1-0.3 GeV
               FluxU_gam(igam,4)=flux_gam(igam,4)+Ferr_gam(igam,4)
               flux_gam(igam,4)=flux_gam(igam,4)*1.e-8
               FluxU_gam(igam,4)=FluxU_gam(igam,4)*1.e-8
               FluxL_gam(igam,4)=FluxL_gam(igam,4)*1.e-8
               call fluxtofdens(slope_gam(igam,1),1.,3.,flux_gam(igam,4),2.,fdens,nudens)
               flux_gam(igam,4)=fdens
               frequency_gam(igam,4)=nudens
               call fluxtofdens(slope_gam(igam,1),1.,3.,FluxU_gam(igam,4),2.,fdens,nudens)
               FluxU_gam(igam,4)=fdens
               call fluxtofdens(slope_gam(igam,1),1.,3.,FluxL_gam(igam,4),2.,fdens,nudens)
               FluxL_gam(igam,4)=fdens
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,5)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,5) !0.1-0.3 GeV
               FluxL_gam(igam,5)=flux_gam(igam,5)-Ferr_gam(igam,5)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,5) !0.1-0.3 GeV
               FluxU_gam(igam,5)=flux_gam(igam,5)+Ferr_gam(igam,5)
               flux_gam(igam,5)=flux_gam(igam,5)*1.e-8
               FluxU_gam(igam,5)=FluxU_gam(igam,5)*1.e-8
               FluxL_gam(igam,5)=FluxL_gam(igam,5)*1.e-8
               call fluxtofdens(slope_gam(igam,1),3.,10.,flux_gam(igam,5),6.,fdens,nudens)
               flux_gam(igam,5)=fdens
               frequency_gam(igam,5)=nudens
               call fluxtofdens(slope_gam(igam,1),3.,10.,FluxU_gam(igam,5),6.,fdens,nudens)
               FluxU_gam(igam,5)=fdens
               call fluxtofdens(slope_gam(igam,1),3.,10.,FluxL_gam(igam,5),6.,fdens,nudens)
               FluxL_gam(igam,5)=fdens
               gam_type(igam)='2AGILE'
            ELSE IF (catalog(1:4) =='fmev') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) poserr_gam(igam)
               poserr_gam(igam)=poserr_gam(igam)*3600.
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,1)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,1) !30-100 MeV
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_gam(igam,2)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_gam(igam,2) !100-300 MeV
               slope_gam(igam,1)=2.0 !!!!!!!!!!change later
               FluxU_gam(igam,1)=flux_gam(igam,1)+Ferr_gam(igam,1)
               FluxL_gam(igam,1)=flux_gam(igam,1)-Ferr_gam(igam,1)
               FluxU_gam(igam,2)=flux_gam(igam,2)+Ferr_gam(igam,2)
               FluxL_gam(igam,2)=flux_gam(igam,2)-Ferr_gam(igam,2)
               call fluxtofdens(slope_gam(igam,1),0.03,0.1,flux_gam(igam,1),0.06,fdens,nudens)
               flux_gam(igam,1)=fdens
               frequency_gam(igam,1)=nudens
               call fluxtofdens(slope_gam(igam,1),0.03,0.1,FluxU_gam(igam,1),0.06,fdens,nudens)
               FluxU_gam(igam,1)=fdens
               call fluxtofdens(slope_gam(igam,1),0.03,0.1,FluxL_gam(igam,1),0.06,fdens,nudens)
               FluxL_gam(igam,1)=fdens
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,flux_gam(igam,2),0.2,fdens,nudens)
               flux_gam(igam,2)=fdens
               frequency_gam(igam,2)=nudens
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,FluxU_gam(igam,2),0.2,fdens,nudens)
               FluxU_gam(igam,2)=fdens
               call fluxtofdens(slope_gam(igam,1),0.1,0.3,FluxL_gam(igam,2),0.2,fdens,nudens)
               FluxL_gam(igam,2)=fdens
               gam_type(igam)='FermiMeV'
            ENDIF
            poserr_gam(igam)=poserr_gam(igam)*1.1
            !write(*,*) catalog,FluxU_gam(igam,4),flux_gam(igam,4),FluxL_gam(igam,4)
         ELSE IF ((catalog(1:6) == 'fmonlc') .or. (catalog(1:7) == 'ftaptlc')) THEN
            iflcuv=iflcuv+1
            if (iflcuv .ne. 1) THEN
               do j=1,iflcuv-1
                  if ((ra_flcuv(j) .eq. ra) .and. (dec_flcuv(j) .eq. dec)) THEN
                     if ((name_a(j) == catalog ) .and. (filen_a(j) .ne. filen)) then
c                write(*,'(4x,"The counterpart",i4,3x,2(f9.5,2x),6(es10.3,2x))')
c     &                  j,ra_flcuv(j),dec_flcuv(j),flux_flcuv(j,1:6)
                     is=ie
                     ie=index(string(is+1:len(string)),' ')+is
                     read(string(is+1:ie-1),'(a)') repflux
c                write(*,'(a,i4,"_",a,2x,2(f9.5,2x),a)') ' Repeated one:',
c     &                   filen,catalog,ra,dec,repflux(1:lenact(repflux))
                     iflcuv=iflcuv-1
                     goto 300
                     endif
                  endif
               enddo
            endif
            ra_flcuv(iflcuv)=ra
            dec_flcuv(iflcuv)=dec
            name_a(iflcuv)=catalog
            filen_a(iflcuv)=filen
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdst_flcuv(iflcuv)
            if (catalog(1:6) == 'fmonlc') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) mjded_flcuv(iflcuv)
            else
               mjded_flcuv(iflcuv)=mjdst_flcuv(iflcuv)+30.
            endif
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_flcuv(iflcuv,1)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_flcuv(iflcuv,1)
            FluxU_flcuv(iflcuv,1)=flux_flcuv(iflcuv,1)+Ferr_flcuv(iflcuv,1)
            FluxL_flcuv(iflcuv,1)=flux_flcuv(iflcuv,1)-Ferr_flcuv(iflcuv,1)
            if (catalog(1:6) == 'fmonlc') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) duration(iflcuv)
               if (log10(duration(iflcuv)) .lt. 5.) THEN
                  iflcuv=iflcuv-1
                  goto 300
               endif
            endif
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_flcuv(iflcuv,2)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_flcuv(iflcuv,2)
            FluxU_flcuv(iflcuv,2)=flux_flcuv(iflcuv,2)+Ferr_flcuv(iflcuv,2)
            FluxL_flcuv(iflcuv,2)=flux_flcuv(iflcuv,2)-Ferr_flcuv(iflcuv,2)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_flcuv(iflcuv,3)
            is=ie
            ie=index(string(is+1:len(string)),',')+is
            if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_flcuv(iflcuv,3)
            FluxU_flcuv(iflcuv,3)=flux_flcuv(iflcuv,3)+Ferr_flcuv(iflcuv,3)
            FluxL_flcuv(iflcuv,3)=flux_flcuv(iflcuv,3)-Ferr_flcuv(iflcuv,3)
            if (catalog(1:6) =='fmonlc') then
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) ts(iflcuv)
               flux_flcuv(iflcuv,4)=flux_flcuv(iflcuv,2)
               FluxU_flcuv(iflcuv,4)=FluxU_flcuv(iflcuv,2)
               FluxL_flcuv(iflcuv,4)=FluxL_flcuv(iflcuv,2)
               flcuv_type(iflcuv)='FMonLC'
            else
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_flcuv(iflcuv,4)
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_flcuv(iflcuv,4)
               FluxU_flcuv(iflcuv,4)=flux_flcuv(iflcuv,4)+Ferr_flcuv(iflcuv,4)
               FluxL_flcuv(iflcuv,4)=flux_flcuv(iflcuv,4)-Ferr_flcuv(iflcuv,4)
               flcuv_type(iflcuv)='FTAptLC'
               ts(iflcuv)=25.
            endif
            if ((FluxL_flcuv(iflcuv,1) .le. 0.) .or. (flux_flcuv(iflcuv,1) .le. 0.)
     $        .or. (ts(iflcuv) .lt. 10.)) then
               if (catalog(1:6) == 'fmonlc') FluxU_flcuv(iflcuv,1)=flux_flcuv(iflcuv,1)
               flux_flcuv(iflcuv,1)=0.
               FluxL_flcuv(iflcuv,1)=0.
               if (FluxU_flcuv(iflcuv,1) .le. 0.) FluxU_flcuv(iflcuv,1)=0.
            endif
            if ((FluxL_flcuv(iflcuv,2) .le. 0.) .or. (flux_flcuv(iflcuv,2) .le. 0.)
     $        .or. (ts(iflcuv) .lt. 10.)) then
                if (catalog(1:6) == 'fmonlc') FluxU_flcuv(iflcuv,2)=flux_flcuv(iflcuv,2)
               flux_flcuv(iflcuv,2)=0.
               FluxL_flcuv(iflcuv,2)=0.
               if (FluxU_flcuv(iflcuv,2) .le. 0. ) FluxU_flcuv(iflcuv,2)=0.
            endif
            if ((FluxL_flcuv(iflcuv,3) .le. 0.) .or. (flux_flcuv(iflcuv,3) .le. 0.)
     $        .or. (ts(iflcuv) .lt. 10.)) then
               if (catalog(1:6) == 'fmonlc') FluxU_flcuv(iflcuv,3)=flux_flcuv(iflcuv,3)
              flux_flcuv(iflcuv,3)=0.
              FluxL_flcuv(iflcuv,3)=0.
              if (FluxU_flcuv(iflcuv,3) .le. 0. ) FluxU_flcuv(iflcuv,3)=0.
            endif
            if ((FluxL_flcuv(iflcuv,4) .le. 0.) .or. (flux_flcuv(iflcuv,4) .le. 0.)
     $        .or. (ts(iflcuv) .lt. 10.)) then
               if (catalog(1:6) == 'fmonlc') FluxU_flcuv(iflcuv,4)=flux_flcuv(iflcuv,4)
              flux_flcuv(iflcuv,4)=0.
              FluxL_flcuv(iflcuv,4)=0.
              if (FluxU_flcuv(iflcuv,4) .le. 0. ) FluxU_flcuv(iflcuv,4)=0.
            endif
            if (i4fgl .gt. 0) then
               do s=1,i4fgl
                  call DIST_SKY(ra_gamslp(s),dec_gamslp(s),ra_flcuv(iflcuv),dec_flcuv(iflcuv),dist)
                  !write(*,*) ra_gamslp(s),dec_gamslp(s),ra_flcuv(iflcuv),dec_flcuv(iflcuv),dist*3600.,poserr_gamslp(s)
                  if (dist*3600. .lt. poserr_gamslp(s)*1.1) then
                     slope_flcuv(iflcuv)=gamslp(s)
                  else
                     slope_flcuv(iflcuv)=2.0
                  endif
               enddo
            else
               slope_flcuv(iflcuv)=2.
            endif
            if (catalog(1:7) == 'ftaptlc') then
               poserr_flcuv(iflcuv)=120.
            else
               poserr_flcuv(iflcuv)=12.
            endif
            !write(*,*) flux_flcuv(iflcuv,1:3),slope_flcuv(iflcuv)
            if (catalog(1:6) == 'fmonlc') then
               engmin=0.1
            else
               engmin=0.3
            endif
            !write(*,*) flux_flcuv(iflcuv,1:4),slope_flcuv(iflcuv),i4fgl
            call fluxtofdens(slope_flcuv(iflcuv),engmin,300.,flux_flcuv(iflcuv,1),1.,fdens,nudens)
            flux_flcuv(iflcuv,1)=fdens
            frequency_flcuv(iflcuv,1)=nudens
            call fluxtofdens(slope_flcuv(iflcuv),engmin,300.,FluxU_flcuv(iflcuv,1),1.,fdens,nudens)
            FluxU_flcuv(iflcuv,1)=fdens
            call fluxtofdens(slope_flcuv(iflcuv),engmin,300.,FluxL_flcuv(iflcuv,1),1.,fdens,nudens)
            FluxL_flcuv(iflcuv,1)=fdens
            if (catalog(1:6) == 'fmonlc') then
               fluxind=flux_flcuv(iflcuv,3)
               Ufluxind=FluxU_flcuv(iflcuv,3)
               Lfluxind=FluxL_flcuv(iflcuv,3)
            else
               fluxind=flux_flcuv(iflcuv,2)
               Ufluxind=FluxU_flcuv(iflcuv,2)
               Lfluxind=FluxL_flcuv(iflcuv,2)
            endif
            call fluxtofdens(slope_flcuv(iflcuv),0.3,1.,fluxind,0.6,fdens,nudens)
            fluxind=flux_flcuv(iflcuv,2)
            flux_flcuv(iflcuv,2)=fdens
            frequency_flcuv(iflcuv,2)=nudens
            call fluxtofdens(slope_flcuv(iflcuv),0.3,1.,Ufluxind,0.6,fdens,nudens)
            Ufluxind=FluxU_flcuv(iflcuv,2)
            FluxU_flcuv(iflcuv,2)=fdens
            call fluxtofdens(slope_flcuv(iflcuv),0.3,1.,Lfluxind,0.6,fdens,nudens)
            Lfluxind=FluxL_flcuv(iflcuv,2)
            FluxL_flcuv(iflcuv,2)=fdens
            if (catalog(1:6) == 'fmonlc') then
               engmax=300.
            else
               engmax=30.
               fluxind=flux_flcuv(iflcuv,3)
               Ufluxind=FluxU_flcuv(iflcuv,3)
               Lfluxind=FluxL_flcuv(iflcuv,3)
            endif
            call fluxtofdens(slope_flcuv(iflcuv),1.,engmax,fluxind,10.,fdens,nudens)
            flux_flcuv(iflcuv,3)=fdens
            frequency_flcuv(iflcuv,3)=nudens
            call fluxtofdens(slope_flcuv(iflcuv),1.,engmax,Ufluxind,10.,fdens,nudens)
            FluxU_flcuv(iflcuv,3)=fdens
            call fluxtofdens(slope_flcuv(iflcuv),1.,engmax,Lfluxind,10.,fdens,nudens)
            FluxL_flcuv(iflcuv,3)=fdens
            if (catalog(1:6) == 'fmonlc') then
               engmin=1.
            else
               engmin=30.
            endif
            call fluxtofdens(slope_flcuv(iflcuv),engmin,300.,flux_flcuv(iflcuv,4),100.,fdens,nudens)
            flux_flcuv(iflcuv,4)=fdens
            frequency_flcuv(iflcuv,4)=nudens
            call fluxtofdens(slope_flcuv(iflcuv),engmin,300.,FluxU_flcuv(iflcuv,4),100.,fdens,nudens)
            FluxU_flcuv(iflcuv,4)=fdens
            call fluxtofdens(slope_flcuv(iflcuv),engmin,300.,FluxL_flcuv(iflcuv,4),100.,fdens,nudens)
            FluxL_flcuv(iflcuv,4)=fdens
            !write(*,*) iflcuv,flux_flcuv(iflcuv,1:4)
         ELSE IF ((catalog(1:5) == 'magic') .or. (catalog(1:7) == 'veritas')) then
            ivhe=ivhe+1
            ra_vhe(ivhe)=ra
            dec_vhe(ivhe)=dec
            filen_v(ivhe)=filen
            If (catalog(1:5) == 'magic') then
               imagic=imagic+1
               magicind(ivhe)=imagic
               poserr_vhe(ivhe)=10. !!!set to 10 arcsec
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) frequency_vhe(ivhe)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_vhe(ivhe)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               is=ie
               ie=index(string(is+1:len(string)),' ')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_vhe(ivhe)
               if (Ferr_vhe(ivhe) .eq. -999.) Ferr_vhe(ivhe)=0.
               FluxU_vhe(ivhe)=flux_vhe(ivhe)+Ferr_vhe(ivhe)
               FluxL_vhe(ivhe)=flux_vhe(ivhe)-Ferr_vhe(ivhe)
               if (Ferr_vhe(ivhe) .gt. flux_vhe(ivhe)) then
                  flux_vhe(ivhe)=0.
                  FluxU_vhe(ivhe)=0.
                  FluxL_vhe(ivhe)=0.
               endif
               mjdstart(ivhe)=mjdavg
               mjdend(ivhe)=mjdavg
               gammatev=frequency_vhe(ivhe)/(2.418e26)
               if ((redshift .gt. 0.) .and. (gammatev .ge. 0.02 )) then ! convert the debl flux
                  idebl=idebl+1
                  eblnn(igam+ivhe)=idebl
                  call ebl_flux(gammatev,gammatev,2.2,redshift,reduction_factor)
                  if (((FluxL_vhe(ivhe) .eq. 0.) .and. (FluxU_vhe(ivhe) .eq. 0.)) .or. ((FluxL_vhe(ivhe) .eq. flux_vhe(ivhe)) .and. (FluxU_vhe(ivhe) .eq. flux_vhe(ivhe))))then
                     flux_debl(idebl,1)=0. !correction of EBL
                     FluxU_debl(idebl,1)=0.
                     FluxL_debl(idebl,1)=0.
                  else
                     frequency_debl(idebl,1)=frequency_vhe(ivhe)
                     flux_debl(idebl,1)=flux_vhe(ivhe)*reduction_factor !correction of EBL
                     FluxU_debl(idebl,1)=FluxU_vhe(ivhe)*reduction_factor
                     FluxL_debl(idebl,1)=FluxL_vhe(ivhe)*reduction_factor
                  endif
               endif
               vhe_type(ivhe)='MAGIC'
            else if (catalog(1:7) == 'veritas') then
               iverit=iverit+1
               veritind(ivhe)=iverit
               poserr_vhe(ivhe)=10. !!!set to 10 arcsec
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) frequency_vhe(ivhe)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) flux_vhe(ivhe)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_vhe(ivhe)
               FluxL_vhe(ivhe)=flux_vhe(ivhe)-Ferr_vhe(ivhe)
               FluxU_vhe(ivhe)=flux_vhe(ivhe)+Ferr_vhe(ivhe)
               is=ie
               ie=index(string(is+1:len(string)),',')+is
               if (is .ne. ie-1) read(string(is+1:ie-1),*) Ferr_vhe(ivhe)
               if (Ferr_vhe(ivhe) .lt. 10000.) then
                  FluxU_vhe(ivhe)=flux_vhe(ivhe)+Ferr_vhe(ivhe)
                  is=ie
                  ie=index(string(is+1:len(string)),',')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdstart(ivhe)
                  is=ie
                  ie=index(string(is+1:len(string)),' ')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdend(ivhe)
               else
                  mjdstart(ivhe)=Ferr_vhe(ivhe)
                  is=ie
                  ie=index(string(is+1:len(string)),' ')+is
                  if (is .ne. ie-1) read(string(is+1:ie-1),*) mjdend(ivhe)
               endif
               flux_vhe(ivhe)=flux_vhe(ivhe)*1.E-4*1.602E-19*1.E7*1.E12*frequency_vhe(ivhe)**2
               FluxU_vhe(ivhe)=FluxU_vhe(ivhe)*1.E-4*1.602E-19*1.E7*1.E12*frequency_vhe(ivhe)**2
               FluxL_vhe(ivhe)=FluxL_vhe(ivhe)*1.E-4*1.602E-19*1.E7*1.E12*frequency_vhe(ivhe)**2
               frequency_vhe(ivhe)=(1.602E-19)*(frequency_vhe(ivhe)*1.e12)/(6.626e-34)
               gammatev=frequency_vhe(ivhe)/(2.418e26)
               if ((redshift .gt. 0.) .and. (gammatev .ge. 0.02 )) then ! convert the debl flux
                  idebl=idebl+1
                  eblnn(igam+ivhe)=idebl
                  call ebl_flux(gammatev,gammatev,2.2,redshift,reduction_factor)
                  if (((FluxL_vhe(ivhe) .eq. 0.) .and. (FluxU_vhe(ivhe) .eq. 0.)) .or. ((FluxL_vhe(ivhe) .eq. flux_vhe(ivhe)) .and. (FluxU_vhe(ivhe) .eq. flux_vhe(ivhe))))then
                     flux_debl(idebl,1)=0. !correction of EBL
                     FluxU_debl(idebl,1)=0.
                     FluxL_debl(idebl,1)=0.
                  else
                     frequency_debl(idebl,1)=frequency_vhe(ivhe)
                     flux_debl(idebl,1)=flux_vhe(ivhe)*reduction_factor !correction of EBL
                     FluxU_debl(idebl,1)=FluxU_vhe(ivhe)*reduction_factor
                     FluxL_debl(idebl,1)=FluxL_vhe(ivhe)*reduction_factor
                  endif
               endif
               vhe_type(ivhe)='VERITAS'
            endif
         ENDIF
      ENDDO
 99   CONTINUE
      write(*,*)"     "
      write(*,*) 'Number of candidates each band:',i4p8,ipccs100,ifar,iir,iusno,iuv,ixray,igam,ivhe,iflcuv
      CLOSE (lu_in)
      open(14,file=output_file2,status='unknown',iostat=ier)
c      open(17,file=output_file3,status='unknown',iostat=ier)
      write(*,*)"     "

c      write(*,*) "Check FLC",frequency_flcuv(iflcuv,1),flux_flcuv(iflcuv,1),mjdst_flcuv(iflcuv),iflcuv
c      write(*,*) frequency_flcuv(iflcuv,3),flux_flcuv(iflcuv,3),duration(iflcuv)

      if (aim == 'sed') then
         isource=1
         ra_source(1)=ra_center
         dec_source(1)=dec_center
      endif
c      maxebl=maxval(eblnn)
c      write(*,*) maxebl,idebl

      do i=1,i4p8
         f4p8part(i)=0
         f4p8like(i)=-100.
         do j=1,isource
            CALL DIST_SKY(ra_source(j),dec_source(j),ra_4p8(i),dec_4p8(i),dist)
            if (dist .le. min_dist_4p8) then
            drop=2
            sigma=sqrt(epos(1,j)**2+poserr_4p8(i)**2)
            if (dist*3600. .gt. sigma) drop=5
            call liklihood(epos(1,j),poserr_4p8(i),dist,flux(1,j),drop,like)
            !write(*,*) "5 GHz",i,j,epos(1,j),poserr_4p8(i),dist*3600.,like,(flux(1,j)/(1.4E9*1.E-26))
            if ((like .gt. f4p8like(i)) .and. (2.*poserr_4p8(i) .gt. dist*3600.)) then
               f4p8like(i)=like
               f4p8part(i)=j
            endif
            endif
         enddo
c         if (f4p8part(i) .eq. 0) write(*,*) "Warning!!!Check 4.8 GHz counterpart."
      enddo
c      if (i4p8 .ne. 0) then
c         write(*,*) "5 GHz",f4p8like(1:i4p8)
c         write(*,*) "5 GHz",f4p8part(1:i4p8)
c      endif

      do i=1,ipccs100
         pccspart(i)=0
         pccslike(i)=-100.
         do j=1,isource
            CALL DIST_SKY(ra_source(j),dec_source(j),ra_pccs100(i),dec_pccs100(i),dist)
c            write(*,*) ra_source(j),dec_source(j),ra_pccs100(i),dec_pccs100(i),dist
            if (dist .le. min_dist_pccs100) then
            drop=2
            sigma=sqrt(epos(1,j)**2+poserr_pccs100(i)**2)
            if (dist*3600. .gt. sigma) drop=5
            call liklihood(epos(1,j),poserr_pccs100(i),dist,flux(1,j),drop,like)
            !write(*,*) "100 GHz",i,j,epos(1,j),poserr_pccs100(i),dist*3600.,like,(flux(1,j)/(1.4E9*1.E-26))
            if ((like .gt. pccslike(i)) .and. ( 2.*poserr_pccs100(i) .gt. dist*3600.)) then
               pccslike(i)=like
               pccspart(i)=j
            endif
            endif
         enddo
c         if (pccspart(i) .eq. 0) write(*,*) "Warning!!!Check 100 GHz counterpart."
      enddo
c      if (ipccs100 .ne. 0) then
c         write(*,*) "100 GHz",pccslike(1:ipccs100)
c         write(*,*) "100 GHz",pccspart(1:ipccs100)
c      endif

      do i=1,ifar
         farpart(i)=0
         farlike(i)=-100.
         do j=1,isource
            CALL DIST_SKY(ra_source(j),dec_source(j),ra_far(i),dec_far(i),dist)
            if (dist .le. min_dist_far) then
            drop=2
            sigma=sqrt(poserr_far(i)**2+epos(1,j)**2)
            if (dist*3600. .gt. sigma) drop=5
            call liklihood(epos(1,j),poserr_far(i),dist,flux(1,j),drop,like)
            !write(*,*) "far IR",i,j,epos(1,j),poserr_far(i),dist*3600.,like,(flux(1,j)/(1.4E9*1.E-26))
            if ((like .gt. farlike(i)) .and. (2.*poserr_far(i) .gt. dist*3600.)) then
               farlike(i)=like
               farpart(i)=j
            endif
            endif
         enddo
c         if (farpart(i) .eq. 0) write(*,*) "Warning!!!Check far-IR counterpart."
      enddo
c      if (ifar .ne. 0) then
c         write(*,*) "far IR",farlike(1:ifar)
c         write(*,*) "far IR",farpart(1:ifar)
c      endif

      do i=1,igam
         gamlike(i)=-100.
         gampart(i)=0
         do j=1,isource
            CALL DIST_SKY(ra_source(j),dec_source(j),ra_gam(i),dec_gam(i),dist)
            if (dist .le. min_dist_gam) then
            !write(*,*) epos(1,j),poserr_gam(i)
            drop=2
            sigma=sqrt(poserr_gam(i)**2+epos(1,j)**2)
            if (dist*3600. .gt. sigma) drop=5
            !write(*,*) dist*3600., poserr_gam(i),drop
            call liklihood(epos(1,j),poserr_gam(i),dist,flux(1,j),drop,like)
            !write(*,*) "gamma-ray",i,j,epos(1,j),poserr_gam(i),dist*3600.,like,(flux(1,j)/(1.4E9*1.E-26))
            if ((like .gt. gamlike(i)) .and. ( 2.*poserr_gam(i) .gt. dist*3600.)) then
               gamlike(i)=like
               gampart(i)=j
            endif
            endif
         enddo
c         if (gampart(i) .eq. 0) write(*,*) "Warning!!!Check Gamma-ray counterpar."
      enddo
c      if (igam .ne. 0) then
c         write(*,*) "gamma-ray",gamlike(1:igam)
c         write(*,*) "gamma-ray",gampart(1:igam)
c      endif

      open(lu_output,file=output_file(1:lenact(output_file)),status='unknown',iostat=ier)
      if (ns .ne. 0) then
         sourceu=ns
         sourcel=ns
      else
         sourcel=1
         sourceu=isource
      endif
      DO j=sourcel,sourceu
         flux_x = 0.
         type_average = 0.
         ix = 0
         ir = 0
         flux_r = 0.
         do i=1,npt(j)
            if ((spec_type(i,j) .gt. 50.) .and. (flux(i,j) .ne. 0.)) then
               flux_x=flux_x+flux(i,j)
               ix=ix+1
            endif
            if ((frequency(i,j) .lt. 1.e10) .and. (flux(i,j) .ne. 0.)) then
               flux_r=flux_r+flux(i,j)
               ir=ir+1
            endif
         enddo
         if (ix .ne. 0) flux_x = flux_x/float(ix)
         if (ir .ne. 0) flux_r = flux_r/float(ir)
         if ((flux_x .ne. 0.) .and. (flux_r .ne. 0.)) then
            arx = 1.-log10(flux_x/flux_r)/log10(2.41e17/1.4e9)
         else
            arx=99.99
         endif
         write(*,*)"     "
         write(*,*) "==============================================="
         write(*,'(i4,2x,"candidate with average radio 1.4 GHz flux density ",f9.3," ,",2x,
     &             "average X-ray 1 keV flux ",es10.3)') j,flux_r/(1.4E9*1.E-26),flux_x
         write(*,'(6x,"average X-ray-radio spectral slope: ",f6.3)') arx
         write(*,'("6x,R.A. , Dec. = ", f9.5,2x,",",f9.5)') ra_source(j),dec_source(j)
         write(*,*)"     "
         alphar = 99.99
         a100x = 99.99
         farirx=99.99
         ialphar = 0
         aalphar = 0.
         ir100found=0
         ifarfound=0
         ilowrfound=0

c         write(*,*) '..................Low frequency Radio....................'
         do i=1,ilowr
            call DIST_SKY(ra_source(j),dec_source(j),ra_lowr(i),dec_lowr(i),dist)
            min_dist=sqrt(epos(1,j)**2+poserr_lowr(i)**2)
            IF ( dist*3600. < max(min_dist,2.)) then
               ilowrfound=ilowrfound+1
               flux_lowrcand(ilowrfound)=flux_lowr(i)
               uflux_lowrcand(ilowrfound)=FluxU_lowr(i)
               lflux_lowrcand(ilowrfound)=FluxL_lowr(i)
               freq_lowrcand(ilowrfound)=frequency_lowr(i)
               ra_lowrcand(ilowrfound)=ra_lowr(i)
               dec_lowrcand(ilowrfound)=dec_lowr(i)
               epos_lowrcand(ilowrfound)=poserr_lowr(i)
               lowrdist(ilowrfound)=real(dist)
               lowrcand_type(ilowrfound)=lowr_type(i)
            endif
         enddo
         IF (ilowrfound == 0) then
            write(*,'('' NO WISH object within '',f5.2,'' arcsec'')') max(epos(1,j)*1.3,2.)
         endif

         write(*,*) '........................5 GHz Radio......................'
         DO i=1,i4p8
            CALL DIST_SKY(ra_source(j),dec_source(j),ra_4p8(i),dec_4p8(i),dist)
            !IF (dist < min_dist_4p8) THEN
            !if (f4p8part(i) .eq. j) then
            f4p8part(i)=0
            min_dist=sqrt(poserr_4p8(i)**2+epos(1,j)**2)
            if (dist*3600. .lt. min_dist) then
               f4p8part(i)=j
               ialphar = ialphar +1
               if ((flag_4p8(i,1) == 'E') .or. (flag_4p8(i,1) == 'e')) write(*,*) "Warning!!!!Radio Extended."
               if ((flux_r .ne. 0.) .and. (flux_4p8(i,1) .ne. 0.))
     &             alphar = 1.-log10(flux_r/flux_4p8(i,1))/log10(1.4E9/4.8e9)
               aalphar = aalphar + alphar
               if ((f4p8_type(i) == 'PMN') .or. (f4p8_type(i) == 'GB87') .or. (f4p8_type(i) == 'GB6')) THEN
                  write(*,'(a,"flux density",2x,f9.3,",",2x,f7.3," arcsec away")')
     &               f4p8_type(i),flux_4p8(i,1)/(frequency_4p8(i,1)*1.E-26),dist*3600.
               else if (f4p8_type(i) == 'ATPMN') then
                  write(*,'(a,"flux density (4.8 8.6 GHz)",2(2x,f9.3),",",2x,f7.3," arcsec away")')
     &               f4p8_type(i),flux_4p8(i,1:2)/(frequency_4p8(i,1:2)*1.E-26),dist*3600.
               else if (f4p8_type(i) == 'AT20G') then
                  write(*,'(a,"flux density (5 8 20 GHz)",3(2x,f9.3),",",2x,f7.3," arcsec away")')
     &               f4p8_type(i),flux_4p8(i,1:3)/(frequency_4p8(i,1:3)*1.E-26),dist*3600.
               else if ((f4p8_type(i) == 'KUEHR') .and. (kuehrind(i) .eq. 1)) then
                  write(*,'(a,"flux density (5 8 20 GHz)",3(2x,f9.3),",",2x,f7.3," arcsec away")')
     &               f4p8_type(i),flux_4p8(i,1:3)/(frequency_4p8(i,1:3)*1.E-26),dist*3600.

               ENDIF
               !!!!!write(*,'(a,6x,''Radio alpha: '',f6.3)') f4p8_type(i),alphar
            ENDIF
         ENDDO
         IF (ialphar .ne. 0) then
            alphar = aalphar/float(ialphar)
            !write(*,'(6x,''Radio alpha: '',f5.2,",",2x,f6.3," arcsec away")') alphar,dist*3600
         else
            write(*,'(" No 5 GHz detection within",f5.0,2x,"arcsec")') min_dist_4p8*3600.
         endif

         write(*,*) '.................100 GHz Radio........................'
         DO i=1,ipccs100
            CALL DIST_SKY(ra_source(j),dec_source(j),ra_pccs100(i),dec_pccs100(i),dist)
            pccconv=1.
            !IF (dist < min_dist_pccs100) THEN
            !if (pccspart(i) .eq. j) then
            !write(*,*) poserr_pccs100(i)*1.3, dist*3600.
            pccspart(i)=0
            min_dist=sqrt(poserr_pccs100(i)**2+epos(1,j)**2)
            if (dist*3600. .lt. min_dist) then
               pccspart(i)=j
               ir100found=ir100found+1
               if (pccs100_type(i) == 'PCCS2') then
                  write(*,'(f4.0," GHz flux density",2x,f9.3,",",2x,f7.3," arcmin away")')
     &               frequency_pccs100(i,1)/1.E9,flux_pccs100(i,1)/(frequency_pccs100(i,1)*1.E-26),dist*60.
                  if (frequency_pccs100(ipccs100,1) .gt. 1.E11) pccconv=(100./143.)**(-0.3+1)
                  if ((flux_x .ne. 0.) .and. (flux_pccs100(i,1) .ne. 0.))
     &             a100x = 1.-log10(flux_pccs100(i,1)*pccconv/flux_x)/log10(1.e11/2.418e17)
                  write(*,'(6x,''100 GHz - X-ray slope: '',f6.3)') a100x
               else if (pccs100_type(i) == 'ALMA') then
                  if (almaind(i) .eq. 1) then
                     write(*,'("ALMA",f7.3," arcmin away")') dist*60.
                  endif
               else
                  write(*,'("Planck PCNT",f7.3," arcmin away")')dist*60.
               endif
            ENDIF
        ENDDO
c        write(*,*) ir100found
        if (ir100found .eq. 0) write(*,'(" No 100 GHz detection within",f5.0,2x,"arcmin")') min_dist_pccs100*60.

        write(*,*) '.................Far Infrared........................'
        DO i=1,ifar
        CALL DIST_SKY(ra_source(j),dec_source(j),ra_far(i),dec_far(i),dist)
        !IF (dist < min_dist_pccs100) THEN
           !if (farpart(i) .eq. j) then
            farpart(i)=0
            min_dist=sqrt(poserr_far(i)**2+epos(1,j)**2)
            if (dist*3600. .lt. min_dist) then
              farpart(i)=j
              ifarfound=ifarfound+1
              write(*,'(a," flux density",2x,f9.3,",",2x,f7.3," arcsec away")') far_type(i),
     &               flux_far(i,1)/(frequency_far(i,1)*1.E-26),dist*3600.
              if ((flux_x .ne. 0.) .and. (flux_far(i,1) .ne. 0.))
     &             farirx = 1.-log10(flux_far(i,1)/flux_x)/log10(frequency_far(i,1)/2.418e17)
              write(*,'(6x,''far Infrared - X-ray slope: '',f6.3)') farirx
           ENDIF
         ENDDO
         !write(*,*) ir100found
         if (ifarfound .eq. 0) write(*,'(" No far IR detection within",f5.0,2x,"arcsec")') min_dist_far*3600.

         write(*,*) '.....................IR.............................'
         iofound = 0
         iuvfound = 0
         iirfound = 0
         ii1=0
         ii2=0
         irdistval=10.
c         write(*,*) irlcid(1),irlcid(2)
         do i=1,iir
            call DIST_SKY(ra_source(j),dec_source(j),ra_ir(i),dec_ir(i),dist)
            min_dist=sqrt(epos(1,j)**2+poserr_ir(i)**2)
            if ((ir_type(i) == 'NEOWISE') .and. (irlc_name(i) /= 'nan')) then
               if (dist*3600 .le. irdistval) then
                  irdistval=dist*3600.
                  if (irlc_name(i) /= irlcid(2)) irlcid(2)=irlc_name(i)
               endif
            endif
            if ((ir_type(i) == 'WISEME') .and. (irlc_name(i) /= 'nan')) then
               if (dist*3600 .le. irdistval) then
                  irdistval=dist*3600.
                  if (irlc_name(i) /= irlcid(1)) irlcid(1)=irlc_name(i)
               endif
            endif
            !write(*,*) matchradius,epos(1,j),poserr_ir(i),ir_type(i),dist*3600.
            IF (( dist*3600. < max(min_dist,2.) )  .and. (ir_type(i) == 'WISE')) THEN
               ii1=ii1+1
               iirfound=iirfound+1
               if (ii1 .eq. 1) then
                  flux_ircand(iirfound,1:4)=flux_ir(i,1:4)
                  uflux_ircand(iirfound,1:4)=FluxU_ir(i,1:4)
                  lflux_ircand(iirfound,1:4)=FluxL_ir(i,1:4)
                  irmag_cand(iirfound,1:4)=irmag(i,1:4)
                  freq_ircand(iirfound,1:4)=frequency_ir(i,1:4)
                  ra_ircand(iirfound)=ra_ir(i)
                  dec_ircand(iirfound)=dec_ir(i)
                  epos_ircand(iirfound)=poserr_ir(i)
                  irdist(iirfound)=real(dist)
                  ircand_type(iirfound)=ir_type(i)
               else
                  iirfound=iirfound-1
                  if (irdist(iirfound) .gt. dist) then
                     flux_ircand(iirfound,1:4)=flux_ir(i,1:4)
                     uflux_ircand(iirfound,1:4)=FluxU_ir(i,1:4)
                     lflux_ircand(iirfound,1:4)=FluxL_ir(i,1:4)
                     irmag_cand(iirfound,1:4)=irmag(i,1:4)
                     freq_ircand(iirfound,1:4)=frequency_ir(i,1:4)
                     ra_ircand(iirfound)=ra_ir(i)
                     dec_ircand(iirfound)=dec_ir(i)
                     epos_ircand(iirfound)=poserr_ir(i)
                     irdist(iirfound)=real(dist)
                     ircand_type(iirfound)=ir_type(i)
                  endif
               endif
            else if (( dist*3600. < max(min_dist,2.)  ) .and. (ir_type(i) == '2MASS')) THEN
               iirfound=iirfound+1
               ii2=ii2+1
               if (ii2 .eq. 1) then
                  flux_ircand(iirfound,1:4)=flux_ir(i,1:4)
                  uflux_ircand(iirfound,1:4)=FluxU_ir(i,1:4)
                  lflux_ircand(iirfound,1:4)=FluxL_ir(i,1:4)
                  irmag_cand(iirfound,1:4)=irmag(i,1:4)
                  freq_ircand(iirfound,1:4)=frequency_ir(i,1:4)
                  ra_ircand(iirfound)=ra_ir(i)
                  dec_ircand(iirfound)=dec_ir(i)
                  epos_ircand(iirfound)=poserr_ir(i)
                  irdist(iirfound)=real(dist)
                  ircand_type(iirfound)=ir_type(i)
               else
                  iirfound=iirfound-1
                  if (irdist(iirfound) .gt. dist) then
                     flux_ircand(iirfound,1:4)=flux_ir(i,1:4)
                     uflux_ircand(iirfound,1:4)=FluxU_ir(i,1:4)
                     lflux_ircand(iirfound,1:4)=FluxL_ir(i,1:4)
                     irmag_cand(iirfound,1:4)=irmag(i,1:4)
                     freq_ircand(iirfound,1:4)=frequency_ir(i,1:4)
                     ra_ircand(iirfound)=ra_ir(i)
                     dec_ircand(iirfound)=dec_ir(i)
                     epos_ircand(iirfound)=poserr_ir(i)
                     irdist(iirfound)=real(dist)
                     ircand_type(iirfound)=ir_type(i)
                  endif
               endif
            endif
         enddo
c         write(*,*) iirfound,ii1,ii2
         IF (iirfound == 0) then
            write(*,'('' NO IR object within '',f5.2,'' arcsec'')') max(epos(1,j)*1.3,2.)
         else if (iirlc .gt. 0) then
            iircheck=iirfound
c            write(*,*) iircheck,iirfound
            do i=1,iircheck
               airx=99.99
               arir=99.99
               if ((flux_x .ne. 0.) .and. (flux_ircand(i,2) .ne. 0.))
     &              airx = 1.-log10(flux_ircand(i,2)/flux_x)/log10(freq_ircand(i,2)/2.418e17) !k band, w2 band
               if ((flux_r .ne. 0.) .and. (flux_ircand(i,2) .ne. 0.))
     &               arir = 1.-log10(flux_r/flux_ircand(i,2))/log10(1.4e9/freq_ircand(i,2))
               write(*,'(a,"IR-X-ray slope: ",f6.3,",",2x,"radio-IR slope: ",f6.3,",",2x,f7.3," arcsec away")')
     &             ircand_type(i),airx,arir,irdist(i)*3600.
cccccccc           check candidate has lc or not
c               write(*,*) "IR pts number",iirlc
               if (ircand_type(i) == 'WISE') then
               do s=1,iirlc
c                  if (((ir_type(indirlc(s)) == 'WISEME') .or. (ir_type(indirlc(s)) == 'NEOWISE'))
c     &                .and. (ircand_type(i) == 'WISE')) then
                  call DIST_SKY(ra_ircand(i),dec_ircand(i),ra_ir(indirlc(s)),dec_ir(indirlc(s)),dist)
                     !write(*,*) i,indirlc(s),ra_ircand(i),dec_ircand(i),ra_ir(indirlc(s)),dec_ir(indirlc(s))
                  if (dist*3600. .lt. 2.) then
c                     if ( ((irlc_name(s) == irlcid(1)).and.(ir_type(indirlc(s)) == 'WISEME')) .or.
c     &                   ((irlc_name(s) == irlcid(2)).and.(ir_type(indirlc(s)) == 'NEOWISE'))) then
                     iirfound=iirfound+1
                     flux_ircand(iirfound,1:4)=flux_ir(indirlc(s),1:4)
                     uflux_ircand(iirfound,1:4)=FluxU_ir(indirlc(s),1:4)
                     lflux_ircand(iirfound,1:4)=FluxL_ir(indirlc(s),1:4)
                     irmag_cand(iirfound,1:4)=irmag(indirlc(s),1:4)
                     freq_ircand(iirfound,1:4)=frequency_ir(indirlc(s),1:4)
                     ra_ircand(iirfound)=ra_ir(indirlc(s))
                     dec_ircand(iirfound)=dec_ir(indirlc(s))
                     ircand_type(iirfound)=ir_type(indirlc(s))
                     mjdst_irlccand(iirfound)=mjdst_irlc(indirlc(s))
                     mjded_irlccand(iirfound)=mjded_irlc(indirlc(s))
c                  endif
c                        write(*,*)
                  endif
c                 endif
               enddo
               endif
            enddo
         endif
c         write(*,*) iirfound

         write(*,*) '.......................Optical........................'
         ii1=0
         ii2=0
         ii3=0
         ii4=0
         ii5=0
         DO i=1,iusno
            CALL DIST_SKY(ra_source(j),dec_source(j),ra_usno(i),dec_usno(i),dist)
            min_dist=sqrt(poserr_usno(i)**2+epos(1,j)**2)
            IF (( dist*3600. < max(min_dist,2.)  ) .and. (opt_type(i) == 'USNO'))THEN
               iofound = iofound+1
               ii1=ii1+1
               if (ii1 .eq. 1) then
                  flux_usnocand(iofound,1:5)=flux_usno(i,1:5)
                  uflux_usnocand(iofound,1:5)=FluxU_usno(i,1:5)
                  lflux_usnocand(iofound,1:5)=FluxL_usno(i,1:5)
                  usnomag_cand(iofound,1:5)=usnomag(i,1:5)
                  freq_usnocand(iofound,1:5)=frequency_usno(i,1:5)
                  ra_usnocand(iofound)=ra_usno(i)
                  dec_usnocand(iofound)=dec_usno(i)
                  epos_usnocand(iofound)=poserr_usno(i)
                  optdist(iofound)=real(dist)
                  optcand_type(iofound)=opt_type(i)
               else
                  iofound=iofound-1
                  if (optdist(iofound) .gt. dist) then
                     flux_usnocand(iofound,1:5)=flux_usno(i,1:5)
                     uflux_usnocand(iofound,1:5)=FluxU_usno(i,1:5)
                     lflux_usnocand(iofound,1:5)=FluxL_usno(i,1:5)
                     usnomag_cand(iofound,1:5)=usnomag(i,1:5)
                     freq_usnocand(iofound,1:5)=frequency_usno(i,1:5)
                     ra_usnocand(iofound)=ra_usno(i)
                     dec_usnocand(iofound)=dec_usno(i)
                     epos_usnocand(iofound)=poserr_usno(i)
                     optdist(iofound)=real(dist)
                     optcand_type(iofound)=opt_type(i)
                  endif
               endif
            else IF (( dist*3600. < max(min_dist,2.)  ) .and. (opt_type(i) == 'SDSS'))THEN
               iofound = iofound+1
               ii2=ii2+1
               if (ii2 .eq. 1) then
                  flux_usnocand(iofound,1:5)=flux_usno(i,1:5)
                  uflux_usnocand(iofound,1:5)=FluxU_usno(i,1:5)
                  lflux_usnocand(iofound,1:5)=FluxL_usno(i,1:5)
                  usnomag_cand(iofound,1:5)=usnomag(i,1:5)
                  freq_usnocand(iofound,1:5)=frequency_usno(i,1:5)
                  ra_usnocand(iofound)=ra_usno(i)
                  dec_usnocand(iofound)=dec_usno(i)
                  epos_usnocand(iofound)=poserr_usno(i)
                  optdist(iofound)=real(dist)
                  optcand_type(iofound)=opt_type(i)
               else
                  iofound=iofound-1
                  if (optdist(iofound) .gt. dist) then
                     flux_usnocand(iofound,1:5)=flux_usno(i,1:5)
                     uflux_usnocand(iofound,1:5)=FluxU_usno(i,1:5)
                     lflux_usnocand(iofound,1:5)=FluxL_usno(i,1:5)
                     usnomag_cand(iofound,1:5)=usnomag(i,1:5)
                     freq_usnocand(iofound,1:5)=frequency_usno(i,1:5)
                     ra_usnocand(iofound)=ra_usno(i)
                     dec_usnocand(iofound)=dec_usno(i)
                     epos_usnocand(iofound)=poserr_usno(i)
                     optdist(iofound)=real(dist)
                     optcand_type(iofound)=opt_type(i)
                  endif
               endif
            else IF (( dist*3600. < max(min_dist,2.)  ) .and. (opt_type(i) == 'HSTGSC'))THEN
               iofound = iofound+1
               ii3=ii3+1
               if (ii3 .eq. 1) then
                  flux_usnocand(iofound,1:5)=flux_usno(i,1:5)
                  uflux_usnocand(iofound,1:5)=FluxU_usno(i,1:5)
                  lflux_usnocand(iofound,1:5)=FluxL_usno(i,1:5)
                  usnomag_cand(iofound,1:5)=usnomag(i,1:5)
                  freq_usnocand(iofound,1:5)=frequency_usno(i,1:5)
                  ra_usnocand(iofound)=ra_usno(i)
                  dec_usnocand(iofound)=dec_usno(i)
                  epos_usnocand(iofound)=poserr_usno(i)
                  optdist(iofound)=real(dist)
                  optcand_type(iofound)=opt_type(i)
               else
                  iofound=iofound-1
                  if (optdist(iofound) .gt. dist) then
                     flux_usnocand(iofound,1:5)=flux_usno(i,1:5)
                     uflux_usnocand(iofound,1:5)=FluxU_usno(i,1:5)
                     lflux_usnocand(iofound,1:5)=FluxL_usno(i,1:5)
                     usnomag_cand(iofound,1:5)=usnomag(i,1:5)
                     freq_usnocand(iofound,1:5)=frequency_usno(i,1:5)
                     ra_usnocand(iofound)=ra_usno(i)
                     dec_usnocand(iofound)=dec_usno(i)
                     epos_usnocand(iofound)=poserr_usno(i)
                     optdist(iofound)=real(dist)
                     optcand_type(iofound)=opt_type(i)
                  endif
               endif
            else IF (( dist*3600. < max(min_dist,2.)  ) .and. (opt_type(i) == 'PanSTARRS'))THEN
               iofound = iofound+1
               ii4=ii4+1
               if (ii4 .eq. 1) then
                  flux_usnocand(iofound,1:5)=flux_usno(i,1:5)
                  uflux_usnocand(iofound,1:5)=FluxU_usno(i,1:5)
                  lflux_usnocand(iofound,1:5)=FluxL_usno(i,1:5)
                  usnomag_cand(iofound,1:5)=usnomag(i,1:5)
                  freq_usnocand(iofound,1:5)=frequency_usno(i,1:5)
                  ra_usnocand(iofound)=ra_usno(i)
                  dec_usnocand(iofound)=dec_usno(i)
                  epos_usnocand(iofound)=poserr_usno(i)
                  optdist(iofound)=real(dist)
                  optcand_type(iofound)=opt_type(i)
               else
                  iofound=iofound-1
                  if (optdist(iofound) .gt. dist) then
                     flux_usnocand(iofound,1:5)=flux_usno(i,1:5)
                     uflux_usnocand(iofound,1:5)=FluxU_usno(i,1:5)
                     lflux_usnocand(iofound,1:5)=FluxL_usno(i,1:5)
                     usnomag_cand(iofound,1:5)=usnomag(i,1:5)
                     freq_usnocand(iofound,1:5)=frequency_usno(i,1:5)
                     ra_usnocand(iofound)=ra_usno(i)
                     dec_usnocand(iofound)=dec_usno(i)
                     epos_usnocand(iofound)=poserr_usno(i)
                     optdist(iofound)=real(dist)
                     optcand_type(iofound)=opt_type(i)
                  endif
               endif
            else IF (( dist*3600. < max(min_dist,2.)  ) .and. (opt_type(i) == 'GAIA'))THEN
               iofound = iofound+1
               ii5=ii5+1
               if (ii5 .eq. 1) then
                  flux_usnocand(iofound,1:5)=flux_usno(i,1:5)
                  uflux_usnocand(iofound,1:5)=FluxU_usno(i,1:5)
                  lflux_usnocand(iofound,1:5)=FluxL_usno(i,1:5)
                  usnomag_cand(iofound,1:5)=usnomag(i,1:5)
                  freq_usnocand(iofound,1:5)=frequency_usno(i,1:5)
                  ra_usnocand(iofound)=ra_usno(i)
                  dec_usnocand(iofound)=dec_usno(i)
                  epos_usnocand(iofound)=poserr_usno(i)
                  optdist(iofound)=real(dist)
                  optcand_type(iofound)=opt_type(i)
               else
                  iofound=iofound-1
                  if (optdist(iofound) .gt. dist) then
                     flux_usnocand(iofound,1:5)=flux_usno(i,1:5)
                     uflux_usnocand(iofound,1:5)=FluxU_usno(i,1:5)
                     lflux_usnocand(iofound,1:5)=FluxL_usno(i,1:5)
                     usnomag_cand(iofound,1:5)=usnomag(i,1:5)
                     freq_usnocand(iofound,1:5)=frequency_usno(i,1:5)
                     ra_usnocand(iofound)=ra_usno(i)
                     dec_usnocand(iofound)=dec_usno(i)
                     epos_usnocand(iofound)=poserr_usno(i)
                     optdist(iofound)=real(dist)
                     optcand_type(iofound)=opt_type(i)
                  endif
               endif
            endif
         ENDDO
         IF (iofound == 0) then
            write(*,'('' NO optical object within '', f5.2,'' arcsec'')') max(epos(1,j)*1.3,2.)
         else
            Do i=1,iofound
               aox=99.99
               aro=99.99
               alpho=99.99
               if ((flux_x .ne. 0.) .and. (flux_usnocand(i,2) .ne. 0.))
     &             aox = 1.-log10(flux_usnocand(i,2)/flux_x)/log10(freq_usnocand(i,2)/2.418e17)
               if ((flux_usnocand(i,2) .ne. 0.) .and. (flux_usnocand(i,4) .ne. 0.))
     &          alpho= 1.-log10(flux_usnocand(i,4)/flux_usnocand(i,2))/log10(freq_usnocand(i,4)/freq_usnocand(i,2))
               if ((flux_r .ne. 0.) .and. (flux_usnocand(i,1) .ne. 0.))
     &         aro = 1.-log10(flux_r/flux_usnocand(i,2))/log10(1.4e9/freq_usnocand(i,2))
               if (optcand_type(i) == 'USNO' ) then
                  write(*,*) optcand_type(i),usnomag_cand(i,2),usnomag_cand(i,4)
               ELSE if (optcand_type(i) == 'GAIA' ) then
                  write(*,*) optcand_type(i),usnomag_cand(i,1)
               ELSE
                  write(*,*) optcand_type(i),usnomag_cand(i,1:5)
               endif
               write(*,'(a,"Optical slope: ",f6.3,",",2x,f7.3," arcsec away")')
     &                 optcand_type(i),alpho,optdist(i)*3600.
            ENDDO
         ENDIF
c         write(*,*) iofound,ii1,ii2,ii3,ii4,ii5
         DO i=1,ioptlc
         CALL DIST_SKY(ra_source(j),dec_source(j),ra_usno(indoptlc(i)),dec_usno(indoptlc(i)),dist)
         !IF (dist < min_dist_pccs100) THEN
            !if (farpart(i) .eq. j) then
             optpart(i)=0
             min_dist=sqrt(poserr_usno(indoptlc(i))**2+epos(1,j)**2)
             if (dist*3600. .lt. max(min_dist,2.)) then
               optpart(i)=j
            ENDIF
          ENDDO


         write(*,*) '.......................UV...........................'
         ii1=0
         ii2=0
         ii3=0
         Do i=1,iuv
            !write(*,*) ra_source(j),dec_source(j),ra_uv(i),dec_uv(i)
            CALL DIST_SKY(ra_source(j),dec_source(j),ra_uv(i),dec_uv(i),dist)
            min_dist=sqrt(poserr_uv(i)**2+epos(1,j)**2)
            !write(*,*) dist*3600,uv_type(i),epos(1,j)*1.3
            IF (( dist*3600. < max(min_dist,2.) ) .and. (uv_type(i) == 'GALEX'))THEN
               iuvfound = iuvfound+1
               ii1=ii1+1
               if (ii1 .eq. 1) then
                  flux_uvcand(iuvfound,1:6)=flux_uv(i,1:6)
                  uflux_uvcand(iuvfound,1:6)=FluxU_uv(i,1:6)
                  lflux_uvcand(iuvfound,1:6)=FluxL_uv(i,1:6)
                  uvmag_cand(iuvfound,1:6)=uvmag(i,1:6)
                  freq_uvcand(iuvfound,1:6)=frequency_uv(i,1:6)
                  ra_uvcand(iuvfound)=ra_uv(i)
                  dec_uvcand(iuvfound)=dec_uv(i)
                  epos_uvcand(iuvfound)=poserr_uv(i)
                  uvdist(iuvfound)=real(dist)
                  uvcand_type(iuvfound)=uv_type(i)
               else
                  iuvfound=iuvfound-1
                  if (uvdist(iuvfound) .gt. dist) then
                     flux_uvcand(iuvfound,1:6)=flux_uv(i,1:6)
                     uflux_uvcand(iuvfound,1:6)=FluxU_uv(i,1:6)
                     lflux_uvcand(iuvfound,1:6)=FluxL_uv(i,1:6)
                     uvmag_cand(iuvfound,1:6)=uvmag(i,1:6)
                     freq_uvcand(iuvfound,1:6)=frequency_uv(i,1:6)
                     ra_uvcand(iuvfound)=ra_uv(i)
                     dec_uvcand(iuvfound)=dec_uv(i)
                     epos_uvcand(iuvfound)=poserr_uv(i)
                     uvdist(iuvfound)=real(dist)
                     uvcand_type(iuvfound)=uv_type(i)
                  endif
               endif
            ELSE IF (( dist*3600. < max(min_dist,2.) ) .and. (uv_type(i) == 'XMMOM'))THEN
               iuvfound = iuvfound+1
               ii2=ii2+1
               if (ii2 .eq. 1) then
                  flux_uvcand(iuvfound,1:6)=flux_uv(i,1:6)
                  uflux_uvcand(iuvfound,1:6)=FluxU_uv(i,1:6)
                  lflux_uvcand(iuvfound,1:6)=FluxL_uv(i,1:6)
                  uvmag_cand(iuvfound,1:6)=uvmag(i,1:6)
                  freq_uvcand(iuvfound,1:6)=frequency_uv(i,1:6)
                  ra_uvcand(iuvfound)=ra_uv(i)
                  dec_uvcand(iuvfound)=dec_uv(i)
                  epos_uvcand(iuvfound)=poserr_uv(i)
                  uvdist(iuvfound)=real(dist)
                  uvcand_type(iuvfound)=uv_type(i)
               else
                  iuvfound=iuvfound-1
                  if (uvdist(iuvfound) .gt. dist) then
                     flux_uvcand(iuvfound,1:6)=flux_uv(i,1:6)
                     uflux_uvcand(iuvfound,1:6)=FluxU_uv(i,1:6)
                     lflux_uvcand(iuvfound,1:6)=FluxL_uv(i,1:6)
                     uvmag_cand(iuvfound,1:6)=uvmag(i,1:6)
                     freq_uvcand(iuvfound,1:6)=frequency_uv(i,1:6)
                     ra_uvcand(iuvfound)=ra_uv(i)
                     dec_uvcand(iuvfound)=dec_uv(i)
                     epos_uvcand(iuvfound)=poserr_uv(i)
                     uvdist(iuvfound)=real(dist)
                     uvcand_type(iuvfound)=uv_type(i)
                  endif
               endif
            ELSE IF (( dist*3600. < max(min_dist,2.) ) .and. (uv_type(i) == 'UVOT'))THEN
               iuvfound = iuvfound+1
               ii3=ii3+1
               flux_uvcand(iuvfound,1:6)=flux_uv(i,1:6)
               uflux_uvcand(iuvfound,1:6)=FluxU_uv(i,1:6)
               lflux_uvcand(iuvfound,1:6)=FluxL_uv(i,1:6)
               uvmag_cand(iuvfound,1:6)=uvmag(i,1:6)
               freq_uvcand(iuvfound,1:6)=frequency_uv(i,1:6)
               ra_uvcand(iuvfound)=ra_uv(i)
               dec_uvcand(iuvfound)=dec_uv(i)
               epos_uvcand(iuvfound)=poserr_uv(i)
               uvdist(iuvfound)=real(dist)
               uvcand_type(iuvfound)=uv_type(i)
               !write(*,*) iuvfound,dec_uvcand(iuvfound),epos_uvcand(iuvfound)
            endif
         enddo
         IF (iuvfound == 0) then
            write(*,'('' NO UV object within '',f5.2,'' arcsec'')') max(epos(1,j)*1.3,2.)
         else
            do i=1,iuvfound
               !write(*,*) i,dec_uvcand(i),epos_uvcand(i)
               auvx=99.99
               aruv=99.99
               alphauv=99.99
               if ((flux_x .ne. 0.) .and. (flux_uvcand(i,4) .ne. 0.))
     &         auvx = 1.-log10(flux_uvcand(i,4)/flux_x)/log10(freq_uvcand(i,4)/2.418e17) !w1 Fuv
               if ((flux_r .ne. 0.) .and. (flux_uvcand(i,4) .ne. 0.))
     &         aruv = 1.-log10(flux_r/flux_uvcand(i,4))/log10(1.4e9/freq_uvcand(i,4))
c               write(*,'(a,"UV-X-ray slope: ",f6.3,",",2x,"radio-UV slope: ",f6.3,",",2x,f7.3," arcsec away")')
c     &              uvcand_type(i),auvx,aruv,uvdist(i)*3600.
c               if (uvcand_type(i) == 'GALEX' ) then
c                  write(*,*) uvcand_type(i),uvmag_cand(i,3),uvmag_cand(i,4)
c               ELSE
c                  write(*,*) uvcand_type(i),uvmag_cand(i,1:6)
c               endif
               if ((uvcand_type(i) == 'UVOT') .or. (uvcand_type(i) == 'XMMOM')) then
                  if ((flux_uvcand(i,1) .ne. 0.) .and. (flux_uvcand(i,6) .ne. 0.))
     &             alphauv = 1.-log10(flux_uvcand(i,1)/flux_uvcand(i,6))/log10(freq_uvcand(i,1)/freq_uvcand(i,6))
c u to w2
c                  write(*,*) 'UV slope',alphauv
               endif
            enddo
         endif
c         write(*,*) iuvfound,ii1,ii2,ii3

         ixxfound=0
         write(*,*) '.......................hard X-ray and XRT spectral data..................'
         do i=1,ixray
            xraypart(i)=0
            call Dist_sky(ra_source(j),dec_source(j),ra_xray(i),dec_xray(i),dist)
            !write(*,*) poserr_xray(i),epos(1,j)
            min_dist=sqrt(poserr_xray(i)**2+epos(1,j)**2)
            !write(*,*) min_dist,dist*3600.
            if (dist*3600. < min_dist ) then !5 arcsec fixed value
               xraypart(i)=j
               ixxfound=ixxfound+1
               if ((xray_type(i) == 'BAT105m')) write(*,'(a,2x,f5.3,2x,"acrmin away")') xray_type(i),dist*60.
               if ((xray_type(i) == 'XRTSPEC') .and. (xrtspind(i) .eq. 1))
     &           write(*,'(a,2x,f7.3,2x,"acrsec away")') xray_type(i),dist*3600.
            endif
         enddo
c         if (ixray == 0 ) write(*,*) NO BAT detection within 8 arcmin
         IF (ixxfound == 0) write(*,'('' NO BAT or XRT detection within '',f7.3,'' arcmin'')') min_dist

         igamfound=0
         write(*,*) '.......................Gamma-ray..................'
c         debl=.false.
         do i=1,igam
            call Dist_sky(ra_source(j),dec_source(j),ra_gam(i),dec_gam(i),dist)
            !if (dist < min_dist_gam) then
!            if (gampart(i) .eq. j) then
            gampart(i)=0
            min_dist=sqrt(poserr_gam(i)**2+epos(1,j)**2)
            !write(*,*) dist*3600.,min_dist*1.1,gam_type(i)
            if (dist*3600. .lt. min_dist*1.1) then
               gampart(i)=j
               igamfound=igamfound+1
               if (igamfound > 20) stop 'Too many Gamma-ray candidate'
               if (gam_type(i) == '2BIGB') then
                  if (bigbind(i) .eq. 1) then
                     write(*,'(a,"photon index: ",f5.3,",",2x,f7.3," arcmin away")') gam_type(i),slope_gam(i,1),dist*60
                  endif
               else if (gam_type(i) == '2FHL') then
c                  debl=.true.
                  write(*,'(a,"photon index: ",f5.3,",",2x,f7.3," arcmin away")') gam_type(i),slope_gam(i,2),dist*60
               else
c                  if ((gam_type(i) == '3FHL') .or. (gam_type(i) == '3FGL') .or. (gam_type(i) == '4FGL')) debl=.true.
                  write(*,'(a,"photon index: ",f5.3,",",2x,f7.3," arcmin away")') gam_type(i),slope_gam(i,1),dist*60
               endif
            else
              eblnn(i)=0
            endif
         enddo
         !write(*,*) igamfound
         IF (igamfound == 0) write(*,'('' NO Gamma-ray detection within '',f5.0,'' arcmin'')') min_dist_gam*60.

         iflcuvfound=0
         do i=1,iflcuv
            flcuvpart(i)=0
            call Dist_sky(ra_source(j),dec_source(j),ra_flcuv(i),dec_flcuv(i),dist)
c            write(*,*) i,dist*3600.,epos(1,j)
ccc            min_dist=sqrt(poserr_flcuv(i)**2+epos(1,j)**2)
                !write(*,*) dist*60.
                if (dist*60. < 30. ) then !5 arcsec fixed value
                flcuvpart(i)=j
                iflcuvfound=iflcuvfound+1
            endif
         enddo

         write(*,*) '.......................TeV..................'
         do i=1,ivhe
            if (filen_v(i) .eq. j) then ! no coordinate, so use the file number
               if (magicind(i) .eq. 1) write(*,*) 'MAGIC source'
               if (veritind(i) .eq. 1) write(*,*) 'VERITAS source'
            else
               eblnn(igam+i)=0
            endif
         enddo
         if (ivhe .eq. 0) then
            write(*,*) 'NO TeV detection within 10 arcmin'
         endif

c plot the sed, and output the list, decide the source type

cDO i =0,5
c  types(i) = 0
cENDDO
cno_found = 0
cDO i = 0,5
c   IF (types(i) > no_found) THEN
c      no_found=types(i)
c      type_average = i
c   ENDIF
cENDDO
         allocate(rrxx_flag(arrsize(2),arrsize(1)),rrxx_ref(arrsize(2),arrsize(1)))

         write(*,*) 'TEST number',eblnn(1:igam)
         !CALL graphic_code (flux_x,flux_radio(k)/const(k),type_average,code)
         !write(lu_output,*) ra_source(k),dec_source(k),code
         write(14,'(i4,2x,a,2(2x,f9.5),2x,i2,2x,f5.3)') j,"matched source",ra_source(j),dec_source(j),typer(j),redshift
         write(14,'(2(a,2x),2(a,1x),5(a,2x))') " Frequency","   nufnu  "," nufnu unc."," nufnu unc.","start time"," end time ","Flag","Catalog        ","Reference"
         write(14,'(6(a,2x))') "    Hz    "," erg/cm2/s","   upper  ","   lower  ","    MJD   ","    MJD   "
         write(14,'(a)') "---------------------------------------------------------------------------------------------------------------------------"
         do i=1,npt(j)
c the refs file
            rrxx_ref(i,j)=iref+1
            refs(iref+1)='TBD'
            do r=1,iref
               if (rrxx_type(i,j) == name_cat(r)) THEN
                  rrxx_ref(i,j)=r
               endif
            enddo
            if ((flux(i,j) .eq. lflux(i,j)) .and. (flux(i,j) .eq. uflux(i,j))) then
               rrxx_flag(i,j)=' UL '
            else if ((lflux(i,j) .eq. 0.) .and. (uflux(i,j) .ne. 0.) ) then
               rrxx_flag(i,j)=' UL '
            else
               rrxx_flag(i,j)=' Det'
            endif
            write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency(i,j),flux(i,j),uflux(i,j),lflux(i,j),
     &       mjdst_rrxx(i,j),mjded_rrxx(i,j),rrxx_flag(i,j),rrxx_type(i,j),refs(rrxx_ref(i,j))
            if (frequency(i,j) .lt. 1.E10) then
               call graphic_code(flux(i,j),11,code)
               write(lu_output,'(f9.5,2x,f9.5,2x,i6,2x,f8.3)') ra_rrxx(i,j),dec_rrxx(i,j),int(code),epos(i,j)
            else if ((frequency(i,j) .eq. 2.418E17) .or. ((rrxx_type(i,j) == 'XRTDEEP') .and.
     &           (ra_rrxx(i,j) .ne. 0.)) .or. ((rrxx_type(i,j) == 'MAXI') .and. (ra_rrxx(i,j) .ne. 0.)))then
               call graphic_code(flux(i,j),81,code)
               write(lu_output,'(f9.5,2x,f9.5,2x,i6,2x,f8.3)') ra_rrxx(i,j),dec_rrxx(i,j),int(code),epos(i,j)
            endif
         enddo

         deallocate(rrxx_flag,rrxx_ref)
c         write(17,'(i4,2x,a,2(2x,f9.5),2x,i2)') j,"matched source",ra_source(j),dec_source(j),typer(j)
c         do i=1,nptlc(j)
c            write(17,'(4(es10.3,2x),2(f10.4,2x),i2)') frequency_lc(i,j),flux_lc(i,j),
c     *        uflux_lc(i,j),lflux_lc(i,j),mjdst_lc(i,j),mjden_lc(i,j),int(lcurve_type(i,j))
c         enddo
         do i=1,ilowrfound
            lowr_ref(i)=iref+1
            refs(iref+1)='TBD'
            do r=1,iref
               if (lowrcand_type(i) == name_cat(r)) THEN
                  lowr_ref(i)=r
               endif
            enddo
            if ((flux_lowrcand(i) .eq. lflux_lowrcand(i)) .and. (flux_lowrcand(i) .eq. uflux_lowrcand(i))) then
               lowr_flag(i)=' UL '
            else if ((lflux_lowrcand(i) .eq. 0.) .and. (uflux_lowrcand(i) .ne. 0.) ) then
               lowr_flag(i)=' UL '
            else
               lowr_flag(i)=' Det'
            endif
            write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') freq_lowrcand(i),flux_lowrcand(i),uflux_lowrcand(i),
     *       lflux_lowrcand(i),mjdavg,mjdavg,lowr_flag(i),lowrcand_type(i),refs(lowr_ref(i))
         enddo
         do i=1,i4p8
            f4p8_ref(i)=iref+1
            refs(iref+1)='TBD'
            do r=1,iref
               if (f4p8_type(i) == name_cat(r)) THEN
                  f4p8_ref(i)=r
               endif
            enddo
            if (f4p8part(i) .eq. j) then
            if ((f4p8_type(i) == 'PMN') .or. (f4p8_type(i) == 'GB87') .or.
     &           (f4p8_type(i) == 'KUEHR') .or. (f4p8_type(i) == 'GB6')) THEN
               if ((flux_4p8(i,1) .eq. FluxL_4p8(i,1)) .and. (flux_4p8(i,1) .eq. FluxU_4p8(i,1))) then
                  f4p8_flag(i,1)=' UL '
               else if ((FluxL_4p8(i,1) .eq. 0.) .and. (FluxU_4p8(i,1) .ne. 0.) ) then
                  f4p8_flag(i,1)=' UL '
               else
                  f4p8_flag(i,1)=' Det'
               endif
               write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_4p8(i,1),flux_4p8(i,1),FluxU_4p8(i,1),
     &             FluxL_4p8(i,1),mjdavg,mjdavg,f4p8_flag(i,1),f4p8_type(i),refs(f4p8_ref(i))
            else if (f4p8_type(i) == 'CRATES') then
               do s=2,2
                  if ((flux_4p8(i,s) .eq. FluxL_4p8(i,s)) .and. (flux_4p8(i,s) .eq. FluxU_4p8(i,s))) then
                     f4p8_flag(i,s)=' UL '
                  else if ((FluxL_4p8(i,s) .eq. 0.) .and. (FluxU_4p8(i,s) .ne. 0.) ) then
                     f4p8_flag(i,s)=' UL '
                  else
                     f4p8_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_4p8(i,s),flux_4p8(i,s),FluxU_4p8(i,s),
     &            FluxL_4p8(i,s),mjdavg,mjdavg,f4p8_flag(i,s),f4p8_type(i),refs(f4p8_ref(i))
               enddo
            else if (f4p8_type(i) == 'ATPMN') then
               do s=1,2
                  if ((flux_4p8(i,s) .eq. FluxL_4p8(i,s)) .and. (flux_4p8(i,s) .eq. FluxU_4p8(i,s))) then
                     f4p8_flag(i,s)=' UL '
                  else if ((FluxL_4p8(i,s) .eq. 0.) .and. (FluxU_4p8(i,s) .ne. 0.) ) then
                     f4p8_flag(i,s)=' UL '
                  else
                     f4p8_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_4p8(i,s),flux_4p8(i,s),FluxU_4p8(i,s),
     &              FluxL_4p8(i,s),mjdavg,mjdavg,f4p8_flag(i,s),f4p8_type(i),refs(f4p8_ref(i))
               enddo
            else if ((f4p8_type(i) == 'AT20G') .or. (f4p8_type(i) == 'NORTH20')
     &               .or. (f4p8_type(i) == 'F357det')) then
               do s=1,3
                  if ((flux_4p8(i,s) .eq. FluxL_4p8(i,s)) .and. (flux_4p8(i,s) .eq. FluxU_4p8(i,s))) then
                     f4p8_flag(i,s)=' UL '
                  else if ((FluxL_4p8(i,s) .eq. 0.) .and. (FluxU_4p8(i,s) .ne. 0.) ) then
                     f4p8_flag(i,s)=' UL '
                  else
                     f4p8_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_4p8(i,s),flux_4p8(i,s),FluxU_4p8(i,s),
     &              FluxL_4p8(i,s),mjdavg,mjdavg,f4p8_flag(i,s),f4p8_type(i),refs(f4p8_ref(i))
               enddo
            endif
            endif
         enddo
         do i=1,ipccs100
            pccs100_ref(i)=iref+1
            refs(iref+1)='TBD'
            do r=1,iref
               if (pccs100_type(i) == name_cat(r)) THEN
                  pccs100_ref(i)=r
               endif
            enddo
            if (pccspart(i) .eq. j) then
               if (pccs100_type(i) == 'PCCS2') then
                  if ((flux_pccs100(i,1) .eq. FluxL_pccs100(i,1)) .and. (flux_pccs100(i,1) .eq. FluxU_pccs100(i,1))) then
                     pccs100_flag(i,1)=' UL '
                  else if ((FluxL_pccs100(i,1) .eq. 0.) .and. (FluxU_pccs100(i,1) .ne. 0.) ) then
                     pccs100_flag(i,1)=' UL '
                  else
                     pccs100_flag(i,1)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,i3,a,2x,a)') frequency_pccs100(i,1),flux_pccs100(i,1),
     &              FluxU_pccs100(i,1),FluxL_pccs100(i,1),mjdavg,mjdavg,pccs100_flag(i,1),
     &              'PCCS2',int(frequency_pccs100(i,1)/1.e9),' GHz',refs(pccs100_ref(i))
               else if (pccs100_type(i) == 'ALMA') then
                  if ((flux_pccs100(i,1) .eq. FluxL_pccs100(i,1)) .and. (flux_pccs100(i,1) .eq. FluxU_pccs100(i,1))) then
                     pccs100_flag(i,1)=' UL '
                  else if ((FluxL_pccs100(i,1) .eq. 0.) .and. (FluxU_pccs100(i,1) .ne. 0.) ) then
                     pccs100_flag(i,1)=' UL '
                  else
                     pccs100_flag(i,1)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_pccs100(i,1),flux_pccs100(i,1),
     &             FluxU_pccs100(i,1),FluxL_pccs100(i,1),mjdst_alma(i),
     &             mjded_alma(i),pccs100_flag(i,1),pccs100_type(i),refs(pccs100_ref(i))
               else
                  do s=1,9
                     if ((flux_pccs100(i,s) .eq. FluxL_pccs100(i,s)) .and. (flux_pccs100(i,s) .eq. FluxU_pccs100(i,s))) then
                        pccs100_flag(i,s)=' UL '
                     else if ((FluxL_pccs100(i,s) .eq. 0.) .and. (FluxU_pccs100(i,s) .ne. 0.) ) then
                        pccs100_flag(i,s)=' UL '
                     else
                        pccs100_flag(i,s)=' Det'
                     endif
                     write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_pccs100(i,s),flux_pccs100(i,s),FluxU_pccs100(i,s),
     &               FluxL_pccs100(i,s),mjdavg,mjdavg,pccs100_flag(i,s),pccs100_type(i),refs(pccs100_ref(i))
                  enddo
               endif
            endif
         enddo
         do i=1,ifar
            far_ref(i)=iref+1
            refs(iref+1)='TBD'
            do r=1,iref
               if (far_type(i) == name_cat(r)) THEN
                  far_ref(i)=r
               endif
            enddo
            if (farpart(i) .eq. j) then
               if (far_type(i) == 'SPIRE') THEN
                  if ((flux_far(i,1) .eq. FluxL_far(i,1)) .and. (flux_far(i,1) .eq. FluxU_far(i,1))) then
                     far_flag(i,1)=' UL '
                  else if ((FluxL_far(i,1) .eq. 0.) .and. (FluxU_far(i,1) .ne. 0.) ) then
                     far_flag(i,1)=' UL '
                  else
                     far_flag(i,1)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,i3,4x,a)') frequency_far(i,1),flux_far(i,1),FluxU_far(i,1),
     &            FluxL_far(i,1),mjdavg,mjdavg,far_flag(i,1),'SPIRE',int(3.E14/frequency_far(i,1)),refs(far_ref(i))
               else if ((far_type(i) == 'IRAS-PSC') .or. (far_type(i) == 'AKARIBSC')) THEN
                  do s=1,4
                     if ((flux_far(i,s) .eq. FluxL_far(i,s)) .and. (flux_far(i,s) .eq. FluxU_far(i,s))) then
                        far_flag(i,s)=' UL '
                     else if ((FluxL_far(i,s) .eq. 0.) .and. (FluxU_far(i,s) .ne. 0.) ) then
                        far_flag(i,s)=' UL '
                     else
                        far_flag(i,s)=' Det'
                     endif
                     write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_far(i,s),flux_far(i,s),FluxU_far(i,s),
     &               FluxL_far(i,s),mjdavg,mjdavg,far_flag(i,s),far_type(i),refs(far_ref(i))
                  enddo
               else
                  do s=1,5
                     if ((flux_far(i,s) .eq. FluxL_far(i,s)) .and. (flux_far(i,s) .eq. FluxU_far(i,s))) then
                        far_flag(i,s)=' UL '
                     else if ((FluxL_far(i,s) .eq. 0.) .and. (FluxU_far(i,s) .ne. 0.) ) then
                        far_flag(i,s)=' UL '
                     else
                        far_flag(i,s)=' Det'
                     endif
                     write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_far(i,s),flux_far(i,s),FluxU_far(i,s),
     &               FluxL_far(i,s),mjdavg,mjdavg,far_flag(i,s),far_type(i),refs(far_ref(i))
                  enddo
               endif
            endif
         enddo
         do i=1,iirfound
            ir_ref(i)=iref+1
            refs(iref+1)='TBD'
            do r=1,iref
               if (ircand_type(i) == name_cat(r)) THEN
                  ir_ref(i)=r
               endif
            enddo
            if (ircand_type(i) == 'WISE') then
               call graphic_code(irmag_cand(i,2),51,code)
               do s=1,4
                  if ((flux_ircand(i,s) .eq. lflux_ircand(i,s)) .and. (flux_ircand(i,s) .eq. uflux_ircand(i,s))) then
                     ir_flag(i,s)=' UL '
                  else if ((lflux_ircand(i,s) .eq. 0.) .and. (uflux_ircand(i,s) .ne. 0.) ) then
                     ir_flag(i,s)=' UL '
                  else
                     ir_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') freq_ircand(i,s),flux_ircand(i,s),uflux_ircand(i,s),
     &             lflux_ircand(i,s),mjdavg,mjdavg,ir_flag(i,s),ircand_type(i),refs(ir_ref(i))
               enddo
            else if (ircand_type(i) == 'WISEME') then
               do s=1,4
                  if ((flux_ircand(i,s) .eq. lflux_ircand(i,s)) .and. (flux_ircand(i,s) .eq. uflux_ircand(i,s))) then
                     ir_flag(i,s)=' UL '
                  else if ((lflux_ircand(i,s) .eq. 0.) .and. (uflux_ircand(i,s) .ne. 0.) ) then
                     ir_flag(i,s)=' UL '
                  else
                     ir_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') freq_ircand(i,s),flux_ircand(i,s),uflux_ircand(i,s),
     &          lflux_ircand(i,s),mjdst_irlccand(i),mjded_irlccand(i),ir_flag(i,s),ircand_type(i),refs(ir_ref(i))
               enddo
            else if (ircand_type(i) == 'NEOWISE') then
               do s=1,2
                  if ((flux_ircand(i,s) .eq. lflux_ircand(i,s)) .and. (flux_ircand(i,s) .eq. uflux_ircand(i,s))) then
                     ir_flag(i,s)=' UL '
                  else if ((lflux_ircand(i,s) .eq. 0.) .and. (uflux_ircand(i,s) .ne. 0.) ) then
                     ir_flag(i,s)=' UL '
                  else
                     ir_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') freq_ircand(i,s),flux_ircand(i,s),uflux_ircand(i,s),
     &   lflux_ircand(i,s),mjdst_irlccand(i),mjded_irlccand(i),ir_flag(i,s),ircand_type(i),refs(ir_ref(i))
                      enddo
            else
               call graphic_code(irmag_cand(i,2),52,code)
               do s=2,4
                  if ((flux_ircand(i,s) .eq. lflux_ircand(i,s)) .and. (flux_ircand(i,s) .eq. uflux_ircand(i,s))) then
                     ir_flag(i,s)=' UL '
                  else if ((lflux_ircand(i,s) .eq. 0.) .and. (uflux_ircand(i,s) .ne. 0.) ) then
                     ir_flag(i,s)=' UL '
                  else
                     ir_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') freq_ircand(1,s),flux_ircand(1,s),uflux_ircand(i,s),
     &             lflux_ircand(i,s),mjdavg,mjdavg,ir_flag(i,s),ircand_type(i),refs(ir_ref(i))
               enddo
            endif
            if ((ircand_type(i) /= 'WISEME') .and. (ircand_type(i) /= 'NEOWISE'))
     &        write(lu_output,'(f9.5,2x,f9.5,2x,i6,2x,f8.3)') ra_ircand(i),dec_ircand(i),int(code),epos_ircand(i)
         enddo
         do i=1,iofound
            opt_ref(i)=iref+1
            refs(iref+1)='TBD'
            do r=1,iref
               if (optcand_type(i) == name_cat(r)) THEN
                  opt_ref(i)=r
               endif
            enddo
            if (optcand_type(i) == 'USNO' ) then
               if ((flux_usnocand(i,2) .eq. lflux_usnocand(i,2)) .and. (flux_usnocand(i,2) .eq. uflux_usnocand(i,2))) then
                  opt_flag(i,2)=' UL '
               else if ((lflux_usnocand(i,2) .eq. 0.) .and. (uflux_usnocand(i,2) .ne. 0.) ) then
                  opt_flag(i,2)=' UL '
               else
                  opt_flag(i,2)=' Det'
               endif
               if ((flux_usnocand(i,4) .eq. lflux_usnocand(i,4)) .and. (flux_usnocand(i,4) .eq. uflux_usnocand(i,4))) then
                  opt_flag(i,4)=' UL '
               else if ((lflux_usnocand(i,4) .eq. 0.) .and. (uflux_usnocand(i,4) .ne. 0.) ) then
                  opt_flag(i,4)=' UL '
               else
                  opt_flag(i,4)=' Det'
               endif
               write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a)') freq_usnocand(i,2),flux_usnocand(i,2),uflux_usnocand(i,2),
     &         lflux_usnocand(i,2),mjdavg,mjdavg,opt_flag(i,2),optcand_type(i),refs(opt_ref(i))
               write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a)') freq_usnocand(i,4),flux_usnocand(i,4),uflux_usnocand(i,4),
     &          lflux_usnocand(i,4),mjdavg,mjdavg,opt_flag(i,4),optcand_type(i),refs(opt_ref(i))
            ELSE if (optcand_type(i) == 'GAIA' ) then
                if ((flux_usnocand(i,1) .eq. lflux_usnocand(i,1)) .and. (flux_usnocand(i,1) .eq. uflux_usnocand(i,1))) then
                   opt_flag(i,1)=' UL '
                else if ((lflux_usnocand(i,1) .eq. 0.) .and. (uflux_usnocand(i,1) .ne. 0.) ) then
                   opt_flag(i,1)=' UL '
                else
                   opt_flag(i,1)=' Det'
                endif
                write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') freq_usnocand(i,1),flux_usnocand(i,1),uflux_usnocand(i,1),
     &           lflux_usnocand(i,1),mjdavg,mjdavg,opt_flag(i,1),optcand_type(i),refs(opt_ref(i))
            ELSE
               do s=1,5
                  if ((flux_usnocand(i,s) .eq. lflux_usnocand(i,s)) .and. (flux_usnocand(i,s) .eq. uflux_usnocand(i,s))) then
                     opt_flag(i,s)=' UL '
                  else if ((lflux_usnocand(i,s) .eq. 0.) .and. (uflux_usnocand(i,s) .ne. 0.) ) then
                     opt_flag(i,s)=' UL '
                  else
                     opt_flag(i,s)=' Det'
                  endif
                write(14,'(4(es10.3,2x),2(f10.4,2x)a,2x,a,2x,a)') freq_usnocand(i,s),flux_usnocand(i,s),uflux_usnocand(i,s),
     &            lflux_usnocand(i,s),mjdavg,mjdavg,opt_flag(i,s),optcand_type(i),refs(opt_ref(i))
               enddo
            endif
            intensity=max(usnomag_cand(i,1),usnomag_cand(i,2),usnomag_cand(i,3),usnomag_cand(i,4),usnomag_cand(i,5))
            if (optcand_type(i) == 'SDSS') call graphic_code(usnomag_cand(i,3),61,code)
            if (optcand_type(i) == 'HSTGSC') call graphic_code(intensity,62,code)
            if (optcand_type(i) == 'USNO') call graphic_code(usnomag_cand(i,2),63,code)
            if (optcand_type(i) == 'PanSTARRS') call graphic_code(usnomag_cand(i,3),61,code)
            if (optcand_type(i) == 'GAIA') call graphic_code(usnomag_cand(i,2),63,code)
        write(lu_output,'(f9.5,2x,f9.5,2x,i6,2x,f8.3)') ra_usnocand(i),dec_usnocand(i),int(code),epos_usnocand(i)
         enddo
         do i =1,ioptlc
           optlc_ref(i)=iref+1
           refs(iref+1)='TBD'
           do r=1,iref
              if (opt_type(indoptlc(i)) == name_cat(r)) THEN
                optlc_ref(i)=r
              endif
           enddo
           if (optpart(i) .eq. j) then
              if ((flux_usno(indoptlc(i),1) .eq. FluxL_usno(indoptlc(i),1)) .and. (flux_usno(indoptlc(i),1) .eq. FluxU_usno(indoptlc(i),1))) then
                 optlc_flag(i)=' UL '
              else if ((FluxL_usno(indoptlc(i),1) .eq. 0.) .and. (FluxU_usno(indoptlc(i),1) .ne. 0.) ) then
                 optlc_flag(i)=' UL '
              else
                 optlc_flag(i)=' Det'
              endif
              write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_usno(indoptlc(i),1),flux_usno(indoptlc(i),1),FluxU_usno(indoptlc(i),1),
     &            FluxL_usno(indoptlc(i),1),mjdst_optlc(indoptlc(i)),mjdst_optlc(indoptlc(i)),optlc_flag(i),opt_type(indoptlc(i)),refs(optlc_ref(i))
           endif
         enddo
         do i=1,iuvfound
            uv_ref(i)=iref+1
            refs(iref+1)='TBD'
            do r=1,iref
               if (uvcand_type(i) == name_cat(r)) THEN
                  uv_ref(i)=r
               endif
            enddo
            if (uvcand_type(i) == 'GALEX') then
               do s=3,4
                  if ((flux_uvcand(i,s) .eq. lflux_uvcand(i,s)) .and. (flux_uvcand(i,s) .eq. uflux_uvcand(i,s))) then
                     uv_flag(i,s)=' UL '
                  else if ((lflux_uvcand(i,s) .eq. 0.) .and. (uflux_uvcand(i,s) .ne. 0.) ) then
                     uv_flag(i,s)=' UL '
                  else
                     uv_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x)a,2x,a,2x,a)') freq_uvcand(i,s),flux_uvcand(i,s),uflux_uvcand(i,s),
     &            lflux_uvcand(i,s),mjdavg,mjdavg,uv_flag(i,s),uvcand_type(i),refs(uv_ref(i))
               enddo
            else
               do s=1,6
                  if ((flux_uvcand(i,s) .eq. lflux_uvcand(i,s)) .and. (flux_uvcand(i,s) .eq. uflux_uvcand(i,s))) then
                     uv_flag(i,s)=' UL '
                  else if ((lflux_uvcand(i,s) .eq. 0.) .and. (uflux_uvcand(i,s) .ne. 0.) ) then
                     uv_flag(i,s)=' UL '
                  else
                     uv_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x)a,2x,a,2x,a)') freq_uvcand(i,s),flux_uvcand(i,s),uflux_uvcand(i,s),
     &            lflux_uvcand(i,s),mjdavg,mjdavg,uv_flag(i,s),uvcand_type(i),refs(uv_ref(i))
               enddo
            endif
            intensity=max(uvmag_cand(i,1),uvmag_cand(i,2),uvmag_cand(i,3),uvmag_cand(i,4),uvmag_cand(i,5))
            call graphic_code(intensity,71,code)
            write(lu_output,'(f9.5,2x,f9.5,2x,i6,2x,f8.3)') ra_uvcand(i),dec_uvcand(i),int(code),epos_uvcand(i)
         enddo
         do i=1,ixray
            xray_ref(i)=iref+1
            refs(iref+1)='TBD'
            do r=1,iref
               if (xray_type(i) == name_cat(r)) THEN
                  xray_ref(i)=r
               endif
            enddo
            if (xraypart(i) .eq. j) then
               if ((xray_type(i) == 'XRTSPEC') .or. (xray_type(i) == 'BEPPOSAX')
     &           .or. (xray_type(i) == 'NuBlazar')) then
                  if ((flux_xray(i,1) .eq. FluxL_xray(i,1)) .and. (flux_xray(i,1) .eq. FluxU_xray(i,1))) then
                     xray_flag(i,1)=' UL '
                  else if ((FluxL_xray(i,1) .eq. 0.) .and. (FluxU_xray(i,1) .ne. 0.) ) then
                     xray_flag(i,1)=' UL '
                  else
                     xray_flag(i,1)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_xray(i,1),flux_xray(i,1),FluxU_xray(i,1),
     &            FluxL_xray(i,1),mjdst_xrt(i),mjded_xrt(i),xray_flag(i,1),xray_type(i),refs(xray_ref(i))
               else if (xray_type(i) == 'BAT105m') then
                  do s=1,2
                     if ((flux_xray(i,s) .eq. FluxL_xray(i,s)) .and. (flux_xray(i,s) .eq. FluxU_xray(i,s))) then
                        xray_flag(i,s)=' UL '
                     else if ((FluxL_xray(i,s) .eq. 0.) .and. (FluxU_xray(i,s) .ne. 0.) ) then
                        xray_flag(i,s)=' UL '
                     else
                        xray_flag(i,s)=' Det'
                     endif
                    write(14,'(4(es10.3,2x),2(f10.4,2x)a,2x,a,2x,a)') frequency_xray(i,s),flux_xray(i,s),FluxU_xray(i,s),
     &               FluxL_xray(i,s),mjdavg,mjdavg,xray_flag(i,s),xray_type(i),refs(xray_ref(i))
                  enddo
               else if (xray_type(i) == 'OULC') then
                  do s=1,2
                     if ((flux_xray(i,s) .eq. FluxL_xray(i,s)) .and. (flux_xray(i,s) .eq. FluxU_xray(i,s))) then
                        xray_flag(i,s)=' UL '
                     else if ((FluxL_xray(i,s) .eq. 0.) .and. (FluxU_xray(i,s) .ne. 0.) ) then
                        xray_flag(i,s)=' UL '
                     else
                        xray_flag(i,s)=' Det'
                     endif
                    write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_xray(i,s),flux_xray(i,s),FluxU_xray(i,s),
     &             FluxL_xray(i,s),mjdst_xrt(i),mjded_xrt(i),xray_flag(i,s),xray_type(i),refs(xray_ref(i))
                  enddo
               else if ((xray_type(i) == 'OUSXG') .or. (xray_type(i) == 'OUSXB')) then
                  do s=1,5
                     if ((flux_xray(i,s) .eq. FluxL_xray(i,s)) .and. (flux_xray(i,s) .eq. FluxU_xray(i,s))) then
                        xray_flag(i,s)=' UL '
                     else if ((FluxL_xray(i,s) .eq. 0.) .and. (FluxU_xray(i,s) .ne. 0.) ) then
                        xray_flag(i,s)=' UL '
                     else
                        xray_flag(i,s)=' Det'
                     endif
                     write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_xray(i,s),flux_xray(i,s),FluxU_xray(i,s),
     &               FluxL_xray(i,s),mjdst_xrt(i),mjded_xrt(i),xray_flag(i,s),xray_type(i),refs(xray_ref(i))
                  enddo
               endif
               if (((xrtspind(i) .lt. 2).and. (xray_type(i) == 'XRTSPEC') ).or. (xray_type(i) == 'BAT105m')) then
                  call graphic_code(flux_xray(i,1),82,code)
                  write(lu_output,'(f9.5,2x,f9.5,2x,i6,2x,f8.3)') ra_xray(i),dec_xray(i),int(code),poserr_xray(i)
               endif
            endif
         enddo
         do i=1,igam
            gam_ref(i)=iref+1
            refs(iref+1)='TBD'
            do r=1,iref
               if (gam_type(i) == name_cat(r)) THEN
                  gam_ref(i)=r
               endif
            enddo
            if (gampart(i) .eq. j) then
            if (gam_type(i) == '2FHL' ) then
               call graphic_code(flux_gam(i,1),92,code)
               do s=1,4
                  if ((flux_gam(i,s) .eq. FluxL_gam(i,s)) .and. (flux_gam(i,s) .eq. FluxU_gam(i,s))) then
                     gam_flag(i,s)=' UL '
                  else if ((FluxL_gam(i,s) .eq. 0.) .and. (FluxU_gam(i,s) .ne. 0.) ) then
                     gam_flag(i,s)=' UL '
                  else
                     gam_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_gam(i,s),flux_gam(i,s),FluxU_gam(i,s),
     &             FluxL_gam(i,s),mjdavg,mjdavg,gam_flag(i,s),gam_type(i),refs(gam_ref(i))
               enddo
            else if (gam_type(i) == '3FGL') then
               call graphic_code(flux_gam(i,1),91,code)
               do s=1,7
                  if ((flux_gam(i,s) .eq. FluxL_gam(i,s)) .and. (flux_gam(i,s) .eq. FluxU_gam(i,s))) then
                     gam_flag(i,s)=' UL '
                  else if ((FluxL_gam(i,s) .eq. 0.) .and. (FluxU_gam(i,s) .ne. 0.) ) then
                     gam_flag(i,s)=' UL '
                  else
                     gam_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_gam(i,s),flux_gam(i,s),
     &             FluxU_gam(i,s),FluxL_gam(i,s),mjdavg,mjdavg,gam_flag(i,s),gam_type(i),refs(gam_ref(i))
               enddo
            else if (gam_type(i) == '3FHL') then
               call graphic_code(flux_gam(i,1),92,code)
               do s=1,6
                  if ((flux_gam(i,s) .eq. FluxL_gam(i,s)) .and. (flux_gam(i,s) .eq. FluxU_gam(i,s))) then
                     gam_flag(i,s)=' UL '
                  else if ((FluxL_gam(i,s) .eq. 0.) .and. (FluxU_gam(i,s) .ne. 0.) ) then
                     gam_flag(i,s)=' UL '
                  else
                     gam_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_gam(i,s),flux_gam(i,s),
     &            FluxU_gam(i,s),FluxL_gam(i,s),mjdavg,mjdavg,gam_flag(i,s),gam_type(i),refs(gam_ref(i))
               enddo
            else if (gam_type(i) == '2BIGB') then
               if (bigbind(i) .eq. 1) call graphic_code(flux_gam(i,1),91,code)
               if ((flux_gam(i,1) .eq. FluxL_gam(i,1)) .and. (flux_gam(i,1) .eq. FluxU_gam(i,1))) then
                  gam_flag(i,1)=' UL '
               else if ((FluxL_gam(i,1) .eq. 0.) .and. (FluxU_gam(i,1) .ne. 0.) ) then
                  gam_flag(i,1)=' UL '
               else
                  gam_flag(i,1)=' Det'
               endif
               write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_gam(i,1),flux_gam(i,1),
     &          FluxU_gam(i,1),FluxL_gam(i,1),mjdavg,mjdavg,gam_flag(i,1),gam_type(i),refs(gam_ref(i))
            else if (gam_type(i) == '4FGL-DR3') then
               call graphic_code(flux_gam(i,1),93,code)
               do s=1,9
                  if ((flux_gam(i,s) .eq. FluxL_gam(i,s)) .and. (flux_gam(i,s) .eq. FluxU_gam(i,s))) then
                     gam_flag(i,s)=' UL '
                  else if ((FluxL_gam(i,s) .eq. 0.) .and. (FluxU_gam(i,s) .ne. 0.) ) then
                     gam_flag(i,s)=' UL '
                  else
                     gam_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_gam(i,s),flux_gam(i,s),
     &            FluxU_gam(i,s),FluxL_gam(i,s),mjdavg,mjdavg,gam_flag(i,s),gam_type(i),refs(gam_ref(i))
               enddo
            else if (gam_type(i) == '2AGILE') then
               call graphic_code(flux_gam(i,1),94,code)
               do s=1,5
                  if ((flux_gam(i,s) .eq. FluxL_gam(i,s)) .and. (flux_gam(i,s) .eq. FluxU_gam(i,s))) then
                     gam_flag(i,s)=' UL '
                  else if ((FluxL_gam(i,s) .eq. 0.) .and. (FluxU_gam(i,s) .ne. 0.) ) then
                     gam_flag(i,s)=' UL '
                  else
                     gam_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_gam(i,s),flux_gam(i,s),
     &             FluxU_gam(i,s),FluxL_gam(i,s),mjdavg,mjdavg,gam_flag(i,s),gam_type(i),refs(gam_ref(i))
               enddo
            else if (gam_type(i) == 'FermiMeV') then
               call graphic_code(flux_gam(i,1),94,code)
               do s=1,2
                  if ((flux_gam(i,s) .eq. FluxL_gam(i,s)) .and. (flux_gam(i,s) .eq. FluxU_gam(i,s))) then
                     gam_flag(i,s)=' UL '
                  else if ((FluxL_gam(i,s) .eq. 0.) .and. (FluxU_gam(i,s) .ne. 0.) ) then
                     gam_flag(i,s)=' UL '
                  else
                     gam_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_gam(i,s),flux_gam(i,s),
     &              FluxU_gam(i,s),FluxL_gam(i,s),mjdavg,mjdavg,gam_flag(i,s),gam_type(i),refs(gam_ref(i))
               enddo
            endif
            if (gam_type(i) == '2BIGB') then
               if (bigbind(i) .eq. 1) write(lu_output,'(f9.5,2x,f9.5,2x,i6,2x,f8.3)') ra_gam(i),dec_gam(i),int(code),poserr_gam(i)
            else
               write(lu_output,'(f9.5,2x,f9.5,2x,i6,2x,f8.3)') ra_gam(i),dec_gam(i),int(code),poserr_gam(i)
            endif
            endif
         enddo
         do i=1,iflcuv
            if (flcuvpart(i) .eq. j) then
               do s=1,4
                  if ((flux_flcuv(i,s) .eq. FluxL_flcuv(i,s)) .and. (flux_flcuv(i,s) .eq. FluxU_flcuv(i,s))) then
                     flcuv_flag(i,s)=' UL '
                  else if ((FluxL_flcuv(i,s) .eq. 0.) .and. (FluxU_flcuv(i,s) .ne. 0.) ) then
                     flcuv_flag(i,s)=' UL '
                  else
                     flcuv_flag(i,s)=' Det'
                  endif
                  write(14,'(4(es10.3,2x),2(f10.4,2x)a,2x,a,"  Fermi LAT Aperture LC analysis")') frequency_flcuv(i,s),
     &             flux_flcuv(i,s),FluxU_flcuv(i,s),FluxL_flcuv(i,s),mjdst_flcuv(i),
     &             mjded_flcuv(i),flcuv_flag(i,s),flcuv_type(i)
               enddo
            endif
         enddo
         do i=1,ivhe
            vhe_ref(i)=iref+1
            refs(iref+1)='TBD'
            do r=1,iref
               if (vhe_type(i) == name_cat(r)) THEN
                  vhe_ref(i)=r
               endif
            enddo
            if (filen_v(i) .eq. j) then
               if ((flux_vhe(i) .eq. FluxL_vhe(i)) .and. (flux_vhe(i) .eq. FluxU_vhe(i))) then
                  vhe_flag(i)=' UL '
               else if ((FluxL_vhe(i) .eq. 0.) .and. (FluxU_vhe(i) .ne. 0.) ) then
                  vhe_flag(i)=' UL '
               else
                  vhe_flag(i)=' Det'
               endif
               write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,a,2x,a)') frequency_vhe(i),flux_vhe(i),FluxU_vhe(i),
     &             FluxL_vhe(i),mjdstart(i),mjdend(i),vhe_flag(i),vhe_type(i),refs(vhe_ref(i))
            endif
         enddo
         do i=1,idebl
            debl=.false.
            do s=1,igam+ivhe
               if (i .eq. eblnn(s)) debl=.true.
            enddo
            !write(*,*) maxebl
c            if (i .gt. maxebl) debl=.true.
            if (debl) then
            do s=1,5
               if ((flux_debl(i,s) .eq. FluxL_debl(i,s)) .and. (flux_debl(i,s) .eq. FluxU_debl(i,s))) then
                  debl_flag(i,s)=' UL '
               else if ((FluxL_debl(i,s) .eq. 0.) .and. (FluxU_debl(i,s) .ne. 0.) ) then
                  debl_flag(i,s)=' UL '
               else
                  debl_flag(i,s)=' Det'
               endif
               !if (frequency_debl(i,s) .gt. 0.)
               write(14,'(4(es10.3,2x),2(f10.4,2x),a,2x,"DEBL        EBL-corrected flux")')
     &         frequency_debl(i,s),flux_debl(i,s),FluxU_debl(i,s),FluxL_debl(i,s),mjdavg,mjdavg,debl_flag(i,s)
            enddo
            endif
         enddo
c         write(*,*) '.......................source type and cataloged..................'
c         do i=1,icat
c            call Dist_sky(ra_source(j),dec_source(j),ra_cat(i),dec_cat(i),dist)
c            if (dist < min_dist_other) then
c               if (type_cat(i) .eq. -4) then
c                  write(*,*) "Warning!!Might be associated to a cluster."
c               else if (type_cat(i) .eq. -2) then
c                  write(*,*) "Known blazar!!"
c               else if (type_cat(i) .eq. -3) then
c                  write(*,*) "Flat radio spectrum source!!"
c               else
c                  write(*,*) "Already in 3HSP!!"
c               endif
c            endif
c         enddo
         write(*,*) '        '
         !if (j .ne. isource ) write(14,*) "===================="
      ENDDO

      deallocate(typer,ra_source,dec_source)
      deallocate(ra_cat,dec_cat)
      deallocate(type_cat,name_cat,refs)
      deallocate(filen_v,vhe_ref,ra_vhe,dec_vhe)
      deallocate(frequency_vhe,flux_vhe,FluxU_vhe,FluxL_vhe)
      deallocate(poserr_vhe,Ferr_vhe,mjdstart,mjdend)
      deallocate(vhe_flag,vhe_type)
      deallocate(veritind,magicind,eblnn,debl_flag)
      deallocate(flux_debl,FluxU_debl,FluxL_debl,frequency_debl)
      deallocate(filen_g,gampart,gam_ref,bigbind)
      deallocate(ra_gam,dec_gam,gamlike,poserr_gam)
      deallocate(frequency_gam,slope_gam,specerr_gam)
      deallocate(flux_gam,FluxL_gam,FluxU_gam,Ferr_gam)
      deallocate(gam_flag,namegam,name_g,gam_type)
      deallocate(ra_xray,dec_xray,xrtspind,xraypart)
      deallocate(slope_xray,poserr_xray,mjdst_xrt,mjded_xrt)
      deallocate(frequency_xray,flux_xray,FluxU_xray,FluxL_xray,Ferr_xray)
      deallocate(name_x,xray_type,filen_x)
      deallocate(xray_ref,xray_flag)
      deallocate(filen_u,uv_type,name_u)
      deallocate(ra_uv,dec_uv,poserr_uv)
      deallocate(uvmag,flux_uv,frequency_uv,FluxU_uv,FluxL_uv,uvmagerr)
      deallocate(uv_ref,uv_flag,uvcand_type)
      deallocate(ra_uvcand,dec_uvcand,uvdist,epos_uvcand)
      deallocate(flux_uvcand,uvmag_cand,freq_uvcand,uflux_uvcand,lflux_uvcand)
      deallocate(ra_flcuv,dec_flcuv,flcuvpart,filen_a)
      deallocate(mjdst_flcuv,mjded_flcuv,ts,duration,slope_flcuv,poserr_flcuv)
      deallocate(frequency_flcuv,flux_flcuv,FluxU_flcuv,FluxL_flcuv,Ferr_flcuv)
      deallocate(flcuv_flag,flcuv_type,name_a)
      deallocate(optpart,filen_o,indoptlc,optlc_ref,optlc_flag,opt_type,name_o)
      deallocate(ra_usno,dec_usno,poserr_usno,mjdst_optlc)
      deallocate(flux_usno,frequency_usno,FluxU_usno,FluxL_usno,usnomagerr,usnomag)
      deallocate(filen_i,ir_ref,indirlc,ra_ir,dec_ir,ra_ircand,dec_ircand)
      deallocate(irdist,poserr_ir,epos_ircand)
      deallocate(mjdst_irlc,mjded_irlc,mjdst_irlccand,mjded_irlccand)
      deallocate(flux_ir,irmag,frequency_ir,FluxU_ir,FluxL_ir,irmagerr)
      deallocate(flux_ircand,irmag_cand,freq_ircand,uflux_ircand,lflux_ircand)
      deallocate(flag_ir,ir_flag,ir_type,ircand_type,name_i,irlc_name)
      deallocate(filen_f,far_ref,farpart)
      deallocate(ra_far,dec_far,farlike,poserr_far)
      deallocate(frequency_far,flux_far,Ferr_far,FluxU_far,FluxL_far)
      deallocate(far_flag,name_f,far_type)
      deallocate(filen_p,pccspart,pccs100_ref,almaind,ra_pccs100,dec_pccs100)
      deallocate(mjdst_alma,mjded_alma,poserr_pccs100,pccslike)
      deallocate(flux_pccs100,frequency_pccs100,flux2_pccs100,snr_pccs100)
      deallocate(Ferr_pccs100,FluxU_pccs100,FluxL_pccs100,Ferr2_pccs100)
      deallocate(pccs100_flag,name_p,pccs100_type,date_alma)
      deallocate(filen_r,f4p8part,f4p8_ref,kuehrind)
      deallocate(ra_4p8,dec_4p8,f4p8like,poserr_4p8)
      deallocate(flux_4p8,frequency_4p8,Ferr_4p8,FluxU_4p8,FluxL_4p8)
      deallocate(flag_4p8,f4p8_flag,f4p8_type,name_r)
      deallocate(filen_l,lowr_ref,ra_lowr,dec_lowr)
      deallocate(frequency_lowr,flux_lowr,Ferr_lowr,FluxU_lowr,FluxL_lowr,poserr_lowr)
      deallocate(lowr_flag,name_l,lowr_type)
      deallocate(ra_lowrcand,ra_usnocand,ra_gamslp)
      deallocate(dec_lowrcand,dec_usnocand,dec_gamslp)
      deallocate(optdist,epos_usnocand,lowrdist,gamslp,poserr_gamslp)
      deallocate(freq_lowrcand,flux_lowrcand,uflux_lowrcand,lflux_lowrcand,epos_lowrcand)
      deallocate(freq_usnocand,flux_usnocand,uflux_usnocand,lflux_usnocand,usnomag_cand)
      deallocate(opt_ref,opt_flag,optcand_type,lowrcand_type)
      deallocate(mjdst_rrxx,mjded_rrxx,frequency,flux,uflux,lflux,epos)
      deallocate(spec_type,ra_rrxx,dec_rrxx,rrxx_type,recordmjd,npt,zsource)

      close(lu_output)
      close(14)
      END
ccccccccc
      SUBROUTINE DIST_SKY(alpha1,delta1,alpha2,delta2,dist)
      IMPLICIT NONE
      REAL*8 dist,alpha1,alpha2,delta1,delta2,costheta
      REAL*8 radian
      radian=57.2957795
      costheta=sin(delta1/radian)*sin(delta2/radian)+
     &         cos(delta1/radian)*cos(delta2/radian)*
     &         cos((alpha1-alpha2)/radian)
      if (costheta .gt. 1.) costheta=1.
      dist=acos(costheta)*radian
c      write(*,*) costheta,dist*3600.
      RETURN
      END
c
      SUBROUTINE graphic_code(intensity,band_type,code)
      IMPLICIT none 
      REAL*4 intensity,code,rfl_max,rfl_min,xfl_min,xfl_max,irfl_min,irfl_max
      REAL*4 uvfl_max,uvfl_min,gfl_min,gfl_max,optfl_max,optfl_min
      INTEGER*4 component,band_type
      rfl_min=0.8*1.4e9*1E-26 ! 0.8 mJy
      rfl_max=8000.*1.4e9*1E-26 ! 8 Jy
      xfl_min = 1.e-16 ! 1.e-16 erg/cm2/s, nufnu
      xfl_max = 5.e-11
      if ((band_type .gt. 10) .and. (band_type .lt. 20)) then
         component=int(alog10(intensity/rfl_min)/alog10(rfl_max/rfl_min)*99.)
      else if ((band_type .gt. 80) .and. (band_type .lt. 90)) then
         component=int(alog10(intensity/xfl_min)/alog10(xfl_max/xfl_min)*99.)
      else if (band_type .eq. 51) then
         irfl_max= 8
         irfl_min= 18
         component=int(alog10(irfl_min/intensity)/alog10(irfl_min/irfl_max)*99.)
      else if (band_type .eq. 52) then
         irfl_max= 7
         irfl_min= 17
         component=int(alog10(irfl_min/intensity)/alog10(irfl_min/irfl_max)*99.)
      else if (band_type .eq. 61) then
         optfl_max= 8
         optfl_min= 25
         component=int(alog10(optfl_min/intensity)/alog10(optfl_min/optfl_max)*99.)
      else if (band_type .eq. 62) then
         optfl_max= 8
         optfl_min= 24
         component=int(alog10(optfl_min/intensity)/alog10(optfl_min/optfl_max)*99.)
      else if (band_type .eq. 63 ) then
         optfl_max= 8
         optfl_min= 23
         component=int(alog10(optfl_min/intensity)/alog10(optfl_min/optfl_max)*99.)
      else if ((band_type .gt. 70) .and. (band_type .lt. 80)) then
         uvfl_max= 14
         uvfl_min= 24
         component=int(alog10(uvfl_min/intensity)/alog10(uvfl_min/uvfl_max)*99.)
      else if ((band_type .eq. 91) .or. (band_type .eq. 93)) then
         gfl_max= 1.E-9
         gfl_min= 1.e-13
         component=int(alog10(intensity/gfl_min)/alog10(gfl_max/gfl_min)*99.)
      else if (band_type .eq. 92) then
         gfl_max= 1.E-9
         gfl_min= 5.e-13
         component=int(alog10(intensity/gfl_min)/alog10(gfl_max/gfl_min)*99.)
      else if (band_type .eq. 94) then
         gfl_max= 8.e-10
         gfl_min= 5.e-12
         component=int(alog10(intensity/gfl_min)/alog10(gfl_max/gfl_min)*99.)
      endif
      IF (component .GE. 99) THEN
         component = 99
      ELSE IF (component .LE. 1) THEN
         component = 1
      ENDIF
      code = band_type*100.+component
      RETURN  
      END
cccc
      subroutine liklihood(poserr1,poserr2,dist,intensity,drop,like)
      integer drop
      real*4 poserr1,poserr2,intensity,like,density,mjy,ratio,sigma,pi
      real*8 dist
      mjy=intensity/(1.4e9*1.E-26)
      pi=3.1415926
      if (mjy .le. 100.) then
         density=10**(2.13654-0.9096*log10(mjy))
      else
         density=10**(4.6326-2.0377*log10(mjy))
      endif
      sigma=sqrt((poserr1*poserr1)+(poserr2*poserr2))
      ratio=2.*(real(dist)*3600./sigma)**drop
      like=exp(-ratio/2.)/(2.*pi*density*(sigma/3600.)**2)
      !write(*,*) 'LIKE',dist*3600.,sigma,drop,density,ratio,exp(-ratio/2.),(2.*pi*density*(sigma/3600.)**2)
      !write(*,*) "the likily hood",like
      end
c
      SUBROUTINE fluxtofdens(gamma,bandl,bandu,flux,gev,fdens,nudens)
      real*4 bandu,bandl,flux,nudens,fdens,conval
c      write(*,*) gamma,flux,gev,bandu,bandl
c      if (gamma .ne. 2.d0 ) then
      conval=(1./(-gamma+1.))*((bandu)**(-gamma+1.)-(bandl)**(-gamma+1.))
c      else
c      write(*,*) '2FHL'
c      conval=log10(bandu/bandl)
c      endif
      fdens=gev*(flux/conval)*((gev)**(-gamma))!!!!To photon flux at gev
      fdens=fdens*1.602E-19*1.E7*(gev*1.E9)
      nudens=(1.602E-19)*(gev*1.E9)/(6.626e-34)
      RETURN
      end
c
      SUBROUTINE fluxtofdens2(alpha,bandl,bandu,flux,kev,fdens,nudens)
      real*4 alpha,bandu,bandl,flux,nudens,fdens,conval,kev!,nuu,nul
      !write(*,*) alpha,flux,kev,bandu,bandl
      if (alpha .ne. 1. ) then
      conval=(1./(-alpha+1.))*((bandu)**(-alpha+1.)-(bandl)**(-alpha+1.))
      else
      conval=log10(bandu/bandl)
      endif
      !write(*,*) kev,(1./conval)*kev*kev**(-alpha)!/(kev*2.418E-12)
      fdens=kev*(flux/conval)*((kev)**(-alpha))!!!!
      nudens=(1.602E-19)*(kev*1.e3)/(6.626e-34)
      RETURN
      end
cccccc
      subroutine date_to_mjd(year,month,date,hour,minute,second,mjd)
      IMPLICIT NONE
      integer*4 year,month,date,hour,minute,second,leap
      real*8 mjd
      leap = (month-14)/12        !In leap years, -1 for Jan, Feb, else 0
      mjd = date-32075+1461*(year+4800+leap)/4+367*(month-2-leap*12)/12-3*((year+4900+leap)/100)/4
      mjd = mjd + (float(hour)/24.) + (float(minute)/24./60.) + (float(second)/24./3600.) - 0.5 -2400000.5
      return
      end
ccccccc
      subroutine mjd_to_date(mjd,year,month,date,hour)
      IMPLICIT NONE
      integer*4 year,month,date,hour,l,n
      integer*8 intjd
      real*8 mjd,frac
      mjd=mjd+2400000.5
      intjd=int8(mjd)
      frac = mjd - intjd + 0.5          !Fractional part of calendar day
      if (frac .gt. 1.0) then
          frac=frac-1.
          intjd=intjd+1
      endif
      hour = int(real(frac))*24.0
      l = int(intjd) + 68569
      n = 4*l / 146097
      l = l - (146097*n + 3) / 4
      year = 4000*(l+1) / 1461001
      l = l - 1461*year / 4 + 31        !1461 = 365.25 * 4
      month = 80*l / 2447
      date = l - 2447*month / 80
      l = month/11
      month = month + 2 - 12*l
      year = 100*(n-49) + year + l
c      write(*,*) year,month,date,hour
      return
      end

c
        SUBROUTINE mag2flux (nh,m_band,filter,flux,frequency)
c
c  converts u,v,i,h,b,r,j,k magnitudes into monochromatic fluxes
c  in units of erg/cm2/s for nufnu vs nu plots 
c
        IMPLICIT none
        REAL*4 nh, flux , av , m_band, a_band, Rv 
        REAL*4 c, frequency, a
        REAL*8 x,aa,bb,c1,c2,dx,px,ebv,y,const,lambda
        CHARACTER*3 filter
c        print *,' nh, m_band, filter ', nh, m_band,filter
c        call upc(filter)
        IF ( (filter(1:3).NE.'U  ') .AND. (filter(1:3).NE.'B  ') .AND.
     &       (filter(1:3).NE.'V  ') .AND. (filter(1:3).NE.'R  ') .AND.
     &       (filter(1:3).NE.'I  ') .AND. (filter(1:3).NE.'J  ') .AND.
     &       (filter(1:3).NE.'H  ') .AND. (filter(1:3).NE.'K  ') .AND.
     &       (filter(1:3).NE.'u  ') .AND. (filter(1:3).NE.'g  ') .AND.
     &       (filter(1:3).NE.'r  ') .AND. (filter(1:3).NE.'i  ') .AND.
     &       (filter(1:3).NE.'z  ') .AND. (filter(1:3).NE.'psg') .AND.
     &       (filter(1:3).NE.'psr') .AND. (filter(1:3).NE.'psi') .AND.
     &       (filter(1:3).NE.'psz') .AND. (filter(1:3).NE.'psy') .AND.
     &       (filter(1:3).NE.'fuv') .and. (filter(1:3).NE.'nuv') .AND.
     &       (filter(1:3).NE.'su ') .and. (filter(1:3).NE.'sv ') .AND.
     &       (filter(1:3).NE.'sb ') .and. (filter(1:3).NE.'sw1') .AND.
     &       (filter(1:3).NE.'sm2') .and. (filter(1:3).NE.'sw2') .AND.
     &       (filter(1:3).NE.'xu ') .and. (filter(1:3).NE.'xv ') .AND.
     &       (filter(1:3).NE.'xb ') .and. (filter(1:3).NE.'xw1') .AND.
     &       (filter(1:3).NE.'xm2') .and. (filter(1:3).NE.'xw2') .AND.
     &       (filter(1:3).NE.'ww1') .and. (filter(1:3).NE.'ww2') .AND.
     &       (filter(1:3).NE.'ww3') .and. (filter(1:3).NE.'ww4') .AND.
     &       (filter(1:3).NE.'bbG') ) THEN
         write (*,*) ' mag2flux: Filter not supported  '
         stop
        ENDIF
c extintion law taken from Cardelli et al. 1989 ApJ 345, 245
c in UV apply the UV relation from Fitzpatrick 1999
        Rv=3.1
        av = Rv*(-0.055+nh*1.987e-22) !!dust map from BH1978, assumed constant gas-to-dust ratio
ccccccccccc check nh=2.29e21*Av
        ebv=av/Rv
        if (av < 0.) av=0.
        if (filter(1:3) == 'U  ') then
           lambda=3600.d0
           const=log10(1810.d0)-23.d0
        else if (filter(1:3) == 'B  ') then
           lambda=4400.d0
           const=log10(4260.d0)-23.d0
        else if (filter(1:3) == 'V  ') then
           lambda=5500.d0
           const=log10(3640.d0)-23.d0
        else if (filter(1:3) == 'R  ') then
           lambda=6400.d0
           const=log10(3080.d0)-23.d0
        else if (filter(1:3) == 'I  ') then
           lambda=7900.d0
           const=log10(2550.d0)-23.d0
        else if (filter(1:3) == 'J  ') then !2MASS
           lambda=12350.d0
           const=log10(1594.d0)-23.d0
        else if (filter(1:3) == 'H  ') then
           lambda=16620.d0
           const=log10(1024.d0)-23.d0
        else if (filter(1:3) == 'K  ') then
           lambda=21590.d0
           const=log10(666.7d0)-23.d0
        else if (filter(1:3) == 'u  ') then !effective wavelength from SDSS, Doi et al. 2010 ApJ 139, 1628
           !m_band=m_band-0.04 !calibrate of the SDSS u band to AB mag
           lambda=3568.d0
           const=log10(3631.d0)-23.d0 !3631 is the 0 mag flux of AB mag system
        else if (filter(1:3) == 'g  ') then
           lambda=4653.d0
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'r  ') then !effective wavelength from SDSS, Doi et al. 2010 ApJ 139, 1628
           lambda=6148.d0
           const=log10(3631.d0)-23.d0 !3631 is the 0 mag flux of AB mag system
        else if (filter(1:3) == 'i  ') then
           lambda=7468.d0
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'z  ') then
           lambda=8863.d0
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'psg') then !tonry et al. 2012
           lambda=4810.d0
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'psr') then
           lambda=6170.d0
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'psi') then
           lambda=7520.d0
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'psz') then
           lambda=8660.d0
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'psy') then
           lambda=9620.d0
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'bbG') then
           lambda=6730.d0 !!!Jordi et al. 2010
           const=log10(2918.d0)-23.d0 !!!!the zero mag. flux are estimated from Vega flux!!!
        else if (filter(1:3) == 'fuv') then
           lambda=1538.6d0
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'nuv') then
           lambda=2315.7d0
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'su ') then
           lambda=3501.d0 !from Poole et al. (2008) effective wavelength
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'sb ') then
           lambda=4329.d0 !from Poole et al. (2008) effective wavelength
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'sv ') then
           lambda=5402.d0 !from Poole et al. (2008) effective wavelength
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'sw1') then
           lambda=2634.d0 !from Poole et al. (2008) effective wavelength
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'sm2') then
           lambda=2231.d0 !from Poole et al. (2008) effective wavelength
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'sw2') then
           lambda=2030. !from Poole et al. (2008) effective wavelength
           const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'xu ') then
        lambda=3440. !from Page et al. (2008) effective wavelength
        const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'xb ') then
        lambda=4500.d0 !from Page et al. (2011) effective wavelength
        const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'xv ') then
        lambda=5430.d0 !from Page et al. (2011) effective wavelength
        const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'xw1') then
        lambda=2910. !from Page et al. (2011) effective wavelength
        const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'xm2') then
        lambda=2310.d0 !from Page et al. (2011) effective wavelength
        const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'xw2') then
        lambda=2120. !from Page et al. (2011) effective wavelength
        const=log10(3631.d0)-23.d0
        else if (filter(1:3) == 'ww1') then
           lambda=34000.d0
           const=log10(309.540d0)-23.d0
        else if (filter(1:3) == 'ww2') then
           lambda=46000.d0
           const=log10(171.787d0)-23.d0
        else if (filter(1:3) == 'ww3') then
           lambda=120000.d0
           const=log10(31.674d0)-23.d0
        else if (filter(1:3) == 'ww4') then
           lambda=220000.d0
           const=log10(8.363d0)-23.d0
        endif
c lambda from Amstrongs to microns
        !write(*,*) filter(1:3),m_band,a_band,const
        x=(10000.d0/lambda)
        if ((x .lt. 1.1d0) .and. (x .ge. 0.3d0)) then
        a_band=(0.574*(real(x)**1.61)-0.527*(real(x)**1.61)/Rv)*av ! the a_lambda
        else if ((x .lt. 3.3d0) .and. (x .ge. 1.1d0)) then
        y=x-1.82d0
        aa=1+(0.17699d0*y)-(0.50447d0*y**2d0)-(0.02427d0*y**3d0)+(0.72085d0*y**4d0)
     &  +(0.01979d0*y**5d0)-(0.77530d0*y**6d0)+(0.32999d0*y**7d0)
        bb=1.41338d0*y+(2.28305d0*y**2d0)+(1.07233d0*y**3d0)-(5.38434d0*y**4d0)
     &  -(0.62251d0*y**5d0)+(5.30260d0*y**6d0)-(2.09002d0*y**7d0)
c        write(*,*) aa,bb,av,Rv
        a_band=(real(aa)+(real(bb)/Rv))*av
        else if ((x .le. 10.d0) .and. (x .ge. 3.3d0)) then
        c2=-0.824+(4.717/Rv)
        c1=2.03-(3.007*c2)
        dx=(x*x)/((((x**2)-(4.596**2))**2)+((x*0.99)**2))
        px=(0.5392*((x-5.9)**2))+(0.05644*((x-5.9)**3))
        if (x .le. 5.9d0) px=0.
        !write(*,*) (c1+c2*x+3.23*dx+0.41*px),c1,c2,dx,px
        a_band=real((c1+c2*x+3.23*dx+0.41*px))*real(ebv)+av
        else
        a_band=0.
        endif
        c=2.9979e10
        a=1.0
        frequency=c/(real(lambda)*1.e-8)
        flux = 10.**(-0.4*(m_band -a_band)+real(const))*frequency
        !if (filter(1:3)=='xw1') write(*,*) 'W1',m_band,a_band,flux
        if (m_band .eq. 0.) flux=0.
        !m_band=0.
        !lambda=0.
        RETURN 
        END

