      PROGRAM aladin_interface
c
c
      IMPLICIT none
      INTEGER*4 ier, lu_in, lu_in_ra, lu_script, lenact,i,ia,ii
      INTEGER*4 lu_in2, lu_out,source_number,length,in,im
      REAL*4 dummy, fov
      REAL*8 ra, dec,ra_center,dec_center,ra2,dec2
      CHARACTER*1 a0,a1,a2
      CHARACTER*200 input_file,input_file2,input_file_ra,output_file,string,script_file
      character*800 stringin
      LOGICAL there,ok


      CALL rdforn(stringin,length)
c      IF ( length.NE.0 ) THEN
      CALL rmvlbk(stringin)
c         write(*,*) string,length
c      write(*,*) stringin,length
      in=index(stringin(1:length),' ')
      input_file_ra=stringin(1:in-1)
      im=index(stringin(in+1:length),' ')+in
c         write(*,*) in,im
      input_file=stringin(in+1:im-1)
      in=im
      im=index(stringin(in+1:length),' ')+in
      input_file2=stringin(in+1:im-1)
c         in=im
c         im=index(string(in+1:length),' ')+in
c         input_file4=string(in+1:im-1)
      in=im
      im=index(stringin(in+1:length),' ')+in
      script_file=stringin(in+1:im-1)
c      write(*,*) in,im
      in=im
      im=index(stringin(in+1:length),' ')+in
c      write(*,*) in,im
      output_file=stringin(in+1:im-1)
c      write(*,*) output_file
      in=im
      read(stringin(in+1:length),*) fov

c      write(*,*) output_file
c      write(*,*) script_file

c      script_file   = '/Users/paologiommi/app/aladin_script.js'
c      input_file_ra ='output1.csv'
c      input_file    ='find_out_temp.txt'
c      input_file2    ='candidates_posix.txt'
c      output_file='vou-aladin.html'

      ok = .TRUE.
      lu_in = 10
      lu_in2 = 11
      lu_in_ra = 12
      lu_script = 13
      lu_out = 21
      a1 = "'"
      a2 = "}"
      a0 = "{"
c      CALL rdforn(string,length)
c      IF ( length.NE.0 ) then
c          CALL rmvlbk(string)
c          READ(string(1:lenact(string)),*)
c     &                fov
c      ENDIF
      INQUIRE (FILE=input_file,EXIST=there)
      IF (.NOT.there) THEN
         write (*,'('' file '',a,'' not found '')')
     &     input_file(1:lenact(input_file))
         STOP
      ENDIF
      open(lu_script,file=script_file,status='old',iostat=ier)
      open(lu_in_ra,file=input_file_ra,status='old',iostat=ier)
      read (lu_in_ra,'(a)') string
      read(string(4:13),*)  ra_center
      read(string(20:29),*) dec_center
      close (lu_in_ra)
      open(lu_in,file=input_file,status='old',iostat=ier)
      open(lu_in2,file=input_file2,status='old',iostat=ier)
      IF (ier.ne.0) THEN
        write (*,*) ' Error ',ier,' opening file ', input_file2
      ENDIF
      open(lu_out,file=output_file,status='unknown',iostat=ier)
      write(lu_out,'(''<HTML>'')')
      write(lu_out,'(''<HEAD>'')')
      write(lu_out,'(''<img src=http://openuniverse.asi.it/images_ou/OU-logo.png height="40" width="213"> '')')
      write(lu_out,'(3x,''<script type="text/javascript" src="https://code.jquery.com/jquery-1.10.1.min.js"></script>'')')
      write(lu_out,'(3x,''<link rel="stylesheet" href="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.css" >'')')
      write(lu_out,'(3x,''<script type="text/javascript">var jqMenu = jQuery.noConflict();</script>'')')
      write(lu_out,'(''</HEAD>'')')
c      write(lu_out,'(3x,''<script type="text/javascript">'')')
c      write(lu_out,'(3x,''var hipsDir=null;</script>'')')
c      write(lu_out,'(''<H1>VOU-Blazars - Aladin interface</H1>'')')
c      write(lu_out,'(''<script type="text/javascript">'')')
c      write(lu_out,'(''hipsDir = location.href;'')')
c      write(lu_out,'(''hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length))+''
c     &                 ,a,''/HIPS_GRBSurvey'',a,'';'')') a1,a1
c      write(lu_out,'(''document.getElementById("hipsBase").innerHTML=hipsDir;'')')
c      write(lu_out,'(''</script>'')')
      write(lu_out,'(''<script type="text/javascript" src="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.js"''
     &               '' charset="utf-8"></script>'')')
       write(lu_out,'(''<fieldset style="margin:15px 15px 30px; padding-bottom:30px; '',
     &                ''padding-top:15px; padding-left:15px; width:auto; margin-bottom:10px;">'')')
       write(lu_out,'(''<legend class='',a,''title'',a,''><b> VOU-Blazars candidates overlay.</b>'',
     &                ''&emsp;&emsp;  Based on "Aladin Lite" developed at CDS, Strasbourg Observatory, France<br></legend>'')')
     &                 a1,a1
      write(lu_out,'(''<div id="infoDiv" style="margin-left:120px" >&nbsp; </div>'')')
      write(lu_out,'(''<div id="aladin-lite-div" style="margin-left:120px;width:80vw;height:80vh;"></div>'')')
      write(lu_out,'(''<script type="text/javascript">'')')
      i = 0
      DO WHILE(ok)
         READ(lu_in,*,end=98) ra,dec,ii
         ia = abs(ii)/10000
         IF ( ( ia > 0 ) .AND. ( (abs(ii)-ia*10000) > 0) ) THEN 
           i = i +1 
           IF (i == 1) THEN 
              IF (ra_center > 99.999) THEN 
              write(lu_out,'(''var aladin = $.aladin("#aladin-lite-div", {survey: '',a,''P/DSS/colored''
     &                  ,a,'', showSimbadPointerControl: true,''
     &                  ''fov:'',f6.3,'',target:'',a,f8.4,1x,f9.4,a,a,'');'')')a1,a1,fov*2./60.,a1,ra_center,dec_center,a1,a2
              ELSE IF (ra_center > 9.999) THEN 
              write(lu_out,'(''var aladin = $.aladin("#aladin-lite-div", {survey: '',a,''P/DSS/colored''
     &                  ,a,'', showSimbadPointerControl: true,''
     &                  ''fov:'',f6.3,'',target:'',a,f7.4,1x,f9.4,a,a,'');'')')a1,a1,fov*2./60.,a1,ra_center,dec_center,a1,a2
              ELSE
              write(lu_out,'(''var aladin = $.aladin("#aladin-lite-div", {survey: '',a,''P/DSS/colored''
     &                  ,a,'', showSimbadPointerControl: true,''
     &                  ''fov:'',f6.3,'',target:'',a,f9.4,1x,f9.4,a,a,'');'')')a1,a1,fov*2./60.,a1,ra_center,dec_center,a1,a2
              ENDIF
              DO WHILE (ok) 
                 read (lu_script,'(a)',end=97) string
                 write(lu_out,'(a)') string(1:lenact(string))
              ENDDO
97            CONTINUE
c              write(lu_out,'(''aladin.createImageSurvey('',a,
c     &                   ''Swift XRT @ OpenUniverse'',a,'', '',a,''Swift XRT @ OpenUniverse'',a,'','')') a1,a1,a1,a1
c              write(lu_out,'(''hipsDir, '',a,''equatorial'',a,'', 7, {imgFormat: '',a,''png'',a,''});'')') a1,a1,a1,a1
              write(lu_out,'(''var overlay = A.graphicOverlay({color: '',a,''cyan'',a,'', lineWidth: 5});'')') a1,a1
              write(lu_out,'(''aladin.addOverlay(overlay);'')')
           ENDIF
           IF (ia == 1) THEN
              write(lu_out,'(''overlay.add(A.circle('',f9.4,1x,'','',f9.4,'', .0027, ''
     &                    ,a,''color: '',a,''orangered'',a,a,'' ));'')') ra,dec,a0,a1,a1,a2
           ELSE IF (ia == 2) THEN
              write(lu_out,'(''overlay.add(A.circle('',f9.4,1x,'','',f9.4,'', .0027, ''
     &                    ,a,''color: '',a,''aqua'',a,a,'' ));'')') ra,dec,a0,a1,a1,a2
           ELSE IF (ia == 3) THEN
              IF ( ii > 0 ) THEN 
                 write(lu_out,'(''overlay.add(A.circle('',f9.4,1x,'','',f9.4,'', .0027, ''
     &                    ,a,''color: '',a,''#1a1aff'',a,a,'' ));'')') ra,dec,a0,a1,a1,a2
              ELSE
                 write(lu_out,'(''overlay.add(A.circle('',f9.4,1x,'','',f9.4,'', .0027, ''
     &                    ,a,''color: '',a,''yellowgreen'',a,a,'' ));'')') ra,dec,a0,a1,a1,a2
              ENDIF
           ELSE IF ( (ia == 8) .OR. (ia == 9) ) THEN
              write(lu_out,'(''overlay.add(A.circle('',f9.4,1x,'','',f9.4,'', .0027, ''
     &                    ,a,''color: '',a,''silver'',a,a,'' ));'')') ra,dec,a0,a1,a1,a2
           ENDIF
         ENDIF
      ENDDO 
 98   CONTINUE
      DO WHILE(ok)
         READ(lu_in2,*,end=99) ra2,dec2,dummy,dummy,source_number
         write(lu_out,'(''cat.addSources([A.source('',f9.4,1x,'','',f9.4,'','',
     &         a,''name: '',a,'' Nr. '',i3,a,a,'')]);'')') ra2,dec2,a0,a1,source_number,a1,a2
      ENDDO
 99   CONTINUE
      write(lu_out,'(''</script>'')')
      write(lu_out,'(''<script type="text/javascript">'')')
      write(lu_out,'(''document.getElementById("hipsBase").innerHTML=hips'')')
      write(lu_out,'(''</script>'')')
      write(lu_out,'(''</HTML>'')')
      END
