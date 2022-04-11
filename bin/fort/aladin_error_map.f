      PROGRAM aladin_error_map
c
c
      IMPLICIT none
      INTEGER*4 ier, lu_in, lu_in_ra, lu_script, lenact,i,ia,ii
      INTEGER*4 lu_in2, lu_out,source_number,im,in,length
      REAL*4 dummy, radius
      REAL*8 ra, dec,ra_center,dec_center,ra2,dec2
      CHARACTER*1 a0,a1,a2
      CHARACTER*200 input_file,input_file2,input_file_ra,output_file,string,script_file
      character*800 stringin
      LOGICAL there,ok

      CALL rdforn(stringin,length)
      CALL rmvlbk(stringin)
      in=index(stringin(1:length),' ')
      input_file_ra=stringin(1:in-1)
      im=index(stringin(in+1:length),' ')+in
      input_file=stringin(in+1:im-1)
      in=im
      im=index(stringin(in+1:length),' ')+in
      script_file=stringin(in+1:im-1)
      in=im
      im=index(stringin(in+1:length),' ')+in
      output_file=stringin(in+1:length)

c      input_file_ra ='phase2'
c      input_file    ='error_map.txt'
c      script_file   = '~/app/aladin_script.js'
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
      INQUIRE (FILE=input_file,EXIST=there)
      IF (.NOT.there) THEN
         write (*,'('' file '',a,'' not found '')')
     &     input_file(1:lenact(input_file))
         STOP
      ENDIF
      open(lu_script,file=script_file,status='old',iostat=ier)
      open(lu_in_ra,file=input_file_ra,status='old',iostat=ier)
      DO WHILE (ok) 
        read(lu_in_ra,'(a)') string 
        IF (string(11:16) == 'Dec. =') goto 199
      ENDDO  
 199  CONTINUE
      read(string(17:28),*) ra_center
      read(string(30:38),*) dec_center
      close (lu_in_ra)
      open(lu_in,file=input_file,status='old',iostat=ier)
      IF (ier.ne.0) THEN
        write (*,*) ' Error ',ier,' opening file ', input_file2
      ENDIF
      open(lu_out,file=output_file,status='unknown',iostat=ier)
      write(lu_out,'(''<HTML>'')')
      write(lu_out,'(''<HEAD>'')')
      write(lu_out,'(3x,''<script type="text/javascript" src="https://code.jquery.com/jquery-1.12.1.min.js" charset="utf-8"></script>'')')
      write(lu_out,'(3x,''<link rel="stylesheet" href="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.css" >'')')
      write(lu_out,'(3x,''<script type="text/javascript">var jqMenu = jQuery.noConflict();</script>'')')
      write(lu_out,'(3x,''<script type="text/javascript">'')')
      write(lu_out,'(3x,''var hipsDir=null;</script>'')')
      write(lu_out,'(''</HEAD>'')')
      write(lu_out,'(''<img src=http://openuniverse.asi.it/images_ou/OU-logo.png height="40" width="213"> '')')
c      write(lu_out,'(''<H1>VOU-Blazars/Error circles map - Aladin interface</H1>'')')
c      write(lu_out,'(''<script type="text/javascript">'')')
c      write(lu_out,'(''hipsDir = location.href;'')')
c      write(lu_out,'(''hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length))+''
c     &                 ,a,''/HIPS_GRBSurvey'',a,'';'')') a1,a1
c      write(lu_out,'(''document.getElementById("hipsBase").innerHTML=hipsDir;'')')
c      write(lu_out,'(''</script>'')')
      write(lu_out,'(''<script type="text/javascript" src="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.js"''
     &               '' charset="utf-8"></script>'')')
       write(lu_out,'(''<fieldset style="margin:50px 100px 30px; padding-bottom:30px; '',
     &                ''padding-top:15px; padding-left:15px; width:82vw; margin-bottom:10px;">'')')
       write(lu_out,'(''<legend class='',a,''title'',a,''><b>VOU-Blazars - SED error circles map</b>&emsp;'',
     &                ''&emsp;-- Based on "Aladin Lite" developed at CDS, Strasbourg Observatory, France<br></legend>'')')
     &                 a1,a1
      write(lu_out,'(''<div id="aladin-lite-div" style="margin-left:20px; width:80vw; height:70vh;"></div>'')')
      write(lu_out,'(''&nbsp; &nbsp; <input id="allwise" type="radio" name="survey" value="P/allWISE/color"><label for="allwise">AllWISE<label>'')')
      write(lu_out,'(''<input id="DSS" type="radio" name="survey" value="P/DSS2/color"><label for="DSS">DSS color<label>'')')
      write(lu_out,'(''<input id="PANSTARRS" type="radio" name="survey" value="P/PanSTARRS/DR1/color-z-zg-g"><label for="P/PanSTARRS/DR1/color">PanSTARRS DR1 color<label>'')')
      write(lu_out,'(''<input id="SDSS9" type="radio" name="survey" value="P/SDSS9/color"><label for="SDSS9">SDSS9 color<label>'')')
      write(lu_out,'(''<input id="GALEX" type="radio" name="survey" value="P/GALEXGR6/AIS/color"><label for="GALEX">GALEXGR6 color<label>'')')
      write(lu_out,'(''<input id="XMM-PN" type="radio" name="survey" value="P/XMM/PN/color"><label for="XMM-ON">XMM-PN color<label>'')')

      write(lu_out,'(''</fieldset>'')')
      write(lu_out,'(''<script type="text/javascript">'')')
      i = 0
      DO WHILE(ok)
         READ(lu_in,*,end=98) ra,dec,ii,radius
         radius= radius/3600.
         ia = abs(ii)/1000
         i = i +1 
         IF (i == 1) THEN 
            IF (ra_center > 99.999) THEN 
            write(lu_out,'(''var aladin = $.aladin("#aladin-lite-div", {survey: '',a,''P/DSS/colored''
     &                  ,a,'', showSimbadPointerControl: true,''
     &                  ''fov:0.05,target:'',a,f8.4,1x,f9.4,a,a,'');'')')a1,a1,a1,ra_center,dec_center,a1,a2
            ELSE IF (ra_center > 9.999) THEN 
            write(lu_out,'(''var aladin = $.aladin("#aladin-lite-div", {survey: '',a,''P/DSS/colored''
     &                  ,a,'', showSimbadPointerControl: true,''
     &                  ''fov:0.05,target:'',a,f7.4,1x,f9.4,a,a,'');'')')a1,a1,a1,ra_center,dec_center,a1,a2
            ELSE
            write(lu_out,'(''var aladin = $.aladin("#aladin-lite-div", {survey: '',a,''P/DSS/colored''
     &                  ,a,'', showSimbadPointerControl: true,''
     &                  ''fov:0.05,target:'',a,f6.4,1x,f9.4,a,a,'');'')')a1,a1,a1,ra_center,dec_center,a1,a2
            ENDIF
            DO WHILE (ok) 
                 read (lu_script,'(a)',end=97) string
                 write(lu_out,'(a)') string(1:lenact(string))
            ENDDO
97          CONTINUE
c            write(lu_out,'(''aladin.createImageSurvey('',a,
c     &                   ''Swift XRT @ OpenUniverse'',a,'', '',a,''Swift XRT @ OpenUniverse'',a,'','')') a1,a1,a1,a1
c            write(lu_out,'(''hipsDir, '',a,''equatorial'',a,'', 7, {imgFormat: '',a,''png'',a,''});'')') a1,a1,a1,a1
            write(lu_out,'(''var overlay = A.graphicOverlay({color: '',a,''cyan'',a,'', lineWidth: 3});'')') a1,a1
            write(lu_out,'(''aladin.addOverlay(overlay);'')')
         ENDIF
         IF (ia == 1) THEN
              write(lu_out,'(''overlay.add(A.circle('',f9.4,1x,'','',f9.4,'', '',f9.5,'', ''
     &                    ,a,''color: '',a,''red'',a,a,'' ));'')') ra,dec,radius,a0,a1,a1,a2
         ELSE IF (ia == 5) THEN
            write(lu_out,'(''overlay.add(A.circle('',f9.4,1x,'','',f9.4,'', '',f9.5,'', ''
     &                    ,a,''color: '',a,''orange'',a,a,'' ));'')') ra,dec,radius,a0,a1,a1,a2
c     &                    ,a,''color: '',a,''aqua'',a,a,'' ));'')') ra,dec,radius,a0,a1,a1,a2
         ELSE IF (ia == 6) THEN
            write(lu_out,'(''overlay.add(A.circle('',f9.4,1x,'','',f9.4,'', '',f9.5,'', ''
     &                    ,a,''color: '',a,''silver'',a,a,'' ));'')') ra,dec,radius,a0,a1,a1,a2
         ELSE IF (ia == 7) THEN
                 write(lu_out,'(''overlay.add(A.circle('',f9.4,1x,'','',f9.4,'', '',f9.5,'', ''
     &                    ,a,''color: '',a,''aquamarine'',a,a,'' ));'')') ra,dec,radius,a0,a1,a1,a2
         ELSE IF (ia == 8) THEN
            write(lu_out,'(''overlay.add(A.circle('',f9.4,1x,'','',f9.4,'', '',f9.5,'', ''
     &                    ,a,''color: '',a,''#1a1aff'',a,a,'' ));'')') ra,dec,radius,a0,a1,a1,a2
         ELSE IF (ia == 9) THEN
            write(lu_out,'(''overlay.add(A.circle('',f9.4,1x,'','',f9.4,'', '',f9.5,'', ''
     &                    ,a,''color: '',a,''BlueViolet'',a,a,'' ));'')') ra,dec,radius,a0,a1,a1,a2
         ENDIF
      ENDDO 

 98   CONTINUE
      DO WHILE(ok)
         READ(lu_in2,*,end=99) ra2,dec2,dummy,dummy,source_number
         write(lu_out,'(''cat.addSources([A.source('',f9.4,1x,'','',f9.4,'','',
     &         a,''name: '',a,'' Nr. '',i3,a,a,'')]);'')') ra2,dec2,a0,a1,source_number,a1,a2
      ENDDO
 99   CONTINUE

      write(lu_out,'(''$("input[name=survey]").change(function() {'')')
      write(lu_out,'(''   aladin.setImageSurvey($(this).val());'')')
      write(lu_out,'(''});'')')

      write(lu_out,'(''</script>'')')
      write(lu_out,'(''<script type="text/javascript">'')')
      write(lu_out,'(''document.getElementById("hipsBase").innerHTML=hips'')')
      write(lu_out,'(''</script>'')')
      write(lu_out,'(''</HTML>'')')
      END
