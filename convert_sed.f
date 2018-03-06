      PROGRAM convert_sed
c
c converts the format of a SED file produced with the find_candidates tool
c to the format required by the ASDC SED tool  
c
c
      IMPLICIT none
      INTEGER*4 ier, lu_in, lu_out, lenact
      REAL*4 mjd, freq, one, err_up,err_lo
      REAL*4 flux,flux_err
      CHARACTER*2 ul
      CHARACTER*80 input_file,output_file,string
      LOGICAL there,ok
      ok = .TRUE. 
      one = 1.0
      !WRITE (*,'('' Enter file name '',$)')
      !READ (*,'(a)') input_file
      !WRITE (*,'('' Enter output file '',$)')
      !READ (*,'(a)') output_file
      input_file='Sed.txt'
      output_file='Out4SedTool.txt'
      mjd = 55000
      lu_in = 10
      lu_out = 11
      INQUIRE (FILE=input_file,EXIST=there)
      IF (.NOT.there) THEN
         write (*,'('' file '',a,'' not found '')')
     &     input_file(1:lenact(input_file))
         STOP
      ENDIF
      open(lu_in,file=input_file,status='old',iostat=ier)
      open(lu_out,file=output_file,status='unknown',iostat=ier)
      IF (ier.ne.0) THEN
        write (*,*) ' Error ',ier,' opening file ', input_file
      ENDIF
      READ(lu_in,'(a)',end=99) string 
      DO WHILE(ok)
         ul = '  '
         READ(lu_in,'(a)',end=99) string 
c         print *,'string ',string(1:lenact(string))
         IF (string(2:2).NE.'=') THEN 
           READ(string(1:lenact(string)),*) freq, flux, err_up, err_lo 
           flux_err = (err_up-err_lo)/2.
           IF (flux.NE.0.) THEN 
             IF (flux_err == 0.) ul='UL'
             write(lu_out,'(e10.4,'' | '',e10.2,'' | '',e10.4,'' | '',
     &                      e10.4,'' | '',f10.2,'' | '',f10.2,'' |'',1x,a,''|'')')
     &                      freq,one,flux,flux_err,mjd,mjd,ul
           ENDIF
         ENDIF
      ENDDO 
 99   CONTINUE
      END
