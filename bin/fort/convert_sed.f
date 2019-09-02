      PROGRAM convert_sed
c
c converts the format of a SED file produced with the find_candidates tool
c to the format required by the ASDC SED tool  
c
c
      IMPLICIT none
      INTEGER*4 ier, lu_in, lu_out, lenact,sfound,rtype,im,in,length
      REAL*4 mjdstart,mjdend, freq, one, err_up,err_lo
      REAL*4 flux,flux_err
      real*8 rra,rdec
      CHARACTER*2 ul
      character*14 stringin
      CHARACTER*80 input_file,output_file
      Character*150 string
      LOGICAL there,ok
      ok = .TRUE. 
      one = 1.0
      !WRITE (*,'('' Enter file name '',$)')
      !READ (*,'(a)') input_file
      !WRITE (*,'('' Enter output file '',$)')
      !READ (*,'(a)') output_file

      call rdforn(string,length)
      call rmvlbk(string)
      in=index(string(1:length),' ')
      input_file=string(1:in-1)
      output_file=string(in+1:length)
c      write(*,*) input_file
c      write(*,*) output_file

      lu_in = 10
      lu_out = 11
      INQUIRE (FILE=input_file,EXIST=there)
      IF (.NOT.there) THEN
         write (*,'('' file '',a,'' not found '')')
     &     input_file(1:lenact(input_file))
         STOP
      ENDIF
      open(lu_in,file=input_file,status='old',iostat=ier)
      IF (ier.ne.0) THEN
        write (*,*) ' Error ',ier,' opening file ', input_file
      ENDIF
      READ(lu_in,'(i4,2x,a,2(2x,f9.5),2x,i2)',end=99) sfound,stringin,rra,rdec,rtype
c      write(*,*) rra,rdec
      read(lu_in,*) string
      read(lu_in,*) string
      read(lu_in,*) string
      open(lu_out,file=output_file,status='unknown',iostat=ier)
      DO WHILE(ok)
         ul = '  '
         READ(lu_in,'(a)',end=99) string 
c         print *,'string ',string(1:lenact(string))
         IF (string(2:2).NE.'=') THEN 
           READ(string(1:lenact(string)),*) freq, flux, err_up, err_lo, mjdstart,mjdend
           flux_err = (err_up-err_lo)/2.
           IF ((flux .NE. 0.) .or. (err_up .ne. 0.)) THEN
             if (flux .lt. 0.) flux = -flux
             IF ((flux_err == 0.) .and. (err_up .ne. 0.)) ul='UL'
             if ((err_up .ne. 0.) .and. (err_lo .eq. 0.)) ul='UL'
             if (flux_err .gt. flux) then
                flux=2.*flux_err
                flux_err=0.
             endif
             write(lu_out,'(f9.5,'' | '',f9.5,'' | '',es10.3,'' | '',es10.3,'' | '',
     &        es10.3,'' | '',es10.3,'' | '',f10.2,'' | '',f10.2,'' |'',1x,a,''|'')')
     &                      rra,rdec,freq,one,flux,flux_err,mjdstart,mjdend,ul
           ENDIF
         ENDIF
      ENDDO 
 99   CONTINUE
      END
