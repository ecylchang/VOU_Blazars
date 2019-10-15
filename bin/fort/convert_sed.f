      PROGRAM convert_sed
c
c converts the format of a SED file produced with the find_candidates tool
c to the format required by the ASDC SED tool  
c
c
      IMPLICIT none
      INTEGER*4 ier, lu_in, lu_out, lenact,sfound,rtype,im,in,length,is,ie
      REAL*4 mjdstart,mjdend, freq, one, err_up,err_lo
      REAL*4 flux,flux_err
      real*8 rra,rdec
      CHARACTER*2 ul
      character*14 stringin,catalog
      CHARACTER*80 input_file,output_file,output_file2,ref1,ref2,ref3,ref4,ref5
      Character*200 string,reff
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
      im=index(string(in+1:length),' ')+in
      output_file=string(in+1:im-1)
      in=im
c      im=index(string(in+1:length),' ')+in
      output_file2=string(in+1:length)

c      write(*,*) input_file
c      write(*,*) output_file2

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
      open(12,file=output_file2,status='unknown',iostat=ier)
c      write(12,'(f10.5,'','',f10.5)') rra,rdec
      write(12,'(a)') "freq. ,flux ,err_flux ,MJD_start ,MJD_end ,catalog,reference"
      DO WHILE(ok)
         ul = '  '
         READ(lu_in,'(a)',end=99) string 
c         print *,'string ',string(1:lenact(string))
         IF (string(2:2).NE.'=') THEN
           READ(string(1:lenact(string)),*) freq, flux, err_up, err_lo, mjdstart,mjdend,catalog
           read(string(85:lenact(string)),'(a)') reff
           ref1=' '
           ref2=' '
           ref3=' '
           ref4=' '
           if (catalog == 'DEBL') reff='3FHL EBL-corrected flux'
           is=index(reff(1:lenact(reff)),',')
           if (is .ne. 0) read(reff(1:is-1),'(a)') ref1
           ie=index(reff(is+1:lenact(reff)),',')+is
           if (is .ne. ie-1) read(reff(is+1:ie-1),'(a)') ref2
           is=ie
           ie=index(reff(is+1:lenact(reff)),',')+is
           if (is .ne. ie-1) read(reff(is+1:ie-1),'(a)') ref3
           is=ie
           ie=index(reff(is+1:lenact(reff)),',')+is
           if (is .ne. ie-1) read(reff(is+1:ie-1),'(a)') ref4
           is=ie
           if (is .ne. ie-1) read(reff(is+1:lenact(reff)),'(a)') ref5
c           write(*,*) ref1(1:lenact(ref1)),ref2(1:lenact(ref2)),ref3(1:lenact(ref3)),ref4(1:lenact(ref4))
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
             write(12,'(es10.3,'','',es10.3,'','',es10.3,'','',f10.2,'','',f10.2,'','',a,'','',a,a,a,a,a)')
     &        freq,flux,flux_err,mjdstart,mjdend,catalog,ref1(1:lenact(ref1)),
     &        ref2(1:lenact(ref2)),ref3(1:lenact(ref3)),ref4(1:lenact(ref4)),ref5(1:lenact(ref5))
           ENDIF
         ENDIF
      ENDDO 
 99   CONTINUE
      END
