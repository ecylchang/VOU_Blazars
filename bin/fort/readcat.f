      program readcat
c this program read the output from vo tool and make a input for find candidates
      implicit none
      integer*4 icat,i,is,ie,ia,it,in,ip,out,iskip
      integer*4 rah,ram,decd,decm,length,ns,lenact
      real*8 radeg,decdeg,ra_center,dec_center,errrad,errmaj,errmin,errang
      real*4 ras,decs,radius,posxerr,posyerr,poserr,posang,major,minor,nh,ra1,ra2,dec1,dec2
      character*1 sign
      character*2 cratecheck
      character*80 catalog,ra,dec,inputlist,test,outputlist,catname
      character*1000 head,value
      character*800 flux,string
      logical ok,there

      call rdforn(string,length)
      call rmvlbk(string)
      in=index(string(1:length),' ')
      inputlist=string(1:in-1)

      out=index(string(in+1:length),' ')+in
      outputlist=string(in+1:out-1)
      iskip=index(outputlist(1:len(outputlist)),'_output')
c      write(*,*) iskip
      read(string(out+1:length),*) ra_center,dec_center,radius,nh,errrad,errmaj,errmin,errang
      !write(*,*) inputlist,outputlist
c      write(*,*) ra_center,dec_center,radius,nh

      open(11,file=inputlist,status='old')
      open(13,file=outputlist)
      write(13,'("RA= ",f9.5,2x,"Dec= ",f9.5,2x,"Searching radius= ",f6.2)')
     &      ra_center,dec_center,radius
      write(13,'("nH= ",es9.3,2x,"Error circle/elliptical= ",4(f6.2,2x))')nh,errrad,errmaj,errmin,errang
      ok=.true.
      do i=1,2000
         icat=0
         read(11,'(a)',end=100) string
c         is=index(string(1:len(string)),'_')
         it=index(string(iskip+1:len(string)),'.')+iskip
         read(string(iskip+1:it-1),'(a)') catalog
c         write(*,*) catalog
         if (inputlist(iskip+1:iskip+12) == "catlist2.txt") then
            is=it
            ie=index(string(is+1:len(string)),'.')+is
c            write(*,*) is,ie
            read(string(is+1:ie-1),*) ns
         else
            ns=0
         endif
         open(12,file=string,status='old')
         read(12,'(a)') head
         do while (ok)
300         continue
            read(12,'(a)',end=200) value
            icat=icat+1
            is=index(value(1:len(value)),',')
            ie=index(value(is+1:len(value)),',')+is
c gamma-ray catalog print name
            if ((catalog(1:it-1) == '3fhl') .or. (catalog(1:it-1) == '3fgl')
     &          .or. (catalog(1:it-1) == 'fermi8yr') .or. (catalog(1:it-1) == 'mst9y')
     &          .or. (catalog(1:it-1) == 'agile') .or. (catalog(1:it-1) == 'fmev')) then
               read(value(1:is-1),'(a)') catname
            endif
c the catalog without source name
            if ((catalog(1:it-1) == 'sumss') .or. (catalog(1:it-1) == 'gb87')
     &          .or. (catalog(1:it-1) == '2mass') .or. (catalog(1:it-1) == 'xrtdeep')
     &          .or. (catalog(1:it-1) == 'panstarrs') .or. (catalog(1:it-1) == 'gaia')
     &          .or. (catalog(1:it-1) == 'xrtspec') .or. (catalog(1:it-1) == 'magic')
     &          .or. (catalog(1:it-1) == 'veritas') ) then
               is=0
               ie=index(value(1:len(value)),',')
            endif
            if ((catalog(1:it-1) == 'crates') .and. (ns .eq. 0)) then !check crates counterpart
               read(value(is+1:ie-1),'(a)') cratecheck
               if ((cratecheck == '1') .or. (cratecheck == '0')) then
                  is=ie
                  ie=index(value(is+1:len(value)),',')+is
               else
                  goto 300
               endif
            endif
            if ((catalog(1:it-1) == 'cma') .or. (catalog(1:it-1) == 'north20')) then
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               if (catalog(1:it-1) == 'north20') then
                  is=ie
                  ie=index(value(is+1:len(value)),',')+is
                  is=ie
                  ie=index(value(is+1:len(value)),',')+is
                  is=ie
                  ie=index(value(is+1:len(value)),',')+is
                  is=ie
                  ie=index(value(is+1:len(value)),',')+is

               endif
            endif
            read(value(is+1:ie-1),'(a)') ra
            ia=index(value(is+1:ie-1),' ')
            if (ia .ne. 0) then
               read(ra(1:2),'(i2)')rah
               read(ra(4:5),'(i2)')ram
               read(ra(7:ie-1),*)ras
               call chra(radeg,rah,ram,ras,0)
            else
               read(ra(1:ie-1),*)radeg
            endif
            is=ie
            ie=index(value(is+1:len(value)),',')+is
            if (catalog(1:it-1) == '3hsp')  ie=is+30
            if (catalog(1:it-1) == '5bzcat') ie=is+30
            if (catalog(1:it-1) == 'zw') ie=is+30
            if (catalog(1:it-1) == 'psz2') ie=is+30
            if (catalog(1:it-1) == 'abell') ie=is+30
            if (catalog(1:it-1) == 'mcxc') ie=is+30
            if (catalog(1:it-1) == 'whl') ie=is+30
            if (catalog(1:it-1) == 'swxcs') ie=is+30
            if (catalog(1:it-1) == 'pulsar') ie=is+30
            if (catalog(1:it-1) == 'f2psr') ie=is+30
            if ((catalog(1:it-1) == 'crates') .and. (ns .eq. 0)) ie=is+30
            if (catalog(1:it-1) == 'north20') ie=is+30
            if (catalog(1:it-1) == 'mst9y') ie=is+30
!the source pos error in front of dec
            if ((catalog(1:it-1) == 'sumss') .or. (catalog(1:it-1) == 'gb6') .or.
     &           (catalog(1:it-1) == 'gb87') .or. (catalog(1:it-1) == 'uvot') .or.
     &           (catalog(1:it-1) == 'gaia') ) then
               read(value(is+1:ie-1),*) posxerr
               if ((catalog(1:it-1) == 'gb87') .or. (catalog(1:it-1) == 'gb6')) posxerr=posxerr*15.
               is=ie
               ie=index(value(is+1:len(value)),',')+is
            endif
            if (catalog(1:it-1) == 'chandra') then
               read(value(is+1:ie-1),*) posang
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               read(value(is+1:ie-1),*) major
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               read(value(is+1:ie-1),*) minor
               posxerr=sqrt(((sin(posang)*major)**2)+((cos(posang)*minor)**2))
               posyerr=sqrt(((cos(posang)*major)**2)+((sin(posang)*minor)**2))
               poserr=max(posxerr,posyerr)
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               !write(*,*) posxerr,posyerr,poserr
            endif
            read(value(is+1:ie-1),'(a)') dec
            if (ia .ne. 0) then
               sign=dec(1:1)
               read(dec(2:3),'(i2)')decd
               read(dec(5:6),'(i2)')decm
               read(dec(8:ie-1),*)decs
               call chdec(decdeg,decd,decm,decs,0)
               if (sign == '-') decdeg=-abs(decdeg)
            else
               read(dec(1:ie-1),*) decdeg
            endif
            if ((catalog(1:it-1) == 'veritas') .or. (catalog(1:it-1) == 'magic')) then
               radeg=ra_center
               decdeg=dec_center
            endif
c read other pos_err
            if ((catalog(1:it-1) == 'nvss') .or.(catalog(1:it-1) == 'hst') .or.
     &            (catalog(1:it-1) == 'wise') .or. (catalog(1:5) == 'spire') .or.
     &            (catalog(1:it-1) == 'panstarrs'))  then
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               if ((catalog(1:it-1) == 'panstarrs')) then
                  is=index(value(is+1:len(value)),',')+is
                  ie=index(value(is+1:len(value)),',')+is
                  posxerr=0.
                  posyerr=0.
               endif
               if (is .ne. ie-1) read(value(is+1:ie-1),*) posxerr
               if (catalog(1:it-1) == 'nvss') posxerr=posxerr*15.
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               if (is .ne. ie-1) read(value(is+1:ie-1),*) posyerr
               poserr=sqrt((posxerr*posxerr)+(posyerr*posyerr))
               !write(*,*) catalog(1:it-1),poserr,posxerr,posyerr
            endif
            if ((catalog(1:it-1) == 'sumss') .or. (catalog(1:it-1) == 'gb6') .or.
     &          (catalog(1:it-1) == 'gb87') .or. (catalog(1:it-1) == 'uvot') .or.
     &          (catalog(1:it-1) == 'gaia')) then
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               read(value(is+1:ie-1),*) posyerr
               poserr=sqrt((posxerr*posxerr)+(posyerr*posyerr))
               !write(*,*) catalog(1:it-1),poserr,posxerr,posyerr
            endif
            if ((catalog(1:it-1) == '2mass') .or. (catalog(1:it-1) == 'fermi8yr'))  then
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               if (is .ne. ie-1) read(value(is+1:ie-1),*) major
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               if (is .ne. ie-1) read(value(is+1:ie-1),*) minor
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               if (is .ne. ie-1) read(value(is+1:ie-1),*) posang
               posxerr=sqrt(((sin(posang)*major)**2)+((cos(posang)*minor)**2))
               posyerr=sqrt(((cos(posang)*major)**2)+((sin(posang)*minor)**2))
               poserr=max(posxerr,posyerr)
               if (catalog(1:it-1) == 'fermi8yr') poserr=poserr*60.
            endif
            if (catalog(1:it-1) == 'cma') then
               is=ie
               ie=index(value(is+1:len(value)),' ')+is
               read(value(is+1:ie-1),*) poserr
            endif
            if (catalog(1:it-1) == 'gaia') poserr=poserr*0.001
c read the flux
            is=ie
            ie=index(value(is+1:len(value)),' ')+is
            read(value(is+1:ie-1),'(a)') flux
            if ((catalog(1:it-1) == '3hsp') .or. (catalog(1:it-1) == '5bzcat') .or.
     &           (catalog(1:it-1) == 'zw') .or. (catalog(1:it-1) == 'psz2') .or.
     &           (catalog(1:it-1) == 'abell') .or. (catalog(1:it-1) == 'mcxc') .or.
     &           (catalog(1:it-1) == 'whl') .or. (catalog(1:it-1) == 'swxcs') .or.
     &      ((catalog(1:it-1) == 'crates') .and. (ns .eq. 0)) .or.
     &       (catalog(1:it-1) == 'pulsar') .or. (catalog(1:it-1) == 'f2psr') .or.
     &       (catalog(1:it-1) == 'mst9y')) then
               ie=index(value(1:len(value)),',')
               read(value(1:ie-1),'(a)') flux
            endif
            if  ((catalog(1:it-1) == 'veritas') .or. (catalog(1:it-1) == 'magic')) then
               read(value(1:ie-1),'(a)') flux
            endif
            if ((catalog(1:5) == '2mass') .or. (catalog(1:it-1) == 'gaia')) then
               is=index(value(is+1:len(value)),',')+is
               ie=index(value(is+1:len(value)),' ')+is
               read(value(is+1:ie-1),'(a)') flux
            endif
            if ((catalog(1:it-1) == 'cma') .or. (catalog(1:it-1) == 'north20')) then
               is=index(value(1:len(value)),',')
               ie=index(value(is+1:len(value)),',')+is
               ie=index(value(ie+1:len(value)),',')+ie
               read(value(is+1:ie-1),'(a)') flux
               if (catalog(1:it-1) == 'north20') then
                  ie=index(value(ie+1:len(value)),',')+ie
                  ie=index(value(ie+1:len(value)),',')+ie
                  ie=index(value(ie+1:len(value)),',')+ie
                  ie=index(value(ie+1:len(value)),',')+ie
                  read(value(is+1:ie-1),'(a)') flux
               endif
            endif
c write the data
            if (catalog(1:7) == 'chandra') then
               write(13,'(i4,",",a,",",2(f9.5,","),f7.3,",",a)')
     &             ns,catalog(1:it-1),decdeg,radeg,poserr,flux(1:ie-1)
            else if ((catalog(1:it-1) == 'nvss') .or. (catalog(1:it-1) == 'sumss') .or.
     &               (catalog(1:it-1) == '2mass') .or. (catalog(1:it-1) == 'gb6') .or.
     &               (catalog(1:it-1) == 'gb87') .or. (catalog(1:it-1) == 'hst') .or.
     &               (catalog(1:it-1) == 'uvot') .or. (catalog(1:it-1) == 'wise') .or.
     &               (catalog(1:5) == 'spire') .or. (catalog(1:it-1) == 'cma') .or.
     &               (catalog(1:it-1) == 'panstarrs') .or. (catalog(1:it-1) == 'gaia') ) then
               write(13,'(i4,",",a,",",2(f9.5,","),f7.3,",",a)')
     &             ns,catalog(1:it-1),radeg,decdeg,poserr,flux(1:ie-1)
            else if ((catalog(1:it-1) == '3fhl') .or. (catalog(1:it-1) == '3fgl') .or.
     &               (catalog(1:it-1) == 'agile') .or. (catalog(1:it-1) == 'fmev')) then
               write(13,'(i4,",",a,",",2(f9.5,","),a,",",a)')
     &         ns,catalog(1:it-1),radeg,decdeg,catname(1:lenact(catname)),flux(1:ie-1)
            else if (catalog(1:it-1) == 'fermi8yr') then
               write(13,'(i4,",",a,",",2(f9.5,","),a,",",f7.3,",",a)')
     &             ns,catalog(1:it-1),radeg,decdeg,catname(1:lenact(catname)),poserr,flux(1:ie-1)
            else
               write(13,'(i4,",",a,",",2(f9.5,","),a)') ns,catalog(1:it-1),radeg,decdeg,flux(1:ie-1)
            endif
         enddo
200      continue
         if (inputlist(iskip+1:iskip+12) == 'catlist2.txt') then
            write(*,'("Candidate nr.",i4,",",2x,a,i4,2x,"point(s)")') ns,catalog(1:it-1),icat
         else
            write(*,'(a,i4,2x,"point(s)")') catalog(1:it-1),icat
         endif
         close(12)
      enddo
100   continue
      close(11)
      close(13)
      end
