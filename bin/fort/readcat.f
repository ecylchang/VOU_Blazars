      program readcat
c this program reads the output from vo tool and prepares an input file for find_candidates
      implicit none
      integer*4 icat,i,is,ie,ia,it,in,out,iskip
      integer*4 rah,ram,decd,decm,length,ns,lenact
      real*8 radeg,decdeg,ra_center,dec_center,errrad,errmaj,errmin,errang
      real*4 ras,decs,radius,posxerr,posyerr,poserr,posang,major,minor,nh,crtflux
      real*4 racsflux,racsflux_err,pa
      character*1 sign
      character*2 cratecheck
      character*200 inputlist,outputlist
      character*80  catalog,ra,dec,catname
      character*400 chandraflux
      character*1000 head,value
      character*800 flux,string
      logical ok

      call rdforn(string,length)
      call rmvlbk(string)
      in=index(string(1:length),' ')
      inputlist=string(1:in-1)

      out=index(string(in+1:length),' ')+in
      outputlist=string(in+1:out-1)
      iskip=index(outputlist(1:len(outputlist)),'_output')
      if (iskip .eq. 0) iskip=index(outputlist(1:len(outputlist)),'/output')
      !print *,'string(out+1:length) ',string(out+1:length)
      read(string(out+1:length),*) ra_center,dec_center,radius,nh,errrad,errmaj,errmin,errang
      !write(*,*) inputlist,outputlist

      open(11,file=inputlist,status='old')
      open(13,file=outputlist)
      write(13,'("RA= ",f9.5,2x,"Dec= ",f9.5,2x,"Searching radius= ",f6.2)')
     &      ra_center,dec_center,radius
      write(13,'("nH= ",es9.3,2x,"Error circle/elliptical= ",4(f6.2,2x))')nh,errrad,errmaj,errmin,errang
      ok=.true.
      do i=1,2000
         icat=0
         read(11,'(a)',end=100) string
         it=index(string(iskip+1:len(string)),'.')+iskip
c  read the catalog name from the input list of catalogs to be checked
         read(string(iskip+1:it-1),'(a)') catalog
         !write(*,*) catalog(1:lenact(catalog)),inputlist(iskip+1:iskip+12)
         if (inputlist(iskip+1:iskip+12) == "catlist2.txt") then
            is=it
            ie=index(string(is+1:len(string)),'.')+is
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
            if ((catalog(1:it-1) == '3fhl') .or. (catalog(1:it-1) == '3fgl') .or. 
     &          (catalog(1:it-1) == '4fgldr3') .or. (catalog(1:it-1) == 'mst12y') .or. 
     &          (catalog(1:it-1) == 'wiseme') .or. (catalog(1:it-1) == '2agile')  .or. 
     &          (catalog(1:it-1) == 'fmev') .or. (catalog(1:it-1) == '4fgldr4')) then
               read(value(1:is-1),'(a)') catname
            endif
c the catalogs without source name
            !print *,'catalog(1:it-1) ',catalog(1:it-1)
            if ((catalog(1:it-1) == 'sumss') .or. (catalog(1:it-1) == 'gb87')
     &          .or. (catalog(1:it-1) == '2bigb') .or. (catalog(1:it-1) == 'smarts')
     &          .or. (catalog(1:it-1) == '2mass') .or. (catalog(1:it-1) == 'xrtdeep')
     &          .or. (catalog(1:it-1) == 'gaia2')  .or. (catalog(1:it-1) == 'panstarrs')  
     &          .or. (catalog(1:it-1) == 'xrtspec') .or. (catalog(1:it-1) == 'magic')
     &          .or. (catalog(1:it-1) == 'veritas') .or. (catalog(1:it-1) == 'mquas')
     &          .or. (catalog(1:it-1) == 'neowise') .or. (catalog(1:it-1) == 'nublazar')
     &          .or. (catalog(1:it-1) == 'bepposax') .or. (catalog(1:it-1) == 'vlssr')
     &          .or. (catalog(1:it-1) == 'ipccs857') .or. (catalog(1:it-1) == 'unwise')) then
c     &          .or. (catalog(1:it-1) == 'unwise')) then
               is=0
               ie=index(value(1:len(value)),',')
            endif
            if (catalog(1:it-1) == 'crates') then !check crates counterpart
               read(value(is+1:ie-1),'(a)') cratecheck
               if (ns .eq. 0) then
                  if ((cratecheck == '1') .or. (cratecheck == '0')) then
                     is=ie
                     ie=index(value(is+1:len(value)),',')+is
                  else
                     goto 300
                  endif
               else
                  is=ie
                  ie=index(value(is+1:len(value)),',')+is
                  read(value(is+1:ie-1),'(a)') crtflux
               endif
               is=ie
               ie=index(value(is+1:len(value)),',')+is
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
            if (catalog(1:it-1) == 'iraspsc') then
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               is=ie
               ie=index(value(is+1:len(value)),',')+is
            endif
            read(value(is+1:ie-1),'(a)') ra
            ia=index(value(is+1:ie-1),' ')
            !print *,' catalog(1:it-1) ra is ie ',catalog(1:it-1),ra,is,ie
            if ((ie-1) > (is+1)) then 
                if (ia .ne. 0) then
                   read(ra(1:2), '(i2)',err=300)rah
                   read(ra(4:5), '(i2)')ram
                   read(ra(7:ie-1),*)ras
                   call chra(radeg,rah,ram,ras,0)
                else
                   read(ra(1:ie-is-1),*)radeg
                endif
            else
                radeg = 0.
            endif
            is=ie
            ie=index(value(is+1:len(value)),',')+is
            if (catalog(1:it-1) == '3hsp')  ie=is+30
            if (catalog(1:it-1) == '4lacdr3')  ie=is+30
            if (catalog(1:it-1) == '5bzcat') ie=is+30
            if (catalog(1:it-1) == 'bros') ie=is+30
            if (catalog(1:it-1) == 'zw') ie=is+30
            if (catalog(1:it-1) == 'psz2') ie=is+30
            if (catalog(1:it-1) == 'abell') ie=is+30
            if (catalog(1:it-1) == 'mcxc') ie=is+30
            if (catalog(1:it-1) == 'whl') ie=is+30
            if (catalog(1:it-1) == 'swxcs') ie=is+30
            if (catalog(1:it-1) == 'pulsar') ie=is+30
            if (catalog(1:it-1) == 'f2psr') ie=is+30
c            if ((catalog(1:it-1) == 'crates') .and. (ns .eq. 0)) ie=is+30
            if (catalog(1:it-1) == 'north20') ie=is+30
            if (catalog(1:it-1) == 'mst12y') ie=is+30
            if (catalog(1:it-1) == 'fgrb') ie=is+30
            if (catalog(1:it-1) == 'f357cat') ie=is+30
!the source pos error in front of dec
            if ((catalog(1:it-1) == 'sumss') .or. (catalog(1:it-1) == 'gb6') .or.
     &      (catalog(1:it-1) == 'gb87') .or. (catalog(1:it-1) == 'uvot') .or.
     &      (catalog(1:it-1) == 'gleam') .or. (catalog(1:it-1) == 'lotss') .or.
     &      (catalog(1:it-1) == 'gaia2') .or. (catalog(1:it-1) == 'tgss150')) then
               read(value(is+1:ie-1),*) posxerr
               if ((catalog(1:it-1) == 'gb87') .or. (catalog(1:it-1) == 'gb6')) posxerr=posxerr*15.
               is=ie
               ie=index(value(is+1:len(value)),',')+is
            endif
            if (catalog(1:it-1) == 'iraspsc') ie=index(value(is+1:len(value)),' ')+is
            read(value(is+1:ie-1),'(a)') dec
            if (ia .ne. 0) then
               sign=dec(1:1)
               read(dec(2:3),'(i2)')decd
               read(dec(5:6),'(i2)')decm
               read(dec(8:ie-1),*)decs
               call chdec(decdeg,decd,decm,decs,0)
               if (sign == '-') decdeg=-abs(decdeg)
            else
               !print *,'Catalog ie is ',catalog(1:lenact(catalog)),ie,is
               !print *,'dec ',dec(1:lenact(dec))
               if (ie > is) then
                  read(dec(1:ie-is-1),*) decdeg
               else
                  goto 185
               endif
            endif
            if (catalog(1:it-1) == 'chandracsc2') then
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               if (value(is+1:is+1) .NE.',') read(value(is+1:ie-1),*) posang
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
               chandraflux=value(is+1:lenact(value))
c               print *,'posang major minor flux ',posang,major,minor,flux
            endif
            if ((catalog(1:it-1) == 'veritas') .or. (catalog(1:it-1) == 'magic')) then
               radeg=ra_center
               decdeg=dec_center
            endif
c read other pos_err
            if ((catalog(1:it-1) == 'nvss') .or.(catalog(1:it-1) == 'hst') .or.
     &          (catalog(1:it-1) == 'wise') .or. (catalog(1:it-1) == 'spire') .or.
     &          (catalog(1:it-1) == 'spire') .or.
     &          (catalog(1:it-1) == 'panstarrs') .or. (catalog(1:it-1) == 'f357det') ) then
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               if (is .ne. ie-1) read(value(is+1:ie-1),*) posxerr
               if (catalog(1:it-1) == 'nvss') posxerr=posxerr*15.
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               if (is .ne. ie-1) read(value(is+1:ie-1),*) posyerr
               poserr=sqrt((posxerr*posxerr)+(posyerr*posyerr))
               if (poserr > 100.) poserr = 100.
               if (catalog(1:it-1) == 'f357det') poserr=poserr*3600
               !write(*,*) catalog(1:it-1),poserr,posxerr,posyerr
            endif
            if (catalog(1:it-1) == 'racs') then
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               read(value(is+1:ie-1),*) racsflux
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               read(value(is+1:ie-1),*) racsflux_err
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               if (is .ne. ie-1) read(value(is+1:ie-1),*) posxerr
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               if (is .ne. ie-1) read(value(is+1:ie-1),*) posyerr
               is=ie
               if (is .ne. ie-1) read(value(is+1:lenact(value)),*) pa
            endif 
            if ((catalog(1:it-1) == 'sumss') .or. (catalog(1:it-1) == 'gb6') .or.
     &      (catalog(1:it-1) == 'gb87') .or. (catalog(1:it-1) == 'uvot') .or.
     &      (catalog(1:it-1) == 'gleam') .or. (catalog(1:it-1) == 'lotss') .or.
     &      (catalog(1:it-1) == 'gaia2') .or. (catalog(1:it-1) == 'tgss150')) then
               is=ie
               ie=index(value(is+1:len(value)),',')+is
               read(value(is+1:ie-1),*) posyerr
               poserr=sqrt((posxerr*posxerr)+(posyerr*posyerr))
               !write(*,*) catalog(1:it-1),poserr,posxerr,posyerr
            endif
            if ((catalog(1:it-1) == '2mass') .or. (catalog(1:it-1) == 'iraspsc')
     &             .or. (catalog(1:it-1) == 'akaribsc') .or. (catalog(1:it-1) == 'vlassql'))  then
               if (catalog(1:it-1) == 'iraspsc') then
                  is=index(value(1:len(value)),',')
               else
                  is=ie
               endif
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
c               if (catalog(1:it-1) == 'fermi8yr') poserr=poserr*60.
            endif
            if (catalog(1:it-1) == 'cma') then
               is=ie
               ie=index(value(is+1:len(value)),' ')+is
               read(value(is+1:ie-1),*) poserr
            endif
            if (catalog(1:it-1) == 'gaia2')poserr=poserr*0.001
c read the flux
            is=ie
            ie=index(value(is+1:len(value)),' ')+is
            if (catalog(1:it-1) == 'mquas') ie=is+40
            read(value(is+1:ie-1),'(a)') flux
            if ((catalog(1:it-1) == '3hsp') .or. (catalog(1:it-1) == '5bzcat') .or.
     &           (catalog(1:it-1) == 'zw') .or. (catalog(1:it-1) == 'psz2') .or.
     &           (catalog(1:it-1) == 'abell') .or. (catalog(1:it-1) == 'mcxc') .or.
     &           (catalog(1:it-1) == 'whl') .or. (catalog(1:it-1) == 'swxcs') .or.
     &           ((catalog(1:it-1) == 'crates') .and. (ns .eq. 0)) .or.
     &           (catalog(1:it-1) == 'pulsar') .or. (catalog(1:it-1) == 'f2psr') .or.
     &           (catalog(1:it-1) == 'mst12y') .or. (catalog(1:it-1) == 'fgrb') .or. 
     &           (catalog(1:it-1) == 'bros') .or. (catalog(1:it-1) == 'f357cat') .or. 
     &           (catalog(1:it-1) == '4lacdr3') ) then
               ie=index(value(1:len(value)),',')
               read(value(1:ie-1),'(a)') flux
            endif
            if  ((catalog(1:it-1) == 'veritas') .or. (catalog(1:it-1) == 'magic')) then
               read(value(1:ie-1),'(a)') flux
            endif
            if ((catalog(1:5) == '2mass') .or. (catalog(1:it-1) == 'gaia2')) then
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
            if (catalog(1:11) == 'chandracsc2') then
               write(13,'(i4,",",a,",",2(f9.5,","),f7.3,",",a)')
     &             ns,catalog(1:it-1),radeg,decdeg,poserr,chandraflux(1:lenact(chandraflux))
            else if ((catalog(1:it-1) == 'nvss') .or. (catalog(1:it-1) == 'sumss') .or.
     &         (catalog(1:it-1) == '2mass') .or. (catalog(1:it-1) == 'gb6') .or.
     &         (catalog(1:it-1) == 'gb87') .or. (catalog(1:it-1) == 'hst') .or.
     &         (catalog(1:it-1) == 'uvot') .or. (catalog(1:it-1) == 'wise') .or.
     &         (catalog(1:5) == 'spire') .or. (catalog(1:it-1) == 'cma') .or.
     &         (catalog(1:it-1) == 'panstarrs') .or. (catalog(1:it-1) == 'gaia2') .or.
     &         (catalog(1:it-1) == 'tgss150') .or. (catalog(1:it-1) == 'gleam') .or.
     &         (catalog(1:it-1) == 'lotss') .or. (catalog(1:it-1) == 'akaribsc') .or.
     &         (catalog(1:it-1) == 'iraspsc') .or. (catalog(1:it-1) == 'vlassql') .or.
     &         (catalog(1:it-1) == 'f357det') ) then
               write(13,'(i4,",",a,",",2(f9.5,","),f9.3,",",a)')
     &             ns,catalog(1:it-1),radeg,decdeg,poserr,flux(1:ie-1)
            else if ((catalog(1:it-1) == '3fhl') .or. (catalog(1:it-1) == '3fgl') .or.
     &               (catalog(1:it-1) == '2agile') .or. (catalog(1:it-1) == 'fmev') .or.
     &               (catalog(1:it-1) == '4fgldr3') .or. (catalog(1:it-1) == 'wiseme') .or. 
     &               (catalog(1:it-1) == '4fgldr4') ) then
               write(13,'(i4,",",a,",",2(f9.5,","),a,",",a)')
     &         ns,catalog(1:it-1),radeg,decdeg,catname(1:lenact(catname)),flux(1:ie-1)
            else if (catalog(1:it-1) == 'racs') then
              write(13,'(i4,",",a,",",2(f9.5,","),f9.1,",",f7.1,",",f6.2,",",f6.2,",",f5.0)')
     &        ns,catalog(1:it-1),radeg,decdeg,racsflux,racsflux_err,posxerr,posyerr,pa
            else if ((catalog(1:it-1) == 'crates') .and. (ns .ne. 0)) then
              write(13,'(i4,",",a,",",2(f9.5,","),a,",",a)')
     &        ns,catalog(1:it-1),radeg,decdeg,crtflux,flux(1:ie-1)
            else
               write(13,'(i4,",",a,",",2(f9.5,","),a)') ns,catalog(1:it-1),radeg,decdeg,flux(1:ie-1)
            endif
185      continue
         enddo
200      continue
         if (inputlist(iskip+1:iskip+12) == 'catlist2.txt') then
            if (icat > 1) then 
               write(*,'("Candidate nr.",i4,",",2x,a,i5,2x,"points")') ns,catalog(1:12),icat
            else
               write(*,'("Candidate nr.",i4,",",2x,a,i5,2x,"point")') ns,catalog(1:12),icat
            endif
         else
            if (icat > 1) then 
               write(*,'(a,i5,2x,"points")') catalog(1:12),icat
            else
               write(*,'(a,i5,2x,"point")') catalog(1:12),icat
            endif
         endif
         close(12)
      enddo
100   continue
      close(11)
      close(13)
      end
