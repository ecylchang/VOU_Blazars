      PROGRAM gnomo_plot_types
c
      IMPLICIT NONE
      INTEGER*4 symbol, i,n,m,in,is,j,iskip,isource,om,iconfig,iarr
      INTEGER*4 no_of_isoalpha, no_of_isodelta, n_points,lenact,arrsize(30)
      INTEGER*4 lu_infile, length, im, ip,lu_out,rah,irm,id,idm
      INTEGER*4 icol1,icol2,icol3,icol4,icol5,icol11,icol14,icol12,icol13
      REAL*4 isoalpha, isodelta, step_delta,cs,xtick,ytick
      REAL*4 xpoly(4),ypoly(4),R,G,B,rasec,decsec,step,cc
      REAL*4 ra1, ra2, dec1, dec2,x_o, y_o, x_err2, y_err2,xx,yy
      REAL*4 raarr(3),decarr(3),ra_center,dec_center,dec_start,ra_o,dec_o
      REAL*4 step_alpha,pos,fmin,fmax,err_radius, radius, radian
      REAL*4 ellipser1, ellipser2, ellipserot,ellipser1_2, ellipser2_2
      REAL*4 ellipserot_2, err_radius_2,ra_err2,dec_err2
      REAL*4 cxray, cradio, cx
      REAL*8 rra, ddec
      CHARACTER*1 sign
      CHARACTER*85 string, newstring,title
      CHARACTER*200 filein,fileout,device,webprograms,array_size
      CHARACTER*800 stringin
      CHARACTER*14 xaxis_label,yaxis_label
      REAL*4 radius_scale
      REAL*4 pro_scale
      LOGICAL ok,there
      integer*4 status

      integer*4,dimension(:),allocatable :: code
      real*4,dimension(:),allocatable :: ra,dec,epos
c      real*4,dimension(:),allocatable :: x1,y1,racat,deccat
      real*4,dimension(:),allocatable :: x,y!,names_ra,names_dec
      real*4,dimension(:),allocatable :: run_alpha,run_dec,x_grid,y_grid,afmin,afmax

      integer*4,dimension(:),allocatable :: s11
      character*4,dimension(:),allocatable :: tcol1,tcol2,tcol3,tcol4,tcol5,tcol11
      real*4,dimension(:),allocatable :: ra_col1,dec_col1,csx1,csr1
      real*4,dimension(:),allocatable :: ra_col2,dec_col2,csx2,csr2
      real*4,dimension(:),allocatable :: ra_col3,dec_col3,csx3,csr3
      real*4,dimension(:),allocatable :: ra_col4,dec_col4,csx4,csr4
      real*4,dimension(:),allocatable :: ra_col5,dec_col5,csx5,csr5
      real*4,dimension(:),allocatable :: ra_col11,dec_col11,csx11,csr11

      real*4,dimension(:),allocatable :: ra_col13,dec_col13,epos_col13,color13,csr13
      integer*4,dimension(:),allocatable :: s12
      real*4,dimension(:),allocatable :: ra_col12,dec_col12,csx12,csr12
      character*4,dimension(:),allocatable :: tcol12
      integer*4,dimension(:),allocatable :: s14
      real*4,dimension(:),allocatable :: ra_col14,dec_col14
      character*4,dimension(:),allocatable :: tcol14

* External references :
      !data symb_size /9.,7.,6.,5.,3.5,2.5,1.6,1.0,0.8,0.4/
      COMMON /ecetype/pro_scale,radius_scale

      i = 0
      ok = .TRUE.
      radius_scale = 60.
      pro_scale = 1.0
      radian = 45.0/atan(1.0)
c
c Get input parameters
c
      CALL rdforn(stringin,length)
      IF ( length.NE.0 ) then
         in = index(stringin(1:length),',')
         filein=stringin(1:in-1)
         iskip=index(filein(1:len(filein)),'_error')
         if (iskip .eq. 0) iskip=index(filein(1:len(filein)),'_find_out')
         if (iskip .eq. 0) iskip=index(filein(1:len(filein)),'_RX')
         if (iskip .eq. 0) iskip=index(filein(1:len(filein)),'/error')
         if (iskip .eq. 0) iskip=index(filein(1:len(filein)),'/find_out')
         if (iskip .eq. 0) iskip=index(filein(1:len(filein)),'/RX')
         im=index(stringin(in+1:length),',')+in
         fileout=stringin(in+1:im-1)
         in=im
         im = index(stringin(in+1:length),',')+in
         device = stringin(in+1:im-1)
         in=im
         im = index(stringin(in+1:length),',')+in
         webprograms=stringin(in+1:im-1)
         ip = index(stringin(im+1:length),',')+im

         INQUIRE (FILE=filein,EXIST=there)
            IF (.NOT.there) THEN
               write (*,'('' file '',a,'' not found '')')
     &        filein(1:lenact(filein))
              STOP
           ENDIF

         read(stringin(ip+1:length),*) ra_center,dec_center,err_radius
     & ,radius,ra_o,dec_o,ellipser1,ellipser2,ellipserot,err_radius_2
     & ,ellipser1_2,ellipser2_2,ellipserot_2,ra_err2,dec_err2
      ELSE
         write(*,
     &'("Enter input_file,device,ra_center,dec_center,error_radius,search_radius",/,
     &             "e.g. gnomo_plot aa.log,test.ps/cps,4.91583,26.04778,40.,60.0")')
         STOP
      ENDIF

      iarr=0
      iconfig=0
      arrsize(1:30)=0
      i = 0
      radius= radius/60.
      ellipserot=-ellipserot
      ellipserot_2=-ellipserot_2
      lu_out = 11
      lu_infile = 12
      array_size=webprograms(1:lenact(webprograms))//'/array_size.cf'
      open(lu_infile,file=filein,status='old')
      open(lu_out,file=fileout,status='unknown')
      open(18,file=array_size,status='old')

      do while(ok)
      read(18,'(a)',end=700) stringin
      !write(*,*) stringin
      if (stringin(1:3) == '---' ) then
         iconfig=iconfig+1
      else
         if (iconfig .eq. 1) then
            iarr=iarr+1
            if (iarr .ge. 3) then
               is=index(stringin(1:len(stringin)),':')
               !write(*,*) 'CHECK',iconfig,iarr,is
               read(stringin(is+1:lenact(stringin)),*) arrsize(iarr-2)
            endif
         else if (iconfig .eq. 5) then
            iarr=iarr+1
            is=index(stringin(1:len(stringin)),':')
            read(stringin(is+1:lenact(stringin)),*) arrsize(iarr-2)
         endif
      endif
      enddo
700   continue
      close(18)

      allocate(code(arrsize(3)),epos(arrsize(3)),ra(arrsize(3)),dec(arrsize(3)))
      allocate(x(arrsize(11)),y(arrsize(11)))
      allocate(run_alpha(arrsize(12)),run_dec(arrsize(12)),x_grid(arrsize(12)),y_grid(arrsize(12)),afmin(arrsize(12)),afmax(arrsize(12)))

      if (filein(iskip+1:iskip+13) == 'error_map.txt') then
      radius=radius/60.
      DO WHILE (ok)
         i = i + 1
         READ (lu_infile,'(a)',END=100) stringin
          READ (stringin(1:lenact(stringin)),*,END=100,ERR=100)
     &            ra(i),dec(i),code(i),epos(i)
      ENDDO
      else
      DO WHILE (ok)
         i = i + 1
         READ (lu_infile,'(a)',END=500) stringin
         READ (stringin(1:lenact(stringin)),*,END=700,ERR=700)
     &            ra(i),dec(i),code(i)
      ENDDO
      endif
 100  CONTINUE
 500  continue
      n_points=i-1
      IF ( n_points.GT.35000 ) THEN
         print *, ' //max no of points (35000) exceeded '
         print *, ' //plotting first 35000 only '
      END IF
      CLOSE (lu_infile)
      symbol = 17
      no_of_isoalpha = 6
      no_of_isodelta = 6
      ra1 = 0.
      ra2 = 360.
      radian = 45.0/atan(1.0)
      dec1 = max(dec_center-radius,-90.)
      dec2 = min(dec_center+radius,+90.)
      dec_start = dec_center - radius

c   draw the grid
      IF ( abs(dec_center)+radius.EQ.90. ) THEN
        ra1 = ra_center - 90.
        ra1 = ra_center + 90.
      END IF
      IF ( abs(dec_center)+radius.LT.90. ) THEN
        ra1 = ra_center - asin(sin(radius/radian)/cos(dec_center/radian)
     &        )*radian
        ra2 = ra_center + asin(sin(radius/radian)/cos(dec_center/radian)
     &        )*radian
      END IF
      rra=ra_center
      ddec=dec_center
      status=0
      cc = cos(dec_center/radian)
      raarr(1)=ra_center
      decarr(1)=dec_center
      raarr(2)=ra_o
      decarr(2)=dec_o
      raarr(3)=ra_err2
      decarr(3)=dec_err2
      CALL gnom_projection(1,ra_center,dec_center,raarr(2),decarr(2),x,y)
      x_o=x(1)
      y_o=y(1)
      CALL pgbegin(0,device,1,1)
      CALL pgscf(1)
      CALL pgvport(0.2,0.8,0.25,0.7)
      CALL pgwindow(-radius*radius_scale/cc,+radius*radius_scale/cc,-radius*radius_scale,+radius*radius_scale)
      R = 0./255.
      G = 72./255.
      B = 189./255.
      CALL PGSCR(4,R,G,B)
      R=235./255.
      G=235./255.
      B=235./255.
      CALL PGSCR(10,R,G,B)
      call pgsci(10)
      xpoly(1)=-radius*radius_scale/cc
      ypoly(1)=-radius*radius_scale
      xpoly(2)=radius*radius_scale/cc
      ypoly(2)=-radius*radius_scale
      xpoly(3)=radius*radius_scale/cc
      ypoly(3)=radius*radius_scale
      xpoly(4)=-radius*radius_scale/cc
      ypoly(4)=radius*radius_scale
      call pgpoly(4,xpoly,ypoly)
      call pgsci(1)
      CALL pgslw(3)
      CALL pgwindow(ra_center+radius/cc,ra_center-radius/cc,dec_center-radius,dec_center+radius)
      step = nint(2.*radius/cc/6.*10.)/10.
      CALL pgtbox('BCNSTDY',step,0,'BCNSTDY',step,0)
      CALL pgwindow(-radius*radius_scale/cc,+radius*radius_scale/cc,-radius*radius_scale,+radius*radius_scale)
      CALL pgbox('CTMI',0.0,0,'CTMI',0.0,0)
      CALL pgmtxt('T',0.2,.08,0.5,'arcmin')
      CALL pgmtxt('R',0.6,.08,0.5,'arcmin')
      call pgsfs(1)
      xtick = 2.*radius/real(no_of_isoalpha)
      ytick = xtick
      call pgsci(1)
      CALL pgsls(4)
      CALL pgscf(1)
      call pgsch(.6)
      icol1=0
      icol2=0
      icol3=0
      icol4=0
      icol5=0
      icol11=0
      icol12=0
      icol13=0
      icol14=0
      cs = 1.
      isource=0
      call pgsci(4)
      CALL pgsls(3)
      CALL pgslw(3)
      IF(ra_center.NE.ra_o.OR.dec_center.NE.dec_o) THEN
           CALL gnom_circle_off (100,dec_center,err_radius,x,y,x_o,y_o)
           CALL pgline(100,x,y)
      ELSE
           CALL gnom_circle (100,dec_center,err_radius,x,y)
           CALL pgline(100,x,y)
      ENDIF
      CALL PGSCI(10)
      R=212./255.
      G=223./255.
      B=234./255.
      CALL PGSCR(10,R,G,B)
      CALL PGPOLY(100,x,y)


      call pgsci(4)
      CALL pgsls(3)
      CALL pgslw(3)
      IF(ra_center.NE.ra_o.OR.dec_center.NE.dec_o) THEN
         CALL gnom_ellipse (100,dec_center,ellipser1,ellipser2,ellipserot,x_o,y_o,x,y)
         CALL pgline(100,x,y)
      ELSE
         CALL gnom_ellipse (100,dec_center,ellipser1,ellipser2,ellipserot,0.0,0.0,x,y)
         CALL pgline(100,x,y)
      ENDIF
      CALL PGSFS(1)
      CALL PGSCI(10)
      CALL PGPOLY(100,x,y)

c     stile linea
      CALL pgsls(5)
c     larghezza linea
      CALL pgslw(3)
c     colore linea
      call pgsci(14)
      call pgscr(14,0.1,0.6,0.3)

      IF(ra_center.NE.ra_err2.OR.dec_center.NE.dec_err2) THEN
      CALL gnom_circle_off (100,dec_center,err_radius_2,x,y,x_err2,y_err2)
      CALL pgline(100,x,y)
      ELSE
      CALL gnom_circle (100,dec_center,err_radius_2,x,y)
      CALL pgline(100,x,y)
      ENDIF
      R=230./255.
      G=248./255.
      B=255./255.
      CALL PGSCR(10,R,G,B)
      call pgsci(10)
      CALL PGPOLY(100,x,y)

c     stile linea
      CALL pgsls(5)
c     larghezza linea
      CALL pgslw(3)
c     colore linea
      call pgsci(14)
      call pgscr(14,0.1,0.6,0.3)

      CALL gnom_projection(1,ra_center,dec_center,raarr(3),decarr(3),x,y)
      x_err2=x(1)
      y_err2=y(1)
      IF(ra_center.NE.ra_err2.OR.dec_center.NE.dec_err2) THEN
      CALL gnom_ellipse (100,dec_center,ellipser1_2,ellipser2_2,ellipserot_2,x_err2,y_err2,x,y)
      CALL pgline(100,x,y)
      ELSE
      CALL gnom_ellipse (100,dec_center,ellipser1_2,ellipser2_2,ellipserot_2,0.0,0.0,x,y)
      CALL pgline(100,x,y)
      ENDIF
      CALL PGSFS(1)
      CALL PGSCI(10)
      call pgsci(10)
      CALL PGPOLY(100,x,y)

      CALL gnom_projection(1,ra_center,dec_center,raarr(1),decarr(1),x,y)
      call pgsci(2)
      call pgsch(1.5)
      call pgpoint(1,x,y,2)

c     now draw the grid
      CALL pgsch(0.6)
      CALL pgsci(1)
      CALL pgslw(3)
      step_alpha = (ra2-ra1)/float(no_of_isoalpha)
      DO n = 0, no_of_isoalpha
         isoalpha = ra_center + (float(n)-no_of_isoalpha/2.)*step_alpha
         pos=.9-float(n+1)/30.-.02
         rra=isoalpha
         if (rra.lt.0.) rra=rra+360.
         string=' '
         call schra(rra,string)
         newstring='   '//string(1:lenact(string))
         DO m = 1, 100
            run_dec(m) = dec_start + float(m-3)*radius/47.
            run_alpha(m) = isoalpha
         END DO
         CALL gnom_projection(100,ra_center,dec_center,run_alpha,run_dec,
     &                       x_grid,y_grid)
         CALL pgline(100,x_grid,y_grid)
      END DO
      step_delta = 2.*radius/float(no_of_isodelta)
      DO n = 0, no_of_isodelta
      isodelta = dec_center + (float(n)-no_of_isodelta/2.)*step_delta
      IF ( abs(isodelta).LE.90. ) THEN
      pos=.6-float(n+1)/30.-.02
      ddec=isodelta
      call schdec(ddec,string)
      newstring='  '//string(1:lenact(string))
      DO m = 1, 100
      run_alpha(m) = ra1 + float(m-3)*(ra2-ra1)/95.
      run_dec(m) = isodelta
      END DO
      CALL gnom_projection(100,ra_center,dec_center,run_alpha,
     &                         run_dec,x_grid,y_grid)
      CALL pgline(100,x_grid,y_grid)
      END IF
      END DO
      rra=ra_center
      ddec=dec_center
      CALL chra(rra,rah,irm,rasec,1)
      CALL chdec(ddec,id,idm,decsec,1)
      sign = '+'
      IF (ddec < 0. ) sign = '-'
      write(title,'(a,I2.2,1x,I2.2,1x,f4.1,a,a,
     &  I2.2,1x,I2.2,1x,f4.1)')
     & 'Image centre R.A.=',rah,irm,rasec,' Dec.=',
     &  sign,abs(id),abs(idm),abs(decsec)
      xaxis_label = 'R.A. (degrees)'
      yaxis_label = 'Dec. (degrees)'
      CALL pgsch(1.2)
      CALL pglabel (xaxis_label, yaxis_label,title)
      CALL gnom_projection(1,ra_center,dec_center,raarr(2),decarr(2),x,y)
c     colore
      CALL pgsch(1.)
      CALL pgsci(2)
      CALL pgpoint(1,x,y,5)
      fmax=afmax(1)
      fmin=afmin(1)
      call pgsci(1)
      CALL pgsls(3)
      CALL pgslw(3)
      IF ( err_radius .GT. 60. ) THEN
      err_radius=60.
      call pgsci(2)
      ELSE
      call pgsci(4)
      ENDIF

      call pgsci(1)
      CALL pgscr(10,.7,.2,.7)
      if (n_points .gt. 35000) n_points=35000

      deallocate(run_alpha,run_dec,x_grid,y_grid,afmin,afmax)

      allocate(ra_col13(arrsize(3)),dec_col13(arrsize(3)),epos_col13(arrsize(3)),color13(arrsize(3)),csr13(arrsize(3)))
      allocate(ra_col12(arrsize(1)+arrsize(2)),dec_col12(arrsize(1)+arrsize(2)),csx12(arrsize(1)+arrsize(2)),csr12(arrsize(1)+arrsize(2)))
      allocate(tcol12(arrsize(1)+arrsize(2)),s12(arrsize(1)+arrsize(2)))
      allocate(ra_col14(arrsize(10)),dec_col14(arrsize(10)))
      allocate(tcol14(arrsize(10)),s14(arrsize(10)))

      allocate(tcol1(arrsize(4)),tcol2(arrsize(5)),tcol3(arrsize(6)),tcol4(arrsize(7)),tcol5(arrsize(8)),tcol11(arrsize(9)))
      allocate(ra_col1(arrsize(4)),dec_col1(arrsize(4)),csx1(arrsize(4)),csr1(arrsize(4)))
      allocate(ra_col2(arrsize(5)),dec_col2(arrsize(5)),csx2(arrsize(5)),csr2(arrsize(5)))
      allocate(ra_col3(arrsize(6)),dec_col3(arrsize(6)),csx3(arrsize(6)),csr3(arrsize(6)))
      allocate(ra_col4(arrsize(7)),dec_col4(arrsize(7)),csx4(arrsize(7)),csr4(arrsize(7)))
      allocate(ra_col5(arrsize(8)),dec_col5(arrsize(8)),csx5(arrsize(8)),csr5(arrsize(8)))
      allocate(ra_col11(arrsize(9)),dec_col11(arrsize(9)),csx11(arrsize(9)),csr11(arrsize(9)),s11(arrsize(9)))

      DO j = 1,n_points
         IF ((code(j) .GT. 10000) .or. (code(j) .LT. -40000)) isource=isource+1
         IF (code(j) .EQ. -9999) isource=isource+1
         if ((code(j) .GT. 10000) .or. (code(j) .LE. -10000)) then
            om = int(code(j)/10000.)
            code(j) = abs(code(j))
            cxray  = int((code(j)-abs(om)*10000.)/100.)
            cradio = code(j)-abs(om)*10000.-cxray*100.
            cs = max(1.0,cradio*8./99.)
            cx = max(1.0,cxray*8./99.)
            IF (om.EQ.1) THEN
               icol1 = icol1 +1
               ra_col1(icol1)=ra(j)
               dec_col1(icol1)=dec(j)
               csr1(icol1)= cs*0.8
               csx1(icol1)= cx*0.7
               write(tcol1(icol1),'(i4)') isource
            ELSE IF (om.EQ.2) THEN
               icol2 = icol2 +1
                  ra_col2(icol2)=ra(j)
                  dec_col2(icol2)=dec(j)
                  csr2(icol2)= cs*0.8
                  csx2(icol2)= cx*0.7
                  write(tcol2(icol2),'(i4)') isource
            ELSE IF (om.EQ.3) THEN
                  icol3 = icol3 +1
                  ra_col3(icol3)=ra(j)
                  dec_col3(icol3)=dec(j)
                  csr3(icol3)= cs*0.8
                  csx3(icol3)= cx*0.7
                  write(tcol3(icol3),'(i4)') isource
            ELSE IF (om.EQ.4) THEN
                  icol4 = icol4 +1
                  ra_col4(icol4)=ra(j)
                  dec_col4(icol4)=dec(j)
                  csr4(icol4)= cs*0.8
                  csx4(icol4)= cx*0.7
                  write(tcol4(icol4),'(i4)') isource
            ELSE IF (om.LT.0) THEN        ! case of catalogued sources
c PG 
               IF (om > -8) THEN !catalog blazars/AGN
                  icol11 = icol11 +1
                  ra_col11(icol11)=ra(j)
                  dec_col11(icol11)=dec(j)
                  if (om .lt. -4) then
                     write(tcol11(icol11),'(i4)') isource
                  else
                     write(tcol11(icol11),'(a)') "    "
                  ENDIF
                  IF (om .lt. -4) om=om+4
                  IF (om .EQ. -1 ) THEN      ! WHSP
                     s11(icol11) = 12
                     csx11(icol11)= 1.6
                  ELSE IF (om .EQ. -4) THEN  ! Cluster of galaxies
                     s11(icol11) = 63
                     csx11(icol11)= 1.4
                  ELSE IF (om .EQ. -2) THEN  ! BZCAT not correspoponding to X-ray source
                     s11(icol11) = -4
                     csx11(icol11)= 0.9
                  ELSE IF (om .EQ. -3) THEN  ! CRATES
                     s11(icol11) = 6
                     csx11(icol11)= 1.2
                     csr11(icol11)= cs*0.8
                     if (cradio .gt. 0.) THEN
c add radio counterparts/ extra radio counpornents, remove crates components
                        icol12=icol12+1
                        icol11=icol11-1
                        ra_col12(icol12)=-ra(j)
                        dec_col12(icol12)=dec(j)
                        s12(icol12) = 17
                        csr12(icol12)= cs*0.8!*0.9
                        csx12(icol12)= 0.
                     endif
                  ENDIF
               ELSE
                  icol12 = icol12 +1 !single radio and X-ray sources
                  ra_col12(icol12)=ra(j)
                  dec_col12(icol12)=dec(j)
                  write(tcol12(icol12),'(i4)') isource
                  IF (om == -8) THEN  ! Simple X-ray  source
                     s12(icol12) = 21
                     csr12(icol12)= 0.
                     csx12(icol12)= cx*0.45
                  ELSE IF (om == -9) THEN  ! Simple radio source
                     s12(icol12) = 17
                     csr12(icol12)= cs*0.55
                     csx12(icol12)= 0.
                  ENDIF
               ENDIF
            ELSE                          ! om = 5: Unknown source type
               icol5 = icol5 +1
               ra_col5(icol5)=ra(j)
               dec_col5(icol5)=dec(j)
               csr5(icol5)= cs
               csx5(icol5)= cx
               write(tcol5(icol5),'(i4)') isource
            ENDIF
         else !((code(i) .GT. 0 ) .and. (code(i) .LT. 10000 )) then
            om = int(code(j)/100.)
            if (om .gt. 0) then !error circle mode
               icol13 = icol13 + 1
               ra_col13(icol13)=ra(j)
               dec_col13(icol13)=dec(j)
               epos_col13(icol13)=epos(j)
               cradio= code(j)-abs(om)*100.
c              cs = max(1.0,cradio*8./99.)
               cs = max(3.0,cradio*8./99.)
               csr13(icol13)= cs*0.8
               color13(icol13)=om

            else !pulsar, GRB, and gamma-ray sources
               icol14 = icol14 + 1
               ra_col14(icol14)=ra(j)
               dec_col14(icol14)=dec(j)
               if ( (om .eq. -11) .or. (om .eq. -12)) then
                  s14(icol14)=7 !for gamma-ray source
               else if (om .eq. -22) then ! for GRB
                  s14(icol14)=10
               else if (om .eq. -33) then ! for 4FGL-radio
                  s14(icol14)=-3
               else if (om .eq. -44) then ! for 4LAC-DR3
                  s14(icol14)=25 
               else if ((om .eq. -70) .or. (om .eq. -77)) then !quasar
                  s14(icol14)=20
               else
                  s14(icol14)=-5 !for pulsar
               endif
               !write(*,*) isource
               if ((om .eq. -99 )) then !for extra pulsar and GRB
                  write(tcol14(icol14),'(i4)') isource
               else
                  write(tcol14(icol14),'(a)') "    "
               ENDIF
               !write(*,*) icol14!,ra_col14(icol14),dec_col14(icol14),tcol14(icol14)
            endif
         endif
      ENDDO
c      write(*,*) icol1,icol2,icol3,icol4,icol5,icol11,icol12,icol13,icol14
c      icol14=0
      !write(*,*) 'CENTER',ra_center,dec_center

      deallocate(ra,dec,code,epos)

      IF (icol1.GT.0) THEN 
        DO j = 1,icol1
           CALL gnom_projection(1,ra_center,dec_center,ra_col1(j),dec_col1(j),x,y)
           R=110./255.
           G=87./255.
           B=36./255.
c           CALL PGSCR(8,R,G,B)
           CALL pgsci(8)
           CALL pgsch(csr1(j))
           CALL pgpoint(1,x,y,17)
           IF (csr1(j).GT.0.8*csx1(j)) CALL pgsci(10)
           CALL pgsch(csx1(j))
           CALL pgpoint(1,x,y,21)
           CALL pgsci(8)
           call pgsch(0.8)
           CALL PGTEXT (X, Y, tcol1(j))
           write (lu_out,'(4(1x,f10.4),a)') ra_col1(j),dec_col1(j),x(1),y(1),tcol1(j)
        ENDDO
      ENDIF
      IF (icol2.GT.0) THEN 
        DO j = 1,icol2
           CALL gnom_projection(1,ra_center,dec_center,ra_col2(j),dec_col2(j),x,y)
           CALL pgsci(5)
           CALL pgsch(csr2(j))
           CALL pgpoint(1,x,y,17)
           IF (csr2(j).GT.0.8*csx2(j)) CALL pgsci(10)
           CALL pgsch(csx2(j))
           CALL pgpoint(1,x,y,21)
           CALL pgsci(5)
           call pgsch(0.8)
           CALL PGTEXT (X, Y, tcol2(j))
           write (lu_out,'(4(1x,f10.4),a)') ra_col2(j),dec_col2(j),x(1),y(1),tcol2(j)
        ENDDO
      ENDIF
      IF (icol3.GT.0) THEN 
        DO j = 1,icol3
           CALL gnom_projection(1,ra_center,dec_center,ra_col3(j),dec_col3(j),x,y)
           CALL pgsci(4)
           CALL pgsch(csr3(j))
           CALL pgpoint(1,x,y,17)
           IF (csr3(j).GT.1.2*csx3(j)) CALL pgsci(10)
           CALL pgsch(csx3(j))
           CALL pgpoint(1,x,y,21)
           CALL pgsci(4)
           call pgsch(0.8)
           CALL PGTEXT (X, Y, tcol3(j))
           write (lu_out,'(4(1x,f10.4),a)') ra_col3(j),dec_col3(j),x(1),y(1),tcol3(j)
        ENDDO
      ENDIF
      IF (icol4.GT.0) THEN 
        DO j = 1,icol4
           CALL gnom_projection(1,ra_center,dec_center,ra_col4(j),dec_col4(j),x,y)
           R=45./255.
           G=73./255.
           B=62./255.
           call pgscr(3,R,G,B)
           CALL pgsci(3)
           CALL pgsch(csr4(j))
           CALL pgpoint(1,x,y,17)
           IF (csr4(j).GT.1.2*csx4(j)) CALL pgsci(10)
           CALL pgsch(csx4(j))
           CALL pgpoint(1,x,y,21)
           CALL pgsci(3)
           call pgsch(0.8)
           CALL PGTEXT (X, Y, tcol4(j))
           write (lu_out,'(4(1x,f10.4),a)') ra_col4(j),dec_col4(j),x(1),y(1),tcol4(j)
        ENDDO
      ENDIF
      IF (icol5.GT.0) THEN 
        DO j = 1,icol5
           CALL gnom_projection(1,ra_center,dec_center,ra_col5(j),dec_col5(j),x,y)
           R=43./255.
           G=32./255.
           B=34./255.
           call pgscr(17,R,G,B)
           CALL pgsci(17)
           CALL pgsch(csr5(j))
           CALL pgpoint(1,x,y,17)
           IF (csr5(j).GT.1.2*csx5(j)) CALL pgsci(10)
           CALL pgsch(csx5(j))
           CALL pgpoint(1,x,y,21)
           CALL pgsci(17)
           call pgsch(0.8)
           CALL PGTEXT (X, Y, tcol5(j))
           write (lu_out,'(4(1x,f10.4),a)') ra_col5(j),dec_col5(j),x(1),y(1),tcol5(j)
        ENDDO
      ENDIF

      IF (icol12.GT.0) THEN
         DO j = 1,icol12
            CALL gnom_projection(1,ra_center,dec_center,abs(ra_col12(j:j)),dec_col12(j),x,y)
            IF (csx12(j) > 0. ) THEN
               CALL pgsch(csx12(j))
               CALL pgsci(4)
            ENDIF
            IF (csr12(j) > 0. ) THEN
               CALL pgsch(csr12(j))
               R = 227./255.
               G = 11./255.
               B = 93./255.
               CALL PGSCR(2,R,G,B)
               CALL pgsci(2)
            ENDIF
            CALL pgpoint(1,x,y,s12(j))
            call pgsch(0.8)
            if ((filein(iskip+1:iskip+17) == 'find_out_temp.txt') .and. (ra_col12(j).gt. 0.)) then
               call pgtext(x,y,tcol12(j))
               write (lu_out,'(4(1x,f10.4),a)') ra_col12(j),dec_col12(j),x(1),y(1),tcol12(j)
            endif
         ENDDO
      ENDIF
      IF (icol11.GT.0) THEN
        DO j = 1,icol11
           CALL gnom_projection(1,ra_center,dec_center,ra_col11(j),dec_col11(j),x,y)
           if (s11(j) .eq. 6) then
              call pgsci(11)
           else
              R = 211./255.
              G = 179./255.
              B = 102./255.
              CALL PGSCR(7,R,G,B)
              call pgsci(7)
           endif
           CALL pgsch(csx11(j))
           CALL pgpoint(1,x,y,s11(j))
           if (tcol11(j) .ne. "    ") then
              call pgsch(0.8)
              CALL PGTEXT (X, Y, tcol11(j))
           write (lu_out,'(4(1x,f10.4),a)') ra_col11(j),dec_col11(j),x(1),y(1),tcol11(j)
           endif
        ENDDO
      ENDIF
      IF (icol13 .GT. 0) THEN
         do j=1,icol13
           CALL gnom_projection(1,ra_center,dec_center,ra_col13(j),dec_col13(j),x,y)
           CALL pgsch(csr13(j))
           if (color13(j) .lt. 20) CALL pgsci(2)
           if ((color13(j) .gt. 50) .and. (color13(j) .lt. 60)) CALL pgsci(8)
           if ((color13(j) .gt. 60) .and. (color13(j) .lt. 70)) CALL pgsci(7)
           if ((color13(j) .gt. 70) .and. (color13(j) .lt. 80)) CALL pgsci(9)
           if ((color13(j) .gt. 80) .and. (color13(j) .lt. 90)) CALL pgsci(4)
           if (color13(j) .gt. 90) CALL pgsci(12)
C - PG 
           R = 211./255.
           G = 179./255.
           B = 102./255.
           CALL PGSCR(7,R,G,B)
           !CALL pgpoint(1,x,y,21)
           R = 255./255.
           G = 150./255.
           B = 0./255.
           CALL PGSCR(8,R,G,B)
C- PG
           CALL pgpoint(1,x,y,21)
           !if (color13(i) .lt. 20) then
              !CALL pgsci(2)
              xx=x(1)
              yy=y(1)
              CALL gnom_circle_off (100,dec_center,epos_col13(j)/60.,x,y,xx,yy)
              CALL pgline(100,x,y)
           !endif
         enddo
      endif
c      write(*,*) 'number of cat. 14 sources',icol14
      if (icol14 .gt. 0) then
         do j=1,icol14
            !write(*,*) ra_col14(j),dec_col14(j),j,s14(j)
            CALL gnom_projection(1,ra_center,dec_center,ra_col14(j),dec_col14(j),x,y)
            if (s14(j) .eq. 7) then
               call pgsch(1.5)
c            else if (s14(j) .eq. 20) then
c               call pgsch(4.)
            else
               call pgsch(1.3)
            endif
            if (s14(j) .eq. 10) then
               call pgsci(3)
            else if (s14(j) .eq. 20) then
               R = 32./255.
               G = 149./255.
               B = 83./255.
               CALL PGSCR(15,R,G,B)
               call pgsci(15)
            else
               call pgsci(12)
            endif
c            write(*,*) j,s14(j),tcol14(j)
            call pgpoint(1,x,y,s14(j))
            if (tcol14(j) .ne. "    ") then
               call pgsch(0.8)
               CALL PGTEXT (X, Y, tcol14(j))
               write (lu_out,'(4(1x,f10.4),a)') ra_col14(j),dec_col14(j),x(1),y(1),tcol14(j)
            endif
         enddo
      endif
c      CALL pgscf(1)

      deallocate(ra_col13,dec_col13,epos_col13,color13,csr13)
      deallocate(ra_col12,dec_col12,csx12,csr12)
      deallocate(tcol12,s12)
      deallocate(ra_col14,dec_col14)
      deallocate(tcol14,s14)

      deallocate(tcol1,tcol2,tcol3,tcol4,tcol5,tcol11)
      deallocate(ra_col1,dec_col1,csx1,csr1)
      deallocate(ra_col2,dec_col2,csx2,csr2)
      deallocate(ra_col3,dec_col3,csx3,csr3)
      deallocate(ra_col4,dec_col4,csx4,csr4)
      deallocate(ra_col5,dec_col5,csx5,csr5)
      deallocate(ra_col11,dec_col11,csx11,csr11,s11)
      deallocate(x,y)

c ------ DC2Truth part
cccc      CALL gnom_projection(n_true,ra_center,dec_center,ratrue,dectrue,x1,y1)
cccc      call pgsci(2)
cccc      call pgsch(2.)
c      CALL pgpoint(n_true,x1,y1,18)
c
c  ----  Catalogues sources part 
c
cccc      CALL gnom_projection(n_cat,ra_center,dec_center,racat,deccat,x1,y1)
cccc      call pgsci(1)
cccc      call pgsch(2.)
c      CALL pgpoint(n_cat,x1,y1,4)

c ------ Names part
cccc      CALL gnom_projection(n_names,ra_center,dec_center,names_ra,names_dec,x1,y1)
cccc      call pgsci(0)
cccc      call pgsch(1.)
c      CALL pgpoint(n_names,x1,y1,1)

c      deallocate(x1,y1,racat,deccat)

      close(lu_out)
      END

**==GNOM_PROJECTION.FOR
 
      SUBROUTINE gnom_projection(npoints,ra_center,dec_center,ra,dec,x,
     &                           y)

      IMPLICIT NONE
      INTEGER npoints, i
      REAL*4 ra_center, dec_center, ra(*), dec(*), x(*), y(*)
      REAL*4 a, f, scal, radian, dd
c      CHARACTER*10 axis_unit
      REAL*4 radius_scale
      REAL*4 pro_scale
c      COMMON /ecetype/pro_scale,radius_scale,axis_unit
      COMMON /ecetype/pro_scale,radius_scale

      scal = 60.
      radian = 45.0/atan(1.0)

      DO i = 1, npoints
        if (ra(i) .lt. 0) ra(i)=-ra(i)
        dd = dec(i)/radian
        a = cos(dd)*cos((ra_center-ra(i))/radian)
        f = scal*radian/(sin(dec_center/radian)*sin(dd)
     &      +a*cos(dec_center/radian))
c        x(i) = -f*cos(dd)*sin((ra(i)-ra_center)/radian) / pro_scale
        x(i) = -f*cos(dd)*sin((ra(i)-ra_center)/radian) / pro_scale / cos(dd)
        y(i) = f*(cos(dec_center/radian)*sin(dd)-a*sin(dec_center/radian 
     &         )) / pro_scale
      END DO
      RETURN
      END

      SUBROUTINE gnom_circle (npoints,dec,radius,x,y)
      IMPLICIT NONE
      INTEGER*4 npoints, i
      REAL*4 alpha, delta_alpha, radian, radius, x(*), y(*),dec
      
      radian = 45.0/atan(1.0)
      delta_alpha = 360./float(npoints-1)
      alpha=-delta_alpha
      DO i =1,npoints
        alpha=alpha+delta_alpha
c        x(i) = radius*cos(alpha/radian)
        x(i) = radius*cos(alpha/radian)/cos(dec/radian)
        y(i) = radius*sin(alpha/radian)
      ENDDO 
      RETURN 
      END

      SUBROUTINE gnom_circle_off (npoints,dec,radius,x,y,xo,yo)
      IMPLICIT NONE
      INTEGER*4 npoints, i
      REAL*4 xo, yo
      REAL*4 alpha, delta_alpha, radian, radius, x(*), y(*),dec
      
      radian = 45.0/atan(1.0)
      delta_alpha = 360./float(npoints-1)
      alpha=-delta_alpha
      DO i =1,npoints
        alpha=alpha+delta_alpha
        x(i) = radius*cos(alpha/radian)/cos(dec/radian) + xo
        y(i) = radius*sin(alpha/radian) + yo
      ENDDO 
      RETURN 
      END


      SUBROUTINE gnom_ellipse (npoints,dec,r1,r2,rot,off1,off2,x,y)
      IMPLICIT NONE
      INTEGER*4 npoints, i
      REAL*4 alpha, delta_alpha, p2, r1, r2, rot, x(*), y(*), x1, y1
      REAL*4 off1, off2, dec, radian
      
      radian = 45.0/atan(1.0)
      p2 = 8.0*atan(1.0)
      delta_alpha = p2/float(npoints-1)
      alpha=0

      DO i =1,npoints
        alpha=alpha+delta_alpha
        x1 = r2*cos(alpha)
        y1 = r1*sin(alpha)

c       rotazione e traslazione
        x(i) = (x1*cos(-p2*rot/360.) - y1*sin(-p2*rot/360.) + off1)/cos(dec/radian)
        y(i) = x1*sin(-p2*rot/360.) + y1*cos(-p2*rot/360.) + off2

      ENDDO 
      RETURN 
      END
