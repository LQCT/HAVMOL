      subroutine incar (words, NWORDS)
      implicit real*8 (a-h,o-z)
      include 'dimmm'
      common/work/item,nitems,istrt(40),inumb(40),iwidth
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),
     4               KR(NATM1),KA(NATM1),KB(NATM1),
     5               PARVAL(NATM13),ZIN(NATM1),X(NATM13)
      logical*1 qbasfil
      character*80 kwd,anl
      character*1 iachar, achar
      character*8 buffa, end
      character*32 file5
      common
     ./optn1/ kwd(48)
     ./optn2/ icar, kb1, iii
     ./files / file5
     ./digits/ iachar(12), achar(13), end, buffa
      common/elemts/ atnam(37)
      character*4 atnam
      character*32 arg
      logical*2 qcomlin
      integer*2 marg,sta,numarg
      common /cl/ qcomlin,numarg,marg(5),sta,arg
      common/thermo/cpt,cpr,cpv,cpe,
     &              h25t,h25r,h25v,h25e,
     &              s25t,s25r,s25v,s25e,
     &              g25t,g25r,g25v,g25e,
     &              ilin,nax,
     &              temper
      character*8 fin /'FINISH  '/, runtype(7), words(*), filebas
      character*1 an

      kb1 = 1
      qbasfil = .false.
      lruntype = 0
      read(7,'(i4)') nat
      read(7,'(a)') anl
      write (kwd(1),1001) 'TITLE'
      write (kwd(2),1001) anl

      k = 2
100   call readc
110   if (item.eq.nitems) goto 120
        call reada(buffa)
        do j=1,NWORDS
          if (buffa.eq.words(j)) then
            goto
     1           (1,  1,  2,  1,  1,
     2            1,  1,  4,  5,  6,
     3            7,  1,  1,  1,  1,
     4            1,  1,  1,  8,  9,
     5            1,  1, 11,  1,  1,
     6            1, 12, 13,  1,  1,
     7            1,  1,  1,120,  1,
     8            1,  1,  1,  1, 14,
     9            1,  1,  1, 15, 16,
     A           17,  1, 19,  1,  1,
     B            1,  1,  1,  1, 25,
     C           26,  1,  1,  1,  1,
     d           20,  1, 20, 20, 20,
     e           20, 20, 20,20,   1), j
          endif
        enddo
        if (buffa.eq.'&' .or. buffa.eq.'#') goto 100
        if (buffa.eq.'BASLIB') goto 21
1       write(10,1000) buffa, ' item not allowed in .CAR input file'
        call honder(75,'INCAR ')

c PRINT 

2       call readi(nprint)
        call wrmemi(k,buffa,nprint)
        goto 110

c CHARGE

4       call readi(ich)
        call wrmemi(k,buffa,ich)
        goto 110

c MULTIPLICITY

5       call readi(imul)
        call wrmemi(k,buffa,imul)
        goto 110

c SCFTYP

6       call reada(buffa)
7       call wrmem2(k,'SCFTYP  ',buffa)
        goto 110

c MAXIT

8       call readi(maxit)
        call wrmemi(k,buffa,maxit)
        goto 110
  
c CONVERGENCE

9       call readi(nconv)
        call wrmemi(k,buffa,nconv)
        goto 110

c RUNTYP

11      call reada(buffa)
20      lruntype = lruntype + 1
        runtype(lruntype) = buffa
        goto 110

c TIME

12      call readf(ttime)
        k = k + 1
        write (kwd(k),1003) buffa,ttime,'H'
        goto 110

c SKIP

13      call wrmem1(k,buffa)
        goto 110

c GRADTOL

14      call readf(exgrad)
        call wrmemf(k,buffa,exgrad)
        goto 110

c STEPSIZE

15      call readf(vibsiz)
        call wrmemf(k,buffa,vibsiz)
        goto 110

c NPOINTS

16      call readi(nvib)
        call wrmemi(k, buffa,nvib)
        goto 110

c STOP

17      fin = '        '
        fin = buffa
        goto 110

c HUCKEL

19      call wrmem2(k,'GUESS   ',buffa)
        goto 110

c Basis set filename

21      call reada(buffa)
        filebas = buffa
        qbasfil = .true.
        goto 110

c     symmetry axis for rotational calculations

25      call readi(nax)
        call wrmemi(k,buffa,nax)
        goto 110

C     reading temperature for thermodynamic properties

26      call readf(temper)
        call wrmemf(k,buffa,temper)
        goto 110

120   continue
      
      call wrmem1(k,'ANGSTROM')
      call wrmem1(k,'ATOMS   ')
      if (.not.qbasfil) filebas = 'MINI1'

      do i=1,nat
        read(7,'(f11.0,2f10.0,i4)') x0,y0,z0,inat
        k = k + 1
        write (kwd(k),1005) atnam(inat), inat, x0, y0, z0
        call wrmem1(k,filebas)
      enddo
      icar = 0
      close(7)
      call wrmem1(k,end)
      if (lruntype.eq.0) then
        lruntype = 1
        runtype(1) = words(61)
      endif
      do i=1,lruntype
        call wrmem2(k,'RUNTYP  ',runtype(i))
        call wrmem1(k,'START   ')
      enddo
      call wrmem1(k,fin)

      if (qcomlin) then
        if (numarg.eq.5) then
          call rcomlin(5,'INCAR ')
          if (arg.eq.'IMMOUT') then
            an='Y'
          else
            an='N'
          endif
        else
          an='Y'
        endif
      else
        write(10,'(/a\)')
     &         ' Do you want to create a <fn>[.IMM] file? (Y/[N]) '
        call upcase(10,an)
      endif
      if (an.eq.'Y') then
        call filen1(1,file5,'.IMM')
        open (21,file=file5,status='UNKNOWN')
        do i=1,k
          write (21,'(1x,a)') kwd(i)
        enddo
        close (21)
      endif
      return

1000  format (/a,1x,a)
1001  format (a,10(1x,a))
1003  format (a,f10.4,1x,a)
1005  format (a4,i5,3f15.5)

      end

      subroutine intern (words, NWORDS)
      implicit real*8 (a-h,o-z)
      include 'dimmm'
      common/work/item,nitems,istrt(40),inumb(40),iwidth
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),
     4               KR(NATM1),KA(NATM1),KB(NATM1),
     5               PARVAL(NATM13),ZIN(NATM1),X(NATM13)
      character*32 arg
      logical*2 qcomlin
      integer*2 marg,sta,numarg
      common /cl/ qcomlin,numarg,marg(5),sta,arg
      logical*1 qbasfil
      character*80 kwd,anl
      common
     ./optn1/ kwd(48)
     ./optn2/ icar, kb1, iii
     ./files / file5
     ./digits/ iachar(12), achar(13), end, buffa
      common/thermo/cpt,cpr,cpv,cpe,
     &              h25t,h25r,h25v,h25e,
     &              s25t,s25r,s25v,s25e,
     &              g25t,g25r,g25v,g25e,
     &              ilin,nax,
     &              temper
      character*32 file5
      parameter (NWC=30)
      character*8 buffa, end,
     &            fin /'FINISH'/, runtype(7), words(*), filebas,
     &            wordsc(NWC)
     1 /'>PRINT'  , '>CHARGE' , '>MULTIPL', '>SCFTYP' , '>OSCF'   ,
     2  '>MAXIT'  , '>CONVERG', '>RUNTYP' , '>TIME'   , '>SKIP'   ,
     3  '>GRADTOL', '>STEPSIZ', '>NPOINTS', '>STOP'   , '>HUCKEL' ,
     4  '>AXIS'   , '>TEMPERA', '        ', '        ', '        ',
     5  '>HFSCF'  , '>CHECK'  , '>INTEGRA', '>SCF'    , '>GRADIEN',
     6  '>GRADONE' ,'>OPTIMIZ', '>FORCE'  , '>PROPERT', '>ROTATIO'/

      common/elemts/ atnam(37)
      character*4 atnam
      character*1 iachar, achar
      character*1 an

      kb1 = 1
      qbasfil = .false.
      lruntype = 0

      k = 0
100   call readc
110   if (item.eq.nitems) goto 120
        call reada(buffa)
        if (buffa.eq.'>') goto 110
        do j=1,NWORDS
          if (buffa.eq.words(j)) then
            goto
     1           (1,  1,  2,  1,  1,
     2            1,  1,  4,  5,  6,
     3            7,  1,  1,  1,  1,
     4            1,  1,  1,  8,  9,
     5            1,  1, 11,  1,  1,
     6            1, 12, 13,  1,  1,
     7            1,  1,  1,120,  1,
     8            1,  1,  1,  1, 14,
     9            1,  1,  1, 15, 16,
     A           17,  1, 19,  1,  1,
     B            1,  1,  1,  1, 25,
     C           26,  1,  1,  1,  1,
     d           20,  1, 20, 20, 20,
     e           20, 20, 20,20,   1), j
          endif
        enddo
        do j=1,NWC
          if (buffa.eq.wordsc(j)) then
            goto
     1           (2,  4,  5,  6,  7,
     2            8,  9, 10, 12, 13,
     3           14, 15, 16, 17, 19,
     4           25, 26,  1,  1,  1,
     5           30, 30, 30, 30, 30,
     6           30, 30, 30, 30, 30),j
          endif
        enddo
        if (buffa.eq.'&' .or. buffa.eq.'#') goto 100
        if (buffa.eq.'BASLIB' .or. buffa.eq.'>BASLIB') goto 21

1       write(10,1000)
     & buffa, ' item not allowed in internal coordinate input file'
        call honder(76,'INTERN')

c PRINT 

2       call readi(nprint)
        call wrmemi(k,words(3),nprint)
        goto 110

c CHARGE

4       call readi(ich)
        call wrmemi(k,words(8),ich)
        goto 110

c MULTIPLICITY

5       call readi(imul)
        call wrmemi(k,words(9),imul)
        goto 110

c SCFTYP

6       call reada(buffa)
7       call wrmem2(k,'SCFTYP  ',buffa)
        goto 110

c MAXIT

8       call readi(maxit)
        call wrmemi(k,words(19),maxit)
        goto 110
  
c CONVERGENCE

9       call readi(nconv)
        call wrmemi(k,words(20),nconv)
        goto 110

c RUNTYP

11      call reada(buffa)
20      lruntype = lruntype + 1
        runtype(lruntype) = buffa
        goto 110

c >RUNTYP

10      call reada(buffa)
        lruntype = lruntype + 1
        runtype(lruntype) = buffa
        goto 110

30      lruntype = lruntype + 1
        runtype(lruntype) = words((NWORDS-20)+(j-20))
        goto 110

c TIME

12      call readf(ttime)
        k = k + 1
        write (kwd(k),1003) words(27),ttime,'H'
        goto 110

c SKIP

13      call wrmem1(k,words(28))
        goto 110

c GRADTOL

14      call readf(exgrad)
        call wrmemf(k,words(40),exgrad)
        goto 110

c STEPSIZE

15      call readf(vibsiz)
        call wrmemf(k,words(44),vibsiz)
        goto 110

c NPOINTS

16      call readi(nvib)
        call wrmemi(k, words(45),nvib)
        goto 110

c STOP

17      fin = '        '
        fin = words(46)
        goto 110

c HUCKEL

19      call wrmem2(k,'GUESS   ',words(48))
        goto 110

c Basis set filename

21      call reada(buffa)
        filebas = buffa
        qbasfil = .true.
        goto 110

c     symmetry axis for rotational calculations

25      call readi(nax)
        call wrmemi(k,words(55),nax)
        goto 110

C     reading temperature for thermodynamic properties

26      call readf(temper)
        call wrmemf(k,words(56),temper)
        goto 110

c TITLE

120   continue

      k = k + 1
      write (kwd(k),1001) 'TITLE'
      read(7,'(a)') anl
      k = k + 1
      write (kwd(k),1001) anl

      call wrmem1(k,'ANGSTROM')
      call wrmem1(k,'ATOMS   ')
      if (.not.qbasfil) filebas = 'MINI1'

      call getgeo

      do i=1,nat
        inat = int(zan(i))
        k = k + 1
        write (kwd(k),1005) atnam(inat), inat, (c(j,i),j=1,3)
        call wrmem1(k,filebas)
      enddo
      icar = 0
      close(7)
      call wrmem1(k,end)
      if (lruntype.eq.0) then
        lruntype = 1
        runtype(1) = words(61)
      endif
      do i=1,lruntype
        call wrmem2(k,'RUNTYP  ',runtype(i))
        call wrmem1(k,'START   ')
      enddo
      call wrmem1(k,fin)

      if (qcomlin) then
        if (numarg.eq.5) then
          call rcomlin(5,'INTERN')
          if (arg.eq.'IMMOUT') then
            an='Y'
          else
            an='N'
          endif
        else
          an='Y'
        endif
      else
        write(10,'(/a\)')
     &         ' Do you want to create a <fn>[.IMM] file? (Y/[N]) '
        call upcase(10,an)
      endif

      if (an.eq.'Y') then
        call filen1(1,file5,'.IMM')
        open (21,file=file5,status='UNKNOWN')
        do i=1,k
          write (21,'(1x,a)') kwd(i)
        enddo
        close (21)
      endif
      return

1000  format (/a,1x,a)
1001  format (a,10(1x,a))
1003  format (a,f10.4,1x,a)
1005  format (a4,i5,3f15.5)

      end

      subroutine getgeo

      implicit real*8 (a-h,o-z)
      include 'dimmm'
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NAA, NBB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),
     4               KR(NATM1),KA(NATM1),KB(NATM1),
     5               PARVAL(NATM13),ZIN(NATM1),X(3,NATM1)
      common
     ./optn2/ icar, kb1, iii
     ./digits/ iachar(12), achar(13), end, buffa
      character*1 iachar, achar
      character*8 end, buffa
      character*6 suber /'GETGEO'/
      logical*1 qnat3, qend
      integer*2 nn(NATM1), na(NATM1), nb(NATM1), nc(NATM1)

      parameter (DEGRAD=0.01745329252d0, CERO=0.d0,
     &           ATETR=109.47d0, ATRIG=120.d0,
     &           A60=60.d0, A180=180.d0, A240=240.d0, A300=300.d0)
      parameter (UNO=1.d0, COS0=0.9999991d0, DEC=.1d0)

c *** input the trial geometry

c     nat    = the number of atoms in the molecule excluding dummies
c     nuq    = the number of atoms used to define the molecular geometry
c              it includes dummy atoms
c     nn(i)  = the atomic number of atom i .
c            = 99, then the i-th atom is a dummy atom used only to
c              simplify the definition of the molecular geometry.
c     zan(i) = the atom number of atom i . no dummy atoms are included
c              in this array as it is used in the scf computations.
c     x(1,i) = the internuclear separation  in angstroms  between atoms
c              na(i) and (i).
c     x(2,i) = the angle nb(i)-na(i)-(i) input in degrees  stored in
c              radians.
c     x(3,i) = the angle between the vectors nc(i)-nb(i) and na(i)-(i)
c              input in degrees - stored in radians.

c     oputput values of c(i,j) are x, y and z coordinates of each atom

c *** note.    much of this    is not included for the first 3 atoms.
c     atom 1   input nn(1) only.
c     atom 2   input nn(2) and x(1,2)  separation between atoms 1 + 2
c     atom 3   input nn(3), x(1,3)     separation between atoms 2 + 3
c              and x(2,3)               angle  atom 1 -atom 2 -atom 3

      nat=0

      i=0

    3 i=i+1

c        icar=2 the IMC file format is used
c            =3 the AMPAC file free format is used
c            =4 the Z matrix free format is used

      if (icar.eq.2) then
        if (i.eq.1) then
          x(1,1)=CERO
          x(2,1)=CERO
          x(3,1)=CERO
          x(2,2)=CERO
          x(3,2)=CERO
          x(3,3)=CERO
          read(7,'(4i3,3f7.0)') nat1,nat2,nat3,kwik,a12,a13,a23
          nn(1) = nat1
          nn(2) = nat2
          x(1,2) = a12
          qnat3 = nat3.ne.0
          if (qnat3) then
            i = 3
            nat=2
            if (nn(1).ne.99) zan(1)=dble(nn(1))
            if (nn(2).ne.99) zan(2)=dble(nn(2))
            nn(i) = nat3
            go to (31,32),kwik
              a23=ATETR
            go to 32
31            a23=ATRIG
32          continue
            x(1,i) = a13
            x(2,i) = a23
          else
            i = 2
            nat=2
            zan(1)=dble(nn(1))
            zan(2)=dble(nn(2))
            qend = .true.
          endif
        else
          read(7,'(a3,2i3,3x,2i3,3f7.0,3i3)')
     &      buffa,nb(i),na(i),nn(i),ilaz,(x(j,i),j=1,3)
          if (buffa.eq.end) then
            qend = .true.
            goto 1000
          else
            call vali (buffa, nc(i), len(buffa))
          endif
          if (nc(i).le.0) then
            qend = .true.
            go to 1000
          endif
          go to (41,42,43,44,45,46,47,48,49),ilaz
            x(2,i)=ATETR
            x(3,i)=CERO
          go to 49
41          x(2,i)=ATETR
            x(3,i)=A60
          go to 49
42          x(2,i)=ATETR
            x(3,i)=ATRIG
          go to 49
43          x(2,i)=ATETR
            x(3,i)=A180
          go to 49
44          x(2,i)=ATETR
            x(3,i)=A240
          go to 49
45          x(2,i)=ATETR
            x(3,i)=A300
          go to 49
46          x(2,i)=ATRIG
            x(3,i)=CERO
          go to 49
47          x(2,i)=ATRIG
            x(3,i)=A180
          go to 49
48          x(2,i)=A180
            x(3,i)=A180
49        continue
        endif
      elseif (icar.eq.3) then
        call ampinp (i,nn(i),x(1,i),na(i),nb(i),nc(i),qend)
      elseif (icar.eq.4) then
        call zinp (i,nn(i),x(1,i),na(i),nb(i),nc(i),qend)
      endif
1000  if (qend) goto 8
      if(na(i).ge.i) then
        nrefm = na(i)
        goto 1001
      elseif (nb(i).ge.i) then
        nrefm = nb(i)
        goto 1001
      elseif (nc(i).ge.i) then
        nrefm = nc(i)
        goto 1001
      endif
      if(nn(i).ne.99) then
        nat=nat+1
        zan(nat)=dble(nn(i))
        jat1=zan(nat)
      endif

c     convert the angles to radians.

      x(2,i)=x(2,i)*DEGRAD
      x(3,i)=x(3,i)*DEGRAD
      go to 3

8     continue
      nuq = i - 1

c     calculo de las coordenadas cartesianas

      c(1,1)=CERO
      c(2,1)=CERO
      c(3,1)=CERO

      c(1,2)=x(1,2)
      c(2,2)=CERO
      c(3,2)=CERO

      if(nat.eq.2) goto 110

      ccos=cos(x(2,3))
      if(na(3).eq.1)then
         c(1,3)=c(1,1)+x(1,3)*ccos
      else
         c(1,3)=c(1,2)-x(1,3)*ccos
      endif
      c(2,3)=x(1,3)*sin(x(2,3))
      c(3,3)=CERO

      if (nat.eq.3) goto 110

      do 100 i=4,nat
         cosa=cos(x(2,i))
         mb=nb(i)
         ma=na(i)
         xb=c(1,mb)-c(1,ma)
         yb=c(2,mb)-c(2,ma)
         zb=c(3,mb)-c(3,ma)
         rba=UNO/sqrt(xb*xb+yb*yb+zb*zb)
         if (abs(cosa).ge.COS0) then
c
c     atoms ma, mb, and (i) are collinear
c
           rba=x(1,i)*rba*cosa
           c(1,i)=c(1,ma)+xb*rba
           c(2,i)=c(2,ma)+yb*rba
           c(3,i)=c(3,ma)+zb*rba
           go to 100
         endif
c
c     the atoms are not collinear
c
   40    mc=nc(i)
         xa=c(1,mc)-c(1,ma)
         ya=c(2,mc)-c(2,ma)
         za=c(3,mc)-c(3,ma)
c
c     rotate about the z-axis to make yb=0, and xb positive.  if xyb is
c     too small, first rotate the y-axis by 90 degrees.
c
         xyb=sqrt(xb*xb+yb*yb)
         k=-1
         if (xyb.gt.DEC) go to 50
         xpa=za
         za=-xa
         xa=xpa
         xpb=zb
         zb=-xb
         xb=xpb
         xyb=sqrt(xb*xb+yb*yb)
         k=+1
c
c     rotate about the y-axis to make zb vanish
c
   50    costh=xb/xyb
         sinth=yb/xyb
         xpa=xa*costh+ya*sinth
         ypa=ya*costh-xa*sinth
         sinph=zb*rba
         cosph=sqrt(abs(UNO-sinph*sinph))
         xqa=xpa*cosph+za*sinph
         zqa=za*cosph-xpa*sinph
c
c     rotate about the x-axis to make za=0, and ya positive.
c
         yza=sqrt(ypa**2+zqa**2)
         if (yza.lt.2.d-2 ) then
            if (yza.lt.1.d-4) goto 70
            write (6,'(/a)')
     &      ' *** CALCULATION ABANDONED AT THIS POINT ***'
            write(6,'(/a/a/a//a,i3)')
     &' Three atoms being used to define the coordinates of a fourth ato
     &m, whose',
     &' bond-angle is not zero or 180 degreees, are in an almost straigh
     &t line. There',
     &' is a high probability that the coordinates of the atom will be i
     &ncorrect.',
     &' The faulty atom is atom number',i
            write(6,'(//a)')
     &' Cartesian coordinates up to faulty atom:'
            write(6,'(/5x,a1,12x,a1,12x,a1,12x,a1)') 'i','x','y','z'
            do 60 j=1,i
   60         write(6,'(i6,f16.5,2f13.5)') j,(c(k,j),k=1,3)
            write(6,'(/6x,a6,i3,a1,i3,a5,i3,a11,f7.4,a)')
     &' Atoms',ma,',',mb,', and',mc,' are within',yza,
     &' angstroms of a straight line.'
            write (*,'(/a)')
     &      ' *** CALCULATION ABANDONED AT THIS POINT ***'
            write(*,'(/a/a/a//a,i3)')
     &' Three atoms being used to define the coordinates of a fourth ato
     &m, whose',
     &' bond-angle is not zero or 180 degreees, are in an almost straigh
     &t line. There',
     &' is a high probability that the coordinates of the atom will be i
     &ncorrect.',
     &' The faulty atom is atom number',i
            call honder(78,suber)
         endif
         coskh=ypa/yza
         sinkh=zqa/yza
         goto 80
   70    continue
c
c   angle too small to be important
c
         coskh=UNO
         sinkh=CERO
   80    continue
c
c     coordinates :-   a=(xqa,yza,0),   b=(rba,0,0),  c=(0,0,0)
c     none are negative.
c     the coordinates of i are evaluated in the new frame.
c
         sina=sin(x(2,i))
         sind=sin(x(3,i))
         cosd=cos(x(3,i))
         xd=x(1,i)*cosa
         yd=x(1,i)*sina*cosd
         zd=x(1,i)*sina*sind
c
c     transform the coordinates back to the original system.
c
         ypd=yd*coskh-zd*sinkh
         zpd=zd*coskh+yd*sinkh
         xpd=xd*cosph-zpd*sinph
         zqd=zpd*cosph+xd*sinph
         xqd=xpd*costh-ypd*sinth
         yqd=ypd*costh+xpd*sinth
         if (k.lt.1) go to 90
         xrd=-zqd
         zqd=xqd
         xqd=xrd
   90    c(1,i)=xqd+c(1,ma)
         c(2,i)=yqd+c(2,ma)
         c(3,i)=zqd+c(3,ma)
  100 continue
c
c *** now remove the translation vectors, if any, from the array <c>
c
  110 continue

      if (nat.ne.nuq) then
        j=0
        do 200 i=1,nat
          if (nn(i).eq.99) go to 200
          j=j+1
          do 190 k=1,3
  190       c(k,j)=c(k,i)
  200   continue
      endif
      return

 1001 write(6,614) nrefm
      write(*,614) nrefm
      call honder(77,suber)

  614 format(//' ***** INPUT GEOMETRY MISTAKE.'/
     &          '       A non defined atom ',i3,' has been taken as refe
     &rence ...')

      end

      subroutine wrmem1 (k,a)
      character*80 kwd
      common
     ./optn1/ kwd(48)
      character*8 a
      k = k + 1
      write (kwd(k),'(a)') a
      return
      end

      subroutine wrmem2 (k,a1,a2)
      character*80 kwd
      common
     ./optn1/ kwd(48)
      character*8 a1,a2
      k = k + 1
      write (kwd(k),'(2(a,1x))') a1,a2
      return
      end

      subroutine wrmemi (k,a1,i1)
      character*80 kwd
      common
     ./optn1/ kwd(48)
      character*8 a1
      k = k + 1
      write (kwd(k),'(a,1x,i3)') a1,i1
      return
      end

      subroutine wrmemf (k,a1,f1)
      character*80 kwd
      common
     ./optn1/ kwd(48)
      real*8 f1
      character*8 a1
      k = k + 1
      write (kwd(k),'(a,1x,f10.4)') a1,f1
      return
      end

      subroutine vali (ia, ival, ncar)
      implicit real*8 (a-h,o-z)
      implicit integer*2 (i-n)
      integer*4 ncar
      common
     ./digits/ iachar(12), achar(13), end, buffa
      character*8 end, buffa
      character achar, iachar, itemp, ia(80)

      ival=0
      limit=12
      ifact=1
      nstrt=ncar
      do 50 i=1,ncar
        itemp=ia(nstrt)
        do 10 j=1,limit
          if(iachar(j).eq.itemp) goto 30
10      continue

30      if(j.lt.11) then
          goto 40
        else
          if(j.eq.11) then
            ival=-ival
            goto 60
          endif
        endif
        goto 20

40      ival=ival+(j-1)*ifact
        ifact=ifact*10
20      nstrt=nstrt-1
50    continue

60    return
      end

      subroutine ampinp (i,nat,a,na,nb,nc,qend)

c input of AMPAC geometry data

      implicit real*8 (a-h,o-z)
      integer*2 nat,na,nb,nc
      common
     ./digits/ iachar(12), achar(13), end, buffa
      character*8 end, buffa
      character*1 iachar, achar

      logical*1 qend
      parameter (CERO=0.d0)
      dimension a(3), l(3)

      qend = .false.
      nat = 0
      do 5 j=1,3
        a(j) = CERO
        l(j) = 0
5     continue
      na = 0
      nb = 0
      nc = 0

      goto (10,20,30), i
        goto 40

c       FIRST CARD of AMPAC files

10    call readc
      call reada (buffa)
      call atnum (n)
      nat = n
      goto 1000

c     SECOND CARD of AMPAC files

20    call readc
      call reada (buffa)
      if (buffa.eq.end) then
        qend = .true.
        goto 1000
      endif
      call atnum (n)
      nat = n
      call readf (a(1))
      call reada (buffa)
      na = 1
      goto 1000

c     THIRD CARD of AMPAC files

30    call readc
      call reada (buffa)
      if (buffa.eq.end) then
        qend = .true.
        goto 1000
      endif
      call atnum (n)
      nat = n
      do j=1,2
        call readf (a(j))
        call reada (buffa)
      enddo
      na = 2
      nb = 1
      goto 1000

c     other CARD of AMPAC files

40    call readc
      call reada (buffa)
      if (buffa.eq.end) then
        qend = .true.
        goto 1000
      endif
      call atnum (n)
      nat = n
      do j=1,3
        call readf (a(j))
        call reada (buffa)
      enddo
      call readi (n)
      na = n
      call readi (n)
      nb = n
      call readi (n)
      nc = n
      if (nat.eq.0) qend = .true.

1000  continue
      return
      end


      subroutine zinp (i,nat,a,na,nb,nc,qend)

c input of Z-matrix geometry data

      implicit real*8 (a-h,o-z)
      integer*2 nat,na,nb,nc
      common
     ./digits/ iachar(12), achar(13), end, buffa
      character*8 end, buffa
      character*1 iachar, achar

      logical*1 qend
      parameter (CERO=0.d0)
      dimension a(3)

      qend = .false.
      nat = 0
      do 5 j=1,3
        a(j) = CERO
5     continue
      na = 0
      nb = 0
      nc = 0

      goto (10,20,30), i
        goto 40

c       FIRST CARD of Z files

10    call readc
      call reada (buffa)
      call atnum (n)
      nat = n
      goto 1000

c     SECOND CARD of Z files

20    call readc
      call reada (buffa)
      if (buffa.eq.end) then
        qend = .true.
        goto 1000
      endif
      call atnum (n)
      nat = n
      call readi (n)
      na = n
      call readf (a(1))
      goto 1000

c     THIRD CARD of Z files

30    call readc
      call reada (buffa)
      if (buffa.eq.end) then
        qend = .true.
        goto 1000
      endif
      call atnum (n)
      nat = n
      call readi (n)
      na = n
      call readf (a(1))
      call readi (n)
      nb = n
      call readf (a(2))
      goto 1000

c     other CARD of Z files

40    call readc
      call reada (buffa)
      if (buffa.eq.end) then
        qend = .true.
        goto 1000
      endif
      call atnum (n)
      nat = n
      call readi (n)
      na = n
      call readf (a(1))
      call readi (n)
      nb = n
      call readf (a(2))
      call readi (n)
      nc = n
      call readf (a(3))
      if (nat.eq.0) qend = .true.

1000  continue
      return
      end

      subroutine atnum (nat)
      implicit real*8 (a-h,o-z)

*     it takes a string 'buffa' which contains the atomic number or the
*     atomic symbol, indistinctly, and returns the atomic number

      common
     ./digits/ iichar(12), achar(13), end, buffa
     ./elemts/ atnam(37)
      character*1 buffa(8), iichar, achar
      character*4 atnam
      character*8 end

      character*2 carac, elem
      logical*1 qnum

      nat = 0
      ifact = 1
      qnum = .false.
      do i=1,12
        if (buffa(1) .eq. iichar(i)) then
          do j=1,8
            if (buffa(j).eq.char(32)) then
              n = j - 1
              exit
            endif
          enddo
          nstrt=n
          do 50 ii=1,n
            itemp=buffa(nstrt)
            do 10 j=1,12
              if(iichar(j).eq.itemp) goto 30
10          continue
            go to 999

30          if(j.lt.11) goto 40
            if(j.eq.11) nat=-nat
            goto 20
40          nat=nat+(j-1)*ifact
            ifact=ifact*10
            nstrt=nstrt-1
50        continue
          qnum = .true.
          exit
        endif
      enddo

      if (.not.qnum) then
        call pack3 (carac,buffa)
        do j=1,37
          call pack3 (elem,atnam(j))
          call upcase (-1, elem)
          if (carac.eq.elem) then
            nat = j
            return
          endif
        enddo
        goto 999
      endif

20    return

999   write (*,1000) (buffa(i),i=1,n)
      write (9,1000) (buffa(i),i=1,n)
1000  format(/' INPUT ERROR IN THE ATOMIC NUMBER TRANSFER :'/
     &        ' Non allowed character in the read string'/
     &        ' => ', 79a1/)
      call honder(79,'ATNUM ')

      end

      subroutine pack3(buff,ibuf)
      character*4 a,buff
      character*1 ibuf(4),lbuf(4)
      equivalence (a,lbuf)

      do 10 i=1,4
10      lbuf(i)=ibuf(i)
      buff=a

      return
      end

      block data io
      implicit real*8 (a-h,o-z)

      common
     ./answer/ si,bl,hyphen,star,plus,comma,equal
     ./digits/ iachar(12), achar(13), end, buffa
     ./aaaaa / ptype(3)
     ./elemts/ atnam(37)
      character si,bl,hyphen,star,plus,comma,equal
      character iachar, achar
      character*8 end, buffa
      character*4 atnam
      character*10 ptype

      data end /'END'/

      data
     &     ptype /'bond lgth.','bond angle' ,'twist ang.'/

      data  si /'Y'/, bl /' '/, hyphen /'-'/, star /'*'/,
     &      plus /'+'/, comma /','/, equal /'='/

      data
     &      achar /'0','1','2','3','4','5','6','7','8','9','+','-','.'/,
     &      iachar /'0','1','2','3','4','5','6','7','8','9','-','+'/

      data atnam/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     1           'Na','Mg','Al','Si','P ','S ','Cl','Ar',
     2           'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co',
     3           'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','  '/

      end

