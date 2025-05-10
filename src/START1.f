      SUBROUTINE START1 (words, NWORDS, types, NTYPES)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      character*3 ext
      CHARACTER*1 IAAA, HY
      CHARACTER*4 NAM, ext2
      CHARACTER*6 SUBER /'START'/
      CHARACTER*8 TITLE,TYPES(*),SCFTYP,RUNTYP,GUESS,CONF
      CHARACTER*8 FIX,STOP,ALTER,NO,OFF,GROUP,YES,ON,CLOSED,OSCF
      CHARACTER*4 BFNAM,ATNAM
      CHARACTER*8 WORDS(*),ANG,BLANK,DEG,BUFF,COMMNT,hours,minut
      character*9 stat
      character*8 COM /'COM'/,COC /'COC'/, car /'CAR'/
      character*32 name, ifilev
      LOGICAL LBASIS, LCM, dumpf
      LOGICAL REVISE
      LOGICAL LADAPT, LSKIP

      PARAMETER (NSH11=NSH1+1)

c isex(1)  = puntero del segmento T (common SYMTR)
c isex(2)  = puntero del segmento ATDAT (common INFOA)
c isex(3)  = puntero del segmento EX (common NSHEL)
c isex(4)  = puntero del segmento VIBSIZ (common OPTNS)
c isex(5)  = puntero del segmento de la matriz S (overlap matrix)
c isex(6)  = puntero del segmento de la matriz H (core matrix)
c isex(7)  = puntero del segmento de la matriz D (densidad alpha o tot.)
c isex(8)  = puntero del segmento V (vectores propios)
c isex(9)  = puntero del segmento de la matriz E (valores propios)
c isex(10) = puntero del segmento de la matriz D (densidad beta,capa a.)
c isex(11) = puntero del segmento HUCKEL (common HUCKEL)
c isex(12) = puntero del segmento AMASS (common ATMWT)
c isex(13) = puntero del segmento EN (common SCFBLK)
c isex(14) = puntero del segmento de la matriz GRD (gradientes en coord)
c isex(15) = puntero del segmento de la matriz V (coord. actuales OPTX)
c isex(16) = puntero del segmento de la matriz HESS (hessiana)
c isex(17) = puntero del segmento X (common OPTMIZ)
c isex(18) = puntero del segmento de la matriz FCM (cont. fuerza)
c isex(19) = puntero del segmento de la matriz L (MO, capa abierta)
c isex(20) = puntero del segmento de la matriz K (capa abierta)
c isex(21) = puntero del segmento GROUP (common SYMTRB)
c isex(22) = puntero del segmento ANAM (common INFOB)
c isex(23) = puntero del segmento D (densidad opcional)
c isex(24) = puntero del segmento TITLE (common OPTNSB)

c Recordar que si se cambia LDA4 debe cambiarse tambien en otras
c subrutinas.

      parameter (
     &           LDA1=216+(NSH11*24/2)+35,
     &           LDA2=4+NATM1+3*NATM1+4+(7*NATM1)/2+NATM13+NATM1+NATM13,
     &           LDA3=4*NPRM1+(7*NSH1)/2+1,
     &           LDA4=51,
     &           LDA11=nsh1,
     &           LDA12=natm1,
     &           LDA13=23,
     &           LDA21=1,
     &           LDA22=NATM1+NBF1+NATM13+NATM1,
     &           LDA24=14)

      character*80 kwd
      common
     ./optn1/ kwd(48)
     ./optn2/ icar, kb1, iii
     ./origin/ itype, jcar, x0, y0, z0
      common/thermo/cpt,cpr,cpv,cpe,
     &              h25t,h25r,h25v,h25e,
     &              s25t,s25r,s25v,s25e,
     &              g25t,g25r,g25v,g25e,
     &              ilin,nax,
     &              temper
      COMMON/DISC/ISEL,ISELR,ISELW,IREP,ICHEK,IPOS(16)
     1           ,IBLKSZ(16)
      COMMON/DISCB/NAM(16)
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMMAX,TTOT
      COMMON/WORK/JREC,JUMP,ISTT(40),INUMB(40),IWIDTH
      COMMON/WORKB/IAAA(80)
      LOGICAL FIXED,LEX,LDAM12,LDAM13,LDAM23,LDIIS
      COMMON/SCFBLK/EN,ETOT,EHF,SH1(2),SH2(2),GAP1(2),GAP2(2),
     1              D12,D13,D23,CANA,CANB,CANC,FIXED,LEX,
     2              LDAM12,LDAM13,LDAM23,LDIIS,
     3              NCYC,ISCHM,LOCK,MAXIT,NCONV,LOKCYC
      COMMON/HUCKEL/HUCKEL(NSH1),ELEMS(3,18),ELEMP(2,18)
      COMMON/GESBLK/ISWOP(2,25),NSWOP,IFILVA,IBLKA1
      CHARACTER*8 BFLAB,PARNAM,ANAM,ATNAME
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),
     4               KR(NATM1),KA(NATM1),KB(NATM1),
     5               PARVAL(NATM13),ZIN(NATM1),X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2
      COMMON/SYMTRY/T(216),INVT(24),ISO(NSH11,24),NT,
     &              NT2,XXX(21),INDEX,NAXIS
      COMMON/SYMTRB/GROUP
      REAL*4 NOPT, MINULL
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      COMMON/SECTOR/NUM3,IBLK3,REVISE,IODA(24),ICLASS(24),
     1              ILEN(24),MAXB,KBLKLA
      LOGICAL LMOS,LNUCL
      COMMON/PRPSEC/ISECD1,ITYPD1,ISECD2,ITYPD2,LMOS,ISECV,
     1  ITYPV,ISECD3,ITYPD3,LNUCL,GX,GY,GZ
      COMMON/ATMWT/AMS(37),AMASS(natm1)
      common/elemts/ atnam(37)
      common/files/file5
      character*32 file5
      DIMENSION ATDAT(LDA2)
      EQUIVALENCE (ATDAT(1),NAT)
      DIMENSION R(NSH1),BFNAM(10)
      save /optn2/

      DATA BFNAM/'s','x','y','z','xx','yy','zz','xy','xz','yz'/
      DATA BLANK, ONE,  BOHR /'        ', 1.0D0, 0.529177D0/
      DATA OSCF,CLOSED/'OSCF','CLOSED'/,HY/'-'/
      DATA FIX,STOP,ALTER/'FIXED','END','ALTER'/
     1    , YES/'YES'/ ON /'ON'/
     2    , NO /'NO'/, OFF /'OFF'/,COMMNT/'NOTE'/
      data hours/'H'/, minut/'M'/
      DATA ANG /'A'/, DEG /'D'/

      IBLOK=0
      ncoorb=0
      KEY=1

C  LBASIS is true if a previous basis already exists and is to be reused

10    LBASIS=.FALSE.

C  READ TITLE AND OPTIONS

20    if (kb1.eq.1) then
        call readc1
      else
        CALL READC
      endif
      IF (JUMP.EQ.0) GO TO 20
      CALL READA(BUFF)
      DO 30 I=1,NWORDS
      IF (BUFF.EQ.WORDS(I)) GO TO
     1      (60,  90, 200, 210, 220,
     2      750, 240, 270, 280, 710,
     3      720, 290, 310, 320, 330,
     4      610, 370, 380, 410, 420,
     5      440,1050, 620, 750, 960,
     6      770,  45, 700, 650,  40,
     7       40, 800, 830,  20, 160,
     8      780, 100, 340,  50, 910,
     9      110, 120, 130, 930, 940,
     a       41,  42, 400,1460,1470,
     b     1450,1480,1490,  42, 660,
     c      670,  20,  20,  20,  20,
     d     1040,1040,1040,1040,1040,
     e     1040,1040,1040,  20,  20),I
30    CONTINUE
      CALL HONDER(1,SUBER)

C  END, FINISH

40    RUNTYP=WORDS(31)
      RETURN

C  END, STOP

41    RUNTYP=WORDS(46)
      RETURN

C  END, continue or new

42    RUNTYP=WORDS(47)
      RETURN

C  TIME

45    CALL READF(TIMLIM)
      call reada(buff)
      if (buff.eq.hours) then
        timlim = timlim * 3600.d0
      elseif (buff.eq.minut) then
        timlim = timlim * 60.d0
      endif
      TIMMAX=TIMLIM
      GO TO 20

C  WIDTH

50    CALL READI(IWIDTH)
      GOTO 20

C  TITLE

60    if (kb1.ne.1) then
        READ(IR,70) (TITLE(I),I=1,10)
      else
        iii = iii + 1
        read(kwd(iii),70) (TITLE(I),I=1,10)
      endif
70    FORMAT(10A8)

76    write(10,'(//80A1)') ' ',(HY,J=1,79)
      WRITE(IW,'(//80A1)') ' ',(HY,J=1,79)
      write(10,80) (TITLE(J),J=1,10)
      WRITE(IW,80) (TITLE(J),J=1,10)
80    FORMAT (1X,9A8,A7)
      write(10,'(80A1//)') ' ',(HY,J=1,79)
      WRITE(IW,'(80A1//)') ' ',(HY,J=1,79)
      GO TO 20

C  SYMMETRY

90    CALL PTGRP(NPRINT,IW)
      GO TO 20

C     EXTRAPOLATE

100   LEX=.TRUE.
      CALL READA(BUFF)
      IF(BUFF.EQ.NO.OR.BUFF.EQ.OFF) LEX=.FALSE.
      GOTO 20

C     DAMP FACTORS D12,D13,D23

110   CALL READF(D12)
      LDAM12=.TRUE.
      GO TO 20

120   CALL READF(D13)
      LDAM13=.TRUE.
      GOTO 20

130   CALL READF(D23)
      LDAM23=.TRUE.
      GO TO 20

C     CANONICAL

340   CALL READF(CANA)
      CALL READF(CANB)
      CALL READF(CANC)
      GO TO 20

C     NCOORB

160   CALL READI(NCOORB)
      GOTO 20

C     PRINT

200   CALL READI(NPRINT)
      GO TO 20

C     ACCURACY

210   CALL READI(ITOL)
      CALL READI(ICUT)
      CALL READI(ITOLI)
      CALL READI(ICUTI)
      IF(ITOLI.EQ.0)ITOLI=ITOL
      IF(ICUTI.EQ.0)ICUTI=ICUT
      GO TO 20

C     NORM

220   CALL READI(NORMF)
      CALL READI(NORMP)
      GO TO 20

C     --------------------   RESTORE -----------------

240   CONTINUE
      IF (JUMP.GT.1) CALL GETFIL(NNNN,IBLKD,LLLL,IFILD)
241   NUM3=IFILD
      IBLK3=IBLKD
      CALL SECINI(IBLKD,IFILD)

C     Get title and options

      CALL SECGET(ISEX(4),4,IBLOK)
      CALL READQ(VIBSIZ,LDA4,IBLOK,IFILD)
      write(ext,'(i3.3)') nruns
      call strings (3,'.',ext,ext2)
      call oufile (ext2)
      TTOT = TTOT1
      CALL SECGET(ISEX(24),24,IBLOK)
      CALL READQC(TITLE,LDA24,IBLOK,IFILD)
      LT = 0
      DO 245 L=1,10
        LTE = LEN(TITLE(L))
        IF (LTE.NE.0) THEN
          LT = LT + LTE
        ELSE
          GO TO 246
        ENDIF
245     CONTINUE
246   CONTINUE

C  Get symmetry data

      CALL SECGET(ISEX(1),1,IBLOK)
      CALL READQ(T,LDA1,IBLOK,IFILD)
      CALL SECGET(ISEX(21),21,IBLOK)
      CALL READQC(GROUP,LDA21,IBLOK,IFILD)

C  Get atom data

      CALL SECGET(ISEX(2),2,IBLOK)
      CALL READQ(ATDAT,LDA2,IBLOK,IFILD)
      CALL SECGET(ISEX(22),22,IBLOK)
      CALL READQC(ANAM,LDA22,IBLOK,IFILD)
      CALL SECGET(ISEX(11),11,IBLOK)
      nsh = nsh1
      CALL READQ(huckel,nsh,IBLOK,IFILD)
      CALL SECGET(ISEX(12),12,IBLOK)
      CALL READQ(amass,lda12,IBLOK,IFILD)

C  Get basis set data

      CALL SECGET(ISEX(3),3,IBLOK)
      CALL READQ(EX,LDA3,IBLOK,IFILD)

C     get /SCFBLK/

      CALL SECGET(ISEX(13),13,IBLOK)
      CALL READQ(EN,LDA13,IBLOK,IFILD)

      WRITE(IW,250) NAM(IFILD)
      write(10,250) NAM(IFILD)
250   FORMAT (//10X, 'JOB RESTORED FROM DUMPFILE <fn>', A4)
      WRITE(IW,251) (HY,K=1,79),(TITLE(K),K=1,10),(HY,K=1,79),NRUNS
      write(10,251) (HY,K=1,79),(TITLE(K),K=1,10),(HY,K=1,79),NRUNS
251   FORMAT (//10X, 'ORIGINAL TITLE WAS  '
     1        //' ',79A1
     2         /' ',9A8,A7
     3         /' ',79A1
     4        //10X, 'RUN NUMBER  ',I4)
      NRUNS=NRUNS+1
*      DO 260 I=1,NAT
*        IZ=DINT(ZAN(I)+0.5D0)
*        IF(IZ.LE.0.OR.IZ.GT.36)IZ=37
*260     AMASS(I)=AMS(IZ)
      GO TO 20

C     CHARGE

270   CALL READI(ICH)
      GO TO 20

C     MULTIPLICITY

280   CALL READI(NMUL(1))
      CALL READI(NMUL(2))
      IF(NMUL(2).EQ.0) NMUL(2)=1
      MUL=NMUL(1)+NMUL(2)-1
      GO TO 20

C     ANGSTROM

290   CFR=1D0/BOHR
      GO TO 20

C     DENSITY

1460  CALL READI(ISECD1)
      CALL READI(ITYPD1)
      IF(jump.EQ.3)GOTO 20
      CALL READI(ISECD2)
      CALL READI(ITYPD2)
      IF(jump.EQ.5)GOTO 20
      CALL READI(ISECD3)
      CALL READI(ITYPD3)
      GOTO 20

C     MOS

1470  LMOS=.TRUE.
      CALL READI(ISECV)
      CALL READI(ITYPV)
      GOTO 20

C   NONUCLEA

1450  LNUCL=.FALSE.
      GOTO 20

C  MOVE ORIGIN

1480  CALL READA(BUFF)
      IF(BUFF.EQ.COM)THEN
        ITYPE=1
        WRITE(10,101)
        WRITE(9,101)
101     FORMAT(/1X,'Origin to be moved to centre of mass'//)
      ELSE IF(BUFF.EQ.COC)THEN
        ITYPE=2
        WRITE(10,102)
        WRITE(9,102)
102     FORMAT(/1X,'Origin to be moved to centre of charge'//)
      ELSE
        CALL REREAD(-1)
        CALL READF(X0)
        CALL READF(Y0)
        CALL READF(Z0)
        WRITE(10,103) X0,  Y0,  Z0
        WRITE(9,103)  X0,  Y0,  Z0
103     FORMAT(/1X,'Origin to be moved to the point',3F16.8)
        x0 = x0*cfr
        y0 = y0*cfr
        z0 = z0*cfr
        ITYPE=3
      ENDIF
      jcar = 0
      call reada (buff)
      if (buff.eq.car) jcar = 1
      if (itype.eq.3) then
        IF(NAT.EQ.0) then
          WRITE(10,1481)
          WRITE(IW,1481)
          call honder(26,suber)
1481      FORMAT(1H ,'You cannot move the origin before reading the',
     1' geometry.'//
     &' * Move the ORIGIN directive in input file to a line after the AT
     &OMS-END or'/
     &' GEOMETRY-END blocks ... ')
        endif
        CALL MOVORG
      endif
      GO TO 20

C  GAUGE

1490  CALL READF(GX)
      CALL READF(GY)
      CALL READF(GZ)
      GO TO 20

C ---------------FILES-------

C     MAINFILE

310   CALL GETFIL(NFILE,IBLKM,MBLKM,IFILM)
      GO TO 20

C     DUMPFILE

320   CALL GETFIL(NNNN,IBLKD,LLLL,IFILD)
      NUM3=IFILD
      IBLK3=IBLKD
      GOTO 20

330   CALL GETFIL(NNNN,IBLKS,LLLL,IFILS)
      GO TO 20

C     GEOMETRY

370   LBASIS=.TRUE.

C     ATOMS

380   CALL ATOMS(LBASIS)
      LBASIS=.FALSE.
      CALL STRUCT (NAT, C, ZAN)
      DO 390 I=1,NAT
      IZ=DINT(ZAN(I)+0.1)
      IF (IZ.LE.0.OR.IZ.GT.36) IZ=37
390   if (amass(i).eq.0.0d0) AMASS(I)=AMS(IZ)
      GO TO 20

C     INITIAL GUESS BY EHT

400   guess = 'HUCKEL'
401   call readf (aminu)
      minull = aminu
      if (minull.eq.0.0) minull = 0.84
      GO TO 20

C     WEIGHTS

610   if (kb1.eq.1) then
        call readc1
      else
        CALL READC
      endif
      CALL READA(BUFF)
      IF(BUFF.EQ.STOP) GO TO 20
      CALL READI(IZ)
      CALL READF(Z)
      anam(iz) = buff
      AMASS(IZ)=Z
      GO TO 610

C     MAXIT

410   CALL READI(MAXIT)
      GO TO 20

C     CONVERGENCE

420   CALL READI(NCONV)
      GO TO 20

C     GUESS

440   CALL READA(GUESS)
      IF(GUESS.EQ.WORDS(7)) GOTO 450
      IF(GUESS.EQ.WORDS(48)) GOTO 401
      GOTO 20

C   -------------- GUESS --------------

C     GUESS RESTORE

450   IF (JUMP .GT. 2) GO TO 460
      IFILVA=IFILD
      IBLKA1=0
      GOTO 20

460   CALL READA4(IFILEV)
      if (len(ifilev).eq.3) then
        do 470 i=1,16
          if(ifilev.eq.nam(i)) then
            ifilva=i
            goto 471
          endif
470     continue
471     name = file5
        call filen1 (1,name,nam(ifilva))
      else
        name = ifilev
      endif
      inquire (file=name, exist=dumpf)
      if (dumpf) then
        stat='AVAILABLE'
      else
        stat='CREATED'
      endif
      write(iw,150) stat,name
      write(10,150) stat,name
150   format (1x, a9,': ',a32)
      close (ifild)
      open(ifild,file=name,status='unknown',access='direct',
     &        form='unformatted',recl=4096)
      ipos(ifild)=0
      call readi(iblka1)
      goto 20

C     SCFTYP

710   CALL READA(SCFTYP)
      GO TO 20

C     OSCF

720   SCFTYP=BUFF
      GO TO 20

C     RUNTYP

620   CALL READA(RUNTYP)
********************************
*      IF (RUNTYP.EQ.'FORCE') THEN
*        WRITE (*,6200)
*        WRITE (9,6200)
*6200  FORMAT(/' *** FATAL ERROR ***'/
*     &' The FORCE option is not available in this MICROMOL version ...'
*     &)
*        CALL HONDER(15,SUBER)
*      ENDIF
********************************
      DO 630 KEY=1,20
      IF (WORDS((NWORDS-20)+KEY) .EQ. RUNTYP) GO TO 640
630   CONTINUE
      CALL HONDER(15,SUBER)
640   CONTINUE
      GOTO 20

C     IMPTPRINT (unknown purpose)

650   CALL READI(I)
      IF (I .EQ. 0) GO TO 20
      IF (I .LT. 1 .OR. I .GT. 15) CALL HONDER(8,SUBER)
      GO TO 650

c     symmetry axis for rotational calculations

660   call readi(nax)
      goto 20

C     reading temperature for thermodynamic properties

670   call readf(temper)
      goto 20

C     SKIP

700   LSKIP=.TRUE.
      GO TO 20

C     SHIFT

750   SHB=0.0D0
      SHD=0.0D0
      CALL READF(SHA)
      IF(JUMP.GE.6) CALL READF(SHB)
      CALL READI(NCYC)
      CALL READF(SHC)
      IF(JUMP.GE.6) CALL READF(SHD)
      CALL READA(BUFF)
      IF(BUFF.EQ.FIX)GOTO 760
      FIXED=.FALSE.
      GAP1(1)=SHA
      GAP1(2)=SHB
      GAP2(1)=SHC
      GAP2(2)=SHD
      SH1(1)=0.5D0
      SH1(2)=0.5D0
      GOTO 20
760   FIXED=.TRUE.
      SH1(1)=SHA
      SH1(2)=SHB
      SH2(1)=SHC
      SH2(2)=SHD
      GO TO 20

C     LOCK

770   LOCK=1
      CALL READI(LOKCYC)
      GO TO 20

C     SWOP

780   NSWOP=0
      CONF=BLANK
790   if (kb1.eq.1) then
        call readc1
      else
        CALL READC
      endif
      CALL READA(BUFF)
      IF (BUFF.EQ.STOP) GO TO 20
      CONF=ALTER
      JREC=0
      NSWOP=NSWOP+1
      CALL READI(ISWOP(1,NSWOP))
      CALL READI(ISWOP(2,NSWOP))
      GO TO 790

C     VARIABLES

800   IF (NVAR .EQ. 0) GO TO 820
      I=0
810   I=I+1
      IF (I+NVAR .GT. NPAR) GO TO 820
      PARNAM(I)=PARNAM(I+NVAR)
      PARVAL(I)=PARVAL(I+NVAR)
      GO TO 810
820   if (kb1.eq.1) then
        call readc1
      else
        CALL READC
      endif
      CALL READA(BUFF)
      IF (BUFF .EQ. COMMNT .OR. BUFF .EQ. BLANK) GO TO 820
      IF (BUFF .EQ. STOP) GO TO 20
      NVAR=NVAR+1
      NPAR=NPAR+1
      NVAL=NPAR
      PARNAM(NPAR)=PARNAM(NVAR)
      PARVAL(NPAR)=PARVAL(NVAR)
      PARNAM(NVAR)=BUFF
      CALL READF(PARVAL(NVAR))
      CALL READA(BUFF)
      IF (BUFF .EQ. ANG) PARVAL(NVAR)=PARVAL(NVAR)/BOHR
      IF (BUFF .EQ. BLANK) PARVAL(NVAR)=PARVAL(NVAR)*CFR
      IF (BUFF .EQ. DEG) PARVAL(NVAR)=PARVAL(NVAR)*CFA
      GO TO 820

C     CONSTANTS

830   NPAR=NVAR
840   if (kb1.eq.1) then
        call readc1
      else
        CALL READC
      endif
      CALL READA(BUFF)
      IF (BUFF .EQ. COMMNT .OR. BUFF .EQ. BLANK) GO TO 840
      IF (BUFF .EQ. STOP) GO TO 20
      NPAR=NPAR+1
      NVAL=NPAR
      PARNAM(NPAR)=BUFF
      CALL READF(PARVAL(NPAR))
      CALL READA(BUFF)
      IF (BUFF .EQ. ANG) PARVAL(NPAR)=PARVAL(NPAR)/BOHR
      IF (BUFF .EQ. BLANK) PARVAL(NPAR)=PARVAL(NPAR)*CFR
      IF (BUFF .EQ. DEG) PARVAL(NPAR)=PARVAL(NPAR)*CFA
      GO TO 840

C     GRADTOL

910   CALL READF(exgrad)
      nopt = exgrad
      IF(NOPT.EQ.0) NOPT=3.0
      GO TO 20

C     STEPSIZE FOR FORCE CONSTANT CALCULATIONS

930   CALL READF(VIBSIZ)
      GO TO 20

C     NUMBER OF POINTS FOR DIFFERENTIATION
C     IN NUMERICAL FORCE CONSTANT CALCULATIONS

940   CALL READI(NVIB)
      GO TO 20

C     DIIS

960   CALL READA(BUFF)
      IF(BUFF.EQ.YES.OR.BUFF.EQ.ON) LDIIS=.TRUE.
      IF(BUFF.EQ.OFF.OR.BUFF.EQ.NO) LDIIS=.FALSE.
      CALL READF(DIISE)
      IF(DIISE.EQ.0.0) DIISE=0.05
      CALL READF(DIISD)
      IF(DIISD.EQ.0.0) DIISD=1.0D-08
      CALL READI(MINVEC)
      IF(MINVEC.EQ.0) MINVEC=4
      GO TO 20

C     HFSCF

1040  RUNTYP=BUFF
      KEY=I-(NWORDS-NTYPES)
      GOTO 20

C    ---------  START      START      START   ------------

1050  IF ((RUNTYP.EQ.TYPES(7) .OR. RUNTYP.EQ.TYPES(8))
     &                       .AND. NCONV.LE.6) THEN
        WRITE (IW,1051)
        WRITE (10,1051)
1051  FORMAT (/' CONVERGENCE EXPONENT IS BEING SET TO 7 TO RUN "OPTIMIZE
     &" OR "FORCE" OPTIONS')
         NCONV=7
      ENDIF
      IF(JUMP.EQ.2) CALL READI(IREST)

C     EXIT FOR SCF

      IF(KEY.EQ.4) GOTO 1350
      WRITE(IW,1060) NAM(IFILM),IBLKM,NAM(IFILD),
     1    IBLKD,NAM(IFILS),IBLKS
1060  FORMAT(///
     1    /1X,'MAIN FILE ON  <fn>',A4,' STARTING AT BLOCK',I4
     2    /1X,'DUMP FILE ON  <fn>',A4,' STARTING AT BLOCK',I4
     3    /1X,'SCRATCH FILE ON   ',A4,' STARTING AT BLOCK',I4//
     4    //' START')
      CALL CHKFIL(IFILD)
      CALL CHKFIL(IFILS)
      WRITE(IW,1070) TIMLIM,NPRINT,ITOL,ICUT,NORMF,NORMP,
     1   IREST,IST,JST,KST,LST,NREC,INTLOC
1070  FORMAT(/' TIME LIMIT = ',F10.3 /
     1    ' OUTPUT OPTION = ',I5/
     2    ' PREFACTOR TOLERANCE FOR INTEGRALS = 1.0D-',I2/
     3    ' INTEGRAL CUTOFF = 1.0D-',I2/
     4    ' OPTION FOR NORMALIZATION OF BASIS FUNCTIONS = ',I5/
     5    ' OPTION FOR UNNORMALIZATION OF PRIMITIVES = ',I5/
     7    ' RESTART OPTION = ',I5//
     8    ' STARTING SHELLS FOR 2E-INTEGRALS = ',4I5/
     9    ' RECORD NUMBER = ',I10/
     A    ' LOCATION OF FIRST INTEGRAL IN THIS RECORD = ',I5)
      IF (IREST.GE.0) GO TO 1190
      N=0
      DO 1090 II=1,NSHELL
        IAT=KATOM(II)
        J=DINT(ZAN(IAT)+0.01D0)
        IF(J.LE.0.OR.J.GE.37)J=37
        MINI=KMIN(II)
        MAXI=KMAX(II)
        DO 1090 I=MINI,MAXI
          N=N+1
          CALL PACK2(BFLAB(N),ATNAM(J),BFNAM(I))
1090      CONTINUE
      IF(NPRINT.EQ.-1)GOTO 1190
      IPASS=1
      FAC=ONE
      WRITE(IW,1100)
1100  FORMAT(//10X,'INTERNUCLEAR DISTANCES ( A.U.)'/10X,30(1H-))
1110  MAX=0
1120  MIN=MAX+1
      MAX=MAX+5
      IF (MAX.GT.NAT) MAX=NAT
      WRITE(IW,1130)
1130  FORMAT (1H )
      WRITE(IW,1140) (ANAM(J),J=MIN,MAX)
1140  FORMAT(17X,5(3X,A8,4X))
      WRITE(IW,1130)
      DO 1170 I=1,NAT
        DO 1150 J=MIN,MAX
1150      R(J)=FAC*RIJ(C,I,J)
        WRITE(IW,1160) I,ANAM(I),(R(J),J=MIN,MAX)
1160    FORMAT(I3,2X,A8,2X,5F13.6)
1170    CONTINUE
      IF (MAX.LT.NAT) GO TO 1120
      IF (IPASS.EQ.2) GO TO 1190
      IPASS=2
      FAC=BOHR
      WRITE(IW,1180)
1180  FORMAT(/10X,'INTERNUCLEAR DISTANCES (ANGS.)'/10X,30(1H-))
      GO TO 1110

1190  CONTINUE
1200  FORMAT(/' TOTAL NUMBER OF SHELLS                      =',I5/
     1        ' TOTAL NUMBER OF BASIS FUNCTIONS             =',I5/
     2        ' NUMBER OF ELECTRONS                         =',I5/
     3        ' CHARGE OF MOLECULE                          =',I5/
     4        ' STATE MULTIPLICITY                          =',I5/
     5        ' NUMBER OF DOUBLY OCCUPIED ORBITALS          =',I5/
     6        ' NO OF SINGLY OCCUPIED ORBITALS (ALPHA)      =',I5/
     7        ' NO OF SINGLY OCCUPIED ORBITALS (BETA)       =',I5/
     8        ' TOTAL NUMBER OF ATOMS                       =',I5)
      IF(SCFTYP.EQ.CLOSED.AND.MUL.NE.1) WRITE(IW,1291)
1291  FORMAT(1H ,'MULTIPLICITY NOT ONE WITH SCF CLOSED')
      NE=-ICH
      DO 1210 I=1,NAT
1210    NE=NE+DINT(ZAN(I)+0.1D0)
      LCM=MOD(NE+MUL,2).EQ. 1
      NA=(NE+MUL-1)/2
      NB=(NE-MUL+1)/2
      IF(SCFTYP.EQ.OSCF) then
        NA=(NE-MUL+1)/2
        NB=NA+(MUL-1)
      endif
      NDOUB=NA
      NSINA=0
      IF(SCFTYP.EQ.OSCF) NSINA=NB-NA
      NSINB=0
      WRITE(IW,1200) NSHELL,NUM,NE,ICH,MUL,NDOUB,NSINA,NSINB,NAT
      write(10,1200) NSHELL,NUM,NE,ICH,MUL,NDOUB,NSINA,NSINB,NAT
1280  IF (LCM) GO TO 1300
      WRITE (IW,1290)
      WRITE (*,1290)
1290  FORMAT (/' Multiplicity inconsistent with number of electrons'
     1    /' JOB ABANDONED'/)
      call honder(18,SUBER)

1300  WRITE(IW,1310) SCFTYP,RUNTYP,GUESS
      write(10,1310) SCFTYP,RUNTYP,GUESS
1310  FORMAT(' SCFTYP = ',A8   / ' RUNTYP = ',A8     /'  GUESS = ',A8)

C    CREATE DIRECTORY AND WRITE MAIN COMMONS IN DUMPFILE "<fn>.ED3"

1340  IF(IREST.EQ.-1)IREST=0

1350  LDA(ISEX(4))=LDA4
      LDA(ISEX(24))=LDA24
      LDA(ISEX(1))=LDA1
      LDA(ISEX(21))=LDA21
      LDA(ISEX(2))=LDA2
      LDA(ISEX(22))=LDA22
      LDA(ISEX(3))=LDA3
      LDA(ISEX(13))=LDA13
      LDA(ISEX(11))=LDA11
      lda(isex(12))=lda12
      CALL SECPUT(ISEX(4),4,LENSEC(LDA4),IBLOK)
      CALL WRT3(VIBSIZ,LDA4,IBLOK,IFILD)
      CALL SECPUT(ISEX(24),24,NRC(LDA24),IBLOK)
      CALL WRT3C(TITLE,LDA24,IBLOK,IFILD)
      CALL SECPUT(ISEX(1),1,LENSEC(LDA1),IBLOK)
      CALL WRT3(T,LDA1,IBLOK,IFILD)
      CALL SECPUT(ISEX(21),21,NRC(LDA21),IBLOK)
      CALL WRT3C(GROUP,LDA21,IBLOK,IFILD)
      CALL SECPUT(ISEX(2),2,LENSEC(LDA2),IBLOK)
      CALL WRT3(ATDAT,LDA2,IBLOK,IFILD)
      CALL SECPUT(ISEX(22),22,NRC(LDA22),IBLOK)
      CALL WRT3C(ANAM,LDA22,IBLOK,IFILD)
      CALL SECPUT(ISEX(3),3,LENSEC(LDA3),IBLOK)
      CALL WRT3(EX,LDA3,IBLOK,IFILD)
      CALL SECPUT(ISEX(13),13,LENSEC(LDA13),IBLOK)
      CALL WRT3(EN,LDA13,IBLOK,IFILD)
      CALL SECPUT(ISEX(11),11,LENSEC(LDA11),IBLOK)
      nsh = nsh1
      CALL WRT3(huckel,nsh,IBLOK,IFILD)
      call secput(isex(12),12,lensec(lda12),iblok)
      call wrt3(amass,lda12,iblok,ifild)
      CALL TIMIT(3)
      CALL REVIND

      RETURN
      END

      subroutine strings (n,carac,buffi,buffo)
      character*1 buffi(*), buffo(32), carac
      buffo(1) = carac
      do i=1,n
        buffo(i+1) = buffi(i)
      enddo
      return
      end

