C  OPTX  A1

      SUBROUTINE OPTX

C  Minimization of the energy by varying the geometry
C  of the molecule. (Adaptation of Andy Komornicki's
C  routine, using Murtagh and Sargent algorithm)

C  This is an updating routine which successively
C  generates the Hessian matrix (i.e. the inverse of the
C  matrix of second derivatives) after a successful one
C  dimensional search.

c  ICODE is a key pointer
c  ICODE = 1 Continue search cycles (update Hessian Matrix)
c  ICODE = 2 Start search cycles (initialize Hessian matrix)
c  ICODE = 3 End of search (being the geometry optimized or not)

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 GRAD1, OPTIMZ, GRAD, RESTR
      LOGICAL RSTART,RESET, INTERN
      character*22 statopt(3) /
     &'CONTINUE SEARCH CYCLES',
     &'START SEARCH CYCLES   ',
     &'END OF SEARCH CYCLES  '/

      COMMON/BIG/HESS(NATM13,NATM13),QJNK(2800)

      LOGICAL LADAPT,LSKIP
      INTEGER P1, P2, P3
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMMAX,TTOT
      CHARACTER*8 ATNAME,ANAM,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(NATM13),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      REAL*4 NOPT
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/OPTMIZ/V(NATM13),G(NATM13),DV(NATM13),FUNC,ALPHA,
     *    V0(NATM13),G0(NATM13),GS0,ALPH
     *    ,FUNC0,NCOORD,NPTS,NSERCH,IUPDAT,ICODE,ISS,BIGG,BIGG0
     *    ,Z(NATM13),Q(NATM13),CC,GNRM

      CHARACTER*32 FILE5
      COMMON
     ./FILES / FILE5
      COMMON/FUNCT/E,EG(NATM13)
      LOGICAL CARF

      COMMON/GESBLK/ISWOP(2,25),NSWOP,IFILVA,IBLKA1

      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TEN=1.0D1)
      DATA GRAD1/'GRADONE'/
      DATA OPTIMZ,GRAD,RESTR /'OPTIMIZE', 'GRADIENT', 'RESTORE'/

      RSTART=.FALSE.
      RESET=.FALSE.
      ISEC15=0
      ISEC16=0
      ISEC4=0
      ISC24=0
      TOLG=TEN**(-NOPT)
      IF(RUNTYP.EQ.OPTIMZ) then
        WRITE(IW,190) CHAR(12),TOLG
        write(10,190) ' ',TOLG
      endif
      INTERN = NVAR.NE.0

C  Check for a restart run

      IF (.NOT.(RUNTYP.EQ.OPTIMZ).OR.IREST2.LE.0) GOTO 10
        CALL SECGET(ISEX(15),15,ISEC15)
        CALL READQ(V,LDA(ISEX(15)),ISEC15,IFILD)

        IF (ICODE.EQ.3) GOTO 120

        RSTART=NSERCH.GE.0
        IF(.NOT.RSTART)GOTO 10
        GUESS=RESTR
        IFILVA=IFILD
        IBLKA1=0
        IF (RSTART .AND. IUPDAT.EQ.0) GO TO 80
        IF (RSTART) GO TO 70


C  Fresh start, evaluate first gradient

10    NSERCH=-1
      BIGG=ZERO
      BIGG0=ZERO
      NPTS=-1
      if (intern) then
        ncoord = nvar
      else
        NCOORD=3*NAT
      endif
      E=ZERO
      FUNC=ZERO
      FUNC0=ZERO
      DO 20 I=1,NCOORD
        EG(I)=ZERO
        DV(I)=ZERO
        V0(I)=ZERO
        G0(I)=ZERO
        if (intern) then
          v(i) = parval(i)
        else
          V(I)=C(I)
        endif
20      G(I)=ZERO

      LDA(ISEX(15))=3*NATM13+2+2*NATM13+8
      CALL SECPUT(ISEX(15),15,LENSEC(LDA(ISEX(15))),ISEC15)
      CALL WRT3(V,LDA(ISEX(15)),ISEC15,IFILD)
      CALL TIMIT(3)
      CALL REVIND

      IREST2=1
      write (*,'(/a)') ' > OPTX is calling VALOPT'
      write (*,'(a)') ' (Searching gradients in the current geometry)'
      CALL VALOPT
      call timit(4)

      IF (INTERN) then
        WRITE (10,1001) (PARNAM(I), PARVAL(I), G(I), I=1,NVAR)
        WRITE (IW,1001) (PARNAM(I), PARVAL(I), G(I), I=1,NVAR)
      endif
      GNRM=DSQRT(VECSUM(G,G,NCOORD))
*      IF(RUNTYP.EQ.OPTIMZ) NPRINT=-5
      IFILVA=IFILD
      IBLKA1=0
      GUESS=RESTR
      DO 30 I=1,NCOORD
        V0(I)=V(I)
30      G0(I)=G(I)
      ALPHA=ONE
      ALPH=ALPHA
      IUPDAT=0
      ICODE=2
      NSERCH=0
      IREST=0
      DO 40 I=1,NCOORD
40      BIGG=DMAX1(DABS(G(I)),BIGG)

      CALL WRT3(V,LDA(ISEX(15)),ISEC15,IFILD)
      CALL SECGET(ISEX(4),4,ISEC4)
      CALL WRT3(VIBSIZ,LDA(ISEX(4)),ISEC4,IFILD)
      CALL SECGET(ISEX(24),24,ISC24)
      CALL WRT3C(TITLE,LDA(ISEX(24)),ISC24,IFILD)
      
      IF(BIGG.LE.TOLG) ICODE=3
      IF(RUNTYP.EQ.GRAD1) ICODE=3
      IF(RUNTYP.EQ.GRAD) ICODE=3
      GO TO 80

C  Start search

50    CONTINUE
*      NPRINT=-5
      ICODE=1
      IUPDAT=1

C  Initialise Hessian matrix

      LENG=LENSEC(NAT*NAT*9)
      CALL SECPUT(ISEX(16),16,LENG,ISEC16)
      LDA(ISEX(16))=NAT*NAT*9
      CALL SETS0(HESS,NCOORD,ISEC16,IFILD)
      CALL TIMIT(3)
      CALL REVIND
      GO TO 65

C  Update Hessian matrix

60    IUPDAT=IUPDAT+1
      CALL SECGET(ISEX(16),16,ISEC16)
      CALL SETS1(HESS,Z,CC,NCOORD,ISEC16,IFILD)

65    ALPHA=ONE
      CALL SECGET(ISEX(16),16,ISEC16)
      CALL MULTRC(HESS,G,DV,NCOORD,ISEC16,IFILD)

C  Symmetrize displacement vector

      IF (.NOT. INTERN) CALL SYMDR(DV)

C  Search for a minimum in the DV direction

70    write (*,'(/a)') ' > OPTX is calling SMSL'
      write (*,'(a)') ' (Searching a minimum in the direction of the dis
     &placement vector)'
      CALL SMSL(TOLG,RSTART,RESET)
      call timit(4)
      IF(NSERCH.LT.50) GO TO 80
      CALL WRT3(VIBSIZ,LDA(ISEX(4)),ISEC4,IFILD)
      CALL WRT3C(TITLE,LDA(ISEX(24)),ISC24,IFILD)
      write(10,130)
      WRITE(IW,130)
      ICODE=3

C  Print updated data

80    CONTINUE
      DFUNC=FUNC-FUNC0
      IF(NPTS.EQ.0) GOTO 100
      WRITE(IW,140) NSERCH,IUPDAT,NPTS,FUNC,GNRM,ALPHA,statopt(ICODE),
     *              FUNC0,DFUNC
      write(10,140) NSERCH,IUPDAT,NPTS,FUNC,GNRM,ALPHA,statopt(ICODE),
     1              FUNC0,DFUNC
140   FORMAT(//5X,'Search  Update  Points     Energy        ',
     *     'Gradient norm      alpha'
     *    /1X, 3I8, F17.8, F15.8, F14.5
     *    //10x,'Status: ',a22
     *    //10X,'Previous energy',F17.8
     *    //10X,'    Convergence', F17.8)
      IF (.not.INTERN) then
        WRITE(IW,150)
        WRITE(10,150)
        DO 90 IAT=1,NAT
          N=3*(IAT-1)
          write(10,160) IAT,ANAM(IAT),ZAN(IAT),
     *                  (V(N+I),I=1,3),(G(N+I),I=1,3)
90        WRITE(IW,160) IAT,ANAM(IAT),ZAN(IAT),
     *                  (V(N+I),I=1,3),(G(N+I),I=1,3)
      else
        WRITE (IW,1001) (PARNAM(I), PARVAL(I), G(I), I=1,NVAR)
        WRITE (10,1001) (PARNAM(I), PARVAL(I), G(I),I=1,NVAR)
1001  FORMAT (1H ,'   Parameter         value          gradient'
     *       /'       name       (bohr or radian)'
     *       /(1H ,6X, A8, 2F16.8))
      endif


      WRITE(IW,180) BIGG,TOLG,BIGG0
      write(10,180) BIGG,TOLG,BIGG0
100   BIGG0=BIGG
      NSERCH=NSERCH+1
      GO TO (60,50,120), ICODE

C  End of run

120   IF (NPTS.GT.0) then
        WRITE (IW,170) NPTS
        write(10,170) NPTS
      endif

C  Printing wavefunction properties of optimized geometry

      IF (RUNTYP.EQ.OPTIMZ) then
        CALL OPTEND

        CALL FILEN1 (1,FILE5,'.CAR')
        INQUIRE (FILE=FILE5, EXIST=CARF)
        OPEN (20,FILE=FILE5,STATUS='UNKNOWN')
        IF (CARF) CALL APPE (20)
        WRITE (20,200) NAT
        WRITE (20,210) TITLE
        WRITE (20,210) 'NOTE ab ','initio o','ptimized',' geometr',
     &                 'y       '
        DO 901 IAT=1,NAT
          N=3*(IAT-1)
901       WRITE(20,220) (V(N+I)*0.529177D0,I=1,3),INT(ZAN(IAT))
        CLOSE (20)
      ENDIF

130   FORMAT(/10X,' ***** MINIMIZATION NOT CONVERGING ***** ')
150   FORMAT(//22X,'Coordinates (Bohr)',14X,'Gradient (Hartree/Bohr)',
     *    /3X,'Atom       Z       x',9X,'y',9X,'z',
     *    13X,'x',9X,'y',9X,'z'/)
160   FORMAT(I2, 1X, A8, 1X, F4.1, 3F10.6, 4X, 3F10.6)
170   FORMAT(//' Number of function evaluations =',I3/)
180   FORMAT(/' Largest component of gradient =',F14.7,
     *     ' ( tolerance = ',F14.7,' )'
     *    /' Previous largest component of gradient =',F14.7)
190   FORMAT(A//6X,14(5H*****)/31X,'GEOMETRY OPTIMIZATION',
     $   /6X,14(5H*****)//6X,'Convergence threshold on gradient = ',
     *     F10.8/)
200   FORMAT (I4)
210   FORMAT (1X,9A8,A7)
220   FORMAT (1X,3F10.5,I4)
      RETURN
      END

      SUBROUTINE APPE (IFILE)
      IMPLICIT INTEGER*2 (I-N)
      INTEGER*4 IFILE
      CHARACTER*1 A
      CHARACTER*11 ANS

      INQUIRE (IFILE,FORM=ANS)
      IF (ANS.EQ.'UNFORMATTED') THEN
        WRITE(10,'(/A/A/)')' ERROR: File is SEQUENTIAL and UNFORMATTED',
     &                     ' It can not be appended by APPE...'
        CALL HONDER (74,'OPTX')
      ENDIF
1     IF (EOF(IFILE)) THEN
        RETURN
      ENDIF
      READ (IFILE,'(A)',END=2) A
      GO TO 1
      RETURN
2     BACKSPACE IFILE
      RETURN
      END

C  SMSL  A1
C   9 NOV 1983
C     COMPILE WITH AUTODBL
      SUBROUTINE SMSL(TOLG,RSTART,RESET)

C  Set up arrays to update Hessian matrix
C  by performing a one dimensional search.

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL RSTART,RESET, INTERN

      LOGICAL LADAPT,LSKIP
      INTEGER P1, P2, P3
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMMAX,TTOT
      CHARACTER*8 PARNAM,BFLAB,ATNAME,ANAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(NATM13),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      DIMENSION ATDAT(1)
      EQUIVALENCE (ATDAT(1),NAT)

      COMMON/OPTMIZ/V(NATM13),G(NATM13),DV(NATM13),FUNC,ALPHA,
     *    V0(NATM13),G0(NATM13),GS0,ALPH,FUNC0,NCOORD,NPTS,NSERCH,
     *    IUPDAT,ICODE,ISS,BIGG,BIGG0,Z(NATM13),Q(NATM13),CC,GNRM

      COMMON/BIG/HESS(NATM13,NATM13),QJNK(2800)

      PARAMETER (ZERO=0.0D0, TWO=2.0D0,
     &           ALPHL=1.0D-02, TEPS=1.0D-05,
     &           DELTP=1.0D-08, TOLC=1.0D-02)


      ISEC2=0
      ISC22=0
      ISEC15=0
      ISEC16=0
      INTERN=(NVAR.NE.0)
      EPS=DELTP
      DELTPM=-DELTP

C  Get restart data if necessary

      CALL SECGET(ISEX(2),2,ISEC2)
      CALL SECGET(ISEX(22),22,ISC22)
      CALL SECGET(ISEX(15),15,ISEC15)
      CALL SECGET(ISEX(16),16,ISEC16)
      IF(.NOT.RSTART) GO TO 20
      CALL READQ(V,LDA(ISEX(15)),ISEC15,IFILD)
      GO TO 70

C  Normal start

20    CONTINUE
      IF(RESET) GO TO 40
      DO 30 I=1,NCOORD
        V0(I)=V(I)
30      G0(I)=G(I)
      FUNC0=FUNC
40    CONTINUE
      RESET=.FALSE.
      GS0=VECSUM(G0,DV,NCOORD)
      ALPH=ALPHA

C  Search along DV direction.

50    DO 60 I=1,NCOORD
60      V(I)=V0(I)-ALPH*DV(I)

C  Save coordinates + displacement vector

      CALL WRT3(V,LDA(ISEX(15)),ISEC15,IFILD)

      IF (INTERN) GO TO 65
      DO 62 I=1,NCOORD
62      C(I)=V(I)
      GO TO 69
65    DO 66 I=1,NVAR
66      PARVAL(I)=V(I)
      CALL POSN(.FALSE.)
      CALL SYMGEN(.FALSE.)
69    CALL WRT3 (ATDAT,LDA(ISEX(2)),ISEC2,IFILD)
      CALL WRT3C(ANAM,LDA(ISEX(22)),ISC22,IFILD)

C  Call function evaluation

70    continue
      write (*,'(/a)') ' .. SMSL is calling VALOPT'
      write (*,'(a)') ' (Searching gradients)'
      CALL VALOPT
      call timit(4)

C  This is no longer a restart job.   Reset RSTART and IREST.

      RSTART=.FALSE.
      IREST=0
      IF(FUNC.GE.FUNC0) GO TO 80
      GNRM=DSQRT(VECSUM(G,G,NCOORD))
      TEST=EPS*GS0*ALPH-TEPS
      DIFF=FUNC0-FUNC
      IF(DIFF.GE.TEST) GO TO 110

C  Quadratic fit to get new step size

80    ALPH=GS0*ALPH**2/(TWO*(GS0*ALPH+FUNC-FUNC0))
      IF(ALPH.GT.ALPHL) GO TO 50
      WRITE(IW,90)
      write(10,90)
90    FORMAT(12X,20(1H*),' ALPHA TOO SMALL ',20(1H*)/
     1 25X,'Reset Hessian matrix to unity'/)
      ICODE=2
      RESET=.TRUE.
      ALPHA=ALPH
      DO 100 I=1,NCOORD
        V(I)=V0(I)
100     G(I)=G0(I)
      RETURN

C  Successful search

110   BIGG=ZERO
      DO 120 I=1,NCOORD
120     BIGG=DMAX1(DABS(G(I)),BIGG)
      IF(BIGG.LE.TOLG) GO TO 150
      ALPHA=ALPH
      DO 130 I=1,NCOORD
130     Q(I)=G(I)-G0(I)
      CALL MULTRC(HESS,Q,Z,NCOORD,ISEC16,IFILD)
      DO 140 I=1,NCOORD
140     Z(I)=-(Z(I)+ALPH*DV(I))
      CC=VECSUM(Q,Z,NCOORD)
      ZTG=VECSUM(Z,G0,NCOORD)
      ZTG=ZTG/CC
      Z2=VECSUM(Z,Z,NCOORD)
      DZ=TOLC*Z2
      IF(DABS(CC).GE.DZ.AND.ZTG.LE.DELTPM) RETURN
      ICODE=2
      RETURN

150   ICODE=3
      ALPHA=ALPH
      WRITE(IW,160)
      write(10,160)
160   FORMAT(///T22,'***** OPTIMIZATION CONVERGED *****'/)
      RETURN
      END

C  VALOPT  A1
      SUBROUTINE VALOPT

C  Calculate HF energy + gradient with respect to nuclear coordinates.

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 OPT
      LOGICAL LADAPT,LSKIP
      INTEGER P1, P2, P3
      CHARACTER*8 ANAM,ATNAME,PARNAM,BFLAB
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13),  ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      LOGICAL FIXED,LEX,LDAM12,LDAM13,LDAM23,LDIIS
      COMMON/SCFBLK/EN,ETOT,EHF,SH1(2),SH2(2),GAP1(2),GAP2(2),
     1              D12,D13,D23,CANA,CANB,CANC,FIXED,LEX,
     2              LDAM12,LDAM13,LDAM23,LDIIS,
     3              NCYC,ISCHM,LOCK,MAXIT,NCONV,LOKCYC

      COMMON/FUNCT/E,EG(NATM13)

      COMMON/OPTMIZ/V(NATM13),G(NATM13),DV(NATM13),
     *    FUNC,ALPHA,V0(NATM13),G0(NATM13),GS0,
     *    ALPH,FUNC0,NCOORD,NPTS,NSERCH,IUPDAT,ICODE
     $    ,ISS,BIGG,BIGG0,Z(NATM13),Q(NATM13),CC,GNRM

      PARAMETER (ZERO=0.0D0, TOLG=1.0D-20)
      DATA OPT/'OPTIMIZE'/

      IBLOK=0
      NC=3*NAT
      NPTS=NPTS+1

      IF(RUNTYP.NE.OPT) GOTO 16
      IF (NSERCH.GT.0) then
        WRITE(IW,9999) NSERCH
        write(10,9999) NSERCH
9999    FORMAT(' Line search',I3)
      endif
      WRITE (IW,9997) NPTS
      WRITE (10,9997) NPTS
9997  FORMAT (35X,19(1H-)/'  Point',I3,25X,'Nuclear coordinates'/
     1 35X,19(1H-)//17X,' Atom',13X,'x',19X,'y',19X,'z'/)
      WRITE (10,9998) (I,ANAM(I),(C(J,I),J=1,3), I=1,NAT)
      WRITE (IW,9998) (I,ANAM(I),(C(J,I),J=1,3), I=1,NAT)
9998  FORMAT (5X,I3,2X,A8,2X,3F20.10)

      CALL STRUCT(NAT, C, ZAN)

16    IF (IREST.le.3) then

C  Calculate HF energy

        write (*,'(/a)') ' . VALOPT is calling HFSCF'
        write (*,'(/a)') ' (Starting SCF procedure in the current geomet
     &ry)'
        CALL HFSCF
        call timit(4)
      else
        CALL SECGET(ISEX(13),13,IBLOK)
        CALL READQ(EN,LDA(ISEX(13)),IBLOK,IFILD)
        WRITE(10,*) ETOT
        WRITE(9,*) ETOT
      endif

      FUNC=ETOT
      E=ETOT
      IF(FUNC.EQ.ZERO) CALL HONDER(77,'VALOPT')
      IF(FUNC.GE.FUNC0) then
        write (9,1001) func, func0
        write (10,1001) func, func0
1001    format (/
     &  ' WARNING:'//
     &  ' New SCF total energy:',t40, f20.12, ' Hartrees'/
     &  ' Former SCF total energy:',t40,f20.12,' Hartrees'//
     &  ' Because energy increased, no gradients are calculated in this 
     &cycle.'/)
        RETURN
      endif

C  Calculate gradient

      DO 20 I=1,NC
   20   EG(I)=ZERO

      write (*,'(/a)') ' . VALOPT is calling HFGRAD'
      write (*,'(a)') ' (Searching gradients in the current geometry)'
      CALL HFGRAD
      call timit(4)

22    IF (NVAR.NE.0) GO TO 30
      DO 25 I=1,NC
        IF(DABS(EG(I)).LT.TOLG) EG(I)=ZERO
25      G(I)=EG(I)
      RETURN

30    CALL POSN(.TRUE.)
      CALL SYMGEN(.TRUE.)
      RETURN
      END

      SUBROUTINE OPTEND

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL LADAPT,LSKIP
      CHARACTER*8 ANAM,ATNAME,PARNAM,BFLAB
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      LOGICAL FIXED,LEX,LDAM12,LDAM13,LDAM23,LDIIS
      COMMON/SCFBLK/EN,ETOT,EHF,SH1(2),SH2(2),GAP1(2),GAP2(2),
     1              D12,D13,D23,CANA,CANB,CANC,FIXED,LEX,
     2              LDAM12,LDAM13,LDAM23,LDIIS,
     3              NCYC,ISCHM,LOCK,MAXIT,NCONV,LOKCYC
      PARAMETER (NDIMQ=NATM13*NATM13+2799, UAEV=27.211611D0)
      COMMON/BIG/Q(NDIMQ),QUNK

      IBLOK=0
      write(10,10) ' ',TITLE
      WRITE(IW,10) CHAR(12),TITLE
10    FORMAT(A//1X,9A8,A7/1X,79('-')/)
      CALL SECGET(ISEX(13),13,IBLOK)
      CALL READQ(EN,LDA(ISEX(13)),IBLOK,IFILD)
      WRITE(IW,20) EN,EN*UAEV,EHF,EHF*UAEV,ETOT,ETOT*UAEV
      write(10,20) EN,EN*UAEV,EHF,EHF*UAEV,ETOT,ETOT*UAEV
20    FORMAT(/T5,'NUCLEAR ENERGY    =',F20.12,' a.u. =',F20.9,' ev'/,
     1        T5,'ELECTRONIC ENERGY =',F20.12,' a.u. =',F20.9,' ev'/,
     2        T5,'TOTAL ENERGY      =',F20.12,' a.u. =',F20.9,' ev')
      NOC=NA
      IPOS=8
50    CONTINUE
      NDIM=NUM
      I1=NX+1
      I2=I1+NUM*NUM
      CALL SECGET(ISEX(IPOS),IPOS,IBLOK)
      CALL READQ(Q(I1),NUM*NCOORB,IBLOK,IFILD)
      IPOS=IPOS+1
      CALL SECGET(ISEX(IPOS),IPOS,IBLOK)
      CALL READQ (Q(I2),NCOORB,IBLOK,IFILD)
      WRITE(IW,70)
      DO 60 I=1,NOC
60      WRITE(IW,80) I,Q(I2+I-1),Q(I2+I-1)*UAEV
70    FORMAT(/,10X,16(1H-),/,10X,'ORBITAL ENERGIES',/,10X,16(1H-)/
     &T32,'au',20X,'ev'/)
80    FORMAT(T18,I5,F20.12,F20.9)
      WRITE(IW,90)
90    FORMAT(/,10X,18(1H-),/,10X,'MOLECULAR ORBITALS',/,10X,18(1H-))
      CALL VECOUT(Q(I1),Q(I2),NUM,NOC,NDIM,IW,BFLAB)
120   CALL STRUCT(NAT,C,ZAN)
      RETURN
      END

      SUBROUTINE POSN(LDERIV)

C  Calculate the positions of the atoms and the derivatives of the
C  positions w.r.t. the variable parameters.

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL LDERIV

      INTEGER P1, P2, P3
      CHARACTER*8 BFLAB,PARNAM,ATNAME,ANAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), X(3,NATM1)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      LOGICAL ACTIVE
      PARAMETER (NQTAP=2800-NATM13/2)
      COMMON/BIG/DXDP(3,NATM1,NATM13),ACTIVE(NATM13),QJNK(NQTAP)
      PARAMETER (ZERO=0.D0)

C   DXDP(K,N,M)   Derivative w.r.t. the Mth parameter of the Kth
C                 cartesian coordinate of the Nth atom.

C   ACTIVE(J)     .TRUE. if any atom position so far depends on
C                 parameter J (set by FNDPAR).

      IF (.NOT. LDERIV) GO TO 20
      DO 10 I=1,3
        DO 10 N=1,NUQ
          DO 10 J=1,NVAL
10          DXDP(I,N,J)=ZERO
      DO 15 J=1,NVAL
15      ACTIVE(J)=.FALSE.

20    DO 1000 N=1,NUQ

C  Branch on KONTYP

        KC=KONTYP(N)
        IF (KC .GT. 1) GO TO 200

C  Cartesian coordinates
100     DO 110 I=1,3
          CALL FNDPAR(N,I,L,SGN,VALUE)
          X(I,N)=ZERO
          IF (L .NE. 0) X(I,N)=VALUE
110       IF (LDERIV .AND. L .NE. 0 .AND. L .LE. NVAR) DXDP(I,N,L)=SGN
        GO TO 1000

C  Polar coordinates

200     IF (KC .GT. 2) GO TO 300
        CALL FNDPAR(N,1,LR,SGNR,R)
        CALL FNDPAR(N,2,LA,SGNTH,THETA)
        CALL FNDPAR(N,3,LB,SGNPH,PHI)

        COSTH=DCOS(THETA)
        SINTH=DSIN(THETA)
        COSPH=DCOS(PHI)
        SINPH=DSIN(PHI)
        X(1,N)=R*SINTH*COSPH
        X(2,N)=R*SINTH*SINPH
        X(3,N)=R*COSTH
        IF (.NOT. LDERIV) GO TO 1000
        IF (LR .EQ. 0 .OR. LR .GT. NVAR) GO TO 210
        DXDP(1,N,LR)=SGNR*SINTH*COSPH
        DXDP(2,N,LR)=SGNR*SINTH*SINPH
        DXDP(3,N,LR)=SGNR*COSTH
210     IF (LA .EQ. 0 .OR. LA .GT. NVAR) GO TO 220
        DXDP(1,N,LA)=R*SGNTH*COSTH*COSPH
        DXDP(2,N,LA)=R*SGNTH*COSTH*SINPH
        DXDP(3,N,LA)=-R*SGNTH*SINTH
220     IF (LB .EQ. 0 .OR. LB .GT. NVAR) GO TO 1000
        DXDP(1,N,LB)=-R*SINTH*SGNPH*SINPH
        DXDP(2,N,LB)=R*SINTH*SGNPH*COSPH
        DXDP(3,N,LB)=ZERO
        GO TO 1000
300     CONTINUE
1000    CONTINUE
C
      RETURN
      END

C  SYMGEN  A1
      SUBROUTINE SYMGEN(LDERIV)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      COMMON /BIG/DXDP(3,NATM1,NATM13),QJNK(2800)
      LOGICAL LDERIV

      INTEGER P1, P2, P3
      CHARACTER*8 ANAM,ATNAME,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5     PARVAL(NATM13), ZIN(NATM1), X(3,NATM1)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)

      CHARACTER*8 GROUP
      COMMON/SYMTRY/T(216),INVT(24),NEWSH(NSH1+1,24),NT,NT2,
     1    XOLD,YOLD,ZOLD,XNEW,YNEW,ZNEW,XP,YP,ZP
     2    ,U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
     3    ,INDEX,NAXIS
      COMMON/SYMTRB/GROUP

      COMMON/OPTMIZ/V(NATM13),G(NATM13),DV(NATM13),FUNC,ALPHA,
     *    V0(NATM13),G0(NATM13),GS0,ALPH
     *    ,FUNC0,NCOORD,NPTS,NSERCH,IUPDAT,ICODE,ISS,BIGG,BIGG0
     *    ,Z(NATM13),Q(NATM13),CC,GNRM

      COMMON/FUNCT/E,EG(3,NATM1)
      PARAMETER (ZERO=0.D0)
      DATA TOL /1.0D-6/

C  Transfer information from the Z-matrix to common /INFOA/, generating
C  symmetry-equivalent centres and deleting dummies.  Basis set
C  information is assumed to be present already.  If LDERIV is set,
C  the derivatives are symmetry-generated also and used to transform the
C  gradient w.r.t. nuclear coordinates (held in EG) to the gradient w.r.t.
C  variable internal coordinates (in G).

      N=0
      NAT=0
      IF (.NOT. LDERIV) GO TO 100
      DO 10 J=1,NVAR
10      G(J)=ZERO

100   N=N+1
      IF (N .GT. NUQ) RETURN

101   IF (ZIN(N) .LT. ZERO) GO TO 100
      NAT=NAT+1
      ZAN(NAT)=ZIN(N)
      ANAM(NAT)=ATNAME(N)
      DO 110 I=1,3
110     C(I,NAT)=X(I,N)

C  Generate equivalent centres

      IF (NT .EQ. 1) GO TO 160
      XS=U1*(X(1,N)-X0)+U2*(X(2,N)-Y0)+U3*(X(3,N)-Z0)
      YS=V1*(X(1,N)-X0)+V2*(X(2,N)-Y0)+V3*(X(3,N)-Z0)
      ZS=W1*(X(1,N)-X0)+W2*(X(2,N)-Y0)+W3*(X(3,N)-Z0)

C  Loop over symmetry-generated points

      DO 150 IT=2,NT
        NN=9*(IT-1)
        XNEW=XS*T(NN+1)+YS*T(NN+2)+ZS*T(NN+3)
        YNEW=XS*T(NN+4)+YS*T(NN+5)+ZS*T(NN+6)
        ZNEW=XS*T(NN+7)+YS*T(NN+8)+ZS*T(NN+9)
        XP=X0+U1*XNEW+V1*YNEW+W1*ZNEW
        YP=Y0+U2*XNEW+V2*YNEW+W2*ZNEW
        ZP=Z0+U3*XNEW+V3*YNEW+W3*ZNEW
        DO 130 IAT=1,NAT
          TEST=(XP-C(1,IAT))**2+(YP-C(2,IAT))**2+(ZP-C(3,IAT))**2
          IF(TEST.LE.TOL) GO TO 150
130       CONTINUE
        NUAT=NAT
        NAT=NAT+1
        C(1,NAT)=C(1,NUAT)
        C(2,NAT)=C(2,NUAT)
        C(3,NAT)=C(3,NUAT)
        C(1,NUAT)=XP
        C(2,NUAT)=YP
        C(3,NUAT)=ZP
        ZAN(NAT)=ZIN(N)
        ANAM(NAT)=ATNAME(N)

        IF (.NOT. LDERIV) GO TO 150

        DO 140 J=1,NVAR
          XD=U1*DXDP(1,N,J)+U2*DXDP(2,N,J)+U3*DXDP(3,N,J)
          YD=V1*DXDP(1,N,J)+V2*DXDP(2,N,J)+V3*DXDP(3,N,J)
          ZD=W1*DXDP(1,N,J)+W2*DXDP(2,N,J)+W3*DXDP(3,N,J)
          XNEW=XD*T(NN+1)+YD*T(NN+2)+ZD*T(NN+3)
          YNEW=XD*T(NN+4)+YD*T(NN+5)+ZD*T(NN+6)
          ZNEW=XD*T(NN+7)+YD*T(NN+8)+ZD*T(NN+9)
          G(J)=G(J) + EG(1,NUAT)*(U1*XNEW+V1*YNEW+W1*ZNEW)
     *              + EG(2,NUAT)*(U2*XNEW+V2*YNEW+W2*ZNEW)
     *              + EG(3,NUAT)*(U3*XNEW+V3*YNEW+W3*ZNEW)
140       CONTINUE

150     CONTINUE

160   IF (.NOT. LDERIV) GO TO 100
      DO 170 J=1,NVAR
        DO 170 I=1,3
170       G(J)=G(J)+EG(I,NAT)*DXDP(I,N,J)

      GO TO 100

      END

C  FNDPAR  A1

      SUBROUTINE FNDPAR(N,I, L,SGN,VALUE)

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      INTEGER P1, P2, P3
      CHARACTER*8 PARNAM,ATNAME,ANAM,BFLAB
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1), K(NATM1,3),
     4    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)

      LOGICAL ACTIVE
      PARAMETER (NQTAP=2800-NATM13/2)
      COMMON /BIG/ DXDP(NATM13,NATM13), ACTIVE(NATM13),QJNK(NQTAP)
      PARAMETER (ZERO=0.D0, ONE=1.D0)

      L=IABS(K(N,I))
      SGN=ONE
      VALUE=ZERO
      IF (L .EQ. 0) RETURN
      SGN=ISIGN(1,K(N,I))
      VALUE=PARVAL(L)*SGN
      ACTIVE(L)=.TRUE.

      RETURN
      END

C     AUTODBL

      SUBROUTINE MULTRC(A,X,Y,NCOORD,IBLK,IFILD)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      DIMENSION X(NCOORD),Y(NCOORD),A(NATM13,NCOORD)
      PARAMETER (ZERO=0.0D0)
      DATA NDIM/NATM13/

      CALL READM(A,NCOORD,NCOORD,NDIM,IBLK,IFILD)
      DO 20 I=1,NCOORD
        DUM=ZERO
        DO 10 J=1,NCOORD
   10     DUM=DUM+A(I,J)*X(J)
   20   Y(I)=DUM

      RETURN
      END

C  SYMDER  A1

      SUBROUTINE SYMDR(DR)

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      INTEGER P1, P2, P3
      CHARACTER*8 ANAM,ATNAME,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), XQ(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      CHARACTER*8 GROUP
      COMMON/SYMTRY/T(216),INVT(24),ISO(NSH1+1,24),NT,NT2,
     1    XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
     2    ,U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
     3    ,INDEX,NAXIS
      COMMON/SYMTRB/GROUP
      DIMENSION DR(3,NATM1)
      PARAMETER (NQTAP=2800-(216+(NATM1*24+24)/2))
      COMMON/BIG/HH(NATM13,NATM13),PTR(3,72),ICT(NATM1,24),
     1           INDIN(24),QJNK(NQTAP)
      PARAMETER (ZERO=0.D0, ONE=1.0D0)

      IF(NT.EQ.1) RETURN

C     ----- CALCULATE TRANFORMATION MATRICES OF COORDINATES. -----

      X=X0+ONE
      Y=Y0
      Z=Z0
      CALL LOCAL(X,Y,Z,XS,YS,ZS)
      XSMAL=XS
      YSMAL=YS
      ZSMAL=ZS
      DO 100 IT=1,NT
        NN=9*(IT-1)
        CALL TRANS(NN)
        N=3*(IT-1)
        NP1 = N+1
        PTR(1,NP1)=XP-X0
        PTR(2,NP1)=YP-Y0
        PTR(3,NP1)=ZP-Z0
  100   CONTINUE
      X=X0
      Y=Y0+ONE
      Z=Z0
      CALL LOCAL(X,Y,Z,XS,YS,ZS)
      XSMAL=XS
      YSMAL=YS
      ZSMAL=ZS
      DO 200 IT=1,NT
        NN=9*(IT-1)
        CALL TRANS(NN)
        N=3*(IT-1)
        NP2 = N+2
        PTR(1,NP2)=XP-X0
        PTR(2,NP2)=YP-Y0
        PTR(3,NP2)=ZP-Z0
  200   CONTINUE
      X=X0
      Y=Y0
      Z=Z0+ONE
      CALL LOCAL(X,Y,Z,XS,YS,ZS)
      XSMAL=XS
      YSMAL=YS
      ZSMAL=ZS
      DO 300 IT=1,NT
        NN=9*(IT-1)
        CALL TRANS(NN)
        N=3*(IT-1)
        NP3 = N+3
        PTR(1,NP3)=XP-X0
        PTR(2,NP3)=YP-Y0
        PTR(3,NP3)=ZP-Z0
  300   CONTINUE

C     ----- SET TRANSFORMATION TABLE: ATOMS VERSUS SYMMETRY OPERATIONS.

      DO 1200 II=1,NSHELL
        IC=KATOM(II)
        DO 1000 IT=1,NT
 1000     INDIN(IT)=ISO(II,IT)
        DO 1100 IT=1,NT
          ID=INDIN(IT)
 1100     ICT(IC,IT)=KATOM(ID)
 1200   CONTINUE

C     ----- SYMMETRYZE DISPLACEMENT VECTOR -----

      DO 2500 IC=1,NAT
        DO 2100 IT=1,NT
          IF(ICT(IC,IT).GT.IC) GO TO 2500
 2100     CONTINUE
        DX=ZERO
        DY=ZERO
        DZ=ZERO
        DO 2200 IT=1,NT
          ICNU=ICT(IC,IT)
          DXP=DR(1,ICNU)
          DYP=DR(2,ICNU)
          DZP=DR(3,ICNU)
          N=3*(IT-1)
          NP1 = N+1
          NP2 = N+2
          NP3 = N+3
          DX=DX+DXP*PTR(1,NP1)+DYP*PTR(2,NP1)+DZP*PTR(3,NP1)
          DY=DY+DXP*PTR(1,NP2)+DYP*PTR(2,NP2)+DZP*PTR(3,NP2)
          DZ=DZ+DXP*PTR(1,NP3)+DYP*PTR(2,NP3)+DZP*PTR(3,NP3)
 2200     CONTINUE
        DR(1,IC)=DX
        DR(2,IC)=DY
        DR(3,IC)=DZ
        DO 2400 IT=1,NT
          ICNU=ICT(IC,IT)
          IF(ICNU.EQ.IC) GO TO 2400
          IF(IT.EQ.NT) GO TO 2350
          IT1=IT+1
          DO 2300 JT=IT1,NT
            IF(ICT(IC,JT).EQ.ICNU) GO TO 2400
 2300       CONTINUE
 2350     CONTINUE
          JT=INVT(IT)
          N=3*(JT-1)
          NP1 = N+1
          NP2 = N+2
          NP3 = N+3
          DXP=DX*PTR(1,NP1)+DY*PTR(2,NP1)+DZ*PTR(3,NP1)
          DYP=DX*PTR(1,NP2)+DY*PTR(2,NP2)+DZ*PTR(3,NP2)
          DZP=DX*PTR(1,NP3)+DY*PTR(2,NP3)+DZ*PTR(3,NP3)
          DR(1,ICNU)=DXP
          DR(2,ICNU)=DYP
          DR(3,ICNU)=DZP
 2400     CONTINUE
 2500   CONTINUE
      DUM=DBLE(NT)
      DO 2600 N=1,NAT
        DO 2600 I=1,3
 2600     DR(I,N)=DR(I,N)/DUM

      RETURN
      END

C  SETS0  A1
      SUBROUTINE SETS0(S,NCOORD,IBLK,IFILD)

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      DIMENSION S(NATM13,NCOORD)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      DATA NDIM/NATM13/

      DO 20 I=1,NCOORD
        DO 10 J=1,NCOORD
   10     S(I,J)=ZERO
   20   S(I,I)=ONE
      CALL WRT3M(S,NCOORD,NCOORD,NDIM,IBLK,IFILD)
      RETURN
      END

C  SETS1  A1
      SUBROUTINE SETS1(S,Z,C,NCOORD,IBLK,IFILD)

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      DIMENSION S(NATM13,NCOORD),Z(NCOORD)
      DATA NDIM/NATM13/
      CALL READM(S,NCOORD,NCOORD,NDIM,IBLK,IFILD)

      DO 10 J=1,NCOORD
        DO 10 I=J,NCOORD
          S(I,J)=S(I,J)+Z(I)*Z(J)/C
   10     S(J,I)=S(I,J)
      CALL WRT3M(S,NCOORD,NCOORD,NDIM,IBLK,IFILD)
      RETURN
      END
