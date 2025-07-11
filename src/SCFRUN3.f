      SUBROUTINE SCFHS
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'

      COMMON/BIG/QQQ(MDM1),QXTRA
      COMMON/BIGG/QQ(MDM1),QQXTRA,Q(NBF1M+NBF1)
      LOGICAL LSWOP,LOKNXT
      LOGICAL NEARLY,CONVGD,LSTOP,EDIS
      LOGICAL FIXED,LEX,LDAM12,LDAM13,LDAM23,LDIIS
      COMMON/SCFBLK/EN,ETOT,EHF,SH1(2),SH2(2),GAP1(2),GAP2(2),
     1              D12,D13,D23,CANA,CANB,CANC,FIXED,LEX,
     2              LDAM12,LDAM13,LDAM23,LDIIS,
     3              NCYC,ISCHM,LOCK,MAXIT,NCONV,LOKCYC
      LOGICAL  LADAPT,LSKIP
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
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMMAX,TTOT
      CHARACTER*8 PARNAM,ATNAM,ANAM,BFLAB
      INTEGER P1,P2,P3
      COMMON/INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(NATM1),C(3,NATM1),
     &             NVAR,NPAR,NVAL,NUQ,
     &             CFR,CFA,KONTYP(NATM1),P1(NATM1),P2(NATM1),P3(NATM1)
     &             ,KR(NATM1),KA(NATM1),KB(NATM1),PARVAL(NATM13)
     &             ,ZIN(NATM1),X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAM(NATM1)
      COMMON/TEMP2/IOSF(50)
      COMMON/MAPPER/IA(NBF1+1)

      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,TEN=1.0D1,M511=511,
     &           M340=340,MVADD=0)

C     High-spin open-shell SCF proceedure from Roothaan
C     but programmed in manner of ATMOL

10    FORMAT(/' ----- NUCLEAR ENERGY ----- = ',F20.12)
20    FORMAT(/,1X,'CONVERGENCE DATA'//
     1 ' MAXIMUM NUMBER OF ITERATIONS = ',I6,/,
     2 ' CONVERGENCE CRITERION        = 1.0D-',I2,//,
     4 '  CYC   TOTAL        ELECTRONIC        DEL(E)    DEL(D) ',
     5 'DEL(T)     RUN   SHA  SHB',/
     6 '        ENERGY',T24,'ENERGY',T69,'TIME'/)
60    FORMAT(//' CLOSED SHELL ***** '/)
70    FORMAT(//' OPEN SHELL ****** '/)
80    FORMAT(/10X,11(1H-),/,10X,'EIGENVALUES',/,10X,11(1H-)/)
90    FORMAT(6X,I5,F20.12)
100   FORMAT(/10X,12(1H-),/,10X,'EIGENVECTORS',/,10X,12(1H-)/)
110   FORMAT(/10X,14(1H-),/,10X,'DENSITY MATRIX',/,10X,14(1H-))
120   FORMAT(I4,F12.5,2F15.7,F10.7,  F6.1,F9.1,2F4.1)
130   FORMAT(/10X,17(1H-),/,10X,'DENSITY CONVERGED',/,10X,17(1H-))
140   FORMAT('... SCF HAS CONVERGED BUT PROGRAM STOPPED',/,
     1 ' FOR TIMLIM BEFORE ULTIMATE CYCLE ...')
150   FORMAT(' ..... SCF HAS NOT CONVERGED ...',
     1 ' ... JOB TO BE RESTARTED ....')
160   FORMAT(/,10X,30(1H-),/,10X,'EXCESSIVE NUMBER OF ITERATIONS',
     1 /,10X,30(1H-))
170   FORMAT(5F15.7)
180   FORMAT(//10X,19(1H-)/10X,'RHF OPEN-SHELL CALCULATION'/
     1 10X,19(1H-))
190   FORMAT(/10X,'FINAL VALUES OF TOTAL ENERGY AND ELECTRONIC'
     1 ,' ENERGY ARE ',/,25X,F16.9,' AND ',F16.9)

      write(10,180)
      WRITE(IW,180)

C     ----- INITIALIZE VARIABLES -----

      TRHMIN=1.0D-16
      TRHMAX=1.0D-08
      AUTO=0.25D0
      ISEC4=0
      ISECS=0
      ISECH=0
      ISECDA=0
      ISECV=0
      ISECDB=0
      ISC24=0
      CALL TRMAT
      DIFF0=ZERO
      LOKNXT=.FALSE.
      DIFF=100.0
      DELDIF=ZERO
      NDIM=NUM
      I1=1
      I2=NX+I1
      I3=1
      I4=NX+I3
      I6=NX+1
      I7=NUM+I6
      J1=I1+NUM*NUM
      J2=I3+NUM*NUM
      I5=2*NX
      write(10,123) MAXDIM,I5
      WRITE(IW,123) MAXDIM,I5
123   FORMAT(/15X,'SCF : core available',I12,', used',I12)
      IF(I5.GT.MAXDIM) write(* ,1001)
      IF(I5.GT.MAXDIM) WRITE(IW,1001)
1001  FORMAT(1H ,'*** INSUFFICIENT CORE ***')

C     ----- SET OCCUPATION NUMBERS -----

      ITYPE=1
      PENAL=ZERO
      NA1=NA+1
      NB1=NB+1
      DAMP=ONE
      DO 200 I=1,NCOORB
        Q(I6+I-1)=ZERO
        Q(I7+I-1)=ZERO
        IF(I.LE.NA)Q(I6+I-1)=TWO
        IF(I.GT.NA.AND.I.LE.NB)Q(I7+I-1)=ONE
200     CONTINUE
      LOK=LOKCYC
      EHFM1=ZERO
      EHFM2=ZERO
      EHF=ZERO
      EHF0=ZERO
      ACURCY=TEN**(-NCONV)
      FACUR=ACURCY*ACURCY
      IF(FACUR.LT.TRHMIN)FACUR=TRHMIN
      EDIFF=ZERO
      ITER=0
      ICA=0
      EDIS=.FALSE.
      LSTOP=.FALSE.
      NEARLY=.FALSE.
      CONVGD=.FALSE.

      CALL SECGET(ISEX(24),24,ISC24)
      CALL SECGET(ISEX(4),4,ISEC4)
      CALL SECGET(ISEX(5),5,ISECS)
      CALL SECGET(ISEX(6),6,ISECH)
      CALL SECGET(ISEX(7),7,ISECDA)
      MTYPE=0
      CALL SECGET(ISEX(8),MTYPE,ISECV)
      CALL SECGET(ISEX(10),10,ISECDB)

      CALL TIMIT(3)
      TIM2=TIM
*      TIM0=TIM

C     ----- NUCLEAR ENERGY

      EN=ENUC(NAT,ZAN,C)
      write(10 ,10) EN
      WRITE(IW,10) EN
      IF(.NOT.FIXED) THEN
        write(10,240)GAP1,NCYC,GAP2
        WRITE(IW,240)GAP1,NCYC,GAP2
      ENDIF
240   FORMAT(/5X,'VARIABLE LEVEL SHIFTERS USED'
     1/5X,'MINIMUM  SEPARATIONS',F6.2,' AND',F6.2,
     2' UNTIL CYCLE',I5,' THEN',F6.2,' AND',F6.2)
      IF(FIXED) THEN
        write(10,250) SH1,NCYC,SH2
        WRITE(IW,250) SH1,NCYC,SH2
      ENDIF
      IF(LOCK.NE.0) THEN
        write(10,260) LOKCYC
        WRITE(IW,260) LOKCYC
      ENDIF
250   FORMAT(/5X,'LEVEL SHIFTERS',F7.2,' AND',F7.2,'  UNTIL CYCLE',
     1I5,'  THEN',F7.2,'  AND',F7.2)
260   FORMAT(/5X,'LOCK OPTION IN FORCE FROM CYCLE ',I6)
      write(10,270)CANA,CANB,CANC
      WRITE(IW,270)CANA,CANB,CANC
270   FORMAT(/5X,'CANONICALISATION FACTORS',3F12.3)
      write(10,275) D12,D13,D23
      WRITE(IW,275) D12,D13,D23
      write(10,20) MAXIT,NCONV
      WRITE(IW,20) MAXIT,NCONV
275   FORMAT(/5X,'DAMP FACTORS -- D12 =',F6.2,' D13 =',F6.2,
     &' D23 =',F6.2)

      IOSF(8)=IOSF(7)+LENSEC(NX)
      IOSF(9)=IOSF(8)+LENSEC(NX)
      IOSF(10)=IOSF(9)+LENSEC(NX)
      IOSF(11)=IOSF(10)+LENSEC(NX)

C     READ INITIAL DENSITY MATRICES

      CALL READQ(Q(1),NX,ISECDA,IFILD)
      CALL READQ(QQQ(I4),NX,ISECDB,IFILD)

      IF(.NOT.LEX)GOTO 285
      CALL WRT3(Q,NX,IOSF(2),IFILS)
      CALL WRT3S(Q,NX,IFILS)
      CALL WRT3S(Q,NX,IFILS)
285   CONTINUE

C     Q(1)=DC (CLOSED)
C     QQQ(I4)=DO  (OPEN)
C     DENSITY MATRICES INCLUDE OCCUPATION NUMBERS IN DEFINITION

C     ----- START SCF PROCEDURE -----

290   CALL SEARCH(IBLKM,IFILM)

C     ----- CONSTRUCT A SKELETON FOCK MATRIX -----

C     Q(1)=DC
C     QQ(I1)=J(R)-K(R)/2 ; R=DC+DO
C     QQQ(I4)=DO
C     QQQ(I3)=K(DO)/2
C     QQ(I2)=WORKSPACE

      CALL HSTHS(Q(1),QQ(I1),QQQ(I4),QQQ(I3),QQ(I2),IA)

C     ----- SYMMETRIZE SKELETON FOCK MATRIX -----

C     QQ(I1)=J-K/2
C     QQ(I2)=WORKSPACE
C     QQQ(I3)=K/2

      IF(.NOT.LADAPT) THEN
        CALL SYMH(QQ(I1),QQ(I2),IA)
        CALL SYMH(QQQ(I3),QQ(I2),IA)
      ENDIF

C     ----- READ IN CORE HAMILTONIAN MATRIX

      ITER=ITER+1
      EHFM2=EHFM1
      EHFM1=EHF0
      EHF0=EHF
      CALL READQ(QQ(I2),NX,ISECH,IFILD)
      DO 300 I=1,NX
        J=I-1
        DUM=QQ(I2+J)+QQ(I1+J)
        QQ(I1+J)=DUM
        QQ(I2+J)=QQ(I2+J)+DUM
        Q(I)=Q(I)+QQQ(I4+J)
300     CONTINUE

C     Q(1)=DT ; QQ(I1)=F=H0+J-K/2 ; QQ(I2)=H0+F ; QQQ(I3)=K/2 ; QQQ(I4)=DO

      EHF=TRACEP(Q(1),QQ(I2),NUM)-TRACEP(QQQ(I4),QQQ(I3),NUM)
      EHF=EHF/TWO

C     STORE K MATRIX AND TOTAL DENSITY MATRIX ON SCRATCHFILE

      CALL WRT3(QQ(I1),NX,IOSF(7),IFILS)
      CALL WRT3(Q(1),NX,IOSF(8),IFILS)
      CALL WRT3(QQQ(I3),NX,IOSF(9),IFILS)

C     READ IN THE PREVIOUS SET OF M.O.'S TO TRANSFORM FOCK MATRIX

      CALL READQ(QQQ(I3),NUM*NCOORB,ISECV,IFILD)

C      TRANSFORM

      CALL FOCKTR(Q(1),QQ(I1),QQQ(I3),QQQ(J2),IA,NCOORB,NUM,NDIM)

C       TRANSFORMED F NOW IN Q(1)
C       STORE ON SCRATCH FILE

      CALL WRT3(Q(1),NX,IOSF(10),IFILS)

C     GET THE K MATRIX

      CALL READQ(QQ(I1),NX,IOSF(9),IFILS)
      CALL FOCKTR(Q(1),QQ(I1),QQQ(I3),QQQ(J2),IA,NCOORB,NUM,NDIM)
330   CALL WRT3(Q(1),NX,IOSF(9),IFILS)

C     TRANSFORMED K NOW IN Q(1)
C     GET THE TRANSFORMED F BACK AND CONSTRUCT TOTAL OPEN-SHELL
C     FOCK MATRIX

      CALL READQ(QQ(I1),NX,IOSF(10),IFILS)
      CALL FFORM(Q(1),QQ(I1),NA,NB,NCOORB,CANA,CANB,CANC,IA,IW,NPRINT)

C     ATTEMPT EXTRAPOLATION

      IF(.NOT.LEX)GOTO 543
      CALL DTR(Q(1),QQ(I1),QQQ(I3),QQQ(J2),NUM,NCOORB,NDIM)
      CALL READQ(QQ(I2),NX,ISECS,IFILD)
      IJ=-1
      DO 111 I=1,NUM
        DO 111 J=1,I
          IJ=IJ+1
          SIJ=QQ(I2+IJ)
          QQQ((I-1)*NUM+J)=SIJ
111       QQQ((J-1)*NUM+I)=SIJ
      CALL FOCKTR(Q(1),QQ(I1),QQQ(1),QQQ(J2),IA,NUM,NUM,NDIM)
      CALL EXTRAP(QQ(I1),QQ(I2),QQQ(1),Q(1),ICA,IW,IFILS,IOSF,
     &   0,EDIS)
      CALL FMOVE(Q(1),QQ(I1),NX)
      CALL READQ(QQQ(I3),NUM*NCOORB,ISECV,IFILD)
      CALL FOCKTR(Q(1),QQ(I1),QQQ(I3),QQQ(J2),IA,NCOORB,NUM,NDIM)

543   THRESH=DIFF*DIFF
      IF(THRESH.GT.TRHMAX)THRESH=TRHMAX
      IF(THRESH.LT.FACUR)THRESH=FACUR
      LSWOP=ITER.GT.1.AND.ITER.LE.LOK.AND.DIFF.LT.AUTO.AND..NOT.LOKNXT
      DIFF0=DIFF
      DIFF=ZERO

C     FIND LEVEL SHIFTERS

      IF(ITER.GE.4)CALL PNLTY(PENAL,DAMP,EHF,EHF0,EHFM1,EHFM2,LSWOP,
     1DELDIF,ICA)
      CALL LEVEL(SHA,Q(1),SH1(1),SH2(1),GAP1(1),GAP2(1),
     1NCYC,ITER,NA,NB,FIXED)
      CALL LEVEL(SHB,Q(1),SH1(2),SH2(2),GAP1(2),GAP2(2),
     1NCYC,ITER,NB,NCOORB,FIXED)
      SHA=SHA+PENAL/TWO
      SHB=SHB+PENAL/TWO

C     ADD LEVEL SHIFTERS

      CALL LSHIFT(Q(1),NCOORB,NA,IA,SHA,ONE)
      CALL LSHIFT(Q(1),NCOORB,NB,IA,SHB,ONE)

C     DAMPING

      IF(LDAM12)CALL DMPER(Q(1),D12,NA+1,NB,1,NA)
      IF(LDAM13)CALL DMPER(Q(1),D13,NB+1,NCOORB,1,NA)
      IF(LDAM23)CALL DMPER(Q(1),D23,NB+1,NCOORB,NA+1,NB)

C     ----- DIAGONALIZE NEW HAMILTONIAN MATRIX -----

*      call hqrii1 (NCOORB,1,NCOORB,-1,Q(1),QQ(J1),NDIM,QQ(I1),.TRUE.)
      CALL LIGEN(Q(1),QQ(I1),QQ(J1),IA,NCOORB,NDIM,ITYPE,THRESH)

C     LOCK OPTION

      IF(ITER.GE.LOK.OR.LOKNXT)CALL LOCKER(QQ(J1),QQ(I1),NCOORB,NCOORB,
     1NDIM)

C     REMOVE LEVEL SHIFTERS

      IF(NA1.GT.NB)GOTO 360
      DO 350 I=NA1,NB
        II=J1+I-1
350     QQ(II)=QQ(II)-SHA
360   IF(NB1.GT.NCOORB)GOTO 380
      DO 370 I=NB1,NCOORB
        II=J1+I-1
370     QQ(II)=QQ(II)-SHA-SHB

C     TRANSFORM VECTORS BACK TO A.O. BASIS

380   IF(NPRINT.NE.5)GOTO 390
      WRITE(IW,100)
      CALL VECOUT(QQ(I1),QQ(J1),NUM,NCOORB,NDIM,IW,BFLAB)
390   CONTINUE
      CALL BACKTR(QQ(I1),QQQ(I3),QQQ(J2),NCOORB,NUM,NDIM)

C     ORDER IF NECESSARY

      LOKNXT=.FALSE.
      IIII=0
      IF(LSWOP)CALL ORDER(QQ(I1),QQ(J1),NUM,NCOORB,
     1NDIM,NB,IIII)
      IF(IIII.NE.0)THEN
        ICA=0
        LOKNXT=.TRUE.
      ENDIF

C     SCHMIDT ORTHOGONALISE EVERY ISCHM CYCLES

      IF((ITER/ISCHM)*ISCHM.NE.ITER)GOTO 410
        CALL READQ(Q(1),NX,ISECS,IFILD)
        CALL ORTHO(QQ(I1),Q(1),QQQ(J2),NUM,NCOORB,NDIM,IA)

410   CALL DMTX(Q(1),QQ(I1),Q(I6),NA,NUM,NDIM)
      CALL DMTX(QQQ(I4),QQ(I1),Q(I7),NB,NUM,NDIM)

C     ----- SAVE MO'S + DENSITY + ORBITAL ENERGIES -----

      IPOS=7
      CALL SAVEMO(Q(1),QQ(I1),QQ(J1),NDIM,IPOS)
      CALL WRT3(QQQ(I4),NX,ISECDB,IFILD)

C     GET PREVIOUS DENSITY MATRIX OFF SCRATCH FILE AND COMPARE

      CALL READQ(QQQ(I3),NX,IOSF(8),IFILS)
      DO 440 I=1,NX
        J=I-1
        DDD=DABS(QQQ(I3+J)-Q(I)-QQQ(I4+J))
440     IF(DDD.GT.DIFF)DIFF=DDD

      DELDIF=DIFF0-DIFF
      EDIFF=EHF-EHF0
      ETOT=EHF+EN
      CALL TIMIT(3)
      TI1=TIM
      DELT=TI1-TIM2
      TIM2=TI1
      write(10,120) ITER,ETOT,EHF,EDIFF,DIFF,DELT,TI1,SHA,SHB
      WRITE(IW,120) ITER,ETOT,EHF,EDIFF,DIFF,DELT,TI1,SHA,SHB
      CONVGD=NEARLY.AND.DIFF.LE.ACURCY
      NEARLY=DIFF.LE.(ACURCY*3.0)
      IF(NEARLY)LOK=ITER
      IF(CONVGD) then
        write(10,130)
        WRITE(IW,130)
        GOTO 480
      endif
      CALL TIMIT(3)

C     ----- EXIT IN CASE OF TIME LIMIT

      IF((TIMLIM-TIM).GT.(DELT+DELT)) GOTO 470
      LSTOP=.TRUE.
      write(*,150)
      WRITE(IW,150)
      GO TO 480
470   IF(ITER.LT.MAXIT) GO TO 290
      write(*,160)
      WRITE(IW,160)
      LSTOP=.TRUE.

C     ----- PRINT ACTUAL VALUE OF ENERGY -----

480   WRITE(IW,190) ETOT,EHF
      write(10,190) ETOT,EHF
      IF(LSTOP) IREST=3
      IF(CONVGD) IREST=0

C      OUTPUT THE M.O.'S , EIGENVALUES AND DENSITY MATRICES

      IF(NPRINT.EQ.-5) GO TO 530
      WRITE(IW,60)
      WRITE(IW,80)
      DO 490 I=1,NA
490   WRITE(IW,90) I,QQ(J1+I-1)
      WRITE(IW,100)
      CALL VECOUT(QQ(I1),QQ(J1),NUM,NA,NDIM,IW,BFLAB)
      IF(NPRINT.NE.-1)WRITE(IW,110)
      IF(NPRINT.NE.-1)CALL FOUT(Q(1),NUM,IW)
      IF(NB.EQ.NA) GO TO 510
      WRITE(IW,70)
      WRITE(IW,80)
      NN = 0
      NLOWER = NA+1
      DO 500 I=NLOWER,NB
        NN = NN+1
500     WRITE(IW,90) NN,QQ(J1+I-1)
      WRITE(IW,100)
      CALL VECOUT(QQ(I1+NA*NUM),QQ(J1+NA),NUM,NB-NA,NDIM,IW,BFLAB)
      IF(NPRINT.NE.-1)WRITE(IW,110)
      IF(NPRINT.NE.-1)CALL FOUT(QQQ(I4),NUM,IW)
510   CONTINUE
      IF (NPRINT.EQ.-1) GO TO 530
      WRITE(IW,520)
520   FORMAT(/'      ----------------'/
     1        '      VIRTUAL ORBITALS'/
     2        '      ----------------')
      CALL VECOUT(QQ(I1+NB*NUM),QQ(J1+NB),NUM,NCOORB-NB,NDIM,IW,
     1            BFLAB)
530   CONTINUE

C     CALCULATE THE LAGRANGIAN FOR USE BY GRADIENT ROUTINES

      CALL FMOVE(QQ(I1),QQQ(I3),NUM*NCOORB)
      CALL READQ (Q(1),NX,IOSF(10),IFILS)
      CALL READQ (QQ(I1),NX,IOSF(9),IFILS)
      LEN=LENSEC(NX)
      LDA(ISEX(19))=NX
      LDA(ISEX(20))=NX
      CALL SECPUT(ISEX(19),19,LEN,ISECL)
      CALL SECPUT(ISEX(20),20,LEN,ISECK)
      CALL TIMIT(0)
      CALL REVIND
      CALL WRT3(QQ(I1),NX,ISECK,IFILD)
      IF(NPRINT.NE.-1.AND.NPRINT.NE.-5)WRITE(IW,575)
575   FORMAT(/1X,'LAGRANGIAN IN M.O. BASIS'/)
      N = -1
      IF (NA.EQ.0) GO TO 580
      DO 570 I = 1,NA
        DO 570 J = 1,I
          N = N+1
570       QQ(I2+N) = Q(N+1)*TWO
580   IF (NB.EQ.NA) GO TO 600
      DO 590 I = NA1,NB
        DO 590 J = 1,I
          N = N+1
          QQ(I2+N) = Q(N+1)*TWO
          IF (J.LE.NA) GO TO 590
          QQ(I2+N)=Q(N+1)-QQ(I1+N)
590       CONTINUE
600   IF(NPRINT.NE.-1.AND.NPRINT.NE.-5)CALL FOUT(QQ(I2),NB,IW)
      CALL DTR (QQ(I2),QQ(I1),QQQ(I3),QQQ(J2),NUM,NB,NDIM)
      DO 610 I = 1,NX
610     QQ(I1+I-1) = -QQ(I1+I-1)
      CALL WRT3 (QQ(I1),NX,ISECL,IFILD)
      IF(.NOT.LEX)GOTO 645
      CALL READQ(Q(1),NX,IOSF(7),IFILS)
      WRITE(IW,630)
630   FORMAT(/5X,'FOCK MATRIX ELEMENTS IN A.O. BASIS'/)
      WRITE(IW,640) (Q(I),I=1,NUM)
640   FORMAT(5X,8F15.3)
645   CALL SECSUM
*      TI=TIM0
      IF(NPRINT.EQ.-5) CALL TIMIT(3)
      IF(NPRINT.NE.-5)CALL TIMIT(1)
      CALL WRT3(VIBSIZ,LDA(ISEX(4)),ISEC4,IFILD)
      IF(LSTOP) GOTO 620
      CALL WHTPS
      RETURN

620   CALL CLENUP
      write(*,'(/a/)') ' MICROMOL INTERRUPTED IN SCFHS'
      write(IW,'(/a/)') ' MICROMOL INTERRUPTED IN SCFHS'
      stop
      END

C    *******************************************

C           FOCK MATRIX ASSEMBLY FOR OPEN-SHELL

      SUBROUTINE FFORM(AK,F,NA,NB,NORB,A,B,C,IA,IW,NPRINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AK(*),IA(*),F(*)

C     Input AK = K(DO)/2 = exchange operator for open shells
C           where DO is the open shell density
C           F = H0 + J(DT)-K(DT)/2 = Fock matrix using
C           total density matrix DT
C    A , B , and C are canonicalisation factors
C    A = B = C = 0 corresponds to default ATMOL choice
C    A = 0 ; B = -1 ; C = 0 corresponds to Roothaan double
C    Fock matrix approach
C    A = 2 ; B = 0 ; C = -2 corresponds to Roothaan single
C    Fock matrix approach

      NX=NORB*(NORB+1)/2
      NA1=NA+1
      NB1=NB+1

C     closed - closed block

      IF(NA.EQ.0)GOTO 3
      DO 2 I=1,NA
        DO 2 J=1,I
          IJ=IA(I)+J
2         AK(IJ)=AK(IJ)*A

C     open - open block

3     IF(NA1.GT.NB)GOTO 5
      DO 4 I=NA1,NB
        DO 4 J=NA1,I
          IJ=IA(I)+J
4         AK(IJ)=AK(IJ)*B

C     virtual - virtual block

5     IF(NB1.GT.NORB)GOTO 7
      DO 6 I=NB1,NORB
        DO 6 J=NB1,I
          IJ=IA(I)+J
6         AK(IJ)=AK(IJ)*C

C     closed - virtual block

7     IF(NA.EQ.0.OR.NB1.GT.NORB)GOTO 9
      DO 8 I=NB1,NORB
        DO 8 J=1,NA
          IJ=IA(I)+J
8         AK(IJ)=0.0d0

C     open - virtual block

9     IF(NA1.GT.NB.OR.NB1.GT.NORB)GOTO 11
      DO 10 I=NB1,NORB
        DO 10 J=NA1,NB
          IJ=IA(I)+J
10        AK(IJ)=-AK(IJ)

C     assemble complete Fock operator.

C       ---------------------------------------------
C       |   F + A*(K/2) |   F+(K/2)    |   F        |
C       ---------------------------------------------
C       |   F + (K/2)   | F + B*(K/2)  | F - (K/2)  |
C       ---------------------------------------------
C       |     F         | F - (K/2)    | F +C*(K/2) |
C       ---------------------------------------------

11    CONTINUE
      DO 12 I=1,NX
12      AK(I)=AK(I)+F(I)
      IF(NPRINT.NE.5)RETURN
      WRITE(IW,13)
*      CALL FOUT(AK,NORB,IW,NPRINT)
      CALL FOUT(AK,NORB,IW)
13    FORMAT(/1X,'COMPLETE FOCK MATRIX IN M.O. BASIS')
      RETURN
      END

C =========================================================
C
C      LEVEL SHIFTERS

      SUBROUTINE DMPER(F,DMP,ILOW,IUP,JLOW,JUP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(*)

      DO 10 I=ILOW,IUP
        DO 10 J=JLOW,JUP
          IJ=I*(I-1)/2+J
10        F(IJ)=F(IJ)*DMP

      RETURN
      END

C     ***************************************

      SUBROUTINE HSTHS(DC,FC,DO,FO,DT,IA)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      COMMON/BLKIN/XP(510),NINT,NJUNK
      DIMENSION IXP(1)
      EQUIVALENCE (IXP(1),XP(1))
      INTEGER P1,P2,P3
      CHARACTER*8 PARNAM,ATNAM,ANAM,BFLAB
      COMMON/INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(NATM1),C(3,NATM1),
     &             NVAR,NPAR,NVAL,NUQ,
     &             CFR,CFA,KONTYP(NATM1),P1(NATM1),P2(NATM1),P3(NATM1)
     &             ,KR(NATM1),KA(NATM1),KB(NATM1),PARVAL(NATM13)
     &             ,ZIN(NATM1),X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAM(NATM1)

      LOGICAL LADAPT,LSKIP
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
      DIMENSION IA(*)
      DIMENSION DC(*),DO(*),FC(*),FO(*),DT(*)
      PARAMETER (HALF=0.5D0, ZERO=0.0D0)

C     DC is the closed shell density matrix
C     DO is the open shell density matrix
C     DT is the total density matrix
C     FC is 2J(DT)-K(DT)
C     FO is the exchange operator for open shells, K(DO)

      DO 10 M=1,NX
        DT(M)=DC(M)+DO(M)
        FC(M)=ZERO
10      FO(M)=ZERO

C**********************
C    SEE COMMENTS IN HSTARO
C***********************

20    CALL GET(XP,NW,IFILM)
      IF(NW.EQ.0) GOTO 30
      CALL PROC2(FC,DT,FO,DO,IA)
      GOTO 20
30    DO 40 I=1,NX
40      FO(I)=-FO(I)
      DO 110 I=1,NX
        FC(I)=FC(I)*HALF
        FO(I)=FO(I)*HALF
110     CONTINUE
      RETURN
      END

C     ******************************************

      SUBROUTINE PROC2(A,P,B,Q,IKY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),Q(*),P(*),IKY(*)
      INTEGER*2 IJIN,KLIN
      real*4 g
      COMMON/BLKIN/G(510),IJIN(510),KLIN(510),MWORD,NJNK
*      COMMON/BLKIN/G(340),IJIN(340),KLIN(340),MWORD,NJNK
      EQUIVALENCE (GG1,GG),(IKYJ,JL),(IKYK,JK),(IKYI,IL)
      PARAMETER (ZERO=0.0D0)

      DO 120 IW=1,MWORD
        GG=G(IW)
        GG2=GG
        GG3=GG
        GG4=GG
        GG5=GG
        GG6=GG
        GG7=GG
        GG8=GG+GG
        GG9=GG8
        N=(IW+IW)+(IW+IW)
        I=IJIN(IW)/256
        K=KLIN(IW)/256
        J=IJIN(IW)-256*I
        L=KLIN(IW)-256*K
        IKYI=IKY(I)
        IJ=IKYI+J
        IK=IKYI+K
        IL=IKYI+L
        IKYK=IKY(K)
        KL=IKYK+L
        IKYJ=IKY(J)
        IF(I.NE.J)GO TO 10
        GG3=ZERO
        GG4=ZERO
        GG6=ZERO
        GG7=ZERO
        GG9=GG
10      IF(IJ.NE.KL)GO TO 20
        GG5=ZERO
        GG6=ZERO
        GG7=ZERO
        GG9=ZERO
20      IF(K.NE.L)GO TO 30
        GG2=ZERO
        GG4=ZERO
        GG7=ZERO
        GG8=GG
30      IF(I.EQ.K)GG1=GG1+GG5
        IF(J-K)50,40,90
40      GG3=GG3+GG6
        GO TO 60
50      GG3=GG6
60      JK=IKYK+J
        IF(J-L)80,70,100
70      GG4=GG4+GG7
        GO TO 100
80      GG4=GG7
        JL=IKY(L)+J
        GO TO 110
90      JK=IKYJ+K
100     JL=IKYJ+L
110     A(IJ)=(GG8+GG8)*P(KL)+A(IJ)
        A(KL)=(GG9+GG9)*P(IJ)+A(KL)
        A(IK)=A(IK)-GG1*P(JL)
        A(JL)=A(JL)-GG4*P(IK)
        A(IL)=A(IL)-GG2*P(JK)
        A(JK)=A(JK)-GG3*P(IL)
        B(IK)=B(IK)-GG1*Q(JL)
        B(JL)=B(JL)-GG4*Q(IK)
        B(IL)=B(IL)-GG2*Q(JK)
        B(JK)=B(JK)-GG3*Q(IL)
120     CONTINUE
      RETURN
      END
