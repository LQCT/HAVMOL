      SUBROUTINE JKDER
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      INTEGER*2 I,K,II,JJ,KK,LL,IT
      INTEGER P1,P2,P3
      LOGICAL DPRES,LAB,LABC,LABCD
      LOGICAL NOFORM,SKIP1,SKIP2,SKIP3,SKIP4,OUT
      CHARACTER*8 ANAM,ATNAME,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL ,NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13),  ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      COMMON/BLKIN/ZZZZZ(510),MWORD,MJNK
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

      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      COMMON/SHLNOS/Q4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1 MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2 NIJ,IJ,KL,IJKL
      CHARACTER*8 CLOSED
      CHARACTER*8 GROUP
      COMMON/SYMTRY/T(216),INVT(24),ISO(NSH1+1,24),NT
     $,NT2,XX(21),INDEX,NAXIS
      COMMON/SYMTRB/GROUP
      COMMON/TGRAD/DGOUT(9)
      COMMON/BIG/AA(100),R(100),X1(100),Y1(100),Z1(100),
     &   DD(400),IJDEN(100),IK(100),
     1  IJX(NATM13),IJY(NATM13),IJZ(NATM13),
     2  KLX(NATM13),KLY(NATM13),KLZ(NATM13),
     3  DIJ(100),DKL(100),IJGT(NATM13),KLGT(NATM13)
     4  ,NOFORM(1296),QBDENS(NBF1M),
     5  AKJNK(88)
      COMMON/BIGG/ABDENS(1297),QQX(1930),QQY(1930),QQZ(1930)
      COMMON/BIGK/QDENS(nbf1m)
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMMAX,TTOT
      COMMON/MAPPER/IA(NBF1+1)
      COMMON/USEFUL/NAT3,NDENS,NFOK,NBSQ,LENB
      COMMON/MISC/TOL(4),CUTOFF,NATOM(4),NPASS,OUT,SKIP1,SKIP2,
     1SKIP3,SKIP4,IJUNK(4)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,INC1,INC2,INC3,INC4,INC5,
     1  LENDD,IOFF,IXI,IYI,IZI,IABD,ISUB(6),IC7,IDA,IFOK

      DIMENSION M0(48),M1(48),M2(48),M3(48)
      DIMENSION QQ(1)
      EQUIVALENCE (QQ(1),ABDENS(1))

      PARAMETER (ZERO=0.0D0)
      DATA CLOSED/'CLOSED'/

      NAT3=NAT*3
      NBSQ=NUM*NUM
      LENB=LENSEC(NX)
      DPRES=.FALSE.
      DO 1 I=1,NSHELL
        IF(KTYPE(I).EQ.3)DPRES=.TRUE.
1       CONTINUE

      ID0=1
      IDA=2570  +  54*24
      IF(DPRES) IDA = 12970  +  193*30
      NDENS=1
      IF(SCFTYP.NE.CLOSED) NDENS=2
      IFOK=IDA+NX*NDENS
      NREQ=IDA+NDENS*NX
      WRITE(IW,123)NREQ
      write(10,123)NREQ
123   FORMAT(///' GRADIENTS: core used',I8,' words')

      KLOC(NSHELL+1)=NUM+1

C     ----- TWO ELECTRON CONTRIBUTION TO THE GRADIENT -----

      NDGOUT=9
      CALL TIMIT(3)
*      TIM0=TIM
      TIM1=TIM
      TIM2=TIM
      DTMAX=ZERO
      CALL DDEBUT(QDENS(1),QBDENS(1))
      MWORD=0
      IF(IST.GT.NSHELL) GO TO 190

C     ----- ISHELL -----

      DO 180 II=IST,NSHELL
        DO 10 IT=1,NT
          ID=ISO(II,IT)
          IF(ID.GT.II) GO TO 180
10        M0(IT)=ID
        ICENI=KATOM(II)

C     ----- JSHELL -----

        J0=JST
        DO 170 JJ=J0,II
          JST=1
          DO 30 IT=1,NT
            ID=M0(IT)
            JD=ISO(JJ,IT)
            IF(JD.GT.II) GO TO 170
            IF(ID.GE.JD) GO TO 20
            ND=ID
            ID=JD
            JD=ND
20          IF(ID.EQ.II.AND.JD.GT.JJ) GO TO 170
            M1(IT)=ID
30          M2(IT)=JD
          LAB=KATOM(JJ).EQ.ICENI

C     STORE INFORMATION ABOUT THE PAIR (IJ)

          CALL DSHELL(1,II,JJ,KK,LL)
          CALL DPRIM
          IF(NIJ.EQ.0)GOTO 170

C     ----- KSHELL -----

          K0=KST
          DO 140 KK=K0,II
            KST=1
            DO 40 IT=1,NT
              KD=ISO(KK,IT)
              IF(KD.GT.II) GO TO 140
40            M3(IT)=KD
            LABC=LAB.AND.KATOM(KK).EQ.ICENI

C     ----- LSHELL -----

            L0=LST
            MAXLL=KK
            IF(KK.EQ.II) MAXLL=JJ
            DO 130 LL=L0,MAXLL
              LST=1
              LABCD=LABC.AND.KATOM(LL).EQ.ICENI
              IF(LABCD)GOTO 130
              N4=0
              DO 70 IT=1,NT
                LD=ISO(LL,IT)
                IF(LD.GT.II) GO TO 130
                KD=M3(IT)
                IF(KD.GE.LD) GO TO 50
                ND=KD
                KD=LD
                LD=ND
50              ID=M1(IT)
                JD=M2(IT)
                IF(ID.NE.II.AND.KD.NE.II) GO TO 70
                IF(KD.LT.ID) GO TO 60
                IF(KD.EQ.ID.AND.LD.LE.JD) GO TO 60
                ND=ID
                ID=KD
                KD=ND
                ND=JD
                JD=LD
                LD=ND
60              IF(JD.LT.JJ) GO TO 70
                IF(JD.GT.JJ) GO TO 130
                IF(KD.LT.KK) GO TO 70
                IF(KD.GT.KK) GO TO 130
                IF(LD.LT.LL) GO TO 70
                IF(LD.GT.LL) GO TO 130
                N4=N4+1
70              CONTINUE

C     ----- CALCULATE Q4 FACTOR FOR THIS GROUP OF SHELLS -----

              Q4=DBLE(FLOAT(NT))/DBLE(FLOAT(N4))

C     ----- CHECK FOR REDUNDANT COMBINATIONS -----

              CALL REDUND(II,JJ,KK,LL,IW)
              IF(NPASS.EQ.0) GO TO 130

C     ----- INITIALIZE DGOUT TO ZERO -----

              DO 80 I=1,NDGOUT
80              DGOUT(I)=ZERO

              CALL DSHELL(2,II,JJ,KK,LL)

C     ----- FORM PRODUCTS OF DENSITY MATRIX ELEMENTS -----

              CALL DAB(QQ(1),QDENS(1),QBDENS(1),IA,II,JJ,
     1                 KK,LL,Q4,SCFTYP)

C     ----MESS ABOUT WITH THE DENSITY MATIX BY ELEMINATING
C     ZERO ELEMENTS

              DO 90 I=1,LENDD
                INC=I
                QQ(INC+LENDD)=QQ(INC)
90              NOFORM(I)=.FALSE.
              ABMAX=1.0D0
              ABMAX=0.0D0
              NN=0
              DO 100 I=1,IJ
                DO 100 K=1,KL
                  NN=NN+1
                  N=IJGT(I)+KLGT(K)
                  AB=DABS(QQ(LENDD+N))
                  IF(AB.LT.CUTOFF)NOFORM(NN)=.TRUE.
                  IF(AB.GT.ABMAX) ABMAX=AB
100               CONTINUE
              NN=0
              IJKL=0
              DO 120 I=1,IJ
                DO 120 K=1,KL
                  NN=NN+1
                  IF(NOFORM(NN))GOTO 120
                  N=IJGT(I)+KLGT(K)
                  IJKL=IJKL+1
                  QQ(IJKL)=QQ(LENDD+N)
120               CONTINUE
              IF(IJKL.EQ.0)GOTO 130
              CALL DGENRL(ABMAX)

C     ----- GENERATE ALL 4 PARTIAL CONTRIBUTIONS TO THE GRADIENT ----

              CALL FORMEG
130           CONTINUE
140         CONTINUE

C     ----- SAVE GRADIENT AND RESTART DATA -----

          CALL TIMIT(3)
          TIM3=TIM-TIM1
          TIM1=TIM
          IF(TIM3.GT.DTMAX) DTMAX=TIM3
          IF((TIMLIM-TIM).GT.(DTMAX*1.8D0+10.0D0)) GOTO 160
          CALL DFINAL(0,II,JJ)
160       CALL DFINAL(2,II,JJ)
170       CONTINUE
        CALL TIMIT(3)
        CALL DFINAL(2,II,II)
        DTIM=TIM-TIM2
        TIM2=TIM
        write(10,150) II,II,II,II,DTIM,TIM
        IF(NPRINT.NE.-5)WRITE(IW,150) II,II,II,II,DTIM,TIM
150   FORMAT(2X,'I,J,K,L =',4(I4,','),
     &       ' CYCLE TIME:',F11.2,' TOTAL TIME:',F11.2)
180     CONTINUE

C     ----- END OF 'SHELL' LOOPS -----

190   CONTINUE
      CALL SYMDE
      CALL DFINAL(1,II,JJ)
*      TI=TIM0
      CALL TIMIT(1)

      RETURN
      END

C=================================================
C        INITIALISING ROUTINE
C================================================
      SUBROUTINE DDEBUT(DA,DB)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      INTEGER P1,P2,P3
      LOGICAL OUT,SKIP(4),IANDJ,KANDL,SAME,NORM
      COMMON/MISC/TOL1,TOL2,TOL3,TOL4,CUTOFF,NATOM(4),NPASS,OUT
     1,SKIP,IANDJ,KANDL,SAME,NORM
      COMMON/INCRS/KLN2(22),IDA,IFOK

      LOGICAL LADAPT,LSKIP
      CHARACTER*8 CLOSED
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
      COMMON/GRAD/DE(NATM13)
      COMMON/MAPPER/IA(NBF1+1)

      COMMON/USEFUL/NAT3,NDENS,NFOK,NBSQ,LENB
      CHARACTER*8 PARNAM,BFLAB,ATNAME,ANAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     4    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      DIMENSION DA(*),DB(*)
      PARAMETER (ONE=1.0D0, TEN=1.0D1, E=2.30258D0)
      DATA CLOSED/'CLOSED'/

      ISECD1=0
      ISECD2=0
      IBLOK=0
      write(10,10)
      WRITE(IW,10)
10    FORMAT(/,10X,22(1H-),/,10X,'GRADIENT OF THE ENERGY',
     1       /,10X,22(1H-)/)

      CALL SECGET(ISEX(7),7,ISECD1)
      CALL READQ(DA,LDA(ISEX(7)),ISECD1,IFILD)
      IF(SCFTYP.EQ.CLOSED) GO TO 35
      CALL SECGET(ISEX(10),10,ISECD2)
      CALL READQ(DB,LDA(ISEX(10)),ISECD2,IFILD)
      DO 20 I=1,NX
      DUMA=DA(I)
      DUMB=DB(I)
      DA(I)=DUMA+DUMB
20    CONTINUE
35    OUT=NPRINT.EQ.-4
      CUTOFF=ONE/(TEN**ICUT)
      TOL1=E*(ITOL+1)
      TOL2=E*(ITOL-1)
      TOL3=ONE/TEN**(ITOL+1)
      TOL4=ONE/TEN**(ITOL-1)
      NORM=NORMF.NE.1.OR.NORMP.NE.1
      CALL SECGET(ISEX(14),14,IBLOK)
      CALL READQ(DE,LDA(ISEX(14)),IBLOK,IFILD)
      RETURN
      END

C==========================================
C     CLOSE DOWN ROUTINE
C========================================
      SUBROUTINE DFINAL(INDEX,II,JJ)

      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      INTEGER P1,P2,P3
      INTEGER*2 II,JJ
      CHARACTER*8 DNAM
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
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMMAX,TTOT

      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     2             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      CHARACTER*8 BFLAB,PARNAM,ATNAME,ANAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     4    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      COMMON/GRAD/DE(NATM13)
      COMMON/FUNCT/EHF,EG(NATM13)
      DIMENSION F(3,NATM1),DNAM(3)
      EQUIVALENCE (EG(1),F(1,1))
      DATA DNAM /'E''X','E''Y','E''Z'/

      IBLOK=0
      IF(INDEX.NE.2) THEN
        IF(INDEX.EQ.1) GO TO 40

C     ----- GET RESTART DATA -----
 
        IF(JJ.EQ.NSHELL) RETURN
        IST=II 
        JST=JJ+1
        IF (JJ.EQ.II) IST=II+1
        IF (JJ.EQ.II) JST=1
        KST=1
        LST=1
        NREC=1
        IREST=5
        INTLOC=1
        write(*,10)
        WRITE(IW,10)
10    FORMAT(//T21,'*** INSUFFICIENT TIME TO CONTINUE ***'
     &//T16,'This job can be restarted to complete evaluation',
     &/T20,'of two-electron contribution to gradient')
        write(10,20) IST,JST,KST,LST
        WRITE(IW,20) IST,JST,KST,LST
20    FORMAT(//T14,'IST, JST, KST, LST ',4I8)
        CALL SECGET(ISEX(14),14,IBLOK)
        CALL WRT3(DE,LDA(ISEX(14)),IBLOK,IFILD)
        CALL SECGET(ISEX(4),4,IBLOK)
        CALL WRT3(VIBSIZ,LDA(ISEX(4)),IBLOK,IFILD)
        CALL SECGET(ISEX(24),24,IBLOK)
        CALL WRT3C(TITLE,LDA(ISEX(24)),IBLOK,IFILD)
        CALL CLENUP
        write(*,'(/a/)') ' MICROMOL INTERRUPTED IN DFINAL'
        write(IW,'(/a/)') ' MICROMOL INTERRUPTED IN DFINAL'
        stop
      ENDIF
      IF(JJ.EQ.NSHELL) RETURN
      ISTR=IST
      JSTR=JST
      KSTR=KST
      LSTR=LST
      IST=II
      JST=JJ+1
      IF (JJ.EQ.II) IST=II+1
      IF(IST.GT.NSHELL)RETURN
      IF (JJ.EQ.II) JST=1
      KST=1
      LST=1
      IREST=5
      CALL SECGET(ISEX(4),4,IBLOK)
      CALL WRT3(VIBSIZ,LDA(ISEX(4)),IBLOK,IFILD)
      CALL SECGET(ISEX(24),24,IBLOK)
      CALL WRT3C(TITLE,LDA(ISEX(24)),IBLOK,IFILD)
      CALL SECGET(ISEX(14),14,IBLOK)
      CALL WRT3(DE,LDA(ISEX(14)),IBLOK,IFILD)
      IST=ISTR
      JST=JSTR
      KST=KSTR
      LST=LSTR
      RETURN

40    NCOORD=3*NAT
      DO 50 I=1,NCOORD
50      EG(I)=DE(I)

C     ----- SAVE 1E+2E-GRADIENT -----

      IREST=0
      IST=1
      JST=1
      KST=1
      LST=1
      NREC=1
      INTLOC=1
      CALL SECGET(ISEX(4),4,IBLOK)
      CALL WRT3(VIBSIZ,LDA(ISEX(4)),IBLOK,IFILD)
      CALL SECGET(ISEX(24),24,IBLOK)
      CALL WRT3C(TITLE,LDA(ISEX(24)),IBLOK,IFILD)
      CALL SECGET(ISEX(14),14,IBLOK)
      CALL WRT3(EG,LDA(ISEX(14)),IBLOK,IFILD)
      MAX=0
60    MIN=MAX+1
      MAX=MAX+5
      IF(MAX.GT.NAT) MAX=NAT
      write(10,90)
      WRITE(IW,90)
      write(10,100) (I,I=MIN,MAX)
      WRITE(IW,100) (I,I=MIN,MAX)
      write(10,90)
      WRITE(IW,90)
      DO 70 N=1,3
        write(10,110) DNAM(N),(F(N,I),I=MIN,MAX)
70      WRITE(IW,110) DNAM(N),(F(N,I),I=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 60
      IF(NPRINT.NE.-5) CALL WHTPS
      WRITE(IW,80)IREST
80    FORMAT(/T12,'End of calculation of the energy gradient, IREST =',
     &I4,/)
90    FORMAT(/)
100   FORMAT(5X,'ATOM',5(6X,I2,7X))
110   FORMAT(7X,A3,5D15.7)
      CALL TIMIT(3)
      CALL REVIND
      RETURN
      END

C==============================================
C    FIND WHICH DERIVATIVES CAN BE OMITTED
C=============================================
      SUBROUTINE REDUND(II,JJ,KK,LL,IW)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      INTEGER*2 II,JJ,KK,LL
      LOGICAL OUT
      LOGICAL SKIP1,SKIP2,SKIP3,SKIP4
      LOGICAL INEJ,INEK,INEL,JNEK,JNEL,KNEL
      COMMON/MISC/TOL(4),CUTOFF,NATOM(4),NPASS,OUT,SKIP1,SKIP2,
     1SKIP3,SKIP4,MJNK(4)

      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      DIMENSION LLL(4)
      EQUIVALENCE (LLL(1),LIT)

      SKIP1=.TRUE.
      SKIP2=.TRUE.
      SKIP3=.TRUE.
      SKIP4=.TRUE.
      NPASS=0
      DO 10 I=1,4
10      NATOM(I)=0
      LIT=KTYPE(II)
      LJT=KTYPE(JJ)
      LKT=KTYPE(KK)
      LLT=KTYPE(LL)
      IAT=KATOM(II)
      JAT=KATOM(JJ)
      KAT=KATOM(KK)
      LAT=KATOM(LL)
      INEJ=IAT.NE.JAT
      INEK=IAT.NE.KAT
      INEL=IAT.NE.LAT
      JNEK=JAT.NE.KAT
      JNEL=JAT.NE.LAT
      KNEL=KAT.NE.LAT
      IF(INEJ) GO TO 50
      IF(INEK) GO TO 30
      IF(INEL) GO TO 20

C      IAT=JAT=KAT=LAT   ----  NO DERIVATIVES

      GO TO 230

C     IAT=JAT=KAT   DERIVATIVE (IJ/KL')

20    SKIP4=.FALSE.
      NATOM(1)=LAT
      NATOM(2)=IAT
      NPASS=1
      GO TO 230
30    IF(KNEL) GO TO 40

C     IAT=JAT  ,  KAT=LAT   DIFFERENTIATE ONE PAIR

      N1=LIT*LJT*(LKT+1)*(LLT+1)
      N2=(LIT+1)*(LJT+1)*LKT*LLT
      IF(N2.LT.N1)GOTO 160

C     IAT=JAT   DERIVATIVES (IJ/K'L) AND (IJ/KL')

40    SKIP3=.FALSE.
      SKIP4=.FALSE.
      NATOM(1)=KAT
      NATOM(2)=LAT
      NATOM(3)=IAT
      NPASS=2
      GO TO 230
50    IF(INEK) GO TO 100
      IF(INEL) GO TO 70

C     IAT=KAT=LAT   DERIVATIVE (IJ'/KL)

60    SKIP2=.FALSE.
      NATOM(1)=JAT
      NATOM(2)=IAT
      NPASS=1
      GO TO 230

C      IAT=KAT    JAT=LAT

70    IF(JNEL)GOTO 90
      IF(II.NE.KK.OR.JJ.NE.LL)GOTO 80
      IF(LJT.GT.LIT)GOTO 60
      GOTO 110
80    N1=(LIT+1)*(LKT+1)*LJT*LLT
      N2=LIT*LKT*(LJT+1)*(LLT+1)
      IF(N1.LT.N2)GOTO 140

C     IAT=KAT   DERIVATIVES (IJ'/KL) AND (IJ/KL')

90    SKIP2=.FALSE.
      SKIP4=.FALSE.
      NATOM(1)=JAT
      NATOM(2)=LAT
      NATOM(3)=IAT
      NPASS=2
      GO TO 230
100   IF(JNEK) GO TO 130
      IF(JNEL) GO TO 120

C     JAT=KAT=LAT    (I'J/KL)

110   SKIP1=.FALSE.
      NATOM(1)=IAT
      NATOM(2)=JAT
      NPASS=1
      GO TO 230

C     ----- JAT = KAT # IAT # LAT -----

120   SKIP1=.FALSE.
      SKIP4=.FALSE.
      NATOM(1)=IAT
      NATOM(2)=LAT
      NATOM(3)=JAT
      NPASS=2
      GO TO 230
130   IF(JNEL) GO TO 150

C      JAT=LAT    DERIVATIVES (I'J/KL) AND (IJ/K'L)

140   SKIP1=.FALSE.
      SKIP3=.FALSE.
      NATOM(1)=IAT
      NATOM(2)=KAT
      NATOM(3)=JAT
      NPASS=2
      GO TO 230
150   IF(KNEL) GO TO 170

C     KAT=LAT   DERIVATIVES (I'J/KL) AND (IJ'/KL)

160   SKIP1=.FALSE.
      SKIP2=.FALSE.
      NATOM(1)=IAT
      NATOM(2)=JAT
      NATOM(3)=KAT
      NPASS=2
      GO TO 230

C     IAT # JAT # KAT # LAT  -- OMIT ONE CENTRE

170   MIN=LIT
      IMIN=1
      DO 180 IPER=2,4
        IF(LLL(IPER).GE.MIN)GOTO 180
        MIN=LLL(IPER)
        IMIN=IPER
180     CONTINUE
      GOTO (190,200,210,220),IMIN
190   NATOM(1)=JAT
      NATOM(2)=KAT
      NATOM(3)=LAT
      NATOM(4)=IAT
      NPASS=3
      SKIP2=.FALSE.
      SKIP3=.FALSE.
      SKIP4=.FALSE.
      GOTO 230
200   NATOM(1)=IAT
      NATOM(2)=KAT
      NATOM(3)=LAT
      NATOM(4)=JAT
      NPASS=3
      SKIP1=.FALSE.
      SKIP3=.FALSE.
      SKIP4=.FALSE.
      GOTO 230
210   NATOM(1)=IAT
      NATOM(2)=JAT
      NATOM(3)=LAT
      NATOM(4)=KAT
      NPASS=3
      SKIP1=.FALSE.
      SKIP2=.FALSE.
      SKIP4=.FALSE.
      GOTO 230
220   SKIP1=.FALSE.
      SKIP2=.FALSE.
      SKIP3=.FALSE.
      NATOM(1)=IAT
      NATOM(2)=JAT
      NATOM(3)=KAT
      NATOM(4)=LAT
      NPASS=3

C     -----

230   IF(.NOT.OUT) RETURN
      WRITE(IW,240) II,JJ,KK,LL,SKIP1,SKIP2,SKIP3,SKIP4,NPASS,
     1(NATOM(I),I=1,4)
240   FORMAT(' II,JJ,KK,LL =',4I3,', SKIP1,2,3,4 =',4L3,
     1 ', NPASS =',I2,' CENTERS =',4I5)
      RETURN
      END

C====================================================
C    TWO PARTICLE DENSITY MATRIX FOR SCF CASES
C=====================================================
      SUBROUTINE DAB(ABDENS,DA,DB,IA,II,JJ,KK,LL,Q4,SCFTYP)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      INTEGER*2 II,JJ,KK,LL,I,J,K,L
      CHARACTER*8 CLOSED,SCFTYP
      INTEGER P1,P2,P3
      LOGICAL IJEQ,KLEQ
      DIMENSION ABDENS(*),DA(*),DB(*),IA(*)

      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),         
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      LOGICAL OPEN
      CHARACTER*8 ATNAME,BFLAB,PARNAM,ANAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,INC,INCL,INCK,INCJ,INCI,NJNK(15)
      PARAMETER (HALF=0.5D0, FOUR=4.0D0)
      DATA CLOSED/'CLOSED'/

      OPEN=SCFTYP.NE.CLOSED
      MINI=KMIN(II)
      MINJ=KMIN(JJ)
      MINK=KMIN(KK)
      MINL=KMIN(LL)
      MAXI=KMAX(II)
      MAXJ=KMAX(JJ)
      MAXK=KMAX(KK)
      MAXL=KMAX(LL)
      LOCI=KLOC(II)-MINI
      LOCJ=KLOC(JJ)-MINJ
      LOCK=KLOC(KK)-MINK
      LOCL=KLOC(LL)-MINL
      NI=1
      DO 40 I=MINI,MAXI
        NJ=NI
        I1=LOCI+I
        II1=IA(I1)
        DO 30 J=MINJ,MAXJ
          NK=NJ
          J1=LOCJ+J
          JJ1=IA(J1)
          MIJ=II1+J1
          IF(J1.GT.I1)MIJ=JJ1+I1
          DIJ=DA(MIJ)*FOUR
          IJEQ=I1.EQ.J1
          DO 20 K=MINK,MAXK
            NL=NK
            K1=LOCK+K
            KK1=IA(K1)
            MIK=II1+K1
            IF(K1.GT.I1)MIK=KK1+I1
            MJK=JJ1+K1
            IF(K1.GT.J1)MJK=KK1+J1
            DIK=DA(MIK)
            DJK=DA(MJK)
            DO 10 L=MINL,MAXL
              NN=NL
              L1=LOCL+L
              LL1=IA(L1)
              MKL=KK1+L1
              IF(L1.GT.K1)MKL=LL1+K1
              MJL=JJ1+L1
              IF(L1.GT.J1)MJL=LL1+J1
              MIL=II1+L1
              IF(L1.GT.I1)MIL=LL1+I1
              KLEQ=K1.EQ.L1
              DFAC=DIJ*DA(MKL)-DIK*DA(MJL)-DJK*DA(MIL)
              IF(OPEN) DFAC=DFAC-DB(MIK)*DB(MJL)-DB(MJK)*DB(MIL)
              IF(IJEQ)DFAC=DFAC*HALF
              IF(KLEQ)DFAC=DFAC*HALF
              ABDENS(NN)=DFAC*Q4
10            NL=NL+INCL
20          NK=NK+INCK
30        NJ=NJ+INCJ
40      NI=NI+INCI

      RETURN
      END

C ========================================
C   TRANSLATIONAL INVARIANCE OF DERIVATIVE
C   USED TO GET MISSING TERMS
C===========================================
      SUBROUTINE FORMEG
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      LOGICAL OUT
      LOGICAL SKIP1,SKIP2,SKIP3,SKIP4
      COMMON/MISC/TOL(4),CUTOFF,NATOM(4),NPASS,OUT,SKIP1,SKIP2,
     &SKIP3,SKIP4,NJNK(4)
      COMMON/GRAD/DE(3,NATM1)
      COMMON/TGRAD/DGOUT(3,3)
      PARAMETER (ZERO=0.0D0)

      DUMX=ZERO
      DUMY=ZERO
      DUMZ=ZERO
      DO 100 IPASS=1,NPASS
        IAT=NATOM(IPASS)
        DUM=DGOUT(1,IPASS)
        DUMX=DUMX+DUM
        DE(1,IAT)=DE(1,IAT)+DUM
        DUM=DGOUT(2,IPASS)
        DUMY=DUMY+DUM
        DE(2,IAT)=DE(2,IAT)+DUM
        DUM=DGOUT(3,IPASS)
        DUMZ=DUMZ+DUM
        DE(3,IAT)=DE(3,IAT)+DUM
  100   CONTINUE
      IAT=NATOM(NPASS+1)
      DE(1,IAT)=DE(1,IAT)-DUMX
      DE(2,IAT)=DE(2,IAT)-DUMY
      DE(3,IAT)=DE(3,IAT)-DUMZ

      RETURN
      END

C  DSHELL     ----   1 NOV 1983
C  AUTODBL
C============================================
C     INFORMATION ABOUT SHELLS
C     EQUIVALENT OF SUBROUINE SHELLS IN INTEGRALS
C==================================================
      SUBROUTINE DSHELL(NELEC,ISH,JSH,KSH,LSH)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      INTEGER*2 I,J,K,L,ISH,JSH,KSH,LSH
      INTEGER P1,P2,P3
      LOGICAL OUT,SKIP(4),IANDJ,KANDL,SAME,NORM
      COMMON/MISC/TOL1,TOL2,TOL3,TOL4,CUTOFF,NATOM(4),NPASS,OUT
     1,SKIP,IANDJ,KANDL,SAME,NORM
      CHARACTER*8 ATNAME,PARNAM,ANAM,BFLAB
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)

      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      COMMON/RTWT/XX,U(NATM1),W(NATM1),NROOTS,JUNK,ZJUNK(NSH1*2)
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1 MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2 NIJ,IJ,KL,IJKL
      LOGICAL LDUM4
      COMMON/BIG/AAA(900),IJDEN(100),IK(100),
     1           IJX(NATM13),IJY(NATM13),IJZ(NATM13),
     2           KLX(NATM13),KLY(NATM13),KLZ(NATM13),
     3           DIJ(100),DKL(100),IJGT(NATM13),KLGT(NATM13),
     4           LDUM4(1296),DUM4(NBF1M),
     5           DUM5(88)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,INC,INCL,INCK,INCJ,INCI,
     &   LENDD,IOFF,IXI,IYI,IZI,IABD,ISUB,IJNK(8)
      COMMON/SHLINF/GA(10),CSA(10),CPA(10),CDA(10),
     1              GB(10),CSB(10),CPB(10),CDB(10),
     2              GC(10),CSC(10),CPC(10),CDC(10),
     3              GD(10),CSD(10),CPD(10),CDD(10),
     4              AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD,
     5              NGA,NGB,NGC,NGD,EXIJ(100),RSMALL

      IF(NELEC.EQ.2) GO TO 40
      IANDJ=ISH.EQ.JSH
      I=KATOM(ISH)
      AX=C(1,I)
      AY=C(2,I)
      AZ=C(3,I)
      I1=KSTART(ISH)
      I2=I1+KNG(ISH)-1
      LIT=KTYPE(ISH)
      MINI=KMIN(ISH)
      MAXI=KMAX(ISH)
      LOCI=KLOC(ISH)-MINI
      NGA=0
      DO 10 I=I1,I2
        NGA=NGA+1
        GA(NGA)=EX(I)
        CSA(NGA)=CS(I)
        CPA(NGA)=CP(I)
10      CDA(NGA)=CD(I)
      J=KATOM(JSH)
      BX=C(1,J)
      BY=C(2,J)
      BZ=C(3,J)
      J1=KSTART(JSH)
      J2=J1+KNG(JSH)-1
      LJT=KTYPE(JSH)
      MINJ=KMIN(JSH)
      MAXJ=KMAX(JSH)
      LOCJ=KLOC(JSH)-MINJ
      NGB=0
      DO 20 J=J1,J2
        NGB=NGB+1
        GB(NGB)=EX(J)
        CSB(NGB)=CS(J)
        CPB(NGB)=CP(J)
20      CDB(NGB)=CD(J)
      RAB=((AX-BX)**2+(AY-BY)**2+(AZ-BZ)**2)
      RETURN

40    KANDL=KSH.EQ.LSH
      SAME=ISH.EQ.KSH.AND.JSH.EQ.LSH
      K=KATOM(KSH)
      CX=C(1,K)
      CY=C(2,K)
      CZ=C(3,K)
      K1=KSTART(KSH)
      K2=K1+KNG(KSH)-1
      LKT=KTYPE(KSH)
      MINK=KMIN(KSH)
      MAXK=KMAX(KSH)
      LOCK=KLOC(KSH)-MINK
      NGC=0
      DO 50 K=K1,K2
        NGC=NGC+1
        GC(NGC)=EX(K)
        CSC(NGC)=CS(K)
        CPC(NGC)=CP(K)
50      CDC(NGC)=CD(K)
      L=KATOM(LSH)
      DX=C(1,L)
      DY=C(2,L)
      DZ=C(3,L)
      L1=KSTART(LSH)
      L2=L1+KNG(LSH)-1
      LLT=KTYPE(LSH)
      MINL=KMIN(LSH)
      MAXL=KMAX(LSH)
      LOCL=KLOC(LSH)-MINL
      NGD=0
      DO 60 L=L1,L2
        NGD=NGD+1
        GD(NGD)=EX(L)
        CSD(NGD)=CS(L)
        CPD(NGD)=CP(L)
60      CDD(NGD)=CD(L)
      NROOTS=(LIT+LJT+LKT+LLT-1)/2
      RCD=((CX-DX)**2+(CY-DY)**2+(CZ-DZ)**2)

C     DETERMINE VARIOUS OFFSETS AND INDEXING ARRAYS

      INCL=1
      INCK=INCL*(MAXL-MINL+1)
      INCJ=INCK*(MAXK-MINK+1)
      INCI=INCJ*(MAXJ-MINJ+1)
      LENDD=INCI*(MAXI-MINI+1)
      IJ=0
      MAX=MAXJ
      DO 30 I=MINI,MAXI
        IF(IANDJ)MAX=I
        DO 30 J=MINJ,MAX
          IJ=IJ+1
          IJGT(IJ)=INCI*(I-MINI)+INCJ*(J-MINJ)+1
30        CONTINUE
      KL=0
      MAX=MAXL
      DO 70 K=MINK,MAXK
        IF(KANDL)MAX=K
        DO 70 L=MINL,MAX
          KL=KL+1
          KLGT(KL)=INCK*(K-MINK)+INCL*(L-MINL)
70        CONTINUE
      IJKL=IJ*KL
      IXI=LENDD+1
      IYI=IXI+LENDD
      IZI=IYI+LENDD
      IABD=LENINT(LENDD*4)
      IOFF=0
      ISUB=LENDD+1
      RETURN
      END

C=============================================
C      STORE INFORMATION REGARDING PAIR IJ
C    EQUIVALENT OF IJPRIM FOR INTEGRAL ROUTINES
C=============================================
      SUBROUTINE DPRIM
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      LOGICAL SKIP(4)
      LOGICAL IANDJ,KANDL,SAME,OUT,NORM
      COMMON/MISC/TOL1,TOL2,TOL3,TOL4,CUTOFF,NATOM(4),NPASS,OUT,SKIP,
     1IANDJ,KANDL,SAME,NORM
      COMMON/SHLINF/AG(10),CSA(10),CPA(10),CDA(10),
     1              BG(10),CSB(10),CPB(10),CDB(10),
     2              CG(10),CSC(10),CPC(10),CDC(10),
     3              DG(10),CSD(10),CPD(10),CDD(10),
     4              XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,
     5              NGA,NGB,NGC,NGD,EXIJ(100),RSMALL
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1 MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2 NIJ,IJ,KL,IJKL
      LOGICAL LDUM4
      COMMON/BIG/A(100),R(100),X1(100),Y1(100),Z1(100),
     &           DD(400),IJDEN(100),IDUM0(100),
     1           IDUM1(3,NATM13),
     2           IDUM2(3,NATM13),
     3           DUM3(200),IDUM3(2,NATM13),
     4           LDUM4(1296),DUM4(NBF1M),
     5           DUM5(88)
      PARAMETER (ONE=1.0D0)

      MAX=MAXJ
      N=0
      NN=0
      DO 50 I=MINI,MAXI
        GO TO (10,10,20,20,10,20,20,20,20,20,
     &         10,20,20,20,20,20,20,20,20,20),I
10        NM=NN
20        NN=NM
        IF(IANDJ) MAX=I
        DO 50 J=MINJ,MAX
          GO TO (30,30,40,40,30,40,40,40,40,40,
     &           30,40,40,40,40,40,40,40,40,40),J
30          NN=NN+1
40          N=N+1
50        IJDEN(N)=NN

C     ----- I PRIMITIVE

      NIJ=0
      JBMAX=NGB
      DO 190 IA=1,NGA
        AI=AG(IA)
        ARRI=AI*RRI
        AXI=AI*XI
        AYI=AI*YI
        AZI=AI*ZI
        CSI=CSA(IA)
        CPI=CPA(IA)
        CDI=CDA(IA)

C      ----- J PRIMITIVE

        IF(IANDJ) JBMAX=IA
        DO 180 JB=1,JBMAX
          AJ=BG(JB)
          AA=AI+AJ
          AINV=ONE/AA
          DUM=AJ*ARRI*AINV
          CSJ=CSB(JB)*AINV
          CPJ=CPB(JB)*AINV
          CDJ=CDB(JB)*AINV
          NM=(NIJ+NIJ)+(NIJ+NIJ)
          NN=NM
          NIJ=NIJ+1
          EXIJ(NIJ)=DEXP(-DUM)
          R(NIJ)=DUM
          A(NIJ)=AA
          X1(NIJ)=(AXI+AJ*XJ)*AINV
          Y1(NIJ)=(AYI+AJ*YJ)*AINV
          Z1(NIJ)=(AZI+AJ*ZJ)*AINV

C      ----- DENSITY FACTOR

          DO 150 I=MINI,MAXI
            IF(IANDJ)MAX=I
            GO TO(60,70,150,150,80,150,150,150,150,150,
     &            80,150,150,150,150,150,150,150,150,150),I
60            DUM1=CSI
            GO TO 90
70            DUM1=CPI
            GO TO 90
80            DUM1=CDI
90          CONTINUE
            DO 140 J=MINJ,MAX
              GO TO(100,110,140,140,120,140,140,140,140,140,
     &              120,140,140,140,140,140,140,140,140,140),J
100             DUM2=DUM1*CSJ
              GO TO 130
110             DUM2=DUM1*CPJ
              GO TO 130
120             DUM2=DUM1*CDJ
130           NN=NN+1
              DD(NN)=DUM2
140           CONTINUE
150         CONTINUE
          IF(.NOT.IANDJ) GO TO 180
          IF(IA.EQ.JB) GO TO 180
          NMP1=NM+1
          NMP2=NM+2
          NMP3=NM+3
          GO TO (170,160,170,170),LIT
160       IF(MINI.EQ.2) GO TO 170
          DD(NMP2)=DD(NMP2)+CSI*CPJ
          DD(NMP3)=DD(NMP3)+DD(NMP3)
170       DD(NMP1)=DD(NMP1)+DD(NMP1)
180       CONTINUE
190     CONTINUE
      IF(NIJ.EQ.0)RETURN
      RSMALL=R(1)
      DO 200 N=1,NIJ
        IF(RSMALL.GT.R(N))RSMALL=R(N)
200     CONTINUE
      IF(RSMALL.GE.TOL1)NIJ=0
      RETURN
      END

C===========================================
C    CONTAINS LOOPS OVER PRIMITIVES
C    HEAVILY MODIFIED VERSION OF ORIGINAL HONDO
C==========================================
      SUBROUTINE DGENRL(ABMAX)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      INTEGER*2 I,K,M, IG,JG,KG,LG
      LOGICAL SKIP,DTRUIJ,DTRUKL
      LOGICAL DOUBLE,IANDJ,KANDL,SAME,OUT,NORM
      LOGICAL NOFORM
      COMMON/SHLINF/AG(10),CSA(10),CPA(10),CDA(10),
     1              BG(10),CSB(10),CPB(10),CDB(10),
     2              CG(10),CSC(10),CPC(10),CDC(10),
     3              DG(10),CSD(10),CPD(10),CDD(10),
     4              XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,
     5              NGA,NGB,NGC,NGD,EXIJ(100),RSMALL
      COMMON/BIG/AA(100),R(100),X1(100),Y1(100),Z1(100),
     &   DD(400),IJDEN(100),IK(100),
     1  IJX(NATM13),IJY(NATM13),IJZ(NATM13),
     2  KLX(NATM13),KLY(NATM13),KLZ(NATM13),
     3  DIJ(100),DKL(100),IJGT(NATM13),KLGT(NATM13),
     4  NOFORM(1296),QBDENS(NBF1M),
     5  ABJNK(88)
      COMMON/BIGG/ABDENS(1297),QQX(1930),QQY(1930),QQZ(1930)
      COMMON/SHLNOS/Q4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1 MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2 NIJ,IJ,KL,IJKL
      COMMON/RTWT/XX,U(NATM1),W(NATM1),NROOTS,JUNK,ZJUNK(NSH1*2)
      COMMON/SETINT/BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00
     1,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL
     2,IN(NATM1),KN(NATM1),NI,NJ,NK,NL,NMAX,MMAX,IJ1,IJ2,KL1,KL2
      COMMON/MISC/TOL1,TOL2,TOL3,TOL4,CUTOFF,NATOM(4),NPASS,OUT,SKIP(4),
     1IANDJ,KANDL,SAME,NORM
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,INC1,INC2,INC3,INC4,INC5,
     1  LENDD,IOFF,IXI,IYI,IZI,IABD,IC1,IC2,IC3,IC4,IC5,IC6,IC7,IDA,IJNK
      COMMON/TGRAD/DGOUT(9)
      DIMENSION IQQ(1),QQ(1)
      EQUIVALENCE (IQQ(1),NOFORM(1)),(QQ(1),ABDENS(1))
      DIMENSION IN1(NATM1)
      PARAMETER (ONE=1.0D0, PT5=0.5D0, PI252=34.986836655250D0)

      IF(IJKL.EQ.0)RETURN
      NI=LIT
      IF(SKIP(1))NI=LIT-1
      NJ=LJT
      IF(SKIP(2))NJ=LJT-1
      NK=LKT
      IF(SKIP(3))NK=LKT-1
      NL=LLT
      IF(SKIP(4))NL=LLT-1
      KLN2=1
      KLN1=NL+1
      IJN2=KLN1*(NK+1)
      IJN1=IJN2*(NJ+1)
      INC1=IJN1*(NI+1)
      IF(MOD(INC1,4).EQ.0)INC1=INC1+1
      IF(NI.GE.NJ) THEN
        IJ1 = IJN1
        IJ2 = IJN2
        XC = XI
        YC = YI
        ZC = ZI
        DXIJ = XI-XJ
        DYIJ = YI-YJ
        DZIJ = ZI-ZJ
      ELSE
        IS=NI
        NI=NJ
        NJ=IS
        IJ1 = IJN2
        IJ2 = IJN1  
        XC = XJ
        YC = YJ
        ZC = ZJ
        DXIJ = XJ-XI
        DYIJ = YJ-YI
        DZIJ = ZJ-ZI
      ENDIF
      IF(NK.GE.NL) THEN
        XD=XK
        YD = YK
        ZD = ZK
        DXKL = XK-XL
        DYKL = YK-YL
        DZKL = ZK-ZL
        KL1=KLN1
        KL2=KLN2
      ELSE
        IS=NL
        NL=NK
        NK=IS
        KL1=KLN2
        KL2=KLN1
        XD=XL
        YD=YL
        ZD=ZL
        DXKL=XL-XK
        DYKL=YL-YK
        DZKL=ZL-ZK
      ENDIF
      NMAX = NI+NJ
      MMAX = NK+NL
      MAX = NMAX+1
      DO 50 I = 1,MAX
        N = I-1
        IF (N .LE. NI) IN1(I) = IJ1*N+1
50      IF (N .GT. NI) IN1(I) = IJ1*NI+IJ2*(N-NI)+1
      MAX=MMAX+1
      DO 60 K=1,MAX
        N=K-1
        IF(N.LE.NK) KN(K)=KL1*N
60      IF(N.GT.NK) KN(K)=KL1*NK+KL2*(N-NK)

C     INDEXING

      CALL INDXA(IJX,IJY,IJZ,IJ,MINI,MAXI,MINJ,MAXJ,IANDJ,IJN1,IJN2,1)
      CALL INDXA(KLX,KLY,KLZ,KL,MINK,MAXK,MINL,MAXL,KANDL,KLN1,KLN2,0)
      NN=0

      DTRUIJ=LIT.GE.3.OR.LJT.GE.3
      DTRUKL=LKT.GE.3.OR.LLT.GE.3

      IC1=1
      IC2=1
      IC3=1
      IC4=INC1+1
      IC5=INC1+1
      IC6=INC1+1
      INC1=INC1*2
      IC7=0

C   ABSORB NORMALISATION FACTORS INTO DENSITY MATRIX

      IF(.NOT.(DTRUIJ.OR.DTRUKL)) GO TO 82
      IF(.NOT.NORM) GO TO 82
      CALL DNORM(ABDENS)
82    CONTINUE

C     K PRIMITIVE

      MAXLG=NGD
      DO 340 KG=1,NGC
        AK=CG(KG)
        BRRK=AK*RRK
        AKXK=AK*XK
        AKYK=AK*YK
        AKZK=AK*ZK
        CSK=CSC(KG)*PI252
        CPK=CPC(KG)*PI252
        CDK=CDC(KG)*PI252

C     ----- L PRIMITIVE

        IF(KANDL)MAXLG=KG
        DO 330 LG=1,MAXLG
          AL=DG(LG)
          B=AK+AL
          BINV=ONE/B
          BBRRK=AL*BRRK*BINV
          IF((BBRRK+RSMALL).GT.TOL1) GO TO 330
          EXKL=DEXP(-BBRRK)
          CSL=CSD(LG)*BINV
          CPL=CPD(LG)*BINV
          CDL=CDD(LG)*BINV
          XB=(AKXK+AL*XL)*BINV
          YB=(AKYK+AL*YL)*BINV
          ZB=(AKZK+AL*ZL)*BINV
          BXBK=B*(XB-XD)
          BYBK=B*(YB-YD)
          BZBK=B*(ZB-ZD)
          BXBI=B*(XB-XC)
          BYBI=B*(YB-YC)
          BZBI=B*(ZB-ZC)

C     ----- DENSITY FACTOR

          DOUBLE=KANDL.AND.KG.NE.LG
          CALL DENFAC(DKL,CSK,CPK,CDK,CSL,CPL,CDL,MINK,MAXK,
     &                MINL,MAXL,KANDL,DOUBLE)
          DKLD=DKL(1)
200       DKLD=DKLD*EXKL
          IF(DABS(DKLD*ABMAX).LT.TOL3)GOTO 330

C     ----- PAIR OF I,J PRIMITIVES

          N=0
          NN=0
          JGMAX=NGB
          DO 320 IG=1,NGA
            AI=AG(IG)
            IF(IANDJ)JGMAX=IG
            DO 320 JG=1,JGMAX
              N=N+1
              IF((BBRRK+R(N)).GE.TOL2)GOTO 320
              AJ=BG(JG)
              DIJD=DD(NN+1)
220           A=AA(N)
              AB=A*B
              AANDB=A+B
              EXPE=DKLD*DIJD*EXIJ(N)/DSQRT(AANDB)
              IF(DABS(EXPE*ABMAX).LT.TOL4)GOTO 320
              RHO=AB/AANDB
              XA=X1(N)
              YA=Y1(N)
              ZA=Z1(N)
              XX=RHO*((XA-XB)**2+(YA-YB)**2+(ZA-ZB)**2)
              AXAK=A*(XA-XD)
              AYAK=A*(YA-YD)
              AZAK=A*(ZA-ZD)
              AXAI=A*(XA-XC)
              AYAI=A*(YA-YC)
              AZAI=A*(ZA-ZC)
              C1X=BXBK+AXAK
              C2X=A*BXBK
              C3X=BXBI+AXAI
              C4X=B*AXAI
              C1Y=BYBK+AYAK
              C2Y=A*BYBK
              C3Y=BYBI+AYAI
              C4Y=B*AYAI
              C1Z=BZBK+AZAK
              C2Z=A*BZBK
              C3Z=BZBI+AZAI
              C4Z=B*AZAI

C     ----- ROOTS AND WEIGHTS FOR QUADRATURE

              CALL ROOTS
              MM=0
              MAX=NMAX+1

C     COMPUTE TWO-ELECTRON  INTEGRALS FOR EACH ROOT

              DO 240 M=1,NROOTS
                U2=U(M)*RHO
                F00=EXPE*W(M)
                DO 230 I=1,MAX
230               IN(I)=IN1(I)+MM
                DUM2=PT5/(AB+U2*AANDB)
                DUM=DUM2+DUM2
                BP01=(A+U2)*DUM2
                B00=U2*DUM2
                B10=(B+U2)*DUM2
                XCP00=(U2*C1X+C2X)*DUM
                XC00 =(U2*C3X+C4X)*DUM
                YCP00=(U2*C1Y+C2Y)*DUM
                YC00 =(U2*C3Y+C4Y)*DUM
                ZCP00=(U2*C1Z+C2Z)*DUM
                ZC00 =(U2*C3Z+C4Z)*DUM
                CALL XYZD(QQX(IC1),QQY(IC2),QQZ(IC3))
240             MM=MM+INC1

C     ----- FORM (I,J//K,L) INTEGRALS

              IOFF=0
              IF(SKIP(1))GOTO 260
              MM=0
              DO 250 M=1,NROOTS
                CALL SUBSD1(QQX(IC1),QQY(IC2),QQZ(IC3),QQX(IC4),
     &                      QQY(IC5),QQZ(IC6),MM,AI)
250             MM=MM+INC1
              CALL DFORM(QQX(IC1),QQY(IC2),QQZ(IC3),QQX(IC4),
     &                   QQY(IC5),QQZ(IC6),ABDENS(1))
              IOFF=IOFF+LENDD
260           IF(SKIP(2))GOTO 280
              MM=0
              DO 270 M=1,NROOTS
                CALL SUBSD2(QQX(IC1),QQY(IC2),QQZ(IC3),QQX(IC4),
     &                      QQY(IC5),QQZ(IC6),MM,AJ)
270             MM=MM+INC1
              CALL DFORM(QQX(IC1),QQY(IC2),QQZ(IC3),QQX(IC4),
     1                   QQY(IC5),QQZ(IC6),ABDENS(1))
              IOFF=IOFF+LENDD
280           IF(SKIP(3))GOTO 300
              MM=0
              DO 290 M=1,NROOTS
                CALL SUBSD3(QQX(IC1),QQY(IC2),QQZ(IC3),QQX(IC4),
     &                      QQY(IC5),QQZ(IC6),MM,AK)
290             MM=MM+INC1
              CALL DFORM(QQX(IC1),QQY(IC2),QQZ(IC3),QQX(IC4),
     1                   QQY(IC5),QQZ(IC6),ABDENS(1))
              IOFF=IOFF+LENDD
300           IF(SKIP(4))GOTO 320
              MM=0
              DO 310 M=1,NROOTS
                CALL SUBSD4(QQX(IC1),QQY(IC2),QQZ(IC3),QQX(IC4),
     &                      QQY(IC5),QQZ(IC6),MM,AL)
310             MM=MM+INC1
              CALL DFORM(QQX(IC1),QQY(IC2),QQZ(IC3),QQX(IC4),
     1                   QQY(IC5),QQZ(IC6),ABDENS(1))

C     END OF LOOPS OVER PRIMITIVES

320           NN=NN+4
330       CONTINUE
340     CONTINUE
      IF(NATOM(1).NE.NATOM(2))GOTO 370
      DGOUT(1)=DGOUT(1)+DGOUT(4)
      DGOUT(2)=DGOUT(2)+DGOUT(5)
      DGOUT(3)=DGOUT(3)+DGOUT(6)
      DGOUT(4)=DGOUT(7)
      DGOUT(5)=DGOUT(8)
      DGOUT(6)=DGOUT(9)
      NPASS=NPASS-1
      NATOM(2)=NATOM(3)
      NATOM(3)=NATOM(4)
      NATOM(4)=0
370   CONTINUE

      RETURN
      END

C===========================================
C    DERIVATIVE SUBSIDIARY INTEGRALS AS LINEAR
C    COMBINATIONS OF ORIGINAL SUBSIDIARY
C    R.D.AMOS OCTOBER 1983 (APPROXIMATELY)
C    OPTIMIZED BY L.A.MONTERO, 1989
C==============================================
      SUBROUTINE SUBSD4(X,Y,Z,XD,YD,ZD,MM,A)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      INTEGER*2 I,J,K,L,N1,N2,N3,N4,N4K,N4KM
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,IJUNK(16)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,JUNK(20)
      DIMENSION X(*),Y(*),Z(*),XD(*),YD(*),ZD(*)
      PARAMETER (TWO=-2.0D0, THREE=-3.0D0)

      A2=A+A
      N1=MM+1

      DO 70 I=1,LIT
        N2=N1
        DO 60 J=1,LJT
          N3=N2
          DO 50 K=1,LKT
            N4=N3
            DO 40 L=1,LLT
              N4K=N4+KLN2
              N4KM=N4-KLN2
              GOTO (10,20,30,35),L
10              XD(N4)=A2*X(N4K)
                YD(N4)=A2*Y(N4K)
                ZD(N4)=A2*Z(N4K)
              GOTO 40
20              XD(N4)=A2*X(N4K)-X(N4KM)
                YD(N4)=A2*Y(N4K)-Y(N4KM)
                ZD(N4)=A2*Z(N4K)-Z(N4KM)
              GOTO 40
30              XD(N4)=A2*X(N4K)+TWO*X(N4KM)
                YD(N4)=A2*Y(N4K)+TWO*Y(N4KM)
                ZD(N4)=A2*Z(N4K)+TWO*Z(N4KM)
              GOTO 40
35              XD(N4)=A2*X(N4K)+THREE*X(N4KM)
                YD(N4)=A2*Y(N4K)+THREE*Y(N4KM)
                ZD(N4)=A2*Z(N4K)+THREE*Z(N4KM)
40            N4=N4+KLN2
50          N3=N3+KLN1
60        N2=N2+IJN2
70      N1=N1+IJN1

      RETURN
      END

      SUBROUTINE SUBSD3(X,Y,Z,XD,YD,ZD,MM,A)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      INTEGER*2 I,J,K,L,N1,N2,N3,N4,N4K,N4KM
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,IJUNK(16)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,JUNK(20)
      DIMENSION X(*),Y(*),Z(*),XD(*),YD(*),ZD(*)
      PARAMETER (TWO=-2.0D0, THREE=-3.0D0)

      A2=A+A
      N1=MM+1

      DO 140 I=1,LIT
        N2=N1
        DO 130 J=1,LJT
          N3=N2
          DO 120 L=1,LLT
            N4=N3
            DO 110 K=1,LKT
              N4K=N4+KLN1
              N4KM=N4-KLN1
              GOTO (80,90,100,105),K
80              XD(N4)=A2*X(N4K)
                YD(N4)=A2*Y(N4K)
                ZD(N4)=A2*Z(N4K)
              GOTO 110
90              XD(N4)=A2*X(N4K)-X(N4KM)
                YD(N4)=A2*Y(N4K)-Y(N4KM)
                ZD(N4)=A2*Z(N4K)-Z(N4KM)
              GOTO 110
100             XD(N4)=A2*X(N4K)+TWO*X(N4KM)
                YD(N4)=A2*Y(N4K)+TWO*Y(N4KM)
                ZD(N4)=A2*Z(N4K)+TWO*Z(N4KM)
              GOTO 110
105             XD(N4)=A2*X(N4K)+THREE*X(N4KM)
                YD(N4)=A2*Y(N4K)+THREE*Y(N4KM)
                ZD(N4)=A2*Z(N4K)+THREE*Z(N4KM)
110           N4=N4K
120         N3=N3+KLN2
130       N2=N2+IJN2
140     N1=N1+IJN1
      RETURN
      END

      SUBROUTINE SUBSD2(X,Y,Z,XD,YD,ZD,MM,A)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      INTEGER*2 I,J,K,L,N1,N2,N3,N4,N4I,N4IM
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,IJUNK(16)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,JUNK(20)
      DIMENSION X(*),Y(*),Z(*),XD(*),YD(*),ZD(*)
      PARAMETER (TWO=-2.0D0, THREE=-3.0D0)

      A2=A+A
      N1=MM+1

      DO 210 I=1,LIT
        N2=N1
        DO 200 K=1,LKT
          N3=N2
          DO 190 L=1,LLT
            N4=N3
            DO 180 J=1,LJT
              N4I=N4+IJN2
              N4IM=N4-IJN2
              GOTO (150,160,170,175),J
150             XD(N4)=A2*X(N4I)
                YD(N4)=A2*Y(N4I)
                ZD(N4)=A2*Z(N4I)
              GOTO 180
160             XD(N4)=A2*X(N4I)-X(N4IM)
                YD(N4)=A2*Y(N4I)-Y(N4IM)
                ZD(N4)=A2*Z(N4I)-Z(N4IM)
              GOTO 180
170             XD(N4)=A2*X(N4I)+TWO*X(N4IM)
                YD(N4)=A2*Y(N4I)+TWO*Y(N4IM)
                ZD(N4)=A2*Z(N4I)+TWO*Z(N4IM)
              GOTO 180
175             XD(N4)=A2*X(N4I)+THREE*X(N4IM)
                YD(N4)=A2*Y(N4I)+THREE*Y(N4IM)
                ZD(N4)=A2*Z(N4I)+THREE*Z(N4IM)
180           N4=N4I
190         N3=N3+KLN2
200       N2=N2+KLN1
210     N1=N1+IJN1
      RETURN
      END

      SUBROUTINE SUBSD1(X,Y,Z,XD,YD,ZD,MM,A)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      INTEGER*2 I,J,K,L,N1,N2,N3,N4,N4I,N4IM
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,IJUNK(16)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,JUNK(20)
      DIMENSION X(*),Y(*),Z(*),XD(*),YD(*),ZD(*)
      PARAMETER (TWO=-2.0D0, THREE=-3.0D0)

      A2=A+A
      N1=MM+1

      DO 280 J=1,LJT
        N2=N1
        DO 270 K=1,LKT
          N3=N2
          DO 260 L=1,LLT
            N4=N3
            DO 250 I=1,LIT
              N4I=N4+IJN1
              N4IM=N4-IJN1
              GOTO (220,230,240,245),I
220             XD(N4)=A2*X(N4I)
                YD(N4)=A2*Y(N4I)
                ZD(N4)=A2*Z(N4I)
              GOTO 250
230             XD(N4)=A2*X(N4I)-X(N4IM)
                YD(N4)=A2*Y(N4I)-Y(N4IM)
                ZD(N4)=A2*Z(N4I)-Z(N4IM)
              GOTO 250
240             XD(N4)=A2*X(N4I)+TWO*X(N4IM)
                YD(N4)=A2*Y(N4I)+TWO*Y(N4IM)
                ZD(N4)=A2*Z(N4I)+TWO*Z(N4IM)
              GOTO 250
245             XD(N4)=A2*X(N4I)+THREE*X(N4IM)
                YD(N4)=A2*Y(N4I)+THREE*Y(N4IM)
                ZD(N4)=A2*Z(N4I)+THREE*Z(N4IM)
250           N4=N4I
260         N3=N3+KLN2
270       N2=N2+KLN1
280     N1=N1+IJN2
      RETURN
      END

C======================================
C    ASSEMBLY OF DERIVATIVE FROM SUBSIDIARIES
C    R.D.AMOS OCTOBER 1983 (APPROXIMATELY)
C =============================================
      SUBROUTINE DFORM(X,Y,Z,XD,YD,ZD,ABDENS)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      LOGICAL NOFORM
      INTEGER*2 I,N,NN,MX,MY,MZ,NR,NOFF1,NOFF2,NOFF3
      COMMON/RTWT/XX,U(NATM1),W(NATM1),NROOTS,JUNK,ZJUNK(NSH1*2)
      COMMON/BIG/AA(100),R(100),X1(100),Y1(100),Z1(100),
     &   DD(400),IJDEN(100),IK(100),
     1  IJX(NATM13),IJY(NATM13),IJZ(NATM13),
     2  KLX(NATM13),KLY(NATM13),KLZ(NATM13),
     3  DIJ(100),DKL(100),IJGT(NATM13),KLGT(NATM13),
     4  NOFORM(1296),DUM4(NBF1M),
     5  DUM5(88)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,INC1,INC2,INC3,INC4,INC5,
     1   LENDD,IOFF,IJUNK(13)
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1 MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2 NIJ,IJ,KL,IJKL
      COMMON/TGRAD/DGOUT(9)
      DIMENSION X(*),Y(*),Z(*),XD(*),YD(*),ZD(*),ABDENS(*)
      PARAMETER (ZERO=0.D0)

      NOFF=IOFF/LENDD
      NOFF1=3*NOFF+1
      NOFF2=NOFF1+1
      NOFF3=NOFF2+1
      I1=INC1
      N=0
      NN=0
      DO 70 I=1,IJ
        DO 70 K=1,KL
          NN=NN+1
          IF (NOFORM(NN)) GOTO 70
          N=N+1
          MX=IJX(I)+KLX(K)
          MY=IJY(I)+KLY(K)
          MZ=IJZ(I)+KLZ(K)
          S1=ZERO
          S2=ZERO
          S3=ZERO
          DO 25 NR=1,NROOTS
            S1=S1+XD(MX)*Y(MY)*Z(MZ)
            S2=S2+X(MX)*YD(MY)*Z(MZ)
            S3=S3+X(MX)*Y(MY)*ZD(MZ)
            MX=I1+MX
            MY=I1+MY
25          MZ=I1+MZ
          ABDENN=ABDENS(N)
          DGOUT(NOFF1)=DGOUT(NOFF1)+S1*ABDENN
          DGOUT(NOFF2)=DGOUT(NOFF2)+S2*ABDENN
          DGOUT(NOFF3)=DGOUT(NOFF3)+S3*ABDENN
70        CONTINUE

      RETURN
      END

C=============================================
C    SUBSIDIARY INTEGRALS
C    FROM ORIGINAL HONDO PROGRAM
C    IMPROVED TO SPEED UP BY L.A.MONTERO, 1989
C=============================================
      SUBROUTINE XYZD(X,Y,Z)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*2 (I-N)
$INCLUDE:'DIMMM'
      LOGICAL N0,N1,M0,M1
      INTEGER*4 I,K,NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,IJKL
      COMMON/SETINT/A(16),I(NATM1),K(NATM1),NIMAX,
     3NJMAX,NKMAX,NLMAX,NMAX,MMAX,IJKL(4)
      DIMENSION X(*),Y(*),Z(*)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)

      BP01=A(1)
      B00=A(2)
      B10=A(3)
      XCP00=A(4)
      XC00=A(5)
      YCP00=A(6)
      YC00=A(7)
      ZCP00=A(8)
      ZC00=A(9)
      F00=A(10)
      DXIJ=A(11)
      DYIJ=A(12)
      DZIJ=A(13)
      DXKL=A(14)
      DYKL=A(15)
      DZKL=A(16)
      IJ1=IJKL(1)
      IJ2=IJKL(2)
      KL1=IJKL(3)
      KL2=IJKL(4)

      N0=NMAX.EQ.0
      N1=NMAX.LE.1
      M0=MMAX.EQ.0
      M1=MMAX.LE.1
      IF(N0)GOTO 320
      IF(M0)GOTO 250

C     ----- I(0,0) -----

      I1=I(1)
      X(I1)=ONE
      Y(I1)=ONE
      Z(I1)=F00

C     ----- I(1,0) -----

      I2=I(2)
      X(I2)=XC00
      Y(I2)=YC00
      Z(I2)=ZC00*F00

C     ----- I(0,1) -----

      K2=K(2)
      I3=I1+K2
      X(I3)=XCP00
      Y(I3)=YCP00
      Z(I3)=ZCP00*F00

C     ----- I(1,1) -----

      I3=I2+K2
      CP10=B00
      X(I3)=XCP00*X(I2)+CP10
      Y(I3)=YCP00*Y(I2)+CP10
      Z(I3)=ZCP00*Z(I2)+CP10*F00
      IF(N1) GO TO 30
      C10=ZERO
      I3=I1
      I4=I2
      DO 20 N=2,NMAX
        C10=C10+B10

C     ----- I(N,0) -----

        I5=I(N+1)
        X(I5)=C10*X(I3)+XC00*X(I4)
        Y(I5)=C10*Y(I3)+YC00*Y(I4)
        Z(I5)=C10*Z(I3)+ZC00*Z(I4)
        CP10=CP10+B00

C     ----- I(N,1) -----

        I3=I5+K2
        X(I3)=XCP00*X(I5)+CP10*X(I4)
        Y(I3)=YCP00*Y(I5)+CP10*Y(I4)
        Z(I3)=ZCP00*Z(I5)+CP10*Z(I4)
10      I3=I4
20      I4=I5
30    IF(M1) GO TO 60
      CP01=ZERO
      C01=B00
      I3=I1
      I4=I1+K2
      DO 50 M=2,MMAX
        CP01=CP01+BP01

C     ----- I(0,M) -----

        I5=I1+K(M+1)
        X(I5)=CP01*X(I3)+XCP00*X(I4)
        Y(I5)=CP01*Y(I3)+YCP00*Y(I4)
        Z(I5)=CP01*Z(I3)+ZCP00*Z(I4)  
        C01=C01+B00

C     ----- I(1,M) -----

        I3=I2+K(M+1)
        X(I3)=XC00*X(I5)+C01*X(I4)
        Y(I3)=YC00*Y(I5)+C01*Y(I4)
        Z(I3)=ZC00*Z(I5)+C01*Z(I4)
40      I3=I4
50      I4=I5
60    IF(N1.OR.M1) GO TO 90

C     ----- I(N,M) -----

      C01=B00
      K3=K2
      DO 80 M=2,MMAX
        K4=K(M+1)
        C01=C01+B00
        I3=I1
        I4=I2
        C10=B10
        DO 70 N=2,NMAX
          I5=I(N+1)
          I5K4=I5+K4
          I3K4=I3+K4
          I4K4=I4+K4
          I4K3=I4+K3
          X(I5K4)=C10*X(I3K4)+XC00*X(I4K4)+C01*X(I4K3)
          Y(I5K4)=C10*Y(I3K4)+YC00*Y(I4K4)+C01*Y(I4K3)
          Z(I5K4)=C10*Z(I3K4)+ZC00*Z(I4K4)+C01*Z(I4K3)
          C10=C10+B10
          I3=I4
70        I4=I5
80      K3=K4
90    IF(NJMAX.EQ.0) GO TO 160

C     ----- I(NI,NJ,M) -----

      M=0
      I5=I(NMAX+1)
100   MIN=NIMAX
      KM=K(M+1)
110   N=NMAX
      I3=I5+KM
120   I4=I(N)+KM
      X(I3)=X(I3)+DXIJ*X(I4)
      Y(I3)=Y(I3)+DYIJ*Y(I4)
      Z(I3)=Z(I3)+DZIJ*Z(I4)
      I3=I4
      N=N-1
      IF(N.GT.MIN) GO TO 120
      MIN=MIN+1
      IF(MIN.LT.NMAX) GO TO 110
      IF(NIMAX.EQ.0) GO TO 150
      I3=IJ2+KM+I1
      DO 140 NJ=1,NJMAX
        I4=I3
        DO 130 NI=1,NIMAX
          I412=I4+IJ1-IJ2
          I42=I4-IJ2
          X(I4)=X(I412)+DXIJ*X(I42)
          Y(I4)=Y(I412)+DYIJ*Y(I42)
          Z(I4)=Z(I412)+DZIJ*Z(I42)
130       I4=I4+IJ1
140     I3=I3+IJ2
150   M=M+1
      IF(M.LE.MMAX) GO TO 100
160   IF(NLMAX.EQ.0) GO TO 240

C     ----- I(NI,NJ,NK,NL) -----

      I5=K(MMAX+1)
      IA=I1
      NI=0
170   NJ=0
      IB=IA
180   MIN=NKMAX
190   M=MMAX
      I3=IB+I5
200   I4=IB+K(M)
      X(I3)=X(I3)+DXKL*X(I4)
      Y(I3)=Y(I3)+DYKL*Y(I4)
      Z(I3)=Z(I3)+DZKL*Z(I4)
      I3=I4
      M=M-1
      IF(M.GT.MIN) GO TO 200
      MIN=MIN+1
      IF(MIN.LT.MMAX) GO TO 190
      IF(NKMAX.EQ.0) GO TO 230
      I3=IB+KL2
      DO 220 NL=1,NLMAX
        I4=I3
        DO 210 NK=1,NKMAX
          I412=I4+KL1-KL2
          I42=I4-KL2
          X(I4)=X(I412)+DXKL*X(I42)
          Y(I4)=Y(I412)+DYKL*Y(I42)
          Z(I4)=Z(I412)+DZKL*Z(I42)
210       I4=I4+KL1
220     I3=I3+KL2
230   NJ=NJ+1
      IB=IB+IJ2
      IF(NJ.LE.NJMAX) GO TO 180
      NI=NI+1
      IA=IA+IJ1
      IF(NI.LE.NIMAX) GO TO 170
 240  RETURN

250   I1 = I(1)
      X(I1) = ONE
      Y(I1) = ONE
      Z(I1) = F00
      IF (N0) GO TO 310

C     ----- I(1,0) -----

      I2 = I(2)
      X(I2) = XC00
      Y(I2) = YC00
      Z(I2) = ZC00*F00
      IF (N1) GO TO 310

C     ----- I(NI,0) -----

      C10 = ZERO
      I3 = I1
      I4 = I2
      DO 260 NI = 2,NMAX
        C10 = C10+B10
        I5 = I(NI+1) 
        X(I5) = C10*X(I3)+XC00*X(I4)
        Y(I5) = C10*Y(I3)+YC00*Y(I4)
        Z(I5) = C10*Z(I3)+ZC00*Z(I4)
        I3 = I4
260     I4 = I5
      IF (NJMAX .EQ. 0) GO TO 310

C     ----- I(NI,NJ,0,0) -----

      I5 = I(NMAX+1)
      MIN = NIMAX
270   NI = NMAX
      I3 = I5
280   I4 = I(NI)
      X(I3) = X(I3)+DXIJ*X(I4)
      Y(I3) = Y(I3)+DYIJ*Y(I4)
      Z(I3) = Z(I3)+DZIJ*Z(I4)
      I3 = I4
      NI = NI-1
      IF (NI .GT. MIN) GO TO 280
      MIN = MIN+1
      IF (MIN .LT. NMAX) GO TO 270
      IF (NIMAX .EQ. 0) GO TO 310
      I3 = IJ2+I1
      DO 300 NJ = 1,NJMAX
        I4 = I3
        DO 290 NI = 1,NIMAX
          I412=I4+IJ1-IJ2
          I42=I4-IJ2
          X(I4) = X(I412)+DXIJ*X(I42)
          Y(I4) = Y(I412)+DYIJ*Y(I42)
          Z(I4) = Z(I412)+DZIJ*Z(I42)
290       I4 = I4+IJ1
300     I3 = I3+IJ2
310   RETURN

320   I1 = I(1)
      X(I1) = ONE
      Y(I1) = ONE
      Z(I1) = F00
      IF (M0) GO TO 380

C     ----- I(0,1) -----

      K2 = K(2)
      I3 = I1+K2
      X(I3) = XCP00
      Y(I3) = YCP00
      Z(I3) = ZCP00*F00
      IF (M1) GO TO 380

C     ----- I(0,NK) -----

      CP01 = ZERO
      I3 = I1
      I4 = I1+K2
      DO 330 NK = 2,MMAX
        CP01 = CP01+BP01
        I5 = I1+K(NK+1)
        X(I5) = CP01*X(I3)+XCP00*X(I4)
        Y(I5) = CP01*Y(I3)+YCP00*Y(I4)
        Z(I5) = CP01*Z(I3)+ZCP00*Z(I4)
        I3 = I4
330     I4 = I5
      IF (NLMAX .EQ. 0) GO TO 380

C     ----- I(0,0,NK,NL) -----

      I5 = K(MMAX+1)
      MIN = NKMAX
340   NK = MMAX
      I3 = I1+I5
350   I4 = I1+K(NK)
      X(I3) = X(I3)+DXKL*X(I4)
      Y(I3) = Y(I3)+DYKL*Y(I4)
      Z(I3) = Z(I3)+DZKL*Z(I4)
      I3 = I4
      NK = NK-1
      IF (NK .GT. MIN) GO TO 350
      MIN = MIN+1
      IF (MIN .LT. MMAX) GO TO 340
      IF (NKMAX .EQ. 0) GO TO 380
      I3 = I1+KL2
      DO 370 NL = 1,NLMAX
        I4 = I3
        DO 360 NK = 1,NKMAX
          I412=I4+KL1-KL2
          I42=I4-KL2
          X(I4) = X(I412)+DXKL*X(I42)
          Y(I4) = Y(I412)+DYKL*Y(I42)
          Z(I4) = Z(I412)+DZKL*Z(I42)
360       I4 = I4+KL1
370     I3 = I3+KL2
380   RETURN
      END

      SUBROUTINE DNORM(ABDENS)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'

      LOGICAL OUT,SKIP(4),IANDJ,KANDL,SAME,NORM
      COMMON/MISC/TOL1,TOL2,TOL3,TOL4,CUTOFF,NATOM(4),NPASS,OUT
     1,SKIP,IANDJ,KANDL,SAME,NORM

      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1 MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2 NIJ,IJ,KL,IJKL
      LOGICAL NOFORM
      COMMON/BIG/AAA(900),IJDEN(100),IK(100),
     1  IJX(NATM13),IJY(NATM13),IJZ(NATM13),
     2  KLX(NATM13),KLY(NATM13),KLZ(NATM13),
     3  DPIJ(100),DPKL(100),IJGT(NATM13),KLGT(NATM13),
     4  NOFORM(1296),DUM4(NBF1M),
     5  DUM5(88)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,INC1,INC2,INC3,INC4,INC5,
     &   LENDD,IJUNK(14)
      DIMENSION DIJ(NATM13),DKL(NATM13),ABDENS(*)
      PARAMETER (ONE=1.0D0)

      ROOT3=DSQRT(3.0D0)
      N=0
      MAX=MAXJ
      DUM1=ONE
      DO 20 I=MINI,MAXI
        IF(I.EQ.8) DUM1=ROOT3
        DUM2=DUM1
        IF(IANDJ) MAX=I
        DO 10 J=MINJ,MAX
          IF(J.EQ.8) DUM2=DUM2*ROOT3
          N=N+1
10        DIJ(N)=DUM2
20      CONTINUE
      N=0
      DUM1=ONE
      MAX=MAXL
      DO 40 K=MINK,MAXK
        IF(K.EQ.8) DUM1=ROOT3
        DUM2=DUM1
        IF(KANDL) MAX=K
        DO 30 L=MINL,MAX 
          IF(L.EQ.8) DUM2=DUM2*ROOT3
          N=N+1
30        DKL(N)=DUM2
40      CONTINUE
      NN=0
      NNN=0
      DO 410 I=1,IJ
        D1=DIJ(I)
        DO 410 K=1,KL
          NN=NN+1
          IF(NOFORM(NN))GOTO 410
          NNN=NNN+1
         ABDENS(NNN)=ABDENS(NNN)*D1*DKL(K)
410      CONTINUE
      RETURN
      END
