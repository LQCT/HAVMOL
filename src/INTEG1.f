      SUBROUTINE STANDV
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      INTEGER P1,P2,P3
      LOGICAL IANDJ,OUT,NORM,DOUBLE
      LOGICAL FIXED,LEX,LDAM12,LDAM13,LDAM23,LDIIS
      COMMON/SCFBLK/EN,ETOT,EHF,SH1(2),SH2(2),GAP1(2),GAP2(2),
     1              D12,D13,D23,CANA,CANB,CANC,FIXED,LEX,
     2              LDAM12,LDAM13,LDAM23,LDIIS,
     3              NCYC,ISCHM,LOCK,MAXIT,NCONV,LOKCYC
      COMMON/TIMEX/TSTART,TI,TX,TIM,TJUNK,TIMMAX,TTOT
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
      LOGICAL REVISE
      COMMON/SECTOR/NUM3,IBLK3,REVISE,IODA(24),ICLASS(24),
     1              ILEN(24),MAXB,KBLKLA
      COMMON/RTWT/XX,U(NATM1),W(NATM1),NROOTS,IJUNK,ZJUNK(NSH1*2)
      CHARACTER*8 PARNAM,BFLAB,ATNAM,ANAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NNP, NE, NA, NB,
     1             ZAN(NATM1),C(3,NATM1),
     &             NVAR,NPAR,NVAL,NUQ,
     &             CFR,CFA,KONTYP(NATM1),P1(NATM1),P2(NATM1),P3(NATM1)
     &             ,KR(NATM1),KA(NATM1),KB(NATM1),PARVAL(NATM13)
     &             ,ZZIN(NATM1),X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAM(NATM1)
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2
      COMMON/BIG/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ,
     1          S(36),G(36),FT(36),
     2          DIJ(100),XIN(100),YIN(100),ZIN(100),
     &          IJX(36),IJY(36),IJZ(36),
     3          Q(2*NBF1M)
      INTEGER*2 I,J,II,JJ,IG,JG,IC,IN,K,MM,NN,NX,NY,NZ
      PARAMETER (ZERO=0.0D0, PT5=0.5D0, ONE=1.0D0, TWO=2.0D0,
     &           THREE=3.0D0, FIVE=5.0D0, SEVEN=7.0D0)
      PARAMETER (PI212=1.1283791670955D0)

      IF(NPRINT.EQ.3) then
        WRITE(IW,10)
10    FORMAT(//,10X,22(1H*)/10X,'ONE-ELECTRON INTEGRALS'/10X,22(1H*)/)
      else
        write(IW,'(/a/)')
     &      ' ..... Calculation of one-electron integrals .....'
        write(10,'(/a/)')
     &      ' ..... Calculation of one-electron integrals .....'
      endif
      IBLOK=0
      II1=1+NNP
      II2=II1+NNP
      IF(II2.GT.MAXDIM)CALL HONDER(24,'STANDV')
      NWORD=NNP
      TOLA=2.30258D0*ITOL
      OUT=NPRINT.EQ.3
      NORM=NORMF.NE.1.OR.NORMP.NE.1
40    DO 50 I=1,NNP
        Q(I)=ZERO
50      Q(II1+I-1)=ZERO
      DO 320 II=1,NSHELL
        I=KATOM(II)
        XI=C(1,I)
        YI=C(2,I)
        ZI=C(3,I)
        I1=KSTART(II)
        I2=I1+KNG(II)-1
        LIT=KTYPE(II)
        MINI=KMIN(II)
        MAXI=KMAX(II)
        LOCI=KLOC(II)-MINI
        DO 320 JJ=1,II
          CUTOFF=10.0D0**(-ICUT)
          TOL=TOLA
          J=KATOM(JJ)
          XJ=C(1,J)
          YJ=C(2,J)
          ZJ=C(3,J)
          J1=KSTART(JJ)
          J2=J1+KNG(JJ)-1
          LJT=KTYPE(JJ)
          MINJ=KMIN(JJ)
          MAXJ=KMAX(JJ)
          LOCJ=KLOC(JJ)-MINJ
          NROOTS=(LIT+LJT)/2
          RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
          IANDJ=II.EQ.JJ
          CALL INDXA(IJX,IJY,IJZ,IJ,MINI,MAXI,MINJ,MAXJ,IANDJ,4,1,1)
          IJ=0
          MAX=MAXJ
          DO 90 I=MINI,MAXI
            IF(IANDJ) MAX=I
            DO 90 J=MINJ,MAX
              IJ=IJ+1
              IF(J.GT.1) GO TO 60
              FT(IJ)=THREE
              GO TO 90
60            CONTINUE
              IF(J.GT.4)GO TO 70
              FT(IJ)=FIVE
              GO TO 90
70            FT(IJ)=SEVEN
90            CONTINUE
          DO 100 I=1,IJ
            S(I)=ZERO
100         G(I)=ZERO
          JGMAX=J2
          DO 270 IG=I1,I2
            AI=EX(IG)
            ARRI=AI*RR
            AXI=AI*XI
            AYI=AI*YI
            AZI=AI*ZI
            CSI=CS(IG)
            CPI=CP(IG)
            CDI=CD(IG)
            IF(IANDJ) JGMAX=IG
            DO 270 JG=J1,JGMAX
              AJ=EX(JG)
              AA=AI+AJ
              AAINV=ONE/AA
              DUM=AJ*ARRI*AAINV
              IF(DUM.GT.TOL) GO TO 270
              FAC=DEXP(-DUM)
              CSJ=CS(JG)*FAC
              CPJ=CP(JG)*FAC
              CDJ=CD(JG)*FAC
              AX=(AXI+AJ*XJ)*AAINV
              AY=(AYI+AJ*YJ)*AAINV
              AZ=(AZI+AJ*ZJ)*AAINV
              DOUBLE=IANDJ.AND.IG.NE.JG
              CALL DENFAN(DIJ,CSI,CPI,CDI,CSJ,CPJ,CDJ,MINI,MAXI,
     1                    MINJ,MAXJ,IANDJ,DOUBLE,NORM)
              T=DSQRT(AA)
              TINV=ONE/T
              T1=-TWO*AJ*AJ*TINV
              T2=-PT5*TINV
              X0=AX
              Y0=AY
              Z0=AZ
              IN=-4
              DO 140 I=1,LIT
                IN=IN+4
                NI=I
                DO 140 J=1,LJT
                  JN=IN+J
                  NJ=J
                  CALL STVINT
                  XIN(JN)=XINT*TINV
                  YIN(JN)=YINT*TINV
                  ZIN(JN)=ZINT*TINV
                  NJ=J+2
                  CALL STVINT
                  XIN(JN+16)=XINT*T1
                  YIN(JN+16)=YINT*T1
                  ZIN(JN+16)=ZINT*T1
                  NJ=J-2
                  IF(NJ.GT.0) GO TO 120
                  XINT=ZERO
                  YINT=ZERO
                  ZINT=ZERO
                  GO TO 130
120               CONTINUE
                  CALL STVINT
130               N=(J-1)*(J-2)
                  DUM=FLOAT(N)*T2
                  XIN(JN+32)=XINT*DUM
                  YIN(JN+32)=YINT*DUM
                  ZIN(JN+32)=ZINT*DUM
140               CONTINUE
              DO 150 I=1,IJ
                NX=IJX(I)
                NY=IJY(I)
                NZ=IJZ(I)
                YZ=YIN(NY)*ZIN(NZ)
                DUM  = YZ*XIN(NX)
                DUM1 = (XIN(NX+16)+XIN(NX+32))*YZ
     1                +(YIN(NY+16)+YIN(NY+32))*XIN(NX)*ZIN(NZ)
     2                +(ZIN(NZ+16)+ZIN(NZ+32))*XIN(NX)*YIN(NY)
                S(I)=S(I)+DIJ(I)*DUM
150             G(I)=G(I)+DIJ(I)*(DUM*AJ*FT(I)+DUM1)
160           CONTINUE
              DUM=PI212/AA
              DO 170 I=1,IJ
170             DIJ(I)=DIJ(I)*DUM
              AAX=AA*AX
              AAY=AA*AY
              AAZ=AA*AZ
190           CONTINUE
              DO 260 IC=1,NAT
                ZNUC=-ZAN(IC)
                CX=C(1,IC)
                CY=C(2,IC)
                CZ=C(3,IC)
                XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
                CALL ROOTS
                MM=0
                DO 230 K=1,NROOTS
                  UU=AA*U(K)
                  WW=W(K)*ZNUC
                  TT=AA+UU
                  T=DSQRT(TT)
                  TINV=ONE/TT
                  X0=(AAX+UU*CX)*TINV
                  Y0=(AAY+UU*CY)*TINV
                  Z0=(AAZ+UU*CZ)*TINV
                  IN=-4+MM
                  DO 220 I=1,LIT
                    IN=IN+4
                    NI=I
                    DO 220 J=1,LJT
                      JN=IN+J
                      NJ=J
                      CALL STVINT
                      XIN(JN)=XINT
                      YIN(JN)=YINT
220                   ZIN(JN)=ZINT*WW
230               MM=MM+16
                DO 250 I=1,IJ
                  NX=IJX(I)
                  NY=IJY(I)
                  NZ=IJZ(I)
                  DUM=ZERO
                  MM=0
                  DO 240 K=1,NROOTS
                    DUM=DUM+XIN(NX+MM)*YIN(NY+MM)*ZIN(NZ+MM)
                    MM=MM+16
 240                CONTINUE
                  FACTOR=DUM*DIJ(I)
                  G(I)=G(I)+FACTOR
250               CONTINUE
260             CONTINUE
270           CONTINUE
          MAX=MAXJ
          NN=0
          DO 310 I=MINI,MAXI
            LI=LOCI+I
            IN=(LI*(LI-1))/2
            IF(IANDJ) MAX=I
            DO 310 J=MINJ,MAX
              LJ=LOCJ+J
              JN=LJ+IN
              NN=NN+1
              Q(JN)=S(NN)
              Q(II1+JN-1)=G(NN)
              IF(DABS(S(NN)).LT.CUTOFF)Q(JN)=ZERO
              IF(DABS(G(NN)).LT.CUTOFF)Q(II1+JN-1)=ZERO
310           CONTINUE
320       CONTINUE
      IF(OUT) CALL PRTHS(Q(1),Q(II1),NUM,IW)
      NBLOK=(NWORD-1)/511+1
      CALL SECPUT(ISEX(5),5,NBLOK,IBLOK)
      CALL WRT3(Q(1),NWORD,IBLOK,IFILD)
      CALL SECPUT(ISEX(6),6,NBLOK,IBLOK)
      CALL WRT3(Q(II1),NWORD,IBLOK,IFILD)
      LDA(ISEX(5))=NWORD
      LDA(ISEX(6))=NWORD
      CALL TIMIT(3)
      CALL REVIND
      WRITE(IW,420)
      write(10,420)
420   FORMAT(/' ..... End of one-electron integral calculation .....'/)
      IF(NPRINT.EQ.-5) CALL TIMIT(3)
      IF(NPRINT.NE.-5) CALL TIMIT(1)
      IF(TIM.LT.TIMLIM) RETURN
      WRITE(IW,430)
      write(*,430)
430   FORMAT(//10X,'*** INSUFFICIENT TIME TO CONTINUE ***'//
     &10X,'This job can be restarted ...')
      CALL CLENUP
      write(*,'(/a/)') ' MICROMOL INTERRUPTED AT SUBR. STANDV'
      write(IW,'(/a/)') ' MICROMOL INTERRUPTED AT SUBR. STANDV'
      STOP
      END

      SUBROUTINE PRTHS(SS,GG,NUM,IW)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      DIMENSION SS(*),GG(*)

10    FORMAT(/10X,13(1H-),/,10X,13HH-CORE MATRIX,/,10X,13(1H-))
20    FORMAT(/10X,14(1H-),/,10X,14HOVERLAP MATRIX,/,10X,14(1H-))
      WRITE(IW,20)
      CALL FOUT(SS,NUM,IW)
      WRITE(IW,10)
      CALL FOUT(GG,NUM,IW)
      RETURN
      END

