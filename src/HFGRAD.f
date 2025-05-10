      SUBROUTINE HFGRAD

C     ----- CALCULATE GRADIENT OF THE HF ENERGY WITH RESPECT
C           TO THE NUCLEAR COORDINATES.

      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GRAD1,GUESS,CONF
      LOGICAL LADAPT,LSKIP
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      COMMON/FUNCT/E,EG(NATM13)
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMMAX,TTOT
      PARAMETER (ZERO=0.0D0)
      DATA GRAD1/'GRADONE'/

      DO 20 IJR=1,NATM13
20      EG(IJR)=ZERO
      CALL TIMIT(3)
      IF((TIMLIM-TIM).LE.3.0D0)GOTO 100

      write(10,8902)
8902  FORMAT(/1H ,'One-electron derivatives')
C     ERLOAD=LOADER(10,'STVDE_.OVL')
      IF(IREST.LE.4) CALL STVDER
      write(10,8903)
8903  FORMAT(/1H ,'One-electron derivatives completed')
      IF(RUNTYP.EQ.GRAD1) RETURN

      write(10,9805)
9805  FORMAT(/1H ,'Two-electron derivatives')
C     ERLOAD=LOADER(10,'JKDER_.OVL')
      IF(IREST.LE.5) CALL JKDER
      write(10,8906)
8906  FORMAT(/1H ,'Two-electron derivatives completed')
      RETURN

100   WRITE(IW,110)
      write(*,110)
110   FORMAT(//10X,'*** INSUFFICIENT TIME TO CONTINUE ***'//
     &10X,'This job can be restarted ...')
      CALL CLENUP
      write(*,'(/a/)') ' MICROMOL INTERRUPTED IN HFGRAD'
      write(IW,'(/a/)') ' MICROMOL INTERRUPTED IN HFGRAD'
      stop
      END

      SUBROUTINE STVDER

C     DRIVING INTEGRAL DERIVATIVE CALCULATIONS

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL LADAPT,LSKIP
      INTEGER P1,P2,P3
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
      CHARACTER*8 ATNAME,ANAM,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMMAX,TTOT

      IBLOK=0
      IBLOKC=0
      IF(IREST.LE.3) IREST=4
      CALL SECGET(ISEX(4),4,IBLOK)
      CALL SECGET(ISEX(24),24,IBLOKC)
      CALL TIMIT(3)
      IF((TIMLIM-TIM).LT.3.0D0)GOTO 10
      CALL HELFEY
      CALL TVDER
      CALL SDER
      CALL WRT3(VIBSIZ,LDA(ISEX(4)),IBLOK,IFILD)
      CALL WRT3C(TITLE,LDA(ISEX(24)),IBLOKC,IFILD)
      RETURN

10    CALL WRT3(VIBSIZ,LDA(ISEX(4)),IBLOK,IFILD)
      CALL WRT3C(TITLE,LDA(ISEX(24)),IBLOKC,IFILD)
      WRITE(IW,20)
      write(*,20)
20    FORMAT(//10X,'*** INSUFFICIENT TIME TO CONTINUE ***'//
     &10X,'This job can be restarted ...')
      CALL CLENUP
      write(*,'(/a/)') ' MICROMOL INTERRUPTED IN STVDER'
      write(IW,'(/a/)') ' MICROMOL INTERRUPTED IN STVDER'
      stop
      END

C   HELFEY  A1
      SUBROUTINE HELFEY

C     CALCULATION OF HELLMANN-FEYNMAN INTEGRALS AND FORCES

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      PARAMETER (NBF1M3=3*NBF1M-1)
      CHARACTER*8 DNAM
      INTEGER P1,P2,P3
      LOGICAL IANDJ,OUT,NORM,DOUBLE,ALLOUT
      COMMON/RTWT/XX,U(NATM1),W(NATM1),NROOTS,NJUNK,ZJUNK(NSH1*2)

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
      CHARACTER*8 ANAM,ATNAME,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NNP, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZZIN(NATM1), X(NATM13)
      COMMON/INFOB/ ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)

      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      COMMON/BIG/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
     1 ,CX,CY,CZ,IJN(100),IJX(NATM13),IJY(NATM13),IJZ(NATM13),
     2 XIN(150),YIN(150),ZIN(150),DFAC(100),
     3 DIJ(100),DUM3(200),HELF(NATM13,NATM1,3),
     4 IDUM4(3,NATM13),DUM4(300),
     5 DUM5(3,NATM1)
      COMMON/BIGG/DDV(NBF1M3),DXTRA,DD(NBF1M)
      COMMON/MAPPER/IA(NBF1+1)
      COMMON/GRAD/DE(3,NATM1)
      COMMON/BLKIN/GOUT(510),NINTH,NJUNK1
      DIMENSION DNAM(3)
      INTEGER*2 I,J,K,II,JJ,IG,JG,IC,MM,ND2R,JN
      PARAMETER (NDIM=5, ZERO=0.0D0, ONE=1.0D0,TWO=2.0D0)
      PARAMETER (PI212=1.1283791670955D0)
      DATA DNAM /'E''X','E''Y','E''Z'/

      CALL ONEPDM(DD(1),DDV(1))

      OUT=NPRINT.EQ.-3
      ALLOUT=NPRINT.EQ.-10
      IF(ALLOUT) WRITE(IW,10)
      IF(ALLOUT)WRITE(IW,11)
      TOL=2.30258D0*ITOL
      DO 20 I=1,3
        DO 20 N=1,NAT
20        DE(I,N)=ZERO
      NINTH=1
      NORM=NORMF.NE.1.OR.NORMP.NE.1

C     ----- ISHELL
      DO 310 II=1,NSHELL
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
C     ----- JSHELL
        DO 310 JJ=1,II
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
          NROOTS=(LIT+LJT+1)/2
          ND2=NDIM*NDIM
          ND2R=ND2*NROOTS
          RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
          IANDJ=II.EQ.JJ
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
          IJ=0
          MAX=MAXJ
          DO 40 I=MINI,MAXI
            IF(IANDJ) MAX=I
            DO 40 J=MINJ,MAX
              IJ=IJ+1
              IJN(IJ)=IA(LOCI+I)+LOCJ+J
              DFAC(IJ)=TWO
40            IF(IANDJ.AND.I.EQ.J) DFAC(IJ)=ONE
          CALL INDXA(IJX,IJY,IJZ,IJ,MINI,MAXI,MINJ,MAXJ,IANDJ,NDIM,1,1)
          DO 60 I=1,IJ
            DO 60 J=1,NAT
              DO 60 K=1,3
60              HELF(I,J,K)=ZERO
C     ----- I PRIMITIVE
          JGMAX=J2
          DO 250 IG=I1,I2
            AI=EX(IG)
            ARRI=AI*RR
            AXI=AI*XI
            AYI=AI*YI
            AZI=AI*ZI
            CSI=CS(IG)
            CPI=CP(IG)
            CDI=CD(IG)
C     ----- J PRIMTIVE
            IF(IANDJ) JGMAX=IG
            DO 250 JG=J1,JGMAX
              AJ=EX(JG)
              AA=AI+AJ
              AAINV=ONE/AA
              DUM=AJ*ARRI*AAINV
              IF(DUM.GT.TOL) GO TO 250
              FAC=DEXP(-DUM)
              CSJ=CS(JG)*FAC
              CPJ=CP(JG)*FAC
              CDJ=CD(JG)*FAC
              AX=(AXI+AJ*XJ)*AAINV
              AY=(AYI+AJ*YJ)*AAINV
              AZ=(AZI+AJ*ZJ)*AAINV
C     ----- DENSITY FACTOR
              DOUBLE=IANDJ.AND.IG.NE.JG
              CALL DENFAN(DIJ,CSI,CPI,CDI,CSJ,CPJ,CDJ,MINI,MAXI,
     &                    MINJ,MAXJ,IANDJ,DOUBLE,NORM)
C     ..... HELLMANN-FEYNMAN TERM .....
              DUM=PI212/AA
              DUM=DUM+DUM
              DO 190 I=1,IJ
190             DIJ(I)=DIJ(I)*DUM
              AAX=AA*AX
              AAY=AA*AY
              AAZ=AA*AZ
              DO 240 IC=1,NAT
                ZNUC=-ZAN(IC)
                CX=C(1,IC)
                CY=C(2,IC)
                CZ=C(3,IC)
                XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
                CALL ROOTS
                MM=0
                DO 210 K=1,NROOTS
                  UU=AA*U(K)
                  WW=W(K)*ZNUC
                  WW=WW*UU
                  TT=AA+UU
                  T=DSQRT(TT)
                  TINV=ONE/TT
                  X0=(AAX+UU*CX)*TINV
                  Y0=(AAY+UU*CY)*TINV
                  Z0=(AAZ+UU*CZ)*TINV
                  IN=-NDIM+MM
                  DO 200 I=1,LIT
                    IN=IN+NDIM
                    NI=I
                    DO 200 J=1,LJT
                      JN=IN+J
                      NJ=J
                      CALL STVINT
                      XIN(JN)=XINT
                      YIN(JN)=YINT
                      ZIN(JN)=ZINT*WW
                      CALL DVINT
                      XIN(JN+ND2R)=XINT
                      YIN(JN+ND2R)=YINT
200                   ZIN(JN+ND2R)=ZINT*WW
210               MM=MM+ND2
                DO 230 I=1,IJ
                  NX=IJX(I)
                  NY=IJY(I)
                  NZ=IJZ(I)
                  DUMX=ZERO
                  DUMY=ZERO
                  DUMZ=ZERO
                  MM=0
                  DO 220 K=1,NROOTS
                    DUMX=DUMX+
     &                   XIN(NX+MM+ND2R)*YIN(NY+MM     )*ZIN(NZ+MM     )
                    DUMY=DUMY+
     &                   XIN(NX+MM     )*YIN(NY+MM+ND2R)*ZIN(NZ+MM     )
                    DUMZ=DUMZ+
     &                   XIN(NX+MM     )*YIN(NY+MM     )*ZIN(NZ+MM+ND2R)
220                 MM=MM+ND2
                  DUM=DIJ(I)
                  DUMX=DUMX*DUM
                  DUMY=DUMY*DUM
                  DUMZ=DUMZ*DUM
                  HELF(I,IC,1)=HELF(I,IC,1)+DUMX
                  HELF(I,IC,2)=HELF(I,IC,2)+DUMY
                  HELF(I,IC,3)=HELF(I,IC,3)+DUMZ
                  DUM=DD(IJN(I))*DFAC(I)
                  DE(1,IC)=DE(1,IC)+DUM*DUMX
                  DE(2,IC)=DE(2,IC)+DUM*DUMY
230               DE(3,IC)=DE(3,IC)+DUM*DUMZ
240             CONTINUE
250           CONTINUE
          IF(.NOT.ALLOUT)GOTO 280
      IIJ=0
      DO 260 I=MINI,MAXI
        MAX=MAXJ
        IF(IANDJ)MAX=I
        IP=LOCI+I
        DO 260 J=MINJ,MAX
          IIJ=IIJ+1
          JP=LOCJ+J
          DO 260 K=1,NAT
260         WRITE(IW,270) IP,JP,K,(HELF(IIJ,K,L),L=1,3),DD(IJN(IIJ))
280       CONTINUE
C     ----- END OF 'PRIMITIVE' LOOPS -----
310       CONTINUE
C     ----- END OF 'SHELL' LOOPS -----

320   IF(.NOT.OUT) GOTO 350
      MAX=0
330   MIN=MAX+1
      MAX=MAX+5
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,360)
      WRITE(IW,370) (I,I=MIN,MAX)
      WRITE(IW,360)
      DO 340 N=1,3
340     WRITE(IW,380) DNAM(N),(DE(N,I),I=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 330
350   CALL NUCDER(NPRINT,IW)
      IF(OUT) CALL WHTPS
      WRITE(IW,390)
10    FORMAT(/,10X,26(1H-),/,10X,'Hellmann-Feynman integrals',
     1 /,10X,26(1H-))
11    FORMAT(//5X,'I',4X,'J',4X,'K',15X,'Ex',18X,'Ey',18X,'Ez',18X,
     & 'Dij')
270   FORMAT(1X,3I5,5X,4F20.8)
360   FORMAT(/)
370   FORMAT(5X,'ATOM',5(6X,I2,6X))
380   FORMAT(7X,A3,5D14.7)
390   FORMAT(/,' ...... End of Hellmann-Feynman force ......',/)
      RETURN
      END

C  SDER  A1
C  AUTODBL VERSION
      SUBROUTINE SDER

C     CALCULATION OF ONE ELECTRON INTEGRAL DERIVATIVES

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL OUT,NORM
      INTEGER P1,P2,P3

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
      CHARACTER*8 ANAM,ATNAME,PARNAM,BFLAB
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)

      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),KTYPE(NSH1),
     2             KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      COMMON/BIG/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
     1 ,CX,CY,CZ,IDUM1(100),IJX(NATM13),IJY(NATM13),IJZ(NATM13),
     2 XIN(75),YIN(75),ZIN(75),XD(75),YD(75),ZD(75),DIJ(100),
     3 SX(100),SY(100),SZ(100),DUM31(NATM13,NATM1,3),
     4 IDUM4(3,NATM13),DUM4(300),
     5 DUM5(3,NATM1)
      COMMON/BIGG/DDX(NBF1M),DDY(NBF1M),DDZ(NBF1M),DD(NBF1M)
      COMMON/MAPPER/IA(NBF1+1)
      COMMON/BLKIN/GOUT(510),NINTH,NJNK
      COMMON/GRAD/GRD(NATM13)
      DIMENSION DDV(1)
      INTEGER*2 I,J,II,JJ,IAT,JAT,IG,JG,JN,NN,N
      EQUIVALENCE(DDV(1),DDX(1))
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      DATA NDIM/5/

C     ----- CALCULATE DERIVATIVES OF THE OVERLAP MATRIX -----

      TOL=2.30258D0*ITOL
      IBLOK=0
      OUT=NPRINT.EQ.-10
      IF(OUT)WRITE(IW,10)
      NORM=NORMF.NE.1.OR.NORMP.NE.1
      LAB1=1
      LAB2=LAB1+NUM*NUM
      CALL LAGRAN (DD(1),DDV(LAB1),DDV(LAB2),NUM)
      LAB1=0
      LAB2=0
*      LAB3=0
      NINTH=1
C     ----- ISHELL
      DO 80 II=1,NSHELL
        IAT=KATOM(II)
        XI=C(1,IAT)
        YI=C(2,IAT)
        ZI=C(3,IAT)
        I1=KSTART(II)
        I2=I1+KNG(II)-1
        LIT=KTYPE(II)
        LIT1=LIT+1
        MINI=KMIN(II)
        MAXI=KMAX(II)
        LOCI=KLOC(II)-MINI
C     ----- JSHELL
        DO 80 JJ=1,II
          JAT=KATOM(JJ)
          IF(IAT.EQ.JAT) GOTO 80
          XJ=C(1,JAT)
          YJ=C(2,JAT)
          ZJ=C(3,JAT)
          J1=KSTART(JJ)
          J2=J1+KNG(JJ)-1
          LJT=KTYPE(JJ)
          MINJ=KMIN(JJ)
          MAXJ=KMAX(JJ)
          LOCJ=KLOC(JJ)-MINJ
          RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
          CALL INDXA (IJX,IJY,IJZ,IJ,MINI,MAXI,MINJ,MAXJ,
     &                .FALSE.,NDIM,1,1)
          DO 20 I=1,IJ
            SX(I)=ZERO
            SY(I)=ZERO
20          SZ(I)=ZERO
C     ----- I PRIMITIVE
          DO 50 IG=I1,I2
            AI=EX(IG)
            ARRI=AI*RR
            AXI=AI*XI
            AYI=AI*YI
            AZI=AI*ZI
            CSI=CS(IG)
            CPI=CP(IG)
            CDI=CD(IG)
C     ----- J PRIMTIVE
            DO 50 JG=J1,J2
              AJ=EX(JG)
              AA=AI+AJ
              AAINV=ONE/AA
              DUM=AJ*ARRI*AAINV
              IF(DUM.GT.TOL) GO TO 50
              FAC=DEXP(-DUM)
              CSJ=CS(JG)*FAC
              CPJ=CP(JG)*FAC
              CDJ=CD(JG)*FAC
              AX=(AXI+AJ*XJ)*AAINV
              AY=(AYI+AJ*YJ)*AAINV
              AZ=(AZI+AJ*ZJ)*AAINV
C     ----- DENSITY FACTOR
              CALL DENFAN(DIJ,CSI,CPI,CDI,CSJ,CPJ,CDJ,MINI,MAXI,
     1                    MINJ,MAXJ,.FALSE.,.FALSE.,NORM)
C     ----- OVERLAP
              T=DSQRT(AA)
              TINV=ONE/T
              X0=AX
              Y0=AY
              Z0=AZ
              IN=-NDIM
              DO 30 I=1,LIT1
                IN=IN+NDIM
                NI=I
                DO 30 J=1,LJT
                  JN=IN+J
                  NJ=J
                  CALL STVINT
                  XIN(JN)=XINT*TINV
                  YIN(JN)=YINT*TINV
30                ZIN(JN)=ZINT*TINV

              CALL ONELD(XIN,YIN,ZIN,XD,YD,ZD,AI,LIT,LJT,1,NDIM)

              DO 40 I=1,IJ
                MX=IJX(I)
                MY=IJY(I)
                MZ=IJZ(I)
                SX(I)=SX(I)+DIJ(I)*XD(MX)*YIN(MY)*ZIN(MZ)
                SY(I)=SY(I)+DIJ(I)*XIN(MX)*YD(MY)*ZIN(MZ)
40              SZ(I)=SZ(I)+DIJ(I)*XIN(MX)*YIN(MY)*ZD(MZ)
C     ----- END OF PRIMITIVE LOOPS -----
50            CONTINUE
C     ----- CALCULATE DERIVATIVES OF OVERLAP MATRIX -----
          N=0
          DO 70 I=MINI,MAXI
            IN=LOCI+I
            DO 60 J=MINJ,MAXJ
              N=N+1
              JN=LOCJ+J
              IF(JN.GT.IN) GO TO 60
              NN=IA(IN)+JN
              DDX(NN)=SX(N)
              DDY(NN)=SY(N)
              DDZ(NN)=SZ(N)
              IF(OUT) WRITE(IW,100) IN,JN,SX(N),SY(N),SZ(N),DD(NN)
60            CONTINUE
70          CONTINUE
80        CONTINUE
90    CONTINUE
C     ----- END OF SHELL LOOPS -----

      CALL GLIMIT
C     ----- FORM THE 1E-GRADIENT -----
      CALL SGRAD(DD,DDX,DDY,DDZ,NPRINT,IW)
101   IST=1
      JST=1
      KST=1
      LST=1
      NREC=1
      INTLOC=1
      IREST=5
      LDA(ISEX(14))=NATM13
      CALL SECPUT(ISEX(14),14,LENSEC(NATM13),IBLOK)
      CALL WRT3(GRD,LDA(ISEX(14)),IBLOK,IFILD)
      CALL TIMIT(1)
      CALL REVIND
10    FORMAT(///5X,'I',4X,'J',15X,'Sx',18X,'Sy',18X,'Sz',18X,'Lij')
100   FORMAT(1X,2I5,5X,4F20.8)
      RETURN
      END

C  SGRAD  M2
      SUBROUTINE SGRAD(DD,DDX,DDY,DDZ,NPRINT,IW)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 DNAM
      INTEGER P1,P2,P3
      LOGICAL OUT
      CHARACTER*8 ATNAME,ANAM,PARNAM,BFLAB
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13),  ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      COMMON/SMALL/EG(3,NATM1),DJNK(132)
      COMMON/MAPPER/IA(NBF1+1)
      COMMON/ALIM/LIMINF(NATM1),LIMSUP(NATM1)
      COMMON/GRAD/DE(3,NATM1)
      DIMENSION DD(*),DDX(*),DDY(*),DDZ(*)
      DIMENSION DNAM(3)
      INTEGER*2 I,J,K,L,KK,N
      PARAMETER (ZERO=0.0D0, TWO=2.0D0)
      DATA DNAM /'E''X','E''Y','E''Z'/

      OUT=NPRINT.EQ.-3.OR.NPRINT.EQ.-10
      EG(1,1)=ZERO
      EG(2,1)=ZERO
      EG(3,1)=ZERO
      DO 20 K=2,NAT
        EG(1,K)=ZERO
        EG(2,K)=ZERO
        EG(3,K)=ZERO
        KINF=LIMINF(K)
        KSUP=LIMSUP(K)
        KM1=K-1
        DO 20 L=1,KM1
          LINF=LIMINF(L)
          LSUP=LIMSUP(L)
          IF(KINF.EQ.0.OR.KSUP.EQ.0) GOTO 20
          IF(LINF.EQ.0.OR.LSUP.EQ.0) GOTO 20
          DO 10 I=KINF,KSUP
            DO 10 J=LINF,LSUP
              IJ=IA(I)+J
              EG(1,K)=EG(1,K)+DDX(IJ)*DD(IJ)
              EG(2,K)=EG(2,K)+DDY(IJ)*DD(IJ)
              EG(3,K)=EG(3,K)+DDZ(IJ)*DD(IJ)
10            CONTINUE
20        CONTINUE
      NAT1=NAT-1
      DO 40 K=1,NAT1
        KINF=LIMINF(K)
        KSUP=LIMSUP(K)
        KP1=K+1
        DO 40 L=KP1,NAT
          LINF=LIMINF(L)
          LSUP=LIMSUP(L)
          IF(KINF.EQ.0.OR.KSUP.EQ.0) GOTO 40
          IF(LINF.EQ.0.OR.LSUP.EQ.0) GOTO 40
          DO 30 I=KINF,KSUP
            DO 30 J=LINF,LSUP
              IJ=IA(J)+I
              EG(1,K)=EG(1,K)-DDX(IJ)*DD(IJ)
              EG(2,K)=EG(2,K)-DDY(IJ)*DD(IJ)
              EG(3,K)=EG(3,K)-DDZ(IJ)*DD(IJ)
30            CONTINUE
40        CONTINUE
      DO 60 K=1,NAT
        DO 60 KK=1,3
60        EG(KK,K)=TWO*EG(KK,K)
C     ----- ADD ALL CONTRIBUTIONS TO 1E-GRADIENT -----
      DO 70 I=1,NAT
        DO 70 J=1,3
70        DE(J,I)=DE(J,I)+EG(J,I)
C     ----- PRINT SECTION -----
      IF(.NOT.OUT) GOTO 170
      WRITE(IW,150)
      MAX=0
80    MIN=MAX+1
      MAX=MAX+5
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,120)
      WRITE(IW,130) (I,I=MIN,MAX)
      WRITE(IW,120)
      DO 90 N=1,3
90      WRITE(IW,140) DNAM(N),(EG(N,I),I=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 80
      WRITE(IW,160)
      MAX=0
100   MIN=MAX+1
      MAX=MAX+5
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,120)
      WRITE(IW,130) (I,I=MIN,MAX)
      WRITE(IW,120)
      DO 110 N=1,3
110     WRITE(IW,140) DNAM(N),(DE(N,I),I=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 100
120   FORMAT(/)
130   FORMAT(5X,'ATOM',4(6X,I2,6X))
140   FORMAT(7X,A3,5D14.7)
150   FORMAT(/,10X,32(1H-),/,10X,'Overlap contribution to gradient',
     1 /,10X,32(1H-))
160   FORMAT(/,10X,22(1H-),/,10X,'One-electron gradient',
     1 /,10X,22(1H-))
170   RETURN
      END

      SUBROUTINE TVDER

C     CALCULATION OF OTHER ONE ELECTRON CONTRIBUTION TO GRADIENTS

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      PARAMETER (NBF1M3=3*NBF1M-1)
      CHARACTER*8 DNAM
      INTEGER P1,P2,P3
      LOGICAL ALLOUT,OUT,NORM

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
      COMMON/RTWT/XX,U(NATM1),W(NATM1),NROOTS,NJNK,XJNK(NSH1*2)
      CHARACTER*8 BFLAB,PARNAM,ATNAME,ANAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NNP, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5     PARVAL(NATM13), ZZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)

      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),KTYPE(NSH1),
     2             KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      CHARACTER*8 GROUP
      COMMON/SYMTRY/TTTT(216),INVT(24),ISO(NSH1+1,24),NT,NT2,XXX(21)
     &,INDEX,NAXIS
      COMMON/SYMTRB/GROUP
      COMMON/BIG/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
     1 ,CX,CY,CZ,IDUM10(100),IDUM12(3,NATM13),
     2 XIN(75),YIN(75),ZIN(75),XAD(75),YAD(75),ZAD(75),DUM2(100),
     2 DELX(25),DELY(25),DELZ(25),DELXA(25),DELYA(25),DELZA(25),
     3 DIJ(100),DUM30(50),DUM31(NATM13,NATM1,3),
     4 IJX(NATM13),IJY(NATM13),IJZ(NATM13),TVX(100),TVY(100),TVZ(100),
     5 TVD(3,NATM1)
      COMMON/BIGG/DDV(NBF1M3),DXTRA,DD(NBF1M)
      COMMON/GRAD/DE(3,NATM1)
      COMMON/BLKIN/GOUT(510),NINTH,NJNK1
      COMMON/MAPPER/IA(NBF1+1)
      INTEGER*2 I,J,K,IC,II,JJ,IAT,JAT,IG,JG,N,NN
      DIMENSION DNAM(3)

      PARAMETER (TWO=2.0D0,HALF=0.5D0,ZERO=0.0D0,ONE=1.0D0)
      PARAMETER (PI212=1.1283791670955D0)
      DATA NDIM/5/
      DATA DNAM /'E''X','E''Y','E''Z'/

      CALL ONEPDM(DD(1),DDV(1))
      OUT=NPRINT.EQ.-3.OR.NPRINT.EQ.-10
      ALLOUT=NPRINT.EQ.-10
      DO 10 J=1,NAT
        DO 10 I=1,3
10        TVD(I,J)=ZERO
      IF(OUT) WRITE(IW,20)
      IF(ALLOUT)WRITE(IW,30)
      TOL=2.30258D0*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
      NINTH=1
C     ----- ISHELL
      DO 190 II=1,NSHELL
        IAT=KATOM(II)
        XI=C(1,IAT)
        YI=C(2,IAT)
        ZI=C(3,IAT)
        I1=KSTART(II)
        I2=I1+KNG(II)-1
        LIT=KTYPE(II)
        LIT1=LIT+1
        MINI=KMIN(II)
        MAXI=KMAX(II)
        LOCI=KLOC(II)-MINI
C     ----- JSHELL
        DO 190 JJ=1,NSHELL
          JAT=KATOM(JJ)
          XJ=C(1,JAT)
          YJ=C(2,JAT)
          ZJ=C(3,JAT)
          J1=KSTART(JJ)
          J2=J1+KNG(JJ)-1
          LJT=KTYPE(JJ)
          MINJ=KMIN(JJ)
          MAXJ=KMAX(JJ)
          LOCJ=KLOC(JJ)-MINJ
          NROOTS=(LIT+LJT+1)/2
          RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2

C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
          CALL INDXA(IJX,IJY,IJZ,IJ,MINI,MAXI,MINJ,MAXJ,.FALSE.,
     1               NDIM,1,1)
          DO 40 I=1,IJ
            TVX(I)=ZERO
            TVY(I)=ZERO
40          TVZ(I)=ZERO

C     ----- I PRIMITIVE

          DO 130 IG=I1,I2
            AI=EX(IG)
            ARRI=AI*RR
            AXI=AI*XI
            AYI=AI*YI
            AZI=AI*ZI
            CSI=CS(IG)
            CPI=CP(IG)
            CDI=CD(IG)

C     ----- J PRIMTIVE

            DO 120 JG=J1,J2
              AJ=EX(JG)
              AA=AI+AJ
              AAINV=ONE/AA
              DUM=AJ*ARRI*AAINV
              IF(DUM.GT.TOL) GO TO 120
              FAC=DEXP(-DUM)
              CSJ=CS(JG)*FAC
              CPJ=CP(JG)*FAC
              CDJ=CD(JG)*FAC
              AX=(AXI+AJ*XJ)*AAINV
              AY=(AYI+AJ*YJ)*AAINV
              AZ=(AZI+AJ*ZJ)*AAINV

C     ----- DENSITY FACTOR

              CALL DENFAN(DIJ,CSI,CPI,CDI,CSJ,CPJ,CDJ,MINI,MAXI,
     1                    MINJ,MAXJ,.FALSE.,.FALSE.,NORM)

              T=DSQRT(AA)
              TINV=ONE/T
              T1=-TWO*AJ*AJ*TINV
              T2=AJ*TINV
              T3=-HALF*TINV
              X0=AX
              Y0=AY
              Z0=AZ
              IN=-NDIM
              DO 50 I=1,LIT1
                IN=IN+NDIM
                NI=I
                DO 50 J=1,LJT
                  JN=IN+J
                  NJ=J
                  CALL STVINT
                  XIN(JN)=XINT*TINV
                  YIN(JN)=YINT*TINV
                  ZIN(JN)=ZINT*TINV
C     ELEMENTS OF DEL-SQUARED
                  DUM=DBLE(FLOAT(J+J-1))*T2
                  DELX(JN)=DUM*XINT
                  DELY(JN)=DUM*YINT
                  DELZ(JN)=DUM*ZINT
                  NJ=J+2
                  CALL STVINT
                  DELX(JN)=DELX(JN)+XINT*T1
                  DELY(JN)=DELY(JN)+YINT*T1
                  DELZ(JN)=DELZ(JN)+ZINT*T1
                  IF(J.LE.2)GOTO 50
                  NJ=J-2
                  CALL STVINT
                  DUM=DBLE(FLOAT((J-1)*(J-2)))*T3
                  DELX(JN)=DELX(JN)+XINT*DUM
                  DELY(JN)=DELY(JN)+YINT*DUM
                  DELZ(JN)=DELZ(JN)+ZINT*DUM
50                CONTINUE
              CALL ONELD(XIN,YIN,ZIN,XAD,YAD,ZAD,AI,LIT,LJT,1,NDIM)
              CALL ONELD(DELX,DELY,DELZ,DELXA,DELYA,DELZA,AI,LIT,LJT,1,
     1                   NDIM )
              DO 60 I=1,IJ
                NX=IJX(I)
                NY=IJY(I)
                NZ=IJZ(I)
                D1=DIJ(I)
                TVX(I)=TVX(I)+D1*(DELXA(NX) * YIN(NY) * ZIN(NZ)
     1                             +XAD(NX) *DELY(NY) * ZIN(NZ)
     2                             +XAD(NX) * YIN(NY) *DELZ(NZ))
                TVY(I)=TVY(I)+D1*(DELX(NX)  * YAD(NY) * ZIN(NZ)
     1                            +XIN(NX) *DELYA(NY) * ZIN(NZ)
     2                            +XIN(NX) *  YAD(NY) *DELZ(NZ))
                TVZ(I)=TVZ(I)+D1*(DELX(NX) * YIN(NY) *  ZAD(NZ)
     1                            +XIN(NX) *DELY(NY) *  ZAD(NZ)
     2                            +XIN(NX) * YIN(NY) *DELZA(NZ))
60              CONTINUE

C     ..... NUCLEAR ATTRACTION

              DUM=PI212/AA
              DO 70 I=1,IJ
70              DIJ(I)=DIJ(I)*DUM
              AAX=AA*AX
              AAY=AA*AY
              AAZ=AA*AZ
              DO 110 IC=1,NAT
                ZNUC=-ZAN(IC)
                CX=C(1,IC)
                CY=C(2,IC)
                CZ=C(3,IC)
                XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
                CALL ROOTS
                DO 100 K=1,NROOTS
                  UU= AA*U(K )
                  WW=W(K)*ZNUC
                  TT=AA+UU
                  T=DSQRT(TT)
                  TINV=ONE/TT
                  X0=(AAX+UU*CX)*TINV
                  Y0=(AAY+UU*CY)*TINV
                  Z0=(AAZ+UU*CZ)*TINV
                  IN=-NDIM
                  DO 80 I=1,LIT1
                    IN=IN+NDIM
                    NI=I
                    DO 80 J=1,LJT
                      JN=IN+J
                      NJ=J
                      CALL STVINT
                      XIN(JN)=XINT
                      YIN(JN)=YINT
80                    ZIN(JN)=ZINT*WW
                  CALL ONELD(XIN,YIN,ZIN,XAD,YAD,ZAD,AI,LIT,LJT,1,NDIM)
                  DO 90 I=1,IJ
                    NX=IJX(I)
                    NY=IJY(I)
                    NZ=IJZ(I)
                    D1=DIJ(I)
                    TVX(I)=TVX(I)+D1*XAD(NX)*YIN(NY)*ZIN(NZ)
                    TVY(I)=TVY(I)+D1*XIN(NX)*YAD(NY)*ZIN(NZ)
                    TVZ(I)=TVZ(I)+D1*XIN(NX)*YIN(NY)*ZAD(NZ)
90                  CONTINUE
100               CONTINUE
110             CONTINUE
C     ----- END OF PRIMITIVE LOOPS -----

120           CONTINUE
130         CONTINUE

          N=0
          DO 180 I=MINI,MAXI
            IN=LOCI+I
            DO 180 J=MINJ,MAXJ
              N=N+1
              JN=LOCJ+J
              IF(IN.LT.JN) THEN
                ID=JN
                JD=IN
              ELSE
                ID=IN
                JD=JN
              ENDIF
160           NN=IA(ID)+JD
              DUM=DD(NN)
              DUM=DUM+DUM
              TVD(1,IAT)=TVD(1,IAT)+DUM*TVX(N)
              TVD(2,IAT)=TVD(2,IAT)+DUM*TVY(N)
              TVD(3,IAT)=TVD(3,IAT)+DUM*TVZ(N)
              IF(ALLOUT)
     &          WRITE(IW,270) IN,JN,TVX(N),TVY(N),TVZ(N),DD(NN)
180           CONTINUE
190       CONTINUE
C     ----- END OF SHELL LOOPS -----

      DO 200 J=1,NAT
        DO 200 I=1,3
200       DE(I,J)=DE(I,J)+TVD(I,J)
      IF(.NOT.OUT) GO TO 230
      MAX=0
210   MIN=MAX+1
      MAX=MAX+5
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,240)
      WRITE(IW,250) (I,I=MIN,MAX)
      WRITE(IW,240)
      DO 220 N=1,3
220     WRITE(IW,260) DNAM(N),(TVD(N,I),I=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 210

230   CONTINUE
280   IF(OUT) CALL WHTPS

20    FORMAT(//,10X,45(1H-),/,10X,'Other one-electron contribution',
     1 ' to the gradient',/,10X,45(1H-))
30    FORMAT(//5X,'I',4X,'J',15X,'Ex',18X,'Ey',18X,'Ez',18X,'Dij')
240   FORMAT(/)
250   FORMAT(5X,'ATOM',5(6X,I2,6X))
260   FORMAT(7X,A3,5D14.7)
270   FORMAT(1X,2I5,5X,4F17.8)
      RETURN
      END

      SUBROUTINE GLIMIT
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      INTEGER P1,P2,P3
      CHARACTER*8 ATNAME,PARNAM,BFLAB,ANAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13),  ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      COMMON/ALIM/LIMINF(NATM1),LIMSUP(NATM1)

      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2
      INTEGER*2 I,LAT,IAT

      DO 11 I=1,NAT
        LIMSUP(I)=0
11      LIMINF(I)=0
      LAT=KATOM(1)
      LIMINF(LAT)=1
      DO 10 I=1,NSHELL
        IAT=KATOM(I)
        IF(LAT.NE.IAT) THEN
          LIMSUP(LAT)=KLOC(I)-1
          LAT=IAT
          LIMINF(LAT)=KLOC(I)
        ENDIF
   10   CONTINUE
      LIMSUP(LAT)=NUM
      RETURN
      END

      SUBROUTINE LAGRAN(DD,V,E,NDIM)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      INTEGER P1,P2,P3
      DIMENSION DD(*),E(*),V(NDIM,*)

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
      CHARACTER*8 OSCF
      CHARACTER*8 PARNAM,ATNAME,ANAM,BFLAB
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13),  ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      PARAMETER (ZERO=0.0D0, TWO=2.0D0)
      DATA OSCF/'OSCF'/

      IF(SCFTYP.EQ.OSCF) GO TO 50
      OCC=TWO
      MTYP=0
      IBLOK=0
      CALL SECGET(ISEX(8),MTYP,IBLOK)
      IBLVEC=IBLOK
      CALL READQ(V,NUM*NCOORB,IBLVEC,IFILD)
      MTYP=0
      CALL SECGET(ISEX(9),MTYP,IBLOK)
      CALL READQ(E,LDA(ISEX(9)),IBLOK,IFILD)
      IJ=0
      DO 20 I=1,NUM
        DO 20 J=1,I
          IJ=IJ+1
          DUM=ZERO
          DO 10 K=1,NA
10          DUM=DUM-E(K)*V(I,K)*V(J,K)
20        DD(IJ)=DUM*OCC
      RETURN

50    IBLOK=0
      MTYP=0
      CALL SECGET(ISEX(19),MTYP,IBLOK)
      CALL READQ(DD,LDA(ISEX(19)),IBLOK,IFILD)
      RETURN
      END

      SUBROUTINE ONELD(X,Y,Z,XD,YD,ZD,A,NI,NJ,ITYPE,NDIM)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      DIMENSION X(NDIM,NDIM),Y(NDIM,NDIM),Z(NDIM,NDIM),
     & XD(NDIM,NDIM),YD(NDIM,NDIM),ZD(NDIM,NDIM)
      PARAMETER (TWO=-2.0D0, THREE=-3.0D0, FOUR=-4.0D0)

      A2=A+A
      GOTO (10,80),ITYPE
10    DO 70 J=1,NJ
        DO 70 I=1,NI
          GOTO (20,30,40,50,60),I
20          XD(J,1)=A2*X(J,2)
            YD(J,1)=A2*Y(J,2)
            ZD(J,1)=A2*Z(J,2)
          GOTO 70
30          XD(J,2)=A2*X(J,3)-X(J,1)
            YD(J,2)=A2*Y(J,3)-Y(J,1)
            ZD(J,2)=A2*Z(J,3)-Z(J,1)
          GOTO 70
40          XD(J,3)=A2*X(J,4)+TWO*X(J,2)
            YD(J,3)=A2*Y(J,4)+TWO*Y(J,2)
            ZD(J,3)=A2*Z(J,4)+TWO*Z(J,2)
          GOTO 70
50          XD(J,4)=A2*X(J,5)+THREE*X(J,3)
            YD(J,4)=A2*Y(J,5)+THREE*Y(J,3)
            ZD(J,4)=A2*Z(J,5)+THREE*Z(J,3)
          GOTO 70
60          XD(J,5)=A2*X(J,6)+FOUR*X(J,4)
            YD(J,5)=A2*Y(J,6)+FOUR*Y(J,4)
            ZD(J,5)=A2*Z(J,6)+FOUR*Z(J,4)
70        CONTINUE
      RETURN

80    CONTINUE
      DO 140 I=1,NI
        DO 140 J=1,NJ
          GOTO (90,100,110,120,130),J
90          XD(1,I)=A2*X(2,I)
            YD(1,I)=A2*Y(2,I)
            ZD(1,I)=A2*Z(2,I)
          GOTO 140
100         XD(2,I)=A2*X(3,I)-X(1,I)
            YD(2,I)=A2*Y(3,I)-Y(1,I)
            ZD(2,I)=A2*Z(3,I)-Z(1,I)
          GOTO 140
110         XD(3,I)=A2*X(4,I)+TWO*X(2,I)
            YD(3,I)=A2*Y(4,I)+TWO*Y(2,I)
            ZD(3,I)=A2*Z(4,I)+TWO*Z(2,I)
          GOTO 140
120         XD(4,I)=A2*X(5,I)+THREE*X(3,I)
            YD(4,I)=A2*Y(5,I)+THREE*Y(3,I)
            ZD(4,I)=A2*Z(5,I)+THREE*Z(3,I)
          GOTO 140
130         XD(5,I)=A2*X(6,I)+FOUR*X(4,I)
            YD(5,I)=A2*Y(6,I)+FOUR*Y(4,I)
            ZD(5,I)=A2*Z(6,I)+FOUR*Z(4,I)
140       CONTINUE
      RETURN
      END
