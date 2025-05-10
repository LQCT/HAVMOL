      SUBROUTINE NUCDER(NPRINT,IW)

C     CALCULATION OF NUCLEAR CONTRIBUTION TO GRADIENT
C     AND HELLMANN-FEYNMAN FORCES

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 DNAM
      INTEGER P1,P2,P3
      LOGICAL OUT
      CHARACTER*8 ANAM,ATNAME,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      COMMON/BIG/DRG(NATM1,NATM1),DNUC(3,NATM1)
      COMMON/GRAD/DE(3,NATM1)
      DIMENSION DNAM(3)
      INTEGER*2 I,K,KK,L,N
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DATA DNAM /'E''X','E''Y','E''Z'/

      OUT=NPRINT.EQ.-3.OR.NPRINT.EQ.-10
      DO 10 K=1,NAT
        DNUC(1,K)=ZERO
        DNUC(2,K)=ZERO
10      DNUC(3,K)=ZERO
      DRG(1,1)=ZERO
      DO 30 K=2,NAT
        DRG(K,K)=ZERO
        K1=K-1
        DO 30 L=1,K1
          RKL=(C(1,K)-C(1,L))**2+(C(2,K)-C(2,L))**2+(C(3,K)-C(3,L))**2
          DRG(K,L)=-ONE/RKL
30        DRG(L,K)=DSQRT(RKL)
      DO 60 KK=1,3
        DO 40 K=2,NAT
          ZAK=ZAN(K)
          KM1=K-1
          DO 40 L=1,KM1
            ZAL=ZAN(L)
            PKL=(C(KK,K)-C(KK,L))/DRG(L,K)
            DNUC(KK,K)=DNUC(KK,K)+PKL*DRG(K,L)*ZAK*ZAL
40          CONTINUE
        NAT1=NAT-1
        DO 50 K=1,NAT1
          ZAK=ZAN(K)
          KP1=K+1
          DO 50 L=KP1,NAT
            ZAL=ZAN(L)
            PKL=(C(KK,K)-C(KK,L))/DRG(K,L)
            DNUC(KK,K)=DNUC(KK,K)+PKL*DRG(L,K)*ZAK*ZAL
50          CONTINUE
60      CONTINUE

C     --- PRINT SECTION -----

      DO 70 K=1,NAT
        DO 70 KK=1,3
70        DE(KK,K)=DE(KK,K)+DNUC(KK,K)
      IF(.NOT.OUT)GOTO 110
      WRITE(IW,80)
      MAX=0
90    MIN=MAX+1
      MAX=MAX+5
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(9,140)
      WRITE(9,150) (I,I=MIN,MAX)
      WRITE(9,140)
      DO 100 N=1,3
100     WRITE(9,160) DNAM(N),(DNUC(N,I),I=MIN,MAX)
      IF(MAX.LT.NAT) GOTO 90
110   WRITE(IW,170)
      write(10,170)
      MAX=0
120   MIN=MAX+1
      MAX=MAX+5
      IF(MAX.GT.NAT) MAX=NAT
      write(10,140)
      WRITE(IW,140)
      write(10,150) (I,I=MIN,MAX)
      WRITE(IW,150) (I,I=MIN,MAX)
      write(10,140)
      WRITE(IW,140)
      DO 130 N=1,3
        write(10,160) DNAM(N),(DE(N,I),I=MIN,MAX)
130     WRITE(IW,160) DNAM(N),(DE(N,I),I=MIN,MAX)
      IF(MAX.LT.NAT) GO TO 120
80    FORMAT(//10X,'Nuclear contribution to gradient'//)
140   FORMAT(/)
150   FORMAT(5X,'ATOM',5(6X,I2,6X))
160   FORMAT(7X,A3,5E14.7)
170   FORMAT(/,10X,22(1H-),/,10X,'Hellmann-Feynman force',
     1 /,10X,22(1H-))
      RETURN
      END
