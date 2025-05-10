C==========================================
C  SYMDE  A1
C========================================
C     SYMMETRISE GRADIENT
C==========================================
      SUBROUTINE SYMDE

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      INTEGER*2 I,IC,II,IT,ICNU,JT,N,NP1,NP2,NP3
      INTEGER P1,P2,P3
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      CHARACTER*8 ANAM,ATNAME,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), XTRA(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      CHARACTER*8 GROUP
      COMMON/SYMTRY/T(216),INVT(24),ISO(NSH1+1,24),NT,NT2,
     &XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
     &,U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
     &,INDEX,NAXIS
      COMMON/SYMTRB/GROUP
      COMMON/GRAD/DE(3,NATM1)
      COMMON/BIG/PTR(3,72),ICT(NATM1,24),INDIN(24),XJUNK(3724)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)

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
        NP1=N+1
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
        NP2=N+2
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
        NP3=N+3
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

C     ----- SYMMETRYZE GRADIENT VECTOR -----

      DO 2500 IC=1,NAT
        DO 2100 IT=1,NT
          IF(ICT(IC,IT).GT.IC) GO TO 2500
 2100     CONTINUE
        DEDX=ZERO
        DEDY=ZERO
        DEDZ=ZERO
        DO 2200 IT=1,NT
          ICNU=ICT(IC,IT)
          DEDXP=DE(1,ICNU)
          DEDYP=DE(2,ICNU)
          DEDZP=DE(3,ICNU)
          N=3*(IT-1)
          NP1=N+1
          NP2=N+2
          NP3=N+3
          DEDX=DEDX+DEDXP*PTR(1,NP1)+DEDYP*PTR(2,NP1)+DEDZP*PTR(3,NP1)
          DEDY=DEDY+DEDXP*PTR(1,NP2)+DEDYP*PTR(2,NP2)+DEDZP*PTR(3,NP2)
          DEDZ=DEDZ+DEDXP*PTR(1,NP3)+DEDYP*PTR(2,NP3)+DEDZP*PTR(3,NP3)
 2200     CONTINUE
        DE(1,IC)=DEDX
        DE(2,IC)=DEDY
        DE(3,IC)=DEDZ
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
          NP1=N+1
          NP2=N+2
          NP3=N+3
          DE(1,ICNU)=DE(1,IC)*PTR(1,NP1)+DE(2,IC)*PTR(2,NP1)
     1                                  +DE(3,IC)*PTR(3,NP1)
          DE(2,ICNU)=DE(1,IC)*PTR(1,NP2)+DE(2,IC)*PTR(2,NP2)
     1                                  +DE(3,IC)*PTR(3,NP2)
          DE(3,ICNU)=DE(1,IC)*PTR(1,NP3)+DE(2,IC)*PTR(2,NP3)
     1                                  +DE(3,IC)*PTR(3,NP3)
 2400     CONTINUE
 2500   CONTINUE
      DUM=DBLE(NT)
      DO 2600 N=1,NAT
        DO 2600 I=1,3
 2600     DE(I,N)=DE(I,N)/DUM
      RETURN
      END
