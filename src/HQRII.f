      SUBROUTINE HQRII1 (N,IV,LV,IORD,AL,E,NVX,V,SCHMID)

C MODULE HQRII1 TO DIAGONALIZE DENSE REAL SYMMETRIC MATRICES
C Copyright by Annik Vivier Bunge, Carlos F. Bunge, Yoshitaka Beppu,
C              Ichizo Ninomiya and Zdenko A. Tomasic, 1986
C Reference: A.V. Bunge, C.F. Bunge, Comput. Chem., 10, 259, 1986
C Addapted for i8086 code microcomputers by L.A. Montero, J.R. Alvarez,
C Universidad de La Habana, Cuba, 1987.
C Actualized for MS-FORTRAN v. 5.1, Universidad de La Habana, Cuba, 1992
C and CSIC, Instituto de Estructura de la Materia, Madrid, Espana, 1992

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*2 (I-N)
      REAL*8 MACHEP
*      CHARACTER*6 SUBER /'HQRII1'/
      INTEGER*4 I,J,KK
     &,N,IV,LV,IORD,NVX
      LOGICAL SCHMID
      PARAMETER (NX=6,
     .           NX1=(NX*NX+NX)/2,
     .           CERO=0.D0,
     .           P5=0.5D0,
     .           UNO=1.D0,
     .           BOUNDR=1.D-6,
     .           MACHEP=5.471592668859014E-020,
     .           GOLDRT=0.618033988749894D0)

C     MACHEP=1.1610630703531D-16 IN i-8086/7/8 KIND PROCESSORS WITH
C     MS-FORTRAN 77 COMPILER, v. 3.xx

C     MACHEP=5.983858066120757D-20 IN i-80286/7 KIND PROCESSORS WITH
C     MS-FORTRAN 77 COMPILER, v. 4.xx

C     MACHEP=5.471592668859014E-020 IN i-486 KIND PROCESSORS WITH
C     MS-FORTRAN 77 COMPILER, v. 5.1

      DIMENSION AL(*),E(*),V(*),
     &          IX(NX),AU(NX1),W1(NX),W2(NX),W3(NX),W4(NX),W5(NX),W6(NX)

      NL = N
      IVL = IV
      LVL = LV
      NVXL = NVX
      IORDL = IORD
      IF (IVL.LE.0) IVL = 1

*      IF ( NL.LT.1
*     &.OR. NL.GT.NX
*     &.OR.(NL.GT.NVXL.AND.LVL.GE.IVL)
*     &.OR. LVL.GT.NL) THEN
*        CALL HONDER(129,SUBER)
*      ENDIF
      IF (NL.EQ.1) THEN
        E(1) = AL(1)
        V(1) = UNO
        RETURN
      ENDIF

      IX(1) = 0
      DO 5 J=2,NL
5       IX(J) = IX(J-1) + J - 1
      NM1 = NL - 1
      NVF =(LVL-1)*NVXL
      IF (NVF.LT.0) NVF = 0

      IF (NL.GT.2) THEN

C     Householder transformation

        NM2 = NL - 2
        DO 70 K=1,NM2
          KP1 = K+1
          W2(K) = AL(K+IX(K))
          SCALE = CERO
          DO 8 J=KP1,NL
8           SCALE = ABS (AL(IX(J)+K))+SCALE
          W1(K) = AL(IX(KP1)+K)
          IF (SCALE.GT.CERO) THEN
            SCALEI = UNO/SCALE
            SUM = CERO
            DO 10 J=KP1,NL
              W2(J) = AL(IX(J)+K)*SCALEI
10            SUM = W2(J)**2 + SUM
            S = SIGN (SQRT (SUM),W2(KP1))
            W1(K) = -S*SCALE
            W2(KP1) = W2(KP1)+S
            AL(IX(KP1)+K) = W2(KP1)*SCALE
            H = W2(KP1)*S
            HUNS = (H*SCALE)*SCALE
            HI = UNO/H
            SUMM = CERO
            DO 15 II = KP1,NL
15            W5(II) = CERO
            DO 50 I = KP1,NL
              IM1 = I-1
              SUM = CERO
              I0 = IX(I)
              W2I = W2(I)
              NREST = MOD(IM1-KP1+1,6)
              DO 20 J = KP1,KP1+NREST-1
                SUM = W2(J)*AL(I0+J) + SUM
20              W5(J) = W5(J) + W2I*AL(I0+J)
              DO 25 J = KP1+NREST,IM1,6
25              SUM = W2(J  )*AL(I0+J  ) + W2(J+1)*AL(I0+J+1)
     .              + W2(J+2)*AL(I0+J+2) + W2(J+3)*AL(I0+J+3)
     .              + W2(J+4)*AL(I0+J+4) + W2(J+5)*AL(I0+J+5)
     .              + SUM
              DO 30 J=KP1+NREST,IM1,6
                W5(J  ) = W5(J  ) + W2I*AL(I0+J  )
                W5(J+1) = W5(J+1) + W2I*AL(I0+J+1)
                W5(J+2) = W5(J+2) + W2I*AL(I0+J+2)
                W5(J+3) = W5(J+3) + W2I*AL(I0+J+3)
                W5(J+4) = W5(J+4) + W2I*AL(I0+J+4)
30              W5(J+5) = W5(J+5) + W2I*AL(I0+J+5)
50            W6(I) = W2I*AL(I0+I) + SUM
            DO 55 I=KP1,NL
              W1(I) = (W5(I)+W6(I))*HI
55            SUMM = W1(I)*W2(I) + SUMM
            U = P5*SUMM*HI
            DO 65 I=KP1,NL
              I0 = IX(I)
              W1(I) = W2(I)*U - W1(I)
              W1I = W1(I)
              W2I = W2(I)
              NREST = MOD(I-KP1+1,6)
              DO 60 J=KP1,KP1+NREST-1
60              AL(I0+J)= W2I*W1(J) + W1I*W2(J) + AL(I0+J)
              DO 65 J=KP1+NREST,I,6
                AL(I0+J  ) = W2I*W1(J  ) + W1I*W2(J  ) + AL(I0+J  )
                AL(I0+J+1) = W2I*W1(J+1) + W1I*W2(J+1) + AL(I0+J+1)
                AL(I0+J+2) = W2I*W1(J+2) + W1I*W2(J+2) + AL(I0+J+2)
                AL(I0+J+3) = W2I*W1(J+3) + W1I*W2(J+3) + AL(I0+J+3)
                AL(I0+J+4) = W2I*W1(J+4) + W1I*W2(J+4) + AL(I0+J+4)
65              AL(I0+J+5) = W2I*W1(J+5) + W1I*W2(J+5) + AL(I0+J+5)
          ELSE
             HUNS = CERO
          ENDIF
70        AL(IX(K)+K) = HUNS
      ENDIF

      NM1NM1 = IX(NM1) + NM1
      NM1N = IX(NL) + NM1
      NN = NM1N + 1
      W2(NM1) = AL(NM1NM1)
      W2(NL) = AL(NN)
      W1(NM1) = AL(NM1N)
      W1(NL) = CERO
      GERSCH = ABS (W2(1)) + ABS (W1(1))
      DO 90 I=1,NM1
90      GERSCH = MAX (ABS (W2(I+1)) + ABS (W1(I)) + ABS (W1(I+1)),
     .  GERSCH)

C     Trap null matrix before it is too late.

*      IF (GERSCH.EQ.CERO) THEN
*        CALL HONDER(33,SUBER)
*      ENDIF

      SUMD = CERO
      SUMCOD = CERO
      DO 100 I=1,NL
        SUMCOD = ABS (W1(I)) + SUMCOD
100     SUMD = ABS (W2(I)) + SUMD
      SCALE = SUMD + SUMCOD
      SCALEI = UNO / SCALE
      DO 105 I = 1,NL
        W1(I) = W1(I)*SCALEI
        W2(I) = W2(I)*SCALEI
        W3(I) = W1(I)
        E(I) = W2(I)
105     V(I+NVF) = E(I)
      EPS = SQRT (MACHEP)
      GERSCH = GERSCH*SCALEI
      DEL = GERSCH*EPS
      DELW5 = GERSCH*MACHEP

C     CHANGE SUGGESTED BY C. BUNGE, OCTOBER 21th., 1987
*      IF (SUM/SUMCOD.GT.DEL) THEN

C     QR method with origin shift.

      DO 160 K=NL,2,-1
110     L = K
120     IF (ABS (W3(L-1)).GT.DEL) THEN
          L = L - 1
          IF (L.GT.1) GO TO 120
        ENDIF
        IF (L.NE.K) THEN
          WW = (E(K-1)+E(K))*P5
          R = E(K)-WW
          Z = WW - SIGN (SQRT (W3(K-1)**2 + R*R),WW)
          EE = E(L) - Z
          E(L) = EE
          FF = W3(L)
          R = SQRT (EE*EE + FF*FF)
          RI = UNO/R
          C = E(L)*RI
          S = W3(L)*RI
          WW = E(L+1) - Z
          E(L) = (FF*C + WW*S)*S + EE + Z
          E(L+1) = C*WW - S*FF
          DO 150 J=L+1,K-1
            R = SQRT (E(J)**2 + W3(J)**2)
            RI = UNO/R
            W3(J-1) = S*R
            EE = E(J)*C
            FF = W3(J)*C
            C = E(J)*RI
            S = W3(J)*RI
            WW = E(J+1) - Z
            E(J) = (FF*C + WW*S)*S + EE + Z
            E(J+1) = C*WW - S*FF
150         CONTINUE
          W3(K-1) = E(K)*S
          E(K) = E(K)*C + Z
          GO TO 110
        ENDIF
160     CONTINUE

C     Straight selection sort of eigenvalues.

      SORTER = UNO
      IF (IORDL.LT.0) SORTER = -UNO
      J = NL
170   L = 1
      II = 1
      LL = 1
      DO 190 I=2,J
        IF ((E(I)-E(L))*SORTER.LE.CERO) THEN
          L = I
        ELSE
          II = I
          LL = L
        ENDIF
190     CONTINUE
      IF (II.NE.LL) THEN
        WW = E(LL)
        E(LL) = E(II)
        E(II) = WW
      ENDIF
        J = II - 1
        IF (J.GT.1) GO TO 170
C     CHANGE SUGGESTED BY C. BUNGE, OCTOBER 21th., 1987
*      ENDIF

      IF (LVL.GE.IVL) THEN

C     Inverse iteration for eigenvectors.

        FN = DBLE (NL)
        EPS1 = SQRT (FN)*EPS
        SEPS = SQRT(EPS)
        EPS2 = (GERSCH*BOUNDR)/(FN*SEPS)
        RN = CERO
        RA = EPS*GOLDRT
        IF = (IVL-2)*NVXL
        DO 312 I = IVL,LVL
          IF = IF + NVXL
          DO 220 J=1,NL
            W3(J) = CERO
            W4(J) = W1(J)
            W5(J) = V(NVF+J) - E(I)
            RN = RN + RA
            IF (RN.GE.EPS) RN = RN - EPS
220         W6(J) = RN
          DO 250 J =1,NM1
            IF (ABS (W5(J)).LE.ABS(W1(J))) THEN
              IF (W1(J).EQ.CERO) W1(J) = DEL
              W2(J) = -W5(J)/W1(J)
              W5(J) = W1(J)
              T = W5(J+1)
              W5(J+1) = W4(J)
              W4(J) = T
              W3(J) = W4(J+1)
              IF (W3(J).EQ.CERO) W3(J) = DEL
              W4(J+1) = CERO
            ELSE
              W2(J) = -W1(J)/W5(J)
            ENDIF
            W4(J+1) = W3(J)*W2(J) + W4(J+1)
250         W5(J+1) = W4(J)*W2(J) + W5(J+1)
          IF (W5(NL).EQ.CERO) W5(NL) = DELW5
          WNM15I = UNO/W5(NM1)
          WN5I = UNO/W5(NL)
          DO 310 ITERE=1,2
            IF (ITERE.NE.1) THEN
              DO 260 J=1,NM1
                IF (W5(J).EQ.W1(J)) THEN
                  T = W6(J)
                  W6(J) = W6(J+1)
                  W6(J+1) = T
                ENDIF
260             W6(J+1) = W6(J)*W2(J) + W6(J+1)
            ENDIF
            W6(NL) = W6(NL)*WN5I
            W6(NM1) = (W6(NM1)-W6(NL)*W4(NM1))*WNM15I
            VN = MAX (ABS (W6(NL)), ABS (W6(NM1)))
            DO 280 K = NM2,1,-1
              W6(K) = (W6(K)-W6(K+1)*W4(K)-W6(K+2)*W3(K))/W5(K)
280           VN =MAX(ABS(W6(K)),VN)
            S = EPS1/VN
            NREST = MOD(NL,6)
            DO 295 J=1,NREST
295           W6(J) = S*W6(J)
            DO 300 J=1+NREST,NL,6
              W6(J  ) = S*W6(J  )
              W6(J+1) = S*W6(J+1)
              W6(J+2) = S*W6(J+2)
              W6(J+3) = S*W6(J+3)
              W6(J+4) = S*W6(J+4)
300           W6(J+5) = S*W6(J+5)
310         CONTINUE
          DO 312 J=1,NL
312         V(IF+J) = W6(J)

C     Use sequential row notation for upper triangular matrix

        IX(1) = 0
        DO 315 J=2,NL
315       IX(J) = IX(J-1) - J + 1 + NL
        IJ = 0
        DO 320 J=1,NL
          DO 320 I=1,J
            IJ = IJ + 1
320         AU(IX(I)+J) = AL(IJ)

C     Back transformation of eigenvectors.

        IG = 1
        IF = (IVL-2)*NVXL
        DO 430  I=IVL,LVL
          IF = IF + NVXL
          DO 322 J=1,NL
322         W6(J) = V(IF+J)
          IM1 = I - 1
          IF (NL.GT.2) THEN
            DO 350 J=1,NM2
              K = NL - J - 1
              K0 = IX(K)
              IF (AU(K0+K).NE.CERO) THEN
                KP1 = K + 1
                SUM = CERO
                NREST = MOD(NL-KP1+1,6)
                DO 325 KK=KP1,KP1+NREST-1
325               SUM = AU(K0+KK)*W6(KK) + SUM
                DO 330 KK=KP1+NREST,NL,6
330               SUM = AU(K0+KK  )*W6(KK  )
     .                + AU(K0+KK+1)*W6(KK+1)
     .                + AU(K0+KK+2)*W6(KK+2)
     .                + AU(K0+KK+3)*W6(KK+3)
     .                + AU(K0+KK+4)*W6(KK+4)
     .                + AU(K0+KK+5)*W6(KK+5)+SUM
                S = -SUM/AU(K0+K)
                DO 335 KK=KP1,KP1+NREST-1
335               W6(KK) = S*AU(K0+KK) + W6(KK)
                DO 340 KK=KP1+NREST,NL,6
                  W6(KK  ) = S*AU(K0+KK  ) + W6(KK  )
                  W6(KK+1) = S*AU(K0+KK+1) + W6(KK+1)
                  W6(KK+2) = S*AU(K0+KK+2) + W6(KK+2)
                  W6(KK+3) = S*AU(K0+KK+3) + W6(KK+3)
                  W6(KK+4) = S*AU(K0+KK+4) + W6(KK+4)
340               W6(KK+5) = S*AU(K0+KK+5) + W6(KK+5)
              ENDIF
350           CONTINUE
          ENDIF
          DO 370 J=IG,I
            JJ = J
370         IF (ABS(E(J)-E(I)).LT.EPS2) GO TO 380
380       IG = JJ
          NREST = MOD(NL,6)

          IF (IG.NE.I.AND.SCHMID) THEN

C         Degenerate eigenvalues. First, orthogonalize.

            KF = (IG-2)*NVXL
            DO 400 K=IG,IM1
              KF = KF + NVXL
              SUM = CERO
              DO 394 J=1,NREST
394             SUM = V(KF+J)*W6(J) + SUM
              DO 395 J=1+NREST,NL,6
395             SUM = V(KF+J  )*W6(J  ) + V(KF+J+1)*W6(J+1)
     .              + V(KF+J+2)*W6(J+2) + V(KF+J+3)*W6(J+3)
     .              + V(KF+J+4)*W6(J+4) + V(KF+J+5)*W6(J+5)
     .              + SUM
              S = -SUM
              DO 398 J=1,NREST
398             W6(J) = S*V(KF+J) + W6(J)
              DO 399 J=1+NREST,NL,6
                W6(J  ) = S*V(KF+J  ) + W6(J  )
                W6(J+1) = S*V(KF+J+1) + W6(J+1)
                W6(J+2) = S*V(KF+J+2) + W6(J+2)
                W6(J+3) = S*V(KF+J+3) + W6(J+3)
                W6(J+4) = S*V(KF+J+4) + W6(J+4)
399             W6(J+5) = S*V(KF+J+5) + W6(J+5)
400           CONTINUE
          ENDIF

C         Normalization.

          SUM = CERO
          DO 405 J=1,NREST
405         SUM = W6(J)**2 + SUM
          DO 410 J=1+NREST,NL,6
410         SUM = W6(J  )**2 + W6(J+1)**2 + W6(J+2)**2
     .          + W6(J+3)**2 + W6(J+4)**2 + W6(J+5)**2 + SUM
          IF (SUM.EQ.CERO) SUM=MACHEP
          S = UNO/SQRT(SUM)
          DO 415 J=1,NREST
415         V(IF+J) = S*W6(J)
          DO 420 J=1+NREST,NL,6
            V(IF+J  ) = S*W6(J  )
            V(IF+J+1) = S*W6(J+1)
            V(IF+J+2) = S*W6(J+2)
            V(IF+J+3) = S*W6(J+3)
            V(IF+J+4) = S*W6(J+4)
420         V(IF+J+5) = S*W6(J+5)
430       CONTINUE
      ENDIF

      DO 450 I=1,NL
450     E(I) = SCALE*E(I)

      RETURN
      END
