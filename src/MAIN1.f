      SUBROUTINE LIGEN(A,Q,E,IKY,N,NDIM,ITYPE,THRESH)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      DIMENSION A(*),E(*),Q(NDIM,N),IKY(*)
      COMMON/SMALL/Y(NBF1),XJNK(105)
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,TENTH=0.1D0)
      PARAMETER (MS1=-1,M1=1)

      DO 20 I=1,N
        DO 10 J=1,N
10        Q(I,J)=ZERO
20      Q(I,I)=ONE

      IF(N.EQ.M1)GOTO 380
30    TE=ZERO
      DO 40 I=2,N
        I1=I+MS1
        IKYI=IKY(I)
        DO 40 J=1,I1
          TEMP=DABS(A(J+IKYI))
          IF(TE.LT.TEMP)TE=TEMP
40    CONTINUE
      
      IF(TE.LT.THRESH)GOTO 200
      TE=TE*TENTH
      DO 190 I=2,N
        I1=I+MS1
        IP1=I+M1
        IKYI=IKY(I)
        ITEST=N-IP1
        II=I+IKYI
        DO 50 IR=1,N
50        Y(IR)=Q(IR,I)
        DO 170 J=1,I1
          IJ=J+IKYI
          VIJ=A(IJ)
          IF(DABS(VIJ) .LT. TE) GO TO 170
          VII=A(II)*HALF
          J1=J+MS1
          JP1=J+M1
          IKYJ=IKY(J)
          JJ=J+IKYJ
          VJJ=A(JJ)*HALF
          TEMP=VII-VJJ
          TEM=DSQRT(TEMP*TEMP+VIJ*VIJ)
          IF(TEMP)60,70,70
60        TEM=-TEM
70        COST=(TEMP+TEM)/VIJ
          SINT=DSQRT(ONE/(ONE+COST*COST))
          COST=COST*SINT
          TEMP=VII+VJJ
          A(II)=TEMP+TEM
          A(JJ)=TEMP-TEM
          A(IJ)=ZERO
          IF(J1)100,100,80
80        CONTINUE
          DO 90 K=1,J1
            JJ=K+IKYJ
            VIJ=A(K+IKYI)
            A(K+IKYI)=VIJ*COST+A(JJ)*SINT
90          A(JJ)=A(JJ)*COST-VIJ*SINT
100       IF(I1 .LT. JP1) GO TO 120
          DO 110 K=JP1,I1
            JJ=IKY(K)+J
            VIJ=A(K+IKYI)
            A(K+IKYI)=VIJ*COST+A(JJ)*SINT
110         A(JJ)=A(JJ)*COST-VIJ*SINT
120       IF(ITEST)150,130,130
130       CONTINUE
          DO 140 K=IP1,N
            IJ=IKY(K)+I
            JJ=J+IKY(K)
            VIJ=A(IJ)
            A(IJ)=VIJ*COST+A(JJ)*SINT
140         A(JJ)=A(JJ)*COST-VIJ*SINT
150       CONTINUE
          DO 160 K=1,N
            VJJ=Q(K,J)
            Q(K,J)=VJJ*COST-Y(K)*SINT
160         Y(K)=Y(K)*COST+VJJ*SINT
170       CONTINUE
      DO 180 IR=1,N
180     Q(IR,I)=Y(IR)
190   CONTINUE
      GOTO 30
200   CONTINUE

      DO 300 I=1,N
300     E(I)=A(IKY(I)+I)
      GOTO(210,250,370),ITYPE

C     SORT EIGENVALUES INTO INCREASING ORDER

210   CONTINUE
      DO 240 MIN=1,N
        JM=MIN
        EM=E(MIN)
        DO 220 J=MIN,N
          IF(E(J).GE.EM)GOTO 220
          EM=E(J)
          JM=J
220       CONTINUE
        IF(JM.EQ.MIN)GOTO 240
        TEMP=E(JM)
        E(JM)=E(MIN)
        E(MIN)=TEMP
        DO 230 J=1,N
          TEMP=Q(J,JM)
          Q(J,JM)=Q(J,MIN)
230       Q(J,MIN)=TEMP
240     CONTINUE
      GOTO 370

C     SORT INTO DECREASING ORDER

250   CONTINUE
      DO 280 MAX=1,N
        JM=MAX
        EM=E(MAX)
        DO 260 J=MAX,N
          IF(E(J).LE.EM)GOTO 260
          JM=J
          EM=E(J)
260       CONTINUE
        IF(JM.EQ.MAX)GOTO 280
        TEMP=E(JM)
        E(JM)=E(MAX)
        E(MAX)=TEMP
        DO 270 J=1,N
          TEMP=Q(J,JM)
          Q(J,JM)=Q(J,MAX)
270       Q(J,MAX)=TEMP
280   CONTINUE

370   CONTINUE
      RETURN

380   E(1)=A(1)
      RETURN
      END

