C     ****************************************
C
C     NUMERICAL FORCE CONSTANT ROUTINES
C
C     HOPEFULLY MACHINE INDEPENDENT.....
C
C     NO CHANGES REQUIRED FOR CRAY
C
C    *****************************************
      SUBROUTINE FGMTRX(IW)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
C
C     ----- SUBROUTINE CONSTRUCTS AN EFFECTIVE WILSON -FG- MATRIX
C           BY INTRODUCING MASS DEPENDENCE INTO THE POTENTIAL FORCE
C           CONSTANT MATRIX.  THE NORMAL MODES AND SPECTROSCOPIC
C           FREQUENCIES ARE THEN EVALUATED.                     -----
C
      INTEGER P1, P2, P3
      CHARACTER*4 CLAB
      CHARACTER*8 ANAM,ATNAME,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(NATM13),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),
     &             PARNAM(NATM13),ATNAME(NATM1)
      COMMON/OPTMIZ/SPACE(4*NATM13+2),NCOORD,IATOM,ICOORD,IVIB
      DIMENSION CLAB(3)
      COMMON/MAPPER/IA(NBF1+1)
      COMMON/BIG/VEC(NATM13,NATM13),
     &           A(NNC1),FC(NNC1),E(NATM13),RM(NATM13)
      common/thermo/cpt,cpr,cpv,cpe,
     &              h25t,h25r,h25v,h25e,
     &              s25t,s25r,s25v,s25e,
     &              g25t,g25r,g25v,g25e,
     &              ilin,nax,
     &              temper
      PARAMETER (ZERO=0.0D0, TFACT=2.6436411815482D+07)
      parameter (
     &           R=8.3144126d0,
     &           A298=298.15d0,
     &           UACM=27.211606d0*8065.479d0,
     &           C2=1.43878645d0,
     &           UNO=1.d0, DOS=2.d0)
      DATA CLAB/'X','Y','Z'/

      NNC=(NCOORD*(NCOORD+1))/2
      DO 10 I=1,NCOORD
        DO 10 J=1,I
10        FC(IA(I)+J)=VEC(I,J)

C     ----- DIAGONALIZE CARTESIAN FORCE MATRIX -----

      DO 20 I=1,NNC
20      A(I)=FC(I)
      NDIM=NATM13
*      call hqrii1 (NCOORD,1,NCOORD,-1,A,E,NDIM,VEC,.TRUE.)
      CALL LIGEN(A,VEC,E,IA,NCOORD,NDIM,1,1.0D-20)
      WRITE(IW,50) CHAR(12)
      write(10 ,50) ' '
      MAX=0
30    MIN=MAX+1
      MAX=MAX+5
      IF(MAX.GT.NCOORD) MAX=NCOORD
      WRITE(IW,60)
      write(10 ,60)
      WRITE(IW,70) (J,J=MIN,MAX)
      write(10 ,70) (J,J=MIN,MAX)
      WRITE(IW,60)
      write(10 ,60)
      WRITE(IW,80) (E(J),J=MIN,MAX)
      write(10 ,80) (E(J),J=MIN,MAX)
      WRITE(IW,60)
      write(10 ,60)
      DO 40 IAT=1,NAT
        I0=3*(IAT-1)
        WRITE(IW,90) IAT,ANAM(IAT),CLAB(1),
     1               (VEC(I0+1,J),J=MIN,MAX)
        write(10 ,90) IAT,ANAM(IAT),CLAB(1),
     1               (VEC(I0+1,J),J=MIN,MAX)
        WRITE(IW,100) CLAB(2),(VEC(I0+2,J),J=MIN,MAX)
        write(10 ,100) CLAB(2),(VEC(I0+3,J),J=MIN,MAX)
        WRITE(IW,100) CLAB(3),(VEC(I0+3,J),J=MIN,MAX)
        write(10 ,100) CLAB(3),(VEC(I0+3,J),J=MIN,MAX)
40      CONTINUE
      IF(MAX.LT.NCOORD) GO TO 30
50    FORMAT(A/,10X,47(1H-),/,10X,'EIGENVECTORS OF CARTESIAN',
     1 ' FORCE CONSTANT MATRIX',/,10X,47(1H-))
60    FORMAT(/)
70    FORMAT(20X,5(4X,I3,5X))
80    FORMAT(20X,5F12.5)
90    FORMAT(I3,3X,A8,2X,A4,5F12.8)
100   FORMAT(16X,A4,5F12.8)

C     ----- CREATE FULL FORCE MATRIX WITH MASS WEIGHTED TRANSFORMATION -

110   CONTINUE
      CALL RAMS(RM,IW)
      DO 120 I=1,NCOORD
        DO 120 J=1,I
          IJ=IA(I)+J
120       A(IJ)=RM(I)*FC(IJ)*RM(J)

C     ----- GET NORMAL MODES AND FREQUENCIES -----

*      call hqrii1 (NCOORD,1,NCOORD,-1,A,E,NDIM,VEC,.TRUE.)
      CALL LIGEN(A,VEC,E,IA,NCOORD,NDIM,1,1.0D-20)

C     ----- CONVERT FREQUENCIES TO INVERSE CM -----

      DO 130 I=1,NCOORD
130     E(I)=DSIGN(1.0D0,E(I))*DSQRT(DABS(TFACT*E(I)))

C     ----- CONSTRUCT MASS WIEGHTED DISPLACEMENT VECTORS
C           AND RENORMALIZE THEM -----

      DO 140 I=1,NCOORD
        DO 140 J=1,NCOORD
140       VEC(J,I)=VEC(J,I)*RM(J)
      DO 160 J=1,NCOORD
        DUM=ZERO
        DO 150 I=1,NCOORD
150       DUM=DUM+VEC(I,J)*VEC(I,J)
        DUM=DSQRT(DUM)
        DO 160 I=1,NCOORD
160       VEC(I,J)=VEC(I,J)/DUM
      WRITE(IW,190)
      write(10 ,190)
      MAX=0
170   MIN=MAX+1
      MAX=MAX+5
      IF(MAX.GT.NCOORD) MAX=NCOORD
      WRITE(IW,60)
      write(10 ,60)
      WRITE(IW,70) (J,J=MIN,MAX)
      write(10,70) (J,J=MIN,MAX)
      WRITE(IW,60)
      write(10 ,60)
      WRITE(IW,80) (E(J),J=MIN,MAX)
      write(10,80) (E(J),J=MIN,MAX)
      WRITE(IW,60)
      write(10 ,60)
      DO 180 IAT=1,NAT
        I0=3*(IAT-1)
        WRITE(IW,90) IAT,ANAM(IAT),CLAB(1),
     1               (VEC(I0+1,J),J=MIN,MAX)
        write(10 ,90) IAT,ANAM(IAT),CLAB(1),
     1               (VEC(I0+1,J),J=MIN,MAX)
        WRITE(IW,100) CLAB(2),(VEC(I0+2,J),J=MIN,MAX)
        write(10 ,100) CLAB(2),(VEC(I0+2,J),J=MIN,MAX)
        WRITE(IW,100) CLAB(3),(VEC(I0+3,J),J=MIN,MAX)
        write(10 ,100) CLAB(3),(VEC(I0+3,J),J=MIN,MAX)
180     CONTINUE
      IF(MAX.LT.NCOORD) GO TO 170
190   FORMAT(/10X,40(1H-)/10X,'NORMAL MODES AND VIBRATIONAL',
     1 ' FREQUENCIES',/,10X,40(1H-))

C     THERMODYNAMIC PROPERTIES

      if (ilin.eq.0) then
        ninic = 7
      else
        ninic = 6
      endif
      cpv = ZERO
      h25v = ZERO
      s25v = ZERO
      g25v = ZERO
      rt = R*temper
      do i=ninic,ncoord
        uf = .9d0*C2*e(i)/temper
        euf = exp(-uf)
        vv1 = (UNO-euf)
        cpv  =  cpv + R*uf**2*euf/(vv1**2)
        vv2 = uf*euf/vv1
        h25v =  h25v + rt*vv2
        s25v =  s25v + R*vv2 - R*log(vv1)
        g25v =  g25v - rt*log(vv1)
      enddo
      write(9,11) temper,cpv,h25v,s25v,g25v
      write(10,11) temper,cpv,h25v,s25v,g25v
11    format (//' THERMODYNAMIC QUANTITIES'/
     &          ' Calculated for the output geometry, as an ideal gas (T
     &=',f10.2,' K)'//
     &         ' Vibrational:'/
     &         ' - Heat capacity [Cp]         =',g16.7,' J/ K mol'/
     &         ' - Enthalpy [Ho(25)-Ho(0)]    =',g16.7,' J/ mol'/
     &         ' - Entropy [So]               =',g16.7,' J/ K mol'/
     &         ' - Free energy [Go(25)-Ho(0)] =',g16.7,' J/ mol')
      if (MUL.gt.1) then
        cpe = ZERO
        h25e = ZERO
        alogm = log(dble(MUL))
        s25e = R*alogm
        g25e = -s25e*temper + h25e
        write(10,13) cpe, h25e, s25e, g25e
        write(iw,13) cpe, h25e, s25e, g25e
13      format (/
     &         ' Electronic (only for ground state):'/
     &         ' - Heat capacity [Cp]         =',g16.7,' J/ K mol'/
     &         ' - Enthalpy [Ho(25)-Ho(0)]    =',g16.7,' J/ mol'/
     &         ' - Entropy [So]               =',g16.7,' J/ K mol'/
     &         ' - Free energy [Go(25)-Ho(0)] =',g16.7,' J/ mol')
      else
        cpe = ZERO
        h25e = ZERO
        s25e = ZERO
        g25e = ZERO
        write(iw,14)
14      format (/' *** Ground state singlets have no electronic contribu
     &tion to thermodynamic'/
     &           ' quantities.')
      endif
      if (cpr.ne.ZERO) then
        write(9,12)
     &     cpv+cpr+cpt+cpe,
     &     h25v+h25r+h25t+h25e,
     &     s25v+s25r+s25t+s25e,
     &     g25v+g25r+g25t+g25e
        write(10,12)
     &     cpv+cpr+cpt+cpe,
     &     h25v+h25r+h25t+h25e,
     &     s25v+s25r+s25t+s25e,
     &     g25v+g25r+g25t+g25e
12      format (/
     &         ' TOTAL:'/
     &         ' - Heat capacity [Cp]         =',g16.7,' J/ K mol'/
     &         ' - Enthalpy [Ho(25)-Ho(0)]    =',g16.7,' J/ mol'/
     &         ' - Entropy [So]               =',g16.7,' J/ K mol'/
     &         ' - Free energy [Go(25)-Ho(0)] =',g16.7,' J/ mol')
      endif
      RETURN
      END
