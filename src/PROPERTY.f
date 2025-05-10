      SUBROUTINE MULKEN

C     MULLIKEN CHARGE ANALYSIS

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL LADAPT,LSKIP
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),KR(NATM1),
     5               KA(NATM1),KB(NATM1),
     4               PARVAL(NATM13),ZIN(NATM1),X(3,NATM1)
      CHARACTER*8  ALAB,BFLAB,PARNAM,ATNAME
      COMMON/INFOB/ALAB(NATM1),BFLAB(NBF1),PARNAM(NATM13),
     1               ATNAME(NATM1)
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),KTYPE(NSH1),
     3             KNG(NSH1),KLOC(NSH1)
     2            ,KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2
      COMMON/SMALL/Y(NATM1*NATM1+1)
      COMMON/MAPPER/IA(NBF1+1)
      COMMON/BIG/Q(NBF1M),Q1(NBF1M),LIMLOW(NATM1),LIMSUP(NATM1),
     2 APOP(NATM1)
      COMMON/BIGG/QQ(MAXDIM),QQ1(NBF1M)
      LOGICAL LMOS,LNUCL
      COMMON/PRPSEC/ISECD1,ITYPD1,ISECD2,ITYPD2,LMOS,ISECV,ITYPV
     1     ,ISECD3,ITYPD3,LNUCL,GX,GY,GZ
      CHARACTER*8 OSCF,CLOSED
      DATA OSCF/'OSCF'/
      DATA CLOSED/'CLOSED'/,ZERO/0.0D0/

1     FORMAT(A//5X,14(5H*****),/,
     1 23X,'ELECTRONIC PROPERTY CALCULATIONS',/,5X,14(5H*****))

      WRITE(IW,1) CHAR(12)
      write(10 ,1) ' '
      IPASS=1
      I1=1
      I2=1
      I3=1
      CALL SECGET(ISEX(5),5,ISECS)
      WRITE(9,'(/a,i2,a)')
     & ' PROPERTIES FROM DENSITY MATRIX ON SECTION ',
     & ISECD1,' IN DUMP FILE'
      CALL SECGET(ISECD1,ITYPD1,IBLD1)
      CALL READQ(Q(1),NX,IBLD1,IFILD)
      IF(SCFTYP.NE.CLOSED) THEN
        CALL SECGET(ISECD2,ITYPD2,IBLD2)
        CALL READQ(Q1(1),NX,IBLD2,IFILD)
        DO 10 I=1,NX
10        Q(I)=Q(I)+Q1(I)
      ENDIF
20    CALL READQ(Q1(1),NX,ISECS,IFILD)
      IF(LMOS)THEN
        CALL SECGET(ISECV,ITYPV,IBLKV)
        CALL READQ(qq(1),NUM*NCOORB,IBLKV,IFILD)
        CALL DTR(Q(1),qq1(1),qq(1),Y,NUM,NCOORB,NUM)
        CALL FMOVE(qq1(1),Q(1),NX)
        CALL DTR(Q1(1),qq1(1),qq(1),Y,NUM,NCOORB,NUM)
        CALL FMOVE(qq1(1),Q1(1),NX)
      ENDIF
      DO 40 I=1,NX
40      Q(I)=Q(I)*Q1(I)
      DO 60 I=1,NUM
        DUM=ZERO
        DO 50 J=1,NUM
          IJ=IA(I)+J
          IF(J.GT.I) IJ=IA(J)+I
50      DUM=DUM+Q(IJ)
60    qq1(I)=DUM
      DO 70 I=1,NAT
        LIMLOW(I)=0
70    LIMSUP(I)=0
      LAT=KATOM(1)
      LIMLOW(LAT)=1
      DO 80 I=1,NSHELL
        IAT=KATOM(I)
        IF(LAT.EQ.IAT) GOTO 80
        LIMSUP(LAT)=KLOC(I)-1
        LAT=IAT
        LIMLOW(LAT)=KLOC(I)
80    CONTINUE
      LIMSUP(LAT)=NUM
      DO 100 I=1,NAT
        APOP(I)=ZERO
        L1=LIMLOW(I)
        M1=LIMSUP(I)
        IF(L1.EQ.0.OR.M1.EQ.0) GOTO 100
        DO 90 J=L1,M1
90        APOP(I)=APOP(I)+qq1(J)
100   CONTINUE
      IF(IPASS.EQ.1) then
        WRITE(IW,110)
        WRITE(10,110)
110     FORMAT(/' CHARGE DENSITIES'/)
      endif
      IF(IPASS.EQ.2) then
        WRITE(10,120)
        IF(IPASS.EQ.2)WRITE(IW,120)
120     FORMAT(/' SPIN DENSITIES'/)
      endif
      WRITE(10,130) (I,BFLAB(I),qq1(I),I=1,NUM)
      WRITE(IW,130) (I,BFLAB(I),qq1(I),I=1,NUM)
130   FORMAT(10X,'BASIS FUNCTION POPULATIONS'
     1//(10X,I5,2X,A8,F12.5))
      WRITE(10,140) (I,ALAB(I),APOP(I),I=1,NAT)
      WRITE(IW,140) (I,ALAB(I),APOP(I),I=1,NAT)
140   FORMAT(/10X,'MULLIKEN GROSS ATOMIC POPULATIONS'
     1//(10X,I5,3x,A8,F12.5))
      IF(IPASS.EQ.2)RETURN

      IF(SCFTYP.EQ.OSCF)GOTO 170
      RETURN
170   CONTINUE
      CALL READQ(Q(1),NX,IBLD2,IFILD)
      IPASS=2
      GOTO 20
      END
C--------------------------------------------------------------------
      SUBROUTINE DIPOLE
C
C----------------------------------------------------------------
C     DIPOLE MOMENTS
C     R D  AMOS
C----------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 CLOSED
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      CHARACTER*8  ANAM,BFLAB,PARNAM,ATNAME
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),KR(NATM1),
     5                KA(NATM1),KB(NATM1),
     4               PARVAL(NATM13),ZIN(NATM1),X(3,NATM1)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),
     1               ATNAME(NATM1)
      LOGICAL LADAPT,LSKIP
      LOGICAL LMOS,LNUCL
      COMMON/SMALL/Y(NATM1*NATM1+1)
      COMMON/PRPSEC/ISECD1,ITYPD1,ISECD2,ITYPD2,LMOS,ISECV,ITYPV
     1  ,ISECD3,ITYPD3,LNUCL,ORIGIN(3)
      COMMON/MAXLEN/MAXQ
      COMMON/BIG/Q(NBF1M),Q1(NBF1M),QXTRA(2*NATM1)
      COMMON/BIGG/QQ(MAXDIM),QQ1(NBF1M)
      COMMON/MAPPER/IA(NBF1+1)
C
      DIMENSION DE(3),DN(3),DT(3)
      DATA ZERO,FAC,CLOSED/0.0D0,2.5417701D0,'CLOSED'/

      IF(ICH.NE.0) WRITE(IW,10)
      IF(ICH.NE.0) WRITE(10,10)
10    FORMAT(//10X,'Warning ---- dipole moment for charged species',
     1' is origin dependent')
C     ----- CALCULATE DIPOLE MOMENT INTEGRALS -----
      WRITE(IW,15) ORIGIN
      WRITE(10,15) ORIGIN
15    FORMAT(//
     &' DIPOLE MOMENT'//
     &' Origin for property evaluation (X,Y,Z in Bohrs):',3f9.5)

      IBLOK=0
      I1=1
      I2=I1
      I3=NX+I2
      IF(ISECD1.EQ.0)ISECD1=ISEX(7)
         CALL SECGET(ISECD1,ITYPD1,IBLOK)
         CALL READQ(Q(1),NX,IBLOK,IFILD)
      IF(SCFTYP.EQ.CLOSED) GOTO 30
      IF(ISECD2.EQ.0)ISECD2=ISEX(10)
      CALL SECGET(ISECD2,ITYPD2,IBLOK)
      CALL READQ(Q1(1),NX,IBLOK,IFILD)
      DO 20 I=1,NX
20    Q(I)=Q(I)+Q1(I)
30    IF(.NOT.LMOS)GOTO 40
      CALL SECGET(ISECV,ITYPV,IBLOK)
      CALL READQ(QQ(1),NUM*NCOORB,IBLOK,IFILD)
      CALL DTR(Q(1),Q1(1),QQ(1),Y,NUM,NUM,NUM)
      CALL FMOVE (Q1(1),Q(1),NX)
40    CONTINUE
      CALL XYZDIP(Q1(1),QQ(1),QQ(I3))
C     ----- ELECTRONIC CONTRIBUTION TO DIPOLE MOMENT -----
      DE(1)=-TRACEP(Q(1),Q1(1),NUM)
      DE(2)=-TRACEP(Q(1),QQ(1),NUM)
      DE(3)=-TRACEP(Q(1),QQ(I3),NUM)
C     ----- NUCLEAR CONTRIBUTION
      DO 70 J=1,3
      DN(J)=ZERO
      IF(.NOT.LNUCL)GOTO 60
      DO 50 I=1,NAT
50    DN(J)=DN(J)+ZAN(I)*(C(J,I)-ORIGIN(J))
60    CONTINUE
70    CONTINUE
      DO 80 J=1,3
80    DT(J)=DE(J)+DN(J)
      DIPOL=DSQRT(DT(1)*DT(1)+DT(2)*DT(2)+DT(3)*DT(3))
      WRITE(IW,90)
      WRITE(10,90)
      WRITE(10,100) DE,DN,DT,DIPOL
      WRITE(IW,100) DE,DN,DT,DIPOL
90    FORMAT(/30X,'In relative atomic units:')
100   FORMAT(/30X,'x',15X,'y',15X,'z'
     1//10X,'Electronic',3(2X,F14.8)
     2/10X,'Nuclear',3X,3(2X,F14.8)
     3/10X,'Total',5X,3(2X,F14.8)
     4//10X,'Dipole Moment',T40,F15.8)
      DIPOL=FAC*DIPOL
      DO 110 I=1,3
      DE(I)=DE(I)*FAC
      DN(I)=DN(I)*FAC
110   DT(I)=DT(I)*FAC
      WRITE(IW,120)
      WRITE(IW,100) DE,DN,DT,DIPOL
      WRITE(10,120)
      WRITE(10,100) DE,DN,DT,DIPOL
120   FORMAT(/30X,'In Debyes')
      RETURN
      END
C------------------------------------------------------------------
C
      SUBROUTINE QPOLE
C-------------------------------------------------------------------
C     QUADRUPOLE MOMENTS     R D AMOS
C------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL LADAPT,LSKIP
      CHARACTER*8 CLOSED
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),KR(NATM1),
     5               KA(NATM1),KB(NATM1),
     4               PARVAL(NATM13),ZIN(NATM1),X(3,NATM1)
      CHARACTER*8  ANAM,BFLAB,PARNAM,ATNAME
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),
     1               ATNAME(NATM1)
      COMMON/MAXLEN/MAXQ
      COMMON/BIG/Q(NBF1M),Q1(NBF1M),QXTRA(2*NATM1)
      COMMON/BIGG/QQ(MAXDIM),QQ1(NBF1M)
      COMMON/BIGK/QQ2(NBF1M)
      COMMON/SMALL/Y(NATM1*NATM1+1)
      LOGICAL LMOS,LNUCL
      COMMON/PRPSEC/ISECD1,ITYPD1,ISECD2,ITYPD2,LMOS,ISECV,ITYPV
     1     ,ISECD3,ITYPD3,LNUCL,GX,GY,GZ
      COMMON/MAPPER/IA(NBF1+1)
C
      DIMENSION SECM(6),SECN(6),SECT(6)
      DATA ZERO,HALF,THRHAL/0.0D0,0.5D0,1.5D0/
      DATA CLOSED/'CLOSED'/
      DATA FAC/1.34504627D0/
C     ----- CALCULATE SECOND MOMENT INTEGRALS -----
      LEN=LENSEC(NX)
      IBLOK=IBLKS
      IPOS=1
      I1=1
      I2=1
      I3=NX+I2
      I4=1
      I5=1
      CALL QUADRU(Q(1),Q1(1),QQ(1),QQ(I3),QQ1(1),QQ2(1))
      CALL WRT3(Q(1),NX,IBLOK,IFILS)
      IBLOK=IBLOK+LEN
      CALL WRT3(Q1(1),NX,IBLOK,IFILS)
      IBLOK=IBLOK+LEN
      CALL WRT3(QQ(1),NX,IBLOK,IFILS)
      IBLOK=IBLOK+LEN
      CALL WRT3(QQ(I3),NX,IBLOK,IFILS)
      IBLOK=IBLOK+LEN
      CALL WRT3(QQ1(1),NX,IBLOK,IFILS)
      IBLOK=IBLOK+LEN
      CALL WRT3(QQ2(1),NX,IBLOK,IFILS)
      IBLOK=IBLOK+LEN
10    CONTINUE
         CALL SECGET(ISECD1,ITYPD1,IBLOK)
         CALL READQ(Q(1),NX,IBLOK,IFILD)
      IF(SCFTYP.EQ.CLOSED) GOTO 30
      CALL SECGET(ISECD2,ITYPD2,IBLOK)
      CALL READQ(Q1(1),NX,IBLOK,IFILD)
      DO 20 I=1,NX
20    Q(I)=Q(I)+Q1(I)
30    IF(.NOT.LMOS)GOTO 40
      CALL SECGET(ISECV,ITYPV,IBLOK)
      CALL READQ(QQ(1),NUM*NCOORB,IBLOK,IFILD)
      CALL DTR(Q(1),Q1(1),QQ(1),Y,NUM,NCOORB,NUM)
      CALL FMOVE(Q1(1),Q(1),NX)
40    CONTINUE
C     ----- ELECTRONIC CONTRIBUTION TO SECOND MOMENT -----
      IBLOK=IBLKS
      DO 50 I=1,6
      CALL READQ(Q1(1),NX,IBLOK,IFILS)
      IBLOK=IBLOK+LEN
      SECM(I)=-TRACEP(Q(1),Q1(1),NUM)
50    CONTINUE
C     ----- NUCLEAR CONTRIBUTION
      DO 60 I=1,6
60    SECN(I)=ZERO
      IF(.NOT.LNUCL)GOTO 80
      DO 70 I=1,NAT
      SECN(1)=SECN(1)+ZAN(I)*(C(1,I)-GX)*(C(1,I)-GX)
      SECN(2)=SECN(2)+ZAN(I)*(C(2,I)-GY)*(C(2,I)-GY)
      SECN(3)=SECN(3)+ZAN(I)*(C(3,I)-GZ)*(C(3,I)-GZ)
      SECN(4)=SECN(4)+ZAN(I)*(C(1,I)-GX)*(C(2,I)-GY)
      SECN(5)=SECN(5)+ZAN(I)*(C(1,I)-GX)*(C(3,I)-GZ)
      SECN(6)=SECN(6)+ZAN(I)*(C(2,I)-GY)*(C(3,I)-GZ)
70    CONTINUE
80    DO 90 I=1,6
      SECT(I)=SECM(I)+SECN(I)
90    CONTINUE
      WRITE(10,100)
      WRITE(IW,100)
100   FORMAT(/5X,'SECOND MOMENTS AND QUADRUPOLE MOMENTS ARE',/,
     1        5x,'CALCULATED RELATIVE TO THE ORIGIN OF THE'/
     &        5x,'COORDINATE FRAME')
      WRITE(10,110)
      WRITE(IW,110)
110   FORMAT(//20X,'SECOND MOMENTS -- IN ATOMIC UNITS')
      WRITE(10,120) SECM,SECN,SECT
      WRITE(IW,120) SECM,SECN,SECT
120   FORMAT(/19X,'XX',9X,'YY',9X,'ZZ',9X,'XY',9X,'XZ',9X,'YZ'
     1//1X,'ELECTRONIC',3X,6(F11.6)
     2/1X,'NUCLEAR',6X,6(F11.6)
     3/1X,'TOTAL',8X,6(F11.6))
      AVM=SECM(1)+SECM(2)+SECM(3)
      AVN=SECN(1)+SECN(2)+SECN(3)
      DO 130 I=1,3
      SECM(I)=THRHAL*SECM(I)-HALF*AVM
130   SECN(I)=THRHAL*SECN(I)-HALF*AVN
      DO 140 I=4,6
      SECM(I)=THRHAL*SECM(I)
140   SECN(I)=THRHAL*SECN(I)
      DO 150 I=1,6
150   SECT(I)=SECN(I)+SECM(I)
      WRITE(10,160)
      WRITE(IW,160)
160   FORMAT(//20X,'QUADRUPOLE MOMENTS -- IN ATOMIC UNITS')
      WRITE(10,120) SECM,SECN,SECT
      WRITE(IW,120) SECM,SECN,SECT
      DO 170 I=1,6
      SECN(I)=SECN(I)*FAC
      SECT(I)=SECT(I)*FAC
170   SECM(I)=SECM(I)*FAC
      WRITE(10,180)
      WRITE(IW,180)
180   FORMAT(/20X,'QUADRUPOLE MOMENTS --  IN E.S.U. UNITS'
     1,' (10**-26 ESU CM**2)')
      WRITE(10,120) SECM,SECN,SECT
      WRITE(IW,120) SECM,SECN,SECT
      RETURN
      END
C--------------------------------------------------------------------
C
C
      SUBROUTINE XYZDIP(XS,YS,ZS)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL IANDJ,NORM,DOUBLE
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NNP, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),KR(NATM1),
     5               KA(NATM1),KB(NATM1),
     4               PARVAL(NATM13),ZZZIN(NATM1),X(3,NATM1)
      CHARACTER*8  ANAM,BFLAB,PARNAM,ATNAME
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),
     1               ATNAME(NATM1)
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),KTYPE(NSH1),
     3             KNG(NSH1),KLOC(NSH1)
     2            ,KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2
      LOGICAL  LADAPT,LSKIP
      COMMON/INTDIP/XINT0,YINT0,ZINT0,XINTX,YINTY,ZINTZ,T,X0,Y0,Z0,
     1 XI,YI,ZI,XJ,YJ,ZJ,NI,NJ,ORIGX,ORIGY,ORIGZ
      COMMON/BIGG/DIJ(100),
     &            XIN(32),YIN(32),ZIN(32),
     &            SX(100),SY(100),SZ(100),
     &            IJX(100),IJY(100),IJZ(100),
     &            QXTRA(MAXDIM+NBF1M-646)
      LOGICAL LMOS,LNUCL
      COMMON/PRPSEC/ISECD1,ITYPD1,ISECD2,ITYPD2,LMOS,ISECV,ITYPV
     1     ,ISECD3,ITYPD3,LNUCL,GX,GY,GZ
      DIMENSION XS(*),YS(*),ZS(*)
      DATA ZERO,ONE/0.0D0,1.0D0/
      ORIGX=GX
      ORIGY=GY
      ORIGZ=GZ
      TOL=2.30258D0*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C     ----- ISHELL
      DO 60 II=1,NSHELL
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
        DO 60 JJ=1,II
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
          RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
          IANDJ=II.EQ.JJ
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
          CALL INDXA(IJX,IJY,IJZ,IJ,MINI,MAXI,MINJ,MAXJ,IANDJ,4,1,1)
          DO 10 I=1,IJ
            SX(I)=ZERO
            SY(I)=ZERO
10        SZ(I)=ZERO
C     ----- I PRIMITIVE
          JGMAX=J2
          DO 40 IG=I1,I2
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
            DO 40 JG=J1,JGMAX
              AJ=EX(JG)
              AA=AI+AJ
              AAINV=ONE/AA
              DUM=AJ*ARRI*AAINV
              IF(DUM.GT.TOL) GO TO 40
              FAC=DEXP(-DUM)
              CSJ=CS(JG)*FAC
              CPJ=CP(JG)*FAC
              CDJ=CD(JG)*FAC
              AX=(AXI+AJ*XJ)*AAINV
              AY=(AYI+AJ*YJ)*AAINV
              AZ=(AZI+AJ*ZJ)*AAINV
C     ----- DENSITY FACTOR
              DOUBLE=IANDJ.AND.IG.NE.JG
              CALL DENFAN
     &                   (DIJ,CSI,CPI,CDI,CSJ,CPJ,CDJ,
     &                    MINI,MAXI,MINJ,MAXJ,IANDJ,DOUBLE,NORM)
C     ----- DIPOLE MOMENT INTEGRALS -----
              T=DSQRT(AA)
              TINV=ONE/T
              X0=AX
              Y0=AY
              Z0=AZ
              IN=-4
              DO 20 I=1,LIT
                IN=IN+4
                NI=I
                DO 20 J=1,LJT
                  JN=IN+J
                  NJ=J
                  CALL DIPINT
                  XIN(JN)=XINT0*TINV
                  YIN(JN)=YINT0*TINV
                  ZIN(JN)=ZINT0*TINV
                  XIN(JN+16)=XINTX*TINV
                  YIN(JN+16)=YINTY*TINV
                  ZIN(JN+16)=ZINTZ*TINV
20              CONTINUE
                DO 30 I=1,IJ
                  NX=IJX(I)
                  NY=IJY(I)
                  NZ=IJZ(I)
                  SX(I)=SX(I)+DIJ(I)*(XIN(NX+16)*YIN(NY   )*ZIN(NZ   ))
                  SY(I)=SY(I)+DIJ(I)*(XIN(NX   )*YIN(NY+16)*ZIN(NZ   ))
30                SZ(I)=SZ(I)+DIJ(I)*(XIN(NX   )*YIN(NY   )*ZIN(NZ+16))
40        CONTINUE
C     ----- SET UP DIPOLE MOMENT MATRICES -----
          MAX=MAXJ
          NN=0
          DO 50 I=MINI,MAXI
            LI=LOCI+I
            IN=(LI*(LI-1))/2
            IF(IANDJ) MAX=I
            DO 50 J=MINJ,MAX
              LJ=LOCJ+J
              JN=LJ+IN
              NN=NN+1
              XS(JN)=SX(NN)
              YS(JN)=SY(NN)
50          ZS(JN)=SZ(NN)
60    CONTINUE
      IF (NPRINT.NE.3) GOTO 105
      WRITE(IW,70)
      CALL FOUT(XS,NUM,IW)
      WRITE(IW,80)
      CALL FOUT(YS,NUM,IW)
      WRITE(IW,90)
      CALL FOUT(ZS,NUM,IW)
70    FORMAT(/,10X,25(1H-),/,10X,'X-DIPOLE MOMENT INTEGRALS',/,
     1 10X,25(1H-))
80    FORMAT(/,10X,25(1H-),/,10X,'Y-DIPOLE MOMENT INTEGRALS',/,
     1 10X,25(1H-))
90    FORMAT(/,10X,25(1H-),/,10X,'Z-DIPOLE MOMENT INTEGRALS',/,
     1 10X,25(1H-))
105   CONTINUE
100   RETURN
      END

      SUBROUTINE QUADRU(XXS,YYS,ZZS,XYS,XZS,YZS)
C--------------------------------------------------------------
C     QUADRUPOLE MOMENT INTEGRALS (A.O. BASIS)    R D AMOS
C-------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL IANDJ,NORM,DOUBLE
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NNP, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),KR(NATM1),
     5               KA(NATM1),KB(NATM1),
     4               PARVAL(NATM13),ZZZIN(NATM1),X(3,NATM1)
      CHARACTER*8  ANAM,BFLAB,PARNAM,ATNAME
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),
     1               ATNAME(NATM1)
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),KTYPE(NSH1),
     3             KNG(NSH1),KLOC(NSH1)
     2            ,KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2
      LOGICAL  LADAPT,LSKIP
      COMMON/INTDIP/XINT0,YINT0,ZINT0,XINTX,YINTY,ZINTZ,T,X0,Y0,Z0,
     1 XI,YI,ZI,XJ,YJ,ZJ,NI,NJ,ORIGX,ORIGY,ORIGZ
      COMMON/BIGG/DIJ(100),
     &            SXX(100),SYY(100),SZZ(100),
     &            SXY(100),SXZ(100),SYZ(100),
     &            XIN(32),YIN(32),ZIN(32),
     &            XXIN(16),YYIN(16),ZZIN(16),
     &            IJX(100),IJY(100),IJZ(100),
     &            QJUNK(MAXDIM+NBF1M-994)
      LOGICAL LMOS,LNUCL
      COMMON/PRPSEC/ISECD1,ITYPD1,ISECD2,ITYPD2,LMOS,ISECV,ITYPV
     1      ,ISECD3,ITYPD3,LNUCL,GX,GY,GZ
      DIMENSION XXS(*),YYS(*),ZZS(*),XYS(*),XZS(*),YZS(*)
      DATA ZERO,ONE/0.0D0,1.0D0/
      ORIGX=GX
      ORIGY=GY
      ORIGZ=GZ
      TOL=2.30258D0*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C     ----- ISHELL
      DO 60 II=1,NSHELL
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
      DO 60 JJ=1,II
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
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      IANDJ=II.EQ.JJ
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
      CALL INDXA(IJX,IJY,IJZ,IJ,MINI,MAXI,MINJ,MAXJ,IANDJ,4,1,1)
      DO 10 I=1,IJ
      SXX(I)=ZERO
      SYY(I)=ZERO
      SXY(I)=ZERO
      SXZ(I)=ZERO
      SYZ(I)=ZERO
10    SZZ(I)=ZERO
C     ----- I PRIMITIVE
      JGMAX=J2
      DO 40 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
C
C
C     ----- J PRIMTIVE
C
      IF(IANDJ) JGMAX=IG
      DO 40 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AAINV=ONE/AA
      DUM=AJ*ARRI*AAINV
      IF(DUM.GT.TOL) GO TO 40
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
     1   MINJ,MAXJ,IANDJ,DOUBLE,NORM)
C     ----- SECOND MOMENT INTEGRALS -----
      T=DSQRT(AA)
      TINV=ONE/T
      X0=AX
      Y0=AY
      Z0=AZ
      IN=-4
      DO 20 I=1,LIT
      IN=IN+4
      NI=I
      DO 20 J=1,LJT
      JN=IN+J
      NJ=J
      CALL DIPINT
      XIN(JN)=XINT0*TINV
      YIN(JN)=YINT0*TINV
      ZIN(JN)=ZINT0*TINV
      XIN(JN+16)=XINTX*TINV
      YIN(JN+16)=YINTY*TINV
      ZIN(JN+16)=ZINTZ*TINV
      CALL QUDINT
      XXIN(JN)=XINTX*TINV
      YYIN(JN)=YINTY*TINV
20    ZZIN(JN)=ZINTZ*TINV
      DO 30 I=1,IJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DD=DIJ(I)
      SXY(I)=SXY(I)+DD*XIN(NX+16)*YIN(NY+16)*ZIN(NZ)
      SXZ(I)=SXZ(I)+DD*XIN(NX+16)*YIN(NY)*ZIN(NZ+16)
      SYZ(I)=SYZ(I)+DD*XIN(NX)*YIN(NY+16)*ZIN(NZ+16)
      SXX(I)=SXX(I)+DD*XXIN(NX)*YIN(NY)*ZIN(NZ)
      SYY(I)=SYY(I)+DD*XIN(NX)*YYIN(NY)*ZIN(NZ)
30    SZZ(I)=SZZ(I)+DD*XIN(NX)*YIN(NY)*ZZIN(NZ)
40    CONTINUE
C     ----- SET UP SECOND MOMENT MATRICES -----
      MAX=MAXJ
      NN=0
      DO 50 I=MINI,MAXI
      LI=LOCI+I
      IN=(LI*(LI-1))/2
      IF(IANDJ) MAX=I
      DO 50 J=MINJ,MAX
      LJ=LOCJ+J
      JN=LJ+IN
      NN=NN+1
      XXS(JN)=SXX(NN)
      YYS(JN)=SYY(NN)
      ZZS(JN)=SZZ(NN)
      XYS(JN)=SXY(NN)
      XZS(JN)=SXZ(NN)
50    YZS(JN)=SYZ(NN)
60    CONTINUE
      IF (NPRINT.NE.3) GOTO 125
      WRITE(IW,70)
      CALL FOUT(XXS,NUM,IW)
      WRITE(IW,80)
      CALL FOUT(YYS,NUM,IW)
      WRITE(IW,90)
      CALL FOUT(ZZS,NUM,IW)
70    FORMAT(/,10X,25(1H-),/,10X,'XX-SECOND MOMENT INTEGRALS',/,
     1 10X,25(1H-))
80    FORMAT(/,10X,25(1H-),/,10X,'YY-SECOND MOMENT INTEGRALS',/,
     1 10X,25(1H-))
90    FORMAT(/,10X,25(1H-),/,10X,'ZZ-SECOND MOMENT INTEGRALS',/,
     1 10X,25(1H-))
      WRITE(IW,100)
      CALL FOUT(XYS,NUM,IW)
100   FORMAT(/' XY-INTEGRALS'/)
      WRITE(IW,110)
      CALL FOUT(XZS,NUM,IW)
110   FORMAT(/' XZ-INTEGRALS'/)
      WRITE(IW,120)
      CALL FOUT(YZS,NUM,IW)
120   FORMAT(/' YZ-INTEGRALS'/)
125   CONTINUE
130   RETURN
      END
C-------------------------------------------------------
      SUBROUTINE DIPINT
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      COMMON/INTDIP/XINT0,YINT0,ZINT0,XINTX,YINTY,ZINTZ,T,X0,Y0,Z0,
     1 XI,YI,ZI,XJ,YJ,ZJ,NI,NJ,GX,GY,GZ
      COMMON/HERMIT/H1,H2(2),H3(3),H4(4),H5(5),H6(6),H7(7)
      COMMON/WERMIT/W1,W2(2),W3(3),W4(4),W5(5),W6(6),W7(7)
      DIMENSION H(28),W(28),MIN(6),MAX(6)
      EQUIVALENCE (H(1),H1),(W(1),W1)
      DATA MIN /1,2,4,7,11,16/
      DATA MAX /1,3,6,10,15,21/
      DATA ZERO /0.0D0/
      XINT0=ZERO
      YINT0=ZERO
      ZINT0=ZERO
      XINTX=ZERO
      YINTY=ZERO
      ZINTZ=ZERO
      NPTS=(NI+NJ-1)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 130 I=IMIN,IMAX
      DUM=W(I)
      PX=DUM
      PY=DUM
      PZ=DUM
      DUM=H(I)/T
      PTX=DUM+X0
      PTY=DUM+Y0
      PTZ=DUM+Z0
      AX=PTX-XI
      AY=PTY-YI
      AZ=PTZ-ZI
      BX=PTX-XJ
      BY=PTY-YJ
      BZ=PTZ-ZJ
      GO TO (50,40,30,20,10),NI
10    PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
20    PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
30    PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
40    PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
50    GO TO (120,110,100,90,80,70,60),NJ
60    PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
70    PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
80    PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
90    PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
100   PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
110   PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
120   CONTINUE
      XINT0=XINT0+PX
      YINT0=YINT0+PY
      ZINT0=ZINT0+PZ
      XINTX=XINTX+PX*(PTX-GX)
      YINTY=YINTY+PY*(PTY-GY)
      ZINTZ=ZINTZ+PZ*(PTZ-GZ)
130   CONTINUE
      RETURN
      END
      SUBROUTINE QUDINT
C------------------------------------------------------------
C     QUADRUPOLE MOMENT SUBSIDIARY INTEGRALS - R D AMOS
C-------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      COMMON/INTDIP/XINT0,YINT0,ZINT0,XINTX,YINTY,ZINTZ,T,X0,Y0,Z0,
     1 XI,YI,ZI,XJ,YJ,ZJ,NI,NJ,GX,GY,GZ
      COMMON/HERMIT/H(28)
      COMMON/WERMIT/W(28)
      DIMENSION MIN(6),MAX(6)
      DATA MIN /1,2,4,7,11,16/
      DATA MAX /1,3,6,10,15,21/
      DATA ZERO /0.0D0/
      XINTX=ZERO
      YINTY=ZERO
      ZINTZ=ZERO
      NPTS=(NI+NJ)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 130 I=IMIN,IMAX
      DUM=W(I)
      PX=DUM
      PY=DUM
      PZ=DUM
      DUM=H(I)/T
      PTX=DUM+X0
      PTY=DUM+Y0
      PTZ=DUM+Z0
      AX=PTX-XI
      AY=PTY-YI
      AZ=PTZ-ZI
      BX=PTX-XJ
      BY=PTY-YJ
      BZ=PTZ-ZJ
      GO TO (50,40,30,20,10),NI
10    PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
20    PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
30    PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
40    PX=PX*AX
      PY=PY*AY
      PZ=PZ*AZ
50    GO TO (120,110,100,90,80,70,60),NJ
60    PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
70    PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
80    PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
90    PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
100   PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
110   PX=PX*BX
      PY=PY*BY
      PZ=PZ*BZ
120   CONTINUE
      XINTX=XINTX+PX*(PTX-GX)*(PTX-GX)
      YINTY=YINTY+PY*(PTY-GY)*(PTY-GY)
      ZINTZ=ZINTZ+PZ*(PTZ-GZ)*(PTZ-GZ)
130   CONTINUE
      RETURN
      END

      SUBROUTINE ELECTR(WRTOUT)
C
C    ELCTRIC FIELDS AT NUCLEI
C    R D AMOS
C
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL IANDJ,WRTOUT,NORM,DOUBLE
      COMMON/RTWT/XX,U(NATM1),W(NATM1),NROOTS
C
      LOGICAL  LADAPT,LSKIP
      CHARACTER*8 DNAM(3)
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NNP, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),KR(NATM1),
     5               KA(NATM1),KB(NATM1),
     4               PARVAL(NATM13),ZZZIN(NATM1),X(3,NATM1)
      CHARACTER*8  ANAM,BFLAB,PARNAM,ATNAME
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),
     1               ATNAME(NATM1)
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),KTYPE(NSH1),
     3             KNG(NSH1),KLOC(NSH1)
     2            ,KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      COMMON/BIG/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     &           NI,NJ,
     &           CX,CY,CZ,
     &           IJN(100),IJX(36),IJY(36),IJZ(36),
     &           XIN(150),YIN(150),ZIN(150),DFAC(100),DIJ(100),
     &           HELF(36,NATM1,3),DD(NBF1M),
     &           QXTRA(NBF1M-771-108*NATM1+2*NBF1)
      COMMON/BIGG/DDV(MAXDIM),DDE(NBF1M)

      COMMON/MAPPER/IA(NBF1+1)
      COMMON/BIGK/DE(3,NATM1),DN(3,NATM1),DJNK(NBF1M-6*NATM1)
C
      COMMON/BLKIN/GOUT(510),NINTH
C
      DATA NDIM/5/
      DATA DNAM/'Ex','Ey','Ez'/
      DATA PI212 /1.1283791670955D0/
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/
C
      WRITE(10,111) (DNAM(I),I=1,3)
      WRITE(IW,111) (DNAM(I),I=1,3)
111   FORMAT(//10X,'ELECTRIC FIELD AT NUCLEI'//14X,'ATOM',
     &  8X,A2,12X,A2,12X,A2/)
      I1=1+NNP
      CALL ONEPDM(DD(1),DDV(1))
C
      TOL=2.30258D0*ITOL
      DO 20 I=1,3
        DO 20 N=1,NAT
20        DE(I,N)=ZERO
      IPOSF1=IBLKS
      NINTH=1
      IF(WRTOUT) CALL SEARCH(IPOSF1,IFILS)
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
60        HELF(I,J,K)=ZERO
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
                CX=C(1,IC)
                CY=C(2,IC)
                CZ=C(3,IC)
                XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
                CALL ROOTS
                MM=0
                DO 210 K=1,NROOTS
                  UU=AA*U(K)
                  WW=W(K)
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
     &                   XIN(NX+MM+ND2R)*YIN(NY+MM   )*ZIN(NZ+MM   )
                    DUMY=DUMY+
     &                   XIN(NX+MM   )*YIN(NY+MM+ND2R)*ZIN(NZ+MM   )
                    DUMZ=DUMZ+
     &                   XIN(NX+MM   )*YIN(NY+MM   )*ZIN(NZ+MM+ND2R)
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
230             DE(3,IC)=DE(3,IC)+DUM*DUMZ
240           CONTINUE
250       CONTINUE

280       IF (.NOT.WRTOUT) GOTO 310
          DO 300 I=1,IJ
            DO 300 J=1,NAT
              DO 290 K=1,3
                GOUT(NINTH)=HELF(I,J,K)
290             NINTH=NINTH+1
              IF(NINTH.LT.511) GOTO 300
              NINTH=NINTH-1
              CALL PUT(GOUT,511,IFILS)
              NINTH=1
              IPOSF1=IPOSF1+1
300       CONTINUE

C     ----- END OF 'PRIMITIVE' LOOPS -----
310   CONTINUE
C     ----- END OF 'SHELL' LOOPS -----

      IF(.NOT.WRTOUT) GOTO 320
      NINTH=NINTH-1
      CALL PUT(GOUT,511,IFILS)
      IPOSF1=IPOSF1+1
320   CALL NUCFIE(DE,NAT,ZAN,C)
      DO 330 I=1,NAT
      WRITE(IW,340) i,ANAM(I),(DE(J,I),J=1,3)
      WRITE(10,340) i,ANAM(I),(DE(J,I),J=1,3)
340   FORMAT(6X,i5,1x,A8,3(F14.7))
330   CONTINUE
      RETURN
      END

      SUBROUTINE NUCFIE(DE,NAT,ZAN,C)
C
C    NUCLEAR COMPONENTS OF FIELDS AT NUCLEI
C
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      DIMENSION  ZAN(NATM1),C(3,NATM1),DE(3,NATM1)
      COMMON/BIG/DRG(NATM1,NATM1),DXTRA(2*NBF1M+2*NBF1-NATM1*NATM1)
      DATA ZERO,ONE/0.0D0,1.0D0/
      DRG(1,1)=ZERO
      DO 30 K=2,NAT
      DRG(K,K)=ZERO
      K1=K-1
      DO 30 L=1,K1
      RKL=ZERO
      DO 20 I=1,3
20    RKL=RKL+(C(I,K)-C(I,L))**2
      DRG(K,L)=ONE/RKL
30    DRG(L,K)=DSQRT(RKL)
      DO 60 KK=1,3
      DO 40 K=2,NAT
      KM1=K-1
      DO 40 L=1,KM1
      ZAL=ZAN(L)
      PKL=(C(KK,K)-C(KK,L))/DRG(L,K)
      DE(KK,K)=DE(KK,K)+PKL*DRG(K,L)*ZAL
40    CONTINUE
      NAT1=NAT-1
      DO 50 K=1,NAT1
      KP1=K+1
      DO 50 L=KP1,NAT
      ZAL=ZAN(L)
      PKL=(C(KK,K)-C(KK,L))/DRG(K,L)
      DE(KK,K)=DE(KK,K)+PKL*DRG(L,K)*ZAL
50    CONTINUE
60    CONTINUE
      RETURN
      END

