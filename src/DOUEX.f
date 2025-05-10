      SUBROUTINE DOUEX(GTHR,NRES)
      IMPLICIT REAL*8 (A-H,O-Z)
C----------------------------
      PARAMETER (nbf1 = 50)
      PARAMETER (NBFT = nbf1*(nbf1+1)/2)
      PARAMETER (NBUF = 50000)
C------------------------------------------
      COMMON/OMEGA/INDEX(NBUF)
      COMMON/GAMM/FELD(NBUF)
      parameter (
      DIMENSION TONEN(NBUF)
      DIMENSION LABOC(nbf1),LABVA(nbf1),LIN(NBUF)
      DIMENSION C(nbf1,nbf1),VAC(nbf1,nbf1),ZAEH(NBUF)
      DIMENSION ONE(NBFT),DENS(NBFT),ONEII(NBFT),ONEAA(NBFT)

      NMAX=NBUF
      INP=5
      IOUT=6
      WRITE(IOUT,866) NRES,GTHR
  866 FORMAT(//10X,'RESTART CONDITION:',I10,
     *       //10X,'THRESHOLD FOR NOMINATOR:',F10.7)
      DO 733 I=1,NBUF
        TONEN(I)=0.0
  733   ZAEH(I)=0.
C-----------------------------------------------------------------------
C          READ NA,NC,NLO,NVAC FROM TAPE 10 GENERAL INFORMATION
C------------------------------------------------------------------------  
      REWIND 10
      READ(10) NA,NC,NLO,NVAC
C----------------------------------------------------------
C          READ FROM TAPE 10 ONE-ELECTRON INTEGRALS
C-------------------------------------------------------
      READ(10) NAME,LLL,(ONE(I),I=1,LLL)
      READ(10) (ONE(I),I=1,LLL)
      READ(10) (ONE(I),I=1,LLL)
   10 FORMAT(/10X,' TEST READ IN:',5X,I4,5X,I8,5X,I4)
  197 FORMAT(//10X,'LABEL FOR OCCUPIED ORBITALS',/,10X,30I3)
  198 FORMAT(//10X,'LABELS FOR VACANT ORBITALS',/,10X,30I3)
C-----------------------------------------------------------
C          READ WAVE FUNCTION FROM FILE 1I.
C--------------------------------------------------------------------
      REWIND 12
      READ(12) NAME,NGG,((C(I,J),I=1,NC),J=1,NLO)
      WRITE(IOUT,664)
  664 FORMAT(/10X,'OCCUPIED ORBITALS',//)
      DO 787 J=1,NLO
        WRITE(IOUT,665) (C(I,J),I=1,NC)
  665   FORMAT(10(F8.5,3X))
        WRITE(IOUT,789)
  789   FORMAT(//)
  787 CONTINUE
C-----------------------------------------------------------------
C          CALCULATE DENSITY MATRIX
C-----------------------------------------------------------------
      NZ=0
      DO 17 I=1,NC
        DO 17 J=1,I
          NZ=NZ+1
          DENS(NZ)=0.
          DO 18 K=1,NLO
   18       DENS(NZ)=DENS(NZ)+C(I,K)*C(J,K)
          DENS(NZ)=DENS(NZ)+DENS(NZ)
   17 CONTINUE

      READ(12) NAME,NGG,((VAC(I,J),I=1,NC),J=1,NVAC)
      WRITE(IOUT,663)
  663 FORMAT(//10X,'VACANT ORBITALS',//)
      DO 796 J=1,NVAC
        WRITE(IOUT,665) (VAC(I,J),I=1,NC)
        WRITE(IOUT,789)
  796 CONTINUE
C-------------------------------------------------------------------
C.......C........PRESET LABOC AND LABVA WITH VALUE 99.
C-----------------------------------------------------------
      DO 601 I=1,NLO
  601   LABOC(I)=99
      DO 602 I=1,NVAC
  602   LABVA(I)=99
C---------------------------------------------------------
C.......C........AUTOMATIC DETERMINATION OF ORBITAL TYPE
C----------------------------------------------------------------
      III=0
      LABOC(1)=0
      DO 603 I=1,NLO
        IF(LABOC(I).LT.99) GOTO 603
        DO 604 J=1,NLO
          S=0.
          IF(LABOC(J).LT.99) GOTO 604
          DO 605 K=1,NC
  605       S=S+ABS(C(K,I))*ABS(C(K,J))
          IF(ABS(S).LT.1.E-8) GOTO 604
          LABOC(J)=III
  604   CONTINUE
        III=III+1
  603 CONTINUE
      III=III-1
C-----------------------------------------------------------------
C.......C........NOW COMPARE VACANT ORBITALS WITH OCCUPIED ONES.
C-------------------------------------------------------------------
      DO 606 I=1,NVAC
        DO 607 J=1,NLO
          S=0.
          DO 608 K=1,NC
  608       S=S+ABS(VAC(K,I))*ABS(C(K,J))
          IF(ABS(S).LT.1.E-8) GOTO 607
          LABVA(I)=LABOC(J)
          GOTO 606
  607   CONTINUE
  606 CONTINUE
C----------------------------------------------------------------
C.......C........CHEK IF ALL TYPES OF VACANT ORBITALS ARE DETERMINED
C---------------------------------------------------------------------
      DO 609 I=1,NVAC
        IF(LABVA(I).LT.99) GOTO 609
        III=III+1
        DO 610 J=1,NVAC
          S=0.
          IF(LABVA(J).LT.99) GOTO 610
          DO 611 K=1,NC
  611       S=S+ABS(VAC(K,I))*ABS(VAC(K,J))
          IF(ABS(S).LT.1.E-8) GOTO 610
          LABVA(J)=III
  610   CONTINUE
  609 CONTINUE
      WRITE(IOUT,197) (LABOC(I),I=1,NLO)
      WRITE(IOUT,198) (LABVA(I),I=1,NVAC)

C------------------------------------------------------
C          END OF SYMMETRY DETERMINATION OF ORBITALS
C--------------------------------------------------------
C          HERE STARTS THE RESTART IN CASE  NRES = 1
C--------------------------------------------------------

      IF(NRES.EQ.99) GOTO 7000
C----------------------------------------------------------------------

      NZ=1
      NVG=0
      REWIND 14
      READ(14) NAME,NN,(FELD(I),I=1,NN)
      READ(14) (INDEX(I),I=1,NN)
      WRITE(IOUT,864) NN
  864 FORMAT(//10X,'NUMBER OF INTEGRAL PACKAGE',I10)
      DO 1010 IC=1,NC
        DO 1020 JC=1,IC
          DO 1030 KC=1,IC
            LAST=KC
            IF(KC.EQ.IC) LAST=JC
            DO 1000 LC=1,LAST
              NVG=NVG+1
              IF(NVG.LE.NN) GOTO 50
              READ(14) NAME,NN,(FELD(I),I=1,NN)
              READ(14) (INDEX(I),I=1,NN)
              NVG=1
50            CONTINUE
              IF(NZ-INDEX(NVG)) 1011,1012,1011

1011          NVG=NVG-1
              NZ=NZ+1
              GOTO 1000

1012          CAA=FELD(NVG)
              NZ=NZ+1
C------------------------------------------------------------------
C.......C........ZAEHLER (I(1)J(2)/K(1)L(2)) - (I(1)J(2)/L(1)K(2))
C.......C........IF SPIN (1) = SPIN(2)
C.......C........OTHERWISE WITHOUT EXCHANGE INTEGRAL
C------------------------------------------------------------------
              NQ=0
              DO 200 I=1,NLO
                CII=C(IC,I)
                CJI=C(JC,I)
                CKI=C(KC,I)
                CLI=C(LC,I)
                DO 207 J=1,I
                  MIJ=LABOC(I)-LABOC(J)
                  CIJ=C(IC,J)
                  CJJ=C(JC,J)
                  CKJ=C(KC,J)
                  CLJ=C(LC,J)
                  CIKJ=CII*CKJ
                  CILJ=CII*CLJ
                  CJKJ=CJI*CKJ
                  CJLJ=CJI*CLJ
                  CKIJ=CKI*CIJ
                  CKJJ=CKI*CJJ
                  CLIJ=CLI*CIJ
                  CLJJ=CLI*CJJ
                  DO 208 K=1,NVAC
                    MIK=LABOC(I)-LABVA(K)
                    MJK=LABOC(J)-LABVA(K)
                    CIK=VAC(IC,K)
                    CJK=VAC(JC,K)
                    CKK=VAC(KC,K)
                    CLK=VAC(LC,K)
                    LEND=NVAC
                    IF(I.EQ.J) LEND=K
                    DO 208 L=1,LEND
                      MIL=LABOC(I)-LABVA(L)
                      MJL=LABOC(J)-LABVA(L)
                      MKL=LABVA(K)-LABVA(L)
                      IF(MIJ.EQ.0.AND.MKL.EQ.0) GOTO 222
                      IF(MIL.EQ.0.AND.MJK.EQ.0) GOTO 222
                      IF(MIK.EQ.0.AND.MJL.EQ.0) GOTO 222
                      GOTO 208
222                   CONTINUE
                      CIL=VAC(IC,L)
                      CJL=VAC(JC,L)
                      CKL=VAC(KC,L)
                      CLL=VAC(LC,L)
                      CJLL=CJK*CLL
                      CJKL=CJK*CKL
                      CILL=CIK*CLL
                      CIKL=CIK*CKL
                      CLJL=CLK*CJL
                      CLIL=CLK*CIL
                      CKJL=CKK*CJL
                      CKIL=CKK*CIL
                      NQ=NQ+1
                      CC1=CAA*(CIKJ*CJLL+CILJ*CJKL+CJKJ*CILL+CJLJ*CIKL
     &                        +CKIJ*CLJL+CKJJ*CLIL+CLIJ*CKJL+CLJJ*CKIL)
                      ZAEH(NQ)=ZAEH(NQ)+CC1
                      IF(I.NE.J.AND.K.NE.L) GOTO 512
                      GOTO 510
512                   CONTINUE
                      NQ=NQ+1
                      ZAEH(NQ)=ZAEH(NQ)-CAA*
     &                         (CIKJ*CLJL+CILJ*CKJL+CJKJ*CLIL+CJLJ*CKIL
     &                         +CKIJ*CJLL+CKJJ*CILL+CLIJ*CJKL+CLJJ*CIKL)
     &                        +CC1
510                   CONTINUE
208               CONTINUE
207             CONTINUE
200           CONTINUE
 1000       CONTINUE
 1030     CONTINUE
 1020   CONTINUE
1010  CONTINUE
C--------------------------------------------------------------
C.......C..SELECT OUT CONFIGURATIONS ACCORDING TO NOMINATOR....
C------------------------------------------------------------------
      NZW=0
      NZER=0
      DO 819 I=1,NQ
        ZAE=ZAEH(I)
        IF(ABS(ZAE).LT.GTHR) GOTO 837
        NZW=NZW+1
        ZAEH(NZW)=ZAE
        GOTO 819
  837   NZER=NZER+1
        LIN(NZER)=I
  819 CONTINUE
      NQ=NZW
      WRITE(IOUT,838) NZER
  838 FORMAT(//5X,'ZERO MATRIX ELEMENTS(NUMBER)',10X,I8)
      NAMO=3200
      WRITE(13) NAMO,NQ,(ZAEH(I),I=1,NQ)
      WRITE(IOUT,10)NAMO,NQ
      NAMO=3300
      WRITE(13)NAMO,NZER,(LIN(I),I=1,NZER)
      WRITE(IOUT,10) NAMO,NZER
 7000 CONTINUE
      IF(NRES.NE.99)GOTO 7200
      REWIND 13
      READ(13) NAMO,NQ,(ZAEH(I),I=1,NQ)
      WRITE(IOUT,10)NAMO,NQ
      READ(13)NAMO,NZER,(LIN(I),I=1,NZER)
      WRITE(IOUT,10)NAMO,NZER
 7200 CONTINUE
C----------------------------------------------------------------------
C     CALCULATE SUM OVER ONE-ELECTRON INTEGRALS FOR EVERY CONFIGURATION
C     CONSTRUCTION OF SPINORBITALS NO LONGER NECESSARY
C     CALCULATE ONE-ELECTRON INTEGRAL
C     FIRST <I/I>
C-------------------------------------------------------------------------
      NZ=0
      DO 11 I=1,NLO
        NZ=NZ+1
        ONEII(NZ)=0.
        DO 11 K=1,NC
          DO 11 L=1,K
            KL=K*(K-1)/2+L
            FAK=C(K,I)*C(L,I)
            IF(L.LT.K)FAK=FAK+FAK
   11       ONEII(NZ)=ONEII(NZ)+FAK*ONE(KL)
      WRITE(IOUT,13)
   13 FORMAT(/5X,'ONE-ELECTRON TERMS: <I/I>',//)
      WRITE(IOUT,15)(ONEII(1),I=1,NZ)
   15 FORMAT(5X,10F10.6)
C---------------------------------------------
C     NOW <A/A>
C---------------------------------------------
      NZ=0
      DO 12 I=1,NVAC
        NZ=NZ+1
        ONEAA(NZ)=0.
        DO 12 K=1,NC
          DO 12 L=1,K
            KL=K*(K-1)/2+L
            FAK=VAC(K,I)*VAC(L,I)
            IF(L.LT.K)FAK=FAK+FAK
   12 ONEAA(NZ)=ONEAA(NZ)+FAK*ONE(KL)
      WRITE(IOUT,14)
   14 FORMAT(/5X,'ONE-ELECTRON TERMS: <A/A>',//)
      WRITE(IOUT,15)(ONEAA(I),I=1,NZ)
C--------------------------------------------------
C     RESTART POINT
C----------------------------------------------------------------------
      IF(NRES.EQ.88)GOTO 9000
C----------------------------------------------------
C     READ INTEGRALS FROM TAPE 14
C-----------------------------------------------------------
      NVG=0
      NZ=1
      NMAX=1540
      REWIND 14
      READ(14) NAME,NN,(FELD(I),I=1,NN)
      READ(14) (INDEX(I),I=1,NN)
      DO 4010 IC=1,NC
        DO 4020 JC=1,IC
          DO 4030 KC=1,IC
            LAST=KC
            IF(KC.EQ.IC) LAST=JC
            DO 4000 LC=1,LAST
              NVG=NVG+1
              IF(NVG.LE.NN) GO TO 55
              READ(14) NAME,NN,(FELD(I),I=1,NN)
              READ(14) (INDEX(I),I=1,NN)
              NVG=1
   55         CONTINUE
              IF(NZ-INDEX(NVG)) 5011,5012,5011

 5011         NVG=NVG-1
              NZ=NZ+1
              GO TO 4000

 5012         CAA=FELD(NVG)
              NZ=NZ+1
              IJ=IC*(IC-1)/2+JC
              IK=IC*(IC-1)/2+KC
              IL=IC*(IC-1)/2+LC
              JK=JC*(JC-1)/2+KC
              IF(JC.LT.KC) JK=KC*(KC-1)/2+JC
              JL=JC*(JC-1)/2+LC
              IF(JC.LT.LC) JL=LC*(LC-1)/2+JC
              KL=KC*(KC-1)/2+LC
              PIJ=DENS(IJ)
              PIK=DENS(IK)
              PIL=DENS(IL)
              PJK=DENS(JK)
              PJL=DENS(JL)
              PKL=DENS(KL)
              NZA=0
              NAB=1
              NQ=0
              DO 250 I=1,NLO
                AI=C(IC,I)
                AJ=C(JC,I)
                AK=C(KC,I)
                AL=C(LC,I)
                DO 251 J=1,I
                  MIJ=LABOC(I)-LABOC(J)
                  BI=C(IC,J)
                  BJ=C(JC,J)
                  BK=C(KC,J)
                  BL=C(LC,J)
                  DO 252 K=1,NVAC
                    MIK=LABOC(I)-LABVA(K)
                    MJK=LABOC(J)-LABVA(K)
                    CI=VAC(IC,K)
                    CJ=VAC(JC,K)
                    CK=VAC(KC,K)
                    CL=VAC(LC,K)
                    LEND=NVAC
                    IF(I.EQ.J) LEND=K
                    DO 550 L=1,LEND
                      MIL=LABOC(I)-LABVA(L)
                      MJL=LABOC(J)-LABVA(L)
                      MKL=LABVA(K)-LABVA(L)
                      IF(MIJ.EQ.0.AND.MKL.EQ.0) GO TO 333
                      IF(MIK.EQ.0.AND.MJL.EQ.0) GO TO 333
                      IF(MIL.EQ.0.AND.MJK.EQ.0) GO TO 333
                      GO TO 550
  333                 CONTINUE
                      NTRIP=1
                      NSING=1
                      IF(I.EQ.J.OR.K.EQ.L) GO TO 450
                      NZA=NZA+1
                      IF(NZA.NE.LIN(NAB)) GO TO 452
                      NAB=NAB+1
                      NTRIP=0
                      GO TO 453
  452                 NTRIP=1
  453                 NZA=NZA+1
                      IF(NZA.NE.LIN(NAB)) GO TO 454
                      NAB=NAB+1
                      NSING=0
                      IF(NTRIP.EQ.0.AND.NSING.EQ.0) GO TO 550
                      GO TO 451
  454                 NSING=1
                      GO TO 451
  450                 NZA=NZA+1
                      IF(NZA.NE.LIN(NAB)) GO TO 451
                      NAB=NAB+1
                      NTRIP=0
                      GO TO 550
  451                 CONTINUE
                      DI=VAC(IC,L)
                      DJ=VAC(JC,L)
                      DK=VAC(KC,L)
                      DL=VAC(LC,L)
                      IF(I.EQ.J.OR.K.EQ.L) GO TO 777
C--------------------------------------------------------------------
C.......C........ CASE: I.NE.J AND. K.NE.L
C-------------------------------------------------------------------
                      IF(NSING.EQ.0) GO TO 777
C------------------------------------------------------------------
C.......C........ CALCULATE KAB+ KIJ - KAJ - KBI
C------------------------------------------------------------------
      EKAJ=(CI*CL*BJ*BK+CI*CK*BJ*BL+CJ*CL*BI*BK+CJ*CK*BI*BL)
      EKIJ=(AI*AL*BJ*BK+AI*AK*BJ*BL+AJ*AL*BI*BK+AJ*AK*BI*BL)
      EKBI=(DI*DL*AJ*AK+DI*DK*AJ*AL+DJ*DL*AI*AK+DJ*DK*AI*AL)
      EKAB=(CI*CL*DJ*DK+CI*CK*DJ*DL+CJ*CL*DI*DK+CJ*CK*DI*DL)
      AUST=2.0*(EKAB+EKIJ-EKAJ-EKBI)
  777 CONTINUE
      EJJ=(CI*CJ+DI*DJ)*(AK*AL+BK*BL)+(CK*CL+DK*DL)*(AI*AJ+BI*BJ)
      EJAB=(CI*CJ*DK*DL+CK*CL*DI*DJ)
      EJIJ=(AI*AJ*BK*BL+AK*AL*BI*BJ)
      ECC=4.0*(EJJ-EJAB-EJIJ)
      EKAI=(CI*CL*AJ*AK+CI*CK*AJ*AL+CJ*CL*AI*AK+CJ*CK*AI*AL)
      EKBJ=(DI*DL*BJ*BK+DI*DK*BJ*BL+DJ*DL*BI*BK+DJ*DK*BI*BL)
      EKK=2.0*(EKAI+EKBJ)
      EEL=ECC-EKK
      COUL=4.0*(PKL*(AI*AJ+BI*BJ-CI*CJ-DI*DJ)+
     *          PIJ*(AK*AL+BK*BL-CK*CL-DK*DL))
      EEL=EEL+COUL
      EXX=PJK*(CI*CL+DI*DL-AI*AL-BI*BL)+
     *    PJL*(CI*CK+DI*DK-AI*AK-BI*BK)+
     *    PIK*(CJ*CL+DJ*DL-AJ*AL-BJ*BL)+
     *    PIL*(CJ*CK+DJ*DK-AJ*AK-BJ*BK)
      EEL=EEL+EXX
      NQ=NQ+1
      TONEN(NQ)=TONEN(NQ)+EEL*CAA
      IF(I.EQ.J.OR.K.EQ.L) GO TO 550
      IF(NSING.EQ.0) GO TO 550
      IF(NTRIP.EQ.0) GO TO 799
      NQ=NQ+1
      TONEN(NQ)=TONEN(NQ)+(EEL+AUST)*CAA
      GO TO 550
  799 TONEN(NQ)=TONEN(NQ)+AUST*CAA
  550               CONTINUE
  252             CONTINUE
  251           CONTINUE
  250         CONTINUE
 4000       CONTINUE
 4030     CONTINUE
 4020   CONTINUE
 4010 CONTINUE
      WRITE(IOUT,710) NQ
  710 FORMAT(//10X,'NUMBER OF DOUBLE-EXCITED CONFIGURATIONS',I10)
C----------------------------------------------------------------------
C.......C........ CALCULATE ENERGY-CONTRIBUTION FOR EACH CONFIGURATION.
C
C.......C.........CALCULATION OF TOTAL ENERGY CONTRIBUTION.
C----------------------------------------------------------------------
      WRITE(IOUT,716)
  716 FORMAT(//19X,'I',6X,'J',6X,'K',6X,'L',6X,'EONE',10X,'NENNER',8X,
     ;'ZAEHLER',8X,'TOTAL ENERGY',//)
      NZA=0
      NAB=1
      NUM=0
      EDOU=0.
      DO 21 I=1,NLO
        DO 27 J=1,I
          FF=1.0
          MIJ=LABOC(I)-LABOC(J)
          DO 25 K=1,NVAC
            MIK=LABOC(I)-LABVA(K)
            MJK=LABOC(J)-LABVA(K)
            LAST=NVAC
            IF(I.EQ.J) LAST=K
            DO 25 L=1,LAST
              MIL=LABOC(I)-LABVA(L)
              MJL=LABOC(J)-LABVA(L)
              MKL=LABVA(K)-LABVA(L)
              IF(MIJ.EQ.0.AND.MKL.EQ.0) GO TO 444
              IF(MIK.EQ.0.AND.MJL.EQ.0) GO TO 444
              IF(MIL.EQ.0.AND.MJK.EQ.0) GO TO 444
              GO TO 25
  444         CONTINUE
              FAK=FF*2.0
              IF(I.EQ.J.AND.K.EQ.L) FAK=FF
              NZA=NZA+1
              IF(NZA.NE.LIN(NAB)) GO TO 23
              NAB=NAB+1
              GO TO 310
   23         NUM=NUM+1
              EONE     =ONEII(I)+ONEII(J)-ONEAA(K)-ONEAA(L)
              TONEN(NUM)=TONEN(NUM)+EONE
              ZW=TONEN(NUM)
              TONEN(NUM)=ZAEH(NUM)**2/TONEN(NUM)
              EDOU=EDOU+TONEN(NUM)*FAK
              WRITE(IOUT,28) NUM,I,J,K,L,EONE,ZW,ZAEH(NUM),TONEN(NUM)
  310         IF(I.NE.J.AND.K.NE.L) GO TO 31
              GO TO 25
   31         NZA=NZA+1
              IF(NZA.NE.LIN(NAB)) GO TO 24
              NAB=NAB+1
              GO TO 25
   24         NUM=NUM+1
              EONE=ONEII(I)+ONEII(J)-ONEAA(K)-ONEAA(L)
              TONEN(NUM)=TONEN(NUM)+EONE
              ZW=TONEN(NUM)
              TONEN(NUM)=ZAEH(NUM)**2/ZW
              EDOU=EDOU+TONEN(NUM)*FAK
              WRITE(IOUT,28) NUM,I,J,K,L,EONE,ZW,ZAEH(NUM),TONEN(NUM)
   25     CONTINUE
   27   CONTINUE
   21 CONTINUE
      WRITE(IOUT,29) NUM
   29 FORMAT(//10X,'TOTAL NUMBER OF CONFIGURATIONS:',I10)
   28 FORMAT(10X,I4,'.',4X,4(I2,4X),4F15.10)
  712 FORMAT(/,3(2X,F15.8,3X,F15.8,10X))
  714 FORMAT  ( 1X,3(1X,F15.8,3X,F15.8,10X))
C---------------------------------------------------
C          END OF ENERGY CONTRIBUTION CALCULATION
C-----------------------------------------------------
      EEV=EDOU*27.21
      EKC=EDOU*627.71
      WRITE(IOUT,730) EDOU,EEV,EKC
  730 FORMAT(///10X,'TOTAL ENERGY CONTRIBUTION OF DOUBLY EXCITED STATES'
     F,2X,F20.10,5X,'(A.U.)',//,60X,F20.10,5X,'(EV)'//,60X,F20.10,5X,'(
     FKCAL)')
 9000 CONTINUE
      END
