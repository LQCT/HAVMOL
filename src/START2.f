      SUBROUTINE ATOMS (LBASIS)

C READS AND EVALUATES ATOMIC PARAMETERS

      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      LOGICAL LBASIS
      LOGICAL  LADAPT,LSKIP,REVISE, bfl, first
      character*6 suber /'ATOMS'/
      CHARACTER*8 ATOM1, atom2
      character*80 kwd
      common
     ./optn1/ kwd(48)
     ./optn2/ icar, kb1, iii
      CHARACTER*4 ALABEL,TYPE
      COMMON/WORK/JREC,JUMP,ISTRT(40),INUMB(40),IWIDTH
      COMMON/BIG/SC,SCC(3),CSINP(NPRM1),CPINP(NPRM1),CDINP(NPRM1)
     &          ,INTYP(NSH1),NS(NSH1),EXX(6),CSS(6),CPP(6),CDD(6)
     &          ,ICT(NSH1,24),KS(NSH1),XJUNK(3333)
      CHARACTER*8 BFLAB,ANAM,ATNAME,PARNAM
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),
     4               KR(NATM1),KA(NATM1),KB(NATM1),
     5               PARVAL(NATM13),ZIN(NATM1),X(3,NATM1)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1)
     3            ,KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2
      COMMON/SYMTRY/T(216),
     &              INVT(24),NEWSH(NSH1+1,24),
     &              NT,NT2,
     &              Xold,Yold,Zold, XNEW,YNEW,ZNEW, XP,YP,ZP,
     &              U1,U2,U3, V1,V2,V3, W1,W2,W3, X0,Y0,Z0,
     &              INDEX,NAXIS
      CHARACTER*8 GROUP
      COMMON/SYMTRB/GROUP
      COMMON/STREAM/IR,IW231
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      real*4 nopt,minull
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR231,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/SECTOR/NUM3,IBLK3,REVISE,IODA(24),ICLASS(24),
     &              ILEN(24),MAXB,KBLKLA
      save /sector/
      COMMON/HUCKEL/HUCKEL(NSH1),ELEMS(3,18),ELEMP(2,18)
      CHARACTER*8 FINI,COMMNT,BLANK,BUFF,LIMIT
      character*32 basis
      DIMENSION NBFS(11),MINF(11),MAXF(11),ALABEL(11),NANGM(11),nhu(11)
      PARAMETER (ZERO=0.0D+00, THREE=3.D0, PI32=5.56832799683170D+00,
     &           DATA=1.0D-10,TOL=1.0D-10)

      DATA FINI/'END'/, COMMNT /'NOTE'/, LIMIT/'LIMIT'/
      DATA ALABEL/'1S','2S','2P','2SP','3S','3P','3SP','3D','S','P','D'/
      DATA NBFS  /  1 ,  1 ,  3 ,   4 ,  1 ,  3 ,   4 ,  6 , 1 , 3 , 6 /
      DATA MINF  /  1 ,  1 ,  2 ,   1 ,  1 ,  2 ,   1 ,  5 , 1 , 2 , 5 /
      DATA MAXF  /  1 ,  1 ,  4 ,   4 ,  1 ,  4 ,   4 , 10 , 1 , 4 ,10 /
      DATA NANGM /  1 ,  1 ,  2 ,   2 ,  1 ,  2 ,   2 ,  3 , 1 , 2 , 3 /
      data nhu   /  1 ,  2 ,  1 ,   1 ,  3 ,  2 ,   2 ,  0 , 1 , 1 , 0 /
      DATA BLANK /'        '/

10    FORMAT(/10X,20(1H*)//10X,'MOLECULAR BASIS SET'//
     1    10X,20(1H*))
20    FORMAT(//1H ,39(2H**)//1X,' ATOM     NUCLEAR',15X,
     1    'COORDINATES',17X,'NUMBER OF'/11X,'CHARGE',7X,'X',13X,
     2    'Y',14X,'Z',7X,'SHELLS'//1H ,39(2H**))
30    FORMAT(/1X,A8,2X,F5.1,3(F12.7,3X),I2)
40    FORMAT(/1H ,78(1H*)//' CONTRACTED PRIMITIVE FUNCTIONS'//
     1    ' Atom Shell Type   Primitive  Exponents   ',
     2    5X,'Contraction coefficients')
50    FORMAT(6X,I3,4X,A4,3X,I3,1X,2F15.6,3H  (,F10.6,2H ))
60    FORMAT(' EXCESSIVE NUMBER OF ATOMS')
70    FORMAT(' EXCESSIVE NUMBER OF SHELLS')
80    FORMAT(' EXCESSIVE NUMBER OF BASIS FUNCTIONS')
90    FORMAT(' EXCESSIVE NUMBER OF CONTRACTED PRIMITIVE FUNCTIONS')
100   FORMAT(' EXCESSIVE CONTRACTION INTO SHELL'/
     &' (more than 10 primitives in a contracted function)')
110   FORMAT(65X,3A4)
120   FORMAT(/1X,5H*****,' TRANSFORMATION TABLE OF ATOMS',
     1    5H*****/30X,' ROWS ARE ATOMS'/30X,
     2    ' COLUMNS ARE SYMMETRY OPERATIONS')
130   FORMAT(//)
140   FORMAT(1X,16A8)
150   FORMAT(1X,A8,15(2X,I2,3X,1H*))
160   FORMAT(1X,16(2X,I2,3X,1H*))
170   FORMAT(/1X,5H*****,' TRANSFORMATION TABLE OF SHELLS',
     1    5H*****/30X,' ROWS ARE SHELLS'/30X,
     2    ' COLUMNS ARE SYMMETRY OPERATIONS')
180   FORMAT(' THE CONTRACTED PRIMITIVE FUNCTIONS HAVE BEEN ',
     1    'UNNORMALIZED')
200   FORMAT(/)
210   FORMAT(/1X,5H*****,' INVERSE TRANSFORMATIONS ',5H*****/)
220   FORMAT(' THE CONTRACTED BASIS FUNCTIONS ARE NOW NORMALIZED TO',
     1    ' UNITY')
230   FORMAT(65X,I3,4X,A4,3X,I3,1X,2F14.6,2H (,F10.6,2H ),
     1    F14.6,2H (,F10.6,2H ))
240   FORMAT(/1X,A8)
931   format (/' THE BASIS FUNCTION LIBRARY ',a,' DOES NOT CONTAIN ',
     &        a,' ATOMS'/
     &' PLEASE, WRITE UP YOUR INPUT VALUES IN THE "<fn>.IMM" FILE FORMAT
     & FOR ATOMS'/)


C     READ IN UNIQUE CENTRES AND ATOMIC BASIS SETS GROUPED
C     IN SHELLS
C     GENERATE ALL NEW CENTERS
C     SET TABLE ( CENTERS VERSUS TRANSFORMATIONS )
C     SET TABLE ( SHELLS VERSUS TRANSFORMATIONS )

c LBASIS comes true when a geometry is changed from a previous calcu-
c lation and is planned to use the same basis set. Is is recommended
c to follow a reaction path.

      IF (LBASIS) NATOLD=NAT
      NAT=0
      NUQ=0
270   IF (.NOT.LBASIS) THEN
        DO 330 I=1,NPRM1
          EX(I)=ZERO
          CS(I)=ZERO
          CP(I)=ZERO
          CD(I)=ZERO
          CSINP(I)=ZERO
          CPINP(I)=ZERO
          CDINP(I)=ZERO
330     CONTINUE
        NSHELL=0
        LOC=0
        NGAUSS=0
      ENDIF

C     ----- READ IN UNIQUE CENTER -----

340   bfl = .false.
      first = .false.

c READAT reads atom identification and coordinates

      CALL READAT(ATOM1)

      IF(ATOM1.EQ.FINI) GO TO 590

      NAT=NAT+1
      IF(NAT.GT.NATM1) GOTO 801
      ANAM(NAT)=ATOM1
      IF (LBASIS .AND. (ZAN(NAT) .NE. ZIN(NUQ))) CALL HONDER(3,suber)
      ZAN(NAT)=ZIN(NUQ)
      DO 350 I=1,3
350     C(I,NAT)=X(I,NUQ)

c If the the former basis set is to be used, the new atom is read now

      IF (LBASIS) GO TO 340

C     ----- READ IN SHELLS -----

      ISHELL=0
      IGAUSS=0
      NS(NAT)=0
      KS(NAT)=NSHELL+1

c Here begins the loop to read each shell

360   CONTINUE
      if (kb1.eq.1) then
        call readc1
      else
        CALL READC
      endif
      CALL READA(BUFF)
      IF (BUFF .EQ. BLANK .OR. BUFF .EQ. COMMNT) GO TO 360
      IF(BUFF.EQ.FINI) GOTO 540
      if (.not.bfl) then
        call pathf (8, buff, '.BFN', basis)
        inquire (file=basis, exist=bfl)
        if (bfl) then
          open(25,file=basis,status='OLD')
          ir = 25
          first = .true.
        endif
      endif

      if (first) then
930     call readc
        call reada (atom2)
        call readf (znuq)
        if (atom2.eq.LIMIT) then
          write(10,931) basis,atom1
          write (iw,931) basis,atom1
          call honder (14,suber)
        endif
        if (zin(nuq).eq.znuq .or. atom1.eq.atom2) then
          call readc
          first = .false.
        else
920       if (atom2.ne.FINI) then
            call readc
            call reada (atom2)
            go to 920
          else
            go to 930
          endif
        endif
      endif

365   if (ishell.ne.0) then
        call readc
        CALL READA(BUFF)
        IF (BUFF .EQ. BLANK .OR. BUFF .EQ. COMMNT) GO TO 365
        IF(BUFF.EQ.FINI) GOTO 540
      endif
      JREC=0
      CALL READI(ISHELL)
      CALL READA4(TYPE)
      CALL READI(IGAUSS)
      CALL READF(sc)

      IF(IGAUSS.GT.10) THEN
        WRITE(IW,100)
        write(10,100)
        CALL HONDER (9,suber)
      ENDIF 

C - ITYP IS THE CODE FOR THE SHELL. IT HAS A FIXED NUMBER OF BASIS
C   FUNCTIONS ACCORDING DATAS /ALABEL/ AND CORRESPONDING /NBFS/

C - LOC ACOUNTS THE TOTAL NUMBER OF BASIS FUNCTIONS IN THE MOLECULE

C - NGAUSS ACOUNTS THE TOTAL AMOUNT OF GAUSSIAN PRIMITIVES

      DO 390 I=1,11
        IF(TYPE.EQ.ALABEL(I)) THEN
          ITYP = I
          GOTO 400
        ENDIF
390   CONTINUE
      CALL HONDER(4,suber)

400   NSHELL=NSHELL+1
      IF(NSHELL.GT.NSH1) GOTO 811

      KMIN(NSHELL)=MINF(ITYP)
      KMAX(NSHELL)=MAXF(ITYP)
      KSTART(NSHELL)=NGAUSS+1
      KATOM(NSHELL)=NAT
      KTY=NANGM(ITYP)
      KTYPE(NSHELL)=KTY
      INTYP(NSHELL)=ITYP
      KNG(NSHELL)=IGAUSS
      KLOC(NSHELL)=LOC+1
      NGAUSS=NGAUSS+IGAUSS
      LOC=LOC+NBFS(ITYP)
      iznat=zin(nat)
      if (kty.eq.1) then
        huckel(nshell) = elems(nhu(ityp),iznat)
      else if (kty.eq.2) then
        huckel(nshell) = elemp(nhu(ityp),iznat)
      else if (kty.ge.3) then
        huckel(nshell) = ZERO
      endif
      if (loc.gt.nbf1) goto 821
      K1=KSTART(NSHELL)
      K2=K1+KNG(NSHELL)-1
      IF(K2.GT.NPRM1)GOTO 831
      NS(NAT)=NS(NAT)+1

C     ----- GENERAL BASIS SET -----

*      IF(ITYP.LE.8) CALL HONDER(5,suber)
      if(ityp.eq.4 .or. ityp.eq.7) call honder(5,SUBER)

      DO 450 K=K1,K2
420     CALL READC
        CALL READA(BUFF)
        IF (BUFF .EQ. BLANK .OR. BUFF .EQ. COMMNT) GO TO 420
        JREC=0
        CALL READI(KDUM)
        CALL READF(EX(K))
        if (sc.ne.ZERO) ex(k) = ex(k)*sc*sc
        CALL READF(C1)
        CALL READF(C2)
440     continue
        IF(TYPE.EQ.ALABEL(1)) CSINP(K)=C1
        IF(TYPE.EQ.ALABEL(2)) CSINP(K)=C1
        IF(TYPE.EQ.ALABEL(3)) CPINP(K)=C1
        IF(TYPE.EQ.ALABEL(4)) CSINP(K)=C1
        IF(TYPE.EQ.ALABEL(4)) CPINP(K)=C2
        IF(TYPE.EQ.ALABEL(5)) CSINP(K)=C1
        IF(TYPE.EQ.ALABEL(6)) CPINP(K)=C1
        IF(TYPE.EQ.ALABEL(7)) CSINP(K)=C1
        IF(TYPE.EQ.ALABEL(7)) CPINP(K)=C2
        IF(TYPE.EQ.ALABEL(8)) CDINP(K)=C1
        IF(TYPE.EQ.ALABEL(9)) CSINP(K)=C1
        IF(TYPE.EQ.ALABEL(10)) CPINP(K)=C1
        IF(TYPE.EQ.ALABEL(11)) CDINP(K)=C1
        CS(K)=CSINP(K)
        CP(K)=CPINP(K)
        CD(K)=CDINP(K)
450     CONTINUE

C     IF(NORMP.NE.1) ... UNNORMALIZATION OF THE PRIMITIVE FUNCTIONS.
C     IF CONTRACTION COEFFICIENTS ARE GIVEN IN TERMS OF NORMALIZED
C     PRIMITIVE FUNCTIONS, CHANGE THEM TO GO WITH UNNORMALIZED
C     PRIMITIVES .

      IF (NORMP.EQ.1) GO TO 500
      DO 490 IG=K1,K2
        EE=EX(IG)+EX(IG)
        FACS=PI32/(EE*DSQRT(EE))
        A4=EE+EE
        FACP=FACS/A4
        FACD=THREE*FACP/A4
        CS(IG)=CS(IG)/DSQRT(FACS)
        CP(IG)=CP(IG)/DSQRT(FACP)
        CD(IG)=CD(IG)/DSQRT(FACD)
490     CONTINUE
500   CONTINUE

C     IF(NORMF.NE.1) NORMALIZE THE CONTRACTED BASIS FUNCTIONS.

      IF(NORMF.EQ.1) GO TO 365
      FACS=ZERO
      FACP=ZERO
      FACD=ZERO
      DO 520 IG=K1,K2
        DO 520 JG=K1,IG
          EE=EX(IG)+EX(JG)
          FAC=EE*DSQRT(EE)
          DUMS=CS(IG)*CS(JG)/FAC
          A4=EE+EE
          DUMP=CP(IG)*CP(JG)/(A4*FAC)
          DUMD=3.0*CD(IG)*CD(JG)/(A4*A4*FAC)
          IF(IG.EQ.JG) GO TO 510
          DUMS=DUMS+DUMS
          DUMP=DUMP+DUMP
          DUMD=DUMD+DUMD
510       FACS=FACS+DUMS
          FACP=FACP+DUMP
          FACD=FACD+DUMD
520       CONTINUE
      DO 530 IG=K1,K2
        IF(FACS.GT.TOL) CS(IG)=CS(IG)/DSQRT(FACS*PI32)
        IF(FACP.GT.TOL) CP(IG)=CP(IG)/DSQRT(FACP*PI32)
        IF(FACD.GT.TOL) CD(IG)=CD(IG)/DSQRT(FACD*PI32)
530     CONTINUE
      GO TO 365

540   if (bfl) then
        ir = 7
        close (25)
      endif

c Reading unique atom basis set is finished
c If no symmetry a new atom is to be read

      IF (NT .EQ. 1) GO TO 340

C     GENERATE EQUIVALENT CENTERS

      XS=U1*(C(1,NAT)-X0)+U2*(C(2,NAT)-Y0)+U3*(C(3,NAT)-Z0)
      YS=V1*(C(1,NAT)-X0)+V2*(C(2,NAT)-Y0)+V3*(C(3,NAT)-Z0)
      ZS=W1*(C(1,NAT)-X0)+W2*(C(2,NAT)-Y0)+W3*(C(3,NAT)-Z0)
      NAT0=NAT

C  Loop over symmetry-generated points

      DO 570 IT=2,NT
        NN=9*(IT-1)
        XNEW=XS*T(NN+1)+YS*T(NN+2)+ZS*T(NN+3)
        YNEW=XS*T(NN+4)+YS*T(NN+5)+ZS*T(NN+6)
        ZNEW=XS*T(NN+7)+YS*T(NN+8)+ZS*T(NN+9)
        XP=X0+U1*XNEW+V1*YNEW+W1*ZNEW
        YP=Y0+U2*XNEW+V2*YNEW+W2*ZNEW
        ZP=Z0+U3*XNEW+V3*YNEW+W3*ZNEW
        DO 550 IAT=1,NAT
          TEST=(XP-C(1,IAT))**2+(YP-C(2,IAT))**2+(ZP-C(3,IAT))**2
          IF(TEST.LE.TOL) GO TO 570
550       CONTINUE
        NUAT=NAT
        NAT=NAT+1
        IF(NAT.GT.NATM1)GOTO 801
        C(1,NAT)=C(1,NUAT)
        C(2,NAT)=C(2,NUAT)
        C(3,NAT)=C(3,NUAT)
        C(1,NUAT)=XP
        C(2,NUAT)=YP
        C(3,NUAT)=ZP
        IF (LBASIS) GO TO 570
        NS(NAT)=NS(NAT0)
        KS(NAT)=KS(NAT-1)+NS(NAT0)
        ANAM(NAT)=ANAM(NAT0)
        ZAN(NAT)=ZAN(NAT0)
        NSHELL=NSHELL+NS(NAT)
        NS1=NS(NAT)
        IF(NS1.EQ.0) GOTO 570
        DO 560 K=1,NS1
          J=KS(NAT0)+K-1
          JJ=KS(NAT)+K-1
          KMIN(JJ)=KMIN(J)
          KMAX(JJ)=KMAX(J)
          KSTART(JJ)=KSTART(J)
          KTYPE(JJ)=KTYPE(J)
          huckel(jj) = huckel(j)
          KNG(JJ)=KNG(J)
          KATOM(JJ)=NAT
          KLOC(JJ)=LOC+1
          INTYP(JJ)=INTYP(J)
          LOC=LOC+NBFS(INTYP(JJ))
560       if (loc.gt.nbf1) goto 821
570     CONTINUE
      GO TO 340

c Reading atoms is finished here

590   IF (LBASIS .AND. NAT .NE. NATOLD) CALL HONDER(3,suber)
      IF (LBASIS) RETURN

C     IF(NORMP.NE.1) THE CONTRACTED BASIS FUNCTIONS HAVE BEEN
C     EXPRESSED IN TERMS OF UNNORMALIZED PRIMITIVE FUNCTIONS.

      IF(NORMP.NE.1) WRITE(IW,180)

C     IF(NORMF.NE.1) THE CONTRACTED BASIS FUNCTIONS HAVE BEEN
C     NORMALIZED TO UNITY

      IF(NORMF.NE.1) WRITE(IW,220)

C     FORM TRANSFORMATION TABLES FOR ATOMS AND SHELLS.

      DO 640 IAT=1,NAT
        CX=C(1,IAT)
        CY=C(2,IAT)
        CZ=C(3,IAT)
        NS1=KS(IAT)-1
        NS2=NS(IAT)
        IF(NS2.EQ.0) GOTO 640
        XS=U1*(CX-X0)+U2*(CY-Y0)+U3*(CZ-Z0)
        YS=V1*(CX-X0)+V2*(CY-Y0)+V3*(CZ-Z0)
        ZS=W1*(CX-X0)+W2*(CY-Y0)+W3*(CZ-Z0)
        DO 630 IT=1,NT
          NN=9*(IT-1)
          XNEW=XS*T(NN+1)+YS*T(NN+2)+ZS*T(NN+3)
          YNEW=XS*T(NN+4)+YS*T(NN+5)+ZS*T(NN+6)
          ZNEW=XS*T(NN+7)+YS*T(NN+8)+ZS*T(NN+9)
          XP=X0+U1*XNEW+V1*YNEW+W1*ZNEW
          YP=Y0+U2*XNEW+V2*YNEW+W2*ZNEW
          ZP=Z0+U3*XNEW+V3*YNEW+W3*ZNEW
          DO 600 I=1,NAT
            TEST=(XP-C(1,I))**2+(YP-C(2,I))**2+(ZP-C(3,I))**2
            IF(TEST.GT.TOL) GO TO 600
            ICTR=I
            GO TO 610
600       CONTINUE
610       ICT(IAT,IT)=ICTR

C  ICT(IAT,IT) is the number of the atom into which atom IAT is taken
C  by symmetry operation IT.

          NS3=KS(ICTR)-1
          DO 620 ISH=1,NS2
            NEWSH(NS1+ISH,IT)=NS3+ISH
620       continue

C  NEWSH(NS,IT) is the number of the shell into which shell NS is taken
C  by symmetry operation IT.

630       CONTINUE
640     CONTINUE

C     NPRINT=1   ... EXTRA PRINTING FOR BASIS SET + SYMMETRY.
C     NPRINT=0   ... NORMAL PRINTING

      IF(NPRINT.EQ.6) GO TO 800
      WRITE(IW,10)
      WRITE(10,20)
      WRITE(IW,20)
      DO 650 IAT=1,NAT
        WRITE(IW,30) ANAM(IAT),ZAN(IAT),C(1,IAT),C(2,IAT),
     1               C(3,IAT),NS(IAT)
        WRITE(10,30) ANAM(IAT),ZAN(IAT),C(1,IAT),C(2,IAT),
     2               C(3,IAT),NS(IAT)
        IF(NS(IAT).EQ.0) GOTO 650
        NS1=KS(IAT)
        NS2=NS1+NS(IAT)-1
        WRITE(10,110) (ALABEL(INTYP(ISH)),ISH=NS1,NS2)
        WRITE(IW,110) (ALABEL(INTYP(ISH)),ISH=NS1,NS2)
650     CONTINUE
      WRITE(IW,40)
      DO 740 IAT=1,NAT
        DO 660 IT=1,NT
          IF(ICT(IAT,IT).GT.IAT) GO TO 740
660       CONTINUE
        WRITE(IW,240) ANAM(IAT)
        NS1=KS(IAT)
        NS2=NS1+NS(IAT)-1
        IF(NS(IAT).EQ.0) GOTO 740
        DO 730 ISH=NS1,NS2
          WRITE(IW,200)
          I1=KSTART(ISH)
          I2=I1+KNG(ISH)-1
          ITYP=INTYP(ISH)
          DO 720 IG=I1,I2
            GO TO (670,670,680,700,670,680,700,690,670,680,690) ,ITYP
670           C1=CS(IG)
              C2=CSINP(IG)
            GO TO 710
680           C1=CP(IG)
              C2=CPINP(IG)
            GO TO 710
690           C1=CD(IG)
              C2=CDINP(IG)
            GOTO 710
700           C1=CS(IG)
              C2=CSINP(IG)
              C3=CP(IG)
              C4=CPINP(IG)
              WRITE(IW,230) ISH,ALABEL(ITYP),IG,EX(IG),C1,C2,C3,C4
            GO TO 720
710           WRITE(IW,50) ISH,ALABEL(ITYP),IG,EX(IG),C1,C2
720         CONTINUE
730       CONTINUE
740     CONTINUE

800   CONTINUE

      IF(NAT.GT.NATM1) GO TO 801
      IF(NSHELL.GT.NSH1) GO TO 811
      IF(LOC.GT.NBF1) GO TO 821
      IF(NGAUSS.GT.NPRM1) GO TO 831

      NUM=LOC
      NX=NUM*(NUM+1)/2
      RETURN
801   WRITE(IW,60)
      CALL HONDER(7,suber)
811   WRITE(IW,70)
      CALL HONDER(10,suber)
821   WRITE(IW,80)
      CALL HONDER(8,suber)
831   WRITE(IW,90)
      CALL HONDER(6,suber)
      END

C     READAT M2

      SUBROUTINE READAT(ATOM)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'

C  Argument:
C  ATOM  Name of new atom

      DIMENSION ACONX(2)
      INTEGER IATCON(3), IAT1,IAT2,IAT3,P
      CHARACTER*8 ATOM
      COMMON/SYMTRY/T(216),
     &              INVT(24),NEWSH(NSH1+1,24),
     &              NT,NT2,
     &              XSMAL,YSMAL,ZSMAL, XNEW,YNEW,ZNEW, XP,YP,ZP,
     &              U1,U2,U3, V1,V2,V3, W1,W2,W3, X0,Y0,Z0,
     &              INDEX,NAXIS
      CHARACTER*8 GROUP
      COMMON/SYMTRB/GROUP
      CHARACTER*8 ANAM,BFLAB,ATNAME,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P(NATM1,3),KR(NATM1),KA(NATM1),KB(NATM1),
     4               PARVAL(NATM13),ZIN(NATM1),X(3,NATM1),
     .       /optn1/ kwd(48)
     ./optn2/ icar, kb1, iii
      character*80 kwd
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      CHARACTER*8 ACONX,COMMNT,FINI,SEP,BLANK,BCONX
      EQUIVALENCE (IAT1,IATCON(1)),(IAT2,IATCON(2)),(IAT3,IATCON(3))

      DATA COMMNT /'NOTE'/
      DATA ACONX /'XYZ','POL'/
      DATA FINI /'END'/, SEP /'NEXT'/, BLANK /'        '/

C  Initialize

      IF (NUQ .EQ. 0) NVAL=NPAR

C  Read atom name and nuclear charge

10    if (kb1.eq.1) then
        call readc1
      else
        CALL READC
      endif
      CALL READA(ATOM)
      IF (ATOM .EQ. BLANK .OR. ATOM .EQ. COMMNT) GO TO 10
      IF (ATOM .EQ. FINI .OR. ATOM .EQ. SEP) RETURN
      NUQ=NUQ+1
      ATNAME(NUQ)=ATOM
      CALL READF(ZNUC)
      ZIN(NUQ)=ZNUC

C  Read connexion type

      CALL READA(BCONX)

C  Possible connexions are:
C          1. Absolute cartesian coordinates        'XYZ' (default)
C          2. Absolute polar coordinates            'POL'

      KCONX=0
      DO 20 K=1,2
20      IF (BCONX .EQ. ACONX(K)) KCONX=K
      KONTYP(NUQ)=KCONX
      IF (KCONX .EQ. 0) CALL REREAD(-1)
      IF (KCONX .EQ. 2) GO TO 30

C  Read in cartesian coordinates

      CALL READP(X(1,NUQ),KR(NUQ),CFR)
      CALL READP(X(2,NUQ),KA(NUQ),CFR)
      CALL READP(X(3,NUQ),KB(NUQ),CFR)
      GO TO 210

C  Read in polar coordinates

30    CALL READP(R,KR(NUQ),CFR)
      CALL READP(THETA,KA(NUQ),CFA)
      CALL READP(PHI,KB(NUQ),CFA)
      X(3,NUQ)=R*DCOS(THETA)
      X(1,NUQ)=R*DSIN(THETA)*DCOS(PHI)
      X(2,NUQ)=R*DSIN(THETA)*DSIN(PHI)
210   IF (ZIN(NUQ) .LT. 0) GO TO 10
      RETURN
      END

      subroutine readc1

      character*80 kwd
      common
     ./optn1/ kwd(48)
     ./optn2/ icar, kb1, ii
      character*1 ia
      logical initem
      common/stream/ir,iw
      common/work/item,nitems,istrt(40),inumb(40),iwidth
      common/workb/ia(80)
     ./answer/ si,bl,hyphen,star,plus,comma,equal
      character si,bl,hyphen,star,plus,comma,equal

c ii is the number of the current line in the internal file
c NITEMS is the amount of words in the line
c ISTRT(i) is the column where word 'i' begins
c INUMB(i) is the amount of characters in word 'i'
c IWIDTH was set to 80 columns in BLOCK DATA

      ii = ii + 1
      nitems=0
      item=0
      initem=.false.

      read (kwd(ii),'(80a1)',end=50) IA

      do 15 j=1,80
        if (ichar(IA(j)).eq.10 .or. ichar(IA(j)).eq.13)
     &    IA(j) = char(32)
        do 15 i=97,122
          ich = ichar(IA(j))
          if (ich .eq. i) IA(j)=CHAR(i-32)
15        continue

      do 40 i=1,iwidth
        if(ia(i).eq.bl .or.
     &     ia(i).eq.equal .or.
     &     ia(i).eq.comma) then
          initem=.false.
        else
          if (initem) then
            inumb(nitems)=inumb(nitems)+1
          else
            nitems=nitems+1
            istrt(nitems)=i
            inumb(nitems)=1
            initem=.true.
          endif
        endif
40      continue
      return

50    write(10,100)
      write (iw,100)
100   format
     & (/' *** MORE THAN EXPECTED LINES IN INTERNAL FILE <KWD> ***'/)
      CALL HONDER (16,'READC1')
      end

C     PTGRP   M2

      SUBROUTINE PTGRP(NPRINT,IW)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      CHARACTER*6 SUBER /'PTGRP'/
      CHARACTER*8 GRP,DRC,FINI,DIRECT,GUFF,BUFF
      COMMON/BLKIN/NIRR,MULT(8,8),IRR(NBF1),AGROUP,XJNK(446)
      COMMON/SYMTRY/T(216),
     &              INVT(24),NEWSH(NSH1+1,24),
     &              NT,NT2,
     &              XSMAL,YSMAL,ZSMAL, XNEW,YNEW,ZNEW, XP,YP,ZP,
     &              U1,U2,U3, V1,V2,V3, W1,W2,W3, X0,Y0,Z0,
     &              INDEX,NAXIS
      CHARACTER*8 GROUP
      COMMON/SYMTRB/GROUP
      COMMON/WORK/JREC,JUMP,ISTRT(40),INUMB(40),IWIDTH
      CHARACTER*1 IA
      COMMON/WORKB/IA(80)

      DIMENSION GRP(16),DRC(2)
      PARAMETER (ZERO=0.0D+00,PT5=0.5D+00,ONE=1.0D+00,TOL=1.0D-10,
     &           PI2=6.28318530717958D+00)

      DATA FINI/'END'/
      DATA GRP/'C1'      ,'CS'      ,'CI'      ,'CN'      ,'S2N'     ,
     &         'CNH'     ,'CNV'     ,'DN'      ,'DNH'     ,'DND'     ,
     &         'CINFV'   ,'DINFH'   ,'T'       ,'TH'      ,'TD'      ,
     &         'O'       /
      DATA DRC/'NORMAL','PARALLEL'/

10    FORMAT(' LINEAR MOLECULE , POINT GROUP IS CINFV OR DINFH ',/,
     1 ' PLEASE USE GROUP CNV OR DNH...')
20    FORMAT(' ILLEGAL POINT GROUP INPUT. STOP..')
30    FORMAT(/,' THE POINT GROUP OF THE MOLECULE IS ...',A8,/,
     1 ' THE ORDER OF THE PRINCIPAL AXIS IS ...',I3)
40    FORMAT(/,' THE ORIGIN OF THE LOCAL FRAME IS AT X =  ',F10.5,
     1 ' Y = ',F10.5,' Z = ',F10.5,/,
     2 ' DIRECTOR COSINES OF THE NEW AXES ',/,34X,3(5X,F10.5),/,34X,
     3 3(5X,F10.5),/,34X,3(5X,F10.5))
50    FORMAT(' ROTATIONS ABOUT PRINCIPAL AXIS')
60    FORMAT(' SIGMA-H FOLLOWED BY ROTATIONS')
70    FORMAT(' C2 FOLLOWED BY ROTATIONS ')
80    FORMAT(' SIGMA-D FOLLOWED BY ROTATIONS')
90    FORMAT(' SIGMA-V FOLLOWED BY ROTATIONS')
100   FORMAT(/,10X,' CENTER OF SYMMETRY AT X = ',F10.5,' Y = ',F10.5,
     1 ' Z = ',F10.5)
110   FORMAT(/,' PLANE OF SYMMETRY DEFINED BY ITS NORMAL U = ',F10.5,
     1 ' V = ',F10.5,' W = ',F10.5)
120   FORMAT(/,10X,3F15.9,/,10X,3F15.9,/,10X,3F15.9)
130   FORMAT(' C2 FOLLOWED BY SIGMA-H FOLLOWED BY ROTATIONS')
140   FORMAT(' SIGMA-D FOLLOWED BY C2 FOLLOWED BY ROTATIONS')
150   FORMAT(' S2N ROTATION FOLLOWED BY ROTATIONS')

160   CALL READC
      CALL READA(GROUP)
      CALL READI(NAXIS)
170   INDEX=17
      DO 180 I=1,16
180     IF(GROUP.EQ.GRP(I)) INDEX=I
      WRITE(IW,30) GROUP,NAXIS
      IF(INDEX.LE.17) GO TO 190
      WRITE(IW,20)
      CALL HONDER(2,suber)
190   IF(INDEX.LE.3) GO TO 410
      IF(INDEX.EQ.11.OR.INDEX.EQ.12) NAXIS=2
      IF(INDEX.EQ.11) INDEX=7
      IF(INDEX.EQ.12) INDEX=9
C     DEFINE LOCAL FRAME
C     READ IN PRINCIPAL AXIS   ( 1 CARD )
C     READ IN X-LOCAL AXIS   ( 1 CARD )
C     DEFAULT OPTION : LOCAL FRAME IDENTICAL TO MASTER FRAME
      CALL READC
      CALL READA(GUFF)
      IF(GUFF.EQ.FINI) GOTO 200
      JREC=0
      CALL READF(X0)
      CALL READF(Y0)
      CALL READF(Z0)
      CALL READF(X1)
      CALL READF(Y1)
      CALL READF(Z1)
      CALL READC
      CALL READA(GUFF)
      IF(GUFF.EQ.FINI) CALL HONDER(2,suber)
      JREC=0
      CALL READF(X2)
      CALL READF(Y2)
      CALL READF(Z2)
      CALL READA(DIRECT)
      RHO=DSQRT((X1-X0)**2+(Y1-Y0)**2+(Z1-Z0)**2)
      CALL READC
      CALL READA(GUFF)
      IF(GUFF.NE.FINI) CALL HONDER(2,suber)
      GOTO 210
200   X0=ZERO
      Y0=ZERO
      Z0=ZERO
      X1=ZERO
      Y1=ZERO
      Y2=ZERO
      Z2=ZERO
      Z1=ONE
      X2=ONE
      DIRECT=DRC(2)
      RHO=ONE
210   CONTINUE
      IF(DIRECT.NE.DRC(1)) DIRECT=DRC(2)
      W1=(X1-X0)/RHO
      W2=(Y1-Y0)/RHO
      W3=(Z1-Z0)/RHO
      WW=W1*W1+W2*W2+W3*W3
      X02=X2-X0
      Y02=Y2-Y0
      Z02=Z2-Z0
      RHO=(W1*X02+W2*Y02+W3*Z02)/WW
      DUM=RHO*W1
      X0=X0+DUM
      X02=X02-DUM
      DUM=RHO*W2
      Y0=Y0+DUM
      Y02=Y02-DUM
      DUM=RHO*W3
      Z0=Z0+DUM
      Z02=Z02-DUM
      UU=(X02*X02+Y02*Y02+Z02*Z02)
      U=DSQRT(UU)
      U1=X02/U
      U2=Y02/U
      U3=Z02/U
      V3=W1*U2-W2*U1
      V2=W3*U1-W1*U3
      V1=W2*U3-W3*U2
      IF(DIRECT.EQ.DRC(2)) GO TO 220
      DUM=U1
      U1=V1
      V1=-DUM
      DUM=U2
      U2=V2
      V2=-DUM
      DUM=U3
      U3=V3
      V3=-DUM
220   CONTINUE
      IF(NPRINT.NE.1) GO TO 230
      WRITE(IW,40) X0,Y0,Z0,U1,V1,W1,U2,V2,W2,U3,V3,W3
230   CONTINUE
      IF(INDEX.GE.13) GO TO 410
C     ROTATION ABOUT PRINCIPAL AXIS
      NN=0
      N=NAXIS
      ALPHA=ZERO
      ALPH=PI2/FLOAT(N)
240   NN=NN+1
      IF(NN.GT.N) GO TO 250
      COSA=DCOS(ALPHA)
      SINA=DSIN(ALPHA)
      I=9*(NN-1)
      T(I+1)=COSA
      T(I+5)=COSA
      T(I+2)=-SINA
      T(I+4)=SINA
      T(I+3)=ZERO
      T(I+6)=ZERO
      T(I+7)=ZERO
      T(I+8)=ZERO
      T(I+9)=ONE
      ALPHA=ALPHA+ALPH
      GO TO 240
C     END OF GROUP 4
250   NT=N
      II=9*NT
      IF(NPRINT.NE.1) GO TO 260
      WRITE(IW,50)
      N1=1
      N2=NAXIS
      CALL PRINT(N1,N2,IW)
260   CONTINUE
      IF(INDEX.EQ.4) GO TO 630
      IF(INDEX.EQ.5) GO TO 500
      IF(INDEX.EQ.7) GO TO 380
      IF(INDEX.NE.6.AND.INDEX.NE.9) GO TO 310
C     SIGMA-H PLANE  EQUATION (Z=0) IN LOCAL FRAME
      NN=0
270   NN=NN+1
      IF(NN.GT.NT) GO TO 290
C     GROUP 6 0R 9
      I=II+9*(NN-1)
      DO 280 J=1,8
280     T(I+J)=T(I+J-II)
      T(I+9)=-T(I+9-II)
      GO TO 270
290   NT=NT+NT
      II=9*NT
      IF(NPRINT.NE.1) GO TO 300
      WRITE(IW,60)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT(N1,N2,IW)
300   CONTINUE
C     END OF GROUP 6
      IF(INDEX.EQ.6) GO TO 630
C     ONE CP2 AXIS IS THE X-AXIS OF THE LOCAL FRAME
C     GROUP 8 , 9 ,10
310   CONTINUE
      NN=0
320   NN=NN+1
      IF(NN.GT.NT) GO TO 330
      I=II+9*(NN-1)
      T(I+1)=T(I+1-II)
      T(I+2)=-T(I+2-II)
      T(I+3)=-T(I+3-II)
      T(I+4)=T(I+4-II)
      T(I+5)=-T(I+5-II)
      T(I+6)=-T(I+6-II)
      T(I+7)=T(I+7-II)
      T(I+8)=-T(I+8-II)
      T(I+9)=-T(I+9-II)
      GO TO 320
330   NT=NT+NT
      II=9*NT
      IF(NPRINT.NE.1) GO TO 340
      WRITE(IW,70)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT(N1,N2,IW)
      IF(INDEX.NE.9) GO TO 340
      WRITE(IW,130)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT(N1,N2,IW)
340   CONTINUE
C     END OF GROUP 8 AND 9
      IF(INDEX.EQ.8.OR.INDEX.EQ.9) GO TO 630
C     DND GROUP . EQUATION OF PLANE SIGMA-D IS :
C     DSIN(ALPH/4)*X-DCOS(ALPH/4)*Y=0
C     THE X-AXIS IS THE CP2 AXIS.
C     GROUP 10
      BETA=PT5*ALPH
      COSA=DCOS(BETA)
      SINA=DSIN(BETA)
      NN=0
350   NN=NN+1
      IF(NN.GT.NT) GO TO 360
      I=II+9*(NN-1)
      T(I+1)=COSA*T(I+1-II) + SINA*T(I+2-II)
      T(I+2)=SINA*T(I+1-II) - COSA*T(I+2-II)
      T(I+3)=     T(I+3-II)
      T(I+4)=COSA*T(I+4-II) + SINA*T(I+5-II)
      T(I+5)=SINA*T(I+4-II) - COSA*T(I+5-II)
      T(I+6)=     T(I+6-II)
      T(I+7)=COSA*T(I+7-II) + SINA*T(I+8-II)
      T(I+8)=SINA*T(I+7-II) - COSA*T(I+8-II)
      T(I+9)=     T(I+9-II)
      GO TO 350
360   NT=NT+NT
      II=9*NT
      IF(NPRINT.NE.1) GO TO 370
      WRITE(IW,80)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT(N1,N2,IW)
      WRITE(IW,140)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT(N1,N2,IW)
370   CONTINUE
C     END OF GROUP 10
      GO TO 630
C     GROUP 7
C     SIGMA-V IS THE (X-Z) PLANE OF LOCAL FRAME
380   CONTINUE
      NN=0
390   NN=NN+1
      IF(NN.GT.NT) GO TO 400
      I=II+9*(NN-1)
      T(I+1)=T(I+1-II)
      T(I+2)=-T(I+2-II)
      T(I+3)=T(I+3-II)
      T(I+4)=T(I+4-II)
      T(I+5)=-T(I+5-II)
      T(I+6)=T(I+6-II)
      T(I+7)=T(I+7-II)
      T(I+8)=-T(I+8-II)
      T(I+9)=T(I+9-II)
      GO TO 390
400   NT=NT+NT
      II=9*NT
C     END OF GROUP 7
      IF(NPRINT.NE.1) GO TO 630
      WRITE(IW,90)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT(N1,N2,IW)
      GO TO 630
410   CONTINUE
      T(1)=ONE
      T(5)=ONE
      T(9)=ONE
      T(2)=ZERO
      T(3)=ZERO
      T(4)=ZERO
      T(6)=ZERO
      T(7)=ZERO
      T(8)=ZERO
      IF(INDEX.EQ.1) GO TO 420
      IF(INDEX.EQ.2) GO TO 430
      IF(INDEX.EQ.3) GO TO 460
      IF(INDEX.EQ.11.OR.INDEX.EQ.12) GO TO 490
      GO TO 530
420   NT=1
      X0=ZERO
      Y0=ZERO
      Z0=ZERO
      U1=ONE
      V2=ONE
      W3=ONE
      U2=ZERO
      U3=ZERO
      V1=ZERO
      V3=ZERO
      W1=ZERO
      W2=ZERO
      CALL READC
      CALL READA(BUFF)
      IF(BUFF.NE.FINI) CALL HONDER(2,suber)
      GO TO 630
C     CS SYMMETRY GROUP
C     THE 3 POINTS 1,2,3 DEFINE SIGMA-H PLANE
430   CONTINUE
      CALL READC
      CALL READA(GUFF)
      IF(GUFF.EQ.FINI) GOTO 440
      JREC=0
      CALL READF(X1)
      CALL READF(Y1)
      CALL READF(Z1)
      CALL READF(X2)
      CALL READF(Y2)
      CALL READF(Z2)
      CALL READC
      CALL READA(GUFF)
      IF(GUFF.EQ.FINI) CALL HONDER(2,suber)
      JREC=0
      CALL READF(X3)
      CALL READF(Y3)
      CALL READF(Z3)
      RHO=(X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2
      CALL READC
      CALL READA(GUFF)
      IF(GUFF.NE.FINI) CALL HONDER(2,suber)
      GOTO 450
C     DEFAULT OPTION : PLANE IS THE (X,Z) PLANE
440   X1=ZERO
      Y1=ZERO
      Z1=ZERO
      X2=ZERO
      Y2=ZERO
      Z2=ONE
      X3=ONE
      Y3=ZERO
      Z3=ZERO
450   CONTINUE
      NT=2
      W1=(Y2-Y1)*(Z3-Z1)-(Y3-Y1)*(Z2-Z1)
      W2=(Z2-Z1)*(X3-X1)-(Z3-Z1)*(X2-X1)
      W3=(X2-X1)*(Y3-Y1)-(X3-X1)*(Y2-Y1)
      RHO=DSQRT(W1*W1+W2*W2+W3*W3)
      W1=W1/RHO
      W2=W2/RHO
      W3=W3/RHO
      U1=X2-X1
      U2=Y2-Y1
      U3=Z2-Z1
      RHO=DSQRT(U1*U1+U2*U2+U3*U3)
      U1=U1/RHO
      U2=U2/RHO
      U3=U3/RHO
      V1=W2*U3-W3*U2
      V2=W3*U1-W1*U3
      V3=W1*U2-W2*U1
      X0=X1
      Y0=Y1
      Z0=Z1
      T(10)=ONE
      T(14)=ONE
      T(18)=-ONE
      T(11)=ZERO
      T(12)=ZERO
      T(13)=ZERO
      T(15)=ZERO
      T(16)=ZERO
      T(17)=ZERO
      IF(NPRINT.NE.1) GO TO 630
      WRITE(IW,110) W1,W2,W3
      WRITE(IW,120) U1,V1,W1,U2,V2,W2,U3,V3,W3
      GO TO 630
C     CI SYMMETRY GROUP
C     CENTER OF INVERSION IS (X0,Y0,Z0)
460   X0=ZERO
      Y0=ZERO
      Z0=ZERO
470   CALL READC
      CALL READA(GUFF)
      IF(GUFF.EQ.FINI) GOTO 480
      JREC=0
      CALL READF(X0)
      CALL READF(Y0)
      CALL READF(Z0)
      CALL READC
      CALL READA(GUFF)
      IF(GUFF.NE.FINI) CALL HONDER(2,suber)
480   IF(NPRINT.EQ.3) WRITE(IW,100) X0,Y0,Z0
      T(10)=-ONE
      T(14)=-ONE
      T(18)=-ONE
      T(11)=ZERO
      T(12)=ZERO
      T(13)=ZERO
      T(15)=ZERO
      T(16)=ZERO
      T(17)=ZERO
      NT=2
      U1=ONE
      V2=ONE
      W3=ONE
      U2=ZERO
      U3=ZERO
      V1=ZERO
      V3=ZERO
      W1=ZERO
      W2=ZERO
      GO TO 630
490   WRITE(IW,10)
      write(10,10)
      call honder(19,suber)

500   NN=0
      BETA=PT5*ALPH
      COSB=DCOS(BETA)
      SINB=DSIN(BETA)
510   NN=NN+1
      IF(NN.GT.NT) GO TO 520
C     S2N GROUP
C     THE PLANE OF SYMMETRY FOR THE IMPROPER ROTATION
C     IS THE (X,Y) PLANE OF THE LOCAL FRAME
      I=II+9*(NN-1)
      T(I+1)= COSB*T(I+1-II)+SINB*T(I+2-II)
      T(I+2)=-SINB*T(I+1-II)+COSB*T(I+2-II)
      T(I+3)=     -T(I+3-II)
      T(I+4)= COSB*T(I+4-II)+SINB*T(I+5-II)
      T(I+5)=-SINB*T(I+4-II)+COSB*T(I+5-II)
      T(I+6)=     -T(I+6-II)
      T(I+7)= COSB*T(I+7-II)+SINB*T(I+8-II)
      T(I+8)=-SINB*T(I+7-II)+COSB*T(I+8-II)
      T(I+9)=     -T(I+9-II)
      GO TO 510
520   NT=NT+NT
      II=9*NT
      IF(NPRINT.NE.1) GO TO 630
      WRITE(IW,150)
      N1=N2+1
      N2=N2+NAXIS
      CALL PRINT(N1,N2,IW)
      GO TO 630
C     T GROUP AND OTHERS CONTAINING A SUBGROUP T :
C     LOCAL X,Y,Z ARE THE C2 AXES
530   DO 540 I=10,36
540     T(I)=ZERO
      T(10)=ONE
      T(23)=ONE
      T(36)=ONE
      T(14)=-ONE
      T(18)=-ONE
      T(19)=-ONE
      T(27)=-ONE
      T(28)=-ONE
      T(32)=-ONE
      DO 550 II=5,12
        I=9*(II-1)
        J=9*(II-5)
        T(I+1)=T(J+7)
        T(I+2)=T(J+8)
        T(I+3)=T(J+9)
        T(I+4)=T(J+1)
        T(I+5)=T(J+2)
        T(I+6)=T(J+3)
        T(I+7)=T(J+4)
        T(I+8)=T(J+5)
550     T(I+9)=T(J+6)
      NT=12
      IF(INDEX.EQ.13) GO TO 630
      IF(INDEX.EQ.14) GO TO 560
      IF(INDEX.EQ.15) GO TO 590
      GO TO 580
C     TH GROUP
C     EXPAND GROUP BY TAKING DIRECT PRODUCT WITH CI
560   I=9*NT
      DO 570 J=1,I
570     T(J+I)=-T(J)
      NT=NT+NT
      GO TO 630
C     OCTAHEDRAL GROUP IS DIRECT PRODUCT OF T AND A C4 ROTATION
C     ABOUT Z AXIS
580   SIGN=-ONE
      GO TO 600
C     TD GROUP IS DIRECT PRODUCT OF T AND A REFLECTION IN A
C     PLANE ( EQUATION OF THE PLANE   X=Y  )
590   SIGN=ONE
600   DO 610 II=13,24
        I=9*(II-1)
        J=9*(II-13)
        T(I+1)=T(J+4)*SIGN
        T(I+2)=T(J+5)*SIGN
        T(I+3)=T(J+6)*SIGN
        T(I+4)=T(J+1)
        T(I+5)=T(J+2)
        T(I+6)=T(J+3)
        T(I+7)=T(J+7)
        T(I+8)=T(J+8)
610     T(I+9)=T(J+9)
      NT=24
630   CONTINUE
C     FIND THE INVERSE TRANSFORMATIONS
      DO 650 ITR=1,NT
        NN=9*(ITR-1)
        DO 640 IT=1,NT
          II=9*(IT-1)
          TEST= T(NN+1)*T(II+1)+T(NN+2)*T(II+4)+T(NN+3)*T(II+7)
     1         +T(NN+4)*T(II+2)+T(NN+5)*T(II+5)+T(NN+6)*T(II+8)
     2         +T(NN+7)*T(II+3)+T(NN+8)*T(II+6)+T(NN+9)*T(II+9)
     3         -3.0D+00
          IF(DABS(TEST).GT.TOL) GO TO 640
          INVT(ITR)=IT
          GO TO 650
640       CONTINUE
650     CONTINUE
      RETURN
      END

      SUBROUTINE READA4(BUFF)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 IA,IBLNK,IBUF
      CHARACTER*4 BUFF
      DIMENSION IBUF(4)
      COMMON/WORK/JREC,JUMP,ISTRT(40),INUMB(40),IWIDTH
      COMMON/WORKB/IA(80)
      DATA IBLNK/' '/

      DO 10 I=1,4
10      IBUF(I)=IBLNK
      JREC=JREC+1
      IF(JREC.GT.JUMP) GOTO 30
      N=INUMB(JREC)
      NSTRT=ISTRT(JREC)
      IF(N.GT.4) N=4
      DO 20 I=1,N
        IBUF(I)=IA(NSTRT)
20      NSTRT=NSTRT+1
30    CALL PACK4(BUFF,IBUF)

      RETURN
      END

      SUBROUTINE PACK4(BUFF,IBUF)
      CHARACTER*4 A,BUFF
      CHARACTER*1 IBUF(4),LBUF(4)
      EQUIVALENCE (A,LBUF)

      DO 10 I=1,4
10    LBUF(I)=IBUF(I)
      BUFF=A

      RETURN
      END

      SUBROUTINE GETFIL(NFILE,IBLK,LBLK,NOTAPE)
      CHARACTER*1 IA
      CHARACTER*4 ITEST
      character*6 suber /'GETFIL'/
      COMMON /STREAM/ IR,IW
     ./files / file5
      COMMON/WORK/JREC,JUMP,INUMB(80),IWIDTH
      COMMON/WORKB/IA(80)
      COMMON/DISC/ISELL(5),IPOS(16),IBLKSZ(16)
      save /disc/
      COMMON/DISCB/NAM(16)
      character*4 nam
      character*32 file5, file

      IF(JUMP.GT.13)JUMP=13
      GOTO(10,40,40,40,50,50,50,60,60,60,70,70,70),JUMP
10    WRITE(IW,20)
      write(10,20) 
20    FORMAT(//1X,'No files specified')
      WRITE(IW,30) IA
      write(10,30) IA
30    FORMAT(/1X,'Last card read '/1X,80A1)
      call honder(20,suber)

40    NFILE=1
      GOTO 80
50    NFILE=2
      GOTO 80
60    NFILE=3
      GOTO 80
70    NFILE=4
80    file = file5
      DO 160 I=1,NFILE
        CALL READA4(ITEST)
        DO 90 ISEL=1,16
          IF(ITEST.EQ.nam(isel))GOTO 110
90      CONTINUE
        WRITE(IW,100)
        write(10,100)
100     FORMAT(//1X,'Illegal file name')
        WRITE(IW,30) IA
        write(10,30) IA
        call honder(21,suber)

110     call filen1(1,file,ITEST)
        IF(ISEL.NE.8)GOTO 130
        WRITE(IW,120)
        write(10,120)
120     FORMAT(//1X,'<fn>.ED7 is reserved for the Scratchfile')
        WRITE(IW,30) IA
        write(10,30) IA
        call honder(22,suber)

130   IF(IPOS(ISEL).GE.0)GOTO 150
      WRITE(IW,140)
      write(10,140)
140   FORMAT(/1X,'File to be restored has not been opened')
      WRITE(IW,30) IA
      write(10,30) IA
      call honder(23,suber)

150   NOTAPE=ISEL
      CALL READI(IIII)
      IF(IIII.EQ.0)IIII=1
      IBLK=IIII
      CALL READI(LBLK)
160   CONTINUE
      RETURN
      END

C     REREAD  M2

      SUBROUTINE REREAD(K)

C  K positive: reset input so as to reread the Kth item on the current
C  card.
C  K negative: move back -K items (e.g. with K=-1 read the last item aga

      CHARACTER*1 IA
      COMMON/WORK/JREC,JUMP,ISTRT(40),INUMB(40),IWIDTH
      COMMON/WORKB/IA(80)

      IF (K .GT. 0) JREC=K-1
      IF (K .LT. 0) JREC=JREC+K

      RETURN
      END

      FUNCTION RIJ(C,I,J)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      DIMENSION C(3,NSH1)
      RIJSQ=(C(1,I)-C(1,J))**2+(C(2,I)-C(2,J))**2+(C(3,I)-C(3,J))**2
      RIJ=DSQRT(RIJSQ)
      RETURN
      END

      SUBROUTINE READP(Q,K,CF)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 BLANK,S
      CHARACTER*8 B,PARNAM,ATNAME,ANAM,BFLAB
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),
     4               KR(NATM1),KA(NATM1),KB(NATM1),
     5               PARVAL(NATM13),ZIN(NATM1),X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)

      DATA S /'-'/, BLANK /'        '/

      KS=-1
10    KS=-KS
      CALL READA(B)
      IF (B .EQ. S) GO TO 10
20    DO 30 I=1,NPAR
        IF (B .EQ. PARNAM(I)) GO TO 50
30    CONTINUE
      CALL REREAD(-1)
      CALL READF(Q)
      IF (Q .EQ. 0D0) GO TO 60
      IF (Q .LT. 0D0) KS=-KS
      Q=DABS(Q)*CF
      N1=NVAR+1
      DO 40 I=N1,NVAL
        IF (Q .EQ. PARVAL(I)) GO TO 50
40      CONTINUE

C  No match - insert new table entry

      NVAL=NVAL+1
      IF (NVAL .GT. NATM13) CALL HONDER(13,'READP ')
      PARNAM(NVAL)=BLANK
      PARVAL(NVAL)=Q
      I=NVAL
50    K=I*KS
      Q=KS*PARVAL(I)
      RETURN

60    K=0
      RETURN
      END

      SUBROUTINE CHKFIL(IFIL)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*4 NAM
      COMMON/DISC/ISEL(5),IPOS(32)
      save /disc/
      COMMON/DISCB/NAM(16)

      IF(IPOS(IFIL).GE.0)RETURN
      WRITE(9,10) NAM(IFIL)
      write(10,10) NAM(IFIL)
10    FORMAT(//10X,'FILE ',A4,' HAS NOT BEEN OPENED')
      call honder(23,'CHKFIL')
      END

      SUBROUTINE PACK2(CHARP,K,L)
      CHARACTER*8 A,CHARP
      CHARACTER*4 I(2),K,L
      EQUIVALENCE (A,I(1))
      I(1)=K
      I(2)=L
      CHARP=A
      RETURN
      END

      SUBROUTINE MOVORG
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
C     MOVE ORIGIN        R D AMOS      1986
c                        L.A. Montero  1992
      CHARACTER*8 PARNAM,ATNAME,ANAM,BFLAB
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),KR(NATM1),
     5               KA(NATM1),KB(NATM1),
     4               PARVAL(NATM13),ZIN(NATM1),XUNK(NATM13)
      common/origin/ itype, jcar, x, y, z
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),
     1               ATNAME(NATM1)
      COMMON/ATMWT/AMS(37),AMASS(NATM1)
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
     ./files / file5
      character*32 file5
      character*17 com_t, coc_t, ttt
      logical*1 carw
      logical carf
      DATA BOHR/0.52917706D0/
      data com_t /' center of mass'/, coc_t /' center of charge'/,
     &       ttt /'geometry moved to'/

      GOTO(21,22,23),ITYPE

c center of mass

21    X=0.0D0
      Y=0.0D0
      Z=0.0D0
      AT=0.0D0
      DO 31 N=1,NAT
        IZ=IDINT(ZAN(N)+0.1D0)
        IF(IZ.LT.1)GOTO 31
        AM=AMS(IZ)
        AT=AT+AM
        X=X+AM*C(1,N)
        Y=Y+AM*C(2,N)
        Z=Z+AM*C(3,N)
31    CONTINUE
      X=X/AT
      Y=Y/AT
      Z=Z/AT
      XA=X*BOHR
      YA=Y*BOHR
      ZA=Z*BOHR
      WRITE(9,67) X,Y,Z,XA,YA,ZA
      WRITE(10,67) X,Y,Z,XA,YA,ZA
67    FORMAT(/1X,'Current position of centre of mass is'
     & /10X,3F16.8,'  atomic units'
     & /10X,3F16.8,'  angstrom')
      GOTO 23

c center of charge

22    X=0.0D0
      Y=0.0D0
      Z=0.0D0
      AT=0.0D0
      DO 32 N=1,NAT
        AM=ZAN(N)
        IF(AM.LT.0.9999D0)GOTO 32
        AT=AT+AM
        X=X+AM*C(1,N)
        Y=Y+AM*C(2,N)
        Z=Z+AM*C(3,N)
32    CONTINUE
      X=X/AT
      Y=Y/AT
      Z=Z/AT
      XA=X*BOHR
      YA=Y*BOHR
      ZA=Z*BOHR
      WRITE(9,66) X,Y,Z,XA,YA,ZA
      WRITE(10,66) X,Y,Z,XA,YA,ZA
66    FORMAT(/1X,'Current position of centre of charge is'
     & /10X,3F16.8,'  atomic units'
     & /10X,3F16.8,'  angstrom')

c to a new center

23    WRITE(9,47)
      WRITE(10,47)
47    FORMAT(//1X,'New geometry in atomic units'//)
      carw = jcar.eq.1
      if (carw) then
        call filen1 (1,file5,'.car')
        inquire (file=file5, exist=carf)
        open(20,file=file5,status='unknown')
        if (carf) call appe (20)
        write (20,200) nat
200     format (i4)
        write (20,210) title
210     format (1x,10a8)
        goto (101,102,103), itype
101       write(20,201) ttt, com_t
201       format (' ',a,a)
        goto 104
102       write(20,201) ttt, coc_t
        goto 104
103       write (20,202) ttt, x*BOHR, y*BOHR, z*BOHR
202       format (1x,a,3f15.5)
104     continue
      endif
      DO 45 N=1,NAT
        C(1,N)=C(1,N)-X
        C(2,N)=C(2,N)-Y
        C(3,N)=C(3,N)-Z
        WRITE(10,46) C(1,N),C(2,N),C(3,N)
        WRITE(9,46) C(1,N),C(2,N),C(3,N)
46      FORMAT(10X,3F16.8)
45    continue
      WRITE(10,48)
      WRITE(9,48)
48    FORMAT(//1X,'in angstrom'//)
      DO 49 N=1,NAT
        XA=C(1,N)*BOHR
        YA=C(2,N)*BOHR
        ZA=C(3,N)*BOHR
        WRITE(10,46) XA,YA,ZA
        WRITE(9,46)XA,YA,ZA
        if (carw)
     &    write(20,'(f11.5,2f10.5,i4)') xa, ya, za, int(zan(n)+.5)
49    CONTINUE
      if (carw) close(20)
      RETURN
      END
