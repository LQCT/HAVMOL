      SUBROUTINE JANDK
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL LAB,LABC,LABCD
      LOGICAL SKIPA,SKIPB,SKIPC,NPSYM
      LOGICAL OUT,IANDJ,KANDL,SAME,NORM,DTRUIJ,DTRUKL
      LOGICAL RESTAR
      LOGICAL REVISE
      LOGICAL  LADAPT,LSKIP
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      INTEGER P1,P2,P3
      DIMENSION MI(24),MJ(24),MK(24),M0(24)
      COMMON/MISC/TOL(4),CUTOFF,ICOUNT,NINTMX,OUT,IANDJ,KANDL,
     1            SAME,NORM,DTRUIJ,DTRUKL,KLUNK(5)
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),KTYPE(NSH1),
     2             KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2              NIJ,IJ,KL,IJKL
      COMMON/SYMTRY/T(216),INVT(24),ISO(NSH1+1,24),NT
     &             ,NT2,XX(21),INDEX,NAXIS
      CHARACTER*8 GROUP
      COMMON/SYMTRB/GROUP
      COMMON/BIG/AAA(900),IJDEN(100),IK(100),
     1           IJX(36),IJY(36),IJZ(36),
     2           KLX(36),KLY(36),KLZ(36),
     3           DIJ(100),DKL(100),
     &           IJGT(36),KLGT(36),
     &           IGT(1296),IX(1296),IY(1296),IZ(1296),
     4           YJUNK(160)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,INC1,INC2,INC3,INC4,INC5,
     1             LENDD,IOFF,IXI,IYI,IZI,IC6,IC7,IC8,IC9,IJUNK(6)
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMMAX,TTOT
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/SECTOR/NUM3,IBLK3,REVISE,IODA(24),ICLASS(24),
     1              ILEN(24),MAXB,KBLKLA
      CHARACTER*8 PARNAM,BFLAB,ATNAM,ANAM
      COMMON/INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(NATM1),C(3,NATM1)
     &            ,NVAR,NPAR,NVAL,NUQ
     &            ,CFR,CFA,KONTYP(NATM1),P1(NATM1),P2(NATM1),P3(NATM1)
     &            ,KR(NATM1),KA(NATM1),KB(NATM1),PARVAL(NATM13)
     &            ,ZIN(NATM1),X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAM(NATM1)
      PARAMETER (ONE=1.D0,ONEH=1.5D0,TWO=2.D0,EIGHT=8.D0,TEN=10.D0)

      write(10,'(/a/)')
     &      ' ....... Two-electron integrals ......'
      write(IW,'(/a/)')
     &      ' ....... Two-electron integrals ......'
      RESTAR=IREST.EQ.1
      CALL TIMIT(3)
      TIM0=TIM
      TIM1=TIM
      TIM2=TIM
      DTMAX=ONE
      CALL DEBUT
      IF(IST.GT.NSHELL) RETURN
      DO 400 II=IST,NSHELL
        J1=1
        DO 20 IT=1,NT
          ID=ISO(II,IT)
          IF(ID.GT.II) GO TO 390
          MI(IT)=ID
20        CONTINUE
        J0=J1
        IF(RESTAR) J0=JST
        KTPI=KTYPE(II)-1
        ICENI=KATOM(II)
        DO 380 JJ=J0,II
          JST=J1
          DO 50 IT=1,NT
            ID=MI(IT)
            JD=ISO(JJ,IT)
            MJ(IT)=JD
            IF(ID.LT.JD) THEN
              ND=ID
              ID=JD
              JD=ND
            ENDIF
            IF(ID-II) 50,40,380
40          IF(JD-JJ) 50,50,380
50          CONTINUE
          K0=J1
          IF(RESTAR) K0=KST
          KTPJ=KTYPE(JJ)-1+KTPI
          LAB=KATOM(JJ).EQ.ICENI
          DO 370 KK=K0,JJ
            KST=J1
            DO 110 IT=1,NT
              ID=MI(IT)
              JD=MJ(IT)
              KD=ISO(KK,IT)
              MK(IT)=KD
60            IF(ID.LT.JD) THEN
                ND=ID
                ID=JD
                JD=ND
              ENDIF
              IF(JD.LT.KD) THEN
                ND=JD
                JD=KD
                KD=ND
                GO TO 60
              ENDIF
              IF(ID-II) 110,90,370
90            IF(JD-JJ) 110,100,370
100           IF(KD-KK) 110,110,370
110           CONTINUE
            KTPK=KTYPE(KK)-1+KTPJ
            LABC=LAB.AND.KATOM(KK).EQ.ICENI
            L0=J1
            IF(RESTAR) L0=LST
            DO 360 LL=L0,KK
              KTPL=KTYPE(LL)-1+KTPK
              LABCD=LABC.AND.KATOM(LL).EQ.ICENI
              IF(LABCD.AND.(((KTPL/2)*2).NE.KTPL))GOTO 360
120           CONTINUE
              RESTAR=.FALSE.
              LST=J1
              N4=0
              DO 210 IT=1,NT
                ID=MI(IT)
                JD=MJ(IT)
                KD=MK(IT)
                LD=ISO(LL,IT)
130             IF(ID.LT.JD) THEN
                  ND=ID
                  ID=JD
                  JD=ND
                ENDIF
140             IF(JD.LT.KD) THEN
                  ND=JD
                  JD=KD
                  KD=ND
                  GO TO 130
                ENDIF
                IF(KD.LT.LD) THEN
                  ND=KD
                  KD=LD
                  LD=ND
                  GO TO 140
                ENDIF
                IF(ID-II) 210,170,360
170             IF(JD-JJ) 210,180,360
180             IF(KD-KK) 210,190,360
190             IF(LD-LL) 210,200,360
200             N4=N4+1
                M0(N4)=IT
210             CONTINUE
              SKIPA=JJ.EQ.KK
              SKIPB=(II.EQ.KK).OR.(JJ.EQ.LL)
              SKIPC=(II.EQ.JJ).OR.(KK.EQ.LL)
              NPSYM=.FALSE.
              IF(SKIPA.OR.SKIPB.OR.SKIPC) GO TO 300
              NPSYM=.TRUE.
              DO 260 M=1,N4
                IT=M0(M)
                IH=MI(IT)
                JH=MJ(IT)
                IF(JH.LE.IH) THEN
                  ID=IH
                  JD=JH
                ELSE
                  ID=JH
                  JD=IH
                ENDIF
                IF(.NOT.SKIPA) SKIPA=(ID.EQ.II.AND.JD.EQ.KK).OR.
     1                               (ID.EQ.JJ.AND.JD.EQ.LL)
                IF(.NOT.SKIPB) SKIPB=(ID.EQ.II.AND.JD.EQ.LL).OR.
     1                               (ID.EQ.JJ.AND.JD.EQ.KK)
                IF(SKIPA.AND.SKIPB) GO TO 270
                KH=MK(IT)
                IF(KH.LE.IH) THEN
                  ID=IH
                  KD=KH
                ELSE
                  ID=KH
                  KD=IH
                ENDIF
                IF(.NOT.SKIPC) SKIPC=(ID.EQ.II.AND.KD.EQ.LL).OR.
     1                               (ID.EQ.JJ.AND.KD.EQ.KK)
                IF(SKIPA.AND.SKIPC) GO TO 280
                IF(SKIPB.AND.SKIPC) GO TO 290
260             CONTINUE
              GO TO 300
270           SKIPC=.TRUE.
              GO TO 300
280           SKIPB=.TRUE.
              GO TO 300
290           SKIPA=.TRUE.
300           CONTINUE
              Q4=DFLOAT(NT)/DFLOAT(N4)
              IEXCH=1
              ISH=II
              JSH=JJ
              KSH=KK
              LSH=LL
              QQ4=Q4
              IF(SKIPA.AND.NPSYM) QQ4=QQ4+Q4
              IF(SKIPB.AND.NPSYM) QQ4=QQ4+Q4
              IOFF=0
              GO TO 330
310           IF(SKIPA) GO TO 320
              IEXCH=2
              ISH=II
              JSH=KK
              KSH=JJ
              LSH=LL
              QQ4=Q4
              IF(SKIPC.AND.NPSYM) QQ4=QQ4+Q4
              GO TO 330
320           IF(SKIPB.OR.SKIPC) GO TO 350
              IEXCH=3
              ISH=II
              JSH=LL
              KSH=JJ
              LSH=KK
              QQ4=Q4
330           CONTINUE
              CALL SHELLS(1,ISH,JSH,KSH,LSH)
              CALL IJPRIM
              CALL SHELLS(2,ISH,JSH,KSH,LSH)
              IF(NIJ.EQ.0)GOTO 340
              CALL GENRAL
340           CONTINUE
              GO TO (310,320,350),IEXCH
350           CONTINUE
              CALL TIMIT(3)
              DT0=TIM-TIM0
              DT2=TIM-TIM2
              TIM2=TIM
              IF(DT2.GT.DTMAX) DTMAX=DT2
              IF((TIMLIM-TIM).GT.(DTMAX*ONEH+EIGHT)) GOTO 360
              CALL FINAL(0,II,JJ,KK,LL)
360           CONTINUE
370         CONTINUE
          IF((TIMLIM-TIM).LT.(DTMAX*TWO+TEN))CALL FINAL(2,II,JJ,KK,LL)
380       CONTINUE
390     CALL TIMIT(3)
        CALL FINAL(2,II,II,II,II)
        DT0=TIM-TIM0
        DT1=TIM-TIM1
        TIM1=TIM
        IF(NPRINT.NE.-5) WRITE(IW,10) II,II,II,II,NREC,ICOUNT,DT1,TIM
        write(10 ,10) II,II,II,II,NREC,ICOUNT,DT1,TIM
400     CONTINUE
      CALL FINAL(1,II,JJ,KK,LL)
      CALL TIMIT(3)
      DTIM=TIM-TIM0
      write(10,410) DTIM,TIM
      IF(NPRINT.NE.-5) WRITE(IW,410) DTIM,TIM
10    FORMAT(' I,J,K,L =',4(I3,','),' NREC= ',I5,' ILC= ',I5
     1 ,' DEL(T) =',F8.2,' T=' ,F10.2)
410   FORMAT(/' TIME FOR INTEGRAL EVALUATION ',F10.3,' SECS.'/
     &        ' TOTAL ELAPSED TIME           ',F10.3,' SECS.'/)
      RETURN
      END

      SUBROUTINE DEBUT
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      LOGICAL  LADAPT,LSKIP,REVISE
      LOGICAL OUT,IANDJ,KANDL,SAME,NORM,DTRUIJ,DTRUKL
      real*4 xx
      INTEGER*2 IJ,KL
      COMMON/BLKIN/XX(510),IJ(510),KL(510),NWORD,NJUNK
*      COMMON/BLKIN/XX(340),IJ(340),KL(340),NWORD,NJUNK
      COMMON/INTPR/QINT(3),VALINT(3),JCINT(16)
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/SECTOR/NUM3,IBLK3,REVISE,IODA(24),ICLASS(24),
     &              ILEN(24),MAXB,KBLKLA
      COMMON/MISC/TOL(4),CUTOFF,ICOUNT,NINTMX,OUT,IANDJ,KANDL
     1          ,SAME,NORM,DTRUIJ,DTRUKL,KLUNK(5)
      PARAMETER (ONE=1.0D0,TEN=10.0D0,E=2.30258D0)

      OUT=NPRINT.EQ.4
      IF(OUT) write(10,10)
      IF(OUT) WRITE(IW,10)
10    FORMAT(//,10X,20(1H*),//,10X,'Two-electron integrals',
     1 //,10X,20(1H*))
      NORM=NORMF.NE.1.OR.NORMP.NE.1
*     NINTMX=340
      NINTMX=510
      JCINT(1)=0
      CUTOFF=ONE/(TEN**ICUT)
      TOL(1)=E*(ITOL+1)
      TOL(2)=E*(ITOL-1)
      TOL(3)=ONE/TEN**(ITOL+1)
      TOL(4)=ONE/TEN**(ITOL-1)
   35 IF(IREST.NE.1) GOTO 60
      ICOUNT=INTLOC
      IPOSM=IBLKM+(NREC-1)
      write(10,40) IFILM,IBLKM,NREC,ICOUNT,IST,JST,KST,LST,IPOSM
      WRITE(IW,40) IFILM,IBLKM,NREC,ICOUNT,IST,JST,KST,LST,IPOSM
40    FORMAT(//10X,'Integrals (re-)started'//10X,'integrals are on',
     1' file  ',I3,'  from block  ',I5
     2/10X,'NREC= ',I4,' ICOUNT= ',I4,' IST= ',I4,' JST= ',I4,' KST= ',
     3I4,' LST= ',I4
     4/10X,'restarting at block',I6)
      CALL SEARCH(IPOSM,IFILM)
      IF(INTLOC.EQ.1) GOTO 50
      CALL GET(XX,NW,IFILM)
      CALL SEARCH(IPOSM,IFILM)
50    CONTINUE
      CALL WHTPS
      RETURN
60    IST=1
      JST=1
      KST=1
      LST=1
      INTLOC=1
      NREC=1
      ICOUNT=1
      IREST=1
      CALL SEARCH(IBLKM,IFILM)
      IPOSM=1
      RETURN
      END

      SUBROUTINE FINAL(INDEX,II,JJ,KK,LL)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      LOGICAL OUT,IANDJ,KANDL,SAME,NORM,DTRUIJ,DTRUKL
      LOGICAL  LADAPT,LSKIP,REVISE
      COMMON/MISC/TOL(4),CUTOFF,ICOUNT,NINTMX,OUT,IANDJ,KANDL
     1           ,SAME,NORM,DTRUIJ,DTRUKL,KLUNK(5)
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/SECTOR/NUM3,IBLK3,REVISE,IODA(24),ICLASS(24),
     &              ILEN(24),MAXB,KBLKLA
      real*4 xx
      INTEGER*2 IOUT,JOUT
      COMMON/BLKIN/XX(510),IOUT(510),JOUT(510),NWORD,NJUNK
*      COMMON/BLKIN/XX(340),IOUT(340),JOUT(340),NWORD,NJUNK
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2              NIJ,IJ,KL,IJKL
      COMMON/INTPR/Q(3),V(3),JC,N1(3),J1(3),J2(3),J3(3),J4(3)
      DATA M0,M511/0,511/

      IBLOK=0
      IF(JC.EQ.0.OR.INDEX.EQ.2) GO TO 20
      WRITE(IW,10)(J1(M),J2(M),J3(M),J4(M),Q(M),N1(M),V(M),M=1,JC)
10    FORMAT(3(4I3,F6.3,I5,E20.12))
20    CONTINUE
      ISTR=IST
      JSTR=JST
      KSTR=KST
      LSTR=LST
      IF(INDEX.EQ.1) GO TO 60
      IST=II
      JST=JJ
      KST=KK
      LST=LL+1
      IF(LST.LE.KK) GO TO 30
      LST=1
      KST=KK+1
      IF(KST.LE.JJ) GO TO 30
      KST=1
      JST=JJ+1
      IF(JST.LE.II) GO TO 30
      JST=1
      IST=II+1
      IF(IST.GT.NSHELL) RETURN

30    INTLOC=ICOUNT
      IREST=1
      CALL SECGET(ISEX(24),24,IBLOK)
      CALL WRT3C(TITLE,LDA(ISEX(24)),IBLOK,IFILD)
      CALL SECGET(ISEX(4),4,IBLOK)
      CALL WRT3(VIBSIZ,LDA(ISEX(4)),IBLOK,IFILD)
      IF(INDEX.EQ.2)GOTO 120
      write(*,40) NREC,INTLOC,IST,JST,KST,LST
      WRITE(IW,40) NREC,INTLOC,IST,JST,KST,LST
40    FORMAT(//10X,'*** INSUFFICIENT TIME TO CONTINUE ***'//
     &' This job must be restarted to complete the evaluation of two-ele
     &ctron integrals'
     2//10X,'NREC= ',I6,' INTLOC= ',I4,' IST= ',I4,' JST= ',I4,
     3' KST= ',I4,' LST= ',I4)
      IF(INTLOC.EQ.1) GOTO 50
      IPOSM=IBLKM+(NREC-1)
      CALL SEARCH(IPOSM,IFILM)
      NWORD=INTLOC-1
      CALL PUT(XX,M511,IFILM)
50    CALL CLENUP
      write(*,'(/a/)') ' MICROMOL INTERRUPTED AT SUBR. FINAL'
      write(IW,'(/a/)') ' MICROMOL INTERRUPTED AT SUBR. FINAL'
      STOP

60    CONTINUE
      NWORD=ICOUNT-1
      IST=1
      JST=1
      KST=1
      LST=1
      IREST=2
      IPOSM=IBLKM+(NREC-1)
      CALL SEARCH(IPOSM,IFILM)
      IF(ICOUNT.EQ.1) GOTO 70
      CALL PUT(XX,M511,IFILM)
      IPOSM=IPOSM+1
70    NWORD=0
      MBLKM=IPOSM
      CALL PUT(XX,M0,IFILM)
      NINT=NINTMX*(NREC-1)+ICOUNT-1
      write(10,100) NINT
      WRITE(IW,100) NINT
      write(10 ,90) NREC,ICOUNT
      WRITE(IW,90) NREC,ICOUNT
      write(10,80) IST,JST,KST,LST
      WRITE(IW,80) IST,JST,KST,LST
80     FORMAT(' IST= ',I4,' JST= ',I4,' KST= ',I4,' LST= ',I4)
      write(10,110) IREST
      WRITE(IW,110)IREST
90     FORMAT(' There are ',I7,' record(s) of two-electron',
     1' integrals,  INTLOC=',I4)
100   FORMAT(/' Total number of two-electron integrals= ',I8)
110   FORMAT(/' End of two-electron integrals, IREST=',I2)
      INTLOC=1
      NREC=1
      ICOUNT=1
      CALL SECGET(ISEX(24),24,IBLOK)
      CALL WRT3C(TITLE,LDA(ISEX(24)),IBLOK,IFILD)
      CALL SECGET(ISEX(4),4,IBLOK)
      CALL WRT3(VIBSIZ,LDA(ISEX(4)),IBLOK,IFILD)
      CALL WHTPS
      RETURN
120   IST=ISTR
      JST=JSTR
      KST=KSTR
      LST=LSTR
      RETURN
      END

      SUBROUTINE SHELLS (NELEC,ISH,JSH,KSH,LSH)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL OUT,IANDJ,KANDL,SAME
      LOGICAL NORM,DTRUIJ,DTRUKL
      INTEGER P1,P2,P3
      CHARACTER*8 PARNAM,ATNAM,ANAM,BFLAB
      COMMON/INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(NATM1),C(3,NATM1),
     &             NVAR,NPAR,NVAL,NUQ,
     &             CFR,CFA,KONTYP(NATM1),P1(NATM1),P2(NATM1),P3(NATM1)
     &             ,KR(NATM1),KA(NATM1),KB(NATM1),PARVAL(NATM13)
     &             ,ZIN(NATM1),X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAM(NATM1)
      COMMON/RTWT/XX,U(NATM1),W(NATM1),NROOTS,JUNK,ZJUNK(NSH1*2)
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,INC1,INCL,INCK,INCJ,INCI,
     1             LENDD,IOFF,IXI,IYI,IZI,INC6,IC7,IC8,IC9,NGM
     2             ,KUNK(5)
      COMMON/BIG/AAA(900),IJDEN(100),IK(100),
     1           IJX(36),IJY(36),IJZ(36),
     2           KLX(36),KLY(36),KLZ(36),
     3           DIJ(100),DKL(100),
     &           IJGT(36),KLGT(36),
     &           IGT(1296),IX(1296),IY(1296),IZ(1296),
     4           YJUNK(160)
      COMMON/BIGG/QQ(mdm1),QJNK,QQQ(nbf1m+nbf1)
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2              NIJ,IJ,KL,IJKL
      COMMON/SHLINF/GA(10),CSA(10),CPA(10),CDA(10),
     1              GB(10),CSB(10),CPB(10),CDB(10),
     2              GC(10),CSC(10),CPC(10),CDC(10),
     3              GD(10),CSD(10),CPD(10),CDD(10),
     4              AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD,
     5              NGA,NGB,NGC,NGD,EXIJ(100),RSMALL
      COMMON/MISC/TOL(4),CUTOFF,ICOUNT,NINTMX,OUT,IANDJ,KANDL,SAME
     1           ,NORM,DTRUIJ,DTRUKL,KLUNK(5)

      IF(NELEC.NE.2) THEN
c ISH and JSH are the numbers of the shells in the system
        IANDJ=ISH.EQ.JSH
        DTRUIJ=.FALSE.
c I and J are the atoms A and B where ISH and JSH are in, respectively
        I=KATOM(ISH)
        AX=C(1,I)
        AY=C(2,I)
        AZ=C(3,I)
c I1 and J1 are the starting number of basis functions in ISH and JSH
c I2 and J2 are the ending number of basis functions in ISH and JSH
        I1=KSTART(ISH)
        I2=I1+KNG(ISH)-1
        LIT=KTYPE(ISH)
        MINI=KMIN(ISH)
        MAXI=KMAX(ISH)
        LOCI=KLOC(ISH)-MINI
        NGA=0
        DO 10 I=I1,I2
          NGA=NGA+1
          GA(NGA)=EX(I)
          CSA(NGA)=CS(I)
          CPA(NGA)=CP(I)
10        CDA(NGA)=CD(I)
        J=KATOM(JSH)
        BX=C(1,J)
        BY=C(2,J)
        BZ=C(3,J)
        J1=KSTART(JSH)
        J2=J1+KNG(JSH)-1
        LJT=KTYPE(JSH)
        MINJ=KMIN(JSH)
        MAXJ=KMAX(JSH)
        LOCJ=KLOC(JSH)-MINJ
        NGB=0
        DO 20 J=J1,J2
          NGB=NGB+1
          GB(NGB)=EX(J)
          CSB(NGB)=CS(J)
          CPB(NGB)=CP(J)
          CDB(NGB)=CD(J)
20        CONTINUE
        RAB=((AX-BX)**2+(AY-BY)**2+(AZ-BZ)**2)
        DTRUIJ=MAXI.GT.4.OR.MAXJ.GT.4
      ELSE
        KANDL=KSH.EQ.LSH
        SAME=ISH.EQ.KSH.AND.JSH.EQ.LSH
        DTRUKL=.FALSE.
        K=KATOM(KSH)
        CX=C(1,K)
        CY=C(2,K)
        CZ=C(3,K)
        K1=KSTART(KSH)
        K2=K1+KNG(KSH)-1
        LKT=KTYPE(KSH)
        MINK=KMIN(KSH)
        MAXK=KMAX(KSH)
        LOCK=KLOC(KSH)-MINK
        NGC=0
        DO 40 K=K1,K2
          NGC=NGC+1
          GC(NGC)=EX(K)
          CSC(NGC)=CS(K)
          CPC(NGC)=CP(K)
40        CDC(NGC)=CD(K)
        L=KATOM(LSH)
        DX=C(1,L)
        DY=C(2,L)
        DZ=C(3,L)
        L1=KSTART(LSH)
        L2=L1+KNG(LSH)-1
        LLT=KTYPE(LSH)
        MINL=KMIN(LSH)
        MAXL=KMAX(LSH) 
        LOCL=KLOC(LSH)-MINL
        NGD=0 
        DO 50 L=L1,L2
          NGD=NGD+1
          GD(NGD)=EX(L)
          CSD(NGD)=CS(L)
          CPD(NGD)=CP(L)
50        CDD(NGD)=CD(L)

        NROOTS=(LIT+LJT+LKT+LLT-2)/2
        RCD=((CX-DX)**2+(CY-DY)**2+(CZ-DZ)**2)
        DTRUKL=MAXK.GT.4.OR.MAXL.GT.4
        NGL=MAXL-MINL+1 
        NGK=MAXK-MINK+1
        NGJ=MAXJ-MINJ+1
        NGI=MAXI-MINI+1
        INCL=1
        INCK=NGL
        INCJ=INCK*NGK
        INCI=INCJ*NGJ
        LENDD=INCI*NGI
70      CONTINUE
        IF(MOD(LENDD,4).EQ.0)LENDD=LENDD+1
        CALL INDXB(IJGT,MINI,MAXI,MINJ,MAXJ,IANDJ,INCI,INCJ,1)
        CALL INDXB(KLGT,MINK,MAXK,MINL,MAXL,KANDL,INCK,INCL,0)
        KLN2=1
        KLN1=LLT
        IJN2=KLN1*LKT
        IJN1=IJN2*LJT
        INC1=IJN1*LIT
        IF(MOD(INC1,4).EQ.0)INC1=INC1+1
        CALL INDXA(IJX,IJY,IJZ,IJ,MINI,MAXI,MINJ,MAXJ,IANDJ,
     1             IJN1,IJN2,1)
        CALL INDXA(KLX,KLY,KLZ,KL,MINK,MAXK,MINL,MAXL,KANDL,
     1             KLN1,KLN2,0)
        MAX=KL
        DO 80 I=1,IJ
          IF(SAME) MAX=I
80        IK(I)=MAX
        IXI=LENDD
        IYI=IXI+LENDD
        IZI=IYI+LENDD
        IJKL=0
        DO 90 I=1,IJ
          MAX=IK(I)
          DO 90 K=1,MAX
            IJKL=IJKL+1
            IGT(IJKL)=IJGT(I)+KLGT(K)
            IGT(IXI+IJKL)=IJX(I)+KLX(K)
            IGT(IYI+IJKL)=IJY(I)+KLY(K)
            IGT(IZI+IJKL)=IJZ(I)+KLZ(K)
90          CONTINUE
100     INC6=1
        IC7=LENDD+1
        IC8=IC7+INC1
        IC9=IC8+INC1
        INC1=INC1*3
        DO 200 I=1,LENDD
200       QQ(I)=0.0D0
        IXI=IXI+1
        IYI=IYI+1
        IZI=IZI+1
      ENDIF
      RETURN
      END

      SUBROUTINE IJPRIM
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL OUT,IANDJ,KANDL,SAME,NORM,
     &DTRUIJ,DTRUKL
      COMMON/MISC/TOL(4),CUTOFF,ICOUNT,NINTMX,OUT,IANDJ,KANDL
     1           ,SAME,NORM,DTRUIJ,DTRUKL,KLUNK(5)
      COMMON/SHLINF/AG(10),CSA(10),CPA(10),CDA(10),
     1              BG(10),CSB(10),CPB(10),CDB(10),
     2              CG(10),CSC(10),CPC(10),CDC(10),
     3              DG(10),CSD(10),CPD(10),CDD(10),
     4              XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,
     5              NGA,NGB,NGC,NGD,EXIJ(100),RSMALL
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2              NIJ,IJ,KL,IJKL
      COMMON/BIG/A(100),R(100),
     &           X1(100),Y1(100),Z1(100),DIJ(400),
     &           IJDEN(100),IDUM0(100),
     1           IDUM10(6,36),ADUM1(200),IDUM11(2,36),
     2           IDUM2(5184),ADUM2(160)
      PARAMETER (ONE=1.0D0)

      MAX=MAXJ
      N=0
      NN=0
      DO 50 I=MINI,MAXI
        GO TO (10,10,20,20,10,20,20,20,20,20),I
10        NM=NN
20        NN=NM
        IF(IANDJ) MAX=I
        DO 50 J=MINJ,MAX
          GO TO (30,30,40,40,30,40,40,40,40,40),J
30          NN=NN+1
40          N=N+1
50        IJDEN(N)=NN
      NIJ=0
      JBMAX=NGB
      DO 190 IA=1,NGA
        AI=AG(IA)
        ARRI=AI*RRI
        AXI=AI*XI
        AYI=AI*YI
        AZI=AI*ZI
        CSI=CSA(IA)
        CPI=CPA(IA)
        CDI=CDA(IA)
        IF(IANDJ) JBMAX=IA
        DO 180 JB=1,JBMAX
          AJ=BG(JB)
          AA=AI+AJ
          AINV=ONE/AA
          DUM=AJ*ARRI*AINV
          CSJ=CSB(JB)*AINV
          CPJ=CPB(JB)*AINV
          CDJ=CDB(JB)*AINV
          NM=(NIJ+NIJ)+(NIJ+NIJ)
          NN=NM
          NIJ=NIJ+1
          EXIJ(NIJ)=DEXP(-DUM)
          R(NIJ)=DUM
          A(NIJ)=AA
          X1(NIJ)=(AXI+AJ*XJ)*AINV
          Y1(NIJ)=(AYI+AJ*YJ)*AINV
          Z1(NIJ)=(AZI+AJ*ZJ)*AINV
          DO 150 I=MINI,MAXI
            IF(IANDJ)MAX=I
            GO TO(60,70,150,150,80,150,150,150,150,150),I
60            DUM1=CSI
            GO TO 90
70            DUM1=CPI
            GOTO 90
80            DUM1=CDI
90          CONTINUE
            DO 140 J=MINJ,MAX
              GO TO(100,110,140,140,120,140,140,140,140,140),J
100             DUM2=DUM1*CSJ
              GO TO 130
110             DUM2=DUM1*CPJ
              GOTO 130
120             DUM2=DUM1*CDJ
130           NN=NN+1
              DIJ(NN)=DUM2
140           CONTINUE
150         CONTINUE
          IF(.NOT.IANDJ) GO TO 180
          IF(IA.EQ.JB) GO TO 180
          GO TO (170,160,170,170),LIT
160       IF(MINI.EQ.2) GO TO 170
          DIJ(NM+2)=DIJ(NM+2)+CSI*CPJ
          DIJ(NM+3)=DIJ(NM+3)+DIJ(NM+3)
170       DIJ(NM+1)=DIJ(NM+1)+DIJ(NM+1)
180       CONTINUE
190     CONTINUE
      IF(NIJ.EQ.0)GOTO 210
      RSMALL=R(1)
      DO 200 N=1,NIJ
        IF(RSMALL.GT.R(N))RSMALL=R(N)
200     CONTINUE
      IF(RSMALL.GE.TOL(1))NIJ=0
210   CONTINUE
      RETURN
      END

      SUBROUTINE GENRAL
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL OUT,IANDJ,KANDL,SAME,DOUBLE,NORM
      LOGICAL DTRUIJ,DTRUKL,DTRU
      COMMON/SHLINF/AG(10),CSA(10),CPA(10),CDA(10),
     1              BG(10),CSB(10),CPB(10),CDB(10),
     2              CG(10),CSC(10),CPC(10),CDC(10),
     3              DG(10),CSD(10),CPD(10),CDD(10),
     4              XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,
     5              NGA,NGB,NGC,NGD,EXIJ(100),RSMALL
      COMMON/BIG/AA(100),R(100),
     &           X1(100),Y1(100),Z1(100),
     &           DD(400),IJDEN(100),IK(100),
     1           IJX(36),IJY(36),IJZ(36),
     2           KLX(36),KLY(36),KLZ(36),
     3           DIJ(100),DKL(100),
     &           IJGT(36),KLGT(36),
     &           IGT(1296),IX(1296),IY(1296),IZ(1296),
     &           YJUNK(160)
      COMMON/BIGG/QQ(mdm1),QJNK,QQQ(nbf1m+nbf1)
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2              NIJ,IJ,KL,IJKL
      COMMON/RTWT/XX,U(NATM1),W(NATM1),NROOTS,JUNK,ZJUNK(NSH1*2)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,INC1,INC2,INC3,INC4,INC5,
     1             LENDD,IOFF,IXI,IYI,IZI,IC6,IC7,IC8,IC9,IJUNK(6)
      COMMON/MISC/TOL(4),CUTOFF,ICOUNT,NINTMX,OUT,IANDJ,KANDL
     1           ,SAME,NORM,DTRUIJ,DTRUKL,KLUNK(5)
      COMMON/SETINT/BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00
     1             ,F00,DXIJ,DYIJ,DZIJ,DXKL,DYKL,DZKL,
     2             IN(NATM1),KN(NATM1)
     3             ,NI,NJ,NK,NL,NMAX,MMAX,IJ1,IJ2,KL1,KL2
      DIMENSION IN1(NATM1),IQQ(1)
      EQUIVALENCE (IQQ,IGT)
      PARAMETER (PI252=34.986836655250D0, ONE=1.0D0, PT5=0.5D0)

      DTRU=DTRUIJ.OR.DTRUKL
      FACTOR=PI252*QQ4
      IF (LIT .GE. LJT) THEN
        NI = LIT-1
        NJ = LJT-1
        IJ1 = IJN1
        IJ2 = IJN2  
        XC = XI 
        YC = YI
        ZC = ZI
        DXIJ = XI-XJ
        DYIJ = YI-YJ
        DZIJ = ZI-ZJ
      ELSE
        NI = LJT-1
        NJ = LIT-1
        IJ1 = IJN2
        IJ2 = IJN1
        XC = XJ
        YC = YJ
        ZC = ZJ
        DXIJ = XJ-XI
        DYIJ = YJ-YI
        DZIJ = ZJ-ZI
      ENDIF
      IF (LKT .GE. LLT) THEN
        NK = LKT-1
        NL = LLT-1
        KL1 = KLN1
        KL2 = KLN2
        XD = XK
        YD = YK
        ZD = ZK
        DXKL = XK-XL
        DYKL = YK-YL
        DZKL = ZK-ZL
      ELSE
        NK = LLT-1
        NL = LKT-1
        KL1 = KLN2
        KL2 = KLN1
        XD = XL
        YD = YL
        ZD = ZL
        DXKL = XL-XK
        DYKL = YL-YK
        DZKL = ZL-ZK
      ENDIF
      NMAX = NI+NJ
      MMAX = NK+NL
      MAX = NMAX+1
      DO 50 I = 1,MAX
        N = I-1
        IF (N .LE. NI) IN1(I) = IJ1*N+1
        IF (N .GT. NI) IN1(I) = IJ1*NI+IJ2*(N-NI)+1
50      CONTINUE
      MAX = MMAX+1
      DO 60 K = 1,MAX
        N = K-1
        IF (N .LE. NK) KN(K) = KL1*N
        IF (N .GT. NK) KN(K) = KL1*NK+KL2*(N-NK)
60      CONTINUE
      LGMAX=NGD
      DO 150 KG=1,NGC
        AK=CG(KG)
        BRRK=AK*RRK
        AKXK=AK*XK
        AKYK=AK*YK
        AKZK=AK*ZK
        CSK=CSC(KG)*FACTOR
        CPK=CPC(KG)*FACTOR
        CDK=CDC(KG)*FACTOR
        IF(KANDL) LGMAX=KG
        DO 140 LG=1,LGMAX
          AL=DG(LG)
          B=AK+AL
          BINV=ONE/B
          BBRRK=AL*BRRK*BINV
          IF((BBRRK+RSMALL).GT.TOL(1)) GO TO 140
          EXKL=DEXP(-BBRRK)
          CSL=CSD(LG)*BINV
          CPL=CPD(LG)*BINV
          CDL=CDD(LG)*BINV
          XB=(AKXK+AL*XL)*BINV
          YB=(AKYK+AL*YL)*BINV
          ZB=(AKZK+AL*ZL)*BINV
          BXBK=B*(XB-XD)
          BYBK=B*(YB-YD)
          BZBK=B*(ZB-ZD)
          BXBI=B*(XB-XC)
          BYBI=B*(YB-YC)
          BZBI=B*(ZB-ZC)
          DOUBLE=KANDL.AND.KG.NE.LG
          CALL DENFAC(DKL,CSK,CPK,CDK,CSL,CPL,CDL,MINK,MAXK,
     1                MINL,MAXL,KANDL,DOUBLE)
          DKLD=DKL(1)
          NN=0
          DKLD=DKLD*EXKL 
          IF(DABS(DKLD).LT.TOL(3))GOTO 140
          DO 130 N=1,NIJ
            DUM=BBRRK+R(N)
            IF(DUM.GT.TOL(2)) GO TO 130
            DIJD=DD(NN+1)
            DIJD=DIJD*DKLD
            A=AA(N)
            AB=A*B
            AANDB=A+B
            EXPE=DIJD*EXIJ(N)/DSQRT(AANDB)
            IF(DABS(EXPE).LT.TOL(4))GOTO 130
            RHO=AB/AANDB
            XA=X1(N)
            YA=Y1(N)
            ZA=Z1(N)
            XX=RHO*((XA-XB)**2+(YA-YB)**2+(ZA-ZB)**2)
            AXAK=A*(XA-XD)
            AYAK=A*(YA-YD)
            AZAK=A*(ZA-ZD)
            AXAI=A*(XA-XC)
            AYAI=A*(YA-YC)
            AZAI=A*(ZA-ZC)
            C1X=BXBK+AXAK
            C2X=A*BXBK
            C3X=BXBI+AXAI
            C4X=B*AXAI
            C1Y=BYBK+AYAK
            C2Y=A*BYBK
            C3Y=BYBI+AYAI
            C4Y=B*AYAI
            C1Z=BZBK+AZAK
            C2Z=A*BZBK
            C3Z=BZBI+AZAI
            C4Z=B*AZAI
            CALL ROOTS
            MM=0
            MAX=NMAX+1
            DO 120 M=1,NROOTS
              U2=U(M)*RHO
              F00=EXPE*W(M)
              DO 110 I=1,MAX
110             IN(I)=IN1(I)+MM
              DUM=ONE/(AB+U2*AANDB)
              DUM2=PT5*DUM
              BP01=(A+U2)*DUM2
              B00=U2*DUM2
              B10=(B+U2)*DUM2
              XCP00=(U2*C1X+C2X)*DUM
              XC00 =(U2*C3X+C4X)*DUM
              YCP00=(U2*C1Y+C2Y)*DUM
              YC00 =(U2*C3Y+C4Y)*DUM
              ZCP00=(U2*C1Z+C2Z)*DUM
              ZC00 =(U2*C3Z+C4Z)*DUM
              CALL XYZE(QQ(IC7),QQ(IC8),QQ(IC9))
120           MM=MM+INC1
            CALL FMAE(QQ(IC6),QQ(IC7),QQ(IC8),
     &                QQ(IC9),IQQ(1),IQQ(IXI),IQQ(IYI),IQQ(IZI))
130         NN=NN+4
140       CONTINUE
150     CONTINUE
      IF(DTRU.AND.NORM)CALL INORM(QQ(IC6))
      CALL QOUT(QQ(IC6))
250   CONTINUE

      RETURN
      END

C  INDXB  A1
C
      SUBROUTINE INDXB(IJGT,MINI,MAXI,
     &                      MINJ,MAXJ,IANDJ,INC1,INC2,INC3)
      LOGICAL IANDJ
      DIMENSION IJGT(*)

      IJ1=0
      JMAX=MAXJ
      DO 1 I=MINI,MAXI
        IF(IANDJ) JMAX=I
        DO 1 J=MINJ,JMAX
          IJ1=IJ1+1
          IJGT(IJ1)=INC1*(I-MINI)+INC2*(J-MINJ)+INC3
1         CONTINUE
      RETURN
      END

      SUBROUTINE XYZE(X,Y,Z)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL N0,N1,M0,M1
      COMMON/SETINT/A(16),I(NATM1),K(NATM1),NIMAX,
     3              NJMAX,NKMAX,NLMAX,NMAX,MMAX,IJKL(4)
      DIMENSION X(*),Y(*),Z(*)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)

      BP01=A(1)
      B00=A(2)
      B10=A(3)
      XCP00=A(4)
      XC00=A(5)
      YCP00=A(6)
      YC00=A(7)
      ZCP00=A(8)
      ZC00=A(9)
      F00=A(10)
      DXIJ=A(11)
      DYIJ=A(12)
      DZIJ=A(13)
      DXKL=A(14)
      DYKL=A(15)
      DZKL=A(16)
      IJ1=IJKL(1)
      IJ2=IJKL(2)
      KL1=IJKL(3)
      KL2=IJKL(4)
      N0=NMAX.EQ.0
      N1=NMAX.LE.1
      M0=MMAX.EQ.0
      M1=MMAX.LE.1
      IF(N0)GOTO 320
      IF(M0)GOTO 250
      I1=I(1)
      X(I1)=ONE
      Y(I1)=ONE
      Z(I1)=F00
      I2=I(2)
      X(I2)=XC00
      Y(I2)=YC00
      Z(I2)=ZC00*F00
      K2=K(2)
      I3=I1+K2
      X(I3)=XCP00
      Y(I3)=YCP00
      Z(I3)=ZCP00*F00
      I3=I2+K2
      CP10=B00
      X(I3)=XCP00*X(I2)+CP10
      Y(I3)=YCP00*Y(I2)+CP10
      Z(I3)=ZCP00*Z(I2)+CP10*F00
      IF(N1) GO TO 30
      C10=ZERO
      I3=I1
      I4=I2
      DO 20 N=2,NMAX
        C10=C10+B10
        I5=I(N+1)
        X(I5)=C10*X(I3)+XC00*X(I4)
        Y(I5)=C10*Y(I3)+YC00*Y(I4)
        Z(I5)=C10*Z(I3)+ZC00*Z(I4)
        CP10=CP10+B00
        I3=I5+K2
        X(I3)=XCP00*X(I5)+CP10*X(I4)
        Y(I3)=YCP00*Y(I5)+CP10*Y(I4)
        Z(I3)=ZCP00*Z(I5)+CP10*Z(I4)
10      I3=I4
20      I4=I5
30    IF(M1) GO TO 60
      CP01=ZERO
      C01=B00
      I3=I1
      I4=I1+K2
      DO 50 M=2,MMAX
        CP01=CP01+BP01
        I5=I1+K(M+1)
        X(I5)=CP01*X(I3)+XCP00*X(I4)
        Y(I5)=CP01*Y(I3)+YCP00*Y(I4)
        Z(I5)=CP01*Z(I3)+ZCP00*Z(I4)
        C01=C01+B00
        I3=I2+K(M+1)
        X(I3)=XC00*X(I5)+C01*X(I4)
        Y(I3)=YC00*Y(I5)+C01*Y(I4)
        Z(I3)=ZC00*Z(I5)+C01*Z(I4)
40      I3=I4
50      I4=I5
60    IF(N1.OR.M1) GO TO 90
      C01=B00
      K3=K2
      DO 80 M=2,MMAX
        K4=K(M+1)
        C01=C01+B00
        I3=I1
        I4=I2
        C10=B10
        DO 70 N=2,NMAX
          I5=I(N+1)
          I5K4=I5+K4
          I3K4=I3+K4
          I4K4=I4+K4
          I4K3=I4+K3
          X(I5K4)=C10*X(I3K4)+XC00*X(I4K4)+C01*X(I4K3)
          Y(I5K4)=C10*Y(I3K4)+YC00*Y(I4K4)+C01*Y(I4K3)
          Z(I5K4)=C10*Z(I3K4)+ZC00*Z(I4K4)+C01*Z(I4K3)
          C10=C10+B10
          I3=I4
70        I4=I5
80      K3=K4
90    IF(NJMAX.EQ.0) GO TO 160
      M=0
      I5=I(NMAX+1)
100   MIN=NIMAX
      KM=K(M+1)
110   N=NMAX
      I3=I5+KM
120   I4=I(N)+KM
      X(I3)=X(I3)+DXIJ*X(I4)
      Y(I3)=Y(I3)+DYIJ*Y(I4)
      Z(I3)=Z(I3)+DZIJ*Z(I4)
      I3=I4
      N=N-1
      IF(N.GT.MIN) GO TO 120
      MIN=MIN+1
      IF(MIN.LT.NMAX) GO TO 110
      IF(NIMAX.EQ.0) GO TO 150
      I3=IJ2+KM+I1
      DO 140 NJ=1,NJMAX
        I4=I3
        DO 130 NI=1,NIMAX
          I4M=I4-IJ2
          I4P=I4M+IJ1
          X(I4)=X(I4P)+DXIJ*X(I4M)
          Y(I4)=Y(I4P)+DYIJ*Y(I4M)
          Z(I4)=Z(I4P)+DZIJ*Z(I4M)
130       I4=I4+IJ1
140     I3=I3+IJ2
150   M=M+1
      IF(M.LE.MMAX) GO TO 100
160   IF(NLMAX.EQ.0) GO TO 240
      I5=K(MMAX+1)
      IA=I1
      NI=0
170   NJ=0
      IB=IA
180   MIN=NKMAX
190   M=MMAX
      I3=IB+I5
200   I4=IB+K(M)
      X(I3)=X(I3)+DXKL*X(I4)
      Y(I3)=Y(I3)+DYKL*Y(I4)
      Z(I3)=Z(I3)+DZKL*Z(I4)
      I3=I4
      M=M-1
      IF(M.GT.MIN) GO TO 200
      MIN=MIN+1
      IF(MIN.LT.MMAX) GO TO 190
      IF(NKMAX.EQ.0) GO TO 230
      I3=IB+KL2
      DO 220 NL=1,NLMAX
        I4=I3
        DO 210 NK=1,NKMAX
          I4M=I4-KL2
          I4P=I4M+KL1
          X(I4)=X(I4P)+DXKL*X(I4M)
          Y(I4)=Y(I4P)+DYKL*Y(I4M)
          Z(I4)=Z(I4P)+DZKL*Z(I4M)
210       I4=I4+KL1
220     I3=I3+KL2
230   NJ=NJ+1
      IB=IB+IJ2
      IF(NJ.LE.NJMAX) GO TO 180
      NI=NI+1
      IA=IA+IJ1
      IF(NI.LE.NIMAX) GO TO 170
240   RETURN
250   I1 = I(1)
      X(I1) = ONE
      Y(I1) = ONE
      Z(I1) = F00
      IF (N0) GO TO 310
      I2 = I(2)
      X(I2) = XC00
      Y(I2) = YC00
      Z(I2) = ZC00*F00
      IF (N1) GO TO 310
      C10 = ZERO
      I3 = I1
      I4 = I2
      DO 260 NI = 2,NMAX
        C10 = C10+B10
        I5 = I(NI+1)
        X(I5) = C10*X(I3)+XC00*X(I4)
        Y(I5) = C10*Y(I3)+YC00*Y(I4)
        Z(I5) = C10*Z(I3)+ZC00*Z(I4)
        I3 = I4
260     I4 = I5
      IF (NJMAX .EQ. 0) GO TO 310
      I5 = I(NMAX+1)
      MIN = NIMAX
270   NI = NMAX
      I3 = I5
280   I4 = I(NI)
      X(I3) = X(I3)+DXIJ*X(I4)
      Y(I3) = Y(I3)+DYIJ*Y(I4)
      Z(I3) = Z(I3)+DZIJ*Z(I4)
      I3 = I4
      NI = NI-1
      IF (NI .GT. MIN) GO TO 280
      MIN = MIN+1
      IF (MIN .LT. NMAX) GO TO 270
      IF (NIMAX .EQ. 0) GO TO 310
      I3 = IJ2+I1
      DO 300 NJ = 1,NJMAX
        I4 = I3
        DO 290 NI = 1,NIMAX
          I4M=I4-IJ2
          I4P=I4M+IJ1
          X(I4) = X(I4P)+DXIJ*X(I4M)
          Y(I4) = Y(I4P)+DYIJ*Y(I4M)
          Z(I4) = Z(I4P)+DZIJ*Z(I4M)
290       I4 = I4+IJ1
300     I3 = I3+IJ2
310   RETURN

320   I1 = I(1)
      X(I1) = ONE
      Y(I1) = ONE
      Z(I1) = F00
      IF (M0) GO TO 380
      K2 = K(2)
      I3 = I1+K2
      X(I3) = XCP00
      Y(I3) = YCP00
      Z(I3) = ZCP00*F00
      IF (M1) GO TO 380
      CP01 = ZERO
      I3 = I1
      I4 = I1+K2
      DO 330 NK = 2,MMAX
        CP01 = CP01+BP01
        I5 = I1+K(NK+1)
        X(I5) = CP01*X(I3)+XCP00*X(I4)
        Y(I5) = CP01*Y(I3)+YCP00*Y(I4)
        Z(I5) = CP01*Z(I3)+ZCP00*Z(I4)
        I3 = I4
330     I4 = I5
      IF (NLMAX .EQ. 0) GO TO 380
      I5 = K(MMAX+1)
      MIN=NKMAX
340   NK = MMAX
      I3 = I1+I5
350   I4 = I1+K(NK)
      X(I3) = X(I3)+DXKL*X(I4)
      Y(I3) = Y(I3)+DYKL*Y(I4)
      Z(I3) = Z(I3)+DZKL*Z(I4)
      I3 = I4
      NK = NK-1
      IF (NK .GT. MIN) GO TO 350
      MIN = MIN+1
      IF (MIN .LT. MMAX) GO TO 340
      IF (NKMAX .EQ. 0) GO TO 380
      I3 = I1+KL2
      DO 370 NL = 1,NLMAX
        I4 = I3
        DO 360 NK = 1,NKMAX
          I4M=I4-KL2
          I4P=I4M+KL1
          X(I4) = X(I4P)+DXKL*X(I4M)
          Y(I4) = Y(I4P)+DYKL*Y(I4M)
          Z(I4) = Z(I4P)+DZKL*Z(I4M)
360       I4 = I4+KL1
370     I3 = I3+KL2

380   RETURN
      END

      SUBROUTINE INORM(ABDENS)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      DIMENSION ABDENS(*)
      LOGICAL OUT,IANDJ,KANDL,SAME,NORM,DTRUIJ,DTRUKL
      COMMON/MISC/TOL(4),CUTOFF,ICOUNT,NINTMX,OUT,IANDJ,KANDL
     &           ,SAME,NORM,DTRUIJ,DTRUKL,KLUNK(5)
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     &              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     &              NIJ,IJ,KL,IJKL
      COMMON/BIG/AAA(900),IJDEN(100),IK(100),
     1           IJX(36),IJY(36),IJZ(36),
     2           KLX(36),KLY(36),KLZ(36),
     3           DIJ(100),DKL(100),
     &           IJGT(36),KLGT(36),
     &           YJUNK(2752)
      PARAMETER (ONE=1.0D0, ROOT3=1.7320508D0, ROOT5=2.236068D0)

      N=0
      MAX=MAXJ
      DUM1=ONE
      DO 20 I=MINI,MAXI
        IF(I.EQ.8) DUM1=ROOT3
        IF(I.EQ.14)DUM1=ROOT5
        IF(I.EQ.20)DUM1=DUM1*ROOT3
        DUM2=DUM1
        IF(IANDJ) MAX=I
        DO 10 J=MINJ,MAX
          IF(J.EQ.8) DUM2=DUM2*ROOT3
          IF(J.EQ.14)DUM2=DUM2*ROOT5
          IF(J.EQ.20)DUM2=DUM2*ROOT3
          N=N+1
  10      DIJ(N)=DUM2
  20    CONTINUE
      N=0
      DUM1=ONE
      MAX=MAXL
      DO 40 K=MINK,MAXK
        IF(K.EQ.8) DUM1=ROOT3
        IF(K.EQ.14)DUM1=ROOT5
        IF(K.EQ.20)DUM1=DUM1*ROOT3
        DUM2=DUM1
        IF(KANDL) MAX=K
        DO 30 L=MINL,MAX
          IF(L.EQ.8) DUM2=DUM2*ROOT3
          IF(L.EQ.14)DUM2=DUM2*ROOT5
          IF(L.EQ.20)DUM2=DUM2*ROOT3
          N=N+1
  30      DKL(N)=DUM2
  40    CONTINUE

      DO 50 I=1,IJ
        D1=DIJ(I)
        N1=IJGT(I)
        MAX=IK(I)
        DO 50 K=1,MAX
          N=N1+KLGT(K)
          ABDENS(N)=ABDENS(N)*D1*DKL(K)
  50      CONTINUE

      RETURN
      END

      SUBROUTINE FMAE(G,X,Y,Z,IGT,IX,IY,IZ)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      COMMON/RTWT/XX,U(NATM1),W(NATM1),NROOTS,JUNK,ZJUNK(NSH1*2)
      COMMON/SHLNOS/QQ4,LIT(17),IJ,KL,IJKL
      COMMON/BIG/AAA(900),IJDEN(100),IK(100),
     1           IJX(36),IJY(36),IJZ(36),
     2           KLX(36),KLY(36),KLZ(36),
     3           DIJ(100),DKL(100),
     &           IJGT(36),KLGT(36),
     &           YJUNK(2752)
      COMMON/INCRS/KLN2,KLN1,IJN2,IJN1,INC,JJUNK(19)
      DIMENSION G(*),X(*),Y(*),Z(*),IGT(*),IX(*),IY(*),IZ(*)

      I1=INC
      DO 120 I=1,IJKL
        N=IGT(I)
        MX=IX(I)
        MY=IY(I)
        MZ=IZ(I)
        S=0.0
        DO 130 NR=1,NROOTS
          S=S+X(MX)*Y(MY)*Z(MZ)
          MX=I1+MX
          MY=I1+MY
130       MZ=I1+MZ
        G(N)=G(N)+S
120     CONTINUE

      RETURN
      END

      SUBROUTINE QOUT(G)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      LOGICAL IANDJ,KANDL,SAME,OUT,LJUNK1,LJUNK2
     &,LJUNK3,DTRUIJ,DTRUKL
      LOGICAL  LADAPT,LSKIP,REVISE
      DIMENSION G(*)
      COMMON/MISC/TOL(4),CUTOFF,ICOUNT,NINTMX,OUT,IANDJ,KANDL,SAME
     &           ,LJUNK1,DTRUIJ,DTRUKL,LJUNK2,LJUNK3,KLUNK(3)
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/SECTOR/NUM3,IBLK3,REVISE,IODA(24),ICLASS(24),
     1              ILEN(24),MAXB,KBLKLA
      COMMON/SHLNOS/QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     1              MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     2              NIJ,IJ,KL,IJKL
      COMMON/BIG/AAA(900),IJDEN(100),IK(100),
     1           IJX(36),IJY(36),IJZ(36),
     2           KLX(36),KLY(36),KLZ(36),
     3           DIJ(100),DKL(100),
     &           IJGT(36),KLGT(36),
     &           YJUNK(2752)
      real*4 xx
      INTEGER*2 IOUT,JOUT
      COMMON/BLKIN/XX(510),IOUT(510),JOUT(510),NWORD,NJUNK
*      COMMON/BLKIN/XX(340),IOUT(340),JOUT(340),NWORD,NJUNK

      IJN=0
      JMAX=MAXJ
      DO 70 I=MINI,MAXI
        IF(IANDJ) JMAX=I
        I1=LOCI+I
        DO 60 J=MINJ,JMAX
          IJN=IJN+1
          N1=IJGT(IJN)
          I2=LOCJ+J
          LMAX=MAXL
          KLN=0
          DO 50 K=MINK,MAXK
            IF(KANDL) LMAX=K
            I3=LOCK+K
            DO 40 L=MINL,LMAX
              KLN=KLN+1
              IF(SAME.AND.KLN.GT.IJN) GO TO 60
              NN=N1+KLGT(KLN)
              VAL=G(NN)
              IF(DABS(VAL).LT.CUTOFF) GO TO 40
              I4=LOCL+L
              LAB1=MAX0(I1,I2)*256+MIN0(I1,I2)
              LAB2=MAX0(I3,I4)*256+MIN0(I3,I4)
              IF (OUT) CALL INTOUT(I1,I2,I3,I4,QQ4,NN,VAL)
              XX(ICOUNT)=VAL
              IOUT(ICOUNT)=MAX0(LAB1,LAB2)
              JOUT(ICOUNT)=MIN0(LAB1,LAB2)
              ICOUNT=ICOUNT+1
              IF(ICOUNT.LE.NINTMX) GO TO 40
              NWORD=NINTMX
              CALL PUT(XX,511,IFILM)
              ICOUNT=1
              NREC=NREC+1
40            CONTINUE
50          CONTINUE
60        CONTINUE
70      CONTINUE

      RETURN
      END

      SUBROUTINE INTOUT(I1,I2,I3,I4,Q4,NN,VAL)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/INTPR/Q(3),V(3),JC,N1(3),J1(3),J2(3),J3(3),J4(3)

      JC=JC+1
      J1(JC)=I1
      J2(JC)=I2
      J3(JC)=I3
      J4(JC)=I4
      Q(JC)=Q4
      N1(JC)=NN
      V(JC)=VAL
      IF(JC.LT.3) GO TO 100
      JC=0
      WRITE(9,500)(J1(M),J2(M),J3(M),J4(M),Q(M),N1(M),V(M),M=1,3)
  500 FORMAT(3(4I3,F6.3,I5,E20.12))
  100 CONTINUE
      RETURN
      END
