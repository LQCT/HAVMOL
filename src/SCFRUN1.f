C  QMAT  A1

      SUBROUTINE QMAT(S,V,E,work,NDIM)

C CREATES THE CANONICAL MOLECULAR ORBITAL SET FROM OVERLAP MATRIX BY
C DIAGONALIZATION AND FURTHER SELECTION OF SIGNIFICANT EIGENVALUES
* S is dimensioned as a one dimensional array with the lower triangle
*   and diagonal terms of the overlap matrix
* V is a two dimensional array dimensioned to allow quite all the
*   maximum core available
* E is a linear array of the maximum number of basis functions allowed
*   dimensioned from V to the tail of the maximum core available

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL  LADAPT,LSKIP
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
      COMMON/MAPPER/IA(NBF1+1)
      CHARACTER*8 BFLAB,PARNAM,ATNAM,ANAM
      COMMON/INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(NATM1),C(3,NATM1),
     &             NVAR,NPAR,NVAL,NUQ,
     &             CFR,CFA,KONTYP(NATM1),P1(NATM1),P2(NATM1),P3(NATM1)
     &             ,KR(NATM1),KA(NATM1),KB(NATM1),PARVAL(NATM13)
     &             ,ZIN(NATM1),X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAM(NATM1)
      DIMENSION S(*),E(*),V(NDIM,*),work(*)
*      PARAMETER (ZERO=0.D0,ONE=1.D0,TOL=1.0D-05)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TOL=1.0D-06)

C     ----- READ IN OVERLAP MATRIX AND DIAGONALIZE IT -----

      IBLOK=0
      ITYPE=5
      if (ncoorb.eq.0) ncoorb=num
      CALL SECGET(ISEX(5),ITYPE,IBLOK)
      CALL READQ(S,NX,IBLOK,IFILD)
      if (nprint.eq.2) then
        DO 20 I=1,NUM
20        E(I)=ZERO
        DO 50 I=1,NUM
          DO 50 J=1,NUM
            IF(J.GT.I) THEN
              IJ=IA(J)+I
            ELSE
              IJ=IA(I)+J
            ENDIF
            V(I,J)=S(IJ)
50          CONTINUE
        WRITE(IW,60)
60      FORMAT(//,5X,14(1H-),/,5X,'OVERLAP MATRIX',/,5X,14(1H-))
        CALL VECOUT(V,E,NUM,NUM,NDIM,IW,BFLAB)
      endif

*      call hqrii1 (NUM,1,NUM,1,S,E,NDIM,V,.false.,work)
      CALL LIGEN(S,V,E,IA,NUM,NDIM,2,1.0D-16)
      if (nprint.eq.2) then
        WRITE(IW,80)
80      FORMAT(//,5X,17(1H-),/,5X,'EIGENVECTORS OF S',/,5X,17(1H-))
        CALL VECOUT(V,E,NUM,NUM,NDIM,IW,BFLAB)
      endif

C     ----- ELIMINATE EIGENVECTORS FOR WHICH EIGENVALUE IS LESS
C           THAN TOL=1.0D-06

      DUM=E(NUM)
      J=0
      K=0
      DO 120 I=1,NUM
        IF(E(I).LT.TOL) GO TO 110
        J=J+1
        E(J)=ONE/DSQRT(E(I))
        DO 100 N=1,NUM
100       V(N,J)=V(N,I)
        GO TO 120
110     K=K+1
120     CONTINUE
      NCOORB=ncoorb-K
      IF(K.NE.0) WRITE(IW,130) DUM,K,TOL,NCOORB
130   FORMAT(//' ..... WARNING .....',//,
     1 ' THE SMALLEST EIGENVALUE OF THE OVERLAP MATRIX IS ',F15.8,/,
     2 ' THERE IS(ARE) ',I5,' EIGENVALUE(S) LESS THAN ',D10.2,/,
     3 ' THE NUMBER OF CANONICAL ORBITALS KEPT IS ',I5)

C     ----- FORM CANONICAL ORTHONORMAL ORBITALS AND STORE IN -----
C     ----- SCRATCH FILE                                     -----

      DO 140 J=1,NCOORB
        DO 140 I=1,NUM
140       V(I,J)=V(I,J)*E(J)
      CALL WRT3(V,NUM*NCOORB,IBLKS,IFILS)

      ITYPE=4
      CALL SECGET(ISEX(4),ITYPE,IBLOK)
      CALL WRT3(VIBSIZ,LDA(ISEX(4)),IBLOK,IFILD)
      CALL SECGET(ISEX(24),24,IBLOK)
      CALL WRT3C(TITLE,LDA(ISEX(24)),IBLOK,IFILD)
      if (nprint.eq.2) then
        WRITE(IW,150)
150     FORMAT(//5X,30(1H-),/,5X,'CANONICAL ORTHONORMAL ORBITALS',/,
     1 5X,30(1H-))
        CALL VECOUT(V,E,NUM,NCOORB,NDIM,IW,BFLAB)
      endif
      RETURN
      END

C  HCORE  A1

      SUBROUTINE HCORE(H,V,E,COORB,X,H0,NDIM)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL OUT
      INTEGER P1,P2,P3
      LOGICAL  LADAPT,LSKIP,iandj
      DIMENSION H(*),H0(*),E(*),X(*),V(NDIM,*),COORB(NDIM,*)
      COMMON/MAPPER/IA(NBF1+1)

      CHARACTER*8 PARNAM,BFLAB,ATNAM,ANAM
      COMMON/INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(NATM1),C(3,NATM1),
     &             NVAR,NPAR,NVAL,NUQ,
     &             CFR,CFA,KONTYP(NATM1),P1(NATM1),P2(NATM1),P3(NATM1)
     &             ,KR(NATM1),KA(NATM1),KB(NATM1),PARVAL(NATM13)
     &             ,ZIN(NATM1),XTRA(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAM(NATM1)
      CHARACTER*8 TITLE,SCFTYP,RUNTYP,GUESS,CONF
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      real*4 nopt,minull
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),
     2             KTYPE(NSH1),KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2
      common/huckel/huckel(nsh1),elems(3,18),elemp(2,18)
      CHARACTER*8 ALTER,BLANK,CLOSED,OSCF
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.0D0)
      DATA ALTER/'ALTER'/, BLANK /'        '/
      DATA CLOSED/'CLOSED'/, OSCF/'OSCF'/

      write(10,1101) NB,SCFTYP
      WRITE(IW,1101) NB,SCFTYP
1101  FORMAT(//' HCORE, NB =',I5,'   SCFTYP = ',A8)

      OUT=NPRINT.EQ.2.OR.NPRINT.EQ.5

C READ OVERLAP MATRIX TEMPORARILY INTO <H> IF EHT GUESS IS USED

      if (minull.ne.0.) then

        iblok=0
        itype=5
        call secget(isex(itype),itype,iblok)
        call readq(h,nx,iblok,ifild)

C     ----- INITIAL GUESS MATRIX BY EHT -----

        do 5 ii=1,nshell
          mini = kmin(ii)
          maxi = kmax(ii)
          loci = kloc(ii) - mini
          do 5 jj=1,ii
            minj = kmin(jj)
            maxj = kmax(jj)
            locj = kloc(jj) - minj
            iandj = ii.eq.jj
            max = maxj
            do 6 i=mini,maxi
              li = loci + i
              in = (li*(li-1))/2
              if (iandj) max = i
              do 6 j=minj,max
                lj = locj + j
                jn = lj + in
                if (li.eq.lj) then
                  h0(jn) = huckel(ii)
                else
                  h0(jn) = minull*h(jn)*(huckel(ii)+huckel(jj))
                endif
6               continue
5           continue
        write(iw,3000)
        write(10,3000)
3000    format(/24X,33(1H-)
     &         /24X,'INITIAL GUESS FROM THE EHT MATRIX'
     &         /24X,33(1H-))
        if (nprint.eq.2) then
          WRITE(IW,601)
601       FORMAT(//5X,16(1H-),/,5X,'GUESS EHT MATRIX',/,5X,16(1H-))
          CALL PRT2 (H0,COORB,X,NUM,IA,NDIM,IW,BFLAB)
        endif

      else
 
        write(10,10)
        WRITE(IW,10)
10    FORMAT(//9X,'INITIAL GUESS ORBITALS GENERATED FROM ONE-ELECTRON',
     1 ' HAMILTONIAN',/)

C READ IN ONE ELECTRON MATRIX INTO <H0>

        IBLOK=0
        ITYPE=6
        CALL SECGET(ISEX(ITYPE),6,IBLOK)
        CALL READQ(H0,NX,IBLOK,IFILD)

        if (nprint.eq.2) then
          WRITE(IW,600)
600       FORMAT(//,5X,19(1H-),/,5X,'ONE ELECTRON MATRIX',/,5X,19(1H-))
          CALL PRT2 (H0,COORB,X,NUM,IA,NDIM,IW,BFLAB)
        endif
      endif

      CALL READM(COORB,NUM,NCOORB,NDIM,IBLKS,IFILS)

C     ----- FOCK TRANSFORMATION H= Q*H0*Q -----

      CALL FOCKTR(H,H0,COORB,X,IA,NCOORB,NUM,NDIM)

C     ----- DIAGONALIZE H MATRIX -----

*      call hqrii1 (NCOORB,1,NCOORB,-1,H,E,NDIM,V,.TRUE.,h0)
      CALL LIGEN(H,V,E,IA,NCOORB,NDIM,1,1.0D-16)

C     ----- BACK-TRANSFORM THE EIGENVECTORS -----

      CALL BACKTR(V,COORB,X,NCOORB,NUM,NDIM)

C     ----- SELECT MO'S AND OCCUPATION NUMBERS
C           CALCULATE DENSITY MATRICES.

      NOC=NA
      OCC=TWO
      DO 20 I = 1,NCOORB
        X(I) = ZERO
        IF (I.LE.NA) X(I)=OCC
20      CONTINUE
      IF(OUT) WRITE(IW,30)
30    FORMAT(/30X,32(1H-)/30X,'Initial guess molecular orbitals'/
     1 30X,32(1H-))
      NO=1
      IPOS=7
50    IF(NOC.EQ.0)GOTO 60
      IF(CONF.EQ.ALTER) CALL NUCONF(V,NO,NUM,NDIM,IW)
      NOCT=NOC
      IF(SCFTYP.EQ.OSCF) NOCT=NB
      write(10,1102) NB,NOCT,SCFTYP
      WRITE(IW,1102) NB,NOCT,SCFTYP
1102  FORMAT(/' IN HCORE, NB =',I5,'  NOCT =',I5,'  SCFTYP = ',A8)
      IF(OUT) CALL VECOUT(V,E,NUM,NOCT,NDIM,IW,BFLAB)
60    CALL DMTX(H,V,X,NOC,NUM,NDIM)

C     ----- SAVE MO'S + DENSITY + ORBITAL ENERGIES -----

      CALL SAVEMO(H,V,E,NDIM,IPOS)
      IF(SCFTYP.EQ.CLOSED.OR.NB.EQ.0) GO TO 100
      DO 70 I=1,NCOORB
        X(I)=ZERO
        IF(I.LE.NB) X(I)=ONE
        IF(I.LE.NA) X(I)=ZERO
70      CONTINUE
      CALL DMTX(H,V,X,NB,NUM,NDIM)
      NBLOK=LENSEC(NX)
      LDA(ISEX(10))=NX
      IBLOK=0
      ITYPE=10
      CALL SECPUT(ISEX(10),ITYPE,NBLOK,IBLOK)
      CALL WRT3(H,NX,IBLOK,IFILD)
100   CONTINUE
      CONF=BLANK
      RETURN
      END

      SUBROUTINE PRT2 (A,B,C,N,IA,NDIM,IW,BFLAB)
      IMPLICIT REAL*8 (A-H,O-Z)
$INCLUDE:'DIMMM'
      CHARACTER*8 BFLAB(NBF1)
      DIMENSION A(*), B(N,*), C(*), IA(*)
      parameter (ZERO=0.d0)

      DO 200 I=1,N
200      C(I)=ZERO
      DO 500 I=1,N
        DO 500 J=1,N
          IF(J.GT.I) THEN
            IJ=IA(J)+I
          ELSE
            IJ=IA(I)+J
          ENDIF
          B(I,J)=A(IJ)
500        CONTINUE
      CALL VECOUT(B,C,N,N,NDIM,IW,BFLAB)
      RETURN
      END

C  RESTMO  A1

      SUBROUTINE RESTMO(D,V,E,X,NDIM)
      IMPLICIT REAL*8(A-H,O-Z)
$INCLUDE:'DIMMM'
      LOGICAL OUT
      CHARACTER*4 NAM
      LOGICAL  LADAPT,LSKIP
      CHARACTER*8 CLOSED,OSCF
      INTEGER P1,P2,P3
      COMMON /DISC/ ISEL,ISELR,ISELW,IREP,ICHEK,IPOS(16),
     1    IBLKSZ(16)
      COMMON/DISCB/NAM(16)

      CHARACTER*8 PARNAM,BFLAB,ATNAM,ANAM
      COMMON/INFOA/NAT,ICH,MUL,NUM,NX,NE,NA,NB,ZAN(NATM1),C(3,NATM1),
     &             NVAR,NPAR,NVAL,NUQ,
     &             CFR,CFA,KONTYP(NATM1),P1(NATM1),P2(NATM1),P3(NATM1)
     &             ,KR(NATM1),KA(NATM1),KB(NATM1),PARVAL(NATM13)
     &             ,ZIN(NATM1),XTRA(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAM(NATM1)
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
      COMMON/MAPPER/IA(NBF1+1)

      COMMON/GESBLK/ISWOP(2,25),NSWOP,IFILVA,IBLKA1
      COMMON/SMALL/Y(NBF1),YJNK(105)

      CHARACTER*8 BLANK,ALTER
      DIMENSION D(*),E(*),X(*),V(NDIM,*)
      PARAMETER (ZERO=0.0D0,TWO=2.0D0,ONE=1.0D0)

      DATA BLANK/'        '/
      DATA ALTER/'ALTER'/
      DATA CLOSED,OSCF/'CLOSED','OSCF'/

      IBLOK=0
      OUT=NPRINT.EQ.2.OR.NPRINT.EQ.5
      ITYPE=0
*      ITYPE=8
      IF(IBLKA1.EQ.0)CALL SECGET(ISEX(8),ITYPE,IBLKA1)
      IBLVEC=IBLKA1
      WRITE(IW,70) NAM(IFILVA)
      OCC=TWO
      DO 5 I=1,NCOORB
        E(I)=ZERO
        IF(I.LE.NB) E(I)=ONE
5       IF(I.LE.NA) E(I)=OCC
      DO 10 I=1,NCOORB
        X(I)=ZERO
10      IF(I.LE.NA) X(I)=OCC
      NOC=NA
      NO=1
      IPOSD=7
      IF(OUT)WRITE(IW,80)
20    CONTINUE
      CALL READQ(V,NUM*NCOORB,IBLVEC,IFILVA)

c     reading overlap matrix into density matrix vector

      ITYPE=5
      CALL SECGET(ISEX(5),ITYPE,IBLOK)
      CALL READQ(D,NX,IBLOK,IFILD)
      CALL ORTHO(V,D,Y,NUM,NCOORB,NDIM,IA)
      IF(CONF.EQ.ALTER) CALL NUCONF(V,NO,NUM,NDIM,IW)
      NOCT=NOC
      IF(SCFTYP.EQ.OSCF) NOCT=NB
      IF(OUT)CALL VECOUT(V,E,NUM,NOCT,NDIM,IW,BFLAB)

c     calculating density matrix 

30    CALL DMTX(D,V,X,NOC,NUM,NDIM)

C     ----- SAVE MO'S + DENSITY + ORBITAL ENERGIES. -----

      CALL SAVEMO(D,V,E,NDIM,IPOSD)
      IF(SCFTYP.EQ.CLOSED.OR.NB.EQ.0) GO TO 60
      DO 40 I=1,NCOORB
        X(I)=ZERO
        IF(I.LE.NB) X(I)=ONE
        IF(I.LE.NA) X(I)=ZERO
40      E(I)=X(I)
      CALL DMTX(D,V,X,NB,NUM,NDIM)
      LEN=LENSEC(NX)
      ITYPE=10
      IBLOK=0
      CALL SECPUT(ISEX(10),ITYPE,LEN,IBLOK)
      LDA(ISEX(10))=NX
      CALL WRT3(D,NX,IBLOK,IFILD)
60    CONF=BLANK
      RETURN

70    FORMAT(//' Initial guess orbitals to be read from dumpfile ', A4/)
80    FORMAT(/30X,32(1H-)/30X,'Initial guess molecular orbitals'/
     1 30X,32(1H-))
      END

      SUBROUTINE NUCONF(V,NO,NUM,NDIM,IW)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      DIMENSION V(NDIM,*)
      LOGICAL L

      COMMON/GESBLK/ISWOP(2,25),NSWOP,IFILVA,IBLKA1

9999  FORMAT(//5X,'SWOP OPTION IN FORCE'//5X,'INITIAL ORBITAL',
     &' GUESS ALTERED AS FOLLOWS'//)
9998  FORMAT(5X,'ORBITALS',I5,'  AND',I5,' INTERCHANGED')

C  NO is +1 if alpha orbitals are being swopped, -1 if beta

      L=.TRUE.
      IF(NSWOP.EQ.0) RETURN
      DO 1 I=1,NSWOP
        I1=ISWOP(1,I)
        IF (I1*NO .LE. 0) GO TO 1
        I1=IABS(I1)
        I2=IABS(ISWOP(2,I))
        IF (L) WRITE(IW,9999)
        WRITE(IW,9998) I1,I2
        DO 2 J=1,NUM
          XX=V(J,I1)
          V(J,I1)=V(J,I2)
2         V(J,I2)=XX
        L=.FALSE.
1       CONTINUE
      RETURN
      END

      SUBROUTINE SAXPY(N,SA,SX,SY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SX(*),SY(*)

      IF(N.LE.0.OR.SA.EQ.0) RETURN
      DO 10 I=1,N
        SY(I)=SA*SX(I)+SY(I)
10      CONTINUE
      RETURN
      END

