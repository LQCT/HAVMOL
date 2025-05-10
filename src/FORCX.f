      subroutine forcx
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'

C     This subroutine calculates the force constant matrix by
C     differencing the gradient of the potential energy; it is
C     an adaptation of A. Komornicki's original routine.

C     Options available are :
C        NPOINTS = 0          Simplex method
C        NPOINTS = 1          Forward difference formula
C        NPOINTS = 2          Central difference formula

C        ( NPOINTS = NVIB  in this subroutine )

      LOGICAL RSTART,SKIP
      INTEGER P1, P2, P3
      CHARACTER*8 RESTRT,PARNAM,BFLAB,ANAME,ATNAME,GROUP
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(NATM13),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13),  ZIN(NATM1), XJ(NATM13)
      COMMON/INFOB/ANAME(NATM1),BFLAB(NBF1),
     &            PARNAM(NATM13),ATNAME(NATM1)

      COMMON/GESBLK/ISWOP(2,25),NSWOP,IFILVA,IBLKA1

      COMMON/OPTMIZ/X(NATM13),G(NATM13),F,X0(NATM13),G0(NATM13),F0
     1,NCOORD,IATOM,ICOORD,IVIB,XJNK(116)
      COMMON/FUNCT/E,EG(NATM13)
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMMAX,TTOT
      COMMON/SYMTRY/T(216),INVT(24),ISO(NSH1+1,24),NT,NT2,
     1    XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
     2    ,U1,U2,U3,V1,V2,V3,W1,W2,W3,XX0,Y0,Z0
     3    ,INDEX,NAXIS
      COMMON/SYMTRB/GROUP
      LOGICAL LADAPT,LSKIP
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
      LOGICAL F2JNK
      COMMON/BIG/FCM(NATM13,NATM13),
     &           F1JNK(228),
     &           I1JNK(NATM1,24),I2JNK(24),
     &           F2JNK(NATM1,NATM1)
      DIMENSION SKIP(NATM1),ALPHA(NATM13),MCOORD(NATM13)
      DIMENSION D(2)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, ONEM=-1.0D0, VIB=1.0D-02)

      DATA RESTRT/'RESTORE'/

10    FORMAT(A//5X,14(5H*****),/,
     1 23X,'FORCE CONSTANT MATRIX CALCULATION',/,5X,14(5H*****))
11    FORMAT(//1X,'-',I1,'- POINT DIFFERENCE FORMULA WITH ',
     1 'STEPSIZE = ',F10.5,' BOHR'/)
12    FORMAT(//7X,'SIMPLEX METHOD WITH STEPSIZE = ',F10.5,
     1 ' BOHR'/)

      ISEC4=0
      ISC24=0
      ISEC17=0
      ISEC18=0
      NCOORD=3*NAT
      LDA(ISEX(17))=2*NATM13+1+2*NATM13+1 + LENINT(4)
      IF(NVIB.GT.2) NVIB=0
      IF(VIBSIZ.EQ.ZERO) VIBSIZ=VIB
      D(1)= VIBSIZ
      D(2)=-VIBSIZ
      CALL SECGET(ISEX(4),4,ISEC4)
      CALL SECGET(ISEX(24),24,ISC24)
      WRITE(IW,10) CHAR(12)
      write(10 ,10) ' '
      IF(NVIB.NE.0) then
        WRITE(IW,11) NVIB,VIBSIZ
        write(10,11) NVIB,VIBSIZ
      else
        WRITE(IW,12) VIBSIZ
        write(10,12) VIBSIZ
      endif

C.... Find which atoms must be displaced and which are related
C     by symmetry and hence do not need to be displaced.

      CALL P1CENT(SKIP)

C.... Check whether this is a restart.

      IF(IREST3.EQ.0)GOTO 20
      CALL SECGET(ISEX(17),17,ISEC17)
      CALL READQ (X,LDA(ISEX(17)),ISEC17,IFILD)
      RSTART=IVIB.GT.0
      IF(.NOT.RSTART) GOTO 20
      IFILVA=IFILD
      IBLKA1=0
      GUESS=RESTRT
      CALL SECGET(ISEX(18),18,ISEC18)
      IF (NVIB.EQ.0) GOTO 200
      IF(IVIB.GT.NVIB) GO TO 150
      GOTO 80

C.... Start a new force-constant calculation:
C           evaluate energy and gradient at the supplied geometry
C           and set force-constant matrix to zero.

20    CONTINUE
      ICOORD=0
      IVIB=0
      IATOM=0
      E=ZERO
      F=ZERO
      F0=ZERO
      DO 30 I=1,NCOORD
        EG(I)=ZERO
        G(I)=ZERO
30      X(I)=C(I)
      IREST3=1
      CALL SECPUT(ISEX(17),17,LENSEC(LDA(ISEX(17))),ISEC17)
      CALL WRT3(X,LDA(ISEX(17)),ISEC17,IFILD)
      CALL TIMIT(3)
      CALL REVIND

C.... Evaluate energy and gradient (except for the simplex method).

      IF(NVIB.NE.0) then
        write (*,'(/a)') ' > FORCX is calling VALFOR'
        write (*,'(a)') ' (Evaluate energy and gradient (except for the
     & simplex method)'
        CALL VALFOR
        call timit(4)
      endif
      F0=F
      DO 40 I=1,NCOORD
        G0(I)=G(I)
        X0(I)=X(I)
40      CONTINUE
      LDA(ISEX(18))=NCOORD*NCOORD
      LENG=LENSEC(NCOORD)
      CALL SECPUT(ISEX(18),18,LENG,ISEC18)
      CALL WRT3(X,LDA(ISEX(17)),ISEC17,IFILD)
      CALL WRT3(VIBSIZ,LDA(ISEX(4)),ISEC4,IFILD)
      CALL WRT3C(TITLE,LDA(ISEX(24)),ISC24,IFILD)
      CALL SETFCM(FCM,NCOORD,NCOORD,G,0,ISEC18,IFILD)
      CALL TIMIT(3)
      CALL REVIND
      IF(NVIB.EQ.0) GOTO 200

C.... Force constants by forward or central difference formula.

      IVIB=1
      IFILVA=IFILD
      IBLKA1=0
      GUESS=RESTRT
50    IF(IVIB.GT.NVIB) GO TO 120
      IATOM=1
60    IF(IATOM.GT.NAT) GO TO 110
      IF(SKIP(IATOM)) GO TO 100
      ICOORD=1
70    IF(ICOORD.GT.3) GO TO 100
80    CONTINUE
      N=3*(IATOM-1)+ICOORD
      DO 90 I=1,NCOORD
        X(I)=X0(I)
90      C(I)=X(I)
      X(N)=X0(N)+D(IVIB)
      C(N)=X(N)
      CALL WRT3(X,LDA(ISEX(17)),ISEC17,IFILD)
      NTSAVE=NT
      NT=1
      write (*,'(/a)') ' > FORCX is calling VALFOR'
      write (*,'(a)') ' (Evaluate energy and gradient)'
      CALL VALFOR
      call timit(4)
      NT=NTSAVE
      CALL SETFCM(FCM,NCOORD,N,G,IVIB,ISEC18,IFILD)
      ICOORD=ICOORD+1
      GO TO 70
100   CONTINUE
      IATOM=IATOM+1
      GO TO 60
110   CONTINUE
      IVIB=IVIB+1
      GO TO 50
120   CONTINUE
      DO 130 IATOM=1,NAT
        IF(.NOT.SKIP(IATOM)) GO TO 140
130   CONTINUE
140   CONTINUE
      ICOORD=1
      CALL WRT3(X,LDA(ISEX(17)),ISEC17,IFILD)
150   CONTINUE
      NDIM=NATM13
      CALL READM(FCM,NCOORD,NCOORD,NDIM,ISEC18,IFILD)

C.... Symmetrize and print out the force-constant matrix :
C        ( A(i,j) + A(j,i) )/2  in lower triangle ;
C        ( A(i,j) - A(j,i) )/2  in upper triangle .

C.... Then calculate normal modes and vibrational frequencies.

160   CONTINUE
      CALL SYMFCM
      RETURN

C.... Force constants by simplex method : evaluate the gradient
C     at the vertices of a simplex in vibrational coordinate space.

200   CONTINUE

C.... Get map from real coords to vibrating coords

      MAP=0
      IATOM=0
210   IATOM = IATOM+1
      IF (IATOM.GT.NAT) GOTO 230
      IF (SKIP(IATOM)) GOTO 210
      N = 3*(IATOM-1)
      DO 220 ICOORD=1,3
220     MCOORD(MAP+ICOORD)=N+ICOORD
      MAP=MAP+3
      GOTO 210

C.... Generate simplex coordinates : number of vertices = MAP

230   DO 250 I=1,MAP
        ALPHA(I) = ONEM/DBLE(FLOAT(I))
        IF (I.EQ.1) GOTO 250
        Z = DSQRT(ONE - ALPHA(I)**2)
        IM = I-1
        DO 240 J=1,IM
240       ALPHA(J)=ALPHA(J)*Z
250     CONTINUE

      IATOM=0
      ICOORD=0
*      NPRINT = -5
      IF(IVIB.GT.1) GO TO 265

C.... First vertex

      IVIB=1
251   CONTINUE
      DO 252  I=1,NCOORD
252     X(I)=X0(I)
      DO 254  I=1,MAP
254     X(MCOORD(I)) = X(MCOORD(I)) + ALPHA(I)*VIBSIZ
      DO 256  I=1,NCOORD
256     C(I)=X(I)
      CALL WRT3 (X,LDA(ISEX(17)),ISEC17,IFILD)
      NTSAVE=NT
      NT=1
      write (*,'(/a)') ' > FORCX is calling VALFOR'
      write (*,'(a)') ' (Evaluate energy and gradient: first vertex)'
      CALL VALFOR
      call timit(4)
      NT=NTSAVE
      CALL SETFCM (FCM,NCOORD,MCOORD(1),G,1,ISEC18,IFILD)

C.... Subsequent vertices :

260   CONTINUE
265   CONTINUE
      IF(IVIB.GT.MAP) GO TO 320
      JVIB=IVIB+1
      DO 270 I=1,NCOORD
270     X(I)=X0(I)
      IFILVA=IFILD
      IBLKA1=0
      GUESS=RESTRT
      X(MCOORD(IVIB)) = X(MCOORD(IVIB)) -
     1      DBLE(FLOAT(IVIB))*ALPHA(IVIB)*VIBSIZ
280   IF (IVIB.GE.MAP) GOTO 300
      DO 290 I=JVIB,MAP
290     X(MCOORD(I))=X(MCOORD(I)) + ALPHA(I)*VIBSIZ
300   CONTINUE
      DO 310 I=1,NCOORD
310     C(I)=X(I)
      CALL WRT3 (X,LDA(ISEX(17)),ISEC17,IFILD)
      NTSAVE=NT
      NT=1
      IVIB=IVIB+1
      write (*,'(/a)') ' > FORCX is calling VALFOR'
      write (*,'(a)') ' (Evaluate energy and gradient: following vertex)
     &'
      CALL VALFOR
      call timit(4)
      NT=NTSAVE

C.... Jump out of loop after MAP vertices

      CALL WRT3 (X,LDA(ISEX(17)),ISEC17,IFILD)
      IF (IVIB.GT.MAP) GOTO 320
      CALL SETFCM(FCM,NCOORD,MCOORD(IVIB),G,1,ISEC18,IFILD)
      GOTO 260

C.... End of gradient evaluations : now construct force-constant matrix.

320   NDIM=NATM13
      CALL READM(FCM,NCOORD,NCOORD,NDIM,ISEC18,IFILD)

C.... Subtract off last point to eliminate G0

      DO 330 I=1,MAP
        DO 330 J=1,NCOORD
330       FCM(J,MCOORD(I))=FCM(J,MCOORD(I)) - G(J)

C.... Backsubstitution to get (unsymmetrised) fcm

      DO 340 I=2,MAP
        DO 340 J=1,NCOORD
340       FCM(J,MCOORD(I))=FCM(J,MCOORD(I))+FCM(J,MCOORD(I-1))
      DO 350 I=2,MAP
        DO 350 J=1,NCOORD
350       FCM(J,MCOORD(I)) = FCM(J,MCOORD(I)) / DBLE(FLOAT(I))
      K=MAP-1
      DO 360 I=1,K
        DO 360 J=1,NCOORD
360       FCM(J,MCOORD(I))=(FCM(J,MCOORD(I))-FCM(J,MCOORD(I+1)))
     &                     /ALPHA(I)
      DO 370 J=1,NCOORD
370     FCM(J,MCOORD(MAP))=FCM(J,MCOORD(MAP))/(ALPHA(MAP)+ONEM)

      NVIB=1
380   CONTINUE
      GOTO 160
      END

*  VALFOR  A1
      SUBROUTINE VALFOR
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'

C     CALCULATE HF ENERGY + GRADIENT WITH RESPECT TO
C     NUCLEAR COORDINATES.

      LOGICAL LSKIP,LADAPT
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
      COMMON/FUNCT/E,EG(NATM13)

      LOGICAL FIXED,LEX,LDAM12,LDAM13,LDAM23,LDIIS
      COMMON/SCFBLK/EN,ETOT,EHF,SH1(2),SH2(2),GAP1(2),GAP2(2),
     1              D12,D13,D23,CANA,CANB,CANC,FIXED,LEX,
     2              LDAM12,LDAM13,LDAM23,LDIIS,
     3              NCYC,ISCHM,LOCK,MAXIT,NCONV,LOKCYC

      INTEGER P1, P2, P3
      CHARACTER*8 ANAME,ATNAME,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(NATM13),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13),  ZIN(NATM1), XJ(NATM13)
      COMMON/INFOB/ANAME(NATM1),BFLAB(NBF1),
     &            PARNAM(NATM13),ATNAME(NATM1)
      COMMON/OPTMIZ/X(NATM13),G(NATM13),F,X0(NATM13),G0(NATM13),F0
     1,NCOORD,IATOM,ICOORD,IVIB,XJNK(116)
      PARAMETER (ZERO=0.0D0, TOLG=1.0D-20)
      character*26 method(3) /
     &                       'Simplex                   ',
     &                       'Forward difference formula',
     &                       'Central difference formula'/
      character*1 xyz(3) /'X', 'Y', 'Z'/

   20 FORMAT(/'   ATOM',I3/
     1 40X,19(1H-)/
     & '   COORD ',a1,30X,'NUCLEAR COORDINATES'/
     2 40X,19(1H-)/
     & '   METHOD: ',a26//
     & 7X,' ATOM',13X,'X',19X,'Y',19X,'Z'/)
   30 FORMAT(1X,I3,2X,A8,2X,3F20.10)

      WRITE(IW,20) IATOM,xyz(ICOORD),method(IVIB+1)
      write(10,20) IATOM,xyz(ICOORD),method(IVIB+1)

      DO 10 I=1,NAT
        N=3*(I-1)
        WRITE(IW,30) I,ANAME(I),(C(N+J),J=1,3)
10      write(10 ,30) I,ANAME(I),(C(N+J),J=1,3)

C     ----- SPECIAL PRINT OPTION IF IATOM.GT.0

      IF(IATOM.GT.0) NPRINT=-5

C     ----- CALCULATE HF ENERGY -----

      IF(IREST.GT.3) GOTO 40
      write (*,'(/a)') ' VALFOR is calling HFSCF'
      write (*,'(a)') ' (Starting SCF procedure)'
      CALL HFSCF
40    DO 50 I=1,NCOORD
   50   EG(I)=ZERO
      write (*,'(/a)') ' VALFOR is calling HFGRAD'
      write (*,'(a)') ' (Calculating gradients)'
      CALL HFGRAD
52    CONTINUE
      DO 60 I=1,NCOORD
        IF(DABS(EG(I)).LT.TOLG) EG(I)=ZERO
        G(I)=EG(I)
   60   X(I)=C(I)
      RETURN
      END

      SUBROUTINE P1CENT(SKIP)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      LOGICAL SKIP

      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),KTYPE(NSH1),
     2             KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      INTEGER P1, P2, P3
      CHARACTER*8 PARNAM,BFLAB,ATNAME,ANAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(NATM13),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13),  ZIN(NATM1), XJ(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),
     &            PARNAM(NATM13),ATNAME(NATM1)
      CHARACTER*8 GROUP
      COMMON/SYMTRY/T(216),INVT(24),ISO(NSH1+1,24),NT,NT2,
     1    XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
     2    ,U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
     3    ,INDEX,NAXIS
      COMMON/SYMTRB/GROUP
      DIMENSION SKIP(NATM1)

      DO 100 IAT=1,NAT
  100   SKIP(IAT)=.TRUE.
      DO 400 II=1,NSHELL
        IAT=KATOM(II)
        DO 300 IT=1,NT
          ID=ISO(II,IT)
          IC=KATOM(ID)
          IF(IC.GT.IAT) GO TO 400
  300     CONTINUE
        SKIP(IAT)=.FALSE.
  400   CONTINUE
      RETURN
      END

C  SETFCM  A1
      SUBROUTINE SETFCM(A,M,N,G,NC,IBLK,IFILD)
c     Store force constant matrix, according the stay and lind of
c     calculation.
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      DIMENSION A(NATM13,M),G(M)
      PARAMETER (NDIM=NATM13, ZERO=0.0D0)

      IF(NC.le.0) then
        DO 10 I=1,M
          DO 10 J=1,M
   10       A(I,J)=ZERO
        CALL WRT3M(A,M,M,NDIM,IBLK,IFILD)
        RETURN
      else
   20   CALL READM(A,M,M,NDIM,IBLK,IFILD)
        IF(NC.le.1) then
          DO 30 I=1,M
   30       A(I,N)=G(I)
          CALL WRT3M(A,M,M,NDIM,IBLK,IFILD)
          RETURN
        endif
   40   CONTINUE
        DO 50 I=1,M
   50     A(I,N)=A(I,N)-G(I)
        CALL WRT3M(A,M,M,NDIM,IBLK,IFILD)
        RETURN
      endif
      END

C  SYMFCM  A1
      SUBROUTINE SYMFCM
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      LOGICAL SKIP
      INTEGER P1,P2,P3
      COMMON/NSHEL/EX(NPRM1),CS(NPRM1),CP(NPRM1),CD(NPRM1),
     1             KSTART(NSH1),KATOM(NSH1),KTYPE(NSH1),
     2             KNG(NSH1),KLOC(NSH1),
     3             KMIN(NSH1),KMAX(NSH1),NSHELL,NSEL2

      CHARACTER*8 PARNAM,BFLAB,ANAM,ATNAME
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(NATM13),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13),  ZIN(NATM1), XJ(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),
     &            PARNAM(NATM13),ATNAME(NATM1)
      CHARACTER*8 GROUP
      COMMON/SYMTRY/T(216),INVT(24),ISO(NSH1+1,24),NT,NT2,
     1    XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
     2    ,U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
     3    ,INDEX,NAXIS
      COMMON/SYMTRB/GROUP
      COMMON/OPTMIZ/X1(NATM13),G1(NATM13),F1,X00(NATM13),G00(NATM13),F00
     &,NCOORD,IATOM,ICOORD,IVIB,XJNK(116)
      COMMON/BIG/A(NATM13,NATM13),
     &           O(3,3),V(3),PTR(3,72),
     &           ICT(NATM1,24),INDIN(24),
     &           SKIP(NATM1,NATM1)
      LOGICAL LADAPT,LSKIP
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
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)

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
        PTR(1,N+1)=XP-X0
        PTR(2,N+1)=YP-Y0
        PTR(3,N+1)=ZP-Z0
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
        PTR(1,N+2)=XP-X0
        PTR(2,N+2)=YP-Y0
        PTR(3,N+2)=ZP-Z0
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
        PTR(1,N+3)=XP-X0
        PTR(2,N+3)=YP-Y0
        PTR(3,N+3)=ZP-Z0
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

C     ----- CHECK WHICH BLOCK OF THE -FCM- HAS BEEN COMPUTED -----

      DO 2000 IAT=1,NAT
        DO 2000 JAT=1,NAT
 2000     SKIP(IAT,JAT)=.FALSE.
      DO 2300 JAT=1,NAT
        DO 2100 IT=1,NT
          IF(ICT(JAT,IT).GT.JAT) GO TO 2300
 2100     CONTINUE
        DO 2200 IAT=1,NAT
 2200     SKIP(IAT,JAT)=.TRUE.
 2300   CONTINUE

C     ----- SYMMETRIZE -FCM- -----

      NCOORD=3*NAT
      IF(NVIB.EQ.2) GO TO 2450
      DO 2400 I=1,NCOORD
        DO 2400 J=1,NCOORD
 2400     A(I,J)=A(I,J)-G00(I)
 2450 DUM=VIBSIZ*NVIB
      DO 2500 I=1,NCOORD
        DO 2500 J=1,NCOORD
 2500     A(I,J)=A(I,J)/DUM
      DO 3600 IAT=1,NAT
        DO 3500 JAT=1,NAT
          IF(SKIP(IAT,JAT)) GO TO 3500
          DO 3400 IT=1,NT
            KAT=ICT(IAT,IT)
            LAT=ICT(JAT,IT)
            IF(.NOT.SKIP(KAT,LAT)) GO TO 3400
            N=3*(IT-1)
            LOCI=3*(IAT-1)
            LOCJ=3*(JAT-1)
            LOCK=3*(KAT-1)
            LOCL=3*(LAT-1)
            DO 2600 K=1,3
              DO 2600 L=1,3
 2600           O(K,L)=A(LOCK+K,LOCL+L)
            DO 2900 K=1,3
              DO 2800 L=1,3
                DUM=ZERO
                DO 2700 M=1,3
 2700             DUM=DUM+O(K,M)*PTR(M,N+L)
 2800           V(L)=DUM
                DO 2900 L=1,3
 2900         O(K,L)=V(L)
            DO 3200 L=1,3
              DO 3100 K=1,3
                DUM=ZERO
                DO 3000 M=1,3
 3000             DUM=DUM+PTR(M,N+K)*O(M,L)
 3100           V(K)=DUM
              DO 3200 K=1,3
 3200           O(K,L)=V(K)
              DO 3300 K=1,3
                DO 3300 L=1,3
 3300             A(LOCI+K,LOCJ+L)=O(K,L)
              SKIP(IAT,JAT)=.TRUE.
              GO TO 3500
 3400       CONTINUE
 3500     CONTINUE
 3600   CONTINUE
      DO 3700 I=1,NCOORD
        DO 3700 J=1,I
          DIFF=(A(I,J)-A(J,I))/TWO
          SUM =(A(I,J)+A(J,I))/TWO
          A(J,I)=DIFF
 3700     A(I,J)=SUM

C     ----- PRINT AND PUNCH FORCE CONSTANT MATRIX -----

      CALL FCMOUT(A,NCOORD,IW)

C     ----- CONVERT TO MASS WEIGHTED COORDINATES AND CALCULATE
C           NORMAL MODES AND VIBRATIONAL FREQUENCIES -----

      CALL FGMTRX(IW)
      RETURN
      END

C  FCMOUT  M2
      SUBROUTINE FCMOUT(FC,NCOORD,IW)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      INTEGER P1, P2, P3
      CHARACTER*4 CLAB
      CHARACTER*8 ANAM,ATNAME,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(NATM13),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13),  ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),
     &            PARNAM(NATM13),ATNAME(NATM1)
      DIMENSION FC(NATM13,NCOORD),CLAB(3)
      DATA CLAB/'X','Y','Z'/

      WRITE(IW,90) CHAR(12)
      MAX=0
   10 MIN=MAX+1
      MAX=MAX+2
      IF(MAX.GT.NAT) MAX=NAT
      WRITE(IW,30)
      WRITE(IW,40) (N,N=MIN,MAX)
      WRITE(IW,30)
      WRITE(IW,50) (ANAM(N),N=MIN,MAX)
      WRITE(IW,30)
      WRITE(IW,60) ((CLAB(M),M=1,3),N=MIN,MAX)
      WRITE(IW,30)
      J0=3*(MIN-1)+1
      J1=3*MAX
      DO 20 IAT=1,NAT
        I0=3*(IAT-1)
        WRITE(IW,70) IAT,ANAM(IAT),CLAB(1),
     1               (FC(I0+1,J),J=J0,J1)
        WRITE(IW,80) CLAB(2) ,(FC(I0+2,J),J=J0,J1)
        WRITE(IW,80) CLAB(3), (FC(I0+3,J),J=J0,J1)
   20   CONTINUE
      IF(MAX.LT.NAT) GOTO 10
      RETURN
   30 FORMAT(/)
   40 FORMAT(20X,2(13X,I3,14X))
   50 FORMAT(20X,2(11X,A8,11X))
   60 FORMAT(20X,2(3X,A4,6X,A4,6X,A4,3X))
   70 FORMAT(1X,I3,2X,A8,1X,A4,1X,6F10.7)
   80 FORMAT(15X,A4,1X,6F10.7)
   90 FORMAT(A//10X,43(1H-),/,10X,'SYMMETRIZED CARTESIAN FORCE',
     1' CONSTANT MATRIX',/,10X,43(1H-))
      END

      SUBROUTINE RAMS(RM,IW)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'

C     ----- CONSTRUCT THE -G- MATRIX.

      INTEGER P1, P2, P3
      CHARACTER*8 ANAM,ATNAME,PARNAM,BFLAB
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(NATM13),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13),  ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),
     &            PARNAM(NATM13),ATNAME(NATM1)
      COMMON/ATMWT/AMS(37),AMASS(natm1)
      DIMENSION RM(NATM13)
      PARAMETER (ONE=1.0D0)

      WRITE(IW,50) CHAR(12)
      write(10 ,50) ' '
      DO 20 IAT=1,NAT
        write(10 ,60) IAT,ANAM(IAT),AMASS(IAT)
20      WRITE(IW,60) IAT,ANAM(IAT),AMASS(IAT)
      K=0
      DO 30 IAT=1,NAT
        DUM=ONE/DSQRT(AMASS(IAT))
        DO 30 J=1,3
          K=K+1
30        RM(K)=DUM
      RETURN
40    FORMAT(8F10.5)
50    FORMAT(A/10X,21(1H-),/,10X,'ATOMIC WEIGHTS (A.U.)',/,
     1 10X,21(1H-),/)
60    FORMAT(I5,5X,A8,2X,F15.5)
      END
