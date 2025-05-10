      SUBROUTINE INIT (ini)

* INITIALIZES ALL DATA TO DEFAULTS

      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      CHARACTER*1 IAAA
      CHARACTER*4 NAM
      CHARACTER*9 STAT
      character*32 NAME,scrat
      character*24 vd
      COMMON/DISC/ISEL,ISELR,ISELW,IREP,ICHEK,IPOS(16),IBLKSZ(16)
      COMMON/DISCB/NAM(16)
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMMAX,TTOT
      COMMON /STREAM/ IREAD,IWRT
      COMMON/WORK/JREC,JUMP,ISTT(40),INUMB(40),IWIDTH
      COMMON/WORKB/IAAA(80)
      LOGICAL FIXED,LEX,LDAM12,LDAM13,LDAM23,LDIIS
      COMMON/SCFBLK/EN,ETOT,EHF,SH1(2),SH2(2),GAP1(2),GAP2(2),
     1              D12,D13,D23,CANA,CANB,CANC,FIXED,LEX,
     2              LDAM12,LDAM13,LDAM23,LDIIS,
     3              NCYC,ISCHM,LOCK,MAXIT,NCONV,LOKCYC
      CHARACTER*8 ANAM,ATNAME,PARNAM,BFLAB
      INTEGER P1, P2, P3
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1               ZAN(NATM1),C(3,NATM1),
     2               NVAR,NPAR,NVAL,NUQ,CFR,CFA,KONTYP(NATM1),
     3               P1(NATM1),P2(NATM1),P3(NATM1),
     4               KR(NATM1),KA(NATM1),KB(NATM1),
     5               PARVAL(NATM13),ZIN(NATM1),X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      COMMON/GESBLK/ISWOP(2,25),NSWOP,IFILVA,IBLKA1
      CHARACTER*8 GROUP
      COMMON/SYMTRY/T(216),INVT(24),ISO(NSH1+1,24),NT,
     &              NT2,XXX(21),INDEX,NAXIS
      COMMON/SYMTRB/GROUP
      LOGICAL LADAPT,LSKIP,DUMPF
      CHARACTER*8 SCFTYP,RUNTYP,GUESS,TITLE,CONF
      COMMON/OPTNSB/TITLE(10),SCFTYP,RUNTYP,GUESS,CONF
      real*4 nopt, minull
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      LOGICAL REVISE
      COMMON/SECTOR/NUM3,IBLK3,REVISE,IODA(24),ICLASS(24),
     1              ILEN(24),MAXB,KBLKLA
      COMMON/MAPPER/IA(NBF1+1)
      COMMON/RTWT/XX,U(NATM1),W(NATM1),NPT,NPT2,DJI(NSH1*2)
      common/files/file5
      character*32 file5
      common/thermo/cpt,cpr,cpv,cpe,
     &              h25t,h25r,h25v,h25e,
     &              s25t,s25r,s25v,s25e,
     &              g25t,g25r,g25v,g25e,
     &              ilin,nax,
     &              temper
     ./origin/ itype, jcar, x0, y0, z0
      character*32 arg
      logical*2 qcl
      integer*2 marg,sta,numarg
      common /cl/ qcl,numarg,marg(5),sta,arg
      dimension ccth(12)
      equivalence (cpt,ccth)
      LOGICAL LMOS,LNUCL
      COMMON/PRPSEC/ISECD1,ITYPD1,ISECD2,ITYPD2,LMOS,ISECV,ITYPV
     1     ,ISECD3,ITYPD3,LNUCL,ORIGIN(3)
      CHARACTER*8 BLANK,GRPC1,HFSCF,HCORE,TET,CLOSED
      DIMENSION ISPACE(1)
      EQUIVALENCE (NPRINT,ISPACE(1))
      save vd
      PARAMETER (PT25=0.25D0, ZERO=0.0D0, ONE=1.0D0, THREE=3.0D0,
     &           PT5=0.5D0, PI=3.14159265358979D0)
      DATA BLANK /'        '/, GRPC1 /'C1'/, TET /'T'/
      DATA HFSCF,HCORE,CLOSED /'HFSCF','HCORE','CLOSED'/

      if (ini.eq.0) TIMLIM=86400.0D0
      TIMMAX=TIMLIM

      do i=1,12
        ccth(i) = ZERO
      enddo
      nax = 0
      ilin = 0
      temper = ZERO

      DO 50 I=1,94
        ISPACE(I)=0
   50   CONTINUE
      DO 10 I=1,NBF1
10      IA(I)=I*(I-1)/2
      NORMF=0
      NORMP=0
      DO 20 I=1,NSH1*2
20      DJI(I)=ONE/(I+I-1)
      DO 30 I=1,10
30      TITLE(I)=BLANK
      itype=0
      IR=7
      IW=9
      IREAD=IR
      IWRT=IW
      if (ini.eq.0) then
        TTGO=TIMLIM
      else
        timmax = ttgo
      endif
*      SAFT=TTGO/100.0D0
*      IF (SAFT.LT.THREE) SAFT=THREE
*      TIMLIM=TTGO-SAFT
      timlim = ttgo
C     DEFAULTS IN OPTIONS BLOCK

      NOPT=3.0
      NVIB=1
      IREST2=0
      IREST3=0
      VIBSIZ=0.001D0
      ITOL=18
      NRUNS=1
      ICUT=10
      IREST=-1
      NPRINT=-1
      IST=1
      JST=1
      KST=1
      LST=1
      INTLOC=1
      NREC=1
      IBLKM=1
      IBLKD=1
      IBLKS=1
      IFILM=1
      IFILD=2
      IFILS=8
      DO 80 I=1,24
80      ISEX(I)=I
      GUESS=HCORE
      RUNTYP=HFSCF
      SCFTYP=CLOSED
      NSWOP=0
      IFILVA=0
      IBLKA1=0
      LADAPT=.FALSE.
      LSKIP=.FALSE.
      MINULL = 0.0

C     DEFAULTS IN BLOCK /INFOA/

      ICH=0
      MUL=1
      DO 90 I=1,NATM1
90      ANAM(I)=BLANK
      DO 100 I=1,NBF1
100     BFLAB(I)=BLANK
      CFR=ONE
      CFA=PI/180.0D0
      NVAR=0
      NPAR=1
      NVAL=1
      PARNAM(1)=TET
      PARVAL(1)=DACOS(-ONE/THREE)

C     SCF BLOCK DEFAULT OPTIONS

      LOCK=0
      LOKCYC=999
      NCONV=5
      MAXIT=50
      NCYC=5
      ISCHM=3
      SH1(1)=PT5
      SH1(2)=PT5
      SH2(1)=PT25
      SH2(2)=PT25
      LDAM12=.FALSE.
      LDAM13=.FALSE.
      LDAM23=.FALSE.
      D12=ONE
      D13=ONE
      D23=ONE
      FIXED=.FALSE.
      GAP1(1)=ONE
      GAP1(2)=ONE
      GAP2(1)=PT5
      GAP2(2)=PT5
      LEX=.TRUE.

C     DEFAULT SYMMETRY

      NT=1
      ISO(1,1)=1
      INDEX=1
      NAXIS=0
      DO 110 I=1,216
110     T(I)=ZERO
      DO 120 I=10,21
120     XXX(I)=ZERO
      T(1)=ONE
      T(5)=ONE
      T(9)=ONE
      XXX(10)=ONE
      XXX(14)=ONE
      XXX(18)=ONE
      GROUP=GRPC1

C     DEFAULT /PRPSEC/

      ISECD1=7
      ISECD2=10
      ISECD3=23
      ISECV=8
      LMOS=.FALSE.
      LNUCL=.TRUE.
      ITYPD1=0
      ITYPD2=0
      ITYPD3=0
      ITYPV=0
      ORIGIN(1)=ZERO
      ORIGIN(2)=ZERO
      ORIGIN(3)=ZERO

      REVISE=.FALSE.
      NUM3=2
      IBLK3=1
      DO 130 I=1,24
130     IODA(I)=-1
      MAXB=8*MAXDIM
      KBLKLA=1
*      WRITE (IW,140)
      WRITE (10,140)
140   FORMAT (/' EXTERNAL FILES:'/)
      DO 160 ISEL=1,16

*     SUBROUTINE ROPEN

        GOTO (1,1,5,5,5,5,5,2,5,5,5,5,5,5,5,5), ISEL

* OPENING DUMP FILES

1         NAME=file5
          call filen1 (1,NAME,NAM(ISEL))
          INQUIRE (FILE=NAME, EXIST=DUMPF)      
          OPEN(ISEL,FILE=NAME,STATUS='UNKNOWN',ACCESS='DIRECT',
     &         FORM='UNFORMATTED',RECL=4096)
          IPOS(ISEL)=0
          IF (DUMPF) THEN
            STAT='AVAILABLE'
          ELSE
            STAT='CREATED'
          ENDIF
*          WRITE (IW,150) STAT,NAME
          WRITE (10,150) STAT,NAME
        GO TO 5

* OPENING SCRATCH FILE (UNIT 7)

2         if (ini.eq.0) then
            if (qcl) then
              if (numarg.gt.3) then
                call rcomlin(4,'INIT  ')
                vd=arg
              else
                vd=char(32)
              endif
            else
              write(10,'(/a/a/a\)')
     & ' If you have a virtual disk with sufficient capacity for the scr
     &atch file,',
     &' enter now the DOS virtual drive and path to write such scratch f
     &ile. The',
     &' default is the current drive => '
              read(10,'(a)') vd
            endif
          endif
          scrat = file5
          call filen2 (vd,scrat)
          call filen1 (1,scrat,nam(isel))
          OPEN (ISEL,FILE=scrat,STATUS='UNKNOWN',
     &          ACCESS='DIRECT',FORM='UNFORMATTED',RECL=4096)
          IPOS(ISEL)=0
*          WRITE (IW,151) 'SCRATCH  : ',scrat
          WRITE (10,151) 'SCRATCH  : ',scrat

5       CONTINUE
160   CONTINUE

      CALL TIMIT(2)

150   FORMAT (1X, A9,': ',A32)
151   FORMAT (/1X, A11,A32)
      RETURN
      END
