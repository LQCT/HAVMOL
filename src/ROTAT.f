      SUBROUTINE  ROTAT
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'

C    WORKS OUT ROTATIONAL CONSTANTS
C    FROM BERKELEY CODES, MODIFIED BY R D AMOS
C    THERMODYNAMIC PROPERTIES BY LUIS MONTERO

      COMMON/VIB104/PARA,WAVE,CONST,CYCL,CONV

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
     4               PARVAL(NATM13),ZIN(NATM1),XX(3,NATM1)
      CHARACTER*8  ANAM,BFLAB,PARNAM,ATNAME
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),
     1               ATNAME(NATM1)
     ./origin/ itype, jcar, x0, y0, z0

      LOGICAL LADAPT,LSKIP

      COMMON/ATMWT/AMS(37),AMASS(NATM1)

      DIMENSION X(NATM1),Y(NATM1),Z(NATM1),W(NATM1)

      DATA PH,BK,AVN / 6.626176D+00 , 1.380662D+00 , 6.022045D+00 /
      DATA A0 / 0.52917706D+00 /
      DATA CL,TEMP / 2.99792458D+00 , 298.15D+00 /
      DATA PI / 3.1415926536D+00 /

    5 FORMAT(//10x,21('*')/10x,'ROTATIONAL PROPERTIES'/10x,21('*')//
     &        2X,'MOLECULAR GEOMETRY (ATOMIC UNITS)'//
     &        6X,4H NO.,11X,2H X,14X,2H Y,14X,2H Z,12X,'AT.W.'/)
    7 FORMAT(2X,I7,5X,4F16.9)
    8 FORMAT(/2X,'MOLECULAR GEOMETRY (ANGSTROMS)'//
     1 6X,4H NO.,11X,2H X,14X,2H Y,14X,2H Z,12X,'AT.W.'/)

C   PHYSICAL CONSTANTS

      PARA=1.0D+00/AVN
      WAVE=1.0D+04*DSQRT(AVN)/(2.0D+00*PI*CL)
      CONST=1.0D+02*(PH*AVN)/(8.0D+00*PI*PI*CL)
      CYCL=1.0D+06*(PH*AVN)/(8.0D+00*PI*PI)
      CONV=1.0D-01*(PH*CL)/(2.0D+00*BK*TEMP)

      NATOM=NAT
      WRITE(10,5)
      WRITE(9,5)
      DO 210 I=1,NATOM
        X(I)=C(1,I)
        Y(I)=C(2,I)
        Z(I)=C(3,I)
        W(I)=AMASS(I)
        WRITE(10,7) I,X(I),Y(I),Z(I),W(I)
        WRITE(9,7) I,X(I),Y(I),Z(I),W(I)
  210 CONTINUE
      WRITE(10,8)
      WRITE(9,8)
      DO 103 I=1,NATOM
        X(I)=X(I)*A0
        Y(I)=Y(I)*A0
        Z(I)=Z(I)*A0
        WRITE(9,7) I,X(I),Y(I),Z(I),W(I)
        WRITE(10,7) I,X(I),Y(I),Z(I),W(I)
103   CONTINUE
      CALL MOMENT(X,Y,Z,W,NATOM)
      if (itype.eq.1) call movorg
      END

      SUBROUTINE MOMENT(X,Y,Z,W,NATOM)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      DIMENSION T(6),PR(3),EG(3,3),PA(3),pl(2)
      DIMENSION X(NATOM),Y(NATOM),Z(NATOM),W(NATOM)
      COMMON/VIB104/PARA,WAVE,CONST,CYCL,CONV
      COMMON/SYMTRY/TT(216),INVT(24),ISO(NSH1+1,24),NT,NT2,
     1              XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
     2             ,U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
     3             ,INDEX,NAXIS
      COMMON/MAPPER/IKY(NBF1+1)
      COMMON/ATMWT/AMS(37),AMASS(NATM1)
      common/thermo/cpt,cpr,cpv,cpe,
     &              h25t,h25r,h25v,h25e,
     &              s25t,s25r,s25v,s25e,
     &              g25t,g25r,g25v,g25e,
     &              ilin,nax,
     &              temper
      parameter (ZERO=0.0D+00,
     &           UNO=1.d0,
     &           DOS=2.d0,
     &           R=8.3144126d0,
     &           PI=3.1415926536D+00,
     &           A15=1.5d0,
     &           A05=.5d0,
     &           A298=298.15d0,
     &           A25=2.5d0,
     &           AM3=1.0d-3,
     &           C2=1.43878645d0,
     &           FAT=1.1650d0,
     &           FAR=2.7106d0,
     &           FAR2=2.1886d0
     &)


    1 FORMAT(/,2X,' CENTER OF MASS'//8X,' CMX',11X,' CMY',11X,' CMZ'//
     1 2X,3(5X,F10.7)/)
    2 FORMAT(/3X,'CARTESIAN COORDINATES REFERENCED TO CENTER OF MASS'//
     &        7X,4H NO.,9X,2H X,14X,2H Y,14X,2H Z,12X,'AT.W.'/)
    3 FORMAT(2X,I7,5X,4F15.7)
    4 FORMAT(/,2X,' INERTIA TENSOR')
    5 FORMAT(/,2X,' PRINCIPAL MOMENTS OF INERTIA'//
     1 4H NO.,11X,11Hin AMU.A**2,7X,10Hin g.cm**2,8X,
     2 9Hin cm**-1,10X,6Hin MHz/34X,8H(10**39)/)
    6 FORMAT(2X,I2,3X,4D18.5)
    7 FORMAT(/,2X,' EIGENVECTORS OF INERTIA TENSOR'/)
C   THE CALCULATION OF CENTER OF MASS

      SUMW=ZERO
      SUMWX=ZERO
      SUMWY=ZERO
      SUMWZ=ZERO
      DO 101 I=1,NATOM
        SUMW=SUMW+W(I)
        SUMWX=SUMWX+W(I)*X(I)
        SUMWY=SUMWY+W(I)*Y(I)
        SUMWZ=SUMWZ+W(I)*Z(I)
  101 CONTINUE
      CMX=SUMWX/SUMW
      CMY=SUMWY/SUMW
      CMZ=SUMWZ/SUMW
      WRITE(10,1) CMX,CMY,CMZ
      WRITE(9,1) CMX,CMY,CMZ

C   THE CARTESIAN COORDINATES W.R.T. CENTER OF MASS

      WRITE(10,2)
      WRITE(9,2)
      DO 102 I=1,NATOM
        X(I)=X(I)-CMX
        Y(I)=Y(I)-CMY
        Z(I)=Z(I)-CMZ
        WRITE(10,3) I,X(I),Y(I),Z(I),W(I)
        WRITE(9,3) I,X(I),Y(I),Z(I),W(I)
  102 CONTINUE

C   THE CALCULATION OF INERTIA TENSOR

      DO 103 I=1,6
  103   T(I)=ZERO
      DO 104 I=1,NATOM
        T(1)=T(1)+W(I)*(Y(I)**2+Z(I)**2)
        T(3)=T(3)+W(I)*(Z(I)**2+X(I)**2)
        T(6)=T(6)+W(I)*(X(I)**2+Y(I)**2)
        T(2)=T(2)-W(I)*X(I)*Y(I)
        T(4)=T(4)-W(I)*Z(I)*X(I)
        T(5)=T(5)-W(I)*Y(I)*Z(I)
  104 CONTINUE
      AIXX=T(1)
      AIYY=T(3)
      AIZZ=T(6)

c     if any cartesian component is too small, the molecule is linear

      if (aixx.le.AM3 .or. aiyy.le.AM3 .or. aizz.le.AM3) ilin=1

      WRITE(9,4)
      CALL FOUT(T,3,9)

C   THE CALCULATION OF PRINCIPAL MOMENTS OF INERTIA

      write(9,5)
      write(10,5)
      IF(ILIN.EQ.0) then

C   NON LINEAR POLYATOMIC

        I=0
        CALL LIGEN(T,EG,PR,IKY,3,3,1,1.0D-16)
        DO 105 I=1,3
          PA(I)=PR(I)*PARA
          PB=CONST/PR(I)
          PC=CYCL/PR(I)
          WRITE(10,6) I,PR(I),PA(I),PB,PC
          WRITE(9,6) I,PR(I),PA(I),PB,PC
  105   CONTINUE
        WRITE(9,7)
        CALL MOUT(EG,3,3,3,9)
      else

C   DIATOMIC AND LINEAR POLYATOMIC MOLECULES

        DO 106 I=1,3
          II=I*(I+1)/2
          pr(i)=ZERO
          PA(I)=ZERO
          PB=ZERO
          PC=ZERO
          IF(T(II).LE.1.0D-3) then
            T(II)=ZERO
          else
            pr(i)=t(ii)
            PA(I)=T(II)*PARA
            pal = pa(i)
            PB=CONST/T(II)
            PC=CYCL/T(II)
          endif
          WRITE(9,6) I,T(II),PA(I),PB,PC
          WRITE(10,6) I,T(II),PA(I),PB,PC
  106   CONTINUE
      endif
      if (temper.eq.ZERO) temper = A298
      write(9,8) temper
      write(10,8) temper
    8 format (//' THERMODYNAMIC QUANTITIES'/
     &          ' Calculated for the output geometry, as an ideal gas (T
     &=',f10.2,' K)'/)

C     TRANSLATIONAL PROPERTIES

        fa1 = R*FAT
      cpt   = A25*R
      h25t  = cpt*temper
        cc1 = A15*R*log(sumw)
        cc2 = A25*R*log(temper)
      s25t  = cc1 + cc2 - fa1
      g25t  = -(temper*cc1 + temper*cc2 - temper*fa1) + h25t
      write(9,10) cpt,h25t,s25t,g25t
      write(10,10) cpt,h25t,s25t,g25t
10    format (/' Translational:'/
     &         ' - Heat capacity [Cp]         =',g16.7,' J/ K mol'/
     &         ' - Enthalpy [Ho(25)-Ho(0)]    =',g16.7,' J/ mol'/
     &         ' - Entropy [So]               =',g16.7,' J/ K mol'/
     &         ' - Free energy [Go(25)-Ho(0)] =',g16.7,' J/ mol')

C     ROTATIONAL PROPERTIES

      if (ilin.eq.0) then

c       non linear polyatomic molecule routine

        if (nax.eq.0) then
          if (naxis.eq.0) then
            sigma=UNO
          else
            sigma=dble(naxis)
          endif
        else
          sigma = dble(nax)
        endif

c       formulae of Knox

        cpr   = A15*R
        h25r  = cpr*temper
         cc1  = R*A05*log(pr(1)*pr(2)*pr(3))
         cc2  = R*A15*log(temper)
        s25r  = cc1 + cc2 - R*log(sigma) - R*FAR
        g25r  = h25r - temper*s25r
      else

c       routine for linear molecules

        k = 0
        do 500 i=1,3
          if (pr(i).eq.ZERO) goto 500
          k = k + 1
          pl(k) = pr(i)
500       continue
        if (cmx.eq.ZERO .and. cmy.eq.ZERO .and. cmz.eq.ZERO) then
          sigma = DOS
        else
          sigma = UNO
        endif                

c       formulae of JANAF

        br = 2.79932d-39/(pal*1.d-39)
        cc = C2*br/temper
        cc1= 2.2222222222d-2*R*cc**2
        cpr  = R + cc1
        rt = R*temper
        h25r = rt - 3.33333333333d-1*cc*rt - cc1*temper
        s25r = R - R*log(cc*sigma) - 1.111111111111d-2*R*cc
        g25r = rt*log(cc*sigma) - 3.333333333333d-1*rt*cc -
     &         1.1111111111111d-2*rt*cc**2

c        formulae of Knox

*        cpr  = DOS*R
*          rt = R*temper
*        h25r = rt
*        s25r = R*log(pl(1)) + R*log(temper) - R*log(sigma) - R*FAR2
*        g25r = h25r - temper*s25r
      endif

      write(9,11) cpr,h25r,s25r,g25r
      write(10,11) cpr,h25r,s25r,g25r
11    format (/' Rotational:'/
     &         ' - Heat capacity [Cp]         =',g16.7,' J/ K mol'/
     &         ' - Enthalpy [Ho(25)-Ho(0)]    =',g16.7,' J/ mol'/
     &         ' - Entropy [So]               =',g16.7,' J/ K mol'/
     &         ' - Free energy [Go(25)-Ho(0)] =',g16.7,' J/ mol')

      RETURN
      END

      SUBROUTINE MOUT(W,NBASIS,NCOL,NDIM,IW)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      DIMENSION W(*)
10    FORMAT(/3X,8I14)
20    FORMAT(/)
30    FORMAT(7X,I3,8F14.7)
      M=1
      N=5
40    IF(NCOL.LT.M)RETURN
      IF(N.GT.NCOL)N=NCOL
      WRITE(IW,10)(I,I=M,N)
      WRITE(IW,20)
      DO 50 J=1,NBASIS
50      WRITE(IW,30)J,(W(J+(I-1)*NDIM),I=M,N)
      M=M+5
      N=N+5
      GOTO 40
      END
