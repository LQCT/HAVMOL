      subroutine upcase (ir,line)
*     READS AND CONVERTS LOWER TO UPPERCASE LETTERS
      character*80 line
      integer*4 ir
*      if (ir.ge.0) read (ir,'(a)',end=50) line
      if (ir.ge.0) read (ir,'(a)') line
      call uc (line)
      return
*50    write(10,'(a,i4)') 'EOF found in input file ',ir
*      call honder(81,'UPCASE')
*      return
      end

      subroutine uc (line)
      implicit integer*2 (i-m)
      character*1 line(80)
      do l=80,1,-1
        if (line(l).ne.char(32)) then
          ncar = l
          exit
        endif
      enddo
      do 15 j=1,ncar
        do 15 i=97,122
          ich = ichar(line(j))
          if (ich .eq. i) line(j)=CHAR(i-32)
          if (ich.eq.13 .or. ich.eq.10) exit
15        continue
      return
      end

      subroutine filen1 (kk,fname,ext)
      implicit integer*2 (i-n)

c     Default extension for DOS filenames. If:

c     kk.eq.0 then fname remains unchanged with the original extension
c                        if it exists, and if it is not, 'ext' is added.
c     kk.gt.0            changes its extension to the one in parameter
c                        'ext', in all cases.
c     kk.lt.0            returns the existing extension in 'ext',
c                        without changes
 
      character*1 fname(32), ext(4), BL/' '/, PT/'.'/
      integer*4 kk

      do 10 i=1,32
        if (fname(i).eq.PT) then
c         the file name has an explicit extension
          if (kk.le.0) go to 50
          exit
        else if (fname(i).eq.BL) then
          exit
        endif
10    continue

      do j=1,4
        fname(i-1+j) = ext(j)
      enddo
      return

50    continue
      if (kk.lt.0) then
        do k=1,4
            ext(k)=fname(i+k-1)
        enddo      
      endif
      return
      end

      BLOCK DATA DATIN
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      CHARACTER*4 NAM
      parameter (ZERO=0.d0)

C     INITIAL VALUES FOR OPTIONS

      COMMON/OPTNS/IDUM(104)

C     DATA FOR INTEGRATION POINTS AND WEIGHTS

      COMMON/RTDATA/WHI(15),WLOW(15),RHI(15),RLOW(15),RFAC(15),
     1AMPS(9),IPOINT(9),MADD(20)

C     ATMOL INPUT

      COMMON/WORK/JREC,JUMP,ISTRT(40),INUMB(40),IWIDTH
      DIMENSION ITAP(81)
      EQUIVALENCE (ITAP,JUMP)
      COMMON/DISCB/NAM(16)
      COMMON/DISC/ISEL,ISELR,ISELW,IREP,ICHECK,IPOS(16),
     1IBLKSZ(16)
      COMMON/HUCKEL/HUCKEL(NSH1),ELEMS(3,18),ELEMP(2,18)
      COMMON/ATMWT/AMS(37),AMASS(natm1)

C  Hni and Wni are the i-th point and weight for n-point Gauss-Hermite
C  quadrature.

      COMMON/HERMIT/H11,H21,H22,H31,H32,H33,H41,H42,H43,H44,
     1    H51,H52,H53,H54,H55,H61,H62,H63,H64,H65,H66,
     2    H71,H72,H73,H74,H75,H76,H77
      COMMON/WERMIT/W11,W21,W22,W31,W32,W33,W41,W42,W43,W44,
     1    W51,W52,W53,W54,W55,W61,W62,W63,W64,W65,W66,
     2    W71,W72,W73,W74,W75,W76,W77

      DATA AMS/1.007825D0,4.002604D0,7.01600D0,9.01218D0,11.009305D0,
     &12.0000000D0,14.003074D0,15.994915D0,18.99840D0,19.992440D0,
     &22.98977D0,23.98504D0,26.98153D0,27.97693D0,30.97376D0,
     &31.97207D0,34.96885D0,39.96238D0,38.96371D0,39.96259D0,44.95592D0,
     &47.94795D0,50.94398D0,51.94051D0,54.93805D0,55.9349D0,58.93319D0,
     &57.9353D0,62.9295D0,63.92915D0,68.9257D0,73.921D0,74.921D0,
     &9.916D0,78.914D0,83.911D0,0.0D0/
      data amass/natm1*0.0d0/

C     SET DATA IN /RTDATA/ USED BY INTEGRATION ROUTINES

      DATA MADD/17,17,18,18,19,19,20,20,21,21,22,22,
     &   23,23,24,24,25,25,26,26/

C     ROOTS AND WEIGHTS GENERATED USING NAG ROUTINE D01BCF
C     RHI CONTAINS SQUARES OF GAUSS-HERMITE ROOTS OF ORDER
C     2,4,6.......18
C     WHI ARE CORRESPONDING WEIGHTS
C     RLOW CONTAINS SQUARES OF GAUSS-LEGENDRE ROOTS OF ORDER
C     2,4,6,......,18
C     WLOW ARE CORRESPONDING WEIGHTS
C     RFAC CONTAINS RATIOS RLOW/RHI

      DATA RHI/
     &0.50000000000000D+00,0.27247448713916D+01,0.27525512860841D+00,
     &0.55253437422632D+01,0.17844927485433D+01,0.19016350919349D+00,
     &0.85886356890120D+01,0.39269635013583D+01,0.13390972881264D+01,
     &0.14530352150332D+00,0.11807189489972D+02,0.64147297336620D+01,
     &0.30859374437175D+01,0.10745620124369D+01,0.11758132021178D+00/
      DATA WHI/
     &0.88622692545276D+00,0.81312835447246D-01,0.80491409000551D+00,
     &0.45300099055090D-02,0.15706732032286D+00,0.72462959522439D+00,
     &0.19960407221138D-03,0.17077983007414D-01,0.20780232581489D+00,
     &0.66114701255824D+00,0.76404328552331D-05,0.13436457467813D-02,
     &0.33874394455481D-01,0.24013861108232D+00,0.61086263373533D+00/
      DATA RLOW/
     &0.33333333333333D+00,0.74155574714581D+00,0.11558710999705D+00,
     &0.86949939491826D+00,0.43719785275109D+00,0.56939115967007D-01,
     &0.92215660849206D+00,0.63467747623464D+00,0.27618431387246D+00,
     &0.33648268067507D-01,0.94849392628837D+00,0.74833462838728D+00,
     &0.46159736149627D+00,0.18783156765245D+00,0.22163568807218D-01/
      DATA WLOW/
     &0.10000000000000D+01,0.34785484513745D+00,0.65214515486255D+00,
     &0.17132449237917D+00,0.36076157304814D+00,0.46791393457269D+00,
     &0.10122853629038D+00,0.22238103445337D+00,0.31370664587789D+00,
     &0.36268378337836D+00,0.66671344308689D-01,0.14945134915058D+00,
     &0.21908636251598D+00,0.26926671931000D+00,0.29552422471475D+00/
      DATA RFAC/
     &0.66666666666667D+00,0.27215603006790D+00,0.41992718021791D+00,
     &0.15736566546394D+00,0.24499839134006D+00,0.29942188282334D+00,
     &0.10736939391571D+00,0.16162041638918D+00,0.20624663818033D+00,
     &0.23157228207121D+00,0.80331896688366D-01,0.11665879303696D+00,
     &0.14958091987121D+00,0.17479825778177D+00,0.18849566212812D+00/
      DATA JREC/-1/,IWIDTH/80/, ITAP/81*0/
      DATA ICHECK,ISEL,ISELR,ISELW,IREP/5*0/
      DATA IPOS/16*-1/
      DATA IBLKSZ/8*0,8*4096/
      DATA NAM/'.ED0','.ED1','.ED2','.ED3',
     &         '.ED4','.ED5','.ED6','.ED7',
     &          'MT0', 'MT1', 'MT2', 'MT3',
     &          'MT4', 'MT5', 'MT6', 'MT7'/
      DATA AMPS/24.0D0,30.0D0,35.0D0,39.0D0,43.0D0,47.0D0,
     $51.0D0,55.0D0,60.0D0/
      DATA IPOINT/0,1,3,6,10,15,21,28,36/
      DATA H11 /0.0D+00/
      DATA W11 /1.77245385090552D+00/

      DATA H21,H22 /-.707106781186548D+00,.707106781186548D+00/
      DATA W21,W22 /.8862269254528D+00,.8862269254528D+00/

      DATA H31,H32,H33 /-1.22474487139159D+00,0.0D+00,
     1    1.22474487139159D+00/
      DATA W31,W32,W33 /0.2954089751509D+00,1.181635900604D+00,
     1    0.2954089751509D+00/

      DATA H41,H42,H43,H44
     1    /-1.65068012388578D+00,-0.524647623275290D+00,
     2    0.524647623275290D+00,1.65068012388578D+00/
      DATA W41,W42,W43,W44 /8.131283544725D-02,8.049140900055D-01,
     1    8.049140900055D-01,8.131283544725D-02/

      DATA H51,H52,H53,H54,H55 /-2.02018287045609D+00,
     1    -0.958572464613819D+00, 0.0D+00,
     2    0.958572464613819D+00, 2.02018287045609D+00 /
      DATA W51,W52,W53,W54,W55 /1.995324205905D-02,3.936193231522
     1D-01, 9.453087204829D-01,3.936193231522D-01,1.995324205905D-02/

      DATA H61,H62,H63,H64,H65,H66 /-2.350604973674D+00,
     1    -1.335849074014D+00,-0.436077411928D+00,0.436077411928D+00,
     2    1.335849074014D+00, 2.350604973674D+00/
      DATA W61,W62,W63,W64,W65,W66 / 4.530009905509D-03,
     1    1.570673203229D-01, 7.246295952244D-01,7.246295952244D-01,
     2    1.570673203229D-01, 4.530009905509D-03/

      DATA H71,H72,H73,H74,H75,H76,H77 /-2.651961356835233D0,
     1    -1.673551628767471D0, -0.816287882858965D0, 0.0D0,
     2    0.816287882858965D0, 1.673551628767471D0, 2.651961356835233D0/
      DATA W71,W72,W73,W74,W75,W76,W77 /9.717812450995D-4,
     1    5.451558281913D-2, 4.256072526101D-1, 8.102646175568D-1,
     2    4.256072526101D-1, 5.451558281913D-2, 9.717812450995D-4 /
      DATA IDUM/104*0/

c Atomic Orbital Eigenvalues from:
c E. Clementi, "Tables of Atomic Functions", suppl. of IBM Journal
c of Research and Development, 9, 2 (1965)

      DATA ELEMS
     &/ -.5d0,                  2*ZERO,  -.91795d0,              2*ZERO,
     & -2.47994d0, -.19259d0,     ZERO, -4.76131d0, -.31095d0,     ZERO,
     & -7.69448d0, -.49411d0,     ZERO,-11.32333d0, -.70429d0,     ZERO,
     &-15.62433d0, -.94239d0,     ZERO,-20.66330d0,-1.24081d0,     ZERO,
     &-26.37451d0,-1.56700d0,     ZERO,-32.76074d0,-1.92239d0,     ZERO,
     &-40.46668d0,-2.79022d0,-.18114d0,-49.02498d0,-3.76463d0,-.25227d0,
     &-58.49929d0,-4.91012d0,-.39269d0,-68.81016d0,-6.15512d0,-.53801d0,
     &-79.98579d0,-7.50438d0,-.69525d0,-92.00999d0,-9.00629d0,-.87912d0,
     &-108.87793d0,-10.60076d0,-1.06552d0,
     &-118.60161d0,-12.31329d0,-9.56285d0/
      DATA ELEMP
     &/               2*ZERO,                 2*ZERO,
     &                2*ZERO,                 2*ZERO,
     5 -.30993d0,       ZERO,  -.43265d0,       ZERO,
     7 -.56548d0,       ZERO,  -.62866d0,       ZERO,
     9 -.72438d0,       ZERO,  -.84201d0,       ZERO,
     1-1.51124d0,       ZERO, -2.27874d0,       ZERO,
     3-3.21741d0,  -.20950d0, -4.25419d0,  -.29587d0,
     5-5.41341d0,  -.39192d0, -6.68963d0,  -.43708d0,
     7-8.06460d0,  -.49991d0, -9.56285d0,  -.58375d0/

      END
      subroutine pathf (ncb, buff, ext, buffo)
      implicit integer*2 (i-n)
      integer*4 ncb
      character*1 buff(ncb), ext(4), buffo(32), file5(32)
      common
     &/files/ file5

      do 5 i=1,ncb
        if (buff(i).eq.char(32)) then
          ncar = i - 1
          goto 6
        endif
5       continue
6     continue

      ibs = 0
      do 10 i=1,32
        if (file5(i).eq.'\') then
          ibs = i
        endif
10      continue

20    if (ibs.ne.0) then

c caso en el que existe un "path" en el fichero inicial
 
        do 30 j=1,ibs
30        buffo(j) = file5(j)
        do 32 l=1,ncar
32        buffo(ibs+l) = buff(l)
        do 34 j=1,ncar
          if (buff(j).eq.'.') then
            return
          endif
34        continue
        do 36 k=1,4
36        buffo(ibs+ncar+k) = ext(k)
          return
      elseif (file5(2).eq.':') then

c caso en el que existe un "drive" en el fichero inicial
 
        buffo(1) = file5(1)
        buffo(2) = file5(2)
        do 42 l=1,ncar
42        buffo(2+l) = buff(l)
        do 44 j=1,ncar
          if (buff(j).eq.'.') then
            return
          endif
44        continue
        do 46 k=1,4
46        buffo(2+ncar+k) = ext(k)
          return

      else

c caso en el que no existe ni "drive" ni "path" en el fichero inicial

        do 52 l=1,ncar
52        buffo(l) = buff(l)
        do 54 j=1,ncar
          if (buff(j).eq.'.') then
            return
          endif
54        continue
        do 56 k=1,4
56        buffo(ncar+k) = ext(k)
          return
      endif

      end

      subroutine filen2 (dd,fname)
      implicit integer*2 (i-n)

C path PARA LOS NOMBRES DE FICHERO DOS

      character*1 fname(*),dd(*),name(32),bl

      bl = char(32)

C SI NO SE ENTRO path SE TERMINA

      if (dd(1).eq.bl) go to 50

C CONTEO DE LOS CARACTERES DEL path Y DEPOSITO EN ND

      do 5 i=2,32
        if (dd(i).eq.bl) then
          nd = i - 1
          go to 6
        endif
5       continue

C CONTEO DE CARACTERES DEL path DEL NOMBRE DEL FICHERO ORIGINAL

6     l = 0
      i = 1
10    if (fname(i).eq.':') then
        l = i
        do 15 j=i+1,32
          if (fname(j).eq.bl) then
            i = j
            go to 20
          endif
15        if (fname(j).eq.'\') l = j
      endif
      if (fname(i).eq.bl) go to 20
      i = i + 1
      if (i.le.32) go to 10

C CONSTRUCCION DEL NUEVO NOMBRE DE FICHERO

20    n = i - 1
      do 25 j=1,nd
25      name(j) = dd(j)
      do 30 j=1,n
30      name(nd+j) = fname(l+j)
      do 40 j=1,nd+(n-l)
40      fname(j) = name(j)

50    return
      end

C PROP1E
      SUBROUTINE PROP1E
C     DRIVING ROUTINE FOR ONE-ELECTRON PROPERTIES
      IMPLICIT REAL*8(A-H,O-Z)
      common
     ./origin/ itype, jcar, x0, y0, z0
      CALL MULKEN
      CALL DIPOLE
      CALL QPOLE
      CALL ELECTR(.FALSE.)
      if (itype.eq.2) call movorg
      RETURN
      END

      SUBROUTINE HONDER (N,SUBER)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 IA,SPACE,Q
      character*6 suber
      COMMON/STREAM/IR,IW
      COMMON/WORK/JREC,JUMP,ISTRT(40),INUMB(40),IWIDTH
      COMMON/WORKB/IA(80)
      DATA    SPACE /' '/, Q /'^'/
C
500   WRITE(IW,502) N,SUBER
      write(*,502) N,SUBER
502   FORMAT(/1X,'------ HONDO error ',I3,' in subroutine ',a)

      WRITE(IW,501)
      write(*,501)
501   FORMAT(/' Current input line is')

      WRITE(IW,1) (IA(I),I=1,79)
      write(*,1) (IA(I),I=1,79)
1     FORMAT(1X,79A1)

      IF (JREC .EQ. 0) K=1
      IF (JREC .GT. 0) K=ISTRT(JREC)
      IF (JREC .GT. JUMP) K=ISTRT(JUMP)+INUMB(JUMP)+1
      WRITE (IW,1001) (SPACE, I=1,K), Q
      WRITE (*,1001) (SPACE, I=1,K), Q
1001  FORMAT (80A1)

      CALL CLENUP
      STOP '*** MICROMOL ABORTED ***'
      END

      SUBROUTINE CLENUP
* Ends the run
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/STREAM/IR,IW
      CALL TIMIT (0)
      CALL REVIND
      CALL SECSUM
      CLOSE (8, STATUS='DELETE')
      RETURN
      END

      SUBROUTINE REVIND
* Recording total elapsed time and parameters for resuming a run.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/OPTNS/TIML1M,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      LOGICAL REVISE
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMLIM,TTOT
      COMMON/SECTOR/NUM3,IBLK3,REVISE,
     1ioda(24),ICLASS(24),ILEN(24),MAXB,KBLKLA
      DIMENSION C(37)
      EQUIVALENCE (ioda,C)
      PARAMETER (M37=37)

      TTOT1 = TTOT
      if (lda(isex(4)).ne.0) then
        CALL SECPUT(ISEX(4),4,LENSEC(LDA(isex(4))),IBLOK)
        CALL WRT3(VIBSIZ,LDA(isex(4)),IBLOK,NUM3)
        IF(.NOT.REVISE) THEN
          CALL SEARCH(IBLK3,NUM3)
          CALL PUT(C,M37,NUM3)
          REVISE=.TRUE.
        ENDIF
      else
        write(*,1000)
        write(iw,1000)
1000    format (/' *** NO DATA AVAILABLE TO RESTART RUN ***'/
     &           '     Please, beguin again ...'/)
      endif
      RETURN
      END

      SUBROUTINE SEARCH(IBLK,IUNIT)
c gives the IBLK value to IPOS of unit IUNIT
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/DISC/ISEL,ISELR,ISELW,IREP,ICHECK,IPOS(16),IBLKSZ(16)

      IPOS(IUNIT)=IBLK
      RETURN
      END

      SUBROUTINE SECSUM
* Output dump file parameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*4 IED
      COMMON/OPTNS/TIMLIM,VIBSIZ,NOPT,NVIB,IREST2,IREST3,
     &            LADAPT,LSKIP,NPRINT,ITOL,ICUT,NORMF,NORMP,
     &            IREST,NREC,INTLOC,IST,JST,KST,LST,IR,IW,
     &            IBLKD,IFILD,IBLKS,IFILS,IBLKM,IFILM,MBLKM,
     &            NOTAPE(3),IBLK(3),LBLK(3),ILOW,NCOORB,NSA,
     &            MINVEC,DIISE,DIISD,MINULL,
     &            ICONV,NRUNS,NMUL(2),IPOSM,LDA(24),ISEX(24),
     &            TTOT1
      LOGICAL REVISE
      COMMON/SECTOR/NUM3,IBLK3,REVISE,
     1ioda(24),ICLASS(24),ILEN(24),MAXB,KBLKLA
      COMMON/DISCB/IED(16)
      character*16 dsect(24)
      data dsect /' Symmetry Data 1',
     &            ' System Data 1  ',
     &            ' Basis Functions',
     &            ' Optns.& Rest. 1',
     &            ' Overlap Matrix ',
     &            ' One Electr. H  ',
     &            ' Density Matrix ',
     &            ' Wave Functions ',
     &            ' Eigenvalues    ',
     &            ' Dens.Mat.OSHELL',
     &            ' Common HUCKEL  ',
     &            ' Common ATMWT   ',
     &            ' SCF Data       ',
     &            ' Gradients      ',
     &            ' Opt.Coord.&Data',
     &            ' Inv.Hessian    ',
     &            ' Force Cnst.Data',
     &            ' Forc.Cnst.Mat. ',
     &            ' Lagrang. OSHELL',
     &            ' K Mat. OSHELL  ',
     &            ' Symmetry Data 2',
     &            ' System Data 2  ',
     &            ' Dens.mat.optn. ',
     &            ' Optns.& Rest. 2'/

      WRITE(9,30)IED(NUM3),IBLK3,KBLKLA,MAXB
      write(10,30)IED(NUM3),IBLK3,KBLKLA,MAXB
30    FORMAT(//' *SUMMARY OF DUMPFILE ON <fn>',A4,' AT BLOCK',I6/2H *,/
     1' *CURRENT LENGTH =',I6,' BLOCKS'/
     2' *MAXIMUM LENGTH =',I6,' BLOCKS'/2H *)
      WRITE (9,31)
31    FORMAT
     &(' *       SECTION         TYPE  BLOCK LENGTH/4096 bytes')
      DO 60 I=1,24
        IF(ioda(I))60,60,40
40      M=IBLK3+ioda(I)
        WRITE(9,50)I,dsect(i),ICLASS(I),M,ILEN(I)
60      CONTINUE
50    FORMAT(2H *,I3,2x,a16,I5,I6,I11)

      WRITE (9,70) TTOT1,TTOT1/3600.D0
      WRITE (*,70) TTOT1,TTOT1/3600.D0
70    FORMAT (/' TOTAL TIME ELAPSED BY THIS SYSTEM =',E16.8,' SECS.'/
     &         '                                   =',F16.5,' HOURS')
      RETURN
      END

C  SIMPLIFIED I/O ROUTINES

      SUBROUTINE FMOVE(A,B,NWORD)
      IMPLICIT REAL*8 (A-H,O-Z)
C  Moves NWORD elements of array A into array B
      DIMENSION A(*),B(*)

      DO 10 I=1,NWORD
        B(I)=A(I)
10      CONTINUE
      RETURN
      END

      SUBROUTINE FMOVEC(A,B,NWORD)
      IMPLICIT REAL*8 (A-H,O-Z)
C  Moves NWORD elements of array A into array B
      CHARACTER*8 A,B
      DIMENSION A(*),B(*)
      DO 10 I=1,NWORD
      B(I)=A(I)
10    CONTINUE
      RETURN
      END

      SUBROUTINE SECINI(IBL,NUM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL REVISE
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMLIM,TTOT
      COMMON/SECTOR/NUM3,IBLK3,REVISE,
     1ioda(24),ICLASS(24),ILEN(24),MAXB,KBLKLA
      DIMENSION C(37)
      EQUIVALENCE (ioda,C)
      DATA M37/37/

      NUM3=NUM
      IBLK3=IBL
      CALL SEARCH(IBLK3,NUM3)
      REVISE=.TRUE.
      CALL GET(C,NW,NUM3)
      IF(NW.NE.M37)CALL HONDER(61,'SECINI')
      RETURN
      END

      SUBROUTINE SECPUT (MPOS,MTYPE,LENGTH,IBLOCK)
* Initializes or actualizes the block number IBLOCK and type MTYPE of
* section MPOS. If initializes, it gives the length LENGTH to section
* MPOS.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL REVISE
      CHARACTER*6 SUBER /'SECPUT'/
      COMMON/STREAM/IR,IW
      COMMON/SECTOR/NUM3,IBLK3,REVISE,
     1ioda(24),ICLASS(24),ILEN(24),MAXB,KBLKLA
      PARAMETER (M1=1)
      IF (MPOS.LT.M1) CALL HONDER(63,SUBER)
      IF(ioda(MPOS).LE.0) THEN
        REVISE=.FALSE.
        M=LENGTH+KBLKLA
        IF(M.GT.MAXB)CALL HONDER(71,SUBER)
        IBLOCK=KBLKLA+IBLK3
        ioda(MPOS)=KBLKLA
        ICLASS(MPOS)=MTYPE
        ILEN(MPOS)=LENGTH
        KBLKLA=M
      ELSE
        IF(LENGTH.GT.ILEN(MPOS)) THEN
          WRITE (10,100) LENGTH, MPOS, ILEN(MPOS)
          WRITE (IW,100) LENGTH, MPOS, ILEN(MPOS)
          CALL HONDER(72,SUBER)
        ENDIF
        IF(MTYPE.NE.ICLASS(MPOS))REVISE=.FALSE.
        ICLASS(MPOS)=MTYPE  
        IBLOCK=ioda(MPOS)+IBLK3
      ENDIF
100   FORMAT (/' ERROR: LENGTH',I5,' IS NOT PREVIOUSLY ASSIGNED TO MPOS'
     &       ,I2/' WHICH CURRENT LENGTHS IS',I5)
      RETURN
      END

      SUBROUTINE WRT3(Q,NWORD,IBLK,NUM3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(*)
      CALL SEARCH(IBLK,NUM3)
      L=511
      J=1
      K=NWORD
10    K=K-511
      IF(K)20,30,30
20    L=K+511
30    CALL PUT(Q(J),L,NUM3)
      J=J+511
      IF(K)40,40,10
40    RETURN
      END

      SUBROUTINE WRT3M(Q,N,M,NMAX,IBLK,NUM3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(NMAX,M)
      COMMON/BLKIN/A(511)
      DATA M1,M511/1,511/
      CALL SEARCH(IBLK,NUM3)
      K=M1
      DO 3 J=1,M
      DO 3 I=1,N
      A(K)=Q(I,J)
      K=K+M1
      IF(K.LE.M511) GOTO 3
      K=M1
      CALL PUT(A,M511,NUM3)
3     CONTINUE
      IF(K.GT.M1) CALL PUT(A,M511,NUM3)
      RETURN
      END

      SUBROUTINE PUT(C,NWORD,IUNIT)
* Writing the buffer NWORD bytes of BUFF into the 511 bytes of file IUNIT
* at IREC=IPOSN
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BUFIN/BUFF(511)
      DIMENSION C(*)
      COMMON/DISC/ISEL,ISELR,ISELW,IREP,ICHECK,IPOS(16),IBLKSZ(16)
      IPOSN=IPOS(IUNIT)
      CALL FMOVE(C,BUFF,NWORD)
      WRITE(IUNIT,REC=IPOSN,ERR=999) (BUFF(I),I=1,511),NWORD
      IPOS(IUNIT)=IPOSN+1
      RETURN
999   write(10,'(a/3(a,i4),a/)')
     &      ' *** Output error to direct access file ***',
     &      ' *** UNIT =',iunit,'; REC =',iposn,'; NWORD =',nword,' ***'
      CALL HONDER(73,'PUT')
      END

      SUBROUTINE WRT3C(Q,NWORD,IBLK,NUM3)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 Q
      DIMENSION Q(*)
      CALL SEARCH(IBLK,NUM3)
      L=127
      J=1
      K=NWORD
10    K=K-127
      IF(K)20,30,30
20    L=K+127
30    CALL PUTC(Q(J),L,NUM3)
      J=J+127
      IF(K)40,40,10
40    RETURN
      END

      SUBROUTINE PUTC(C,NWORD,IUNIT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 C,CUFF
      COMMON/CUFIN/CUFF(127)
      DIMENSION C(*)
      COMMON/DISC/ISEL,ISELR,ISELW,IREP,ICHECK,IPOS(16),IBLKSZ(16)
      IPOSN=IPOS(IUNIT)
      CALL FMOVEC(C,CUFF,NWORD)
      WRITE(IUNIT,REC=IPOSN,ERR=999) (CUFF(I),I=1,127),NWORD
      IPOS(IUNIT)=IPOSN+1
      RETURN
999   write(10,'(a/3(a,i4),a/)')
     &      ' *** Output error to direct access file ***',
     &      ' *** UNIT =',iunit,'; REC =',iposn,'; NWORD =',nword,' ***'
      CALL HONDER(73,'PUTC')
      END

      SUBROUTINE SECGET (MPOS,MTYPE,IBLOCK)
* Caracterizes the type MTYPE and block number IBLOCK of section MPOS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL REVISE
      CHARACTER*6 SUBER /'SECGET'/
      COMMON/SECTOR/NUM3,IBLK3,REVISE,
     1ioda(24),ICLASS(24),ILEN(24),MAXB,KBLKLA
      PARAMETER (M1=1)
10    IF(MPOS.LT.M1)CALL HONDER(63,SUBER)
      IF(ioda(MPOS).LT.M1)CALL HONDER(64,SUBER)
      IF(MTYPE.LT.M1)MTYPE=ICLASS(MPOS)
      IF(ICLASS(MPOS).NE.MTYPE)CALL HONDER(65,SUBER)
20    IBLOCK=ioda(MPOS)+IBLK3
      RETURN
      END

      SUBROUTINE READQ(Q,NWORD,IBLK,NUM3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(*)

      CALL SEARCH(IBLK,NUM3)
      J=1
10    CONTINUE
      CALL GET(Q(J),L,NUM3)
      J=J+L
      IF(J.LE.NWORD)GO TO 10
      RETURN
      END

      SUBROUTINE READS(Q,NWORD,NUM3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(*)

      J=1
10    CONTINUE
      CALL GET(Q(J),L,NUM3)
      J=J+L
      IF(J.LE.NWORD)GO TO 10
      RETURN
      END

      SUBROUTINE READM(Q,N,M,NMAX,IBLK,NUM3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(NMAX,M)
      COMMON/BLKIN/A(511)
      DATA M1,M511/1,511/

      CALL SEARCH(IBLK,NUM3)
      CALL GET(A,NW,NUM3)
      K=M1
      DO 2 J=1,M
        DO 2 I=1,N
          Q(I,J)=A(K)
          K=K+M1
          IF(K.LE.M511) GOTO 2
          K=M1
          CALL GET(A,NW,NUM3)
2         CONTINUE
      RETURN
      END

      SUBROUTINE GET(C,NWORD,IUNIT)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BUFIN/BUFF(511)
      DIMENSION C(*)
      COMMON/DISC/ISEL,ISELR,ISELW,IREP,ICHECK,IPOS(16),IBLKSZ(16)
      IPOSN=IPOS(IUNIT)
      READ(IUNIT,REC=IPOSN) (BUFF(I),I=1,511),NWORD
      CALL FMOVE(BUFF,C,NWORD)
      IPOS(IUNIT)=IPOSN+1
      RETURN
      END

C   EXTRA SUBROUTINES EXPLICITLY FOR CHARACTER VARIABLES
C   S. M. COLWELL JAN '86

      SUBROUTINE READQC(Q,NWORD,IBLK,NUM3)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 Q
      DIMENSION Q(*)
      CALL SEARCH(IBLK,NUM3)
      J=1
10    CONTINUE
      CALL GETC(Q(J),L,NUM3)
      J=J+L
      IF(J.LE.NWORD)GO TO 10
      RETURN
      END

      SUBROUTINE GETC(C,NWORD,IUNIT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 C,CUFF
      COMMON/CUFIN/CUFF(127)
      DIMENSION C(*)
      COMMON/DISC/ISEL,ISELR,ISELW,IREP,ICHECK,IPOS(16),IBLKSZ(16)
      IPOSN=IPOS(IUNIT)
      READ(IUNIT,REC=IPOSN) (CUFF(I),I=1,127),NWORD
      CALL FMOVEC(CUFF,C,NWORD)
      IPOS(IUNIT)=IPOSN+1
      RETURN
      END

C     FOUT  M2

      SUBROUTINE FOUT(F,NUMSCF,IW)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'dimmm'
      DIMENSION F(*)
      COMMON/MAPPER/IA(NBF1+1)
      DIMENSION N(8)

10    FORMAT(/)
20    FORMAT(/5X,5(6X,I3,6X))
30    FORMAT(I5,5(F15.8))
      MAX=5
      IMAX=0
40    IMIN=IMAX+1
      IMAX=IMAX+MAX
      IF(IMAX.GT.NUMSCF) IMAX=NUMSCF
      WRITE(IW,20) (I,I=IMIN,IMAX)
      WRITE(IW,10)
      DO 60 J=1,NUMSCF
        NN=0
        DO 50 I=IMIN,IMAX
          NN=NN+1
50        N(NN)=IA(MAX0(I,J))+MIN0(I,J)
        WRITE(IW,30) J,(F(N(K)),K=1,NN)
60      CONTINUE
      IF(IMAX.LT.NUMSCF) GO TO 40
      RETURN
      END

      SUBROUTINE WHTPS
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*4 NAM
      COMMON /STREAM/ IR,IW
      COMMON/DISC/ISEL,ISELR ,ISELW,IREP,ICHEK,IPOS(16),IBLKSZ(16)
      COMMON/DISCB/NAM(16)

C...  PART OF THE BASIC ATMOL ROUTINE SYSTEM

      DO 20 K=1,16
      IF(IPOS(K))20,30,30
   20 CONTINUE
      RETURN
   30 WRITE(IW,40)
   40 FORMAT(/12X,'POSITIONED'/' F.EXT.      AT BLOCK')
      DO 60 I=K,16
      IF(IPOS(I))60,50,50
   50 WRITE(IW,70)NAM(I),IPOS(I)
   60 CONTINUE
      WRITE(IW,70)
   70 FORMAT(1X,A4,8X,I5)
      RETURN
      END

C  INDXA  A1

      SUBROUTINE INDXA(IJX,IJY,IJZ,IJ,MINI,MAXI,
     &MINJ,MAXJ,IANDJ,INC1,INC2,INC3)
      LOGICAL IANDJ
      DIMENSION IJX(*),IJY(*),IJZ(*)
     &         ,IX(10),IY(10),IZ(10)
      DATA IX/0,1,0,0,2,0,0,1,1,0/
      DATA IY/0,0,1,0,0,2,0,1,0,1/
      DATA IZ/0,0,0,1,0,0,2,0,1,1/

      IJ=0
      JMAX=MAXJ
      DO 1 I=MINI,MAXI
        NX=IX(I)*INC1+INC3
        NY=IY(I)*INC1+INC3
        NZ=IZ(I)*INC1+INC3
        IF(IANDJ) JMAX=I
        DO 1 J=MINJ,JMAX
          IJ=IJ+1
          IJX(IJ)=NX+IX(J)*INC2
          IJY(IJ)=NY+IY(J)*INC2
          IJZ(IJ)=NZ+IZ(J)*INC2
1         CONTINUE
      RETURN
      END

      SUBROUTINE DENFAN(DIJ,CSI,CPI,CDI,CSJ,CPJ,CDJ,
     1    MINI,MAXI,MINJ,MAXJ,IANDJ,DOUBLE,NORM)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL IANDJ,DOUBLE,NORM
      DIMENSION DIJ(*)
      parameter (SQRT3=1.732050807568877d0)

      MAX=MAXJ
      NN=0
      DO 140 I=MINI,MAXI
        GOTO (10,20,60,60,30,60,60,40,60,60),I
10        DUM1=CSI
        GO TO 60
20        DUM1=CPI
        GO TO 60
30        DUM1=CDI
        GO TO 60
40        IF(NORM) DUM1=DUM1*SQRT3
        GO TO 60
60      IF(IANDJ) MAX=I
        DO 140 J=MINJ,MAX
          GO TO (70,90,130,130,100,130,130,110,130,130),J
70          DUM2=DUM1*CSJ
            IF(.NOT.DOUBLE) GO TO 130
            IF(I.GT.1) GO TO 80
            DUM2=DUM2+DUM2
          GO TO 130
80          DUM2=DUM2+CSI*CPJ
          GO TO 130
90          DUM2=DUM1*CPJ
            IF(DOUBLE) DUM2=DUM2+DUM2
          GO TO 130
100         DUM2=DUM1*CDJ
            IF(DOUBLE) DUM2=DUM2+DUM2
          GO TO 130
110         IF(NORM)   DUM2=DUM2*SQRT3
          GO TO 130
130       NN=NN+1
140       DIJ(NN)=DUM2
      RETURN
      END

      SUBROUTINE STVINT
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/BIG/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     &           NI,NJ,
     &           CX,CY,CZ
      COMMON/HERMIT/H(28)
      COMMON/WERMIT/W(28)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      INTEGER*2 IMIN, IMAX, I
      DIMENSION MIN(7),MAX(7)
      DATA MIN /1,2,4, 7,11,16,22/
      DATA MAX /1,3,6,10,15,21,28/

      XINT=ZERO
      YINT=ZERO
      ZINT=ZERO
      NPTS=(NI+NJ)/2
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 13 I=IMIN,IMAX
        PX=ONE
        PY=ONE
        PZ=ONE
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
        GO TO (5,15,14,4,3,2,1),NI
    1     PX=   AX
          PY=   AY
          PZ=   AZ
    2       PX=PX*AX
            PY=PY*AY
            PZ=PZ*AZ
    3         PX=PX*AX
              PY=PY*AY
              PZ=PZ*AZ
    4           PX=PX*AX
                PY=PY*AY
                PZ=PZ*AZ
14                PX=PX*AX
                  PY=PY*AY
                  PZ=PZ*AZ
15                  PX=PX*AX
                    PY=PY*AY
                    PZ=PZ*AZ
    5                 GO TO (12,11,10,9,8,7,6),NJ
    6                   PX=PX*BX
                        PY=PY*BY
                        PZ=PZ*BZ
    7                     PX=PX*BX
                          PY=PY*BY
                          PZ=PZ*BZ
    8                       PX=PX*BX
                            PY=PY*BY
                            PZ=PZ*BZ
    9                         PX=PX*BX
                              PY=PY*BY
                              PZ=PZ*BZ
   10                           PX=PX*BX
                                PY=PY*BY
                                PZ=PZ*BZ
   11                             PX=PX*BX
                                  PY=PY*BY
                                  PZ=PZ*BZ
   12                               DUM=W(I)
                                    XINT=XINT+DUM*PX
                                    YINT=YINT+DUM*PY
                                    ZINT=ZINT+DUM*PZ
   13                               CONTINUE
      RETURN
      END

      SUBROUTINE ROOTS
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      LOGICAL JUMP
      DIMENSION F(NSH1*2)
      COMMON/RTDATA/WHI(15),WLOW(15),RHI(15),RLOW(15),RFAC(15),
     &   AMPS(9),IPOINT(9),MADD(20)
      DIMENSION  ROOT(NATM1),WT(NATM1),
     &           BETA(NATM1),GAMMA(NATM1),V(NATM1),D(NATM1)
      COMMON/RTWT/X,ROOT,WT,NPT,NPT2,DJI(NSH1*2)
      PARAMETER (EPS=1.0D-12, ALIM=1.0D-11, ZERO=0.D0, ONE=1.D0)
      INTEGER*2 M1,M2,IFLAG,IPT,IFPT,MM,M24,IBAS,JBAS,K,IBMM,JBMM,
     &          KM1,KM2,IFLAGK,IFLAGT, I
      PARAMETER (M1=1, M2=2, M24=24)
      DATA MS2/-2/
      DATA D(1)/1.0D0/

      NPT2=NPT-2
      MPRIME=NPT+NPT-1
      MPRIM2=MPRIME-1
      IFLAG=IPOINT(NPT)
      XMAX=AMPS(NPT)
      IF(X.LE.ALIM) THEN
        DO 20 IPT=1,NPT
          IFPT = IFLAG+IPT
          ROOT(IPT)=RLOW(IFPT)
20        WT(IPT)=WLOW(IFPT)
      ELSE
        IF(X.LT.XMAX) THEN
40        CALL FFUN(X,F,NPT,DJI,MADD)
          A=F(2)/F(1)
          IF(NPT.LE.M1) THEN
            ROOT(1)=A
            WT(1)=F(1)
            GO TO 190
          ENDIF
          A1=A
          DO 60 MM=1,MPRIM2
60          F(MM+M24)=F(MM+M2)-A*F(MM+M1)
          IBAS=0
          JBAS=M24
          MAX=MPRIM2
          DO 80 K=2,NPT
            CFAC=F(JBAS+M1)
            AX=F(JBAS+M2)/CFAC
            BET=A-AX
            KM1 = K - M1
            BETA(KM1)=BET
            A=AX
            GAM=-CFAC/F(IBAS+M1)
            GAMMA(KM1)=GAM
            IF(K.EQ.NPT) GO TO 90
            MAX=MAX+MS2
            DO 70 MM=1,MAX
              IBMM = IBAS+MM
              JBMM = JBAS+MM
70            F(IBMM)=BET*F(JBMM+M1)+GAM*F(IBMM+M2)+F(JBMM+M2)
            IBAS=JBAS
80          JBAS=M24-JBAS
90        GAM=ZERO
          DO 100 K=1,NPT
            IFLAGK = IFLAG+K
            R=RLOW(IFLAGK)/(RFAC(IFLAGK)*X+ONE)
            IF(X.GT.17.0D0)R=RHI(IFLAGK)/X
            ROOT(K)=R
100         GAM=GAM+R
          EPSI=EPS*AX
          DO 160 MM=1,NPT
            RRR=ROOT(MM)
            R=RRR*AX/GAM
            JUMP=.FALSE.
            DO 140 ICOUNT=1,50
              V(1)=R-A1
              TOP=R+BETA(1)
              DD=V(1)+TOP
              D(2)=DD
              VV=V(1)*TOP+GAMMA(1)
              IF(NPT2)130,130,110
110             V(2)=VV
                DO 120 K=3,NPT
                  KM1 = K-M1
                  KM2 = K-M2
                  TOP=R+BETA(KM1)
                  DD=TOP*DD+GAMMA(KM1)*D(KM2)+VV
                  D(K)=DD
                  VV=TOP*VV+GAMMA(KM1)*V(KM2)
120               V(K)=VV
130           SCALE=VV/DD
              R=R-SCALE
              IF(JUMP) GOTO 150
              IF(DABS(SCALE).LT.EPSI) JUMP=.TRUE.
140           CONTINUE
C      WRITE(9,123)X,NPT
123   FORMAT(1X,'NO CONVERGENCE',F15.6,I6)
150         ROOT(MM)=R
            WT(MM)=CFAC/(DD*V(NPT-M1))
            AX=AX-R
160         GAM=GAM-RRR
        ELSE
          TOP=ONE/X
          SCALE=DSQRT(TOP)
          DO 180 IPT=1,NPT
            IFLAGT = IFLAG+IPT
            ROOT(IPT)=RHI(IFLAGT)*TOP
180         WT(IPT)=WHI(IFLAGT)*SCALE
        ENDIF
      ENDIF
190   CONTINUE
      DO 200 I=1,NPT
        RR=ROOT(I)
200     ROOT(I)=RR/(ONE-RR)
      RETURN
      END

      SUBROUTINE FFUN(X,F,NPT,DJI,MADD)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*2 (I-N)
      INTEGER*4 NPT,MADD
      DIMENSION F(*),MADD(*),DJI(*)
      DIMENSION XMAX(20)
      PARAMETER (ZERO=0.D0)
      DATA XMAX/24.0D0,29.0D0,32.0D0,35.0D0,37.0D0,40.0D0,42.0D0,
     &   44.0D0,47.0D0,49.0D0,51.0D0,53.0D0,55.0D0,57.0D0,58.0D0,
     &   59.0D0,59.0D0,59.0D0,59.0D0,59.0D0/

      M=NPT+NPT-1
      FACT=4.8D0+0.4D0*NPT
      IF(X.LT.FACT) THEN
        IF((X-1.0D-10).LT.ZERO) THEN
          J=M+1
          DO 140 I=1,J
140         F(I)=DJI(I)
          RETURN
        ELSE
          IF((X-1.0D-5).LT.ZERO) THEN
            J=M+1
            DO 170 I=1,J
170           F(I)=DJI(I)-X*DJI(I+1)
            RETURN
          ELSE
            X2=X+X
            I=X2+MADD(M)
            J=I
            X3=DEXP(-X)
            F(J+1)=ZERO
            DO 190 K=1,I
              F(J)=(F(J+1)*X2+X3)*DJI(J)
190           J=J-1
          ENDIF
        ENDIF
        RETURN
      ENDIF
      X2=0.5D0/X
      XINV=X2+X2
      F1=DSQRT(1.5707963267949D0*X2)
      IF (X.GE.XMAX(M)) THEN
        IF (M.LE.0) RETURN
        F(1)=F1
        DO 110 I=1,M
          F(I+1)=F(I)*X2
110       X2=X2+XINV
        RETURN
      ENDIF
      X3=DEXP(-X)
      IF(X.GT.21.6D0) THEN
        F1=F1-X3*X2
        GOTO 50
      ENDIF
      IF(X.GT.18.2D0) THEN
        F1=F1-X3*X2*(1.0D0-X2)
        GO TO 50
      ENDIF
      IF(X.GT.(12.0D0+0.1D0*NPT)) THEN
        F1= (( 1.9623264149430D-1*XINV-4.9695241464490D-1)*XINV
     1                                -6.0156581186481D-5)*X3 + F1
        GOTO 50
      ENDIF
      IF(X.GT.(9.2D0+0.2D0*NPT)) THEN
        F1= (((-1.8784686463512D-1*XINV+2.2991849164985D-1)*XINV
     1        -4.9893752514047D-1)*XINV-2.1916512131607D-5)*X3
     2      + F1
        GO TO 50
      ENDIF
      F1=F1+((((((4.6897511375022D-1*XINV-6.9955602298985D-1)*XINV
     1+5.3689283271887D-1)*XINV-3.2883030418398D-1)*XINV
     2+2.4645596956002D-1)*XINV-4.9984072848436D-1)*XINV
     3-3.1501078774085D-6)*X3
50    F(1)=F1
      IF (M.LE.0) RETURN
      X23=X2*X3
      DO 70 I=1,M
        F(I+1)=F(I)*X2-X23
70      X2=X2+XINV
      RETURN

      END

      SUBROUTINE DENFAC(DKL,CSK,CPK,CDK,CSL,CPL,CDL,MINK,
     & MAXK,MINL,MAXL,KANDL,DOUBLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL KANDL,DOUBLE
      DIMENSION DKL(*)

      N=0
      MAX=MAXL
      DO 100 K=MINK,MAXK
        IF(KANDL)MAX=K
        GO TO(10,20,40,40,30,40,40,40,40,40),K
10        DUM1=CSK
        GO TO 40
20        DUM1=CPK
        GO TO 40
30        DUM1=CDK
40      DO 100 L=MINL,MAX
          GO TO(50,70,90,90,80,90,90,90,90,90),L
50          DUM2=DUM1*CSL
            IF(.NOT.DOUBLE) GO TO 90
            IF(K.GT.1) GO TO 60
            DUM2=DUM2+DUM2
            GO TO 90
60          DUM2=DUM2+CSK*CPL
          GO TO 90
70          DUM2=DUM1*CPL
            IF(DOUBLE) DUM2=DUM2+DUM2
          GO TO 90
80          DUM2=DUM1*CDL
            IF(DOUBLE) DUM2=DUM2+DUM2
90        N=N+1
100       DKL(N)=DUM2
      RETURN
      END

C  LOCAL  A1
C     COMPILE WITH AUTODBL
      SUBROUTINE LOCAL(X,Y,Z,XS,YS,ZS)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 GROUP
      COMMON/SYMTRY/T(216),INVT(24),NEWSH(NSH1+1,24),NT,NT2,
     &XOLD,YOLD,ZOLD,XNEW,YNEW,ZNEW,XP,YP,ZP,
     &U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
     &,INDEX,NAXIS
      COMMON/SYMTRB/GROUP

      XX0 = X-X0
      YY0 = Y-Y0
      ZZ0 = Z-Z0
      XS=U1*(XX0)+U2*(YY0)+U3*(ZZ0)
      YS=V1*(XX0)+V2*(YY0)+V3*(ZZ0)
      ZS=W1*(XX0)+W2*(YY0)+W3*(ZZ0)
      RETURN
      END

C  TRANS  A1
      SUBROUTINE TRANS(NN)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      CHARACTER*8 GROUP
      COMMON/SYMTRY/T(216),INVT(24),NEWSH(NSH1+1,24),NT,NT2,
     &XOLD,YOLD,ZOLD,XNEW,YNEW,ZNEW,XP,YP,ZP,
     &U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
     &,INDEX,NAXIS
      COMMON/SYMTRB/GROUP

      XNEW=XOLD*T(NN+1)+YOLD*T(NN+2)+ZOLD*T(NN+3)
      YNEW=XOLD*T(NN+4)+YOLD*T(NN+5)+ZOLD*T(NN+6)
      ZNEW=XOLD*T(NN+7)+YOLD*T(NN+8)+ZOLD*T(NN+9)
      XP=X0+U1*XNEW+V1*YNEW+W1*ZNEW
      YP=Y0+U2*XNEW+V2*YNEW+W2*ZNEW
      ZP=Z0+U3*XNEW+V3*YNEW+W3*ZNEW
      RETURN
      END

      SUBROUTINE STRUCT(N, C, AZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      INTEGER P
      LOGICAL  LADAPT,LSKIP,REVISE
      DIMENSION COV(36),C(3,NATM1),AZ(NATM1)
      COMMON/SMALL/BL(NATM1,NATM1),V(8),M(4,8)

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
      COMMON/SECTOR/NUM3,IBLK3,REVISE,IODA(24),ICLASS(24),
     &              ILEN(24),MAXB,KBLKLA
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0)
      PARAMETER (PI=3.14159265358979D0)

*  Covalent radii (from "Table of Periodic Properties of the Elements",
*  Sargent-Welch Sci. Corp., 1968) for atoms up to argon (in a.u.).

      DATA COV /0.605D0, 1.758D0,
     1 2.33D0,1.70D0,1.55D0,1.46D0,1.42D0,1.38D0,1.36D0,1.34D0,
     2 2.91D0,2.57D0,2.23D0,2.10D0,2.00D0,1.93D0,1.87D0,1.85D0,
     3 3.84D0,3.29D0,2.72D0,2.49D0,2.31D0,2.23D0,2.21D0,2.21D0,
     4 2.19D0,2.17D0,2.21D0,2.36D0,2.38D0,2.31D0,2.27D0,2.19D0,
     5 2.16D0,2.12D0/

      NM1=N-1
      IF(NM1.LE.0)RETURN

* SET BOND LENGTH MATRIX <BL> AND MAKE IT NEGATIVE WHEN NO BOND OCCURS

      DO 10 I=1,NM1
        IZ=IDINT(AZ(I)+0.1D0)
        COVI=TWO
        IF (IZ .GT. 0 .AND. IZ .LE. 36) COVI=COV(IZ)
        BL(I,I)=ZERO
        DO 10 J=I+1,N
          BL(I,J)=DSQRT(
     1    (C(1,J)-C(1,I))**2 + (C(2,J)-C(2,I))**2 + (C(3,J)-C(3,I))**2)
          JZ=IDINT(AZ(J)+0.1D0)
          COVJ=TWO
          IF (JZ .GT. 0 .AND. JZ .LE. 36) COVJ=COV(JZ)
          IF (BL(I,J) .GT. (COVI+COVJ)*1.2) BL(I,J)=-BL(I,J)
          BL(J,I)=BL(I,J)
10        CONTINUE

      BL(N,N)=ZERO

* PRINT BOND LENGTHS

      WRITE (IW,20)
20    FORMAT (/' BOND LENGTHS IN BOHR [ANGSTROM]')

      MAX=5
      P=0
      DO 40 I=1,NM1
        DO 40 J=I+1,N
          IF (BL(I,J) .LE. ZERO) GO TO 40
          P=P+1
          IF (P .GT. MAX) THEN
            CALL PRTBUF(MAX,2, M, V)
            P=1
          ENDIF
          M(1,P)=I
          M(2,P)=J
          V(P)=BL(I,J)
40        CONTINUE

      IF (P .GT. 0) CALL PRTBUF(P,2, M, V)
      IF (P .EQ. 0) WRITE (IW,50)
50    FORMAT ('      NONE')

      WRITE (IW,60)
60    FORMAT (/' BOND ANGLES')

      MAX=3
      P=0
      DO 100 J=1,NM1
        DO 90 I=1,N
          IF (BL(I,J) .LE. ZERO .OR. J .EQ. I) GO TO 90
          JP1=J+1
          DO 80 K=JP1,N
            IF (BL(I,K) .LE. ZERO .OR. K .EQ. I) GO TO 80
            P=P+1
            IF (P .LE. MAX) GO TO 70
            CALL PRTBUF(MAX,3, M, V)
            P=1
70          M(1,P)=J
            M(2,P)=I
            M(3,P)=K
            COSINE=
     &      (BL(I,J)**2+BL(I,K)**2-BL(J,K)**2)/(TWO*BL(I,J)*BL(I,K))
            IF(DABS(COSINE).GT.ONE) COSINE=DSIGN(ONE,COSINE)
            V(P)=(180.0D0/PI)*DACOS(COSINE)
80          CONTINUE
90        CONTINUE
100     CONTINUE

      IF (P .GT. 0) CALL PRTBUF(P,3, M, V)
      IF (P .EQ. 0) WRITE (IW,50)
      WRITE (IW,110)
110   FORMAT (/' DIHEDRAL ANGLES')

      MAX=3
      P=0
      DO 160 I=1,NM1
        IP1=I+1
        DO 150 J=IP1,N
          IF (BL(I,J) .LE. ZERO) GO TO 150
C  I and J are bonded.  Construct unit vector from I to J.
          A=BL(I,J)
          XIJ=(C(1,J)-C(1,I))/A
          YIJ=(C(2,J)-C(2,I))/A
          ZIJ=(C(3,J)-C(3,I))/A
          DO 140 K=1,N
            IF (BL(I,K).LE.ZERO .OR. K.EQ.I .OR. K.EQ.J) GO TO 140
C  I and K are bonded.  Construct unit vector in JIK plane, normal
C  to IJ.
            XIK=C(1,K)-C(1,I)
            YIK=C(2,K)-C(2,I)
            ZIK=C(3,K)-C(3,I)
            DOT=XIK*XIJ + YIK*YIJ + ZIK*ZIJ
            XIK=XIK-DOT*XIJ
            YIK=YIK-DOT*YIJ
            ZIK=ZIK-DOT*ZIJ
            A=DSQRT(XIK**2+YIK**2+ZIK**2)
            IF (A .LT. 1.0D-4) GO TO 140
            XIK=XIK/A
            YIK=YIK/A
            ZIK=ZIK/A
            DO 130 L=1,N
              IF (BL(J,L).LE.ZERO .OR. L.EQ.I .OR. L.EQ.K) GO TO 130
C  J and L are bonded.  Construct unit vector in IJL plane, normal
C  to IJ.
              XJL=C(1,L)-C(1,J)
              YJL=C(2,L)-C(2,J)
              ZJL=C(3,L)-C(3,J)
              DOT=XJL*XIJ + YJL*YIJ + ZJL*ZIJ
              XJL=XJL-DOT*XIJ
              YJL=YJL-DOT*YIJ
              ZJL=ZJL-DOT*ZIJ
              A=DSQRT(XJL**2+YJL**2+ZJL**2)
              IF (A .LT. 1.0D-4) GO TO 140
              XJL=XJL/A
              YJL=YJL/A
              ZJL=ZJL/A
              P=P+1
              IF (P .LE. MAX) GO TO 120
              CALL PRTBUF(MAX,4, M, V)
              P=1
120           M(1,P)=K
              M(2,P)=I
              M(3,P)=J
              M(4,P)=L
              COSINE=XIK*XJL+YIK*YJL+ZIK*ZJL
              SINE=XIJ*(YJL*ZIK-ZJL*YIK)
     1             + YIJ*(ZJL*XIK-XJL*ZIK) + ZIJ*(XJL*YIK-YJL*XIK)
              V(P)=(180.0D0/PI)*DATAN2(SINE,COSINE)
130           CONTINUE
140         CONTINUE
150       CONTINUE
160     CONTINUE

      IF (P .GT. 0) CALL PRTBUF(P,4, M, V)
      IF (P .EQ. 0) WRITE (IW,50)

      WRITE (IW,170)
170   FORMAT (/)

      RETURN
      END

      SUBROUTINE PRTBUF(N,K, M, V)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL  LADAPT,LSKIP
      CHARACTER*1 BRA,KET
      LOGICAL REVISE
      DIMENSION V(8),M(4,8)

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
      DATA BRA,KET/'[',']'/

      DO 10 I=1,N
        DO 10 J=2,K
10        IF (M(J,I) .LT. 10) M(J,I)=-M(J,I)

      GO TO (20,20,60,80), K

20    WRITE (IW,30) (M(1,I), M(2,I), V(I), I=1,N)
30    FORMAT (/5(I2, '-', I2, F9.5, 2X))
      DO 40 I=1,N
40      V(I)=V(I)*0.529177D0
      WRITE (IW,50) (BRA, V(I), KET, I=1,N)
50    FORMAT (5(5X, A1, F7.5, A1,2X))
      RETURN

60    WRITE (IW,70) (M(1,I), M(2,I), M(3,I), V(I), I=1,N)
70    FORMAT(/1X,5(I2,2('-',I2),F12.6,2X),I2,2('-',I2),F12.6)
      RETURN

80    WRITE (IW,90) (M(1,I), M(2,I), M(3,I), M(4,I), V(I), I=1,N)
90    FORMAT (/1X, 5(I2, 3('-', I2), F12.6, 2X))
      RETURN
      END

      FUNCTION LENINT(N)
C     LENINT is the number of REALS corresponding to N INTEGERS
      LENINT=(N+1)/2
CCRAY LENINT=N
      RETURN
      END

      FUNCTION VECSUM(A,B,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*)
      VECSUM=0.0D0
      DO 10 I=1,N
        VECSUM=VECSUM+A(I)*B(I)
10      CONTINUE
      RETURN
      END

C  VECOUT  A1

      SUBROUTINE VECOUT(V,E,N,M,NDIM,IW,BFLAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 BFLAB
      DIMENSION BFLAB(*)
      DIMENSION V(NDIM,M),E(*)

      MAX=0
   10 MIN=MAX+1
      MAX=MAX+5
      IF(MAX.GT.M) MAX=M
      WRITE(IW,9999)
      WRITE(IW,9996) (E(J),J=MIN,MAX)
      WRITE(IW,9999)
      WRITE(IW,9998) (J,J=MIN,MAX)
      WRITE(IW,9999)
      DO 20 I=1,N
20    WRITE(IW,9997) I,BFLAB(I),(V(I,J),J=MIN,MAX)
      IF(MAX.LT.M) GO TO 10
9999  FORMAT(/)
9998  FORMAT(15X,5(5X,I3,5X))
9997  FORMAT(I5,2X,A8,7F13.8)
9996  FORMAT(15X,5F13.8)
      RETURN
      END

C     PRINT  M2

      SUBROUTINE PRINT(N1,N2,IW)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      COMMON/SYMTRY/T(216),INVT(24),ISO(NSH1+1,24),NT
     &,NT2,XX(21),INDEX,NAXIS
      DIMENSION NN(48)

      IMAX=N1-1
10    IMIN=IMAX+1
      IMAX=IMAX+4
      IF(IMAX.GT.N2) IMAX=N2
      NJ=9*N1-8
      DO 30 J=1,3
      NI=0
      DO 20 I=IMIN,IMAX
      NN(I)=NJ+NI
20    NI=NI+9
      WRITE(IW,40) (T(NN(I)),T(NN(I)+1),T(NN(I)+2),I=IMIN,IMAX)
30    NJ=NJ+3
      WRITE(IW,50)
      IF(IMAX.LT.N2) GO TO 10
40    FORMAT(4X,4(3F10.5,2H *))
50    FORMAT(/)
      RETURN
      END

      FUNCTION LENSEC(N)
      LENSEC=(N-1)/511+1
      RETURN
      END

      FUNCTION NRC(N)
      NRC=(N-1)/127+1
      RETURN
      END

      SUBROUTINE DTR (DMO, DAO, COORB, T,NBASIS,NCOORB,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)

C  Transform density matrix from MO basis to AO basis

C  DMO density matrix in MO basis (triangular form)
C  DAO array for density matrix in AO basis (triangular form)
C  COORB MO coefficients
C  T workspace
C  NBASIS number of basis functions

      DIMENSION DMO(*), DAO(*), COORB(NDIM,*), T(*)

      IJ=0
      DO 70 J=1,NBASIS
        T(1)=DMO(1)*COORB(J,1)
        KK=1
        DO 40 K=2,NCOORB
          X=0.D0
          KM1=K-1
          DO 30 L=1,KM1
            X=X+DMO(KK+L)*COORB(J,L)
30          T(L)=T(L)+DMO(KK+L)*COORB(J,K)
          T(K)=X+DMO(KK+K)*COORB(J,K)
40        KK=KK+K
        DO 60 I=1,J
          X=0.D0
          DO 50 K=1,NCOORB
50          X=X+COORB(I,K)*T(K)
          IJ=IJ+1
60        DAO(IJ)=X
70      CONTINUE

      RETURN
      END

      SUBROUTINE ONEPDM(DA,DB)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      INTEGER P1,P2,P3
      CHARACTER*8 ANAM,ATNAME,BFLAB,PARNAM
      COMMON /INFOA/ NAT, ICH, MUL, NUM, NX, NE, NA, NB,
     1    ZAN(NATM1), C(3,NATM1),
     2    NVAR, NPAR, NVAL, NUQ, CFR, CFA, KONTYP(NATM1),
     3    P1(NATM1), P2(NATM1), P3(NATM1),
     4    KR(NATM1), KA(NATM1), KB(NATM1),
     5    PARVAL(NATM13), ZIN(NATM1), X(NATM13)
      COMMON/INFOB/ANAM(NATM1),BFLAB(NBF1),PARNAM(NATM13),ATNAME(NATM1)
      LOGICAL LADAPT,LSKIP
      CHARACTER*8 CLOSED
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
      DIMENSION DA(*),DB(*)
      DATA CLOSED/'CLOSED'/

      MTYP=0
      MTYPO=0
      IBLOK=0
      IBLOKO=0
      CALL SECGET(ISEX(7),MTYP,IBLOK)
      CALL READQ(DA,NX,IBLOK,IFILD)
      IF(SCFTYP.EQ.CLOSED) RETURN
      CALL SECGET(ISEX(10),MTYPO,IBLOKO)
      CALL READQ(DB,NX,IBLOKO,IFILD)
      DO 100 I=1,NX
100   DA(I)=DA(I)+DB(I)
      END

C  DVINT  A1
C     MUST BE COMPILED WITH AUTODBL

      SUBROUTINE DVINT

C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----

      IMPLICIT REAL*8(A-H,O-Z)
      include 'dimmm'
      COMMON/BIG/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     &           NI,NJ,
     &           CX,CY,CZ
      COMMON/HERMIT/H1,H2(2),H3(3),H4(4),H5(5),H6(6),H7(7)
      COMMON/WERMIT/W1,W2(2),W3(3),W4(4),W5(5),W6(6),W7(7)
      DIMENSION H(21),W(21),MIN(6),MAX(6)
      EQUIVALENCE (H(1),H1),(W(1),W1)
      PARAMETER (ZERO=0.0D0)
      INTEGER*2 I,NPTS
      DATA MIN /1,2,4,7,11,16/
      DATA MAX /1,3,6,10,15,21/

      XINT=ZERO
      YINT=ZERO
      ZINT=ZERO
      NPTS=(NI+NJ-1)/2+1
      IMIN=MIN(NPTS)
      IMAX=MAX(NPTS)
      DO 13 I=IMIN,IMAX
        DUM=H(I)/T
        PTX=DUM+X0
        PTY=DUM+Y0
        PTZ=DUM+Z0
        PX=PTX-CX
        PY=PTY-CY
        PZ=PTZ-CZ
        AX=PTX-XI
        AY=PTY-YI
        AZ=PTZ-ZI
        BX=PTX-XJ
        BY=PTY-YJ
        BZ=PTZ-ZJ

        GO TO (5,4,3,2,1),NI
    1     PX=PX*AX
          PY=PY*AY
          PZ=PZ*AZ
    2       PX=PX*AX
            PY=PY*AY
            PZ=PZ*AZ
    3         PX=PX*AX
              PY=PY*AY
              PZ=PZ*AZ
    4           PX=PX*AX
                PY=PY*AY
                PZ=PZ*AZ
    5           GO TO (12,11,10,9,8,7,6),NJ
    6             PX=PX*BX
                  PY=PY*BY
                  PZ=PZ*BZ
    7               PX=PX*BX
                    PY=PY*BY
                    PZ=PZ*BZ
    8                 PX=PX*BX
                      PY=PY*BY
                      PZ=PZ*BZ
    9                   PX=PX*BX
                        PY=PY*BY
                        PZ=PZ*BZ
   10                     PX=PX*BX
                          PY=PY*BY
                          PZ=PZ*BZ
   11                       PX=PX*BX
                            PY=PY*BY
                            PZ=PZ*BZ
   12                         DUM=W(I)
                              XINT=XINT+DUM*PX
                              YINT=YINT+DUM*PY
                              ZINT=ZINT+DUM*PZ
   13   CONTINUE
      RETURN
      END

      SUBROUTINE READC

      CHARACTER*1 IA
      INTEGER*2 I
      LOGICAL INITEM
      COMMON/STREAM/IR,IW
      COMMON/WORK/ITEM,NITEMS,ISTRT(40),INUMB(40),IWIDTH
     ./answer/ si,bl,hyphen,star,plus,comma,equal
      character si,bl,hyphen,star,plus,comma,equal
      COMMON/WORKB/IA(80)

c NITEMS is the amount of words in the line (JUMP in other subroutines)
c ISTRT(i) is the column where word 'i' begins
c INUMB(i) is the amount of characters in word 'i'
c IWIDTH was set to 80 columns in BLOCK DATA

      NITEMS=0
      ITEM=0
      INITEM=.FALSE.
      call upcase (ir,ia)
      DO 40 I=1,IWIDTH
        IF(IA(I).EQ.BL .or.
     &     IA(I).EQ.EQUAL .or.
     &     IA(I).EQ.COMMA) then
          INITEM=.FALSE.
        else
          IF (INITEM) then
            INUMB(NITEMS)=INUMB(NITEMS)+1
          else
            NITEMS=NITEMS+1
            ISTRT(NITEMS)=I
            INUMB(NITEMS)=1
            INITEM=.TRUE.
          endif
        endif
40      CONTINUE
      RETURN
      END

C     READF  M2

      SUBROUTINE READF(BUFF)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 ITEMP,IA
      COMMON/WORK/JREC,JUMP,ISTRT(40),INUMB(40),IWIDTH
     ./digits/ iachar(12), achar(13), end, buffa
      COMMON/WORKB/IA(80)
      character iachar, achar
      character*8 end, buffa

      BUFF=0.0D0
      JREC=JREC+1
      IF(JREC.GT.JUMP) GOTO 80
      LIMIT=13
      IFACT2=0
      N=INUMB(JREC)
      FACT=1.0D0
      IST=ISTRT(JREC)
      NSTRT=IST+N-1
      DO 60 I=1,N
        ITEMP=IA(NSTRT)
        DO 10 J=1,LIMIT
          IF(ACHAR(J).EQ.ITEMP) GOTO 30
10        CONTINUE
20      CALL HONDER(68,'READF ')
30      IF(J.LT.11) GOTO 50
        IF(J.LT.13) GOTO 40
        IFACT2=I-1
        LIMIT=12
        GOTO 60
40      IF(NSTRT.NE.IST) GOTO 20
        IF(J.EQ.12) BUFF=-BUFF
        GOTO 70
50      BUFF=BUFF+(J-1)*FACT
        FACT=FACT*10.0D0
60      NSTRT=NSTRT-1
70    BUFF=BUFF/10.D0**IFACT2

80    RETURN
      END

C     READI  M2

      SUBROUTINE READI(JBUF)
      CHARACTER*1 ITEMP,IA
      COMMON/WORK/JREC,JUMP,ISTRT(40),INUMB(40),IWIDTH
     ./digits/ iachar(12), achar(13), end, buffa
      COMMON/WORKB/IA(80)
      character iachar, achar
      character*8 end, buffa

      JBUF=0
      JREC=JREC+1
      IF(JREC.GT.JUMP) GOTO 60
      N=INUMB(JREC)
      IFACT=1
      IST=ISTRT(JREC)
      NSTRT=IST+N-1
      DO 50 I=1,N
        ITEMP=IA(NSTRT)
        DO 10 J=1,12
          IF(IACHAR(J).EQ.ITEMP) GOTO 30
10        CONTINUE
20      CALL HONDER(69,'READI ')
30      IF(J.LT.11) GOTO 40
        IF(NSTRT.NE.IST) GOTO 20
        IF(J.EQ.11) JBUF=-JBUF
        GOTO 60
40      JBUF=JBUF+(J-1)*IFACT
        IFACT=IFACT*10
50      NSTRT=NSTRT-1

60    RETURN
      END

C     READA  M2

      SUBROUTINE READA (BUFF)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 IA,IBLNK,IBUF
      CHARACTER*8 BUFF
      DIMENSION IBUF(8)
      COMMON/WORK/JREC,JUMP,ISTRT(40),INUMB(40),IWIDTH
      COMMON/WORKB/IA(80)
      DATA IBLNK/' '/

      DO 10 I=1,8
10      IBUF(I)=IBLNK
      JREC=JREC+1
      IF(JREC.LE.JUMP) THEN
        N=INUMB(JREC)
        NSTRT=ISTRT(JREC)
        IF(N.GT.8) N=8
        DO 20 I=1,N
          IBUF(I)=IA(NSTRT)
20        NSTRT=NSTRT+1
      ENDIF
30    CALL PACK8(BUFF,IBUF)
      RETURN
      END

      SUBROUTINE PACK8(BUFF,IBUF)
      CHARACTER*8 A,BUFF
      CHARACTER*1 IBUF(8),LBUF(8)
      EQUIVALENCE (A,LBUF)

      DO 10 I=1,8
10      LBUF(I)=IBUF(I)
      BUFF=A

      RETURN
      END

