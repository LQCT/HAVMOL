C  MAIN PROGRAM ---- CALLS IOFILES AND MASTER

C  INITIAL DATA IN BLOCKS DATA CINIT AND DATIN
C  S.M. COLWELL JAN 1985
C  L.A. MONTERO NOV 1988, FEB 1989, MAY 1989, OCT. 1989, FEB. 1990, MAY 1991
C  April 1992, October 1992, January 1993, June 1993

      include 'screen.inc'
      PROGRAM HAVMOL
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'dimmm'
      include 'screenp.inc'

      call clearscreen($gclearscreen)
      open (10, file='CON')
      call comlin
      CALL MACHON (10)
      CALL IOFILES (iclos)
      CALL MASTER (iclos)
      STOP '*** HAVMOL terminated ***'
      END

      subroutine comlin
      character*32 arg
      logical*2 qcl
      integer*2 marg,sta,numarg
      common /cl/ qcl,numarg,marg(5),sta,arg

      numarg = NARGS()
      qcl = numarg.gt.1
      numarg = numarg - 1
      do i=1,5
        marg(i) = i
      enddo
      return
      end

      subroutine rcomlin (n,subpr)
      integer*4 n
      character*6 subpr
      character*32 arg,blanco/'                                '/
      logical*2 qcl
      integer*2 marg,sta,numarg
      common /cl/ qcl,numarg,marg(5),sta,arg
      arg=blanco
      call getarg (marg(n),arg,sta)
      if (sta.lt.0) call honder(80,subpr)
      call upcase(-1,arg)
      return
      end

      SUBROUTINE MACHON (II)
      include 'dimmm'
      WRITE (II ,70)
70    FORMAT (/
     &37x,'HAVMOL'//
     &13x,'Ab Initio Hartree-Fock Molecular Orbital Calculations'//
     & 5x,'Version 1.1 by Luis A. Montero, University of Havana, Havana,
     & June 1993'/
     &16X,'after the IBM PC MICROMOL, Cambridge, Sept. 1987'/
     &26x,'by R.D. Amos and S.M. Colwell'/
     &)
      write (II,80) NATM1, NSH1, NBF1, NPRM1
80    FORMAT (t21,'Maximum number of atoms:          ',i4
     &/       t21,'Maximum number of shells:         ',i4
     &/       t21,'Maximum number of basis functions:',i4
     &/       t21,'Maximum number of primitives:     ',i4/)
      RETURN
      END

      SUBROUTINE IOFILES (iclos)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'screenp.inc'
C     OPENING FILES IN DOS
      character*80 kwd,anl,anl0
      character*32 file5,file6
      character*6 suber /'IOFILE'/
      common
     ./files / file5
     ./optn1/ kwd(48)
     ./optn2/ icar, kb1, iii
      common mes,idi,ian,iho,min,ise,ics
      integer*2 mes,idi,ian,iho,min,ise,ics
      character*32 arg
      logical*2 qcl
      integer*2 marg,sta,numarg
      common /cl/ qcl,numarg,marg(5),sta,arg
      character*4 ext1,ext
      character*5 dd
      character*1 an
      save /cl/

      iclos = 0
      KK = 0
      kb1 = 0
      iii = 0
      dd = 'main'
      irst = 0
      CALL GETDAT (IAN,MES,IDI)
      CALL GETTIM (IHO,MIN,ISE,ICS)
      WRITE (10,150) MES,IDI,IAN,IHO,MIN,ISE,ICS

      if (qcl) then
        call rcomlin (1,suber)
        if (arg.eq.'RESTORE') then
          an='Y'
        else
          an='N'
        endif
      else
        write(10,'(/a\)')
     &' Are you trying to restore a previous calculation? (Y/[N]) '
        call upcase (10,an)
      endif

c     restart procedure

      if (an.eq.'Y') then
        if (qcl) then
          call rcomlin (2,suber)
          if (arg.eq.'INTERNAL') then
            an='Y'
            qcl = .false.
          else
            an='N'
          endif
        else
          write(10,'(/a\)')
     &  ' Do you want to enter restart commands from the keyboard? (Y/[N
     &]) '
          call upcase (10,an)
        endif

        if (an.eq.'Y') then
          kb1 = 1
          write (kwd(1),'(a)') 'RESTORE'
          write(10,'(/a/a/a)')
     &    ' Enter up to 47 lines of HAVMOL commands (including the STA
     &RT <n> and FINISH',
     &    ' or STOP lines, RESTORE must NOT be written). SYMMETRY direct
     &ives can not be',
     &    ' changed.'
          anl0 = 'RESTORE'
          do 15 i2=2,48
            write(10,'(a\)') ' => '
            call upcase (10,anl)
            write (kwd(i2),'(a)') anl
            if (anl.eq.'STOP' .or.
     &          anl.eq.'FINISH' .or.
     &         (anl.eq.'END' .and.
     &          index(anl0,'R').eq.4 .and.
     &          index(anl0,'T').eq.2))
     &        go to 16
            anl0 = anl
15          continue
16        dd = 'dummy'
        endif

        ext1 = '.IM1'
        irst = 1
      else

c     option input files

        if (qcl) then
          call rcomlin(2,suber)
          call filen1(-1,arg,ext)
          if (ext.eq.'.CAR') then
            an = '2'
          elseif (ext.eq.'.IMC') then
            an = '3'
          elseif (ext.eq.'.AMI') then
            an = '4'
          elseif (ext.eq.'.ZMT') then
            an = '5'
          else
            an = '1'
          endif
        else
          write(10,'(/a//a/a/a/a/a/a\)')
     &' Enter the code for a kind of main input file:',
     &'   HAVMOL (.IMM) ......................... [1]',
     &'   Cartesian coordinate (CAR) ............. 2',
     &'   Internal coordinate (IMC) .............. 3',
     &'   AMPAC (AMI) ............................ 4',
     &'   Z matrix (ZMT) ......................... 5',
     &' => '
          call upcase (10,an)
        endif
        if (an.eq.'2') then
          ext1 = '.CAR'
          icar = 1
        else if (an.eq.'3') then
          ext1 = '.IMC'
          icar = 2
        else if (an.eq.'4') then
          ext1 = '.AMI'
          icar = 3
        else if (an.eq.'5') then
          ext1 = '.ZMT'
          icar = 4
        else
          ext1 = '.IMM'
          icar = 0
        endif
      endif

c     Input file name

      if (qcl) then
        call rcomlin (2,suber)
        file5 = arg
      else
        write(10,'(/a,a,a,a,a\)')
     &  ' Enter the ',dd,' input file name <fn>[',ext1,'] => '
        read(10,'(a)') file5
      endif
      call filen1 (0,file5,ext1)
      if (kb1.ne.1) open(7,file=file5,status='old')
      if (irst.eq.0) call oufile ('.OMM')

      if (qcl) then
        if (numarg.gt.2) then
          call rcomlin (3,suber)
          if (arg.eq.'NOSCREENOUT') then
            an='Y'
          else
            an='N'
          endif
        else
          an='N'
        endif
      else
        write(10,'(/a\)')
     &' Do you want to supress the detailed run footprint on screen? (Y/
     &[N]) '
        call upcase (10,an)
      endif
      if (an.eq.'Y') iclos = 1

c     clearing the screen

      call clearscreen($gclearscreen)

  150 FORMAT (/' Date : ',I2,'/',I2.2,'/',I4,T60,'Time : ',I2,
     &2(':',I2.2),'.',I2.2/)
  170 FORMAT(/' Enter the main output file name. You can use the default
     . <fn> of input'/' file by striking <CR>. If you enter CON, the out
     .put will be written in the'/' screen, and if PRN, to the installed
     . printer.'/' [<fn>[',a,']] => '\)
      RETURN
      END

C     MASTER
C     CALLING ROUTINE FOR ALL SECTIONS OF THE PROGRAM

      SUBROUTINE MASTER (iclos)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 TYPE(15), TITLE,SCFTYP,RUNTYP,GUESS,CONF

      LOGICAL  LADAPT,LSKIP
      LOGICAL REVISE
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
      common/thermo/cpt,cpr,cpv,cpe,
     &              h25t,h25r,h25v,h25e,
     &              s25t,s25r,s25v,s25e,
     &              g25t,g25r,g25v,g25e,
     &              ilin,nax,
     &              temper
      parameter (ZERO=0.d0)
      character*1 an
      character*6 suber /'MASTER'/
      character*32 arg
      logical*2 qcl
      integer*2 marg,sta,numarg
      common /cl/ qcl,numarg,marg(5),sta,arg
      save ini, irot, TYPE
      DATA TYPE/
     1  'HFSCF'   ,'FINISH'  ,'INTEGRAL','SCF'     ,'GRADIENT',
     2  'GRADONE' ,'OPTIMIZE','FORCE'   ,'STOP'    ,'CONTINUE',
     3  'PROPERTY','ROTATION','NEW'     ,'CHECK'   ,'        '/

      ini = 0
      irot = 0
200   CALL INIT (ini)

      if (iclos.gt.0) then
        close (10)
        open(10,file='NUL')
      endif

20    CONTINUE
      write (*,'(/a)') ' >> MASTER is calling START'
      write (*,'(a)') ' (Lookig into main input file for data)'
      CALL START
      DO 30 I=1,15
        IF(RUNTYP.EQ.TYPE(I))GOTO 40
30    CONTINUE
      CALL HONDER(15,suber)

40    GOTO ( 50, 60, 80, 90,100,
     &      100,120,140,160,160,
     &      150,155,160, 60, 50) I

50    write (*,'(/a)') ' >> MASTER is calling HFSCF'
      write (*,'(a)') ' (Starting SCF procedure)'
      CALL HFSCF
      GOTO 20
80    write (*,'(/a)') ' >> MASTER is calling INTEG'
      write (*,'(a)') ' (Calculating electron integrals)'
      CALL INTEG
      GOTO 20
90    write (*,'(/a)') ' >> MASTER is calling SCFRUN'
      write (*,'(a)') ' (Calculating SCF energy)'
      CALL SCFRUN
      GOTO 20
100   IF(IREST.GT.3) GOTO 110
      write (*,'(/a)') ' >> MASTER is calling HFSCF'
      write (*,'(a)') ' (Starting SCF procedure)'
      CALL HFSCF
110   write (*,'(/a)') ' >> MASTER is calling HFGRAD'
      write (*,'(a)') ' (Calculating energy gradients)'
      CALL HFGRAD
      GO TO 20

C  GEOMETRY OPTIMIZATION

120   write (*,'(/a)') ' >> MASTER is calling OPTX'
      write (*,'(a)') ' (Calculating an optimized geometry)'
      CALL OPTX
      write (*,'(/a)') ' >> MASTER is calling ROTAT'
      write (*,'(a)') ' (Calculating rotational properties)'
      CALL ROTAT
      irot = 1
      GO TO 20

C  NUMERICAL FORCE CONSTANT EVALUATION

140   continue
      if (irot.eq.0) then
        write (*,'(/a)') ' >> MASTER is calling ROTAT'
        CALL ROTAT
        write (*,'(a)') ' (Calculating rotational properties)'
        irot = 1
      endif
      write (*,'(/a)') ' >> MASTER is calling FORCX'
      write (*,'(a)') ' (Calculating force constants and vibrational ene
     &rgies)'
      CALL FORCX
      GO TO 20

150   CONTINUE
      write (*,'(/a)') ' >> MASTER is calling PROP1E'
      write (*,'(a)') ' (Calculating electronic properties)'
      CALL PROP1E
155   continue
      if (irot.eq.0) then
        write (*,'(/a)') ' >> MASTER is calling ROTAT'
        write (*,'(a)') ' (Calculating rotational properties)'
        CALL ROTAT
        irot = 1
      endif
      GO TO 20

160   CALL CLENUP
      DO 6 ISEL=1,16
        GOTO (2,2,6,6,6,6,6,6,6,6,6,6,6,6,6,6), ISEL
* CLOSING DUMP FILES
2       CLOSE (ISEL)
6     CONTINUE

      if (runtyp.eq.type(10) .or. runtyp.eq.type(13)) then
        write (iw,71) char(12)
        write(10,71) char(12)
        ini = 1
        goto 200
      endif

        GO TO 151

60    CONTINUE
      CALL CLENUP
      DO 5 ISEL=1,16
        GOTO (1,1,5,5,5,5,5,5,5,5,5,5,5,5,5,5), ISEL
* CLOSING DUMP FILES
1       CLOSE (ISEL,STATUS='DELETE')
5       CONTINUE

151   WRITE(IW,70)
      write(10,70)
70    FORMAT(///1X,79('-')/34X,'ALL FINISHED'/1X,79('-'))
71    FORMAT(///1X,79('-')/36X,'FINISHED'/1X,79('-')/a)
      RETURN
      END

      subroutine start
      implicit real*8 (a-h,o-z)
      include 'dimmm'
      common
     ./optn2/ icar, kb1, iii
      character*8 words, opts
      parameter (NWORDS=80, NTYPES=20)
      DIMENSION WORDS(NWORDS)
      DIMENSION OPTS(NWORDS-NTYPES),TYPES(NTYPES)
      EQUIVALENCE (WORDS(1),OPTS(1)),(WORDS(NWORDS-NTYPES+1),TYPES(1))
      DATA OPTS
     1 / 'TITLE'   , 'SYMMETRY', 'PRINT'   , 'ACCURACY', 'NORM'    ,
     2   'LEVEL'   , 'RESTORE' , 'CHARGE'  , 'MULTIPLI', 'SCFTYP'  ,
     3   'OSCF'    , 'ANGSTROM', 'MAINFILE', 'DUMPFILE', 'SCRATCHF',
     4   'WEIGHTS' , 'GEOMETRY', 'ATOMS'   , 'MAXIT'   , 'CONVERGE',
     5   'GUESS'   , 'START'   , 'RUNTYP'  , 'SHIFT'   , 'DIIS'    ,
     6   'LOCK'    , 'TIME'    , 'SKIP'    , 'IMPTPRIN', 'END'     ,
     7   'FINISH'  , 'VARIABLE', 'CONSTANT', 'NOTE'    , 'NCOORB'  ,
     8   'SWAP'    , 'EXTRAPOL', 'CANONICA', 'WIDTH'   , 'GRADTOL' ,
     9   'DAMP12'  , 'DAMP13'  , 'DAMP23'  , 'STEPSIZE', 'NPOINTS' ,
     A   'STOP'    , 'CONTINUE', 'HUCKEL'  , 'DENSITY' , 'MOS'     ,
     B   'NONUCLEA', 'ORIGIN'  , 'GAUGE'   , 'NEW'     , 'AXIS'    ,
     C   'TEMPERAT',4*'        '/
      DATA TYPES
     $/  'HFSCF'   , 'CHECK'   , 'INTEGRAL', 'SCF'     , 'GRADIENT',
     $   'GRADONE' , 'OPTIMIZE', 'FORCE'   , 'PROPERTY', 'ROTATION',
     $10*'       '/

      goto (1,2,2,2), icar
      goto 10
1       call incar (words, NWORDS)
      goto 10
2       call intern (words, NWORDS)
10    call start1 (words, NWORDS, types, NTYPES)
      return
      end

      subroutine oufile (ext)
      character*4 ext
      character*32 file5, file6
      common
     ./files / file5
      common mes,idi,ian,iho,min,ise,ics
      integer*2 mes,idi,ian,iho,min,ise,ics

      file6 = file5
      call upcase(-1,file6)
      call filen1 (1,file6,ext)
      open (9,file=file6,status='unknown')
      write(9,150) mes,idi,ian,iho,min,ise,ics
      call machon (9)
  150 FORMAT (/' Date : ',I2,'/',I2.2,'/',I4/
     &         ' Time : ',I2,2(':',I2.2),'.',I2.2/)
      return
      end
