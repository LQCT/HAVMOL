      SUBROUTINE TIMIT(INDEX)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 INDEX
      INTEGER*2 IH,IM,IS,IHS,IA,ME,ID,N1,ID1,IAP1,MEP1,IDP1
      COMMON/TIMEX/TSTART,TI,TX,TIM,TTGO,TIMLIM,TTOT
      PARAMETER (N1=1, SES=60.D0, CENT=0.01D0, CERO=0.D0, OM=86400.D0)
      LOGICAL SETDAT, QD
      EXTERNAL SETDAT
      SAVE

* IF INDEX
*         =0 TIME IS CHECKED AND A BALANCE AND DATE ARE PRINTED
*         =1 TIME IS CHECKED AND A BALANCE IS ONLY PRINTED IN OUTPUT
*         =2 ALL TIME INDEX ARE SET TO 0 AND TSTART TAKES THE AMOUNT OF
*            SECONDS FROM 0:00 h TODAY
*         =3 TIME IS CHECKED
*         =4 TIME IS CHECKED AND A BALANCE IS ONLY PRINTED IN THE SCREEN
* TSTART is the amount of seconds elapsed from 00:00:00.00 today
* TIM is the amount of seconds from start, today

      CALL GETTIM(IH,IM,IS,IHS)
      IF (INDEX.EQ.2) THEN
        DIH=DBLE(IH)
        DIM=DBLE(IM)
        DIS=DBLE(IS)
        DIHS=DBLE(IHS)
        TSTART=(DIH*SES+DIM)*SES+DIS+CENT*DIHS
        TIM=CERO
        TX=CERO
        TI=CERO
      ELSE
        DIH=DBLE(IH)
        DIM=DBLE(IM)
        DIS=DBLE(IS)
        DIHS=DBLE(IHS)
        TIM1=(DIH*SES+DIM)*SES+DIS+CENT*DIHS-TSTART
        IF(TIM1.LT.TIM) THEN
          TSTART=TSTART-OM
          TIM1=TIM1+OM
          CALL GETDAT (IA,ME,ID)
          IDP1 = ID + N1
          MEP1 = ME + N1
          IAP1 = IA + N1
          IF (ID1.EQ.ID) THEN
            IF (SETDAT(IA,ME,IDP1)) THEN
              ID1 = IDP1
              GO TO 1
            ELSEIF (SETDAT(IA,MEP1,N1)) THEN
              ID1 = N1
              GO TO 1
            ELSE
              ID1 = N1
              QD = SETDAT(IAP1,N1,N1)
            ENDIF
          ENDIF
        ENDIF
1       TIM=TIM1
        TX=TIM-TI
        if (index.ne.4) TI=TIM
        TI=TIM
        TTGO=TIMLIM-TIM
        TTOT = TTOT + TX
      ENDIF
      IF(INDEX.LE.1) THEN
        write(9,1001)
        write(9,1002) tx
        write(9,1003) tim,ttgo
        write(*,1001)
        write(*,1003) tim,ttgo
        if (index.eq.0) then
          call getdat (ia,me,id)
          write (9,150) me,id,ia,ih,im,is,ihs
          write(10,150) me,id,ia,ih,im,is,ihs
          id1 = id
        endif
      ENDIF
150   FORMAT (//' Date : ',I2,'/',I2.2,'/',I4,T60,'Time : ',I2,2(':',I2.
     .2),'.',I2.2/)
1001  FORMAT(//1X,79('-'))
1002  format(24X,'Last step used',T41,F12.2,' s')
1003  format(24X,'Total used',T41,F12.2,' s'
     &/24X,'Remaining ',T41,F12.2,' s'
     &/1X,79('-'))
      RETURN
      END
