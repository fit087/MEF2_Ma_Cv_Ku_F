C    -------------------
      SUBROUTINE SYSTIM
     .           ( VECTOR , TOTCPU , ELAPSD , CTIME , CDATE , IRC)
C    -------------------
C
C
C     RETORNA TEMPOS E DATA DO SISTEMA
C
C
      CHARACTER*10 CTIME , CDATE
c      CHARACTER*24 buffer1, buffer2

      REAL*8      VECTOR , TOTCPU , ELAPSD
c      REAL        tarray(2)
c	integer*4 iticks 
	integer*2   iyr,imon,iday,ihr,imin,isec,i100th

c ---------------------------------------------------------------------
c                1234567890

        CTIME = '          '
        CDATE = '          '

        
c     Lahey 5.00

c        call time (buffer1)
c        do i = 1, 8
c                CTIME(i:i) = buffer1(i:i)
c        enddo
c
c        call date (buffer2)
c        do i = 1, 9
c                CDATE(i:i) = buffer2(i:i)
c        enddo
        
c     Microsoft Fortran (PowerStation, etc.)        
        
        call getdat (iyr,imon,iday)
        write(CDATE(1:10),'(i2.2,1h/,i2.2,3h/19,i2.2)') iday,imon,iyr
        do i = 1,10
          if (CDATE(i:i) .eq. ' ') CDATE(i:i) = '0'
        enddo
        
        call gettim (ihr,imin,isec,i100th)
        write(CTIME(1:8),'(i2.2,1h:,i2.2,1h:,i2.2)') ihr,imin,isec
        do i = 1,8
          if (CTIME(i:i) .eq. ' ') CTIME(i:i) = '0'
        enddo
                
        elapsd = (ihr*3600 + imin*60 + isec) + dfloat(i100th)*1.d-2
        totcpu = elapsd
        vector = elapsd 
        
        irc = 0
	
c	Sun - Vax library

c        ELAPSD  = etime (tarray)
c        TOTCPU  = tarray(1)
c        VECTOR  = tarray(2)

c 	Lahey 5.00

c	call timer (iticks)
c
c        ELAPSD  = iticks * 1.e-2
c        TOTCPU  = elapsd
c        VECTOR  = elapsd


      RETURN
C
      END
C
C ----------------------------------------------------------------------
C
C     CONTROLE DOS TEMPOS DE EXECUCAO E DE CPU
C
C ----------------------------------------------------------------------
C
C     J.ALVES - 900401 - COPPE/UFRJ
C
C      UTLILIZACAO :
C      -------------
C
C        1) AO INICIAR O PROGRAMA CHAME A ROTINA "TIME-TABLE" PARA
C           A INICIALIZACAO DAS TABELAS DE CONTROLE DE TEMPOS.
C
C        2) A PARTIR DAI, PARA CADA ROTINA (OU TRECHO) DO PROGRAMA
C           CHAME O ENTRY-POINT "BOT" COM O PARAMETRO "TSKID" (CHAR*6)
C           IDENTIFICANDO O NOME DA ROTINA (SIMPLES IDENTIFICACAO).
C
C        3) AO FINAL DO TRECHO CHAME O ENTRY-POINT "EOT" COM O PARAMETRO
C           "TSKID".
C
C      OPERACAO :
C      ----------
C
C           O PRESENTE MODULO CONSTROI UMA TABELA DE TEMPOS LIQUIDOS DE
C         EXECUCAO (ELAPSED) E DE CPU PARA CADA ROTINA (OU TRECHO) DO
C         PROGRAMA, MANTENDO UM REGISTRO (LOG) DE SUA ATIVIDADE.
C
C           O MODULO REQUER A UTILIZACAO DE ROTINAS DE TEMPO ESPECIFICAS
C         DO AMBIENTE COMPUTACIONAL. NO CASO, A CHAMADA PARA TAIS ROTI-
C         NAS FOI ESTABELECIDA NA ROTIMA "SYSTIM"
C
C           O LOG E GERADO NA UNIDADE LOGICA "ILOG" QUE DEVERA ESTAR
C         PREVIAMENTE DEFINIDA E CONECTADA.
C
C           AS "TASKS" (ROTINAS OU TRECHOS DO PROGARMA) SAO ORGANIZADAS
C         SEGUNDO OS SEUS NIVEIS HIERARQUICOS, SENDO "0" A DE NIVEL MAIS
C         ALTO, "1" A DE NIVEL IMEDIATAMENTE INFERIOR E ASSIM POR DIANTE
C         TAL ESTRUTURA CONSTITUI UMA "STACK LAST-IN-FIRST-OUT".
C
C           O NUMERO DE TASKS E LIMITADO PELO PARAMTERO "MTASK" E A PRO-
C         FUNDIDADE MAXIMA DA STACK PELO PARAMETRO "MLEVEL".
C
C           AO ENCERRAR-SE A TASK "0" A TABELA GERAL DE ATIVIDADE SERA
C         IMPRESSA NA UNIDADE "ILOG".
C
C ---------------------------------------------------- (JLDA 900401) ---
C
C    -----------------------
      SUBROUTINE TIME_TABLE
C    -----------------------

      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER ( MTASK  = 30   ,     ! MAXIMUM NO. OF TASKS
     .            MLEVEL =  5 )       ! MAXIMUM NO. OF NESTED TASKS
                                      !                (STACK DEPTH)

      CHARACTER*6 TASKID , TSKID
      CHARACTER*9 STATUS
      CHARACTER*3 FLAG
      CHARACTER*1 CR
      CHARACTER*10 CDATE
      CHARACTER*10 CTIME
c      REAL*8      TOTAL,SCALAR,VECTOR

      COMMON / TASKS / LEVEL               , ! CURRENT LEVEL
     .                 NTASK               , ! CURRENT TASK
     .                 IRCACC              , ! RETURN CODE FROM CLOCKX
     .                 IALIGN              , ! DOUBLE WORD ALIGNMENT
     .                 DVCPU   ( 0:MTASK ) , ! VIRTUAL CPU TIME
     .                 DTCPU   ( 0:MTASK ) , ! TOTAL   CPU TIME
     .                 DTELP   ( 0:MTASK ) , ! ELAPSED     TIME
     .                 DVCP0               ,
     .                 DTCP0               ,
     .                 DELP0               ,
     .                 IACTIVE ( 0:MLEVEL) , ! STACK OF ACTIVE TASKS
     .                 TASKID  ( 0:MTASK )   ! TASK IDENTIFICATOR

      include 'tapes.h'

      LEVEL  = 0
      NTASK  = 0
      IRCACC = 0

      RETURN

C    -----------
      ENTRY BOT   ( TSKID )
C    -----------

      TASKID  (NTASK) = TSKID                ! STORE CALLED TASK ID
      IACTIVE (LEVEL) = NTASK                ! PUT CURRENT TASK ON STACK

      CALL SYSTIM ( DVCPU (NTASK) ,          ! GET TIMING INFO FOR
     .              DTCPU (NTASK) ,          !     CURRENT TASK AT
     .              DTELP (NTASK) ,          !     ITS BEGINNING
     .              CTIME          ,          !
     .              CDATE     , IRC     )     !

      IRCACC = IRCACC + IRC

      STATUS = ' STARTED '                   ! SET STATUS

      IF ( LEVEL .EQ. 0 ) THEN               ! CHECK CONDITION FOR
           CR   = '1'                        !
           FLAG = 'BOJ'                      !       "BEGIN OF JOB"
           DVCP0 = DVCPU (0)
           DTCP0 = DTCPU (0)
           DELP0 = DTELP (0)
      ELSE                                   ! OR
           CR   = ' '                        !
           FLAG = 'BOT'                      !       "BEGIN OF TASK"
      ENDIF

      WRITE (ILOG,6910) CR,CTIME,FLAG,LEVEL,TSKID,STATUS

      NTASK   = NTASK + 1                    ! INCREMENT TASK COUNTER
      LEVEL   = LEVEL + 1                    ! SET NEXT STACK POSITION

      RETURN

C    -----------
      ENTRY EOT  ( TSKID )
C    -----------

      LEVEL   = LEVEL - 1                    ! DECREASE STACK DEPTH
      ITASK   = IACTIVE ( LEVEL )            ! RETRIEVE CURRENT TASK PTR

      CALL SYSTIM ( RVCPU ,                  ! GET TIMING INFO FOR
     .              RTCPU ,                  !     CURRENT TASK AT
     .              RTELP ,                  !     ITS END
     .              CTIME  ,                  !
     .              CDATE  , IRC )            !

      IRCACC = IRCACC + IRC

      STATUS = ' ENDED   '                   ! SET STATUS TO 'ENDED'

      DVCPU (ITASK) = RVCPU - DVCPU ( ITASK )   ! COMPUTE TIMING VALUES
      DTCPU (ITASK) = RTCPU - DTCPU ( ITASK )   !         FOR CURRENT
      DTELP (ITASK) = RTELP - DTELP ( ITASK )   !         TASK

      VCPU = DVCPU ( ITASK )
      TCPU = DTCPU ( ITASK )
      TELP = DTELP ( ITASK )

      DO ILEVEL = 0 , LEVEL-1
         IPTASK = IACTIVE (ILEVEL)
         DVCPU ( IPTASK ) = DVCPU ( IPTASK ) + DVCPU ( ITASK )
         DTCPU ( IPTASK ) = DTCPU ( IPTASK ) + DTCPU ( ITASK )
         DTELP ( IPTASK ) = DTELP ( IPTASK ) + DTELP ( ITASK )
      ENDDO

C     LEVEL = LEVEL - 1                      ! DECREASE QUEUE DEPTH

      IF ( LEVEL .EQ. 0 ) THEN               ! CHECK CONDITION FOR
           CR   = ' '                        !
           FLAG = 'EOJ'                      !      "END OF JOB"
           VCP0 = ( RVCPU - DVCP0 )
           TCP0 = ( RTCPU - DTCP0 )
           ELP0 = ( RTELP - DELP0 )
      ELSE                                   ! OR
           CR   = ' '                        !
           FLAG = 'EOT'                      !      "END OF TASK"
      ENDIF

      WRITE (ILOG,6920) CR,CTIME,FLAG,LEVEL,TSKID,STATUS
     .                 ,TCPU,VCPU,TELP

      IF ( LEVEL .EQ. 0 ) THEN                 ! CHECK IF CURRENT
           WRITE (ILOG,6940)                   ! TASK IS THE MAIN TASK
           DO ITASK = 0 , NTASK-1              ! IF SO
                VCPU = DVCPU(ITASK)
                TCPU = DTCPU(ITASK)            !    OUTPUT TIMING TABLE
                TELP = DTELP(ITASK)            !
                PT   = 1.D0
                PJ   = 1.D0
                PV   = 1.D0
                IF ( TELP.NE.0.) THEN
                   PT   =  TCPU/TELP * 100.D0
                ENDIF
                IF (TCP0 .NE. 0 ) THEN
                   PJ   =  TCPU/TCP0 * 100.D0
                ENDIF
                IF (TCPU .NE. 0 ) THEN
                  PV   =  VCPU/TCPU * 100.D0
                ENDIF
                                               !            FOR
                WRITE (ILOG,6950)              !
     .                    TASKID(ITASK)    ,   !         ALL TASKS     .                    TCPU             ,   !
     .                    VCPU             ,   !
     .                    TELP             ,   !
     .                    PV               ,
     .                    PT               ,
     .                    PJ
           ENDDO                               !
           VCEP = 0.0
           TCEP = 0.0
           IF (ELP0 .NE. 0 ) TCEP = TCP0/ELP0 * 100.
           IF (TCP0 .NE. 0 ) VCEP = VCP0/TCP0 * 100.
           WRITE (ILOG,6960) TCP0,VCP0,ELP0,VCEP,TCEP , IRCACC

      ENDIF                                    !


      RETURN

 6910 FORMAT (A,A,' (',A,') - ',I2,' ',A,A)
 6920 FORMAT (A,A,' (',A,') - ',I2,' ',A,A,
     .                         ' T',F10.3,
     .                         ' V',F10.3,
     .                         ' E',F10.3   )
 6940 FORMAT (
     .'1-------------------------------------------------------------',
     .'---------',/,
     .'  TASK-ID      TOTAL      SYSTEM      ELAPSED    SYS     PCPU',
     .'    PJOB',/,
     .' -------------------------------------------------------------',
     .'---------',/)
 6950 FORMAT (' ',2X,A6,3(2X,F10.3),3(2X,F6.2) )
 6960 FORMAT (/,
     .' -------------------------------------------------------------',
     .'---------',///,
     .' -------------------------------------------------------------',
     .'---------',//,
     .    ' ',2X,'TOTAIS',//,
     .    ' ',2X,'CPU TOTAL            ',F15.3,' SEG',/,
     .    ' ',2X,'CPU VETORIAL         ',F15.3,' SEG',/,
     .    ' ',2X,'ELAPSED              ',F15.3,' SEG',//,
     .    ' ',2X,'PORCENTAGEM VETORIAL ',F15.3,' %  ',/,
     .    ' ',2X,'PORCENTAGEM CPU      ',F15.3,' %  ',//,
     .    ' ',2X,'RETURN CODE CPUTIM   ',I15  ,' (ACUMULADO)',//,
     .' -------------------------------------------------------------',
     .'---------' )

      END


