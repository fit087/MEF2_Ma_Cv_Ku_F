C    ------------------
      PROGRAM TRIEL_A
C    ------------------

C
C **********************************************************************
C *               Triangulo Elasticidade 2D                            *
C *               elemento tria3                                       *
C *                                                                    *
C **********************************************************************
C COPPE-PEC/UFRJ/1995 JLDA
C COPPE/LAMCE 2005 CES + JRPG + NOGJ

      CALL TIME_TABLE        !controla o tempo de execução
      CALL IOMNGR (0)       !preparas os arquivos para a leitura, prepara a memoria e todos os arquivos
      CALL BOT ('triel ')

	call setdb
	call contrl
	
	call inmesh
	call loads       !le os carregamentos
	
	call matriz
!	call ldriver        !resolve.deve mexerse para adicionar o transiente
      !call ldriver_newmark
	
      CALL EOT ('triel ')
      CALL IOMNGR (1)

      STOP
      END
