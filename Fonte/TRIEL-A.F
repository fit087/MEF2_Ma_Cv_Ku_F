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

      CALL TIME_TABLE
      CALL IOMNGR (0)
      CALL BOT ('triel ')

	call setdb
	call contrl
	
	call inmesh
	call loads
	
	call matriz
	call ldriver
	
      CALL EOT ('triel ')
      CALL IOMNGR (1)

      STOP
      END
