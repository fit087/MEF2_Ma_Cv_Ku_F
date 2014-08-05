      SUBROUTINE PROFIL2 ( MAXA,MHT,LM,NEQ,NEL,NDT,NWK )
C.......................................................................
C.                                                                     .
C.    ROTINA PARA COMPUTAR OS PARAMETROS DO SKYLINE                    .
C.                                                                     .
C.    output maxa,mht,nwk    
c 
c
c     input lm,neq,nel,nd
c
c PARAMETROS :                                                     .
C.                                                                     .
C.    IEN    --> TABELA DAS INCIDENCIAS                                .
C.    ID     --> MAPA DAS EQUACOES CORRESPONDENTES AOS NOS             .
C.    MAXA   --> PONTEIRO DA DIAGONAL DA MATRIZ SKYLINE                .
C.    MHT    --> ALTURAS EFETIVAS DE COLUNA                            .
C.    LM     --> MAPA DOS DOF DOS ELEMENTOS                            .
C.                                                                     .
C.    NEQ    --> NUMERO DE EQUACOES                                    .
C.    NEL    --> NUMERO DE ELEMENTOS                                   .
C.    NNOEL  --> NUMERO DE NOS POR ELEMENTO                            .
C.    NGL    --> NUMERO DE DOF POR NO                                  .
C.    ND     --> NUMERO DE DOF POR ELEMENTO ( NNOEL*NGL )              .
C.    NWK    --> NUMERO DE TERMOS DA MATRIZ SKYLINE                    .
C.                                                                     .
C.                                                                     .
C.......................................................................
C
	include 'tapes.h'

      INTEGER          MAXA ( 1 )                  ,
     .                 MHT  ( NEQ )                ,
     .                 LM   ( nel, NDT )

C --- CALCULA AS ALTURAS MAXIMAS DE COLUNA

	call bot ('PROFIL')
      CALL IZERO ( MHT , NEQ )

      DO 50 KEL = 1 , NEL

	LS = 1000000

	DO 30 I = 1 , NDT
	  IF ( LM(kel,I) .gt. 0 )   THEN
	    IF ( LM(kel,I) .LT. LS )    LS = LM(kel,I)
	  ENDIF
   30   CONTINUE

	DO 40 I = 1 , NDT
	  IF ( LM(kel,I) .gt. 0 )   THEN
	    ME = LM(kel,I) - LS
	    IF ( ME .GT. MHT(LM(kel,I)) )   MHT(LM(kel,I)) = ME
	  ENDIF
   40   CONTINUE

   50 CONTINUE

C --- CALCULA OS ENDERECOS DA DIAGONAL PRINCIPAL

      NEQ1 = NEQ + 1

      CALL IZERO ( MAXA , NEQ1 )

      MAXA (1) = 1
      MAXA (2) = 2
      MK       = 0

      IF ( NEQ .GT. 1 )   THEN

	  DO 60 KEQ = 2 , NEQ
	    IF ( MHT ( KEQ ) .GT. MK )  MK = MHT ( KEQ )
	    MAXA ( KEQ+1 ) = MAXA ( KEQ ) + MHT ( KEQ ) + 1
   60     CONTINUE

      ENDIF

      MK  = MK + 1
      NWK = MAXA ( NEQ1 ) - MAXA ( 1 )
      MM  = NWK / NEQ

      WRITE ( IOUT,1000 )   NEQ , NWK , MK , MM

	call eot ('PROFIL')

      RETURN

C --- FORMATOS

 1000 FORMAT(///,
     .' Informacoes do Sistema de Equacoes:',/,
     .t10,'Numero de Equacoes        (neq) = ',I10,/,
     .t10,'Coeficientes Armazenados  (nwk) = ',I10,/,
     .t10,'Meia-Banda Maxima         (mk ) = ',I10,/,
     .t10,'Meia Banda Media          (mm ) = ',I10 )

      END



	SUBROUTINE DCLEAR ( A , N )
      DOUBLE PRECISION A(*)
      DO 001 I = 1 , N
  001    A(I) = 0.D00
      RETURN
      END


	SUBROUTINE IZERO ( IA , N )
      INTEGER IA(*)
      DO 001 I = 1 , N
  001   IA(I) = 0
      RETURN
      END


	SUBROUTINE DCLEAR2 ( A , N , M)
      DOUBLE PRECISION A
	DIMENSION        A(N,M)
      
	DO I = 1 , N
	   DO J = 1 , M
	      
		  A(I,J) = 0.D00
	   
	   ENDDO
	ENDDO

      RETURN
      END