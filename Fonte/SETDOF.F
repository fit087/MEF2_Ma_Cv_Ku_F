	subroutine setdof (idof,ngl)
	dimension idof(6)

	NGL = 0
	DO KDOF = 1,6
		IF ( IDOF(KDOF) .EQ. 0 ) THEN
			NGL        = NGL + 1
			IDOF(KDOF) = NGL
		ELSE
			IDOF(KDOF) = 0
		ENDIF
	enddo

	return
	end




