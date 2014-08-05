	subroutine ptload (f0,id,idof,nloadc,numnp,neq,ngl)
	implicit real*8 (a-h,o-z)
	include 'tapes.h'

	dimension       f0      (0:neq)         ,
     .                  id      (ngl,numnp)     ,
     .                  idof    (6)             ,
     .                  force   (6)

	call bot ('ptload')

	write (iout,101) 

	do i = 1, nloadc

		read  (iin ,100)  n,force
		write (iout,200)  n,force

		if (n.gt.numnp.or.n.le.0) then
			write (iout,300) n
			stop
		endif

		do ii = 1, 6
			idf = idof(ii)
			if (idf.ne.0) then
			    ieq = id(idf,n)
			    f0(ieq) = f0(ieq) + force(ii)
			endif
		enddo
		f0(0) = 0.d0

	enddo
	
	call eot ('ptload')
	return

  100   format (i5,6f10.0)
  101	format (/,' Cargas Nodais:',/,
     .		  ' no #',t40,'vetor de carga',/)
  200   format (' ',i5,1p,6e15.7)
  300   format (' *** (PTLOAD) Erro de Dados: no invalido ',
     .		'(',i5,')')

	end
	