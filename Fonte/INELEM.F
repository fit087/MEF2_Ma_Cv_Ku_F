	subroutine inelem (id,incid,mtype,lm,nume,nnoel,ndt,ngl)
	include	'tapes.h'
	dimension	incid	(nume,1)	,
     .			lm	    (nume,1)	,
     .			mtype	(nume)		,
     .			id	    (ngl,1)
 
	call bot ('inelem')

	write (iout,101)

	do i = 1, nume

		read  (iin,100)  n,(incid(n,ii),ii=1,nnoel),mtype(n)
		write (iout,200) n,(incid(n,ii),ii=1,nnoel),mtype(n)

		if (n.gt.nume.or.n.lt.1) then
			write (iout,300) n
			stop
		endif

	enddo

	kk = 0

	do j = 1, nnoel
		do igl = 1, ngl
			kk = kk + 1
			do i = 1, nume
				lm(i,kk) = id(igl,incid(i,j))
			enddo
		enddo
	enddo	

	call eot ('inelem')
	return

  100	format (6i5)

  101	format (/,' Dados dos Elementos',//,
     .		' elem #',t16,'incid',t51,'mat #')
  200	format (' ',i5,5x,4i5,5x,i5)
  300	format (' *** (INELEM) Erro de Dados: número de elemento',
     .		' inválido (',i5,')')

	end




