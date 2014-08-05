	subroutine joints (id,x,y,z)
	implicit real*8 (a-h,o-z)

      include 'cntl.h'
	include	'tapes.h'
	dimension         id      (ngl,numnp)     ,
     .                  x       (numnp)         ,
     .                  y       (numnp)         ,
     &                  z       (numnp)         ,
     .                  inc     (6)             

	call bot ('joints')

	write (iout,101)

	do i = 1, numnp

		read  (iin ,100) n,inc,x(n),y(n),z(n)
		write (iout,200) n,inc,x(n),y(n),z(n)

		if ( n.gt.numnp.or.n.lt.1.or.n.ne.i) then
			write(iout,300) n,i
			stop
		endif

		do igl = 1, 6
			jgl = idof(igl)
			if (jgl.ne.0) id(jgl,n) = inc(igl)
		enddo

	enddo

	call eot ('joints')
	return

  100   format (7i5,3f10.0)
  101	format (/,' Node Data ',/
     .		,' node #',t26,'id',t66,'coordinates',/)	
  200   format (' ',i5,5x,6i5,5x,1p,3e15.7)
  300   format (' *** (JOINTS) Data Error : node .gt. numnp (',
     .            i5,i5')')

	end




