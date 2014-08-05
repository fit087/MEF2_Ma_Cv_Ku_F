	subroutine loads
	implicit real*8 (a-h,o-z)

	include         'common.h'
      	include         'cntl.h'
	include         'tapes.h' 	
	
	call bot ('loads ')

            
	call locate ('id      ',nid,nr,nc)
	call define ('f0      ',nf0,neq+1,1)
	
c     cargas nodais	

	call ptload (ia(nf0),ia(nid),idof,ncarg,numnp,neq,ngl)
      
      call eot ('loads ')
      
	end




