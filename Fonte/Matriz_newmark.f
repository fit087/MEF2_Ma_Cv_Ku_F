    subroutine matriz
	implicit real*8 (a-h,o-z)

	include         'common.h'
	include         'cntl.h'
	include         'tapes.h'
	include         'pointers.h'

	call bot ('matriz')

	call defini ('maxa     ',nmaxa,neq+1,1)
	call defini ('mht      ',nmht,neq,1)
	call locate ('lm       ',nlm,nume,ndt)

	call profil2 (ia(nmaxa),ia(nmht),ia(nlm),neq,nume,ndt,nwk)

	call define ('stiff    ',nstiff,nwk,1)
    
    call define ('massa    ',nmassa,nwk,1)    
!	call dclear (ia(nstiff),nwk)

      call assemb
      call assemb_massa

	call eot ('matriz')

	return
	end




