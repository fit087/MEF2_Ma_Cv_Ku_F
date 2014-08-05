	subroutine ldriver
	implicit real*8 (a-h,o-z)

	include 'common.h'
	include 'cntl.h'
	include 'pointers.h'


	call bot ('ldrivr')

c	call pre_plas

	call locate ('stiff   ',nstiff,nwk,nc)
	call locate ('maxa    ',nmaxa ,nr ,nc)
	call locate ('f0      ',nf0   ,nr ,nc)


	call locate ('id      ',nid,numnp,ngl)

	call define ('disp    ',ndisp ,neq+1,1)
	call define ('fi      ',nfi   ,neq+1,1)
	call define ('f       ',nf    ,neq+1,1)

	call define ('du      ',ndu   ,neq+1,1)

	call dclear (ia(ndisp),neq+1)

	call solver (ia(nstiff),ia(nmaxa),ia(ndisp),ia(nf),ia(nf0),
     &            ia(nfi),ia(ndu),nwk,ia(nid))

	call eot ('ldrivr')

	return
	end
