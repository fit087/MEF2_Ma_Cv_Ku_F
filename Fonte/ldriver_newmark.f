    subroutine ldriver_newmark
	implicit real*8 (a-h,o-z)

	include 'common.h'
	include 'cntl.h'
	include 'pointers.h'


	call bot ('ldrivr')

!   call pre_plas

!     Cria um vetorzão para almacenar todos os dados

	call locate ('stiff   ',nstiff,nwk,nc)
	call locate ('maxa    ',nmaxa ,nr ,nc)
	call locate ('f0      ',nf0   ,nr ,nc)


	call locate ('id      ',nid,numnp,ngl)

	call define ('disp    ',ndisp ,neq+1,1)
	call define ('fi      ',nfi   ,neq+1,1)
	call define ('f       ',nf    ,neq+1,1)
    
    call define ('upred     ',nupred,neq+1,1)
    call define ('vpred      ',nvpred,neq+1,1)
    
    call define ('ucorr      ',nucorr,neq+1,1)
    call define ('vcorr      ',nvcorr,neq+1,1)


	call define ('du      ',ndu   ,neq+1,1)

	call dclear (ia(ndisp),neq+1)
    
    
    nsteps  = int(timef/dt)
    nflag   = nsteps/nfile
    
    
    
	call solver_newmark (ia(nstiff),ia(nmassa),ia(nmaxa),ia(ndisp),ia(nf),ia(nf0),
     &            ia(nfi),ia(ndu),nwk,ia(nid),ia(nupred),ia(nvpred),ia(nucorr),
     &           ia(nvpred),nsteps,nflag,timef,dt,betan,gamman)
     
     !deve adicionar ia(ndamp) se houvese
    
 !
	!call solver (ia(nstiff),ia(nmaxa),ia(ndisp),ia(nf),ia(nf0),
 !    &            ia(nfi),ia(ndu),nwk,ia(nid))

	call eot ('ldrivr')

	return
	end
