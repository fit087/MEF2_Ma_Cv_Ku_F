	subroutine ldriver_newmark
	implicit real*8 (a-h,o-z)

	include 'common.h'
	include 'cntl.h'
	include 'pointers.h'
      include 'time.h'


	call bot ('ldrivr')

!	call pre_plas
      
      call locate ('x       ',nx ,numnp,nc )
	call locate ('y       ',ny ,numnp,nc )
	call locate ('z       ',nz ,numnp,nc )
      call locate ('incid   ', ninc  , nume   ,nnoel)
      call locate ('prop    ', nprop , nummat ,ncprop)
      call locate ('mtype   ', nmtype, nume   ,nc)
      call locate ('lm      ', nlm   , nume   ,ndt)
      
	call locate ('stiff   ',nstiff,nwk,nc)
	call locate ('maxa    ',nmaxa ,nr ,nc)
	call locate ('f0      ',nf0   ,nr ,nc)


	call locate ('id      ',nid,numnp,ngl)
      
	call define ('fi      ',nfi   ,neq+1,1)
	call define ('f       ',nf    ,neq+1,1)
    
      call define ('upred   ',nupred,neq+1,1)
      call define ('vpred   ',nvpred,neq+1,1)
    
      call define ('ucorr   ',nucorr,neq+1,1)
      call define ('vcorr   ',nvcorr,neq+1,1)

	call define ('du      ',ndu   ,neq+1,1)
    
         nsteps = int(timef/dt)
         nflag  = int(nsteps/nfile)

	call solver_newmark (ia(nx),ia(ny),ia(nz),
     &            ia(ninc),ia(nprop),ia(nmtype),ia(nlm),ia(nstiff),
     &            ia(nmaxa),ia(nf),ia(nf0),ia(nfi),ia(ndu),
     &            nwk,ia(nid),ia(nupred),ia(nvpred),ia(nucorr),
     &            ia(nvcorr),nsteps,nflag,timef,dt,betan,gamman,
     &            damp1,damp2)

	call eot ('ldrivr')

	return
	end
