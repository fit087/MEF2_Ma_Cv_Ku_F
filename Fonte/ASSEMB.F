      subroutine assemb
      
c     chama rotina de elemento para montagem da matriz global linear      
      
      implicit real*8 (a-h,o-z)
      include 'common.h'
      common  /mesh/  nnoel,ndt
      call bot ('assemb')

      call locate ('x       ', nx    , numnp  ,nc)
      call locate ('y       ', ny    , numnp  ,nc)
      call locate ('incid   ', ninc  , nume   ,nnoel)
      call locate ('prop    ', nprop , nummat ,ncprop)
      call locate ('mtype   ', nmtype, nume   ,nc)
      call locate ('lm      ', nlm   , nume   ,ndt)
      call locate ('stiff   ', nstiff, nwk    ,nc)
      call locate ('maxa    ', nmaxa , neq1   ,nc) 
      

      call trie2d  (ia(nstiff),ia(nmaxa),ia(nlm),
     .              ia(nx),ia(ny),ia(ninc),
     .              ia(nmtype),ia(nprop),
     .              numnp, nume, nummat, nwk, ndt, neq1)

  
	call eot ('assemb')
      return
      end



