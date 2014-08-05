      subroutine newmark (stiff,massa,f0,upred,vpred,an,nwk,betan,gamman,dt,xke)
      
      implicit real*8 (a-h,o-z)
      
      dimension stiff (nwk), massa (nwk), xke (nwk), f0(0:neq), 
     &                 upred (0:neq), vpred(0:neq), an (0:neq)
      
      DO i=1,nwk
          xke(i) = massa(i) + betan*dt**2*stiff(nwk)
      ENDDO
      
      call COLSOL (stiff,fk
      
      call COLSOL (xke,an,maxa,neq,nwk,neq1,1,iplt)
      call COLSOL (xke,an,maxa,neq,nwk,neq1,2,iplt)
      
      
      endsubroutine
      