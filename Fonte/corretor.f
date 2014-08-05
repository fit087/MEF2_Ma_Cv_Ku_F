      subroutine corretor(neq,dt,betan,gamman,upred,vpred,
     &                    an,ucorr,vcorr)
      
      implicit real*8 (a-h,o-z)
      
      dimension  ucorr (0:neq)    ,
     &           vcorr (0:neq)    ,
     &           upred (0:neq)    ,
     &           vpred (0:neq)    ,
     &           an    (0:neq)
      
      DO ieq=1,neq
          ucorr(ieq)=upred(ieq)+betan*dt**2*an(ieq)
          vcorr(ieq)=vpred(ieq)+gamman*dt*an(ieq)
      ENDDO
      
      end