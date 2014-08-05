      subroutine preditor(neq,dt,betan,gamman,ucorr,vcorr,
     &                    an,upred,vpred)
      
      implicit real*8 (a-h,o-z)
      
      dimension  ucorr (0:neq)    ,
     &           vcorr (0:neq)    ,
     &           upred (0:neq)    ,
     &           vpred (0:neq)    ,
     &           an    (0:neq)
      
      DO ieq=1,neq
          upred(ieq)=ucorr(ieq)+vcorr(ieq)*dt+(0.5d0*dt**2)*(1.d0-2.d0*
     &               betan)*an(ieq)
          vpred(ieq)=vcorr(ieq)+dt*(1.d0-gamman)*an(ieq)
      ENDDO
      
      upred(0)=0.d0   !notaçao cientifica .d0 = E+00
      vpred(0)=0.d0
      
      end