      subroutine newmark(stiff,massa,nwk,betan,gamman,dt,xke)
      
      implicit real*8 (a-h,o-z)
      
      dimension stiff(nwk),massa(nwk),xke(nwk),f0(0:neq),
      & upred(0:neq),an(0:neq)
      
      do i=1,nwk
          xke(i)=massa(i)+betan*dt**2*stiff(nwk)
          an(i)=f0(i)-stiff(i)*upred              !vetor for�a de newmark
          
      enddo
      
      
            do i=1,nwk
                        do n=
                
                          an(i)=f0(i)-stiff(i)*upred              !vetor for�a de newmark

                
                      enddo

                                    !stiff n�o � da mesma dimens�o de f0?
      
      
      
      call COLSOL(xke,an,maxa,neq,nwk,neq1,1,iplt)!escalonando     an � for�a
      call COLSOL(xke,an,maxa,neq,nwk,neq1,1,iplt)!resolvendo      an cuspe an acelera��o
      
      
      
      end subroutine