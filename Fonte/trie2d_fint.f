      subroutine trie2d_fint (stiff,fp,up,maxa,lm,x,y,z,incid,
     &                          mtype,prop,numnp,nume,nummat,nwk,nd,
     &                          neq1,neqp,nnoel,ncprop,
     &                          damp1,damp2,betan,gamman,dt,neq,
     &                          fi,upred,vpred)
      
      implicit real*8 (a-h,o-z)
      
      include 'tapes.h'
      
      dimension       stiff  (nwk)           ,
     &                fp     (0:neq)         ,
     &                up     (0:neq)         ,   !Declara que o vetor inicia em 0
     &                fi     (0:neq)         ,
     &                upred  (0:neq)         ,
     &                vpred  (0:neq)         , 
     &                lm     (nume,nd)
      
      dimension       incid  (nume,4)        ,
     &                x      (numnp)         ,
     &                y      (numnp)         ,
     &                z      (numnp)         ,
     &                mtype  (nume)          , 
     &                prop   (nummat,ncprop)
      
      dimension       ske      (21)  ,
     &                sme      (21)  ,
     &                lmaux    (6)   ,
     &                sk       (6,6) ,
     &                sm       (6,6) ,
     &                uepred   (6)   ,
     &                vepred   (6)   ,
     &                fe1_pred (6)   ,
     &                fe2_pred (6)   ,
     &                fe_pred  (6) 
      
      rewind ielmnt

	do iel = 1, nume
          
          ske(:) = 0.d0
          sme(:) = 0.d0

	    x21 = x (incid(iel,2)) - x (incid(iel,1))
	    y12 = y (incid(iel,1)) - y (incid(iel,2))

	    x13 = x (incid(iel,1)) - x (incid(iel,3))
	    y31 = y (incid(iel,3)) - y (incid(iel,1))

	    x32 = x (incid(iel,3)) - x (incid(iel,2))
	    y23 = y (incid(iel,2)) - y (incid(iel,3))

	    area2 = x (incid(iel,1)) * y23
     .          + x (incid(iel,2)) * y31
     .          + x (incid(iel,3)) * y12


	    if ( area2 .le. 0.d0 ) then
              write (iout,300) iel
	        stop
	    endif

c       estado plano de deformacoes (thick = 1.) 
c
c       area2    = 2 * area
c       volume   = 0.5 * area2 * thick
c       det(jac) = area2
c       vjac2    = volume / det(jacobiano)**2
c
c                   | y23  y31  y12 |
c       B = 1/area2 |               |  = dNi/dXj  (i=1,2,3 ; j=1,2)
c                   | x32  x31  x21 |


          thic = prop ( mtype(iel),3 )

	    vjac2 = 0.5d0 / area2 * thic

c       matriz constitutiva D

	    dd1 = prop ( mtype(iel),6 )
	    dd2 = prop ( mtype(iel),7 )
	    dd3 = prop ( mtype(iel),8 ) 
	  
c       matriz de rigidez de elemento  int(BtDB.dV)

	    ske (1)  = vjac2 * ( y23*dd1*y23 + x32*dd3*x32 )
	    ske (2)  = vjac2 * ( y23*dd2*x32 + x32*dd3*y23 )
	    ske (4)  = vjac2 * ( y23*dd1*y31 + x32*dd3*x13 )
	    ske (7)  = vjac2 * ( y23*dd2*x13 + x32*dd3*y31 )
	    ske (11) = vjac2 * ( y23*dd1*y12 + x32*dd3*x21 )
	    ske (16) = vjac2 * ( y23*dd2*x21 + x32*dd3*y12 )

	    ske (3)  = vjac2 * ( x32*dd1*x32 + y23*dd3*y23 )
	    ske (5)  = vjac2 * ( x32*dd2*y31 + y23*dd3*x13 )
	    ske (8)  = vjac2 * ( x32*dd1*x13 + y23*dd3*y31 )
	    ske (12) = vjac2 * ( x32*dd2*y12 + y23*dd3*x21 )
	    ske (17) = vjac2 * ( x32*dd1*x21 + y23*dd3*y12 )

	    ske (6)  = vjac2 * ( y31*dd1*y31 + x13*dd3*x13 )
	    ske (9)  = vjac2 * ( y31*dd2*x13 + x13*dd3*y31 )
	    ske (13) = vjac2 * ( y31*dd1*y12 + x13*dd3*x21 )
	    ske (18) = vjac2 * ( y31*dd2*x21 + x13*dd3*y12 )

	    ske (10) = vjac2 * ( x13*dd1*x13 + y31*dd3*y31 )
	    ske (14) = vjac2 * ( x13*dd2*y12 + y31*dd3*x21 )
	    ske (19) = vjac2 * ( x13*dd1*x21 + y31*dd3*y12 )

	    ske (15) = vjac2 * ( y12*dd1*y12 + x21*dd3*x21 )
	    ske (20) = vjac2 * ( y12*dd2*x21 + x21*dd3*y12 )

	    ske (21) = vjac2 * ( x21*dd1*x21 + y12*dd3*y12 )
          
c       massa do elemento (area*thic*rho)

          thic = prop ( mtype(iel),3 )
	    dens = prop ( mtype(iel),4 )
          
          emass = 0.5d0 * area2 * thic * dens !m esta implicito que é enteiro
          emass_3 = emass/3.d0
	  
c       matriz de massa de elemento

	    sme (1)  = emass_3
          sme (3)  = emass_3
          sme (6)  = emass_3
          sme(10)  = emass_3
          sme(15)  = emass_3
          sme(21)  = emass_3
 
c      montando matriz quadrada          
          
          kk = 0                             
          do jj = 1 ,nd
              do ii = 1 ,jj
                kk = kk + 1
                sk(ii,jj) = ske(kk)
                sk(jj,ii) = sk(ii,jj)
                
                sm(ii,jj) = sme(kk)
                sm(jj,ii) = sm(ii,jj)
              enddo
          enddo 
          
c      mapeando deslocamentos e velocidades preditos

          do ind =1, nd
	        ieq = lm(iel,ind) 
              if(ieq.gt.0) then
	            uepred(ind) = upred(ieq)
                  vepred(ind) = vpred(ieq)
	        endif
          enddo

c     calculando forca interna

        do ind = 1, nd
           do jnd = 1, nd
              fe1_pred(ind)= fe1_pred(ind) + sk(ind,jnd) * ( uepred(jnd)
     &                       + gamman*vepred(jnd) )
               
           enddo
        enddo
   
   
        do ind = 1, nd
           do jnd = 1, nd
              fe2_pred(ind)= fe2_pred(ind) + betan*sm(ind,jnd)
     &                       *vepred(jnd)
           enddo
        enddo   
   
        fe_pred(:) = fe1_pred(:) + fe2_pred(:)
        
c     assemblando froca interna
        
        do ind = 1, nd                   !assemblando
          ieq = lm(iel,ind)
          if (ieq.gt.0) then
            fi(ieq) = fi(ieq) + fe_pred(ind)  
          endif
        enddo
          
      enddo
          
          
  300 format (' *** (TRIE2D) Area nao positiva p/ o elemento (',i5,')')
      
      end