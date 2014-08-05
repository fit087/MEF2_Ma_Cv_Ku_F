	subroutine trie2d (stiff,maxa,lm,x,y,incid,mtype,prop,
     .                   numnp,nume,nummat,nwk,ndt,neq1)
	implicit real*8 (a-h,o-z)
	include 'tapes.h'

c       matriz de rigidez - triangulo linear - elasticidade plana (epd)
c
c       ske =   s1      s2      s4      s7      s11     s16
c                       s3      s5      s8      s12     s17     
c                               s6      s9      s13     s18
c                                       s10     s14     s19
c                       simm.                   s15     s20
c                                                       s21
c

	dimension         stiff   (nwk)           ,
     .                  lm      (nume,ndt)       ,
     .                  maxa    (neq1)

	dimension         x       (numnp)          ,
     .                  y       (numnp)          ,
     .                  incid   (nume,4)        ,
     .                  mtype   (nume)          ,
     .                  prop    (nummat,8)       ,
     .                  ske     (21)            ,
     .                  lmaux   (6)              
     


	call bot ('trie2d')

	rewind ielmnt

	do iel = 1, nume

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

	    do j = 1, ndt
	        lmaux(j) = lm (iel,j)
          enddo 
	          
	    call addban2 (stiff,maxa,ske,lmaux,ndt)
	  
	
	                                   
c	    kk = 0                             
c          do jj=1,6
c          do ii=1,jj
c              kk = kk + 1
c              sx(ii,jj) = ske(kk)
c              sx(jj,ii) = sx(ii,jj)
c          enddo
c          enddo 	 
c          
c          write (iout,'(a,t10,6i12)') 'lm',lmaux
c                                                      
c          write (iout,fm1) 'sx',((sx(ii,jj),jj=1,nd),
c     .                                        ii=1,nd)
c	    do ii=1,6
c	       ieq = lm(iel,ii)
c	        if (ieq.gt.0) then
c	        do jj = 1,6
c	            jeq = lm(iel,jj)
c	            if (jeq.gt.0) then
c	            sg(ieq,jeq) = sg(ieq,jeq) + sx (ii,jj)
c	            endif
c	         enddo
c	       endif
c	     enddo          
c	   
c	 
c         write (iout,fm2) 'sg',((sg(ii,jj),jj=1,neq),
c     .                                        ii=1,neq)
	     
	    write (ielmnt) ske    

	enddo

	rewind ielmnt

	call eot ('trie2d')
c	deallocate (sg,stat=ierr)
	return

  300 format (' *** (TRIE2D) Area nao positiva p/ o elemento (',i5,')')

	end




 