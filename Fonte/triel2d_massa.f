	subroutine trie2d_massa (massa,maxa,lm,x,y,incid,mtype,prop,
     .                         numnp,nume,nummat,nwk,ndt,neq1)
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

	dimension         massa   (nwk)           ,
     .                  lm      (nume,ndt)       ,
     .                  maxa    (neq1)

	dimension         x       (numnp)          ,
     .                  y       (numnp)          ,
     .                  incid   (nume,4)        ,
     .                  mtype   (nume)          ,
     .                  prop    (nummat,8)       ,
     .                  ske     (21)            ,
     .                  lmaux   (6)              
     


	call bot ('trie2d_massa')

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

c       massa do elemento (area*thic*rho)

          thic = prop ( mtype(iel),3 )
	    dens = prop ( mtype(iel),4 )
          
          mass = 0.5d0 * area2 * thic * dens
          mass_3 = mass/3.d0
	  
c       matriz de massa de elemento

	    ske (1)  = mass_3
          ske (3)  = mass_3
          ske (6)  = mass_3
          ske(10)  = mass_3
          ske(15)  = mass_3
          ske(21)  = mass_3

	    do j = 1, ndt
	        lmaux(j) = lm (iel,j)
          enddo 
	          
	    call addban2 (massa,maxa,ske,lmaux,ndt)
	  
         !do i=1,nwk
         !    damp(i) = alphad*stiff(i) + betad*massa(i)
         !enddo
	
	                                   
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




 