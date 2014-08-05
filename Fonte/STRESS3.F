c       Calculo das Tensoes para tria3

	subroutine stress3 (lm,x,y,incid,mtype,prop,
     .                    numnp,nume,nummat,nwk,ndt,disp,neq,
     &                    tx,ty,txy,tz,s1,s2,
     .                    ex,ey,exy,ez,
     &                    neplot,ngeplot,nnplot,ngnplot,keyopt)
      
	implicit real*8 (a-h,o-z)
	include 'tapes.h'
	


	dimension      	  b       (4,8)            ,
     .                  lm      (nume,ndt)       ,
     .	              disp    (0:neq)          ,
     .                  Tx      (nume)           ,
     .                  Ty      (nume)           ,
     .                  Txy      (nume)          ,
     .                  Ex       (nume)          ,
     .                  Ey       (nume)          ,
     .                  Exy      (nume)          ,
     .                  S1      (nume)           ,
     .                  S2      (nume)           ,
     .                  Ez      (nume)           ,
     .                  Tz      (nume)           

	dimension         x       (numnp)          ,
     .                  y       (numnp)          ,
     .                  incid   (nume,4)         ,
     .                  mtype   (nume)           ,
     .                  prop    (nummat,13)        
     
   
	call bot ('stres3')
      
	 EN=0.d0
      
       write (iout,200) 
        
      s_xz=0.d0
      s_yz=0.d0
      e_xz=0.d0
      e_yz=0.d0      
      
	do iel = 1, nume
          
          no1 = incid(iel,1)
          no2 = incid(iel,2)
          no3 = incid(iel,3)
          
          y23 = y(no2) - y(no3)
          y31 = y(no3) - y(no1)
          y12 = y(no1) - y(no2)
          
          x32 = x(no3) - x(no2)
          x13 = x(no1) - x(no3)
          x21 = x(no2) - x(no1)
         
          area2 = x (no1) * y23
     .          + x (no2) * y31
     .          + x (no3) * y12

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

c       matriz constitutiva D

	    dd1 = prop ( mtype(iel),6 ) 
	    dd2 = prop ( mtype(iel),7 ) 
	    dd3 = prop ( mtype(iel),8 ) 
	  
	  poiss  = prop ( mtype(iel),2 )
   	  eyoung = prop ( mtype(iel),1 )

c       matriz B
	          
      Area2_inv= 1.d0 /  area2

	b11 = y23 * area2_inv
	b12 = y31 * area2_inv
	b13 = y12 * area2_inv
	
	b21 = x32 * area2_inv
	b22 = x13 * area2_inv
	b23 = x21 * area2_inv

c      Deslocamentos ui do elemento iel 

	u1= disp ( lm( iel , 1) )
	v1= disp ( lm( iel , 2) )

	u2= disp ( lm( iel , 3) )
	v2= disp ( lm( iel , 4) )
	
	u3= disp ( lm( iel , 5) )
	v3= disp ( lm( iel , 6) )

c     deformacoes infinitesimais
          
          e_xx   =   b11*u1 + b12*u2 + b13*u3
          e_yy   =   b21*v1 + b22*v2 + b23*v3
          e_xy   =   b21*u1 + b22*u2 + b23*u3
     &           +   b11*v1 + b12*v2 + b13*v3
          e_zz   =   0.d0
       
c      Tensoes Tx , Ty , Txy

          s_xx = dd1 * e_xx + dd2 * e_yy
          s_yy = dd2 * e_xx + dd1 * e_yy
          s_xy = dd3 * e_xy
          s_zz = 0.d0
          
          if (keyopt.eq.1) s_zz =   poiss * (s_xx + s_yy)            ! tensao z para EPD
          if (keyopt.eq.2) e_zz = - poiss * (s_xx + s_yy) / eyoung   ! deformaçao z para EPT
 
          tx  (iel) = s_xx
          ty  (iel) = s_yy
          txy (iel) = s_xy
          tz  (iel) = s_zz
          
          ex  (iel) = e_xx
          ey  (iel) = e_yy
          exy (iel) = e_xy
          ez  (iel) = e_zz
          
c      sxx syy szz sxy sxz syz
          
          write (iout,300) iel, s_xx,s_yy,s_zz,s_xy,s_xz,s_yz
          write (iout,400)      e_xx,e_yy,e_zz,e_xy,e_xz,e_yz
          
!	Tx(iel)= area2_inv * (  y23*dd1*u1 + x32*dd2*v1 + y31*dd1*u2 + 
!     &   	                    x13*dd2*v2 + y12*dd1*u3 + x21*dd2*v3 )
!	
!	
!	Ty(iel)= area2_inv * (  y23*dd2*u1 + x32*dd1*v1 + y31*dd2*u2 + 
!     &   	                    x13*dd1*v2 + y12*dd2*u3 + x21*dd1*v3 )
!
!	
!	Txy(iel)= area2_inv * ( x32*dd3*u1 + y23*dd3*v1 + x13*dd3*u2 + 
!     &   	                    y31*dd2*v2 + x21*dd1*u3 + y12*dd2*v3 )
!
!
!
!	Ex(iel) = area2_inv * ( y23 * u1  + y31 * u2 + y12 * u3 )
!	Ey(iel) = area2_inv * ( x32 * v1  + x13 * v2 + x21 * v3 )
!	Exy(iel)= area2_inv * ( x32 * v1  + x13 * v2 + x21 * v3 + 
!     &	                    y23 * u1  + y31 * u2 + y12 * u3 )
!


c      Calculo da Energia

!      EN=EN + ( area2_inv * Tx(iel) * ( u1 * y23 + u2 * y31 + u3 * y12)
!     .	  +   area2_inv * Ty(iel) * ( v1 * x32 + v2 * x13 + v3 * x21)
!     .      +   area2_inv * Txy(iel) * (  v1 * y23 + v2 * y31 + v3 * y12
!     .      +   u1 * x32 + u2 * x13 + u3 * x21 ) )


c      Tensoes Principais

      S1(iel) = ( Tx(iel) + Ty(iel) )*0.5d0 + 
     .sqrt( ( (Tx(iel)-Ty(iel))* 0.5d0 )*((Tx(iel) - Ty(iel))* 0.5d0) +
     .       Txy(iel)*Txy(iel))
	
      S2(iel) = ( Tx(iel) + Ty(iel) )*0.5d0 - 
     .sqrt( ( (Tx(iel)-Ty(iel))* 0.5d0 )*((Tx(iel) - Ty(iel))* 0.5d0) +
     .       Txy(iel)*Txy(iel) )

	
	  
!	if (keyopt.eq.1) then ! Caso EPD
!
!	Tz(iel)= poiss * ( Tx(iel) + Ty(iel) )
!
!      elseif (keyopt.eq.2) then ! Caso EPT
!          
!      Tz(iel) = 0.d0
!
!	Ez(iel)= - (poiss / E)*( Tx(iel) + Ty(iel) )
!	
!
!      endif  

      enddo
      

      

c     Plota Tensao no Elemento definido

	if (ngeplot.eq.1) then

	write (iesaida,*)Tx  (lm( neplot , 1))
	
	write (iesaida,*)Ty  (lm( neplot , 1))

	write (iesaida,*)Txy (lm( neplot , 1))

	else if (ngeplot.eq.2) then

	write (iesaida,*)Tx  (lm( neplot , 2))
	
	write (iesaida,*)Ty  (lm( neplot , 2))

	write (iesaida,*)Txy (lm( neplot , 2))

	endif


	write(iout,*)
	write(iout,*)"Norma da Energia = ",EN

	rewind ielmnt

	call eot ('stres3')
	return
      
c                                     sxx syy szz sxy sxz syz
      
200     format (////,
     &          ' Tensoes      iel  sxx syy szz sxy sxz syz',/,
     &          ' Deformacoes       exx eyy ezz exy exz eyz',/)
300       format (i5,t10,6(1pe15.8,1x) )
400       format (   t10,6(1pe15.8,1x)/)
          
 	end




 