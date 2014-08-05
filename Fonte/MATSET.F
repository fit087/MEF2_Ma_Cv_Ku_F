	subroutine matset (prop,nummat,ncprop,keyopt)
	implicit real*8 (a-h,o-z)

c     leitura de propriedades do material

      real*8  prop(nummat,ncprop)
	include 'tapes.h'

	call bot ('matset')
	
	conv = 4.d0*datan(1.d0)/180.d0
	
	do i = 1, nummat

      read  (iin,100)  n,eyoung,poiss,thic,dens
      
      if (keyopt.eq.1) thic = 1.0d0
		  
      bulk = eyoung/(3.d0*(1.d0-2.d0*poiss))
      
      write (iout,101)   n,eyoung,poiss,thic,dens
      
      if (n.gt.nummat.or.n.le.0) then
         write (iout,300) n
         stop
	endif

c     armazena propriedades do material
		
		prop(n,1) = eyoung
		prop(n,2) = poiss
		prop(n,3) = thic
		prop(n,4) = dens
 		prop(n,5) = 0.d0

c	tensor de elasticidade plana (EPD) CASE 1  ( d1 d2 d3 )

	if (keyopt.eq.1) then 
          
		prop(n,6) = eyoung * (1.d0 - poiss)
     .			/ ((1.d0 + poiss) * (1.d0 - 2.d0 * poiss))
		prop(n,7) = eyoung  * poiss 
     .   / ((1.d0 + poiss) * (1.d0 - 2.d0 * poiss))
		prop(n,8) = 0.5d0 * eyoung / (1.d0 + poiss)
	
c     		
		 
          prop(n,9)  = 0.d0
	    prop(n,10) = 0.d0
	    prop(n,11) = 0.d0              !   dsin(phi*conv)
	    prop(n,12) = 0.d0
	    prop(n,13) = 0.d0              !   cohes * dcos(phi*conv)  

c     calculo de D para Estado Plano de Tensao (EPT) CASE 2
	
	elseif (keyopt.eq.2) then
	
		prop(n,6) = eyoung / (1.d0 - poiss**2)
		prop(n,7) = eyoung * poiss / (1.d0 - poiss**2)
		prop(n,8) = 0.5d0 * eyoung * (1.d0 - poiss) / (1.d0 - poiss**2)
	
	endif


	enddo
     

	call eot ('matset')
	return

  100   format (i5,4f10.0)
  101	format (/,' Propriedades do Material: mat #',i5,//,
     &           t10,'Modulo de Elasticidade   (eyoung ) = ',1pe15.8,/,
     &           t10,'Coeficiente de Poisson   (poisson) = ',1pe15.8,/,
     &           t10,'Espessura                (thic   ) = ',1pe15.8,/,
     &           t10,'Massa específica         (dens   ) = ',1pe15.8,/)
     
     
  200 format (' ',i8,1p,9e15.7)
  300	format (' *** (MATSET) Erro de Dados: grupo de material ',
     .		'inválido (',i8,')')
	end





