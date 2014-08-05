	subroutine equat (id)
	
	parameter   (   is_fixed        =   1   ,
     .                is_free         =   0   ,
     .                is_prescribed   =  -1  )

      include     'tapes.h'
      include     'cntl.h'       
      
	dimension 	 id	(ngl,numnp)  ,
     .            idaux(6)        
      character*5 dofnam(6)
      
      data dofnam /'    u','    v','    w','   ru','  rv','   rw'/

c     numeracao das equacoes
c     	     
c     neq  = numero de equacoes livres
c     neqp = numero de equacoes prescritas (nao nulas)	   
      
      call bot ('equat ')
	     
	neq  = 0
	neqp = 0
	
	do i = 1, numnp
	do j = 1, ngl

c		select case( id(j,i) )
		   
c		case (is_fixed)
		if (id(j,i).eq.is_fixed) then
          
          		id(j,i) = 0
              
c		case (is_free)  
		elseif (id(j,i).eq.is_free) then            		   
          
			neq = neq + 1
			id(j,i) = neq
			      
c		case (is_prescribed)
		elseif (id(j,i).eq.is_prescribed) then
          
			neqp = neqp + 1
			id(j,i) = -neqp
              
c		end select
		endif
          
	enddo
	enddo 
	
	write (iout,100) neq,neqp,dofnam
      
  100 format (//,' Mapa de Equacoes:',//,
     .           ' equações livres     (neq ):',5x,i5,/,
     .           ' equações prescritas (neqp):',5x,i5,//
     .           t5,'   no',5x,6a5,//)
      
  200 format (t5,i5,5x,6i5)
  
      do igl = 1, 6
       idaux(igl) = 0
      enddo
       
      do no = 1, numnp
          
          do j = 1,6
                igl = idof(j)
                if (igl.ne.0) idaux(j) = id(igl,no)
          enddo
          
          write (iout,200) no, idaux
          
      enddo      
                    
      call eot ('equat ')
	return
	end




