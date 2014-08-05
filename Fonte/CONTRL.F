	subroutine contrl 
	implicit real*8 (a-h,o-z)
	
	include     'cntl.h'
	include     'tapes.h'

	call bot ('contrl')

	read (iin,100) hed
	read (iin,120) numnp,nume,nummat,ncarg,keyopt
	
	view3d  = .false.
	ensight = .false.
	xgraph  = .false.
	desloc  = .false.
	tensao  = .false.  

	ixno  = 0
	ixgl  = 0
	fxv3d = 1.d0
	fyv3d = 1.d0
	
	read (iin,'(5l1)')      view3d,ensight,xgraph,displac,stress 
      read (iin,'(2f10.0)')   fxv3d,fyv3d

	call setdof (idof,ngl)

	write (iout,200) hed
	
	if (keyopt.eq.1) then

	write(iout,*)"Estado Plano de Deformacao"

	elseif (keyopt.eq.2) then
	
	write(iout,*)"Estado Plano de Tensao"

	endif
	
	write (iout,220) idof,numnp,nume,nummat,ncarg,keyopt
	write (iout,230) view3d,ensight,xgraph,desloc,tensao
 	write (iout,240) ixno,ixgl
	write (iout,250) fxv3d,fyv3d
      
	call eot ('contrl')

      return

  100 format (a)
  110 format (6i1)
  120 format (16i5)
  130 format (2i8)

  200 format (' Title : ',/,' ',a)
  220 format (/,' Dados de Controle ',/,
     &            '  Graus de Liberdade       (idof)   : ',6i8,/,
     &            '  Numero de Nos            (numnp)  : ',i8,/,
     &            '  Numero de Elementos      (nume)   : ',i8,/,
     &            '  Numero de Materiais      (nummat) : ',i8,/,
     &            '  Numero de Nos Carregados (ncarg)  : ',i8,/
     &            '  Opcao (EPD=1; EPT=2)     (keyopt) : ',i8,// )     
  230 format (/,' Geracao de Saidas ',/,
     &            '  View3D         : ',l1,/,
     &            '  Ensight        : ',l1,/,
     &            '  Xgraph         : ',l1,/,
     &            '  Deslocamentos  : ',l1,/,
     &            '  Tensoes        : ',l1,//)
  240 format (/,' Saida Xgraph (no, g.l.) ',/,
     &            '  no             : ',i5,/,
     &            '  g.liberdade    : ',i5,//)
  250 format (/,' Fator de Deformada View3D ',/,
     &            '  fxv3d          : ',f10.5,/,
     &            '  fyv3d          : ',f10.5,//)
	end
