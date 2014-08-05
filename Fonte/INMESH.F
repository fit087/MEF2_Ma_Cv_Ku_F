	subroutine inmesh
	
c     leitura da malha
	
	implicit real*8 (a-h,o-z)

	include 'common.h'
	include 'pointers.h'
	include 'cntl.h'

	call bot ('inmesh')
      
      nnoel = 3

	call defini ('id      ',nid,numnp,ngl)
	call define ('x       ',nx ,numnp,1  )
	call define ('y       ',ny ,numnp,1  )
	call define ('z       ',nz ,numnp,1  )

	ndt = nnoel * ngl

c     coordenadas e condicoes de contorno

	call joints (ia(nid),ia(nx),ia(ny),ia(nz))

c     numeracao das equacoes
	
	call equat  (ia(nid))

c     propriedades dos materiais

      ncprop = 13      
	call define ('prop    ',nprop,nummat,ncprop)   

      call matset (ia(nprop),nummat,ncprop,keyopt)

c	nd    = nnoel * ngl



	call defini ('incid   ',ninc,nume,nnoel )
	call defini ('mtype   ',nmtype,nume,1 )
	call defini ('lm      ',nlm,nume,ndt  )

c     definicao dos elementos

	call inelem (ia(nid),ia(ninc),ia(nmtype),ia(nlm),
     .			nume,nnoel,ndt,ngl)

	incre = 0
	call ensg_geo (ia(nx),ia(ny),ia(nz),ia(ninc),incre)

	call eot ('inmesh')
	return
	end

     


