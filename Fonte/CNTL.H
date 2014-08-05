      character*80 	hed
      dimension 	idof(6)

      logical conver  ,
     &        view3d  ,
     &        ensight ,
     &        xgraph  ,
     &        desloc  ,
     &        tensao  ,
     &        flgout

      common /cntl/	idof,method,numnp,nume,nummat,nedges,ncarg,
     &			    ngl,neq,neqp,ndt,maxite,
     &              nnoel,nincre,conver,en,
     &              view3d,ensight,xgraph,desloc,tensao,flgout,
     &              ixno,ixgl,fxv3d,fyv3d,hed,ntria,nquad,keyopt,
     &              neplot,ngeplot,nnplot,ngnplot
	 
	 

