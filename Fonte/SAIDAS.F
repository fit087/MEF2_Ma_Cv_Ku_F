      subroutine ensg_geo (x,y,z,incid,istep)
      
      implicit real*8 (a-h,o-z)
      
      include 'tapes.h'
      include  'cntl.h'
      
      dimension   x       (numnp)     ,
     &            y       (numnp)     ,
     &            z       (numnp)     ,
     &            incid   (nume,nnoel)
     
      
      character*70 filename
                       
c     arquivo Ensight para a geometria do modelo
      
      write (filename(1:3),'(i3.3)') istep
      
      filename = filer(1:iblank-1)//'.geo'//filename(1:3)
      
      write(*,'(a)') ' Gravando '//filename(1:50)
                                  
      open (unit=iplt,file=filename,form='formatted')                 
                       
      write (iplt,'(a)') 'EPD'
      write (iplt,'(a)') 'Ensight Geometry COPPE/UFRJ'
      write (iplt,'(a)') 'node id given'
      write (iplt,'(a)') 'element id given'
      write (iplt,'(a,/,i8)') 'coordinates',numnp
      
      write (iplt,'(i8,1p,3e12.5)') (i,x(i),y(i),z(i),i=1,numnp)
      
      write (iplt,'(a,i8)') 'part',1
      write (iplt,'(a)')    'malha'
      
      write (iplt,'(a,/,i8)') 'tria3',nume
      
      write (iplt,'(4i8)') (i,(incid(i,j),j=1,nnoel),i=1,nume)
      
      close (unit=iplt)

      return
      end


      subroutine ensg_esca (tipo,esca,istep,nelx)
      
      implicit real*8 (a-h,o-z)
      
      include 'tapes.h'
	include  'cntl.h'
      
      dimension   esca     (nelx) 
      
      character*70 filename
      character*3  tipo

c     arquivo Ensight para escalares nodais

      write (filename(1:3),'(i3.3)') istep

      filename = filer(1:iblank-1)//'.'//tipo//filename(1:3)

      write(*,'(a)') ' Gravando '//filename(1:50)
      
      open (unit=iplt,file=filename,form='formatted')
      
      write (iplt,'(a,i5,1x,a)') 'Ensight Escalar passo ',istep

      write (iplt,'(a/a)') 'part 1','tria3' 
      write (iplt,'(1p,6e12.5)') (esca(i),i=1,nelx)

      close (unit=iplt)
      
      return
      end


      subroutine ensg_vec (tipo,disp,id,istep)
      
      implicit real*8 (a-h,o-z)
      
      include 'tapes.h'
      include  'cntl.h'

      dimension   disp (0:neq)	,
     &		  id   (ngl,numnp) 
      
      character*70 filename
      character*3  tipo

c     arquivo Ensight para vetores

      write (filename(1:3),'(i3.3)') istep

      filename = filer(1:iblank-1)//'.'//tipo//filename(1:3)

      write(*,'(a)') ' Gravando '//filename(1:50)

      
      open (unit=iplt,file=filename,form='formatted')
      
      write (iplt,'(a,i5,1x,a)') 'Ensight Vector passo ',istep

      write (iplt,'(1p,6e12.5)') 
     &        (disp(id(1,i)),disp(id(2,i)),0.d0,i=1,numnp)

      close (unit=iplt)
      
      return
      end


      subroutine ensg_vec2 (tipo,vec,istep,idvec,nvec)
      
      implicit real*8 (a-h,o-z)
      
      include 'tapes.h'
      include  'cntl.h'

      dimension   vec (idvec,3)	
        
      character*70 filename
      character*3  tipo

ccc     arquivo Ensight para vetores

      write (filename(1:3),'(i3.3)') istep

      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:3)

      write(*,'(a)') ' Gravando '//filename(1:50)

      
      open (unit=iplt,file=filename,form='formatted')
      
      write (iplt,'(a,i5,1x,a)') 'Ensight Vector passo ',istep
      write (iplt,'(1p,6e12.5)') (vec(i,1),
     &                            vec(i,2),
     &                            vec(i,3),
     &				                i=1,nvec)
      
      close (unit=iplt)
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine ensg_par (xyz_par,istep,idpar,numpar)
      
      implicit real*8 (a-h,o-z)
      
      include 'tapes.h'
      
      dimension   xyz_par       (idpar,3)         
      
      character*70 filename
                       
c     arquivo Ensight para particulas
      
      write (filename(1:3),'(i3.3)') istep
      
      write(filename,'(a)') filer(1:iblank-1)//'.p'//filename(1:3)
      
      write(*,'(a)') ' Gravando '//filename(1:50)
                                  
      open (unit=iplt,file=filename,form='formatted')                 
                       
      write (iplt,'(a)') 'Ensight Particles (BETH) COPPE/UFRJ'
      write (iplt,'(a,/,i8)') 'particle coordinates',numpar
      
      write (iplt,'(i8,1p,3e12.5)') (i,(xyz_par(i,j),j=1,3),i=1,numpar)
      
        
      close (unit=iplt)

      return
      end


CHICAO 14/01/2002
      SUBROUTINE ENSG_SCA1 (TIPO,SCA,ISTEP,NUMNP, npinballs)
      
      implicit real*8 (a-h,o-z)
      
      include 'tapes.h'
      
      dimension   sca     (numnp) 
      
      character*70 filename
      character*1  tipo

c     arquivo Ensight para escalares nodais

      write (filename(1:3),'(i3.3)') istep

      write(filename,'(a)') filer(1:iblank-1)//'.'//tipo//filename(1:3)

      write(*,'(a)') ' Gravando '//filename(1:50)

      
      open (unit=iplt,file=filename,form='formatted')
      
      write (iplt,'(a,i5)') 'Ensight Scalar passo ',istep
      write (iplt,'(1p,6e12.5)') (sca(i),i=1,npinballs)
      
      close (unit=iplt)
      
      RETURN
      END
CHICAO
      

      SUBROUTINE ENSG_tens (TIPO,sx,sy,sz,txy,ISTEP,nume)
c
c
c      sxx syy szz sxy sxz syz
c
c      
      implicit real*8 (a-h,o-z)
      
      include 'tapes.h'
      
      dimension   sx(1), sy(1), sz(1), txy(1) 
      
      character*70 filename
      character*3  tipo

c     arquivo Ensight para tensores por elemento

      write (filename(1:3),'(i3.3)') istep

      filename = filer(1:iblank-1)//'.'//tipo//filename(1:3)

      write(*,'(a)') ' Gravando '//filename(1:50)

      
      open (unit=iplt,file=filename,form='formatted')
      
      write (iplt,'(a,i5)') 'Ensight Tensor passo ',istep
      write (iplt,'(a/a)') 'part 1','tria3'
           
      write (iplt,'(1p,6e12.5)') 
     &            (sx(i),sy(i),sz(i),txy(i),0.d0,0.d0,i=1,nume)
      
      
c      sxx syy szz txy txz tyz

      
      close (unit=iplt)
      
      RETURN
      END
