      subroutine solver_newmark  (x,y,z,incid,prop,mtype,lm,
     &             stiff,maxa,f,f0,fi,
     &             du,nwk,id,upred,vpred,ucorr,vcorr,
     &             nsteps,nflag,timef,dt,betan,gamman,damp1,damp2)
c

c
c     jlda	feb-1994
c	resumed	jan-1995/2014
c

      implicit real*8 (a-h,o-z)
      include 'common.h'
      include 'cntl.h'
      include 'tapes.h'
      include 'pointers.h'
                                                  
      dimension       x       (numnp)         ,
     &                y       (numnp)         ,
     &                z       (numnp)         ,
     &                incid   (nume,nnoel)    ,
     &                prop    (nummat,8)      ,
     &                mtype   (nume)          ,
     &                lm      (nume,ndt)      , 
     &                stiff   (nwk)           ,
     &                f       (0:neq)         ,
     &                f0      (0:neq)         ,
     &                fi      (0:neq)         ,
     &                du      (0:neq)         ,
     &                upred   (0:neq)         ,
     &                vpred   (0:neq)         ,
     &                ucorr   (0:neq)         ,
     &                vcorr   (0:neq)         ,
     &                an      (0:neq)         , 
     &                id	(ngl,numnp) 
      
      Dimension
     .                  tx      (nume)        ,
     .                  ty      (nume)        ,
     .                  txy      (nume)       ,
     .                  ex       (nume)       ,
     .                  ey       (nume)       ,
     .                  exy      (nume)       ,
     .                  s1       (nume)       ,
     .                  s2       (nume)       ,
     .                  s3       (nume)       ,
     .                  ez      (nume)        ,
     .                  tz      (nume)                  
     
      character*1   flgs
      logical       fimarq /.false./

      call bot ('solver') 

      neq1 = neq + 1
      icount = 0
      ucorr(:) = 0.d0
      vcorr(:) = 0.d0
      an(:) = 0.d0
      
      call effectmass_trie2d (stiff,fp,up,maxa,lm,x,y,z,incid,mtype,
     &                 prop,numnp,nume,nummat,nwk,ndt,neq+1,
     &                 neqp,nnoel,8,damp1,damp2,betan,gamman,
     &                 dt,neq,ndaux)

      DO istep = 1,nsteps
          
          call preditor(neq,dt,betan,gamman,ucorr,vcorr,an,upred,vpred)
          call trie2d_fint (stiff,fp,up,maxa,lm,x,y,z,incid,
     &                      mtype,prop,numnp,nume,nummat,nwk,ndt,
     &                      neq+1,neqp,nnoel,8,
     &                      damp1,damp2,betan,gamman,dt,neq,
     &                      fi,upred,vpred)
          
          DO ieq=1,neq
              an(ieq)=f0(ieq)-f(ieq)-fi(ieq)
          ENDDO
          
          call COLSOL(stiff,an,maxa,neq,nwk,neq+1,2,ierror)
          
          call corretor(neq,dt,betan,gamman,upred,vpred,an,ucorr,vcorr,)
          
          
          
          IF (mod(istep,nflag).eq.0) THEN
              icount = icount + 1
              CALL ENSG_VEC('dis',ucorr,id,icount)
          ENDIF
          
      ENDDO
            
	call stress3 (ia(nlm), ia(nx), ia(ny), ia(ninc), ia(nmtype), ia(nprop),
     .              numnp,nume,nummat,nwk,ndt,f0,neq,
     &              tx,ty,txy,tz,s1,s2,
     &              ex,ey,exy,ez,
     .              neplot,ngeplot,nnplot,ngnplot,s1,s2,keyopt)
	

	CALL ENSG_VEC('dis',f0,id,0)


c      Saida das Tensoes

	call ensg_esca ('txx',tx,0,nume)
	call ensg_esca ('tyy',ty,0,nume)
	call ensg_esca ('txy',txy,0,nume)
	call ensg_esca ('tzz',tz,0,nume)
      
c      Saida dos E
	
	call ensg_esca ('exx',ex,0,nume)
	call ensg_esca ('eyy',ey,0,nume)
	call ensg_esca ('exy',exy,0,nume)
	call ensg_esca ('ezz',ez,0,nume)

c      Saida das Tensoes Principais
 
      call ensg_esca ('S11',s1,0,nume)
      call ensg_esca ('S22',s2,0,nume)
      call ensg_esca ('S33',0.d0,0,nume)

      call ENSG_tens ('str',tx,ty,tz,txy,0,nume)

      call eot ('solver')
      return

  100 format (///,' Iteracao Newton-Raphson Modificada ',/,
     .            ' Numero de Incrementos (nincre) : ',i5,/)
  200 format ( /, ' Incremento # ',i5,//,
     .            '  iter        rel.rnorm           rel.unorm ',
     .            '        ns  kp',/)
  300 format ( '  ',i5,1p,2e20.8,1x,2i5)
  400 format (    '   Convergencia alcancada ',//)

  500 format ( /, ' *** Convergencia nao alcancada          *** ',/,
     &            ' *** Incremento ',i5,' em ',i5,' iteracoes *** ',/,
     .            ' *** Execucao interrompida               *** ',/)

      end


!PREDITOR
!DO ieq=1,neq
!      upred = ucorr + vcorr*dt + (0.5d0 * dt**2) * (1.d0-2.d0*betan) * an
!      vpred = vcorr + dt * (1.d0-gamman) * an
!ENDDO
!
!
!CORRETOR
!DO ieq=1,neq
!      ucorr = upred + betan * dt**2 * an
!      vcorr = vpred + gamman * dt * an
!ENDDO