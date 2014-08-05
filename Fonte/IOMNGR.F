C    -------------------
      SUBROUTINE IOMNGR  ( INT )
C    -------------------
C
	include	'tapes.h'
	common  /iolist/        ntm,ntr,nin,not,nsp,nfl,nt7,nt8
	
	character*70         filein     ,
     .                     fileesaida ,
     .                     filensaida ,
     .                     fileout    ,
     .                     fileplt    ,
     .                     filelog    ,
     .                     filecon    ,
     .                     filedbg    ,
     .                     filedef    ,
     .                     filexgf

      integer*4            count	

      IF ( INT .EQ. 0 ) THEN 
      
c     checa linha de comandos por parametros
      
         count = nargs()
         if (count.gt.1) then
         
            call getarg (1,filer,iblank)
            iblank = iblank + 1
            ibacks = index (filer,'\',.true.)
            filein = filer(1:iblank-1)//'.dat'
         
         else   

C --- DEFINICAO DAS UNIDADES LOGICAS PARA I/O

	      write (*,'(a)') ' Arquivo de DADOS   ?'
	      read  (*,'(a)')  filein
	
	      write (*,'(a)') ' Radical das SAIDAS ?'
            read  (*,'(a)')  filer
	
            iblank=index(filer,' ')
            ibacks=index(filer,'\',.true.)

      endif                      

c     monta os nomes dos arquivos

          fileout     =   filer(1:iblank-1)//'.lis'
	    fileesaida  =   filer(1:iblank-1)//'.esai'
	    filensaida  =   filer(1:iblank-1)//'.nsai'
          fileplt     =   filer(1:iblank-1)//'.plt'
	    filelog     =   filer(1:iblank-1)//'.log'
          filecon     =   filer(1:iblank-1)//'.con'
          filedbg     =   filer(1:iblank-1)//'.dbg' 
          filedef     =   filer(1:iblank-1)//'.def'      
          filexgf     =   filer(1:iblank-1)//'.xgf'
      
      write (*,*) fileout
	write (*,*) fileesaida
      write (*,*) filensaida
      write (*,*) fileplt
      write (*,*) filelog
      write (*,*) filecon
      write (*,*) filedbg
      write (*,*) filedef
      write (*,*) filexgf
	

      IWORK1 = 0
      IWORK2 = 0
      IWORK3 = 0
      IELDAT = 0

      IIN    = 1
      IOUT   = 2
      IELMNT = 3
      IARQ   = 4
      ILOAD  = 12
      

      ICON   = 10
      ILOG   = 11
      IPLT   = 20
      IDEF   = 21
      IASF   = 22
      IDBG   = 30
      
 	ntm = ilog
	not = iout

C --------------------------------------------------------------------


      ibuff = 1

      OPEN (UNIT=IIN ,FILE=filein,
     &        BLANK='ZERO',mode='READ',status='OLD')
c     OPEN (UNIT=IOUT     ,FILE=fileout,fileopt='BUFFER=10'  )
      OPEN (UNIT=IOUT     ,FILE=fileout,blocksize=ibuff  ) 
      OPEN (UNIT=IPLT     ,FILE=fileplt  )
	OPEN (UNIT=INSAIDA  ,FILE=filensaida  )
	OPEN (UNIT=IESAIDA  ,FILE=fileesaida  )
      
	open (unit=ielmnt,status='SCRATCH',form='BINARY')
      open (unit=iload ,status='SCRATCH',form='BINARY')
c     open (unit=ielmnt,status='SCRATCH',form='UNFORMATTED')
      
      OPEN (UNIT=ICON     ,FILE=filecon)
      OPEN (UNIT=ILOG     ,FILE=filelog,blocksize=ibuff)
c     OPEN (UNIT=ILOG,FILE=filelog,fileopt='BUFFER=10')
      OPEN (UNIT=IDBG     ,FILE=filedbg,blocksize=ibuff  )
      OPEN (UNIT=IDEF     ,FILE=filedef  )
      OPEN (UNIT=IASF     ,FILE=filexgf  )

      ELSE

		CLOSE (UNIT=IIN   )
		CLOSE (UNIT=IOUT  )
		CLOSE (UNIT=ICON  )
		CLOSE (UNIT=IDBG  )
		CLOSE (UNIT=IPLT  )
		CLOSE (UNIT=IDEF  )
		CLOSE (UNIT=IASF  )
		close (unit=ielmnt)
		close (unit=iload)
	    close (unit=INSAIDA)
	    close (unit=IESAIDA)

      ENDIF

      RETURN
      END
