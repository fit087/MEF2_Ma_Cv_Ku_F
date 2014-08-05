	subroutine setdb

c       estabelece as caracteristicas da base de dados

	include		'common.h'
c
c       common.h :
c	definicao de "parameter       (maxmem)"
c	definicao de "common          mtot,np,ia(maxmem)"
c
	character*3     type
	common /dbsys/  numa,next,idir,ip(3),type(3)

c
c       mtot  = memoria total na area de trabalho
c       ia    = area de trabalho (arranjos+tabelas)
c       ip(1) = no. de bytes para inteiros
c       ip(2) = no. de bytes para reais
c       ip(3) = no. de bytes para caracter
c       numa  = no. de arranjos alocados
c       next  = proxima posicao livre no common
c       idir  = posicao do inicio da tabela de arranjos
c               sera disposta ao final do arranjo ia
c
	mtot  = maxmem
	idir  = maxmem
	numa  = 0
	next  = 1

	ip(1) = 4
	ip(2) = 8
	ip(3) = 1

	type(1) = 'int'
	type(2) = 'dbl'
	type(3) = 'chr'

	return
	end
c
	subroutine define (name,na,nr,nc)
	character*1 name(8)
      include 'common.h'
	np = 2
	call defin (name,na,nr,nc)
	return
	end
c
	subroutine defini (name,na,nr,nc)
	character*1 name(8)
      include 'common.h'
	np = 1
	call defin (name,na,nr,nc)
	return
	end
c
	subroutine definh (name,na,nr,nc)
	character*1 name(8)
      include 'common.h'
	np = 3
	call defin (name,na,nr,nc)
	return
	end
c
	subroutine defin (name,na,nr,nc)
	character*1 name(8)
	character*3 type
      include 'common.h'
	common  /dbsys/ numa,next,idir,ip(3),type(3)
	common  /iolist/        ntm,ntr,nin,not,nsp,nfl,nt7,nt8
	
c       name = nome do arranjo (8 caracteres)
c       na   = pointer para o arranjo no commom 
c       nr   = numero de linhas
c       nc   = numero de colunas
c       mtot = fim do diretorio
c       numa = numero de arranjos na base de dados
c       next = proxima posicao livre no arranjo 
c       idir = inicio do diretorio de arranjos no common
c       ip   = numero de bytes contidos no tipo de dado
c       lenr = numero de bytes no registro logico
c       np   = tipo do dado
c               1 inteiro
c               2 real
c               3 caracter
c
c       definicao dos diretorios
c       idir(1,n) = nome do arranjo - iname (8 caracteres)
c       idir(9,n) = numero de linhas    - nr
c       idir(10,n) = numero de colunas  - nc
c       idir(11,n) = tipo de dado       - np
c       idir(12,n) = posicao na memoria - na
c                 = -1 arquivo sequencial em disco
c                 = -2 arquivo direto em disco
c       idir(13,n) = tamanho do arranjo
c       idir(14,n)= 0 se estiver em memoria
c
c       definicao do diretorio para arquivos diretos
c       idir(9,n) = numero de inteiros
c       idir(10,n) = numero de palavras reais
c       idir(11,n) = numero de bytes
c       idir(12,n) = numero de registros logicos
c       idir(13,n) = numero do registro logico
c       idir(14,n) = lun se estiver na unidade logica lun
c
	if ( (np.ne.1) .and. (np.ne.3) ) np = 2
	nsize = ( nr*nc*ip(np)-1)/(ip(1)*2)
	nsize = nsize*2 + 2
	na = next
	next = next + nsize
c
	numa = numa + 1
	idir = idir - 14
	i    = idir
c
	if ( i.lt.next ) then
		i = next - i + mtot - 1
		write (ntm,2000) name,i, mtot 
		write (not,2000) name,i, mtot
		stop
	endif
c
	call icon (name,ia(i))
	ia(i+8) = nr
	ia(i+9) = nc
	ia(i+10) = np
	ia(i+11) = na
	ia(i+12) = nsize
	ia(i+13) = 0
	write (ntm,100) name,type(np),nr,nc,na,nsize
c
	return

  100   format (' (DEFIN ) ',8a1,a,' nr:',i6,' nc:',i6,
     .          ' at:',i6,' size:',i6)
 2000   format (' (DEFIN ) Memoria Esgotada em DEFIN para: ',8a1,
     .        /,'          Dimensao Necessaria  = ',i7,
     .        /,'          Dimensao Disponivel  = ',i7)
	end
c
	subroutine defdir (name,nr,nc,istr)
	character*1 name(8)
	character*3 type
      include 'common.h'
	common  /dbsys/ numa,next,idir,ip(3),type(3)
	common  /iolist/        ntm,ntr,nin,not,nsp,nfl,nt7,nt8

	if (np.eq.0) np = 2

	numa = numa + 1
	idir = idir - 14
	i = idir

	if ( i.lt.next ) then
		i = next - i + mtot - 1
		write (ntm,2000) name,i,mtot
		write (not,2000) name,i,mtot
		stop
	endif

	call icon (name,ia(i))
	ia(i+8) = nr
	ia(i+9) = nc
	ia(i+10) = np 
	ia(i+11) = istr
	ia(i+12) = 0
	ia(i+13) = 0
	write (ntm,100) name,type(np),nr,nc,istr
	return

  100   format (' (DEFDIR) ',8a1,a,' nr:',i6,' nc:',i6,' at:',i6)

 2000   format (' (DEFDIR) Memoria esgotada em DEFDIR para: ',8a1,
     .        /,'          Dimensao necessaria  = ',i7,
     .        /,'          Dimensao disponivel  = ',i7)

	end
c
	subroutine locate (name,na,nr,nc)

	character*1     name(8)
	integer         iname(8)
	character*3     type

      include 'common.h'
	common  /dbsys/         numa,next,idir,ip(3),type(3)
	common  /iolist/        ntm,ntr,nin,not,nsp,nfl,nt7,nt8

	na = 0
	call icon (name,iname)
	i = ifind (iname,0)
	if ( i.ne.0 ) then
		nr = ia(i+8)
		nc = ia(i+9)
		np = ia(i+10)
		na = ia(i+11)
		write (ntm,100) name,type(np),nr,nc,na,ia(i+12)
	else
		write (ntm,110) name
		write (not,110) name
		stop
	endif
	return
  100   format (' (LOCATE) ',8a1,a,' nr:',i6,' nc:',i6,
     .          ' at:',i6,' size:',i6)
  110   format (' (LOCATE) ',8a1,' NAO FOI ENCONTRADO.')
	end
c
	subroutine icon (name,iname)

c       converte caracter para inteiro

	character*1 name(8)
	dimension   iname(8)
	do i = 1,8
		iname(i) = ichar(name(i))
	enddo
	return
	end
c
	integer function ifind (iname,lun)
	dimension       iname(8)
	character*3 type
      include 'common.h'
	common  /dbsys/ numa,next,idir,ip(3),type(3)
	
c       procura diretorio do arranjo (=0 nao esta em memoria)

	i = idir
	do 100 n = 1, numa

		if (lun.ne.ia(i+13)) go to 100

		if (iname(1).ne.ia(i  )) go to 100
		if (iname(2).ne.ia(i+1)) go to 100
		if (iname(3).ne.ia(i+2)) go to 100
		if (iname(4).ne.ia(i+3)) go to 100
		if (iname(5).ne.ia(i+4)) go to 100
		if (iname(6).ne.ia(i+5)) go to 100
		if (iname(7).ne.ia(i+6)) go to 100
		if (iname(8).eq.ia(i+7)) go to 200

  100   i = i + 14
	i = 0
  200   ifind = i

	return
	end