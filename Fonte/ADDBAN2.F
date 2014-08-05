      subroutine addban2    (a,maxa,s,lm,ndt)
c  
c
c **********************************************************************
c  r o t i n a    a d d b a n 2
c
c	para montagem da matriz de rigidez dos elementos na
c	matriz de rigidez global compactada (skyline), triangulo superior
c
c	a matriz de elemento a ser montada deve estar armazenada por 
c	colunas, apenas os termos acima da diagonal inclusive:
c
c	    s1      s2      s4	     s7	     s11	s16	    s22     s29
c		        s3      s5	     s8	     s12	s17	    s23     s30
c			            s6	     s9	     s13	s18	    s24     s31
c				                s10     s14	    s19	    s25     s32
c					                    s15     s20	    s26     s33
c						                        s21	    s27     s34
c							                            s28     s35
c                                                                 s36
c
c
c **********************************************************************
c
      implicit real*8 (a-h,o-z)        
      include 'tapes.h'
      dimension a(1),maxa(1),s(1),lm(1)

	do j = 1, ndt

		   jeq = lm(j)
		   if (jeq.gt.0) then

		   do i = 1, j

			      ieq = lm(i)
			      if (ieq.gt.0) then

			          is = j*(j-1)/2 + i

			          if (jeq.ge.ieq) then
			              kk = maxa(jeq) + jeq-ieq 
c                      write (iout,100) i,j  , ieq,jeq 
			              
			          else
			              kk = maxa(ieq) + ieq-jeq 
c                      write (iout,100) i,j  , jeq,ieq 

			          endif
			          
			          a(kk) = a(kk) + s(is) 

              
c  100         format ('local(',2i2,') -> global(',2i2,') pos ',2i5)                       

			      endif
			      
		   enddo

		   endif
	
	enddo

     	return
      end




