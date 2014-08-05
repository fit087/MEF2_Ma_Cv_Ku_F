      SUBROUTINE A1 ( A , V , MAXA , N )

      DOUBLE PRECISION    A(*) , V(*)
      INTEGER             MAXA(*) , N , I

      DO 10 I = 1 , N
   10   A(I) = A(I) / V ( MAXA(N-I) )

      RETURN
      END

      SUBROUTINE A2 ( A , V , DA , N )

      DOUBLE PRECISION    A(*) , V(*) , DA
      INTEGER             N , I

      DO 10 I = 1 , N
   10   V( N-I ) = V( N-I ) + DA*A(I)

      RETURN
      END

      DOUBLE PRECISION FUNCTION BD ( A , V , MAXA , N )
      DOUBLE PRECISION A(*) , V(*)
      INTEGER          MAXA(*) , N , I

      BD = 0.D0

      DO 10 I = 1 , N
   10   BD = BD + A(I)*A(I) / V ( MAXA(N-I) )

      RETURN
      END

      DOUBLE PRECISION FUNCTION CD ( A , V , N )
      DOUBLE PRECISION A(*) , V(*)
      INTEGER          N , I

      CD = 0.D0

      DO 10 I = 1 , N
   10   CD = CD + A(I)*V(N-I)

      RETURN
      END
      
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
	dtemp = dtemp + dx(ix)*dy(iy)
	ix = ix + incx
	iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
	dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
	dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end