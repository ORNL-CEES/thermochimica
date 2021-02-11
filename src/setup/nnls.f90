MODULE precision
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15, 60)
END MODULE precision




!     SUBROUTINE nnls(a, m, n, b, x, rnorm, w, indx, mode)
!
!  Algorithm NNLS: NONNEGATIVE LEAST SQUARES
!
!  The original version of this code was developed by
!  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
!  1973 JUN 15, and published in the book
!  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
!  Revised FEB 1995 to accompany reprinting of the book by SIAM.
!
!  This translation into Fortran 90 by Alan Miller, February 1997
!  Latest revision - 15 April 1997

!  N.B. The following call arguments have been removed:
!       mda, zz
!
!  GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
!  N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM
!
!                   A * X = B  SUBJECT TO X  >=  0
!  ------------------------------------------------------------------
!                  Subroutine Arguments
!
!  A(), M, N   ON ENTRY, A() CONTAINS THE M BY N MATRIX, A.
!              ON EXIT, A() CONTAINS THE PRODUCT MATRIX, Q*A , WHERE Q IS AN
!              M x M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY THIS SUBROUTINE.
!  B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CONTAINS Q*B.
!  X()     ON ENTRY X() NEED NOT BE INITIALIZED.
!          ON EXIT X() WILL CONTAIN THE SOLUTION VECTOR.
!  RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE RESIDUAL VECTOR.
!  W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN THE DUAL
!          SOLUTION VECTOR.   W WILL SATISFY W(I) = 0. FOR ALL I IN SET P
!          AND W(I) <= 0. FOR ALL I IN SET Z
!  INDX()  AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
!          ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS P AND Z
!          AS FOLLOWS..
!              INDX(1)   THRU INDX(NSETP) = SET P.
!              INDX(IZ1) THRU INDX(IZ2)   = SET Z.
!              IZ1 = NSETP + 1 = NPP1
!              IZ2 = N
!  MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS.
!          1   THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
!          2   THE DIMENSIONS OF THE PROBLEM ARE BAD.
!              EITHER M <= 0 OR N <= 0.
!          3   ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS.
!
!  ------------------------------------------------------------------
SUBROUTINE nnls (a, m, n, b, x, rnorm, w, indx, mode)
!     ------------------------------------------------------------------
USE precision
IMPLICIT NONE
INTEGER, INTENT(IN)       :: m, n
INTEGER, INTENT(OUT)      :: mode
INTEGER, INTENT(OUT), dimension(n)     :: indx
REAL(8), INTENT(IN OUT), dimension(m,n) :: a
REAL(8), INTENT(IN OUT), dimension(m) :: b
REAL(8), INTENT(OUT), dimension(n)    :: x
REAL(8), INTENT(OUT)    :: rnorm
REAL(8), INTENT(OUT), dimension(n)    :: w

INTERFACE
  SUBROUTINE g1(a, b, cterm, sterm, sig)
    USE precision
    IMPLICIT NONE
    REAL (dp), INTENT(IN)  :: a, b
    REAL (dp), INTENT(OUT) :: cterm, sterm, sig
  END SUBROUTINE g1

  SUBROUTINE h12(mode, lpivot, l1, m, u, up, c, ice, icv, ncv)
    USE precision
    IMPLICIT NONE
    INTEGER, INTENT(IN)                     :: mode, lpivot, l1, m, ice, icv, &
                                               ncv
    REAL (dp), DIMENSION(:), INTENT(IN OUT) :: u, c
    REAL (dp), INTENT(IN OUT)               :: up
  END SUBROUTINE h12
END INTERFACE

! Local variables

INTEGER                 :: i, ii, ip, iter, itmax, iz, iz1, iz2, izmax,   &
                           j, jj, jz, l, mda, npp1, nsetp
REAL (dp), DIMENSION(m) :: zz
REAL (dp), DIMENSION(1) :: dummy
REAL (dp)               :: alpha, asave, cc, factor = 0.01_dp, sm, &
                           ss, t, temp, two = 2.0_dp, unorm, up, wmax,    &
                           zero = 0.0_dp, ztest
!     ------------------------------------------------------------------
mode = 1
IF (m <= 0 .OR. n <= 0) THEN
  mode = 2
  RETURN
END IF
iter = 0
itmax = 3*n

!                    INITIALIZE THE ARRAYS indx() AND X().

DO i = 1,n
  x(i) = zero
  indx(i) = i
END DO

iz2 = n
iz1 = 1
nsetp = 0
npp1 = 1
!                             ******  MAIN LOOP BEGINS HERE  ******
!                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
!                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.

30 IF (iz1 > iz2 .OR. nsetp >= m) GO TO 350

!         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().

DO iz = iz1,iz2
  j = indx(iz)
  w(j) = DOT_PRODUCT(a(npp1:m,j), b(npp1:m))
END DO

!                                   FIND LARGEST POSITIVE W(J).
60 wmax = zero
DO iz = iz1,iz2
  j = indx(iz)
  IF (w(j) > wmax) THEN
    wmax = w(j)
    izmax = iz
  END IF
END DO

!             IF WMAX  <=  0. GO TO TERMINATION.
!             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.

IF (wmax <= zero) GO TO 350
iz = izmax
j = indx(iz)

!     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.
!     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID
!     NEAR LINEAR DEPENDENCE.

asave = a(npp1,j)
CALL h12 (1, npp1, npp1+1, m, a(:,j), up, dummy, 1, 1, 0)
unorm = zero
IF (nsetp  /=  0) THEN
  unorm = SUM( a(1:nsetp,j)**2 )
END IF
unorm = SQRT(unorm)
IF (unorm + ABS(a(npp1,j))*factor - unorm  >  zero) THEN

!        COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
!        AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).

  zz(1:m) = b(1:m)
  CALL h12 (2, npp1, npp1+1, m, a(:,j), up, zz, 1, 1, 1)
  ztest = zz(npp1)/a(npp1,j)

!                                     SEE IF ZTEST IS POSITIVE

  IF (ztest > zero) GO TO 140
END IF

!     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.
!     RESTORE A(NPP1,J), SET W(J) = 0., AND LOOP BACK TO TEST DUAL
!     COEFFS AGAIN.

a(npp1,j) = asave
w(j) = zero
GO TO 60

!     THE INDEX  J = indx(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
!     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER
!     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN
!     COL J,  SET W(J) = 0.

140 b(1:m) = zz(1:m)

indx(iz) = indx(iz1)
indx(iz1) = j
iz1 = iz1+1
nsetp = npp1
npp1 = npp1+1

mda = SIZE(a,1)
IF (iz1  <=  iz2) THEN
  DO jz = iz1,iz2
    jj = indx(jz)
    CALL h12 (2, nsetp, npp1, m, a(:,j), up, a(:,jj), 1, mda, 1)
  END DO
END IF

IF (nsetp /= m) THEN
  a(npp1:m,j) = zero
END IF

w(j) = zero
!                                SOLVE THE TRIANGULAR SYSTEM.
!                                STORE THE SOLUTION TEMPORARILY IN ZZ().
CALL solve_triangular(zz)

!                       ******  SECONDARY LOOP BEGINS HERE ******

!                          ITERATION COUNTER.

210 iter = iter+1
IF (iter > itmax) THEN
  mode = 3
  WRITE (*,'(/a)') ' NNLS quitting on iteration count.'
  GO TO 350
END IF

!                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.
!                                  IF NOT COMPUTE ALPHA.

alpha = two
DO ip = 1,nsetp
  l = indx(ip)
  IF (zz(ip)  <=  zero) THEN
    t = -x(l)/(zz(ip)-x(l))
    IF (alpha > t) THEN
      alpha = t
      jj = ip
    END IF
  END IF
END DO

!          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL
!          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.

IF (alpha == two) GO TO 330

!          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO
!          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.

DO ip = 1,nsetp
  l = indx(ip)
  x(l) = x(l) + alpha*(zz(ip)-x(l))
END DO

!        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I
!        FROM SET P TO SET Z.

i = indx(jj)
260 x(i) = zero

IF (jj /= nsetp) THEN
  jj = jj+1
  DO j = jj,nsetp
    ii = indx(j)
    indx(j-1) = ii
    CALL g1 (a(j-1,ii), a(j,ii), cc, ss, a(j-1,ii))
    a(j,ii) = zero
    DO l = 1,n
      IF (l /= ii) THEN

!                 Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))

        temp = a(j-1,l)
        a(j-1,l) = cc*temp + ss*a(j,l)
        a(j,l)   = -ss*temp + cc*a(j,l)
      END IF
    END DO

!                 Apply procedure G2 (CC,SS,B(J-1),B(J))

    temp = b(j-1)
    b(j-1) = cc*temp + ss*b(j)
    b(j)   = -ss*temp + cc*b(j)
  END DO
END IF

npp1 = nsetp
nsetp = nsetp-1
iz1 = iz1-1
indx(iz1) = i

!        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
!        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
!        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY
!        THAT ARE NONPOSITIVE WILL BE SET TO ZERO
!        AND MOVED FROM SET P TO SET Z.

DO jj = 1,nsetp
  i = indx(jj)
  IF (x(i) <= zero) GO TO 260
END DO

!         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.

zz(1:m) = b(1:m)
CALL solve_triangular(zz)
GO TO 210
!                      ******  END OF SECONDARY LOOP  ******

330 DO ip = 1,nsetp
  i = indx(ip)
  x(i) = zz(ip)
END DO
!        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.
GO TO 30

!                        ******  END OF MAIN LOOP  ******

!                        COME TO HERE FOR TERMINATION.
!                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.

350 sm = zero
IF (npp1 <= m) THEN
  sm = SUM( b(npp1:m)**2 )
ELSE
  w(1:n) = zero
END IF
rnorm = SQRT(sm)
RETURN

CONTAINS

SUBROUTINE solve_triangular(zz)

!     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE
!     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ().

REAL (dp), INTENT(IN OUT) :: zz(:)

DO l = 1,nsetp
  ip = nsetp+1-l
  IF (l  /=  1) zz(1:ip) = zz(1:ip) - a(1:ip,jj)*zz(ip+1)
  jj = indx(ip)
  zz(ip) = zz(ip) / a(ip,jj)
END DO

RETURN
END SUBROUTINE solve_triangular

END SUBROUTINE nnls



SUBROUTINE g1(a, b, cterm, sterm, sig)

!     COMPUTE ORTHOGONAL ROTATION MATRIX..

!  The original version of this code was developed by
!  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
!  1973 JUN 12, and published in the book
!  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
!  Revised FEB 1995 to accompany reprinting of the book by SIAM.

!     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))
!                        (-S,C)         (-S,C)(B)   (   0          )
!     COMPUTE SIG = SQRT(A**2+B**2)
!        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT
!        SIG MAY BE IN THE SAME LOCATION AS A OR B .
!     ------------------------------------------------------------------
USE precision
IMPLICIT NONE
REAL (dp), INTENT(IN)  :: a, b
REAL (dp), INTENT(OUT) :: cterm, sterm, sig

!     Local variables
REAL (dp) :: one = 1.0D0, xr, yr, zero = 0.0D0
!     ------------------------------------------------------------------
IF (ABS(a) > ABS(b)) THEN
  xr = b / a
  yr = SQRT(one + xr**2)
  cterm = SIGN(one/yr, a)
  sterm = cterm * xr
  sig = ABS(a) * yr
  RETURN
END IF

IF (b /= zero) THEN
  xr = a / b
  yr = SQRT(one + xr**2)
  sterm = SIGN(one/yr, b)
  cterm = sterm * xr
  sig = ABS(b) * yr
  RETURN
END IF

!      SIG = ZERO
cterm = zero
sterm = one
RETURN
END SUBROUTINE g1



!     SUBROUTINE h12 (mode, lpivot, l1, m, u, up, c, ice, icv, ncv)

!  CONSTRUCTION AND/OR APPLICATION OF A SINGLE
!  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B

!  The original version of this code was developed by
!  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
!  1973 JUN 12, and published in the book
!  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
!  Revised FEB 1995 to accompany reprinting of the book by SIAM.
!     ------------------------------------------------------------------
!                     Subroutine Arguments

!     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a
!            Householder transformation, or Algorithm H2 to apply a
!            previously constructed transformation.
!     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT.
!     L1,M   IF L1  <=  M   THE TRANSFORMATION WILL BE CONSTRUCTED TO
!            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M
!            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
!     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot
!            vector.  IUE is the storage increment between elements.
!            On exit when MODE = 1, U() and UP contain quantities
!            defining the vector U of the Householder transformation.
!            on entry with MODE = 2, U() and UP should contain
!            quantities previously computed with MODE = 1.  These will
!            not be modified during the entry with MODE = 2.
!     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH
!            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE
!            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED.
!            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS.
!     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().
!     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().
!     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV  <=  0
!            NO OPERATIONS WILL BE DONE ON C().
!     ------------------------------------------------------------------
SUBROUTINE h12(mode, lpivot, l1, m, u, up, c, ice, icv, ncv)
!     ------------------------------------------------------------------

USE precision
IMPLICIT NONE
INTEGER, INTENT(IN)                     :: mode, lpivot, l1, m, ice, icv, ncv
REAL (dp), DIMENSION(:), INTENT(IN OUT) :: u, c
REAL (dp), INTENT(IN OUT)               :: up

!  Local variables
INTEGER          :: i, i2, i3, i4, incr, j
REAL (dp)        :: b, cl, clinv, one = 1.0D0, sm
!     ------------------------------------------------------------------
IF (0 >= lpivot .OR. lpivot >= l1 .OR. l1 > m) RETURN
cl = ABS(u(lpivot))
IF (mode /= 2) THEN
!                            ****** CONSTRUCT THE TRANSFORMATION. ******
  DO j = l1, m
    cl = MAX(ABS(u(j)),cl)
  END DO
  IF (cl <= 0) RETURN
  clinv = one / cl
  sm = (u(lpivot)*clinv) ** 2 + SUM( (u(l1:m)*clinv)**2 )
  cl = cl * SQRT(sm)
  IF (u(lpivot) > 0) THEN
    cl = -cl
  END IF
  up = u(lpivot) - cl
  u(lpivot) = cl
ELSE
!            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******

  IF (cl <= 0) RETURN
END IF
IF (ncv <= 0) RETURN

b = up * u(lpivot)
!                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.

IF (b < 0) THEN
  b = one / b
  i2 = 1 - icv + ice * (lpivot-1)
  incr = ice * (l1-lpivot)
  DO j = 1, ncv
    i2 = i2 + icv
    i3 = i2 + incr
    i4 = i3
    sm = c(i2) * up
    DO i = l1, m
      sm = sm + c(i3) * u(i)
      i3 = i3 + ice
    END DO
    IF (sm /= 0) THEN
      sm = sm * b
      c(i2) = c(i2) + sm * up
      DO i = l1, m
        c(i4) = c(i4) + sm * u(i)
        i4 = i4 + ice
      END DO
    END IF
  END DO ! j = 1, ncv
END IF

RETURN
END SUBROUTINE h12


!================ A sample test program =====================

! PROGRAM test_nnls
! Fit a mixture of exponentials by NNLS.
! The fitting of a sum of exponentials:

!   Y = Sum A(i).exp(-b(i).t)

! to a set of values of Y is very ill-conditioned problem in non-linear
! least squares.   The parameters to be fitted are the amplitudes A(i)
! and the b(i)'s.   Often the number of terms is unknown.   Sometimes the
! discrete A(i)'s are replaced by a continuous function and the sum replaced
! with an integral.

! A common practice is to fit the sum with the b(i)'s held fixed and
! then constrain the amplitudes A(i) to be positive.   That is done in
! this example.   We will fit:

!  Y = Sum A(i).exp(-t/t0(i))

! where the t0(i)'s increase by a factor of sqrt(2) from 2.0 to 1024.
! The `data' will be generated with 3 exponential terms:

!  Y = 100 * ( exp(-t/5) + exp(-t/50) + exp(-t/500) ) + noise

! The times (t) will be 10, 20, ..., 1000.

! The nearest fitted t0'2 to each of these t0's are:
!   5 between   4   and   5.66
!  50 between  32   and  45.3
! 500 between 362.0 and 512

! Test example added 29 June 1998

! USE precision
! IMPLICIT NONE
!
! REAL (dp) :: x(100,19), y(100), t0(19), arg, arglimit, noise(100), t, rnorm, &
!              b(19), w(19)
! INTEGER   :: i, indx(19), j, mode
!
! INTERFACE
!   SUBROUTINE nnls (a, m, n, b, x, rnorm, w, indx, mode)
!     USE precision
!     IMPLICIT NONE
!     INTEGER, INTENT(IN)       :: m, n
!     INTEGER, INTENT(OUT)      :: indx(:), mode
!     REAL (dp), INTENT(IN OUT) :: a(:,:), b(:)
!     REAL (dp), INTENT(OUT)    :: x(:), rnorm, w(:)
!   END SUBROUTINE nnls
! END INTERFACE
!
! t0(1) = 2.0
! t0(2) = t0(1) * SQRT(2.0_dp)
! DO i = 3, 19
!   t0(i) = 2.0_dp * t0(i-2)
! END DO
!
! ! Calculate the X-matrix, avoiding underflow
!
! arglimit = LOG( 2.0_dp * TINY(1.0_dp) )
! DO j = 1, 19
!   DO i = 1, 100
!     arg = -10._dp * i / t0(j)
!     IF (arg > arglimit) THEN
!       x(i,j) = EXP(arg)
!     ELSE
!       x(i:100,j) = 0.0_dp
!       EXIT
!     END IF
!   END DO
! END DO
!
! ! Generate Y's with uniformly-distributed noise.
!
! CALL RANDOM_NUMBER( noise )
! DO i = 1, 100
!   t = 10._dp * i
!   y(i) = 100._dp * ( exp(-t/5._dp) + exp(-t/50._dp) + exp(-t/500._dp) ) +  &
!          noise(i) - 0.5_dp
! END DO
!
! ! Now call NNLS to do the fitting.
!
! CALL nnls (x, 100, 19, y, b, rnorm, w, indx, mode)
!
! SELECT CASE (mode)
!   CASE (1)
!     DO i = 1, 19
!       IF (b(i) > 0.0_dp) THEN
!         WRITE(*, '(a, f9.1, a, f9.1)')  &
!                  ' Time constant: ', t0(i), '  Fitted amplitude = ', b(i)
!       END IF
!     END DO
!     WRITE(*, '(a, 19i3)') ' Array INDX =', indx
!     WRITE(*, '(a/ (" ", 10f8.0))') ' Alternate solution:', w
!     WRITE(*, '(a, f9.2)') ' rnorm = ', rnorm
!   CASE (2)
!     WRITE(*, *) 'Error in input argument 2 or 3'
!   CASE (3)
!     WRITE(*, *) 'Failed to converge'
! END SELECT
!
! STOP
! END PROGRAM test_nnls
