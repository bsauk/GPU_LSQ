PROGRAM test8

!     Test handling of a case with more variables than cases.
!     X1, .., X4 from random number generator; Y = X1 + X2 + X3 + X4
!     X5 = X1 + X2
!     X6 =      X2 + X3
!     X7 =           X3 + X4
!     X8 = X1           + X4   Hence Y = (X5 + X6 + X7 + X8)/2

!     ncases = 5

!     The matrix R should look like:
!                (1)   1   2   3   4   5   6   7
!                     (1)  8   9  10  11  12  13
!                         (1) 14  15  16  17  18
!                             (1) 19  20  21  22
!                                 (1) 23  24  25
!                                     (1) 26  27
!                                         (1) 28
!                                             (1)
!     where (1) indicates an implicit 1.0 and the other numbers are the
!     locations within array R.   Elements 26-28 should be zero.
!
!--------------------------------------------------------------------------

USE lsq
IMPLICIT NONE

INTEGER, PARAMETER     :: np = 8
INTEGER                :: ier, case, i, list(np), new
INTEGER (KIND(123456)) :: ix, iy, iz
REAL (KIND(0.E0))      :: rand
REAL ( lsq_kind )      :: wt = 1.0, x(8), y, beta(np)
LOGICAL                :: lindep(np)

common /randc/ ix, iy, iz

ix = 777
iy = 777
iz = 777

CALL startup(np, .false.)

WRITE(*, *)'Example with 8 variables but only 5 observations'
WRITE(*, *)'Variables X1 to X4 are random, Y = X1 + X2 + X3 + X4'
WRITE(*, *)'   X5 = X1 + X2          '
WRITE(*, *)'   X6 =      X2 + X3     '
WRITE(*, *)'   X7 =           X3 + X4'
WRITE(*, *)'   X8 = X1           + X4'
WRITE(*, *)

!     Generate 5 lines of data

DO case = 1, 5
  DO i = 1, 4
    x(i) = rand()
  END DO ! i = 1, 4
  x(5) = x(1) + x(2)
  x(6) = x(2) + x(3)
  x(7) = x(3) + x(4)
  x(8) = x(4) + x(1)
  y = x(5) + x(7)
  CALL includ(wt, x, y)
END DO ! case = 1, 5

!     Now look at the factorization

    WRITE(*, 900) d
900 format(' D = ', 8f9.5)
    WRITE(*, 910) r(1:7)
910 format(' r = '/ 1x, 'Row 1: (1) ', 7f9.5)
    WRITE(*, 920) r(8:13)
920 format(1x, 'Row 2:', 9x, ' (1) ', 6f9.5)
    WRITE(*, 930) r(14:18)
930 format(1x, 'Row 3:', 18x, ' (1) ', 5f9.5)
    WRITE(*, 940) r(19:22)
940 format(1x, 'Row 4:', 27x, ' (1) ', 4f9.5)
    WRITE(*, 950) r(23:25)
950 format(1x, 'Row 5:', 36x, ' (1) ', 3f9.5)
    WRITE(*, 960) r(26:28)
960 format(1x, 'Row 6:', 45x, ' (1) ', 2f9.5/ 1x, 'Row 7:', 54x, ' (1) ', f9.5)

CALL tolset
CALL sing(lindep, ier)
WRITE(*, *)' IER = ', ier
WRITE(*, *)' LINDEP = ', lindep
CALL ss

!     Calculate regression coefficients for first 4 variables.

    CALL regcf(beta, 4, ier)
    WRITE(*, 970) (beta(i),i=1,4)
970 FORMAT(' Regression coefficients:'/ 1x, 4f10.5//)

!     Re-order variables as X5 .. X8, X1 .. X4
!     N.B.  X5 + X7 = X6 + X8

WRITE(*, *)'New order of variables is X5 .. X8, X1 .. X4'
WRITE(*, *)

new = 5
DO i = 1, 8
  list(i) = new
  IF (i .eq. 4) THEN
    new = 1
  else
    new = new + 1
  end if
END DO ! i = 1, 8
CALL reordr(list, 4, 1, ier)

!     Now look at the factorization again

WRITE(*, 900) d
WRITE(*, 910) r(1:7)
WRITE(*, 920) r(8:13)
WRITE(*, 930) r(14:18)
WRITE(*, 940) r(19:22)
WRITE(*, 950) r(23:25)
WRITE(*, 960) r(26:28)
CALL sing(lindep, ier)
WRITE(*, *)' IER = ', ier
WRITE(*, *)' LINDEP = ', lindep

!     Calculate regression coefficients for first 4 variables.

CALL regcf(beta, 4, ier)
WRITE(*, 970) beta(1:4)

END PROGRAM test8

include 'whran.f90'
