PROGRAM test3

!     Test vmove for treatment of singularities.

USE lsq
IMPLICIT NONE

REAL (KIND(0.0))       :: rand
REAL ( lsq_kind )      :: x(7), y, beta(7), one = 1.0
INTEGER                :: i, ier
INTEGER (KIND(123456)) :: ix, iy, iz
LOGICAL                :: lindep(7)
CHARACTER (LEN=1)      :: key
common /randc/ ix, iy, iz

ix = 777
iy = 777
iz = 777

!     Initialize orthogonal reduction.

CALL startup(7, .false.)     ! No constant to be fitted

!     Generate data such that:
!             x2 = x4 - x5
!             x6 = x1 - x3
!             y  = x1 + x3 + x4 + x5 + x7

do i = 1, 12
  x(1) = rand()
  x(3) = rand()
  x(4) = rand()
  x(5) = rand()
  x(7) = rand()
  x(2) = x(4) - x(5)
  x(6) = x(1) - x(3)
  y = x(1) + x(3) + x(4) + x(5) + x(7)
  CALL includ(one, x, y)
END DO
    WRITE(*, *)'As output from includ:'
    WRITE(*, 900) d, r, rhs
900 FORMAT(' d:'/1x,7f11.6/                                         &
           ' r:'/1x,6f11.6/ 12x,5f11.6/ 23x,4f11.6/ 34x,3f11.6/     &
                 45x,2f11.6/ 56x,f11.6/                             &
           ' rhs:'/1x,7f11.6//)

!   Set up arrays TOL & RSS.

    CALL tolset
    CALL ss
    WRITE(*, 920) rss
920 FORMAT(' RSS:'/ 1x, 7f11.5/)

WRITE(*, *)'Press ENTER to continue'
READ(*, '(a)') key

!   Use SING to set near zeroes to zero.

CALL sing(lindep, ier)
WRITE(*, *)'After being processed by routine SING:'
IF (ier == 0) THEN
  WRITE(*, *)'QR-factorization is not singular'
ELSE
  DO i = 1, 7
    IF (lindep(i)) THEN
      WRITE(*, *) 'Variable', i, ' is exactly linearly related to earlier variables'
    END IF
  END DO ! i = 1, 7
END IF ! (ier == 0)

WRITE(*, *)
WRITE(*, 900) d, r, rhs

!   Swap rows 4 & 5 and rows 6 & 7.

    CALL vmove(4, 5, ier)
    CALL vmove(6, 7, ier)
    WRITE(*, *)'After interchange of variables 4 & 5, and 6 & 7'
    WRITE(*, 900) d, r, rhs
    CALL regcf(beta, 7, ier)
    WRITE(*, 910) beta
910 FORMAT(' Regression coefficients:'/ 1x,7f11.6)
END PROGRAM test3

include 'whran.f90'

