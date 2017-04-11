PROGRAM test2

!     Test vmove using the Longley data.

USE lsq
IMPLICIT NONE

REAL ( lsq_kind ) :: x(16,6), y(16), xrow(7), yelem, beta(7), one = 1.D0
INTEGER           :: nvar, i, j, ifault
LOGICAL           :: fit_const, lindep(7)
CHARACTER         :: text*20, key

WRITE(*, *)'Using the Longley data'
WRITE(*, *)

OPEN(9, file='longley.dat', status='old')
READ(9, '(a20)') text
READ(9, '(a20)') text
DO i = 1, 16
  READ(9, *) (x(i,j),j=1,6), y(i)
END DO

nvar = 6
fit_const = .true.
CALL startup(nvar, fit_const)

!     Use includ to form the orthogonal reduction.

DO i = 1, 16
  xrow(1) = one
  DO j = 1, 6
    xrow(j+1) = x(i,j)
  END DO
  yelem = y(i)
  CALL includ(one, xrow, yelem)
END DO

!     CALL sing to detect singularities.

CALL tolset
CALL sing(lindep, ifault)
IF (ifault == 0) THEN
  WRITE(*, *)'QR-factorization is not singular'
ELSE
  DO i = 1, nvar
    IF (lindep(i)) THEN
      WRITE(*, *) 'Variable', i, ' is exactly linearly related to earlier variables'
    END IF
  END DO ! i = 1, nvar
END IF ! (ifault == 0)

WRITE(*, *)'sserr = ', sserr
CALL ss
WRITE(*, *)

!     Calculate regression coefficients, using vmove to cycle through
!     the order of the variables.

    DO i = 1, 7
      CALL regcf(beta, 7, ifault)
      WRITE(*, 900) beta
900   FORMAT(' Regn. coeffs.:'/ 1x, 5g15.7/ 1x, 2g15.7)
      WRITE(*, 910) d, r, rhs
910   FORMAT(' d: ', 5g13.5/1x,2g13.5 &
           / ' r:'/ 1x,6g13.5/ 14x,5g13.5/ 27x,4g13.5/    &
          40x,3g13.5/ 53x,2g13.5/ 66x,g13.5/              &
           ' rhs: '/1x, 5g13.5/ 1x, 2g13.5//)

WRITE(*, *)'Press ENTER to continue'
READ(*, '(a)') key

!     Move variable in 1st position to the last.

      CALL vmove(2, 7, ifault)
    END DO

END PROGRAM test2

