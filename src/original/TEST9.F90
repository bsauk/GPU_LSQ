PROGRAM test9

!     The data set AKE_B253.dat is from p.253 of the book:
!     Bjorck, A. (1996) `Numerical methods for least squares problems',
!     Publisher: SIAM   ISBN 0-89871-360-9
!     It is very ill-conditioned but NOT singular.
!     The regression coefficients should all be 1.0 exactly.

USE lsq
IMPLICIT NONE

REAL ( lsq_kind ) :: x(20), y, beta(20), one = 1.D0, error(20), proj(20)
INTEGER           :: nvar, i, j, ifault
LOGICAL           :: fit_const, lindep(20)
CHARACTER         :: text*20, key

WRITE(*, *)'Using data file ake_b253.dat for a problem in book of Ake Bjorck'
WRITE(*, *)

OPEN(9, file='ake_b253.dat', status='old')
READ(9, *) text

nvar = 20
fit_const = .FALSE.
CALL startup(nvar, fit_const)

!     Use includ to form the orthogonal reduction.

DO i = 1, 20
  READ(9, *) x, y
  CALL includ(one, x, y)
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

DO i = 1, 20
  CALL regcf(beta, 20, ifault)
  error = beta - one
  proj = SQRT(d) * rhs
  WRITE(*, 900) (beta(j), error(j), d(j), proj(j), j=1,20)
  900 FORMAT(' Regn.coeffs.   Error     Row mult.  Y-projection',   &
             / (1x, f10.7, 2x, 3g12.4))
  WRITE(*, *)'Press ENTER to continue'
  READ(*, '(a)') key

!     Move variable in 1st position to the last.

  CALL vmove(1, 20, ifault)
  CALL sing(lindep, ifault)
  IF (ifault /= 0) WRITE(*, *)'Singularity detected'
END DO

END PROGRAM test9

