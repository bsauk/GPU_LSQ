PROGRAM test1

  !     Test vmove.
  !     The LS fit to the data is:  Y = 1 + X + X**2 + X**3
  !     i.e. all of the regression coefficients should be exactly 1.0.
  !     An extra variable equal to X + X**2 is inserted to give a singularity.

  USE lsq
  IMPLICIT NONE

  REAL ( dp ) :: x(11,4), y(11), xrow(5), yelem, beta(5), one = 1.D0
  INTEGER           :: nvar, i, j, ifault
  LOGICAL           :: fit_const, lindep(5)
  CHARACTER (LEN=1) :: key

  DATA y/65647., 70638., 75889., 81399., 87169., 93202., 99503., &
       106079., 112939., 120094., 127557./

  WRITE(*, *)'Fitting nasty cubic'
  WRITE(*, *)'1st 4 regression coefficients should equal 1.0'
  WRITE(*, *)'Last variable = X + X^2, to introduce a deliberate singularity'
  WRITE(*, *)

  DO i = 1, 11
     x(i,1) = i + 39
     x(i,2) = x(i,1)**2
     x(i,3) = x(i,1) * x(i,2)
     x(i,4) = x(i,1) + x(i,2)
  END DO

  !     Use includ to form the orthogonal reduction.

  nvar = 4
  fit_const = .true.
  CALL startup(nvar, fit_const)

  DO i = 1, 11
     xrow(1) = one
     DO j = 1, 4
        xrow(j+1) = x(i,j)
     END DO
     yelem = y(i)
     CALL includ(one, xrow, yelem)
  END DO

  !     CALL tolset to set tolerances, then sing to detect singularities.

  CALL tolset
  CALL sing(lindep, ifault)
  WRITE(*, *)'From routine SING, IFAULT =', ifault
  WRITE(*, *)'Array LINDEP:', lindep
  WRITE(*, *)
  WRITE(*, *)'sserr = ', sserr, '   Should be 286.000'

  !     Calculate residual sums of squares (RSS).

  CALL ss

  !     Calculate regression coefficients, using vmove to cycle through
  !     the order of the variables.

  DO i = 1, 6
     CALL regcf(beta, 5, ifault)
     WRITE(*, 920) vorder
920  FORMAT(1x, 'Variable order:'/ 1x, i10, 4i15)
     WRITE(*, 900) beta
900  FORMAT(' Regn. coeffs.:'/ 4x, 5g15.7)
     WRITE(*, 910) d, r, rhs, rss
910  FORMAT(' d: ', 5g15.7/ ' r: ', 4g15.7/19x,3g15.7/34x,2g15.7/ &
          49x,g15.7/ ' rhs: '/4x, 5g15.7/ ' rss: '/4x, 5g15.7//)

     WRITE(*, *)'Press ENTER to continue'
     READ(*, '(a)') key

     !     Move variable in 1st position to the last.

     CALL vmove(1, 5, ifault)
  END DO
END PROGRAM test1


