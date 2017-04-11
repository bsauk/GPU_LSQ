PROGRAM test5

!     Test of subroutine COV.

USE lsq
IMPLICIT NONE

REAL ( lsq_kind ) :: rho = 0.6D0, var, covmat(21), sterr(6),     &
                     one = 1.0, zero = 0.0
INTEGER           :: pos, row, col, ier

CALL startup(7, .false.)     ! No constant

pos = 1
DO row = 1, 6
  d(row) = one
  rhs(row) = one
  DO col = row+1, 7
    IF (col .eq. row+1) THEN
      r(pos) = rho
    ELSE
      r(pos) = zero
    END IF
    pos = pos + 1
  END DO ! col = row+1, 7
END DO ! row = 1, 6

sserr = one                  ! This would normally be calculated by INCLUD
nobs = 20                    ! Ditto
CALL ss

CALL cov(6, var, covmat, 21, sterr, ier)
IF (ier == 0) THEN
  WRITE(*, 900) covmat
  900 FORMAT(1x, 6f12.6/ 13x, 5f12.6/ 25x, 4f12.6/ 37x, 3f12.6/  &
            49x, 2f12.6/ 61x, f12.6)
ELSE
  WRITE(*, *)'Error', ier, ' returned from routine COV'
END IF

END PROGRAM test5
