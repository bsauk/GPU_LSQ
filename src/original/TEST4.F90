PROGRAM test4

!     Test of subroutine INV.

USE lsq
IMPLICIT NONE

REAL ( lsq_kind ) :: rinv(15), rho = 0.6D0, zero = 0.0
INTEGER           :: pos, row, col

CALL startup(7, .false.)

pos = 1
DO row = 1, 6
  DO col = row+1, 7
    IF (col .eq. row+1) THEN
      r(pos) = rho
    ELSE
      r(pos) = zero
    END IF
    pos = pos + 1
  END DO ! col = row+1, 7
END DO ! row = 1, 6

CALL inv(6, rinv)
WRITE(*, 900) rinv
900 FORMAT(1x, 5g14.6/ 15x, 4g14.6/ 29x, 3g14.6/ 43x, 2g14.6/ 57x, g14.6)

END PROGRAM test4
