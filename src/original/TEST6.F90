PROGRAM test6

!     Use the STEAM data from Draper & Smith to test PARTIAL_CORR.

USE lsq
IMPLICIT NONE

INTEGER           :: ier, case, col
REAL ( lsq_kind ) :: x(0:9), y, cormat(36), ycorr(0:9), one = 1.0
CHARACTER         :: text*80

OPEN(9, file='steam.dat', status='old')
CALL startup(9, .true.)

WRITE(*, *)'Using the STEAM data from Draper & Smith'
WRITE(*, *)'Variable names - dependent variable is last'
READ(9, '(a)') text
WRITE(*, '(1x, a)') text
WRITE(*, *)

DO case = 1, 25
  READ(9, *) (x(col),col=1,9), y
  x(0) = one
  CALL includ(one, x, y)
END DO ! case = 1, 25

WRITE(*, *)'ier, sserr = ', ier, sserr

CALL partial_corr(1, cormat, 36, ycorr, ier)
if (ier .ne. 0) WRITE(*, *)'ier = ', ier
WRITE(*, 900) ycorr(1:9)
900 FORMAT(' Correlations with the dependent variable:'/ 1x, 9f8.4/)
WRITE(*, 910) cormat
910 FORMAT(' Correlations amongst the predictors:'/ &
             1x, '1.0', 3x, 8f8.4/ &
             9x, '1.0', 3x, 7f8.4/ &
            17x, '1.0', 3x, 6f8.4/ &
            25x, '1.0', 3x, 5f8.4/ &
            33x, '1.0', 3x, 4f8.4/ &
            41x, '1.0', 3x, 3f8.4/ &
            49x, '1.0', 3x, 2f8.4/ &
            57x, '1.0', 3x, f8.4/ &
            65x, '1.0'/)
WRITE(*, *)'Compare with the table on page 616 of 2nd edition of D & S'


END PROGRAM test6

