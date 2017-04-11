PROGRAM test7

!      Test the calculation of the diagonal elements of the hat matrix
!      using the Cloud seeding example on page 4 of Cook & Weisberg.
!      The data are contained in the file: woodley.dat
!      The elements of the diagonal of the hat matrix are on page 127.
!      The model is basically that given in (1.1.3) on page 5, but with
!      the cube root transformation of variable P.

USE lsq
IMPLICIT NONE

INTEGER           :: ier, row, in = 8, col, ncases = 24, np = 11
REAL ( lsq_kind ) :: x(24,6), y(24), xrow(11), hii(24), wt = 1.0, total,   &
                     third, one = 1.0, zero = 0.0
LOGICAL           :: seeded

third = one / 3

WRITE(*, *)'Using the cloud-seeding data from page 4 of Cook & Weisberg'
WRITE(*, *)

!      Read in the data.

OPEN(in, file='woodley.dat', status='old')
DO row = 1, ncases
  READ(in, *)(x(row, col), col=1,6), y(row)
END DO

CALL startup(np, .false.)
xrow(1) = one
DO row = 1, ncases
  xrow(3) = x(row, 2)
  xrow(4) = x(row, 3)
  xrow(5) = x(row, 4)
  xrow(6) = x(row, 5)**third
  xrow(7) = x(row, 6)
  seeded = (x(row, 1) .eq. one)
  IF (seeded) THEN
    xrow(2) = one
    xrow(8) = x(row, 3)
    xrow(9) = x(row, 4)
    xrow(10) = xrow(6)
    xrow(11) = x(row, 6)
  ELSE
    xrow(2) = zero
    xrow(8) = zero
    xrow(9) = zero
    xrow(10) = zero
    xrow(11) = zero
  END IF
  CALL includ(wt, xrow, y(row))
END DO ! row = 1, ncases

CALL tolset
total = zero
DO row = 1, ncases
  xrow(3) = x(row, 2)
  xrow(4) = x(row, 3)
  xrow(5) = x(row, 4)
  xrow(6) = x(row, 5)**third
  xrow(7) = x(row, 6)
  seeded = (x(row, 1) .eq. one)
  IF (seeded) THEN
    xrow(2) = one
    xrow(8) = x(row, 3)
    xrow(9) = x(row, 4)
    xrow(10) = xrow(6)
    xrow(11) = x(row, 6)
  ELSE
    xrow(2) = zero
    xrow(8) = zero
    xrow(9) = zero
    xrow(10) = zero
    xrow(11) = zero
  END IF
  CALL hdiag(xrow, np, hii(row), ier)
  total = total + hii(row)
END DO ! row = 1, ncases

WRITE(*, *)'Diagonal elements of Hat matrix:'
WRITE(*, '(1x, 6f10.4)') hii
WRITE(*, *)
WRITE(*, *)'Total of diagonal elements = ', total, '   Should be exactly 11'
WRITE(*, *)
WRITE(*, *)'The diagonal elements should agree with those given on page 127,'
WRITE(*, *)'of Cook & Weisbergs book'

END PROGRAM test7
