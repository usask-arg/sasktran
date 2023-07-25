! ---------------------------------------------------------------------
! Fortran-90 versions of machine-constant routines R1MACH, D1MACH, I1MACH
!
! Usage:  In whatever routine needed, put the statement
!         USE RDI1MACH_f90, ONLY : R1MACH, D1MACH, I1MACH
!         (list only the x1MACH routines you actually use after ONLY)
! ---------------------------------------------------------------------


MODULE RDI1MACH_f90

! ---------------------------------------------------------------------
!  {R,D,I}1MACH revisited: no more uncommenting DATA statements
!  
!  The widely-used original routines were modified to use Fortran-90 
!  intrinsic functions.  This was not completely possible with I1MACH, 
!  which returns some parameters (logical unit numbers of standard
!  input, standard output, and standard error) that may require
!  user customization. 
!  
!  Presented at the IFIP WG 2.5 International Workshop on 
!  "Current Directions in Numerical Software and High Performance 
!  Computing", 19 - 20 October 1995, Kyoto, Japan. 
!  
!  David Gay (dmg@bell-labs.com)
!  Eric Grosse (ehg@bell-labs.coM)
!  Bell Laboratories
!  700 Mountain Avenue
!  Murray Hill, New Jersey 07974-0636
!  USA 
!  
!  References:
!  
!  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
!    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
!    August 1996.
!
!  http://www.nsc.liu.se/~boein/ifip/kyoto/workshop-info/proceedings/einarsson/d1mach.html
!
! ---------------------------------------------------------------------

PRIVATE
PUBLIC :: R1MACH, D1MACH, I1MACH

CONTAINS

REAL FUNCTION R1MACH (I)
!
!   R1MACH can be used to obtain machine-dependent parameters for
!   single precision numbers.  The results for various values of I are:
!
!   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   R1MACH(3) = B**(-T), the smallest relative spacing.
!   R1MACH(4) = B**(1-T), the largest relative spacing.
!   R1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!***REFERENCES  
!
!  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
!     ACM Transactions on Mathematical Software 4, 177-188 (1978).
!
!  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
!     Hill, NJ. submitted to ACM Transactions on Mathematical Software,
!     August 1996. 
!
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   960329  Modified for Fortran 90 (BE after suggestions by Eric Grosse)      
! --------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: I
REAL :: B, X = 1.0

B = RADIX(X)

SELECT CASE (I)
  CASE (1)
    R1MACH = B**(MINEXPONENT(X)-1) ! the smallest positive magnitude.
  CASE (2)
    R1MACH = HUGE(X)               ! the largest magnitude.
  CASE (3)
    R1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
  CASE (4)
    R1MACH = B**(1-DIGITS(X))      ! the largest relative spacing.
  CASE (5)
    R1MACH = LOG10(B)
  CASE DEFAULT
    CALL ERRMSG('R1MACH - input argument out of bounds', .TRUE.)
END SELECT

RETURN
END FUNCTION R1MACH

DOUBLE PRECISION FUNCTION D1MACH (I)
!
!   D1MACH can be used to obtain machine-dependent parameters for
!   double precision numbers.  The results for various values of I are:
!
!   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH(3) = B**(-T), the smallest relative spacing.
!   D1MACH(4) = B**(1-T), the largest relative spacing.
!   D1MACH(5) = LOG10(B)
!
!   Assume double precision numbers are represented in the T-digit,
!   base-B form
!
!        sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!***REFERENCES  
!
!  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
!     ACM Transactions on Mathematical Software 4, 177-188 (1978).
!
!  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
!    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
!    August 1996. 
!
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   960329  Modified for Fortran 90 (BE after suggestions by Eric Grosse)      
! --------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: I
DOUBLE PRECISION :: B, X = 1.D0

B = RADIX(X)

SELECT CASE (I)
  CASE (1)
    D1MACH = B**(MINEXPONENT(X)-1) ! the smallest positive magnitude.
  CASE (2)
    D1MACH = HUGE(X)               ! the largest magnitude.
  CASE (3)
    D1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
  CASE (4)
    D1MACH = B**(1-DIGITS(X))      ! the largest relative spacing.
  CASE (5)
    D1MACH = LOG10(B)
  CASE DEFAULT
    CALL ERRMSG('D1MACH - input arg out of bounds', .TRUE.)
END SELECT

RETURN
END FUNCTION D1MACH

INTEGER FUNCTION I1MACH (I)
!
!   I1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  The results for various values of I are:
!
!   I/O unit numbers (**MAY REQUIRE USER CUSTOMIZATION**):
!     I1MACH( 1) = the standard input unit.
!     I1MACH( 2) = the standard output unit.
!     I1MACH( 3) = the standard punch unit (obsolete, will cause error)
!     I1MACH( 4) = the standard error message unit.
!                  (the error message unit is usually 0 in UNIX systems)
!
!   Words:
!     I1MACH( 5) = the number of bits per integer storage unit.
!     I1MACH( 6) = the number of characters per integer storage unit.
!                  (obsolete, will cause an error)
!
!   Integers:
!     assume integers are represented in the S-digit, base-A form
!
!          sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!     where 0 <= X(I) < A for I=0,...,S-1.
!
!     I1MACH( 7) = A, the base.
!     I1MACH( 8) = S, the number of base-A digits.
!     I1MACH( 9) = A**S - 1, the largest magnitude.
!
!   Floating-Point Numbers:
!     Assume floating-point numbers are represented in the T-digit,
!     base-B form
!                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!     where 0 <= X(I) .LT. B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
!
!     I1MACH(10) = B, the base.
!
!   Single-Precision:
!     I1MACH(11) = T, the number of base-B digits.
!     I1MACH(12) = EMIN, the smallest exponent E.
!     I1MACH(13) = EMAX, the largest exponent E.
!
!   Double-Precision:
!     I1MACH(14) = T, the number of base-B digits.
!     I1MACH(15) = EMIN, the smallest exponent E.
!     I1MACH(16) = EMAX, the largest exponent E.
!
!***REFERENCES  
!
!  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
!     ACM Transactions on Mathematical Software 4, 177-188 (1978).
!
!  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
!    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
!    August 1996. 
!
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   960411  Modified for Fortran 90 (BE after suggestions by Eric Grosse)    
! --------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: I
REAL :: X_single  = 1.0
DOUBLE PRECISION :: X_double = 1.D0

SELECT CASE (I)
  CASE (1)
    I1MACH = 5 ! Input unit
  CASE (2)
    I1MACH = 6 ! Output unit
  CASE (3)
    CALL ERRMSG('I1MACH: input arg = 3 is obsolete', .TRUE.)
  CASE (4)
    I1MACH = 0 ! Error message unit
  CASE (5)
    I1MACH = BIT_SIZE(I)
  CASE (6)
    CALL ERRMSG('I1MACH: input arg = 6 is obsolete', .TRUE.)
  CASE (7)
    I1MACH = RADIX(1)
  CASE (8)
    I1MACH = BIT_SIZE(I) - 1
  CASE (9)
    I1MACH = HUGE(1)
  CASE (10)
    I1MACH = RADIX(X_single)
  CASE (11)
    I1MACH = DIGITS(X_single)
  CASE (12)
    I1MACH = MINEXPONENT(X_single)
  CASE (13)
    I1MACH = MAXEXPONENT(X_single)
  CASE (14)
    I1MACH = DIGITS(X_double)
  CASE (15)
    I1MACH = MINEXPONENT(X_double)
  CASE (16)
    I1MACH = MAXEXPONENT(X_double) 
  CASE DEFAULT
    CALL ERRMSG('I1MACH: input argument out of bounds', .TRUE.)
END SELECT

RETURN
END FUNCTION I1MACH

SUBROUTINE ERRMSG( MESSAG, FATAL )
!
!    Print a warning or error message to the default logical unit;
!    abort if error is fatal, otherwise return to calling program.
!    (user can insert call to symbolic dump routine; see below)
!    NOTE:  FATAL is an optional argument; omitting it simply
!           prints MESSAG with no limit on number of messages
!
IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
LOGICAL, OPTIONAL, INTENT(IN) :: FATAL
CHARACTER(LEN=*), INTENT(IN) :: MESSAG
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
INTEGER, PARAMETER :: MAXMSG = 100
INTEGER :: NUMMSG = 0
LOGICAL :: MSGLIM = .FALSE.
!-----------------------------------------------

IF (.NOT.PRESENT(FATAL)) THEN

   WRITE ( *, '(/,2A,/)' ) ' ****** ERROR *****  ', MESSAG
   RETURN

END IF

IF (FATAL) THEN

   WRITE ( *, '(/,2A,/)' ) ' ****** ERROR *****  ', MESSAG

!     (insert call to symbolic dump routine here if desired)

   STOP

END IF

NUMMSG = NUMMSG + 1
IF (MSGLIM) RETURN

IF (NUMMSG <= MAXMSG) THEN

   WRITE ( *, '(/,2A,/)' ) ' ****** WARNING *****  ', MESSAG

ELSE

   MSGLIM = .TRUE.
   WRITE ( *, '(//,A,//)' ) &
' ****** TOO MANY WARNING MESSAGES -- will no longer be printed *******'

END IF

RETURN
END SUBROUTINE ERRMSG

END MODULE RDI1MACH_f90

