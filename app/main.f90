! *************************************************************************************************
!
!                               _______       ______  ____  _   ________
!                              /_  __/ |     / / __ \/ __ \/ | / /_  __/
!                               / /  | | /| / / / / / /_/ /  |/ / / /
!                              / /   | |/ |/ / /_/ / ____/ /|  / / /
!                             /_/    |__/|__/\____/_/   /_/ |_/ /_/
!
!
! MIT License
!
! (C) Federico Perini, 2023
!     A Modern Fortran port of the TWOPNT code for Boundary Value Problems.
!
!     Original code written by Dr. Joseph F. Grcar, Sandia National Laboratories, Livermore, CA.
!     VERSION 3.29 OF APRIL 1998
!
!     References:
!         J. F. Grcar, "The Twopnt Program for Boundary Value Problems, Sandia National Labs
!         Report SAND91-8230, Livermore, California, April 1992.  Reprinted February 1996.
!
! *************************************************************************************************
! SAMPLE MAIN PROGRAM FOR SIMULATING SWIRLING FLOWS
! *************************************************************************************************
program TWMAIN
    use twopnt_core
    use iso_fortran_env, only: output_unit
    implicit none

    CHARACTER NAME*16, REPORT*16, SIGNAL*16, VERSIO*80
    REAL(BK) A, ABOVE, BELOW, BUFFER, CONDIT, F, F0, G, G0, H, K, LAMBDA, &
             MU, OMEGA, RHO, RWORK, STRIDE, T, T0, TMAX, TZERO, U, U0, &
             VALUE, WMAX, X, ZMAX

    INTEGER :: ASIZE, COMPS, GROUPA, GROUPB, ISIZE, IWORK, J, LENGTH, N, NAMES, PIVOT, PMAX, &
               POINTS, RSIZE, TEXT
    LOGICAL :: ACTIVE, ERROR, MARK, RETURN, TIME

    character(len=*), parameter :: ID = 'TWMAIN:  '
    integer, parameter :: TEXT = OUTPUT_UNIT
    integer, parameter :: ISIZE  = 5000, RSIZE = 100000
    integer, parameter :: COMPS  = 5
    integer, parameter :: GROUPA = 0
    integer, parameter :: GROUPB = 0
    integer, parameter :: PMAX   = 200
    integer, parameter :: ASIZE = (6 * COMPS - 1) * COMPS * PMAX
    integer, parameter :: NAMES = COMPS + GROUPA + GROUPB

      DIMENSION
     +   A(ASIZE), ABOVE(COMPS), ACTIVE(COMPS), BELOW(COMPS),
     +   BUFFER(COMPS * PMAX), F(PMAX), F0(PMAX), G(PMAX), G0(PMAX),
     +   H(PMAX), IWORK(ISIZE), K(PMAX), LAMBDA(PMAX), MARK(PMAX),
     +   MU(PMAX), NAME(NAMES), PIVOT(COMPS * PMAX), RHO(PMAX),
     +   RWORK(RSIZE), T(PMAX), T0(PMAX), U(COMPS, PMAX), U0(COMPS *
     +   PMAX), X(PMAX)

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  OPEN FILES.

C      OPEN (FILE = 'twopnt.out', STATUS='UNKNOWN',
C     1      FORM = 'FORMATTED', UNIT = TEXT)

C///////////////////////////////////////////////////////////////////////
C
C     SET PROBLEM PARAMETERS.
C
C///////////////////////////////////////////////////////////////////////

C     ROTATION RATE AT POROUS DISK
      OMEGA = 4.0 * 3.14159

C     GAS TEMPERATURE AT POROUS DISK
      TMAX = 1000.0

C     GAS TEMPERATURE AT SOLID DISK
      TZERO = 300.0

C     FLOW AT POROUS DISK
      WMAX = - 2.0

C     DISTANCE BETWEEN DISKS
      ZMAX = 5.0

C///////////////////////////////////////////////////////////////////////
C
C     SET TWOPNT CONTROLS.
C
C///////////////////////////////////////////////////////////////////////

C///  CHOOSE THE INITIAL GRID SIZE.

      POINTS = 6
      N = GROUPA + COMPS * POINTS + GROUPB

C///  SPECIFY THE CONTROLS.

      CALL TWSETL (ERROR, TEXT, 'ADAPT', .TRUE.)
      IF (ERROR) GO TO 9001

      CALL TWSETI (ERROR, TEXT, 'LEVELD', 1)
      IF (ERROR) GO TO 9002

      CALL TWSETI (ERROR, TEXT, 'LEVELM', 1)
      IF (ERROR) GO TO 9002

      CALL TWSETI (ERROR, TEXT, 'STEPS1', 50)
      IF (ERROR) GO TO 9002

      CALL TWSETI (ERROR, TEXT, 'STEPS2', 50)
      IF (ERROR) GO TO 9002

      VALUE = 1.0E-3
      CALL TWSETR (ERROR, TEXT, 'STRID0', VALUE)
      IF (ERROR) GO TO 9003

      VALUE = 3.16
      CALL TWSETR (ERROR, TEXT, 'TINC', VALUE)
      IF (ERROR) GO TO 9003

      VALUE = 0.1
      CALL TWSETR (ERROR, TEXT, 'TOLER1', VALUE)
      IF (ERROR) GO TO 9003

      VALUE = 0.1
      CALL TWSETR (ERROR, TEXT, 'TOLER2', VALUE)
      IF (ERROR) GO TO 9003

C*****RELAXED TOLERANCES
C      VALUE = 1.0E-9
C      CALL TWSETR (ERROR, TEXT, 'SSABS', VALUE)
C      IF (ERROR) GO TO 9003
C
C      VALUE = 1.0E-9
C      CALL TWSETR (ERROR, TEXT, 'TDABS', VALUE)
C      IF (ERROR) GO TO 9003
C
C      VALUE = 1.0E-4
C      CALL TWSETR (ERROR, TEXT, 'SSREL', VALUE)
C      IF (ERROR) GO TO 9003
C
C      VALUE = 1.0E-4
C      CALL TWSETR (ERROR, TEXT, 'TDREL', VALUE)
C      IF (ERROR) GO TO 9003
C*****END RELAXED TOLERANCES

C///////////////////////////////////////////////////////////////////////
C
C     INITIALIZE ARRAYS FOR TWOPNT.
C
C///////////////////////////////////////////////////////////////////////

C///  FORM THE INITIAL GRID.

      DO 0100 J = 1, POINTS
         X(J) = (REAL (J - 1) / REAL (POINTS - 1)) * ZMAX
0100  CONTINUE

C///  CHOOSE GUESSES FOR UNKNOWNS.

      DO 0200 J = 1, POINTS
C        F
         U(1, J) = 0.0

C        G
         U(2, J) = OMEGA * X(J) / ZMAX

C        H
         U(3, J) = WMAX * (1.0 - X(J) / ZMAX)

C        LAMBDA
         U(4, J) = 0.0

C        T
         U(5, J) = (1.0 - X(J) / ZMAX) * TZERO + (X(J) / ZMAX) * TMAX
0200  CONTINUE

C///  ASSIGN LIMITS FOR UNKNOWNS.

      ABOVE(1) = 4.0
      BELOW(1) = - 4.0
      ABOVE(2) = 1.0E4
      BELOW(2) = - 1.0E4
      ABOVE(3) = 1.0E4
      BELOW(3) = - 1.0E4
      ABOVE(4) = 1.0E4
      BELOW(4) = - 1.0E4
      ABOVE(5) = 2.0 * TMAX
      BELOW(5) = 0.5 * TZERO

C///  ASSIGN NAMES FOR UNKNOWNS.

      NAME(1) = 'F'
      NAME(2) = 'G'
      NAME(3) = 'H'
      NAME(4) = 'LAMBDA'
      NAME(5) = 'T'

C///  CHOOSE UNKNOWNS TO EXAMINE FOR GRID ADAPTION.

      ACTIVE(1) = .TRUE.
      ACTIVE(2) = .TRUE.
      ACTIVE(3) = .TRUE.
      ACTIVE(4) = .FALSE.
      ACTIVE(5) = .TRUE.

C///////////////////////////////////////////////////////////////////////
C
C     CALL TWOPNT.
C
C///////////////////////////////////////////////////////////////////////

C*****PRECISION > DOUBLE
      VERSIO = 'DOUBLE PRECISION VERSION 3.22'
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      VERSIO = 'SINGLE PRECISION VERSION 3.22'
C*****END PRECISION > SINGLE

C     SUBROUTINE TWOPNT
C    +  (ERROR, TEXT, VERSIO,
C    +   ABOVE, ACTIVE, BELOW, BUFFER, COMPS, CONDIT, GROUPA, GROUPB,
C    +   ISIZE, IWORK, MARK, NAME, NAMES, PMAX, POINTS, REPORT, RSIZE,
C    +   RWORK, SIGNAL, STRIDE, TIME, U, X)

      SIGNAL = ' '
0300  CONTINUE
      CALL TWOPNT
     +  (ERROR, TEXT, VERSIO,
     +   ABOVE, ACTIVE, BELOW, BUFFER, COMPS, CONDIT, GROUPA, GROUPB,
     +   ISIZE, IWORK, MARK, NAME, NAMES, PMAX, POINTS, REPORT, RSIZE,
     +   RWORK, SIGNAL, STRIDE, TIME, U, X)
      IF (ERROR) GO TO 9004

C///////////////////////////////////////////////////////////////////////
C
C     SERVICE REQUESTS FROM TWOPNT.
C
C///////////////////////////////////////////////////////////////////////

C///  EVALUATE THE FUNCTION.

      IF (SIGNAL .EQ. 'RESIDUAL') THEN

C     SUBROUTINE TWFUNC
C    +  (ERROR, TEXT,
C    +   BUFFER, F, F0, G, G0, H, K, LAMBDA, MU, OMEGA, POINTS, RHO,
C    +   STRIDE, T, T0, TIME, TMAX, TZERO, U0, WMAX, X)

      CALL TWFUNC
     +  (ERROR, TEXT,
     +   BUFFER, F, F0, G, G0, H, K, LAMBDA, MU, OMEGA, POINTS, RHO,
     +   STRIDE, T, T0, TIME, TMAX, TZERO, U0, WMAX, X)
      IF (ERROR) GO TO 9005

C///  EVALUATE AND FACTOR THE JACOBIAN.

      ELSE IF (SIGNAL .EQ. 'PREPARE') THEN

      RETURN = .FALSE.
0400  CONTINUE

C     SUBROUTINE TWPREP
C    +  (ERROR, TEXT,
C    +   A, ASIZE, BUFFER, COMPS, CONDIT, GROUPA, GROUPB, PIVOT, POINTS,
C    +   RETURN)

      CALL TWPREP
     +  (ERROR, TEXT,
     +   A, ASIZE, BUFFER, COMPS, CONDIT, GROUPA, GROUPB, PIVOT, POINTS,
     +   RETURN)
      IF (ERROR) GO TO 9006

      IF (RETURN) THEN

C     SUBROUTINE TWFUNC
C    +  (ERROR, TEXT,
C    +   BUFFER, F, F0, G, G0, H, K, LAMBDA, MU, OMEGA, POINTS, RHO,
C    +   STRIDE, T, T0, TIME, TMAX, TZERO, U0, WMAX, X)

      CALL TWFUNC
     +  (ERROR, TEXT,
     +   BUFFER, F, F0, G, G0, H, K, LAMBDA, MU, OMEGA, POINTS, RHO,
     +   STRIDE, T, T0, TIME, TMAX, TZERO, U0, WMAX, X)
      IF (ERROR) GO TO 9005

      GO TO 0400
      END IF

C///  SHOW THE SOLUTION.

C     SUBROUTINE TWSHOW
C    +  (ERROR, TEXT,
C    +   BUFFER, COMPS, GRID, GROUPA, GROUPB, POINTS, X)

      ELSE IF (SIGNAL .EQ. 'SHOW') THEN
         CALL TWSHOW
     +     (ERROR, TEXT, BUFFER, COMPS, .TRUE., GROUPA, GROUPB,
     +      POINTS, X)
         IF (ERROR) GO TO 9007

C///  SOLVE THE LINEAR EQUATIONS.

      ELSE IF (SIGNAL .EQ. 'SOLVE') THEN

C     SUBROUTINE TWSOLV
C    +  (ERROR, TEXT,
C    +   A, ASIZE, BUFFER, COMPS, GROUPA, GROUPB, PIVOT, POINTS)

      CALL TWSOLV
     +  (ERROR, TEXT,
     +   A, ASIZE, BUFFER, COMPS, GROUPA, GROUPB, PIVOT, POINTS)
      IF (ERROR) GO TO 9008

C///  RETAIN THE SOLUTION FOR TIME INTEGRATION.

      ELSE IF (SIGNAL .EQ. 'RETAIN') THEN
         CALL TWCOPY (N, BUFFER, U0)

C///  UPDATE THE GRID.

      ELSE IF (SIGNAL .EQ. 'UPDATE') THEN
         N = GROUPA + COMPS * POINTS + GROUPB

C///  BOTTOM OF THE BLOCK TO SERVICE REQUESTS.

      END IF
      IF (SIGNAL .NE. ' ') GO TO 0300

C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  CHECK FOR SUCCESS.

      ERROR = REPORT .EQ. 'NONE FOUND'
      IF (ERROR) GO TO 9009

C///  WRITE A SUMMARY.

      WRITE (TEXT, 10001) ID, U(4, 1), OMEGA, TZERO, TMAX, WMAX

C///  WRITE THE DATA FOR PLOTTING.

C     SUBROUTINE TWPLOT
C    +  (ERROR, TEXT,
C    +   COMPS, DATA, POINTS, U, X)

C      CALL TWPLOT
C     +  (ERROR, TEXT,
C     +   COMPS, 16, POINTS, U, X)
      IF (ERROR) GO TO 9010



C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      GO TO 99999

9001  IF (0 .LT. TEXT) WRITE (TEXT, 99001) ID
      GO TO 99999

9002  IF (0 .LT. TEXT) WRITE (TEXT, 99002) ID
      GO TO 99999

9003  IF (0 .LT. TEXT) WRITE (TEXT, 99003) ID
      GO TO 99999

9004  IF (0 .LT. TEXT) WRITE (TEXT, 99004) ID
      GO TO 99999

9005  IF (0 .LT. TEXT) WRITE (TEXT, 99005) ID
      GO TO 99999

9006  IF (0 .LT. TEXT) WRITE (TEXT, 99006) ID
      GO TO 99999

9007  IF (0 .LT. TEXT) WRITE (TEXT, 99007) ID
      GO TO 99999

9008  IF (0 .LT. TEXT) WRITE (TEXT, 99008) ID
      GO TO 99999

9009  IF (0 .LT. TEXT) THEN
         CALL twlast(LENGTH, REPORT)
         WRITE (TEXT, 99009) ID, REPORT(1:LENGTH)
      END IF
      GO TO 99999

9010  IF (0 .LT. TEXT) WRITE (TEXT, 99010) ID
      GO TO 99999


    return

    ! INFORMATIVE MESSAGES.
    10001 FORMAT(/1X, A9, 1P, E10.2, 0P, '  LAMBDA' &
                 /10X, F10.5, '  OMEGA' &
                 /10X, F10.5, '  TZERO' &
                 /10X, F10.5, '  TMAX'  &
                 /10X, F10.5, '  WMAX')

    ! ERROR MESSAGES.
    99001 FORMAT(/1X, A9, 'ERROR.  TWSETL FAILS.')
    99002 FORMAT(/1X, A9, 'ERROR.  TWSETI FAILS.')
    99003 FORMAT(/1X, A9, 'ERROR.  TWSETR FAILS.')
    99004 FORMAT(/1X, A9, 'ERROR.  TWOPNT FAILS.')
    99005 FORMAT(/1X, A9, 'ERROR.  RESID FAILS.')
    99006 FORMAT(/1X, A9, 'ERROR.  TWPREP FAILS.')
    99007 FORMAT(/1X, A9, 'ERROR.  TWSHOW FAILS.')
    99008 FORMAT(/1X, A9, 'ERROR.  TWSOLV FAILS.')
    99009 FORMAT(/1X, A9, 'ERROR.  TWOPNT DOES NOT SOLVE THE PROBLEM.'//10X, '   REPORT:  ', A)
    99010 FORMAT(/1X, A9, 'ERROR.  TWPLOT FAILS.')

    contains

      ! Compute density, viscosity and thermal conductivity of Argon at given T
      elemental subroutine AR_TRANSPORT(T,RHO,MU,K)
          real(RK), intent(in)  :: T
          real(RK), intent(out) :: RHO,MU,K

          ! PRESSURE AT 1 STANDARD ATMOSPHERE
          real(RK), PARAMETER :: P = 1013250.0_RK

          ! MOLECULAR WEIGHT OF ARGON
          real(RK), PARAMETER :: W = 39.948_RK

          ! UNIVERSAL GAS CONSTANT
          real(RK), PARAMETER :: R = 83140000.0_RK

          ! SPECIFIC HEAT OF ARGON AT CONSTANT PRESSURE [ERGS / (GM * K)]
          real(RK), PARAMETER :: CP = R * 2.5_RK / W

          real(RK) :: A

          RHO = (P * W) / (R * T)
          A   = LOG(T)
          MU  = EXP(((0.0121673_RK * A - 0.284023_RK) * A + 2.85205_RK) * A - 17.6539_RK)
          K   = EXP(((0.0121673_RK * A - 0.284023_RK) * A + 2.85205_RK) * A - 1.78377_RK)

      end subroutine AR_TRANSPORT

      ! SAMPLE RESIDUAL FUNCTION FOR SIMULATING SWIRLING FLOWS.
      subroutine TWFUNC(ERROR, TEXT, &
                        BUFFER, F, F0, G, G0, H, K, LAMBDA, MU, OMEGA, POINTS, RHO, &
                        STRIDE, T, T0, TIME, TMAX, TZERO, U0, WMAX, X)

          character(*), parameter :: id = 'TWFUNC:  '
          integer     , parameter :: COMPS = 5 ! Number of variables

          integer, intent(in)  :: POINTS ! Number of 1D variables
          integer, intent(in)  :: TEXT   ! output unit (0 for NONE)
          logical, intent(out) :: ERROR
          logical, intent(in)  :: TIME

          real(RK) :: A0, A1, A2, B0, B1, B2, C0, C1, C2, CP, OMEGA, P, R, STRIDE, TEMP, &
                      TMAX, TZERO, W, WMAX
          real(RK), dimension(POINTS) :: F, F0, G, G0, H, K, LAMBDA, MU, RHO, T, T0, X
          real(RK), dimension(COMPS, POINTS) :: BUFFER, U0(COMPS, POINTS)
          INTEGER :: J,
          INTRINSIC :: EXP, LOG

          ! Check the arguments
          ERROR = .NOT. POINTS>=2
          if (ERROR) then
             if (TEXT>0) WRITE (TEXT, 99001) ID, POINTS
             return
          end if

          ! UNPACK THE VARIABLES.
          unroll_variables: DO J = 1, POINTS
             F(J)      = BUFFER(1, J)
             G(J)      = BUFFER(2, J)
             H(J)      = BUFFER(3, J)
             LAMBDA(J) = BUFFER(4, J)
             T(J)      = BUFFER(5, J)
             F0(J)     = U0(1, J)
             G0(J)     = U0(2, J)
             T0(J)     = U0(5, J)
          end do unroll_variables

          ! CALCULATE DENSITIES, VISCOSITIES AND THERMAL CONDUCTIVITIES.
          call AR_TRANSPORT(T,RHO,MU,K)

          ! F, G AND T EQUATIONS 1) 2) 3)
          ! RHO * (dF/dt + F ** 2 - G ** 2 + F' H) - (MU F')' + P = 0
          ! RHO * (dG/dt + 2 F G + G' H) - (MU G')' = 0
          ! RHO * CP * (dT/dt + H T') - (K T')' = 0

          ! BOUNDARY CONDITIONS on points 1 and POINTS
          BUFFER(1, 1)      = F(1) - zero
          BUFFER(1, POINTS) = F(POINTS) - zero

          BUFFER(2, 1)      = G(1) - OMEGA
          BUFFER(2, POINTS) = G(POINTS) - zero

          BUFFER(5, 1)      = T(1) - TZERO
          BUFFER(5, POINTS) = T(POINTS) - TMAX

          ! EQUATIONS.
          equations: DO J = 2, POINTS - 1
             TEMP = ((X(J - 1) - X(J)) * (X(J - 1) - X(J + 1)))
             A0 = (X(J) - X(J + 1)) / TEMP
             B0 = (MU(J) + MU(J - 1)) / TEMP
             C0 = (K(J) + K(J - 1)) / TEMP

             TEMP = ((X(J) - X(J + 1)) * (X(J - 1) - X(J + 1)))
             A2 = (X(J) - X(J - 1)) / TEMP
             B2 = (MU(J + 1) + MU(J)) / TEMP
             C2 = (K(J + 1) + K(J)) / TEMP

             A1 = - (A0 + A2)
             B1 = - (B0 + B2)
             C1 = - (C0 + C2)

             BUFFER(1, J) = RHO(J) * (F(J) ** 2 - G(J) ** 2 + H(J) &
                          * (A0 * F(J - 1) + A1 * F(J) + A2 * F(J + 1))) &
                          - (B0 * F(J - 1) + B1 * F(J) + B2 * F(J + 1)) &
                          + LAMBDA(J)

             BUFFER(2, J) = RHO(J) * (2.0 * F(J) * G(J) + H(J) &
                          * (A0 * G(J - 1) + A1 * G(J) + A2 * G(J + 1))) &
                          - (B0 * G(J - 1) + B1 * G(J) + B2 * G(J + 1))

             BUFFER(5, J) = RHO(J) * CP * H(J) &
                          * (A0 * T(J - 1) + A1 * T(J) + A2 * T(J + 1)) &
                          - (C0 * T(J - 1) + C1 * T(J) + C2 * T(J + 1))

             IF (TIME) THEN
                 BUFFER(1, J) = BUFFER(1, J) + RHO(J) * (F(J) - F0(J)) / STRIDE
                 BUFFER(2, J) = BUFFER(2, J) + RHO(J) * (G(J) - G0(J)) / STRIDE
                 BUFFER(5, J) = BUFFER(5, J) + RHO(J) * CP * (T(J) - T0(J)) / STRIDE
             END IF
          end do equations

          ! H AND P EQUATIONS 4) 5)
          !
          ! 2 F + H' + H * (LOG RHO)' = 0
          !
          ! P CONSTANT

          ! BOUNDARY CONDITION.
          BUFFER(3,      1) = H(1) - zero
          BUFFER(4, POINTS) = H(POINTS) - WMAX

          ! EQUATIONS.
          do J = 2, POINTS
             BUFFER(3, J) = F(J) + F(J - 1) &
                          + (H(J) - H(J - 1)) / (X(J) - X(J - 1)) &
                          + 0.5 * (H(J) + H(J - 1)) &
                          * (LOG (RHO(J)) - LOG (RHO(J - 1))) / (X(J) - X(J - 1))

             BUFFER(4, J - 1) = LAMBDA(J) - LAMBDA(J - 1)
          end do

          return

          ! ERROR MESSAGES.
          99001 FORMAT(/1X, A9, 'ERROR.  THERE MUST BE AT LEAST TWO POINTS.'//10X, I10, '  POINTS')

      end subroutine TWFUNC



end program TWMAIN
