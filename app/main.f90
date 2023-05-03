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
!         J. F. Grcar, "The Twopnt Program for Boundary Value Problems", Sandia National Labs
!         Report SAND91-8230, Livermore, California, April 1992.  Reprinted February 1996.
!
! *************************************************************************************************
! SAMPLE MAIN PROGRAM FOR SIMULATING SWIRLING FLOWS
! *************************************************************************************************
program TWMAIN
    use twopnt_core
    use iso_fortran_env, only: output_unit
    implicit none

    ! Parameters
    character(len=*), parameter :: ID = 'TWMAIN:  '
    integer, parameter :: TEXT = OUTPUT_UNIT
    integer, parameter :: COMPS  = 5
    integer, parameter :: GROUPA = 0
    integer, parameter :: GROUPB = 0
    integer, parameter :: PMAX   = 200
    integer, parameter :: ASIZE = (6 * COMPS - 1) * COMPS * PMAX
    integer, parameter :: NAMES = COMPS + GROUPA + GROUPB
    logical, parameter :: RELAXED_TOLERANCES = .false.

    character(len=16) :: NAME(NAMES), REPORT, SIGNAL
    character(len=80) :: VERSIO

    type(twcom) :: settings
    type(twsize) :: sizes
    type(twwork) :: work
    type(twfunctions) :: problem
    real(RK) :: A(ASIZE), ABOVE(COMPS), BELOW(COMPS)
    real(RK), dimension(COMPS*PMAX) :: BUFFER
    real(RK), dimension(COMPS,PMAX) :: U,U0
    real(RK), dimension(PMAX) :: F,F0,G,G0,H,K,LAMBDA,MU,RHO,T,T0,X
    real(RK) :: CONDIT,STRIDE
    logical  :: ACTIVE(COMPS), MARK(PMAX)
    integer  :: PIVOT(COMPS*PMAX)
    INTEGER :: J, LENGTH, N
    LOGICAL :: ERROR, TIME

    ! PROLOGUE. OPEN FILES.
!   OPEN (FILE = 'twopnt.out', STATUS='UNKNOWN',FORM = 'FORMATTED', UNIT = TEXT)

    ! *** PROBLEM PARAMETERS. ***

    ! ROTATION RATE AT POROUS DISK
    real(RK), parameter :: OMEGA = 4.0_RK * pi

    ! GAS TEMPERATURE AT POROUS DISK
    real(RK), parameter :: TMAX = 1000.0_RK

    ! GAS TEMPERATURE AT SOLID DISK
    real(RK), parameter :: TZERO = 300.0_RK

    ! FLOW AT POROUS DISK
    real(RK), parameter :: WMAX = - 2.0_RK

    ! DISTANCE BETWEEN DISKS
    real(RK), parameter :: ZMAX = 5.0_RK

    ! *** SET TWOPNT CONTROLS. ***

    ! CHOOSE THE INITIAL GRID SIZE.
    sizes = twsize(GROUPA,COMPS,6,PMAX,GROUPB)
    N = sizes%N()

    ! SPECIFY THE CONTROLS.
    call settings%set(ERROR, TEXT, 'ADAPT', .TRUE.)
    call settings%set(ERROR, TEXT, 'LEVELD', 1)
    call settings%set(ERROR, TEXT, 'LEVELM', 1)
    call settings%set(ERROR, TEXT, 'STEPS1', 50)
    call settings%set(ERROR, TEXT, 'STEPS2', 50)
    call settings%set(ERROR, TEXT, 'STRID0', 1.0e-3_RK)
    call settings%set(ERROR, TEXT, 'TINC', 3.16_RK)
    call settings%set(ERROR, TEXT, 'TOLER1', 0.1_RK)
    call settings%set(ERROR, TEXT, 'TOLER2', 0.1_RK)
    if (RELAXED_TOLERANCES) then
        call settings%set(ERROR, TEXT, 'SSABS', 1.0e-9_RK)
        call settings%set(ERROR, TEXT, 'TDABS', 1.0e-9_RK)
        call settings%set(ERROR, TEXT, 'SSREL', 1.0e-4_RK)
        call settings%set(ERROR, TEXT, 'TDREL', 1.0e-4_RK)
    endif

    ! INITIALIZE ARRAYS FOR TWOPNT.

    ! FORM THE INITIAL GRID.
    FORALL(J=1:sizes%POINTS) X(J) = ZMAX*(REAL(J-1,RK)/REAL(sizes%POINTS-1,RK))

    ! CHOOSE GUESSES FOR THE UNKNOWNS.
    init_unknowns: DO J = 1, sizes%POINTS

         ! F
         U(1, J) = zero

         ! G
         U(2, J) = OMEGA * X(J) / ZMAX

         ! H
         U(3, J) = WMAX * (one - X(J) / ZMAX)

         ! LAMBDA
         U(4, J) = zero

         ! T
         U(5, J) = (one - X(J) / ZMAX) * TZERO + (X(J) / ZMAX) * TMAX

    end do init_unknowns

    ! ASSIGN LIMITS FOR THE UNKNOWNS.
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

    ! ASSIGN NAMES FOR THE UNKNOWNS.
    NAME(1) = 'F'
    NAME(2) = 'G'
    NAME(3) = 'H'
    NAME(4) = 'LAMBDA'
    NAME(5) = 'T'

    ! CHOOSE UNKNOWNS TO EXAMINE FOR GRID ADAPTION.
    ACTIVE(1) = .TRUE.
    ACTIVE(2) = .TRUE.
    ACTIVE(3) = .TRUE.
    ACTIVE(4) = .FALSE.
    ACTIVE(5) = .TRUE.

    ! Initialize functions
    problem%save_sol => savesol
    problem%fun      => residual

    ! CALL TWOPNT.
    VERSIO = 'DOUBLE PRECISION VERSION 3.22'
    SIGNAL = ' '

    ITERATE: DO

          ! Call driver
          CALL TWOPNT(SETTINGS, ERROR, TEXT, VERSIO, sizes, ABOVE, ACTIVE, BELOW, BUFFER, CONDIT, &
                      WORK, MARK, NAME, NAMES, REPORT, SIGNAL, STRIDE, TIME, U, X, problem)

          IF (ERROR) GO TO 9004

          ! SERVICE REQUESTS FROM TWOPNT.
          select case (SIGNAL)

             case ('RESIDUAL')

                ! Evaluate residual
                call residual(error,text,sizes%points,time,stride,x,buffer)
                IF (ERROR) GO TO 9005

             case ('PREPARE')

                ! EVALUATE AND FACTOR THE JACOBIAN
                CALL problem%prep(ERROR, TEXT, A, ASIZE, BUFFER, sizes, CONDIT, PIVOT, time, stride, x)
                IF (ERROR) GO TO 9006

             case ('SHOW')

                 ! SHOW THE SOLUTION.
                 CALL TWSHOW(ERROR, TEXT, BUFFER, sizes, .TRUE., X)
                 IF (ERROR) GO TO 9007

             case ('SOLVE')

                 ! SOLVE THE LINEAR EQUATIONS.
                 CALL TWSOLV(ERROR, TEXT, A, ASIZE, BUFFER, sizes, PIVOT)
                 IF (ERROR) GO TO 9008

             case ('UPDATE')

                 ! UPDATE THE GRID.
                 N = sizes%N()

             case (' ')

                 ! Iteration finished
                 exit ITERATE

          end select

    END DO ITERATE

    ! CHECK FOR SUCCESS.
    ERROR = REPORT == 'NONE FOUND'
    IF (ERROR) GO TO 9009

    ! WRITE A SUMMARY.
    WRITE (TEXT, 10001) ID, U(4, 1), OMEGA, TZERO, TMAX, WMAX

    ! SUCCESSFUL RETURN.
    return

    ! ERROR HANDLING.
    9004  IF (TEXT>0) WRITE (TEXT, 99004) ID; return
    9005  IF (TEXT>0) WRITE (TEXT, 99005) ID; return
    9006  IF (TEXT>0) WRITE (TEXT, 99006) ID; return
    9007  IF (TEXT>0) WRITE (TEXT, 99007) ID; return
    9008  IF (TEXT>0) WRITE (TEXT, 99008) ID; return
    9009  IF (TEXT>0) THEN
             CALL twlast(LENGTH, REPORT)
             WRITE (TEXT, 99009) ID, REPORT(1:LENGTH)
          END IF
    return

    ! INFORMATIVE MESSAGES.
    10001 FORMAT(/1X, A9, 1P, E10.2, 0P, '  LAMBDA' &
                 /10X, F10.5, '  OMEGA' &
                 /10X, F10.5, '  TZERO' &
                 /10X, F10.5, '  TMAX'  &
                 /10X, F10.5, '  WMAX')

    ! ERROR MESSAGES.
    99004 FORMAT(/1X, A9, 'ERROR.  TWOPNT FAILS.')
    99005 FORMAT(/1X, A9, 'ERROR.  RESID FAILS.')
    99006 FORMAT(/1X, A9, 'ERROR.  TWPREP FAILS.')
    99007 FORMAT(/1X, A9, 'ERROR.  TWSHOW FAILS.')
    99008 FORMAT(/1X, A9, 'ERROR.  TWSOLV FAILS.')
    99009 FORMAT(/1X, A9, 'ERROR.  TWOPNT DOES NOT SOLVE THE PROBLEM.'//10X, '   REPORT:  ', A)

    contains

      !  Save function interface
      subroutine savesol(vars,buffer)
         type(twsize), intent(in) :: vars
         real(RK),     intent(in) :: buffer(vars%N())


         integer :: N

         N = vars%N()

         ! RETAIN THE SOLUTION FOR TIME INTEGRATION.
         CALL TWCOPY (N,BUFFER,U0)

      end subroutine savesol

      !  Residual function interface:
      subroutine residual(error,text,points,time,stride,x,buffer)
         logical, intent(out) :: error  ! .true. if something went wrong
         integer, intent(in)  :: text   ! output unit (0 for NONE)
         integer, intent(in)  :: points ! Number of 1D variables
         logical, intent(in)  :: time   ! Time-resolved or steady state
         real(RK), intent(in) :: stride !
         real(RK), intent(in) :: x(*)   ! dimensioned >=PMAX, contains the grid
         real(RK), intent(inout) :: buffer(*) ! on input: contains the approximate solution

         CALL TWFUNC(ERROR, TEXT, &
                     BUFFER, F, F0, G, G0, H, K, LAMBDA, MU, OMEGA, POINTS, RHO, &
                     STRIDE, T, T0, TIME, TMAX, TZERO, U0, WMAX, X)

      end subroutine residual

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

          real(RK) :: A0, A1, A2, B0, B1, B2, C0, C1, C2, OMEGA, STRIDE, TEMP, TMAX, TZERO, WMAX
          real(RK), dimension(POINTS) :: F, F0, G, G0, H, K, LAMBDA, MU, RHO, T, T0, X
          real(RK), dimension(COMPS, POINTS) :: BUFFER, U0(COMPS, POINTS)
          INTEGER :: J
          INTRINSIC :: EXP, LOG

          ! PRESSURE AT 1 STANDARD ATMOSPHERE
          real(RK), PARAMETER :: P = 1013250.0_RK

          ! MOLECULAR WEIGHT OF ARGON
          real(RK), PARAMETER :: W = 39.948_RK

          ! UNIVERSAL GAS CONSTANT
          real(RK), PARAMETER :: R = 83140000.0_RK

          ! SPECIFIC HEAT OF ARGON AT CONSTANT PRESSURE [ERGS / (GM * K)]
          real(RK), PARAMETER :: CP = R * 2.5_RK / W

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
