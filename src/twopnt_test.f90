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
!
! *************************************************************************************************
! SAMPLE BVP IMPLEMENTATION FOR SIMULATING SWIRLING FLOWS
! *************************************************************************************************
module twopnt_test
    use iso_fortran_env, only: real64
    use twopnt_core
    implicit none
    private

    !> Derive a BVP class for this problem
    type, public, extends(TwoPntBVProblem) :: SwirlingFlow

        ! ROTATION RATE AT POROUS DISK
        real(RK) :: OMEGA = 4.0_RK * pi

        ! GAS TEMPERATURE AT POROUS DISK
        real(RK) :: TMAX = 1000.0_RK

        ! GAS TEMPERATURE AT SOLID DISK
        real(RK) :: TZERO = 300.0_RK

        ! FLOW AT POROUS DISK
        real(RK) :: WMAX = - 2.0_RK

        ! DISTANCE BETWEEN DISKS
        real(RK) :: ZMAX = 5.0_RK

        integer  :: COMPS  = 5
        integer  :: GROUPA = 0
        integer  :: GROUPB = 0
        integer  :: PMAX   = 200
        logical  :: RELAXED_TOLERANCES = .false.

        ! Problem fixed variables and buffers
        integer :: N
        real(RK), allocatable :: ABOVE(:)
        real(RK), allocatable :: BELOW(:)
        real(RK), allocatable :: U(:,:),U0(:,:)
        real(RK), allocatable, dimension(:) :: F,F0,G,G0,H,K,LAMBDA,MU,RHO,T,T0
        logical , allocatable, dimension(:) :: ACTIVE, MARK

        contains

            ! Initialize problem
            procedure :: new => swirling_new

            ! Print summary
            procedure :: summary => swirl_summary


            ! Problem functions
            procedure :: save_sol    => savesol
            procedure :: fun         => residual
            procedure :: update_grid => on_grid_update

    end type SwirlingFlow

    contains

      ! Initialize swirling flow problem
      subroutine swirling_new(this,error,text)
         class(SwirlingFlow), intent(inout) :: this
         logical, intent(out) :: error
         integer, intent(in) :: text

         integer :: J

         allocate(this%F(this%PMAX))
         allocate(this%F0(this%PMAX))
         allocate(this%G(this%PMAX))
         allocate(this%G0(this%PMAX))
         allocate(this%H(this%PMAX))
         allocate(this%K(this%PMAX))
         allocate(this%LAMBDA(this%PMAX))
         allocate(this%MU(this%PMAX))
         allocate(this%RHO(this%PMAX))
         allocate(this%T(this%PMAX))
         allocate(this%T0(this%PMAX))
         allocate(this%MARK(this%PMAX))
         allocate(this%U(this%COMPS,this%PMAX),this%U0(this%COMPS,this%PMAX))
         allocate(this%ACTIVE(this%COMPS),source=.true.)
         allocate(this%BELOW(this%COMPS),this%ABOVE(this%COMPS))

         associate(ACTIVE=>this%ACTIVE,BELOW=>this%BELOW,ABOVE=>this%ABOVE,settings=>this%setup,&
                   domain=>this%domain,U=>this%U,TZERO=>this%TZERO,TMAX=>this%TMAX,OMEGA=>this%OMEGA,&
                   WMAX=>this%WMAX,ZMAX=>this%ZMAX)

         ERROR = .false.

         ! CHOOSE UNKNOWNS TO EXAMINE FOR GRID ADAPTION.
         ACTIVE(4) = .FALSE.

         ! ASSIGN LIMITS FOR THE UNKNOWNS.
         BELOW(1) = - 4.0;            ABOVE(1) = 4.0
         BELOW(2) = - 1.0E4;          ABOVE(2) = 1.0E4
         BELOW(3) = - 1.0E4;          ABOVE(3) = 1.0E4
         BELOW(4) = - 1.0E4;          ABOVE(4) = 1.0E4
         BELOW(5) = 0.5 * TZERO;      ABOVE(5) = 2.0 * TMAX

         ! *** SET TWOPNT CONTROLS. ***

         ! ASSIGN INITIAL PROBLEM SIZES AND NAMES FOR THE UNKNOWNS.
         call domain%new(ERROR, TEXT,this%GROUPA,this%COMPS,6,this%PMAX,this%GROUPB,ABOVE,BELOW, &
                        [character(6) :: 'F','G','H','LAMBDA','T'], &
                        XRANGE=[zero,this%ZMAX], ACTIVE=ACTIVE)
         if (error) return

         ! CHOOSE THE INITIAL GRID SIZE.
         this%N = domain%N()

         ! SPECIFY THE CONTROLS.
         call settings%set(ERROR, TEXT, 'ADAPT', .true.)
         call settings%set(ERROR, TEXT, 'LEVELD', 1)
         call settings%set(ERROR, TEXT, 'LEVELM', 1)
         call settings%set(ERROR, TEXT, 'STEPS1', 50)
         call settings%set(ERROR, TEXT, 'STEPS2', 50)
         call settings%set(ERROR, TEXT, 'STRID0', 1.0e-3_RK)
         call settings%set(ERROR, TEXT, 'TINC', 3.16_RK)
         call settings%set(ERROR, TEXT, 'TOLER1', 0.1_RK)
         call settings%set(ERROR, TEXT, 'TOLER2', 0.1_RK)
         if (this%RELAXED_TOLERANCES) then
             call settings%set(ERROR, TEXT, 'SSABS', 1.0e-9_RK)
             call settings%set(ERROR, TEXT, 'TDABS', 1.0e-9_RK)
             call settings%set(ERROR, TEXT, 'SSREL', 1.0e-4_RK)
             call settings%set(ERROR, TEXT, 'TDREL', 1.0e-4_RK)
         endif


         ! CHOOSE GUESSES FOR THE UNKNOWNS.
         init_unknowns: DO J = 1, domain%POINTS
             U(1, J) = zero ! F
             U(2, J) = OMEGA * domain%X(J) / ZMAX ! G
             U(3, J) = WMAX * (one - domain%X(J) / ZMAX) ! H
             U(4, J) = zero ! LAMBDA
             U(5, J) = (one - domain%X(J) / ZMAX) * TZERO + (domain%X(J) / ZMAX) * TMAX ! T
         end do init_unknowns


         endassociate

      end subroutine swirling_new

      !  Save function interface
      subroutine savesol(this, vars,buffer)
         class(SwirlingFlow), intent(inout) :: this
         type(TwoPntBVPDomain), intent(in) :: vars
         real(RK),     intent(in) :: buffer(vars%N())

         integer :: N

         N = vars%N()

         ! RETAIN THE SOLUTION FOR TIME INTEGRATION.
         CALL TWCOPY (N,BUFFER,this%U0)

      end subroutine savesol

      !  Residual function interface:
      subroutine residual(this,error,text,points,time,stride,x,buffer)
         class(SwirlingFlow), intent(inout) :: this
         logical, intent(out) :: error  ! .true. if something went wrong
         integer, intent(in)  :: text   ! output unit (0 for NONE)
         integer, intent(in)  :: points ! Number of 1D variables
         logical, intent(in)  :: time   ! Time-resolved or steady state
         real(RK), intent(in) :: stride !
         real(RK), intent(in) :: x(*)   ! dimensioned >=PMAX, contains the grid
         real(RK), intent(inout) :: buffer(*) ! on input: contains the approximate solution

         CALL TWFUNC(ERROR, TEXT, &
                     BUFFER, this%F, this%F0, this%G, this%G0, this%H, this%K, this%LAMBDA, this%MU, &
                     this%OMEGA, points, this%RHO, STRIDE, this%T, this%T0, TIME, this%TMAX, this%TZERO, &
                     this%U0, this%WMAX, X)

      end subroutine residual

      ! Grid update function interface
      subroutine on_grid_update(this,error,vars,u)
         class(SwirlingFlow), intent(inout) :: this
         logical, intent(out)     :: error
         type(TwoPntBVPDomain), intent(in) :: vars
         real(RK), intent(inout)  :: u(:)

         this%N = vars%N()
         error = this%N<=0

      end subroutine on_grid_update

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

      ! WRITE A SUMMARY.
      subroutine swirl_summary(this,text)
          class(SwirlingFlow), intent(in) :: this
          integer, intent(in) :: text

          character(len=*), parameter :: ID = 'TWMAIN:  '

          WRITE (TEXT, 1) ID, this%U(4, 1), this%OMEGA, this%TZERO, this%TMAX, this%WMAX

          ! INFORMATIVE MESSAGES.
          1 FORMAT(/1X, A9, 1P, E10.2, 0P, '  LAMBDA' &
                  /10X, F10.5, '  OMEGA' &
                  /10X, F10.5, '  TZERO' &
                  /10X, F10.5, '  TMAX'  &
                  /10X, F10.5, '  WMAX')

      end subroutine swirl_summary

end module twopnt_test

