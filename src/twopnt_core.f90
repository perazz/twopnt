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
module twopnt_core
    use iso_fortran_env, only: real64
    implicit none
    private

    public :: twlast,twcopy,twcom,twsolv,twshow,twopnt

    integer, parameter, public :: RK = real64

    ! Numeric constants
    real(RK), parameter, public :: zero       = 0.0_RK
    real(RK), parameter, public :: half       = 0.5_RK
    real(RK), parameter, public :: one        = 1.0_RK
    real(RK), parameter, public :: pi         = acos(-1.0_RK)
    real(RK), parameter, public :: minute     = 60.0_RK
    real(RK), parameter, public :: hour       = 3600.0_RK
    real(RK), parameter, public :: hundred    = 100.0_RK
    real(RK), parameter, public :: smallnum04 = 1.0e-4_RK

    ! Machine epsilon and the absolute and relative perturbations.
    real(RK), parameter, public :: eps   = epsilon(0.0_RK)
    real(RK), parameter, public :: absol = sqrt(eps)
    real(RK), parameter, public :: relat = sqrt(eps)

    ! Error codes
    integer,  parameter, public :: qnull = 0
    integer,  parameter, public :: qbnds = 1
    integer,  parameter, public :: qdvrg = 2

    ! Task codes
    integer,  parameter, public :: qgrid  =  1
    integer,  parameter, public :: qtimst =  2
    integer,  parameter, public :: qsearc =  3
    integer,  parameter, public :: qrefin =  4
    integer,  parameter, public :: qfunct =  5
    integer,  parameter, public :: qjacob =  6
    integer,  parameter, public :: qsolve =  7
    integer,  parameter, public :: qother =  8
    integer,  parameter, public :: qtotal =  9
    integer,  parameter, public :: qentry = 10
    integer,  parameter, public :: qexit  = 11

    character(len=*), parameter, public :: qname(11) = ['    GRID','  EVOLVE','  SEARCH','  REFINE',&
                                                        'FUNCTION','JACOBIAN','   SOLVE','   OTHER',&
                                                        '   TOTAL','   START','    EXIT']

    ! Maximum number of grids attempted
    integer,  parameter, public :: gmax = 100

    logical, parameter :: DEBUG = .true.
    integer, parameter :: CONTRL_MAX_LEN = 40
    integer, parameter :: MAX_ERROR_LINES = 20
    integer, parameter :: MAX_DECAY_ITERATIONS = 5
    integer, parameter :: DEFAULT_MAX_POINTS = 200

    ! Supported versions
    character(len=8), parameter :: vnmbr(*) = [character(len=8) :: '3.18', '3.19', '3.20', '3.21', &
                                                                   '3.22', '3.23', '3.24', '3.25', &
                                                                   '3.26', '3.27', '3.28', '3.29']
    integer, parameter :: vnmbrs = size(vnmbr)

    ! Settings structure (formerly a common block)
    integer, parameter :: cntrls = 22

    ! TWOPNT settings
    type, public :: twcom

        ! Adaptive grid size
        logical  :: adapt  = .false.
        integer  :: leveld = 1
        integer  :: levelm = 1
        logical  :: padd   = .false.
        integer  :: ipadd  = 0
        real(RK) :: ssabs  = 1.0e-9_RK
        integer  :: ssage  = 10
        real(RK) :: ssrel  = 1.0e-6_RK
        logical  :: steady = .true.
        integer  :: steps0 = 0
        integer  :: steps1 = 200
        integer  :: steps2 = 100
        real(RK) :: strid0 = 1.0e-4_RK
        real(RK) :: tdabs  = 1.0e-9_RK
        integer  :: tdage  = 20
        real(RK) :: tdec   = 3.1623_RK
        real(RK) :: tdrel  = 1.0e-6_RK
        real(RK) :: tinc   = 10.0_RK
        real(RK) :: tmax   = 1.0e-2_RK
        real(RK) :: tmin   = 1.0e-20_RK
        real(RK) :: toler0 = 1.0e-9_RK
        real(RK) :: toler1 = 0.2_RK
        real(RK) :: toler2 = 0.2_RK

        contains

           procedure :: init => twinit
           procedure, private :: twsetr
           procedure, private :: twseti
           procedure, private :: twsetl
           generic :: set => twsetr,twseti,twsetl

    end type twcom

    ! TWOPNT work arrays
    type, public :: twwork

        integer, allocatable :: vary (:)   ! (PMAX)
        integer, allocatable :: vary1(:)   ! (PMAX)
        integer, allocatable :: vary2(:)   ! (PMAX)

        real(RK), allocatable :: above(:)  ! (GROUPA + COMPS * PMAX + GROUPB)
        real(RK), allocatable :: below(:)  ! (GROUPA + COMPS * PMAX + GROUPB)
        real(RK), allocatable :: ratio1(:) ! PMAX
        real(RK), allocatable :: ratio2(:) ! PMAX
        real(RK), allocatable :: s0(:)     ! (GROUPA + COMPS * PMAX + GROUPB)
        real(RK), allocatable :: s1(:)     ! (GROUPA + COMPS * PMAX + GROUPB)
        real(RK), allocatable :: usave(:)  ! (GROUPA + COMPS * PMAX + GROUPB)
        real(RK), allocatable :: vsave(:)  ! (GROUPA + COMPS * PMAX + GROUPB)
        real(RK), allocatable :: v1(:)     ! (GROUPA + COMPS * PMAX + GROUPB)
        real(RK), allocatable :: xsave(:)  ! PMAX
        real(RK), allocatable :: y0(:)     ! (GROUPA + COMPS * PMAX + GROUPB)
        real(RK), allocatable :: y1(:)     ! (GROUPA + COMPS * PMAX + GROUPB)

        contains

            procedure :: init => partition_working_space
            procedure :: load_bounds => expand_bounds

    end type twwork

    ! TWOPNT variables
    type, public :: twsize

        ! Group A variables
        integer :: groupa = 0

        ! Unknowns on each grid point
        integer :: comps  = 0

        ! Number of grid points
        integer :: points = 0

        ! Max number of grid points
        integer :: pmax = DEFAULT_MAX_POINTS

        ! Group B variables
        integer :: groupb = 0

        contains

           ! Number of unknowns to be solved on the current grid
           procedure, non_overridable :: N    => twsize_N

           ! Max number of unknowns (storage size)
           procedure, non_overridable :: NMAX => twsize_max

           ! Number of independent variables
           procedure, non_overridable :: NVAR => twsize_NVARS

           ! Check variable sizes
           procedure, non_overridable :: check => twsize_checks

           ! Indices of relevant sizes
           procedure, non_overridable :: idx_B    => twsize_idx_B ! groupb indices
           procedure, non_overridable, private :: twsize_idx_comp_all ! grid component indices
           procedure, non_overridable, private :: twsize_idx_comp_one
           generic :: idx_comp => twsize_idx_comp_all,twsize_idx_comp_one

    end type twsize

    ! TWOPNT problem functions
    type, public :: twfunctions

        ! These functions need to be implemented by the user
        procedure(twopnt_save), nopass, pointer :: save_sol => null()
        procedure(twopnt_residual), nopass, pointer :: fun => null()

        contains

           procedure, nopass :: show => twshow

           ! Provide dense algebra implementation of Jacobian handling
           procedure, pass(this) :: prep => twprep

    end type twfunctions

    ! TWOPNT statistics arrays
    type, public :: twstat

        integer :: grid = 0 ! grid number
        integer :: step = 0 ! step number

        real(RK) :: detail(gmax,qtotal)  = zero
        real(RK) :: timer(qtotal)        = zero
        real(RK) :: total(qtotal)        = zero
        integer  :: event(gmax,qtotal)   = 0
        integer  :: gsize(gmax)          = 0

        contains

           procedure :: new      => init_stats
           procedure :: new_grid => stats_new_grid
           procedure :: tick => tick_stats
           procedure :: tock => tock_stats
           procedure :: print_stats

    end type twstat

    ! TWOPNT state structure
    type, public :: twstate

       ! Present task: ENTRY
       integer :: qtask = qentry

       ! Allow further time evolution
       logical :: allow = .true.

       ! Solution found
       logical :: found = .true.

    end type twstate

    interface realloc
        module procedure realloc_int
        module procedure realloc_real
    end interface realloc

    ! Procedure interfaces
    abstract interface

       ! Residual function calculation:
       subroutine twopnt_residual(error,text,points,time,stride,x,buffer)
          import RK
          implicit none
          logical, intent(out) :: error  ! .true. if something went wrong
          integer, intent(in)  :: text   ! output unit (0 for NONE)
          integer, intent(in)  :: points ! Number of 1D variables
          logical, intent(in)  :: time   ! Time-resolved or steady state
          real(RK), intent(in) :: stride !
          real(RK), intent(in) :: x(*)   ! dimensioned >=PMAX, contains the grid
          real(RK), intent(inout) :: buffer(*) ! on input: contains the approximate solution
                                               ! on output:
       end subroutine twopnt_residual

       ! Pass solution at the beginning of the timestep to the solver
       subroutine twopnt_save(vars,buffer)
          import RK, twsize
          implicit none
          ! Initial solution passed back to the user
          type(twsize), intent(in) :: vars
          real(RK), intent(in) :: buffer(vars%N())
       end subroutine twopnt_save

    end interface


    contains

       ! Indices of groupb variables
       pure function twsize_idx_B(this) result(igroupb)
          class(twsize), intent(in) :: this
          integer, allocatable :: igroupb(:)
          integer :: j

          if (this%groupb<=0) then
             allocate(igroupb(0))
          else
             igroupb = this%groupa + this%comps*this%points +[(j,j=1,this%groupb)]
          end if

       end function twsize_idx_B

       ! Indices of all grid variables of component comp (one per grid point)
       elemental integer function twsize_idx_comp_one(this,comp,point) result(icomp)
          class(twsize), intent(in) :: this
          integer,       intent(in) :: comp,point
          icomp = this%groupa + comp + this%comps*(point - 1)
       end function twsize_idx_comp_one

       ! Indices of all grid variables of component comp (one per grid point)
       pure function twsize_idx_comp_all(this,comp) result(icomp)
          class(twsize), intent(in) :: this
          integer,       intent(in) :: comp
          integer, allocatable :: icomp(:)

          integer :: point

          if (this%comps<=0 .or. comp<=0 .or. comp>this%comps .or. this%points<=0) then
             allocate(icomp(0))
          else
             icomp = this%groupa + comp + [(this%comps*(point - 1), point=1,this%points)]
          end if

       end function twsize_idx_comp_all

       ! Check variable sizes
       subroutine twsize_checks(this,error,id,text)
          class(twsize), intent(in)  :: this
          logical,       intent(out) :: error
          integer,       intent(in)  :: text
          character(*),  intent(in)  :: id

          ! Check no negatives
          error = .not. all([this%comps,this%points,this%groupa,this%groupb]>=0)
          sizes: if (error) then
              if (text>0) write (text, 1) id, this%comps, this%points, this%groupa, this%groupb
              return
          endif sizes

          ! Check that there is at least 1 point if there is 1 component
          error = .not. ((this%comps>0) .eqv. (this%points>0))
          unknowns: if (error) then
              if (text>0) write (text, 2) id, this%comps, this%points
              return
          end if unknowns

          ! Check total size
          error = .not. (this%N() > 0)
          if (error) then
              if (text>0) write (text, 3) id, this%comps, this%points, this%groupa, this%groupb, &
                                               this%N()
              return
          end if

          ! Check storage size
          too_many_points: if (.not.(this%points<=this%pmax)) then
              if (text>0) write (text, 4) id, this%points, this%pmax
              return
          endif too_many_points

          1 format(/1X, a9, 'ERROR.  NUMBERS OF ALL TYPES OF UNKNOWNS MUST BE AT' &
                  /10X, 'LEAST ZERO.' &
                 //10X, i10, '  COMPS, COMPONENTS' &
                  /10X, i10, '  POINTS' &
                  /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                  /10X, i10, '  GROUPB, GROUP B UNKNOWNS')
          2 format(/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE' &
                  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE.' &
                 //10X, i10, '  COMPS, COMPONENTS' &
                  /10X, i10, '  POINTS')
          3 format(/1X, a9, 'ERROR.  TOTAL UNKNOWNS MUST BE POSITIVE.' &
                 //10X, i10, '  COMPS, COMPONENTS' &
                  /10X, i10, '  POINTS' &
                  /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                  /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                  /10X, i10, '  TOTAL NUMBER')
          4 format(/1X, a9, 'ERROR.  THERE ARE TOO MANY POINTS.' &
                 //10X, i10, '  POINTS' &
                  /10X, i10, '  PMAX, LIMIT ON POINTS')

       end subroutine twsize_checks

       ! Return total number of variables
       elemental integer function twsize_NVARS(this) result(N)
          class(twsize), intent(in) :: this
          N = this%GROUPA + this%COMPS + this%GROUPB
       end function twsize_NVARS

       ! Return total number of unknowns
       elemental integer function twsize_N(this) result(N)
          class(twsize), intent(in) :: this
          N = this%GROUPA + this%COMPS * this%POINTS + this%GROUPB
       end function twsize_N

       ! Return max allocation size
       elemental integer function twsize_max(this) result(N)
          class(twsize), intent(in) :: this
          N = this%GROUPA + this%COMPS * this%PMAX + this%GROUPB
       end function twsize_max

       ! Initialize the control structure
       elemental subroutine twinit (this)
          class(twcom), intent(inout) :: this

          character(len=9), parameter :: id = 'TWINIT:  '

          ! SET THE CONTROLS.
          this%adapt  = .false.
          this%leveld = 1
          this%levelm = 1
          this%padd   = .false.
          this%ipadd  = 0
          this%ssabs  = 1.0e-9_RK
          this%ssage  = 10
          this%ssrel  = 1.0e-6_RK
          this%steady = .true.
          this%steps0 = 0
          this%steps1 = 200
          this%steps2 = 100
          this%strid0 = 1.0e-4_RK
          this%tdabs  = 1.0e-9_RK
          this%tdage  = 20
          this%tdec   = 3.1623_RK
          this%tdrel  = 1.0e-6_RK
          this%tinc   = 10.0_RK
          this%tmax   = 1.0e-2_RK
          this%tmin   = 1.0e-20_RK
          this%toler0 = 1.0e-9_RK
          this%toler1 = 0.2_RK
          this%toler2 = 0.2_RK

          return
      end subroutine twinit

      ! Set a control that takes a real value
      subroutine twsetr(this, error, text, contrl, value)
          class(twcom), intent(inout) :: this
          logical,          intent(out) :: error
          integer,          intent(in)  :: text  ! output unit
          character(len=*), intent(in)  :: contrl
          real(RK),         intent(in)  :: value

          ! Local variables
          logical  :: found

          ! SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
          logical, parameter :: mess = .false.

          character(len=*), parameter :: id = 'TWSETR:  '

          ! WRITE ALL MESSAGES.
          print_only: if (text>0 .and. mess) then
               write (text, 1) id
               write (text, 2) id, twtrim(contrl,CONTRL_MAX_LEN)
               write (text, 3) id, twtrim(contrl,CONTRL_MAX_LEN)
               write (text, 5) id, twtrim(contrl,CONTRL_MAX_LEN)
               stop
          end if print_only

          found = .false.

          ! Set the controls.
          select case (contrl)
             case ('SSSABS')
                 found = .true.
                 this%ssabs = value
             case ('SSREL')
                 found = .true.
                 this%ssrel = value
             case ('STRID0')
                 found = .true.
                 this%strid0 = value
             case ('TDABS')
                 found = .true.
                 this%tdabs = value
             case ('TDEC')
                 found = .true.
                 this%tdec = value
             case ('TDREL')
                 found = .true.
                 this%tdrel = value
             case ('TINC')
                 found = .true.
                 this%tinc = value
             case ('TMAX')
                 found = .true.
                 this%tmax = value
             case ('TMIN')
                 found = .true.
                 this%tmin = value
             case ('TOLER0')
                 found = .true.
                 this%toler0 = value
             case ('TOLER1')
                 found = .true.
                 this%toler1 = value
             case ('TOLER2')
                 found = .true.
                 this%toler2 = value
             case ('ADAPT','STEADY')
                 error = .true.
                 if (text>0) write (text, 2) id, twtrim(contrl,CONTRL_MAX_LEN)
                 return
             case ('LEVELD','LEVELM','PADD','SSAGE','STEPS0','STEPS1','STEPS2','TDAGE')
                 error = .true.
                 if (text>0) write (text, 3) id, twtrim(contrl,CONTRL_MAX_LEN)
                 return
          end select

          error = .not. found
          if (error .and. text>0) write (text, 5) id, twtrim(contrl,CONTRL_MAX_LEN)

          return

          ! Error messages.
          1 format(/1X, a9, 'ERROR.  TWINIT FAILS.')
          2 format(/1X, a9, 'ERROR.  THE CONTROL TAKES A LOGICAL VALUE WHICH' &
                  /10X, 'MUST BE SET USING TWSETL.' &
                 //10X, '     CONTROL:  ', a)
          3 format(/1X, a9, 'ERROR.  THE CONTROL TAKES AN INTEGER VALUE WHICH' &
                  /10X, 'MUST BE SET USING TWSETI.' &
                 //10X, '     CONTROL:  ', a)
          5 format(/1X, a9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.' &
                 //10X, '     CONTROL:  ', a)

      end subroutine twsetr

      ! Set a control that takes an integer value,
      subroutine twseti(this, error, text, contrl, value)
          class(twcom), intent(inout) :: this
          logical,          intent(out) :: error
          integer,          intent(in)  :: text  ! output unit
          character(len=*), intent(in)  :: contrl
          integer,          intent(in)  :: value

          ! Local variables
          logical  :: found

          ! SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
          logical, parameter :: mess = .false.

          character(len=*), parameter :: id = 'TWSETI:  '

          ! WRITE ALL MESSAGES.
          print_only: if (text>0 .and. mess) then
               write (text, 1) id
               write (text, 2) id, twtrim(contrl,CONTRL_MAX_LEN)
               write (text, 3) id, twtrim(contrl,CONTRL_MAX_LEN)
               write (text, 5) id, twtrim(contrl,CONTRL_MAX_LEN)
               stop
          end if print_only

          found = .false.

          ! Set the controls.
          select case (contrl)
             case ('LEVELD')
                 found = .true.
                 this%leveld = value
             case ('LEVELM')
                 found = .true.
                 this%levelm = value
             case ('PADD')
                 found = .true.
                 this%padd = .true.
                 this%ipadd = value
             case ('SSAGE')
                 found = .true.
                 this%ssage = value
             case('STEPS0')
                 found = .true.
                 this%steps0 = value
             case('STEPS1')
                 found = .true.
                 this%steps1 = value
             case('STEPS2')
                 found = .true.
                 this%steps2 = value
             case ('TDAGE')
                 found = .true.
                 this%tdage = value
             case ('ADAPT','STEADY')
                 error = .true.
                 if (text>0) write (text, 2) id, twtrim(contrl,CONTRL_MAX_LEN)
                 return
             case ('SSABS','SSREL','STRID0','TDABS','TDEC','TDREL','TINC','TMAX','TMIN',&
                   'TOLER0','TOLER1','TOLER2')
                 error = .true.
                 if (text>0) write (text, 3) id, twtrim(contrl,CONTRL_MAX_LEN)
                 return
          end select

          error = .not. found
          if (error .and. text>0) write (text, 5) id, twtrim(contrl,CONTRL_MAX_LEN)

          return

          ! Error messages.
          1 format(/1X, a9, 'ERROR.  TWINIT FAILS.')
          2 format(/1X, a9, 'ERROR.  THE CONTROL TAKES A LOGICAL VALUE WHICH' &
                  /10X, 'MUST BE SET USING TWSETL.' &
                 //10X, '     CONTROL:  ', a)
          3 format(/1X, a9, 'ERROR.  THE CONTROL TAKES A REAL VALUE WHICH MUST BE' &
                  /10X, 'SET USING TWSETR.' &
                 //10X, '     CONTROL:  ', a)
          5 format(/1X, a9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.' &
                 //10X, '     CONTROL:  ', a)

      end subroutine twseti

      ! Set a control that takes a logical value
      subroutine twsetl(this, error, text, contrl, value)
          class(twcom), intent(inout) :: this
          logical,          intent(out) :: error
          integer,          intent(in)  :: text  ! output unit
          character(len=*), intent(in)  :: contrl
          logical,          intent(in)  :: value


          ! Local variables
          logical  :: found

          ! SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
          logical, parameter :: mess = .false.

          character(len=*), parameter :: id = 'TWSETL:  '

          ! WRITE ALL MESSAGES.
          print_only: if (text>0 .and. mess) then
               write (text, 1) id
               write (text, 2) id, twtrim(contrl,CONTRL_MAX_LEN)
               write (text, 3) id, twtrim(contrl,CONTRL_MAX_LEN)
               write (text, 5) id, twtrim(contrl,CONTRL_MAX_LEN)
               stop
          end if print_only

          found = .false.

          ! Set the controls.
          select case (contrl)
             case ('ADAPT')
                 found = .true.
                 this%adapt = value
             case ('STEADY')
                 found = .true.
                 this%steady = value
             case ('LEVELD','LEVELM','PADD','SSAGE','STEPS0','STEPS1','STEPS2','TDAGE')
                 error = .true.
                 if (text>0) write (text, 2) id, twtrim(contrl,CONTRL_MAX_LEN)
                 return
             case ('SSABS','SSREL','STRID0','TDABS','TDEC','TDREL','TINC','TMAX','TMIN',&
                   'TOLER0','TOLER1','TOLER2')
                 error = .true.
                 if (text>0) write (text, 3) id, twtrim(contrl,CONTRL_MAX_LEN)
                 return
          end select

          error = .not. found
          if (error .and. text>0) write (text, 5) id, twtrim(contrl,CONTRL_MAX_LEN)

          return

          ! Error messages.
          1 format(/1X, a9, 'ERROR.  TWINIT FAILS.')
          2 format(/1X, a9, 'ERROR.  THE CONTROL TAKES AN INTEGER VALUE WHICH' &
                  /10X, 'MUST BE SET USING TWSETI.' &
                 //10X, '     CONTROL:  ', a)
          3 format(/1X, a9, 'ERROR.  THE CONTROL TAKES A REAL VALUE WHICH MUST BE' &
                  /10X, 'SET USING TWSETR.' &
                 //10X, '     CONTROL:  ', a)
          5 format(/1X, a9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.' &
                 //10X, '     CONTROL:  ', a)

      end subroutine twsetl

      ! *******************************************************************************************************
      ! STRINGS
      ! *******************************************************************************************************

      ! FIND THE LAST NONBLANK CHARACTER IN A STRING.
      pure subroutine twlast(length, string)
         character(len=*), intent(in) :: string
         integer, intent(out) :: length

         intrinsic :: len_trim

         length = len_trim(string)
         if (length==0) length = 1
         return
      end subroutine twlast

      !> Trim a string to max length, add trailing
      pure function twtrim(string,maxlen) result(trimmed)
          character(len=*), intent(in) :: string
          integer, intent(in) :: maxlen
          character(len=:), allocatable :: trimmed

          integer :: length

          call twlast(length,string)
          if (length<=maxlen) then
              trimmed = string
          else
              length = maxlen
              trimmed = string(1:maxlen-3)//'...'
          end if
      end function twtrim

      ! Write a common logarithm to a character string
      subroutine twlogr (string, value)
          character(len=*), intent(inout) :: string
          real(RK), intent(in) :: value

          intrinsic :: len, log10

          if (len(string)>=6) then
              if (value<zero) then
                  string = ' '
              elseif (value==zero) then
                  string = '  ZERO'
              else
                  write(string,'(F6.2)') log10(value)
              end if
          else
              string = '******'
          end if
      end subroutine twlogr

      ! SQUEEZE LEADING BLANKS AND MULTIPLE BLANKS FROM A CHARACTER
      ! STRING.  RETURN THE LENGTH OF THE SQUEEZED STRING.
      subroutine twsqez (length, string)
          integer, intent(out) :: length
          character(len=*), intent(inout) :: string
          !implicit complex (a - z)

          character :: char
          integer :: j
          intrinsic :: len
          logical :: blank

          ! SQUEEZE THE STRING.
          length = 0
          blank = .true.
          do j = 1, len (string)
             char = string (j : j)
             if (.not. blank .or. char .ne. ' ') then
                blank = char .eq. ' '
                length = length + 1
                string (length : length) = char
             end if
          end do

          ! ADJUST THE LENGTH AND PAD THE STRING.
          if (length>0) then
             if (string (length : length) .eq. ' ') length = length - 1
             if (length < len (string)) string (length + 1 : ) = ' '
          else
             length = 1
          end if

          return
      end subroutine twsqez

      ! *******************************************************************************************************
      ! MATH
      ! *******************************************************************************************************

      subroutine twshow(error, text, buffer, vars, grid, x)
          logical,      intent(out) :: error
          integer,      intent(in)  :: text
          type(twsize), intent(in) :: vars
          logical,      intent(in)  :: grid
          real(RK),     intent(in) :: buffer(vars%N()), x(*)

          ! Local variables
          character(len=80) :: string, title(6)
          integer :: cols, comp, count, first, groups, j, last, length, point
          intrinsic :: min, merge
          character(*), parameter :: id = 'TWSHOW:  '

          ! Check the arguments.
          call vars%check(error,id,text)
          if (error) return

          !///  COUNT THE GROUPS.
          groups = 0
          if (0 < vars%groupa) groups = groups + 1
          if (0 < vars%groupb) groups = groups + 1
          if (vars%comps>0 .and. vars%points>0) groups = groups + 1

          !///  CHOOSE NUMBER OF DATA COLUMNS.
          cols = merge(5,6,grid)

          !///////////////////////////////////////////////////////////////////////
          !
          !
          !
          !///////////////////////////////////////////////////////////////////////
          print_data: if (text>0) then

              ! (2) PRINT THE GROUPED DATA.
              if (0 < vars%groupa) then
                 if (1 < groups) write (text, 11) 'GROUP A UNKNOWNS'
                 write (text, 12) (j, buffer(j), j = 1, vars%groupa)
              end if

              if (0 < vars%groupb) then
                 if (1 < groups) write (text, 11) 'GROUP B UNKNOWNS'
                 write (text, 12) &
                    (j, buffer(vars%groupa + vars%comps*vars%points + j), j = 1, vars%groupb)
              end if

              ! (2) PRINT THE COMPONENTS AT POINTS.
              if (vars%comps>0 .and. vars%points>0) then

                 if (1 < groups) write (text, 11) 'COMPONENTS AT POINTS'

                 components: do first = 1, vars%comps, cols
                    count = 0
                    last = min (first + cols - 1, vars%comps)
                    do comp = first, last
                        count = count + 1
                        title(count) = ' '
                        write (string, '(A5, I5)') 'COMP ', comp
                        call twsqez (length, string)
                        title(count) (11 - length : 10) = string
                    end do

                    if (grid) then
                       write (text, 13) 'GRID POINT', (title(j), j = 1, count)
                    else
                       write (text, 13) (title(j), j = 1, count)
                    end if

                    if (count == cols) then
                       if (grid) then
                          write (text, 14) (point, x(point), (buffer(vars%idx_comp(comp,point)), &
                             comp = first, last), point = 1, vars%points)
                       else
                          write (text, 15) (point, (buffer(vars%idx_comp(comp,point)), &
                             comp = first, last), point = 1, vars%points)
                       end if
                    else
                       do point = 1, vars%points
                          if (grid) then
                             write (text, 14) point, x(point), (buffer(vars%idx_comp(comp,point)), &
                                comp = first, last)
                          else
                             write (text, 15) point, (buffer(vars%idx_comp(comp,point)), comp = first, last)
                          end if
                       end do
                    end if
                 end do components
              end if

          end if print_data


          ! INFORMATIVE MESSAGES.
          11 format(/10X, a)
          12 format(/(10X, 4(i3, '> ', 1pe10.3)))
          13 format(/14X, 6(1X, a10))
          14 format(10X, 0p, i3, '>', f11.6, 1p, 5E11.3)
          15 format(10X, 0p, i3, '>', 1p, 6E11.3)

      end subroutine twshow

      ! COPY ONE VECTOR TO ANOTHER.
      pure subroutine twcopy (n, from, to)
          integer , intent(in)  :: n
          real(RK), intent(in)  :: from(*)
          real(RK), intent(out) ::  to(*)
          if (n>0) to(1:n) = from(1:n)
          return
      end subroutine twcopy

      ! Compute the max-norm of a vector
      pure real(RK) function twnorm (n,x)
          integer,  intent(in) :: n
          real(RK), intent(in) :: x(n)
          intrinsic :: abs,maxval
          twnorm = zero
          if (n>0) twnorm = maxval(abs(x),1)
      end function twnorm

      ! SOLVE A SYSTEM OF LINEAR EQUATIONS USING THE MATRIX PREPARED BY TWPREP.
      subroutine twsolv(error, text, a, asize, buffer, vars, pivot)

          integer     , intent(in) :: asize,text
          type(twsize), intent(in) :: vars
          integer     , intent(in) :: pivot(vars%N())
          real(RK)    , intent(in) :: a(asize)
          real(RK)    , intent(inout) :: buffer(vars%N())
          logical     , intent(out) :: error

          ! Local variables
          integer   :: n,width
          intrinsic :: max
          character(*), parameter :: id = 'TWSOLV:  '

          ! SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
          logical, parameter :: mess = .false.

          !***** (1) PROLOGUE *****

          ! WRITE MESSAGES only.
          if (mess .and. text>0) then
              write (text, 1) id,vars%comps,vars%points,vars%groupa,vars%groupb,n
              write (text, 2) id,vars%comps,vars%points,vars%groupa,vars%groupb,n,width,(3*width+2)*n,asize
              stop
          end if

          ! CHECK THE ARGUMENTS.
          n = vars%N()
          call vars%check(error,id,text)
          if (error) return

          width = vars%comps+max(vars%comps,vars%groupa,vars%groupb) - 1
          error = .not. ((3 * width + 2) * n <= asize)
          if (error) then
              if (text>0) write (text,2) id, vars%comps,vars%points,vars%groupa,vars%groupb, n, width, &
                                         (3*width+2)*n, asize
              return
          end if

          !***** (2) SCALE AND SOLVE THE EQUATIONS. *****
          buffer(1:n) = buffer(1:n) * a(1:n)
          call twgbsl(a(n + 1), 3 * width + 1, n, width, width, pivot, buffer)

          return

          !Error messages.
          1 format (/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE' &
                   /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES' &
                   /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS' &
                   /10X, 'MUST BE POSITIVE.' &
                  //10X, i10, '  COMPS, COMPONENTS' &
                   /10X, i10, '  POINTS' &
                   /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                   /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                   /10X, i10, '  TOTAL UNKNOWNS')

          2 format (/1X, a9, 'ERROR.  THE MATRIX SPACE IS TOO SMALL.' &
                  //10X, i10, '  COMPS, COMPONENTS' &
                   /10X, i10, '  POINTS' &
                   /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                   /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                   /10X, i10, '  MATRIX ORDER' &
                   /10X, i10, '  STRICT HALF BANDWIDTH' &
                  //10X, i10, '  SPACE EXPECTED' &
                   /10X, i10, '  ASIZE, PROVIDED')

          return
      end subroutine twsolv


      ! SOLVE A SYSTEM OF LINEAR EQUATIONS FOR TWSOLV.  BASED ON _GBSL
      ! FROM THE LINPACK LIBRARY.
      subroutine twgbsl (abd, lda, n, lower, upper, pivot, b)
          real(RK), intent(in) :: abd(lda,*)
          real(RK), intent(inout) :: b(*)
          integer , intent(in) :: pivot(*)

          real(RK)  :: t
          integer   :: j, jdiag, k, l, la, lb, lda, lm, lower, n, upper
          intrinsic :: min

          jdiag = upper + lower + 1

          if (0 < lower) then
             do k = 1, n - 1
                  l = pivot(k)
                  t = b(l)
                  if (l /= k) then
                     b(l) = b(k)
                     b(k) = t
                  end if

                  lm = min (lower, n - k)
                  do j = 1, lm
                     b(k+j) = t * abd(jdiag + j, k) + b(k + j)
                  end do
              end do
          end if

          do k = n, 1, - 1
             b(k) = b(k) / abd (jdiag, k)
             lm = min (k, jdiag) - 1
             la = jdiag - lm - 1
             lb = k - lm - 1
             t = - b(k)
             do j = 1, lm
                 b(lb + j) = t * abd(la + j, k) + b(lb + j)
             end do
          end do

      return
      end subroutine twgbsl

      ! Factor a banded matrix and estimate the reciprocal of its condition number.
      ! Based on _GBCO from the LINPACK library.
      subroutine twgbco (a, lda, n, lower, upper, pivot, rcond, z)
          integer, intent(in)     :: lda, n, lower, upper
          integer, intent(inout)  :: pivot(n)
          real(RK), intent(inout) :: a(lda,n),z(n)
          real(RK), intent(out)   :: rcond

          real(RK)  :: dsum, anorm, ek, s, sm, t, wk, wkm, ynorm
          integer   :: first, info, j, jdiag, ju, k, last, mm
          intrinsic :: abs, max, min, sign, sum

          jdiag = lower + upper + 1

          ! Compute the 1-norm of A
          anorm = zero
          do k = 1, n
             first = max(lower + 1, jdiag + 1 - k)
             last  = min(jdiag + lower, jdiag + n - k)
             anorm = max(anorm, sum(abs(a(first:last,k))))
          end do

          ! Factor A
          call twgbfa (a, lda, n, lower, upper, pivot, info)

          ! Solve transpose(U) * W = E
          ek = one
          z  = zero

          ju = 0
          solve_1: do k = 1, n
             if (z(k) /= zero) ek = sign(ek, -z(k))

             if (abs (ek - z(k)) > abs (a(jdiag, k))) then
                s  = abs (a(jdiag, k)) / abs (ek - z(k))
                ek = s * ek
                z  = s*z
             end if

             wk  = ek - z(k)
             wkm = - ek - z(k)
             s   = abs (wk)
             sm  = abs (wkm)
             if (a(jdiag, k) /= zero) then
                wk  = wk / a(jdiag, k)
                wkm = wkm / a(jdiag, k)
             else
                wk  = one
                wkm = one
             end if

             ju = min (max (ju, upper + pivot(k)), n)
             mm = jdiag
             if (k + 1 <= ju) then
                do j = k + 1, ju
                   mm = mm - 1
                   sm = sm + abs (z(j) + wkm * a(mm, j))
                   z(j) = z(j) + wk * a(mm, j)
                   s = s + abs (z(j))
                end do

                if (s < sm) then
                   t = wkm - wk
                   wk = wkm
                   mm = jdiag
                   do j = k + 1, ju
                      mm = mm - 1
                      z(j) = z(j) + t * a(mm, j)
                   end do
                end if
             end if

             z(k) = wk

          end do solve_1

          s = one/sum(abs(z))
          z = s*z

          ! Solve transpose(L) * Y = W
          solve_2: do k = n, 1, - 1
             dsum = zero
             do j = 1, min (lower, n - k)
                dsum = dsum + dble (a(jdiag + j, k)) * dble (z(k + j))
             end do
             z(k) = z(k) + dsum

             if (one < abs (z(k))) then
                s = one / abs (z(k))
                z = s*z
             end if

             j    = pivot(k)
             t    = z(j)
             z(j) = z(k)
             z(k) = t
          end do solve_2

          s = one / sum(abs(z))
          z = s*z
          ynorm = one

          ! Solve L * V = Y
          solve_3: do k = 1, n
             j = pivot(k)
             t = z(j)
             z(j) = z(k)
             z(k) = t

             do j = 1, min (lower, n - k)
                z(k + j) = t * a(jdiag + j, k) + z(k + j)
             end do

             if (one < abs (z(k))) then
                s = one / abs (z(k))
                z = s*z
                ynorm = s * ynorm
             end if

          end do solve_3

          s = one / sum(abs(z))
          z = s*z
          ynorm = s * ynorm

          ! Solve U * Z = W
          solve_4: do k = n, 1, - 1
             if (abs (z(k)) > abs (a(jdiag, k))) then
                s = abs (a(jdiag, k)) / abs (z(k))
                z = s*z
                ynorm = s*ynorm
             end if

             if (a(jdiag, k) /= zero) then
                z(k) = z(k) / a(jdiag, k)
             else
                z(k) = one
             end if

             t = - z(k)
             do j = 1, min (k, jdiag) - 1
                z(k - j) = t * a(jdiag - j, k) + z(k - j)
             end do
          end do solve_4

          s = one / sum(abs(z))
          z = s*z
          ynorm = s * ynorm

          ! FORM RCOND
          rcond = merge(ynorm/anorm,zero,anorm/=zero)

          return
      end subroutine twgbco

      ! Factor a banded matrix for twgbco. Based on _GBFA from the LINPACK library.
      subroutine twgbfa (a, lda, n, lower, upper, pivot, info)
          real(RK), intent(inout) :: a(lda, n)
          integer , intent(inout) :: pivot(n)
          integer , intent(in)    :: lda, n, lower, upper
          integer , intent(out)   :: info

          ! Local variables
          real(RK)  :: t, value
          integer   :: i, j, jk, k, pdiag, pjk, save
          intrinsic :: abs, max, min

          ! Initialize packed row position of the diagonal
          pdiag = lower + upper + 1

          ! Info stores the value of a zero diagonal
          info = 0

          ! Loop over columns
          columns: do k = 1, n

             ! Initialize fill-in space
             a(:lower, k) = zero

             ! Loop over the previous columns
             previous: do j = max (1, k - lower - upper), k - 1
                pjk = packed(pivot(j), k)
                jk  = packed(j, k)
                t   = a(pjk, k)
                if (pjk /= jk) then
                   a(pjk, k) = a(jk, k)
                   a(jk, k) = t
                end if

                if (t /= zero) &
                forall (i = 1: min (lower, n - j)) a(jk + i, k) = t * a(pdiag + i, j) + a(jk + i, k)

             end do previous

             ! Find the pivot
             save = pdiag
             value = abs (a(pdiag, k))
             do i = pdiag + 1, pdiag + min (lower, n - k)
                if (value < abs (a(i, k))) then
                   save = i
                   value = abs (a(i, k))
                end if
             end do
             pivot(k) = true (save, k)

             ! Interchange if necessary
             if (save /= pdiag) then
                t = a(save, k)
                a(save, k) = a(pdiag, k)
                a(pdiag, k) = t
             end if

             ! Scale the lower column
             if (a(save, k) /= zero) then
                t = -one/ a(pdiag, k)
                forall (i = pdiag + 1: pdiag + min (lower, n - k)) a(i, k) = t * a(i, k)
             else
                info = k
             end if

          end do columns

          ! The final column is trivial
          pivot(n) = n
          if (a(pdiag, n) == zero) info = n

          return

          contains

              ! Statement functions

              ! Packed row position of entry J in column K
              elemental integer function packed(j, k)
                 integer, intent(in) :: j, k
                 packed = j - k + pdiag
              end function packed

              ! True row position of entry J packed in column K
              elemental integer function true(j, k)
                 integer, intent(in) :: j, k
                 true = j - pdiag + k
              end function true

      end subroutine twgbfa

      ! Evaluate a block tridiagonal Jacobian matrix by one-sided finite differences and reverse
      ! communication, pack the matrix into the LINPACK banded form, scale the rows, and factor the
      ! matrix using LINPACK's SGBCO
      subroutine twprep(this, error, text, a, asize, buffer, vars, condit, pivot, time, stride, x)
          class(twfunctions), intent(inout) :: this
          logical,      intent(out)   :: error
          integer,      intent(in)    :: text ! output unit
          integer,      intent(in)    :: asize
          type(twsize), intent(in)    :: vars
          integer,      intent(inout) :: pivot (vars%N())
          real(RK),     intent(inout) :: buffer(vars%N())
          real(RK),     intent(inout) :: a(asize), condit
          logical ,     intent(in)    :: time
          real(RK),     intent(in)    :: stride
          real(RK),     intent(in)    :: x(:) ! mesh


          real(RK) :: delta, temp
          integer :: block, blocks, cfirst, clast, col, count, diag, j, lda, length, n, offset, &
                     rfirst, rlast, row, skip, width
          intrinsic :: abs, int, max, min, mod, sqrt
          logical :: found

          ! Parameters
          character(len=*), parameter :: id = 'TWPREP:  '

          ! Check that the residual function is present.
          if (.not.associated(this%fun)) then
             error = .true.
             if (text>0) write (text, 1) id
             return
          end if

          ! CHECK THE ARGUMENTS.
          n = vars%N()
          call vars%check(error,id,text)
          if (error) return

          width = vars%comps + max(vars%comps,vars%groupa,vars%groupb) - 1
          error = .not. ((3 * width + 2) * n <= asize)
          if (error) then
             if (text>0) write (text, 2) id, vars%comps,vars%points,vars%groupa,vars%groupb,n,width,(3*width+2) * n, asize
             return
          end if

          ! Initialize counters and pointers

          ! Main diagonal row in the packing which places diagonals in rows
          diag = 2 * width + 1

          ! Packed row dimension
          lda  = 3 * width + 1
          skip = 2 * width + 1

          ! BLOCKS AND BLOCK SIZES
          ! Temporarily store block sizes and pointers in array "pivot"
          blocks = 0
          if (0 < vars%groupa) then
             blocks = blocks + 1
             pivot(blocks) = vars%groupa
          end if

          do j = 1, vars%points
             blocks = blocks + 1
             pivot(blocks) = vars%comps
          end do

          if (0 < vars%groupb) then
             blocks = blocks + 1
             pivot(blocks) = vars%groupb
          end if

          ! ***** (2) INITIALIZE THE COLUMNS OF THE MATRIX *****

          ! Store evaluation vector
          a(:n) = buffer(:n)

          ! Clear matrix
          a(n+1:n+lda*n) = zero

          ! EVALUATE THE FUNCTION AT THE UNPERTURBED X.
          call this%fun(error,text,vars%points,time,stride,x,buffer)
          if (error) return

          ! Place function values into the matrix.
          clast = 0
          do block = 1, blocks
             cfirst = clast + 1
             clast  = clast + pivot(block)

             if (block > 1) then
                rfirst = cfirst - pivot(block - 1)
             else
                rfirst = cfirst
             end if

             if (block < blocks) then
                rlast = clast + pivot(block + 1)
             else
                rlast = clast
             end if

             do col = cfirst, clast
                offset = n + diag - col + lda * (col - 1)
                do row = rfirst, rlast
                   a(offset + row) = buffer(row)
                end do
             end do
          end do

          ! ***** (3) Form the columns of the matrix. *****

          column_groups: do

              found = .false.

              ! Restore the evaluation vector
              buffer(:n) = a(:n)

              ! Perturb vector at independent positions.
              block = 1
              cfirst = 1
              perturb_vector: do while (block<=blocks)
                   if (0 < pivot(block)) then
                       found = .true.
                       col = cfirst - 1 + pivot(block)
                       delta = relat * a(col) + sign(absol,a(col))
                       buffer(col) = buffer(col) + delta
                       count = 3
                   else
                       count = 1
                   end if

                   do j = 1, count
                      if (block == 1 .and. 0 < vars%groupa) then
                         cfirst = cfirst + vars%groupa
                      else if (block == blocks .and. 0 < vars%groupb) then
                         cfirst = cfirst + vars%groupb
                      else
                         cfirst = cfirst + vars%comps
                      end if
                      block = block + 1
                   end do
              end do perturb_vector

              if (.not. found) exit column_groups

              ! EVALUATE THE FUNCTION AT THE PERTURBED VALUES.
              call this%fun(error,text,vars%points,time,stride,x,buffer)
              if (error) return

              ! DIFFERENCE TO FORM THE COLUMNS OF THE JACOBIAN MATRIX.
              block  = 1
              cfirst = 1
              form_columns: do while (block<=blocks)
                 if (0 < pivot(block)) then
                    col = cfirst - 1 + pivot(block)
                    pivot(block) = pivot(block) - 1

                    delta  = relat * a(col) + sign(absol,a(col))
                    temp   = one / delta
                    offset = n + diag - col + lda * (col - 1)

                    if (block == 1 .and. 0 < vars%groupa) then
                       clast = cfirst + vars%groupa - 1
                    else if (block == blocks .and. 0 < vars%groupb) then
                       clast = cfirst + vars%groupb - 1
                    else
                       clast = cfirst + vars%comps - 1
                    end if

                    if (1 < block) then
                       if (block == 2 .and. 0 < vars%groupa) then
                          rfirst = cfirst - vars%groupa
                       else
                          rfirst = cfirst - vars%comps
                       end if
                    else
                       rfirst = cfirst
                    end if

                    if (block < blocks) then
                       if (block == blocks - 1 .and. 0 < vars%groupb) then
                          rlast = clast + vars%groupb
                       else
                          rlast = clast + vars%comps
                       end if
                    else
                       rlast = clast
                    end if

                    forall(row=rfirst:rlast) a(offset+row) = (buffer(row) - a(offset+row))*temp

                    count = 3
                 else
                    count = 1
                 end if

                 do j = 1, count
                    if (block == 1 .and. 0 < vars%groupa) then
                       cfirst = cfirst + vars%groupa
                    else if (block == blocks .and. 0 < vars%groupb) then
                       cfirst = cfirst + vars%groupb
                    else
                       cfirst = cfirst + vars%comps
                    end if
                    block = block + 1
                 end do

              end do form_columns

          end do column_groups

          ! ***** (4) CHECK FOR ZERO COLUMNS. *****
          call count_zero_columns(a,n,diag,lda,width,count)
          error = .not. (count == 0)
          if (error) then
              call print_invalid_rowscols(id,text,vars,a,count,.false.)
              return
          endif

          ! ***** (5) SCALE THE ROWS. *****
          call scale_rows(a,n,diag,lda,width,count)
          error = .not. (count == 0)
          if (error) then
             call print_invalid_rowscols(id,text,vars,a,count,.true.)
             return
          end if

          ! ***** (6) FACTOR THE MATRIX.
          call twgbco(a(n+1), lda, n, width, width, pivot, condit, buffer)
          error = condit == zero
          if (error) then
              if (text>0) write (text, 3) id
              return
          end if
          condit = one/condit

          return

          ! Formats section
          1 format(/1X, a9, 'ERROR.  THE PROBLEM FUNCTION IS UNDEFINED.')
          2 format(/1X, a9, 'ERROR.  THE MATRIX SPACE IS TOO SMALL.' &
                 //10X, i10, '  COMPS, COMPONENTS' &
                  /10X, i10, '  POINTS' &
                  /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                  /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                  /10X, i10, '  MATRIX ORDER' &
                  /10X, i10, '  STRICT HALF BANDWIDTH' &
                 //10X, i10, '  SPACE REQUIRED' &
                  /10X, i10, '  ASIZE, PROVIDED')
          3 format(/1X, a9, 'ERROR.  THE JACOBIAN MATRIX IS SINGULAR.')

      end subroutine twprep

      !> Sum columns, put result in a(1:n), return count of zero columns
      pure subroutine count_zero_columns(a,n,diag,lda,width,count)
          real(RK), intent(inout) :: a(*)
          integer , intent(in)    :: n,diag,lda,width
          integer, intent(out)    :: count

          integer :: col,row,offset
          real(RK) :: col_sum

          count = 0
          do col = 1, n
             offset = n + diag - col + lda * (col - 1)
             col_sum = zero
             do row = max (col - width, 1), min (col + width, n)
                col_sum = col_sum + abs (a(offset + row))
             end do
             a(col) = col_sum
             if (col_sum == zero) count = count + 1
          end do

      end subroutine count_zero_columns

      !>
      pure subroutine scale_rows(a,n,diag,lda,width,count)
          real(RK), intent(inout) :: a(*)
          integer , intent(in)    :: n,diag,lda,width
          integer, intent(out)    :: count

          integer :: col,row,offset
          real(RK) :: col_sum,temp
          intrinsic :: min, abs, max

          count = 0
          rows: do row = 1, n
             offset = n + diag + row
             col_sum = zero
             do col = max (row - width, 1), min (row + width, n)
                col_sum = col_sum + abs (a(offset - col + lda * (col - 1)))
             end do

             if (col_sum == zero) then
                count = count + 1
                a(row) = col_sum
             else
                temp = one / col_sum
                a(row) = temp

                do col = max (row - width, 1), min (row + width, n)
                     a(offset - col + lda * (col - 1)) = &
                     a(offset - col + lda * (col - 1)) * temp
                end do
             endif
          end do rows

      end subroutine scale_rows

      subroutine print_invalid_rowscols(id,text,vars,a,invalid,rows)
         character(*), intent(in) :: id
         integer, intent(in) :: text
         type(twsize), intent(in) :: vars
         real(RK), intent(in) :: a(*)
         integer,  intent(in) :: invalid
         logical, intent(in) :: rows

         integer :: count,j,length
         character(len=80) :: string

         if (text==0) return

         if (rows) then
            write (text, 4) id, vars%comps, vars%points, vars%groupa, vars%groupb, vars%N(), invalid
         else
            write (text, 3) id, vars%comps, vars%points, vars%groupa, vars%groupb, vars%N(), invalid
         endif

         count = 0
         loop_matrix: do j = 1, vars%N()
            if (a(j) == zero) then
               count = count + 1
               if (count <= MAX_ERROR_LINES) then
                  if (j <= vars%groupa) then
                     write (string, '(A, I10)') 'GROUP A ', j

                  else if (j <= vars%groupa + vars%comps * vars%points) then
                     write (string, '(A, I10, A, I10)') &
                        ' COMPONENT ', mod (j - vars%groupa - 1, vars%comps) + 1, &
                        ' AT POINT ', int ((j - vars%groupa - 1) / vars%comps) + 1

                  else
                     write (string, '(A, I10)') &
                        'GROUP B ', j - vars%groupa - vars%comps * vars%points
                  end if
                  call twsqez (length, string)
                  write (text, 1) string (1 : length)
               end if
            end if
         end do loop_matrix
         if (MAX_ERROR_LINES < count) write (text, 2)

         1 format(10X, a)
         2 format(10X, '... MORE')
         3 format(/1X, a9, 'ERROR.  SOME COLUMNS ARE ZERO.' &
                //10X, i10, '  COMPS, COMPONENTS' &
                 /10X, i10, '  POINTS' &
                 /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                 /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                 /10X, i10, '  TOTAL COLUMNS' &
                 /10X, i10, '  ZERO COLUMNS' &
                //10X, 'UNKNOWNS WITH ZERO COLUMNS:'/)

         4 format(/1X, a9, 'ERROR.  SOME ROWS ARE ZERO.' &
                //10X, i10, '  COMPS, COMPONENTS' &
                 /10X, i10, '  POINTS' &
                 /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                 /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                 /10X, i10, '  TOTAL ROWS' &
                 /10X, i10, '  ZERO ROWS' &
                //10X, 'ZERO ROWS:'/)

      end subroutine print_invalid_rowscols

      ! *******************************************************************************************************
      ! UTILITIES
      ! *******************************************************************************************************

      ! Obtain computing time in seconds.
      subroutine twtime(timer)
         real(RK), intent(out) :: timer
         real :: temp

         call cpu_time(temp)
         timer = real(temp,RK)

      end subroutine twtime

      ! Obtain elapset computing time in seconds.
      subroutine twlaps (timer)
         real(RK), intent(inout) :: timer

         real(RK) :: temp

         call twtime(temp)
         timer = temp - timer

      end subroutine twlaps

      ! Perform time evolution
      subroutine evolve(error, text, above, below, buffer, vars, condit, desire, &
                        leveld, levelm, name, names, report, s0, s1, signal, &
                        step, steps2, strid0, stride, success, tdabs, tdage, tdec, &
                        tdrel, time, tinc, tmax, tmin, v0, v1, vsave, y0, y1, ynorm, x, functions)

      integer,          intent(in)    :: text
      type(twsize),     intent(in)    :: vars
      logical,          intent(out)   :: error
      logical,          intent(out)   :: time
      real(RK),         intent(in)    :: strid0
      real(RK),         intent(in)    :: tinc      ! dt increment factor
      real(RK),         intent(in)    :: tmin,tmax ! dt bounds
      real(RK),         intent(in)    :: tdec,tdrel,tdabs
      integer,          intent(in)    :: desire ! Desired number of timesteps
      integer,          intent(in)    :: levelm,leveld,steps2,tdage
      real(RK),         intent(inout) :: ynorm
      integer,          intent(inout) :: step
      integer,          intent(in)    :: names
      integer,          intent(out)   :: report
      character(len=*), intent(in)    :: name(names)
      character(len=*), intent(inout) :: signal
      real(RK), dimension(vars%N()), intent(in)    :: above,below,x
      real(RK), dimension(vars%N()), intent(inout) :: buffer,s0,s1,v0,v1,y0,y1,vsave
      type(twfunctions), intent(in)   :: functions

      real(RK)  :: change,condit,csave,dummy,high,low,stride
      integer   :: age,agej,count,first,last,length,number,route,xrepor
      intrinsic :: log10, max, min
      logical   :: exist,jacob,success,xsucce
      character(len=80) :: cword,jword,remark,yword

      character(len=*), parameter :: id = 'EVOLVE:  '

      ! SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.
      save

      ! Initialization.
      time    = .false. ! Turn off reverse communication flags.
      error   = .false. ! Turn off all completion status flags.
      report  = qnull
      success = .false.

      ! IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.
      if (signal /= ' ') then
         go to (1010, 1020, 1060, 1080, 1090) route
         error = .true.
         if (text>0) write (text, 21) id, route
         return
      end if

      ! Check the arguments
      call vars%check(error,id,text)
      if (error) return

      error = .not. (0 < desire)
      if (error) then
         if (text>0) write (text, 23) id, desire
         return
      end if

      error = .not. (tdec>=one .and. tinc>=one)
      if (error) then
         if (text>0) write (text, 24) id, tdec, tinc
         return
      end if

      ! stride bounds OK
      error = .not. (tmin>zero .and. tmax>=tmin)
      if (error) then
         if (text>0) write (text, 25) id, tmin, tmax
         return
      end if

      ! strid0 in bounds
      error = .not. (strid0>=tmin .and. tmax>=strid0)
      if (error) then
         if (text>0) write (text, 26) id, tmin, strid0, tmax
         return
      end if

      error = .not. step>=0
      if (error) then
         if (text>0) write (text, 27) id, step
         return
      end if

      error = tinc>one .and. .not. steps2>0
      if (error) then
         if (text>0) write (text, 28) id, steps2
         return
      end if

      ! ***** Time evolution. *****
      save_initial_solution: if (step<=0) then
         stride = strid0
         age    = 0

         ! Send latest solution to the problem handler
         buffer = v0
         call functions%save_sol(vars,buffer)

      endif save_initial_solution
      1010 continue

      ! Initialize timestepping
      exist = .false.
      first = step
      last  = step + desire

      ! Print header and initial function
      if (levelm>0 .and. text>0) then
         buffer = v0
         time = .false.
         call functions%fun(error,text,vars%points,time,stride,x,buffer)
         1020 continue
         ynorm = twnorm(vars%N(), buffer)
         call print_evolve_header(text,id,levelm,step,ynorm,stride)
      endif

      time_integration: do while (step < last)

          low  = tmin
          high = tmax

          ! Increase dt if possible
          if (age == steps2 .and. stride < high .and. tinc>one) then
             age    = 0
             exist  = .false.
             low    = stride * tdec
             stride = min (high, stride * tinc)
             if (levelm>1.and.text>0)            write (text, 10) id, step, yword, log10(stride)
          else
             if (levelm>1.and.text>0.and.step>0) write (text, 9) id, step, yword, log10(stride)
          end if

          xsucce = .false.
          change = zero

          newton_search: do while (change==zero .or. .not.xsucce)

              ! STORE THE LATEST SOLUTION SHOULD THE SEARCH FAIL
              call twcopy (vars%N(), v0, vsave)

              count = 0
              csave = zero
              jword = ' '
              jacob = .false.

              1060  continue
              if (jacob) then
                 exist = .true.
                 count = count + 1
                 csave = max(condit, csave)
                 if (csave == zero) then
                    write (jword, '(I3, 3X, A6)') count, '    NA'
                 else
                    write (jword, '(I3, 3X, F6.2)') count, log10 (csave)
                 end if
              end if

              call search(error, text, above, agej, below, buffer, vars, condit, exist, &
                          leveld - 1, levelm - 1, name, names, xrepor, &
                          s0, s1, signal, number, xsucce, v0, v1, tdabs, tdage, tdrel, y0, dummy, y1,&
                          x, functions)
              if (error) then
                 if (text>0) write (text, 29) id
                 return
              end if

              if (signal /= ' ') then
                 jacob = signal == 'PREPARE'
                 time  = .true.
        !        GO TO 1060 WHEN ROUTE = 3
                 route = 3
                 return
              end if

              ! Newton search unsuccessful
              if (.not. xsucce) then
                 if (levelm==1 .and. text>0) then
                    if (xrepor == qbnds) then
                       length = 6
                       remark = 'BOUNDS'
                    else if (xrepor == qdvrg) then
                       length = 7
                       remark = 'DIVERGE'
                    else
                       length = 1
                       remark = ' '
                    end if
                    write (text, 1) step + 1, log10 (stride), number, jword, remark (1:length)
                 end if

                 ! Retry with decreased strinde, if possible
                 decrease_dt: if (low<stride .and. one<tdec) then
                    age = 0
                    call twcopy (vars%N(), vsave, v0)
                    exist = .false.
                    high = stride / tinc
                    stride = max (low, stride / tdec)
                    if (levelm>1 .and. text>0) write (text, 11) id, step, yword, log10 (stride)
                    cycle newton_search
                 else
                    ! FAILURE
                    exit time_integration
                 end if decrease_dt

              end if

              ! If the solution is not changing and we can still increase stride, do it. Otherwise, failure.
              buffer = v0-vsave
              change = twnorm(vars%N(),buffer)
              call twlogr(cword, change)

              unchanged_solution: if (change == zero) then
                 if (levelm==1 .and. text>0) write (text, 4) step+1,'  ZERO',log10(stride),number,jword

                 increase_dt: if (tinc>one .and. stride < high) then
                    age = 0
                    exist = .false.
                    low = stride*tdec
                    stride = min (high, stride * tinc)
                    if (levelm>1 .and. text>0) write (text, 12) id, step, yword, log10 (stride)
                    cycle newton_search
                 else
                    exit time_integration
                 end if increase_dt

              end if unchanged_solution

          end do newton_search

          ! New timestep
          age  = age + 1
          step = step + 1

          ! Save latest solution for use by the function
          ! GO TO 1080 WHEN ROUTE = 4
          buffer = v0
          call functions%save_sol(vars,buffer)
          1080 continue

          ! Summary
          print_summary: if (levelm>0 .and. text>0) then
             buffer = v0
             time = .false.
             call functions%fun(error,text,vars%points,time,stride,x,buffer)
             1090 continue
             ynorm  = twnorm(vars%N(), buffer)
             call twlogr (yword, ynorm)
             if (levelm==1) write (text, 5) step,yword,cword,log10(stride),number,jword
          end if print_summary

      end do time_integration

      ! Epilogue.
      if (levelm>0 .and. text>0) then

         if (step == first) then
            write (text, 6) id
         else if (step == last) then
            if (levelm>1) then; write (text, 13) id, step, yword
            else;               write (text, 7) id; endif
         else
            if (levelm>1) then; write (text, 14) id, step, yword
            else;               write (text, 8) id; endif
         end if

         if (first<last .and. leveld==1) then
            write (text, 15) id
            call functions%show(error,text,v0,vars,.true.,x)
         end if
      end if

      ! SET THE COMPLETION STATUS FLAGS.
      success = step>first
      if (step<last) report = xrepor

      return

      ! Informative messages
      1 format(10X, i6, 21X, f6.2, 3X, i5, 3X, a12, 3X, a)
      4 format(10X, i6, 12X, a6, 3X, f6.2, 3X, i5, 3X, a12)
      5 format(10X, i6, 2(3X, a6), 3X, f6.2, 3X, i5, 3X, a12)
      6 format(/1X, a9, 'FAILURE.  NO TIME EVOLUTION.')
      7 format(/1X, a9, 'SUCCESS.  TIME EVOLUTION COMPLETED.')
      8 format(/1X, a9, 'PARTIAL SUCCESS.  TIME EVOLUTION INCOMPLETE.')
      9 format(/1X, a9, 'CONTINUE TIME EVOLUTION.' &
                 //10X, i10, '  LATEST TIME POINT' &
                  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE' &
                  /10X, f10.2, '  LOG10 STRIDE TO NEXT TIME POINT' &
                 //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')
      10 format(/1X, a9, 'CONTINUE TIME EVOLUTION WITH INCREASED STRIDE.' &
                 //10X, i10, '  LATEST TIME POINT' &
                  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE' &
                  /10X, f10.2, '  LOG10 INCREASED STRIDE TO NEXT TIME POINT' &
                 //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')
      11 format(/1X, a9, 'RETRY THE STEP WITH A DECREASED TIME STRIDE.' &
                 //10X, i10, '  LATEST TIME POINT' &
                  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE' &
                  /10X, f10.2, '  LOG10 DECREASED STRIDE TO NEXT TIME POINT' &
                 //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE, AGAIN.')
      12 format(/1X, a9, 'THE SOLUTION DID NOT CHANGE.  RETRYING THE STEP' &
                  /10X, 'WITH AN INCREASED TIME STRIDE.' &
                 //10X, i10, '  LATEST TIME POINT' &
                  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE' &
                  /10X, f10.2, '  LOG10 INCREASED STRIDE TO NEXT TIME POINT' &
                 //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE, AGAIN.')
      13 format(/1X, a9, 'SUCCESS.  TIME EVOLUTION COMPLETED.' &
                 //10X, i10, '  LAST TIME POINT' &
                  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE')
      14 format(/1X, a9, 'PARTIAL SUCCESS.  TIME EVOLUTION INCOMPLETE.' &
                 //10X, i10, '  LAST TIME POINT' &
                  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE')
      15 format(/1X, a9, 'THE LATEST SOLUTION:')

      ! Error messages
      21 format(/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
                 //10X, i10, '  ROUTE')
      23 format(/1X, a9, 'ERROR.  THE NUMBER OF TIME STEPS MUST BE POSITIVE.' &
                 //10X, i10, '  STEPS0 OR STEPS1, DESIRED NUMBER OF STEPS')
      24 format(/1X, a9, 'ERROR.  THE FACTORS FOR CHANGING THE TIME STRIDE' &
                  /10X, 'MUST BE NO SMALLER THAN 1.' &
                 //10X, 1p, e10.2, '  TDEC, DECREASE FACTOR', &
                  /10X, 1p, e10.2, '  TINC, INCREASE FACTOR')
      25 format(/1X, a9, 'ERROR.  THE BOUNDS ON THE TIME STRIDE ARE OUT OF' &
                  /10X, 'ORDER.' &
                 //10X, 1p, e10.2, '  TMIN, SHORTEST STRIDE' &
                  /10X, 1p, e10.2, '  TMAX, LONGEST STRIDE')
      26 format(/1X, a9, 'ERROR.  THE INITIAL TIME STRIDE MUST LIE BETWEEN' &
                  /10X, 'THE LOWER AND UPPER BOUNDS.' &
                 //10X, 1p, e10.2, '  TMIN, SHORTEST STRIDE' &
                  /10X, 1p, e10.2, '  STRID0, INITIAL STRIDE' &
                  /10X, 1p, e10.2, '  TMAX, LONGEST STRIDE')
      27 format(/1X, a9, 'ERROR.  THE COUNT OF TIME STEPS MUST BE ZERO OR' &
                  /10X, 'POSITIVE.' &
                 //10X, i10, '  STEP')
      28 format(/1X, a9, 'ERROR.  THE TIME STEPS BEFORE STRIDE INCREASES' &
                  /10X, 'MUST BE POSITIVE.' &
                 //10X, i10, '  STEPS2, TIME STEPS BEFORE STRIDE INCREASES')
      29 format(/1X, a9, 'ERROR.  SEARCH FAILS.')

      end subroutine evolve

      ! Perform the damped, modified Newton's search
      subroutine search(error, text, above, age, below, buffer, vars, condit, exist, &
                        leveld, levelm, name, names, report, s0, s1, signal, steps, &
                        success, v0, v1, xxabs, xxage, xxrel, y0, y0norm, y1, x, functions)

          type(twsize)    , intent(in)    :: vars
          integer         , intent(out)   :: report
          logical         , intent(out)   :: error
          logical         , intent(in)    :: exist ! Do we have a Jacobian already
          integer         , intent(in)    :: text
          real(RK)        , intent(in)    :: xxabs,xxrel ! settings
          integer         , intent(in)    :: xxage
          integer         , intent(in)    :: names
          character(len=*), intent(in)    :: name(names)
          character(len=*), intent(inout) :: signal
          real(RK), dimension(vars%N()), intent(in)    :: above,below,x
          real(RK), dimension(vars%N()), intent(inout) :: buffer,s0,s1,v0,v1,y0,y1
          type(twfunctions), intent(in)   :: functions

          real(RK) :: abs0,abs1, condit, deltab, deltad, rel0, rel1, s0norm, s1norm, &
                      value, y0norm, y1norm
          integer  :: age, entry, expone, leveld, levelm, number, route, steps
          intrinsic :: abs, int, log10, max, min, mod
          logical   :: force,success,converged,update_jac

          character(len=*), parameter :: id = 'SEARCH:  '

          ! SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.
          save

          ! *** Initialization. ***

          ! Turn off all completion flags
          error   = .false.
          report  = qnull
          success = .false.

          ! IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.
          if (signal /= ' ') then
             go to (2020, 2040, 2050, 2140, 2150) route
             error = .true.
             if (text>0) write (text, 5) id, route
             return
          end if

          ! Check number of variables
          call vars%check(error,id,text)
          if (error) return

          ! Check variable names
          error = .not. (names== 1 .or. names == vars%NVAR())
          number_of_names: if (error) then
              if (text>0) write (text,7) id,names,vars%comps,vars%groupa,vars%groupb,vars%NVAR()
              return
          end if number_of_names

          ! Check variable bounds
          error = any(.not.below<above)
          if (error) then
             !call print_invalid_bounds(id,text,name,groupa,groupb,comps,points,below,above)
             return
          end if

          ! Check the unknowns are valid
          error = any(.not.(below<=v0 .and. v0<=above))
          if (error) then
             call print_invalid_ranges(id,text,name,vars,below,above,v0)
             return
          end if

          ! Check tolerances
          error = .not. (zero<=xxabs .and. zero<=xxrel)
          if (error) then
             if (text>0) write (text,8) id, xxabs, xxrel
             return
          end if

          ! Check Jacobian update age
          error = .not. xxage>0
          if (error) then
             if (text>0) write (text, 9) id, xxage
             return
          end if

          ! PRINT THE HEADER.
          if (levelm>=1 .and. text>0) call print_search_header(text,id)

          !///////////////////////////////////////////////////////////////////////
          !
          !     SIR ISAAC NEWTON'S ALGORITHM.
          !
          !///////////////////////////////////////////////////////////////////////

          ! Number of steps
          number     = 0
          update_jac = .false.
          converged  = .false.

          newton_iterations: do while (.not.converged)

              ! Evaluate Jacobian at v0. Re-evaluate y0=F(v0) in case F changes when the Jacobian does.
              ! Solve J*s0 = v0; Evaluate relative and absolute errors
              update_jacobian: if (age>=xxage &                      ! Jacobian is too old
                             .or. (age>0 .and. update_jac) &         ! Jacobian not new, probably inaccurate
                             .or. (number==0 .and. .not.exist)) then ! First initialization

                  call twcopy(vars%N(),v0,buffer)
                  signal = 'PREPARE'

                  ! GO TO 2020 WHEN ROUTE = 1
                  route = 1
                  steps = number ! Copy the protected local variable
                  return

                  2020 continue
                  signal = ' '
                  age = 0

                  ! Turn off Jacobian request
                  update_jac = .false.

              endif update_jacobian

              ! EVALUATE Y0 := F(V0).  SOLVE J S0 = Y0.  EVAUATE ABS0 AND REL0.
              update_F: if (age==0 .or. number==0) then

                  buffer = v0
                  signal = 'RESIDUAL'
            !     GO TO 2040 WHEN ROUTE = 2
                  route = 2
                  steps = number ! Copy the protected local variable
                  return

                  2040 continue
                  signal = ' '
                  y0     = buffer
                  y0norm = twnorm(vars%N(),y0)

                  buffer = y0
                  signal = 'SOLVE'
            !     GO TO 2050 WHEN ROUTE = 3
                  route = 3
                  steps = number ! Copy the protected local variable
                  return
            2050  continue
                  signal = ' '
                  s0     = buffer
                  s0norm = twnorm(vars%N(),s0)

                  ! Check for success
                  call check_convergence(xxrel,xxabs,v0,s0,abs0,rel0,success)
                  if (success) exit newton_iterations

              endif update_F

              ! CHOOSE DELTAB.
              call newton_damping(v0,s0,above,below,deltab,force,entry,value)

              error = deltab < zero
              if (error) then
                  if (text>0) write (text,10) id, deltab
                  return
              end if

              !///  0 < DELTAB?
              if (.not. (zero < deltab)) then

                 ! If deltab becomes negative after some iterations, try a recovery by updating
                 ! the Jacobian.
                 update_jac = age>0
                 if (update_jac) cycle newton_iterations

                 if (levelm>0 .and. text>0) then
                    call print_newt_summary(text,number,y0norm,s0norm,abs0,rel0,deltab,deltad,condit)
                    call print_invalid_ranges(id,text,name,vars,below,above,v0,s0)
                 end if

                 report  = qbnds
                 success = .false.
                 steps = number ! Copy the protected local variable
                 return
              end if

              ! Exponential decay of the damping parameter.
              deltad = one ! Current
              expone = 0   ! Number of exponential decay iterations
              s1norm = huge(zero)

              ! Perform a simple backtracking iteration
              decay_iterations: do while (expone<=MAX_DECAY_ITERATIONS .and. .not.s1norm<=s0norm)

                  ! V1 := V0 - DELTAB DELTAD S0.
                  v1 = v0 - (deltab*deltad)*s0

                  ! KEEP V1 IN BOUNDS DESPITE ROUNDING ERROR.
                  v1 = min(max(v1,below),above)
                  if (expone==0 .and. force) v1(entry) = value

                  ! EVALUATE Y1 := F(V1)
                  buffer = v1
                  signal = 'RESIDUAL'
            !     GO TO 2140 WHEN ROUTE = 4
                  route = 4
                  steps = number ! Copy the protected local variable
                  return

            2140  continue
                  signal = ' '
                  y1     = buffer
                  y1norm = twnorm (vars%N(), y1)

                  ! SOLVE J*S1 = Y1
                  buffer = y1
                  signal = 'SOLVE'
            !     GO TO 2150 WHEN ROUTE = 5
                  route = 5
                  steps = number ! Copy the protected local variable
                  return

            2150  continue
                  signal = ' '
                  s1     = buffer
                  s1norm = twnorm (vars%N(), s1)

                  ! Check convergence (abs1, rel1)
                  call check_convergence(xxrel,xxabs,v1,s1,abs1,rel1,converged)

                  ! Check progress
                  deltad = half * deltad
                  expone = expone + 1

              end do decay_iterations

              is_diverging: if (expone>MAX_DECAY_ITERATIONS) then

                 ! Check if we can try restarting this iteration with an updated Jacobian
                 update_jac = age>0
                 if (update_jac) cycle newton_iterations

                 ! Failed too many times.
                 if (levelm>0 .and. text>0) then
                    call print_newt_summary(text,number,y0norm,s0norm,abs0,rel0,deltab,deltad,condit)
                    write (text, 1) id
                 end if

                 report  = qdvrg
                 success = .false.
                 steps   = number ! Copy the protected local variable
                 return

              end if is_diverging

              ! Print summary.
              if (levelm>0 .and. text>0) &
              call print_newt_summary(text,number,y0norm,s0norm,abs0,rel0,deltab,deltad,condit)

              ! Advance step
              age    = age + 1
              number = number + 1
              s0     = s1; s0norm = s1norm
              v0     = v1; y0norm = y1norm
              y0     = y1
              abs0   = abs1
              rel0   = rel1

          end do newton_iterations

          ! SUCCESS!

          ! Print summary.
          if (levelm>0 .and. text>0) then
             call print_newt_summary(text,number,y0norm,s0norm,abs0,rel0,deltab,deltad,condit)

             if (leveld>0) then
                ! Ask to display the final solution
                write (text, 3) id
                call functions%show(error,text,v0,vars,.true.,x)
             else
                write (text, 4) id
             end if
          end if

          signal  = ' '
          success = .true.
          steps   = number ! Copy the protected local variable
          return

          ! Informative messages.
          1 format(/1X, a9, 'FAILURE.  THE SEARCH DIVERGES.')
          3 format(/1X, a9, 'SUCCESS.  THE SOLUTION:')
          4 format(/1X, a9, 'SUCCESS.')

          ! Error messages.
          5 format(/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
                 //10X, i10, '  ROUTE')
          7 format(/1X, a9, 'ERROR.  THE NUMBER OF NAMES IS WRONG.' &
                 //10X, i10, '  NAMES' &
                 //10X, i10, '  COMPS, COMPONENTS' &
                  /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                  /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                  /10X, i10, '  TOTAL NUMBER')
          8 format(/1X, a9, 'ERROR.  THE BOUNDS FOR THE ABSOLUTE AND RELATIVE' &
                  /10X, 'CONVERGENCE TESTS MUST BE ZERO OR POSITIVE.' &
                 //10X, 1p, e10.2, '  SSABS OR TDABS, ABSOLUTE ERROR' &
                  /10X, 1p, e10.2, '  SSREL OR TDREL, RELATIVE ERROR')
          9 format(/1X, a9, 'ERROR.  THE RETIREMENT AGE OF THE JACOBIAN MATRIX' &
                  /10X, 'MUST BE POSITIVE.' &
                 //10X, i10, '  SSAGE OR TDAGE, MATRIX RETIREMENT AGE')
         10 format(/1X, a9, 'ERROR.  THE DAMPING COEFFICIENT FOR STAYING' &
                  /10X, 'IN BOUNDS IS NEGATIVE.' &
                 //10X, 1p, e10.2, '  DELTA B')

      end subroutine search

      pure subroutine check_convergence(RTOL,ATOL,v0,s0,abs0,rel0,converged)
          real(RK), intent(in)  :: RTOL,ATOL
          real(RK), intent(in)  :: v0(:),s0(size(v0))
          real(RK), intent(out) :: abs0,rel0 ! relative, absolute accuracy
          logical , intent(out) :: converged

          real(RK) :: sj,vj
          integer :: j
          intrinsic :: max, abs, size

          abs0 = zero
          rel0 = zero
          do j = 1, size(v0)
             sj = abs(v0(j) - (v0(j) - s0(j)))
             vj = abs(v0(j))
             if (sj>RTOL*vj) abs0 = max (abs0, sj)
             if (sj>ATOL .and. vj>zero) rel0 = max (rel0, sj / vj)
          end do

          converged = rel0<=RTOL .and. abs0<=ATOL

      end subroutine check_convergence

      ! Deltab is the largest damping coefficient in [0,1] that keeps v1 within bounds.
      ! If v1 belongs on the boundary, then provisions are made to force it there despite
      ! rounding error.
      pure subroutine newton_damping(v0,s0,above,below,deltab,force,entry,value)
         real(RK), intent(out) :: deltab !> Damped Newton parameter
         logical , intent(out) :: force  !>
         integer , intent(out) :: entry
         real(RK), intent(out) :: value
         real(RK), intent(in)  :: v0(:)
         real(RK), intent(in), dimension(size(v0)) :: s0,above,below

         integer :: j
         real(RK) :: temp

         ! Initialize no damping
         deltab = one
         entry  = 0
         force  = .false.
         value  = zero

         do j = 1, size(v0)
            if (s0(j)>max(zero, v0(j) - below(j))) then
                temp = (v0(j) - below(j)) / s0(j)
                if (temp<deltab) then
                   deltab = temp
                   entry = j
                   force = .true.
                   value = below(j)
                end if
             else if (s0(j) < min (zero, v0(j) - above(j))) then
                temp = (v0(j) - above(j)) / s0(j)
                if (temp<deltab) then
                   deltab = temp
                   entry = j
                   force = .true.
                   value = above(j)
                end if
             end if
         end do

      end subroutine newton_damping

      subroutine print_newt_summary(text,number,y0nrm,s0nrm,eabs,erel,db,dd,condit)
          integer, intent(in) :: text
          integer, intent(in) :: number ! of iterations
          real(RK), intent(in) :: y0nrm,s0nrm
          real(RK), intent(in) :: eabs,erel
          real(RK), intent(in) :: db,dd
          real(RK), optional, intent(in) :: condit

          character(len=16) :: column(7)

          if (text==0) return

          column = ' '

          call twlogr (column(1), y0nrm)
          call twlogr (column(3), s0nrm)
          call twlogr (column(4), eabs)
          call twlogr (column(5), erel)
          if (db /= one) call twlogr (column(6), db)
          if (dd /= one) call twlogr (column(7), dd)
          if (present(condit)) call twlogr (column(2),condit)

          write (text, 1) number, column

          1 format(10X, i6, 3(3X, a6), 2(3X, a6, 2X, a6))

      end subroutine print_newt_summary

      subroutine print_evolve_header(text,id,levelm,step,ynorm,stride)
          integer, intent(in) :: text,levelm,step
          real(RK), intent(in) :: ynorm,stride
          character(*), intent(in) :: id

          character(len=80) :: header(2,3),yword

          !               123456789_123456789_123456789_123456789_1234
          !               123456   123456   123456   123456   12345
          header(1, 1) = '  TIME   LOG10                      NEWTON S'
          header(1, 2) = ' POINT   ------------------------   --------'
          header(1, 3) = 'NUMBER   NORM F   CHANGE   STRIDE   STEPS   '

          !               123456789_123456789_1
          !               123   123456   123456
          header(2, 1) = 'EARCH                '
          header(2, 2) = '---------------------'
          header(2, 3) = 'J''S   COND J   REMARK'

          call twlogr(yword, ynorm)

          if (levelm==1) then
            if (step==0) then
               write (text, 11) id, header, step, yword
            else
               write (text, 12) id, header, step, yword
            end if
         else if (levelm>1 .and. step==0) then
            if (step == 0) then
               write (text, 21) id, step, yword, log10(stride)
            else
               write (text, 22) id, step, yword, log10(stride)
            end if
         end if

         11 format(/1X, a9, 'BEGIN TIME EVOLUTION.' &
               /3(/10X, a44, a21) &
                  /10X, i6, 3X, a6)
         12 format(/1X, a9, 'CONTINUE TIME EVOLUTION.' &
               /3(/10X, a44, a21) &
                  /10X, i6, 3X, a6)
         21 format(/1X, a9, 'BEGIN TIME EVOLUTION.' &
                 //10X, i10, '  LATEST TIME POINT' &
                  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE' &
                  /10X, f10.2, '  LOG10 STRIDE TO NEXT TIME POINT' &
                 //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')
         22 format(/1X, a9, 'CONTINUE TIME EVOLUTION.' &
                 //10X, i10, '  LATEST TIME POINT' &
                  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE' &
                  /10X, f10.2, '  LOG10 STRIDE TO NEXT TIME POINT' &
                 //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')

      end subroutine print_evolve_header


      subroutine print_search_header(iunit,id)
          integer, intent(in) :: iunit
          character(*), intent(in) :: id

          character(len=80) :: header(3,2)

          ! PRINT THE HEADER.
          !               123456789_123456789_123456789_123456789_1234
          !               123456   123456   123456   123456   123456
          header(1, 1) = '         LOG10                              '
          header(2, 1) = '  SLTN   -----------------------------------'
          header(3, 1) = 'NUMBER   NORM F   COND J   NORM S      ABS A'
          header(1, 2) = '                       '
          header(2, 2) = '-----------------------'
          header(3, 2) = 'ND REL    DELTA B AND D'
          if (iunit>0) write(iunit,1) id, header

          1 format(/1X, a9, 'SOLVE NONLINEAR, NONDIFFERENTIAL EQUATIONS.'/4(/10X,a44,a23)/)

      end subroutine print_search_header

      ! TWOPNT driver.
      subroutine twopnt(setup, error, text, versio, vars, &
                        above, active, below, buffer, condit,  &
                        work, mark, name, names, report, signal, stride, time, u, x, &
                        functions)

      type(twcom) , intent(inout) :: setup
      logical     , intent(out)   :: error
      integer     , intent(in)    :: text
      character(*), intent(in)    :: versio
      type(twsize), intent(inout) :: vars
      character(*), intent(inout) :: signal,report
      integer     , intent(in)    :: names

      character(*), intent(inout) :: name(names) ! Names of the variables
      real(RK)    , intent(inout), dimension(vars%NVAR()) :: above,below
      logical     , intent(inout) :: active(*),mark(*)
      real(RK)    , intent(inout) :: buffer(vars%NMAX())
      real(RK)    , intent(inout) :: condit
      type(twwork), intent(inout) :: work
      real(RK)    , intent(inout) :: u(vars%NMAX())
      real(RK)    , intent(inout) :: x(:)
      type(twfunctions), intent(in) :: functions

      ! Local variables
      character(*), parameter :: id = 'TWOPNT:  '

      real(RK) ::  maxcon, ratio(2), stride, temp, ynorm
      type(twstat) :: stats
      integer :: age, comps, desire, j, jacobs, label, length, nsteps, psave, qtask, qtype, return, &
                 route, steps, xrepor
      intrinsic :: max
      logical :: allow, exist, found, satisf, time

      character(len=80) :: column(3),header(6),string

      ! SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.
      save

      !***** INITIALIZATION. *****

      time = .false. ! Turn off all reverse communication flags

      ! If this is a return call, continue where the program had paused.
      if (signal /= ' ') then
         go to (9912, 9922, 9932) route
         error = .true.

         if (text>0) write (text, 01) id, route
         return
      end if

      !***** ENTRY BLOCK.  INITIALIZE A NEW PROBLEM. *****

      ! Turn off all status reports.
      error = .false.
      report = ' '

      ! Check the version.
      call check_version(versio,text,error); if (error) return

      ! Additional settings initialization
      if (.not.setup%padd) setup%ipadd = vars%pmax

      ! Print entry banner
      string = vnmbr(vnmbrs)
      if (setup%levelm>0 .and. text>0) write (text, 10001) id, precision_flag(), trim(string)

      ! CHECK THE ARGUMENTS.
      error = .not. (setup%leveld <= setup%levelm)
      printing_levels: if (error) then
          if (text>0) write (text, 5) id, setup%leveld, setup%levelm
          return
      end if printing_levels

      ! Check variable sizes
      call vars%check(error,id,text); if (error) return

      error = .not. (names == 1 .or. names == vars%NVAR())
      number_of_names: if (error) then
          if (text>0) write (text, 9) id,names,vars%comps,vars%groupa,vars%groupb,vars%NVAR()
          return
      end if number_of_names

      ! Check variable bounds
      error = any(.not.below<above)
      if (error) then
          call print_invalid_bounds(id,text,name,vars,below,above)
          return
      end if

      ! PARTITION THE INTEGER WORK SPACE.
      call work%init(text,vars,error)
      if (error) return

      ! ONE-TIME INITIALIZATION.
      allow = .true.  ! Allow further time evolution
      qtask = qentry  ! Present task: ENTRY
      found = .true.  ! SOLUTION FLAG

      ! Init solver statistics
      call stats%new(vars%points)

      ! EXPAND THE BOUNDS.
      call work%load_bounds(above,below,vars)

      ! SAVE THE INITIAL SOLUTION.
      psave = vars%points
      if (setup%adapt .and. vars%points>0) call twcopy(vars%points, x, work%xsave)
      call twcopy (vars%N(), u, work%usave)

!     GO TO 1090 WHEN RETURN = 1
      return = 1
      ! Save the last solution
      call twcopy (vars%N(), u, buffer)
      signal = 'SAVE'
      ! GO TO 9912 WHEN ROUTE = 1
      route = 1
      return


1090  continue

!///  PRINT LEVELS 11, 21, AND 22.
      if (setup%leveld>0 .and. text>0) then
         write (text, 10002) id, 'INITIAL GUESS:'
         ! GO TO 1100 WHEN RETURN = 2
         return = 2
         go to 9921
      end if
1100  continue

!///  PRINT LEVEL 10 AND 11.

!                  123456789_123456789_123456789_1234
!                  12345678   123456  123456   123456
      header(1) = '            LOG10   LOG10         '
      header(2) = '    TASK   NORM F  COND J   REMARK'

      if (setup%levelm == 1 .and. text>0) then
         if (setup%leveld>0) write (text, 10002) id, 'SOLVE THE PROBLEM.'
         write (text, 10003) (header(j), j = 1, 2)
         ! GO TO 1110 WHEN LABEL = 1
         label = 1
         go to 7010
      end if
      1110  continue

      !///////////////////////////////////////////////////////////////////////
      !
      !     DECISION BLOCK.  THE PREVIOUS TASK DETERMINES THE NEXT.
      !
      !///////////////////////////////////////////////////////////////////////

      new_task: do

          ! ENTRY WAS THE PREVIOUS TASK.
          qtask = twopnt_next_task(setup,stats,qtask,found,satisf,allow,report,desire,error)
          if (error) then
              if (text>0) write (text, 14) id
              return
          end if

          ! Branch to the next task.
          select case (qtask)
              case (qexit) ! *** EXIT BLOCK. ***

                  ! Complete statistics for the last grid
                  call stats%tock(qgrid)

                  ! Restore the solution.
                  if (report /= ' ') then
                     ! BE CAREFUL NOT TO ASSIGN A VALUE TO A PARAMETER
                     if (vars%points /= psave) vars%points = psave
                     if (setup%adapt .and. vars%points>0) call twcopy(vars%points, work%xsave, x)
                     call twcopy(vars%N(), work%usave, u)
                  end if

                  ! PRINT LEVEL 11 OR 21.

                  ! SAVE THE STATUS REPORTS DURING REVERSE COMMUNICATION
                  string = report

                  if (setup%leveld == 1 .and. text>0) then
                     write (text, 10002) id, 'FINAL SOLUTION:'
                     ! GO TO 3020 WHEN RETURN = 3
                     return = 3
                     go to 9921
                  end if
            3020  continue

                  ! RESTORE THE STATUS REPORTS AFTER REVERSE COMMUNICATION
                  report = string

                  ! Complete the total time statistics
                  call stats%tock(qtotal)

                  ! TOP OF THE REPORT BLOCK.
                  print_reports: if (setup%levelm>0 .and. text>0) then

                     call stats%print_stats(text,setup%adapt)

                      ! Report the completion status.
                      if (setup%levelm>0) then
                         if (report == ' ') then
                            write (text, 10010) id
                         else if (report == 'NO SPACE') then
                            write (string, '(I10)') vars%points
                            write (text, 10011) id, trim(string), ratio, setup%toler1, setup%toler2
                         else if (report == 'NOT SOLVED') then
                            write (text, 10012) id
                         else if (report == 'SOME SOLVED') then
                            write (string, '(I10)') vars%points
                            write (text, 10013) id, trim(string), ratio, setup%toler1, setup%toler2
                         else
                             error = .true.
                             if (text>0) write (text, 16) id
                             return
                         end if
                      end if

                  end if print_reports

                  return

              case (qsearc) ! *** SEARCH BLOCK. ***

                  ! INITIALIZE STATISTICS ON ENTRY TO THE SEARCH BLOCK.
                  call stats%tick(qsearc)
                  jacobs = 0
                  maxcon = zero

                  ! PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE SEARCH BLOCK.
                  if (setup%levelm>1 .and. text>0) write (text, 10014) id

                  ! PREPARE TO CALL SEARCH.

                  ! Save the solution should the search fail
                  call twcopy (vars%N(), u, work%vsave)
                  exist = .false.

                  ! CALL SEARCH.
                  age = 0
            4020  continue

                  call search(error, text, &
                             work%above, age, work%below, buffer, vars, condit, &
                             exist, setup%leveld - 1, setup%levelm - 1, name, names, &
                             xrepor, work%s0, work%s1, signal, nsteps, found, &
                             u, work%v1, setup%ssabs, setup%ssage, setup%ssrel, work%y0, ynorm, &
                             work%y1, x, functions)
                  if (error) then
                      if (text>0) write (text, 17) id
                      return
                  end if

                  ! PASS REQUESTS FROM SEARCH TO THE CALLER.
                  if (signal /= ' ') then
                     ! GO TO 4020 WHEN RETURN = 4
                     return = 4
                     go to 9931
                  end if

                  ! REACT TO THE COMPLETION OF SEARCH.
                  save_or_restore: if (found) then

                     psave = vars%points
                     if (setup%adapt .and. vars%points>0) call twcopy (vars%points, x, work%xsave)
                     call twcopy(vars%N(), u, work%usave)

                     ! GO TO 4030 WHEN RETURN = 5
                     return = 5
                     ! Save the last solution
                     call twcopy (vars%N(), u, buffer)
                     signal = 'SAVE'
                     ! GO TO 9912 WHEN ROUTE = 1
                     route = 1
                     return

                     4030  continue

                  else save_or_restore
                     ! RESTORE THE SOLUTION
                     call twcopy(vars%N(), work%vsave, u)
                  end if save_or_restore

                  ! COMPLETE STATISTICS FOR THE SEARCH BLOCK.
                  call stats%tock(qsearc)

                  ! PRINT LEVEL 10 OR 11 ON EXIT FROM THE SEARCH BLOCK.
                  if (setup%levelm == 1 .and. text>0) then
                     ! GO TO 4040 WHEN LABEL = 2
                     label = 2
                     exit new_task
                  end if
            4040  continue

            !///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE SEARCH BLOCK.

                  if (setup%levelm>1) then
                     if (found) then
                        if (text>0) write (text, 10015) id
                     else
                        if (text>0) write (text, 10016) id
                     end if
                  end if

            !///  BOTTOM OF THE SEARCH BLOCK.

                  cycle new_task


              case (qrefin) ! *** REFINE BLOCK. ***

                  ! INITIALIZE STATISTICS ON ENTRY TO THE REFINE BLOCK.
                  call stats%tick(qrefin)

                  ! PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE REFINE BLOCK.
                  if (setup%levelm>1 .and. text>0) write (text, 10017) id

                  ! PREPARE TO CALL REFINE.

                  ! Save group B values (will be shifted by the new grid size)
                  if (vars%groupb>0) work%vsave(:vars%groupb) = u(vars%idx_B())
                  exist = .false.

                  ! CALL REFINE.
                  5030  continue
                  ! Refine only modifies comps*points
                  call refine(error, text, active, &
                              buffer(vars%groupa+1), vars%comps, setup%leveld - 1, setup%levelm - 1, mark, &
                              found, setup%ipadd, vars%pmax, vars%points, ratio, work%ratio1, work%ratio2,   &
                              signal, satisf, setup%toler0, setup%toler1, setup%toler2, u(vars%groupa + 1), &
                              work%vary1, work%vary2, work%vary, x, functions)
                  if (error) then
                      if (text>0) write (text, 18) id
                      return
                  end if

                  ! SERVICE REQUESTS FROM REFINE: PASS REQUESTS TO THE CALLER.
                  if (signal /= ' ') then
                     ! INSERT THE GROUP A AND B UNKNOWNS
                     buffer(1:vars%groupa) = u(1:vars%groupa)
                     buffer(vars%idx_B()) = work%vsave(1:vars%groupb)

                     ! GO TO 5030 WHEN RETURN = 6
                     return = 6
                     go to 9931
                  end if

                  ! REACT TO THE COMPLETION OF REFINE.
                  refine_found: if (found) then

                     ! Initialize statistics for the new grid
                     call stats%new_grid(vars%points)

                     ! Insert the group B values
                     if (vars%groupb>0) u(vars%idx_B()) = work%vsave(:vars%groupb)

                     ! Expand bounds to new grid size
                     call work%load_bounds(above,below,vars)

                     ! SAVE THE LATEST SOLUTION
                     ! GO TO 5100 WHEN RETURN = 7
                     return = 7
                     call twcopy (vars%N(),u,buffer)
                     signal = 'SAVE'
                     ! GO TO 9912 WHEN ROUTE = 1
                     route = 1
                     return
                     5100 continue

                  endif refine_found

                  ! COMPLETE STATISTICS FOR THE REFINE BLOCK.
                  call stats%tock(qrefin)

                  ! PRINT LEVEL 10 OR 11 ON EXIT FROM THE REFINE BLOCK.
                  if (setup%levelm == 1 .and. text>0) then
                     write (text, '()')
                     ! GO TO 5120 WHEN LABEL = 3
                     label = 3
                     exit new_task
                  end if
                  5120 continue

                  ! PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE REFINE BLOCK.
                  if (setup%levelm>1) then
                     if (found) then
                        if (text>0) write (text, 10018) id
                     else
                        if (text>0) write (text, 10019) id
                     end if
                  end if

                  ! BOTTOM OF THE REFINE BLOCK.
                  cycle new_task

              case (qtimst) ! *** EVOLVE BLOCK. ***

                  ! INITIALIZE STATISTICS ON ENTRY TO THE EVOLVE BLOCK.
                  call stats%tick(qtimst)
                  jacobs = 0
                  maxcon = zero
                  steps = stats%step

                  ! PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE EVOLVE BLOCK.
                  if (setup%levelm>1 .and. text>0) write (text, 10020) id

                  ! CALL EVOLVE.
                  6020  continue
                  call evolve(error, text, work%above, work%below, &
                         buffer, vars, condit, desire, setup%leveld - 1, &
                         setup%levelm - 1, name, names, xrepor, work%s0, work%s1, signal, &
                         stats%step, setup%steps2, setup%strid0, stride, found, setup%tdabs, &
                         setup%tdage, setup%tdec, setup%tdrel, time, setup%tinc, setup%tmax, &
                         setup%tmin, u, work%v1, work%vsave, work%y0, work%y1, ynorm, x, functions)
                  if (error) then
                      if (text>0) write (text, 19) id
                      return
                  end if

                  ! PASS REQUESTS FROM EVOLVE TO THE CALLER.
                  if (signal /= ' ') then
                     ! GO TO 6020 WHEN RETURN = 8
                     return = 8
                     go to 9931
                  end if

                  ! REACT TO THE COMPLETION OF EVOLVE.
                  if (found) then
                     ! SAVE THE LATEST SOLUTION
                     ! GO TO 6030 WHEN RETURN = 9
                     return = 9

                     ! Save the last solution
                     call twcopy (vars%N(), u, buffer)
                     signal = 'SAVE'
                     ! GO TO 9912 WHEN ROUTE = 1
                     route = 1
                     return

                  end if
                  6030  continue

                  ! ALLOW FURTHER TIME EVOLUTION.
                  allow = xrepor == qnull

                  ! COMPLETE STATISTICS FOR THE EVOLVE BLOCK.
                  call stats%tock(qtimst)

                  steps = stats%step - steps

                  ! PRINT LEVEL 10 OR 11 ON EXIT FROM THE EVOLVE BLOCK.
                  if (setup%levelm==1 .and. text>0) then
                     ! GO TO 6040 WHEN LABEL = 4
                     label = 4
                     exit new_task
                  end if
                  6040  continue

                  ! PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE EVOLVE BLOCK.
                  if (setup%levelm>1) then
                     if (found) then
                        if (text>0) write (text, 10021) id
                     else
                        if (text>0) write (text, 10022) id
                     end if
                  end if

                  ! BOTTOM OF THE EVOLVE BLOCK.

              case default
                  error = .true.
                  if (text>0) write (text, 15) id
                  return
          end select

      end do new_task

!///////////////////////////////////////////////////////////////////////
!
!     BLOCK TO PRINT LOG LINES.
!
!///////////////////////////////////////////////////////////////////////

      7010 continue

      column(:) = ' '
      string = ' '

      ! COLUMN 1: NAME OF THE TASK
      column(1) = trim(qname(qtask))

      ! COLUMN 2: NORM OF THE STEADY STATE FUNCTION
      if (found) then
         ! GO TO 7030 WHEN RETURN = 10
         return = 10

         ! EVALUATE THE STEADY STATE FUNCTION.
         call stats%tick(qfunct)
         call twcopy (vars%N(), from=u, to=buffer)
         time = .false.
         call functions%fun(error,text,vars%points,time,stride,x,buffer)
         call stats%tock(qfunct,event=.true.)

         go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030) return
         error = .true.
         if (text>0) write (text, 21) id, return
         return

         7030 continue
         temp = twnorm (vars%N(), buffer)
         call twlogr (column(2),temp)
      endif

      ! COLUMN 3: LARGEST CONDITION NUMBER
      if (any(qtask==[qsearc,qtimst]) .and. maxcon/=zero) call twlogr (column(3), maxcon)

      ! REMARK
      if (qtask == qsearc) then
         if (xrepor == qdvrg) then
            string = 'DIVERGING'
         else if (xrepor == qnull) then
            if (nsteps == 1) then
               write (string, '(I10, A)') nsteps, ' SEARCH STEP'
            else
               write (string, '(I10, A)') nsteps, ' SEARCH STEPS'
            end if
         else if (xrepor == qbnds) then
            string = 'GOING OUT OF BOUNDS'
         else
            string = '?'
         end if
      else if (qtask == qtimst) then
         if (xrepor == qbnds .or. xrepor == qdvrg .or. xrepor == qnull) then
            write (string, '(I10, A, 1P, E10.1, A)') steps, ' TIME STEPS, ', stride, ' LAST STRIDE'
         else
            string = '?'
         end if
      else if (qtask == qentry .and. setup%adapt) then
         write (string, '(I10, A)') vars%points, ' GRID POINTS'
      else if (qtask == qrefin) then
         if (found) then
            write (string, '(F10.2, A, F10.2, A, I10, A)') &
               ratio(1), ' AND ', ratio(2), ' RATIOS, ', vars%points, ' GRID POINTS'
         else
            write (string, '(F10.2, A, F10.2, A)') &
               ratio(1), ' AND ', ratio(2), ' RATIOS'
         end if
      end if

      call twsqez (length, string)
      if (text>0) write (text, 10023) column, string (1 : length)

      go to (1110, 4040, 5120, 6040) label
      error = .true.
      if (text>0) write (text, 20) id, label
      return

!///////////////////////////////////////////////////////////////////////
!
!     REQUEST REVERSE COMMUNICATION.
!
!///////////////////////////////////////////////////////////////////////

!///  SAVE THE SOLUTION.

9912  continue
      signal = ' '

      go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030) return
      error = .true.
      if (text>0) write (text, 21) id, return
      return

!///  PRINT THE LATEST SOLUTION.

9921  continue

      call twcopy (vars%N(),u,buffer)
      call functions%show(error,text,buffer,vars,.true.,x)
      9922 continue
      signal = ' '

      go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030) return
      error = .true.
      if (text>0) write (text, 21) id, return
      return

!///  PASS REQUESTS FROM SEARCH, REFINE, OR EVOLVE TO THE CALLER.

9931  continue

      ! THIS MUST BE SAVED TO GATHER STATISTICS AT REENTRY.
      ! THE REVERSE COMMUNICATION FLAGS WILL NOT BE SAVED BECAUSE THEY ARE CLEARED AT EVERY ENTRY.
      qtype = identify_request(signal)

      ! COUNT THE JACOBIANS
      if (qtype == qjacob) jacobs = jacobs + 1
      call stats%tick(qtype)

      !  GO TO 9932 WHEN ROUTE = 3
      route = 3
      return
9932  continue

!     SAVE THE CONDITION NUMBER
      if (qtype == qjacob) maxcon = max (maxcon, condit)

      call stats%tock(qtype,event=.true.)

      go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030) return
      error = .true.
      if (text>0) write (text, 21) id, return
      return

          ! Error messages
           1 format(/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
                  //10X, i10, '  ROUTE')
           5 format(/1X, a9, 'ERROR.  THE PRINTING LEVELS ARE OUT OF ORDER.' &
                   /10X, 'LEVELD CANNOT EXCEED LEVELM.' &
                  //10X, i10, '  LEVELD, FOR SOLUTIONS' &
                   /10X, i10, '  LEVELM, FOR MESSAGES')
           9 format(/1X, a9, 'ERROR.  THE NUMBER OF NAMES IS WRONG.' &
                  //10X, i10, '  NAMES' &
                  //10X, i10, '  COMPS, COMPONENTS' &
                   /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                   /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                   /10X, i10, '  TOTAL NUMBER')
          14 format(/1X, a9, 'ERROR.  NEITHER THE INITIAL TIME EVOLUTION NOR THE' &
                   /10X, 'SEARCH FOR THE STEADY STATE IS ALLOWED.')
          15 format(/1X, a9, 'ERROR.  UNKNOWN TASK.')
          16 format(/1X, a9, 'ERROR.  UNKNOWN REPORT CODE.')
          17 format(/1X, a9, 'ERROR.  SEARCH FAILS.')
          18 format(/1X, a9, 'ERROR.  REFINE FAILS.')
          19 format(/1X, a9, 'ERROR.  EVOLVE FAILS.')
          20 format(/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
                  //10X, i10, '  LABEL')
          21 format(/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
                  //10X, i10, '  RETURN')

       ! Informative messages
       10001 format(/1X, a9, a, ' (TWO POINT BOUNDARY VALUE PROBLEM) SOLVER,' &
                    /10X, 'VERSION ', a,' OF APRIL 1998 BY DR. JOSEPH F. GRCAR.')
       10002 format(/1X, a9, a)
       10003 format(3(/10X, a35)/)
       10010 format(/1X, a9, 'SUCCESS.  PROBLEM SOLVED.')
       10011 format(/1X, a9, 'FAILURE.  A SOLUTION WAS FOUND FOR A GRID WITH ', a &
                   /10X, 'POINTS, BUT ONE OR BOTH RATIOS ARE TOO LARGE.' &
                  //22X, '   RATIO 1     RATIO 2' &
                  //10X, '     FOUND', 2F12.2 &
                   /10X, '   DESIRED', 2F12.2 &
                  //10X, 'A LARGER GRID COULD NOT BE FORMED.')
       10012 format(/1X, a9, 'FAILURE.  NO SOLUTION WAS FOUND.')
       10013 format(/1X, a9, 'FAILURE.  A SOLUTION WAS FOUND FOR A GRID WITH ', a &
                   /10X, 'POINTS, BUT ONE OR BOTH RATIOS ARE TOO LARGE.' &
                  //22X, '   RATIO 1     RATIO 2' &
                  //10X, '     FOUND', 2F12.2 &
                   /10X, '   DESIRED', 2F12.2 &
                  //10X, 'A SOLUTION COULD NOT BE FOUND FOR A LARGER GRID.')
       10014 format(/1X, a9, 'CALLING SEARCH TO SOLVE THE STEADY STATE PROBLEM.')
       10015 format(/1X, a9, 'SEARCH FOUND THE STEADY STATE.')
       10016 format(/1X, a9, 'SEARCH DID NOT FIND THE STEADY STATE.')
       10017 format(/1X, a9, 'CALLING REFINE TO PRODUCE A NEW GRID.')
       10018 format(/1X, a9, 'REFINE SELECTED A NEW GRID.')
       10019 format(/1X, a9, 'REFINE DID NOT SELECT A NEW GRID.')
       10020 format(/1X, a9, 'CALLING EVOLVE TO PERFORM TIME EVOLUTION.')
       10021 format(/1X, a9, 'EVOLVE PERFORMED A TIME EVOLUTION.')
       10022 format(/1X, a9, 'EVOLVE DID NOT PERFORM A TIME EVOLUTION.')
       10023 format(10X, a8, 3X, a6, 2X, a6, 3X, a)


      end subroutine twopnt

      ! Perform automatic grid selection
      subroutine refine(error, text, active, buffer, comps, leveld, levelm, mark, newx, &
                        padd, pmax, points, ratio, ratio1, ratio2, signal, success, toler0, &
                        toler1, toler2, u, vary1, vary2, weight, x, functions)

          integer,  intent(in)     :: text,leveld,levelm
          integer,  intent(in)     :: pmax,padd
          integer,  intent(in)     :: comps
          integer,  intent(inout)  :: points
          logical,  intent(inout)  :: error, newx, success
          logical,  intent(inout)  :: active(comps)
          logical,  intent(inout)  :: mark(pmax)
          real(RK), intent(inout)  :: ratio1(pmax),ratio2(pmax),ratio(2),x(pmax)
          real(RK), intent(inout)  :: buffer(comps*pmax),u(comps,pmax)
          real(RK), intent(in)     :: toler0,toler1,toler2
          integer,  intent(inout), dimension(pmax) :: vary1,vary2,weight
          character(len=*), intent(inout) :: signal
          type(twfunctions), intent(in) :: functions

          character(len=*), parameter :: id = 'REFINE:  '

          character(len=80) :: word
          real(RK) :: differ,left,length,lower,maxmag,mean,range,right,temp,temp1,temp2,upper
          integer  :: act,counted,former,itemp,j,k,least,more,most,new,old,route,signif,total
          intrinsic :: abs, max, min, count, minval, maxval

          ! Save local variables during returns for reverse communication.
          save

          ! ***** prologue. *****

          ! every-time initialization: turn off all completion status flags.
          error   = .false.
          newx    = .false.
          success = .false.

          ! if this is a return call, then continue where the program paused.
          if (signal/=' ') then
             go to (4060, 5040) route
             error = .true.
             ! invalid route
             if (text>0) write (text, 101) id, route
             return
          end if

          ! Levelm printing.
          if (0<levelm .and. text>0) write (text, 1) id

          ! Check the arguments.
          error = .not. (1<=comps .and. 2<=points)
          if (error) then
               if (text>0) write (text, 102) id, comps, points
               return
          end if

          error = .not. (0<=padd)
          if (error) then
               if (text>0) write (text, 103) id, padd
               return
          endif

          error = .not. (points<=pmax)
          if (error) then
               if (text>0) write (text, 104) id, points, pmax
               return
          end if

          ! Check there is at least one variable that affects grid adaption
          counted = count(active)
          error = .not. (counted>=1)
          if (error) then
               if (text>0) write (text, 105) id
               return
          end if

          error = .not. toler0>=zero
          if (error) then
              if (text>0) write (text, 106) id, toler0
              return
          end if

          ! Check tolerances in [0,1]
          error = .not. (zero <= toler1 .and. toler1 <= one &
                  .and.  zero <= toler2 .and. toler2 <= one)
          if (error) then
              if (text>0) write (text, 107) id, toler1, toler2
              return
          end if

          ! Check monotonic
          counted = count(x(1:points-1)<x(2:points))
          error = .not. (counted==0 .or. counted==points-1)
          if (error) then
              if (text>0) write (text, 108) id
              return
          end if

          ! at each interval, count the active, significant components that vary too greatly.
          act    = 0  ! number of active components
          signif = 0  ! number of significant components: max(u)-min(u)>=tol*max(|u|)
          mark(:points) = .false.
          ratio1(:points) = zero
          ratio2(:points) = zero
          vary1 (:points) = 0
          vary2 (:points) = 0

          ! top of the loop over the components.
          active_components: do j = 1, comps

             if (.not.active(j)) cycle active_components

             act = act + 1

             ! find range and maximum magnitude of this component.
             lower = u(j, 1)
             upper = u(j, 1)
             do k = 2, points
                lower = min (lower, u(j,k))
                upper = max (upper, u(j,k))
             end do
             range  = upper - lower
             maxmag = max(abs(lower), abs(upper))

             ! decide whether the component is significant.
             if (.not. abs(range)>toler0*max(one,maxmag)) cycle active_components

             ! this is a significant component.
             signif = signif + 1

             ! at each interval, see whether the component'S CHANGE EXCEEDS SOME
             ! fraction of the component'S GLOBAL CHANGE.
             max_du: do k = 1, points - 1
                  differ = abs (u(j,k+1)-u(j,k))
                  if (zero<range) ratio1(k) = max (ratio1(k), differ / range)
                  if (toler1 * range<differ) vary1(k) = vary1(k) + 1
             end do max_du

             ! find the global change of the component'S DERIVATIVE.
             temp  = grad(u,x,comp=j,point=1)
             lower = temp
             upper = temp
             max_grad: do k = 2, points - 1
                 temp  = grad(u,x,comp=j,point=k)
                 lower = min(lower, temp)
                 upper = max(upper, temp)
             end do max_grad
             range = upper - lower

             ! at each interior point, see whether the derivative'S CHANGE
             ! exceeds some fraction of the derivative'S GLOBAL CHANGE.
             right =  grad(u,x,comp=j,point=1)
             do k = 2, points - 1
                 left = right
                 right = grad(u,x,comp=j,point=k)
                 differ = abs (left - right)
                 if (zero<range) ratio2(k) = max (ratio2(k), differ / range)
                 if (toler2 * range < differ) vary2(k) = vary2(k) + 1
             end do

          end do active_components

          ! save the maximum ratios.
          ratio(1) = max(zero,maxval(ratio1(1:points-1),1))
          ratio(2) = max(zero,maxval(ratio2(2:points-1),1))

          ! ***** select the intervals to halve. *****

          ! weight the intervals in which variations that are too large occur.
          most = 0
          amr_intervals: do k = 1, points - 1
             weight(k) = vary1(k)
             if (1<k) weight(k) = weight(k) + vary2(k)
             if (k<points - 1) weight(k) = weight(k) + vary2(k + 1)
             if (0<weight(k)) most = most + 1
          end do amr_intervals

          ! sort the weights using interchange sort.
          do k = 1, points - 1
             do j = k + 1, points - 1
                if (weight(j)>weight(k)) then
                   itemp = weight(j)
                   weight(j) = weight(k)
                   weight(k) = itemp
                end if
             end do
             if (weight(k) == 0) exit
          end do

          ! find the least weight of intervals to halve.
          more = max (0, min (most, padd, pmax - points))
          if (more>0) then
             least = weight(more)
          else
             least = 1 + weight(1)
          end if

          ! reconstruct the weights.
          do k = 1, points - 1
             weight(k) = vary1(k)
             if (k>1)        weight(k) = weight(k) + vary2(k)
             if (k<points-1) weight(k) = weight(k) + vary2(k + 1)
          end do

          ! mark the intervals to halve.
          counted = 0
          to_be_halved: do k = 1, points - 1
             if (counted<more .and. least <= weight(k)) then
                counted = counted + 1
                mark(k) = .true.
             end if
          end do to_be_halved

          ! hack  if one point is marked, mark them all
          !      if (counted>0) then
          !         counted = 0
          !         do k = 1, points - 1
          !            counted = counted + 1
          !            mark(k) = .true.
          !         enddo
          !      endif

          more = counted

          ! ***** halve the intervals, if any. *****

          ! total number of points in the new grid.
          total = points + more

          add_points: if (more>0) then

              counted = 0
              length  = abs(x(points) - x(1))
              check_degenerate: do k = 1, points - 1
                 if (mark(k)) then
                    mean = half*(x(k)+x(k+1))
                    ! Check this interval is not degenerate
                    if (.not. ((x(k)  <mean .and. mean<x(k+1)) .or. &
                               (x(k+1)<mean .and. mean<x(k)))) counted = counted + 1
                 end if
              end do check_degenerate
              error = counted>0
              if (error) then
                  if (text>0) write (text, 109) id
                  return
              end if

              ! add the new points, interpolate x and the bounds.
              new = total
              new_points: do old = points, 2, - 1

                 ! Copy right boundary
                 x(new)   = x(old)
                 u(:,new) = u(:,old)

                 new = new - 1

                 ! Interpolate solution and location
                 if (mark(old-1)) then
                    x(new)   = half*(x(old)+x(old-1))
                    u(:,new) = half*(u(:,old)+u(:,old-1))
                    new = new - 1
                 end if
              end do new_points

              ! mark the new points.
              new = total
              mark_new_points: do old = points, 2, - 1
                 mark(new) = .false.
                 new = new - 1
                 if (mark(old-1)) then
                    mark(new) = .true.
                    new = new - 1
                 end if
              end do mark_new_points
              mark(new) = .false.

              ! update the number of points.
              former = points
              points = total

              ! allow the user to update the solution.
              call twcopy (comps*points, u, buffer)

              ! Request to update the grid
              signal = 'UPDATE'
              ! go to 4060 when route = 1
              route = 1
              return

              ! Restart after grid update.
              ! DO NOT MOVE THIS out of (more>0) block
              4060  continue
              signal = ' '
              call twcopy (comps*points,buffer,u)

          end if add_points

          ! ***** epilogue. *****

          ! print summary
          if (levelm>0 .and. text>0) then

             temp1 = maxval(ratio1(1:former-1),1)
             temp2 = maxval(ratio2(2:former-1),1)

             if (signif == 0) then
                write (text, 2) id
             else

                write (text, 3) temp1, temp2, toler1, toler2
                if (most == 0) then
                    write (text, 4) id
                else if (more == 0) then
                    write (text, 5) id
                else
                    write (text, 6)

                    old = 0
                    do k = 1, points
                       if (.not. mark(k)) then
                          old = old + 1
                          if (1<k) then
                             if (vary1(old-1)/=zero) then
                                write (word, '(F4.2, I4)') ratio1(old - 1), vary1(old-1)
                             else
                                write (word, '(F4.2, I4)') ratio1(old - 1)
                             end if
                             if (mark(k - 1)) then
                                write (text, 7) k - 1, x(k - 1), word
                             else
                                write (text, 8) word
                             end if
                          end if

                          if (k>1 .and. k<points) then
                             if (vary2(old)/=0) then
                                 write (word, '(F4.2, I4)') ratio2(old), vary2(old)
                             else
                                 write (word, '(F4.2, I4)') ratio2(old)
                             end if
                             write (text, 9) k, x(k), word
                          else
                             write (text, 9) k, x(k)
                          end if
                       end if
                    end do
                end if
             end if

             if (leveld>0 .and. more>0) then
                write (text, 10) id
                call twcopy(comps*points,from=u,to=buffer)
                signal = 'SHOW'
                ! go to 5040 when route = 2
                route = 2
                return
             end if
          end if

          ! Restart after SHOW request
          5040 continue

          ! set the completion status flags.
          signal  = ' '
          newx    = more > 0
          success = most == 0

          return

          ! Formats section.
             ! Informative messages.
             1 format(/1X, a9, 'SELECT A GRID.')
             2 format(/1X, a9, 'SUCCESS.  THE GRID IS ADEQUATE BECAUSE ALL ACTIVE' &
                      /10X, 'COMPONENTS ARE INSIGNIFICANT.')
             !               123456789-   1234567   1234567
             3 format(/15X, '             RATIO 1   RATIO 2' &
                      /15X, '             -------   -------' &
                      /15X, '    ACTUAL', 2F10.3 &
                      /15X, '   DESIRED', 2F10.3)
             4 format(/1X, a9, 'SUCCESS.  THE GRID IS ADEQUATE.')
             5 format(/1X, a9, 'FAILURE.  MORE POINTS ARE NEEDED BUT NONE CAN BE ADDED.')
             !               123456   123456789-123456   12345678   12345678
             6 format(/10X, 'THE NEW GRID (* MARKS NEW POINTS):' &
                     //10X, '                             LARGEST RATIOS AND'&
                      /10X, ' INDEX         GRID POINT      NUMBER TOO LARGE'&
                      /10X, '------   ----------------   -------------------'&
                      /10X, '                             RATIO 1    RATIO 2')
             7 format(10X, i6, '*  ', 1p, e16.9, 0p, 3X, a8)
             8 format(38X, a8)
             9 format(10X, i6, '   ', 1p, e16.9, 0p, 14X, a8)
            10 format(/1X, a9, 'THE SOLUTION GUESS FOR THE NEW GRID:')

            ! Error messages.
            101 format(/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
                     //10X, i10, '  ROUTE')
            102 format(/1X, a9, 'ERROR.  THERE MUST BE AT LEAST ONE COMPONENT AND AT' &
                      /10X, 'LEAST TWO POINTS.' &
                     //10X, i10, '  COMPS, COMPONENTS' &
                      /10X, i10, '  POINTS')
            103 format(/1X, a9, 'ERROR.  THE LIMIT ON POINTS ADDED TO A GRID MUST BE' &
                      /10X, 'ZERO OR POSITIVE.'&
                     //10X, i10, '  PADD, LIMIT ON ADDED POINTS')
            104 format(/1X, a9, 'ERROR.  POINTS IS OUT OF RANGE.' &
                     //10X, i10, '  POINTS'&
                      /10X, i10, '  PMAX, LIMIT ON POINTS')
            105 format(/1X, a9, 'ERROR.  THERE ARE NO ACTIVE COMPONENTS.')
            106 format(/1X, a9, 'ERROR.  THE BOUNDS ON MAGNITUDE AND RELATIVE CHANGE' &
                      /10X, 'OF MAGNITUDE FOR INSIGNIFICANT COMPONENTS MUST BE'&
                      /10X, 'POSITIVE.'&
                     //10X, 1p, e10.2, '  TOLER0, SIGNIFICANCE LEVEL')
            107 format(/1X, a9, 'ERROR.  THE BOUNDS ON RELATIVE CHANGES IN MAGNITUDE'&
                      /10X, 'AND ANGLE MUST LIE BETWEEN 0 AND 1.'&
                     //10X, 1p, e10.2, '  TOLER1'&
                      /10X, 1p, e10.2, '  TOLER2')
            108 format(/1X, a9, 'ERROR.  THE GRID IS NOT ORDERED.')
            109 format(/1X, a9, 'ERROR.  SOME INTERVALS IN THE GRID ARE TOO SHORT.'&
                      /10X, 'THE NEW GRID WOULD NOT BE ORDERED.')

      return
      end subroutine refine

      !> Return gradient of the solution vector at given location and point
      pure real(RK) function grad(u,x,comp,point)
          real(RK), intent(in) :: u(:,:),x(:)
          integer , intent(in) :: comp  ! ID of the unknown
          integer , intent(in) :: point ! ID of the finite-difference node (<=points-1)
          grad = (u(comp,point+1)-u(comp,point)) / (x(point+1)-x(point))
      end function grad

      !> Return precision flag based on the current real kind
      function precision_flag()
          character(len=16) :: precision_flag

          if (precision(zero)==precision(0.0)) then
             precision_flag = 'SINGLE PRECISION'
          elseif (precision(zero)==precision(0.d0)) then
             precision_flag = 'DOUBLE PRECISION'
          else
             precision_flag = 'UNKNWN PRECISION'
          endif

      end function precision_flag

      !> Check that the requested version is supported
      subroutine check_version(version_label,text,error)
         character(len=*), intent(in) :: version_label
         integer, intent(in) :: text
         logical, intent(out) :: error

         integer :: j
         logical :: match_found
         character(*), parameter :: id = 'TWOPNT:  '

         ! Check all compatible versions for a match
         match_found = .false.
         do j = 1, vnmbrs
            match_found = version_label == precision_flag()//' VERSION ' // vnmbr(j)
            if (match_found) exit
         end do

         error = .not.match_found

         if (error .and. text>0) then
            write (text, 1) id, trim(version_label),precision_flag(),vnmbr(vnmbrs)
            do j = vnmbrs - 1, 1, - 1
                write (text,'(10X,A,A)')' CAN REPLACE:  '//precision_flag()//' VERSION ',vnmbr(j)
            end do
         end if

         return

         1 format(/1X, a9, 'ERROR.  THE CALLING PROGRAM EXPECTS A VERSION OF' &
                 /10X, 'TWOPNT NOT COMPATIBLE WITH THIS VERSION.' &
                //10X, '     EXPECTS:  ', a &
                //10X, 'THIS VERSION:  ',a,' VERSION ', a)

      end subroutine check_version

      ! Print computer time
      function print_time(seconds) result(string)
          real(RK), intent(in) :: seconds
          character(len=:), allocatable :: string

          character(len=80) :: buffer

          if (hour <= seconds) then
             write (buffer, '(F10.2, A)') seconds / hour, ' HOURS'
          else if (minute <= seconds) then
             write (buffer, '(F10.2, A)') seconds / minute, ' MINUTES'
          else
             write (buffer, '(F10.2, A)') seconds, ' SECONDS'
          end if

          string = trim(buffer)

      end function print_time

      subroutine partition_working_space(this,text,vars,error)
          class(twwork), intent(inout) :: this
          integer,       intent(in)    :: text
          type(twsize),  intent(in)    :: vars
          logical,       intent(out)   :: error

          character(*), parameter :: id = 'TWOPNT:  '
          integer :: N

          N = vars%NMAX()

          call realloc(this%vary,vars%pmax,error);    if (error) goto 1
          call realloc(this%vary1,vars%pmax,error);   if (error) goto 1
          call realloc(this%vary2,vars%pmax,error);   if (error) goto 1
          call realloc(this%above,N,error);           if (error) goto 1
          call realloc(this%below,N,error);           if (error) goto 1
          call realloc(this%ratio1,vars%pmax,error);  if (error) goto 1
          call realloc(this%ratio2,vars%pmax,error);  if (error) goto 1
          call realloc(this%s0,N,error);              if (error) goto 1
          call realloc(this%s1,N,error);              if (error) goto 1
          call realloc(this%usave,N,error);           if (error) goto 1
          call realloc(this%vsave,N,error);           if (error) goto 1
          call realloc(this%v1,N,error);              if (error) goto 1
          call realloc(this%xsave,vars%pmax,error);   if (error) goto 1
          call realloc(this%y0,N,error);              if (error) goto 1
          call realloc(this%y1,N,error);              if (error) goto 1

          ! Success!
          return

          ! Error control
          1 if (text>0) write (text, 10) id; return

          10 format(/1X, a9, 'ERROR.  TWGRAB FAILS.')

      end subroutine partition_working_space

      pure subroutine expand_bounds(this,above,below,vars)
         class(twwork), intent(inout) :: this
         real(RK),      intent(in)    :: above(:),below(:)
         type(twsize),  intent(in)    :: vars

         integer :: ptr,j,k

         ! Input arrays: bounds for groups and components only
         ! Working array: bounds for each point
         ! +--------+-------+--------+
         ! | groupa | comps | groupb |
         ! +--------+-------+--------+

         ! Working array: bounds for each point
         ! +--------+------------------------------------+--------+
         ! | groupa | point1 | point2 | point.. | points | groupb |
         ! |        |          comps*points              |        |
         ! +--------+------------------------------------+--------+

         ! EXPAND THE BOUNDS.
         ptr = 1
         do j = 1, vars%groupa
             this%above(ptr) = above(j)
             this%below(ptr) = below(j)
             ptr = ptr + 1
         end do

         do k = 1, vars%points
             do j = 1, vars%comps
                this%above(ptr) = above(vars%groupa + j)
                this%below(ptr) = below(vars%groupa + j)
                ptr = ptr + 1
             end do
         end do

         do j = 1, vars%groupb
             this%above(ptr) = above(vars%groupa + vars%comps + j)
             this%below(ptr) = below(vars%groupa + vars%comps + j)
             ptr = ptr + 1
         end do

      end subroutine expand_bounds

      pure subroutine realloc_real(array,min_size,error)
         real(RK), allocatable, intent(inout) :: array(:)
         integer , intent(in)  :: min_size
         real(RK), allocatable :: tmp(:)
         logical , intent(out) :: error
         integer :: stat

         if (allocated(array)) then
             if (size(array)<min_size) then
                allocate(tmp(min_size),stat=stat)
                call move_alloc(from=tmp,to=array)
             else
                stat = 0
             end if
         else
             allocate(array(min_size),stat=stat)
         end if

         error = stat/=0

      end subroutine realloc_real

      pure subroutine realloc_int(array,min_size,error)
         integer , allocatable, intent(inout) :: array(:)
         integer , intent(in)  :: min_size
         integer , allocatable :: tmp(:)
         logical , intent(out) :: error
         integer :: stat

         if (allocated(array)) then
             if (size(array)<min_size) then
                allocate(tmp(min_size),stat=stat)
                call move_alloc(from=tmp,to=array)
             else
                stat = 0
             end if
         else
             allocate(array(min_size),stat=stat)
         end if

         error = stat/=0

      end subroutine realloc_int

      subroutine init_stats(this,points)
         class(twstat), intent(inout) :: this
         integer      , intent(in)    :: points

          ! Init statistics arrays
          this%total = zero
          this%detail = zero
          this%event = 0
          this%gsize = 0

          ! Init total time counter
          call twtime(this%timer(qtotal))

          ! Initialize first grid
          this%grid = 1
          this%gsize(this%grid) = points
          call twtime(this%timer(qgrid))

          ! Initialize time step number
          this%step = 0

      end subroutine init_stats

      subroutine tock_stats(this,task,event)
          class(twstat), intent(inout) :: this
          integer, intent(in) :: task
          logical, optional, intent(in) :: event

          call twlaps(this%timer(task))
          this%total(task) = this%total(task) + this%timer(task)
          if (this%grid <= gmax) this%detail(this%grid, task) = this%detail(this%grid, task) &
                                                              + this%timer(task)
          ! Task-specific finalization
          select case (task)
             case (qgrid)
                 this%detail(this%grid, qother) = this%detail(this%grid,task) &
                                                - sum(this%detail(this%grid,[qfunct,qjacob,qsolve]))
             case (qtotal)
                 this%total(qother) = this%total(task)-sum(this%total([qfunct,qjacob,qsolve]))
          end select

          if (this%grid<=gmax .and. present(event)) then
              if (event) this%event(this%grid, task) = this%event(this%grid, task) + 1
          end if

      end subroutine tock_stats

      subroutine tick_stats(this,task)
         class(twstat), intent(inout) :: this
         integer, intent(in) :: task
         call twtime(this%timer(task))
      end subroutine tick_stats

      subroutine stats_new_grid(this,points)
         class(twstat), intent(inout) :: this
         integer, intent(in) :: points

         ! Complete statistics for the last grid
         call this%tock(qgrid)

         ! Initialize statistics for the new grid
         this%grid = this%grid + 1
         if (this%grid <= gmax) then
            call twtime (this%timer(qgrid))
            this%gsize(this%grid) = points
         end if

      end subroutine stats_new_grid

      ! Identify the request
      integer function identify_request(signal) result(qtype)
         character(len=*), intent(in) :: signal

         select case (signal)
            case ('RESIDUAL'); qtype = qfunct
            case ('PREPARE');  qtype = qjacob
            case ('SOLVE');    qtype = qsolve
            case default;      qtype = qother
         end select

      end function identify_request

      ! Print the identifier for this variable's unknown class
      elemental subroutine unknown_class_id(i,groupa,comps,clss,localid)
          integer, intent(in) :: i,groupa,comps
          character, intent(out) :: clss
          integer  , intent(out) :: localid
          if (i<=groupa) then
             clss    = 'A'     ! GROUP A UNKNOWNS
             localid = i
          elseif (i<=groupa+comps) then
             clss    = 'C'     ! COMPONENTS AT POINTS
             localid = i-groupa
          else
             clss    = 'B'     ! GROUP B UNKNOWNS
             localid = i-groupa-comps
          end if
      end subroutine unknown_class_id

      ! Print the identifier for this variable's unknown class
      pure function unknown_class_name(i,groupa,comps) result(s)
          integer, intent(in) :: i,groupa,comps
          character(len=:), allocatable :: s
          if (i<=groupa) then
             s = 'GROUP A UNKNOWNS'
          elseif (i<=groupa+comps) then
             s = 'COMPONENTS AT POINTS'
          else
             s = 'GROUP B UNKNOWNS'
          end if
      end function unknown_class_name

      ! Detailed error message for a case with invalid variable bounds
      subroutine print_invalid_bounds(id,text,name,vars,below,above)
         integer     , intent(in) :: text
         type(twsize), intent(in) :: vars
         real(RK),     intent(in) :: below(vars%N())
         real(RK),     intent(in) :: above(vars%N())
         character(len=*), intent(in) :: id,name(:)

         integer :: j,counter,len1,len2,local,length
         intrinsic :: count
         character(len=80) :: ctemp1,ctemp2,string
         character :: c

         if (text==0) return

         counter = count(.not.below<above)

         ! Header
         write (text,1) id, vars%groupa, vars%groupb, vars%comps, vars%NVAR(), counter

         counter = 0
         loop_bounds: do j = 1, vars%NVAR()

            if (below(j) < above(j)) cycle loop_bounds

            counter = counter + 1
            if (counter > MAX_ERROR_LINES) exit loop_bounds

            ! Use variable names, if provided
            if (size(name) == vars%NVAR()) then
                ctemp1 = name(j)
            else
                ctemp1 = ' '
            end if
            call twsqez(len1,ctemp1)

            call unknown_class_id(j,vars%groupa,vars%comps,c,local)
            write(ctemp2,2) c,local
            call twsqez (len2, ctemp2)

            ! Limit to 40 columns
            if (ctemp1 == ' ') then
               string = ctemp2
               length = len2
            else if (len1 + 2 + len2 <= CONTRL_MAX_LEN) then
               string = ctemp1 (1 : len1) // '  ' // ctemp2
               length = len1 + 2 + len2
            else if (len1 + 1 + len2 <= CONTRL_MAX_LEN) then
               string = ctemp1 (1 : len1) // ' ' // ctemp2
               length = len1 + 1 + len2
            else
               len1   = CONTRL_MAX_LEN - len2 - 4
               string = ctemp1 (1 : len1) // '... ' // ctemp2
               length = CONTRL_MAX_LEN
            end if

            write (text, 3) below(j), above(j), string(1:length)

         end do loop_bounds

         if (MAX_ERROR_LINES<counter) write (text, 4)

         ! Formats section
         1 format(/1X, a9, 'ERROR.  THE LOWER AND UPPER BOUNDS ON SOME UNKNOWNS' &
                 /10X, 'ARE OUT OF ORDER.' &
                //10X, i10, '  GROUP A UNKNOWNS (A)' &
                 /10X, i10, '  GROUP B UNKNOWNS (B)' &
                 /10X, i10, '  COMPONENTS AT POINTS (C)' &
                 /10X, i10, '  TOTAL TYPES OF UNKNOWNS' &
                 /10X, i10, '  NUMBER OF BOUNDS OUT OF ORDER' &
                //10X, '     LOWER       UPPER' &
                 /10X, '     BOUND       BOUND   UNKNOWN'/)

         2 format('(',a,1x,i10,')')
         3 format(10X, 1p, e10.2, 2X, e10.2, 3X, a)
         4 format(10X, '  ... MORE')

      end subroutine print_invalid_bounds

      subroutine print_invalid_ranges(id,text,name,vars,below,above,v0,s0)
         integer     , intent(in) :: text
         type(twsize), intent(in) :: vars
         real(RK),     intent(in) :: below (vars%N())
         real(RK),     intent(in) :: above (vars%N())
         real(RK),     intent(in) :: v0    (vars%N())
         real(RK),     intent(in), optional :: s0(vars%N())
         character(len=*), intent(in) :: id,name(:)

         integer :: i,j,counter,len1,len2,length
         logical :: search,verified(vars%N())
         character(len=80) :: ctemp1,ctemp2,string

         if (text==0) return

         search = present(s0)

         if (search) then
             verified = (below==v0.and.s0>zero).or.(above==v0.and.s0<zero)
         else
             verified = (below<=v0.and.v0<=above)
         end if

         counter = count(.not.verified)

         if (search) then
            write (text,7) id
         else
            write (text,1) id, vars%groupa,vars%groupb,vars%comps,vars%points,vars%N(),counter
         end if

         counter = 0
         loop_bounds: do j = 1, vars%N()

            if (verified(j)) cycle loop_bounds

            counter = counter + 1
            if (counter > MAX_ERROR_LINES) exit loop_bounds

            if (j <= vars%groupa) then
               i = j
            else if (j <= vars%groupa + vars%comps*vars%points) then
               i = vars%groupa + mod (j - vars%groupa - 1, vars%comps) + 1
            else
               i = j - vars%groupa - vars%comps*vars%points
            end if

            ! Use variable names, if provided
            if (size(name) == vars%NVAR()) then
                ctemp1 = name(i)
            else
                ctemp1 = ' '
            end if
            call twsqez(len1,ctemp1)

            if (j <= vars%groupa) then
                write (ctemp2, 2) 'A', i
            else if (j <= vars%groupa + vars%comps*vars%points) then
                write (ctemp2, 3) 'C', i, 'P', int((j-vars%groupa-1)/vars%comps) + 1
            else
                write (ctemp2, 2) 'B', i
            end if
            call twsqez (len2, ctemp2)

            ! Limit to 40 columns
            if (ctemp1 == ' ') then
               string = ctemp2
               length = len2
            else if (len1 + 2 + len2 <= 30) then
               string = ctemp1 (1 : len1) // '  ' // ctemp2
               length = len1 + 2 + len2
            else if (len1 + 1 + len2 <= 30) then
               string = ctemp1 (1 : len1) // ' ' // ctemp2
               length = len1 + 1 + len2
            else
               len1   = 30 - len2 - 4
               string = ctemp1 (1 : len1) // '... ' // ctemp2
               length = 30
            end if

            if (search) then
                write (text, 6) merge('LOWER','UPPER',below(j)==v0(j)),v0(j),string(1:length)
            else
                write (text, 5) below(j), v0(j), above(j), string(1:length)
            end if

         end do loop_bounds

         if (MAX_ERROR_LINES<counter) write (text, 4)
         return

         ! Formats section
         1 format(/1X, a9, 'ERROR.  THE GUESSES FOR SOME UNKNOWNS ARE OUT OF' &
                 /10X, 'BOUNDS.' &
                //10X, i10, '  GROUP A UNKNOWNS (A)' &
                 /10X, i10, '  GROUP B UNKNOWNS (B)' &
                 /10X, i10, '  COMPONENTS AT POINTS (C)' &
                 /10X, i10, '  POINTS (P)' &
                 /10X, i10, '  TOTAL UNKNOWNS' &
                 /10X, i10, '  NUMBER OUT OF BOUNDS' &
                //10X, '     LOWER                   UPPER' &
                 /10X, '     BOUND       VALUE       BOUND   UNKNOWN'/)
         2 format('(',a,1x,i10,')')
         3 format('(',a,1x,i10,1x,a,1x,i10,')')
         4 format(10X, '  ... MORE')
         5 format(10X, 1p, e10.2, 2X, e10.2, 2X, e10.2, 3X, a)
         6 format(10X, a5, 2X, 1p, e10.2, 3X, a)

         7 format(/1X, a9, 'FAILURE.  THE SEARCH FOR THE FOLLOWING UNKNOWNS GOES' &
                 /10X, 'OUT OF BOUNDS.' &
                //10X, 'BOUND       VALUE   UNKNOWN'/)

      end subroutine print_invalid_ranges


      subroutine print_stats(this,text,adapt)
         class(twstat), intent(inout) :: this
         integer, intent(in) :: text
         logical, intent(in) :: adapt

         integer :: j,k
         real(RK) :: temp
         character(len=80) :: header(6)
         character(*), parameter :: id = 'TWOPNT: '

         associate(total=>this%total,detail=>this%detail,gsize=>this%gsize,grid=>this%grid,&
                   event=>this%event)

         time_output: if (total(qtotal)>zero) then

             ! Report total computer time
             write (text, 10004) id, print_time(total(qtotal))

             ! REPORT PERCENT OF TOTAL COMPUTER TIME.
             temp = hundred/total(qtotal)

             adaptive_grid: if (adapt) then

                !                  123456789_123456789_123456789_12345678
                !                  123456  123456  123456 123456 123456
                      header(1) = '                TASK                  '
                      header(3) = '  GRID    GRID  --------------------  '
                      header(5) = 'POINTS  TOTALS  EVOLVE SEARCH REFINE  '

                !                  123456789_123456789_1234567
                !                  123456 123456 123456 123456
                      header(2) = 'SUBTASK                    '
                      header(4) = '---------------------------'
                      header(6) = 'EVAL F PREP J  SOLVE  OTHER'

                      write (text, 10005) header,(gsize(j),(temp*detail(j,k), k=1,8), j=1,grid)
                      if (grid>1)    write (text, 10006) (temp * total(k), k = 2, 8)
                      if (grid>gmax) write (text, 10007)

                ! REPORT AVERAGE COMPUTER TIME.

                !                  123456789_123456789_123456789_1234567
                !                  123456   1234567  1234567  1234567
                      header(1) = '         AVERAGE SECONDS             '
                      header(3) = '  GRID   -------------------------   '
                      header(5) = 'POINTS    EVAL F   PREP J    SOLVE   '


                !                  123456789_123456789_12345
                !                  1234567  1234567  1234567
                      header(2) = 'NUMBER OF SUBTASKS       '
                      header(4) = '-------------------------'
                      header(6) = ' EVAL F   PREP J    SOLVE'

                      write (text, 10009) header, &
                                          (gsize(j), (detail(j,k) / event(j, k), k = 5,7), &
                                          (event(j, k), k = 5,7), j = 1, grid)

             else adaptive_grid

                !                  123456789_123456789_123456789_123456789_123456789_1
                !                  123456   123456   123456   123456   123456   123456
                      header(1) = 'SUBTASK                             TASK           '
                      header(2) = '---------------------------------   ---------------'
                      header(3) = 'EVAL F   PREP J    SOLVE    OTHER   EVOLVE   SEARCH'

                      write (text, 10008) (header(j), j=1,3), '  % OF TOTAL', &
                         (temp * total(k), k = 5, 8), (temp * total(k), k = 2, 3), &
                         'MEAN SECONDS', (detail(1, k) / event(1, k), k = 5, 7), &
                         '    QUANTITY', (event(1, k), k = 5, 7)

             end if adaptive_grid

         end if time_output

         endassociate

         10004 format(/1X, a9, a, ' TOTAL COMPUTER TIME (SEE BREAKDOWN BELOW).')
         10005 format(/10X, 'PERCENT OF TOTAL COMPUTER TIME FOR VARIOUS TASKS:' &
                      /3(/10X, a38, a27) &
                     //(10X, i6, 2X, f6.1, 1X, 3(1X, f6.1), 1X, 4(1X, f6.1)))
         10006 format(/12X, 'TASK TOTALS:', 1X, 3(1X, f6.1), 1X, 4(1X, f6.1))
         10007 format(/10X, 'SOME GRIDS ARE OMITTED, BUT THE TOTALS ARE FOR ALL.')
         10008 format(3(/24X, a51) &
                     //10X, a12, f8.1, 5F9.1 &
                      /10X, a12, f8.3, 2F9.3 &
                      /10X, a12, i8, 2i9)
         10009 format(/10X, 'AVERAGE COMPUTER TIMES FOR, AND NUMBERS OF, SUBTASKS:' &
                     /3(/10X, a37, a25) &
                     //(10X, i6, 3X, f7.3, 2X, f7.3, 2X, f7.3, 1X, 3(2X, i7)))


      end subroutine print_stats

      ! Decision block: the previous task determines the next
      integer function twopnt_next_task(setup,stats,old_task,found,satisf,allow,report,desire,error) result(qtask)
          type(twcom),  intent(in)    :: setup
          type(twstat), intent(inout)  :: stats
          integer,      intent(in)     :: old_task
          logical,      intent(in)     :: found,satisf
          logical,      intent(inout)  :: allow
          character(*), intent(inout)  :: report
          integer,      intent(inout)  :: desire
          logical,      intent(out)    :: error

          error = .false.
          qtask = old_task

          select case (old_task)

             case (qentry)
                 if (setup%steps0>0) then
                    qtask  = qtimst
                    desire = setup%steps0
                 else if (setup%steady) then
                    qtask  = qsearc
                 else
                    error  = .true.
                    return
                 end if

             case (qsearc) ! SEARCH WAS THE PREVIOUS TASK.

                 if (found) then
                    if (setup%adapt) then
                       qtask = qrefin
                    else
                       qtask  = qexit
                       report = ' '
                    end if
                 else
                    if (allow .and. setup%steps1>0) then
                       qtask  = qtimst
                       desire = setup%steps1
                    else
                       qtask = qexit
                       if (stats%grid>1) then
                          report = 'SOME SOLVED'
                       else
                          report = 'NOT SOLVED'
                       end if
                    end if
                 end if

             case (qrefin) ! REFINE WAS THE PREVIOUS TASK.

                 if (found) then
                    stats%step = 0
                    qtask = qsearc
                    allow = .true.
                 else
                    qtask = qexit
                    if (satisf) then
                       report = ' '
                    else
                       report = 'NO SPACE'
                    end if
                 end if

             case (qtimst) ! EVOLVE WAS THE PREVIOUS TASK.

                 if (found) then
                    if (setup%steady) then
                       qtask = qsearc
                    else
                       qtask = qexit
                       report = ' '
                    end if
                 else
                    qtask = qexit
                    if (stats%grid>1) then
                       report = 'SOME SOLVED'
                    else
                       report = 'NOT SOLVED'
                    end if
                 end if

          end select

      end function twopnt_next_task

end module twopnt_core

