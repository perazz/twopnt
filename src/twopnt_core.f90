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
!
! *************************************************************************************************
module twopnt_core
    use iso_fortran_env, only: real64
    implicit none
    private

    public :: twlast,twcopy,twcom,twsolv,twshow,twprep,twopnt

    integer, parameter, public :: RK = real64

    ! Numeric constants
    real(RK), parameter, public :: zero = 0.0_RK
    real(RK), parameter, public :: one  = 1.0_RK
    real(RK), parameter, public :: pi   = acos(-1.0_RK)

    ! Machine epsilon and the absolute and relative perturbations.
    real(RK), parameter, public :: eps   = epsilon(0.0_RK)
    real(RK), parameter, public :: absol = sqrt(eps)
    real(RK), parameter, public :: relat = sqrt(eps)

    ! Error codes
    integer,  parameter, public :: qnull = 0
    integer,  parameter, public :: qbnds = 1
    integer,  parameter, public :: qdvrg = 2

    ! LOCATION OF DATA IN ARRAYS DETAIL, EVENT, TIMER, AND TOTAL.  THE
!     LOCATIONS ARE CHOSEN TO SIMPLIFY WRITE STATEMENTS.  DETAIL USES
!     ONLY 1 : 8, EVENT USES ONLY 5 : 8, TIMER USES 1 : 9, AND TOTAL
!     USES ONLY 2 : 9.  IN ADDITION, 2, 3, 4, 10, AND 11 ARE USED AS
!     MNEMONIC VALUES FOR QTASK.

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

    logical, parameter :: DEBUG = .true.
    integer, parameter :: CONTRL_MAX_LEN = 40

    ! Settings structure (formerly a common block)
    integer, parameter :: cntrls = 22

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

    contains

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

      subroutine twshow(error, text, buffer, comps, grid, groupa, groupb, points, x)

          logical, intent(out) :: error
          integer, intent(in)  :: text,comps,groupa,groupb,points
          logical, intent(in)  :: grid
          real(RK), intent(in) :: buffer(groupa + comps * points + groupb), x(*)

          ! Local variables
          character(len=80) :: string, title(6)
          integer :: cols, comp, count, first, groups, j, last, length, point
          intrinsic :: min
          character(*), parameter :: id = 'TWSHOW:  '

          ! SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
          logical, parameter :: mess = .false.

          !///////////////////////////////////////////////////////////////////////
          !
          !     (1) PROLOGUE.
          !
          !///////////////////////////////////////////////////////////////////////

          !///  WRITE ALL MESSAGES.
          if (mess .and. 0 < text) then
              write (text, 1) id, comps, points, groupa, groupb, groupa + comps * points + groupb
              stop
          end if

          !///  CHECK THE ARGUMENTS.
          error = .not. (((0 < comps) .eqv. (0 < points)) .and. &
                           0 <= comps .and. 0 <= points .and. 0 <= groupa .and. &
                           0 <= groupb .and. 0 < groupa + comps * points + groupb)
          if (error) then
              if (0 < text) write (text, 1) id, comps, points, groupa, groupb, &
                                            groupa + comps * points + groupb
              return
          end if

          !///  COUNT THE GROUPS.
          groups = 0
          if (0 < groupa) groups = groups + 1
          if (0 < groupb) groups = groups + 1
          if (0 < comps .and. 0 < points) groups = groups + 1

          !///  CHOOSE NUMBER OF DATA COLUMNS.
          cols = merge(5,6,grid)

          !///////////////////////////////////////////////////////////////////////
          !
          !
          !
          !///////////////////////////////////////////////////////////////////////
          print_data: if (0 < text) then

              ! (2) PRINT THE GROUPED DATA.
              if (0 < groupa) then
                 if (1 < groups) write (text, 11) 'GROUP A UNKNOWNS'
                 write (text, 12) (j, buffer(j), j = 1, groupa)
              end if

              if (0 < groupb) then
                 if (1 < groups) write (text, 11) 'GROUP B UNKNOWNS'
                 write (text, 12) &
                    (j, buffer(groupa + comps * points + j), j = 1, groupb)
              end if

              ! (2) PRINT THE COMPONENTS AT POINTS.
              if (0 < comps .and. 0 < points) then

                 if (1 < groups) write (text, 11) 'COMPONENTS AT POINTS'

                 components: do first = 1, comps, cols
                    count = 0
                    last = min (first + cols - 1, comps)
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
                          write (text, 14) (point, x(point), &
                             (buffer(groupa + comp + comps * (point - 1)), &
                             comp = first, last), point = 1, points)
                       else
                          write (text, 15) (point, &
                             (buffer(groupa + comp + comps * (point - 1)), &
                             comp = first, last), point = 1, points)
                       end if
                    else
                       do point = 1, points
                          if (grid) then
                             write (text, 14) point, x(point), &
                                (buffer(groupa + comp + comps * (point - 1)), &
                                comp = first, last)
                          else
                             write (text, 15) point, &
                                (buffer(groupa + comp + comps * (point - 1)), &
                                comp = first, last)
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

          ! ERROR MESSAGES.
          1 format(/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE' &
                  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES' &
                  /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS' &
                  /10X, 'MUST BE POSITIVE.' &
                 //10X, i10, '  COMPS, COMPONENTS' &
                  /10X, i10, '  POINTS' &
                  /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                  /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                  /10X, i10, '  TOTAL UNKNOWNS')

      end subroutine twshow

      ! COPY ONE VECTOR TO ANOTHER.
      pure subroutine twcopy (n, x, y)
          integer , intent(in)  :: n
          real(RK), intent(in)  :: x(n)
          real(RK), intent(out) ::  y(n)
          y = x
          return
      end subroutine twcopy

      ! Compute the max-norm of a vector
      pure subroutine twnorm (n, value, x)
          integer,  intent(in) :: n
          real(RK), intent(in) :: x(n)
          real(RK), intent(out) :: value

          intrinsic :: abs, max

          value = zero
          if (n>0) value = maxval(abs(x),1)

      end subroutine twnorm

      ! SOLVE A SYSTEM OF LINEAR EQUATIONS USING THE MATRIX PREPARED BY TWPREP.
      subroutine twsolv(error, text, a, asize, buffer, comps, groupa, groupb, pivot, points)

          integer , intent(in) :: asize,text,comps,groupa,groupb,points
          integer , intent(in) :: pivot(groupa + comps * points + groupb)
          double precision, intent(in) :: a(asize)
          double precision, intent(inout) :: buffer(groupa + comps * points + groupb)
          logical , intent(out) :: error

          ! Local variables
          integer   :: n,width
          intrinsic :: max
          character(*), parameter :: id = 'TWSOLV:  '

          ! SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
          logical, parameter :: mess = .false.

          !***** (1) PROLOGUE *****

          ! WRITE MESSAGES only.
          if (mess .and. text>0) then
              write (text, 1) id, comps, points, groupa, groupb, n
              write (text, 2) id, comps, points, groupa, groupb, n, width, (3*width+2)*n, asize
              stop
          end if

          ! CHECK THE ARGUMENTS.
          n = groupa + comps * points + groupb
          error = .not. (((0 < comps) .eqv. (0 < points)) .and. &
                           0 <= comps .and. 0 <= points .and. 0 <= groupa .and. &
                           0 <= groupb .and. 0 < n)
          if (error) then
              if (text>0) write (text, 1) id, comps, points, groupa, groupb, n
              return
          end if

          width = comps + max (comps, groupa, groupb) - 1
          error = .not. ((3 * width + 2) * n <= asize)
          if (error) then
              if (text>0) write (text,2) id, comps, points, groupa, groupb, n, width, &
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
          double precision, intent(in) :: abd(lda,*)
          double precision, intent(inout) :: b(*)
          integer , intent(in) :: pivot(*)

          double precision  :: t
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
      subroutine twprep(error, text, a, asize, buffer, comps, condit, groupa, groupb, pivot, points, return_call)

          integer,  intent(in)    :: text ! output unit
          integer,  intent(in)    :: asize,groupa,groupb,comps,points
          integer,  intent(inout) :: pivot (groupa + comps * points + groupb)
          real(RK), intent(inout) :: buffer(groupa + comps * points + groupb)
          real(RK), intent(inout) :: a(asize), condit
          logical , intent(inout) :: return_call

          real(RK) :: delta, temp
          integer :: block, blocks, cfirst, clast, col, count, diag, j, lda, length, n, offset, &
                     rfirst, rlast, route, row, skip, width
          intrinsic :: abs, int, max, min, mod, sqrt
          logical :: error, found
          character(len=80) :: string

          ! Parameters
          integer, parameter :: lines = 20
          character(len=*), parameter :: id = 'TWPREP:  '

          ! SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
          logical, parameter :: mess = .false.

          save

          ! ***** (1) PROLOGUE *****

          ! Every-time initialization

          ! IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.
          if (return_call) then
             return_call = .false.
             go to (2030, 3050) route
             error = .true.
             go to 9001
          endif

          ! CHECK THE ARGUMENTS.
          n = groupa + comps * points + groupb
          error = .not. (((0 < comps) .eqv. (0 < points)) .and. &
                  0 <= comps .and. 0 <= points .and. 0 <= groupa .and. &
                  0 <= groupb .and. 0 < n)
          if (error) go to 9002

          width = comps + max (comps, groupa, groupb) - 1
          error = .not. ((3 * width + 2) * n <= asize)
          if (error) go to 9003

          ! WRITE ALL MESSAGES.
          if (mess .and. 0 < text) then
             route = 0
             go to 9001
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
          if (0 < groupa) then
             blocks = blocks + 1
             pivot(blocks) = groupa
          end if

          do j = 1, points
             blocks = blocks + 1
             pivot(blocks) = comps
          end do

          if (0 < groupb) then
             blocks = blocks + 1
             pivot(blocks) = groupb
          end if

          ! ***** (2) INITIALIZE THE COLUMNS OF THE MATRIX *****

          ! Store evaluation vector
          a(:n) = buffer(:n)

          ! Clear matrix
          a(n+1:n+lda*n) = zero

          ! EVALUATE THE FUNCTION AT THE UNPERTURBED X.

    !     GO TO 2030 WHEN ROUTE = 1
          route = 1
          return_call = .true.
          go to 99999
    2030  continue

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
                      if (block == 1 .and. 0 < groupa) then
                         cfirst = cfirst + groupa
                      else if (block == blocks .and. 0 < groupb) then
                         cfirst = cfirst + groupb
                      else
                         cfirst = cfirst + comps
                      end if
                      block = block + 1
                   end do
              end do perturb_vector

              if (.not. found) exit column_groups

              ! EVALUATE THE FUNCTION AT THE PERTURBED VALUES.
              ! GO TO 3050 WHEN ROUTE = 2
              route = 2
              return_call = .true.
              go to 99999
              3050  continue

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

                    if (block == 1 .and. 0 < groupa) then
                       clast = cfirst + groupa - 1
                    else if (block == blocks .and. 0 < groupb) then
                       clast = cfirst + groupb - 1
                    else
                       clast = cfirst + comps - 1
                    end if

                    if (1 < block) then
                       if (block == 2 .and. 0 < groupa) then
                          rfirst = cfirst - groupa
                       else
                          rfirst = cfirst - comps
                       end if
                    else
                       rfirst = cfirst
                    end if

                    if (block < blocks) then
                       if (block == blocks - 1 .and. 0 < groupb) then
                          rlast = clast + groupb
                       else
                          rlast = clast + comps
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
                    if (block == 1 .and. 0 < groupa) then
                       cfirst = cfirst + groupa
                    else if (block == blocks .and. 0 < groupb) then
                       cfirst = cfirst + groupb
                    else
                       cfirst = cfirst + comps
                    end if
                    block = block + 1
                 end do

              end do form_columns

          end do column_groups

          ! ***** (4) CHECK FOR ZERO COLUMNS. *****
          call count_zero_columns(a,n,diag,lda,width,count)
          error = .not. (count == 0)
          if (error) go to 9004

          ! ***** (5) SCALE THE ROWS. *****
          call scale_rows(a,n,diag,lda,width,count)
          error = .not. (count == 0)
          if (error) go to 9005

          ! ***** (6) FACTOR THE MATRIX.
          call twgbco(a(n+1), lda, n, width, width, pivot, condit, buffer)
          error = condit == zero
          if (error) go to 9006
          condit = one/condit

    !///////////////////////////////////////////////////////////////////////
    !
    !     INFORMATIVE MESSAGES.
    !
    !///////////////////////////////////////////////////////////////////////

    80001 format(10X, a)
    80002 format(10X, '... MORE')

    !///////////////////////////////////////////////////////////////////////
    !
    !     ERROR MESSAGES.
    !
    !///////////////////////////////////////////////////////////////////////

          go to 99999

    9001  if (0 < text) write (text, 99001) id, route
          if (.not. mess) go to 99999

    9002  if (0 < text) write (text, 99002) id, &
             comps, points, groupa, groupb, n
          if (.not. mess) go to 99999

    9003  if (0 < text) write (text, 99003) id, &
             comps, points, groupa, groupb, n, width, &
             (3 * width + 2) * n, asize
          if (.not. mess) go to 99999

    9004  if (0 < text) then
             write (text, 99004) id, comps, points, groupa, groupb, &
                groupa + comps * points + groupb, count
             count = 0
             do 8010 j = 1, groupa + comps * points + groupb
                if (a(j) == 0.0 .or. mess) then
                   count = count + 1
                   if (count <= lines) then
                      if (j <= groupa) then
                         write (string, '(A, I10)') 'GROUP A ', j
                      else if (j <= groupa + comps * points) then
                         write (string, '(A, I10, A, I10)') &
                            ' COMPONENT ', mod (j - groupa - 1, comps) + 1, &
                            ' AT POINT ', int ((j - groupa - 1) / comps) + 1
                      else
                         write (string, '(A, I10)') &
                            'GROUP B ', j - groupa - comps * points
                      end if
                      call twsqez (length, string)
                      write (text, 80001) string (1 : length)
                   end if
                end if
    8010     continue
             if (lines < count) write (text, 80002)
          end if
          if (.not. mess) go to 99999

    9005  if (0 < text) then
             write (text, 99005) id, comps, points, groupa, groupb, &
                groupa + comps * points + groupb, count
             count = 0
             do 8020 j = 1, groupa + comps * points + groupb
                if (a(j) == 0.0 .or. mess) then
                   count = count + 1
                   if (count <= lines) then
                      if (j <= groupa) then
                         write (string, '(A, I10)') 'GROUP A ', j
                      else if (j <= groupa + comps * points) then
                         write (string, '(A, I10, A, I10)') &
                            ' COMPONENT ', mod (j - groupa - 1, comps) + 1, &
                            ' AT POINT ', int ((j - groupa - 1) / comps) + 1
                      else
                         write (string, '(A, I10)') &
                            'GROUP B ', j - groupa - comps * points
                      end if
                      call twsqez (length, string)
                      write (text, 80001) string (1 : length)
                   end if
                end if
    8020     continue
             if (lines < count) write (text, 80002)
          end if
          if (.not. mess) go to 99999

        9006  if (0 < text) write (text, 99006) id
              if (.not. mess) go to 99999

        ! Formats section
        99001 format(/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
                   //10X, i10, '  ROUTE')

        99002 format(/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE' &
                    /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES' &
                    /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS' &
                    /10X, 'MUST BE POSITIVE.' &
                   //10X, i10, '  COMPS, COMPONENTS' &
                    /10X, i10, '  POINTS' &
                    /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                    /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                    /10X, i10, '  TOTAL UNKNOWNS')

        99003 format(/1X, a9, 'ERROR.  THE MATRIX SPACE IS TOO SMALL.' &
                   //10X, i10, '  COMPS, COMPONENTS' &
                    /10X, i10, '  POINTS' &
                    /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                    /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                    /10X, i10, '  MATRIX ORDER' &
                    /10X, i10, '  STRICT HALF BANDWIDTH' &
                   //10X, i10, '  SPACE REQUIRED' &
                    /10X, i10, '  ASIZE, PROVIDED')

        99004 format(/1X, a9, 'ERROR.  SOME COLUMNS ARE ZERO.' &
                   //10X, i10, '  COMPS, COMPONENTS' &
                    /10X, i10, '  POINTS' &
                    /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                    /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                    /10X, i10, '  TOTAL COLUMNS' &
                    /10X, i10, '  ZERO COLUMNS' &
                   //10X, 'UNKNOWNS WITH ZERO COLUMNS:'/)

        99005 format(/1X, a9, 'ERROR.  SOME ROWS ARE ZERO.' &
                   //10X, i10, '  COMPS, COMPONENTS' &
                    /10X, i10, '  POINTS' &
                    /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                    /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                    /10X, i10, '  TOTAL ROWS' &
                    /10X, i10, '  ZERO ROWS' &
                   //10X, 'ZERO ROWS:'/)

        99006 format(/1X, a9, 'ERROR.  THE JACOBIAN MATRIX IS SINGULAR.')

        stop
        99999 continue
        return
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

      ! Reserve space in an array
      subroutine twgrab (error, last, first, number)
          integer, intent(inout) :: first,last
          integer, intent(in)    :: number
          logical, intent(out)   :: error

          intrinsic :: max

          ! Check the arguments.
          error = .not. (0 <= last);   if (error) return
          error = .not. (0 <= number); if (error) return

          ! Grab the space.
          first = last + 1
          last  = last + max (1, number)

      end subroutine twgrab


      subroutine evolve &
        (error, text, &
         above, below, buffer, comps, condit, desire, groupa, groupb, &
         leveld, levelm, name, names, points, report, s0, s1, signal, &
         step, steps2, strid0, stride, succes, tdabs, tdage, tdec, &
         tdrel, time, tinc, tmax, tmin, v0, v1, vsave, y0, y1, ynorm)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     EVOLVE
!
!     PERFORM TIME EVOLUTION.
!
!///////////////////////////////////////////////////////////////////////

      !implicit complex (a - p, r - z), integer (q)
      character &
         cword*80, header*80, id*9, jword*80, name*(*), remark*80, &
         signal*(*), yword*80
!**** PRECISION > DOUBLE
      double precision    above, below, buffer, change, condit, csave, dummy, high, low, &
         s0, s1, strid0, stride, tdabs, tdec, tdrel, tinc, tmax, tmin, &
         v0, v1, vsave, y0, y1, ynorm
      integer &
         age, agej, comps, count, desire, first, groupa, groupb, j, &
         last, length, leveld, levelm, names, number, points, qbnds, &
         qdvrg, qnull, report, route, step, steps2, tdage, text, &
         xrepor
      intrinsic &
         log10, max, min
      logical &
         error, exist, jacob, mess, succes, time, xsucce

      parameter (id = 'EVOLVE:  ')

!     REPORT CODES
      parameter (qnull = 0, qbnds = 1, qdvrg = 2)

      dimension &
         above(groupa + comps * points + groupb), &
         below(groupa + comps * points + groupb), &
         buffer(groupa + comps * points + groupb), header(2, 3), &
         name(names), &
         s0(groupa + comps * points + groupb), &
         s1(groupa + comps * points + groupb), &
         v0(groupa + comps * points + groupb), &
         v1(groupa + comps * points + groupb), &
         vsave(groupa + comps * points + groupb), &
         y0(groupa + comps * points + groupb), &
         y1(groupa + comps * points + groupb)

!///  SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.

      save

!///////////////////////////////////////////////////////////////////////
!
!     PROLOGUE.
!
!///////////////////////////////////////////////////////////////////////

!///  INITIALIZE.

!     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

!     TURN OF REVERSE COMMUNICATION FLAGS.
      time = .false.

!     TURN OFF ALL COMPLETION STATUS FLAGS.
      error = .false.
      report = qnull
      succes = .false.

!///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      if (signal /= ' ') then
         go to (1010, 1020, 1060, 1080, 1090, 2020) route
         error = .true.
         go to 9001
      end if

!///  CHECK THE ARGUMENTS.

      error = .not. (((0 < comps) .eqv. (0 < points)) .and. &
         0 <= comps .and. 0 <= points .and. 0 <= groupa .and. &
         0 <= groupb .and. 0 < groupa + comps * points + groupb)
      if (error) go to 9002

      error = .not. (0 < desire)
      if (error) go to 9003

      error = .not. (1.0 <= tdec .and. 1.0 <= tinc)
      if (error) go to 9004

      error = .not. (0.0 < tmin .and. tmin <= tmax)
      if (error) go to 9005

      error = .not. (tmin <= strid0 .and. strid0 <= tmax)
      if (error) go to 9006

      error = .not. (0 <= step)
      if (error) go to 9007

      error = 1.0 < tinc .and. .not. 0 < steps2
      if (error) go to 9008

!///  WRITE ALL MESSAGES.

!                     123456789_123456789_123456789_123456789_1234
!                     123456   123456   123456   123456   12345
      header(1, 1) = '  TIME   LOG10                      NEWTON S'
      header(1, 2) = ' POINT   ------------------------   --------'
      header(1, 3) = 'NUMBER   NORM F   CHANGE   STRIDE   STEPS   '

!                     123456789_123456789_1
!                     123   123456   123456
      header(2, 1) = 'EARCH                '
      header(2, 2) = '---------------------'
      header(2, 3) = 'J''S   COND J   REMARK'

      if (mess .and. 0 < text) then
         route = 0
         ynorm = 1.0E-4
         call twlogr (yword, ynorm)

         write (text, 10001) id, header, step, yword
         write (text, 20001) id, step, yword, log10 (strid0)
         write (text, 10002) id, header, step, yword
         write (text, 20002) id, step, yword, log10 (strid0)
         write (text, 20003) id, step, yword, log10 (strid0)
         write (text, 10006) id
         write (text, 10008) id
         write (text, 20007) id, step, yword
         write (text, 20004) id, step, yword, log10 (strid0)
         write (text, 10007) id
         write (text, 20006) id, step, yword
         write (text, 20008) id
         write (text, 20005) id, step, yword, log10 (strid0)

         go to 9001
      end if

!///////////////////////////////////////////////////////////////////////
!
!     TIME EVOLUTION.
!
!///////////////////////////////////////////////////////////////////////

!///  0 < M?

      if (0 < step) go to 1010
         stride = strid0
         age = 0

!        RETAIN THE LATEST SOLUTION FOR USE BY THE FUNCTION
         call twcopy (groupa + comps * points + groupb, v0, buffer)
         signal = 'RETAIN'
!        GO TO 1010 WHEN ROUTE = 1
         route = 1
         go to 99999
1010  continue
      signal = ' '

!///  FIRST := STEP, LAST := STEP + DESIRE.

      exist = .false.
      first = step
      last = step + desire

!///  PRINT.

      if (.not. (0 < levelm .and. 0 < text)) go to 1030
         call twcopy (groupa + comps * points + groupb, v0, buffer)
         signal = 'RESIDUAL'
         time = .false.
!        GO TO 1020 WHEN ROUTE = 2
         route = 2
         go to 99999
1020     continue
         signal = ' '
         call twnorm (groupa + comps * points + groupb, ynorm, buffer)
         call twlogr (yword, ynorm)

         if (1 == levelm) then
            if (step == 0) then
               write (text, 10001) id, header, step, yword
            else
               write (text, 10002) id, header, step, yword
            end if
         else if (1 < levelm .and. 0 == step) then
            if (step == 0) then
               write (text, 20001) id, step, yword, log10 (stride)
            else
               write (text, 20002) id, step, yword, log10 (stride)
            end if
         end if
1030  continue

!///  LOW := TMIN, HIGH := TMAX.

1040  continue

      low = tmin
      high = tmax

!///  IF AGE = STEPS2 AND STRIDE < HIGH AND 1 < TINC, THEN INCREASE
!///  STRIDE.

      if (age == steps2 .and. stride < high .and. 1.0 < tinc) &
         then
         age = 0
         exist = .false.
         low = stride * tdec
         stride = min (high, stride * tinc)
         if (1 < levelm .and. 0 < text) &
            write (text, 20003) id, step, yword, log10 (stride)
      else
         if (1 < levelm .and. 0 < text .and. 0 < step) &
            write (text, 20002) id, step, yword, log10 (stride)
      end if

!///  NEWTON SEARCH.

1050  continue

!     STORE THE LATEST SOLUTION SHOULD THE SEARCH FAIL
      call twcopy (groupa + comps * points + groupb, v0, vsave)

      count = 0
      csave = 0.0
      jword = ' '
      jacob = .false.

1060  continue

      if (jacob) then
         exist = .true.
         count = count + 1
         csave = max (condit, csave)
         if (csave == 0.0) then
            write (jword, '(I3, 3X, A6)') count, '    NA'
         else
            write (jword, '(I3, 3X, F6.2)') count, log10 (csave)
         end if
      end if

!     SUBROUTINE SEARCH &
!       (ERROR, TEXT, &
!        ABOVE, AGE, BELOW, BUFFER, COMPS, CONDIT, EXIST, GROUPA, &
!        GROUPB, LEVELD, LEVELM, NAME, NAMES, POINTS, REPORT, S0, S1, &
!        SIGNAL, STEPS, SUCCES, V0, V1, XXABS, XXAGE, XXREL, Y0, Y0NORM, &
!        Y1)

      call search &
        (error, text, &
         above, agej, below, buffer, comps, condit, exist, groupa, &
         groupb, leveld - 1, levelm - 1, name, names, points, xrepor, &
         s0, s1, signal, number, xsucce, v0, v1, tdabs, tdage, tdrel, &
         y0, dummy, y1)
      if (error) go to 9009

      if (signal /= ' ') then
         jacob = signal == 'PREPARE'
         time = .true.
!        GO TO 1060 WHEN ROUTE = 3
         route = 3
         go to 99999
      end if

!///  UNSUCCESSFUL?

      if (.not. xsucce) then
         if (1 == levelm .and. 0 < text) then
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
            write (text, 10003) step + 1, log10 (stride), number, jword, &
               remark (1 : length)
         end if

!///  IF ALSO LOW < STRIDE AND 1 < TDEC, THEN DECREASE STRIDE.

         if (low < stride .and. 1.0 < tdec) then
            age = 0
            call twcopy (groupa + comps * points + groupb, vsave, v0)
            exist = .false.
            high = stride / tinc
            stride = max (low, stride / tdec)
            if (1 < levelm .and. 0 < text) write (text, 20004) &
               id, step, yword, log10 (stride)
            go to 1050
         end if

!///  OTHERWISE END, FAILURE.

         go to 2010
      end if

!///  IF NO CHANGE AND STRIDE < HIGH AND 1.0 < TINC, THEN
!///  INCREASE STRIDE.  OTHERWISE END, FAILURE.

      do 1070 j = 1, groupa + comps * points + groupb
         buffer(j) = v0(j) - vsave(j)
1070  continue
      call twnorm (groupa + comps * points + groupb, change, buffer)
      call twlogr (cword, change)

      if (change == 0.0) then
         if (1 == levelm .and. 0 < text) then
            write (text, 10004) &
               step + 1, '  ZERO', log10 (stride), number, jword
         end if

         if (1.0 < tinc .and. stride < high) then
            age = 0
            exist = .false.
            low = stride * tdec
            stride = min (high, stride * tinc)
            if (1 < levelm .and. 0 < text) &
               write (text, 20005) id, step, yword, log10 (stride)
            go to 1050
         end if
         go to 2010
      end if

!///  AGE := AGE + 1, M := M + 1.

      age = age + 1
      step = step + 1

!     RETAIN THE LATEST SOLUTION FOR USE BY THE FUNCTION.
!     GO TO 1080 WHEN ROUTE = 4
      route = 4
      call twcopy (groupa + comps * points + groupb, v0, buffer)
      signal = 'RETAIN'
      go to 99999
1080  continue
      signal = ' '

!///  PRINT.

      if (.not. (0 < levelm .and. 0 < text)) go to 1100
         call twcopy (groupa + comps * points + groupb, v0, buffer)
         signal = 'RESIDUAL'
         time = .false.
!        GO TO 1090 WHEN ROUTE = 5
         route = 5
         go to 99999
1090     continue
         signal = ' '
         call twnorm (groupa + comps * points + groupb, ynorm, buffer)
         call twlogr (yword, ynorm)

         if (1 == levelm) write (text, 10005) &
            step, yword, cword, log10 (stride), number, jword
1100  continue

!///  M < LAST?

      if (step < last) go to 1040

!///////////////////////////////////////////////////////////////////////
!
!     EPILOGUE.
!
!///////////////////////////////////////////////////////////////////////

2010  continue

!///  PRINT.

      if (0 < levelm .and. 0 < text) then
         if (1 == levelm) then
            if (step == first) then
               write (text, 10006) id
            else if (step == last) then
               write (text, 10007) id
            else
               write (text, 10008) id
            end if
         else if (1 < levelm) then
            if (step == first) then
               write (text, 10006) id
            else if (step == last) then
               write (text, 20006) id, step, yword
            else
               write (text, 20007) id, step, yword
            end if
         end if

         if (first < last .and. 1 == leveld) then
            write (text, 20008) id
            call twcopy (groupa + comps * points + groupb, v0, buffer)
            signal = 'SHOW'
!           GO TO 2020 WHEN ROUTE = 6
            route = 6
            go to 99999
         end if
      end if

2020  continue
      signal = ' '

!///  SET THE COMPLETION STATUS FLAGS.

      succes = first < step
      if (step < last) report = xrepor

!///////////////////////////////////////////////////////////////////////
!
!     INFORMATIVE MESSAGES.
!
!///////////////////////////////////////////////////////////////////////

10001 format &
        (/1X, a9, 'BEGIN TIME EVOLUTION.' &
        /3(/10X, a44, a21) &
        /10X, i6, 3X, a6)

10002 format &
        (/1X, a9, 'CONTINUE TIME EVOLUTION.' &
        /3(/10X, a44, a21) &
        /10X, i6, 3X, a6)

10003 format &
        (10X, i6, 21X, f6.2, 3X, i5, 3X, a12, 3X, a)

10004 format &
        (10X, i6, 12X, a6, 3X, f6.2, 3X, i5, 3X, a12)

10005 format &
        (10X, i6, 2(3X, a6), 3X, f6.2, 3X, i5, 3X, a12)

10006 format &
        (/1X, a9, 'FAILURE.  NO TIME EVOLUTION.')

10007 format &
        (/1X, a9, 'SUCCESS.  TIME EVOLUTION COMPLETED.')

10008 format &
        (/1X, a9, 'PARTIAL SUCCESS.  TIME EVOLUTION INCOMPLETE.')

20001 format &
        (/1X, a9, 'BEGIN TIME EVOLUTION.' &
       //10X, i10, '  LATEST TIME POINT' &
        /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE' &
        /10X, f10.2, '  LOG10 STRIDE TO NEXT TIME POINT' &
       //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')

20002 format &
        (/1X, a9, 'CONTINUE TIME EVOLUTION.' &
       //10X, i10, '  LATEST TIME POINT' &
        /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE' &
        /10X, f10.2, '  LOG10 STRIDE TO NEXT TIME POINT' &
       //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')

20003 format &
        (/1X, a9, 'CONTINUE TIME EVOLUTION WITH INCREASED STRIDE.' &
       //10X, i10, '  LATEST TIME POINT' &
        /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE' &
        /10X, f10.2, '  LOG10 INCREASED STRIDE TO NEXT TIME POINT' &
       //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')

20004 format &
        (/1X, a9, 'RETRY THE STEP WITH A DECREASED TIME STRIDE.' &
       //10X, i10, '  LATEST TIME POINT' &
        /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE' &
        /10X, f10.2, '  LOG10 DECREASED STRIDE TO NEXT TIME POINT' &
       //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE, AGAIN.')

20005 format &
        (/1X, a9, 'THE SOLUTION DID NOT CHANGE.  RETRYING THE STEP' &
        /10X, 'WITH AN INCREASED TIME STRIDE.' &
       //10X, i10, '  LATEST TIME POINT' &
        /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE' &
        /10X, f10.2, '  LOG10 INCREASED STRIDE TO NEXT TIME POINT' &
       //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE, AGAIN.')

20006 format &
        (/1X, a9, 'SUCCESS.  TIME EVOLUTION COMPLETED.' &
       //10X, i10, '  LAST TIME POINT' &
        /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE')

20007 format &
        (/1X, a9, 'PARTIAL SUCCESS.  TIME EVOLUTION INCOMPLETE.' &
       //10X, i10, '  LAST TIME POINT' &
        /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE')

20008 format &
       (/1X, a9, 'THE LATEST SOLUTION:')

!///////////////////////////////////////////////////////////////////////
!
!     ERROR MESSAGES.
!
!///////////////////////////////////////////////////////////////////////

      go to 99999

9001  if (0 < text) write (text, 99001) id, route
      if (.not. mess) go to 99999

9002  if (0 < text) write (text, 99002) id, &
         comps, points, groupa, groupb, groupa + comps * points + groupb
      if (.not. mess) go to 99999

9003  if (0 < text) write (text, 99003) id, desire
      if (.not. mess) go to 99999

9004  if (0 < text) write (text, 99004) id, tdec, tinc
      if (.not. mess) go to 99999

9005  if (0 < text) write (text, 99005) id, tmin, tmax
      if (.not. mess) go to 99999

9006  if (0 < text) write (text, 99006) id, tmin, strid0, tmax
      if (.not. mess) go to 99999

9007  if (0 < text) write (text, 99007) id, step
      if (.not. mess) go to 99999

9008  if (0 < text) write (text, 99008) id, steps2
      if (.not. mess) go to 99999

9009  if (0 < text) write (text, 99009) id
      if (.not. mess) go to 99999

99001 format &
        (/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
       //10X, i10, '  ROUTE')

99002 format &
        (/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE' &
        /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES' &
        /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS' &
        /10X, 'MUST BE POSITIVE.' &
       //10X, i10, '  COMPS, COMPONENTS' &
        /10X, i10, '  POINTS' &
        /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
        /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
        /10X, i10, '  TOTAL UNKNOWNS')

99003 format &
        (/1X, a9, 'ERROR.  THE NUMBER OF TIME STEPS MUST BE POSITIVE.' &
       //10X, i10, '  STEPS0 OR STEPS1, DESIRED NUMBER OF STEPS')

99004 format &
        (/1X, a9, 'ERROR.  THE FACTORS FOR CHANGING THE TIME STRIDE' &
        /10X, 'MUST BE NO SMALLER THAN 1.' &
       //10X, 1p, e10.2, '  TDEC, DECREASE FACTOR', &
        /10X, 1p, e10.2, '  TINC, INCREASE FACTOR')

99005 format &
        (/1X, a9, 'ERROR.  THE BOUNDS ON THE TIME STRIDE ARE OUT OF' &
        /10X, 'ORDER.' &
       //10X, 1p, e10.2, '  TMIN, SHORTEST STRIDE' &
        /10X, 1p, e10.2, '  TMAX, LONGEST STRIDE')

99006 format &
        (/1X, a9, 'ERROR.  THE INITIAL TIME STRIDE MUST LIE BETWEEN' &
        /10X, 'THE LOWER AND UPPER BOUNDS.' &
       //10X, 1p, e10.2, '  TMIN, SHORTEST STRIDE' &
        /10X, 1p, e10.2, '  STRID0, INITIAL STRIDE' &
        /10X, 1p, e10.2, '  TMAX, LONGEST STRIDE')

99007 format &
        (/1X, a9, 'ERROR.  THE COUNT OF TIME STEPS MUST BE ZERO OR' &
        /10X, 'POSITIVE.' &
       //10X, i10, '  STEP')

99008 format &
        (/1X, a9, 'ERROR.  THE TIME STEPS BEFORE STRIDE INCREASES' &
        /10X, 'MUST BE POSITIVE.' &
       //10X, i10, '  STEPS2, TIME STEPS BEFORE STRIDE INCREASES')

99009 format &
        (/1X, a9, 'ERROR.  SEARCH FAILS.')

!///  EXIT.

      stop
99999 continue
      return
      end
      subroutine search &
        (error, text, &
         above, age, below, buffer, comps, condit, exist, groupa, &
         groupb, leveld, levelm, name, names, points, report, s0, s1, &
         signal, steps, succes, v0, v1, xxabs, xxage, xxrel, y0, y0norm, &
         y1)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     SEARCH
!
!     PERFORM THE DAMPED, MODIFIED NEWTON'S SEARCH.
!
!///////////////////////////////////////////////////////////////////////

      !implicit complex (a - z)

      character &
         column*16, ctemp1*80, ctemp2*80, header*80, id*9, name*(*), &
         signal*(*), string*80
!**** PRECISION > DOUBLE
      double precision    above, abs0, abs1, below, buffer, condit, deltab, deltad, rel0, &
         rel1, s0, s0norm, s1, s1norm, sj, temp, v0, v1, value, vj, &
         xxabs, xxrel, y0, y0norm, y1, y1norm, zero
      integer &
         age, comps, count, entry, expone, groupa, groupb, i, j, k, &
         len1, len2, length, leveld, levelm, lines, names, number, &
         points, qbnds, qdvrg, qnull, report, route, steps, text, xxage
      intrinsic &
         abs, int, log10, max, min, mod
      logical &
         error, exist, force, mess, succes

      parameter (id = 'SEARCH:  ')
      parameter (lines = 20)
      parameter (zero = 0.0)

!     REPORT CODES
      parameter (qnull = 0, qbnds = 1, qdvrg = 2)

      dimension &
         above(groupa + comps * points + groupb), &
         below(groupa + comps * points + groupb), &
         buffer(groupa + comps * points + groupb), column(7), &
         header(3, 2), name(names), &
         s0(groupa + comps * points + groupb), &
         s1(groupa + comps * points + groupb), &
         v0(groupa + comps * points + groupb), &
         v1(groupa + comps * points + groupb), &
         y0(groupa + comps * points + groupb), &
         y1(groupa + comps * points + groupb)

!///  SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.

      save

!///////////////////////////////////////////////////////////////////////
!
!     PROLOGUE.
!
!///////////////////////////////////////////////////////////////////////

!///  EVERY-TIME INITIALIZATION.

!     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

!     TURN OFF ALL COMPLETION STATUS FLAGS.
      error = .false.
      report = qnull
      succes = .false.

!///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      if (signal /= ' ') then
         go to (2020, 2040, 2050, 2140, 2150, 2180) route
         error = .true.
         go to 9001
      end if

!///  ONE-TIME INITIALIZATION.

      number = 0

!///  CHECK THE ARGUMENTS.

      error = .not. (((0 < comps) .eqv. (0 < points)) .and. &
         0 <= comps .and. 0 <= points .and. 0 <= groupa .and. &
         0 <= groupb .and. 0 < groupa + comps * points + groupb)
      if (error) go to 9002

      error = .not. (names == 1 .or. &
         names == groupa + comps + groupb)
      if (error) go to 9003

      count = 0
      do j = 1, groupa + comps * points + groupb
         if (.not. (below(j) < above(j))) count = count + 1
      end do
      error = count /= 0
      if (error) go to 9004

      count = 0
      do j = 1, groupa + comps * points + groupb
         if (.not. (below(j) <= v0(j) .and. v0(j) <= above(j))) count = count + 1
      end do
      error = count /= 0
      if (error) go to 9005

      error = .not. (zero<=xxabs .and. zero<=xxrel)
      if (error) go to 9006

      error = .not. (0 < xxage)
      if (error) go to 9007

!///  WRITE ALL MESSAGES.

      if (mess .and. 0 < text) then
         route = 0

         write (text, 10003) id
         write (text, 10002) id
         count = 0
         do 1030 j = 1, groupa + comps * points + groupb
            count = count + 1
            if (count <= lines) then
               if (j <= groupa) then
                  i = j
               else if (j <= groupa + comps * points) then
                  i = groupa + mod (j - groupa - 1, comps) + 1
               else
                  i = j - groupa - comps * points
               end if

               if (names == comps + groupa + groupb) then
                  ctemp1 = name(i)
               else
                  ctemp1 = ' '
               end if
               call twsqez (len1, ctemp1)

               if (j <= groupa) then
                  write (ctemp2, 80001) 'A', i
               else if (j <= groupa + comps * points) then
                  write (ctemp2, 80002) 'C', i, &
                     'P', int ((j - groupa - 1) / comps) + 1
               else
                  write (ctemp2, 80001) 'B', i
               end if
               call twsqez (len2, ctemp2)

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
                  len1 = 30 - len2 - 4
                  string = ctemp1 (1 : len1) // '... ' // ctemp2
                  length = 30
               end if

               write (text, 80003) 'LOWER', v0(j), string (1 : length)
            end if
1030     continue
         if (lines < count) write (text, 80004)
         write (text, 10001) id
         write (text, 10006) id
         write (text, 10005) id

         go to 9001
      end if

!///  PRINT THE HEADER.

!                     123456789_123456789_123456789_123456789_1234
!                     123456   123456   123456   123456   123456
      header(1, 1) = '         LOG10                              '
      header(2, 1) = '  SLTN   -----------------------------------'
      header(3, 1) = 'NUMBER   NORM F   COND J   NORM S      ABS A'

!                     123456789_123456789_123
!                     123456   123456  123456
      header(1, 2) = '                       '
      header(2, 2) = '-----------------------'
      header(3, 2) = 'ND REL    DELTA B AND D'

      if (levelm >= 1 .or. mess) then
         if (0 < text) write (text, 10001) &
            id, ((header(j, k), k = 1, 2), j = 1, 3)
      end if

!///////////////////////////////////////////////////////////////////////
!
!     SIR ISSAC NEWTON'S ALGORITHM.
!
!///////////////////////////////////////////////////////////////////////

!///  J EXIST?

      if (.not. exist) go to 2010

!///  AGE < XXAGE?

      if (age < xxage) go to 2030

!///  EVALUATE J AT V0.  RE-EVALUATE Y0 := F(V0) IN CASE F CHANGES WHEN
!///  J DOES.  SOLVE J S0 = Y0.  EVAUATE ABS0 AND REL0.

2010  continue

      call twcopy (groupa + comps * points + groupb, v0, buffer)
      signal = 'PREPARE'
!     GO TO 2020 WHEN ROUTE = 1
      route = 1
      go to 99999
2020  continue
      signal = ' '
      age = 0

!     JACOBIAN EVALUATION SHOULD RETURN A NEW RESIDUAL TOO.

      if (0 < levelm .and. 0 < text) then
         if (0.0 < condit) then
            write (column(2), '(F6.2)') log10 (condit)
         else
            column(2) = '    NA'
         end if
      end if

!///  EVALUATE Y0 := F(V0).  SOLVE J S0 = Y0.  EVAUATE ABS0 AND REL0.

2030  continue

      call twcopy (groupa + comps * points + groupb, v0, buffer)
      signal = 'RESIDUAL'
!     GO TO 2040 WHEN ROUTE = 2
      route = 2
      go to 99999
2040  continue
      signal = ' '
      call twcopy (groupa + comps * points + groupb, buffer, y0)
      call twnorm (groupa + comps * points + groupb, y0norm, y0)

      call twcopy (groupa + comps * points + groupb, y0, buffer)
      signal = 'SOLVE'
!     GO TO 2050 WHEN ROUTE = 3
      route = 3
      go to 99999
2050  continue
      signal = ' '
      call twcopy (groupa + comps * points + groupb, buffer, s0)
      call twnorm (groupa + comps * points + groupb, s0norm, s0)

      abs0 = 0.0
      rel0 = 0.0
      do 2060 j = 1, groupa + comps * points + groupb
         sj = abs (v0(j) - (v0(j) - s0(j)))
         vj = abs (v0(j))
         if (xxrel * vj < sj) abs0 = max (abs0, sj)
         if (xxabs < sj .and. 0.0 < vj) &
            rel0 = max (rel0, sj / vj)
2060  continue

!///  CHECK FOR SUCCESS.

      if (abs0 <= xxabs .and. rel0 <= xxrel) go to 2170

!///  CHOOSE DELTAB.

2070  continue

!     DELTAB IS THE LARGEST DAMPING COEFFICIENT BETWEEN 0 AND 1 THAT
!     KEEPS V1 WITHIN BOUNDS.  IF V1 BELONGS ON THE BOUNDARY, THEN
!     PROVISIONS ARE MADE TO FORCE IT THERE DESPITE ROUNDING ERROR.

      deltab = 1.0
      force = .false.
      do 2080 j = 1, groupa + comps * points + groupb
         if (s0(j) > max (zero, v0(j) - below(j))) then
            temp = (v0(j) - below(j)) / s0(j)
            if (temp < deltab) then
               deltab = temp
               entry = j
               force = .true.
               value = below(j)
            end if
         else if (s0(j) < min (zero, v0(j) - above(j))) then
            temp = (v0(j) - above(j)) / s0(j)
            if (temp < deltab) then
               deltab = temp
               entry = j
               force = .true.
               value = above(j)
            end if
         end if
2080  continue

      error = deltab < 0.0
      if (error) go to 9008

!///  0 < DELTAB?

      if (.not. (0.0 < deltab)) then
         if (0 < age) go to 2010

         if (0 < levelm .and. 0 < text) then
            call twlogr (column(1), y0norm)
            call twlogr (column(3), s0norm)
            call twlogr (column(4), abs0)
            call twlogr (column(5), rel0)
            column(6) = ' '
            if (deltab /= 1.0) call twlogr (column(6), deltab)
            column(7) = ' '
            if (deltad /= 1.0) call twlogr (column(7), deltad)
            write (text, 10004) number, column
            write (text, 10002) id

            count = 0
            do 2090 j = 1, groupa + comps * points + groupb
               if ((below(j) == v0(j) .and. 0.0 < s0(j)) .or. &
                  (v0(j) == above(j) .and. s0(j) < 0.0)) then
                  count = count + 1
                  if (count <= lines) then
                     if (j <= groupa) then
                        i = j
                     else if (j <= groupa + comps * points) then
                        i = groupa + mod (j - groupa - 1, comps) + 1
                     else
                        i = j - groupa - comps * points
                     end if

                     if (names == comps + groupa + groupb) then
                        ctemp1 = name(i)
                     else
                        ctemp1 = ' '
                     end if
                     call twsqez (len1, ctemp1)

                     if (j <= groupa) then
                        write (ctemp2, 80001) 'A', i
                     else if (j <= groupa + comps * points) then
                        write (ctemp2, 80002) 'C', i, &
                           'P', int ((j - groupa - 1) / comps) + 1
                     else
                        write (ctemp2, 80001) 'B', i
                     end if
                     call twsqez (len2, ctemp2)

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
                        len1 = 30 - len2 - 4
                        string = ctemp1 (1 : len1) // '... ' // ctemp2
                        length = 30
                     end if

                     if (below(j) == v0(j)) then
                        write (text, 80003) &
                           'LOWER', v0(j), string (1 : length)
                     else
                        write (text, 80003) &
                           'UPPER', v0(j), string (1 : length)
                     end if
                  end if
               end if
2090        continue
            if (lines < count) write (text, 80005)
         end if

         report = qbnds
         succes = .false.
         go to 99999
      end if

!///  DELTAD := 1.

      deltad = 1.0
      expone = 0

!///  V1 := V0 - DELTAB DELTAD S0.  EVALUATE Y1 := F(V1).  SOLVE
!///  J S1 = Y1.  EVALUATE ABS1 AND REL1.

2100  continue

      temp = deltab * deltad
      do 2110 j = 1, groupa + comps * points + groupb
         v1(j) = v0(j) - temp * s0(j)
2110  continue

!     KEEP V1 IN BOUNDS DESPITE ROUNDING ERROR.

      do 2120 j = 1, groupa + comps * points + groupb
         v1(j) = min (v1(j), above(j))
2120  continue
      do 2130 j = 1, groupa + comps * points + groupb
         v1(j) = max (v1(j), below(j))
2130  continue
      if (expone == 0 .and. force) v1(entry) = value

      call twcopy (groupa + comps * points + groupb, v1, buffer)
      signal = 'RESIDUAL'
!     GO TO 2140 WHEN ROUTE = 4
      route = 4
      go to 99999
2140  continue
      signal = ' '
      call twcopy (groupa + comps * points + groupb, buffer, y1)
      call twnorm (groupa + comps * points + groupb, y1norm, y1)

      call twcopy (groupa + comps * points + groupb, y1, buffer)
      signal = 'SOLVE'
!     GO TO 2150 WHEN ROUTE = 5
      route = 5
      go to 99999
2150  continue
      signal = ' '
      call twcopy (groupa + comps * points + groupb, buffer, s1)
      call twnorm (groupa + comps * points + groupb, s1norm, s1)

      abs1 = 0.0
      rel1 = 0.0
      do 2160 j = 1, groupa + comps * points + groupb
         sj = abs (v1(j) - (v1(j) - s1(j)))
         vj = abs (v1(j))
         if (xxrel * vj < sj) abs1 = max (abs1, sj)
         if (xxabs < sj .and. 0.0 < vj) &
            rel1 = max (rel1, sj / vj)
2160  continue

!///  NORM S1 < OR = NORM S0?

      if (s1norm <= s0norm) then
      else
         deltad = 0.5 * deltad
         expone = expone + 1
         if (expone <= 5) go to 2100
            if (0 < age) go to 2010
               if (0 < levelm .and. 0 < text) then
                  call twlogr (column(1), y0norm)
                  call twlogr (column(3), s0norm)
                  call twlogr (column(4), abs0)
                  call twlogr (column(5), rel0)
                  column(6) = ' '
                  if (deltab /= 1.0) call twlogr (column(6), deltab)
                  column(7) = ' '
                  if (deltad /= 1.0) call twlogr (column(7), deltad)
                  write (text, 10004) number, column
                  write (text, 10003) id
               end if
               report = qdvrg
               succes = .false.
               go to 99999
      end if

!///  PRINT.

      if (0 < levelm .and. 0 < text) then
         call twlogr (column(1), y0norm)
         call twlogr (column(3), s0norm)
         call twlogr (column(4), abs0)
         call twlogr (column(5), rel0)
         column(6) = ' '
         if (deltab /= 1.0) call twlogr (column(6), deltab)
         column(7) = ' '
         if (deltad /= 1.0) call twlogr (column(7), deltad)
         write (text, 10004) number, column
         column(2) = ' '
      end if

!///  S0 := S1, U := V1, Y0 := Y1, AGE := AGE + 1.

      age = age + 1
      number = number + 1
      call twcopy (groupa + comps * points + groupb, s1, s0)
      call twcopy (groupa + comps * points + groupb, v1, v0)
      call twcopy (groupa + comps * points + groupb, y1, y0)
      s0norm = s1norm
      y0norm = y1norm
      abs0 = abs1
      rel0 = rel1

!///  S0 SMALL VS V0?

      if (.not. (abs0 <= xxabs .and. rel0 <= xxrel)) then
         if (age < xxage) go to 2070
         go to 2010
      end if

!///  SUCCESS.

2170  continue

!///  PRINT.

      if (0 < levelm .and. 0 < text) then
         call twlogr (column(1), y0norm)
         call twlogr (column(3), s0norm)
         call twlogr (column(4), abs0)
         call twlogr (column(5), rel0)
         column(6) = ' '
         column(7) = ' '
         if (0 < leveld) then
            write (text, 10004) number, column
            write (text, 10005) id
            signal = 'SHOW'
            call twcopy (groupa + comps * points + groupb, v0, buffer)
!           GO TO 2180 WHEN ROUTE = 6
            route = 6
            go to 99999
         else
            write (text, 10004) number, column
            write (text, 10006) id
         end if
      end if

2180  continue
      signal = ' '

      succes = .true.

!///////////////////////////////////////////////////////////////////////
!
!     INFORMATIVE MESSAGES.
!
!///////////////////////////////////////////////////////////////////////

10001 format &
        (/1X, a9, 'SOLVE NONLINEAR, NONDIFFERENTIAL EQUATIONS.' &
        /4(/10X, a44, a23)/)

10002 format &
       (/1X, a9, 'FAILURE.  THE SEARCH FOR THE FOLLOWING UNKNOWNS GOES' &
        /10X, 'OUT OF BOUNDS.' &
       //10X, 'BOUND       VALUE   UNKNOWN' &
        /)

10003 format &
        (/1X, a9, 'FAILURE.  THE SEARCH DIVERGES.')

10004 format &
        (10X, i6, 3(3X, a6), 2(3X, a6, 2X, a6))

10005 format &
        (/1X, a9, 'SUCCESS.  THE SOLUTION:')

10006 format &
        (/1X, a9, 'SUCCESS.')

80001 format &
         ('(', a, ' ', i10, ')')

80002 format &
        ('(', a, ' ', i10, ' ', a, ' ', i10, ')')

80003 format &
        (10X, a5, 2X, 1p, e10.2, 3X, a)

80004 format &
        (30X, '... MORE')

80005 format &
        (10X, '  ... MORE')

80006 format &
        (10X, 1p, e10.2, 2X, e10.2, 3X, a)

80007 format &
        (10X, 1p, e10.2, 2X, e10.2, 2X, e10.2, 3X, a)

!///////////////////////////////////////////////////////////////////////
!
!     ERROR MESSAGES.
!
!///////////////////////////////////////////////////////////////////////

      go to 99999

9001  if (0 < text) write (text, 99001) id, route
      if (.not. mess) go to 99999

9002  if (0 < text) write (text, 99002) id, &
         comps, points, groupa, groupb, groupa + comps * points + groupb
      if (.not. mess) go to 99999

9003  if (0 < text) write (text, 99003) id, &
         names, comps, groupa, groupb, groupa + comps + groupb
      if (.not. mess) go to 99999

9004  if (0 < text) then
         write (text, 99004) id, &
            groupa, groupb, comps, groupa + comps + groupb, count
         count = 0
         do 8010 j = 1, groupa + comps + groupb
            if (.not. (below(j) < above(j)) .or. mess) then
               count = count + 1
               if (count <= lines) then
                  if (names == comps + groupa + groupb) then
                     ctemp1 = name(j)
                  else
                     ctemp1 = ' '
                  end if
                  call twsqez (len1, ctemp1)

                  if (j <= groupa) then
                     write (ctemp2, 80001) 'A', j
                  else if (j <= groupa + comps) then
                     write (ctemp2, 80001) 'C', j - groupa
                  else
                     write (ctemp2, 80001) 'B', j - groupa - comps
                  end if
                  call twsqez (len2, ctemp2)

                  if (ctemp1 == ' ') then
                     string = ctemp2
                     length = len2
                  else if (len1 + 2 + len2 <= 40) then
                     string = ctemp1 (1 : len1) // '  ' // ctemp2
                     length = len1 + 2 + len2
                  else if (len1 + 1 + len2 <= 40) then
                     string = ctemp1 (1 : len1) // ' ' // ctemp2
                     length = len1 + 1 + len2
                  else
                     len1 = 40 - len2 - 4
                     string = ctemp1 (1 : len1) // '... ' // ctemp2
                     length = 40
                  end if

                  write (text, 80006) &
                     below(j), above(j), string (1 : length)
               end if
            end if
8010     continue
         if (lines < count) write (text, 80005)
      end if
      if (.not. mess) go to 99999

9005  if (0 < text) then
         write (text, 99005) id, groupa, groupb, comps, points, &
            groupa + comps * points + groupb, count
         count = 0
         do 8020 j = 1, groupa + comps * points + groupb
            if (.not. (below(j) <= v0(j) .and. v0(j) <= above(j)) &
               .or. mess) then
               count = count + 1
               if (count <= lines) then
                  if (j <= groupa) then
                     i = j
                  else if (j <= groupa + comps * points) then
                     i = groupa + mod (j - groupa - 1, comps) + 1
                  else
                     i = j - groupa - comps * points
                  end if

                  if (names == comps + groupa + groupb) then
                     ctemp1 = name(i)
                  else
                     ctemp1 = ' '
                  end if
                  call twsqez (len1, ctemp1)

                  if (j <= groupa) then
                     write (ctemp2, 80001) 'A', i
                  else if (j <= groupa + comps * points) then
                     write (ctemp2, 80002) 'C', i, &
                        'P', int ((j - groupa - 1) / comps) + 1
                  else
                     write (ctemp2, 80001) 'B', i
                  end if
                  call twsqez (len2, ctemp2)

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
                     len1 = 30 - len2 - 4
                     string = ctemp1 (1 : len1) // '... ' // ctemp2
                     length = 30
                  end if

                  write (text, 80007) &
                     below(j), v0(j), above(j), string (1 : length)
               end if
            end if
8020     continue
         if (lines < count) write (text, 80005)
      end if
      if (.not. mess) go to 99999

9006  if (0 < text) write (text, 99006) id, xxabs, xxrel
      if (.not. mess) go to 99999

9007  if (0 < text) write (text, 99007) id, xxage
      if (.not. mess) go to 99999

9008  if (0 < text) write (text, 99008) id, deltab
      if (.not. mess) go to 99999

99001 format &
        (/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
       //10X, i10, '  ROUTE')

99002 format &
        (/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE' &
        /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES' &
        /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS' &
        /10X, 'MUST BE POSITIVE.' &
       //10X, i10, '  COMPS, COMPONENTS' &
        /10X, i10, '  POINTS' &
        /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
        /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
        /10X, i10, '  TOTAL UNKNOWNS')

99003 format &
        (/1X, a9, 'ERROR.  THE NUMBER OF NAMES IS WRONG.' &
       //10X, i10, '  NAMES' &
       //10X, i10, '  COMPS, COMPONENTS' &
        /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
        /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
        /10X, i10, '  TOTAL NUMBER')

99004 format &
        (/1X, a9, 'ERROR.  THE LOWER AND UPPER BOUNDS ON SOME UNKNOWNS' &
        /10X, 'ARE OUT OF ORDER.' &
       //10X, i10, '  GROUP A UNKNOWNS (A)' &
        /10X, i10, '  GROUP B UNKNOWNS (B)' &
        /10X, i10, '  COMPONENTS AT POINTS (C)' &
        /10X, i10, '  TOTAL TYPES OF UNKNOWNS' &
        /10X, i10, '  NUMBER OF BOUNDS OUT OF ORDER' &
       //10X, '     LOWER       UPPER' &
        /10X, '     BOUND       BOUND   UNKNOWN' &
        /)

99005 format &
        (/1X, a9, 'ERROR.  THE GUESSES FOR SOME UNKNOWNS ARE OUT OF' &
        /10X, 'BOUNDS.' &
       //10X, i10, '  GROUP A UNKNOWNS (A)' &
        /10X, i10, '  GROUP B UNKNOWNS (B)' &
        /10X, i10, '  COMPONENTS AT POINTS (C)' &
        /10X, i10, '  POINTS (P)' &
        /10X, i10, '  TOTAL UNKNOWNS' &
        /10X, i10, '  NUMBER OUT OF BOUNDS' &
       //10X, '     LOWER                   UPPER' &
        /10X, '     BOUND       VALUE       BOUND   UNKNOWN' &
        /)

99006 format &
        (/1X, a9, 'ERROR.  THE BOUNDS FOR THE ABSOLUTE AND RELATIVE' &
        /10X, 'CONVERGENCE TESTS MUST BE ZERO OR POSITIVE.' &
       //10X, 1p, e10.2, '  SSABS OR TDABS, ABSOLUTE ERROR' &
        /10X, 1p, e10.2, '  SSREL OR TDREL, RELATIVE ERROR')

99007 format &
        (/1X, a9, 'ERROR.  THE RETIREMENT AGE OF THE JACOBIAN MATRIX' &
        /10X, 'MUST BE POSITIVE.' &
       //10X, i10, '  SSAGE OR TDAGE, MATRIX RETIREMENT AGE')

99008 format &
        (/1X, a9, 'ERROR.  THE DAMPING COEFFICIENT FOR STAYING' &
        /10X, 'IN BOUNDS IS NEGATIVE.' &
       //10X, 1p, e10.2, '  DELTA B')

!///  EXIT.

      stop
99999 continue

!     COPY THE PROTECTED LOCAL VARIABLE
      steps = number

      return
      end subroutine search


      ! TWOPNT driver.
      subroutine twopnt(setup, error, text, versio, &
                        above, active, below, buffer, comps, condit, groupa, groupb, &
                        isize, iwork, mark, name, names, pmax, points, report, rsize, &
                        rwork, signal, stride, time, u, x)

      type(twcom), intent(inout) :: setup

      character &
         column*80, ctemp1*80, ctemp2*80, header*80, name*(*), &
         report*(*), signal*(*), string*80, versio*(*), vnmbr*8
      real(RK) ::  above, below, buffer, condit, detail, maxcon, ratio, rwork, stride, temp, timer, total, u, x, ynorm
      integer &
         age, cntrls, comps, count, desire, event, grid, groupa, &
         groupb, ilast, isize, iwork, j, jacobs, k, label, len1, &
         len2, length, names, nsteps, pmax, &
         points, psave, qabove, qbelow, qrat1, qrat2, &
         qs0, qs1, qtask, qtype, qusave, &
         qv1, qvary, qvary1, qvary2, qvsave, qxsave, qy0, qy1, return, &
         rlast, route, rsize, size, step, steps, text, xrepor
      intrinsic :: max
      logical &
         active, allow, error, exist, first, flag, found, &
         mark, satisf, time

      character(*), parameter :: id = 'TWOPNT:  '
      integer,      parameter :: gmax = 100
      integer,      parameter :: lines = 20
      integer,      parameter :: vnmbrs = 12

      ! SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      logical,      parameter :: mess = .false.

      dimension &
         above(groupa + comps + groupb), active(*), below(groupa + comps &
         + groupb), buffer(groupa + comps * pmax + groupb), column(3), &
         detail(gmax, qtotal), event(gmax, qtotal), header(6), &
          iwork(isize), mark(*), &
         name(names), ratio(2), rwork(rsize), &
         size(gmax), timer(qtotal), total(qtotal), u(groupa + comps * &
         pmax + groupb), vnmbr(vnmbrs), x(*)


!///  SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.

      save

!///////////////////////////////////////////////////////////////////////
!
!     PROLOGUE.
!
!///////////////////////////////////////////////////////////////////////

!///  EVERY-TIME INITIALIZATION.

!     TURN OFF ALL REVERSE COMMUNICATION FLAGS.
      time = .false.

!///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      if (signal /= ' ') then
         go to (9912, 9922, 9932, 9942) route
         error = .true.
         go to 9001
      end if

!///////////////////////////////////////////////////////////////////////
!
!     ENTRY BLOCK.  INITIALIZE A NEW PROBLEM.
!
!///////////////////////////////////////////////////////////////////////

!///  TURN OFF ALL STATUS REPORTS.

      error = .false.
      report = ' '

!///  WRITE ALL MESSAGES.

      if (mess .and. 0 < text) then
         label = 0
         return = 0
         route = 0

         write (text, 10004) id, '???'
         write (text, 10020) id
         write (text, 10017) id
         write (text, 10014) id
         write (text, 10001) id, 'DOUBLE PRECISION', twtrim(vnmbr(vnmbrs),CONTRL_MAX_LEN)
         write (text, 10022) id
         write (text, 10021) id
         write (text, 10011) id, '???', ratio, setup%toler1, setup%toler2
         write (text, 10013) id, '???', ratio, setup%toler1, setup%toler2
         write (text, 10012) id
         write (text, 10002) id, 'FINAL SOLUTION:'
         write (text, 10002) id, 'INITIAL GUESS:'
         write (text, 10019) id
         write (text, 10018) id
         write (text, 10016) id
         write (text, 10015) id
         write (text, 10001) id, 'SINGLE PRECISION', twtrim(vnmbr(vnmbrs),CONTRL_MAX_LEN)
         write (text, 10002) id, 'SOLVE THE PROBLEM.'
         write (text, 10010) id

         go to 9001
      end if

!///  CHECK THE VERSION.

      data vnmbr &
         / '3.18', '3.19', '3.20', '3.21', '3.22', '3.23', '3.24', &
           '3.25', '3.26', '3.27', '3.28', '3.29' /

      flag = .false.
      do j = 1, vnmbrs
         flag = flag .or. versio == 'DOUBLE PRECISION VERSION ' // vnmbr(j)
      end do
      error = .not. flag
      if (error) go to 9002

      ! Additional settings initialization
      if (.not.setup%padd) setup%ipadd = pmax


      error = .not. (count == cntrls)
      if (error) go to 9004

!///  PRINT THE ENTRY BANNER AT ALL PRINT LEVELS.

      string = vnmbr(vnmbrs)
      call twlast (length, string)
      if ((0 < setup%levelm .or. mess) .and. 0 < text) &
         write (text, 10001) id, 'DOUBLE PRECISION', string (1:length)

!///  CHECK THE ARGUMENTS.

      error = .not. (setup%leveld <= setup%levelm)
      if (error) go to 9005

      error = .not. (0 <= comps .and. 0 <= points .and. &
         0 <= groupa .and. 0 <= groupb)
      if (error) go to 9006

      error = .not. ((0 < comps) .eqv. (0 < points))
      if (error) go to 9007

      error = .not. (0 < groupa + comps * points + groupb)
      if (error) go to 9008

      error = .not. (names == 1 .or. &
         names == groupa + comps + groupb)
      if (error) go to 9009

      error = .not. (points <= pmax)
      if (error) go to 9010

      count = 0
      do 1020 j = 1, groupa + comps + groupb
         if (.not. (below(j) < above(j))) count = count + 1
1020  continue
      error = count /= 0
      if (error) go to 9011

!///  PARTITION THE INTEGER WORK SPACE.

!     SUBROUTINE TWGRAB (ERROR, LAST, FIRST, NUMBER)

      ilast = 0

!     VARY(PMAX)
      call twgrab (error, ilast, qvary, pmax)
      if (error) go to 9012

!     VARY1(PMAX)
      call twgrab (error, ilast, qvary1, pmax)
      if (error) go to 9012

!     VARY2(PMAX)
      call twgrab (error, ilast, qvary2, pmax)
      if (error) go to 9012

!///  PARTITION THE REAL WORK SPACE.

!     SUBROUTINE TWGRAB (ERROR, LAST, FIRST, NUMBER)

      rlast = 0

!     ABOVE(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qabove, groupa + comps * pmax + groupb)
      if (error) go to 9012

!     BELOW(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qbelow, groupa + comps * pmax + groupb)
      if (error) go to 9012

!     RATIO1(PMAX)
      call twgrab (error, rlast, qrat1, pmax)
      if (error) go to 9012

!     RATIO2(PMAX)
      call twgrab (error, rlast, qrat2, pmax)
      if (error) go to 9012

!     S0(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qs0, groupa + comps * pmax + groupb)
      if (error) go to 9012

!     S1(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qs1, groupa + comps * pmax + groupb)
      if (error) go to 9012

!     USAVE(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qusave, groupa + comps * pmax + groupb)
      if (error) go to 9012

!     VSAVE(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qvsave, groupa + comps * pmax + groupb)
      if (error) go to 9012

!     V1(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qv1, groupa + comps * pmax + groupb)
      if (error) go to 9012

!     XSAVE(PMAX)
      call twgrab (error, rlast, qxsave, pmax)
      if (error) go to 9012

!     Y0(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qy0, groupa + comps * pmax + groupb)
      if (error) go to 9012

!     Y1(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qy1, groupa + comps * pmax + groupb)
      if (error) go to 9012

!///  CHECK THE WORK SPACES' SIZES.

      error = .not. (ilast <= isize .and. rlast <= rsize)
      if (error) go to 9013

!///  ONE-TIME INITIALIZATION.

!     ALLOW FURTHER TIME EVOLUTION
      allow = .true.

!     PRESENT TASK
      qtask = qentry

!     STATISTICS ARRAYS
      do 1040 k = 1, qtotal
         total(k) = 0.0
         do 1030 j = 1, gmax
            detail(j, k) = 0.0
            event(j, k) = 0
1030     continue
1040  continue

!     TOTAL TIME STATISTIC
      call twtime (timer(qtotal))

!     GRID POINTER AND STATISTICS FOR THE FIRST GRID
      grid = 1
      size(grid) = points
      call twtime (timer(qgrid))

!     TIME STEP NUMBER
      step = 0

!     SOLUTION FLAG
      found = .true.

!///  EXPAND THE BOUNDS.

      count = 0

      do 1050 j = 1, groupa
         rwork(qabove + count) = above(j)
         rwork(qbelow + count) = below(j)
         count = count + 1
1050  continue

      do 1070 k = 1, points
         do 1060 j = 1, comps
            rwork(qabove + count) = above(groupa + j)
            rwork(qbelow + count) = below(groupa + j)
            count = count + 1
1060     continue
1070  continue

      do 1080 j = 1, groupb
         rwork(qabove + count) = above(groupa + comps + j)
         rwork(qbelow + count) = below(groupa + comps + j)
         count = count + 1
1080  continue

!///  SAVE THE INITIAL SOLUTION.

      psave = points
      if (setup%adapt .and. 0 < points) &
         call twcopy (points, x, rwork(qxsave))
      call twcopy (groupa + comps * points + groupb, u, rwork(qusave))

!     GO TO 1090 WHEN RETURN = 1
      return = 1
      go to 9911
1090  continue

!///  PRINT LEVELS 11, 21, AND 22.

      if (0 < setup%leveld .and. 0 < text) then
         write (text, 10002) id, 'INITIAL GUESS:'
!        GO TO 1100 WHEN RETURN = 2
         return = 2
         go to 9921
      end if
1100  continue

!///  PRINT LEVEL 10 AND 11.

!                  123456789_123456789_123456789_1234
!                  12345678   123456  123456   123456
      header(1) = '            LOG10   LOG10         '
      header(2) = '    TASK   NORM F  COND J   REMARK'

      if (setup%levelm == 1 .and. 0 < text) then
         if (0 < setup%leveld) write (text, 10002) id, &
            'SOLVE THE PROBLEM.'
         write (text, 10003) (header(j), j = 1, 2)
!        GO TO 1110 WHEN LABEL = 1
         label = 1
         go to 7010
      end if
1110  continue

!///////////////////////////////////////////////////////////////////////
!
!     DECISION BLOCK.  THE PREVIOUS TASK DETERMINES THE NEXT.
!
!///////////////////////////////////////////////////////////////////////

2010  continue

!///  ENTRY WAS THE PREVIOUS TASK.

      if (qtask == qentry) then
         if (0 < setup%steps0) then
            qtask = qtimst
            desire = setup%steps0
         else if (setup%steady) then
            qtask = qsearc
         else
            error = .true.
            go to 9014
         end if

!///  SEARCH WAS THE PREVIOUS TASK.

      else if (qtask == qsearc) then
         if (found) then
            if (setup%adapt) then
               qtask = qrefin
            else
               qtask = qexit
               report = ' '
            end if
         else
            if (allow .and. 0 < setup%steps1) then
               qtask = qtimst
               desire = setup%steps1
            else
               qtask = qexit
               if (1 < grid) then
                  report = 'SOME SOLVED'
               else
                  report = 'NOT SOLVED'
               end if
            end if
         end if

!///  REFINE WAS THE PREVIOUS TASK.

      else if (qtask == qrefin) then
         if (found) then
            step = 0
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

!///  EVOLVE WAS THE PREVIOUS TASK.

      else if (qtask == qtimst) then
         if (found) then
            if (setup%steady) then
               qtask = qsearc
            else
               qtask = qexit
               report = ' '
            end if
         else
            qtask = qexit
            if (1 < grid) then
               report = 'SOME SOLVED'
            else
               report = 'NOT SOLVED'
            end if
         end if
      end if

!///  BRANCH TO THE NEXT TASK.

      if (qtask == qexit) go to 3010
      if (qtask == qsearc) go to 4010
      if (qtask == qrefin) go to 5010
      if (qtask == qtimst) go to 6010
      error = .true.
      go to 9015

!///////////////////////////////////////////////////////////////////////
!
!     EXIT BLOCK.
!
!///////////////////////////////////////////////////////////////////////

3010  continue

!///  COMPLETE STATISTICS FOR THE LAST GRID.

      call twlaps (timer(qgrid))
      if (grid <= gmax) then
         detail(grid, qgrid) = timer(qgrid)
         detail(grid, qother) &
            = detail(grid, qgrid) - (detail(grid, qfunct) &
            + detail(grid, qjacob) + detail(grid, qsolve))
      end if

!///  RESTORE THE SOLUTION.

      if (report /= ' ') then
!        BE CAREFUL NOT TO ASSIGN A VALUE TO A PARAMETER
         if (points /= psave) points = psave
         if (setup%adapt .and. 0 < points) &
            call twcopy (points, rwork(qxsave), x)
         call twcopy &
            (groupa + comps * points + groupb, rwork(qusave), u)
      end if

!///  PRINT LEVEL 11 OR 21.

!     SAVE THE STATUS REPORTS DURING REVERSE COMMUNICATION
      string = report

      if (setup%leveld == 1 .and. 0 < text) then
         write (text, 10002) id, 'FINAL SOLUTION:'
!        GO TO 3020 WHEN RETURN = 3
         return = 3
         go to 9921
      end if
3020  continue

!     RESTORE THE STATUS REPORTS AFTER REVERSE COMMUNICATION
      report = string

!///  COMPLETE THE TOTAL TIME STATISTICS.

      call twlaps (timer(qtotal))
      total(qtotal) = timer(qtotal)
      total(qother) = total(qtotal) &
         - (total(qfunct) + total(qjacob) + total(qsolve))

!///  TOP OF THE REPORT BLOCK.

      if (0 < setup%levelm .and. 0 < text) then
         if (0.0 < total(qtotal)) then

!///  REPORT TOTAL COMPUTER TIME.

      temp = total(qtotal)
      if (3600.0 <= temp) then
         write (string, '(F10.2, A)') temp / 3600.0, ' HOURS'
      else if (60.0 <= temp) then
         write (string, '(F10.2, A)') temp / 60.0, ' MINUTES'
      else
         write (string, '(F10.2, A)') temp, ' SECONDS'
      end if

      call twsqez (length, string)
      write (text, 10004) id, string (1 : length)

!///  REPORT PERCENT OF TOTAL COMPUTER TIME.

      temp = 100.0 / total(qtotal)
      if (setup%adapt) then

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

      write (text, 10005) header, &
         (size(j), (temp * detail(j, k), k = 1, 8), j = 1, grid)
      if (1 < grid) write (text, 10006) (temp * total(k), k = 2, 8)
      if (gmax < grid) write (text, 10007)

      else

!                  123456789_123456789_123456789_123456789_123456789_1
!                  123456   123456   123456   123456   123456   123456
      header(1) = 'SUBTASK                             TASK           '
      header(2) = '---------------------------------   ---------------'
      header(3) = 'EVAL F   PREP J    SOLVE    OTHER   EVOLVE   SEARCH'

      write (text, 10008) &
         (header(j), j = 1, 3), '  % OF TOTAL', &
         (temp * total(k), k = 5, 8), (temp * total(k), k = 2, 3), &
         'MEAN SECONDS', (detail(1, k) / event(1, k), k = 5, 7), &
         '    QUANTITY', (event(1, k), k = 5, 7)

      end if

!///  REPORT AVERAGE COMPUTER TIME.

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

      if (setup%adapt) write (text, 10009) header, &
         (size(j), (detail(j, k) / event(j, k), k = 5, 7), &
         (event(j, k), k = 5, 7), j = 1, grid)

      end if

!///  REPORT THE COMPLETION STATUS.

      if (0 < setup%levelm) then
         if (report == ' ') then
            write (text, 10010) id
         else if (report == 'NO SPACE') then
            write (string, '(I10)') points
            call twsqez (length, string)
            write (text, 10011) &
               id, string (1 : length), ratio, setup%toler1, setup%toler2
         else if (report == 'NOT SOLVED') then
            write (text, 10012) id
         else if (report == 'SOME SOLVED') then
            write (string, '(I10)') points
            call twsqez (length, string)
            write (text, 10013) &
               id, string (1 : length), ratio, setup%toler1, setup%toler2
         else
            error = .true.
            go to 9016
         end if
      end if

!///  BOTTOM OF THE REPORT BLOCK.

      end if

!///  BOTTOM OF THE EXIT BLOCK.

      return

!///////////////////////////////////////////////////////////////////////
!
!     SEARCH BLOCK.
!
!///////////////////////////////////////////////////////////////////////

4010  continue

!///  INITIALIZE STATISTICS ON ENTRY TO THE SEARCH BLOCK.

      call twtime (timer(qsearc))
      first = .true.
      jacobs = 0
      maxcon = 0.0

!///  PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE SEARCH BLOCK.

      if (1 < setup%levelm) then
         if (0 < text) write (text, 10014) id
      end if

!///  PREPARE TO CALL SEARCH.

!     SAVE THE SOLUTION SHOULD THE SEARCH FAIL
      call twcopy (groupa + comps * points + groupb, u, rwork(qvsave))

      exist = .false.

!///  CALL SEARCH.

      age = 0
4020  continue

!     SUBROUTINE SEARCH &
!       (ERROR, TEXT, &
!        ABOVE, AGE, BELOW, BUFFER, COMPS, CONDIT, EXIST, GROUPA, &
!        GROUPB, setup%leveld, setup%levelm, NAME, NAMES, POINTS, REPORT, S0, S1, &
!        SIGNAL, STEPS, SUCCES, V0, V1, XXABS, XXAGE, XXREL, Y0, Y0NORM, &
!        Y1)

      call search &
        (error, text, &
         rwork(qabove), age, rwork(qbelow), buffer, comps, condit, &
         exist, groupa, groupb, setup%leveld - 1, setup%levelm - 1, name, names, &
         points, xrepor, rwork(qs0), rwork(qs1), signal, nsteps, found, &
         u, rwork(qv1), setup%ssabs, setup%ssage, setup%ssrel, rwork(qy0), ynorm, &
         rwork(qy1))
      if (error) go to 9017

!///  PASS REQUESTS FROM SEARCH TO THE CALLER.

      if (signal /= ' ') then
!        GO TO 4020 WHEN RETURN = 4
         return = 4
         go to 9931
      end if

!///  REACT TO THE COMPLETION OF SEARCH.

      if (found) then
!        SAVE THE LATEST SOLUTION

         psave = points
         if (setup%adapt .and. 0 < points) &
            call twcopy (points, x, rwork(qxsave))
         call twcopy &
            (groupa + comps * points + groupb, u, rwork(qusave))

!        GO TO 4030 WHEN RETURN = 5
         return = 5
         go to 9911
      else
!        RESTORE THE SOLUTION
         call twcopy &
            (groupa + comps * points + groupb, rwork(qvsave), u)
      end if
4030  continue

!///  COMPLETE STATISTICS FOR THE SEARCH BLOCK.

      call twlaps (timer(qsearc))
      total(qsearc) = total(qsearc) + timer(qsearc)
      if (grid <= gmax) &
         detail(grid, qsearc) = detail(grid, qsearc) + timer(qsearc)

!///  PRINT LEVEL 10 OR 11 ON EXIT FROM THE SEARCH BLOCK.

      if (setup%levelm == 1 .and. 0 < text) then
!        GO TO 4040 WHEN LABEL = 2
         label = 2
         go to 7010
      end if
4040  continue

!///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE SEARCH BLOCK.

      if (1 < setup%levelm) then
         if (found) then
            if (0 < text) write (text, 10015) id
         else
            if (0 < text) write (text, 10016) id
         end if
      end if

!///  BOTTOM OF THE SEARCH BLOCK.

      go to 2010

!///////////////////////////////////////////////////////////////////////
!
!     REFINE BLOCK.
!
!///////////////////////////////////////////////////////////////////////

5010  continue

!///  INITIALIZE STATISTICS ON ENTRY TO THE REFINE BLOCK.

      call twtime (timer(qrefin))

!///  PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE REFINE BLOCK.

      if (1 < setup%levelm) then
         if (0 < text) write (text, 10017) id
      end if

!///  PREPARE TO CALL REFINE.

!     SAVE THE GROUP B VALUES
      do 5020 j = 1, groupb
         rwork(qvsave - 1 + j) = u(groupa + comps * points + j)
5020  continue

      exist = .false.

!///  CALL REFINE.

5030  continue

!     SUBROUTINE REFINE &
!       (ERROR, TEXT, &
!        ACTIVE, BUFFER, COMPS, setup%leveld, setup%levelm, MARK, NEWX, PADD, PMAX, &
!        POINTS, RATIO, RATIO1, RATIO2, SIGNAL, SUCCES, TOLER0, TOLER1, &
!        TOLER2, U, VARY1, VARY2, WEIGHT, X)

      call refine &
        (error, text, &
         active, &
         buffer(groupa + 1), comps, setup%leveld - 1, setup%levelm - 1, mark, &
         found, setup%ipadd, pmax, points, ratio, rwork(qrat1), rwork(qrat2), &
         signal, satisf, setup%toler0, setup%toler1, setup%toler2, u(groupa + 1), &
         iwork(qvary1), iwork(qvary2), iwork(qvary), x)
      if (error) go to 9018

!///  SERVICE REQUESTS FROM REFINE: PASS REQUESTS TO THE CALLER.

      if (signal /= ' ') then
!        INSERT THE GROUP A AND B UNKNOWNS
         buffer(1:groupa) = u(1:groupa)
         do j = 1, groupb
            buffer(groupa + comps * points + j) = rwork(qvsave - 1 + j)
         end do

!        GO TO 5030 WHEN RETURN = 6
         return = 6
         go to 9931
      end if

!///  REACT TO THE COMPLETION OF REFINE.

      if (.not. found) go to 5110

!        COMPLETE STATISTICS FOR THE OLD GRID
         call twlaps (timer(qgrid))
         if (grid <= gmax) then
            detail(grid, qgrid) = timer(qgrid)
            detail(grid, qother) &
               = detail(grid, qgrid) - (detail(grid, qfunct) &
               + detail(grid, qjacob) + detail(grid, qsolve))
         end if

!        INITIALIZE STATISTICS FOR THE NEW GRID
         grid = grid + 1
         if (grid <= gmax) then
            call twtime (timer(qgrid))
            size(grid) = points
         end if

!        INSERT THE GROUP B VALUES
         do 5060 j = 1, groupb
            u(groupa + comps * points + j) = rwork(qvsave - 1 + j)
5060     continue

!        EXPAND THE BOUNDS
         count = groupa

         do 5080 k = 1, points
            do 5070 j = 1, comps
               rwork(qabove + count) = above(groupa + j)
               rwork(qbelow + count) = below(groupa + j)
               count = count + 1
5070        continue
5080     continue

         do 5090 j = 1, groupb
            rwork(qabove + count) = above(groupa + comps + j)
            rwork(qbelow + count) = below(groupa + comps + j)
            count = count + 1
5090     continue

!        SAVE THE LATEST SOLUTION
!        GO TO 5100 WHEN RETURN = 7
         return = 7
         go to 9911
5100     continue

5110  continue

!///  COMPLETE STATISTICS FOR THE REFINE BLOCK.

      call twlaps (timer(qrefin))
      total(qrefin) = total(qrefin) + timer(qrefin)
      if (grid <= gmax) &
         detail(grid, qrefin) = detail(grid, qrefin) + timer(qrefin)

!///  PRINT LEVEL 10 OR 11 ON EXIT FROM THE REFINE BLOCK.

      if (setup%levelm == 1 .and. 0 < text) then
         write (text, '()')
!        GO TO 5120 WHEN LABEL = 3
         label = 3
         go to 7010
      end if
5120  continue

!///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE REFINE BLOCK.

      if (1 < setup%levelm) then
         if (found) then
            if (0 < text) write (text, 10018) id
         else
            if (0 < text) write (text, 10019) id
         end if
      end if

!///  BOTTOM OF THE REFINE BLOCK.

      go to 2010

!///////////////////////////////////////////////////////////////////////
!
!     EVOLVE BLOCK.
!
!///////////////////////////////////////////////////////////////////////

6010  continue

!///  INITIALIZE STATISTICS ON ENTRY TO THE EVOLVE BLOCK.

      call twtime (timer(qtimst))
      first = .true.
      jacobs = 0
      maxcon = 0.0
      steps = step

!///  PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE EVOLVE BLOCK.

      if (1 < setup%levelm) then
         if (0 < text) write (text, 10020) id
      end if

!///  CALL EVOLVE.

6020  continue

!     SUBROUTINE EVOLVE &
!       (ERROR, TEXT, &
!        ABOVE, BELOW, BUFFER, COMPS, CONDIT, DESIRE, GROUPA, GROUPB, &
!        setup%leveld, setup%levelm, NAME, NAMES, POINTS, REPORT, S0, S1, SIGNAL, &
!        STEP, STEPS2, STRID0, STRIDE, SUCCES, TDABS, TDAGE, TDEC, &
!        TDREL, TIME, TINC, TMAX, TMIN, V0, V1, VSAVE, Y0, Y1, YNORM)

      call evolve &
        (error, text, &
         rwork(qabove), rwork(qbelow), buffer, comps, condit, desire, &
         groupa, groupb, setup%leveld - 1, setup%levelm - 1, name, names, points, &
         xrepor, rwork(qs0), rwork(qs1), signal, step, setup%steps2, setup%strid0, &
         stride, found, setup%tdabs, setup%tdage, setup%tdec, setup%tdrel, time, setup%tinc, setup%tmax, &
         setup%tmin, u, rwork(qv1), rwork(qvsave), rwork(qy0), rwork(qy1), &
         ynorm)
      if (error) go to 9019

!///  PASS REQUESTS FROM EVOLVE TO THE CALLER.

      if (signal /= ' ') then
!        GO TO 6020 WHEN RETURN = 8
         return = 8
         go to 9931
      end if

!///  REACT TO THE COMPLETION OF EVOLVE.

      if (found) then
!        SAVE THE LATEST SOLUTION
!        GO TO 6030 WHEN RETURN = 9
         return = 9
         go to 9911
      end if
6030  continue

!///  ALLOW FURTHER TIME EVOLUTION.

      allow = xrepor == qnull

!///  COMPLETE STATISTICS FOR THE EVOLVE BLOCK.

      call twlaps (timer(qtimst))
      total(qtimst) = total(qtimst) + timer(qtimst)
      if (grid <= gmax) &
         detail(grid, qtimst) = detail(grid, qtimst) + timer(qtimst)
      steps = step - steps

!///  PRINT LEVEL 10 OR 11 ON EXIT FROM THE EVOLVE BLOCK.

      if (setup%levelm == 1 .and. 0 < text) then
!        GO TO 6040 WHEN LABEL = 4
         label = 4
         go to 7010
      end if
6040  continue

!///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE EVOLVE BLOCK.

      if (1 < setup%levelm) then
         if (found) then
            if (0 < text) write (text, 10021) id
         else
            if (0 < text) write (text, 10022) id
         end if
      end if

!///  BOTTOM OF THE EVOLVE BLOCK.

      go to 2010

!///////////////////////////////////////////////////////////////////////
!
!     BLOCK TO PRINT LOG LINES.
!
!///////////////////////////////////////////////////////////////////////

7010  continue

      do 7020 j = 1, 3
         column(j) = ' '
7020  continue

      string = ' '

!     COLUMN 1: NAME OF THE TASK
      if (qtask == qentry) column(1) = '   START'
      if (qtask == qsearc) column(1) = '  SEARCH'
      if (qtask == qrefin) column(1) = '  REFINE'
      if (qtask == qtimst) column(1) = '  EVOLVE'

!     COLUMN 2: NORM OF THE STEADY STATE FUNCTION
      if (.not. found) go to 7040
!        GO TO 7030 WHEN RETURN = 10
         return = 10
         go to 9941
7030     continue
         call twnorm (groupa + comps * points + groupb, temp, buffer)
         call twlogr (column(2), temp)
7040  continue

!     COLUMN 3: LARGEST CONDITION NUMBER
      if (qtask == qsearc .or. qtask == qtimst) then
         if (maxcon /= 0.0) call twlogr (column(3), maxcon)
      end if

!     REMARK
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
         if (xrepor == qbnds .or. xrepor == qdvrg .or. &
            xrepor == qnull) then
            write (string, '(I10, A, 1P, E10.1, A)') &
               steps, ' TIME STEPS, ', stride, ' LAST STRIDE'
         else
            string = '?'
         end if
      else if (qtask == qentry .and. setup%adapt) then
         write (string, '(I10, A)') points, ' GRID POINTS'
      else if (qtask == qrefin) then
         if (found) then
            write (string, '(F10.2, A, F10.2, A, I10, A)') &
               ratio(1), ' AND ', ratio(2), ' RATIOS, ', points, &
               ' GRID POINTS'
         else
            write (string, '(F10.2, A, F10.2, A)') &
               ratio(1), ' AND ', ratio(2), ' RATIOS'
         end if
      end if

      call twsqez (length, string)
      if (0 < text) write (text, 10023) column, string (1 : length)

      go to (1110, 4040, 5120, 6040) label
      error = .true.
      go to 9020

!///////////////////////////////////////////////////////////////////////
!
!     REQUEST REVERSE COMMUNICATION.
!
!///////////////////////////////////////////////////////////////////////

!///  SAVE THE SOLUTION.

9911  continue

      call twcopy (groupa + comps * points + groupb, u, buffer)
      signal = 'SAVE'
!     GO TO 9912 WHEN ROUTE = 1
      route = 1
      return
9912  continue
      signal = ' '

      go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030) &
         return
      error = .true.
      go to 9021

!///  PRINT THE LATEST SOLUTION.

9921  continue

      call twcopy (groupa + comps * points + groupb, u, buffer)
      signal = 'SHOW'
!     GO TO 9922 WHEN ROUTE = 2
      route = 2
      return
9922  continue
      signal = ' '

      go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030) &
         return
      error = .true.
      go to 9021

!///  PASS REQUESTS FROM SEARCH, REFINE, OR EVOLVE TO THE CALLER.

9931  continue

!     IDENTIFY THE REQUEST.  THIS MUST BE SAVED TO GATHER STATISTICS
!     AT REENTRY.  THE REVERSE COMMUNICATION FLAGS WILL NOT BE SAVED
!     BECAUSE THEY ARE CLEARED AT EVERY ENTRY.
      if (signal == 'RESIDUAL') then
         qtype = qfunct
      else if (signal == 'PREPARE') then
         qtype = qjacob
      else if (signal == 'SOLVE') then
         qtype = qsolve
      else
         qtype = qother
      end if

!     COUNT THE JACOBIANS
      if (qtype == qjacob) jacobs = jacobs + 1

      call twtime (timer(qtype))

!     GO TO 9932 WHEN ROUTE = 3
      route = 3
      return
9932  continue

!     SAVE THE CONDITION NUMBER
      if (qtype == qjacob) maxcon = max (maxcon, condit)

      call twlaps (timer(qtype))
      total(qtype) = total(qtype) + timer(qtype)
      if (grid <= gmax) then
         detail(grid, qtype) = detail(grid, qtype) + timer(qtype)
         event(grid, qtype) = event(grid, qtype) + 1
      end if

      go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030) &
         return
      error = .true.
      go to 9021

!///  EVALUATE THE STEADY STATE FUNCTION.

9941  continue
      call twtime (timer(qfunct))

      call twcopy (groupa + comps * points + groupb, u, buffer)
      signal = 'RESIDUAL'
      time = .false.
!     GO TO 9942 WHEN ROUTE = 4
      route = 4
      return
9942  continue
      signal = ' '

      call twlaps (timer(qfunct))
      total(qfunct) = total(qfunct) + timer(qfunct)
      if (grid <= gmax) then
         detail(grid, qfunct) = detail(grid, qfunct) + timer(qfunct)
         event(grid, qfunct) = event(grid, qfunct) + 1
      end if

      go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030) &
         return
      error = .true.
      go to 9021

!///////////////////////////////////////////////////////////////////////
!
!     INFORMATIVE MESSAGES.
!
!///////////////////////////////////////////////////////////////////////

10001 format(/1X, a9, a, ' (TWO POINT BOUNDARY VALUE PROBLEM) SOLVER,' &
             /10X, 'VERSION ', a, &
             ' OF APRIL 1998 BY DR. JOSEPH F. GRCAR.')
10002 format(/1X, a9, a)
10003 format(3(/10X, a35)/)
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
80001 format('(', a, ' ', i10, ')')
80002 format(10X, 1p, e10.2, 2X, e10.2, 3X, a)
80003 format(10X, '  ... MORE')

!///////////////////////////////////////////////////////////////////////
!
!     ERROR MESSAGES.
!
!///////////////////////////////////////////////////////////////////////

!     return

9001  if (0 < text) write (text, 01) id, route
      if (.not. mess) return

9002  if (0 < text) then
         call twlast (length, versio)
         write (text, 02) id, versio (1 : length), vnmbr(vnmbrs)
         do j = vnmbrs - 1, 1, - 1
            write (text, '(10X, A, A)') &
               ' CAN REPLACE:  double precision VERSION ', vnmbr(j)
         end do
      end if
      if (.not. mess) return

9004  if (0 < text) write (text, 04) id, cntrls, count
      if (.not. mess) return

9005  if (0 < text) write (text, 05) id, setup%leveld, setup%levelm
      if (.not. mess) return

9006  if (0 < text) write (text, 06) id, &
         comps, points, groupa, groupb
      if (.not. mess) return

9007  if (0 < text) write (text, 07) id, comps, points
      if (.not. mess) return

9008  if (0 < text) write (text, 08) id, &
         comps, points, groupa, groupb, groupa + comps * points + groupb
      if (.not. mess) return

9009  if (0 < text) write (text, 09) id, &
         names, comps, groupa, groupb, groupa + comps + groupb
      if (.not. mess) return

9010  if (0 < text) write (text, 10) id, points, pmax
      if (.not. mess) return

9011  if (0 < text) then
         write (text, 11) id, &
            groupa, groupb, comps, groupa + comps + groupb, count
         count = 0
         do 8010 j = 1, groupa + comps + groupb
            if (.not. (below(j) < above(j)) .or. mess) then
               count = count + 1
               if (count <= lines) then
                  if (names == comps + groupa + groupb) then
                     ctemp1 = name(j)
                  else
                     ctemp1 = ' '
                  end if
                  call twsqez (len1, ctemp1)

                  if (j <= groupa) then
                     write (ctemp2, 80001) 'A', j
                  else if (j <= groupa + comps) then
                     write (ctemp2, 80001) 'C', j - groupa
                  else
                     write (ctemp2, 80001) 'B', j - groupa - comps
                  end if
                  call twsqez (len2, ctemp2)

                  if (ctemp1 == ' ') then
                     string = ctemp2
                     length = len2
                  else if (len1 + 2 + len2 <= 40) then
                     string = ctemp1 (1 : len1) // '  ' // ctemp2
                     length = len1 + 2 + len2
                  else if (len1 + 1 + len2 <= 40) then
                     string = ctemp1 (1 : len1) // ' ' // ctemp2
                     length = len1 + 1 + len2
                  else
                     len1 = 40 - len2 - 4
                     string = ctemp1 (1 : len1) // '... ' // ctemp2
                     length = 40
                  end if

                  write (text, 80002) &
                     below(j), above(j), string (1 : length)
               end if
            end if
8010     continue
         if (lines < count) write (text, 80003)
      end if
      if (.not. mess) return

9012  if (0 < text) write (text, 12) id
      if (.not. mess) return

9013  if (0 < text) write (text, 13) id, &
         isize, rsize, ilast, rlast
      if (.not. mess) return

9014  if (0 < text) write (text, 14) id
      if (.not. mess) return

9015  if (0 < text) write (text, 15) id
      if (.not. mess) return

9016  if (0 < text) write (text, 16) id
      if (.not. mess) return

9017  if (0 < text) write (text, 17) id
      if (.not. mess) return

9018  if (0 < text) write (text, 18) id
      if (.not. mess) return

9019  if (0 < text) write (text, 19) id
      if (.not. mess) return

9020  if (0 < text) write (text, 20) id, label
      if (.not. mess) return

9021  if (0 < text) write (text, 21) id, return
      if (.not. mess) return

      stop

        ! Error messages
        01 format(/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
                   //10X, i10, '  ROUTE')
        02 format(/1X, a9, 'ERROR.  THE CALLING PROGRAM EXPECTS A VERSION OF' &
                    /10X, 'TWOPNT NOT COMPATIBLE WITH THIS VERSION.' &
                   //10X, '     EXPECTS:  ', a &
                   //10X, 'THIS VERSION:  double precision VERSION ', a)
        04 format(/1X, a9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.' &
                   //10X, i10, '  CONTROLS' &
                    /10X, i10, '  COUNTED')
        05 format(/1X, a9, 'ERROR.  THE PRINTING LEVELS ARE OUT OF ORDER.' &
                    /10X, 'LEVELD CANNOT EXCEED LEVELM.' &
                   //10X, i10, '  LEVELD, FOR SOLUTIONS' &
                    /10X, i10, '  LEVELM, FOR MESSAGES')
        06 format(/1X, a9, 'ERROR.  NUMBERS OF ALL TYPES OF UNKNOWNS MUST BE AT' &
                    /10X, 'LEAST ZERO.' &
                   //10X, i10, '  COMPS, COMPONENTS' &
                    /10X, i10, '  POINTS' &
                    /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                    /10X, i10, '  GROUPB, GROUP B UNKNOWNS')
        07 format(/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE' &
                    /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE.' &
                   //10X, i10, '  COMPS, COMPONENTS' &
                    /10X, i10, '  POINTS')
        08 format(/1X, a9, 'ERROR.  TOTAL UNKNOWNS MUST BE POSITIVE.' &
                   //10X, i10, '  COMPS, COMPONENTS' &
                    /10X, i10, '  POINTS' &
                    /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                    /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                    /10X, i10, '  TOTAL NUMBER')
        09 format(/1X, a9, 'ERROR.  THE NUMBER OF NAMES IS WRONG.' &
                   //10X, i10, '  NAMES' &
                   //10X, i10, '  COMPS, COMPONENTS' &
                    /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
                    /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
                    /10X, i10, '  TOTAL NUMBER')
        10 format(/1X, a9, 'ERROR.  THERE ARE TOO MANY POINTS.' &
                   //10X, i10, '  POINTS' &
                    /10X, i10, '  PMAX, LIMIT ON POINTS')
        11 format(/1X, a9, 'ERROR.  THE LOWER AND UPPER BOUNDS ON SOME UNKNOWNS' &
                    /10X, 'ARE OUT OF ORDER.' &
                   //10X, i10, '  GROUP A UNKNOWNS (A)' &
                    /10X, i10, '  GROUP B UNKNOWNS (B)' &
                    /10X, i10, '  COMPONENTS AT POINTS (C)' &
                    /10X, i10, '  TOTAL TYPES OF UNKNOWNS' &
                    /10X, i10, '  NUMBER OF BOUNDS OUT OF ORDER' &
                   //10X, '     LOWER       UPPER' &
                    /10X, '     BOUND       BOUND   UNKNOWN'/)
        12 format(/1X, a9, 'ERROR.  TWGRAB FAILS.')
        13 format(/1X, a9, 'ERROR.  ONE OR BOTH WORK SPACES ARE TOO SMALL.' &
                   //25X, '   INTEGER        REAL' &
                   //10X, ' PRESENT SIZE', 2i12 &
                    /10X, 'REQUIRED SIZE', 2i12)
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



      end subroutine twopnt

      ! Perform automatic grid selection
      subroutine refine(error, text, &
                        active, buffer, comps, leveld, levelm, mark, newx, padd, pmax, &
                        points, ratio, ratio1, ratio2, signal, succes, toler0, toler1, &
                        toler2, u, vary1, vary2, weight, x)

      character(len=*), parameter :: id = 'REFINE:  '

      ! set true to print examples of all messages.
      logical, parameter :: mess = .false.

      character signal*(*), word*80
      double precision :: buffer, differ, left, length, lower, maxmag, mean, range, &
        ratio, ratio1, ratio2, right, temp, temp1, temp2, toler0, toler1, toler2, u, upper, x
      integer :: act, comps, count, former, itemp, j, k, least, leveld, levelm, &
        more, most, new, old, padd, pmax, points, route, signif, text, total, vary1, vary2, weight
      intrinsic :: abs, max, min
      logical   :: active, error, mark, newx, succes

      dimension active(comps), buffer(comps * pmax), mark(pmax), ratio(2), ratio1(pmax), ratio2(pmax), &
                u(comps, pmax), vary1(pmax),vary2(pmax), weight(pmax), x(pmax)

!///  save local variables during returns for reverse communication.

      save

!///////////////////////////////////////////////////////////////////////
!
!     prologue.
!
!///////////////////////////////////////////////////////////////////////

!///  every-time initialization.


!     turn off all completion status flags.
      error  = .false.
      newx   = .false.
      succes = .false.

!///  if this is a return call, then continue where the program paused.

      if (signal/=' ') then
         go to (4060, 5040) route
         error = .true.

         if (0<text) write (text, 101) id, route
         return
      end if

      ! write all messages.
      if (mess .and. 0<text) then
         route = 0
         write (text, 1) id
         write (text, 2) id
         write (text, 4) id
         write (text, 5) id
         write (text, 10) id
         write (text, 101) id, route
         write (text, 102) id, comps, points
         write (text, 103) id, padd
         write (text, 104) id, points, pmax
         write (text, 105) id
         write (text, 106) id, toler0
         write (text, 107) id, toler1, toler2
         write (text, 108) id
         write (text, 109) id
         stop
      end if

      ! Levelm printing.
      if (0<levelm .and. 0<text) write (text, 1) id

      ! Check the arguments.
      error = .not. (1 <= comps .and. 2 <= points)
      if (error) then
           if (0<text) write (text, 102) id, comps, points
           return
      end if

      error = .not. (0 <= padd)
      if (error) then
           if (0<text) write (text, 103) id, padd
           return
      endif

      error = .not. (points <= pmax)
      if (error) then
           if (0<text) write (text, 104) id, points, pmax
           return
      end if

      count = 0
      do 1010 j = 1, comps
         if (active(j)) count = count + 1
1010  continue
      error = .not. (1 <= count)
      if (error) then
           if (0<text) write (text, 105) id
           return
      end if

      error = .not. toler0>=zero
      if (error) then
          if (0<text) write (text, 106) id, toler0
          return
      end if

      error = .not. (zero <= toler1 .and. toler1 <= one &
              .and.  zero <= toler2 .and. toler2 <= one)
      if (error) then
          if (0<text) write (text, 107) id, toler1, toler2
          return
      end if

      count = 0
      do 1020 k = 1, points - 1
         if (x(k)<x(k + 1)) count = count + 1
1020  continue
      error = .not. (count == 0 .or. count == points - 1)
      if (error) then
          if (0<text) write (text, 108) id
          return
      end if

!///////////////////////////////////////////////////////////////////////
!
!     at each interval, count the active, significant components that
!     vary too greatly.
!
!///////////////////////////////////////////////////////////////////////

      act = 0
      signif = 0

      do 2010 k = 1, points
         mark(k) = .false.
         ratio1(k) = zero
         ratio2(k) = zero
         vary1(k) = 0
         vary2(k) = 0
2010  continue

!///  top of the loop over the components.

      do 2060 j = 1, comps
         if (active(j)) then
            act = act + 1

!///  find the range and maximum magnitude of the component.

            lower = u(j, 1)
            upper = u(j, 1)
            do 2020 k = 2, points
               lower = min (lower, u(j, k))
               upper = max (upper, u(j, k))
2020        continue
            range = upper - lower
            maxmag = max (abs (lower), abs (upper))

!///  decide whether the component is significant.

            if (abs (range)>max (toler0, toler0 * maxmag)) then
               signif = signif + 1

!///  at each interval, see whether the component'S CHANGE EXCEEDS SOME
!///  fraction of the component'S GLOBAL CHANGE.

               do 2030 k = 1, points - 1
                  differ = abs (u(j, k + 1) - u(j, k))
                  if (zero<range) ratio1(k) = max (ratio1(k), differ / range)
                  if (toler1 * range<differ) vary1(k) = vary1(k) + 1
2030           continue

!///  find the global change of the component'S DERIVATIVE.

               temp = (u(j, 2) - u(j, 1)) / (x(2) - x(1))
               lower = temp
               upper = temp
               do 2040 k = 2, points - 1
                  temp  = (u(j, k+1) - u(j, k)) / (x(k+1) - x(k))
                  lower = min(lower, temp)
                  upper = max(upper, temp)
2040           continue
               range = upper - lower

!///  at each interior point, see whether the derivative'S CHANGE
!///  exceeds some fraction of the derivative'S GLOBAL CHANGE.

               right = (u(j, 2) - u(j, 1)) / (x(2) - x(1))
               do 2050 k = 2, points - 1
                  left = right
                  right = (u(j, k + 1) - u(j, k)) / (x(k + 1) - x(k))
                  differ = abs (left - right)
                  if (zero<range) ratio2(k) = max (ratio2(k), differ / range)
                  if (toler2 * range < differ) vary2(k) = vary2(k) + 1
2050           continue

!///  bottom of the loop over the components.

            end if
         end if
2060  continue

      ! save the maximum ratios.
      ratio(1) = max(zero,maxval(ratio1(1:points-1),1))
      ratio(2) = max(zero,maxval(ratio2(2:points-1),1))

!///////////////////////////////////////////////////////////////////////
!
!     select the intervals to halve.
!
!///////////////////////////////////////////////////////////////////////

!///  weight the intervals in which variations that are too large occur.

      most = 0
      do 3010 k = 1, points - 1
         weight(k) = vary1(k)
         if (1<k) weight(k) = weight(k) + vary2(k)
         if (k<points - 1) weight(k) = weight(k) + vary2(k + 1)
         if (0<weight(k)) most = most + 1
3010  continue

!///  sort the weights.

      do 3030 k = 1, points - 1
         do 3020 j = k + 1, points - 1
            if (weight(j)>weight(k)) then
               itemp = weight(j)
               weight(j) = weight(k)
               weight(k) = itemp
            end if
3020     continue
         if (weight(k) == 0) go to 3040
3030  continue
3040  continue

!///  find the least weight of intervals to halve.

      more = max (0, min (most, padd, pmax - points))
      if (0<more) then
         least = weight(more)
      else
         least = 1 + weight(1)
      end if

!///  reconstruct the weights.

      do 3050 k = 1, points - 1
         weight(k) = vary1(k)
         if (1<k) weight(k) = weight(k) + vary2(k)
         if (k<points - 1) weight(k) = weight(k) + vary2(k + 1)
3050  continue

!///  mark the intervals to halve.

      count = 0
      do 3060 k = 1, points - 1
         if (count<more .and. least <= weight(k)) then
            count = count + 1
            mark(k) = .true.
         end if
3060  continue

! hack  if one point is marked, mark them all
!      if (count>0) then
!         count = 0
!         do k = 1, points - 1
!            count = count + 1
!            mark(k) = .true.
!         enddo
!      endif

      more = count

!///////////////////////////////////////////////////////////////////////
!
!     halve the intervals, if any.
!
!///////////////////////////////////////////////////////////////////////

!///  form the total number of points in the new grid.

      total = points + more
      if (0 == more) go to 4070

!///  top of the block to create the new grid.  check the order.

      count = 0
      length = abs (x(points) - x(1))
      do 4010 k = 1, points - 1
         if (mark(k)) then
            mean = 0.5 * (x(k) + x(k + 1))
            if (.not. ((x(k)<mean .and. mean<x(k + 1)) .or. (x(k + 1)<mean .and. mean<x(k)))) &
               count = count + 1
         end if
4010  continue
      error = .not. (count == 0)
      if (error) then
          if (0<text) write (text, 109) id
          return
      end if

!///  add the new points.  interpolate x and the bounds.

      new = total
      do 4040 old = points, 2, - 1
         x(new) = x(old)
         u(:,new) = u(:,old)
         do 4020 j = 1, comps
            u(j, new) = u(j, old)
4020     continue
         new = new - 1

         if (mark(old - 1)) then
            x(new) = 0.5 * (x(old) + x(old - 1))
            do 4030 j = 1, comps
               u(j, new) = 0.5 * (u(j, old) + u(j, old - 1))
4030        continue
            new = new - 1
         end if
4040  continue

!///  mark the new points.

      new = total
      do 4050 old = points, 2, - 1
         mark(new) = .false.
         new = new - 1
         if (mark(old - 1)) then
            mark(new) = .true.
            new = new - 1
         end if
4050  continue
      mark(new) = .false.

!///  update the number of points.

      former = points
      points = total

!///  allow the user to update the solution.

      call twcopy (comps * points, u, buffer)
      signal = 'UPDATE'
!     go to 4060 when route = 1
      route = 1
      return
4060  continue
      signal = ' '
      call twcopy (comps * points, buffer, u)

!///  bottom of the block to create a new grid.

4070  continue

!///////////////////////////////////////////////////////////////////////
!
!     epilogue.
!
!///////////////////////////////////////////////////////////////////////

!///  print.

      if (0<levelm .and. 0<text) then
         temp1 = ratio1(1)
         do 5010 k = 2, former - 1
            temp1 = max (temp1, ratio1(k))
5010     continue

         temp2 = ratio2(2)
         do 5020 k = 3, former - 1
            temp2 = max (temp2, ratio2(k))
5020     continue

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
               do 5030 k = 1, points
                  if (.not. mark(k)) then
                     old = old + 1
!                     dx = zerod0
                     if (1<k) then
!                        dx = x(k)-x(k-1)
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

                     if (1<k .and. k<points) then
                        if (vary2(old)/=zero) then
                           write (word, '(F4.2, I4)') ratio2(old), vary2(old)
                        else
                           write (word, '(F4.2, I4)') ratio2(old)
                        end if
                        write (text, 9) k, x(k), word
                     else
                        write (text, 9) k, x(k)
                     end if
                  end if
5030           continue
            end if
         end if

         if (0<leveld .and. 0<more) then
            write (text, 10) id
            call twcopy (comps * points, u, buffer)
            signal = 'SHOW'
!           go to 5040 when route = 2
            route = 2
            return
         end if
      end if
5040  continue
      signal = ' '

!///  set the completion status flags.

      newx = 0<more
      succes = 0 == most

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

      stop

      return
      end subroutine refine


end module twopnt_core

