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
module twopnt
    use iso_fortran_env, only: real64
    implicit none (type, external)
    private

    integer, parameter :: RK = real64

    logical, parameter :: DEBUG = .true.

    contains

      ! *******************************************************************************************************
      ! MATH
      ! *******************************************************************************************************

      ! SOLVE A SYSTEM OF LINEAR EQUATIONS USING THE MATRIX PREPARED BY TWPREP.
      subroutine twsolv(error, text, a, asize, buffer, comps, groupa, groupb, pivot, points)

          integer , intent(in) :: asize,text,comps,groupa,groupb,points
          integer , intent(in) :: pivot(groupa + comps * points + groupb)
          real(RK), intent(in) :: a(asize)
          real(RK), intent(inout) :: buffer(groupa + comps * points + groupb)
          logical , intent(out) :: error

          ! Local variables
          integer   :: n,width
          intrinsic :: max
          character(*), parameter :: id = 'TWSOLV:  '

          ! SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
          logical, parameter :: mess = .false.

          !///////////////////////////////////////////////////////////////////////
          !
          !     (1) PROLOGUE.
          !
          !///////////////////////////////////////////////////////////////////////

          !///  WRITE MESSAGES only.
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

          !///////////////////////////////////////////////////////////////////////
          !
          !     (2) SCALE AND SOLVE THE EQUATIONS.
          !
          !///////////////////////////////////////////////////////////////////////

          buffer(1:n) = buffer(1:n) * a(1:n)

          call twgbsl(a(n + 1), 3 * width + 1, n, width, width, pivot, buffer)

          !///////////////////////////////////////////////////////////////////////
          !
          !     ERROR MESSAGES.
          !
          !///////////////////////////////////////////////////////////////////////
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

          !///  EXIT.
          stop

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

end module twopnt

