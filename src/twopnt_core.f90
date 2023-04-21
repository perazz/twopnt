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
          double precision, intent(in) :: a(asize)
          double precision, intent(inout) :: buffer(groupa + comps * points + groupb)
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

      ! *******************************************************************************************************
      ! UTILITIES
      ! *******************************************************************************************************

      ! Obtain computing time in seconds.
      subroutine twtime(timer)
         double precision, intent(out) :: timer
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

end module twopnt_core

!     cvs $revision: 1.1.1.1 $ reposited $date: 2006/05/26 19:09:34 $

!///////////////////////////////////////////////////////////////////////

!
!///////////////////////////////////////////////////////////////////////

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

      implicit complex (a - p, r - z), integer (q)
      character &
         cword*80, header*80, id*9, jword*80, name*(*), remark*80, &
         signal*(*), yword*80
!**** PRECISION > DOUBLE
      double precision    above, below, buffer, change, condit, csave, dummy, high, low, &
         s0, s1, strid0, stride, tdabs, tdec, tdrel, tinc, tmax, tmin, &
         v0, v1, vsave, y0, y1, ynorm
      external &
         search, twcopy, twlogr, twnorm
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
         name(names), s0(groupa + comps * points + groupb), &
         v0(groupa + comps * points + groupb), &
         v1(groupa + comps * points + groupb), &
         vsave(groupa + comps * points + groupb), &
         y0(groupa + comps * points + groupb)

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

      implicit complex (a - z)

      character &
         column*16, ctemp1*80, ctemp2*80, header*80, id*9, name*(*), &
         signal*(*), string*80
!**** PRECISION > DOUBLE
      double precision    above, abs0, abs1, below, buffer, condit, deltab, deltad, rel0, &
         rel1, s0, s0norm, s1, s1norm, sj, temp, v0, v1, value, vj, &
         xxabs, xxrel, y0, y0norm, y1, y1norm, zero
      external &
         twcopy, twlogr, twnorm, twsqez
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
      do 1010 j = 1, groupa + comps * points + groupb
         if (.not. (below(j) < above(j))) count = count + 1
1010  continue
      error = count /= 0
      if (error) go to 9004

      count = 0
      do 1020 j = 1, groupa + comps * points + groupb
         if (.not. (below(j) <= v0(j) .and. v0(j) <= above(j))) &
            count = count + 1
1020  continue
      error = count /= 0
      if (error) go to 9005

      error = .not. (0.0 <= xxabs .and. 0.0 <= xxrel)
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
      end
      subroutine twcopy (n, x, y)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWCOPY
!
!     COPY ONE VECTOR TO ANOTHER.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      integer j, n
!**** PRECISION > DOUBLE
      double precision    x, y

      dimension x(n), y(n)

      do 0100 j = 1, n
         y(j) = x(j)
0100  continue

      return
      end subroutine twcopy
      subroutine tweps (eps)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWEPS
!
!     FIND MACHINE EPSILON.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
!**** PRECISION > DOUBLE
      double precision    eps, value
      integer scrtch
!**** MACHINE EPSILON > COMPUTED
!      LOGICAL
!     +   SAME &
!**** end MACHINE EPSILON > COMPUTED

      parameter (scrtch = 98)

!///  IEEE STANDARD &
!***  PRECISION > DOUBLE
      value = 1.1102230246251565D-16
!**** END PRECISION > DOUBLE &
!**** PRECISION > SINGLE
!      VALUE = 5.9604645E-08 &
!**** END PRECISION > SINGLE &

!***  MACHINE EPSILON > IEEE STANDARD
      eps = value
!**** END MACHINE EPSILON > IEEE STANDARD

!///  COMPUTED &

!***  MACHINE EPSILON > COMPUTED
!      OPEN (ACCESS = 'SEQUENTIAL', FORM = 'UNFORMATTED',
!     +   STATUS = 'SCRATCH', UNIT = SCRTCH)
!
!      EPS = 1
!1010  CONTINUE
!      EPS = 0.5 * EPS
!
!      VALUE = 1 + EPS
!
!      REWIND (SCRTCH)
!      WRITE (SCRTCH) VALUE
!
!      REWIND (SCRTCH)
!      READ (SCRTCH) VALUE
!
!      SAME = 1 == VALUE
!
!      IF (.NOT. SAME) GO TO 1010
!
!      CLOSE (UNIT = SCRTCH) &
!**** END MACHINE EPSILON > COMPUTED

!///  EXIT.

      return
      end subroutine tweps
      subroutine twgbco (a, lda, n, lower, upper, pivot, rcond, z)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWGBCO
!
!     FACTOR A BANDED MATRIX AND ESTIMATE THE RECIPROCAL OF ITS
!     CONDITION NUMBER.  BASED ON _GBCO FROM THE LINPACK LIBRARY.
!
!///////////////////////////////////////////////////////////////////////

      double precision dsum
!**** PRECISION > DOUBLE
      double precision a, anorm, ek, rcond, s, sm, sum, t, wk, wkm, ynorm, z
      external &
         twgbfa
      integer &
         first, info, j, jdiag, ju, k, last, lda, lower, mm, n, pivot, &
         upper
      intrinsic &
         abs, dble, max, min, sign

      dimension &
         a(lda,n), pivot(n), z(n)

      jdiag = lower + upper + 1

!///  COMPUTE THE 1-NORM OF A

      anorm = 0.0
      do 1020 k = 1, n
         first = max (lower + 1, jdiag + 1 - k)
         last = min (jdiag + lower, jdiag + n - k)
         sum = 0.0
         do 1010 j = first, last
            sum = sum + abs (a(j, k))
1010     continue
         anorm = max (anorm, sum)
1020  continue

!///  FACTOR A

      call twgbfa (a, lda, n, lower, upper, pivot, info)

!///  SOLVE TRANSPOSE(U) * W = E

      ek = 1.0
      do 2010 j = 1, n
         z(j) = 0.0
2010  continue

      ju = 0
      do 2050 k = 1, n
         if (z(k) /= 0.0) ek = sign (ek, - z(k))

         if (abs (ek - z(k)) > abs (a(jdiag, k))) then
            s = abs (a(jdiag, k)) / abs (ek - z(k))
            ek = s * ek
            do 2020 j = 1, n
               z(j) = s * z(j)
2020        continue
         end if

         wk = ek - z(k)
         wkm = - ek - z(k)
         s = abs (wk)
         sm = abs (wkm)
         if (a(jdiag, k) /= 0.0) then
            wk = wk / a(jdiag, k)
            wkm = wkm / a(jdiag, k)
         else
            wk = 1.0
            wkm = 1.0
         end if

         ju = min (max (ju, upper + pivot(k)), n)
         mm = jdiag
         if (k + 1 <= ju) then
            do 2030 j = k + 1, ju
               mm = mm - 1
               sm = sm + abs (z(j) + wkm * a(mm, j))
               z(j) = z(j) + wk * a(mm, j)
               s = s + abs (z(j))
2030        continue

            if (s < sm) then
               t = wkm - wk
               wk = wkm
               mm = jdiag
               do 2040 j = k + 1, ju
                  mm = mm - 1
                  z(j) = z(j) + t * a(mm, j)
2040           continue
            end if
         end if

         z(k) = wk
2050  continue

      sum = 0.0
      do 2060 j = 1, n
         sum = sum + abs (z(j))
2060  continue
      s = 1.0 / sum

      do 2070 j = 1, n
         z(j) = s * z(j)
2070  continue

!///  SOLVE TRANSPOSE(L) * Y = W

      do 3030 k = n, 1, - 1
         dsum = 0.0
         do 3010 j = 1, min (lower, n - k)
            dsum = dsum + dble (a(jdiag + j, k)) * dble (z(k + j))
3010     continue
         z(k) = z(k) + dsum

         if (1.0 < abs (z(k))) then
            s = 1.0 / abs (z(k))
            do 3020 j = 1, n
               z(j) = s * z(j)
3020        continue
         end if

         j = pivot(k)
         t = z(j)
         z(j) = z(k)
         z(k) = t
3030  continue

      sum = 0.0
      do 3040 j = 1, n
         sum = sum + abs (z(j))
3040  continue
      s = 1.0 / sum

      do 3050 j = 1, n
         z(j) = s * z(j)
3050  continue

      ynorm = 1.0

!///  SOLVE L * V = Y

      do 4030 k = 1, n
         j = pivot(k)
         t = z(j)
         z(j) = z(k)
         z(k) = t

         do 4010 j = 1, min (lower, n - k)
            z(k + j) = t * a(jdiag + j, k) + z(k + j)
4010     continue

         if (1.0 < abs (z(k))) then
            s = 1.0 / abs (z(k))
            do 4020 j = 1, n
               z(j) = s * z(j)
4020        continue

            ynorm = s * ynorm
         end if
4030  continue

      sum = 0.0
      do 4040 j = 1, n
         sum = sum + abs (z(j))
4040  continue
      s = 1.0 / sum

      do 4050 j = 1, n
         z(j) = s * z(j)
4050  continue

      ynorm = s * ynorm

!///  SOLVE U * Z = W

      do 5030 k = n, 1, - 1
         if (abs (z(k)) > abs (a(jdiag, k))) then
            s = abs (a(jdiag, k)) / abs (z(k))
            do 5010 j = 1, n
               z(j) = s * z(j)
5010        continue
            ynorm = s*ynorm
         end if

         if (a(jdiag, k) /= 0.0) then
            z(k) = z(k) / a(jdiag, k)
         else
            z(k) = 1.0
         end if

         t = - z(k)
         do 5020 j = 1, min (k, jdiag) - 1
            z(k - j) = t * a(jdiag - j, k) + z(k - j)
5020     continue
5030  continue

      sum = 0.0
      do 5040 j = 1, n
         sum = sum + abs (z(j))
5040  continue
      s = 1.0 / sum

      do 5050 j = 1, n
         z(j) = s * z(j)
5050  continue

      ynorm = s * ynorm

!///  FORM RCOND

      if (anorm /= 0.0) then
         rcond = ynorm / anorm
      else
         rcond = 0.0
      end if

!///  EXIT

      return
      end subroutine twgbco
      subroutine twgbfa (a, lda, n, lower, upper, pivot, info)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWGBFA
!
!     FACTOR A BANDED MATRIX FOR TWGBCO. BASED ON _GBFA FROM THE LINPACK
!     LIBRARY.
!
!/////////////////////////////////////////////////////////////////////// &

!***  PRECISION > DOUBLE
      double precision    a, t, value
      integer &
         i, info, j, jk, k, lda, lower, n, pack, pdiag, pivot, pjk, &
         save, true, upper
      intrinsic &
         abs, max, min

      dimension &
         a(lda, n), pivot(n)

!///  STATEMENT FUNTIONS

!     PACKED ROW POSITION OF ENTRY J IN COLUMN K
      pack (j, k) = j - k + pdiag

!     TRUE ROW POSITION OF ENTRY J PACKED IN COLUMN K
      true (j, k) = j - pdiag + k

!///  INITIALIZE

!     PACKED ROW POSITION OF THE DIAGONAL
      pdiag = lower + upper + 1

      info = 0

!///  TOP OF THE LOOP OVER COLUMNS

      do 2060 k = 1, n

!///  INITIALIZE THE FILL-IN SPACE

         do 2010 i = 1, lower
            a(i, k) = 0.0
2010     continue

!///  LOOP OVER THE PREVIOUS COLUMNS

         do 2030 j = max (1, k - lower - upper), k - 1
            pjk = pack (pivot(j), k)
            jk = pack (j, k)
            t = a(pjk, k)
            if (pjk /= jk) then
               a(pjk, k) = a(jk, k)
               a(jk, k) = t
            end if

            if (t /= 0.0) then
               do 2020 i = 1, min (lower, n - j)
                  a(jk + i, k) = t * a(pdiag + i, j) + a(jk + i, k)
2020           continue
            end if
2030     continue

!///  FIND THE PIVOT

         save = pdiag
         value = abs (a(pdiag, k))
         do 2040 i = pdiag + 1, pdiag + min (lower, n - k)
            if (value < abs (a(i, k))) then
               save = i
               value = abs (a(i, k))
            end if
2040     continue
         pivot(k) = true (save, k)

!///  INTERCHANGE IF NECESSARY

         if (save /= pdiag) then
            t = a(save, k)
            a(save, k) = a(pdiag, k)
            a(pdiag, k) = t
         end if

!///  SCALE THE LOWER COLUMN

         if (a(save, k) /= 0.0) then
            t = - 1.0 / a(pdiag, k)
            do 2050 i = pdiag + 1, pdiag + min (lower, n - k)
               a(i, k) = t * a(i, k)
2050        continue
         else
            info = k
         end if

!///  BOTTOM OF THE LOOP OVER COLUMNS

2060  continue

!///  THE FINAL COLUMN IS TRIVIAL

      pivot(n) = n
      if (a(pdiag, n) == 0.0) info = n

!///  EXIT

      return
      end subroutine twgbfa

      subroutine twgrab (error, last, first, number)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWGRAB
!
!     RESERVE SPACE IN AN ARRAY.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

      integer &
         first, last, number
      intrinsic &
         max
      logical &
         error

!///  CHECK THE ARGUMENTS.

      error = .not. (0 <= last)
      if (error) go to 99999

      error = .not. (0 <= number)
      if (error) go to 99999

!///  GRAB THE SPACE.

      first = last + 1
      last = last + max (1, number)

!///  EXIT.

99999 continue
      return
      end subroutine twgrab
      subroutine twinit (error, text, force)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWINIT
!
!     INITIALIZE THE CONTROLS.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character id*9
!**** PRECISION > DOUBLE
      double precision rvalue
      integer cntrls, count, ivalue, text
      logical error, first, force, lvalue, mess

      parameter (id = 'TWINIT:  ')
      parameter (cntrls = 22)

      dimension ivalue(cntrls), lvalue(cntrls), rvalue(cntrls)

      common / twcomi / ivalue
      common / twcoml / lvalue
      common / twcomr / rvalue

!     THE GNU F77 COMPILER REQUIRES THE SAVE TO PRECEED THE DATA

      save first

      data first / .true. /

!///  WRITE ALL MESSAGES.

!     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

      if (mess .and. 0 < text) go to 9001

!///  TOP OF THE BLOCK TO SET THE CONTROLS.

      if (first .or. force) then
         first = .false.

!///  SET THE CONTROLS.

      count = 0

!     ADAPT

      count = count + 1
      lvalue(count) = .false.

!     LEVELD

      count = count + 1
      ivalue(count) = 1

!     LEVELM

      count = count + 1
      ivalue(count) = 1

!     PADD

      count = count + 1
      lvalue(count) = .false.

!     SSABS

      count = count + 1
      rvalue(count) = 1.0E-9

!     SSAGE

      count = count + 1
      ivalue(count) = 10

!     SSREL

      count = count + 1
      rvalue(count) = 1.0E-6

!     STEADY

      count = count + 1
      lvalue(count) = .true.

!     STEPS0

      count = count + 1
      ivalue(count) = 0

!     STEPS1

      count = count + 1
      ivalue(count) = 200

!     STEPS2

      count = count + 1
      ivalue(count) = 100

!     STRID0

      count = count + 1
      rvalue(count) = 1.0E-4

!     TDABS

      count = count + 1
      rvalue(count) = 1.0E-9

!     TDAGE

      count = count + 1
      ivalue(count) = 20

!     TDEC

      count = count + 1
      rvalue(count) = 3.1623

!     TDREL

      count = count + 1
      rvalue(count) = 1.0E-6

!     TINC

      count = count + 1
      rvalue(count) = 10.0

!     TMAX

      count = count + 1
      rvalue(count) = 1.0E-2

!     TMIN

      count = count + 1
      rvalue(count) = 1.0E-20

!     TOLER0

      count = count + 1
      rvalue(count) = 1.0E-9

!     TOLER1

      count = count + 1
      rvalue(count) = 0.2

!     TOLER2

      count = count + 1
      rvalue(count) = 0.2

!///  BOTTOM OF THE BLOCK TO SET THE CONTROLS.

         error = .not. (count == cntrls)
         if (error) go to 9001
      end if

!///  ERROR MESSAGES.

      go to 99999

9001  if (0 < text) write (text, 99001) id, cntrls, count
      if (.not. mess) go to 99999

99001 format &
        (/1X, a9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.' &
       //10X, i10, '  CONTROLS' &
        /10X, i10, '  COUNTED')

!///  EXIT.

      stop
99999 continue
      return
      end subroutine twinit
      subroutine twlaps (timer)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWLAPS
!
!     OBTAIN ELAPSED COMPUTING TIME IN SECONDS.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      external twtime
!**** PRECISION > DOUBLE
      double precision    temp, timer

      call twtime (temp)
      timer = temp - timer

      return
      end subroutine twlaps
      subroutine twlast (length, string)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWLAST
!
!     FIND THE LAST NONBLANK CHARACTER IN A STRING.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character &
         string*(*)
      integer &
         j, length
      intrinsic &
         len

      do 0100 j = len (string), 1, - 1
         if (string(j : j) /= ' ') then
            length = j
            go to 0200
         end if
0100  continue
      length = 1
0200  continue

      return
      end subroutine twlast
      subroutine twlogr (string, value)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWLOGR
!
!     WRITE A COMMON LOGARITHM TO A CHARACTER STRING.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

      character string*(*)
!**** PRECISION > DOUBLE
      double precision    value
      intrinsic len, log10

      if (6 <= len (string)) then
         if (value < 0.0) then
            string = ' '
         else if (value == 0.0) then
            string = '  ZERO'
         else
            write (string, '(F6.2)') log10 (value)
         end if
      else
         string = '******'
      end if

      return
      end subroutine twlogr
      subroutine twnorm (n, value, x)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWNORM
!
!     COMPUTE THE MAX-NORM OF A VECTOR.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

!***  PRECISION > DOUBLE
      double precision    value, x
      integer &
         j, n
      intrinsic &
         abs, max

      dimension x(n)

      value = 0.0
      do 0100 j = 1, n
         value = max (value, abs (x(j)))
0100  continue

      return
      end subroutine twnorm
      subroutine twopnt &
        (error, text, versio, &
         above, active, below, buffer, comps, condit, groupa, groupb, &
         isize, iwork, mark, name, names, pmax, points, report, rsize, &
         rwork, signal, stride, time, u, x)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWOPNT
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character &
         column*80, ctemp1*80, ctemp2*80, header*80, id*9, name*(*), &
         report*(*), signal*(*), string*80, versio*(*), vnmbr*8
!**** PRECISION > DOUBLE
      double precision    above, below, buffer, condit, detail, maxcon, ratio, rvalue, &
         rwork, ssabs, ssrel, strid0, stride, tdabs, tdec, tdrel, temp, &
         timer, tinc, tmax, tmin, toler0, toler1, toler2, total, u, x, &
         ynorm
      external &
         evolve, refine, search, twcopy, twgrab, twlaps, twlast, twlogr, &
         twnorm, twsqez, twtime, twinit
      integer &
         age, cntrls, comps, count, desire, event, gmax, grid, groupa, &
         groupb, ilast, isize, ivalue, iwork, j, jacobs, k, label, len1, &
         len2, length, leveld, levelm, lines, names, nsteps, padd, pmax, &
         points, psave, qabove, qbelow, qbnds, qdvrg, qentry, qexit, &
         qfunct, qgrid, qjacob, qnull, qother, qrat1, qrat2, qrefin, &
         qs0, qs1, qsearc, qsolve, qtask, qtimst, qtotal, qtype, qusave, &
         qv1, qvary, qvary1, qvary2, qvsave, qxsave, qy0, qy1, return, &
         rlast, route, rsize, size, ssage, step, steps, steps0, steps1, &
         steps2, tdage, text, vnmbrs, xrepor
      intrinsic :: max
      logical &
         active, adapt, allow, error, exist, first, flag, found, lvalue, &
         mark, mess, satisf, steady, time

      parameter (id = 'TWOPNT:  ')
      parameter (cntrls = 22)
      parameter (gmax = 100)
      parameter (lines = 20)
      parameter (vnmbrs = 12)

!     REPORT CODES
      parameter (qnull = 0, qbnds = 1, qdvrg = 2)

!     LOCATION OF DATA IN ARRAYS DETAIL, EVENT, TIMER, AND TOTAL.  THE
!     LOCATIONS ARE CHOSEN TO SIMPLIFY WRITE STATEMENTS.  DETAIL USES
!     ONLY 1 : 8, EVENT USES ONLY 5 : 8, TIMER USES 1 : 9, AND TOTAL
!     USES ONLY 2 : 9.  IN ADDITION, 2, 3, 4, 10, AND 11 ARE USED AS
!     MNEMONIC VALUES FOR QTASK.
      parameter &
        (qgrid  =  1, &
         qtimst =  2, &
         qsearc =  3, &
         qrefin =  4, &
         qfunct =  5, &
         qjacob =  6, &
         qsolve =  7, &
         qother =  8, &
         qtotal =  9, &
         qentry = 10, &
         qexit  = 11)

      dimension &
         above(groupa + comps + groupb), active(*), below(groupa + comps &
         + groupb), buffer(groupa + comps * pmax + groupb), column(3), &
         detail(gmax, qtotal), event(gmax, qtotal), header(6), &
         ivalue(cntrls), iwork(isize), lvalue(cntrls), mark(*), &
         name(names), ratio(2), rvalue(cntrls), rwork(rsize), &
         size(gmax), timer(qtotal), total(qtotal), u(groupa + comps * &
         pmax + groupb), vnmbr(vnmbrs), x(*)

      common / twcomi / ivalue
      common / twcoml / lvalue
      common / twcomr / rvalue

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
         string = vnmbr(vnmbrs)
         call twlast (length, string)
         write (text, 10001) id, 'double precision', string (1 : length)
         write (text, 10022) id
         write (text, 10021) id
         write (text, 10011) id, '???', ratio, toler1, toler2
         write (text, 10013) id, '???', ratio, toler1, toler2
         write (text, 10012) id
         write (text, 10002) id, 'FINAL SOLUTION:'
         write (text, 10002) id, 'INITIAL GUESS:'
         write (text, 10019) id
         write (text, 10018) id
         write (text, 10016) id
         write (text, 10015) id
         string = vnmbr(vnmbrs)
         call twlast (length, string)
         write (text, 10001) id, 'SINGLE PRECISION', string (1 : length)
         write (text, 10002) id, 'SOLVE THE PROBLEM.'
         write (text, 10010) id

         go to 9001
      end if

!///  CHECK THE VERSION.

      data vnmbr &
         / '3.18', '3.19', '3.20', '3.21', '3.22', '3.23', '3.24', &
           '3.25', '3.26', '3.27', '3.28', '3.29' /

      flag = .false.
      do 1010 j = 1, vnmbrs
         flag = flag .or. versio == 'double precision VERSION ' // vnmbr(j)
1010  continue
      error = .not. flag
      if (error) go to 9002

!///  SET THE CONTROLS.

!     SUBROUTINE TWINIT (ERROR, TEXT, FORCE)

      call twinit (error, text, .false.)
      if (error) go to 9003

      count = 0

      count = count + 1
      adapt = lvalue(count)

      count = count + 1
      leveld = ivalue(count)

      count = count + 1
      levelm = ivalue(count)

      count = count + 1
      if (lvalue(count)) then
         padd = ivalue(count)
      else
         padd = pmax
      end if

      count = count + 1
      ssabs = rvalue(count)

      count = count + 1
      ssage = ivalue(count)

      count = count + 1
      ssrel = rvalue(count)

      count = count + 1
      steady = lvalue(count)

      count = count + 1
      steps0 = ivalue(count)

      count = count + 1
      steps1 = ivalue(count)

      count = count + 1
      steps2 = ivalue(count)

      count = count + 1
      strid0 = rvalue(count)

      count = count + 1
      tdabs = rvalue(count)

      count = count + 1
      tdage = ivalue(count)

      count = count + 1
      tdec = rvalue(count)

      count = count + 1
      tdrel = rvalue(count)

      count = count + 1
      tinc = rvalue(count)

      count = count + 1
      tmax = rvalue(count)

      count = count + 1
      tmin = rvalue(count)

      count = count + 1
      toler0 = rvalue(count)

      count = count + 1
      toler1 = rvalue(count)

      count = count + 1
      toler2 = rvalue(count)

      error = .not. (count == cntrls)
      if (error) go to 9004

!///  PRINT THE ENTRY BANNER AT ALL PRINT LEVELS.

      string = vnmbr(vnmbrs)
      call twlast (length, string)
      if ((0 < levelm .or. mess) .and. 0 < text) &
         write (text, 10001) id, 'double precision', string (1 : length)

!///  CHECK THE ARGUMENTS.

      error = .not. (leveld <= levelm)
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
      if (adapt .and. 0 < points) &
         call twcopy (points, x, rwork(qxsave))
      call twcopy (groupa + comps * points + groupb, u, rwork(qusave))

!     GO TO 1090 WHEN RETURN = 1
      return = 1
      go to 9911
1090  continue

!///  PRINT LEVELS 11, 21, AND 22.

      if (0 < leveld .and. 0 < text) then
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

      if (levelm == 1 .and. 0 < text) then
         if (0 < leveld) write (text, 10002) id, &
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
         if (0 < steps0) then
            qtask = qtimst
            desire = steps0
         else if (steady) then
            qtask = qsearc
         else
            error = .true.
            go to 9014
         end if

!///  SEARCH WAS THE PREVIOUS TASK.

      else if (qtask == qsearc) then
         if (found) then
            if (adapt) then
               qtask = qrefin
            else
               qtask = qexit
               report = ' '
            end if
         else
            if (allow .and. 0 < steps1) then
               qtask = qtimst
               desire = steps1
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
            if (steady) then
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
         if (adapt .and. 0 < points) &
            call twcopy (points, rwork(qxsave), x)
         call twcopy &
            (groupa + comps * points + groupb, rwork(qusave), u)
      end if

!///  PRINT LEVEL 11 OR 21.

!     SAVE THE STATUS REPORTS DURING REVERSE COMMUNICATION
      string = report

      if (leveld == 1 .and. 0 < text) then
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

      if (0 < levelm .and. 0 < text) then
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
      if (adapt) then

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

      if (adapt) write (text, 10009) header, &
         (size(j), (detail(j, k) / event(j, k), k = 5, 7), &
         (event(j, k), k = 5, 7), j = 1, grid)

      end if

!///  REPORT THE COMPLETION STATUS.

      if (0 < levelm) then
         if (report == ' ') then
            write (text, 10010) id
         else if (report == 'NO SPACE') then
            write (string, '(I10)') points
            call twsqez (length, string)
            write (text, 10011) &
               id, string (1 : length), ratio, toler1, toler2
         else if (report == 'NOT SOLVED') then
            write (text, 10012) id
         else if (report == 'SOME SOLVED') then
            write (string, '(I10)') points
            call twsqez (length, string)
            write (text, 10013) &
               id, string (1 : length), ratio, toler1, toler2
         else
            error = .true.
            go to 9016
         end if
      end if

!///  BOTTOM OF THE REPORT BLOCK.

      end if

!///  BOTTOM OF THE EXIT BLOCK.

      go to 99999

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

      if (1 < levelm) then
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
!        GROUPB, LEVELD, LEVELM, NAME, NAMES, POINTS, REPORT, S0, S1, &
!        SIGNAL, STEPS, SUCCES, V0, V1, XXABS, XXAGE, XXREL, Y0, Y0NORM, &
!        Y1)

      call search &
        (error, text, &
         rwork(qabove), age, rwork(qbelow), buffer, comps, condit, &
         exist, groupa, groupb, leveld - 1, levelm - 1, name, names, &
         points, xrepor, rwork(qs0), rwork(qs1), signal, nsteps, found, &
         u, rwork(qv1), ssabs, ssage, ssrel, rwork(qy0), ynorm, &
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
         if (adapt .and. 0 < points) &
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

      if (levelm == 1 .and. 0 < text) then
!        GO TO 4040 WHEN LABEL = 2
         label = 2
         go to 7010
      end if
4040  continue

!///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE SEARCH BLOCK.

      if (1 < levelm) then
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

      if (1 < levelm) then
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
!        ACTIVE, BUFFER, COMPS, LEVELD, LEVELM, MARK, NEWX, PADD, PMAX, &
!        POINTS, RATIO, RATIO1, RATIO2, SIGNAL, SUCCES, TOLER0, TOLER1, &
!        TOLER2, U, VARY1, VARY2, WEIGHT, X)

      call refine &
        (error, text, &
         active, &
         buffer(groupa + 1), comps, leveld - 1, levelm - 1, mark, &
         found, padd, pmax, points, ratio, rwork(qrat1), rwork(qrat2), &
         signal, satisf, toler0, toler1, toler2, u(groupa + 1), &
         iwork(qvary1), iwork(qvary2), iwork(qvary), x)
      if (error) go to 9018

!///  SERVICE REQUESTS FROM REFINE: PASS REQUESTS TO THE CALLER.

      if (signal /= ' ') then
!        INSERT THE GROUP A AND B UNKNOWNS
         do 5040 j = 1, groupa
            buffer(j) = u(j)
5040     continue
         do 5050 j = 1, groupb
            buffer(groupa + comps * points + j) = rwork(qvsave - 1 + j)
5050     continue

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

      if (levelm == 1 .and. 0 < text) then
         write (text, '()')
!        GO TO 5120 WHEN LABEL = 3
         label = 3
         go to 7010
      end if
5120  continue

!///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE REFINE BLOCK.

      if (1 < levelm) then
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

      if (1 < levelm) then
         if (0 < text) write (text, 10020) id
      end if

!///  CALL EVOLVE.

6020  continue

!     SUBROUTINE EVOLVE &
!       (ERROR, TEXT, &
!        ABOVE, BELOW, BUFFER, COMPS, CONDIT, DESIRE, GROUPA, GROUPB, &
!        LEVELD, LEVELM, NAME, NAMES, POINTS, REPORT, S0, S1, SIGNAL, &
!        STEP, STEPS2, STRID0, STRIDE, SUCCES, TDABS, TDAGE, TDEC, &
!        TDREL, TIME, TINC, TMAX, TMIN, V0, V1, VSAVE, Y0, Y1, YNORM)

      call evolve &
        (error, text, &
         rwork(qabove), rwork(qbelow), buffer, comps, condit, desire, &
         groupa, groupb, leveld - 1, levelm - 1, name, names, points, &
         xrepor, rwork(qs0), rwork(qs1), signal, step, steps2, strid0, &
         stride, found, tdabs, tdage, tdec, tdrel, time, tinc, tmax, &
         tmin, u, rwork(qv1), rwork(qvsave), rwork(qy0), rwork(qy1), &
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

      if (levelm == 1 .and. 0 < text) then
!        GO TO 6040 WHEN LABEL = 4
         label = 4
         go to 7010
      end if
6040  continue

!///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE EVOLVE BLOCK.

      if (1 < levelm) then
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
      else if (qtask == qentry .and. adapt) then
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
      go to 99999
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
      go to 99999
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
      go to 99999
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
      go to 99999
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

10001 format &
         (/1X, a9, a, ' (TWO POINT BOUNDARY VALUE PROBLEM) SOLVER,' &
         /10X, 'VERSION ', a, &
         ' OF APRIL 1998 BY DR. JOSEPH F. GRCAR.')

10002 format &
         (/1X, a9, a)

10003 format &
         (3(/10X, a35)/)

10004 format &
        (/1X, a9, a, ' TOTAL COMPUTER TIME (SEE BREAKDOWN BELOW).')

10005 format &
        (/10X, 'PERCENT OF TOTAL COMPUTER TIME FOR VARIOUS TASKS:' &
         /3(/10X, a38, a27) &
        //(10X, i6, 2X, f6.1, 1X, 3(1X, f6.1), 1X, 4(1X, f6.1)))

10006 format &
        (/12X, 'TASK TOTALS:', 1X, 3(1X, f6.1), 1X, 4(1X, f6.1))

10007 format &
        (/10X, 'SOME GRIDS ARE OMITTED, BUT THE TOTALS ARE FOR ALL.')

10008 format &
        (3(/24X, a51) &
        //10X, a12, f8.1, 5F9.1 &
         /10X, a12, f8.3, 2F9.3 &
         /10X, a12, i8, 2i9)

10009 format &
        (/10X, 'AVERAGE COMPUTER TIMES FOR, AND NUMBERS OF, SUBTASKS:' &
         /3(/10X, a37, a25) &
        //(10X, i6, 3X, f7.3, 2X, f7.3, 2X, f7.3, 1X, 3(2X, i7)))

10010 format &
        (/1X, a9, 'SUCCESS.  PROBLEM SOLVED.')

10011 format &
        (/1X, a9, 'FAILURE.  A SOLUTION WAS FOUND FOR A GRID WITH ', a &
        /10X, 'POINTS, BUT ONE OR BOTH RATIOS ARE TOO LARGE.' &
        //22X, '   RATIO 1     RATIO 2' &
        //10X, '     FOUND', 2F12.2 &
         /10X, '   DESIRED', 2F12.2 &
        //10X, 'A LARGER GRID COULD NOT BE FORMED.')

10012 format &
        (/1X, a9, 'FAILURE.  NO SOLUTION WAS FOUND.')

10013 format &
        (/1X, a9, 'FAILURE.  A SOLUTION WAS FOUND FOR A GRID WITH ', a &
        /10X, 'POINTS, BUT ONE OR BOTH RATIOS ARE TOO LARGE.' &
        //22X, '   RATIO 1     RATIO 2' &
        //10X, '     FOUND', 2F12.2 &
         /10X, '   DESIRED', 2F12.2 &
        //10X, 'A SOLUTION COULD NOT BE FOUND FOR A LARGER GRID.')

10014 format &
        (/1X, a9, 'CALLING SEARCH TO SOLVE THE STEADY STATE PROBLEM.')

10015 format &
        (/1X, a9, 'SEARCH FOUND THE STEADY STATE.')

10016 format &
        (/1X, a9, 'SEARCH DID NOT FIND THE STEADY STATE.')

10017 format &
        (/1X, a9, 'CALLING REFINE TO PRODUCE A NEW GRID.')

10018 format &
        (/1X, a9, 'REFINE SELECTED A NEW GRID.')

10019 format &
        (/1X, a9, 'REFINE DID NOT SELECT A NEW GRID.')

10020 format &
        (/1X, a9, 'CALLING EVOLVE TO PERFORM TIME EVOLUTION.')

10021 format &
        (/1X, a9, 'EVOLVE PERFORMED A TIME EVOLUTION.')

10022 format &
        (/1X, a9, 'EVOLVE DID NOT PERFORM A TIME EVOLUTION.')

10023 format &
         (10X, a8, 3X, a6, 2X, a6, 3X, a)

80001 format &
         ('(', a, ' ', i10, ')')

80002 format &
         (10X, 1p, e10.2, 2X, e10.2, 3X, a)

80003 format &
         (10X, '  ... MORE')

!///////////////////////////////////////////////////////////////////////
!
!     ERROR MESSAGES.
!
!///////////////////////////////////////////////////////////////////////

!     GO TO 99999

9001  if (0 < text) write (text, 99001) id, route
      if (.not. mess) go to 99999

9002  if (0 < text) then
         call twlast (length, versio)
         write (text, 99002) id, versio (1 : length), vnmbr(vnmbrs)
         do 9901 j = vnmbrs - 1, 1, - 1
            write (text, '(10X, A, A)') &
               ' CAN REPLACE:  double precision VERSION ', vnmbr(j)
9901     continue
      end if
      if (.not. mess) go to 99999

9003  if (0 < text) write (text, 99003) id
      if (.not. mess) go to 99999

9004  if (0 < text) write (text, 99004) id, cntrls, count
      if (.not. mess) go to 99999

9005  if (0 < text) write (text, 99005) id, leveld, levelm
      if (.not. mess) go to 99999

9006  if (0 < text) write (text, 99006) id, &
         comps, points, groupa, groupb
      if (.not. mess) go to 99999

9007  if (0 < text) write (text, 99007) id, comps, points
      if (.not. mess) go to 99999

9008  if (0 < text) write (text, 99008) id, &
         comps, points, groupa, groupb, groupa + comps * points + groupb
      if (.not. mess) go to 99999

9009  if (0 < text) write (text, 99009) id, &
         names, comps, groupa, groupb, groupa + comps + groupb
      if (.not. mess) go to 99999

9010  if (0 < text) write (text, 99010) id, points, pmax
      if (.not. mess) go to 99999

9011  if (0 < text) then
         write (text, 99011) id, &
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
      if (.not. mess) go to 99999

9012  if (0 < text) write (text, 99012) id
      if (.not. mess) go to 99999

9013  if (0 < text) write (text, 99013) id, &
         isize, rsize, ilast, rlast
      if (.not. mess) go to 99999

9014  if (0 < text) write (text, 99014) id
      if (.not. mess) go to 99999

9015  if (0 < text) write (text, 99015) id
      if (.not. mess) go to 99999

9016  if (0 < text) write (text, 99016) id
      if (.not. mess) go to 99999

9017  if (0 < text) write (text, 99017) id
      if (.not. mess) go to 99999

9018  if (0 < text) write (text, 99018) id
      if (.not. mess) go to 99999

9019  if (0 < text) write (text, 99019) id
      if (.not. mess) go to 99999

9020  if (0 < text) write (text, 99020) id, label
      if (.not. mess) go to 99999

9021  if (0 < text) write (text, 99021) id, return
      if (.not. mess) go to 99999

99001 format &
        (/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
       //10X, i10, '  ROUTE')

99002 format &
        (/1X, a9, 'ERROR.  THE CALLING PROGRAM EXPECTS A VERSION OF' &
        /10X, 'TWOPNT NOT COMPATIBLE WITH THIS VERSION.' &
       //10X, '     EXPECTS:  ', a &
       //10X, 'THIS VERSION:  double precision VERSION ', a)

99003 format &
        (/1X, a9, 'ERROR.  TWINIT FAILS.')

99004 format &
        (/1X, a9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.' &
       //10X, i10, '  CONTROLS' &
        /10X, i10, '  COUNTED')

99005 format &
        (/1X, a9, 'ERROR.  THE PRINTING LEVELS ARE OUT OF ORDER.' &
        /10X, 'LEVELD CANNOT EXCEED LEVELM.' &
       //10X, i10, '  LEVELD, FOR SOLUTIONS' &
        /10X, i10, '  LEVELM, FOR MESSAGES')

99006 format &
        (/1X, a9, 'ERROR.  NUMBERS OF ALL TYPES OF UNKNOWNS MUST BE AT' &
        /10X, 'LEAST ZERO.' &
       //10X, i10, '  COMPS, COMPONENTS' &
        /10X, i10, '  POINTS' &
        /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
        /10X, i10, '  GROUPB, GROUP B UNKNOWNS')

99007 format &
        (/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE' &
        /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE.' &
       //10X, i10, '  COMPS, COMPONENTS' &
        /10X, i10, '  POINTS')

99008 format &
        (/1X, a9, 'ERROR.  TOTAL UNKNOWNS MUST BE POSITIVE.' &
       //10X, i10, '  COMPS, COMPONENTS' &
        /10X, i10, '  POINTS' &
        /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
        /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
        /10X, i10, '  TOTAL NUMBER')

99009 format &
        (/1X, a9, 'ERROR.  THE NUMBER OF NAMES IS WRONG.' &
       //10X, i10, '  NAMES' &
       //10X, i10, '  COMPS, COMPONENTS' &
        /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
        /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
        /10X, i10, '  TOTAL NUMBER')

99010 format &
        (/1X, a9, 'ERROR.  THERE ARE TOO MANY POINTS.' &
       //10X, i10, '  POINTS' &
        /10X, i10, '  PMAX, LIMIT ON POINTS')

99011 format &
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

99012 format &
        (/1X, a9, 'ERROR.  TWGRAB FAILS.')

99013 format &
        (/1X, a9, 'ERROR.  ONE OR BOTH WORK SPACES ARE TOO SMALL.' &
       //25X, '   INTEGER        REAL' &
       //10X, ' PRESENT SIZE', 2i12 &
        /10X, 'REQUIRED SIZE', 2i12)

99014 format &
        (/1X, a9, 'ERROR.  NEITHER THE INITIAL TIME EVOLUTION NOR THE' &
        /10X, 'SEARCH FOR THE STEADY STATE IS ALLOWED.')

99015 format &
        (/1X, a9, 'ERROR.  UNKNOWN TASK.')

99016 format &
        (/1X, a9, 'ERROR.  UNKNOWN REPORT CODE.')

99017 format &
        (/1X, a9, 'ERROR.  SEARCH FAILS.')

99018 format &
        (/1X, a9, 'ERROR.  REFINE FAILS.')

99019 format &
        (/1X, a9, 'ERROR.  EVOLVE FAILS.')

99020 format &
        (/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
       //10X, i10, '  LABEL')

99021 format &
        (/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
       //10X, i10, '  RETURN')

!///  EXIT.

      stop
99999 continue
      return
      end
      subroutine twprep &
        (error, text, &
         a, asize, buffer, comps, condit, groupa, groupb, pivot, points, &
         return)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWPREP
!
!     EVALUATE A BLOCK TRIDIAGONAL JACOBIAN MATRIX BY ONE-SIDED FINITE
!     DIFFERENCES AND REVERSE COMMUNICATION, PACK THE MATRIX INTO THE
!     LINPACK BANDED FORM, SCALE THE ROWS, AND FACTOR THE MATRIX USING
!     LINPACK'S SGBCO.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

      character id*9, string*80
      double precision    a, absol, buffer, condit, delta, eps, relat, sum, temp
      external &
         tweps, twgbco, twsqez
      integer &
         asize, block, blocks, cfirst, clast, col, comps, count, diag, &
         groupa, groupb, j, lda, length, lines, n, offset, pivot, &
         points, rfirst, rlast, route, row, skip, text, width
      intrinsic &
         abs, int, max, min, mod, sqrt
      logical &
         error, found, mess, return

      parameter (id = 'TWPREP:  ')
      parameter (lines = 20)

      dimension &
         a(asize), pivot(groupa + comps * points + groupb), &
         buffer(groupa + comps * points + groupb)

      save

!///////////////////////////////////////////////////////////////////////
!
!     (1) PROLOGUE.
!
!///////////////////////////////////////////////////////////////////////

!///  EVERY-TIME INITIALIZATION.

!     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

!///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      if (return) then
         return = .false.
         go to (2030, 3050) route
         error = .true.
         go to 9001
      endif

!///  CHECK THE ARGUMENTS.

      n = groupa + comps * points + groupb
      error = .not. (((0 < comps) .eqv. (0 < points)) .and. &
         0 <= comps .and. 0 <= points .and. 0 <= groupa .and. &
         0 <= groupb .and. 0 < n)
      if (error) go to 9002

      width = comps + max (comps, groupa, groupb) - 1
      error = .not. ((3 * width + 2) * n <= asize)
      if (error) go to 9003

!///  WRITE ALL MESSAGES.

      if (mess .and. 0 < text) then
         route = 0
         go to 9001
      end if

!///  FORM MACHINE EPSILON AND THE ABSOLUTE AND RELATIVE PERTURBATIONS.

      call tweps (eps)
      absol = sqrt (eps)
      relat = sqrt (eps)

!///  INITIALIZE COUNTERS AND POINTERS.

!     MAIN DIAGONAL ROW IN THE PACKING WHICH PLACES DIAGONALS IN ROWS

      diag = 2 * width + 1

!     PACKED ROW DIMENSION

      lda = 3 * width + 1
      skip = 2 * width + 1

!     BLOCKS AND BLOCK SIZES
!     ARRAY PIVOT HOLDS BLOCK SIZES AND POINTERS TEMPORARILY

      blocks = 0
      if (0 < groupa) then
         blocks = blocks + 1
         pivot(blocks) = groupa
      end if

      do 1020 j = 1, points
         blocks = blocks + 1
         pivot(blocks) = comps
1020  continue

      if (0 < groupb) then
         blocks = blocks + 1
         pivot(blocks) = groupb
      end if

!///////////////////////////////////////////////////////////////////////
!
!     (2) INITIALIZE THE COLUMNS OF THE MATRIX.
!
!///////////////////////////////////////////////////////////////////////

!///  STORE THE EVALUATION VECTOR.

      do 2010 j = 1, n
         a(j) = buffer(j)
2010  continue

!///  CLEAR THE MATRIX.

       do 2020 j = n + 1, (3 * width + 2) * n
          a(j) = 0.0
2020  continue

!///  EVALUATE THE FUNCTION AT THE UNPERTURBED X.

!     GO TO 2030 WHEN ROUTE = 1
      route = 1
      return = .true.
      go to 99999
2030  continue

!///  PLACE THE FUNCTION VALUES IN THE MATRIX.

      clast = 0
      do 2060 block = 1, blocks
         cfirst = clast + 1
         clast = clast + pivot(block)

         if (1 < block) then
            rfirst = cfirst - pivot(block - 1)
         else
            rfirst = cfirst
         end if

         if (block < blocks) then
            rlast = clast + pivot(block + 1)
         else
            rlast = clast
         end if

         do 2050 col = cfirst, clast
            offset = n + diag - col + lda * (col - 1)
            do 2040 row = rfirst, rlast
               a(offset + row) = buffer(row)
2040        continue
2050     continue
2060  continue

!///////////////////////////////////////////////////////////////////////
!
!     (3) FORM THE COLUMNS OF THE MATRIX.
!
!///////////////////////////////////////////////////////////////////////

!///  TOP OF THE LOOP OVER GROUPS OF COLUMNS.

3010  continue
      found = .false.

!///  RESTORE THE EVALUATION VECTOR.

      do 3020 j = 1, n
         buffer(j) = a(j)
3020  continue

!///  PERTURB THE VECTOR AT INDEPENDENT POSITIONS.

      block = 1
      cfirst = 1
3030  continue
         if (0 < pivot(block)) then
            found = .true.
            col = cfirst - 1 + pivot(block)
            if (0 <= a(col)) then
               delta = relat * a(col) + absol
            else
               delta = relat * a(col) - absol
            end if
            buffer(col) = buffer(col) + delta

            count = 3
         else
            count = 1
         end if

         do 3040 j = 1, count
            if (block == 1 .and. 0 < groupa) then
               cfirst = cfirst + groupa
            else if (block == blocks .and. 0 < groupb) then
               cfirst = cfirst + groupb
            else
               cfirst = cfirst + comps
            end if
            block = block + 1
3040     continue

      if (block <= blocks) go to 3030

!///  EXIT OF THE LOOP OVER GROUPS OF COLUMNS.

      if (.not. found) go to 3090

!///  EVALUATE THE FUNCTION AT THE PERTURBED VALUES.

!     GO TO 3050 WHEN ROUTE = 2
      route = 2
      return = .true.
      go to 99999
3050  continue

!///  DIFFERENCE TO FORM THE COLUMNS OF THE JACOBIAN MATRIX.

      block = 1
      cfirst = 1
3060  continue
         if (0 < pivot(block)) then
            col = cfirst - 1 + pivot(block)
            pivot(block) = pivot(block) - 1

            if (0 <= a(col)) then
               delta = relat * a(col) + absol
            else
               delta = relat * a(col) - absol
            end if
            temp = 1.0 / delta
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

            do 3070 row = rfirst, rlast
               a(offset + row) = (buffer(row) - a(offset + row)) * temp
3070        continue

            count = 3
         else
            count = 1
         end if

         do 3080 j = 1, count
            if (block == 1 .and. 0 < groupa) then
               cfirst = cfirst + groupa
            else if (block == blocks .and. 0 < groupb) then
               cfirst = cfirst + groupb
            else
               cfirst = cfirst + comps
            end if
            block = block + 1
3080     continue

      if (block <= blocks) go to 3060

!///  BOTTOM OF THE LOOP OVER GROUPS OF COLUMNS.

      go to 3010
3090  continue

!///////////////////////////////////////////////////////////////////////
!
!     (4) CHECK FOR ZERO COLUMNS.
!
!///////////////////////////////////////////////////////////////////////

      count = 0
      do 4020 col = 1, n
         offset = n + diag - col + lda * (col - 1)
         sum = 0.0
         do 4010 row = max (col - width, 1), min (col + width, n)
            sum = sum + abs (a(offset + row))
4010     continue
         a(col) = sum

         if (sum == 0.0) count = count + 1
4020  continue

      error = .not. (count == 0)
      if (error) go to 9004

!///////////////////////////////////////////////////////////////////////
!
!     (5) SCALE THE ROWS.
!
!///////////////////////////////////////////////////////////////////////

      count = 0
      do 5030 row = 1, n
         offset = n + diag + row
         sum = 0.0
         do 5010 col = max (row - width, 1), min (row + width, n)
            sum = sum + abs (a(offset - col + lda * (col - 1)))
5010     continue

         if (sum == 0.0) then
            count = count + 1
            a(row) = sum
         else
            temp = 1.0 / sum
            a(row) = temp

            do 5020 col = max (row - width, 1), min (row + width, n)
               a(offset - col + lda * (col - 1)) &
                  = a(offset - col + lda * (col - 1)) * temp
5020        continue
         endif
5030  continue

      error = .not. (count == 0)
      if (error) go to 9005

!///////////////////////////////////////////////////////////////////////
!
!     (6) FACTOR THE MATRIX.
!
!///////////////////////////////////////////////////////////////////////

      call twgbco &
        (a(n + 1), lda, n, width, width, pivot, condit, buffer)

      error = condit == 0.0
      if (error) go to 9006

      condit = 1.0 / condit

!///////////////////////////////////////////////////////////////////////
!
!     INFORMATIVE MESSAGES.
!
!///////////////////////////////////////////////////////////////////////

80001 format &
        (10X, a)

80002 format &
        (10X, '... MORE')

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
        (/1X, a9, 'ERROR.  THE MATRIX SPACE IS TOO SMALL.' &
       //10X, i10, '  COMPS, COMPONENTS' &
        /10X, i10, '  POINTS' &
        /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
        /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
        /10X, i10, '  MATRIX ORDER' &
        /10X, i10, '  STRICT HALF BANDWIDTH' &
       //10X, i10, '  SPACE REQUIRED' &
        /10X, i10, '  ASIZE, PROVIDED')

99004 format &
        (/1X, a9, 'ERROR.  SOME COLUMNS ARE ZERO.' &
       //10X, i10, '  COMPS, COMPONENTS' &
        /10X, i10, '  POINTS' &
        /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
        /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
        /10X, i10, '  TOTAL COLUMNS' &
        /10X, i10, '  ZERO COLUMNS' &
       //10X, 'UNKNOWNS WITH ZERO COLUMNS:' &
        /)

99005 format &
        (/1X, a9, 'ERROR.  SOME ROWS ARE ZERO.' &
       //10X, i10, '  COMPS, COMPONENTS' &
        /10X, i10, '  POINTS' &
        /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
        /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
        /10X, i10, '  TOTAL ROWS' &
        /10X, i10, '  ZERO ROWS' &
       //10X, 'ZERO ROWS:' &
        /)

99006 format &
        (/1X, a9, 'ERROR.  THE JACOBIAN MATRIX IS SINGULAR.')

!///  EXIT.

      stop
99999 continue
      return
      end
      subroutine twseti (error, text, contrl, value)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWSETI
!
!     SET A CONTROL THAT TAKES AN INTEGER VALUE.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character &
         contrl*(*), id*9, string*80
      external &
         twinit, twlast
      integer &
         cntrls, count, ivalue, length, text, value
      logical &
         error, found, mess, lvalue

      parameter (id = 'TWSETI:  ')
      parameter (cntrls = 22)

      dimension ivalue(cntrls), lvalue(cntrls)

      common / twcomi / ivalue
      common / twcoml / lvalue

!///  WRITE ALL MESSAGES.

!     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

      if (mess .and. 0 < text) go to 9001

!///  INITIALIZE THE CONTROLS.

      call twinit (error, text, .false.)
      if (error) go to 9001

!///  SET THE CONTROLS.

      count = 0
      found = .false.

!     ADAPT

      count = count + 1
      if (contrl == 'ADAPT') then
         error = .true.
         go to 9002
      end if

!     LEVELD

      count = count + 1
      if (contrl == 'LEVELD') then
         found = .true.
         ivalue(count) = value
      end if

!     LEVELM

      count = count + 1
      if (contrl == 'LEVELM') then
         found = .true.
         ivalue(count) = value
      end if

!     PADD

      count = count + 1
      if (contrl == 'PADD') then
         found = .true.
         ivalue(count) = value
         lvalue(count) = .true.
      end if

!     SSABS

      count = count + 1
      if (contrl == 'SSABS') then
         error = .true.
         go to 9003
      end if

!     SSAGE

      count = count + 1
      if (contrl == 'SSAGE') then
         found = .true.
         ivalue(count) = value
      end if

!     SSREL

      count = count + 1
      if (contrl == 'SSREL') then
         error = .true.
         go to 9003
      end if

!     STEADY

      count = count + 1
      if (contrl == 'STEADY') then
         error = .true.
         go to 9002
      end if

!     STEPS0

      count = count + 1
      if (contrl == 'STEPS0') then
         found = .true.
         ivalue(count) = value
      end if

!     STEPS1

      count = count + 1
      if (contrl == 'STEPS1') then
         found = .true.
         ivalue(count) = value
      end if

!     STEPS2

      count = count + 1
      if (contrl == 'STEPS2') then
         found = .true.
         ivalue(count) = value
      end if

!     STRID0

      count = count + 1
      if (contrl == 'STRID0') then
         error = .true.
         go to 9003
      end if

!     TDABS

      count = count + 1
      if (contrl == 'TDABS') then
         error = .true.
         go to 9003
      end if

!     TDAGE

      count = count + 1
      if (contrl == 'TDAGE') then
         found = .true.
         ivalue(count) = value
      end if

!     TDEC

      count = count + 1
      if (contrl == 'TDEC') then
         error = .true.
         go to 9003
      end if

!     TDREL

      count = count + 1
      if (contrl == 'TDREL') then
         error = .true.
         go to 9003
      end if

!     TINC

      count = count + 1
      if (contrl == 'TINC') then
         error = .true.
         go to 9003
      end if

!     TMAX

      count = count + 1
      if (contrl == 'TMAX') then
         error = .true.
         go to 9003
      end if

!     TMIN

      count = count + 1
      if (contrl == 'TMIN') then
         error = .true.
         go to 9003
      end if

!     TOLER0

      count = count + 1
      if (contrl == 'TOLER0') then
         error = .true.
         go to 9003
      end if

!     TOLER1

      count = count + 1
      if (contrl == 'TOLER1') then
         error = .true.
         go to 9003
      end if

!     TOLER2

      count = count + 1
      if (contrl == 'TOLER2') then
         error = .true.
         go to 9003
      end if

      error = .not. (count == cntrls)
      if (error) go to 9004

      error = .not. found
      if (error) go to 9005

!///  ERROR MESSAGES.

      go to 99999

9001  if (0 < text) write (text, 99001) id
      if (.not. mess) go to 99999

9002  if (0 < text) then
         call twlast (length, contrl)
         if (length <= 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99002) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9003  if (0 < text) then
         call twlast (length, contrl)
         if (length <= 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99003) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9004  if (0 < text) write (text, 99004) id, cntrls, count
      if (.not. mess) go to 99999

9005  if (0 < text) then
         call twlast (length, contrl)
         if (length <= 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99005) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

99001 format &
        (/1X, a9, 'ERROR.  TWINIT FAILS.')

99002 format &
        (/1X, a9, 'ERROR.  THE CONTROL TAKES A LOGICAL VALUE WHICH' &
        /10X, 'MUST BE SET USING TWSETL.' &
       //10X, '     CONTROL:  ', a)

99003 format &
        (/1X, a9, 'ERROR.  THE CONTROL TAKES A REAL VALUE WHICH MUST BE' &
        /10X, 'SET USING TWSETR.' &
       //10X, '     CONTROL:  ', a)

99004 format &
        (/1X, a9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.' &
       //10X, i10, '  CONTROLS' &
        /10X, i10, '  COUNTED')

99005 format &
        (/1X, a9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.' &
       //10X, '     CONTROL:  ', a)

!///  EXIT.

      stop
99999 continue
      return
      end subroutine twseti
      subroutine twsetl (error, text, contrl, value)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWSETL
!
!     SET A CONTROL THAT TAKES A LOGICAL VALUE.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character &
         contrl*(*), id*9, string*80
      external &
         twinit, twlast
      integer &
         cntrls, count, length, text
      logical &
         error, found, lvalue, mess, value

      parameter (id = 'TWSETL:  ')
      parameter (cntrls = 22)

      dimension lvalue(cntrls)

      common / twcoml / lvalue

!///  WRITE ALL MESSAGES.

!     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

      if (mess .and. 0 < text) go to 9001

!///  INITIALIZE THE CONTROLS.

      call twinit (error, text, .false.)
      if (error) go to 9001

!///  SET THE CONTROLS.

      count = 0
      found = .false.

!     ADAPT

      count = count + 1
      if (contrl == 'ADAPT') then
         found = .true.
         lvalue(count)= value
      end if

!     LEVELD

      count = count + 1
      if (contrl == 'LEVELD') then
         error = .true.
         go to 9002
      end if

!     LEVELM

      count = count + 1
      if (contrl == 'LEVELM') then
         error = .true.
         go to 9002
      end if

!     PADD

      count = count + 1
      if (contrl == 'PADD') then
         error = .true.
         go to 9002
      end if

!     SSABS

      count = count + 1
      if (contrl == 'SSABS') then
         error = .true.
         go to 9003
      end if

!     SSAGE

      count = count + 1
      if (contrl == 'SSAGE') then
         error = .true.
         go to 9002
      end if

!     SSREL

      count = count + 1
      if (contrl == 'SSREL') then
         error = .true.
         go to 9003
      end if

!     STEADY

      count = count + 1
      if (contrl == 'STEADY') then
         found = .true.
         lvalue(count)= value
      end if

!     STEPS0

      count = count + 1
      if (contrl == 'STEPS0') then
         error = .true.
         go to 9002
      end if

!     STEPS1

      count = count + 1
      if (contrl == 'STEPS1') then
         error = .true.
         go to 9002
      end if

!     STEPS2

      count = count + 1
      if (contrl == 'STEPS2') then
         error = .true.
         go to 9002
      end if

!     STRID0

      count = count + 1
      if (contrl == 'STRID0') then
         error = .true.
         go to 9003
      end if

!     TDABS

      count = count + 1
      if (contrl == 'TDABS') then
         error = .true.
         go to 9003
      end if

!     TDAGE

      count = count + 1
      if (contrl == 'TDAGE') then
         error = .true.
         go to 9002
      end if

!     TDEC

      count = count + 1
      if (contrl == 'TDEC') then
         error = .true.
         go to 9003
      end if

!     TDREL

      count = count + 1
      if (contrl == 'TDREL') then
         error = .true.
         go to 9003
      end if

!     TINC

      count = count + 1
      if (contrl == 'TINC') then
         error = .true.
         go to 9003
      end if

!     TMAX

      count = count + 1
      if (contrl == 'TMAX') then
         error = .true.
         go to 9003
      end if

!     TMIN

      count = count + 1
      if (contrl == 'TMIN') then
         error = .true.
         go to 9003
      end if

!     TOLER0

      count = count + 1
      if (contrl == 'TOLER0') then
         error = .true.
         go to 9003
      end if

!     TOLER1

      count = count + 1
      if (contrl == 'TOLER1') then
         error = .true.
         go to 9003
      end if

!     TOLER2

      count = count + 1
      if (contrl == 'TOLER2') then
         error = .true.
         go to 9003
      end if

      error = .not. (count == cntrls)
      if (error) go to 9004

      error = .not. found
      if (error) go to 9005

!///  ERROR MESSAGES.

      go to 99999

9001  if (0 < text) write (text, 99001) id
      if (.not. mess) go to 99999

9002  if (0 < text) then
         call twlast (length, contrl)
         if (length <= 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99002) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9003  if (0 < text) then
         call twlast (length, contrl)
         if (length <= 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99003) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9004  if (0 < text) write (text, 99004) id, cntrls, count
      if (.not. mess) go to 99999

9005  if (0 < text) then
         call twlast (length, contrl)
         if (length <= 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99005) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

99001 format &
        (/1X, a9, 'ERROR.  TWINIT FAILS.')

99002 format &
        (/1X, a9, 'ERROR.  THE CONTROL TAKES AN INTEGER VALUE WHICH' &
        /10X, 'MUST BE SET USING TWSETI.' &
       //10X, '     CONTROL:  ', a)

99003 format &
        (/1X, a9, 'ERROR.  THE CONTROL TAKES A REAL VALUE WHICH MUST BE' &
        /10X, 'SET USING TWSETR.' &
       //10X, '     CONTROL:  ', a)

99004 format &
        (/1X, a9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.' &
       //10X, i10, '  CONTROLS' &
        /10X, i10, '  COUNTED')

99005 format &
        (/1X, a9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.' &
       //10X, '     CONTROL:  ', a)

!///  EXIT.

      stop
99999 continue
      return
      end subroutine twsetl
      subroutine twsetr (error, text, contrl, value)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWSETR
!
!     SET A CONTROL THAT TAKES A REAL VALUE.
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character &
         contrl*(*), id*9, string*80
      double precision    rvalue, value
      external &
         twinit, twlast
      integer &
         cntrls, count, length, text
      logical &
         error, found, mess

      parameter (id = 'TWSETR:  ')
      parameter (cntrls = 22)

      dimension rvalue(cntrls)

      common / twcomr / rvalue

!///  WRITE ALL MESSAGES.

!     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

      if (mess .and. 0 < text) go to 9001

!///  INITIALIZE THE CONTROLS.

      call twinit (error, text, .false.)
      if (error) go to 9001

!///  SET THE CONTROLS.

      count = 0
      found = .false.

!     ADAPT

      count = count + 1
      if (contrl == 'ADAPT') then
         error = .true.
         go to 9002
      end if

!     LEVELD

      count = count + 1
      if (contrl == 'LEVELD') then
         error = .true.
         go to 9003
      end if

!     LEVELM

      count = count + 1
      if (contrl == 'LEVELM') then
         error = .true.
         go to 9003
      end if

!     PADD

      count = count + 1
      if (contrl == 'PADD') then
         error = .true.
         go to 9003
      end if

!     SSABS

      count = count + 1
      if (contrl == 'SSABS') then
         found = .true.
         rvalue(count) = value
      end if

!     SSAGE

      count = count + 1
      if (contrl == 'SSAGE') then
         error = .true.
         go to 9003
      end if

!     SSREL

      count = count + 1
      if (contrl == 'SSREL') then
         found = .true.
         rvalue(count) = value
      end if

!     STEADY

      count = count + 1
      if (contrl == 'STEADY') then
         error = .true.
         go to 9002
      end if

!     STEPS0

      count = count + 1
      if (contrl == 'STEPS0') then
         error = .true.
         go to 9003
      end if

!     STEPS1

      count = count + 1
      if (contrl == 'STEPS1') then
         error = .true.
         go to 9003
      end if

!     STEPS2

      count = count + 1
      if (contrl == 'STEPS2') then
         error = .true.
         go to 9003
      end if

!     STRID0

      count = count + 1
      if (contrl == 'STRID0') then
         found = .true.
         rvalue(count) = value
      end if

!     TDABS

      count = count + 1
      if (contrl == 'TDABS') then
         found = .true.
         rvalue(count) = value
      end if

!     TDAGE

      count = count + 1
      if (contrl == 'TDAGE') then
         error = .true.
         go to 9003
      end if

!     TDEC

      count = count + 1
      if (contrl == 'TDEC') then
         found = .true.
         rvalue(count) = value
      end if

!     TDREL

      count = count + 1
      if (contrl == 'TDREL') then
         found = .true.
         rvalue(count) = value
      end if

!     TINC

      count = count + 1
      if (contrl == 'TINC') then
         found = .true.
         rvalue(count) = value
      end if

!     TMAX

      count = count + 1
      if (contrl == 'TMAX') then
         found = .true.
         rvalue(count) = value
      end if

!     TMIN

      count = count + 1
      if (contrl == 'TMIN') then
         found = .true.
         rvalue(count) = value
      end if

!     TOLER0

      count = count + 1
      if (contrl == 'TOLER0') then
         found = .true.
         rvalue(count) = value
      end if

!     TOLER1

      count = count + 1
      if (contrl == 'TOLER1') then
         found = .true.
         rvalue(count) = value
      end if

!     TOLER2

      count = count + 1
      if (contrl == 'TOLER2') then
         found = .true.
         rvalue(count) = value
      end if

      error = .not. (count == cntrls)
      if (error) go to 9004

      error = .not. found
      if (error) go to 9005

!///  ERROR MESSAGES.

      go to 99999

9001  if (0 < text) write (text, 99001) id
      if (.not. mess) go to 99999

9002  if (0 < text) then
         call twlast (length, contrl)
         if (length <= 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99002) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9003  if (0 < text) then
         call twlast (length, contrl)
         if (length <= 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99003) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9004  if (0 < text) write (text, 99004) id, cntrls, count
      if (.not. mess) go to 99999

9005  if (0 < text) then
         call twlast (length, contrl)
         if (length <= 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99005) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

99001 format &
        (/1X, a9, 'ERROR.  TWINIT FAILS.')

99002 format &
        (/1X, a9, 'ERROR.  THE CONTROL TAKES A LOGICAL VALUE WHICH' &
        /10X, 'MUST BE SET USING TWSETL.' &
       //10X, '     CONTROL:  ', a)

99003 format &
        (/1X, a9, 'ERROR.  THE CONTROL TAKES AN INTEGER VALUE WHICH' &
        /10X, 'MUST BE SET USING TWSETI.' &
       //10X, '     CONTROL:  ', a)

99004 format &
        (/1X, a9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.' &
       //10X, i10, '  CONTROLS' &
        /10X, i10, '  COUNTED')

99005 format &
        (/1X, a9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.' &
       //10X, '     CONTROL:  ', a)

!///  EXIT.

      stop
99999 continue
      return
      end subroutine twsetr
      subroutine twshow &
        (error, text, &
         buffer, comps, grid, groupa, groupb, points, x)

!///////////////////////////////////////////////////////////////////////
!
!     T W O P N T
!
!     TWSHOW
!
!///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

      character id*9, string*80, title*80
      double precision    buffer, x
      external  twsqez
      integer &
         cols, comp, comps, count, first, groupa, groupb, groups, j, &
         last, length, point, points, text
      intrinsic &
         min
      logical &
         error, grid, mess

      parameter (id = 'TWSHOW:  ')

      dimension &
         buffer(groupa + comps * points + groupb), title(6), x(*)

!///////////////////////////////////////////////////////////////////////
!
!     (1) PROLOGUE.
!
!///////////////////////////////////////////////////////////////////////

!///  WRITE ALL MESSAGES.

!     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

      if (mess .and. 0 < text) go to 9001

!///  CHECK THE ARGUMENTS.

      error = .not. (((0 < comps) .eqv. (0 < points)) .and. &
         0 <= comps .and. 0 <= points .and. 0 <= groupa .and. &
         0 <= groupb .and. 0 < groupa + comps * points + groupb)
      if (error) go to 9001

!///  COUNT THE GROUPS.

      groups = 0
      if (0 < groupa) groups = groups + 1
      if (0 < groupb) groups = groups + 1
      if (0 < comps .and. 0 < points) groups = groups + 1

!///  CHOOSE NUMBER OF DATA COLUMNS.

      if (grid) then
         cols = 5
      else
         cols = 6
      end if

!///////////////////////////////////////////////////////////////////////
!
!     (2) PRINT THE GROUPED DATA.
!
!///////////////////////////////////////////////////////////////////////

      if (0 < text) then

      if (0 < groupa) then
         if (1 < groups) write (text, 10001) 'GROUP A UNKNOWNS'
         write (text, 10002) (j, buffer(j), j = 1, groupa)
      end if

      if (0 < groupb) then
         if (1 < groups) write (text, 10001) 'GROUP B UNKNOWNS'
         write (text, 10002) &
            (j, buffer(groupa + comps * points + j), j = 1, groupb)
      end if

!///////////////////////////////////////////////////////////////////////
!
!     (2) PRINT THE COMPONENTS AT POINTS.
!
!///////////////////////////////////////////////////////////////////////

      if (0 < comps .and. 0 < points) then
         if (1 < groups) write (text, 10001) 'COMPONENTS AT POINTS'

         do 2030 first = 1, comps, cols
            count = 0
            last = min (first + cols - 1, comps)
            do 2010 comp = first, last
               count = count + 1
               title(count) = ' '
               write (string, '(A5, I5)') 'COMP ', comp
               call twsqez (length, string)
               title(count) (11 - length : 10) = string
2010        continue

            if (grid) then
               write (text, 10003) &
                  'GRID POINT', (title(j), j = 1, count)
            else
               write (text, 10003) (title(j), j = 1, count)
            end if

            if (count == cols) then
               if (grid) then
                  write (text, 10004) (point, x(point), &
                     (buffer(groupa + comp + comps * (point - 1)), &
                     comp = first, last), point = 1, points)
               else
                  write (text, 10005) (point, &
                     (buffer(groupa + comp + comps * (point - 1)), &
                     comp = first, last), point = 1, points)
               end if
            else
               do 2020 point = 1, points
                  if (grid) then
                     write (text, 10004) point, x(point), &
                        (buffer(groupa + comp + comps * (point - 1)), &
                        comp = first, last)
                  else
                     write (text, 10005) point, &
                        (buffer(groupa + comp + comps * (point - 1)), &
                        comp = first, last)
                  end if
2020           continue
            end if
2030     continue
      end if

      end if

!///////////////////////////////////////////////////////////////////////
!
!     INFORMATIVE MESSAGES.
!
!///////////////////////////////////////////////////////////////////////

10001 format &
        (/10X, a)

10002 format &
       (/(10X, 4(i3, '> ', 1pe10.3)))

10003 format &
        (/14X, 6(1X, a10))

10004 format &
        (10X, 0p, i3, '>', f11.6, 1p, 5E11.3)

10005 format &
        (10X, 0p, i3, '>', 1p, 6E11.3)

!///////////////////////////////////////////////////////////////////////
!
!     ERROR MESSAGES.
!
!///////////////////////////////////////////////////////////////////////

      go to 99999

9001  if (0 < text) write (text, 99001) id, &
         comps, points, groupa, groupb, groupa + comps * points + groupb
      if (.not. mess) go to 99999

99001 format &
        (/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE' &
        /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES' &
        /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS' &
        /10X, 'MUST BE POSITIVE.' &
       //10X, i10, '  COMPS, COMPONENTS' &
        /10X, i10, '  POINTS' &
        /10X, i10, '  GROUPA, GROUP A UNKNOWNS' &
        /10X, i10, '  GROUPB, GROUP B UNKNOWNS' &
        /10X, i10, '  TOTAL UNKNOWNS')
      if (.not. mess) go to 99999

!///  EXIT.

      stop
99999 continue
      return
      end




