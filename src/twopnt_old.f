c     cvs $revision: 1.1.1.1 $ reposited $date: 2006/05/26 19:09:34 $

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     VERSION 3.29 OF APRIL 1998
C
C     THE TWOPNT PROGRAM FOR BOUNDARY VALUE PROBLEMS
C
C     WRITTEN BY: DR. JOSEPH F. GRCAR
C                 SANDIA NATIONAL LABORATORIES
C                 MAIL STOP 9051
C                 LIVERMORE, CALIFORNIA  94551-0969  USA
C
C                 (925) 294-2662
C                 (FTS) 234-2662
C
C                 na.grcar@na-net.ornl.gov
C                 sepp@california.sandia.gov
C
C///////////////////////////////////////////////////////////////////////
C
C     DOCUMENTATION:
C
C     J. F. Grcar, "The Twopnt Program for Boundary Value Problems,"
C     Sandia National Laboratories Report SAND91-8230, Livermore,
C     California, April 1992.  Reprinted February 1996.
C
C///////////////////////////////////////////////////////////////////////
C
C     CHANGES FROM THE PREVIOUS VERSION:
C
C     1) PUT CHANGE BLOCK AROUND DECLARATION OF SAME IN TWEPS.
C
C     2) ALTER AREA CODE.
C
C///////////////////////////////////////////////////////////////////////

      subroutine evolve
     +  (error, text,
     +   above, below, buffer, comps, condit, desire, groupa, groupb,
     +   leveld, levelm, name, names, points, report, s0, s1, signal,
     +   step, steps2, strid0, stride, succes, tdabs, tdage, tdec,
     +   tdrel, time, tinc, tmax, tmin, v0, v1, vsave, y0, y1, ynorm)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     EVOLVE
C
C     PERFORM TIME EVOLUTION.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - p, r - z), integer (q)
      character
     +   cword*80, header*80, id*9, jword*80, name*(*), remark*80,
     +   signal*(*), yword*80
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   above, below, buffer, change, condit, csave, dummy, high, low,
     +   s0, s1, strid0, stride, tdabs, tdec, tdrel, tinc, tmax, tmin,
     +   v0, v1, vsave, y0, y1, ynorm
      external
     +   search, twcopy, twlogr, twnorm
      integer
     +   age, agej, comps, count, desire, first, groupa, groupb, j,
     +   last, length, leveld, levelm, names, number, points, qbnds,
     +   qdvrg, qnull, report, route, step, steps2, tdage, text,
     +   xrepor
      intrinsic
     +   log10, max, min
      logical
     +   error, exist, jacob, mess, succes, time, xsucce

      parameter (id = 'EVOLVE:  ')

C     REPORT CODES
      parameter (qnull = 0, qbnds = 1, qdvrg = 2)

      dimension
     +   above(groupa + comps * points + groupb),
     +   below(groupa + comps * points + groupb),
     +   buffer(groupa + comps * points + groupb), header(2, 3),
     +   name(names), s0(groupa + comps * points + groupb),
     +   v0(groupa + comps * points + groupb),
     +   v1(groupa + comps * points + groupb),
     +   vsave(groupa + comps * points + groupb),
     +   y0(groupa + comps * points + groupb)

C///  SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.

      save

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  INITIALIZE.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

C     TURN OF REVERSE COMMUNICATION FLAGS.
      time = .false.

C     TURN OFF ALL COMPLETION STATUS FLAGS.
      error = .false.
      report = qnull
      succes = .false.

C///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      if (signal .ne. ' ') then
         go to (1010, 1020, 1060, 1080, 1090, 2020) route
         error = .true.
         go to 9001
      end if

C///  CHECK THE ARGUMENTS.

      error = .not. (((0 .lt. comps) .eqv. (0 .lt. points)) .and.
     +   0 .le. comps .and. 0 .le. points .and. 0 .le. groupa .and.
     +   0 .le. groupb .and. 0 .lt. groupa + comps * points + groupb)
      if (error) go to 9002

      error = .not. (0 .lt. desire)
      if (error) go to 9003

      error = .not. (1.0 .le. tdec .and. 1.0 .le. tinc)
      if (error) go to 9004

      error = .not. (0.0 .lt. tmin .and. tmin .le. tmax)
      if (error) go to 9005

      error = .not. (tmin .le. strid0 .and. strid0 .le. tmax)
      if (error) go to 9006

      error = .not. (0 .le. step)
      if (error) go to 9007

      error = 1.0 .lt. tinc .and. .not. 0 .lt. steps2
      if (error) go to 9008

C///  WRITE ALL MESSAGES.

C                     123456789_123456789_123456789_123456789_1234
C                     123456   123456   123456   123456   12345
      header(1, 1) = '  TIME   LOG10                      NEWTON S'
      header(1, 2) = ' POINT   ------------------------   --------'
      header(1, 3) = 'NUMBER   NORM F   CHANGE   STRIDE   STEPS   '

C                     123456789_123456789_1
C                     123   123456   123456
      header(2, 1) = 'EARCH                '
      header(2, 2) = '---------------------'
      header(2, 3) = 'J''S   COND J   REMARK'

      if (mess .and. 0 .lt. text) then
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

C///////////////////////////////////////////////////////////////////////
C
C     TIME EVOLUTION.
C
C///////////////////////////////////////////////////////////////////////

C///  0 < M?

      if (0 .lt. step) go to 1010
         stride = strid0
         age = 0

C        RETAIN THE LATEST SOLUTION FOR USE BY THE FUNCTION
         call twcopy (groupa + comps * points + groupb, v0, buffer)
         signal = 'RETAIN'
C        GO TO 1010 WHEN ROUTE = 1
         route = 1
         go to 99999
1010  continue
      signal = ' '

C///  FIRST := STEP, LAST := STEP + DESIRE.

      exist = .false.
      first = step
      last = step + desire

C///  PRINT.

      if (.not. (0 .lt. levelm .and. 0 .lt. text)) go to 1030
         call twcopy (groupa + comps * points + groupb, v0, buffer)
         signal = 'RESIDUAL'
         time = .false.
C        GO TO 1020 WHEN ROUTE = 2
         route = 2
         go to 99999
1020     continue
         signal = ' '
         call twnorm (groupa + comps * points + groupb, ynorm, buffer)
         call twlogr (yword, ynorm)

         if (1 .eq. levelm) then
            if (step .eq. 0) then
               write (text, 10001) id, header, step, yword
            else
               write (text, 10002) id, header, step, yword
            end if
         else if (1 .lt. levelm .and. 0 .eq. step) then
            if (step .eq. 0) then
               write (text, 20001) id, step, yword, log10 (stride)
            else
               write (text, 20002) id, step, yword, log10 (stride)
            end if
         end if
1030  continue

C///  LOW := TMIN, HIGH := TMAX.

1040  continue

      low = tmin
      high = tmax

C///  IF AGE = STEPS2 AND STRIDE < HIGH AND 1 < TINC, THEN INCREASE
C///  STRIDE.

      if (age .eq. steps2 .and. stride .lt. high .and. 1.0 .lt. tinc)
     +   then
         age = 0
         exist = .false.
         low = stride * tdec
         stride = min (high, stride * tinc)
         if (1 .lt. levelm .and. 0 .lt. text)
     +      write (text, 20003) id, step, yword, log10 (stride)
      else
         if (1 .lt. levelm .and. 0 .lt. text .and. 0 .lt. step)
     +      write (text, 20002) id, step, yword, log10 (stride)
      end if

C///  NEWTON SEARCH.

1050  continue

C     STORE THE LATEST SOLUTION SHOULD THE SEARCH FAIL
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
         if (csave .eq. 0.0) then
            write (jword, '(I3, 3X, A6)') count, '    NA'
         else
            write (jword, '(I3, 3X, F6.2)') count, log10 (csave)
         end if
      end if

C     SUBROUTINE SEARCH
C    +  (ERROR, TEXT,
C    +   ABOVE, AGE, BELOW, BUFFER, COMPS, CONDIT, EXIST, GROUPA,
C    +   GROUPB, LEVELD, LEVELM, NAME, NAMES, POINTS, REPORT, S0, S1,
C    +   SIGNAL, STEPS, SUCCES, V0, V1, XXABS, XXAGE, XXREL, Y0, Y0NORM,
C    +   Y1)

      call search
     +  (error, text,
     +   above, agej, below, buffer, comps, condit, exist, groupa,
     +   groupb, leveld - 1, levelm - 1, name, names, points, xrepor,
     +   s0, s1, signal, number, xsucce, v0, v1, tdabs, tdage, tdrel,
     +   y0, dummy, y1)
      if (error) go to 9009

      if (signal .ne. ' ') then
         jacob = signal .eq. 'PREPARE'
         time = .true.
C        GO TO 1060 WHEN ROUTE = 3
         route = 3
         go to 99999
      end if

C///  UNSUCCESSFUL?

      if (.not. xsucce) then
         if (1 .eq. levelm .and. 0 .lt. text) then
            if (xrepor .eq. qbnds) then
               length = 6
               remark = 'BOUNDS'
            else if (xrepor .eq. qdvrg) then
               length = 7
               remark = 'DIVERGE'
            else
               length = 1
               remark = ' '
            end if
            write (text, 10003) step + 1, log10 (stride), number, jword,
     +         remark (1 : length)
         end if

C///  IF ALSO LOW < STRIDE AND 1 < TDEC, THEN DECREASE STRIDE.

         if (low .lt. stride .and. 1.0 .lt. tdec) then
            age = 0
            call twcopy (groupa + comps * points + groupb, vsave, v0)
            exist = .false.
            high = stride / tinc
            stride = max (low, stride / tdec)
            if (1 .lt. levelm .and. 0 .lt. text) write (text, 20004)
     +         id, step, yword, log10 (stride)
            go to 1050
         end if

C///  OTHERWISE END, FAILURE.

         go to 2010
      end if

C///  IF NO CHANGE AND STRIDE .LT. HIGH AND 1.0 .LT. TINC, THEN
C///  INCREASE STRIDE.  OTHERWISE END, FAILURE.

      do 1070 j = 1, groupa + comps * points + groupb
         buffer(j) = v0(j) - vsave(j)
1070  continue
      call twnorm (groupa + comps * points + groupb, change, buffer)
      call twlogr (cword, change)

      if (change .eq. 0.0) then
         if (1 .eq. levelm .and. 0 .lt. text) then
            write (text, 10004)
     +         step + 1, '  ZERO', log10 (stride), number, jword
         end if

         if (1.0 .lt. tinc .and. stride .lt. high) then
            age = 0
            exist = .false.
            low = stride * tdec
            stride = min (high, stride * tinc)
            if (1 .lt. levelm .and. 0 .lt. text)
     +         write (text, 20005) id, step, yword, log10 (stride)
            go to 1050
         end if
         go to 2010
      end if

C///  AGE := AGE + 1, M := M + 1.

      age = age + 1
      step = step + 1

C     RETAIN THE LATEST SOLUTION FOR USE BY THE FUNCTION.
C     GO TO 1080 WHEN ROUTE = 4
      route = 4
      call twcopy (groupa + comps * points + groupb, v0, buffer)
      signal = 'RETAIN'
      go to 99999
1080  continue
      signal = ' '

C///  PRINT.

      if (.not. (0 .lt. levelm .and. 0 .lt. text)) go to 1100
         call twcopy (groupa + comps * points + groupb, v0, buffer)
         signal = 'RESIDUAL'
         time = .false.
C        GO TO 1090 WHEN ROUTE = 5
         route = 5
         go to 99999
1090     continue
         signal = ' '
         call twnorm (groupa + comps * points + groupb, ynorm, buffer)
         call twlogr (yword, ynorm)

         if (1 .eq. levelm) write (text, 10005)
     +      step, yword, cword, log10 (stride), number, jword
1100  continue

C///  M < LAST?

      if (step .lt. last) go to 1040

C///////////////////////////////////////////////////////////////////////
C
C     EPILOGUE.
C
C///////////////////////////////////////////////////////////////////////

2010  continue

C///  PRINT.

      if (0 .lt. levelm .and. 0 .lt. text) then
         if (1 .eq. levelm) then
            if (step .eq. first) then
               write (text, 10006) id
            else if (step .eq. last) then
               write (text, 10007) id
            else
               write (text, 10008) id
            end if
         else if (1 .lt. levelm) then
            if (step .eq. first) then
               write (text, 10006) id
            else if (step .eq. last) then
               write (text, 20006) id, step, yword
            else
               write (text, 20007) id, step, yword
            end if
         end if

         if (first .lt. last .and. 1 .eq. leveld) then
            write (text, 20008) id
            call twcopy (groupa + comps * points + groupb, v0, buffer)
            signal = 'SHOW'
C           GO TO 2020 WHEN ROUTE = 6
            route = 6
            go to 99999
         end if
      end if

2020  continue
      signal = ' '

C///  SET THE COMPLETION STATUS FLAGS.

      succes = first .lt. step
      if (step .lt. last) report = xrepor

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 format
     +  (/1X, a9, 'BEGIN TIME EVOLUTION.'
     +  /3(/10X, a44, a21)
     +  /10X, i6, 3X, a6)

10002 format
     +  (/1X, a9, 'CONTINUE TIME EVOLUTION.'
     +  /3(/10X, a44, a21)
     +  /10X, i6, 3X, a6)

10003 format
     +  (10X, i6, 21X, f6.2, 3X, i5, 3X, a12, 3X, a)

10004 format
     +  (10X, i6, 12X, a6, 3X, f6.2, 3X, i5, 3X, a12)

10005 format
     +  (10X, i6, 2(3X, a6), 3X, f6.2, 3X, i5, 3X, a12)

10006 format
     +  (/1X, a9, 'FAILURE.  NO TIME EVOLUTION.')

10007 format
     +  (/1X, a9, 'SUCCESS.  TIME EVOLUTION COMPLETED.')

10008 format
     +  (/1X, a9, 'PARTIAL SUCCESS.  TIME EVOLUTION INCOMPLETE.')

20001 format
     +  (/1X, a9, 'BEGIN TIME EVOLUTION.'
     + //10X, i10, '  LATEST TIME POINT'
     +  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE'
     +  /10X, f10.2, '  LOG10 STRIDE TO NEXT TIME POINT'
     + //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')

20002 format
     +  (/1X, a9, 'CONTINUE TIME EVOLUTION.'
     + //10X, i10, '  LATEST TIME POINT'
     +  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE'
     +  /10X, f10.2, '  LOG10 STRIDE TO NEXT TIME POINT'
     + //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')

20003 format
     +  (/1X, a9, 'CONTINUE TIME EVOLUTION WITH INCREASED STRIDE.'
     + //10X, i10, '  LATEST TIME POINT'
     +  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE'
     +  /10X, f10.2, '  LOG10 INCREASED STRIDE TO NEXT TIME POINT'
     + //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE.')

20004 format
     +  (/1X, a9, 'RETRY THE STEP WITH A DECREASED TIME STRIDE.'
     + //10X, i10, '  LATEST TIME POINT'
     +  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE'
     +  /10X, f10.2, '  LOG10 DECREASED STRIDE TO NEXT TIME POINT'
     + //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE, AGAIN.')

20005 format
     +  (/1X, a9, 'THE SOLUTION DID NOT CHANGE.  RETRYING THE STEP'
     +  /10X, 'WITH AN INCREASED TIME STRIDE.'
     + //10X, i10, '  LATEST TIME POINT'
     +  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE'
     +  /10X, f10.2, '  LOG10 INCREASED STRIDE TO NEXT TIME POINT'
     + //10X, 'SEARCHING FOR THE NEXT TRANSIENT STATE, AGAIN.')

20006 format
     +  (/1X, a9, 'SUCCESS.  TIME EVOLUTION COMPLETED.'
     + //10X, i10, '  LAST TIME POINT'
     +  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE')

20007 format
     +  (/1X, a9, 'PARTIAL SUCCESS.  TIME EVOLUTION INCOMPLETE.'
     + //10X, i10, '  LAST TIME POINT'
     +  /14X, a6, '  LOG10 STEADY STATE RESIDUAL HERE')

20008 format
     + (/1X, a9, 'THE LATEST SOLUTION:')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      go to 99999

9001  if (0 .lt. text) write (text, 99001) id, route
      if (.not. mess) go to 99999

9002  if (0 .lt. text) write (text, 99002) id,
     +   comps, points, groupa, groupb, groupa + comps * points + groupb
      if (.not. mess) go to 99999

9003  if (0 .lt. text) write (text, 99003) id, desire
      if (.not. mess) go to 99999

9004  if (0 .lt. text) write (text, 99004) id, tdec, tinc
      if (.not. mess) go to 99999

9005  if (0 .lt. text) write (text, 99005) id, tmin, tmax
      if (.not. mess) go to 99999

9006  if (0 .lt. text) write (text, 99006) id, tmin, strid0, tmax
      if (.not. mess) go to 99999

9007  if (0 .lt. text) write (text, 99007) id, step
      if (.not. mess) go to 99999

9008  if (0 .lt. text) write (text, 99008) id, steps2
      if (.not. mess) go to 99999

9009  if (0 .lt. text) write (text, 99009) id
      if (.not. mess) go to 99999

99001 format
     +  (/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, i10, '  ROUTE')

99002 format
     +  (/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES'
     +  /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  TOTAL UNKNOWNS')

99003 format
     +  (/1X, a9, 'ERROR.  THE NUMBER OF TIME STEPS MUST BE POSITIVE.'
     + //10X, i10, '  STEPS0 OR STEPS1, DESIRED NUMBER OF STEPS')

99004 format
     +  (/1X, a9, 'ERROR.  THE FACTORS FOR CHANGING THE TIME STRIDE'
     +  /10X, 'MUST BE NO SMALLER THAN 1.'
     + //10X, 1p, e10.2, '  TDEC, DECREASE FACTOR',
     +  /10X, 1p, e10.2, '  TINC, INCREASE FACTOR')

99005 format
     +  (/1X, a9, 'ERROR.  THE BOUNDS ON THE TIME STRIDE ARE OUT OF'
     +  /10X, 'ORDER.'
     + //10X, 1p, e10.2, '  TMIN, SHORTEST STRIDE'
     +  /10X, 1p, e10.2, '  TMAX, LONGEST STRIDE')

99006 format
     +  (/1X, a9, 'ERROR.  THE INITIAL TIME STRIDE MUST LIE BETWEEN'
     +  /10X, 'THE LOWER AND UPPER BOUNDS.'
     + //10X, 1p, e10.2, '  TMIN, SHORTEST STRIDE'
     +  /10X, 1p, e10.2, '  STRID0, INITIAL STRIDE'
     +  /10X, 1p, e10.2, '  TMAX, LONGEST STRIDE')

99007 format
     +  (/1X, a9, 'ERROR.  THE COUNT OF TIME STEPS MUST BE ZERO OR'
     +  /10X, 'POSITIVE.'
     + //10X, i10, '  STEP')

99008 format
     +  (/1X, a9, 'ERROR.  THE TIME STEPS BEFORE STRIDE INCREASES'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, i10, '  STEPS2, TIME STEPS BEFORE STRIDE INCREASES')

99009 format
     +  (/1X, a9, 'ERROR.  SEARCH FAILS.')

C///  EXIT.

      stop
99999 continue
      return
      end
      subroutine search
     +  (error, text,
     +   above, age, below, buffer, comps, condit, exist, groupa,
     +   groupb, leveld, levelm, name, names, points, report, s0, s1,
     +   signal, steps, succes, v0, v1, xxabs, xxage, xxrel, y0, y0norm,
     +   y1)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     SEARCH
C
C     PERFORM THE DAMPED, MODIFIED NEWTON'S SEARCH.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

      character
     +   column*16, ctemp1*80, ctemp2*80, header*80, id*9, name*(*),
     +   signal*(*), string*80
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   above, abs0, abs1, below, buffer, condit, deltab, deltad, rel0,
     +   rel1, s0, s0norm, s1, s1norm, sj, temp, v0, v1, value, vj,
     +   xxabs, xxrel, y0, y0norm, y1, y1norm, zero
      external
     +   twcopy, twlogr, twnorm, twsqez
      integer
     +   age, comps, count, entry, expone, groupa, groupb, i, j, k,
     +   len1, len2, length, leveld, levelm, lines, names, number,
     +   points, qbnds, qdvrg, qnull, report, route, steps, text, xxage
      intrinsic
     +   abs, int, log10, max, min, mod
      logical
     +   error, exist, force, mess, succes

      parameter (id = 'SEARCH:  ')
      parameter (lines = 20)
      parameter (zero = 0.0)

C     REPORT CODES
      parameter (qnull = 0, qbnds = 1, qdvrg = 2)

      dimension
     +   above(groupa + comps * points + groupb),
     +   below(groupa + comps * points + groupb),
     +   buffer(groupa + comps * points + groupb), column(7),
     +   header(3, 2), name(names),
     +   s0(groupa + comps * points + groupb),
     +   s1(groupa + comps * points + groupb),
     +   v0(groupa + comps * points + groupb),
     +   v1(groupa + comps * points + groupb),
     +   y0(groupa + comps * points + groupb),
     +   y1(groupa + comps * points + groupb)

C///  SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.

      save

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  EVERY-TIME INITIALIZATION.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

C     TURN OFF ALL COMPLETION STATUS FLAGS.
      error = .false.
      report = qnull
      succes = .false.

C///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      if (signal .ne. ' ') then
         go to (2020, 2040, 2050, 2140, 2150, 2180) route
         error = .true.
         go to 9001
      end if

C///  ONE-TIME INITIALIZATION.

      number = 0

C///  CHECK THE ARGUMENTS.

      error = .not. (((0 .lt. comps) .eqv. (0 .lt. points)) .and.
     +   0 .le. comps .and. 0 .le. points .and. 0 .le. groupa .and.
     +   0 .le. groupb .and. 0 .lt. groupa + comps * points + groupb)
      if (error) go to 9002

      error = .not. (names .eq. 1 .or.
     +   names .eq. groupa + comps + groupb)
      if (error) go to 9003

      count = 0
      do 1010 j = 1, groupa + comps * points + groupb
         if (.not. (below(j) .lt. above(j))) count = count + 1
1010  continue
      error = count .ne. 0
      if (error) go to 9004

      count = 0
      do 1020 j = 1, groupa + comps * points + groupb
         if (.not. (below(j) .le. v0(j) .and. v0(j) .le. above(j)))
     +      count = count + 1
1020  continue
      error = count .ne. 0
      if (error) go to 9005

      error = .not. (0.0 .le. xxabs .and. 0.0 .le. xxrel)
      if (error) go to 9006

      error = .not. (0 .lt. xxage)
      if (error) go to 9007

C///  WRITE ALL MESSAGES.

      if (mess .and. 0 .lt. text) then
         route = 0

         write (text, 10003) id
         write (text, 10002) id
         count = 0
         do 1030 j = 1, groupa + comps * points + groupb
            count = count + 1
            if (count .le. lines) then
               if (j .le. groupa) then
                  i = j
               else if (j .le. groupa + comps * points) then
                  i = groupa + mod (j - groupa - 1, comps) + 1
               else
                  i = j - groupa - comps * points
               end if

               if (names .eq. comps + groupa + groupb) then
                  ctemp1 = name(i)
               else
                  ctemp1 = ' '
               end if
               call twsqez (len1, ctemp1)

               if (j .le. groupa) then
                  write (ctemp2, 80001) 'A', i
               else if (j .le. groupa + comps * points) then
                  write (ctemp2, 80002) 'C', i,
     +               'P', int ((j - groupa - 1) / comps) + 1
               else
                  write (ctemp2, 80001) 'B', i
               end if
               call twsqez (len2, ctemp2)

               if (ctemp1 .eq. ' ') then
                  string = ctemp2
                  length = len2
               else if (len1 + 2 + len2 .le. 30) then
                  string = ctemp1 (1 : len1) // '  ' // ctemp2
                  length = len1 + 2 + len2
               else if (len1 + 1 + len2 .le. 30) then
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
         if (lines .lt. count) write (text, 80004)
         write (text, 10001) id
         write (text, 10006) id
         write (text, 10005) id

         go to 9001
      end if

C///  PRINT THE HEADER.

C                     123456789_123456789_123456789_123456789_1234
C                     123456   123456   123456   123456   123456
      header(1, 1) = '         LOG10                              '
      header(2, 1) = '  SLTN   -----------------------------------'
      header(3, 1) = 'NUMBER   NORM F   COND J   NORM S      ABS A'

C                     123456789_123456789_123
C                     123456   123456  123456
      header(1, 2) = '                       '
      header(2, 2) = '-----------------------'
      header(3, 2) = 'ND REL    DELTA B AND D'

      if (levelm .ge. 1 .or. mess) then
         if (0 .lt. text) write (text, 10001)
     +      id, ((header(j, k), k = 1, 2), j = 1, 3)
      end if

C///////////////////////////////////////////////////////////////////////
C
C     SIR ISSAC NEWTON'S ALGORITHM.
C
C///////////////////////////////////////////////////////////////////////

C///  J EXIST?

      if (.not. exist) go to 2010

C///  AGE < XXAGE?

      if (age .lt. xxage) go to 2030

C///  EVALUATE J AT V0.  RE-EVALUATE Y0 := F(V0) IN CASE F CHANGES WHEN
C///  J DOES.  SOLVE J S0 = Y0.  EVAUATE ABS0 AND REL0.

2010  continue

      call twcopy (groupa + comps * points + groupb, v0, buffer)
      signal = 'PREPARE'
C     GO TO 2020 WHEN ROUTE = 1
      route = 1
      go to 99999
2020  continue
      signal = ' '
      age = 0

C     JACOBIAN EVALUATION SHOULD RETURN A NEW RESIDUAL TOO.

      if (0 .lt. levelm .and. 0 .lt. text) then
         if (0.0 .lt. condit) then
            write (column(2), '(F6.2)') log10 (condit)
         else
            column(2) = '    NA'
         end if
      end if

C///  EVALUATE Y0 := F(V0).  SOLVE J S0 = Y0.  EVAUATE ABS0 AND REL0.

2030  continue

      call twcopy (groupa + comps * points + groupb, v0, buffer)
      signal = 'RESIDUAL'
C     GO TO 2040 WHEN ROUTE = 2
      route = 2
      go to 99999
2040  continue
      signal = ' '
      call twcopy (groupa + comps * points + groupb, buffer, y0)
      call twnorm (groupa + comps * points + groupb, y0norm, y0)

      call twcopy (groupa + comps * points + groupb, y0, buffer)
      signal = 'SOLVE'
C     GO TO 2050 WHEN ROUTE = 3
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
         if (xxrel * vj .lt. sj) abs0 = max (abs0, sj)
         if (xxabs .lt. sj .and. 0.0 .lt. vj)
     +      rel0 = max (rel0, sj / vj)
2060  continue

C///  CHECK FOR SUCCESS.

      if (abs0 .le. xxabs .and. rel0 .le. xxrel) go to 2170

C///  CHOOSE DELTAB.

2070  continue

C     DELTAB IS THE LARGEST DAMPING COEFFICIENT BETWEEN 0 AND 1 THAT
C     KEEPS V1 WITHIN BOUNDS.  IF V1 BELONGS ON THE BOUNDARY, THEN
C     PROVISIONS ARE MADE TO FORCE IT THERE DESPITE ROUNDING ERROR.

      deltab = 1.0
      force = .false.
      do 2080 j = 1, groupa + comps * points + groupb
         if (s0(j) .gt. max (zero, v0(j) - below(j))) then
            temp = (v0(j) - below(j)) / s0(j)
            if (temp .lt. deltab) then
               deltab = temp
               entry = j
               force = .true.
               value = below(j)
            end if
         else if (s0(j) .lt. min (zero, v0(j) - above(j))) then
            temp = (v0(j) - above(j)) / s0(j)
            if (temp .lt. deltab) then
               deltab = temp
               entry = j
               force = .true.
               value = above(j)
            end if
         end if
2080  continue

      error = deltab .lt. 0.0
      if (error) go to 9008

C///  0 < DELTAB?

      if (.not. (0.0 .lt. deltab)) then
         if (0 .lt. age) go to 2010

         if (0 .lt. levelm .and. 0 .lt. text) then
            call twlogr (column(1), y0norm)
            call twlogr (column(3), s0norm)
            call twlogr (column(4), abs0)
            call twlogr (column(5), rel0)
            column(6) = ' '
            if (deltab .ne. 1.0) call twlogr (column(6), deltab)
            column(7) = ' '
            if (deltad .ne. 1.0) call twlogr (column(7), deltad)
            write (text, 10004) number, column
            write (text, 10002) id

            count = 0
            do 2090 j = 1, groupa + comps * points + groupb
               if ((below(j) .eq. v0(j) .and. 0.0 .lt. s0(j)) .or.
     +            (v0(j) .eq. above(j) .and. s0(j) .lt. 0.0)) then
                  count = count + 1
                  if (count .le. lines) then
                     if (j .le. groupa) then
                        i = j
                     else if (j .le. groupa + comps * points) then
                        i = groupa + mod (j - groupa - 1, comps) + 1
                     else
                        i = j - groupa - comps * points
                     end if

                     if (names .eq. comps + groupa + groupb) then
                        ctemp1 = name(i)
                     else
                        ctemp1 = ' '
                     end if
                     call twsqez (len1, ctemp1)

                     if (j .le. groupa) then
                        write (ctemp2, 80001) 'A', i
                     else if (j .le. groupa + comps * points) then
                        write (ctemp2, 80002) 'C', i,
     +                     'P', int ((j - groupa - 1) / comps) + 1
                     else
                        write (ctemp2, 80001) 'B', i
                     end if
                     call twsqez (len2, ctemp2)

                     if (ctemp1 .eq. ' ') then
                        string = ctemp2
                        length = len2
                     else if (len1 + 2 + len2 .le. 30) then
                        string = ctemp1 (1 : len1) // '  ' // ctemp2
                        length = len1 + 2 + len2
                     else if (len1 + 1 + len2 .le. 30) then
                        string = ctemp1 (1 : len1) // ' ' // ctemp2
                        length = len1 + 1 + len2
                     else
                        len1 = 30 - len2 - 4
                        string = ctemp1 (1 : len1) // '... ' // ctemp2
                        length = 30
                     end if

                     if (below(j) .eq. v0(j)) then
                        write (text, 80003)
     +                     'LOWER', v0(j), string (1 : length)
                     else
                        write (text, 80003)
     +                     'UPPER', v0(j), string (1 : length)
                     end if
                  end if
               end if
2090        continue
            if (lines .lt. count) write (text, 80005)
         end if

         report = qbnds
         succes = .false.
         go to 99999
      end if

C///  DELTAD := 1.

      deltad = 1.0
      expone = 0

C///  V1 := V0 - DELTAB DELTAD S0.  EVALUATE Y1 := F(V1).  SOLVE
C///  J S1 = Y1.  EVALUATE ABS1 AND REL1.

2100  continue

      temp = deltab * deltad
      do 2110 j = 1, groupa + comps * points + groupb
         v1(j) = v0(j) - temp * s0(j)
2110  continue

C     KEEP V1 IN BOUNDS DESPITE ROUNDING ERROR.

      do 2120 j = 1, groupa + comps * points + groupb
         v1(j) = min (v1(j), above(j))
2120  continue
      do 2130 j = 1, groupa + comps * points + groupb
         v1(j) = max (v1(j), below(j))
2130  continue
      if (expone .eq. 0 .and. force) v1(entry) = value

      call twcopy (groupa + comps * points + groupb, v1, buffer)
      signal = 'RESIDUAL'
C     GO TO 2140 WHEN ROUTE = 4
      route = 4
      go to 99999
2140  continue
      signal = ' '
      call twcopy (groupa + comps * points + groupb, buffer, y1)
      call twnorm (groupa + comps * points + groupb, y1norm, y1)

      call twcopy (groupa + comps * points + groupb, y1, buffer)
      signal = 'SOLVE'
C     GO TO 2150 WHEN ROUTE = 5
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
         if (xxrel * vj .lt. sj) abs1 = max (abs1, sj)
         if (xxabs .lt. sj .and. 0.0 .lt. vj)
     +      rel1 = max (rel1, sj / vj)
2160  continue

C///  NORM S1 < OR = NORM S0?

      if (s1norm .le. s0norm) then
      else
         deltad = 0.5 * deltad
         expone = expone + 1
         if (expone .le. 5) go to 2100
            if (0 .lt. age) go to 2010
               if (0 .lt. levelm .and. 0 .lt. text) then
                  call twlogr (column(1), y0norm)
                  call twlogr (column(3), s0norm)
                  call twlogr (column(4), abs0)
                  call twlogr (column(5), rel0)
                  column(6) = ' '
                  if (deltab .ne. 1.0) call twlogr (column(6), deltab)
                  column(7) = ' '
                  if (deltad .ne. 1.0) call twlogr (column(7), deltad)
                  write (text, 10004) number, column
                  write (text, 10003) id
               end if
               report = qdvrg
               succes = .false.
               go to 99999
      end if

C///  PRINT.

      if (0 .lt. levelm .and. 0 .lt. text) then
         call twlogr (column(1), y0norm)
         call twlogr (column(3), s0norm)
         call twlogr (column(4), abs0)
         call twlogr (column(5), rel0)
         column(6) = ' '
         if (deltab .ne. 1.0) call twlogr (column(6), deltab)
         column(7) = ' '
         if (deltad .ne. 1.0) call twlogr (column(7), deltad)
         write (text, 10004) number, column
         column(2) = ' '
      end if

C///  S0 := S1, U := V1, Y0 := Y1, AGE := AGE + 1.

      age = age + 1
      number = number + 1
      call twcopy (groupa + comps * points + groupb, s1, s0)
      call twcopy (groupa + comps * points + groupb, v1, v0)
      call twcopy (groupa + comps * points + groupb, y1, y0)
      s0norm = s1norm
      y0norm = y1norm
      abs0 = abs1
      rel0 = rel1

C///  S0 SMALL VS V0?

      if (.not. (abs0 .le. xxabs .and. rel0 .le. xxrel)) then
         if (age .lt. xxage) go to 2070
         go to 2010
      end if

C///  SUCCESS.

2170  continue

C///  PRINT.

      if (0 .lt. levelm .and. 0 .lt. text) then
         call twlogr (column(1), y0norm)
         call twlogr (column(3), s0norm)
         call twlogr (column(4), abs0)
         call twlogr (column(5), rel0)
         column(6) = ' '
         column(7) = ' '
         if (0 .lt. leveld) then
            write (text, 10004) number, column
            write (text, 10005) id
            signal = 'SHOW'
            call twcopy (groupa + comps * points + groupb, v0, buffer)
C           GO TO 2180 WHEN ROUTE = 6
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

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 format
     +  (/1X, a9, 'SOLVE NONLINEAR, NONDIFFERENTIAL EQUATIONS.'
     +  /4(/10X, a44, a23)/)

10002 format
     + (/1X, a9, 'FAILURE.  THE SEARCH FOR THE FOLLOWING UNKNOWNS GOES'
     +  /10X, 'OUT OF BOUNDS.'
C              12345  123456789_
     + //10X, 'BOUND       VALUE   UNKNOWN'
     +  /)

10003 format
     +  (/1X, a9, 'FAILURE.  THE SEARCH DIVERGES.')

10004 format
     +  (10X, i6, 3(3X, a6), 2(3X, a6, 2X, a6))

10005 format
     +  (/1X, a9, 'SUCCESS.  THE SOLUTION:')

10006 format
     +  (/1X, a9, 'SUCCESS.')

80001 format
     +   ('(', a, ' ', i10, ')')

80002 format
     +  ('(', a, ' ', i10, ' ', a, ' ', i10, ')')

80003 format
     +  (10X, a5, 2X, 1p, e10.2, 3X, a)

80004 format
     +  (30X, '... MORE')

80005 format
     +  (10X, '  ... MORE')

80006 format
     +  (10X, 1p, e10.2, 2X, e10.2, 3X, a)

80007 format
     +  (10X, 1p, e10.2, 2X, e10.2, 2X, e10.2, 3X, a)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      go to 99999

9001  if (0 .lt. text) write (text, 99001) id, route
      if (.not. mess) go to 99999

9002  if (0 .lt. text) write (text, 99002) id,
     +   comps, points, groupa, groupb, groupa + comps * points + groupb
      if (.not. mess) go to 99999

9003  if (0 .lt. text) write (text, 99003) id,
     +   names, comps, groupa, groupb, groupa + comps + groupb
      if (.not. mess) go to 99999

9004  if (0 .lt. text) then
         write (text, 99004) id,
     +      groupa, groupb, comps, groupa + comps + groupb, count
         count = 0
         do 8010 j = 1, groupa + comps + groupb
            if (.not. (below(j) .lt. above(j)) .or. mess) then
               count = count + 1
               if (count .le. lines) then
                  if (names .eq. comps + groupa + groupb) then
                     ctemp1 = name(j)
                  else
                     ctemp1 = ' '
                  end if
                  call twsqez (len1, ctemp1)

                  if (j .le. groupa) then
                     write (ctemp2, 80001) 'A', j
                  else if (j .le. groupa + comps) then
                     write (ctemp2, 80001) 'C', j - groupa
                  else
                     write (ctemp2, 80001) 'B', j - groupa - comps
                  end if
                  call twsqez (len2, ctemp2)

                  if (ctemp1 .eq. ' ') then
                     string = ctemp2
                     length = len2
                  else if (len1 + 2 + len2 .le. 40) then
                     string = ctemp1 (1 : len1) // '  ' // ctemp2
                     length = len1 + 2 + len2
                  else if (len1 + 1 + len2 .le. 40) then
                     string = ctemp1 (1 : len1) // ' ' // ctemp2
                     length = len1 + 1 + len2
                  else
                     len1 = 40 - len2 - 4
                     string = ctemp1 (1 : len1) // '... ' // ctemp2
                     length = 40
                  end if

                  write (text, 80006)
     +               below(j), above(j), string (1 : length)
               end if
            end if
8010     continue
         if (lines .lt. count) write (text, 80005)
      end if
      if (.not. mess) go to 99999

9005  if (0 .lt. text) then
         write (text, 99005) id, groupa, groupb, comps, points,
     +      groupa + comps * points + groupb, count
         count = 0
         do 8020 j = 1, groupa + comps * points + groupb
            if (.not. (below(j) .le. v0(j) .and. v0(j) .le. above(j))
     +         .or. mess) then
               count = count + 1
               if (count .le. lines) then
                  if (j .le. groupa) then
                     i = j
                  else if (j .le. groupa + comps * points) then
                     i = groupa + mod (j - groupa - 1, comps) + 1
                  else
                     i = j - groupa - comps * points
                  end if

                  if (names .eq. comps + groupa + groupb) then
                     ctemp1 = name(i)
                  else
                     ctemp1 = ' '
                  end if
                  call twsqez (len1, ctemp1)

                  if (j .le. groupa) then
                     write (ctemp2, 80001) 'A', i
                  else if (j .le. groupa + comps * points) then
                     write (ctemp2, 80002) 'C', i,
     +                  'P', int ((j - groupa - 1) / comps) + 1
                  else
                     write (ctemp2, 80001) 'B', i
                  end if
                  call twsqez (len2, ctemp2)

                  if (ctemp1 .eq. ' ') then
                     string = ctemp2
                     length = len2
                  else if (len1 + 2 + len2 .le. 30) then
                     string = ctemp1 (1 : len1) // '  ' // ctemp2
                     length = len1 + 2 + len2
                  else if (len1 + 1 + len2 .le. 30) then
                     string = ctemp1 (1 : len1) // ' ' // ctemp2
                     length = len1 + 1 + len2
                  else
                     len1 = 30 - len2 - 4
                     string = ctemp1 (1 : len1) // '... ' // ctemp2
                     length = 30
                  end if

                  write (text, 80007)
     +               below(j), v0(j), above(j), string (1 : length)
               end if
            end if
8020     continue
         if (lines .lt. count) write (text, 80005)
      end if
      if (.not. mess) go to 99999

9006  if (0 .lt. text) write (text, 99006) id, xxabs, xxrel
      if (.not. mess) go to 99999

9007  if (0 .lt. text) write (text, 99007) id, xxage
      if (.not. mess) go to 99999

9008  if (0 .lt. text) write (text, 99008) id, deltab
      if (.not. mess) go to 99999

99001 format
     +  (/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, i10, '  ROUTE')

99002 format
     +  (/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES'
     +  /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  TOTAL UNKNOWNS')

99003 format
     +  (/1X, a9, 'ERROR.  THE NUMBER OF NAMES IS WRONG.'
     + //10X, i10, '  NAMES'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  TOTAL NUMBER')

99004 format
     +  (/1X, a9, 'ERROR.  THE LOWER AND UPPER BOUNDS ON SOME UNKNOWNS'
     +  /10X, 'ARE OUT OF ORDER.'
     + //10X, i10, '  GROUP A UNKNOWNS (A)'
     +  /10X, i10, '  GROUP B UNKNOWNS (B)'
     +  /10X, i10, '  COMPONENTS AT POINTS (C)'
     +  /10X, i10, '  TOTAL TYPES OF UNKNOWNS'
     +  /10X, i10, '  NUMBER OF BOUNDS OUT OF ORDER'
C              123456789_  123456789_
     + //10X, '     LOWER       UPPER'
     +  /10X, '     BOUND       BOUND   UNKNOWN'
     +  /)

99005 format
     +  (/1X, a9, 'ERROR.  THE GUESSES FOR SOME UNKNOWNS ARE OUT OF'
     +  /10X, 'BOUNDS.'
     + //10X, i10, '  GROUP A UNKNOWNS (A)'
     +  /10X, i10, '  GROUP B UNKNOWNS (B)'
     +  /10X, i10, '  COMPONENTS AT POINTS (C)'
     +  /10X, i10, '  POINTS (P)'
     +  /10X, i10, '  TOTAL UNKNOWNS'
     +  /10X, i10, '  NUMBER OUT OF BOUNDS'
C              123456789_  123456789_  123456789_
     + //10X, '     LOWER                   UPPER'
     +  /10X, '     BOUND       VALUE       BOUND   UNKNOWN'
     +  /)

99006 format
     +  (/1X, a9, 'ERROR.  THE BOUNDS FOR THE ABSOLUTE AND RELATIVE'
     +  /10X, 'CONVERGENCE TESTS MUST BE ZERO OR POSITIVE.'
     + //10X, 1p, e10.2, '  SSABS OR TDABS, ABSOLUTE ERROR'
     +  /10X, 1p, e10.2, '  SSREL OR TDREL, RELATIVE ERROR')

99007 format
     +  (/1X, a9, 'ERROR.  THE RETIREMENT AGE OF THE JACOBIAN MATRIX'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, i10, '  SSAGE OR TDAGE, MATRIX RETIREMENT AGE')

99008 format
     +  (/1X, a9, 'ERROR.  THE DAMPING COEFFICIENT FOR STAYING'
     +  /10X, 'IN BOUNDS IS NEGATIVE.'
     + //10X, 1p, e10.2, '  DELTA B')

C///  EXIT.

      stop
99999 continue

C     COPY THE PROTECTED LOCAL VARIABLE
      steps = number

      return
      end
      subroutine twcopy (n, x, y)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWCOPY
C
C     COPY ONE VECTOR TO ANOTHER.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      integer j, n
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   x, y

      dimension x(n), y(n)

      do 0100 j = 1, n
         y(j) = x(j)
0100  continue

      return
      end
      subroutine tweps (eps)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWEPS
C
C     FIND MACHINE EPSILON.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   eps, value
      integer
     +   scrtch
C*****MACHINE EPSILON > COMPUTED
C      LOGICAL
C     +   SAME
C*****end MACHINE EPSILON > COMPUTED

      parameter (scrtch = 98)

C///  IEEE STANDARD

C*****PRECISION > DOUBLE
      value = 1.1102230246251565D-16
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      VALUE = 5.9604645E-08
C*****END PRECISION > SINGLE

C*****MACHINE EPSILON > IEEE STANDARD
      eps = value
C*****END MACHINE EPSILON > IEEE STANDARD

C///  COMPUTED

C*****MACHINE EPSILON > COMPUTED
C      OPEN (ACCESS = 'SEQUENTIAL', FORM = 'UNFORMATTED',
C     +   STATUS = 'SCRATCH', UNIT = SCRTCH)
C
C      EPS = 1
C1010  CONTINUE
C      EPS = 0.5 * EPS
C
C      VALUE = 1 + EPS
C
C      REWIND (SCRTCH)
C      WRITE (SCRTCH) VALUE
C
C      REWIND (SCRTCH)
C      READ (SCRTCH) VALUE
C
C      SAME = 1 .EQ. VALUE
C
C      IF (.NOT. SAME) GO TO 1010
C
C      CLOSE (UNIT = SCRTCH)
C*****END MACHINE EPSILON > COMPUTED

C///  EXIT.

      return
      end
      subroutine twgbco (a, lda, n, lower, upper, pivot, rcond, z)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWGBCO
C
C     FACTOR A BANDED MATRIX AND ESTIMATE THE RECIPROCAL OF ITS
C     CONDITION NUMBER.  BASED ON _GBCO FROM THE LINPACK LIBRARY.
C
C///////////////////////////////////////////////////////////////////////

      double precision
     +   dsum
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   a, anorm, ek, rcond, s, sm, sum, t, wk, wkm, ynorm, z
      external
     +   twgbfa
      integer
     +   first, info, j, jdiag, ju, k, last, lda, lower, mm, n, pivot,
     +   upper
      intrinsic
     +   abs, dble, max, min, sign

      dimension
     +   a(lda,n), pivot(n), z(n)

      jdiag = lower + upper + 1

C///  COMPUTE THE 1-NORM OF A

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

C///  FACTOR A

      call twgbfa (a, lda, n, lower, upper, pivot, info)

C///  SOLVE TRANSPOSE(U) * W = E

      ek = 1.0
      do 2010 j = 1, n
         z(j) = 0.0
2010  continue

      ju = 0
      do 2050 k = 1, n
         if (z(k) .ne. 0.0) ek = sign (ek, - z(k))

         if (abs (ek - z(k)) .gt. abs (a(jdiag, k))) then
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
         if (a(jdiag, k) .ne. 0.0) then
            wk = wk / a(jdiag, k)
            wkm = wkm / a(jdiag, k)
         else
            wk = 1.0
            wkm = 1.0
         end if

         ju = min (max (ju, upper + pivot(k)), n)
         mm = jdiag
         if (k + 1 .le. ju) then
            do 2030 j = k + 1, ju
               mm = mm - 1
               sm = sm + abs (z(j) + wkm * a(mm, j))
               z(j) = z(j) + wk * a(mm, j)
               s = s + abs (z(j))
2030        continue

            if (s .lt. sm) then
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

C///  SOLVE TRANSPOSE(L) * Y = W

      do 3030 k = n, 1, - 1
         dsum = 0.0
         do 3010 j = 1, min (lower, n - k)
            dsum = dsum + dble (a(jdiag + j, k)) * dble (z(k + j))
3010     continue
         z(k) = z(k) + dsum

         if (1.0 .lt. abs (z(k))) then
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

C///  SOLVE L * V = Y

      do 4030 k = 1, n
         j = pivot(k)
         t = z(j)
         z(j) = z(k)
         z(k) = t

         do 4010 j = 1, min (lower, n - k)
            z(k + j) = t * a(jdiag + j, k) + z(k + j)
4010     continue

         if (1.0 .lt. abs (z(k))) then
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

C///  SOLVE U * Z = W

      do 5030 k = n, 1, - 1
         if (abs (z(k)) .gt. abs (a(jdiag, k))) then
            s = abs (a(jdiag, k)) / abs (z(k))
            do 5010 j = 1, n
               z(j) = s * z(j)
5010        continue
            ynorm = s*ynorm
         end if

         if (a(jdiag, k) .ne. 0.0) then
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

C///  FORM RCOND

      if (anorm .ne. 0.0) then
         rcond = ynorm / anorm
      else
         rcond = 0.0
      end if

C///  EXIT

      return
      end
      subroutine twgbfa (a, lda, n, lower, upper, pivot, info)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWGBFA
C
C     FACTOR A BANDED MATRIX FOR TWGBCO. BASED ON _GBFA FROM THE LINPACK
C     LIBRARY.
C
C///////////////////////////////////////////////////////////////////////

C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   a, t, value
      integer
     +   i, info, j, jk, k, lda, lower, n, pack, pdiag, pivot, pjk,
     +   save, true, upper
      intrinsic
     +   abs, max, min

      dimension
     +   a(lda, n), pivot(n)

C///  STATEMENT FUNTIONS

C     PACKED ROW POSITION OF ENTRY J IN COLUMN K
      pack (j, k) = j - k + pdiag

C     TRUE ROW POSITION OF ENTRY J PACKED IN COLUMN K
      true (j, k) = j - pdiag + k

C///  INITIALIZE

C     PACKED ROW POSITION OF THE DIAGONAL
      pdiag = lower + upper + 1

      info = 0

C///  TOP OF THE LOOP OVER COLUMNS

      do 2060 k = 1, n

C///  INITIALIZE THE FILL-IN SPACE

         do 2010 i = 1, lower
            a(i, k) = 0.0
2010     continue

C///  LOOP OVER THE PREVIOUS COLUMNS

         do 2030 j = max (1, k - lower - upper), k - 1
            pjk = pack (pivot(j), k)
            jk = pack (j, k)
            t = a(pjk, k)
            if (pjk .ne. jk) then
               a(pjk, k) = a(jk, k)
               a(jk, k) = t
            end if

            if (t .ne. 0.0) then
               do 2020 i = 1, min (lower, n - j)
                  a(jk + i, k) = t * a(pdiag + i, j) + a(jk + i, k)
2020           continue
            end if
2030     continue

C///  FIND THE PIVOT

         save = pdiag
         value = abs (a(pdiag, k))
         do 2040 i = pdiag + 1, pdiag + min (lower, n - k)
            if (value .lt. abs (a(i, k))) then
               save = i
               value = abs (a(i, k))
            end if
2040     continue
         pivot(k) = true (save, k)

C///  INTERCHANGE IF NECESSARY

         if (save .ne. pdiag) then
            t = a(save, k)
            a(save, k) = a(pdiag, k)
            a(pdiag, k) = t
         end if

C///  SCALE THE LOWER COLUMN

         if (a(save, k) .ne. 0.0) then
            t = - 1.0 / a(pdiag, k)
            do 2050 i = pdiag + 1, pdiag + min (lower, n - k)
               a(i, k) = t * a(i, k)
2050        continue
         else
            info = k
         end if

C///  BOTTOM OF THE LOOP OVER COLUMNS

2060  continue

C///  THE FINAL COLUMN IS TRIVIAL

      pivot(n) = n
      if (a(pdiag, n) .eq. 0.0) info = n

C///  EXIT

      return
      end
      subroutine twgbsl (abd, lda, n, lower, upper, pivot, b)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWGBSL
C
C     SOLVE A SYSTEM OF LINEAR EQUATIONS FOR TWSOLV.  BASED ON _GBSL
C     FROM THE LINPACK LIBRARY.
C
C///////////////////////////////////////////////////////////////////////

C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   abd, b, t
      integer
     +   j, jdiag, k, l, la, lb, lda, lm, lower, n, pivot, upper
      intrinsic
     +   min

      dimension
     +   abd(lda,*), b(*), pivot(*)

      jdiag = upper + lower + 1

      if (0 .lt. lower) then
         do 1020 k = 1, n - 1
            l = pivot(k)
            t = b(l)
            if (l .ne. k) then
               b(l) = b(k)
               b(k) = t
            end if

            lm = min (lower, n - k)
            do 1010 j = 1, lm
               b(k + j) = t * abd(jdiag + j, k) + b(k + j)
1010        continue
1020     continue
      end if

      do 1040 k = n, 1, - 1
         b(k) = b(k) / abd (jdiag, k)
         lm = min (k, jdiag) - 1
         la = jdiag - lm - 1
         lb = k - lm - 1
         t = - b(k)
         do 1030 j = 1, lm
            b(lb + j) = t * abd(la + j, k) + b(lb + j)
1030     continue
1040  continue

      return
      end
      subroutine twgrab (error, last, first, number)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWGRAB
C
C     RESERVE SPACE IN AN ARRAY.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

      integer
     +   first, last, number
      intrinsic
     +   max
      logical
     +   error

C///  CHECK THE ARGUMENTS.

      error = .not. (0 .le. last)
      if (error) go to 99999

      error = .not. (0 .le. number)
      if (error) go to 99999

C///  GRAB THE SPACE.

      first = last + 1
      last = last + max (1, number)

C///  EXIT.

99999 continue
      return
      end
      subroutine twinit (error, text, force)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWINIT
C
C     INITIALIZE THE CONTROLS.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character
     +   id*9
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   rvalue
      integer
     +   cntrls, count, ivalue, text
      logical
     +   error, first, force, lvalue, mess

      parameter (id = 'TWINIT:  ')
      parameter (cntrls = 22)

      dimension ivalue(cntrls), lvalue(cntrls), rvalue(cntrls)

      common / twcomi / ivalue
      common / twcoml / lvalue
      common / twcomr / rvalue

C     THE GNU F77 COMPILER REQUIRES THE SAVE TO PRECEED THE DATA

      save first

      data first / .true. /

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

      if (mess .and. 0 .lt. text) go to 9001

C///  TOP OF THE BLOCK TO SET THE CONTROLS.

      if (first .or. force) then
         first = .false.

C///  SET THE CONTROLS.

      count = 0

C     ADAPT

      count = count + 1
      lvalue(count) = .false.

C     LEVELD

      count = count + 1
      ivalue(count) = 1

C     LEVELM

      count = count + 1
      ivalue(count) = 1

C     PADD

      count = count + 1
      lvalue(count) = .false.

C     SSABS

      count = count + 1
      rvalue(count) = 1.0E-9

C     SSAGE

      count = count + 1
      ivalue(count) = 10

C     SSREL

      count = count + 1
      rvalue(count) = 1.0E-6

C     STEADY

      count = count + 1
      lvalue(count) = .true.

C     STEPS0

      count = count + 1
      ivalue(count) = 0

C     STEPS1

      count = count + 1
      ivalue(count) = 200

C     STEPS2

      count = count + 1
      ivalue(count) = 100

C     STRID0

      count = count + 1
      rvalue(count) = 1.0E-4

C     TDABS

      count = count + 1
      rvalue(count) = 1.0E-9

C     TDAGE

      count = count + 1
      ivalue(count) = 20

C     TDEC

      count = count + 1
      rvalue(count) = 3.1623

C     TDREL

      count = count + 1
      rvalue(count) = 1.0E-6

C     TINC

      count = count + 1
      rvalue(count) = 10.0

C     TMAX

      count = count + 1
      rvalue(count) = 1.0E-2

C     TMIN

      count = count + 1
      rvalue(count) = 1.0E-20

C     TOLER0

      count = count + 1
      rvalue(count) = 1.0E-9

C     TOLER1

      count = count + 1
      rvalue(count) = 0.2

C     TOLER2

      count = count + 1
      rvalue(count) = 0.2

C///  BOTTOM OF THE BLOCK TO SET THE CONTROLS.

         error = .not. (count .eq. cntrls)
         if (error) go to 9001
      end if

C///  ERROR MESSAGES.

      go to 99999

9001  if (0 .lt. text) write (text, 99001) id, cntrls, count
      if (.not. mess) go to 99999

99001 format
     +  (/1X, a9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.'
     + //10X, i10, '  CONTROLS'
     +  /10X, i10, '  COUNTED')

C///  EXIT.

      stop
99999 continue
      return
      end
      subroutine twlaps (timer)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWLAPS
C
C     OBTAIN ELAPSED COMPUTING TIME IN SECONDS.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      external twtime
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   temp, timer

      call twtime (temp)
      timer = temp - timer

      return
      end
      subroutine twlast (length, string)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWLAST
C
C     FIND THE LAST NONBLANK CHARACTER IN A STRING.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character
     +   string*(*)
      integer
     +   j, length
      intrinsic
     +   len

      do 0100 j = len (string), 1, - 1
         if (string(j : j) .ne. ' ') then
            length = j
            go to 0200
         end if
0100  continue
      length = 1
0200  continue

      return
      end
      subroutine twlogr (string, value)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWLOGR
C
C     WRITE A COMMON LOGARITHM TO A CHARACTER STRING.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

      character
     +   string*(*)
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   value
      intrinsic
     +   len, log10

      if (6 .le. len (string)) then
         if (value .lt. 0.0) then
            string = ' '
         else if (value .eq. 0.0) then
            string = '  ZERO'
         else
            write (string, '(F6.2)') log10 (value)
         end if
      else
         string = '******'
      end if

      return
      end
      subroutine twnorm (n, value, x)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWNORM
C
C     COMPUTE THE MAX-NORM OF A VECTOR.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   value, x
      integer
     +   j, n
      intrinsic
     +   abs, max

      dimension x(n)

      value = 0.0
      do 0100 j = 1, n
         value = max (value, abs (x(j)))
0100  continue

      return
      end
      subroutine twopnt
     +  (error, text, versio,
     +   above, active, below, buffer, comps, condit, groupa, groupb,
     +   isize, iwork, mark, name, names, pmax, points, report, rsize,
     +   rwork, signal, stride, time, u, x)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWOPNT
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character
     +   column*80, ctemp1*80, ctemp2*80, header*80, id*9, name*(*),
     +   report*(*), signal*(*), string*80, versio*(*), vnmbr*8
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   above, below, buffer, condit, detail, maxcon, ratio, rvalue,
     +   rwork, ssabs, ssrel, strid0, stride, tdabs, tdec, tdrel, temp,
     +   timer, tinc, tmax, tmin, toler0, toler1, toler2, total, u, x,
     +   ynorm
      external
     +   evolve, refine, search, twcopy, twgrab, twlaps, twlast, twlogr,
     +   twnorm, twsqez, twtime, twinit
      integer
     +   age, cntrls, comps, count, desire, event, gmax, grid, groupa,
     +   groupb, ilast, isize, ivalue, iwork, j, jacobs, k, label, len1,
     +   len2, length, leveld, levelm, lines, names, nsteps, padd, pmax,
     +   points, psave, qabove, qbelow, qbnds, qdvrg, qentry, qexit,
     +   qfunct, qgrid, qjacob, qnull, qother, qrat1, qrat2, qrefin,
     +   qs0, qs1, qsearc, qsolve, qtask, qtimst, qtotal, qtype, qusave,
     +   qv1, qvary, qvary1, qvary2, qvsave, qxsave, qy0, qy1, return,
     +   rlast, route, rsize, size, ssage, step, steps, steps0, steps1,
     +   steps2, tdage, text, vnmbrs, xrepor
      intrinsic
     +   max
      logical
     +   active, adapt, allow, error, exist, first, flag, found, lvalue,
     +   mark, mess, satisf, steady, time

      parameter (id = 'TWOPNT:  ')
      parameter (cntrls = 22)
      parameter (gmax = 100)
      parameter (lines = 20)
      parameter (vnmbrs = 12)

C     REPORT CODES
      parameter (qnull = 0, qbnds = 1, qdvrg = 2)

C     LOCATION OF DATA IN ARRAYS DETAIL, EVENT, TIMER, AND TOTAL.  THE
C     LOCATIONS ARE CHOSEN TO SIMPLIFY WRITE STATEMENTS.  DETAIL USES
C     ONLY 1 : 8, EVENT USES ONLY 5 : 8, TIMER USES 1 : 9, AND TOTAL
C     USES ONLY 2 : 9.  IN ADDITION, 2, 3, 4, 10, AND 11 ARE USED AS
C     MNEMONIC VALUES FOR QTASK.
      parameter
     +  (qgrid  =  1,
     +   qtimst =  2,
     +   qsearc =  3,
     +   qrefin =  4,
     +   qfunct =  5,
     +   qjacob =  6,
     +   qsolve =  7,
     +   qother =  8,
     +   qtotal =  9,
     +   qentry = 10,
     +   qexit  = 11)

      dimension
     +   above(groupa + comps + groupb), active(*), below(groupa + comps
     +   + groupb), buffer(groupa + comps * pmax + groupb), column(3),
     +   detail(gmax, qtotal), event(gmax, qtotal), header(6),
     +   ivalue(cntrls), iwork(isize), lvalue(cntrls), mark(*),
     +   name(names), ratio(2), rvalue(cntrls), rwork(rsize),
     +   size(gmax), timer(qtotal), total(qtotal), u(groupa + comps *
     +   pmax + groupb), vnmbr(vnmbrs), x(*)

      common / twcomi / ivalue
      common / twcoml / lvalue
      common / twcomr / rvalue

C///  SAVE LOCAL VALUES DURING RETURNS FOR REVERSE COMMUNCIATION.

      save

C///////////////////////////////////////////////////////////////////////
C
C     PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  EVERY-TIME INITIALIZATION.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

C     TURN OFF ALL REVERSE COMMUNICATION FLAGS.
      time = .false.

C///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      if (signal .ne. ' ') then
         go to (9912, 9922, 9932, 9942) route
         error = .true.
         go to 9001
      end if

C///////////////////////////////////////////////////////////////////////
C
C     ENTRY BLOCK.  INITIALIZE A NEW PROBLEM.
C
C///////////////////////////////////////////////////////////////////////

C///  TURN OFF ALL STATUS REPORTS.

      error = .false.
      report = ' '

C///  WRITE ALL MESSAGES.

      if (mess .and. 0 .lt. text) then
         label = 0
         return = 0
         route = 0

         write (text, 10004) id, '???'
         write (text, 10020) id
         write (text, 10017) id
         write (text, 10014) id
         string = vnmbr(vnmbrs)
         call twlast (length, string)
         write (text, 10001) id, 'DOUBLE PRECISION', string (1 : length)
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

C///  CHECK THE VERSION.

      data vnmbr
     +   / '3.18', '3.19', '3.20', '3.21', '3.22', '3.23', '3.24',
     +     '3.25', '3.26', '3.27', '3.28', '3.29' /

      flag = .false.
      do 1010 j = 1, vnmbrs
         flag = flag .or.
C*****PRECISION > DOUBLE
     +      versio .eq. 'DOUBLE PRECISION VERSION ' // vnmbr(j)
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C     +      VERSIO .EQ. 'SINGLE PRECISION VERSION ' // VNMBR(J)
C*****END PRECISION > SINGLE
1010  continue
      error = .not. flag
      if (error) go to 9002

C///  SET THE CONTROLS.

C     SUBROUTINE TWINIT (ERROR, TEXT, FORCE)

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

      error = .not. (count .eq. cntrls)
      if (error) go to 9004

C///  PRINT THE ENTRY BANNER AT ALL PRINT LEVELS.

      string = vnmbr(vnmbrs)
      call twlast (length, string)
      if ((0 .lt. levelm .or. mess) .and. 0 .lt. text)
     +   write (text, 10001) id,
C*****PRECISION > DOUBLE
     +   'DOUBLE PRECISION', string (1 : length)
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C     +   'SINGLE PRECISION', STRING (1 : LENGTH)
C*****END PRECISION > SINGLE

C///  CHECK THE ARGUMENTS.

      error = .not. (leveld .le. levelm)
      if (error) go to 9005

      error = .not. (0 .le. comps .and. 0 .le. points .and.
     +   0 .le. groupa .and. 0 .le. groupb)
      if (error) go to 9006

      error = .not. ((0 .lt. comps) .eqv. (0 .lt. points))
      if (error) go to 9007

      error = .not. (0 .lt. groupa + comps * points + groupb)
      if (error) go to 9008

      error = .not. (names .eq. 1 .or.
     +   names .eq. groupa + comps + groupb)
      if (error) go to 9009

      error = .not. (points .le. pmax)
      if (error) go to 9010

      count = 0
      do 1020 j = 1, groupa + comps + groupb
         if (.not. (below(j) .lt. above(j))) count = count + 1
1020  continue
      error = count .ne. 0
      if (error) go to 9011

C///  PARTITION THE INTEGER WORK SPACE.

C     SUBROUTINE TWGRAB (ERROR, LAST, FIRST, NUMBER)

      ilast = 0

C     VARY(PMAX)
      call twgrab (error, ilast, qvary, pmax)
      if (error) go to 9012

C     VARY1(PMAX)
      call twgrab (error, ilast, qvary1, pmax)
      if (error) go to 9012

C     VARY2(PMAX)
      call twgrab (error, ilast, qvary2, pmax)
      if (error) go to 9012

C///  PARTITION THE REAL WORK SPACE.

C     SUBROUTINE TWGRAB (ERROR, LAST, FIRST, NUMBER)

      rlast = 0

C     ABOVE(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qabove, groupa + comps * pmax + groupb)
      if (error) go to 9012

C     BELOW(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qbelow, groupa + comps * pmax + groupb)
      if (error) go to 9012

C     RATIO1(PMAX)
      call twgrab (error, rlast, qrat1, pmax)
      if (error) go to 9012

C     RATIO2(PMAX)
      call twgrab (error, rlast, qrat2, pmax)
      if (error) go to 9012

C     S0(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qs0, groupa + comps * pmax + groupb)
      if (error) go to 9012

C     S1(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qs1, groupa + comps * pmax + groupb)
      if (error) go to 9012

C     USAVE(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qusave, groupa + comps * pmax + groupb)
      if (error) go to 9012

C     VSAVE(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qvsave, groupa + comps * pmax + groupb)
      if (error) go to 9012

C     V1(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qv1, groupa + comps * pmax + groupb)
      if (error) go to 9012

C     XSAVE(PMAX)
      call twgrab (error, rlast, qxsave, pmax)
      if (error) go to 9012

C     Y0(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qy0, groupa + comps * pmax + groupb)
      if (error) go to 9012

C     Y1(GROUPA + COMPS * PMAX + GROUPB)
      call twgrab (error, rlast, qy1, groupa + comps * pmax + groupb)
      if (error) go to 9012

C///  CHECK THE WORK SPACES' SIZES.

      error = .not. (ilast .le. isize .and. rlast .le. rsize)
      if (error) go to 9013

C///  ONE-TIME INITIALIZATION.

C     ALLOW FURTHER TIME EVOLUTION
      allow = .true.

C     PRESENT TASK
      qtask = qentry

C     STATISTICS ARRAYS
      do 1040 k = 1, qtotal
         total(k) = 0.0
         do 1030 j = 1, gmax
            detail(j, k) = 0.0
            event(j, k) = 0
1030     continue
1040  continue

C     TOTAL TIME STATISTIC
      call twtime (timer(qtotal))

C     GRID POINTER AND STATISTICS FOR THE FIRST GRID
      grid = 1
      size(grid) = points
      call twtime (timer(qgrid))

C     TIME STEP NUMBER
      step = 0

C     SOLUTION FLAG
      found = .true.

C///  EXPAND THE BOUNDS.

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

C///  SAVE THE INITIAL SOLUTION.

      psave = points
      if (adapt .and. 0 .lt. points)
     +   call twcopy (points, x, rwork(qxsave))
      call twcopy (groupa + comps * points + groupb, u, rwork(qusave))

C     GO TO 1090 WHEN RETURN = 1
      return = 1
      go to 9911
1090  continue

C///  PRINT LEVELS 11, 21, AND 22.

      if (0 .lt. leveld .and. 0 .lt. text) then
         write (text, 10002) id, 'INITIAL GUESS:'
C        GO TO 1100 WHEN RETURN = 2
         return = 2
         go to 9921
      end if
1100  continue

C///  PRINT LEVEL 10 AND 11.

C                  123456789_123456789_123456789_1234
C                  12345678   123456  123456   123456
      header(1) = '            LOG10   LOG10         '
      header(2) = '    TASK   NORM F  COND J   REMARK'

      if (levelm .eq. 1 .and. 0 .lt. text) then
         if (0 .lt. leveld) write (text, 10002) id,
     +      'SOLVE THE PROBLEM.'
         write (text, 10003) (header(j), j = 1, 2)
C        GO TO 1110 WHEN LABEL = 1
         label = 1
         go to 7010
      end if
1110  continue

C///////////////////////////////////////////////////////////////////////
C
C     DECISION BLOCK.  THE PREVIOUS TASK DETERMINES THE NEXT.
C
C///////////////////////////////////////////////////////////////////////

2010  continue

C///  ENTRY WAS THE PREVIOUS TASK.

      if (qtask .eq. qentry) then
         if (0 .lt. steps0) then
            qtask = qtimst
            desire = steps0
         else if (steady) then
            qtask = qsearc
         else
            error = .true.
            go to 9014
         end if

C///  SEARCH WAS THE PREVIOUS TASK.

      else if (qtask .eq. qsearc) then
         if (found) then
            if (adapt) then
               qtask = qrefin
            else
               qtask = qexit
               report = ' '
            end if
         else
            if (allow .and. 0 .lt. steps1) then
               qtask = qtimst
               desire = steps1
            else
               qtask = qexit
               if (1 .lt. grid) then
                  report = 'SOME SOLVED'
               else
                  report = 'NOT SOLVED'
               end if
            end if
         end if

C///  REFINE WAS THE PREVIOUS TASK.

      else if (qtask .eq. qrefin) then
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

C///  EVOLVE WAS THE PREVIOUS TASK.

      else if (qtask .eq. qtimst) then
         if (found) then
            if (steady) then
               qtask = qsearc
            else
               qtask = qexit
               report = ' '
            end if
         else
            qtask = qexit
            if (1 .lt. grid) then
               report = 'SOME SOLVED'
            else
               report = 'NOT SOLVED'
            end if
         end if
      end if

C///  BRANCH TO THE NEXT TASK.

      if (qtask .eq. qexit) go to 3010
      if (qtask .eq. qsearc) go to 4010
      if (qtask .eq. qrefin) go to 5010
      if (qtask .eq. qtimst) go to 6010
      error = .true.
      go to 9015

C///////////////////////////////////////////////////////////////////////
C
C     EXIT BLOCK.
C
C///////////////////////////////////////////////////////////////////////

3010  continue

C///  COMPLETE STATISTICS FOR THE LAST GRID.

      call twlaps (timer(qgrid))
      if (grid .le. gmax) then
         detail(grid, qgrid) = timer(qgrid)
         detail(grid, qother)
     +      = detail(grid, qgrid) - (detail(grid, qfunct)
     +      + detail(grid, qjacob) + detail(grid, qsolve))
      end if

C///  RESTORE THE SOLUTION.

      if (report .ne. ' ') then
C        BE CAREFUL NOT TO ASSIGN A VALUE TO A PARAMETER
         if (points .ne. psave) points = psave
         if (adapt .and. 0 .lt. points)
     +      call twcopy (points, rwork(qxsave), x)
         call twcopy
     +      (groupa + comps * points + groupb, rwork(qusave), u)
      end if

C///  PRINT LEVEL 11 OR 21.

C     SAVE THE STATUS REPORTS DURING REVERSE COMMUNICATION
      string = report

      if (leveld .eq. 1 .and. 0 .lt. text) then
         write (text, 10002) id, 'FINAL SOLUTION:'
C        GO TO 3020 WHEN RETURN = 3
         return = 3
         go to 9921
      end if
3020  continue

C     RESTORE THE STATUS REPORTS AFTER REVERSE COMMUNICATION
      report = string

C///  COMPLETE THE TOTAL TIME STATISTICS.

      call twlaps (timer(qtotal))
      total(qtotal) = timer(qtotal)
      total(qother) = total(qtotal)
     +   - (total(qfunct) + total(qjacob) + total(qsolve))

C///  TOP OF THE REPORT BLOCK.

      if (0 .lt. levelm .and. 0 .lt. text) then
         if (0.0 .lt. total(qtotal)) then

C///  REPORT TOTAL COMPUTER TIME.

      temp = total(qtotal)
      if (3600.0 .le. temp) then
         write (string, '(F10.2, A)') temp / 3600.0, ' HOURS'
      else if (60.0 .le. temp) then
         write (string, '(F10.2, A)') temp / 60.0, ' MINUTES'
      else
         write (string, '(F10.2, A)') temp, ' SECONDS'
      end if

      call twsqez (length, string)
      write (text, 10004) id, string (1 : length)

C///  REPORT PERCENT OF TOTAL COMPUTER TIME.

      temp = 100.0 / total(qtotal)
      if (adapt) then

C                  123456789_123456789_123456789_12345678
C                  123456  123456  123456 123456 123456
      header(1) = '                TASK                  '
      header(3) = '  GRID    GRID  --------------------  '
      header(5) = 'POINTS  TOTALS  EVOLVE SEARCH REFINE  '

C                  123456789_123456789_1234567
C                  123456 123456 123456 123456
      header(2) = 'SUBTASK                    '
      header(4) = '---------------------------'
      header(6) = 'EVAL F PREP J  SOLVE  OTHER'

      write (text, 10005) header,
     +   (size(j), (temp * detail(j, k), k = 1, 8), j = 1, grid)
      if (1 .lt. grid) write (text, 10006) (temp * total(k), k = 2, 8)
      if (gmax .lt. grid) write (text, 10007)

      else

C                  123456789_123456789_123456789_123456789_123456789_1
C                  123456   123456   123456   123456   123456   123456
      header(1) = 'SUBTASK                             TASK           '
      header(2) = '---------------------------------   ---------------'
      header(3) = 'EVAL F   PREP J    SOLVE    OTHER   EVOLVE   SEARCH'

      write (text, 10008)
     +   (header(j), j = 1, 3), '  % OF TOTAL',
     +   (temp * total(k), k = 5, 8), (temp * total(k), k = 2, 3),
     +   'MEAN SECONDS', (detail(1, k) / event(1, k), k = 5, 7),
     +   '    QUANTITY', (event(1, k), k = 5, 7)

      end if

C///  REPORT AVERAGE COMPUTER TIME.

C                  123456789_123456789_123456789_1234567
C                  123456   1234567  1234567  1234567
      header(1) = '         AVERAGE SECONDS             '
      header(3) = '  GRID   -------------------------   '
      header(5) = 'POINTS    EVAL F   PREP J    SOLVE   '


C                  123456789_123456789_12345
C                  1234567  1234567  1234567
      header(2) = 'NUMBER OF SUBTASKS       '
      header(4) = '-------------------------'
      header(6) = ' EVAL F   PREP J    SOLVE'

      if (adapt) write (text, 10009) header,
     +   (size(j), (detail(j, k) / event(j, k), k = 5, 7),
     +   (event(j, k), k = 5, 7), j = 1, grid)

      end if

C///  REPORT THE COMPLETION STATUS.

      if (0 .lt. levelm) then
         if (report .eq. ' ') then
            write (text, 10010) id
         else if (report .eq. 'NO SPACE') then
            write (string, '(I10)') points
            call twsqez (length, string)
            write (text, 10011)
     +         id, string (1 : length), ratio, toler1, toler2
         else if (report .eq. 'NOT SOLVED') then
            write (text, 10012) id
         else if (report .eq. 'SOME SOLVED') then
            write (string, '(I10)') points
            call twsqez (length, string)
            write (text, 10013)
     +         id, string (1 : length), ratio, toler1, toler2
         else
            error = .true.
            go to 9016
         end if
      end if

C///  BOTTOM OF THE REPORT BLOCK.

      end if

C///  BOTTOM OF THE EXIT BLOCK.

      go to 99999

C///////////////////////////////////////////////////////////////////////
C
C     SEARCH BLOCK.
C
C///////////////////////////////////////////////////////////////////////

4010  continue

C///  INITIALIZE STATISTICS ON ENTRY TO THE SEARCH BLOCK.

      call twtime (timer(qsearc))
      first = .true.
      jacobs = 0
      maxcon = 0.0

C///  PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE SEARCH BLOCK.

      if (1 .lt. levelm) then
         if (0 .lt. text) write (text, 10014) id
      end if

C///  PREPARE TO CALL SEARCH.

C     SAVE THE SOLUTION SHOULD THE SEARCH FAIL
      call twcopy (groupa + comps * points + groupb, u, rwork(qvsave))

      exist = .false.

C///  CALL SEARCH.

      age = 0
4020  continue

C     SUBROUTINE SEARCH
C    +  (ERROR, TEXT,
C    +   ABOVE, AGE, BELOW, BUFFER, COMPS, CONDIT, EXIST, GROUPA,
C    +   GROUPB, LEVELD, LEVELM, NAME, NAMES, POINTS, REPORT, S0, S1,
C    +   SIGNAL, STEPS, SUCCES, V0, V1, XXABS, XXAGE, XXREL, Y0, Y0NORM,
C    +   Y1)

      call search
     +  (error, text,
     +   rwork(qabove), age, rwork(qbelow), buffer, comps, condit,
     +   exist, groupa, groupb, leveld - 1, levelm - 1, name, names,
     +   points, xrepor, rwork(qs0), rwork(qs1), signal, nsteps, found,
     +   u, rwork(qv1), ssabs, ssage, ssrel, rwork(qy0), ynorm,
     +   rwork(qy1))
      if (error) go to 9017

C///  PASS REQUESTS FROM SEARCH TO THE CALLER.

      if (signal .ne. ' ') then
C        GO TO 4020 WHEN RETURN = 4
         return = 4
         go to 9931
      end if

C///  REACT TO THE COMPLETION OF SEARCH.

      if (found) then
C        SAVE THE LATEST SOLUTION

         psave = points
         if (adapt .and. 0 .lt. points)
     +      call twcopy (points, x, rwork(qxsave))
         call twcopy
     +      (groupa + comps * points + groupb, u, rwork(qusave))

C        GO TO 4030 WHEN RETURN = 5
         return = 5
         go to 9911
      else
C        RESTORE THE SOLUTION
         call twcopy
     +      (groupa + comps * points + groupb, rwork(qvsave), u)
      end if
4030  continue

C///  COMPLETE STATISTICS FOR THE SEARCH BLOCK.

      call twlaps (timer(qsearc))
      total(qsearc) = total(qsearc) + timer(qsearc)
      if (grid .le. gmax)
     +   detail(grid, qsearc) = detail(grid, qsearc) + timer(qsearc)

C///  PRINT LEVEL 10 OR 11 ON EXIT FROM THE SEARCH BLOCK.

      if (levelm .eq. 1 .and. 0 .lt. text) then
C        GO TO 4040 WHEN LABEL = 2
         label = 2
         go to 7010
      end if
4040  continue

C///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE SEARCH BLOCK.

      if (1 .lt. levelm) then
         if (found) then
            if (0 .lt. text) write (text, 10015) id
         else
            if (0 .lt. text) write (text, 10016) id
         end if
      end if

C///  BOTTOM OF THE SEARCH BLOCK.

      go to 2010

C///////////////////////////////////////////////////////////////////////
C
C     REFINE BLOCK.
C
C///////////////////////////////////////////////////////////////////////

5010  continue

C///  INITIALIZE STATISTICS ON ENTRY TO THE REFINE BLOCK.

      call twtime (timer(qrefin))

C///  PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE REFINE BLOCK.

      if (1 .lt. levelm) then
         if (0 .lt. text) write (text, 10017) id
      end if

C///  PREPARE TO CALL REFINE.

C     SAVE THE GROUP B VALUES
      do 5020 j = 1, groupb
         rwork(qvsave - 1 + j) = u(groupa + comps * points + j)
5020  continue

      exist = .false.

C///  CALL REFINE.

5030  continue

C     SUBROUTINE REFINE
C    +  (ERROR, TEXT,
C    +   ACTIVE, BUFFER, COMPS, LEVELD, LEVELM, MARK, NEWX, PADD, PMAX,
C    +   POINTS, RATIO, RATIO1, RATIO2, SIGNAL, SUCCES, TOLER0, TOLER1,
C    +   TOLER2, U, VARY1, VARY2, WEIGHT, X)

      call refine
     +  (error, text,
     +   active,
     +   buffer(groupa + 1), comps, leveld - 1, levelm - 1, mark,
     +   found, padd, pmax, points, ratio, rwork(qrat1), rwork(qrat2),
     +   signal, satisf, toler0, toler1, toler2, u(groupa + 1),
     +   iwork(qvary1), iwork(qvary2), iwork(qvary), x)
      if (error) go to 9018

C///  SERVICE REQUESTS FROM REFINE: PASS REQUESTS TO THE CALLER.

      if (signal .ne. ' ') then
C        INSERT THE GROUP A AND B UNKNOWNS
         do 5040 j = 1, groupa
            buffer(j) = u(j)
5040     continue
         do 5050 j = 1, groupb
            buffer(groupa + comps * points + j) = rwork(qvsave - 1 + j)
5050     continue

C        GO TO 5030 WHEN RETURN = 6
         return = 6
         go to 9931
      end if

C///  REACT TO THE COMPLETION OF REFINE.

      if (.not. found) go to 5110

C        COMPLETE STATISTICS FOR THE OLD GRID
         call twlaps (timer(qgrid))
         if (grid .le. gmax) then
            detail(grid, qgrid) = timer(qgrid)
            detail(grid, qother)
     +         = detail(grid, qgrid) - (detail(grid, qfunct)
     +         + detail(grid, qjacob) + detail(grid, qsolve))
         end if

C        INITIALIZE STATISTICS FOR THE NEW GRID
         grid = grid + 1
         if (grid .le. gmax) then
            call twtime (timer(qgrid))
            size(grid) = points
         end if

C        INSERT THE GROUP B VALUES
         do 5060 j = 1, groupb
            u(groupa + comps * points + j) = rwork(qvsave - 1 + j)
5060     continue

C        EXPAND THE BOUNDS
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

C        SAVE THE LATEST SOLUTION
C        GO TO 5100 WHEN RETURN = 7
         return = 7
         go to 9911
5100     continue

5110  continue

C///  COMPLETE STATISTICS FOR THE REFINE BLOCK.

      call twlaps (timer(qrefin))
      total(qrefin) = total(qrefin) + timer(qrefin)
      if (grid .le. gmax)
     +   detail(grid, qrefin) = detail(grid, qrefin) + timer(qrefin)

C///  PRINT LEVEL 10 OR 11 ON EXIT FROM THE REFINE BLOCK.

      if (levelm .eq. 1 .and. 0 .lt. text) then
         write (text, '()')
C        GO TO 5120 WHEN LABEL = 3
         label = 3
         go to 7010
      end if
5120  continue

C///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE REFINE BLOCK.

      if (1 .lt. levelm) then
         if (found) then
            if (0 .lt. text) write (text, 10018) id
         else
            if (0 .lt. text) write (text, 10019) id
         end if
      end if

C///  BOTTOM OF THE REFINE BLOCK.

      go to 2010

C///////////////////////////////////////////////////////////////////////
C
C     EVOLVE BLOCK.
C
C///////////////////////////////////////////////////////////////////////

6010  continue

C///  INITIALIZE STATISTICS ON ENTRY TO THE EVOLVE BLOCK.

      call twtime (timer(qtimst))
      first = .true.
      jacobs = 0
      maxcon = 0.0
      steps = step

C///  PRINT LEVEL 20, 21, OR 22 ON ENTRY TO THE EVOLVE BLOCK.

      if (1 .lt. levelm) then
         if (0 .lt. text) write (text, 10020) id
      end if

C///  CALL EVOLVE.

6020  continue

C     SUBROUTINE EVOLVE
C    +  (ERROR, TEXT,
C    +   ABOVE, BELOW, BUFFER, COMPS, CONDIT, DESIRE, GROUPA, GROUPB,
C    +   LEVELD, LEVELM, NAME, NAMES, POINTS, REPORT, S0, S1, SIGNAL,
C    +   STEP, STEPS2, STRID0, STRIDE, SUCCES, TDABS, TDAGE, TDEC,
C    +   TDREL, TIME, TINC, TMAX, TMIN, V0, V1, VSAVE, Y0, Y1, YNORM)

      call evolve
     +  (error, text,
     +   rwork(qabove), rwork(qbelow), buffer, comps, condit, desire,
     +   groupa, groupb, leveld - 1, levelm - 1, name, names, points,
     +   xrepor, rwork(qs0), rwork(qs1), signal, step, steps2, strid0,
     +   stride, found, tdabs, tdage, tdec, tdrel, time, tinc, tmax,
     +   tmin, u, rwork(qv1), rwork(qvsave), rwork(qy0), rwork(qy1),
     +   ynorm)
      if (error) go to 9019

C///  PASS REQUESTS FROM EVOLVE TO THE CALLER.

      if (signal .ne. ' ') then
C        GO TO 6020 WHEN RETURN = 8
         return = 8
         go to 9931
      end if

C///  REACT TO THE COMPLETION OF EVOLVE.

      if (found) then
C        SAVE THE LATEST SOLUTION
C        GO TO 6030 WHEN RETURN = 9
         return = 9
         go to 9911
      end if
6030  continue

C///  ALLOW FURTHER TIME EVOLUTION.

      allow = xrepor .eq. qnull

C///  COMPLETE STATISTICS FOR THE EVOLVE BLOCK.

      call twlaps (timer(qtimst))
      total(qtimst) = total(qtimst) + timer(qtimst)
      if (grid .le. gmax)
     +   detail(grid, qtimst) = detail(grid, qtimst) + timer(qtimst)
      steps = step - steps

C///  PRINT LEVEL 10 OR 11 ON EXIT FROM THE EVOLVE BLOCK.

      if (levelm .eq. 1 .and. 0 .lt. text) then
C        GO TO 6040 WHEN LABEL = 4
         label = 4
         go to 7010
      end if
6040  continue

C///  PRINT LEVEL 20, 21, OR 22 ON EXIT FROM THE EVOLVE BLOCK.

      if (1 .lt. levelm) then
         if (found) then
            if (0 .lt. text) write (text, 10021) id
         else
            if (0 .lt. text) write (text, 10022) id
         end if
      end if

C///  BOTTOM OF THE EVOLVE BLOCK.

      go to 2010

C///////////////////////////////////////////////////////////////////////
C
C     BLOCK TO PRINT LOG LINES.
C
C///////////////////////////////////////////////////////////////////////

7010  continue

      do 7020 j = 1, 3
         column(j) = ' '
7020  continue

      string = ' '

C     COLUMN 1: NAME OF THE TASK
      if (qtask .eq. qentry) column(1) = '   START'
      if (qtask .eq. qsearc) column(1) = '  SEARCH'
      if (qtask .eq. qrefin) column(1) = '  REFINE'
      if (qtask .eq. qtimst) column(1) = '  EVOLVE'

C     COLUMN 2: NORM OF THE STEADY STATE FUNCTION
      if (.not. found) go to 7040
C        GO TO 7030 WHEN RETURN = 10
         return = 10
         go to 9941
7030     continue
         call twnorm (groupa + comps * points + groupb, temp, buffer)
         call twlogr (column(2), temp)
7040  continue

C     COLUMN 3: LARGEST CONDITION NUMBER
      if (qtask .eq. qsearc .or. qtask .eq. qtimst) then
         if (maxcon .ne. 0.0) call twlogr (column(3), maxcon)
      end if

C     REMARK
      if (qtask .eq. qsearc) then
         if (xrepor .eq. qdvrg) then
            string = 'DIVERGING'
         else if (xrepor .eq. qnull) then
            if (nsteps .eq. 1) then
               write (string, '(I10, A)') nsteps, ' SEARCH STEP'
            else
               write (string, '(I10, A)') nsteps, ' SEARCH STEPS'
            end if
         else if (xrepor .eq. qbnds) then
            string = 'GOING OUT OF BOUNDS'
         else
            string = '?'
         end if
      else if (qtask .eq. qtimst) then
         if (xrepor .eq. qbnds .or. xrepor .eq. qdvrg .or.
     +      xrepor .eq. qnull) then
            write (string, '(I10, A, 1P, E10.1, A)')
     +         steps, ' TIME STEPS, ', stride, ' LAST STRIDE'
         else
            string = '?'
         end if
      else if (qtask .eq. qentry .and. adapt) then
         write (string, '(I10, A)') points, ' GRID POINTS'
      else if (qtask .eq. qrefin) then
         if (found) then
            write (string, '(F10.2, A, F10.2, A, I10, A)')
     +         ratio(1), ' AND ', ratio(2), ' RATIOS, ', points,
     +         ' GRID POINTS'
         else
            write (string, '(F10.2, A, F10.2, A)')
     +         ratio(1), ' AND ', ratio(2), ' RATIOS'
         end if
      end if

      call twsqez (length, string)
      if (0 .lt. text) write (text, 10023) column, string (1 : length)

      go to (1110, 4040, 5120, 6040) label
      error = .true.
      go to 9020

C///////////////////////////////////////////////////////////////////////
C
C     REQUEST REVERSE COMMUNICATION.
C
C///////////////////////////////////////////////////////////////////////

C///  SAVE THE SOLUTION.

9911  continue

      call twcopy (groupa + comps * points + groupb, u, buffer)
      signal = 'SAVE'
C     GO TO 9912 WHEN ROUTE = 1
      route = 1
      go to 99999
9912  continue
      signal = ' '

      go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030)
     +   return
      error = .true.
      go to 9021

C///  PRINT THE LATEST SOLUTION.

9921  continue

      call twcopy (groupa + comps * points + groupb, u, buffer)
      signal = 'SHOW'
C     GO TO 9922 WHEN ROUTE = 2
      route = 2
      go to 99999
9922  continue
      signal = ' '

      go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030)
     +   return
      error = .true.
      go to 9021

C///  PASS REQUESTS FROM SEARCH, REFINE, OR EVOLVE TO THE CALLER.

9931  continue

C     IDENTIFY THE REQUEST.  THIS MUST BE SAVED TO GATHER STATISTICS
C     AT REENTRY.  THE REVERSE COMMUNICATION FLAGS WILL NOT BE SAVED
C     BECAUSE THEY ARE CLEARED AT EVERY ENTRY.
      if (signal .eq. 'RESIDUAL') then
         qtype = qfunct
      else if (signal .eq. 'PREPARE') then
         qtype = qjacob
      else if (signal .eq. 'SOLVE') then
         qtype = qsolve
      else
         qtype = qother
      end if

C     COUNT THE JACOBIANS
      if (qtype .eq. qjacob) jacobs = jacobs + 1

      call twtime (timer(qtype))

C     GO TO 9932 WHEN ROUTE = 3
      route = 3
      go to 99999
9932  continue

C     SAVE THE CONDITION NUMBER
      if (qtype .eq. qjacob) maxcon = max (maxcon, condit)

      call twlaps (timer(qtype))
      total(qtype) = total(qtype) + timer(qtype)
      if (grid .le. gmax) then
         detail(grid, qtype) = detail(grid, qtype) + timer(qtype)
         event(grid, qtype) = event(grid, qtype) + 1
      end if

      go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030)
     +   return
      error = .true.
      go to 9021

C///  EVALUATE THE STEADY STATE FUNCTION.

9941  continue
      call twtime (timer(qfunct))

      call twcopy (groupa + comps * points + groupb, u, buffer)
      signal = 'RESIDUAL'
      time = .false.
C     GO TO 9942 WHEN ROUTE = 4
      route = 4
      go to 99999
9942  continue
      signal = ' '

      call twlaps (timer(qfunct))
      total(qfunct) = total(qfunct) + timer(qfunct)
      if (grid .le. gmax) then
         detail(grid, qfunct) = detail(grid, qfunct) + timer(qfunct)
         event(grid, qfunct) = event(grid, qfunct) + 1
      end if

      go to (1090, 1100, 3020, 4020, 4030, 5030, 5100, 6020, 6030, 7030)
     +   return
      error = .true.
      go to 9021

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 format
     +   (/1X, a9, a, ' (TWO POINT BOUNDARY VALUE PROBLEM) SOLVER,'
     +   /10X, 'VERSION ', a,
     +   ' OF APRIL 1998 BY DR. JOSEPH F. GRCAR.')

10002 format
     +   (/1X, a9, a)

10003 format
     +   (3(/10X, a35)/)

10004 format
     +  (/1X, a9, a, ' TOTAL COMPUTER TIME (SEE BREAKDOWN BELOW).')

10005 format
     +  (/10X, 'PERCENT OF TOTAL COMPUTER TIME FOR VARIOUS TASKS:'
     +   /3(/10X, a38, a27)
     +  //(10X, i6, 2X, f6.1, 1X, 3(1X, f6.1), 1X, 4(1X, f6.1)))

10006 format
     +  (/12X, 'TASK TOTALS:', 1X, 3(1X, f6.1), 1X, 4(1X, f6.1))

10007 format
     +  (/10X, 'SOME GRIDS ARE OMITTED, BUT THE TOTALS ARE FOR ALL.')

10008 format
     +  (3(/24X, a51)
     +  //10X, a12, f8.1, 5F9.1
     +   /10X, a12, f8.3, 2F9.3
     +   /10X, a12, i8, 2i9)

10009 format
     +  (/10X, 'AVERAGE COMPUTER TIMES FOR, AND NUMBERS OF, SUBTASKS:'
     +   /3(/10X, a37, a25)
     +  //(10X, i6, 3X, f7.3, 2X, f7.3, 2X, f7.3, 1X, 3(2X, i7)))

10010 format
     +  (/1X, a9, 'SUCCESS.  PROBLEM SOLVED.')

10011 format
     +  (/1X, a9, 'FAILURE.  A SOLUTION WAS FOUND FOR A GRID WITH ', a
     +  /10X, 'POINTS, BUT ONE OR BOTH RATIOS ARE TOO LARGE.'
C               123456789_  123456789_
     +  //22X, '   RATIO 1     RATIO 2'
     +  //10X, '     FOUND', 2F12.2
     +   /10X, '   DESIRED', 2F12.2
     +  //10X, 'A LARGER GRID COULD NOT BE FORMED.')

10012 format
     +  (/1X, a9, 'FAILURE.  NO SOLUTION WAS FOUND.')

10013 format
     +  (/1X, a9, 'FAILURE.  A SOLUTION WAS FOUND FOR A GRID WITH ', a
     +  /10X, 'POINTS, BUT ONE OR BOTH RATIOS ARE TOO LARGE.'
C               123456789_  123456789_
     +  //22X, '   RATIO 1     RATIO 2'
     +  //10X, '     FOUND', 2F12.2
     +   /10X, '   DESIRED', 2F12.2
     +  //10X, 'A SOLUTION COULD NOT BE FOUND FOR A LARGER GRID.')

10014 format
     +  (/1X, a9, 'CALLING SEARCH TO SOLVE THE STEADY STATE PROBLEM.')

10015 format
     +  (/1X, a9, 'SEARCH FOUND THE STEADY STATE.')

10016 format
     +  (/1X, a9, 'SEARCH DID NOT FIND THE STEADY STATE.')

10017 format
     +  (/1X, a9, 'CALLING REFINE TO PRODUCE A NEW GRID.')

10018 format
     +  (/1X, a9, 'REFINE SELECTED A NEW GRID.')

10019 format
     +  (/1X, a9, 'REFINE DID NOT SELECT A NEW GRID.')

10020 format
     +  (/1X, a9, 'CALLING EVOLVE TO PERFORM TIME EVOLUTION.')

10021 format
     +  (/1X, a9, 'EVOLVE PERFORMED A TIME EVOLUTION.')

10022 format
     +  (/1X, a9, 'EVOLVE DID NOT PERFORM A TIME EVOLUTION.')

10023 format
     +   (10X, a8, 3X, a6, 2X, a6, 3X, a)

80001 format
     +   ('(', a, ' ', i10, ')')

80002 format
     +   (10X, 1p, e10.2, 2X, e10.2, 3X, a)

80003 format
     +   (10X, '  ... MORE')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

C     GO TO 99999

9001  if (0 .lt. text) write (text, 99001) id, route
      if (.not. mess) go to 99999

9002  if (0 .lt. text) then
         call twlast (length, versio)
         write (text, 99002) id, versio (1 : length), vnmbr(vnmbrs)
         do 9901 j = vnmbrs - 1, 1, - 1
            write (text, '(10X, A, A)')
C*****PRECISION > DOUBLE
     +         ' CAN REPLACE:  DOUBLE PRECISION VERSION ', vnmbr(j)
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C     +         ' CAN REPLACE:  SINGLE PRECISION VERSION ', VNMBR(J)
C*****END PRECISION > SINGLE
9901     continue
      end if
      if (.not. mess) go to 99999

9003  if (0 .lt. text) write (text, 99003) id
      if (.not. mess) go to 99999

9004  if (0 .lt. text) write (text, 99004) id, cntrls, count
      if (.not. mess) go to 99999

9005  if (0 .lt. text) write (text, 99005) id, leveld, levelm
      if (.not. mess) go to 99999

9006  if (0 .lt. text) write (text, 99006) id,
     +   comps, points, groupa, groupb
      if (.not. mess) go to 99999

9007  if (0 .lt. text) write (text, 99007) id, comps, points
      if (.not. mess) go to 99999

9008  if (0 .lt. text) write (text, 99008) id,
     +   comps, points, groupa, groupb, groupa + comps * points + groupb
      if (.not. mess) go to 99999

9009  if (0 .lt. text) write (text, 99009) id,
     +   names, comps, groupa, groupb, groupa + comps + groupb
      if (.not. mess) go to 99999

9010  if (0 .lt. text) write (text, 99010) id, points, pmax
      if (.not. mess) go to 99999

9011  if (0 .lt. text) then
         write (text, 99011) id,
     +      groupa, groupb, comps, groupa + comps + groupb, count
         count = 0
         do 8010 j = 1, groupa + comps + groupb
            if (.not. (below(j) .lt. above(j)) .or. mess) then
               count = count + 1
               if (count .le. lines) then
                  if (names .eq. comps + groupa + groupb) then
                     ctemp1 = name(j)
                  else
                     ctemp1 = ' '
                  end if
                  call twsqez (len1, ctemp1)

                  if (j .le. groupa) then
                     write (ctemp2, 80001) 'A', j
                  else if (j .le. groupa + comps) then
                     write (ctemp2, 80001) 'C', j - groupa
                  else
                     write (ctemp2, 80001) 'B', j - groupa - comps
                  end if
                  call twsqez (len2, ctemp2)

                  if (ctemp1 .eq. ' ') then
                     string = ctemp2
                     length = len2
                  else if (len1 + 2 + len2 .le. 40) then
                     string = ctemp1 (1 : len1) // '  ' // ctemp2
                     length = len1 + 2 + len2
                  else if (len1 + 1 + len2 .le. 40) then
                     string = ctemp1 (1 : len1) // ' ' // ctemp2
                     length = len1 + 1 + len2
                  else
                     len1 = 40 - len2 - 4
                     string = ctemp1 (1 : len1) // '... ' // ctemp2
                     length = 40
                  end if

                  write (text, 80002)
     +               below(j), above(j), string (1 : length)
               end if
            end if
8010     continue
         if (lines .lt. count) write (text, 80003)
      end if
      if (.not. mess) go to 99999

9012  if (0 .lt. text) write (text, 99012) id
      if (.not. mess) go to 99999

9013  if (0 .lt. text) write (text, 99013) id,
     +   isize, rsize, ilast, rlast
      if (.not. mess) go to 99999

9014  if (0 .lt. text) write (text, 99014) id
      if (.not. mess) go to 99999

9015  if (0 .lt. text) write (text, 99015) id
      if (.not. mess) go to 99999

9016  if (0 .lt. text) write (text, 99016) id
      if (.not. mess) go to 99999

9017  if (0 .lt. text) write (text, 99017) id
      if (.not. mess) go to 99999

9018  if (0 .lt. text) write (text, 99018) id
      if (.not. mess) go to 99999

9019  if (0 .lt. text) write (text, 99019) id
      if (.not. mess) go to 99999

9020  if (0 .lt. text) write (text, 99020) id, label
      if (.not. mess) go to 99999

9021  if (0 .lt. text) write (text, 99021) id, return
      if (.not. mess) go to 99999

99001 format
     +  (/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, i10, '  ROUTE')

99002 format
     +  (/1X, a9, 'ERROR.  THE CALLING PROGRAM EXPECTS A VERSION OF'
     +  /10X, 'TWOPNT NOT COMPATIBLE WITH THIS VERSION.'
     + //10X, '     EXPECTS:  ', a
C*****PRECISION > DOUBLE
     + //10X, 'THIS VERSION:  DOUBLE PRECISION VERSION ', a)
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C     + //10X, 'THIS VERSION:  SINGLE PRECISION VERSION ', A)
C*****END PRECISION > SINGLE

99003 format
     +  (/1X, a9, 'ERROR.  TWINIT FAILS.')

99004 format
     +  (/1X, a9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.'
     + //10X, i10, '  CONTROLS'
     +  /10X, i10, '  COUNTED')

99005 format
     +  (/1X, a9, 'ERROR.  THE PRINTING LEVELS ARE OUT OF ORDER.'
     +  /10X, 'LEVELD CANNOT EXCEED LEVELM.'
     + //10X, i10, '  LEVELD, FOR SOLUTIONS'
     +  /10X, i10, '  LEVELM, FOR MESSAGES')

99006 format
     +  (/1X, a9, 'ERROR.  NUMBERS OF ALL TYPES OF UNKNOWNS MUST BE AT'
     +  /10X, 'LEAST ZERO.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS')

99007 format
     +  (/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS')

99008 format
     +  (/1X, a9, 'ERROR.  TOTAL UNKNOWNS MUST BE POSITIVE.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  TOTAL NUMBER')

99009 format
     +  (/1X, a9, 'ERROR.  THE NUMBER OF NAMES IS WRONG.'
     + //10X, i10, '  NAMES'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  TOTAL NUMBER')

99010 format
     +  (/1X, a9, 'ERROR.  THERE ARE TOO MANY POINTS.'
     + //10X, i10, '  POINTS'
     +  /10X, i10, '  PMAX, LIMIT ON POINTS')

99011 format
     +  (/1X, a9, 'ERROR.  THE LOWER AND UPPER BOUNDS ON SOME UNKNOWNS'
     +  /10X, 'ARE OUT OF ORDER.'
     + //10X, i10, '  GROUP A UNKNOWNS (A)'
     +  /10X, i10, '  GROUP B UNKNOWNS (B)'
     +  /10X, i10, '  COMPONENTS AT POINTS (C)'
     +  /10X, i10, '  TOTAL TYPES OF UNKNOWNS'
     +  /10X, i10, '  NUMBER OF BOUNDS OUT OF ORDER'
C              123456789_  123456789_
     + //10X, '     LOWER       UPPER'
     +  /10X, '     BOUND       BOUND   UNKNOWN'
     +  /)

99012 format
     +  (/1X, a9, 'ERROR.  TWGRAB FAILS.')

99013 format
     +  (/1X, a9, 'ERROR.  ONE OR BOTH WORK SPACES ARE TOO SMALL.'
C              123456789_  123456789_
     + //25X, '   INTEGER        REAL'
C              123456789_123
     + //10X, ' PRESENT SIZE', 2i12
     +  /10X, 'REQUIRED SIZE', 2i12)

99014 format
     +  (/1X, a9, 'ERROR.  NEITHER THE INITIAL TIME EVOLUTION NOR THE'
     +  /10X, 'SEARCH FOR THE STEADY STATE IS ALLOWED.')

99015 format
     +  (/1X, a9, 'ERROR.  UNKNOWN TASK.')

99016 format
     +  (/1X, a9, 'ERROR.  UNKNOWN REPORT CODE.')

99017 format
     +  (/1X, a9, 'ERROR.  SEARCH FAILS.')

99018 format
     +  (/1X, a9, 'ERROR.  REFINE FAILS.')

99019 format
     +  (/1X, a9, 'ERROR.  EVOLVE FAILS.')

99020 format
     +  (/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, i10, '  LABEL')

99021 format
     +  (/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, i10, '  RETURN')

C///  EXIT.

      stop
99999 continue
      return
      end
      subroutine twprep
     +  (error, text,
     +   a, asize, buffer, comps, condit, groupa, groupb, pivot, points,
     +   return)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWPREP
C
C     EVALUATE A BLOCK TRIDIAGONAL JACOBIAN MATRIX BY ONE-SIDED FINITE
C     DIFFERENCES AND REVERSE COMMUNICATION, PACK THE MATRIX INTO THE
C     LINPACK BANDED FORM, SCALE THE ROWS, AND FACTOR THE MATRIX USING
C     LINPACK'S SGBCO.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

      character
     +   id*9, string*80
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   a, absol, buffer, condit, delta, eps, relat, sum, temp
      external
     +   tweps, twgbco, twsqez
      integer
     +   asize, block, blocks, cfirst, clast, col, comps, count, diag,
     +   groupa, groupb, j, lda, length, lines, n, offset, pivot,
     +   points, rfirst, rlast, route, row, skip, text, width
      intrinsic
     +   abs, int, max, min, mod, sqrt
      logical
     +   error, found, mess, return

      parameter (id = 'TWPREP:  ')
      parameter (lines = 20)

      dimension
     +   a(asize), pivot(groupa + comps * points + groupb),
     +   buffer(groupa + comps * points + groupb)

      save

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  EVERY-TIME INITIALIZATION.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

C///  IF THIS IS A RETURN CALL, THEN CONTINUE WHERE THE PROGRAM PAUSED.

      if (return) then
         return = .false.
         go to (2030, 3050) route
         error = .true.
         go to 9001
      endif

C///  CHECK THE ARGUMENTS.

      n = groupa + comps * points + groupb
      error = .not. (((0 .lt. comps) .eqv. (0 .lt. points)) .and.
     +   0 .le. comps .and. 0 .le. points .and. 0 .le. groupa .and.
     +   0 .le. groupb .and. 0 .lt. n)
      if (error) go to 9002

      width = comps + max (comps, groupa, groupb) - 1
      error = .not. ((3 * width + 2) * n .le. asize)
      if (error) go to 9003

C///  WRITE ALL MESSAGES.

      if (mess .and. 0 .lt. text) then
         route = 0
         go to 9001
      end if

C///  FORM MACHINE EPSILON AND THE ABSOLUTE AND RELATIVE PERTURBATIONS.

      call tweps (eps)
      absol = sqrt (eps)
      relat = sqrt (eps)

C///  INITIALIZE COUNTERS AND POINTERS.

C     MAIN DIAGONAL ROW IN THE PACKING WHICH PLACES DIAGONALS IN ROWS

      diag = 2 * width + 1

C     PACKED ROW DIMENSION

      lda = 3 * width + 1
      skip = 2 * width + 1

C     BLOCKS AND BLOCK SIZES
C     ARRAY PIVOT HOLDS BLOCK SIZES AND POINTERS TEMPORARILY

      blocks = 0
      if (0 .lt. groupa) then
         blocks = blocks + 1
         pivot(blocks) = groupa
      end if

      do 1020 j = 1, points
         blocks = blocks + 1
         pivot(blocks) = comps
1020  continue

      if (0 .lt. groupb) then
         blocks = blocks + 1
         pivot(blocks) = groupb
      end if

C///////////////////////////////////////////////////////////////////////
C
C     (2) INITIALIZE THE COLUMNS OF THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

C///  STORE THE EVALUATION VECTOR.

      do 2010 j = 1, n
         a(j) = buffer(j)
2010  continue

C///  CLEAR THE MATRIX.

       do 2020 j = n + 1, (3 * width + 2) * n
          a(j) = 0.0
2020  continue

C///  EVALUATE THE FUNCTION AT THE UNPERTURBED X.

C     GO TO 2030 WHEN ROUTE = 1
      route = 1
      return = .true.
      go to 99999
2030  continue

C///  PLACE THE FUNCTION VALUES IN THE MATRIX.

      clast = 0
      do 2060 block = 1, blocks
         cfirst = clast + 1
         clast = clast + pivot(block)

         if (1 .lt. block) then
            rfirst = cfirst - pivot(block - 1)
         else
            rfirst = cfirst
         end if

         if (block .lt. blocks) then
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

C///////////////////////////////////////////////////////////////////////
C
C     (3) FORM THE COLUMNS OF THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

C///  TOP OF THE LOOP OVER GROUPS OF COLUMNS.

3010  continue
      found = .false.

C///  RESTORE THE EVALUATION VECTOR.

      do 3020 j = 1, n
         buffer(j) = a(j)
3020  continue

C///  PERTURB THE VECTOR AT INDEPENDENT POSITIONS.

      block = 1
      cfirst = 1
3030  continue
         if (0 .lt. pivot(block)) then
            found = .true.
            col = cfirst - 1 + pivot(block)
            if (0 .le. a(col)) then
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
            if (block .eq. 1 .and. 0 .lt. groupa) then
               cfirst = cfirst + groupa
            else if (block .eq. blocks .and. 0 .lt. groupb) then
               cfirst = cfirst + groupb
            else
               cfirst = cfirst + comps
            end if
            block = block + 1
3040     continue

      if (block .le. blocks) go to 3030

C///  EXIT OF THE LOOP OVER GROUPS OF COLUMNS.

      if (.not. found) go to 3090

C///  EVALUATE THE FUNCTION AT THE PERTURBED VALUES.

C     GO TO 3050 WHEN ROUTE = 2
      route = 2
      return = .true.
      go to 99999
3050  continue

C///  DIFFERENCE TO FORM THE COLUMNS OF THE JACOBIAN MATRIX.

      block = 1
      cfirst = 1
3060  continue
         if (0 .lt. pivot(block)) then
            col = cfirst - 1 + pivot(block)
            pivot(block) = pivot(block) - 1

            if (0 .le. a(col)) then
               delta = relat * a(col) + absol
            else
               delta = relat * a(col) - absol
            end if
            temp = 1.0 / delta
            offset = n + diag - col + lda * (col - 1)

            if (block .eq. 1 .and. 0 .lt. groupa) then
               clast = cfirst + groupa - 1
            else if (block .eq. blocks .and. 0 .lt. groupb) then
               clast = cfirst + groupb - 1
            else
               clast = cfirst + comps - 1
            end if

            if (1 .lt. block) then
               if (block .eq. 2 .and. 0 .lt. groupa) then
                  rfirst = cfirst - groupa
               else
                  rfirst = cfirst - comps
               end if
            else
               rfirst = cfirst
            end if

            if (block .lt. blocks) then
               if (block .eq. blocks - 1 .and. 0 .lt. groupb) then
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
            if (block .eq. 1 .and. 0 .lt. groupa) then
               cfirst = cfirst + groupa
            else if (block .eq. blocks .and. 0 .lt. groupb) then
               cfirst = cfirst + groupb
            else
               cfirst = cfirst + comps
            end if
            block = block + 1
3080     continue

      if (block .le. blocks) go to 3060

C///  BOTTOM OF THE LOOP OVER GROUPS OF COLUMNS.

      go to 3010
3090  continue

C///////////////////////////////////////////////////////////////////////
C
C     (4) CHECK FOR ZERO COLUMNS.
C
C///////////////////////////////////////////////////////////////////////

      count = 0
      do 4020 col = 1, n
         offset = n + diag - col + lda * (col - 1)
         sum = 0.0
         do 4010 row = max (col - width, 1), min (col + width, n)
            sum = sum + abs (a(offset + row))
4010     continue
         a(col) = sum

         if (sum .eq. 0.0) count = count + 1
4020  continue

      error = .not. (count .eq. 0)
      if (error) go to 9004

C///////////////////////////////////////////////////////////////////////
C
C     (5) SCALE THE ROWS.
C
C///////////////////////////////////////////////////////////////////////

      count = 0
      do 5030 row = 1, n
         offset = n + diag + row
         sum = 0.0
         do 5010 col = max (row - width, 1), min (row + width, n)
            sum = sum + abs (a(offset - col + lda * (col - 1)))
5010     continue

         if (sum .eq. 0.0) then
            count = count + 1
            a(row) = sum
         else
            temp = 1.0 / sum
            a(row) = temp

            do 5020 col = max (row - width, 1), min (row + width, n)
               a(offset - col + lda * (col - 1))
     +            = a(offset - col + lda * (col - 1)) * temp
5020        continue
         endif
5030  continue

      error = .not. (count .eq. 0)
      if (error) go to 9005

C///////////////////////////////////////////////////////////////////////
C
C     (6) FACTOR THE MATRIX.
C
C///////////////////////////////////////////////////////////////////////

      call twgbco
     +  (a(n + 1), lda, n, width, width, pivot, condit, buffer)

      error = condit .eq. 0.0
      if (error) go to 9006

      condit = 1.0 / condit

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

80001 format
     +  (10X, a)

80002 format
     +  (10X, '... MORE')

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      go to 99999

9001  if (0 .lt. text) write (text, 99001) id, route
      if (.not. mess) go to 99999

9002  if (0 .lt. text) write (text, 99002) id,
     +   comps, points, groupa, groupb, n
      if (.not. mess) go to 99999

9003  if (0 .lt. text) write (text, 99003) id,
     +   comps, points, groupa, groupb, n, width,
     +   (3 * width + 2) * n, asize
      if (.not. mess) go to 99999

9004  if (0 .lt. text) then
         write (text, 99004) id, comps, points, groupa, groupb,
     +      groupa + comps * points + groupb, count
         count = 0
         do 8010 j = 1, groupa + comps * points + groupb
            if (a(j) .eq. 0.0 .or. mess) then
               count = count + 1
               if (count .le. lines) then
                  if (j .le. groupa) then
                     write (string, '(A, I10)') 'GROUP A ', j
                  else if (j .le. groupa + comps * points) then
                     write (string, '(A, I10, A, I10)')
     +                  ' COMPONENT ', mod (j - groupa - 1, comps) + 1,
     +                  ' AT POINT ', int ((j - groupa - 1) / comps) + 1
                  else
                     write (string, '(A, I10)')
     +                  'GROUP B ', j - groupa - comps * points
                  end if
                  call twsqez (length, string)
                  write (text, 80001) string (1 : length)
               end if
            end if
8010     continue
         if (lines .lt. count) write (text, 80002)
      end if
      if (.not. mess) go to 99999

9005  if (0 .lt. text) then
         write (text, 99005) id, comps, points, groupa, groupb,
     +      groupa + comps * points + groupb, count
         count = 0
         do 8020 j = 1, groupa + comps * points + groupb
            if (a(j) .eq. 0.0 .or. mess) then
               count = count + 1
               if (count .le. lines) then
                  if (j .le. groupa) then
                     write (string, '(A, I10)') 'GROUP A ', j
                  else if (j .le. groupa + comps * points) then
                     write (string, '(A, I10, A, I10)')
     +                  ' COMPONENT ', mod (j - groupa - 1, comps) + 1,
     +                  ' AT POINT ', int ((j - groupa - 1) / comps) + 1
                  else
                     write (string, '(A, I10)')
     +                  'GROUP B ', j - groupa - comps * points
                  end if
                  call twsqez (length, string)
                  write (text, 80001) string (1 : length)
               end if
            end if
8020     continue
         if (lines .lt. count) write (text, 80002)
      end if
      if (.not. mess) go to 99999

9006  if (0 .lt. text) write (text, 99006) id
      if (.not. mess) go to 99999

99001 format
     +  (/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.'
     + //10X, i10, '  ROUTE')

99002 format
     +  (/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES'
     +  /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  TOTAL UNKNOWNS')

99003 format
     +  (/1X, a9, 'ERROR.  THE MATRIX SPACE IS TOO SMALL.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  MATRIX ORDER'
     +  /10X, i10, '  STRICT HALF BANDWIDTH'
     + //10X, i10, '  SPACE REQUIRED'
     +  /10X, i10, '  ASIZE, PROVIDED')

99004 format
     +  (/1X, a9, 'ERROR.  SOME COLUMNS ARE ZERO.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  TOTAL COLUMNS'
     +  /10X, i10, '  ZERO COLUMNS'
     + //10X, 'UNKNOWNS WITH ZERO COLUMNS:'
     +  /)

99005 format
     +  (/1X, a9, 'ERROR.  SOME ROWS ARE ZERO.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  TOTAL ROWS'
     +  /10X, i10, '  ZERO ROWS'
     + //10X, 'ZERO ROWS:'
     +  /)

99006 format
     +  (/1X, a9, 'ERROR.  THE JACOBIAN MATRIX IS SINGULAR.')

C///  EXIT.

      stop
99999 continue
      return
      end
      subroutine twseti (error, text, contrl, value)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSETI
C
C     SET A CONTROL THAT TAKES AN INTEGER VALUE.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character
     +   contrl*(*), id*9, string*80
      external
     +   twinit, twlast
      integer
     +   cntrls, count, ivalue, length, text, value
      logical
     +   error, found, mess, lvalue

      parameter (id = 'TWSETI:  ')
      parameter (cntrls = 22)

      dimension ivalue(cntrls), lvalue(cntrls)

      common / twcomi / ivalue
      common / twcoml / lvalue

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

      if (mess .and. 0 .lt. text) go to 9001

C///  INITIALIZE THE CONTROLS.

      call twinit (error, text, .false.)
      if (error) go to 9001

C///  SET THE CONTROLS.

      count = 0
      found = .false.

C     ADAPT

      count = count + 1
      if (contrl .eq. 'ADAPT') then
         error = .true.
         go to 9002
      end if

C     LEVELD

      count = count + 1
      if (contrl .eq. 'LEVELD') then
         found = .true.
         ivalue(count) = value
      end if

C     LEVELM

      count = count + 1
      if (contrl .eq. 'LEVELM') then
         found = .true.
         ivalue(count) = value
      end if

C     PADD

      count = count + 1
      if (contrl .eq. 'PADD') then
         found = .true.
         ivalue(count) = value
         lvalue(count) = .true.
      end if

C     SSABS

      count = count + 1
      if (contrl .eq. 'SSABS') then
         error = .true.
         go to 9003
      end if

C     SSAGE

      count = count + 1
      if (contrl .eq. 'SSAGE') then
         found = .true.
         ivalue(count) = value
      end if

C     SSREL

      count = count + 1
      if (contrl .eq. 'SSREL') then
         error = .true.
         go to 9003
      end if

C     STEADY

      count = count + 1
      if (contrl .eq. 'STEADY') then
         error = .true.
         go to 9002
      end if

C     STEPS0

      count = count + 1
      if (contrl .eq. 'STEPS0') then
         found = .true.
         ivalue(count) = value
      end if

C     STEPS1

      count = count + 1
      if (contrl .eq. 'STEPS1') then
         found = .true.
         ivalue(count) = value
      end if

C     STEPS2

      count = count + 1
      if (contrl .eq. 'STEPS2') then
         found = .true.
         ivalue(count) = value
      end if

C     STRID0

      count = count + 1
      if (contrl .eq. 'STRID0') then
         error = .true.
         go to 9003
      end if

C     TDABS

      count = count + 1
      if (contrl .eq. 'TDABS') then
         error = .true.
         go to 9003
      end if

C     TDAGE

      count = count + 1
      if (contrl .eq. 'TDAGE') then
         found = .true.
         ivalue(count) = value
      end if

C     TDEC

      count = count + 1
      if (contrl .eq. 'TDEC') then
         error = .true.
         go to 9003
      end if

C     TDREL

      count = count + 1
      if (contrl .eq. 'TDREL') then
         error = .true.
         go to 9003
      end if

C     TINC

      count = count + 1
      if (contrl .eq. 'TINC') then
         error = .true.
         go to 9003
      end if

C     TMAX

      count = count + 1
      if (contrl .eq. 'TMAX') then
         error = .true.
         go to 9003
      end if

C     TMIN

      count = count + 1
      if (contrl .eq. 'TMIN') then
         error = .true.
         go to 9003
      end if

C     TOLER0

      count = count + 1
      if (contrl .eq. 'TOLER0') then
         error = .true.
         go to 9003
      end if

C     TOLER1

      count = count + 1
      if (contrl .eq. 'TOLER1') then
         error = .true.
         go to 9003
      end if

C     TOLER2

      count = count + 1
      if (contrl .eq. 'TOLER2') then
         error = .true.
         go to 9003
      end if

      error = .not. (count .eq. cntrls)
      if (error) go to 9004

      error = .not. found
      if (error) go to 9005

C///  ERROR MESSAGES.

      go to 99999

9001  if (0 .lt. text) write (text, 99001) id
      if (.not. mess) go to 99999

9002  if (0 .lt. text) then
         call twlast (length, contrl)
         if (length .le. 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99002) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9003  if (0 .lt. text) then
         call twlast (length, contrl)
         if (length .le. 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99003) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9004  if (0 .lt. text) write (text, 99004) id, cntrls, count
      if (.not. mess) go to 99999

9005  if (0 .lt. text) then
         call twlast (length, contrl)
         if (length .le. 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99005) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

99001 format
     +  (/1X, a9, 'ERROR.  TWINIT FAILS.')

99002 format
     +  (/1X, a9, 'ERROR.  THE CONTROL TAKES A LOGICAL VALUE WHICH'
     +  /10X, 'MUST BE SET USING TWSETL.'
     + //10X, '     CONTROL:  ', a)

99003 format
     +  (/1X, a9, 'ERROR.  THE CONTROL TAKES A REAL VALUE WHICH MUST BE'
     +  /10X, 'SET USING TWSETR.'
     + //10X, '     CONTROL:  ', a)

99004 format
     +  (/1X, a9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.'
     + //10X, i10, '  CONTROLS'
     +  /10X, i10, '  COUNTED')

99005 format
     +  (/1X, a9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.'
     + //10X, '     CONTROL:  ', a)

C///  EXIT.

      stop
99999 continue
      return
      end
      subroutine twsetl (error, text, contrl, value)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSETL
C
C     SET A CONTROL THAT TAKES A LOGICAL VALUE.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character
     +   contrl*(*), id*9, string*80
      external
     +   twinit, twlast
      integer
     +   cntrls, count, length, text
      logical
     +   error, found, lvalue, mess, value

      parameter (id = 'TWSETL:  ')
      parameter (cntrls = 22)

      dimension lvalue(cntrls)

      common / twcoml / lvalue

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

      if (mess .and. 0 .lt. text) go to 9001

C///  INITIALIZE THE CONTROLS.

      call twinit (error, text, .false.)
      if (error) go to 9001

C///  SET THE CONTROLS.

      count = 0
      found = .false.

C     ADAPT

      count = count + 1
      if (contrl .eq. 'ADAPT') then
         found = .true.
         lvalue(count)= value
      end if

C     LEVELD

      count = count + 1
      if (contrl .eq. 'LEVELD') then
         error = .true.
         go to 9002
      end if

C     LEVELM

      count = count + 1
      if (contrl .eq. 'LEVELM') then
         error = .true.
         go to 9002
      end if

C     PADD

      count = count + 1
      if (contrl .eq. 'PADD') then
         error = .true.
         go to 9002
      end if

C     SSABS

      count = count + 1
      if (contrl .eq. 'SSABS') then
         error = .true.
         go to 9003
      end if

C     SSAGE

      count = count + 1
      if (contrl .eq. 'SSAGE') then
         error = .true.
         go to 9002
      end if

C     SSREL

      count = count + 1
      if (contrl .eq. 'SSREL') then
         error = .true.
         go to 9003
      end if

C     STEADY

      count = count + 1
      if (contrl .eq. 'STEADY') then
         found = .true.
         lvalue(count)= value
      end if

C     STEPS0

      count = count + 1
      if (contrl .eq. 'STEPS0') then
         error = .true.
         go to 9002
      end if

C     STEPS1

      count = count + 1
      if (contrl .eq. 'STEPS1') then
         error = .true.
         go to 9002
      end if

C     STEPS2

      count = count + 1
      if (contrl .eq. 'STEPS2') then
         error = .true.
         go to 9002
      end if

C     STRID0

      count = count + 1
      if (contrl .eq. 'STRID0') then
         error = .true.
         go to 9003
      end if

C     TDABS

      count = count + 1
      if (contrl .eq. 'TDABS') then
         error = .true.
         go to 9003
      end if

C     TDAGE

      count = count + 1
      if (contrl .eq. 'TDAGE') then
         error = .true.
         go to 9002
      end if

C     TDEC

      count = count + 1
      if (contrl .eq. 'TDEC') then
         error = .true.
         go to 9003
      end if

C     TDREL

      count = count + 1
      if (contrl .eq. 'TDREL') then
         error = .true.
         go to 9003
      end if

C     TINC

      count = count + 1
      if (contrl .eq. 'TINC') then
         error = .true.
         go to 9003
      end if

C     TMAX

      count = count + 1
      if (contrl .eq. 'TMAX') then
         error = .true.
         go to 9003
      end if

C     TMIN

      count = count + 1
      if (contrl .eq. 'TMIN') then
         error = .true.
         go to 9003
      end if

C     TOLER0

      count = count + 1
      if (contrl .eq. 'TOLER0') then
         error = .true.
         go to 9003
      end if

C     TOLER1

      count = count + 1
      if (contrl .eq. 'TOLER1') then
         error = .true.
         go to 9003
      end if

C     TOLER2

      count = count + 1
      if (contrl .eq. 'TOLER2') then
         error = .true.
         go to 9003
      end if

      error = .not. (count .eq. cntrls)
      if (error) go to 9004

      error = .not. found
      if (error) go to 9005

C///  ERROR MESSAGES.

      go to 99999

9001  if (0 .lt. text) write (text, 99001) id
      if (.not. mess) go to 99999

9002  if (0 .lt. text) then
         call twlast (length, contrl)
         if (length .le. 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99002) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9003  if (0 .lt. text) then
         call twlast (length, contrl)
         if (length .le. 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99003) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9004  if (0 .lt. text) write (text, 99004) id, cntrls, count
      if (.not. mess) go to 99999

9005  if (0 .lt. text) then
         call twlast (length, contrl)
         if (length .le. 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99005) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

99001 format
     +  (/1X, a9, 'ERROR.  TWINIT FAILS.')

99002 format
     +  (/1X, a9, 'ERROR.  THE CONTROL TAKES AN INTEGER VALUE WHICH'
     +  /10X, 'MUST BE SET USING TWSETI.'
     + //10X, '     CONTROL:  ', a)

99003 format
     +  (/1X, a9, 'ERROR.  THE CONTROL TAKES A REAL VALUE WHICH MUST BE'
     +  /10X, 'SET USING TWSETR.'
     + //10X, '     CONTROL:  ', a)

99004 format
     +  (/1X, a9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.'
     + //10X, i10, '  CONTROLS'
     +  /10X, i10, '  COUNTED')

99005 format
     +  (/1X, a9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.'
     + //10X, '     CONTROL:  ', a)

C///  EXIT.

      stop
99999 continue
      return
      end
      subroutine twsetr (error, text, contrl, value)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSETR
C
C     SET A CONTROL THAT TAKES A REAL VALUE.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character
     +   contrl*(*), id*9, string*80
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   rvalue, value
      external
     +   twinit, twlast
      integer
     +   cntrls, count, length, text
      logical
     +   error, found, mess

      parameter (id = 'TWSETR:  ')
      parameter (cntrls = 22)

      dimension rvalue(cntrls)

      common / twcomr / rvalue

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

      if (mess .and. 0 .lt. text) go to 9001

C///  INITIALIZE THE CONTROLS.

      call twinit (error, text, .false.)
      if (error) go to 9001

C///  SET THE CONTROLS.

      count = 0
      found = .false.

C     ADAPT

      count = count + 1
      if (contrl .eq. 'ADAPT') then
         error = .true.
         go to 9002
      end if

C     LEVELD

      count = count + 1
      if (contrl .eq. 'LEVELD') then
         error = .true.
         go to 9003
      end if

C     LEVELM

      count = count + 1
      if (contrl .eq. 'LEVELM') then
         error = .true.
         go to 9003
      end if

C     PADD

      count = count + 1
      if (contrl .eq. 'PADD') then
         error = .true.
         go to 9003
      end if

C     SSABS

      count = count + 1
      if (contrl .eq. 'SSABS') then
         found = .true.
         rvalue(count) = value
      end if

C     SSAGE

      count = count + 1
      if (contrl .eq. 'SSAGE') then
         error = .true.
         go to 9003
      end if

C     SSREL

      count = count + 1
      if (contrl .eq. 'SSREL') then
         found = .true.
         rvalue(count) = value
      end if

C     STEADY

      count = count + 1
      if (contrl .eq. 'STEADY') then
         error = .true.
         go to 9002
      end if

C     STEPS0

      count = count + 1
      if (contrl .eq. 'STEPS0') then
         error = .true.
         go to 9003
      end if

C     STEPS1

      count = count + 1
      if (contrl .eq. 'STEPS1') then
         error = .true.
         go to 9003
      end if

C     STEPS2

      count = count + 1
      if (contrl .eq. 'STEPS2') then
         error = .true.
         go to 9003
      end if

C     STRID0

      count = count + 1
      if (contrl .eq. 'STRID0') then
         found = .true.
         rvalue(count) = value
      end if

C     TDABS

      count = count + 1
      if (contrl .eq. 'TDABS') then
         found = .true.
         rvalue(count) = value
      end if

C     TDAGE

      count = count + 1
      if (contrl .eq. 'TDAGE') then
         error = .true.
         go to 9003
      end if

C     TDEC

      count = count + 1
      if (contrl .eq. 'TDEC') then
         found = .true.
         rvalue(count) = value
      end if

C     TDREL

      count = count + 1
      if (contrl .eq. 'TDREL') then
         found = .true.
         rvalue(count) = value
      end if

C     TINC

      count = count + 1
      if (contrl .eq. 'TINC') then
         found = .true.
         rvalue(count) = value
      end if

C     TMAX

      count = count + 1
      if (contrl .eq. 'TMAX') then
         found = .true.
         rvalue(count) = value
      end if

C     TMIN

      count = count + 1
      if (contrl .eq. 'TMIN') then
         found = .true.
         rvalue(count) = value
      end if

C     TOLER0

      count = count + 1
      if (contrl .eq. 'TOLER0') then
         found = .true.
         rvalue(count) = value
      end if

C     TOLER1

      count = count + 1
      if (contrl .eq. 'TOLER1') then
         found = .true.
         rvalue(count) = value
      end if

C     TOLER2

      count = count + 1
      if (contrl .eq. 'TOLER2') then
         found = .true.
         rvalue(count) = value
      end if

      error = .not. (count .eq. cntrls)
      if (error) go to 9004

      error = .not. found
      if (error) go to 9005

C///  ERROR MESSAGES.

      go to 99999

9001  if (0 .lt. text) write (text, 99001) id
      if (.not. mess) go to 99999

9002  if (0 .lt. text) then
         call twlast (length, contrl)
         if (length .le. 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99002) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9003  if (0 .lt. text) then
         call twlast (length, contrl)
         if (length .le. 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99003) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

9004  if (0 .lt. text) write (text, 99004) id, cntrls, count
      if (.not. mess) go to 99999

9005  if (0 .lt. text) then
         call twlast (length, contrl)
         if (length .le. 40) then
            string = contrl
         else
            length = 40
            string = contrl (1 : 37) // '...'
         end if
         write (text, 99005) id, string (1 : length)
      end if
      if (.not. mess) go to 99999

99001 format
     +  (/1X, a9, 'ERROR.  TWINIT FAILS.')

99002 format
     +  (/1X, a9, 'ERROR.  THE CONTROL TAKES A LOGICAL VALUE WHICH'
     +  /10X, 'MUST BE SET USING TWSETL.'
     + //10X, '     CONTROL:  ', a)

99003 format
     +  (/1X, a9, 'ERROR.  THE CONTROL TAKES AN INTEGER VALUE WHICH'
     +  /10X, 'MUST BE SET USING TWSETI.'
     + //10X, '     CONTROL:  ', a)

99004 format
     +  (/1X, a9, 'ERROR.  THE NUMBER OF CONTROLS IS INCONSISTENT.'
     + //10X, i10, '  CONTROLS'
     +  /10X, i10, '  COUNTED')

99005 format
     +  (/1X, a9, 'ERROR.  THE CONTROL IS NOT RECOGNIZED.'
     + //10X, '     CONTROL:  ', a)

C///  EXIT.

      stop
99999 continue
      return
      end
      subroutine twshow
     +  (error, text,
     +   buffer, comps, grid, groupa, groupb, points, x)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSHOW
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

      character
     +   id*9, string*80, title*80
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   buffer, x
      external
     +   twsqez
      integer
     +   cols, comp, comps, count, first, groupa, groupb, groups, j,
     +   last, length, point, points, text
      intrinsic
     +   min
      logical
     +   error, grid, mess

      parameter (id = 'TWSHOW:  ')

      dimension
     +   buffer(groupa + comps * points + groupb), title(6), x(*)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

      if (mess .and. 0 .lt. text) go to 9001

C///  CHECK THE ARGUMENTS.

      error = .not. (((0 .lt. comps) .eqv. (0 .lt. points)) .and.
     +   0 .le. comps .and. 0 .le. points .and. 0 .le. groupa .and.
     +   0 .le. groupb .and. 0 .lt. groupa + comps * points + groupb)
      if (error) go to 9001

C///  COUNT THE GROUPS.

      groups = 0
      if (0 .lt. groupa) groups = groups + 1
      if (0 .lt. groupb) groups = groups + 1
      if (0 .lt. comps .and. 0 .lt. points) groups = groups + 1

C///  CHOOSE NUMBER OF DATA COLUMNS.

      if (grid) then
         cols = 5
      else
         cols = 6
      end if

C///////////////////////////////////////////////////////////////////////
C
C     (2) PRINT THE GROUPED DATA.
C
C///////////////////////////////////////////////////////////////////////

      if (0 .lt. text) then

      if (0 .lt. groupa) then
         if (1 .lt. groups) write (text, 10001) 'GROUP A UNKNOWNS'
         write (text, 10002) (j, buffer(j), j = 1, groupa)
      end if

      if (0 .lt. groupb) then
         if (1 .lt. groups) write (text, 10001) 'GROUP B UNKNOWNS'
         write (text, 10002)
     +      (j, buffer(groupa + comps * points + j), j = 1, groupb)
      end if

C///////////////////////////////////////////////////////////////////////
C
C     (2) PRINT THE COMPONENTS AT POINTS.
C
C///////////////////////////////////////////////////////////////////////

      if (0 .lt. comps .and. 0 .lt. points) then
         if (1 .lt. groups) write (text, 10001) 'COMPONENTS AT POINTS'

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
               write (text, 10003)
     +            'GRID POINT', (title(j), j = 1, count)
            else
               write (text, 10003) (title(j), j = 1, count)
            end if

            if (count .eq. cols) then
               if (grid) then
                  write (text, 10004) (point, x(point),
     +               (buffer(groupa + comp + comps * (point - 1)),
     +               comp = first, last), point = 1, points)
               else
                  write (text, 10005) (point,
     +               (buffer(groupa + comp + comps * (point - 1)),
     +               comp = first, last), point = 1, points)
               end if
            else
               do 2020 point = 1, points
                  if (grid) then
                     write (text, 10004) point, x(point),
     +                  (buffer(groupa + comp + comps * (point - 1)),
     +                  comp = first, last)
                  else
                     write (text, 10005) point,
     +                  (buffer(groupa + comp + comps * (point - 1)),
     +                  comp = first, last)
                  end if
2020           continue
            end if
2030     continue
      end if

      end if

C///////////////////////////////////////////////////////////////////////
C
C     INFORMATIVE MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

10001 format
     +  (/10X, a)

10002 format
     + (/(10X, 4(i3, '> ', 1pe10.3)))

10003 format
     +  (/14X, 6(1X, a10))

10004 format
     +  (10X, 0p, i3, '>', f11.6, 1p, 5E11.3)

10005 format
     +  (10X, 0p, i3, '>', 1p, 6E11.3)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      go to 99999

9001  if (0 .lt. text) write (text, 99001) id,
     +   comps, points, groupa, groupb, groupa + comps * points + groupb
      if (.not. mess) go to 99999

99001 format
     +  (/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES'
     +  /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  TOTAL UNKNOWNS')
      if (.not. mess) go to 99999

C///  EXIT.

      stop
99999 continue
      return
      end
      subroutine twsolv
     +  (error, text,
     +   a, asize, buffer, comps, groupa, groupb, pivot, points)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSOLV
C
C     SOLVE A SYSTEM OF LINEAR EQUATIONS USING THE MATRIX PREPARED BY
C     TWPREP.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)

      character
     +   id*9
C*****PRECISION > DOUBLE
      double precision
C*****END PRECISION > DOUBLE
C*****PRECISION > SINGLE
C      REAL
C*****END PRECISION > SINGLE
     +   a, buffer
      external
     +   twgbsl
      integer
     +   asize, comps, groupa, groupb, j, n, pivot, points, text, width
      intrinsic
     +   max
      logical
     +   error, mess

      parameter (id = 'TWSOLV:  ')

      dimension
     +   a(asize), buffer(groupa + comps * points + groupb),
     +   pivot(groupa + comps * points + groupb)

C///////////////////////////////////////////////////////////////////////
C
C     (1) PROLOGUE.
C
C///////////////////////////////////////////////////////////////////////

C///  WRITE ALL MESSAGES.

C     SET TRUE TO PRINT EXAMPLES OF ALL MESSAGES.
      mess = .false.

      if (mess .and. 0 .lt. text) go to 9001

C///  CHECK THE ARGUMENTS.

      n = groupa + comps * points + groupb
      error = .not. (((0 .lt. comps) .eqv. (0 .lt. points)) .and.
     +   0 .le. comps .and. 0 .le. points .and. 0 .le. groupa .and.
     +   0 .le. groupb .and. 0 .lt. n)
      if (error) go to 9001

      width = comps + max (comps, groupa, groupb) - 1
      error = .not. ((3 * width + 2) * n .le. asize)
      if (error) go to 9002

C///////////////////////////////////////////////////////////////////////
C
C     (2) SCALE AND SOLVE THE EQUATIONS.
C
C///////////////////////////////////////////////////////////////////////

      do 2010 j = 1, n
         buffer(j) = buffer(j) * a(j)
2010  continue

      call twgbsl
     +  (a(n + 1), 3 * width + 1, n, width, width, pivot, buffer)

C///////////////////////////////////////////////////////////////////////
C
C     ERROR MESSAGES.
C
C///////////////////////////////////////////////////////////////////////

      go to 99999

9001  if (0 .lt. text) write (text, 99001) id,
     +   comps, points, groupa, groupb, n
      if (.not. mess) go to 99999

9002  if (0 .lt. text) write (text, 99002) id,
     +   comps, points, groupa, groupb, n, width,
     +   (3 * width + 2) * n, asize
      if (.not. mess) go to 99999

99001 format
     +  (/1X, a9, 'ERROR.  NUMBERS OF COMPONENTS AND POINTS MUST BE'
     +  /10X, 'EITHER BOTH ZERO OR BOTH POSITIVE, NUMBERS OF ALL TYPES'
     +  /10X, 'OF UNKNOWNS MUST BE AT LEAST ZERO, AND TOTAL UNKNOWNS'
     +  /10X, 'MUST BE POSITIVE.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  TOTAL UNKNOWNS')

99002 format
     +  (/1X, a9, 'ERROR.  THE MATRIX SPACE IS TOO SMALL.'
     + //10X, i10, '  COMPS, COMPONENTS'
     +  /10X, i10, '  POINTS'
     +  /10X, i10, '  GROUPA, GROUP A UNKNOWNS'
     +  /10X, i10, '  GROUPB, GROUP B UNKNOWNS'
     +  /10X, i10, '  MATRIX ORDER'
     +  /10X, i10, '  STRICT HALF BANDWIDTH'
     + //10X, i10, '  SPACE EXPECTED'
     +  /10X, i10, '  ASIZE, PROVIDED')

C///  EXIT.

      stop
99999 continue
      return
      end
      subroutine twsqez (length, string)

C///////////////////////////////////////////////////////////////////////
C
C     T W O P N T
C
C     TWSQEZ
C
C     SQUEEZE LEADING BLANKS AND MULTIPLE BLANKS FROM A CHARACTER
C     STRING.  RETURN THE LENGTH OF THE SQUEEZED STRING.
C
C///////////////////////////////////////////////////////////////////////

      implicit complex (a - z)
      character
     +   char*1, string*(*)
      integer
     +   j, length
      intrinsic
     +   len
      logical
     +   blank

C///  SQUEEZE THE STRING.

      length = 0
      blank = .true.
      do 0100 j = 1, len (string)
         char = string (j : j)
         if (.not. blank .or. char .ne. ' ') then
            blank = char .eq. ' '
            length = length + 1
            string (length : length) = char
         end if
0100  continue

C///  ADJUST THE LENGTH AND PAD THE STRING.

      if (0 .lt. length) then
         if (string (length : length) .eq. ' ') length = length - 1
         if (length .lt. len (string)) string (length + 1 : ) = ' '
      else
         length = 1
      end if

C///  EXIT.

      return
      end

