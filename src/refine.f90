      ! Perform automatic grid selection
      subroutine refine(error, text, &
                        active, buffer, comps, leveld, levelm, mark, newx, padd, pmax, &
                        points, ratio, ratio1, ratio2, signal, succes, toler0, toler1, &
                        toler2, u, vary1, vary2, weight, x)

      use twopnt_core, only: twcopy
      implicit none

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

         if (0<text) write (text, 99001) id, route
         return
      end if

      ! write all messages.
      if (mess .and. 0<text) then
         route = 0
         write (text, 10001) id
         write (text, 10002) id
         write (text, 10004) id
         write (text, 10005) id
         write (text, 10010) id
         write (text, 99001) id, route
         write (text, 99002) id, comps, points
         write (text, 99003) id, padd
         write (text, 99004) id, points, pmax
         write (text, 99005) id
         write (text, 99006) id, toler0
         write (text, 99007) id, toler1, toler2
         write (text, 99008) id
         write (text, 99009) id
         stop
      end if

      ! Levelm printing.
      if (0<levelm .and. 0<text) write (text, 10001) id

      ! Check the arguments.
      error = .not. (1 <= comps .and. 2 <= points)
      if (error) then
           if (0<text) write (text, 99002) id, comps, points
           return
      end if

      error = .not. (0 <= padd)
      if (error) then
           if (0<text) write (text, 99003) id, padd
           return
      endif

      error = .not. (points <= pmax)
      if (error) then
           if (0<text) write (text, 99004) id, points, pmax
           return
      end if

      count = 0
      do 1010 j = 1, comps
         if (active(j)) count = count + 1
1010  continue
      error = .not. (1 <= count)
      if (error) then
           if (0<text) write (text, 99005) id
           return
      end if

      error = .not. toler0>=zero
      if (error) then
          if (0<text) write (text, 99006) id, toler0
          return
      end if

      error = .not. (zero <= toler1 .and. toler1 <= one &
              .and.  zero <= toler2 .and. toler2 <= one)
      if (error) then
          if (0<text) write (text, 99007) id, toler1, toler2
          return
      end if

      count = 0
      do 1020 k = 1, points - 1
         if (x(k)<x(k + 1)) count = count + 1
1020  continue
      error = .not. (count == 0 .or. count == points - 1)
      if (error) then
          if (0<text) write (text, 99008) id
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
          if (0<text) write (text, 99009) id
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
            write (text, 10002) id
         else
            write (text, 10003) temp1, temp2, toler1, toler2
            if (most == 0) then
               write (text, 10004) id
            else if (more == 0) then
               write (text, 10005) id
            else
               write (text, 10006)

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
                           write (text, 10007) k - 1, x(k - 1), word
                        else
                           write (text, 10008) word
                        end if
                     end if

                     if (1<k .and. k<points) then
                        if (vary2(old)/=zero) then
                           write (word, '(F4.2, I4)') ratio2(old), vary2(old)
                        else
                           write (word, '(F4.2, I4)') ratio2(old)
                        end if
                        write (text, 10009) k, x(k), word
                     else
                        write (text, 10009) k, x(k)
                     end if
                  end if
5030           continue
            end if
         end if

         if (0<leveld .and. 0<more) then
            write (text, 10010) id
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
            10001 format(/1X, a9, 'SELECT A GRID.')
            10002 format(/1X, a9, 'SUCCESS.  THE GRID IS ADEQUATE BECAUSE ALL ACTIVE' &
                        /10X, 'COMPONENTS ARE INSIGNIFICANT.')
            !                   123456789-   1234567   1234567
            10003 format(/15X, '             RATIO 1   RATIO 2' &
                         /15X, '             -------   -------' &
                         /15X, '    ACTUAL', 2F10.3 &
                         /15X, '   DESIRED', 2F10.3)
            10004 format(/1X, a9, 'SUCCESS.  THE GRID IS ADEQUATE.')
            10005 format(/1X, a9, 'FAILURE.  MORE POINTS ARE NEEDED BUT NONE CAN BE ADDED.')
            !                   123456   123456789-123456   12345678   12345678
            10006 format(/10X, 'THE NEW GRID (* MARKS NEW POINTS):' &
                        //10X, '                             LARGEST RATIOS AND'&
                         /10X, ' INDEX         GRID POINT      NUMBER TOO LARGE'&
                         /10X, '------   ----------------   -------------------'&
                         /10X, '                             RATIO 1    RATIO 2')
            10007 format(10X, i6, '*  ', 1p, e16.9, 0p, 3X, a8)
            10008 format(38X, a8)
            10009 format(10X, i6, '   ', 1p, e16.9, 0p, 14X, a8)
            10010 format(/1X, a9, 'THE SOLUTION GUESS FOR THE NEW GRID:')

            ! Error messages.
            99001 format(/1X, a9, 'ERROR.  THE COMPUTED GOTO IS OUT OF RANGE.' &
                       //10X, i10, '  ROUTE')
            99002 format(/1X, a9, 'ERROR.  THERE MUST BE AT LEAST ONE COMPONENT AND AT' &
                        /10X, 'LEAST TWO POINTS.' &
                       //10X, i10, '  COMPS, COMPONENTS' &
                        /10X, i10, '  POINTS')
            99003 format(/1X, a9, 'ERROR.  THE LIMIT ON POINTS ADDED TO A GRID MUST BE' &
                        /10X, 'ZERO OR POSITIVE.'&
                       //10X, i10, '  PADD, LIMIT ON ADDED POINTS')
            99004 format(/1X, a9, 'ERROR.  POINTS IS OUT OF RANGE.' &
                       //10X, i10, '  POINTS'&
                        /10X, i10, '  PMAX, LIMIT ON POINTS')
            99005 format(/1X, a9, 'ERROR.  THERE ARE NO ACTIVE COMPONENTS.')
            99006 format(/1X, a9, 'ERROR.  THE BOUNDS ON MAGNITUDE AND RELATIVE CHANGE' &
                        /10X, 'OF MAGNITUDE FOR INSIGNIFICANT COMPONENTS MUST BE'&
                        /10X, 'POSITIVE.'&
                       //10X, 1p, e10.2, '  TOLER0, SIGNIFICANCE LEVEL')
            99007 format(/1X, a9, 'ERROR.  THE BOUNDS ON RELATIVE CHANGES IN MAGNITUDE'&
                        /10X, 'AND ANGLE MUST LIE BETWEEN 0 AND 1.'&
                       //10X, 1p, e10.2, '  TOLER1'&
                        /10X, 1p, e10.2, '  TOLER2')
            99008 format(/1X, a9, 'ERROR.  THE GRID IS NOT ORDERED.')
            99009 format(/1X, a9, 'ERROR.  SOME INTERVALS IN THE GRID ARE TOO SHORT.'&
                        /10X, 'THE NEW GRID WOULD NOT BE ORDERED.')

      stop

      return
      end subroutine refine
