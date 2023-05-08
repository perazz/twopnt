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
    use twopnt_test
    use iso_fortran_env, only: output_unit
    implicit none

    type(SwirlingFlow) :: problem

    ! Parameters

    character(len=16) :: REPORT
    integer, parameter :: TEXT = output_unit
    logical :: ERROR

    ! Initialize problem object
    call problem%new(ERROR,TEXT)
    print *, 'erro1r=',error
    if (ERROR) return

    !    ! PROLOGUE. OPEN FILES.
    !    !OPEN (FILE = 'twopnt.out', STATUS='UNKNOWN',FORM = 'FORMATTED', UNIT = TEXT)

    ! Call driver
    call problem%run(ERROR, TEXT, REPORT, problem%U)
    print *, 'error2=',error

    ! Print a summary
    call problem%summary(TEXT)


    ! close(TEXT)

    return

end program TWMAIN
