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
    character(len=16) :: reason
    integer, parameter :: text = output_unit
    logical :: error

    ! Initialize problem object
    call problem%new(error,text)
    if (error) then
        write(text,1)
        return
    endif

    ! Call driver
    call problem%run(error, text, reason, problem%U)

    ! Print a summary
    call problem%summary(text)

    1 format(/1X, a9, 'ERROR.  CANNOT INITIALIZE SWIRL FLOW PROBLEM.')

end program TWMAIN
