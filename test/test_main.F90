!-------------------------------------------------------------------------------------------------------------
!
!> \file    test_main.F90
!> \brief   Main test driver for Thermochimica using test-drive framework
!> \details Runs all test suites and reports results
!
!-------------------------------------------------------------------------------------------------------------

program test_main
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_error_handling, only : collect_error_handling
    use test_systems, only : collect_systems
    implicit none

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    ! Define all test suites
    testsuites = [ &
        new_testsuite("error_handling", collect_error_handling), &
        new_testsuite("systems", collect_systems) &
        ]

    ! Print header
    write(error_unit, '(A)') ""
    write(error_unit, '(A)') "=========================================="
    write(error_unit, '(A)') "  Thermochimica Test Suite (test-drive)"
    write(error_unit, '(A)') "=========================================="
    write(error_unit, '(A)') ""

    ! Run each test suite
    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
        write(error_unit, '(A)') ""
    end do

    ! Report final results
    write(error_unit, '(A)') "=========================================="
    if (stat > 0) then
        write(error_unit, '(A, I0, A)') "FAILED: ", stat, " test(s) failed!"
        write(error_unit, '(A)') "=========================================="
        error stop
    else
        write(error_unit, '(A)') "SUCCESS: All tests passed!"
        write(error_unit, '(A)') "=========================================="
    end if

end program test_main
