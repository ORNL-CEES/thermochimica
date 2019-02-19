
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo14.F90
    !> \brief   Unit test - maximum solution phases.
    !> \author  M. Poschmann
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    02/15/2019    M. Poschmann        Original
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this unit test is to ensure that Thermochimica can read a data file with many
    !! solution phases (42).
    !
    !-------------------------------------------------------------------------------------------------------------

program ThermoTest14

    USE ModuleThermoIO

    implicit none

    ! Initialize variables:
    cThermoFileName         = DATA_DIRECTORY // 'test14.dat'

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    if (INFOThermo == 0) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo14: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The unit test failed.
        print *, 'TestThermo14: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Reset Thermochimica:
    call ResetThermo

end program ThermoTest14
