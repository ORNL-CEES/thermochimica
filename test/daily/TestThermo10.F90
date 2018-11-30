
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo10.F90
    !> \brief   Unit test - pressure NAN.
    !> \author  M.H.A. Piro, B.W.N. Fitzpatrick
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
    !    02/07/2012    M.H.A. Piro         Original code
    !    28/11/2018    M. Poschmann        Created new test with C-O database
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this unit test is to ensure that Thermochimica does not proceed when the mass is
    !! out of range.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo10

    USE ModuleThermoIO

    implicit none

    ! Initialize variables:
    dTemperature            = 300D0
    dPressure               = -1D0
    dPressure               = sqrt(dPressure)
    dElementMass            = 1D0
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = DATA_DIRECTORY // 'C-O.dat'

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if (INFOThermo == 2) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo10: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The unit test failed.
        print *, 'TestThermo10: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo10
