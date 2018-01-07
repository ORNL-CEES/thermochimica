
program ThermoTest07

    !-------------------------------------------------------------------------------------------------------------
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer        Description of change
    !    ----          ----------        ---------------------
    !    02/07/2012    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this unit test is to ensure that Thermochimica does not proceed if the pressure of the
    ! system is a NAN.  Note that a dummy variable is added to only this test.
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO

    implicit none

    integer :: STATUS = 0

    real(8):: dTemp

    ! Initialize variables:
    dTemperature            = 300D0
    dPressure               = 0D0
    dElementMass            = 1D0
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/Example1.dat'
    dTemp                   = 0D0

    ! Make dPressure a NAN:
    dPressure = dPressure / dTemp

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if (INFOThermo == 2) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo07: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo07: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program ThermoTest07
