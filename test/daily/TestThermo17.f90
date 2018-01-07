
program TestThermo17

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
    !    02/10/2012    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The sole purpose of this unit test is to test the parser for a data-file that containing the maximum number
    ! of species supported by FactSage.
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Initialize variables:
    dTemperature            = 1000D0
    dPressure               = 1D0
    dElementMass            = 1D0
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/Example14.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if (INFOThermo == 17) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo17: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo17: FAIL <---', INFOThermo
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo17
