
program TestThermo15

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
    !    02/09/2012    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this unit test is to test both the parser and Thermochimica for a system comprised of the
    ! maximum number of elements that can be handled by FactSage.
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
    cThermoFileName         = '../data/Example15.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if (INFOThermo == 0) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo15: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo15: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo15
