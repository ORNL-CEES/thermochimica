
program TestThermo05

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
    ! The purpose of this unit test is to ensure that Thermochimica does not proceed when the input units are
    ! specified but not in supported units.  The mass unit 'tonnes' is not supported by Thermochimica.
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO

    implicit none

    integer :: STATUS = 0

    ! Initialize variables:
    dTemperature            = 300D0
    dPressure               = 1D0
    dElementMass            = 1D0
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'tonnes'
    cThermoFileName         = '../data/Example1.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if (INFOThermo == 4) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo05: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo05: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo05
