
program TestThermo18

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
    !    02/11/2012    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this unit test is to test the parser for a data-file that contains phases with magnetic
    ! contributions to the standard molar Gibbs energy terms and convergence.
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Initialize variables:
    dTemperature            = 1000D0
    dPressure               = 1D0
    dElementMass(12)        = 1D0
    dElementMass(29)        = 6D0
    dElementMass(8)         = 3D0
    dElementMass(26)        = 3D0
    dElementMass(6)         = 1D0

    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/Example16.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if ((INFOThermo == 0).AND.(dMolFraction(9) > 0.61D0).AND.(dMolFraction(9) < 0.63D0)) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo18: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo18: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo18
