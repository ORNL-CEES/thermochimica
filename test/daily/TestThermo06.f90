
program ThermoTest06

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
    ! The purpose of this unit test is to ensure that Thermochimica does not proceed when the number of elements
    ! in the system is below a minimum value (2).  Note that the data-file represents the uranium and oxygen
    ! system, but the input specifies 0 mol of oxygen and 1 mol of uranium.
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO

    implicit none

    integer :: STATUS = 0

    ! Initialize variables:
    dTemperature            = 300D0
    dPressure               = 1D0
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/Example3.dat'

    dElementMass            = 0D0   ! All elements
    dElementMass(92)        = 1D0   ! Uranium

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if (INFOThermo == 5) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo06: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo06: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program ThermoTest06
