

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
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    02/07/2012    M.H.A. Piro         Original code
    !    11/07/2018    B.W.N. Fitzpatrick  Changed to a C-O database
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this unit test is to ensure that Thermochimica does not proceed when the input units are
    ! not specified.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo03

    USE ModuleThermoIO

    implicit none

    ! Initialize variables:
    dTemperature            = 300D0
    dPressure               = 1D0
    dElementMass            = 1D0
    !cThermoInputUnits(1)    = 'K'
    !cThermoInputUnits(2)    = 'atm'
    !cThermoInputUnits(3)    = 'moles'
    cThermoFileName         = DATA_DIRECTORY // 'C-O.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if (INFOThermo == 4) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo03: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo03: FAIL <---'
    end if

    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo03
