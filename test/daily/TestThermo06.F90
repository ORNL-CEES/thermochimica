
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo06.F90
    !> \brief   Unit test - no mass specified.
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
    !    05/07/2018    B.W.N. Fitzpatrick  Modified purpose of code (to see if unspecified mass exists gracefully)
    !
    ! Purpose:
    ! ========
    !> \details To verify that Thermochimica gracefully exits when mass is not specified.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo06

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    ! Initialize variables:
    dTemperature            = 2000D0
    dPressure               = 1D0
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = DATA_DIRECTORY // 'C-O.dat'

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if (INFOThermo == 5) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo06: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The unit test failed.
        print *, 'TestThermo06: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo06
