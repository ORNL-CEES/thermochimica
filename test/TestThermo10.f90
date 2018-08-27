
program ThermoTest10

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
    ! The purpose of this unit test is to ensure that Thermochimica does not proceed if the pressure of the 
    ! system is a NAN.  Note that a dummy variable is added to this test. 
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO

    implicit none

    real(8):: dTemp

    ! Initialize variables:
    dTemperature            = 300D0
    dPressure               = 0D0
    dElementMass            = 1D0
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/C-O.dat'
    dTemp                   = 0D0

    ! Make dPressure a NAN:
    dPressure = dPressure / dTemp

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)
        
    ! Call Thermochimica:
    call Thermochimica
                
    if (INFOThermo == 2) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo10: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo10: FAIL <---'
    end if
        
    ! Reset Thermochimica:
    call ResetThermo
    
end program ThermoTest10
