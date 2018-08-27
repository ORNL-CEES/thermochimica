
program ThermoTest12

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
    ! The purpose of this unit test is to ensure that Thermochimica does not proceed if the number of chemical  
    ! elements is out of range.
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO

    implicit none

    real(8):: dTemp

    ! Initialize variables:
    dTemperature            = 300D0
    dPressure               = 1D0
    dElementMass            = 1D0
    ! Add one to few chemical elements
    dElementMass(6)         = 0



    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/C-O.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)
        
    ! Call Thermochimica:
    call Thermochimica
                
    if (INFOThermo == 5) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo12: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo12: FAIL <---', INFOThermo
    end if
        
    ! Reset Thermochimica:
    call ResetThermo
    
end program ThermoTest12
