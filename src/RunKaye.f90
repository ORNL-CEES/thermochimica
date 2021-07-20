subroutine runKaye

    USE ModuleThermoIO

    implicit none

    ! Specify units:
    call SetStandardUnits()

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 400D0
    dElementMass(42)       = 0.9D0            ! Mo
    dElementMass(44)       = 0.1D0            ! Ru

    ! Specify output mode:
    iPrintResultsMode     = 2

    ! Parse the ChemSage data-file:
    call ParseKayeDataFile()

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    ! Perform post-processing of results:
    if (iPrintResultsMode > 0)  call PrintResults

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end subroutine runKaye
