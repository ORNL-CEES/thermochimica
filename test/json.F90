program json

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'Kaye_NobleMetals.dat'

    ! Specify values:
    dTemperature          = 300D0
    dPressure             = 1.0D0
    dElementMass          = 0D0
    dElementMass(46)      = 0.9D0                              ! Pd
    dElementMass(44)      = 0.1D0                              ! Ru
    dElementMass(43)      = 0.0D0                              ! Tc
    dElementMass(42)      = 1.9D0                              ! Mo

    ! Specify output and debug modes:
    iPrintResultsMode     = 2
    lDebugMode            = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    ! Perform post-processing of results:
    if (iPrintResultsMode > 0)  call PrintResults

    call WriteJSON

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program json
