program heatCapacityLiF

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'gasOnlyExample.dat'

    ! Specify values:
    dTemperature          = 1950D0
    dPressure             = 1.0D0
    dElementMass          = 0D0
    dElementMass(3)       = 1D0                             ! Li
    dElementMass(4)       = 1D0                             ! Be
    dElementMass(9)       = 3D0                             ! F

    ! Specify output and debug modes:
    iPrintResultsMode     = 2
    lDebugMode            = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    if (iPrintResultsMode > 0)  call PrintResults
    call HeatCapacity
    print *, 'Heat Capacity [J/K]: ', dHeatCapacity

    call ResetThermoAll

end program heatCapacityLiF
