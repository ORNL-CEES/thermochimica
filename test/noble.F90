program noble

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    integer   :: nConPhasesIn, nSolnPhasesIn
    integer, allocatable   :: iAssemblageIn(:)
    real(8), allocatable   :: dMolesPhaseIn(:), dMolFractionIn(:)
    real(8)   :: dGibbsEnergySysOut

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'Kaye_NobleMetals.dat'

    ! Specify values:
    dTemperature          = 1650D0
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

    ! Test GibbsEnergy

    allocate(iAssemblageIn(nElements),dMolesPhaseIn(nElements),dMolFractionIn(nSpecies))

    nConPhasesIn = nConPhases
    nSolnPhasesIn = nSolnPhases
    iAssemblageIn = iAssemblage
    dMolesPhaseIn = dMolesPhase
    dMolFractionIn = dMolFraction

    call GibbsEnergy(nConPhasesIn, nSolnPhasesIn, iAssemblageIn, dMolesPhaseIn, dMolFractionIn, dGibbsEnergySysOut)

    call ResetThermo

    dTemperature          = 300D0
    dPressure             = 1.0D0
    dElementMass          = 0D0
    dElementMass(46)      = 0.9D0                              ! Pd
    dElementMass(44)      = 0.1D0                              ! Ru
    dElementMass(43)      = 0.0D0                              ! Tc
    dElementMass(42)      = 1.9D0                              ! Mo

    if (INFOThermo == 0)        call Thermochimica

    ! Perform post-processing of results:
    if (iPrintResultsMode > 0)  call PrintResults

    call GibbsEnergy(nConPhasesIn, nSolnPhasesIn, iAssemblageIn, dMolesPhaseIn, dMolFractionIn, dGibbsEnergySysOut)

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program noble
