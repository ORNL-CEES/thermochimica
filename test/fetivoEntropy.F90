program fetivoEntropy

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS

    implicit none

    integer   ::  nConPhasesIn, nSolnPhasesIn
    integer, allocatable   :: iAssemblageIn(:)
    real(8), allocatable   ::  dMolesPhaseIn(:), dMolFractionIn(:)
    real(8)   :: dGibbsEnergySysOut

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'FeTiVO-noMix.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2000D0
    dElementMass(8)        = 1.35D0         ! O
    dElementMass(26)       = 1D0            ! Fe

    ! Specify output mode:
    iPrintResultsMode     = 2

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    ! Perform post-processing of results:
    if (iPrintResultsMode > 0)  call PrintResults


    allocate(iAssemblageIn(nElements),dMolesPhaseIn(nElements),dMolFractionIn(nSpecies))

    nConPhasesIn = nConPhases
    nSolnPhasesIn = nSolnPhases
    iAssemblageIn = iAssemblage
    dMolesPhaseIn = dMolesPhase
    dMolFractionIn = dMolFraction

    call GibbsEnergy(nConPhasesIn, nSolnPhasesIn, iAssemblageIn, dMolesPhaseIn,dMolFractionIn,dGibbsEnergySysOut)

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program fetivoEntropy
