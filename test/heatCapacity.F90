program heatCapacity

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    integer   :: nConPhasesIn, nSolnPhasesIn, i
    integer, allocatable   :: iAssemblageIn(:)
    real(8), allocatable   :: dMolesPhaseIn(:), dMolFractionIn(:)
    real(8)   :: dGibbsEnergySysOut, dHeatCapacity

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'Pu_U_O_gasOnly.dat'

    ! Specify values:
    dTemperature          = 4100D0
    dPressure             = 1.0D0
    dElementMass          = 0D0
    dElementMass(8)       = 2.01D0                             ! O
    dElementMass(92)      = 1.0D0                              ! U

    ! Specify output and debug modes:
    iPrintResultsMode     = 2
    lDebugMode            = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    ! Perform post-processing of results:
    if (iPrintResultsMode > 0)  call PrintResults

    dHeatCapacity = 0D0;
    do i = 1, nSpecies
      print *, i, cSpeciesName(i), dStdHeatCapacity(i), dMolesSpecies(i)
      dHeatCapacity = dHeatCapacity + dStdHeatCapacity(i) * dMolesSpecies(i)
    end do

    ! dHeatCapacity = dHeatCapacity

    print *, dHeatCapacity

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program heatCapacity
