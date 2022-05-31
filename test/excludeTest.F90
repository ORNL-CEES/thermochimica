program exclude

    USE ModuleThermoIO

    implicit none


    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'Kaye_NobleMetals.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 400D0
    dElementMass(42)       = 0.9D0            ! Mo
    dElementMass(44)       = 0.1D0            ! Ru

    ! Specify output mode:
    iPrintResultsMode     = 2

    ! Exclude some phases
    nPhasesExcluded = 2
    cPhasesExcluded(1) = 'BCCN'
    cPhasesExcluded(2) = 'Pd11Mo9_s1(s)'
    ! nPhasesExcludedExcept = 0
    ! cPhasesExcludedExcept(1) = 'BCCN'
    ! cPhasesExcludedExcept(2) = 'HCPN'

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    ! Perform post-processing of results:
    if (iPrintResultsMode > 0)  call PrintResults

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program exclude
