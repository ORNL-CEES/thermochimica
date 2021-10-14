program TestThermo62

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    integer   :: i,j,k
    real(8) :: gibbsCheck, p1check, p2check, s1check
    logical :: subqPass, gasPass

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY //'NaCl-AlCl3.dat'

    ! Specify values:
    dTemperature          = 2000D0
    dPressure             = 1.0D0
    dElementMass          = 0D0
    dElementMass(17)      = 2D0                              ! Cl
    dElementMass(13)      = 1D0                              ! Al

    gibbsCheck = -9.64834D+05
    p1check    = 0.80307D0
    s1check    = 3.8514D-03
    p2check    = 0.60797D0

    ! Specify output and debug modes:
    iPrintResultsMode     = 2
    lDebugMode            = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    subqPass = .FALSE.
    gasPass = .FALSE.
    if (INFOThermo == 0) then
        if (DABS((dGibbsEnergySys - gibbsCheck)/gibbsCheck) < 1D-3) then
            do i = 1, nSolnPhases
                k = nElements + 1 - i
                j = -iAssemblage(k)
                if (cSolnPhaseName(j) == 'MSsoln') then
                    if (DABS((dMolesPhase(k)-p1check)/p1check) < 1D-3) then
                        if (DABS((dMolFraction(nSpeciesPhase(j-1)+1)-s1check)/s1check) < 1D-3) then
                            subqPass = .TRUE.
                        end if
                    end if
                else if (cSolnPhaseName(j) == 'gas_ideal') then
                    if (DABS((dMolesPhase(k)-p2check)/p2check) < 1D-3) then
                        gasPass = .TRUE.
                    end if
                end if
            end do
        end if
    end if

    if (subqPass .AND. gasPass) then
        ! The test passed:
        print *, 'TestThermo62: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo62: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo62
