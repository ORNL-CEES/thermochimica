program TestThermo63

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    integer   :: i,j,k
    real(8) :: gibbsCheck, p1check, p2check, s1check
    logical :: subqPass, solPass

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY //'ClAlNa.dat'

    ! Specify values:
    dTemperature          = 1000D0
    dPressure             = 1.0D0
    dElementMass          = 0D0
    dElementMass(17)      = 3D0                              ! Cl
    dElementMass(11)      = 1D0                              ! Na
    dElementMass(13)      = 1D0                              ! Al

    gibbsCheck = -1.17685D+06
    p1check    = 4.2176D0
    s1check    = 0.23035D0
    p2check    = 4.4748D-2

    ! Specify output and debug modes:
    iPrintResultsMode     = 2
    lDebugMode            = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    subqPass = .FALSE.
    solPass = .FALSE.
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
                end if
            end do
            do i = 1, nConPhases
                j = iAssemblage(i)
                if (ADJUSTL(TRIM(cSpeciesName(j))) == 'NaCl_S1(s)') then
                    if (DABS((dMolesPhase(i)-p2check)/p2check) < 1D-3) then
                        solPass = .TRUE.
                    end if
                end if
            end do
        end if
    end if

    if (subqPass .AND. solPass) then
        ! The test passed:
        print *, 'TestThermo63: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo63: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo63
