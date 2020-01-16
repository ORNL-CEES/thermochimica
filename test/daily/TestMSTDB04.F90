program mstdb

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleParseCS

    implicit none

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'LiBeUPuThF_T1.dat'


    ! Specify values:
    dTemperature          = 1500D0
    dPressure             = 1.0D0
    dElementMass          = 0D0

    dElementMass(9)       = 3.75D0
    dElementMass(92)      = 1D0

    ! Specify output and debug modes:
    iPrintResultsMode     = 0
    lDebugMode            = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    if (iPrintResultsMode > 0)  call PrintResults

    ! Check results:
    if (INFOThermo == 0) then
        if (&
        ((DABS(dMolFraction(3) - 0.55153)/0.55153) < 1D-3) .AND. &
        ((DABS(dMolFraction(4) - 5.1534D-2)/5.1534D-2) < 1D-3) .AND. &
        ((DABS(dMolFraction(5) - 0.39693)/0.39693) < 1D-3) .AND. &
        ((DABS(dGibbsEnergySys - (-2.18912D6))/(-2.18912D6)) < 1D-3))  then
            ! The test passed:
            print *, 'TestMSTDB04: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestMSTDB04: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestMSTDB04: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    call ResetThermoAll

end program mstdb
