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
    dTemperature          = 800D0
    dPressure             = 1.0D0
    dElementMass          = 0D0

    dElementMass(9)       = 13.5D0
    dElementMass(92)      = 1D0
    dElementMass(94)      = 1D0
    dElementMass(90)      = 1D0
    dElementMass(4)       = 1D0
    dElementMass(3)       = 1D0

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
        ((DABS(dMolFraction(13) - 3.1291D-02)/3.1291D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(14) - 0.29198D00)/0.29198D00) < 1D-3) .AND. &
        ((DABS(dMolFraction(15) - 1.1543D-04)/1.1543D-04) < 1D-3) .AND. &
        ((DABS(dMolFraction(16) - 4.1756D-03)/4.1756D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(17) - 9.8813D-03)/9.8813D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(18) - 6.8670D-04)/6.8670D-04) < 1D-3) .AND. &
        ((DABS(dMolFraction(19) - 0.35729D00)/0.35729D00) < 1D-3) .AND. &
        ((DABS(dMolFraction(20) - 5.5945D-03)/5.5945D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(21) - 8.1356D-03)/8.1356D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(22) - 2.3646D-02)/2.3646D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(23) - 6.2663D-02)/6.2663D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(24) - 1.1940D-03)/1.1940D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(25) - 1.5495D-02)/1.5495D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(26) - 0.12382D00)/0.12382D00) < 1D-3) .AND. &
        ((DABS(dMolFraction(27) - 2.9097D-03)/2.9097D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(28) - 1.2802D-02)/1.2802D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(29) - 5.9238D-03)/5.9238D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(30) - 3.2676D-02)/3.2676D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(31) - 6.3078D-04)/6.3078D-04) < 1D-3) .AND. &
        ((DABS(dMolFraction(32) - 3.4841D-03)/3.4841D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(33) - 5.6062D-03)/5.6062D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(34) - 0.65393D00)/0.65393D00) < 1D-3) .AND. &
        ((DABS(dMolFraction(35) - 0.34607D00)/0.34607D00) < 1D-3) .AND. &
        ((DABS(dElementPotential(1)*dTemperature*dIdealConstant - &
        (-3.604762D05))/(-3.604762D05)) < 1D-3) .AND. &
        ((DABS(dElementPotential(1)*dTemperature*dIdealConstant - &
        (-2.783790D05))/(-2.783790D05)) < 1D-3) .AND. &
        ((DABS(dElementPotential(1)*dTemperature*dIdealConstant - &
        (-4.434106D05))/(-4.434106D05)) < 1D-3) .AND. &
        ((DABS(dElementPotential(1)*dTemperature*dIdealConstant - &
        (-4.517702D05))/(-4.517702D05)) < 1D-3) .AND. &
        ((DABS(dElementPotential(1)*dTemperature*dIdealConstant - &
        (-1.873157D05))/(-1.873157D05)) < 1D-3) .AND. &
        ((DABS(dElementPotential(1)*dTemperature*dIdealConstant - &
        (-2.225840D05))/(-2.225840D05)) < 1D-3) .AND. &
        ((DABS(dGibbsEnergySys - (-7.59106D06))/(-7.59106D06)) < 1D-3))  then
            ! The test passed:
            print *, 'TestMSTDB02: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestMSTDB02: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestMSTDB02: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    call ResetThermoAll

end program mstdb
