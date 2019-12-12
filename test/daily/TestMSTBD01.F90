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

    dElementMass(9)       = 13.5D0
    dElementMass(92)      = 1D0
    dElementMass(94)      = 1D0
    dElementMass(90)      = 1D0
    dElementMass(4)       = 1D0
    dElementMass(3)       = 1D0

    ! Specify output and debug modes:
    iPrintResultsMode     = 2
    lDebugMode            = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        if (&
        ((DABS(dMolFraction(13) - 8.0590D-03)/8.0590D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(14) - 3.5118D-02)/3.5118D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(15) - 1.1534D-02)/1.1534D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(16) - 3.6175D-02)/3.6175D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(17) - 1.8828D-02)/1.8828D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(18) - 6.9123D-02)/6.9123D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(19) - 5.9242D-02)/5.9242D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(20) - 2.8986D-02)/2.8986D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(21) - 2.5245D-02)/2.5245D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(22) - 4.4147D-02)/4.4147D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(23) - 6.9456D-02)/6.9456D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(24) - 3.7824D-02)/3.7824D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(25) - 1.9373D-02)/1.9373D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(26) - 4.9538D-02)/4.9538D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(27) - 3.1646D-02)/3.1646D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(28) - 5.6092D-02)/5.6092D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(29) - 4.1913D-02)/4.1913D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(30) - 9.2999D-02)/9.2999D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(31) - 6.1080D-02)/6.1080D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(32) - 0.12900D0)/0.12900D0) < 1D-3) .AND. &
        ((DABS(dMolFraction(33) - 7.4626D-02)/7.4626D-02) < 1D-3) .AND. &
        ((DABS(dGibbsEnergySys - (-8.47555D6))/(-8.47555D6)) < 1D-3))  then
            ! The test passed:
            print *, 'TestMSTDB01: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestMSTDB01: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestMSTDB01: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    call ResetThermoAll

end program mstdb
