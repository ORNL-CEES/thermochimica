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
    dTemperature          = 1800D0
    dPressure             = 1.0D0
    dElementMass          = 0D0

    dElementMass(9)       = 3.175D0
    dElementMass(92)      = 0.05D0
    dElementMass(94)      = 0D0
    dElementMass(90)      = 0D0
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
        ((DABS(dMolFraction(1)  - 3.0653D-04)/3.0653D-04) < 1D-3) .AND. &
        ((DABS(dMolFraction(2)  - 0.11258D00)/0.11258D00) < 1D-3) .AND. &
        ((DABS(dMolFraction(3)  - 7.8118D-01)/7.8118D-01) < 1D-3) .AND. &
        ((DABS(dMolFraction(4)  - 9.3920D-02)/9.3920D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(5)  - 5.1978D-13)/5.1978D-13) < 1D-3) .AND. &
        ((DABS(dMolFraction(6)  - 1.2012D-02)/1.2012D-02) < 1D-3) .AND. &
        
        ((DABS(dMolFraction(7)  - 0.32446D00)/0.32446D00) < 1D-3) .AND. &
        ((DABS(dMolFraction(8)  - 6.8356D-02)/6.8356D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(9)  - 2.4117D-04)/2.4117D-04) < 1D-3) .AND. &
        ((DABS(dMolFraction(10) - 1.1394D-03)/1.1394D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(11) - 0.51831D00)/0.51831D00) < 1D-3) .AND. &
        ((DABS(dMolFraction(12) - 3.0169D-02)/3.0169D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(13) - 3.8912D-03)/3.8912D-03) < 1D-3) .AND. &
        ((DABS(dMolFraction(14) - 3.7166D-02)/3.7166D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(15) - 1.5083D-02)/1.5083D-02) < 1D-3) .AND. &
        ((DABS(dMolFraction(16) - 1.1878D-03)/1.1878D-03) < 1D-3) .AND. &
        ((DABS(dGibbsEnergySys - (-8.47555D6))/(-8.47555D6)) < 1D-3))  then
            ! The test passed:
            print *, 'TestMSTDB03: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestMSTDB03: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestMSTDB03: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    call ResetThermoAll

end program mstdb
