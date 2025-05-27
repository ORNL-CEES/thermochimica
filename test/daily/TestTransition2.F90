program TestTransition1

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleTesting

    implicit none

    ! Init variables
    real(8) :: dTempMin, dTempMax, dTempTolerance
    real(8), dimension(10) :: dPhaseTransitionTemp, dTestTransitionTemp
    integer :: i, iTransitions, iTestTransitions
    logical :: lPass

    i = 1
    iTransitions = 0
    dPhaseTransitionTemp = 0D0
    dTempTolerance       = 1D0

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'MoPdRuRhTc-Kaye.dat'
    
    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1000
    dElementMass(44)       = 0.55 ! Mo
    dElementMass(42)       = 0.45 ! Ru

    ! Init test values
    dTempMin              = 1000
    dTempMax              = 2500
    iTestTransitions      = 4
    dTestTransitionTemp   = [1415.93D0, 1824.70D0, 2295.85D0, 2354.94D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0]
    lPass                 = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)
    call Thermochimica ! or I can init hardcode nElements

    !Test call
    call PhaseTransition(dTempMin, dTempMax, dTempTolerance, dPhaseTransitionTemp, iTransitions)
    call checkTransitionTest(dPhaseTransitionTemp, dTestTransitionTemp, 1D0, lPass)

    if (lPass) then
        ! The test passed:
        print *, 'TestTransition2: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestTransition2: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestTransition1