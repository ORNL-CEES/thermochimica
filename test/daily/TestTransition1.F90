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
    dTempTolerance       = 1D-1

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'AlMg-Liang.dat'
    
    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1000
    dElementMass(13)       = 0.5 ! Al
    dElementMass(12)       = 0.5 ! Mg

    ! Init test values
    dTempMin              = 300
    dTempMax              = 1000
    iTestTransitions      = 4
    dTestTransitionTemp   = [523.24D0, 667.39D0, 729.97D0, 735.02D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0]
    lPass                 = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)
    call Thermochimica ! or I can init hardcode nElements

    !Test call
    call PhaseTransition(dTempMin, dTempMax, dTempTolerance, dPhaseTransitionTemp, iTransitions)
    call checkTransitionTest(dPhaseTransitionTemp, dTestTransitionTemp, 1D0, lPass)
    do i = 1, size(dPhaseTransitionTemp)
        print *, dPhaseTransitionTemp(i)
    end do 
    if (lPass) then
        ! The test passed:
        print *, 'TestTransition1: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestTransition1: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestTransition1