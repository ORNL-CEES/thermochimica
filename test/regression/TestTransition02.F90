    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestTransition02.F90
    !> \brief   Testing phase transition subroutine
    !> \author  A.E.F. Fitzsimmons
    !
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    05/26/2024    A.E.F. Fitzsimmons   Original Code
    !
    ! Purpose:
    ! ========
    !> \details Testing new Transition function
    !!
    !
    !-------------------------------------------------------------------------------------------------------------
program TestTransition02

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleTesting

    implicit none

    ! Init variables
    real(8) :: dTempMin, dTempMax, dTempTolerance
    real(8), dimension(10) :: dPhaseTransitionTemp, dTestTransitionTemp
    integer :: iTransitions, iTestTransitions
    logical :: lPass
    
    dTempTolerance       = 1D-1

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'MoPdRuRhTc-Kaye.dat'
    
    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1000
    dElementMass(44)       = 0.55 ! Ru
    dElementMass(42)       = 0.45 ! Mo

    ! Init test values
    dTempMin              = 1000
    dTempMax              = 2500
    iTestTransitions      = 4
    dTestTransitionTemp   = [1415.93D0, 1824.70D0, 2295.85D0, 2354.94D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0]
    lPass                 = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)
    call Thermochimica 

    !Test call
    call PhaseTransition(dTempMin, dTempMax, dTempTolerance, dPhaseTransitionTemp, iTransitions)
    call checkTransitionTest(dPhaseTransitionTemp, dTestTransitionTemp, 1D0, lPass)

    if (lPass) then
        ! The test passed:
        print *, 'TestTransition02: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestTransition02: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestTransition02