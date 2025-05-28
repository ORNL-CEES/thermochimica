    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestTransition6.F90
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
program TestTransition6

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
    cThermoFileName        = DATA_DIRECTORY // 'CsTe-1.dat'
    
    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1000
    dElementMass(52)       = 0.25 ! Te
    dElementMass(55)       = 0.75 ! Cs

    ! Init test values
    dTempMin              = 300
    dTempMax              = 1500
    iTestTransitions      = 3
    dTestTransitionTemp   = [301.33D0, 850.82D0, 1027.51D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0]
    lPass                 = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)
    call Thermochimica

    ! Test call
    call PhaseTransition(dTempMin, dTempMax, dTempTolerance, dPhaseTransitionTemp, iTransitions)
    call checkTransitionTest(dPhaseTransitionTemp, dTestTransitionTemp, 1D0, lPass)


    if (lPass) then
        ! The test passed:
        print *, 'TestTransition6: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestTransition6: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestTransition6