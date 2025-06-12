    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestTransition4.F90
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
program TestTransition4

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleTesting

    implicit none

    ! Init variables
    real(8) :: dTempMin, dTempMax, dTempTolerance
    real(8), dimension(10) :: dPhaseTransitionTemp, dTestTransitionTemp
    integer :: iTransitions, iTestTransitions
    logical :: lPass

    dTempTolerance       = 1D0

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'ZrH-Dupin.dat'
    
    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1000
    dElementMass(40)       = 0.97 ! Zr
    dElementMass(1)        = 0.03 ! H

    ! Init test values
    dTempMin              = 300
    dTempMax              = 1300
    iTestTransitions      = 3
    dTestTransitionTemp   = [735.39D0, 1064.22D0, 1103.52D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0]
    lPass                 = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)
    call Thermochimica 

    ! Test call
    call PhaseTransition(dTempMin, dTempMax, dTempTolerance, dPhaseTransitionTemp, iTransitions)
    call checkTransitionTest(dPhaseTransitionTemp, dTestTransitionTemp, 1D0, lPass)

    if (lPass) then
        ! The test passed:
        print *, 'TestTransition4: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestTransition4: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestTransition4