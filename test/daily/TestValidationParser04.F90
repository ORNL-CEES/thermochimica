    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ValidationParserTest01.F90
    !> \brief   Testing validation parser
    !> \author  A.E.F. Fitzsimmons
    !
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    07/15/2025    A.E.F. Fitzsimmons   Original Code
    !
    ! Purpose:
    ! ========
    !> \details Debug call test to many elements.
    !!
    !
    !-------------------------------------------------------------------------------------------------------------
program TestValidationParser04
    USE ModuleTesting
    USE ModuleThermoIO

    implicit none
    real(8), dimension(15) :: dPhaseTransitionTemp
    logical :: lPass

    lPass = .FALSE.
    dPhaseTransitionTemp = 0D0
    nMaxPhaseTransition = SIZE(dPhaseTransitionTemp)

    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    dPressure             = 1D0
    dTemperature          = 1000D0
    cThermoFileName       = DATA_DIRECTORY // "MoPdRuRhTc-Kaye.dat" 

    call ParseCSDataFile(cThermoFileName)
    call ParseValidation(DATA_DIRECTORY // 'IO_CSV-01.csv', lPass)
    
    if (lPass) then
        ! The test passed:
        print *, 'TestValidationParser04: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestValidationParser04: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestValidationParser04
