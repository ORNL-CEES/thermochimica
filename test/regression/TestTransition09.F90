    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestTransition09.F90
    !> \brief   Testing validation parser
    !> \author  A.E.F. Fitzsimmons
    !
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    07/14/2025    A.E.F. Fitzsimmons   Original Code
    !
    ! Purpose:
    ! ========
    !> \details Testing new validation parser, using verification measurements.
    !!
    !
    !-------------------------------------------------------------------------------------------------------------
program TestTransition09
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
    call ParseCSVFile(CSV_DIRECTORY // 'Kaye_Pd-Tc.csv', lPass)
    
    if (lPass) then
        ! The test passed:
        print *, 'TestTransition09: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestTransition09: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestTransition09
