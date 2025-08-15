    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestValidation01.F90
    !> \brief   Testing validation parser
    !> \author  A.E.F. Fitzsimmons
    !
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    07/30/2025    A.E.F. Fitzsimmons   Original Code
    !
    ! Purpose:
    ! ========
    !> \details Validation of LiF-CsF system. Basked on Lipkina Project
    !!
    !
    !-------------------------------------------------------------------------------------------------------------
program TestValidation01
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
    cThermoFileName       = DATA_DIRECTORY // "Lipkina_LiFCsF.dat" 

    call ParseCSDataFile(cThermoFileName)
    call ParseCSVFile(CSV_DIRECTORY // 'MSTDBTC_LiF-CsF.csv', lPass)
    
    !call ThermoDebug

    ! Test for wrong number of elements
    if (lPass) then
        ! The test passed:
        print *, 'TestValidation01: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestValidation01: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestValidation01
