
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo30.F90
    !> \brief   Spot test - W-Au-Ar-O.
    !> \author  M.H.A. Piro, B.W.N. Fitzpatrick
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer           Description of change
    !    ----          ----------           ---------------------
    !    05/14/2013    M.H.A. Piro          Original code
    !    31/08/2018    B.W.N. Fitzpatrick   Change system to fictive RKMP model
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a fictive system labelled W-Au-Ar-O.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo30

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'W-Au-Ar-O_01.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1455D0
    dElementMass(74)       = 1.95D0        ! W
    dElementMass(79)       = 1D0           ! Au
    dElementMass(18)       = 2D0           ! Ar
    dElementMass(8)        = 10D0          ! O

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        if ((DABS(dGibbsEnergySys - (-4.620D5))/((-4.620D5))) < 1D-3) then
            ! The test passed:
            print *, 'TestThermo30: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestThermo30: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestThermo30: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    !call ThermoDebug

    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo30
