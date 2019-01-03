
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo43.F90
    !> \brief   Spot test - 400K with 40% Pd, 60% Ru.
    !> \author  M.H.A. Piro, B.W.N. Fitzpatrick
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    05/14/2013    M.H.A. Piro         Original code
    !    08/31/2018    B.W.N. Fitzpatrick  Modification to use Kaye's Pd-Ru-Tc-Mo system
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for the Pd-Ru-Tc-Mo system at 400K with 40% Pd, 60% Ru.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo43

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'Kaye_NobleMetals.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 400D0
    dElementMass(46)       = 0.4D0        ! Pd
    dElementMass(44)       = 0.6D0        ! Ru

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        if (((DABS(dMolFraction(3) - 0.9945D0)/0.9945D0) < 1D-3).AND. &
        ((DABS(dMolFraction(4) - 5.4695D-3)/5.4695D-3) < 1D-3).AND. &
        ((DABS(dGibbsEnergySys - (-1.338D4))/(1.338D4)) < 1D-3))  then
            ! The test passed:
            print *, 'TestThermo43: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestThermo43: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestThermo43: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo43
