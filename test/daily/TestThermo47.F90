
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo47.F90
    !> \brief   Spot test - 1973K with 30% Mo, 40% Pd, 30% Ru.
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
    !! results for the Pd-Ru-Tc-Mo system at 1973K with 30% Mo, 40% Pd, 30% Ru.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo47

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
    dTemperature           = 1973D0
    dElementMass(42)       = 0.3D0        ! Mo
    dElementMass(46)       = 0.4D0        ! Pd
    dElementMass(44)       = 0.3D0        ! Ru

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        if (((DABS(dMolFraction(5) - 0.2475622D0)/0.2475622D0) < 1D-3).AND. &
        ((DABS(dMolFraction(6) - 0.6205065D0)/0.6205065D0) < 1D-3).AND. &
        ((DABS(dMolFraction(7) - 0.131931D0)/0.131931D0) < 1D-3).AND. &
        ((DABS(dGibbsEnergySys - (-1.38528D5))/(-1.38528D5)) < 1D-3))  then
            ! The test passed:
            print *, 'TestThermo47: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestThermo47: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestThermo47: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

! Reset Thermochimica:
call ResetThermo

end program TestThermo47
