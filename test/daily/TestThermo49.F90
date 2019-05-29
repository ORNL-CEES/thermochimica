
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo49.F90
    !> \brief   Spot test - 1000K with 40% Pd, 60% Tc.
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
    !! results for the Pd-Ru-Tc-Mo system at 1000K with 40% Pd, 60% Tc.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo49

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
    dTemperature           = 1900D0
    dElementMass           = 0D0
    dElementMass(43)       = 0.125D0        ! Tc
    dElementMass(46)       = 0.874D0        ! Pd

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        if (((DABS(dMolFraction(15) - 0.90358D0)/0.90358D0) < 1D-3).AND. &
        ((DABS(dMolFraction(4) - 0.18541D0)/0.18541D0) < 1D-3).AND. &
        ((DABS(dGibbsEnergySys - (-1.28092D5))/(-1.28092D5)) < 1D-3))  then
            ! The test passed:
            print *, 'TestThermo49: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestThermo49: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestThermo49: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

! Reset Thermochimica:
call ResetThermo

end program TestThermo49
