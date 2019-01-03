
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo33.F90
    !> \brief   Spot test - W-Au-Ar-Ne-O, 900K.
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
    !    08/31/2018    B.W.N. Fitzpatrick  Modification to use a fictive system
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a fictive system labelled W-Au-Ar-Ne-O at 900K.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo33

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'W-Au-Ar-Ne-O_04.dat'

    ! Specify values:
    dPressure              = 2D0
    dTemperature           = 900D0
    dElementMass(74)       = 20D0        ! W
    dElementMass(79)       = 2D0         ! Au
    dElementMass(18)       = 7D0         ! Ar
    dElementMass(8)        = 5D0         ! O
    dElementMass(10)       = 1D0         ! Ne

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        if ((DABS(dMolFraction(1) - 6.9028D-11)/6.903D-11 < 1D-3).AND. &
        (DABS(dMolFraction(2) - 0.5758)/0.5758 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (3.220D6))/(3.220D6) < 1D-3))  then
            ! The test passed:
            print *, 'TestThermo33: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestThermo33: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestThermo33: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    !call ThermoDebug

    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo33
