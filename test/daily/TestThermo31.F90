
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo31.F90
    !> \brief   Spot test - W-Au-Ar-O_02.
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
    !    08/31/2018    B.W.N. Fitzpatrick   Change to a fictive database
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a fictive W-Au-Ar-O_02 system.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo31

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'W-Au-Ar-O_02.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1000D0
    dElementMass(74)       = 1D0        ! W
    dElementMass(79)       = 3D0        ! Au
    dElementMass(18)       = 5D0        ! Ar
    dElementMass(8)        = 2D0        ! O

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    !iPrintResultsMode = 2
    !call PrintResults

    ! Check results:
    if (INFOThermo == 0) then
        if ((DABS(dGibbsEnergySys - (6.769E5))/((6.769D5))) < 1D-3) then
            ! The test passed:
            print *, 'TestThermo31: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestThermo31: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestThermo31: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    !call ThermoDebug

    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo31
