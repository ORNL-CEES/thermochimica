
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo32.F90
    !> \brief   Spot test - W-Au-Ar-Ne-O, 2452K.
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
    !! results for a fictive system labelled W-Au-Ar-Ne-O at 2452K.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo32

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'W-Au-Ar-Ne-O_03.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2452D0
    dElementMass(74)       = 1.95D0        ! W
    dElementMass(79)       = 1D0           ! Au
    dElementMass(18)       = 2D0           ! Ar
    dElementMass(8)        = 10D0          ! O
    dElementMass(10)       = 10D0          ! Ne

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    !iPrintResultsMode = 2
    !call PrintResults

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dGibbsEnergySys - (1.672D7))/((1.672D7))) < 1D-3) then
            ! The test passed:
            print *, 'TestThermo32: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestThermo32: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestThermo32: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo32
