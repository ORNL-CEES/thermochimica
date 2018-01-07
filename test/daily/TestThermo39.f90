

    !-------------------------------------------------------------------------------------------------------------
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer        Description of change
    !    ----          ----------        ---------------------
    !    09/06/2015    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that the U-Zr binary system
    ! is correctly computed.  This data-file was kindly provided by T.M. Besmann.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo39

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName       = '../data/UPMZ_MHP.dat'

    ! Specify values:
    dTemperature            = 942D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(92)        = 0.628D0
    dElementMass(40)        = 0.372D0

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(21) - 4.8495D-3)/4.8495D-3 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-5.79316D4))/(-5.79316D4) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo39: PASS'
        else
            ! The test failed.
            print *, 'TestThermo39: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo39: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo39
