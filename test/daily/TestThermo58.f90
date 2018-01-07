

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
    !    09/19/2015    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that Thermochimica computes the correct results for
    ! a specific point on the Cs-I phase diagram.  This particular test should yield a very narrow miscibility
    ! gap (Liquid-liquid).
    !
    ! The data-file has been kindly provided by J.-C. Dumas (CEA).  At the time that this test was created,
    ! a publication was not yet available.  Once a publicaation has been made in the open-literature, I will
    ! be sure to provide an appropriate reference.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo58

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Specify units:
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/CsI_RKMP_MHP.dat'

    ! Initialize variables:
    dTemperature            = 650D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(55)        = 0.89D0
    dElementMass(53)        = 0.11D0

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if (((DABS(dMolFraction(8) - 1.5854D-2)/1.5854D-2 < 1D-3).OR. &
            (DABS(dMolFraction(11) - 1.5854D-2)/1.5854D-2 < 1D-3)).AND. &
            (DABS(dGibbsEnergySys - (-9.71654D4))/(-9.71654D4) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo58: PASS'
        else
            ! The test failed.
            print *, 'TestThermo58: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo58: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo58
