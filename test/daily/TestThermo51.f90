

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
    ! The purpose of this application test is to ensure that the Pu-U-Mo-Zr
    ! quaternary system is correctly computed.  This data-file was kindly
    ! provided by T.M. Besmann, which was converted from the TAF-ID database.
    ! At equilibrium for this particular temperature, pressure and composition,
    ! the BCC and C15 Laves phases should be stable.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo51

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Specify units:
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/UPMZ_MHP.dat'

    ! Specify values:
    dTemperature            = 1000D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(94)        = 1D0
    dElementMass(92)        = 1D0
    dElementMass(42)        = 1D0
    dElementMass(40)        = 1D0

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(14) - 0.21131D0)/0.21131D0 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-2.73704D5))/(-2.73704D5) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo51: PASS'
        else
            ! The test failed.
            print *, 'TestThermo51: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo51: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo51
