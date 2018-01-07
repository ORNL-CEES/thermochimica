

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
    ! the BCC phase should be stable.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo54

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
    dTemperature            = 301D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(94)        = 0.49631070492427776D0
    dElementMass(92)        = 0.76628926756619464D0
    dElementMass(42)        = 0.11009266551580255D0
    dElementMass(40)        = 0.72518063921759246D0

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(69) - 0.41623D0)/0.41623D0 < 1D-3).AND. &
            (DABS(dMolFraction(42) - 0.3147D0)/0.3147D0 < 1D-3).AND. &
            (DABS(dGibbsEnergySys - (-3.97487D4))/(-3.97487D4) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo54: PASS'
        else
            ! The test failed.
            print *, 'TestThermo54: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo54: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo54
