

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
    !    09/13/2015    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that the Pu-Mo binary system
    ! is correctly computed.  This data-file was kindly provided by T.M. Besmann.
    ! This particular application test considers the Pu-Mo system at a particular
    ! temperature and composition that should correspond to only the liquid and
    ! BCC phases.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo50

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
    dTemperature            = 915D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(94)        = 0.99D0     ! Pu
    dElementMass(42)        = 0.01D0     ! Mo

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(5) - 0.99D0)/0.99D0 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-6.94934D4))/(-6.94934D4) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo50: PASS'
        else
            ! The test failed.
            print *, 'TestThermo50: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo50: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo50
