

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
    ! This particular application test considers the U-Mo system at a particular
    ! temperature and composition that should correspond to a bcc-bcc
    ! miscibility gap.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo41

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
    dTemperature            = 1500D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(92)        = 0.5D0     ! U
    dElementMass(42)        = 0.5D0     ! Mo

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(6) - 0.42311D0)/0.42311D0 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-1.03462D5))/(-1.03462D5) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo41: PASS'
        else
            ! The test failed.
            print *, 'TestThermo41: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo41: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo41
