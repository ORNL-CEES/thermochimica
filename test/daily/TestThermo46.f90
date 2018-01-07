

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
    ! The purpose of this application test is to ensure that the Mo-Zr binary system
    ! is correctly computed.  This data-file was kindly provided by T.M. Besmann.
    ! This particular application test considers the Mo-Zr system at a particular
    ! temperature and composition that should correspond to only the C15
    ! Laves phase.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo46

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
    dTemperature            = 1800D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(40)        = 0.35D0     ! Zr
    dElementMass(42)        = 0.65D0     ! Mo

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(18) - 5.8037D-4)/5.8037D-4 < 1D-3).AND. &
        (DABS(dMolFraction(16) - (0.95057D0))/(0.95057D0) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo46: PASS'
        else
            ! The test failed.
            print *, 'TestThermo46: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo46: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo46
