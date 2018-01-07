

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
    !    08/21/2015    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that the U-La-O is computed correctly. The fluorite
    ! oxide should be the only phase that is stable.  More details about this treatment are provided below:
    !
    !   D. Shin and T.M. Besmann, "Thermodynamic modeling of the (U,La)O_{2\pm x} solid solution phase,"
    !   Journal of Nuclear Materials, 433 (2013) 227-232.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo36

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = '../data/ULaO.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 783D0
    dElementMass           = 0d0

    dElementMass(92)       = 3.6971D3          ! U
    dElementMass(57)       = 2.5200D1          ! La
    dElementMass(8)        = 7.8030D3          ! O

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(25) - 6.118D-3)/6.118D-3 < 1D-3).AND. &
        (DABS(dMolFraction(17) - 0.71745D0)/0.71745D0 < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo36: PASS'
        else
            ! The test failed.
            print *, 'TestThermo36: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo36: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo36
