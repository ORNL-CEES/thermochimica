

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
    ! oxide is stable in equilibrium with U and La metallic phases.
    ! More details about this treatment are provided below:
    !
    !   D. Shin and T.M. Besmann, "Thermodynamic modeling of the (U,La)O_{2\pm x} solid solution phase,"
    !   Journal of Nuclear Materials, 433 (2013) 227-232.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo37

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
    dTemperature           = 1000D0
    dElementMass           = 0D0

    dElementMass(92)       = 1D0        ! U
    dElementMass(57)       = 0.05D0     ! La
    dElementMass(8)        = 1.95D0     ! O

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(25) - 1.5754D-03)/1.5754D-03 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-1.17525D6))/(-1.17525D6) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo37: PASS'
        else
            ! The test failed.
            print *, 'TestThermo37: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo37: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo37
