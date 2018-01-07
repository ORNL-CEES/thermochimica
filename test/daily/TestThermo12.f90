
program TestThermo12

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
    !    02/07/2012    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this unit test is to ensure that Thermochimica converges for a system representing the
    ! noble metal phases, which are also known as the "white phase" in irradiated nuclear fuel.  This particular
    ! unit test is taken from the Verification Manual.  The mole fraction of Mo in the body centred cubic phase
    ! is checked in addition to INFOThermo.
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Initialize variables:
    dTemperature            = 850D0
    dPressure               = 1D0
    dElementMass(42)        = 100D0     ! Mo
    dElementMass(43)        = 10D0      ! Tc
    dElementMass(44)        = 1D0       ! Ru
    dElementMass(45)        = 0.1D0     ! Rh
    dElementMass(46)        = 0.01D0    ! Pd
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/Example6b.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if ((INFOThermo == 0).AND.(dMolFraction(12) >= 0.9D0).AND.(dMolFraction(12) <= 0.91D0)) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo12: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo12: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo12
