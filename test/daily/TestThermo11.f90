
program TestThermo11

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
    ! The purpose of this unit test is to ensure that Thermochimica converges for a system representing combustion
    ! of methane in air (O-N-C-H system).  This particular unit test is taken from the Verification Manual.  The
    ! thermochemical activity of O2(g) is checked in addition to INFOThermo.
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Initialize variables:
    dTemperature            = 2000D0
    dPressure               = 1D0
    dElementMass(1)         = 4D0       ! H
    dElementMass(6)         = 1D0       ! C
    dElementMass(7)         = 15.04D0   ! N
    dElementMass(8)         = 4D0       ! O
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/Example12.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if ((INFOThermo == 0).AND.(dMolFraction(39) >= 1.5D-3).AND.(dMolFraction(39) <= 1.7D-3)) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo11: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo11: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo11
