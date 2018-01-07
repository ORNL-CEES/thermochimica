
program TestThermo20

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
    !    08/30/2012    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that Thermochimica computes the correct results for a
    ! scenario involving a miscibility gap under difficult conditions in the Zr-Nb system.  The reason why this
    ! partiuclar test can be difficult is because initially the system is predicted to be homogenous and the
    ! mole fractions of both constituents are 0.5.  Finding the miscibility gap depends strongly on the initial
    ! estimates of the mole fractions of the "other phase", which could be returned to the same value.
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Initialize variables:
    dTemperature            = 1200D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(40)        = 1D0       ! Zr
    dElementMass(41)        = 1D0       ! Nb

    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/Example20.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if ((INFOThermo == 0).AND.(DABS(dMolFraction(3) - 0.572740520982999D0) < 1D-3).AND. &
        (DABS(dMolFraction(5) - 0.235023883219681D0) < 1D-3)) then
        ! The test passed:
        print *, 'TestThermo20: PASS'
    else
        ! The test failed.
        print *, 'TestThermo20: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo20
