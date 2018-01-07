

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
    !    01/16/2013    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that Thermochimica computes the correct results for a
    ! heterogeneous system containing a solution phase represented by a Redlich-Kister-Muggiano-Polynomial (RKMP)
    ! model with binary, ternary and quaternary parameters.  The mole fractions of two constituents are checked in
    ! addition to the integral Gibbs energy of the system.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo24

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0

    ! Specify units:
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/RKMP04.dat'

    ! Specify values:
    dTemperature            = 327.4D0
    dPressure               = 1D0

    dElementMass(92)        = 947D0
    dElementMass(94)        = 731D0
    dElementMass(8)         = 497D0
    dElementMass(54)        = 374D0
    dElementMass(36)        = 421D0

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica


    ! NOTE:
    ! =====
    ! For some strange reason, FactSage reports a different value for the integral Gibbs energy
    ! for this test problem; however, the molar quantities for all phases and the mole fractions
    ! of all constituents are precise to the fifth significant figure.  I did a hand calculation
    ! of this problem and it is consistent with Thermochimica.  Very strange.

    ! Check results:
    if ((INFOThermo == 0).AND.(DABS(dMolFraction(2) - 0.2778391D0)/0.2778391D0 < 1D-3).AND. &
        (DABS(dMolFraction(5) - 0.456534277D0)/0.456534277D0 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - 699432566D0)/699432566D0 < 1D-3)) then
        ! The test passed:
        print *, 'TestThermo24: PASS'
    else
        ! The test failed.
        print *, 'TestThermo24: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo24
