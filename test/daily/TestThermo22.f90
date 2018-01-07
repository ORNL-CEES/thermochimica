

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
    ! homogenous system containing a solution phase represented by a Redlich-Kister-Muggiano-Polynomial (RKMP)
    ! model with a ternary parameter.  The mole fractions of two constituents are checked in
    ! addition to the integral Gibbs energy of the system.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo22

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0

    ! Specify units:
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/RKMP02.dat'

    ! Specify values:
    dTemperature            = 500D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(92)        = 1D0
    dElementMass(94)        = 2D0
    dElementMass(8)         = 3D0
    dElementMass(54)        = 4D0


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if ((INFOThermo == 0).AND.(DABS(dMolFraction(2) - 0.2D0) < 1D-4).AND. &
        (DABS(dMolFraction(4) - 0.4D0) < 1D-4).AND.(DABS(dGibbsEnergySys - 1.21371E+06) < 10D0)) then
        ! The test passed:
        print *, 'TestThermo22: PASS'
    else
        ! The test failed.
        print *, 'TestThermo22: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo22
