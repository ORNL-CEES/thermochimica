

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
    ! model with a zeroth order binary parameter.  The mole fractions of two constituents are checked in
    ! addition to the integral Gibbs energy of the system.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo21

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0

    ! Specify units:
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/RKMP01.dat'

    ! Specify values:
    dTemperature            = 1000D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(92)        = 0.45D0
    dElementMass(94)        = 0.05D0
    dElementMass(8)         = 0.2D0
    dElementMass(54)        = 0.3d0


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if ((INFOThermo == 0).AND.(DABS(dMolFraction(2) - 5.0000E-02) < 1D-4).AND. &
        (DABS(dMolFraction(4) - 0.30000D0) < 1D-4).AND.(DABS(dGibbsEnergySys - 1.06831E+05) < 10D0)) then
        ! The test passed:
        print *, 'TestThermo21: PASS'
    else
        ! The test failed.
        print *, 'TestThermo21: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo21
