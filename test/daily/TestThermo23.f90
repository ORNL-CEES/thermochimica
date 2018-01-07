

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
    !    07/17/2013    M.H.A. Piro       Modified data-file and updated expected results.
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that Thermochimica computes the correct results for a
    ! homogenous system containing a solution phase represented by a Redlich-Kister-Muggiano-Polynomial (RKMP)
    ! model with a quaternary parameter.  The mole fractions of two constituents are checked in
    ! addition to the integral Gibbs energy of the system.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo23

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Specify units:
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/RKMP05.dat'

    ! Specify input values:
    dTemperature            = 1496.3D0
    dPressure               = 1D0
    dElementMass(92)        = 99D0
    dElementMass(94)        = 56D0
    dElementMass(8)         = 96D0
    dElementMass(54)        = 74D0
    dElementMass(36)        = 36D0


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then

        if ((DABS(dMolFraction(2) - 0.15512D0)/0.15512D0 < 1D-3).AND. &
            (DABS(dMolFraction(5) - 9.9723D-02)/9.9723D-02 < 1D-3).AND.(DABS(dGibbsEnergySys - &
            (-15758228.9057433D0))/(-15758228.9057433D0) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo23: PASS'
        else
            ! The test failed.
            print *, 'TestThermo23: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo23: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo23
