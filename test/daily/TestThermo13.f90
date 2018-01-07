
program TestThermo13

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
    !    02/09/2012    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this unit test is to ensure that Thermochimica converges for a system representing
    ! irradiated nuclear fuel.  Thermochemical data is from the RMC Fuel Thermochemical Treatment (RFTT) and
    ! this exact scenario is validated by Corcoran et al (reference below) by oxidation experiments.  This
    ! particular unit test is taken from the Verification Manual.  The thermochemical activity of O2(g)
    ! is checked in addition to INFOThermo.
    !
    !   E.C. Corcoran, B.J. Lewis, W.T. Thompson, J. Mouris, Z. He, Controlled Oxidation Experiments of
    !   Simulated Irradiated UO2 Fuel in Relation to Thermochemical Modelling, Journal of Nuclear Materials,
    !   414 (2011) 73-82.
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Initialize variables:
    dTemperature            = 1173.15D0
    dPressure               = 1D0
    dElementMass(94)        = 0D0
    dElementMass(93)        = 0D0
    dElementMass(92)        = 0.56849D-2
    dElementMass(60)        = 0.18974d-6
    dElementMass(59)        = 0D0
    dElementMass(58)        = 0.75956d-7
    dElementMass(57)        = 0.3831d-7
    dElementMass(56)        = 0.99634d-7
    dElementMass(55)        = 0.12584d-8
    dElementMass(54)        = 0D0
    dElementMass(53)        = 0D0
    dElementMass(52)        = 0D0
    dElementMass(46)        = 0.10003d-4
    dElementMass(45)        = 0.73875d-6
    dElementMass(44)        = 0.42121d-5
    dElementMass(43)        = 0D0
    dElementMass(42)        = 0.14263d-4
    dElementMass(40)        = 0.11667d-6
    dElementMass(39)        = 0D0
    dElementMass(38)        = 0.86762d-7
    dElementMass(37)        = 0D0
    dElementMass(8)         = 1.14467207D-2
    dElementMass(1)         = 0D0
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/Example7c.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if ((INFOThermo == 0).AND.(dMolFraction(2) >= 6D-13).AND.(dMolFraction(2) <= 6.5D-13)) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo13: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo13: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo13
