

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
    !    07/29/2014    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that the subroutine computes the correct (O/M)f
    ! of the UO2 solid solution phase.  The mass of each chemical element corresponds to a burnup of 100 GWd/t(U)
    ! of PWR nuclear fuel, as provided by T.M. Besmann (note: I think K.T. Clarno did the isotopic calculations,
    ! but I may be wrong).  This application test tests the CompOtoMRatio subroutine.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo35

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0

    integer :: INFO
    character(8):: cPhaseIn
    real(8):: dOtoMRatio


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = '../data/DBV6_TMB_modified.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2000D0
    dElementMass           = 0d0

    dElementMass(94)       = 5.7050D01         ! Pu
    dElementMass(92)       = 3.6971D3          ! U
    dElementMass(60)       = 7.8800D1          ! Nd
    dElementMass(59)       = 2.1840D1          ! Pr
    dElementMass(58)       = 5.4720D1          ! Ce
    dElementMass(57)       = 2.5200D1          ! La
    dElementMass(56)       = 3.5550D1          ! Ba
    dElementMass(55)       = 5.7380D1          ! Cs
    dElementMass(54)       = 1.3150D2          ! Xe
    dElementMass(53)       = 5.3120D0          ! I
    dElementMass(52)       = 1.2290D1          ! Te
    dElementMass(46)       = 6.7630D1          ! Pd
    dElementMass(45)       = 7.5310D0          ! Rh
    dElementMass(44)       = 8.9420D1          ! Ru
    dElementMass(43)       = 1.8530D1          ! Tc
    dElementMass(42)       = 1.0080D2          ! Mo
    dElementMass(40)       = 9.8930D1          ! Zr
    dElementMass(39)       = 1.1550D1          ! Y
    dElementMass(38)       = 2.1870D1          ! Sr
    dElementMass(37)       = 9.3000D0          ! Rb
    dElementMass(8)        = 8.4030D3           ! O


    ! Initialize variables for the CompOtoMRatio subroutien:
    INFO       = 0
    dOtoMRatio = 0D0
    cPhaseIn = 'O2ZRU_C'

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Compute the oxygen to metal ratio:
    call CompOtoMRatio(cPhaseIn, dOtoMRatio, INFO)

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(147) - 0.85312D0)/0.85312D0 < 1D-3).AND. &
        (DABS(dMolFraction(258) - 0.22088D0)/0.22088D0 < 1D-3).AND. &
        (DABS(dOtoMRatio - 2.007125D0)/(2.007125D0) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo35: PASS'
        else
            ! The test failed.
            print *, 'TestThermo35: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo35: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo35
