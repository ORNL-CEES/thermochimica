
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo33.F90
    !> \brief   Spot test - W-Au-Ar-Ne-O, 900K.
    !> \author  M.H.A. Piro, B.W.N. Fitzpatrick
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    05/14/2013    M.H.A. Piro         Original code
    !    08/31/2018    B.W.N. Fitzpatrick  Modification to use a fictive system
    !    04/17/2024    A.E.F. Fitzsimmons  Naming convention change
    !    08/27/2024    A.E.F. Fitzsimmons  SQA, standardizing tests
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a fictive system labelled W-Au-Ar-Ne-O at 900K.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo18

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none
    !Gibbs Energy result via Factsage
    real(8) :: gibbscheck
    real(8) :: molFractionTest(2)

    !Init test values
    gibbscheck = 3.06480D+06
    molFractionTest = [0.75306881663786374D0,3.0917444033201544D-002]


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'WAuArNeO-2.dat'

    ! Specify values:
    dPressure              = 2D0
    dTemperature           = 900D0
    dElementMass(74)       = 20D0        ! W
    dElementMass(79)       = 2D0         ! Au
    dElementMass(18)       = 7D0         ! Ar
    dElementMass(8)        = 5D0         ! O
    dElementMass(10)       = 1D0         ! Ne

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        if ((DABS(dMolFraction(1) - molFractionTest(1))/molFractionTest(1) < 1D-3).AND. &
        (DABS(dMolFraction(9) - molFractionTest(2))/molFractionTest(2) < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (gibbscheck))/(gibbscheck) < 1D-3))  then
            ! The test passed:
            print *, 'TestThermo18: PASS'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(0)
        else
            ! The test failed.
            print *, 'TestThermo18: FAIL <---'
            ! Reset Thermochimica:
            call ResetThermo
            call EXIT(1)
        end if
    else
        ! The test failed.
        print *, 'TestThermo18: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    !call ThermoDebug

    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo18
