
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo54.F90
    !> \brief   Spot test - 1973K with 10% Mo, 30% Pd, 60% Ru.
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
    !    08/31/2018    B.W.N. Fitzpatrick  Modification to use Kaye's Pd-Ru-Tc-Mo system
    !    05/06/2024    A.E.F. Fitzsimmons  Naming convention change
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for the Pd-Ru-Tc-Mo system at 1973K with 10% Mo, 30% Pd, 60% Ru.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo33

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleTesting

    implicit none

    ! Init variables
    logical :: lPass
    real(8) :: dGibbsCheck, dHeatCapacityCheck
    integer :: nSpeciesTest
    integer, allocatable :: iSpeciesIndexTest(:)
    real(8), allocatable :: dMolFractionTest(:)

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1800D0
    dElementMass(43)       = 0.01D0        ! Tc
    dElementMass(46)       = 0.09D0        ! Pd
    dElementMass(42)       = 0.9D0        ! Mo

    ! Init test values
    dGibbsCheck            = -1.05595D05
    dHeatCapacityCheck     = 40.0333
    nSpeciesTest           = 6
    iSpeciesIndexTest      = [10, 11, 12, 19, 20, 21] !Tc, Pd, Mo, Tc, Mo, Pd
    dMolFractionTest       = [6.8001D-04, 9.3597D-01, 5.4310D-02, 1.3208D-02, 4.7934D-01, 5.1949D-01]
    lPass                  = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    call HeatCapacity

    ! Execute the test for mole fractions, gibbs energy and heat capacity
    call testMolFraction(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    ! Deallocation
    deallocate(iSpeciesIndexTest, dMolFractionTest)

    if (lPass) then
        ! The test passed:
        print *, 'TestThermo33: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo33: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo33
