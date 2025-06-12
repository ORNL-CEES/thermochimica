
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo65.F90
    !> \brief   Spot test - 400K with 25% Mo, 25% Ru, 25% Pd, 25% Tc.
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
    !    10/01/2021    M. Poschmann         Original code
    !    05/06/2024    A.E.F. Fitzsimmons   Naming convention change
    !    10/28/2024    A.E.F. Fitzsimmons   SQA Remodle
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for the Pd-Ru-Tc-Mo system at 400K with 25% Mo, 25% Ru, 25% Pd, 25% Tc. The database is
    !! modified such that it contains only the BCC phase, but with ternary miscibility gap possible.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo44

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
    cThermoFileName        = DATA_DIRECTORY // 'ternaryMiscibility-Kaye.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 400D0
    dElementMass(42)       = 1D0        ! Mo
    dElementMass(43)       = 1D0        ! Tc
    dElementMass(44)       = 1D0        ! Ru
    dElementMass(46)       = 1D0        ! Pd

    ! Init test values
    dGibbsCheck            = -4.48928D04
    dHeatCapacityCheck     = 118.263
    nSpeciesTest           = 6
    iSpeciesIndexTest      = [1, 2, 3, 4, 10, 11] !Mo, Pd, Tc, Ru, Pd, Tc
    dMolFractionTest       = [2.1655D-04, 4.0808D-03, 1.1254D-01, 8.8316D-01, 9.2991D-01, 6.6589D-02]
    lPass                  = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    call HeatCapacity

    ! Execute the test for mole fractions, gibbs energy and heat capacity
    call testProperties(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    ! Deallocation
    deallocate(iSpeciesIndexTest, dMolFractionTest)

    if (lPass) then
        ! The test passed:
        print *, 'TestThermo44: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo44: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo44
