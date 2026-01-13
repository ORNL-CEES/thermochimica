
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo86.F90
    !> \brief   SUBM mistmatch coefficients and charges.
    !> \author  M.H.A. Piro, M. Poschmann
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
    !    11/08/2022    M. Poschmann        SUBM Test Case
    !    04/17/2024    A.E.F. Fitzsimmons  Naming convention change
    !    10/28/2024    A.E.F. Fitzsimmons  SQA Remodle
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results a SUBM phase where not all charges and coefficients are equal, i.e. with Zr(4+) + O(2-) -> ZrO2.
    !!  Two constituents on first sublattice are tested.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo60

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
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'ZrFeKClNaFOLi.dat'

    ! Specify values:
    dTemperature          = 1000
    dPressure             = 1.0D0
    dElementMass          = 0D0

    dElementMass(3)       = 2D0                              ! Li
    dElementMass(8)       = 3D0                              ! O
    dElementMass(40)      = 1D0                              ! Zr

    ! Init test values
    dGibbsCheck           = -1.01677D05
    dHeatCapacityCheck    = 3189.08
    nSpeciesTest          = 2
    iSpeciesIndexTest     = [1, 2] !Li2O, ZrO2
    dMolFractionTest      = [4.9999D-01, 5.0000D-01]
    lPass                 = .FALSE.

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
        print *, 'TestThermo60: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo60: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo60
