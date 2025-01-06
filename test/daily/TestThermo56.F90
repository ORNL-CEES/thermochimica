
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo82.F90
    !> \brief   SUBM second sublattice binary mixing.
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
    !    11/04/2022    M. Poschmann        SUBM Test Case
    !    04/17/2024    A.E.F. Fitzsimmons  Naming convention change
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results a SUBM phase with a binary excess mixing term on the second sublattice.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo56

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

    dElementMass(3)       = 1.2D0                              ! Li
    dElementMass(17)      = 0.5D0                              ! Cl
    dElementMass(9)       = 0.3D0                              ! F
    dElementMass(8)       = 0.2D0                              ! O

    ! Init test values
    dGibbsCheck           = -4.52110D04
    dHeatCapacityCheck    = -4.92203D-03
    nSpeciesTest          = 3
    iSpeciesIndexTest     = [1, 2, 3] !LiF, LiCl, Li2O
    dMolFractionTest      = [2.9999D-01, 5.0000D-01, 2.0000D-01]
    lPass                 = .FALSE.

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
        print *, 'TestThermo56: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo56: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo56
