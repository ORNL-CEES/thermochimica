
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo85.F90
    !> \brief   SUBM all mixing.
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
    !    10/28/2024    A.E.F. Fitzsimmons  SQA Remodle
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results a SUBM phase with all excess mixing terms.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo59

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

    dElementMass(3)       = 1D0                              ! Li
    dElementMass(11)      = 2D0                              ! Na
    dElementMass(17)      = 2D0                              ! Cl
    dElementMass(9)       = 12D0                             ! F
    dElementMass(26)      = 3D0                              ! Fe
    dElementMass(8)       = 1D0                              ! O
    dElementMass(19)      = 4D0                              ! K

    ! Init test values
    dGibbsCheck           = -5.64331D05
    dHeatCapacityCheck    = 2219.86
    nSpeciesTest          = 5
    iSpeciesIndexTest     = [1, 2, 4, 11, 12] !LiF, LiCl, NaCl, KCl, K2O
    dMolFractionTest      = [8.0000D-02, 1.3333D-02, 2.6666D-02, 5.3333D-02, 2.6666D-02]
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
        print *, 'TestThermo59: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo59: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo59
