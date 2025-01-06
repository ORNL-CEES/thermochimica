
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo88.F90
    !> \brief   SUBM mismatch coefficients and charges.
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
    !    11/11/2022    M. Poschmann        SUBM Test Case
    !    04/17/2024    A.E.F. Fitzsimmons  Naming convention change
    !    10/28/2024    A.E.F. Fitzsimmons  SQA Remodle
    !    
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results a SUBM phase where not all charges and coefficients are equal, i.e. with Zr(4+) + O(2-) -> ZrO2.
    !!  Two constituents on both sublattices are tested.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo62

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

    dElementMass(3)       = 3D0                              ! Li
    dElementMass(11)      = 0D0                              ! Na
    dElementMass(17)      = 0D0                              ! Cl
    dElementMass(9)       = 5D0                              ! F
    dElementMass(26)      = 0D0                              ! Fe
    dElementMass(8)       = 1D0                              ! O
    dElementMass(19)      = 0D0                              ! K
    dElementMass(40)      = 1D0                              ! Zr

    ! Init test values
    dGibbsCheck           = -1.90265D05
    dHeatCapacityCheck    = 3.72529E-06
    nSpeciesTest          = 3
    iSpeciesIndexTest     = [1, 2, 3, 4] !Pd, Pd, Ru
    dMolFractionTest      = [5.9999D-01, 1.2000D-01, 1.9999D-01, 8.0000D-02]
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
        print *, 'TestThermo62: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo62: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo62
