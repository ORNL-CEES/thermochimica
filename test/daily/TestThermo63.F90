
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo89.F90
    !> \brief   SUBM worst case scenario.
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
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results a SUBM phase where everything and the kitchen sink is thrown in.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo63

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
    dElementMass(8)       = 11D0                             ! O
    dElementMass(19)      = 4D0                              ! K
    dElementMass(40)      = 5D0                              ! Zr

    ! Init test values
    dGibbsCheck           = -7.61628D05
    dHeatCapacityCheck    = 1198.68
    nSpeciesTest          = 7
    iSpeciesIndexTest     = [1, 2, 3, 7, 8, 10, 11] !Pd, Pd, Ru
    dMolFractionTest      = [2.7906D-02, 4.6511D-03, 5.5813D-02, 8.3720D-02, 1.3953D-02, 1.1162D-01, 1.8604D-02]
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
        print *, 'TestThermo63: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo63: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo63
