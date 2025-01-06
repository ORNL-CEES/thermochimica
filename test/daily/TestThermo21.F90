
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo42.F90
    !> \brief   Spot test - 2310K with 22% Mo, 78% Tc.
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
    !! results for the Pd-Ru-Tc-Mo system at 2310K with 22% Mo, 78% Tc.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo21

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
    dTemperature           = 2310D0
    dElementMass(42)       = 0.22D0        ! Mo
    dElementMass(43)       = 0.78D0        ! Tc

    ! Init test values
    dGibbsCheck            = -1.55880D05
    dHeatCapacityCheck     = 111.147
    nSpeciesTest           = 3
    iSpeciesIndexTest      = [2, 7, 8] !Mo2, Mo, Tc
    dMolFractionTest       = [3.081021D-13, 7.50680D-01, 3.14203D-01]
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
        print *, 'TestThermo21: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo21: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo21
