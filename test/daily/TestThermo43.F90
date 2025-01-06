
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo64.F90
    !> \brief   Spot test - Nb-Zr-O-H 600 K.
    !> \author  M. Poschmann
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
    !    08/31/2021    M. Poschmann         Original code
    !    10/28/2024    A.E.F. Fitzsimmons   SQA Remodle
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that ternary mixing SUBL cases with only
    !!  one specified coefficient are handled correctly.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo43

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS
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
    cThermoFileName        = DATA_DIRECTORY // 'ZIRC-test64.dat'

    ! Specify values:
    dPressure              = 100D0
    dTemperature           = 600D0
    dElementMass(41)       = 1D0            ! Nb
    dElementMass(40)       = 1D0            ! Zr
    dElementMass(8)        = 1D0            ! O
    dElementMass(1)        = 0.1D0          ! H

    ! Init test values
    dGibbsCheck            = -5.24838D05
    dHeatCapacityCheck     = 72.479
    nSpeciesTest           = 7
    iSpeciesIndexTest      = [1, 3, 5, 7, 10, 15, 17] !H, H2O, HO, HZr, NbO2, Zr2, ZrO2
    dMolFractionTest       = [1.1114D-20, 3.8973D-51, 5.8381D-46, 3.0155D-47, 7.4787D-48, 2.7307D-79, 7.9320D-43]
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
        print *, 'TestThermo43: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo43: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo43
