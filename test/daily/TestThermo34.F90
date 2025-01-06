
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo55.F90
    !> \brief   Spot test - Ni-Cr-Fe-H 300 K.
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
    !    05/21/2020    M. Poschmann        Original code
    !    08/27/2024    A.E.F. Fitzsimmons  SQA, standardizing tests
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a system including magnetic solution species and excess magnetic terms.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo34

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
    cThermoFileName        = DATA_DIRECTORY // 'ZIRC-noSUBI.dat'

    ! Specify values:
    dPressure              = 1000D0
    dTemperature           = 300D0
    dElementMass(23)       = 0.1D0          ! V
    dElementMass(24)       = 1D0            ! Cr
    dElementMass(26)       = 2D0            ! Fe
    dElementMass(28)       = 1D0            ! Ni
    dElementMass(50)       = 0.1D0          ! Sn
    dElementMass(1)        = 1D0            ! H

    ! Init test values
    dGibbsCheck            = -4.96674D04
    dHeatCapacityCheck     = 131.162
    nSpeciesTest           = 8
    iSpeciesIndexTest      = [1, 3, 4, 8, 9, 11, 39, 60] !Cr, Fe, H, Ni, Sn, V, Cr:Fe, V:Ni
    dMolFractionTest       = [1.4724D-64, 4.1279D-68, 1.3186D-37, 5.2694D-71, 1.9488D-55, 4.4569D-89, 8.1559D-02, 5.0192D-04]
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
        print *, 'TestThermo34: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo34: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo34
