
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo59.F90
    !> \brief   Spot test - Fe-V-O 2000 K.
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
    !    08/10/2020    M. Poschmann         Original code
    !    04/17/2024    A.E.F. Fitzsimmons   Naming convention change
    !    09/17/2024    A.E.F. Fitzsimmons   Remodle SQA
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a system including SUBQ solution phase with compound constituent.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo38

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS
    USE ModuleTesting

    implicit none

    !Init variables
    logical :: lPass
    real(8) :: dGibbsCheck, dHeatCapacityCheck
    integer :: nSpeciesTest
    integer, allocatable :: iSpeciesIndexTest(:)
    real(8), allocatable :: dMolFractionTest(:)


    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'FeTiVO.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2000D0
    dElementMass(8)        = 2D0              ! O
    dElementMass(23)       = 0.5D0            ! V
    dElementMass(26)       = 0.5D0            ! Fe

    !Init test values
    dGibbsCheck            = -8.92127D05
    dHeatCapacityCheck     = 133.351
    nSpeciesTest           = 4
    iSpeciesIndexTest      = [2, 3, 6, 10] !V, Fe, FeO, Fe[3+]-Fe[2+]-O-O
    dMolFractionTest       = [3.8002D-17, 1.3197D-08, 9.0072D-07, 0.15646D0]
    lPass                  = .FALSE.
    
    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    call HeatCapacity

    !Execute the test for mole fractions, gibbs energy and heat capacity
    call testMolFraction(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    deallocate(iSpeciesIndexTest, dMolFractionTest)

    if (lPass) then
        ! The test passed:
        print *, 'TestThermo38: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo38: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo38
