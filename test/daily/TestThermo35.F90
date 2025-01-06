
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo56.F90
    !> \brief   Spot test - Fe-Cu-C 1400 K.
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
    !    06/30/2020    M. Poschmann        Original code
    !    06/05/2024    A.E.F. Fitzsimmons  Naming convention change
    !    08/27/2024    A.E.F. Fitzsimmons  SQA, standardizing tests
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a system including magnetic solution species and excess magnetic terms. Also tests
    !! G-type asymmetric (and symmetric) excess mixing energy implementation for SUBG.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo35

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
    cThermoFileName       = DATA_DIRECTORY // 'CuFeC-Kang.dat'

    ! Specify values:
    dTemperature          = 1400D0
    dPressure             = 1.0D0
    dElementMass          = 0D0
    dElementMass(6)       = 1.0D0                              ! C
    dElementMass(26)      = 1.0D0                              ! Fe
    dElementMass(29)      = 1.0D0                              ! Cu

    !Init test values
    dGibbsCheck           = -1.73325D05
    dHeatCapacityCheck    = 1.05080D02
    nSpeciesTest          = 3
    iSpeciesIndexTest     = [1, 6, 13] !C-C-Va-Va, Cu-Fe-Va-Va, Fe
    dMolFractionTest      = [1.8840D-09, 4.4977D-02, 0.88662D0]
    lPass                 = .FALSE.
    
    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica
    call HeatCapacity
    
    !Execute the test for mole fractions, gibbs energy and heat capacity
    call testMolFraction(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    deallocate(iSpeciesIndexTest, dMolFractionTest)

    if (lPass) then
        ! The test passed:
        print *, 'TestThermo35: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo35: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo35
