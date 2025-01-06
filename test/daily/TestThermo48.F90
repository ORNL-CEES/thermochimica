
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo73.F90
    !> \brief   Spot test - 2500K with 33% Ca, 33% Mn, 33% S.
    !> \author  M.H.A. Piro, B.A.T. Breeden
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
    !    01/10/2022    B.A.T. Breeden      Modification to use Dupin's Zirc Data base with SUBI
    !    04/17/2024    A.E.F. Fitzsimmons  Naming convention change
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results for the open literature Ca-Mn-S assessment file at 2500K with 1 mol of Ca, 1 mol of Mn, and
    !!  1 mol of S. It also tests mixing term Case #1, # 2, and #3 of the SUBI phase with a miscibility gap
    !!  present.
    !!  The DAT file was pulled from the following article. However, modifications may have been made
    !!  from the original version: D. Dilner, "Thermodynamic description of the Fe–Mn–Ca–Mg–S system," Calphad,
    !!  vol. 53, pp. 55-61, 2016.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo48

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleTesting

    implicit none

    !Init variables
    logical :: lPass
    real(8) :: dGibbsCheck, dHeatCapacityCheck
    integer :: nSpeciesTest
    integer, allocatable :: iSpeciesIndexTest(:)
    real(8), allocatable :: dMolFractionTest(:)

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'CaMnS.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2500D0
    dElementMass(20)       = 1D0          ! Ca
    dElementMass(25)       = 1D0          ! Mn
    dElementMass(16)       = 1D0          ! S

    !Init test values
    dGibbsCheck            = -9.10619D05
    dHeatCapacityCheck     = -5.32116D03
    nSpeciesTest           = 3
    iSpeciesIndexTest      = [1, 19, 20] !S, Ca:S, Mn:S
    dMolFractionTest       = [2.8964D-06, 5.4737D-03, 9.9452D-01]
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
        print *, 'TestThermo48: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo48: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo48
