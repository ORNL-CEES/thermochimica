
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo75.F90
    !> \brief   Spot test - 1500K with 0.5 Ca - 0.2 Mn - 0.4 Fe.
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
    !    02/24/2022    B.A.T. Breeden      SUBI Test Case
    !    04/17/2024    A.E.F. Fitzsimmons  Naming convention change
    !    08/27/2024    A.E.F. Fitzsimmons  Remodel
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results for the open literature Ca-Mn-S assessment file at 1500K with 0.5 mol of Ca, 0.2 mol of Mn, and
    !!  0.4 mol of Fe. It also tests mixing term Case #8 of the SUBI phase.
    !!
    !!  The DAT file was pulled from the following article. However, modifications may have been made
    !!  from the original version: D. Dilner, "Thermodynamic description of the Fe–Mn–Ca–Mg–S system," Calphad,
    !!  vol. 53, pp. 55-61, 2016.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo50

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
    cThermoFileName        = DATA_DIRECTORY // 'FeMnCaS-2.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1500D0
    dElementMass(26)       = 0.4D0          ! Fe
    dElementMass(25)       = 0.2D0          ! Mn
    dElementMass(20)       = 0.5D0          ! Ca

    ! Init test values
    dGibbsCheck            = -1.16594D05
    dHeatCapacityCheck     = 42.6874
    nSpeciesTest           = 3
    iSpeciesIndexTest      = [1, 2, 3] !Ca+2:Va, Fe+2:Va, Mn+2:Va
    dMolFractionTest       = [4.5454D-01, 3.6363D-01, 1.8181D-01]
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
        print *, 'TestThermo50: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo50: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo50
