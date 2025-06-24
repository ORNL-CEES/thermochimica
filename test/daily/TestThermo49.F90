
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo74.F90
    !> \brief   Spot test - 1900K with , 2 mol Fe, 4 mol Mn, 3 mol S.
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
    !    08/27/2024    A.E.F. Fitzsimmons  Remodel
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

program TestThermo49
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
    cThermoFileName        = DATA_DIRECTORY // 'FeMnCaS-1.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1900D0
    dElementMass(26)       = 2D0          ! Fe
    dElementMass(25)       = 4D0          ! Mn
    dElementMass(16)       = 3D0          ! S

    ! Init test values
    dGibbsCheck            = -1.94495D06
    dHeatCapacityCheck     = 128.00
    nSpeciesTest           = 6
    iSpeciesIndexTest      = [10, 11, 12, 15, 16, 17] !Fe+2:S-2, Fe+2:Va, Mn+2:S-2, Fe+2:S-2, Fe+2:Va. Mn+2:S-2
    dMolFractionTest       = [5.5296D-05, 7.7304D-01, 1.6228D-05, 1.1106D-01, 2.1224D-01, 2.3246D-01]
    lPass                  = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    call HeatCapacity
    
    ! Execute the test for mole fractions, gibbs energy and heat capacity
    call testProperties(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    ! Deallocation
    deallocate(iSpeciesIndexTest, dMolFractionTest)

    if (lPass) then
        ! The test passed:
        print *, 'TestThermo49: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo49: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo49
