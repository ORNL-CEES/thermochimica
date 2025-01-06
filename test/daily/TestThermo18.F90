
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo33.F90
    !> \brief   Spot test - W-Au-Ar-Ne-O, 900K.
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
    !    08/31/2018    B.W.N. Fitzpatrick  Modification to use a fictive system
    !    04/17/2024    A.E.F. Fitzsimmons  Naming convention change
    !    08/27/2024    A.E.F. Fitzsimmons  SQA, standardizing tests
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a fictive system labelled W-Au-Ar-Ne-O at 900K.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo18

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
    cThermoFileName        = DATA_DIRECTORY // 'WAuArNeO-2.dat'

    ! Specify values:
    dPressure              = 2D0
    dTemperature           = 900D0
    dElementMass(74)       = 20D0        ! W
    dElementMass(79)       = 2D0         ! Au
    dElementMass(18)       = 7D0         ! Ar
    dElementMass(8)        = 5D0         ! O
    dElementMass(10)       = 1D0         ! Ne

    ! Init test values
    dGibbsCheck            = 3.06480D06
    dHeatCapacityCheck     = 332.449
    nSpeciesTest           = 4
    iSpeciesIndexTest      = [1, 5, 3, 9] !Au, Au, O, Ne
    dMolFractionTest       = [0.75306D0, 2.75309D-59, 0.24693D0, 3.09174D-02]
    lPass                  = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    call HeatCapacity

    !Execute the test for mole fractions, gibbs energy and heat capacity
    call testMolFraction(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    ! Deallocation
    deallocate(iSpeciesIndexTest, dMolFractionTest)
    
    ! Check results:
    if (lPass) then
        ! The test passed:
        print *, 'TestThermo18: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo18: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if


    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo18
