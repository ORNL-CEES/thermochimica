
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo32.F90
    !> \brief   Spot test - W-Au-Ar-Ne-O, 2452K.
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
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a fictive system labelled W-Au-Ar-Ne-O at 2452K.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo17

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
    cThermoFileName        = DATA_DIRECTORY // 'WAuArNeO-1.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2452D0
    dElementMass(74)       = 1.95D0        ! W
    dElementMass(79)       = 1D0           ! Au
    dElementMass(18)       = 2D0           ! Ar
    dElementMass(8)        = 10D0          ! O
    dElementMass(10)       = 10D0          ! Ne

    ! Init test values
    dGibbsCheck            = 1.67160D07
    dHeatCapacityCheck     = 603.816
    nSpeciesTest           = 5
    iSpeciesIndexTest      = [1, 2, 3, 4, 5] !Au, W, O, Ar, Ne
    dMolFractionTest       = [4.00801D-02, 7.81563D-02, 0.40080D0, 8.01603D-02, 0.40080D0]
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
        print *, 'TestThermo17: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo17: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo17
