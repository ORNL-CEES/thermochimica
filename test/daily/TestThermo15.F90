
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo30.F90
    !> \brief   Spot test - W-Au-Ar-O.
    !> \author  M.H.A. Piro, B.W.N. Fitzpatrick
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer           Description of change
    !    ----          ----------           ---------------------
    !    05/14/2013    M.H.A. Piro          Original code
    !    31/08/2018    B.W.N. Fitzpatrick   Change system to fictive RKMP model
    !    04/17/2024    A.E.F. Fitzsimmons   Naming convention change
    !    08/27/2024    A.E.F. Fitzsimmons   SQA, standardizing tests
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a fictive system labelled W-Au-Ar-O.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo15

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
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
    cThermoFileName        = DATA_DIRECTORY // 'WAuArO-1.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1455D0
    dElementMass(74)       = 1.95D0        ! W
    dElementMass(79)       = 1D0           ! Au
    dElementMass(18)       = 2D0           ! Ar
    dElementMass(8)        = 10D0          ! O
    
    ! Init test values
    dGibbsCheck            = -4.62037D05
    dHeatCapacityCheck     = 360D0
    nSpeciesTest           = 4
    iSpeciesIndexTest      = [1, 2, 3, 4] !Au, W, O, Ar
    dMolFractionTest       = [0.66896D-01, 0.13043D0, 0.66889D0, 0.13377D0 ]
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

    ! Check results:
    if (lPass) then
        ! The test passed:
        print *, 'TestThermo15: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo15: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Reset Thermochimica:
    call ResetThermo

end program TestThermo15
