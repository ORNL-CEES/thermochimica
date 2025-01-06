
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo89.F90
    !> \brief   SUBM worst case scenario.
    !> \author  M.H.A. Piro, M. Poschmann
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
    !    11/11/2022    M. Poschmann        SUBM Test Case
    !    04/17/2024    A.E.F Fitzsimmons   CsI data bug test
    !    05/06/2024    A.E.F Fitzsimmons   Naming convention update 
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results of CsI at low pressure. This test is consistent with the test in the tutorial manual.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo64

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
    cInputUnitTemperature = 'C'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'CsI-Pham.dat'

    ! Specify values:
    dTemperature          = 400
    dPressure             = 1D-5
    dElementMass(53)      = 1D0                              ! I
    dElementMass(55)      = 1D0                              ! cs

    ! Init test values
    dGibbsCheck           = -4.41869D05
    dHeatCapacityCheck    = 69.8320
    nSpeciesTest          = 4
    iSpeciesIndexTest     = [1, 3, 5, 6] !Cs, I, CsI, Cs2I2
    dMolFractionTest      = [4.2170D-18, 1.0711D-02, 2.5818D-03, 1.9349D-05]
    lPass                 = .FALSE.

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
        print *, 'TestThermo64: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo64: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo64
