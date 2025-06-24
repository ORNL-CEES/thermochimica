
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo61.F90
    !> \brief   Tricky SUBL vacancy-vacancy
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
    !    24/06/2021    M. Poschmann         Original code
    !    04/17/2024    A.E.F. Fitzsimmons   Naming convention change
    !    09/17/2024    A.E.F. Fitzsimmons   Remodle SQA
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica includes vacancy-vacancy
    !! species in the SUBL model if appropriate.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo40

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
    dTemperature           = 1000D0
    dElementMass(8)        = 2D0              ! O
    dElementMass(22)       = 0.5D0            ! Ti
    
    !Init test values
    dGibbsCheck            = -6.24557D05
    dHeatCapacityCheck     = 5.48677D01
    nSpeciesTest           = 2
    iSpeciesIndexTest      = [3, 7, 9] !Gas: TiO2, Rutilesoln: TiO2, Ti[4+]
    dMolFractionTest       = [1.3533D-23, 9.9772D-01, 2.2805D-03]
    lPass                  = .FALSE.
    
    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    call HeatCapacity

    !Execute the test for mole fractions, gibbs energy and heat capacity
    call testProperties(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    deallocate(iSpeciesIndexTest, dMolFractionTest)

    if (lPass) then
        ! The test passed:
        print *, 'TestThermo40: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo40: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo40
