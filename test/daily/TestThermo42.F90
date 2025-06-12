program TestThermo42

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
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
    cThermoFileName       = DATA_DIRECTORY //'ClAlNa.dat'

    ! Specify values:
    dTemperature          = 1000D0
    dPressure             = 1.0D0
    dElementMass(17)      = 3D0                              ! Cl
    dElementMass(11)      = 1D0                              ! Na
    dElementMass(13)      = 1D0                              ! Al

    !Init test values
    dGibbsCheck           = -1.17685D06
    dHeatCapacityCheck    = 3.10958D02
    nSpeciesTest          = 3
    iSpeciesIndexTest     = [8, 10, 12] !Al-Al-Cl-Cl, Na-Al-Cl-Cl, Al-Al2-Cl-Cl
    dMolFractionTest      = [3.1217D-03, 0.54214D0, 6.6187D-04]
    lPass                 = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica
    call HeatCapacity

    !Execute the test for mole fractions, gibbs energy and heat capacity
    call testProperties(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    ! Deallocation
    deallocate(iSpeciesIndexTest, dMolFractionTest)

    if (lPass) then
        ! The test passed:
        print *, 'TestThermo42: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo42: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo42
