program TestThermo41

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
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY //'ClAlNa.dat'

    ! Specify values:
    dTemperature          = 2000D0
    dPressure             = 1.0D0
    dElementMass(17)      = 2D0                              ! Cl
    dElementMass(13)      = 1D0                              ! Al

   !Init test values
    dGibbsCheck            = -9.64834D05
    dHeatCapacityCheck     = 6.84308D01
    nSpeciesTest           = 3
    iSpeciesIndexTest      = [2, 4, 8] !Al2Cl6, Al2-Al2-Cl-Cl, Al-Al2-Va-Va
    dMolFractionTest       = [1.5843D-04, 5.3213D-03, 0.36809D0]
    lPass                  = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    call HeatCapacity
    
    !Execute the test for mole fractions, gibbs energy and heat capacity
    call testProperties(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    ! Deallocation
    deallocate(iSpeciesIndexTest, dMolFractionTest)

    if (lPass) then
        ! The test passed:
        print *, 'TestThermo41: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo41: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo41
