program mass

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: i,j,k,ierr
    real(8) :: dPhaseMass, dSpeciesMass, dTotalMass, dFraction

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'Kaye_NobleMetals.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1800D0
    dElementMass(43)       = 0.01D0        ! Tc
    dElementMass(46)       = 0.09D0        ! Pd
    dElementMass(42)       = 0.9D0        ! Mo
    iPrintResultsMode = 2

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)


    ! Call Thermochimica:
    call Thermochimica

    call PrintResults

    do i = 1, nSolnPhases
      k = -iAssemblage(nElements + 1 - i)
      call getPhaseMass(k, dPhaseMass, ierr)
      call getPhaseMassFraction(k, dFraction, ierr)
      print *, cSolnPhaseName(k), dPhaseMass, dFraction
      do j = nSpeciesPhase(k - 1) + 1, nSpeciesPhase(k)
        call getSpeciesMass(j, dSpeciesMass, ierr)
        call getSpeciesMassFraction(j, dFraction, ierr)
        print *, cSpeciesName(j), dSpeciesMass, dFraction
      end do
    end do

    call getTotalMass(dTotalMass, ierr)
    print *, dTotalMass

    ! Reset Thermochimica:
    call ResetThermo

end program mass
