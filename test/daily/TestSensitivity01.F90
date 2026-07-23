program TestSensitivity01

  USE ModuleThermo
  USE ModuleThermoIO
  USE ModuleGEMSolver
  USE ModuleEquilibriumSensitivity

  implicit none

  integer :: INFO, restoreInfo, speciesIndex, elementColumn, bestElementColumn, elementAtomicNumber
  integer :: i, j, k, parameter, nSolSpecies
  integer :: nTokenBase, nTokenPlus, nTokenMinus
  real(8) :: rcond, implicitDerivative, exactDerivative, coarseDerivative, step, baseAmount
  real(8) :: logPlusHalf, logMinusHalf, coarseError, fineError, largestDerivative
  real(8) :: balanceDerivative, expectedDerivative
  real(8), dimension(:), allocatable :: baseMolFraction, baseMolesSpecies, baseChemicalPotential
  real(8), dimension(:), allocatable :: baseDrivingForce, delta
  real(8), dimension(:), allocatable :: sensitivityNextMolFraction, sensitivityNextMolesSpecies
  real(8), dimension(:), allocatable :: sensitivityNextChemicalPotential, sensitivityNextDrivingForce
  integer, dimension(256) :: tokenBase, tokenPlus, tokenMinus
  real(8) :: logPlus, logMinus
  logical :: passed

  bestElementColumn = 0
  speciesIndex = 0
  elementAtomicNumber = 0
  call SolveBase
  passed = INFOThermo == 0
  if (.NOT. passed) call Finish

  baseMolFraction = dMolFraction
  baseMolesSpecies = dMolesSpecies
  baseChemicalPotential = dChemicalPotential
  baseDrivingForce = dDrivingForceSoln
  call GetEquilibriumAssemblageToken(tokenBase,SIZE(tokenBase),nTokenBase,INFO)
  passed = passed .AND. INFO == 0

  call ComputeEquilibriumSensitivity(INFO,rcond)
  passed = passed .AND. INFO == 0 .AND. rcond >= 1D-12
  if (.NOT. passed) call Finish
  nSolSpecies = SIZE(iEquilibriumSensitivitySpecies)
  do parameter = 1, nEquilibriumSensitivityParameters
    do j = 1, nElements
      balanceDerivative = 0D0
      do i = 1, nSolSpecies
        k = iEquilibriumSensitivitySpecies(i)
        balanceDerivative = balanceDerivative + dStoichSpecies(k,j) / &
          DFLOAT(iParticlesPerMole(k)) * EXP(dEquilibriumSensitivityLogMoles(i)) * &
          dEquilibriumSensitivity(i,parameter)
      end do
      do i = 1, nConPhases
        k = iAssemblage(i)
        balanceDerivative = balanceDerivative + dStoichSpecies(k,j) * &
          EXP(dEquilibriumSensitivityLogMoles(nSolSpecies+i)) * &
          dEquilibriumSensitivity(nSolSpecies+i,parameter)
      end do
      expectedDerivative = 0D0
      if (parameter == j + 2) expectedDerivative = 1D0
      passed = passed .AND. ABS(balanceDerivative-expectedDerivative) <= 1D-7
    end do
  end do
  if (.NOT. passed) call Finish
  largestDerivative = 0D0
  do elementColumn = 1, nElements
    if (MAXVAL(ABS(dEquilibriumSensitivity(1:SIZE(iEquilibriumSensitivitySpecies),2+elementColumn))) > &
        largestDerivative) then
      largestDerivative = &
        MAXVAL(ABS(dEquilibriumSensitivity(1:SIZE(iEquilibriumSensitivitySpecies),2+elementColumn)))
      speciesIndex = MAXLOC(ABS(dEquilibriumSensitivity( &
        1:SIZE(iEquilibriumSensitivitySpecies),2+elementColumn)),DIM=1)
      bestElementColumn = elementColumn
    end if
  end do
  passed = passed .AND. largestDerivative > 0D0
  if (.NOT. passed) call Finish
  select case (TRIM(cElementName(bestElementColumn)))
  case ('Mo')
    elementAtomicNumber = 42
  case ('Ru')
    elementAtomicNumber = 44
  case default
    passed = .FALSE.
  end select
  if (.NOT. passed) call Finish
  baseAmount = dElementMass(elementAtomicNumber)
  speciesIndex = iEquilibriumSensitivitySpecies(speciesIndex)
  implicitDerivative = dEquilibriumSensitivity( &
    FINDLOC(iEquilibriumSensitivitySpecies,speciesIndex,DIM=1),2+bestElementColumn)

  allocate(delta(nEquilibriumSensitivityParameters))
  delta = 0D0
  delta(1) = 0.1D0
  call ApplyEquilibriumSensitivity(delta,SIZE(delta),INFO)
  passed = passed .AND. INFO == 0
  call RestoreEquilibriumSensitivity(restoreInfo)
  passed = passed .AND. restoreInfo == 0
  passed = passed .AND. MAXVAL(ABS(dMolFraction-baseMolFraction)) <= 1D-14
  passed = passed .AND. MAXVAL(ABS(dMolesSpecies-baseMolesSpecies)) <= 1D-14
  passed = passed .AND. MAXVAL(ABS(dChemicalPotential-baseChemicalPotential)) <= 1D-12
  passed = passed .AND. MAXVAL(ABS(dDrivingForceSoln-baseDrivingForce)) <= 1D-12

  step = 1D-4
  call SolveAtElementAmount(elementAtomicNumber,baseAmount+step,speciesIndex, &
                            logPlus,tokenPlus,nTokenPlus,INFO)
  passed = passed .AND. INFO == 0
  call SolveAtElementAmount(elementAtomicNumber,baseAmount-step,speciesIndex, &
                            logMinus,tokenMinus,nTokenMinus,INFO)
  passed = passed .AND. INFO == 0
  passed = passed .AND. nTokenPlus == nTokenBase .AND. nTokenMinus == nTokenBase
  if (passed) then
    passed = ALL(tokenPlus(1:nTokenBase) == tokenBase(1:nTokenBase)) .AND. &
             ALL(tokenMinus(1:nTokenBase) == tokenBase(1:nTokenBase))
    coarseDerivative = (logPlus-logMinus)/(2D0*step)
    call SolveAtElementAmount(elementAtomicNumber,baseAmount+step/2D0,speciesIndex, &
                              logPlusHalf,tokenPlus,nTokenPlus,INFO)
    passed = passed .AND. INFO == 0
    call SolveAtElementAmount(elementAtomicNumber,baseAmount-step/2D0,speciesIndex, &
                              logMinusHalf,tokenMinus,nTokenMinus,INFO)
    passed = passed .AND. INFO == 0
    passed = passed .AND. nTokenPlus == nTokenBase .AND. nTokenMinus == nTokenBase
    if (passed) passed = ALL(tokenPlus(1:nTokenBase) == tokenBase(1:nTokenBase)) .AND. &
                         ALL(tokenMinus(1:nTokenBase) == tokenBase(1:nTokenBase))
    exactDerivative = (logPlusHalf-logMinusHalf)/step
    passed = passed .AND. ABS(coarseDerivative-exactDerivative) <= &
             1D-3*MAX(1D0,ABS(exactDerivative))
    passed = passed .AND. ABS(implicitDerivative-exactDerivative) <= &
             1D-3*MAX(1D0,ABS(exactDerivative))
    coarseError = ABS(logPlus-(LOG(baseMolesSpecies(speciesIndex))+implicitDerivative*step))
    fineError = ABS(logPlusHalf-(LOG(baseMolesSpecies(speciesIndex))+implicitDerivative*step/2D0))
    if (coarseError > 1D-10 .AND. fineError > 1D-12) &
      passed = passed .AND. coarseError/fineError >= 2.5D0 .AND. coarseError/fineError <= 6D0
  end if

  ! Reproduce the gas-dominant Li-F state whose inactive SUBQ phase used to be
  ! subminimized by sensitivity evaluation and contaminate the following solve.
  call SolveExcessFluorine
  passed = passed .AND. INFOThermo == 0
  if (.NOT. passed) call Finish
  baseMolFraction = dMolFraction
  baseMolesSpecies = dMolesSpecies
  baseChemicalPotential = dChemicalPotential
  baseDrivingForce = dDrivingForceSoln
  call SaveReinitData
  call ComputeEquilibriumSensitivity(INFO,rcond)
  passed = passed .AND. INFO == 0
  passed = passed .AND. MAXVAL(ABS(dMolFraction-baseMolFraction)) <= 1D-14
  passed = passed .AND. MAXVAL(ABS(dMolesSpecies-baseMolesSpecies)) <= 1D-14
  passed = passed .AND. MAXVAL(ABS(dChemicalPotential-baseChemicalPotential)) <= 1D-12
  passed = passed .AND. MAXVAL(ABS(dDrivingForceSoln-baseDrivingForce)) <= 1D-12
  if (.NOT. passed) call Finish
  if (allocated(delta)) deallocate(delta)
  allocate(delta(nEquilibriumSensitivityParameters))
  delta = 0D0
  delta(1) = 0.1D0
  call ApplyEquilibriumSensitivity(delta,SIZE(delta),INFO)
  passed = passed .AND. INFO == 0
  call RestoreEquilibriumSensitivity(restoreInfo)
  passed = passed .AND. restoreInfo == 0
  passed = passed .AND. MAXVAL(ABS(dMolFraction-baseMolFraction)) <= 1D-14
  passed = passed .AND. MAXVAL(ABS(dMolesSpecies-baseMolesSpecies)) <= 1D-14
  passed = passed .AND. MAXVAL(ABS(dChemicalPotential-baseChemicalPotential)) <= 1D-12
  passed = passed .AND. MAXVAL(ABS(dDrivingForceSoln-baseDrivingForce)) <= 1D-12
  call SolveNextExcessFluorine
  passed = passed .AND. INFOThermo == 0
  if (.NOT. passed) call Finish
  sensitivityNextMolFraction = dMolFraction
  sensitivityNextMolesSpecies = dMolesSpecies
  sensitivityNextChemicalPotential = dChemicalPotential
  sensitivityNextDrivingForce = dDrivingForceSoln

  call SolveExcessFluorine
  passed = passed .AND. INFOThermo == 0
  if (.NOT. passed) call Finish
  call SaveReinitData
  call SolveNextExcessFluorine
  passed = passed .AND. INFOThermo == 0
  passed = passed .AND. MAXVAL(ABS(dMolFraction-sensitivityNextMolFraction)) <= 1D-12
  passed = passed .AND. MAXVAL(ABS(dMolesSpecies-sensitivityNextMolesSpecies)) <= 1D-12
  passed = passed .AND. MAXVAL(ABS(dChemicalPotential-sensitivityNextChemicalPotential)) <= 1D-10
  passed = passed .AND. MAXVAL(ABS(dDrivingForceSoln-sensitivityNextDrivingForce)) <= 1D-10

  ! Independently qualify an active ideal-gas model against exact GEM differences.
  call SolveIdealGas
  passed = passed .AND. INFOThermo == 0
  if (.NOT. passed) call Finish
  baseMolesSpecies = dMolesSpecies
  call GetEquilibriumAssemblageToken(tokenBase,SIZE(tokenBase),nTokenBase,INFO)
  passed = passed .AND. INFO == 0
  call ComputeEquilibriumSensitivity(INFO,rcond)
  passed = passed .AND. INFO == 0 .AND. rcond >= 1D-12
  if (.NOT. passed) call Finish
  nSolSpecies = SIZE(iEquilibriumSensitivitySpecies)
  do parameter = 1, nEquilibriumSensitivityParameters
    do j = 1, nElements
      balanceDerivative = 0D0
      do i = 1, nSolSpecies
        k = iEquilibriumSensitivitySpecies(i)
        balanceDerivative = balanceDerivative + dStoichSpecies(k,j) / &
          DFLOAT(iParticlesPerMole(k)) * EXP(dEquilibriumSensitivityLogMoles(i)) * &
          dEquilibriumSensitivity(i,parameter)
      end do
      expectedDerivative = 0D0
      if (parameter == j + 2) expectedDerivative = 1D0
      passed = passed .AND. ABS(balanceDerivative-expectedDerivative) <= 1D-7
    end do
  end do
  if (.NOT. passed) call Finish
  largestDerivative = 0D0
  do elementColumn = 1, nElements
    if (MAXVAL(ABS(dEquilibriumSensitivity(1:nSolSpecies,2+elementColumn))) > &
        largestDerivative) then
      largestDerivative = MAXVAL(ABS(dEquilibriumSensitivity(1:nSolSpecies,2+elementColumn)))
      speciesIndex = MAXLOC(ABS(dEquilibriumSensitivity(1:nSolSpecies,2+elementColumn)),DIM=1)
      bestElementColumn = elementColumn
    end if
  end do
  select case (TRIM(cElementName(bestElementColumn)))
  case ('C')
    elementAtomicNumber = 6
  case ('O')
    elementAtomicNumber = 8
  case default
    passed = .FALSE.
  end select
  if (.NOT. passed) call Finish
  baseAmount = dElementMass(elementAtomicNumber)
  speciesIndex = iEquilibriumSensitivitySpecies(speciesIndex)
  implicitDerivative = dEquilibriumSensitivity( &
    FINDLOC(iEquilibriumSensitivitySpecies,speciesIndex,DIM=1),2+bestElementColumn)
  step = 1D-4
  call SolveIdealGasAtElementAmount(elementAtomicNumber,baseAmount+step,speciesIndex, &
                                    logPlus,tokenPlus,nTokenPlus,INFO)
  passed = passed .AND. INFO == 0
  call SolveIdealGasAtElementAmount(elementAtomicNumber,baseAmount-step,speciesIndex, &
                                    logMinus,tokenMinus,nTokenMinus,INFO)
  passed = passed .AND. INFO == 0
  coarseDerivative = (logPlus-logMinus)/(2D0*step)
  call SolveIdealGasAtElementAmount(elementAtomicNumber,baseAmount+step/2D0,speciesIndex, &
                                    logPlusHalf,tokenPlus,nTokenPlus,INFO)
  passed = passed .AND. INFO == 0
  call SolveIdealGasAtElementAmount(elementAtomicNumber,baseAmount-step/2D0,speciesIndex, &
                                    logMinusHalf,tokenMinus,nTokenMinus,INFO)
  passed = passed .AND. INFO == 0
  passed = passed .AND. nTokenPlus == nTokenBase .AND. nTokenMinus == nTokenBase
  if (passed) passed = ALL(tokenPlus(1:nTokenBase) == tokenBase(1:nTokenBase)) .AND. &
                       ALL(tokenMinus(1:nTokenBase) == tokenBase(1:nTokenBase))
  exactDerivative = (logPlusHalf-logMinusHalf)/step
  passed = passed .AND. ABS(coarseDerivative-exactDerivative) <= &
           1D-3*MAX(1D0,ABS(exactDerivative))
  passed = passed .AND. ABS(implicitDerivative-exactDerivative) <= &
           1D-3*MAX(1D0,ABS(exactDerivative))

  call Finish

contains

  subroutine SolveBase
    dElementMass = 0D0
    dElementMass(42) = 0.8D0
    dElementMass(44) = 0.2D0
    dTemperature = 2250D0
    dPressure = 1D0
    cInputUnitTemperature = 'K'
    cInputUnitPressure = 'atm'
    cInputUnitMass = 'moles'
    cThermoFileName = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'
    iPrintResultsMode = 0
    call ParseCSDataFile(cThermoFileName)
    call Thermochimica
  end subroutine SolveBase

  subroutine SolveAtElementAmount(element, amount, species, logMoles, token, nToken, solveInfo)
    real(8), intent(in) :: amount
    integer, intent(in) :: element, species
    real(8), intent(out) :: logMoles
    integer, dimension(:), intent(out) :: token
    integer, intent(out) :: nToken, solveInfo

    call ResetThermoAll
    dElementMass = 0D0
    dElementMass(42) = 0.8D0
    dElementMass(44) = 0.2D0
    dElementMass(element) = amount
    dTemperature = 2250D0
    dPressure = 1D0
    cInputUnitTemperature = 'K'
    cInputUnitPressure = 'atm'
    cInputUnitMass = 'moles'
    cThermoFileName = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'
    iPrintResultsMode = 0
    call ParseCSDataFile(cThermoFileName)
    call Thermochimica
    solveInfo = INFOThermo
    if (solveInfo == 0) then
      logMoles = LOG(dMolesSpecies(species))
      call GetEquilibriumAssemblageToken(token,SIZE(token),nToken,solveInfo)
    end if
  end subroutine SolveAtElementAmount

  subroutine SolveExcessFluorine
    call ResetThermoAll
    dElementMass = 0D0
    dElementMass(3) = 1.71906013D1
    dElementMass(9) = 3.87715483D1
    dTemperature = 932.5D0
    dPressure = 207500D0
    cInputUnitTemperature = 'K'
    cInputUnitPressure = 'Pa'
    cInputUnitMass = 'moles'
    cThermoFileName = DATA_DIRECTORY // &
      '../../../examples/MSRE/MSTDB-TC_V3.0_Fluorides_No_Functions_8-2.dat'
    iPrintResultsMode = 0
    call ParseCSDataFile(cThermoFileName)
    call Thermochimica
  end subroutine SolveExcessFluorine

  subroutine SolveNextExcessFluorine
    dElementMass = 0D0
    dElementMass(3) = 1.71906013D1 * 1.0001D0
    dElementMass(9) = 3.87715483D1 * 0.9999D0
    dTemperature = 932.6D0
    dPressure = 207510D0
    lReinitRequested = .TRUE.
    call Thermochimica
  end subroutine SolveNextExcessFluorine

  subroutine SolveIdealGas
    call ResetThermoAll
    dElementMass = 0D0
    dElementMass(6) = 1D0
    dElementMass(8) = 2D0
    dTemperature = 4000D0
    dPressure = 1D0
    cInputUnitTemperature = 'K'
    cInputUnitPressure = 'atm'
    cInputUnitMass = 'moles'
    cThermoFileName = DATA_DIRECTORY // 'CO.dat'
    iPrintResultsMode = 0
    call ParseCSDataFile(cThermoFileName)
    call Thermochimica
  end subroutine SolveIdealGas

  subroutine SolveIdealGasAtElementAmount(element, amount, species, logMoles, token, nToken, solveInfo)
    real(8), intent(in) :: amount
    integer, intent(in) :: element, species
    real(8), intent(out) :: logMoles
    integer, dimension(:), intent(out) :: token
    integer, intent(out) :: nToken, solveInfo

    call ResetThermoAll
    dElementMass = 0D0
    dElementMass(6) = 1D0
    dElementMass(8) = 2D0
    dElementMass(element) = amount
    dTemperature = 4000D0
    dPressure = 1D0
    cInputUnitTemperature = 'K'
    cInputUnitPressure = 'atm'
    cInputUnitMass = 'moles'
    cThermoFileName = DATA_DIRECTORY // 'CO.dat'
    iPrintResultsMode = 0
    call ParseCSDataFile(cThermoFileName)
    call Thermochimica
    solveInfo = INFOThermo
    if (solveInfo == 0) then
      logMoles = LOG(dMolesSpecies(species))
      call GetEquilibriumAssemblageToken(token,SIZE(token),nToken,solveInfo)
    end if
  end subroutine SolveIdealGasAtElementAmount

  subroutine Finish
    if (passed) then
      print *, 'TestSensitivity01: PASS'
      call ResetThermoAll
      call EXIT(0)
    else
      print *, 'TestSensitivity01: FAIL <---'
      call ResetThermoAll
      call EXIT(1)
    end if
  end subroutine Finish

end program TestSensitivity01
