!-------------------------------------------------------------------------------------------------------------
! Compute and apply first-order sensitivities of a converged, fixed phase assemblage.
!-------------------------------------------------------------------------------------------------------------

subroutine ComputeEquilibriumSensitivity(INFO, dRcondOut)

  USE ModuleThermo
  USE ModuleThermoIO, ONLY: dTemperature, dPressure, INFOThermo
  USE ModuleGEMSolver
  USE ModuleEquilibriumSensitivity

  implicit none

  integer, intent(out) :: INFO
  real(8), intent(out) :: dRcondOut
  integer :: i, j, k, l, iSol, nSolSpecies, nVar, nSensitivityParams, lapackInfo, restoreInfo
  integer, dimension(:), allocatable :: IPIV, IWORK
  real(8) :: h, anorm, rawAnorm, residualNorm
  real(8), external :: DLANGE
  real(8), dimension(:), allocatable :: x, rBase, rPlus, rMinus, rPlusHalf, rMinusHalf
  real(8), dimension(:), allocatable :: work, rowScale, columnScale
  real(8), dimension(:,:), allocatable :: KKT, RHS, originalKKT, originalRHS, rawFactor, denominator

  call ResetEquilibriumSensitivity
  INFO = 0
  dRcondOut = 0D0
  if (INFOThermo /= 0 .OR. .NOT. lConverged) then
    INFO = 1
    return
  end if

  nSolSpecies = 0
  do i = 1, nSolnPhases
    k = -iAssemblage(nElements - i + 1)
    do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
      if (dMolFraction(j) > 1.66D-24) nSolSpecies = nSolSpecies + 1
    end do
  end do
  nVar = nSolSpecies + nConPhases + nElements
  nSensitivityParams = nElements + 2
  if (nVar <= 0) then
    INFO = 2
    return
  end if

  call CaptureEquilibriumSensitivitySnapshot

  allocate(iEquilibriumSensitivitySpecies(nSolSpecies))
  allocate(dEquilibriumSensitivityElements(nElements))
  allocate(dEquilibriumSensitivityLogMoles(nSolSpecies + nConPhases))
  allocate(dEquilibriumSensitivityElementPotential(nElements))
  allocate(dEquilibriumSensitivity(nVar,nSensitivityParams))
  allocate(x(nVar), rBase(nVar), rPlus(nVar), rMinus(nVar), rPlusHalf(nVar), rMinusHalf(nVar))
  allocate(rowScale(nVar), columnScale(nVar))
  allocate(KKT(nVar,nVar), originalKKT(nVar,nVar), rawFactor(nVar,nVar))
  allocate(RHS(nVar,nSensitivityParams), originalRHS(nVar,nSensitivityParams), denominator(nVar,nSensitivityParams))
  allocate(IPIV(nVar), IWORK(nVar), work(4*nVar))

  iSol = 0
  do i = 1, nSolnPhases
    k = -iAssemblage(nElements - i + 1)
    do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
      if (dMolFraction(j) <= 1.66D-24) cycle
      iSol = iSol + 1
      iEquilibriumSensitivitySpecies(iSol) = j
      x(iSol) = LOG(MAX(dMolesSpecies(j), 1D-300))
    end do
  end do
  do i = 1, nConPhases
    x(nSolSpecies+i) = LOG(MAX(dMolesPhase(i), 1D-300))
  end do
  x(nSolSpecies+nConPhases+1:nVar) = dElementPotential

  call EquilibriumSensitivityResidual(x,nSolSpecies,nVar,rBase)
  dEquilibriumSensitivityBaseResidual = 0D0
  if (nSolSpecies+nConPhases > 0) &
    dEquilibriumSensitivityBaseResidual = &
      MAXVAL(ABS(rBase(1:nSolSpecies+nConPhases)))
  do j = 1, nElements
    dEquilibriumSensitivityBaseResidual = MAX(dEquilibriumSensitivityBaseResidual, &
      ABS(rBase(nSolSpecies+nConPhases+j))/MAX(1D0,ABS(dMolesElement(j))))
  end do
  if (.NOT. (dEquilibriumSensitivityBaseResidual >= 0D0) .OR. &
      dEquilibriumSensitivityBaseResidual > 10D0*MAX(dTolerance(1),dTolerance(5))) then
    INFO = 5
    call RestoreEquilibriumSensitivitySnapshot(restoreInfo)
    deallocate(x,rBase,rPlus,rMinus,rPlusHalf,rMinusHalf,rowScale,columnScale,KKT,originalKKT, &
               rawFactor,RHS,originalRHS,denominator,IPIV,IWORK,work)
    call ResetEquilibriumSensitivity
    return
  end if

  dEquilibriumSensitivityTemperature = dTemperature
  dEquilibriumSensitivityPressure = dPressure
  dEquilibriumSensitivityElements = dMolesElement
  dEquilibriumSensitivityLogMoles = x(1:nSolSpecies+nConPhases)
  dEquilibriumSensitivityElementPotential = dElementPotential
  if (nSolSpecies > 0) then
    dEquilibriumSensitivityMinActiveFraction = &
      MINVAL(dMolFraction(iEquilibriumSensitivitySpecies))
    dEquilibriumSensitivityMinActiveMoles = &
      MINVAL(dMolesSpecies(iEquilibriumSensitivitySpecies))
  else if (nConPhases > 0) then
    dEquilibriumSensitivityMinActiveFraction = 1D0
    dEquilibriumSensitivityMinActiveMoles = MINVAL(dMolesPhase(1:nConPhases))
  end if
  nEquilibriumSensitivityParameters = nSensitivityParams
  nEquilibriumSensitivityVariables = nVar

  h = EPSILON(1D0)**(1D0/3D0)
  KKT = 0D0

  ! Only the model-dependent constituent columns require numerical differentiation.
  do j = 1, nSolSpecies
    x(j) = x(j) + h
    call EquilibriumSensitivityResidual(x, nSolSpecies, nVar, rPlus)
    x(j) = x(j) - 2D0*h
    call EquilibriumSensitivityResidual(x, nSolSpecies, nVar, rMinus)
    x(j) = x(j) + h
    x(j) = x(j) + h/2D0
    call EquilibriumSensitivityResidual(x, nSolSpecies, nVar, rPlusHalf)
    x(j) = x(j) - h
    call EquilibriumSensitivityResidual(x, nSolSpecies, nVar, rMinusHalf)
    x(j) = x(j) + h/2D0
    KKT(:,j) = (rPlusHalf-rMinusHalf)/h
  end do

  ! Pure-phase log amounts enter only the elemental balances.
  do i = 1, nConPhases
    j = nSolSpecies + i
    k = iAssemblage(i)
    do l = 1, nElements
      KKT(nSolSpecies+nConPhases+l,j) = dStoichSpecies(k,l) * EXP(x(j))
    end do
  end do

  ! Element potentials enter only the stationarity equations.
  do j = 1, nElements
    l = nSolSpecies + nConPhases + j
    do i = 1, nSolSpecies
      k = iEquilibriumSensitivitySpecies(i)
      KKT(i,l) = -dStoichSpecies(k,j) / DFLOAT(iParticlesPerMole(k))
    end do
    do i = 1, nConPhases
      k = iAssemblage(i)
      KKT(nSolSpecies+i,l) = -dStoichSpecies(k,j)
    end do
  end do

  RHS = 0D0
  call EquilibriumSensitivityParameterDerivative(x, nSolSpecies, nVar, 1, h, rPlus, rMinus)
  call EquilibriumSensitivityParameterDerivative(x, nSolSpecies, nVar, 1, h/2D0, rPlusHalf, rMinusHalf)
  RHS(:,1) = -(rPlusHalf-rMinusHalf)/(h*dEquilibriumSensitivityTemperature)
  call EquilibriumSensitivityParameterDerivative(x, nSolSpecies, nVar, 2, h, rPlus, rMinus)
  call EquilibriumSensitivityParameterDerivative(x, nSolSpecies, nVar, 2, h/2D0, rPlusHalf, rMinusHalf)
  RHS(:,2) = -(rPlusHalf-rMinusHalf)/(h*dEquilibriumSensitivityPressure)
  do j = 1, nElements
    RHS(nSolSpecies+nConPhases+j,2+j) = 1D0
  end do

  call RestoreEquilibriumSensitivitySnapshot(restoreInfo)
  originalKKT = KKT
  originalRHS = RHS
  rawFactor = KKT
  rawAnorm = DLANGE('1',nVar,nVar,rawFactor,nVar,work)
  call DGETRF(nVar,nVar,rawFactor,nVar,IPIV,lapackInfo)
  if (lapackInfo == 0) call DGECON('1',nVar,rawFactor,nVar,rawAnorm,dEquilibriumSensitivityRawRcond, &
                                   work,IWORK,lapackInfo)
  do i = 1, nVar
    rowScale(i) = MAXVAL(ABS(KKT(i,:)))
    if (rowScale(i) > 0D0) then
      rowScale(i) = 1D0/rowScale(i)
      KKT(i,:) = rowScale(i)*KKT(i,:)
      RHS(i,:) = rowScale(i)*RHS(i,:)
    end if
  end do
  do j = 1, nVar
    columnScale(j) = MAXVAL(ABS(KKT(:,j)))
    if (columnScale(j) > 0D0) then
      columnScale(j) = 1D0/columnScale(j)
      KKT(:,j) = columnScale(j)*KKT(:,j)
    end if
  end do
  anorm = DLANGE('1',nVar,nVar,KKT,nVar,work)
  call DGETRF(nVar,nVar,KKT,nVar,IPIV,lapackInfo)
  if (lapackInfo == 0) call DGECON('1',nVar,KKT,nVar,anorm,dEquilibriumSensitivityRcond,work,IWORK,lapackInfo)
  if (lapackInfo == 0 .AND. dEquilibriumSensitivityRcond >= 1D-12) &
    call DGETRS('N',nVar,nSensitivityParams,KKT,nVar,IPIV,RHS,nVar,lapackInfo)

  if (lapackInfo == 0 .AND. dEquilibriumSensitivityRcond >= 1D-12) then
    do i = 1, nVar
      RHS(i,:) = columnScale(i)*RHS(i,:)
    end do
  end if

  if (lapackInfo /= 0 .OR. dEquilibriumSensitivityRcond < 1D-12 .OR. ANY(RHS /= RHS)) then
    INFO = 3
  else
    denominator = MATMUL(ABS(originalKKT),ABS(RHS)) + ABS(originalRHS)
    residualNorm = MAXVAL(ABS(MATMUL(originalKKT,RHS) - originalRHS))
    dEquilibriumSensitivityLinearResidual = residualNorm
    dEquilibriumSensitivityBackwardError = MAXVAL(ABS(MATMUL(originalKKT,RHS) - originalRHS) / &
                                                   MAX(denominator,TINY(1D0)))
    if (residualNorm > 1D-7 * MAX(1D0,MAXVAL(ABS(originalRHS)))) then
      INFO = 4
    else
      dEquilibriumSensitivity = RHS
      lEquilibriumSensitivityAvailable = .TRUE.
    end if
  end if
  dRcondOut = dEquilibriumSensitivityRcond

  call RestoreEquilibriumSensitivitySnapshot(restoreInfo)
  if (restoreInfo /= 0) INFO = 6
  deallocate(x,rBase,rPlus,rMinus,rPlusHalf,rMinusHalf,rowScale,columnScale,KKT,originalKKT, &
             rawFactor,RHS,originalRHS,denominator,IPIV,IWORK,work)
  if (INFO /= 0) call ResetEquilibriumSensitivity

end subroutine ComputeEquilibriumSensitivity


subroutine EquilibriumSensitivityResidual(x, nSolSpecies, nVar, residual)

  USE ModuleThermo
  USE ModuleGEMSolver
  USE ModuleEquilibriumSensitivity

  implicit none

  integer, intent(in) :: nSolSpecies, nVar
  real(8), dimension(nVar), intent(in) :: x
  real(8), dimension(nVar), intent(out) :: residual
  integer :: i, j, k
  real(8) :: value

  do i = 1, nSolSpecies
    dMolesSpecies(iEquilibriumSensitivitySpecies(i)) = EXP(x(i))
  end do
  do i = 1, nSolnPhases
    k = -iAssemblage(nElements-i+1)
    dMolesPhase(nElements-i+1) = SUM(dMolesSpecies(nSpeciesPhase(k-1)+1:nSpeciesPhase(k)))
  end do
  do i = 1, nConPhases
    dMolesPhase(i) = EXP(x(nSolSpecies+i))
  end do
  dElementPotential = x(nSolSpecies+nConPhases+1:nSolSpecies+nConPhases+nElements)
  call CompChemicalPotential(.FALSE.)

  do i = 1, nSolSpecies
    k = iEquilibriumSensitivitySpecies(i)
    value = 0D0
    do j = 1, nElements
      value = value + dElementPotential(j)*dStoichSpecies(k,j)
    end do
    residual(i) = dChemicalPotential(k) - value/DFLOAT(iParticlesPerMole(k))
  end do
  do i = 1, nConPhases
    k = iAssemblage(i)
    value = dStdGibbsEnergy(k)
    do j = 1, nElements
      value = value - dElementPotential(j)*dStoichSpecies(k,j)
    end do
    residual(nSolSpecies+i) = value
  end do
  do j = 1, nElements
    value = -dMolesElement(j)
    do i = 1, nSolSpecies
      k = iEquilibriumSensitivitySpecies(i)
      value = value + dMolesSpecies(k)*dStoichSpecies(k,j)/DFLOAT(iParticlesPerMole(k))
    end do
    do i = 1, nConPhases
      k = iAssemblage(i)
      value = value + dMolesPhase(i)*dStoichSpecies(k,j)
    end do
    residual(nSolSpecies+nConPhases+j) = value
  end do

end subroutine EquilibriumSensitivityResidual


subroutine EquilibriumSensitivityParameterDerivative(x, nSolSpecies, nVar, parameter, h, rPlus, rMinus)

  USE ModuleThermoIO, ONLY: dTemperature, dPressure
  USE ModuleEquilibriumSensitivity

  implicit none

  integer, intent(in) :: nSolSpecies, nVar, parameter
  real(8), intent(in) :: h
  real(8), dimension(nVar), intent(in) :: x
  real(8), dimension(nVar), intent(out) :: rPlus, rMinus

  dTemperature = dEquilibriumSensitivityTemperature
  dPressure = dEquilibriumSensitivityPressure
  if (parameter == 1) then
    dTemperature = dEquilibriumSensitivityTemperature*(1D0+h)
  else
    dPressure = dEquilibriumSensitivityPressure*(1D0+h)
  end if
  call CompThermoData
  call EquilibriumSensitivityResidual(x,nSolSpecies,nVar,rPlus)

  dTemperature = dEquilibriumSensitivityTemperature
  dPressure = dEquilibriumSensitivityPressure
  if (parameter == 1) then
    dTemperature = dEquilibriumSensitivityTemperature*(1D0-h)
  else
    dPressure = dEquilibriumSensitivityPressure*(1D0-h)
  end if
  call CompThermoData
  call EquilibriumSensitivityResidual(x,nSolSpecies,nVar,rMinus)

end subroutine EquilibriumSensitivityParameterDerivative


subroutine ApplyEquilibriumSensitivity(delta, nSensitivityParams, INFO)

  USE ModuleThermo
  USE ModuleThermoIO, ONLY: dTemperature, dPressure
  USE ModuleGEMSolver
  USE ModuleEquilibriumSensitivity

  implicit none

  integer, intent(in) :: nSensitivityParams
  integer, intent(out) :: INFO
  real(8), dimension(nSensitivityParams), intent(in) :: delta
  integer :: nSolSpecies, restoreInfo
  real(8), dimension(:), allocatable :: x, residual

  INFO = 0
  if (.NOT. lEquilibriumSensitivityAvailable .OR. &
      nSensitivityParams /= nEquilibriumSensitivityParameters) then
    INFO = 1
    return
  end if
  call RestoreEquilibriumSensitivitySnapshot(restoreInfo)
  if (restoreInfo /= 0) then
    INFO = 3
    return
  end if
  nSolSpecies = SIZE(iEquilibriumSensitivitySpecies)
  allocate(x(nEquilibriumSensitivityVariables),residual(nEquilibriumSensitivityVariables))
  x(1:nSolSpecies+nConPhases) = dEquilibriumSensitivityLogMoles + &
    MATMUL(dEquilibriumSensitivity(1:nSolSpecies+nConPhases,:),delta)
  x(nSolSpecies+nConPhases+1:) = dEquilibriumSensitivityElementPotential + &
    MATMUL(dEquilibriumSensitivity(nSolSpecies+nConPhases+1:,:),delta)

  dTemperature = dEquilibriumSensitivityTemperature + delta(1)
  dPressure = dEquilibriumSensitivityPressure + delta(2)
  dMolesElement = dEquilibriumSensitivityElements + delta(3:)
  if (dTemperature <= 0D0 .OR. dPressure <= 0D0 .OR. ANY(dMolesElement < 0D0)) then
    INFO = 2
    deallocate(x,residual)
    call RestoreEquilibriumSensitivitySnapshot(restoreInfo)
    return
  end if
  if (ANY(x(1:nSolSpecies+nConPhases) > LOG(HUGE(1D0))) .OR. &
      ANY(x(1:nSolSpecies+nConPhases) < LOG(TINY(1D0)))) then
    INFO = 2
    deallocate(x,residual)
    call RestoreEquilibriumSensitivitySnapshot(restoreInfo)
    return
  end if
  call CompThermoData
  call EquilibriumSensitivityResidual(x,nSolSpecies,nEquilibriumSensitivityVariables,residual)
  if ((nSolSpecies > 0 .AND. &
       ANY(dMolFraction(iEquilibriumSensitivitySpecies) <= 1.66D-24)) .OR. &
      MAXVAL(ABS(residual(nSolSpecies+nConPhases+1:))) > &
        1D-7*MAX(1D0,MAXVAL(ABS(dMolesElement)))) then
    INFO = 4
    call RestoreEquilibriumSensitivitySnapshot(restoreInfo)
  end if
  deallocate(x,residual)

end subroutine ApplyEquilibriumSensitivity


subroutine RestoreEquilibriumSensitivity(INFO)

  USE ModuleEquilibriumSensitivity

  implicit none

  integer, intent(out) :: INFO

  call RestoreEquilibriumSensitivitySnapshot(INFO)

end subroutine RestoreEquilibriumSensitivity


subroutine RestoreEquilibriumSensitivityBase

  USE ModuleEquilibriumSensitivity

  implicit none

  integer :: INFO

  call RestoreEquilibriumSensitivitySnapshot(INFO)

end subroutine RestoreEquilibriumSensitivityBase


subroutine GetEquilibriumSensitivityInputs(inputs, nParam, INFO)

  USE ModuleEquilibriumSensitivity

  implicit none

  integer, intent(in) :: nParam
  integer, intent(out) :: INFO
  real(8), dimension(nParam), intent(out) :: inputs

  INFO = 0
  if (.NOT. lEquilibriumSensitivityAvailable .OR. nParam /= nEquilibriumSensitivityParameters) then
    INFO = 1
    inputs = 0D0
    return
  end if
  inputs(1) = dEquilibriumSensitivityTemperature
  inputs(2) = dEquilibriumSensitivityPressure
  inputs(3:) = dEquilibriumSensitivityElements

end subroutine GetEquilibriumSensitivityInputs


subroutine GetEquilibriumSensitivityQualification(qualified, INFO)

  USE ModuleThermo
  USE ModuleThermoIO, ONLY: INFOThermo

  implicit none

  integer, intent(out) :: qualified, INFO
  integer :: k

  qualified = 0
  INFO = 0
  if (INFOThermo /= 0) then
    INFO = 1
    return
  end if
  qualified = 1
  do k = 1, nSolnPhasesSys
    select case (cSolnPhaseType(k))
    case ('IDMX','QKTO')
      continue
    case default
      qualified = 0
    end select
  end do

end subroutine GetEquilibriumSensitivityQualification


subroutine GetEquilibriumSensitivityDiagnostics(diagnostics, nDiagnostics, qualified, INFO)

  USE ModuleThermo
  USE ModuleGEMSolver
  USE ModuleEquilibriumSensitivity

  implicit none

  integer, intent(in) :: nDiagnostics
  integer, intent(out) :: qualified, INFO
  real(8), dimension(nDiagnostics), intent(out) :: diagnostics
  integer :: qualificationInfo

  diagnostics = 0D0
  qualified = 0
  INFO = 0
  if (.NOT. lEquilibriumSensitivityAvailable .OR. nDiagnostics < 7) then
    INFO = 1
    return
  end if
  diagnostics(1) = dEquilibriumSensitivityRcond
  diagnostics(2) = dEquilibriumSensitivityRawRcond
  diagnostics(3) = dEquilibriumSensitivityBaseResidual
  diagnostics(4) = dEquilibriumSensitivityLinearResidual
  diagnostics(5) = dEquilibriumSensitivityBackwardError
  diagnostics(6) = dEquilibriumSensitivityMinActiveFraction
  diagnostics(7) = dEquilibriumSensitivityMinActiveMoles

  ! Models are enabled only after their state-isolation and finite-difference regressions exist.
  call GetEquilibriumSensitivityQualification(qualified,qualificationInfo)
  if (qualificationInfo /= 0) qualified = 0
  if (dEquilibriumSensitivityMinActiveFraction > 0D0 .AND. &
      dEquilibriumSensitivityMinActiveFraction <= 1D-18) qualified = 0
  if (dEquilibriumSensitivityMinActiveMoles > 0D0 .AND. &
      dEquilibriumSensitivityMinActiveMoles <= 1D-20) qualified = 0

end subroutine GetEquilibriumSensitivityDiagnostics


subroutine GetEquilibriumAssemblageToken(token, nTokenCapacity, nToken, INFO)

  USE ModuleThermo
  USE ModuleGEMSolver

  implicit none

  integer, intent(in) :: nTokenCapacity
  integer, intent(out) :: nToken, INFO
  integer, dimension(nTokenCapacity), intent(out) :: token
  integer :: i, j, k, modelCode

  token = 0
  nToken = 0
  INFO = 0
  call AppendToken(1)
  call AppendToken(nConPhases)
  call AppendToken(nSolnPhases)
  do i = 1, nConPhases
    call AppendToken(iAssemblage(i))
  end do
  do i = 1, nSolnPhases
    k = -iAssemblage(nElements-i+1)
    select case (cSolnPhaseType(k))
    case ('IDMX'); modelCode = 1
    case ('QKTO'); modelCode = 2
    case ('RKMP'); modelCode = 3
    case ('RKMPM'); modelCode = 4
    case ('SUBL'); modelCode = 5
    case ('SUBLM'); modelCode = 6
    case ('SUBI'); modelCode = 7
    case ('SUBM'); modelCode = 8
    case ('SUBG'); modelCode = 9
    case ('SUBQ'); modelCode = 10
    case default; modelCode = 0
    end select
    call AppendToken(k)
    call AppendToken(modelCode)
    call AppendToken(MERGE(1,0,lMiscibility(k)))
    call AppendToken(COUNT(dMolFraction(nSpeciesPhase(k-1)+1:nSpeciesPhase(k)) > 1.66D-24))
    do j = nSpeciesPhase(k-1)+1, nSpeciesPhase(k)
      if (dMolFraction(j) > 1.66D-24) call AppendToken(j)
    end do
  end do

contains

  subroutine AppendToken(value)
    integer, intent(in) :: value
    nToken = nToken + 1
    if (nToken > nTokenCapacity) then
      INFO = 1
    else
      token(nToken) = value
    end if
  end subroutine AppendToken

end subroutine GetEquilibriumAssemblageToken


subroutine GetEquilibriumSensitivityState(logAmounts, jacobian, nAmountCapacity, &
                                           nParamCapacity, nAmounts, nParam, INFO)

  USE ModuleEquilibriumSensitivity

  implicit none

  integer, intent(in) :: nAmountCapacity, nParamCapacity
  integer, intent(out) :: nAmounts, nParam, INFO
  real(8), dimension(nAmountCapacity), intent(out) :: logAmounts
  real(8), dimension(nAmountCapacity,nParamCapacity), intent(out) :: jacobian

  logAmounts = 0D0
  jacobian = 0D0
  nAmounts = 0
  nParam = 0
  INFO = 0
  if (.NOT. lEquilibriumSensitivityAvailable) then
    INFO = 1
    return
  end if
  nAmounts = SIZE(dEquilibriumSensitivityLogMoles)
  nParam = nEquilibriumSensitivityParameters
  if (nAmounts > nAmountCapacity .OR. nParam > nParamCapacity) then
    INFO = 2
    return
  end if
  logAmounts(1:nAmounts) = dEquilibriumSensitivityLogMoles
  jacobian(1:nAmounts,1:nParam) = dEquilibriumSensitivity(1:nAmounts,1:nParam)

end subroutine GetEquilibriumSensitivityState
