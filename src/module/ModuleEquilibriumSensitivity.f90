! Fixed-assemblage equilibrium sensitivities used by coupled applications.
module ModuleEquilibriumSensitivity

  implicit none

  SAVE

  logical :: lEquilibriumSensitivityAvailable = .FALSE.
  integer :: nEquilibriumSensitivityParameters = 0
  integer :: nEquilibriumSensitivityVariables = 0
  real(8) :: dEquilibriumSensitivityRcond = 0D0
  real(8) :: dEquilibriumSensitivityRawRcond = 0D0
  real(8) :: dEquilibriumSensitivityBaseResidual = 0D0
  real(8) :: dEquilibriumSensitivityLinearResidual = 0D0
  real(8) :: dEquilibriumSensitivityBackwardError = 0D0
  real(8) :: dEquilibriumSensitivityMinActiveFraction = 0D0
  real(8) :: dEquilibriumSensitivityMinActiveMoles = 0D0
  real(8) :: dEquilibriumSensitivityTemperature = 0D0
  real(8) :: dEquilibriumSensitivityPressure = 0D0

  integer, dimension(:), allocatable :: iEquilibriumSensitivitySpecies
  real(8), dimension(:), allocatable :: dEquilibriumSensitivityElements
  real(8), dimension(:), allocatable :: dEquilibriumSensitivityLogMoles
  real(8), dimension(:), allocatable :: dEquilibriumSensitivityElementPotential
  real(8), dimension(:,:), allocatable :: dEquilibriumSensitivity

  ! The sensitivity evaluator uses Thermochimica's module state as workspace.  Keep an exact
  ! snapshot so numerical perturbations and trial evaluations cannot seed a later GEM solve.
  integer :: iEquilibriumSensitivityInfoThermo = 0
  real(8) :: dEquilibriumSensitivitySnapshotTemperature = 0D0
  real(8) :: dEquilibriumSensitivitySnapshotPressure = 0D0
  real(8) :: dEquilibriumSensitivitySnapshotGibbsEnergySys = 0D0
  real(8), dimension(:), allocatable :: dSnapshotStdGibbsEnergy, dSnapshotGibbsSolnPhase
  real(8), dimension(:), allocatable :: dSnapshotMolesSpecies, dSnapshotMagGibbsEnergy
  real(8), dimension(:), allocatable :: dSnapshotChemicalPotential, dSnapshotExcessGibbsParam
  real(8), dimension(:), allocatable :: dSnapshotElementPotential, dSnapshotMolesPhase
  real(8), dimension(:), allocatable :: dSnapshotMolesElement, dSnapshotMolFraction
  real(8), dimension(:), allocatable :: dSnapshotSumMolFractionSoln, dSnapshotDrivingForceSoln
  real(8), dimension(:), allocatable :: dSnapshotPartialExcessGibbs
  real(8), dimension(:,:), allocatable :: dSnapshotEffStoichSolnPhase
  real(8), dimension(:,:,:), allocatable :: dSnapshotSiteFraction

contains

  subroutine ResetEquilibriumSensitivity

    implicit none

    lEquilibriumSensitivityAvailable = .FALSE.
    nEquilibriumSensitivityParameters = 0
    nEquilibriumSensitivityVariables = 0
    dEquilibriumSensitivityRcond = 0D0
    dEquilibriumSensitivityRawRcond = 0D0
    dEquilibriumSensitivityBaseResidual = 0D0
    dEquilibriumSensitivityLinearResidual = 0D0
    dEquilibriumSensitivityBackwardError = 0D0
    dEquilibriumSensitivityMinActiveFraction = 0D0
    dEquilibriumSensitivityMinActiveMoles = 0D0

    if (allocated(iEquilibriumSensitivitySpecies)) deallocate(iEquilibriumSensitivitySpecies)
    if (allocated(dEquilibriumSensitivityElements)) deallocate(dEquilibriumSensitivityElements)
    if (allocated(dEquilibriumSensitivityLogMoles)) deallocate(dEquilibriumSensitivityLogMoles)
    if (allocated(dEquilibriumSensitivityElementPotential)) deallocate(dEquilibriumSensitivityElementPotential)
    if (allocated(dEquilibriumSensitivity)) deallocate(dEquilibriumSensitivity)

    call ResetEquilibriumSensitivitySnapshot

  end subroutine ResetEquilibriumSensitivity

  subroutine ResetEquilibriumSensitivitySnapshot

    implicit none

    if (allocated(dSnapshotStdGibbsEnergy)) deallocate(dSnapshotStdGibbsEnergy)
    if (allocated(dSnapshotGibbsSolnPhase)) deallocate(dSnapshotGibbsSolnPhase)
    if (allocated(dSnapshotMolesSpecies)) deallocate(dSnapshotMolesSpecies)
    if (allocated(dSnapshotMagGibbsEnergy)) deallocate(dSnapshotMagGibbsEnergy)
    if (allocated(dSnapshotChemicalPotential)) deallocate(dSnapshotChemicalPotential)
    if (allocated(dSnapshotExcessGibbsParam)) deallocate(dSnapshotExcessGibbsParam)
    if (allocated(dSnapshotElementPotential)) deallocate(dSnapshotElementPotential)
    if (allocated(dSnapshotMolesPhase)) deallocate(dSnapshotMolesPhase)
    if (allocated(dSnapshotMolesElement)) deallocate(dSnapshotMolesElement)
    if (allocated(dSnapshotMolFraction)) deallocate(dSnapshotMolFraction)
    if (allocated(dSnapshotSumMolFractionSoln)) deallocate(dSnapshotSumMolFractionSoln)
    if (allocated(dSnapshotDrivingForceSoln)) deallocate(dSnapshotDrivingForceSoln)
    if (allocated(dSnapshotPartialExcessGibbs)) deallocate(dSnapshotPartialExcessGibbs)
    if (allocated(dSnapshotEffStoichSolnPhase)) deallocate(dSnapshotEffStoichSolnPhase)
    if (allocated(dSnapshotSiteFraction)) deallocate(dSnapshotSiteFraction)

  end subroutine ResetEquilibriumSensitivitySnapshot

  subroutine CaptureEquilibriumSensitivitySnapshot

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: dTemperature, dPressure, dGibbsEnergySys, INFOThermo
    USE ModuleGEMSolver

    implicit none

    call ResetEquilibriumSensitivitySnapshot
    dEquilibriumSensitivitySnapshotTemperature = dTemperature
    dEquilibriumSensitivitySnapshotPressure = dPressure
    dEquilibriumSensitivitySnapshotGibbsEnergySys = dGibbsEnergySys
    iEquilibriumSensitivityInfoThermo = INFOThermo
    if (allocated(dStdGibbsEnergy)) dSnapshotStdGibbsEnergy = dStdGibbsEnergy
    if (allocated(dGibbsSolnPhase)) dSnapshotGibbsSolnPhase = dGibbsSolnPhase
    if (allocated(dMolesSpecies)) dSnapshotMolesSpecies = dMolesSpecies
    if (allocated(dMagGibbsEnergy)) dSnapshotMagGibbsEnergy = dMagGibbsEnergy
    if (allocated(dChemicalPotential)) dSnapshotChemicalPotential = dChemicalPotential
    if (allocated(dExcessGibbsParam)) dSnapshotExcessGibbsParam = dExcessGibbsParam
    if (allocated(dElementPotential)) dSnapshotElementPotential = dElementPotential
    if (allocated(dMolesPhase)) dSnapshotMolesPhase = dMolesPhase
    if (allocated(dMolesElement)) dSnapshotMolesElement = dMolesElement
    if (allocated(dMolFraction)) dSnapshotMolFraction = dMolFraction
    if (allocated(dSumMolFractionSoln)) dSnapshotSumMolFractionSoln = dSumMolFractionSoln
    if (allocated(dDrivingForceSoln)) dSnapshotDrivingForceSoln = dDrivingForceSoln
    if (allocated(dPartialExcessGibbs)) dSnapshotPartialExcessGibbs = dPartialExcessGibbs
    if (allocated(dEffStoichSolnPhase)) dSnapshotEffStoichSolnPhase = dEffStoichSolnPhase
    if (allocated(dSiteFraction)) dSnapshotSiteFraction = dSiteFraction

  end subroutine CaptureEquilibriumSensitivitySnapshot

  subroutine RestoreEquilibriumSensitivitySnapshot(INFO)

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: dTemperature, dPressure, dGibbsEnergySys, INFOThermo
    USE ModuleGEMSolver

    implicit none

    integer, intent(out) :: INFO

    INFO = 0
    if (.NOT. allocated(dSnapshotMolesSpecies)) then
      INFO = 1
      return
    end if
    dTemperature = dEquilibriumSensitivitySnapshotTemperature
    dPressure = dEquilibriumSensitivitySnapshotPressure
    dGibbsEnergySys = dEquilibriumSensitivitySnapshotGibbsEnergySys
    INFOThermo = iEquilibriumSensitivityInfoThermo
    if (allocated(dSnapshotStdGibbsEnergy)) dStdGibbsEnergy = dSnapshotStdGibbsEnergy
    if (allocated(dSnapshotGibbsSolnPhase)) dGibbsSolnPhase = dSnapshotGibbsSolnPhase
    if (allocated(dSnapshotMolesSpecies)) dMolesSpecies = dSnapshotMolesSpecies
    if (allocated(dSnapshotMagGibbsEnergy)) dMagGibbsEnergy = dSnapshotMagGibbsEnergy
    if (allocated(dSnapshotChemicalPotential)) dChemicalPotential = dSnapshotChemicalPotential
    if (allocated(dSnapshotExcessGibbsParam)) dExcessGibbsParam = dSnapshotExcessGibbsParam
    if (allocated(dSnapshotElementPotential)) dElementPotential = dSnapshotElementPotential
    if (allocated(dSnapshotMolesPhase)) dMolesPhase = dSnapshotMolesPhase
    if (allocated(dSnapshotMolesElement)) dMolesElement = dSnapshotMolesElement
    if (allocated(dSnapshotMolFraction)) dMolFraction = dSnapshotMolFraction
    if (allocated(dSnapshotSumMolFractionSoln)) dSumMolFractionSoln = dSnapshotSumMolFractionSoln
    if (allocated(dSnapshotDrivingForceSoln)) dDrivingForceSoln = dSnapshotDrivingForceSoln
    if (allocated(dSnapshotPartialExcessGibbs)) dPartialExcessGibbs = dSnapshotPartialExcessGibbs
    if (allocated(dSnapshotEffStoichSolnPhase)) dEffStoichSolnPhase = dSnapshotEffStoichSolnPhase
    if (allocated(dSnapshotSiteFraction)) dSiteFraction = dSnapshotSiteFraction

  end subroutine RestoreEquilibriumSensitivitySnapshot

end module ModuleEquilibriumSensitivity
