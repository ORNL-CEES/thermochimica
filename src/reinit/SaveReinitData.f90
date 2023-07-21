
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    SaveReinitData.f90
    !> \brief   Save data for reiniting calculation from previous results.
    !> \author  M. Poschmann
    !> \sa      LoadReinitData.f90
    !
    !
    ! References:
    ! ===========
    !
    ! Revisions:
    ! ==========
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   30/11/2018      M. Poschmann         Create file.
    !
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this subroutine is to save all pertinent data such that a new call to Thermochimica
    !! may be reinited from that data.
    !
    ! Pertinent variables:
    ! ====================
    !> \param   lReinit        A logical indicating whether reinit data is available.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine SaveReinitData

  USE ModuleThermo
  USE ModuleReinit
  USE ModuleThermoIO

  implicit none

  ! If error has occurred, do not save
  if (INFOThermo /= 0) return
  ! If a calculation has not run, do not save
  if (.NOT. allocated(iAssemblage)) return
  ! If no solution phases, do not save
  ! (this seems to cause errors for new conditions where solution phase should be stable)
  if (nSolnPhases == 0) return

  ! Initialize storage variables if not allocated already at the right size

  ! nElements
  if (allocated(dMolesPhase_Old) .AND. sizeof(dMolesPhase_Old) < nElements) &
    deallocate(dMolesPhase_Old)
  if (.NOT. allocated(dMolesPhase_Old)) &
    allocate(dMolesPhase_Old(nElements))

  if (allocated(dElementPotential_Old) .AND. sizeof(dElementPotential_Old) < nElements) &
    deallocate(dElementPotential_Old)
  if (.NOT. allocated(dElementPotential_Old)) &
    allocate(dElementPotential_Old(nElements))

  if (allocated(iAssemblage_Old) .AND. sizeof(iAssemblage_Old) < nElements) &
    deallocate(iAssemblage_Old)
  if (.NOT. allocated(iAssemblage_Old)) &
    allocate(iAssemblage_Old(nElements))

  !nSpecies
  if (allocated(dChemicalPotential_Old) .AND. sizeof(dChemicalPotential_Old) < nSpecies) &
    deallocate(dChemicalPotential_Old)
  if (.NOT. allocated(dChemicalPotential_Old)) &
    allocate(dChemicalPotential_Old(nSpecies))

  if (allocated(dMolFraction_Old) .AND. sizeof(dMolFraction_Old) < nSpecies) &
    deallocate(dMolFraction_Old)
  if (.NOT. allocated(dMolFraction_Old)) &
    allocate(dMolFraction_Old(nSpecies))

  ! Save old chemical potential data
  dChemicalPotential_Old  = dChemicalPotential
  dElementPotential_Old   = dElementPotential
  ! Save old phase data
  iAssemblage_Old         = iAssemblage
  dMolesPhase_Old         = dMolesPhase !/ dNormalizeInput
  dMolFraction_Old        = dMolFraction
  ! Save which elements were included
  iElementsUsed_Old       = min(ceiling(dElementMass),1)

  ! Set reinit data flag to true
  lReinitAvailable = .TRUE.
  return

end subroutine SaveReinitData
