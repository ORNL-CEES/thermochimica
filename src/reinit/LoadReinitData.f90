
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    LoadReinitData.f90
    !> \brief   Load data for reiniting calculation from previous results.
    !> \author  M. Poschmann
    !> \sa      SaveReinitData.f90
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


subroutine LoadReinitData

  USE ModuleThermo
  USE ModuleReinit
  USE ModuleThermoIO

  implicit none

  integer                              :: i, j
  integer,     dimension(0:168)        :: iElementsUsed


  lReinitLoaded = .FALSE.

  ! Initialize storage variables if not allocated already
  if (.NOT. lReinitAvailable) then
    print *, 'Reinit requested but data not available'
    return
  endif

  ! Check that the number of elements hasn't changed. If it has, don't load data.
  i = SIZE(iAssemblage_Old)
  if (i /= nElements) return
  ! Also check that the elements involved themselves are the same.
  iElementsUsed = min(ceiling(dElementMass),1)
  do j = 0, (nElementsPT)
    if (iElementsUsed(j) /= iElementsUsed_Old(j)) return
  enddo

  ! Save old chemical potential data
  if (.NOT. allocated(dChemicalPotential)) allocate(dChemicalPotential(nSpecies))
  if (.NOT. allocated(dElementPotential)) allocate(dElementPotential(nElements))
  dChemicalPotential  = dChemicalPotential_Old
  dElementPotential   = dElementPotential_Old
  ! Save old phase data
  if (.NOT. allocated(iAssemblage)) allocate(iAssemblage(nElements))
  if (.NOT. allocated(dMolesPhase)) allocate(dMolesPhase(nElements))
  if (.NOT. allocated(iAssemblageBest)) allocate(iAssemblageBest(nElements))
  if (.NOT. allocated(dMolesPhaseBest)) allocate(dMolesPhaseBest(nElements))
  if (.NOT. allocated(dMolFraction)) allocate(dMolFraction(nSpecies))
  if (.NOT. allocated(dMolFractionBest)) allocate(dMolFractionBest(nSpecies))
  iAssemblage         = iAssemblage_Old
  dMolesPhase         = dMolesPhase_Old * dNormalizeInput
  dMolFraction        = dMolFraction_Old

  lReinitLoaded = .TRUE.

end subroutine LoadReinitData
