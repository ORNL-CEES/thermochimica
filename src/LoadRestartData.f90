
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    LoadRestartData.f90
    !> \brief   Load data for restarting calculation from previous results.
    !> \author  M. Poschmann
    !> \sa      SaveRestartData.f90
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
    !! may be restarted from that data.
    !
    ! Pertinent variables:
    ! ====================
    !> \param   lRestart        A logical indicating whether restart data is available.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine LoadRestartData

  USE ModuleThermo
  USE ModuleRestart
  USE ModuleThermoIO

  implicit none

  integer                              :: i, j
  integer,     dimension(0:118)        :: iElementsUsed

  ! Initialize storage variables if not allocated already
  if (.NOT. lRestartAvailable) then
    print *, 'Restart requested but data not available'
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
  if (.NOT. allocated(dMolFraction)) allocate(dMolFraction(nSpecies))
  iAssemblage         = iAssemblage_Old
  dMolesPhase         = dMolesPhase_Old
  dMolFraction        = dMolFraction_Old

  lRestartLoaded = .TRUE.

end subroutine LoadRestartData
