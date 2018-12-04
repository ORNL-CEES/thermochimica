
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    SaveRestartData.f90
    !> \brief   Save data for restarting calculation from previous results.
    !> \author  M. Poschmann
    !> \sa      LoadRestartData.f90
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


subroutine SaveRestartData

  USE ModuleThermo
  USE ModuleRestart
  USE ModuleThermoIO

  implicit none

  ! Initialize storage variables if not allocated already
  if (.NOT. lRestartAvailable) then
    allocate(dMolesPhase_Old(nElements),dChemicalPotential_Old(nSpecies),dElementPotential_Old(nElements),&
    dMolFraction_Old(nSpecies))
    allocate(iAssemblage_Old(nElements))
  endif

  ! Save old chemical potential data
  dChemicalPotential_Old  = dChemicalPotential
  dElementPotential_Old   = dElementPotential
  ! Save old phase data
  iAssemblage_Old         = iAssemblage
  dMolesPhase_Old         = dMolesPhase
  dMolFraction_Old        = dMolFraction
  ! Save which elements were included
  iElementsUsed_Old       = min(ceiling(dElementMass),1)

  ! Set restart data flag to true
  lRestartAvailable = .TRUE.
  return

end subroutine SaveRestartData
