
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
  if (.NOT. lRestart) then
    allocate(dMolesPhase_Old(nElements),dChemicalPotential_Old(nSpecies))
    allocate(iPhase_Old(nSpecies),iAssemblage_Old(nElements))
    allocate(cSolnPhaseType_Old(nSolnPhasesSys),cSolnPhaseName_Old(nSolnPhasesSys))
  endif


  ! Save old chemical potential data
  dChemicalPotential_Old  = dChemicalPotential
  ! Save old phase data
  cSolnPhaseName_Old      = cSolnPhaseName
  cSolnPhaseType_Old      = cSolnPhaseType
  iAssemblage_Old         = iAssemblage
  iPhase_Old              = iPhase
  dMolesPhase_Old         = dMolesPhase

  ! Set restart data flag to true
  lRestart = .TRUE.
  return

end subroutine SaveRestartData
