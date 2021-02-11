
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

  ! Initialize storage variables if not allocated already
  if (.NOT. lReinitAvailable) then
    allocate(dMolesPhase_Old(nElements),dChemicalPotential_Old(nSpecies),dElementPotential_Old(nElements),&
    dMolFraction_Old(nSpecies))
    allocate(iAssemblage_Old(nElements))
  endif

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
