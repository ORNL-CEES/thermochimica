
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ResetReinit.f90
    !> \brief   Deallocate allocatable variables used by the ModuleReinit.f90
    !> \author  M. Poschmann
    !> \date    Nov. 30, 2018
    !> \sa      ModuleReinit.f90
    !> \sa      ResetThermo.f90
    !> \sa      ResetAll.f90
    !
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer         Description of change
    !    ----          ----------         ---------------------
    !    30/11/2018    M. Poschmann        Original code
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this subroutine is to attempt to gracefully exit Thermochimica.  Allocatable
    !! arrays related to reiniting are deallocated and memory is stored for output to external packages.
    !
    ! Pertinent variables:
    ! ====================
    ! INFO                  An error is returned if deallocation is unsuccessful.
    ! INFOThermo            An integer scalar identifying whether the program exits successfully or if
    !                       it encounters an error.  A description for each error is given in ThermoDebug.f90.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ResetReinit

    USE ModuleReinit
    USE ModuleThermoIO, ONLY: INFOThermo, lReinitAvailable, lReinitLoaded

    implicit none

    integer :: i, INFO

    lReinitAvailable = .FALSE.
    lReinitLoaded = .FALSE.

    ! Initialize variables:
    i = 0

    if (allocated(dMolesPhase_Old)) deallocate(dMolesPhase_Old, STAT = INFO)
    i = i + INFO
    if (allocated(dChemicalPotential_Old)) deallocate(dChemicalPotential_Old, STAT = INFO)
    i = i + INFO
    if (allocated(dMolFraction_Old)) deallocate(dMolFraction_Old, STAT = INFO)
    i = i + INFO
    if (allocated(dElementPotential_Old)) deallocate(dElementPotential_Old, STAT = INFO)
    i = i + INFO
    if (allocated(iAssemblage_Old)) deallocate(iAssemblage_Old, STAT = INFO)
    i = i + INFO

    ! Return an INFOThermo if deallocation of any of the allocatable variables failed:
    if (i > 0) then
        INFOThermo = 15
    end if

    return

end subroutine ResetReinit
