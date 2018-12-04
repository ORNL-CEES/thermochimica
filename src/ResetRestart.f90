
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ResetRestart.f90
    !> \brief   Deallocate allocatable variables used by the ModuleRestart.f90
    !> \author  M. Poschmann
    !> \date    Nov. 30, 2018
    !> \sa      ModuleRestart.f90
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
    !! arrays related to restarting are deallocated and memory is stored for output to external packages.
    !
    ! Pertinent variables:
    ! ====================
    ! INFO                  An error is returned if deallocation is unsuccessful.
    ! INFOThermo            An integer scalar identifying whether the program exits successfully or if
    !                       it encounters an error.  A description for each error is given in ThermoDebug.f90.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ResetRestart

    USE ModuleRestart
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer::   i, INFO

    ! Initialize variables:
    i = 0

    if (allocated(dMolesPhase_Old)) then
        ! Deallocate integer arrays from ModuleThermo:
        deallocate (dMolesPhase_Old,dChemicalPotential_Old,dMolFraction_Old,dElementPotential_Old, STAT = INFO)
        i = i + INFO
    end if

    if (allocated(iAssemblage_Old)) then
        ! Deallocate real arrays from ModuleThermo:
        deallocate (iAssemblage_Old,STAT = INFO)
        i = i + INFO
    end if

    if (allocated(cSolnPhaseType_Old)) then
        ! Deallocate character arrays from ModuleThermo:
        deallocate (cSolnPhaseType_Old,cSolnPhaseName_Old, STAT = INFO)
        i = i + INFO
    end if

    ! Return an INFOThermo if deallocation of any of the allocatable variables failed:
    if (i > 0) then
        INFOThermo = 15
    end if

    return

end subroutine ResetRestart
