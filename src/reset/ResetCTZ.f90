
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ResetCTZ.f90
    !> \brief   Deallocate allocatable variables used by the ModuleCTZ.f90
    !> \author  M. Poschmann
    !> \date    October. 06, 2021
    !> \sa      ModuleCTZ.f90
    !> \sa      ResetThermo.f90
    !> \sa      ResetAll.f90
    !
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer         Description of change
    !    ----          ----------         ---------------------
    !    06/10/2021    M. Poschmann        Original code
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this subroutine is to attempt to gracefully exit Thermochimica.  Allocatable
    !! arrays related to common tangent zone are deallocated and memory is stored for output to external packages.
    !
    ! Pertinent variables:
    ! ====================
    ! INFO                  An error is returned if deallocation is unsuccessful.
    ! INFOThermo            An integer scalar identifying whether the program exits successfully or if
    !                       it encounters an error.  A description for each error is given in ThermoDebug.f90.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ResetCTZ

    USE ModuleCTZ
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer :: i, INFO

    ! Initialize variables:
    i = 0

    if (allocated(assemblageHistory)) deallocate(assemblageHistory, STAT = INFO)
    i = i + INFO
    if (allocated(assemblageTlimits)) deallocate(assemblageTlimits, STAT = INFO)
    i = i + INFO
    if (allocated(stoichHistory)) deallocate(stoichHistory, STAT = INFO)
    i = i + INFO
    if (allocated(elementHistory)) deallocate(elementHistory, STAT = INFO)
    i = i + INFO

    lCtzInit = .FALSE.

    ! Return an INFOThermo if deallocation of any of the allocatable variables failed:
    if (i > 0) then
        INFOThermo = 15
    end if

    return

end subroutine ResetCTZ
