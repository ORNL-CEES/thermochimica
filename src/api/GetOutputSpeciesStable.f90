
!-------------------------------------------------------------------------------
!
!> \file    GetOutputSpeciesStable.f90
!> \brief   Get a logical array of all species identifying whether it is
!!           stable (T) or not (F).
!> \author  M.H.A. Piro
!> \date    Sept. 28, 2016
!
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   09/28/2016      M.H.A. Piro         Original code
!
!
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to compute and return a
!!   logical array representing the entire species list that indicates
!!   whether a particular species is stable (T) at equilibrium or not (F).
!!   Note that the species are ordered by successive solution phases, then
!!   followed by stoichiometric phases.
!
!
! Pertinent variables:
! ====================
!
!!                                       site fraction of said constituent.
!> \param[out]    INFO                  An integer scalar indicating a successful
!!                                       exit (== 0) or an error (/= 0).
!!
!
!-------------------------------------------------------------------------------


subroutine GetOutputSpeciesStable

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer :: j, k, iFirst, iLast


    ! Allocate memory for variable:
    if (allocated(lSpeciesStable)) deallocate(lSpeciesStable)
    allocate(lSpeciesStable(nSpecies))

    ! Initialize variable:
    lSpeciesStable = .FALSE.

    ! Loop through all stable solution phases:
    do k = 1, nSolnPhases
        ! Compute absolute solution phase index:
        j = -iAssemblage(nElements - k + 1)

        ! Determine first and last species in this phase:
        iFirst = nSpeciesPhase(j - 1) + 1
        iLast  = nSpeciesPhase(j)

        ! Assign TRUE value to all species in this phase:
        lSpeciesStable(iFirst:iLast) = .TRUE.
    end do

    ! Loop through all stable pure condensed phases:
    do k = 1, nConPhases
        ! Compute absolute phase index:
        j = iAssemblage(k)

        ! Assign TRUE value to this phase:
        lSpeciesStable(j) = .TRUE.
    end do

    return

end subroutine GetOutputSpeciesStable
