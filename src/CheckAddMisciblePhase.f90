
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckAddMisciblePhase.f90
    !> \brief   Check the
    !> \author  M.H.A. Piro
    !> \date    Oct. 21, 2012
    !> \sa      AddSolnPhase.f90
    !> \sa      SwapSolnForPureConPhase.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   10/21/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check if the phase that is being
    !! added is miscible and if one of the other corresponding phases is
    !! currently not predicted to be stable.  If so, this subroutine checks if
    !! the absolute index of the phase that is to be added is greater than the other
    !! phase, in which case it swaps the two phases.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in] iPhaseAddAbs  Absolute index of the miscible solution phase that is to be added to the system.
    !
    ! iPhaseAddAbs    - Absolute index of solution phase that is being added.
    ! iPhaseOtherAbs  - Absolute index of solution phase that has a miscibility
    !                    gap with iPhaseAddAbs (default set to zero).    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckAddMisciblePhaseIndex(iPhaseAddAbs)

    USE ModuleThermo
    USE ModuleGEMSolver, ONLY: lSolnPhases

    implicit none

    integer:: i, j, k, l, iPhaseAddAbs, iPhaseOtherAbs


    ! Initialize variables:
    iPhaseOtherAbs = 0

    ! Only proceed if there is at least one other solution phase:
    if (nSolnPhases > 0) then

        ! First, see if any of the other solution phases currently assumed to be
        ! stable corresponds to this miscible phase:
        LOOP_A: do j = 1, nSolnPhasesSys

            ! Cycle if this is the same phase:
            if (j == iPhaseAddAbs) cycle LOOP_A

            ! Check if this phase is NOT stable and corresponds to iPhaseAddAbs:
            if ((cSolnPhaseName(j) == cSolnPhaseName(iPhaseAddAbs)).AND.(lSolnPhases(j) .EQV. .FALSE.)) then
                iPhaseOtherAbs = j
                exit LOOP_A
            end if

        end do LOOP_A

        ! If the absolute index of the solution phase to be added is greater than the corresponding phase
        ! that is predicted to be unstable, then swap the two phases:
        if ((iPhaseOtherAbs < iPhaseAddAbs).AND.(iPhaseOtherAbs /= 0)) then

            ! Store indices of solution species for both phases:
            i = nSpeciesPhase(iPhaseAddAbs-1) + 1
            j = nSpeciesPhase(iPhaseAddAbs)
            k = nSpeciesPhase(iPhaseOtherAbs-1) + 1
            l = nSpeciesPhase(iPhaseOtherAbs)

            ! Swap mole fractions, moles and chemical potential terms:
            dMolFraction(k:l)       = dMolFraction(i:j)
            dMolesSpecies(k:l)      = dMolesSpecies(i:j)
            dChemicalPotential(k:l) = dChemicalPotential(i:j)

            ! Reset the moles of solution species for the other phase:
            dMolesSpecies(i:j) = 0D0

            ! Swap the indices:
            iPhaseAddAbs = iPhaseOtherAbs

        end if

    end if

    return

end subroutine CheckAddMisciblePhaseIndex
