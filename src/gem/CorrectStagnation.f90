
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CorrectStagnation.f90
    !> \brief   Correct a stagnant system by attempting to remove a phase.
    !> \author  M.H.A. Piro
    !> \date    Mar. 15, 2013
    !> \sa      CheckPhaseAssemblage.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/15/2013      M.H.A. Piro         Original code.
    !   03/17/2013      M.H.A. Piro         Try removing pure condesed phases after trying to remove solution
    !                                        phases.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to try removing a solution phase from the system in an attempt
    !! to regain numerical stability.  This subroutine should only be called as a last resort.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! Global variables
    ! ----------------
    !
    ! nSolnPhases       An integer scalar representing the number of stable solution phases.
    ! dMolesPhase       A double real vector represending the number of moles of each stable phase in the system.
    !
    ! Local variables
    ! ---------------
    !
    ! i, j              Integer scalars; temporary index variables.
    ! iTempVec          An integer vector used by SortPick (output) representing the order of phases sorted by
    !                    the number of moles.
    ! dTempVec          A double real vector used by SortPick (input) that is to be sorted.
    ! lPhasePass        A logical variable indicating whether a particular phase assemblased has passed (TRUE)
    !                    or failed (FALSE).
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CorrectStagnation

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModulePhaseConstraints

    implicit none

    integer                       :: i, j, k
    integer, dimension(nElements) :: iTempVec
    real(8), dimension(nElements) :: dTempVec
    logical                       :: lPhasePass, lSwapLater


    ! Initialize variables:
    iTempVec   = 0
    dTempVec   = -1D6
    lPhasePass = .FALSE.
    lSwapLater = .FALSE.

    ! Check to see if a solution phase should be added to the system:
    if (iterGlobal /= iterLast) call CheckSolnPhaseAdd

    ! Proceed if the phase assemblage did not change:
    if (iterGlobal /= iterLast) then

        ! Construct temporary real vector:
        do i = 1, nSolnPhases
            j = nElements - i + 1
            dTempVec(i) = -dMolesPhase(j)
        end do

        ! Sort the indices of all solution phases in the system in ascending order of the number of moles
        ! of each phase.  The order of phases is stored in the temporary integer vector iTempVec:
        call SortPick(nElements, dTempVec, iTempVec)

        ! Loop through solution phases in an attempt to remove a phase.  Do this in ascending order:
        LOOP_SolnRem: do i = 1, nSolnPhases

            ! Store phase index:
            j = iTempVec(i)
            if (j > nSolnPhases) cycle LOOP_SolnRem

            if (nPhaseConstraints > 0) then
                k = -iAssemblage(nElements - j + 1)
                if (k > 0) then
                    if (lPhaseConstrainedSoln(k)) cycle LOOP_SolnRem
                end if
            end if

            ! Try removing this solution phase:
            call RemSolnPhase(j,lPhasePass)

            ! Exit if the new phase assemblage has passed:
            if (lPhasePass) exit LOOP_SolnRem

        end do LOOP_SolnRem

        ! If a solution phase cannot be removed from the system, then try removing a pure condensed phase:
        if (.NOT.(lPhasePass)) then

            iTempVec = 0
            dTempVec = -1D6

            ! Construct temporary real vector:
            dTempVec(1:nConPhases) = -dMolesPhase(1:nConPhases)

            ! Sort phases:
            call SortPick(nElements, dTempVec, iTempVec)

            ! Loop through stable pure condensed phases:
            LOOP_PureConRem: do i = 1, nConPhases

                j = iTempVec(i)

                if (nPhaseConstraints > 0) then
                    if (lPhaseConstrainedCon(iAssemblage(j))) cycle LOOP_PureConRem
                end if

                ! Try removing this pure condensed phase
                call RemPureConPhase(j,lSwapLater,lPhasePass)

                ! Exit if the new phase assemblage has passed:
                if (lPhasePass) exit LOOP_PureConRem

            end do LOOP_PureConRem

        end if

    end if

    return

end subroutine CorrectStagnation
