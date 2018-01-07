
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckPureConPhaseRem.f90
    !> \brief   Check whether a pure condensed phase should be removed.
    !> \author  M.H.A. Piro
    !> \date    May 23, 2012
    !> \sa      CheckPhaseAssemblage.f90
    !> \sa      RemPureConPhase.f90
    !> \sa      RemPureConAddSolnPhase.f90
    !> \sa      SwapPureConPhase.f90
    !> \sa      ShuffleAssemblage.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/05/2012      M.H.A. Piro         Original code.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   05/23/2012      M.H.A. Piro         Allow the ability to loop through all phases in case if the most
    !                                        negative cannot be removed, another can be removed.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check whether a pure condensed phase should be removed from
    !! the system.  The condition for a pure condensed phase to be removed from the system is when the number of
    !! moles of that phase is below a specified tolerance and the change to the number of moles of the phase
    !! is below an arbitrarily small value.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nConPhases            The number of pure condensed phases in the assemblage.
    ! nSolnPhases           The number of solution phases in the assemblage.
    ! nElements             The number of elements in the system.
    ! iPhaseChange          An integer scalar representing the index of a phase that will be added to, or
    !                        withdrawn from, the current estimated phase assemblage.
    ! iPhaseTypeOut         An integer scalar representing the phase type to be shuffled by the ShuffleAssemblage
    !                        subroutine.  This subroutine shuffles the assemblage to determine the best candidate
    !                        phase to be withdrawn from the assemblage.
    !                        iPhaseTypeOut = 0 for a pure condensed phase,
    !                        iPhaseTypeOut = 1 for a solution phase.
    ! lPhasePass            A logical variable indicating whether the new estimated phase assemblage passed
    !                        (.TRUE.) or failed (.FALSE.).
    ! lSwapLater            A logical variable indicating whether a phase should be swapped when lPhasPass
    !                        indicates that the phase cannot be added directly.
    ! dMolesPhase           A double real vector representing the number of moles of each phase.
    ! dPGEFunctionNorm      A double real scalar representing the norm of the functional vector in the PGESolver.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckPureConPhaseRem

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                       :: i, iPhaseChange, iPhaseTypeOut, iMaxDrivingForce
    integer,dimension(nConPhases) :: iTempVec
    real(8)                       :: dMaxDrivingForce, dTemp
    real(8),dimension(nConPhases) :: dTempVec
    logical                       :: lSwapLater, lPhasePass


    ! Initialize variables:
    dTempVec(1:nConPhases) = -dMolesPhase(1:nConPhases)
    iTempVec               = 0

    ! Return if the smallest value of dMolesPhase is already above the tolerance:
    if (MINVAL(dMolesPhase(1:nConPhases)) > dTolerance(7)) return

    ! Sort the indices of all pure condensed phases in the ascending order of the number of moles:
    call SortPick(nConPhases,dTempVec,iTempVec)

    ! Check the iteration history when there are multiple pure condensed phases that are negative:
    call CheckPureConPhaseRemHistory(iTempVec)

    ! Loop through all pure condensed phases currently predicted to be stable:
    LOOP_PureConPhases: do i = 1, nConPhases

        ! Initialize variables:
        lSwapLater   = .FALSE.
        lPhasePass   = .FALSE.

        ! Determine which phase should be removed:
        iPhaseChange = iTempVec(i)
        dTemp        = dMolesPhase(iPhaseChange) - dMolesPhaseLast(iPhaseChange)

        if (dMolesPhase(iPhaseChange) < -1D3) dTemp = -1D0
        if (iterGlobal - iterLast >= 20)      dTemp = -1D0

        ! Check if the # of moles of this phase is less than tolerance for two consecutive iterations and that
        ! it is decreasing by more than 1%:
        IF_RemPureConPhase: if ((dMolesPhase(iPhaseChange) < dTolerance(7)).AND.&
            (dMolesPhaseLast(iPhaseChange) < dTolerance(7)).AND.(dTemp <= 0.01D0*DABS(dMolesPhase(iPhaseChange)))) then

            ! If the pure condensed phase that is to be removed was the last phase that was swapped, then
            ! revert the system back to the last successful phase assemblage:
            if ((iPureConSwap == iAssemblage(iPhaseChange)).AND.(iterLast == iterSwap)) then

                ! Revert the system to the last succesful phase assemblage:
                call RevertSystem(iterSwap)

                ! Reset iPureConSwap:
                iPureConSwap = 0

                ! Exit
                exit LOOP_PureConPhases

            end if

            ! Try removing this pure condensed phase from the phase assemblage:
            call RemPureConPhase(iPhaseChange,lSwapLater,lPhasePass)

            ! Exit if the phase assemblage has passed:
            if (lPhasePass .EQV. .TRUE.) exit LOOP_PureConPhases

            ! If the pure condensed phase could not be removed, but a pure condensed phase should be added,
            ! try swapping the two.
            IF_SwapPureConPhase: if ((lSwapLater .EQV. .TRUE.).AND.(lPhasePass .EQV. .FALSE.)) then

                ! Compute the driving force for all pure condensed phases:
                call CompDrivingForce(iMaxDrivingForce,dMaxDrivingForce)

                ! Check if a pure condensed phase can be added to the system:
                if (dMaxDrivingForce < dTolerance(4)) then

                    ! TO DO: I think that it would be better to create a new subroutine that specifically removes
                    ! iPhaseChange and specifically removes iMaxDrivingForce.  NOTE: This is not the same as
                    ! SwapPureConPhase.  If this were done, then a call to ShuffleAssemblage is unnecessary.

                    iPhaseTypeOut = 0

                    ! Re-organize the phase assemblage so that the best candidate is first:
                    call ShuffleAssemblage(iMaxDrivingForce,iPhaseTypeOut)

                    ! Try swapping this phase for another pure condensed phase in the assemblage:
                    call SwapPureConPhase(iMaxDrivingForce,lSwapLater,lPhasePass)

                    ! Exit if the phase assemblage has passed:
                    if (lPhasePass .EQV. .TRUE.) exit LOOP_PureConPhases

                end if

                if (nSolnPhases > 0) then
                    ! It may be possible that a solution phase needs to be removed first.
                    if (MINVAL(dMolesPhase(nElements - nSolnPhases + 1: nElements)) < dTolerance(7)) &
                        exit LOOP_PureConPhases
                end if

                ! If a pure condensed phase cannot be swapped for another pure condense phase, try adding
                ! swapping it for a solution phase:
                if (lPhasePass .EQV. .FALSE.) call RemPureConAddSolnPhase(lPhasePass)

                ! Exit if the phase assemblage has passed:
                if (lPhasePass .EQV. .TRUE.) exit LOOP_PureConPhases

            end if IF_SwapPureConPhase

        end if IF_RemPureConPhase

    end do  LOOP_PureConPhases

    return

end subroutine CheckPureConPhaseRem
