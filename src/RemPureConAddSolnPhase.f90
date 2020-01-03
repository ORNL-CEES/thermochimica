
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    RemPureConAddSolnPhase.f90
    !> \brief   Simultaneously remove a pure condensed phase and add a solution phase.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckPureConPhaseRem.f90
    !> \sa      AddSolnPhase.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/25/2012      M.H.A. Piro         Original code
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to remove a particular pure condensed phase and add a
    !! solution phase.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> param[out]   lPhasePass      A logical variable indicating whether the new estimated phase assemblage
    !!                               passed (.TRUE.) or failed (.FALSE.).
    !
    !               nConPhases      The number of pure condensed phases in the assemblage
    !               iPhaseRem       Index of pure condensed phase to be removed to the system.
    !               iPhaseAdd       Index of solution phase to be added from the system.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine RemPureConAddSolnPhase(lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::   i, iPhaseAdd, iPhaseRem
    real(8)::   dTemp
    logical::   lSwapLater, lPhasePass


    ! Initialize variables:
    iPhaseAdd  = 0
    iPhaseRem  = MAXVAL(MINLOC(dMolesPhase(1:nConPhases)))
    iPhaseRem  = MAX(iPhaseRem,1)
    lPhasePass = .FALSE.

    ! First, check to see if there are any solution phases that could be added to the system:
    LOOP_SolnPhaseSys: do i = 1, nSolnPhasesSys

        ! Skip this phase if it is already predicted to be stable:
        if (lSolnPhases(i)) cycle LOOP_SolnPhaseSys

        ! Skip this phase if it is not the first "phase" in a phase with a miscibility gap:
        if (lMiscibility(i)) cycle LOOP_SolnPhaseSys

        ! Compute the mole fractions of all constituents in this solution phase:
        call CompMolFraction(i)

        ! This solution phase is not already part of the assemblage.  Check if it should be added:
        if (dSumMolFractionSoln(i) >= 1D0) then
            iPhaseAdd = i
            exit LOOP_SolnPhaseSys
        end if
    end do LOOP_SolnPhaseSys

    ! Continue only if there is at least one solution phase that can be added to the system:
    if (iPhaseAdd > 0) then

        ! Remove the pure condensed phase manually:
        iConPhaseLast           = iAssemblage(iPhaseRem)
        iAssemblage(iPhaseRem)  = iAssemblage(nConPhases)
        iAssemblage(nConPhases) = 0
        dTemp                   = dMolesPhase(iPhaseRem)
        dMolesPhase(iPhaseRem)  = dMolesPhase(nConPhases)
        dMolesPhase(nConPhases) = 0D0
        nConPhases              = nConPhases - 1
        dMolesPhase             = dMolesPhase * 0.95D0

        ! Add the solution phase:
        call AddSolnPhase(iPhaseAdd,lSwapLater,lPhasePass)

        if (lPhasePass) then
            ! The new phase assemblage has passed.
            iterLast                = iterGlobal
            lSolnPhases(iPhaseAdd)  = .TRUE.
        else
            ! The new phase assemblage did not pass.
            nConPhases              = nConPhases + 1
            iAssemblage(nConPhases) = iConPhaseLast
            dMolesPhase(nConPhases) = dTemp
        end if

    end if

    return

end subroutine RemPureConAddSolnPhase
