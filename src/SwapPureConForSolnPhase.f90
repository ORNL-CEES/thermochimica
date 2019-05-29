
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    SwapPureConForSolnPhase.f90
    !> \brief   Swap a pure condensed phase for a solution phase.
    !> \author  M.H.A. Piro
    !> \date    May 2, 2012
    !> \sa      CheckPureConPhaseAdd.f90
    !> \sa      AddPureConPhase.f90
    !> \sa      SwapPureConPhase.f90
    !> \sa      CheckPhaseChange.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/20/2012      M.H.A. Piro         Original code
    !   02/10/2012      M.H.A. Piro         Added the capability to test a phase assemblage that is now entirely
    !                                        comprised of pure condensed phases.
    !   02/13/2012      M.H.A. Piro         Check if iPhaseChange is already part of the assemblage.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   05/02/2012      M.H.A. Piro         If the system is no longer comprised of solution phases, but it cannot
    !                                        be represented by only pure condensed phases, the number of moles of
    !                                        pure condensed phases should be returned to their previous values.
    !   08/21/2012      M.H.A. Piro         Updating the iteration history check.
    !   03/27/2013      M.H.A. Piro         If a pure condensed phase cannot be swapped for a particular solution
    !                                        phase due to a failure reported by CheckPhaseChange, then that
    !                                        particular solution phase is put at end of the listing.  This changes
    !                                        the order that solution phases are swapped.  Therefore, it was possible
    !                                        for a solution phase to be swapped, failed, then the same phase will
    !                                        be swapped in the succeeding cycle, and a different phase will be
    !                                        overlooked.  Correction: store the iAssemblage and dMolesPhase
    !                                        vectors and correct the entire vector.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to add a partiuclar pure condnesed phase and remove a solution
    !! phase.  The particular solution phase that will be removed will be determined in this subroutine.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iPhaseChange    An integer scalar representing the absolute index of a pure condensed phase
    !!                               that is to be added to the system.
    !> \param[out]  lPhasePass      A logical variable indicating whether the new phase assemblage has passed.
    !
    ! nSolnPhases           An integer scalar representing the number of solution phases in the assemblage.
    ! iAssemblage           Integer vector containing the indices of all phases currently estimated to contribute
    !                        to the equilibrium phase assemblage.
    ! iSolnPhaseLast        Index of the last solution phase to be added to, withdrawn from, or exchanged, from
    !                        the estimated phase assemblage.
    ! dMolesPhase           A double real vector representing the number of moles of each phase.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine SwapPureConForSolnPhase(iPhaseChange,lPhasePass)

    USE ModuleThermoIO, ONLY: INFOTHERmo
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                      :: i, j, k, iPhaseChange, INFO, iterBack
    integer,dimension(nElements) :: iAssemblageTest
    real(8)                      :: dTemp
    real(8),dimension(nElements) :: dTempVec
    logical                      :: lPhasePass, lSwapLater


    ! Initialize variables:
    lPhasePass      = .FALSE.
    dMolesPhaseLast = dMolesPhase

    ! Ensure that this phase is not already part of the assemblage:
    do i = 1,nConPhases
        if (iAssemblage(i) == iPhaseChange) return
    end do

    ! Loop through all solution phases in the current estimated phase assemblgae to see which one should
    ! be swapped:
    LOOP_SolnPhases: do i = 1, nSolnPhases

        ! Check iteration history
        j = nElements - i + 1
        k = nElements - nSolnPhases + 1

        iAssemblageTest(1:nElements)  = iAssemblage(1:nElements)
        iAssemblageTest(j)            = iAssemblage(k)
        iAssemblageTest(k)            = 0
        iAssemblageTest(nConPhases+1) = iPhaseChange
        lSwapLater                    = .FALSE.

        if (iterGlobal < 1000) then
            !iterBack = 200
            iterBack = 500
        else
            iterBack = 1000
        end if

        ! Check the iteration history:
        call CheckIterHistory(iAssemblageTest,iterBack,lSwapLater)

        ! Skip to this next solution phase if the phase assemblage has been previously considered:
        if (lSwapLater .EQV. .TRUE.) cycle LOOP_SolnPhases

        ! Store the info for the solution phase to be removed to temporary variables:
        iAssemblageTest         = iAssemblage
        dTempVec                = dMolesPhase
        j                       = nElements - i + 1
        dTemp                   = dMolesPhase(j)
        iSolnPhaseLast          = -iAssemblage(j)
        k                       = nElements - nSolnPhases + 1
        dMolesPhase(j)          = dMolesPhase(k)
        iAssemblage(j)          = iAssemblage(k)
        iAssemblage(k)          = 0
        dMolesPhase(k)          = 0D0
        nConPhases              = nConPhases + 1
        nSolnPhases             = nSolnPhases - 1
        iAssemblage(nConPhases) = iPhaseChange

        ! Check that this phase change is acceptable:
        call CheckPhaseChange(lPhasePass,INFO)

        if (nSolnPhases == 0) then

            ! Check the system when there are only pure condensed phases:
            call CheckSysOnlyPureConPhases

            lRevertSystem = .FALSE.

            !if (lPhasePass .EQV. .TRUE.) then
            if (lConverged .EQV. .TRUE.) then
                ! The system is comprised of only pure condensed phases.
                lConverged = .TRUE.
                exit LOOP_SolnPhases
            else
                ! The system cannot be represented by only pure condensed phases.  Return mole numbers to their
                ! previous values:
                INFOThermo  = 0
                dMolesPhase = dMolesPhaseLast
            end if

        end if

        ! Check if the new phase assemblage has passed:
        if (lPhasePass .EQV. .TRUE.) then
            ! This phase assemblage can be considered.
            iterLastSoln = iterGlobal
            iterLast     = iterGlobal
            iterSwap     = iterGlobal
            iPureConSwap = iPhaseChange
            iSolnSwap    = 0
            lSolnPhases(iSolnPhaseLast) = .FALSE.
            exit LOOP_SolnPhases
        else
            ! This phase assemblage cannot be considered.  Revert back to the previous assemblage:
            nConPhases   = nConPhases - 1
            nSolnPhases  = nSolnPhases + 1
            iAssemblage  = iAssemblageTest
            dMolesPhase  = dTempVec
            cycle LOOP_SolnPhases
        end if

    end do LOOP_SolnPhases

    return

end subroutine SwapPureConForSolnPhase