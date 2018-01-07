
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckSolnPhaseAdd.f90
    !> \brief   Check if a solution phase should be added to the system.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckPhaseAssemblage.f90
    !> \sa      SortPick.f90
    !> \sa      AddSolnPhase.f90
    !> \sa      ShuffleAssemblage.f90
    !> \sa      SwapSolnForPureConPhase.f90
    !> \sa      SwapSolnPhase.f90
    !> \sa      CompMolFraction.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   05/04/2012      M.H.A. Piro         Original code.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   08/22/2012      M.H.A. Piro         Apply check for a miscibility gap.
    !   08/31/2012      M.H.A. Piro         Call CompMolFractionQKTO after CheckMiscibilityGap is called to
    !                                        reduce the residuals of the chemical potential terms.  This
    !                                        procedure promotes convergence.
    !   09/23/2012      M.H.A. Piro         A solution phase is no longer added if x_{\lambda} > 1, rather when
    !                                        the driving force is negative.  Although these statements are equivalent,
    !                                        it is computationally advantageous to work with the driving force instead
    !                                        of the sum of mole fractions.
    !   09/24/2012      M.H.A. Piro         ShuffleAssemblage now shuffles the entire assemblage and indicates
    !                                        whether a pure condensed or solution phase should be swapped first.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check whether a solution phase should be added to the current
    !! estimated assemblage of stable phases in the system.  The condition for a solution phase to be added
    !! to the system is when the sum of hypothetical mole fractions of all constituents within a solution phase
    !! exceeds unity (within tolerance).  A solution phase can be added to the system when the total number of
    !! phases is less than the number of elements in the system; however, it must replace an existing phase when
    !! the number of phases equals the number of elements (to prevent contradicting Gibbs' Phase Rule).
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nConPhases            The number of pure condensed phases in the assemblage.
    ! nSolnPhases           The number of solution phases in the assemblage.
    ! nSolnPhasesSys        The number of solution phases in the database.
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
    ! lSwapCheck            A logical variable indicating whether the system will allow a solution phase swap
    !                        to be checked.
    ! dSumMolFractionSoln   A double real vector representing the sum of mole fractions in each solution phase.
    ! dGEMFunctionNorm      A double real scalar representing the norm of the functional vector in the PGESolver.
    ! lSolnPhase            A logical vector indicating whether a particular solution phase is stable (TRUE) or
    !                        unstable (FALSE).
    ! lMiscibility          A logical vector indicating whether a particular solution phase contains a miscibility
    !                        gap (TRUE) or not (FALSE).
    ! dDrivingForceSoln     A double real vector representing the driving force of each solution phase in the
    !                        database.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckSolnPhaseAdd

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                           :: i, j
    integer,dimension(nSolnPhasesSys) :: iTempVec
    logical                           :: lSwapLater, lPhasePass, lSwapCheck


    ! Initialize variables:
    lSwapLater = .FALSE.
    lSwapCheck = .FALSE.
    lPhasePass = .FALSE.
    iTempVec   = 0

    if (iterGlobal - iterlast < 5) return

    ! Determine whether a solution phase will be checked for a swap:
    !if ((dGEMFunctionNorm <= 1D-3).OR.(iterGlobal - iterlast >= 50)) lSwapCheck = .TRUE.
    if ((dGEMFunctionNorm <= dTolerance(11)).OR.(iterGlobal - iterlast >= 50)) lSwapCheck = .TRUE.

    ! Sort the indices of all solution phases in the system in descending order of the driving
    ! force of each phase.  The order of phases is stored in the temporary integer vector iTempVec:
    call SortPick(nSolnPhasesSys, dDrivingForceSoln, iTempVec)

    ! Loop through all solutions phases in the database:
    LOOP_SolnPhaseAdd: do j = 1, nSolnPhasesSys

        ! Absolute solution phase index (defined by iTempVec, which is sorted by dDrivingForceSoln):
        ! (The order of iTempVec is reversed  because SortPick sorts in ascending order, not descending)
        i = iTempVec(nSolnPhasesSys - j + 1)

        ! Skip this phase if it is already part of the assemblage:
        if (lSolnPhases(i) .EQV. .TRUE.) cycle LOOP_SolnPhaseAdd

        ! Only add a miscible phase if the functional norm is below a certain value:
       ! if ((lMiscibility(i) .EQV. .TRUE.).AND.(dGEMFunctionNorm > 1D-4)) cycle LOOP_SolnPhaseAdd
        if ((lMiscibility(i) .EQV. .TRUE.).AND.(dGEMFunctionNorm > 1D-4).AND.(iterGlobal - iterLast < 300)) cycle LOOP_SolnPhaseAdd

        ! Check if a solution phase should be added:
        !IF_SolnPhaseAdd: if ((dDrivingForceSoln(i) < dTolerance(4)).AND. &
         !   (nConPhases + nSolnPhases < nElements - nChargedConstraints)) then
        IF_SolnPhaseAdd: if ((dDrivingForceSoln(i) < dTolerance(4)).AND. &
            (nConPhases + nSolnPhases < nElements)) then

            ! Try adding this solution phase to the assemblage:
            call AddSolnPhase(i,lSwapLater,lPhasePass)

            ! If the phase assemblage changed, then return control to the main solver:
            if (lPhasePass .EQV. .TRUE.) exit LOOP_SolnPhaseAdd

            ! The solution phase could not be added the system directly, but it might be possible to
            ! swap another phase that is currently in the phase assmeblage for this one:
            if ((lSwapLater .EQV. .TRUE.).AND.(lSwapCheck .EQV. .TRUE.)) then

                ! Check if this solution phase can be swapped with another phase:
                call CheckSolnPhaseSwap(i,lPhasePass)

                ! Exit if this phase assemblage is appropriate:
                if (lPhasePass .EQV. .TRUE.) exit LOOP_SolnPhaseAdd

            end if

        elseif ((dDrivingForceSoln(i) < dTolerance(4)).AND.(nConPhases + nSolnPhases == nElements).AND. &
        !elseif ((dDrivingForceSoln(i) < dTolerance(4)).AND.(nConPhases + nSolnPhases == nElements - nChargedConstraints).AND. &
            (lSwapCheck .EQV. .TRUE.).AND.(iterGlobal - iterlast > 10)) then

            ! Check if this solution phase can be swapped with another phase:
            call CheckSolnPhaseSwap(i,lPhasePass)

            ! Exit if this phase assemblage is appropriate:
            if (lPhasePass .EQV. .TRUE.) exit LOOP_SolnPhaseAdd

        end if  IF_SolnPhaseAdd

    end do LOOP_SolnPhaseAdd

    return

end subroutine CheckSolnPhaseAdd


!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to check whether a solution phase
    ! should swap for a pure condensed phase or another solution phase.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! i             Absolute index of phase to be added to the system.
    ! lPhasePass    A logical scalar indicating whether the new phase assemblage
    !                has passed (TRUE) or not (FALSE).
    ! iPhaseTypeOut An integer scalar representing the phase type to be
    !                withdrawn from the system to accomodate a new phase.  This
    !                variable is equal to zero if a pure condensed phase is to
    !                be removed or it is equal to one if a solution phase is to
    !                be removed.
    ! iPhaseChange  An integer scalar representing the phase to be added to the
    !                system (used to be consistent with variable naming
    !                convention in the corresponding subroutines.
    !
    !---------------------------------------------------------------------------


subroutine CheckSolnPhaseSwap(i,lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::   i, iPhaseTypeOut, iPhaseChange
    logical::   lPhasePass


    ! Initialize variables:
    iPhaseChange  = -i

    ! Shuffle the phase assemblage so that the best candidate to be swapped is first:
    call ShuffleAssemblage(iPhaseChange,iPhaseTypeOut)

    ! Check if a pure condensed phase or a solution phase should be swapped first:
    if (iPhaseTypeOut == 0) then
        ! A pure condensed phase should be swapped first.

        ! Swap a solution phase for a pure condensed phase:
        call SwapSolnForPureConPhase(iPhaseChange,lPhasePass)

        ! Swap a solution phase for another solution phase:
        if (iterGlobal /= iterLast)  call SwapSolnPhase(i,lPhasePass)

    elseif (iPhaseTypeOut == 1) then
        ! A solution phase should be swapped first.

        ! Swap a solution phase for another solution phase:
        call SwapSolnPhase(i,lPhasePass)

        ! Swap a solution phase for a pure condensed phase:
        if (iterGlobal /= iterLast) call SwapSolnForPureConPhase(iPhaseChange,lPhasePass)

    end if

end subroutine  CheckSolnPhaseSwap


    !---------------------------------------------------------------------------
    !                       END - CheckSolnPhaseAdd.f90
    !---------------------------------------------------------------------------
