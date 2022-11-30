
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckPhaseAssemblage.f90
    !> \brief   Check whether the estimated phase assemblage needs to be modified.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      GEMSolver.f90
    !> \sa      CheckPureConPhaseRem.f90
    !> \sa      CheckPureConPhaseAdd.f90
    !> \sa      CheckSolnPhaseRem.f90
    !> \sa      CheckSolnPhaseAdd.f90
    !> \sa      CheckStagnation.f90
    !> \sa      RevertSystem.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code.
    !   07/31/2011      M.H.A. Piro         Clean up code: remove unnecessary variables, update variable names.
    !   09/14/2011      M.H.A. Piro         Modified strategy for adding solution phases to the assemblage and
    !                                        when a phase should swap for another phase.
    !   09/22/2011      M.H.A. Piro         Added the capability to revert back to the last previous phase
    !                                        assemblage when a solution phase should be removed, but can't because
    !                                        it will yield an inappropriate Jacobian.
    !   09/23/2011      M.H.A. Piro         Added the SortAssemblage subroutine to provide a good first estimate
    !                                        of a pure condensed phase to be replaced by either a new phase.
    !                                        Care of Delta Flight 238.
    !   10/26/2011      M.H.A. Piro         Clean up code: modules, fragment code into multiple subroutines.
    !   01/05/2012      M.H.A. Piro         Created the SwapSolnPhase subroutine.
    !   01/19/2012      M.H.A. Piro         Added the capability to swap a pure condensed phase for a solution
    !                                        phase.
    !   02/08/2012      M.H.A. Piro         Modified algorithm for adding a solution phase: the order that
    !                                        solution phases are considered to be added to the system are in
    !                                        descending order of the sum of their hypothetical mole fractions.
    !   02/10/2012      M.H.A. Piro         Added the capability to test a phase assemblage that is now entirely
    !                                        comprised of pure condensed phases.
    !   02/11/2012      M.H.A. Piro         Added the condition that the change to the number of moles of a pure
    !                                        condensed phase must be non-positive to remove this phase from the
    !                                        system.
    !   02/28/2012      M.H.A. Piro         If a pure condensed phase could not be removed from the system but a
    !                                        different pure condensed phase should be added, then try swapping the
    !                                        two.
    !   04/05/2012      M.H.A. Piro         Clean up code by creating four separate subroutines:
    !                                        CheckPureConPhaseRem, CheckPureConPhaseAdd, CheckSolnPhaseRem,
    !                                        CheckSolnPhaseAdd
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   06/11/2012      M.H.A. Piro         Change order of calls: check if either phase should be removed first,
    !                                        then adding a phase as opposed to pure con remove, pure con add,
    !                                        soln rem, soln add.
    !   07/04/2012      M.H.A. Piro         Only allow a phase to be added to the system when the system is
    !                                        predicted to be stagnant (happy Independence Day).
    !   09/23/2012      M.H.A. Piro         Check whether CheckPureConPhaseAdd should be called before
    !                                        CheckSolnPhaseAdd or vice versa.
    !   09/24/2012      M.H.A. Piro         Do not update the driving force of a solution phase in this subroutine
    !                                        right after it has been calculated in CheckConvergence.
    !   10/19/2012      M.H.A. Piro         If the system has not changed in 50 iterations and the functional
    !                                        norm is greater than 1000, then the system is reverted to the last
    !                                        phase assemblage.  If the last phase change was reverted, then
    !                                        you would be in a continuous loop.  Revert to the first assemblage.
    !   03/18/2013      M.H.A. Piro         Previously, the RevertSystem subroutine is called if the functional
    !                                        norm is above 1000 and 50 iterations have passed since the phase
    !                                        assemblage has changed.  An additional constraint is imposed to this
    !                                        scenario that requires the functional norm to be changing by less than
    !                                        1%.
    !   04/9/2013       M.H.A. Piro         If a previously called subroutine says that the system needs to be
    !                                        reverted, but the system has never changed, nothing will happen.
    !                                        Now, the CorrectStagnation subroutine is called for this scenario
    !                                        which tries to remove a pure condensed phase from the system.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check whether a phase should be added to, removed from, or
    !! swap another phase currently in the estimated phase assemblage and to perform any necessary action.  The
    !! conditions for adding or removing a phase from the assemblage follow:
    !!
    !! <ol>
    !! <li> Remove a pure condensed phase when the number of moles of that phase is negative, </li>
    !! <li> Remove a solution phase when the number of moles of that phase is less than a specified tolernace,
    !!       </li>
    !! <li> Add a pure condensed phase when the driving force is negative, </li>
    !! <li> Add a solution phase when the sum of all mole fractions within the phase is greater than unity.
    !!       An equivalent statement is that the driving force of this phase is negative. </li>
    !! </ol>
    !!
    !! To prevent contradicting the Gibbs Phase Rule (i.e., the maximum number of phases is equal to the number of
    !! elements when temperature and pressure are constant), a phase replaces another phase when the number of
    !! phases already in the system equals the number  of elements.  The phase assemblage is not checked every
    !! iteration to give the numerical solution a chance to converge.
    !!
    !! The main subroutines used by this solver are summarized below:
    !! <table border="1" width="800">
    !! <tr>
    !!    <td> <b> File name </td> <td> Description </b> </td>
    !! </tr>
    !! <tr>
    !!    <td> RevertSystem.f90 </td>
    !!    <td> Revert the system to a particular phase assemblage corresponding to a partiuclar iteration.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckPureConPhaseRem.f90 </td>
    !!    <td> Check if a pure condensed phase should be removed.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckSolnPhaseRem.f90 </td>
    !!    <td> Check if a solution phase should be removed.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckStagnation.f90 </td>
    !!    <td> Check if the system is stagnant. Specifically, this subroutine counts the # of solution phases
    !!          that have molar quantities that have changed by more than a specified amount.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckPureConPhaseAdd.f90 </td>
    !!    <td> Check if a pure condensed phase should be added to the system.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckSolnPhaseAdd.f90 </td>
    !!    <td> Check if a solution phase should be added to the system.  </td>
    !! </tr>
    !! </table>
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! iterGlobal            An integer scalar representing the global iteration count.
    ! iterLast              An integer scalar representing the last global iteration that resulted in a change
    !                        in the estimated phase assemblage.  This variable will be updated when the phase
    !                        assemblage changes.
    ! iterRevert            An integer scalar representing the last global iteration that resulted in reverting
    !                        the system to a previously considered phase assemblage.
    ! iterHistory           An integer matrix representing the indices of phases in the system that have been
    !                        previously considered in the iteration history.
    ! iAssemblage           An integer vector representing the indices of phases currently estimated to be stable.
    ! dMolesPhaseChange     A double real scalar representing the tolerance for the relative change of the number
    !                        of moles of a solution phase.
    ! dMaxDrivingForce      A double real scalar representing the maximum driving force of all pure condensed
    !                        phases.
    ! dDrivingForceSoln     A double real vector representing the driving force of each solution phase.
    ! dMaxChange            A double real scalar indicating the maximum change of the number of moles of a
    !                        solution phase.
    ! iMaxDrivingForce      An integer scalar representing the index of the pure condensed phase associated with
    !                        dMaxDrivingForce (default = 0).
    ! nPhasesCheck          An integer scalar representing the number of solution phases that have molar
    !                        quantities that have changed by a fraction represented by dMolesPhaseChange.
    !                        For example, if dMolesPhaseChange = 0.01, then nPhasesCheck represents the number of
    !                        solution phases that have values of dMolesPhase that have changed by >= 1%.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckPhaseAssemblage

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer      :: nPhasesCheck, j, iMaxDrivingForce
    real(8)      :: dMolesPhaseChange, dMaxChange, dMaxDrivingForce
    logical      :: lAddPhase


    ! Initialize variables:
    nPhasesCheck      = 0
    dMaxChange        = 0D0
    dMolesPhaseChange = 0.01D0
    lAddPhase         = .FALSE.

    if ((iterGlobal <  30).AND.(iterGlobal - iterLast < 2)) return
    if ((iterGlobal >= 30).AND.(iterGlobal - iterLast < 5)) return

    ! Store the indices of phases in the current estimated phase assemblage in a iteration history matrix:
    iterHistory(1:nElements,iterGlobal) = iAssemblage(1:nElements)

    ! If the system has become stagnant and is very far from convergence, revert the system:
    if ((iterGlobal - iterLast > 50).AND.(dGEMFunctionNorm > 1D3)) then
        if (dGEMFunctionNormLast == 0) then
            call RevertSystem(2)
        else if ((dGEMFunctionNorm/dGEMFunctionNormLast > 0.99D0).AND.(dGEMFunctionNormLast - dGEMFunctionNorm < 5D0)) then 
            call RevertSystem(2)
        end if
    end if

    ! If the system has become stagnant and the direction vector has been signficandly dampened, revert:
    if ((MAXVAL(DABS(dUpdateVar)) > 1D15).AND.(iterGlobal - iterLast >= 150)) lRevertSystem = .TRUE.

    ! If the system has become stagnant and the direction vector is massive, but hte system has never changed,
    ! the system cannot be reverted.  Try removing phases from the system:
    if ((MAXVAL(DABS(dUpdateVar)) > 1D15).AND.(iterGlobal >= 150).AND.(iterLast == 0)) call CorrectStagnation

    ! If a previously called subroutine requires that the system be reverted, but the system has never changed
    ! (i.e., the system cannot be reverted), then try removing a pure condensed phase from the system.
    if ((lRevertSystem).AND.(iterLast == 0)) call CorrectStagnation

    ! If the norm of the direction vector is massive and the last time the phase assemblage changed a
    ! pair of phases were swapped, try reverting the system to the last successful iteration:
    if ((MAXVAL(DABS(dUpdateVar)) > 1D15).AND.(iterLast == iterSwap).AND.(iterGlobal - iterLast >= 30)) call RevertSystem(iterLast)

    ! If the system has not changed in 500 iterations and the functional norm has not decreased by 1%, revert the system:
    if ((iterGlobal - iterLast >= 500).AND.(dGEMFunctionNorm > 0.99D0 * dGEMFunctionNormLast).AND. &
        (dGEMFunctionNorm > 0.01D0)) call RevertSystem(2)

    ! Check if the system is to be reverted to a previously considered phase assemblage:
    if (lRevertSystem) then

        if ((iterGlobal /= iterLast).AND.(nConPhases > 0))  call CheckPureConPhaseRem

        if ((iterGlobal /= iterLast).AND.(nSolnPhases > 0)) call CheckSolnPhaseRem

        if (iterGlobal /= iterLast) call RevertSystem(2)

    end if

    ! If the system has not been reverted, proceed to checking if the phase assemblage should be changed:
    IF_CheckRevert: if (iterGlobal /= iterRevert) then

        ! 1) Check if a pure condensed phase should be removed:
        if ((iterGlobal /= iterLast).AND.(nConPhases > 0))  call CheckPureConPhaseRem

        ! 2) Check if a solution phase should be removed:
        if ((iterGlobal /= iterLast).AND.(nSolnPhases > 0)) call CheckSolnPhaseRem

        ! 3) Check if a phase should be added to the system:
        IF_CheckPhaseAdd: if ((iterGlobal /= iterLast).AND.((iterGlobal - iterLast > 200).OR. &
            (dGEMFunctionNorm < dTolerance(13)))) then

            ! Check if the system is stagnant (specifically, count the number of solution phases that are
            ! changing by more than a specified value):
            call CheckStagnation(dMolesPhaseChange,dMaxChange,nPhasesCheck)

            ! Check if a phase should be added/swapped to the system only when the functional norm is
            ! below a certain tolerance AND a solution phase is not being pushed out (dMaxChange <= 0.95)
            ! AND the system is stagnant:
            IF_CheckStagnant: if ((dMaxChange <= 0.95D0).AND.((nPhasesCheck == 0).OR.(iterGlobal - iterLast >= 30))) then

                ! Compute the driving force for all pure condensed phases:
                call CompDrivingForce(iMaxDrivingForce,dMaxDrivingForce)

                if (lDebugMode) print *, 'condensed phase driving forces:', iMaxDrivingForce, dMaxDrivingForce

                ! Loop through all metastable solution phases to compute the mole fractions of constituents
                ! and the driving force of each solution phase:
                do j = 1, nSolnPhasesSys

                    ! The determination of the driving force in the CheckConvergence subroutine is more rigorous
                    ! than in CompMolFraction.  If the current driving force is negative and the functional norm
                    ! is very small, then cycle to the next solution phase and use the current driving force.
                    if ((iterLastMiscGapCheck == iterGlobal - 1).AND.(dDrivingForceSoln(j) < 0D0)) cycle

                    ! Only compute the mole fractions of constituents belonging to unstable phases:
                    if (.NOT.(lSolnPhases(j))) then

                        if (lMiscibility(j)) then

                            ! This phase may contain a miscibility gap.  Periodically check for a miscibility gap:
                            if (iterGlobal - iterLastMiscGapCheck >= 50) then

                                ! Check for a miscibility gap:
                                call CheckMiscibilityGap(j,lAddPhase)
                            else
                                call CompMolFraction(j)
                            end if

                        else
                            ! This is phase does not contain a miscibility gap.  Compute the mole fractions:
                            call CompMolFraction(j)

                        end if
                    end if

                end do

                if (lDebugMode) print *, 'solution phase driving forces:', MINLOC(dDrivingForceSoln), MINVAL(dDrivingForceSoln)


                ! Determine whether a pure condensed phase should be considered first or a solution phase.
                ! This decision is based on whether the lowest driving force of a pure condensed phase is less
                ! than that of all solution phases:
                IF_AddOrder: if (dMaxDrivingForce < MINVAL(dDrivingForceSoln)) then

                    ! 4A) Check if a pure condensed phase should be added/swapped:
                    if (iterGlobal /= iterLast) call CheckPureConPhaseAdd(iMaxDrivingForce, dMaxDrivingForce)

                    ! 5A) Check if a solution phase should be added/swapped:
                    if (iterGlobal /= iterLast) call CheckSolnPhaseAdd

                else

                    ! 4B) Check if a solution phase should be added/swapped:
                    if (iterGlobal /= iterLast) call CheckSolnPhaseAdd

                    ! 5B) Check if a pure condensed phase should be added/swapped:
                    if (iterGlobal /= iterLast) call CheckPureConPhaseAdd(iMaxDrivingForce, dMaxDrivingForce)

                end if IF_AddOrder

            end if IF_CheckStagnant

        end if IF_CheckPhaseAdd

    end if IF_CheckRevert

    ! If the phase assemblage has changed, relax the functional norm for the following iteration:
    if (iterGlobal == iterLast) dGEMFunctionNorm = 10D0 * dGEMFunctionNorm

    return

end subroutine CheckPhaseAssemblage
