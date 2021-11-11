
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckIterHistory.f90
    !> \brief   Check the iteration history to see if a particular phase assemblage has previously been
    !!           considered.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      AddPureConPhase.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   02/13/2012      M.H.A. Piro         Original code
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   05/03/2012      M.H.A. Piro         Previously, a check was made of the SUM of differences of integer
    !                                        vectors iAssemblage and iterHistory.  It is possible that the
    !                                        vectors are different, but the SUM is still zero.  A check was
    !                                        added in the event that this sum is zero and then a more precise
    !                                        comparison is made.
    !   07/04/2012      M.H.A. Piro         If the functional norm is below a specified tolernace, do not
    !                                        check the iteration history.
    !   07/12/2012      M.H.A. Piro         Allow for an entire phase assemblage to be provided as input.
    !                                        This makes this subroutine more general and it can be used for
    !                                        changing the phase assemblage in a number of ways.
    !   07/25/2012      M.H.A. Piro         Allow for an addition check to see order that phases are added
    !                                        to, removed from, or exchanged in the system.
    !   08/17/2012      M.H.A. Piro         Fix bug: the check that looks at the previous phase assemblage
    !                                        exits if the phase assemblage is the same, but it didn't cycle
    !                                        LOOP_History if this was false.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check whether a particular phase assemblage has been
    !! considered in the iteration history.  An additional test is performed that checks the order that a phase
    !! is added to, removed from, or swapped from the system.  For example, a system may be comprised of
    !! phases A, B, C and D at iteration w, phase D is removed at iteration x, and phase B is prematurely
    !! removed at iteration y.  At iteration z, one determines that phase B should be added to the system,
    !! which would result in an identical phase assemblage at iteration x.  This test recognizes that the
    !! phase assemblages considered at iteration x and z are identical, but the changes to the system that
    !! results in these assemblages are different.
    !!
    !! The logical variable lSwapLater is returned indicating whether this phase assemblage can be considered
    !! or not.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iAssemblageTest An integer vector representing the phase assemblage to be tested.
    !> \param[in]   iterBack        An integer scalar indicating how many iteratoins to go back in history.
    !> \param[out]  lSwapLater      A logical variable indicating whether a phase should be swapped later on.
    !
    ! iAssemblage           Integer vector containing the indices of phases estimated to be part of the
    !                       equilibrium phase assemblage.
    ! dGEMFunctionNorm      A double real scalar representing the functional norm in the Gibbs Energy Minimization
    !                        (GEM) solver.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckIterHistory(iAssemblageTest,iterBack,lSwapLater)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::                       i, j, k, l, iterBack, iterFirst
    integer,dimension(nElements)::  iAssemblageTest, iAssemblageLast
    logical::                       lSwapLater

    ! Initialize variables:
    lSwapLater      = .FALSE.
    iAssemblageLast = 0

    ! Record the last successful phase assemblage from iterGlobal:
    LOOP_Last: do i = iterGlobal, 1, -1
        if ((SUM(iterHistory(1:nElements,i) - iAssemblageTest(1:nElements)) /= 0) .AND. &
            (SUM(iterHistory(1:nElements,i)) /= 0)) then
            ! Store the last successful phase assemblage:
            iAssemblageLast(1:nElements) = iterHistory(1:nElements,i)
            exit LOOP_Last
        end if
    end do LOOP_Last

    ! Check iteration history:
    iterFirst = MAX (50, iterGlobal - iterBack)

    LOOP_History: do j = (iterGlobal - 1), iterFirst, -1
        ! Quick check:
        IF_Quick_Check: if (SUM(iAssemblageTest(1:nElements) - iterHistory(1:nElements,j)) == 0) then
            ! Double check that this phase assemblage is consistent:
            LOOP_Outer: do k = 1, nElements        ! Loop through current phase assemblage.
                LOOP_Inner: do l = 1, nElements     ! Loop through previous phase assemblage.
                    ! The phase indices are consistent; move on to the next phase in the current phase assemblage:
                    if (iterHistory(l,j) == iAssemblageTest(k)) cycle LOOP_Outer
                end do LOOP_Inner
                ! There is a discrepancy between the two phase assemblages.
                ! Move on to the next assemblage in the iteraiton history:
                cycle LOOP_History
            end do LOOP_Outer

            ! CONSIDER THIS TEMPORARY TO CONSIDER:
            if (iterGlobal <= 1000) then
                lSwapLater = .TRUE.
                exit LOOP_History
            end if

            ! The current phase assemblage is the same as the assemblage corresponding to iteration j.
            ! Now, compare the phase assemblages before j and iterGlobal.
            IF_Last_Assemblage: if (SUM(iAssemblageLast) /= 0) then
                k = MAX(1,iterFirst - 10)
                LOOP_A: do i = j, k, -1
                    ! Cycle if the row is full of zeros:
                    if ((iterHistory(1,i) == 0).AND.(iterHistory(nElements,i) == 0)) cycle LOOP_A
                    ! Check for a change in the estimated phase assemblage:
                    if (SUM(iterHistory(1:nElements,i) - iterHistory(1:nElements,j)) /= 0) then
                        ! The last successful phase assemblage before iteration j is at iteration i.
                        ! Test if this is the same as the phase assemblage iAssemblageLast:
                        if (SUM(iAssemblageLast(1:nElements) - iterHistory(1:nElements,i)) == 0) then
                            ! This phase assemblage has been previously considered.  Swap later.
                            lSwapLater = .TRUE.
                            exit LOOP_History
                        else
                            ! This phase assemblage has not been previously considered.  Do not swap later.
                            cycle LOOP_History
                        end if
                    end if
                end do LOOP_A
            end if IF_Last_Assemblage
        end if IF_Quick_Check
    end do LOOP_History

    return

end subroutine CheckIterHistory
