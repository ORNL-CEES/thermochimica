
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckPureConPhaseRemHistory.f90
    !> \brief   Check the iteration history when there are multiple pure condensed phases that should be
    !!           removed from the system.
    !> \author  M.H.A. Piro
    !> \date    June 11, 2013
    !> \sa      CheckPureConPhaseRem.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   06/11/2013      M.H.A. Piro         Original code.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check the iteration history when multiple pure condensed
    !! phases are negative.  In this situation, there are multiple candidate phase assemblages that may be
    !! considered; however, Thermochimica only allows one phase to be withdrawn from the system at a time.
    !! To avoid cycling issues, the iteration history is checked instead of just taking the phase
    !! with the most negative mass.  In other words, the system might keep trying to remove the same phase
    !! over and over again.
    !!
    !! First, a check is performed to confirm whether there are multiple pure condensed phases with negative
    !! molar quantities.  Second, the iteration history is checked for each candidate phase assmeblage by
    !! withdrawing one of the negative phases.  Third, if the phase assemblage has been previously considered,
    !! this phase is placed at the end of the iTempVec vector.  Thus,
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[inout] iTempVec  A temporary integer vector representing the order of pure condensed phases to be
    !!                        considered for removal.
    ! nConPhases            The number of pure condensed phases in the assemblage.
    ! nSolnPhases           The number of solution phases in the assemblage.
    ! nElements             The number of elements in the system.
    ! nNegPureConPhases     An integer scalar representing the number of negative pure condensed phases.
    ! iterBack              An integer scalar representing the number of iterations that should be checked in
    !                        the history.
    ! iAssemblageTest       A temporary integer vector used to represent the indices of phases in the phase
    !                        assemblage.
    ! lSwapLater            A logical scalar used within the CheckIterHistory subroutine indicating whether
    !                        a phase should be swapped later (TRUE) or not (FALSE).
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckPureConPhaseRemHistory(iTempVec)

    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                       :: i, j, k, nNegPureConPhases, iterBack
    integer,dimension(nConPhases) :: iTempVec
    integer,dimension(nElements)  :: iAssemblageTest
    logical                       :: lSwapLater


    ! Initialize variables:
    iterBack          = 500
    nNegPureConPhases = 0
    lSwapLater        = .FALSE.

    ! Count the number of pure condensed phases that are negative:
    do i = 1, nConPhases
        if (dMolesPhase(i) <= 0D0) nNegPureConPhases = nNegPureConPhases + 1
    end do

    ! Proceed if there is more than one pure condensed phase that has a negative mass:
    IF_MultiNegPureCon: if (nNegPureConPhases > 1) then

        ! Loop through the pure condensed phases that are negative:
        LOOP_NegPureCon: do i = 1, nNegPureConPhases

            ! Reset the iAssemblageTest vector:
            iAssemblageTest = iAssemblage

            ! Store the phase index (previously sorted);
            j = iTempVec(i)

            ! Error check: the value in the iTempVec vector must be positive:
            if (j == 0) then
                ! Report an error and exit:
                INFOThermo = 35
                exit LOOP_NegPureCon
            end if

            ! Remove this pure condensed phase from the temporary integer vector:
            iAssemblageTest(j)          = iAssemblage(nConPhases)

            ! Reset the previous entry:
            iAssemblageTest(nConPhases) = 0

            ! Reduce the number of pure condensed phases by one:
            nConPhases = nConPhases - 1

            ! Check whether this particular phase assemblage has been previously considered:
            call CheckIterHistory(iAssemblageTest,iterBack,lSwapLater)

            ! Reset the number of pure condensed phases after performing iteration history check:
            nConPhases = nConPhases + 1

            ! If the phase assmeblage has already been considered, then shuffle the iTempVec vector:
            if (lSwapLater) then

                ! Shuffle the vector:
                k                    = iTempVec(nConPhases)
                iTempVec(nConPhases) = j
                iTempVec(i)          = k

                ! Reinitialize the lSwapLater variable:
                lSwapLater           = .FALSE.

                ! Exit the loop (I know that the entire vector should be shuffled, but this is a very
                ! rare occurance):
                exit LOOP_NegPureCon

            end if

        end do LOOP_NegPureCon

    end if IF_MultiNegPureCon

    return

end subroutine CheckPureConPhaseRemHistory
