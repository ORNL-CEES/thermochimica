
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ShuffleAssemblage.f90
    !> \brief   Shuffle the phase assemblage in the order that is most favorable for phase exchange.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      SortPick.f90
    !> \sa      CheckPureConPhaseAdd.f90
    !> \sa      CheckPureConPhaseRem.f90
    !> \sa      CheckSolnPhaseAdd.f90
    !> \sa      GetNewAssemblage.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   09/22/2011      M.H.A. Piro         Original code
    !   10/21/2011      M.H.A. Piro         Clean up code: modules
    !   01/19/2012      M.H.A. Piro         Added the capability of sorting solution phases.
    !   02/14/2012      M.H.A. Piro         Sort the assemblage in terms of the Euclidean norm.
    !   03/06/2012      M.H.A. Piro         Added the capability to shuffle the assemblage when the new
    !                                       phase and the phase to be removed are both solution phases.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   09/24/2012      M.H.A. Piro         Sort the entire phase assemblage, not just one type of phase.
    !                                        The advantage of this is that another subroutine can determine
    !                                        whether it is best to swap a pure condensed or solution phase first.
    !   09/29/2012      M.H.A. Piro         If the system contains a pair of miscible phases, then set the
    !                                        Euclidean Norm to an arbitrarily large value to place these
    !                                        phases to the back of the list.
    !   11/4/2012       M.H.A. Piro         Fixed typo in Euclidean norm vector.  Coefficients that do not
    !                                        correspond to a phase should have a value of 1000 instead of 0,
    !                                        otherwise the sorting routine will place it at the top.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to shuffle the current estimated phase assemblage such that
    !! the order of phases in iAssemblage are in oder of atomic similarity to a new phase that is to be added
    !! to the system.  The Gibbs Phase Rule dictates that the maximum number of phases that can coexist cannot
    !! exceed the number of system components (taken here as chemical elements).  Therefore, if the number of
    !! phases currently expected to be stable is equal to the number of elements, another phase must be
    !! withdrawn from the system to accomodate this new phase.
    !!
    !! The principle of this technique is to quantify the atomic similarity between the new phase and all other
    !! phases in the system.  The phase with the most similar atomic constituency is the most likely best
    !! candidate. The principle of this technique is to compute the Euclidean norm vector in nElements dimensional
    !! space, where one point is the stoichiometry of the new phase to be introduced and each other point
    !! represents the other phases in the current phase assemblage.  The iAssemblage vector is reorganized
    !! in descending order of the Euclidean Norm.
    !!
    !! The principle of this technique are discussed in greater detail in the following literature:
    !! - M.H.A. Piro and S. Simunovic, "Performance Enhancing Algorithms for Computing Thermodynamic
    !!   Equilibria," CALPHAD, 39 (2012) 104-110.
    !!
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iNewPhase       Index of the new phase to be introduced.  If it is positive, it is a pure
    !!                               condensed phase.  If it is negative, it is a solution phase.
    !> \param[in]   iPhaseTypeOut   An integer idenitifying the type of phase with the smallest Euclidean Norm.
    !!                               This can be used to determine whether a pure condensed phase should be
    !!                               swapped first or a solution phase.
    !!
    !!                               iPhaseTypeOut = 0: pure condensed phase;
    !!                               iPhaseTypeOut = 1: solution phase.
    !
    ! iAssemblage           An integer vector containing the indices of phases in the assemblage
    !                        (1:nConphases represent pure condensed phases and
    !                        (nElements-nSolnPhases:nElements) represent solution phases.
    ! nElements             The number of elements in the system.
    ! nConPhases            The number of pure condensed phases in the assemblage
    ! nSolnPhases           The number of solution phases in the assemblage
    ! nSolnPhasesSys        The number of solution phases in the system
    ! dMolesPhase           The number of moles of a phase.  These are directly mapped to phases in iAssemblage
    ! dEuclideanNorm        A double vector representing the Euclidean norm between the stoichiometry of the new
    !                       phase each existing pure condensed phase in the assemblage.
    ! dAtomFractionSpecies  Atomic fraction of a particular element (column) for a particular species (row).
    ! dEffStoichSolnPhase   The effective stoichiometry of a particular element (column) in a particular solution
    !                       phase (row).
    ! iVec                  An integer vector used for internal operations for the sorting routine.
    ! iTempVec              An integer vector storing the phase assemblage at the beginning of the calculation
    !                        to allow the system to be reverted in case if ShuffleAssemblage fails.
    ! dTempVec              A double real temporary vector storing the number of moles of all phases at the
    !                        beginning of the calculation to allow the system to be reverted in case if
    !                        ShuffleAssemblage fails.
    ! iMisciblePhaseID      An integer scalar representing the relative phase index of a phase containing
    !                        a miscibility gap.  If there aren't any phases with a miscibility gap, it is set
    !                        to zero.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ShuffleAssemblage(iNewPhase,iPhaseTypeOut)

    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                      :: i, j, k, l, iNewPhase, iPhaseTypeOut, iAssemblageSum, iMisciblePhaseID
    integer,dimension(nElements) :: iVec, iVecB, iTempVec
    real(8)                      :: dTemp
    real(8),dimension(nElements) :: dEuclideanNorm, dTempVec


    ! Initialize variables:
    iVec             = 0
    iMisciblePhaseID = 0
    dEuclideanNorm   = 1D3

    ! Store the current phase assemblage and mole vector in case if a failure is encountered:
    iTempVec       = iAssemblage
    iAssemblageSum = SUM(iAssemblage)
    dTempVec       = dMolesPhase

    ! Compute the effective stoichiometry of each solution phase predicted to be stable:
    do i = 1, nSolnPhases
        k = -iAssemblage(nElements - i + 1)     ! Absolute solution phase index

        call CompStoichSolnPhase(k)

    end do

    ! If the new phase is a solution phase, compute the effective stoichiometry of that phase:
    if (iNewPhase < 0) call CompStoichSolnPhase(-iNewPhase)

    ! Return if an error was encoutered:
    if (INFOThermo /= 0) return

    ! Check the type of phase:
    IF_CheckType: if (iNewPhase > 0) then
        ! The new phase is a pure condensed phase.

        ! Compute the Euclidean norm for pure condensed phases:
        do i = 1, nConPhases
            dEuclideanNorm(i) = 0D0
            do l = 1, nElements
                dTemp             = dAtomFractionSpecies(iAssemblage(i),l) - dAtomFractionSpecies(iNewPhase,l)
                dEuclideanNorm(i) = dEuclideanNorm(i) + (dTemp)**2
            end do
        end do

        ! Compute the Euclidean norm for solution phases:
        do i = 1, nSolnPhases
            j = nElements - i + 1   ! Relative solution phase index
            k = -iAssemblage(j)     ! Absolute solution phase index
            dEuclideanNorm(j) = 0D0

            do l = 1, nElements
                dTemp             = dAtomFractionSpecies(iNewPhase,l) - dEffStoichSolnPhase(k,l)
                dEuclideanNorm(j) = dEuclideanNorm(j) + (dTemp)**2
            end do

            ! If this is a miscibility gap, store the phase index:
            if (lMiscibility(k)) iMisciblePhaseID = j

        end do

    elseif (iNewPhase < 0) then
        ! The new phase is a solution phase.

        ! Compute the Euclidean norm for pure condensed phases:
        do i = 1, nConPhases
            dEuclideanNorm(i) = 0D0
            do j = 1, nElements
                dTemp             = dAtomFractionSpecies(iAssemblage(i),j) - dEffStoichSolnPhase(-iNewPhase,j)
                dEuclideanNorm(i) = dEuclideanNorm(i) + (dTemp)**2
            end do
        end do

        ! Compute the Euclidean norm for solution phases:
        do i = 1, nSolnPhases
            j = nElements - i + 1   ! Relative solution phase index
            k = -iAssemblage(j)     ! Absolute solution phase index
            dEuclideanNorm(j) = 0D0

            do l = 1, nElements
                k                 = -iAssemblage(j)
                dTemp             = dEffStoichSolnPhase(k,l) - dEffStoichSolnPhase(-iNewPhase,l)
                dEuclideanNorm(j) = dEuclideanNorm(j) + (dTemp)**2
            end do

            ! If this is a miscibility gap, store the phase index:
            if (lMiscibility(k)) iMisciblePhaseID = j

        end do

    end if IF_CheckType

    ! If the Euclidean norm is extremely high, return control to the previous subroutine:
    if (MAXVAL(dEuclideanNorm) > 1D5) return

    ! Check if a pair of miscible phases are predicted to be stable.  If so, set the Euclidean norm to an
    ! arbitrarily high value to place these phases at the back of the list:
    if (iMisciblePhaseID /= 0) then

        do i = 1, nSolnPhases
            j = nElements - i + 1
            k = -iAssemblage(j)

            ! Set the Euclidean norm to an arbitrarily high value for all miscibile phases:
            if (j == iMisciblePhaseID) dEuclideanNorm(j) = 1000D0

            if (cSolnPhaseName(k) == cSolnPhaseName(-iAssemblage(iMisciblePhaseID))) dEuclideanNorm(j) = 1000D0

        end do

    end if

    ! Swap the phase with the lowest Euclidean norm for the first phase in the assemblage:
    IF_Euclid: if (iNewPhase /= 0) then

        dEuclideanNorm = -dEuclideanNorm

        ! Sort the Euclidean Norm vector:
        call SortPick(nElements, dEuclideanNorm, iVec)

        ! Determine what type of phase has the smallest Euclidean norm:
        if (iAssemblage(iVec(1)) < 0) then
            ! A solution phase has the smallest Euclidean norm.
            iPhaseTypeOut = 1
        else
            ! A pure condensed phase has the smallest Euclidean norm.
            iPhaseTypeOut = 0
        end if

        ! Store temporary variables:
        iVecB          = iAssemblage
        dEuclideanNorm = dMolesPhase

        ! Reinitialize variables:
        k           = 0
        l           = 0
        iAssemblage = 0
        dMolesPhase = 0d0

        ! Shuffle the phase assemblage:
        do i = 1, nElements
            j = iVec(i)
            if (iVecB(j) > 0) then
                l = l + 1
                iAssemblage(l) = iVecB(j)
                dMolesPhase(l) = dEuclideanNorm(j)
            elseif (iVecB(j) < 0) then
                k = k + 1
                iAssemblage(nElements - k + 1) = iVecB(j)
                dMolesPhase(nElements - k + 1) = dEuclideanNorm(j)
            end if
        end do
    end if IF_Euclid

    ! Check to make sure the phase assemblage contains the same selection of phases, otherwise revert:
    j = SUM(iAssemblage)

    if (iAssemblageSUM /= j) then
        iAssemblage = iTempVec
        dMolesPhase = dTempVec
    end if

    return

end subroutine ShuffleAssemblage
