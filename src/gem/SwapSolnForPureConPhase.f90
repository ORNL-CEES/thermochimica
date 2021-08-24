
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    SwapSolnForPureConPhase.f90
    !> \brief   Swap a particular solution phase for a pure condensed phase.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckSolnPhaseAdd.f90
    !> \sa      CompMolSolnPhase.f90
    !> \sa      CheckPhaseChange.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/20/2012      M.H.A. Piro         Original code
    !   02/28/2012      M.H.A. Piro         Check iteration history.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   05/25/2012      M.H.A. Piro         Apply a minimum number of moles of solution species.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to add a particular solution phase and remove a pure condnesed
    !! phase.  This subroutine will determine which pure condensed phase is to be removed from the system.
    !! The iteration history is checked to ensure that this particular phase assemblage has not been previously
    !! considered.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iPhaseChange    An integer scalar representing the index of a solution phase that is to be
    !!                               introduced to the system.
    !> \param[out]  lPhasePass      A logical variable indicating whether a particular phase assemblage is
    !!                               appopriate for testing (TRUE) or not (FALSE).
    !
    ! nConPhases                    An integer scalar representing the number of pure condensed phases in the
    !                                assemblage.
    ! nSolnPhases                   An integer scalar representing the number of solution phases in the
    !                                assemblage.
    ! iAssemblage                   An integer vector containing the indices of all phases currently estimated
    !                                to contribute to the equilibrium phase assemblage.
    ! iSolnPhaseLast                Index of the last solution phase to be added to, withdrawn from, or
    !                                exchanged, from the estimated phase assemblage.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine SwapSolnForPureConPhase(iPhaseChange,lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                       :: i, j, k, iPhaseChange, INFO, iterBack
    integer, dimension(nElements) :: iAssemblageTest, iTempVec
    real(8), dimension(nElements) :: dTempVec
    logical                       :: lPhasePass, lSwapLater, lCompEverything


    ! Initialize variables:
    lPhasePass            = .FALSE.
    lSwapLater            = .FALSE.
    lCompEverything       = .FALSE.
    iTempVec(1:nElements) = iAssemblage(1:nElements)

    ! Check if this phase contains a miscibility gap:
    if (lMiscibility(-iPhaseChange)) then

        ! Loop through phases currently predicted to be stable and see if it coresponds to iPhaseChange:
        LOOP_A: do i = 1, nSolnPhases
            j = -iAssemblage(nElements - i + 1)

            if (cSolnPhaseName(j) == cSolnPhaseName(-iPhaseChange)) then
                ! Phase i and -iPhaseChange represent a miscibility gap:
                lSwapLater = .TRUE.
                exit LOOP_A
            end if

        end do LOOP_A

        ! If phase -iPhaseChange
        if (.NOT.(lSwapLater)) call CompMolFraction(-iPhaseChange)

        ! Reinitialize variables:
        lSwapLater = .FALSE.
    else

        ! Compute the mole fractions of solution species in this phase:
        call CompMolFraction(-iPhaseChange)

    end if

    ! If the current phase has a miscibility gap, check the index numbers:
    k = -iPhaseChange
    if (lMiscibility(k)) call CheckAddMisciblePhaseIndex(k)
    iPhaseChange = -k

    ! Compute the stoichiometry of this phase:
    call CompStoichSolnPhase(-iPhaseChange)

    LOOP_ConPhase: do i = 1, nConPhases

        ! Check iteration history:
        iAssemblageTest(1:nElements)             = iTempVec(1:nElements)
        iAssemblageTest(i)                       = iTempVec(nConPhases)
        iAssemblageTest(nConPhases)              = 0
        iAssemblageTest(nElements - nSolnPhases) = iPhaseChange

        if (iterGlobal < 1000) then
            iterBack = 200
        else
            iterBack = 1000
        end if

        ! Check whether this particular phase assemblage has been previously considered:
        call CheckIterHistory(iAssemblageTest,iterBack,lSwapLater)

        if (lSwapLater) cycle LOOP_ConPhase

        ! Store the info for the pure condensed phase to be removed to temporary variables:
        dTempVec        = dMolesPhase
        iAssemblageTest = iAssemblage
        iConPhaseLast   = iAssemblage(i)

        ! Move the info for the last pure conndesed phase to the location being removed:
        dMolesPhase(i) = dMolesPhase(nConPhases)
        iAssemblage(i) = iAssemblage(nConPhases)

        ! Reset information for the last pure condensed phase in the assemblage:
        dMolesPhase(nConPhases) = 0D0
        iAssemblage(nConPhases) = 0

        ! Add the new solution phase:
        nConPhases     = nConPhases - 1
        nSolnPhases    = nSolnPhases + 1
        j              = nElements - nSolnPhases + 1
        iAssemblage(j) = iPhaseChange

        ! Compute the number of moles of each solution phase and establish the Jacobian constraint vector
        !call CompMolSolnPhase
        call CompMolAllSolnPhases

        ! Compute the number of moles of solution species:
        k = nElements - nSolnPhases + 1
        do j = nSpeciesPhase(-iPhaseChange-1) + 1, nSpeciesPhase(-iPhaseChange)
            dMolesSpecies(j) = dMolFraction(j) * dMolesPhase(k) / dSumMolFractionSoln(-iPhaseChange)
            dMolesSpecies(j) = DMAX1(dMolesSpecies(j), dTolerance(8))
        end do

        ! Compute the chemical potentials:
        call CompChemicalPotential(lCompEverything)

        ! Check that this phase change is acceptable:
        call CheckPhaseChange(lPhasePass,INFO)

        lRevertSystem = .FALSE.

        if (lPhasePass) then
            ! This phase assemblage can be considered.
            iterLast     = iterGlobal
            iterSwap     = iterGlobal
            iSolnSwap    = -iPhaseChange
            iPureConSwap = 0
            lPhasePass   = .TRUE.

            dDrivingForceSoln(-iPhaseChange) = 0D0
            lSolnPhases(-iPhaseChange)       = .TRUE.

            exit LOOP_ConPhase
        else
            ! This phase assemblage cannot be considered.  Revert back to the previous assemblage:
            nConPhases  = nConPhases + 1
            nSolnPhases = nSolnPhases - 1
            iAssemblage = iAssemblageTest
            dMolesPhase = dTempVec
            cycle LOOP_ConPhase
        end if

    end do LOOP_ConPhase

    return

end subroutine SwapSolnForPureConPhase
