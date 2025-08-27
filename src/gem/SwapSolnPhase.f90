
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    SwapSolnPhase.f90
    !> \brief   Swap a particular solution phase for another solution phase.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckSolnPhaseAdd.f90
    !> \sa      SwapSolnForPureConPhase.f90
    !> \sa      CompMolSolnPhase.f90
    !> \sa      CheckPhaseChange.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/05/2012      M.H.A. Piro         Original code
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   05/09/2012      M.H.A. Piro         Correct indexing error
    !   07/12/2012      M.H.A. Piro         Add iteration history check.
    !   08/22/2012      M.H.A. Piro         Do not allow one phase to swap another if they both belong to the
    !                                        same phase with a miscibility gap.
    !   09/29/2012      M.H.A. Piro         Fix indexing error in dMolesPhase when computed dMolesSpecies.
    !   10/19/2012      M.H.A. Piro         Check if the phase to be removed has a corresponding phase with
    !                                        a miscibility gap, and if so, make sure that the phase with the
    !                                        highest absolute index is removed.
    !   10/03/2015      M.H.A. Piro         iterBack depends on iterGlobal.  The motivation is that if the
    !                                        iteration count gets really high, it might have to go further back.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to attempt to swap a particular solution phase for another
    !! solution phase in the current estimated phase assemblage.  The logical variable lPhasePass returns a
    !! value of FALSE if the phase in question cannot swap for any other solution phase in the system.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iPhaseChange        An integer scalar representing the index of a solution phase that is
    !!                                   to be introduced to the system.
    !> \param[out]  lPhasePass          A logical variable indicatin whether the phase assemblage is appropriate
    !!                                   for consdieration (TRUE) or not (FALSE).
    !
    ! nElements             An integer scalar representing the number of elements in the system.
    ! nSolnPhases           An integer scalar representing the number of solution phases in the assemblage.
    ! iAssemblage           Integer vector containing the indices of all phases currently estimated to contribute
    !                       to the equilibrium phase assemblage.
    ! iPhaseChange
    ! iSolnPhaseLast        Index of the last solution phase to be added to, withdrawn from, or exchanged, from
    !                       the estimated phase assemblage.
    ! dMolesPhase           A double real vector representing the number of moles of each phase.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine SwapSolnPhase(iPhaseChange,lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                       :: i, j, k, l, iPhaseChange, INFO, iterBack, nVar
    integer, dimension(nElements) :: iAssemblageTest, iAssemblageTemp
    real(8)                       :: dTemp
    real(8), dimension(nElements) :: dMolesPhaseTemp
    real(8), dimension(nSpecies)  :: dMolFractionTemp, dMolesSpeciesTemp
    logical                       :: lPhasePass, lSwapLater, lCompEverything

    ! Initialize variables:
    lCompEverything = .FALSE.
    lPhasePass      = .FALSE.
    lSwapLater      = .FALSE.

    if (iterGlobal > 1500) then
        iterBack = 1000
    else
        iterBack = 200
    end if

    ! CONSIDER: THE FOLLOWING CALL TO COMPMOLFRACTION IS MOST LIKELY UNNECSSARY:

    ! Compute the mole fractions of solution species in this phase:
    call CompMolFraction(iPhaseChange)

    ! Compute the stoichiometry of this phase:
    call CompStoichSolnPhase(iPhaseChange)

    dMolesPhaseTemp   = dMolesPhase
    dMolFractionTemp  = dMolFraction
    dMolesSpeciesTemp = dMolesSpecies
    iAssemblageTemp   = iAssemblage

    ! Loop through all solution phases in the current estimated phase assemblgae to see which one should
    ! be swapped:
    LOOP_SolnPhase: do i = 1, nSolnPhases
        j = -iAssemblage(nElements - i + 1)
        ! Check if this phase has a miscibility gap, and if so, make sure that it does not swap itself:
        if ((lMiscibility(iPhaseChange)).OR.(lMiscibility(j))) then
            ! Either the phase that is to be added to the system or the current phase is not the first
            ! "phase" that contains a miscibility gap.
            ! Check if they belong to the same "phase", and if so, skip to the next phase:
            if (cSolnPhaseName(j) == cSolnPhaseName(iPhaseChange)) cycle LOOP_SolnPhase
        end if

        ! Check if this phase has a corresponding phase with a miscibility gap, and if so, make sure that
        ! the phase with the highest absolute index is used:
        LOOP_CheckMiscible: do k = 1, nSolnPhases
            ! Cycle if they are the same phase:
            if (k == i) cycle LOOP_CheckMiscible
            ! Store the absolute phase index of k:
            l = -iAssemblage(nElements-k+1)
            ! Check if the pair of phases contain a miscibility gap:
            if (cSolnPhaseName(l) == cSolnPhaseName(j)) then
                ! If the absolute phase index of j is less than l, cycle to the next phase:
                if (j < l) cycle LOOP_SolnPhase
            end if
        end do LOOP_CheckMiscible

        ! Check if this phase assemblage has been previously considered:
        if (iterGlobal > 60) then
            j                  = nElements - i + 1
            iAssemblageTest    = iAssemblage
            iAssemblageTest(j) = -iPhaseChange
            ! Check whether this particular phase assemblage has been previously considered:
            call CheckIterHistory(iAssemblageTest,iterBack,lSwapLater)
            ! This phase assemblage has been considered.  Move on to the next phase:
            if (lSwapLater) cycle LOOP_SolnPhase
        end if

        k              = nElements - i + 1
        dTemp          = dMolesPhase(k)
        iSolnPhaseLast = -iAssemblage(k)
        iAssemblage(k) = -iPhaseChange

        ! Compute the number of moles of each solution phase and establish the Jacobian constraint vector
        call CompMolAllSolnPhases

        ! Compute the number of moles of solution species:
        k = nElements - i + 1
        do j = nSpeciesPhase(iPhaseChange-1) + 1, nSpeciesPhase(iPhaseChange)
            dMolesSpecies(j) = dMolFraction(j) * dMolesPhase(k)
        end do

        ! Compute the chemical potentials:
        call CompChemicalPotential(lCompEverything)

        ! Check that this phase change is acceptable:
        nVar = nElements + nConPhases + nSolnPhases
        call ResizeGEMWorkspace(nVar)
        call CheckPhaseChange(lPhasePass,INFO)

        lRevertSystem = .FALSE.

        if (lPhasePass) then
            ! This phase assemblage can be considered.
            iterLastSoln = iterGlobal
            iterLast     = iterGlobal
            iterSwap     = iterGlobal
            iSolnSwap    = iPhaseChange
            iPureConSwap = 0

            dDrivingForceSoln(iPhaseChange) = 0D0
            lSolnPhases(iPhaseChange)       = .TRUE.
            lSolnPhases(iSolnPhaseLast)     = .FALSE.
            exit LOOP_SolnPhase
        else
            ! This phase assemblage cannot be considered.  Revert back to the previous assemblage:
            dMolesPhase   = dMolesPhaseTemp
            dMolFraction  = dMolFractionTemp
            dMolesSpecies = dMolesSpeciesTemp
            iAssemblage   = iAssemblageTemp
            dPartialExcessGibbs = dPartialExcessGibbsLast
            cycle LOOP_SolnPhase
        end if

    end do LOOP_SolnPhase

    return

end subroutine SwapSolnPhase
