
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    RemSolnPhase.f90
    !> \brief   Remove a solution phase from the system.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckSolnPhaseRem.f90
    !> \sa      CheckPhaseChange.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code
    !   10/26/2011      M.H.A. Piro         Clean up code: modules, fragment code into multiple subroutines.
    !   02/01/2012      M.H.A. Piro         Corrected a typo: when a solution phase is removed but is then
    !                                       determined that it cannot be removed, the solution phase index was
    !                                       not returned to its original value.
    !   04/26/2012      M.H.A. Piro         Convert to Gibbs energy minimization solver.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to remove a solution phase from the estimated phase assemblage.
    !! The CheckPhaseChange.f90 subroutine is called to ensure that the new phase assemblage is appropriate for
    !! further consideration.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iPhaseChange        An integer scalar representing the solution phase index corresponding to
    !!                                   the iAssemblage and dMolesPhase vectors.
    !> \param[out]  lPhasePass          A logical variable indicating whether the new phase assemblage has passed.
    !
    ! iterLast                          An integer scalar representing the global iteration number when the
    !                                    phase assemblage changed.
    ! nConPhases                        An integer scalar representing the number of pure condensed phases
    !                                    currently predicted to be stable at equilibrium.
    ! nSolnPhases                       An integer scalar representing the number of solution phases currently
    !                                    predicted to be stable at equilibrium.
    ! iAssemblage                       An integer vector representing the indices of phases currently predicted
    !                                    to be stable at equilibrium.
    ! dMolesPhase                       A double real vector representing the number of moles of phases
    !                                    predicted to be stable at equilibrium.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine RemSolnPhase(iPhaseChange,lPhasePass)

    USE ModuleThermoIO, ONLY: INFOTHermo
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::                      i, j, k, iPhaseChange, INFO, iSolnPhaseLastOld, nVar
    real(8)::                      dTemp
    real(8),dimension(nElements):: dTempVec
    logical::                      lPhasePass


    ! Initialize variables:
    lPhasePass = .FALSE.
    dTempVec   = dMolesPhase

    ! Remove this solution phase:
    j               = nElements - nSolnPhases  + 1
    k               = nElements - iPhaseChange + 1
    dTemp           = dMolesPhase(k)
    dMolesPhase(k)  = dMolesPhase(j)
    dMolesPhase(j)  = 0d0
    iSolnPhaseLastOld = iSolnPhaseLast
    iSolnPhaseLast  = -iAssemblage(k)           ! Absolute index of the solution phase that is to be removed.
    iAssemblage(k)  = iAssemblage(j)
    iAssemblage(j)  = 0

    nSolnPhases = 0
    do i =1, SIZE(iAssemblage)
        if (iAssemblage(i) < 0) nSolnPhases = nSolnPhases + 1
    end do

    if (nSolnPhases == 0) then
        ! Check the system if there are only pure condensed phases and no solution phases:
        call CheckSysOnlyPureConPhases

        if (INFOThermo == 27) INFOThermo = 0

        if (lConverged) then
            ! The system can be represented by only pure condensed phases.
            lPhasePass = .TRUE.
        else
            ! The system has not converged.  Return the num
            dMolesPhase = dTempVec
        end if
    else
        ! Check to make sure that the new phase assemblage is valid:
        j = MAX(1,nConPhases)
        LOOP_PhaseCheck: do i = 1, j

            nVar = nElements + nConPhases + nSolnPhases
            call ResizeGEMWorkspace(nVar)
            call CheckPhaseChange(lPhasePass,INFO)

            if ((INFO > nElements + nSolnPhases) .AND. (nConPhases > 0)) then
                ! A pure condensed phase should be removed.
                k                       = INFO - nElements - nSolnPhases
                iConPhaseLast           = iAssemblage(k)
                iAssemblage(k)          = iAssemblage(nConPhases)
                dMolesPhase(k)          = dMolesPhase(nConPhases)
                iAssemblage(nConPhases) = 0
                dMolesPhase(nConPhases) = 0D0
                nConPhases              = nConPhases - 1
            else
                ! This phase assemblage is acceptable.
                exit LOOP_PhaseCheck
            end if
        end do LOOP_PhaseCheck
    end if

    if (.NOT.(lConverged)) then
        if (lPhasePass) then
            ! The new phase assemblage has passed.
            lSolnPhases(iSolnPhaseLast) = .FALSE.
            iterLastSoln    = iterGlobal
            iterLast        = iterGlobal

            call CheckRemMisciblePhase(iSolnPhaseLast)

            ! Compute the number of moles of all phases:
            !call CompMolAllSolnPhases

        else
            ! The phase in question cannot be removed. Return variables to their previous values:
            nSolnPhases     = nSolnPhases + 1
            j               = nElements - nSolnPhases  + 1
            k               = nElements - iPhaseChange + 1
            iAssemblage(j)  = iAssemblage(k)
            iAssemblage(k)  = -iSolnPhaseLast
            dMolesPhase(j)  = dMolesPhase(k)
            dMolesPhase(k)  = dTemp
            iSolnPhaseLast  = iSolnPhaseLastOld

            if (iterGlobal - iterLast >= 40) lRevertSystem = .TRUE.

        end if
    else
        ! Placeholder: the system is comprised of only pure condensed phases and the system
        ! has converged.  Do nothing.
    end if

    return

end subroutine RemSolnPhase




!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to check if the phase that is being
    ! removed is miscible and if one of the other corresponding phases is
    ! currently predicted to be stable.  If so, this subroutine checks if the
    ! absolute index of the phase that is to be removed is less than the other
    ! phase, in which case it swaps the two phases.  The motivation for doing
    ! this is to have a consistent record of the phase assemblage in the
    ! variable iterHistory.
    !
    !
    ! iPhaseRemoveAbs - Absolute index of solution phase that is being removed.
    ! iPhaseOtherAbs  - Absolute index of solution phase that has a miscibility
    !                    gap with iPhaseRemoveAbs (default set to zero).
    ! iPhaseOtherRel  - Relative index of solution phase corresponding to
    !                    iPhaseOtherAbs.
    !
    !---------------------------------------------------------------------------

subroutine CheckRemMisciblePhase(iPhaseRemoveAbs)

    USE ModuleThermo
    USE ModuleGEMSolver, ONLY: lSolnPhases

    implicit none

    integer:: i, j, k, l, iPhaseRemoveAbs, iPhaseOtherAbs, iPhaseOtherRel


    ! Initialize variables:
    iPhaseOtherAbs = 0
    iPhaseOtherRel = 0

    ! Only proceed if there is at least one other solution phase:
    if (nSolnPhases > 0) then

        ! First, see if any of the other solution phases currently assumed to be
        ! stable correspond to this miscible phase:
        LOOP_A: do j = 1, nSolnPhases
            k = -iAssemblage(nElements-j+1)

            if (cSolnPhaseName(k) == cSolnPhaseName(iPhaseRemoveAbs)) then
                iPhaseOtherRel = j
                iPhaseOtherAbs = k
                exit LOOP_A
            end if

        end do LOOP_A

        if (iPhaseOtherRel /= 0) then
            if (iPhaseOtherAbs > iPhaseRemoveAbs) then
                ! Store indices of solution species for both phases:
                i = nSpeciesPhase(iPhaseOtherAbs-1) + 1
                j = nSpeciesPhase(iPhaseOtherAbs)
                k = nSpeciesPhase(iPhaseRemoveAbs-1) + 1
                l = nSpeciesPhase(iPhaseRemoveAbs)

                ! Swap mole fractions, moles and chemical potential terms:
                dMolFraction(k:l)       = dMolFraction(i:j)
                dMolesSpecies(k:l)      = dMolesSpecies(i:j)
                dChemicalPotential(k:l) = dChemicalPotential(i:j)

                ! Reset the moles of solution species for the other phase:
                dMolesSpecies(i:j) = 0D0

                ! Swap the phase indices in iAssemblage:
                iAssemblage(nElements-iPhaseOtherRel+1) = -iPhaseRemoveAbs

                ! Swap the logical variables for lSolnPhases:
                lSolnPhases(iPhaseRemoveAbs) = .TRUE.
                lSolnPhases(iPhaseOtherAbs)  = .FALSE.

            end if

        end if

    end if

end subroutine CheckRemMisciblePhase


    !---------------------------------------------------------------------------
    !                       END - RemSolnPhase.f90
    !---------------------------------------------------------------------------
