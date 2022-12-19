
!-------------------------------------------------------------------------------------------------------------
    !
    !> \file    RemPureConPhase.f90
    !> \brief   Remove a pure condensed phase from the system.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckPureConPhaseRem.f90
    !> \sa      CheckPhaseChange.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code
    !   10/25/2011      M.H.A. Piro         Clean up code: modules, fragment code into multiple subroutines.
    !   03/25/2012      M.H.A. Piro         Added the capability of checking the iteration history when there
    !                                       is more than one pure condensed phase that needs to be removed.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   05/30/2012      M.H.A. Piro         Store the previous phase assemblage info in case if a pure con phase
    !                                        is initially removed and one (or more) additional phases are removed,
    !                                        but only the first phase is returned to the system.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to remove a pure condensed phase from the estimated phase
    !! assemblage.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iPhaseChange    An integer scalar representing the index of the pure condensed phase to be
    !!                               removed from the system.
    !> \param[out]  lSwapLater      A logical scalar indicating whether this phase should be swapped later.
    !> \param[out]  lPhasePass      A logical scalar indicating whether the new phase assemblage has passed.
    !
!-------------------------------------------------------------------------------------------------------------


subroutine RemPureConPhase(iPhaseChange,lSwapLater,lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::                        i, j, k, iPhaseChange, INFO, nConPhasesLast
    integer,dimension(nConPhases)::  iTempVec
    real(8)::                        dTemp
    real(8),dimension(nConPhases)::  dTempVec
    logical::                        lSwapLater, lPhasePass


    ! Initialize variables:
    lSwapLater = .FALSE.
    lPhasePass = .FALSE.

    ! Store the phase assemblage and number of moles in a temporary vector:
    iTempVec(1:nConPhases) = iAssemblage(1:nConPhases)
    dTempVec(1:nConPhases) = dMolesPhase(1:nConPhases)
    nConPhasesLast         = nConPhases

    ! If iPhaseChange was the last one added, give the system a chance to converge:
    if ((iConPhaseLast == iAssemblage(iPhaseChange)).AND.(iterGlobal - iterLastCon < iterStep).AND.(.NOT. lConverged)) return
    print *, 'rem con ', iAssemblage
    print *, 'rem con ', dMolesPhase

    ! Remove the pure condensed phase corresponding to iPhaseChange:
    iConPhaseLast             = iAssemblage(iPhaseChange)
    iAssemblage(iPhaseChange) = iAssemblage(nConPhases)
    iAssemblage(nConPhases)   = 0
    dTemp                     = dMolesPhase(iPhaseChange)
    dMolesPhase(iPhaseChange) = dMolesPhase(nConPhases)
    dMolesPhase(nConPhases)   = 0D0
    nConPhases                = nConPhases - 1
    if (.NOT. lConverged) dMolesPhase = dMolesPhase * 0.95D0
    print *, 'rem con ', iAssemblage

    ! Check that this phase change is acceptable:
    k = MAX(1, nConPhases)
    do i = 1, k

        call CheckPhaseChange(lPhasePass,INFO)

        if ((INFO > nElements + nSolnPhases) .AND. (nConPhases > 0)) then
            ! A pure condensed phase should be removed.
            j = INFO - nElements - nSolnPhases
            iAssemblage(j)          = iAssemblage(nConPhases)
            dMolesPhase(j)          = dMolesPhase(nConPhases)
            iAssemblage(nConPhases) = 0
            dMolesPhase(nConPhases) = 0D0
            nConPhases              = nConPhases - 1
            iterLast                = iterGlobal
        else
            ! This phase assemblage is appropriate.
            exit
        end if
    end do

    if (lPhasePass) then
        ! The new phase assemblage is acceptable.
        iterLastCon               = iterGlobal
        iterLast                  = iterGlobal
        dMolesPhase               = DABS(dMolesPhase)

        ! Compute the number of moles of all phases:
        call CompMolAllSolnPhases
    else
        ! The phase in question cannot be removed.  Add it back to the assemblage:
        nConPhases                = nConPhasesLast
        iAssemblage(1:nConPhases) = iTempVec(1:nConPhases)
        dMolesPhase(1:nConPhases) = dTempVec(1:nConPhases)
        lSwapLater                = .TRUE.
    end if

    return

end subroutine RemPureConPhase
