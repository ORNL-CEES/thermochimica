
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    RemSolnAddPureConPhase.f90
    !> \brief   Simultaneously remove a particular solution phase and add a particular pure condensed phase.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckSolnPhaseRem.f90
    !> \sa      RemSolnPhase.f90
    !> \sa      CheckSysOnlyPureConPhases.f90
    !> \sa      CheckConvergence.f90
    !> \sa      CheckPhaseChange.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   02/29/2012      M.H.A. Piro         Original code
    !   03/19/2012      M.H.A. Piro         Corrected indexing problem when the phase assemblage fails and is
    !                                       returned to the previous value.  iPhaseRem wasn't used correctly.
    !   04/03/2012      M.H.A. Piro         The CheckSysOnlyPureConPhases subroutine was called if there were
    !                                       zero solution phases present, but the CheckConvergence subroutine
    !                                       needed to be called as well.
    !   04/26/2012      M.H.A. Piro         Convert to Gibbs energy minimization solver and dOxygen.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to remove a particular solution phase and add a particular
    !! pure condensed phase.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iPhaseAdd   An integer scalar representing the index of pure condensed phase to be added
    !!                           to the system.
    !> \param[in]   iPhaseRem   An integer scalar representing the ndex of solution phase to be removed from
    !!                           the system.
    !> \param[out]  lPhasePass  A logical variable indicating whether the new estimated phase assemblage passed
    !!                           (.TRUE.) or failed (.FALSE.).
    !
    ! nConPhases                The number of pure condensed phases in the assemblage
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine RemSolnAddPureConPhase(iPhaseAdd,iPhaseRem,lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::   j, iPhaseAdd, iPhaseRem, INFO, iTemp
    real(8)::   dTemp
    logical::   lPhasePass


    ! First, make sure that this phase is not already part of the assemblage.
    do j = 1, nConPhases
        if (iAssemblage(j) == iPhaseAdd) return
    end do

    ! Initialize variables:
    INFO       = 0
    lPhasePass = .FALSE.

    ! Remove solution phase:
    j                         = nElements - nSolnPhases + 1
    iTemp                     = iAssemblage(iPhaseRem)
    iAssemblage(iPhaseRem)    = iAssemblage(j)
    iAssemblage(j)            = 0
    dTemp                     = dMolesPhase(iPhaseRem)
    dMolesPhase(iPhaseRem)    = dMolesPhase(j)
    dMolesPhase(j)            = 0D0
    nSolnPhases               = nSolnPhases - 1

    ! Add pure condensed phase:
    iAssemblage(nConPhases+1) = iPhaseAdd
    nConPhases                = nConPhases + 1

    if (nSolnPhases == 0) then
        ! Check the system when there are only pure condensed phases:
        call CheckSysOnlyPureConPhases

        if (lConverged .EQV. .TRUE.) lPhasePass = .TRUE.

    else
        ! Check that this phase change is acceptable:
        call CheckPhaseChange(lPhasePass,INFO)

        lRevertSystem = .FALSE.

    end if

    if (lConverged .EQV. .FALSE.) then

        if (lPhasePass .EQV. .TRUE.) then
            ! This phase assemblage can be considered.
            iterLastSoln            = iterGlobal
            iterLast                = iterGlobal
            lSolnPhases(-iTemp)     = .FALSE.
        else
            ! This phase assemblage cannot be considered.  Revert back to the previous assemblage:
            iAssemblage(nConPhases) = 0
            nConPhases              = nConPhases - 1
            nSolnPhases             = nSolnPhases + 1
            j                       = nElements - nSolnPhases + 1
            iAssemblage(j)          = iAssemblage(iPhaseRem)
            dMolesPhase(j)          = dMolesPhase(iPhaseRem)
            iAssemblage(iPhaseRem)  = iTemp
            dMolesPhase(iPhaseRem)  = dTemp
        end if

    end if

    return

end subroutine RemSolnAddPureConPhase