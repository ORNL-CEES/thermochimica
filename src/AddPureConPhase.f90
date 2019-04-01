
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    AddPureConPhase.f90
    !> \brief   Add a pure condensed phase to the system.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckPureConPhaseAdd.f90
    !> \sa      CheckIterHistory.f90
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
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   07/12/2012      M.H.A. Piro         Updated the CheckIterationHistory subroutine to be more general.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to add a pure condensed phase to the estimated phase
    !!  assemblage.  The new phase assemblage is tested to ensure that it is appropriate for furhter consideration.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iPhaseChange    An integer scalar representing the pure condensed phase that is to be added
    !!                               to the system.
    !> \param[out]  lSwapLater      A logical scalar indicating whether the phase should be swapped with another
    !!                               phase later on in another subroutine.
    !> \param[out]  lPhasePass      A local scalar indicating whether the new estimated phase assemblage passed
    !!                               (.TRUE.) or failed (.FALSE.).
    !
    ! nConPhases            The number of pure condensed phases in the assemblage.
    ! iAssemblage           Integer vector containing the indices of phases estimated to be part of the
    !                        equilibrium phase assemblage.
    ! iAssemblageTest       An integer vector representing the phase assemblage to be tested.
    ! iterBack              An integer scalar representing how many iterations to go back in the iteration
    !                        history.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine AddPureConPhase(iPhaseChange,lSwapLater,lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                       :: i, iPhaseChange, iterBack, INFO
    integer, dimension(nElements) :: iAssemblageTest
    logical                       :: lPhasePass, lSwapLater


    ! Initialize variables:
    lPhasePass = .FALSE.
    lSwapLater = .FALSE.

    ! Ensure that this phase is not already part of the assemblage:
    do i = 1, nConPhases
        if (iAssemblage(i) == iPhaseChange) return
    end do

    ! If iPhaseChange was the last one removed, give the system a chance to converge:
    if ((iConPhaseLast == iPhaseChange).AND.(iterGlobal - iterLastCon < 3*iterStep)) return

    ! Check iteration history:
    if (iterGlobal > 50) then

        iAssemblageTest(1:nElements)  = iAssemblage(1:nElements)
        iAssemblageTest(nConPhases+1) = iPhaseChange
        if (iterGlobal < 1000) then
            iterBack = 200
        else
            iterBack = 1000
        end if

        ! Check whether this particular phase assemblage has been previously considered:
        call CheckIterHistory(iAssemblageTest,iterBack,lSwapLater)

    end if

    ! Proceed if the iteration history check passed:
    if (lSwapLater .EQV. .FALSE.) then

        ! Add the new pure condensed phase to the assemblage:
        nConPhases                  = nConPhases + 1
        iAssemblage(nConPhases)     = iPhaseChange

        ! Check to make sure that the phase change is acceptable:
        call CheckPhaseChange(lPhasePass,INFO)

        if (lPhasePass .EQV. .FALSE.) then
            ! The phase in question cannot be added. Revert the system:
            iAssemblage(nConPhases) = 0
            dMolesPhase(nConPhases) = 0D0
            nConPhases              = nConPhases - 1
            lSwapLater              = .TRUE.
        else
            ! This phase will be added to the assemblage.
            iConPhaseLast           = iPhaseChange
            iterLastCon             = iterGlobal
            iterLast                = iterGlobal
        end if

    end if

    return

end subroutine AddPureConPhase