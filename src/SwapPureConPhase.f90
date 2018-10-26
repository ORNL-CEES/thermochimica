
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    SwapPureConPhase.f90
    !> \brief   Swap a pure condensed phase for another pure condensed phase.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckPureConPhaseAdd.f90
    !> \sa      AddPureConPhase.f90
    !> \sa      SwapPureConForSolnPhase.f90
    !> \sa      CheckPhaseChange.f90
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code
    !   10/28/2011      M.H.A. Piro         Clean up code: modules, fragment code into multiple subroutines.
    !   02/13/2012      M.H.A. Piro         Check if iPhaseChange is already part of the assemblage.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   08/01/2012      M.H.A. Piro         Fix iteration history check.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to swap one pure condensed phase for another pure condensed 
    !! phase in the estimated phase assemblage .
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iPhaseChange        An integer scalar representing the absolute phase index of a pure 
    !!                                   condensed phase.
    !> \param[out]  lPhasePass          A logical variable indicating whether the phase assemblage has passed or
    !!                                   failed.
    !> \param[out]  lSwapLater          A logical variable indicating whether the phase should be swapped later.
    !
    ! nConPhases                        The number of pure condensed phases in the assemblage.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine SwapPureConPhase(iPhaseChange,lSwapLater,lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                       :: i, iPhaseChange, INFO, iterBack
    integer, dimension(nElements) :: iAssemblageTest
    real(8), dimension(nElements) :: dTempVec
    logical                       :: lSwapLater, lPhasePass   


    ! Initialize variables:
    lPhasePass = .FALSE.
    iterBack   = 300

    ! Ensure that this phase is not already part of the assemblage:
    do i = 1,nConPhases
        if (iAssemblage(i) == iPhaseChange) return
    end do

    LOOP_ConPhase: do i = 1, nConPhases
                        
        ! Check the iteration history:
        if (iterGlobal > 50) then
            iAssemblageTest    = iAssemblage
            iAssemblageTest(i) = iPhaseChange
            
            ! Check whether this particular phase assemblage has been previously considered:
            call CheckIterHistory(iAssemblageTest,iterBack,lSwapLater)
            
            ! If this phase assemblage has been previously considered, skip to the next phase:
            if (lSwapLater .EQV. .TRUE.) cycle LOOP_ConPhase
            
        end if
        
        ! Store temporary variables:
        dTempVec        = dMolesPhase
        iAssemblageTest = iAssemblage
        iConPhaseLast   = iAssemblage(i)
        iAssemblage(i)  = iPhaseChange
        
        ! Check that this phase change is acceptable:                    
        call CheckPhaseChange(lPhasePass,INFO)
        
        if (lPhasePass .EQV. .FALSE.) then
            ! The phase in question cannot be swapped.  Revert back to the previous assemblage:                
            iAssemblage = iAssemblageTest
            dMolesPhase = dTempVec
            lSwapLater  = .TRUE.
            cycle LOOP_ConPhase
        else
            ! The new phase assemblage is acceptable.
            iterLastCon    = iterGlobal
            iterLast       = iterGlobal
            iterSwap       = iterGlobal
            iPureConSwap   = iPhaseChange
            iSolnSwap      = 0 
            dMolesPhase(i) = 0D0
            exit LOOP_ConPhase
        end if
        
    end do LOOP_ConPhase
    
    return

end subroutine SwapPureConPhase