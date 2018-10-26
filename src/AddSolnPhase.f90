
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    AddSolnPhase.f90
    !> \brief   Add a solution phase to the system.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckSolnPhaseAdd.f90
    !> \sa      CompMolSolnPhase.f90
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code
    !   10/26/2011      M.H.A. Piro         Clean up code: modules, fragment code into multiple subroutines.
    !   03/24/2012      M.H.A. Piro         Added iteration history check when pure condensed phases are removed
    !                                       after the call to CheckPhaseChange.f90.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   07/30/2012      M.H.A. Piro         Do not proceed if the number of moles of any phase is negative.
    !   09/05/2012      M.H.A. Piro         Only call CompMolFraction if the phase is not miscible.  The 
    !                                        reason why this is done is because the subroutine CompMolFraction
    !                                        computes the mole fractions based on the element potentials and 
    !                                        the partial excess Gibbs energy of mixing from the last 
    !                                        iteration; whereas a phase with a miscibility gap should use the 
    !                                        mole fractions computed from Subminimization.
    !   09/07/2012      M.H.A. Piro         Remove call to compute the effective stoichiometry of the solution
    !                                        solution phase that is added because it is calculated in the call
    !                                        to CompMolAllSolnPhases.
    !   09/21/2012      M.H.A. Piro         Added the CheckAddMisciblePhase subroutine.
    !   09/26/2015      M.H.A. Piro         Added a call to CheckMiscibilityGap to ensure that the phase has 
    !                                        a composition that corresponds to a global minimum rather than
    !                                        a local minima.
    !   09/27/2015      M.H.A. Piro         Apply a different value for iterBack when the global iteration 
    !                                        count is above 1500.  The motivation is that sometimes the
    !                                        correct phase assemblage may be stable, but getting to equilibrium
    !                                        depends on the order that phases are added to the system.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to add a solution phase to the estimated phase assemblage.  
    !! The new phase assemblage is tested to ensure that it is appropriate and a logical variable lPhasePass is 
    !! returned.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iPhaseChange        An integer scalar representing the absolute index of a solution phase
    !!                                   to be added to the system.
    !> \param[out]  lSwapLater          A logical variable indicating whether a phase should be swapped for  
    !!                                   another phase when this particular phase cannot be added directly to
    !!                                   the system.
    !> \param[out]  lPhasePass          A logical variable indicating whether the new estimated phase assemblage 
    !!                                   passed (.TRUE.) or failed (.FALSE.).
    !
    ! nConPhases                        The number of pure condensed phases in the assemblage.
    ! iAssemblage                       An integer vector containing the indices of phases estimated to be part
    !!                                   of the equilibrium phase assemblage.
    ! dMolesPhase                       A double real vector containing the number of moles of each phase.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine AddSolnPhase(iPhaseChange,lSwapLater,lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::                       i, j, k, iPhaseChange, INFO, nConPhasesLast, iterBack
    integer,dimension(nElements)::  iAssemblageTest, iAssemblageLast
    real(8),dimension(nElements)::  dTempVec
    logical::                       lPhasePass, lSwapLater, lCompEverything


    ! Initialize variables:
    iAssemblageLast = iAssemblage
    nConPhasesLast  = nConPhases 
    dTempVec        = dMolesPhase
    lSwapLater      = .FALSE.
    lPhasePass      = .FALSE.
    lCompEverything = .FALSE.

    if (iterGlobal > 1500) then
        iterBack = 1000
    else
        iterBack = 200
    end if

    ! If this phase was the last phase to be removed, then give the system a chance to converge:
    if ((iPhaseChange == iSolnPhaseLast).AND.(iterLast == iterLastSoln).AND. & 
        (iterGlobal - iterLastSoln <= iterStep*2)) return
        
    ! Do not try to add a solution phase if any pure condensed phases are negative:
    if ((MINVAL(dMolesPhase(1:nConPhases)) < 0D0).AND.(iterGlobal - iterLast < 20)) return

    ! If the current phase has a miscibility gap, check the index numbers:
    if (lMiscibility(iPhaseChange) .EQV. .TRUE.) then

        call CheckAddMisciblePhaseIndex(iPhaseChange)

    else

        ! Ensure that the mole fractions of this solution phase correspond to a global minimum, rather than a local minima:
        call CheckMiscibilityGap(iPhaseChange,lSwapLater)
        lSwapLater = .FALSE.

    end if

    ! Check if this phase assemblage has been previously considered:
    if (iterGlobal > 60) then
    
        j                  = nElements - nSolnPhases
        iAssemblageTest    = iAssemblage 
        iAssemblageTest(j) = -iPhaseChange

        ! Check whether this particular phase assemblage has been previously considered:
        call CheckIterHistory(iAssemblageTest,iterBack,lSwapLater)

        ! This phase assemblage has been considered.  Move on to the next phase:
        if (lSwapLater .EQV. .TRUE.) return
        
    end if
    
    ! Add this solution phase to the assemblage:
    nSolnPhases                              = nSolnPhases + 1
    iAssemblage(nElements - nSolnPhases + 1) = -iPhaseChange
    
    ! Compute the number of moles of all solution phases:
    call CompMolSolnPhase
    
    ! Absolute solution phase index:
    k = nElements - nSolnPhases + 1

    ! Compute the number of moles of solution species:
    do i = nSpeciesPhase(iPhaseChange-1) + 1, nSpeciesPhase(iPhaseChange)
        dMolesSpecies(i) = dMolFraction(i) * dMolesPhase(k) 
    end do
    
    ! Make sure that the new phase assemblage yields an appropriate Hessian matrix:
    j = MAX(1,nConPhases)
    
    LOOP_CheckAssemblage: do i = 1, j
    
        ! Check that the phase change is acceptable:
        call CheckPhaseChange(lPhasePass,INFO)
                
        IF_CheckINFO: if (INFO > nElements + nSolnPhases) then
            
            ! Remove a pure condensed phase:
            k                       = INFO - nElements - nSolnPhases
            iAssemblage(k)          = iAssemblage(nConPhases)
            iAssemblage(nConPhases) = 0
            dMolesPhase(k)          = dMolesPhase(nConPhases)
            dMolesPhase(nConPhases) = 0D0
            nConPhases              = nConPhases - 1
            
            ! Check if this phase assemblage has been previously considered:
            if (iterGlobal > 60) then
                
                iAssemblageTest    = iAssemblage 

                ! Check whether this particular phase assemblage has been previously considered:
                call CheckIterHistory(iAssemblageTest,iterBack,lSwapLater)
            
                ! This phase assemblage has been considered.  Move on to the next phase:
                if (lSwapLater .EQV. .TRUE.) then
                    lPhasePass = .FALSE.
                    exit LOOP_CheckAssemblage
                end if 
            end if
        else
            ! This phase assemblage is acceptable.
            exit LOOP_CheckAssemblage
            
        end if IF_CheckINFO
        
    end do LOOP_CheckAssemblage
             
    ! Check if this new phase assemblage is acceptable:
    if (lPhasePass .EQV. .FALSE.) then
        ! The phase in question cannot be removed.  Revert the system:
        nSolnPhases         = nSolnPhases - 1        
        iAssemblage         = iAssemblageLast
        dMolesPhase         = dTempVec
        nConPhases          = nConPhasesLast
        lSwapLater          = .TRUE.
        
        ! The following is (probably) unnecessary:
        dPartialExcessGibbs = dPartialExcessGibbsLast
    else
        ! The new phase assemblage is appropriate for testing.
        iterLastSoln   = iterGlobal
        iterlast       = iterGlobal
        iSolnPhaseLast = iPhaseChange
        lPhasePass     = .TRUE.
        
        dDrivingForceSoln(iPhaseChange) = 0D0
        lSolnPhases(iPhaseChange)       = .TRUE.   
    end if
    
    return
    
end subroutine AddSolnPhase
