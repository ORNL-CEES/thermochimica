
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckSolnPhaseRem.f90
    !> \brief   Check whether a solution phase needs to be removed from the system.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckPhaseAssemblage.f90
    !> \sa      RemSolnPhase.f90
    !> \sa      RemSolnAddPureConPhase.f90
    !> \sa      CompDrivingForce.f90
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/05/2012      M.H.A. Piro         Original code.
    !   04/26/2012      M.H.A. Piro         Convert to Gibbs energy minimization solver.
    !   05/02/2012      M.H.A. Piro         If a solution phase swapped for another phase and is then 
    !                                        immediately withdrawn from the system, revert to the previously 
    !                                        successful iteration.
    !   08/15/2012      M.H.A. Piro         If a solution phase cannot be removed directly from the system
    !                                        and it cannot be swapped for a pure condensed phase or another
    !                                        solution phase, then try removing each pure condensed phase in 
    !                                        the system.
    !   10/04/2012      M.H.A. Piro         Undo the last change and instead revert the system to the last 
    !                                        successful phase assemblage.
    !   11/07/2012      M.H.A. Piro         If a solution phase cannot be removed directly from the system
    !                                        and it cannot be swapped for a pure condensed phase or another
    !                                        solution phase, check if the last time the phase assemblgae 
    !                                        changed corresponded to a swap.  Try swapping the two phases.
    !   11/24/2012      M.H.A. Piro         If a solution phase cannot be removed directly from the system
    !                                        and a pure condensed phase was recently swapped, revert to 
    !                                        the last successful iteration.
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check whether a solution phase should be removed from the 
    !! system.  A solution phase should be removed from the sytem if the number of moles of that phase is 
    !! less than a certain tolernace (specified in InitThermo.f90).    It may be possible that this phase may
    !! not be removed directly from the system because it would result in a singularity.  Provisions are in 
    !! place to check if a pure condensed phase should take the place of this solution phase.  If this does not
    !! work, then a check is performed to see if a different solution phase should take the place of this phase.
    !! If this does not work and the number of moles of this solution phase is sufficiently small, then pure
    !! condensed phases are systematically removed from the system in a last attempt to converge.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nSolnPhases           The number of solution phases in the assemblage.
    ! nSolnPhasesSys        The number of solution phases in the database.
    ! nElements             The number of elements in the system.
    ! lPhasePass            A logical variable indicating whether the new estimated phase assemblage passed 
    !                        (.TRUE.) or failed (.FALSE.).
    ! dChemicalPotential    A double real vector representing the chemical potential of each species.  To be
    !                        precise, this is defined as the difference between the standard molar Gibbs energy 
    !                        and the chemical potential defined by the element potentials (represented in 
    !                        dimensionless units and per formula mass).
    ! dMolesPhase           A double real vector representing the number of moles of each phase.
    ! dSumMolFractionSoln   A double real vector representing the sum of mole fractions in each solution phase.
    ! dPGEFunctionNorm      A double real scalar representing the norm of the functional vector in the PGESolver.
    ! dTolerance            A double real vector representing numerical tolerances defined in InitThermo.f90.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckSolnPhaseRem

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::   i, j, k, iMaxDrivingForce
    real(8)::   dMaxDrivingForce
    logical::   lPhasePass, lSwapLater
    
    
    ! Initialize variables:
    iMaxDrivingForce = 0
    dMaxDrivingForce = 0D0
    lPhasePass       = .FALSE.
    lSwapLater       = .FALSE.
    
    ! Loop through all solution phases currently predicted to be stable at equilibrium:
    LOOP_SolnRem: do i = 1, nSolnPhases
        
        ! Index of phase in dMolesPhase/iAssemblage vectors (not absolute index):
        j = nElements - i + 1   
        
        ! Check if the number of moles of this solution phase is less than a specified value:
        IF_RemSolnPhase: if (dMolesPhase(j) < dTolerance(7)) then
                                
            ! If the solution phase that is to be removed was the last phase that was swapped, then
            ! revert the system back to the last successful phase assemblage:
            if ((iSolnSwap == -iAssemblage(nElements - i + 1)).AND.(iterLast == iterSwap)) then
                
                ! Revert system:
                call RevertSystem(iterSwap)
                
                ! Reset iSolnSwap:
                iSolnSwap = 0
                
                ! Exit after the system has reverted:
                exit LOOP_SolnRem
    
            end if
                            
            ! Try removing this solution phase from the system:
            call RemSolnPhase(i,lPhasePass)  
                                                                                                            
            ! Exit if the phase assemblage has passed:
            if (lPhasePass .EQV. .TRUE.) exit LOOP_SolnRem
                          
            ! It may be possible that a pure condensed phase needs to be added to the system.  Try swapping 
            ! these two phases.
                          
            ! Compute the driving force for all pure condensed phases:
            call CompDrivingForce(iMaxDrivingForce,dMaxDrivingForce)

            ! Check if a pure condensed phase should be swapped:
            if ((iMaxDrivingForce /= 0).AND.(dMaxDrivingForce < dTolerance(4))) then
                ! The index of the pure condensed phase that is to be added is equal to iMaxDrivingForce if
                ! it is not zero.
                                                                                                            
                ! Try removing this solution phase and adding a pure condensed phase:
                call RemSolnAddPureConPhase(iMaxDrivingForce,j,lPhasePass)
                                                                                     
                ! If the new phase assemblage has passed, then exit:
                if (lPhasePass .EQV. .TRUE.) exit LOOP_SolnRem
                            
            end if

            ! It may be possible that a solution phase needs to be removed from the system, but it cannot be 
            ! removed directly for appropriate reasons, while simultaneously a different solution phase
            ! should be added to the system.  Try swapping these two phases:
                        
            ! Compute the "hypothetical" mole fractions:
            LOOP_SwapSolnPhaseSpecific: do k = 1, nSolnPhasesSys
                                                                                            
                ! Make sure that this phase doesn't try to replace itself with a miscible phase:
                if (cSolnPhaseName(k) == cSolnPhaseName(-iAssemblage(nElements-i+1))) cycle LOOP_SwapSolnPhaseSpecific
            
                ! Compute the mole fractions and the sum of mole fractions for solution phases that are currently
                ! not estimated to be stable at equilibrium:
                if (lSolnPhases(k) .EQV. .FALSE.) call CompMolFraction(k)
          
                ! Consider swapping these two phases if the driving force of phase k is less than tolerance:
                if ((lSolnPhases(k) .EQV. .FALSE.).AND.(dDrivingForceSoln(k) < dTolerance(4))) then
                    
                    ! Try swapping these two solution phases.
                    call SwapSolnPhaseSpecific(k,i,lPhasePass)
                                               
                    ! Swapping these two solution phases was successful.
                    if (lPhasePass .EQV. .TRUE.) exit LOOP_SolnRem
                
                end if
                                                
            end do LOOP_SwapSolnPhaseSpecific
            
            ! If a solution phase was recently swapped for another phase from the system, try swapping
            ! this phase with the phase that is trying to be removed:
            if ((iSolnPhaseLast /= 0).AND.(iterSwap == iterLast)) then
                
                ! Make sure that the solution phase corresponding to iSolnPhaseLast is not already predicted to
                ! be stable:
                lSwapLater = .FALSE.

                LOOP_B: do k = 1, nSolnPhases 
                    if (-iAssemblage(nElements - k + 1) == iSolnPhaseLast) then 
                        lSwapLater = .TRUE. 
                        exit LOOP_B
                    end if
                end do LOOP_B
                
                ! The phase corresponding to iSolnPhaseLast is not currently predicted to be stable.  Try 
                ! swapping that phase with the phase that is trying to be removed from the system.
                if (lSwapLater .EQV. .FALSE.) then
                
                    k = iSolnPhaseLast
                
                    ! Swap the two phases:
                    call SwapSolnPhaseSpecific(k, i, lPhasePass)
                                                            
                    ! Exit if the new phase assemblage has passed:
                    if (lPhasePass .EQV. .TRUE.) exit LOOP_SolnRem
                
                else
                    lSwapLater = .FALSE.
                end if
            
            end if
            
            ! If a pure condensed phase was recently swapped for another phase from the system, try
            ! reverting to that phase assemblage:
            if ((iterSwap == iterLast).AND.(iPureConSwap /= 0)) then
            
                ! Revert to the last phase assemblage:
                call RevertSystem(iterLast)
                
                ! Exit from loop:
                exit LOOP_SolnRem
            
            end if
            
            ! If a solution phase is extremely small and it cannot be directly removed or swapped with a pure 
            ! condensed phase or another solution phase, try removing a pure condensed phase (last resort):
            if ((dMolesPhase(j) < dTolerance(14)).AND.(nConPhases > 0)) then
            
                ! Loop through all pure condensed phases in the system
                do j = 1, nConPhases
                    
                    ! Try removing this pure condensed phase:
                    call RemPureConPhase(j,lSwapLater,lPhasePass)
                    
                    ! If the new phase assemblage can be removed, exit:
                    if (lPhasePass .EQV. .TRUE.) exit LOOP_SolnRem
                    
                end do
             
             end if
    
        end if IF_RemSolnPhase
        
    end do LOOP_SolnRem
    
    return
                
end subroutine CheckSolnPhaseRem