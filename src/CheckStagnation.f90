
    !---------------------------------------------------------------------------
    !
    !> \file    CheckStagnation.f90
    !> \brief   Check if the system is stagnant.
    !> \author  M.H.A. Piro
    !> \date    July 4, 2012
    !> \sa      CheckPhaseAssemblage.f90
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !   07/04/2012      M.H.A. Piro     Original code (happy Independence Day)
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check whether the system
    !! has become stagnant by determining the number of solution phases that
    !! have molar quantities that changed by a specified value.  For example,
    !! this subroutine will return the number of solution phases that have
    !! changed by more than 5%. 
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]  dMolesPhaseChange    A double real scalar representing the 
    !!                                   relative change of the number of moles
    !!                                   of a solution phase (e.g., = 0.05).
    !> \param[out] nPhasesCheck         An integer scalar representing the 
    !!                                   number of solution phases that have
    !!                                   molar quantities that have changed more
    !!                                   than the specified amount.
    !> \param[out] dMaxChange           A double real scalar representing the 
    !!                                   maximum relative change of dMolesPhase.
    !
    !---------------------------------------------------------------------------



subroutine CheckStagnation(dMolesPhaseChange,dMaxChange,nPhasesCheck)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none
    
    integer :: nPhasesCheck, j, k
    real(8) :: dMolesPhaseChange, dTemp, dMaxChange
    
    
    ! Initialize variables:
    nPhasesCheck = 0
    dMaxChange   = 0D0
    
    ! Check to make sure that dMolesPhaseChange is a reasonable value:
    if (dMolesPhaseChange > 0D0) then
    
        ! Loop through all solution phases:
        do j = 1, nSolnPhases
        
            ! Determine relative solution phase index:
            k = nElements - j + 1
        
            ! Compute the relative change of the number of moles of this solution phase:
            dTemp = DABS(dMolesPhase(k) - dMolesPhaseLast(k)) / dMolesPhaseLast(k)
        
            ! Count the number of phases that have significant changes to their molar quantities:
            if (dTemp >= dMolesPhaseChange) nPhasesCheck = nPhasesCheck + 1

            ! Compute maximum fractional change:
            dMaxChange = DMAX1(dTemp, dMaxChange)
                
        end do
        
    end if
    
    return
    
end subroutine CheckStagnation
