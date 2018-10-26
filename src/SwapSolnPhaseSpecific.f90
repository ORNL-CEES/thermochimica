
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    SwapSolnPhaseSpecific.f90
    !> \brief   Swap one specific solution phase for another specific solution phase.
    !> \author  M.H.A. Piro
    !> \date    May 8, 2012
    !> \sa      CheckSolnPhaseRem.f90
    !> \sa      SwapSolnPhase.f90
    !> \sa      CompMolSolnPhase.f90
    !> \sa      CompStoichSolnPhase.f90
    !> \sa      CompChemicalPotential.f90
    !> \sa      CheckPhaseChange.f90
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   05/08/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to attempt to swap a particular solution phase for another
    !! specific solution phase in the current estimated phase assemblage.  The logical variable lPhasePass 
    !! returns a value of FALSE if the phase in question cannot swap for any other solution phase in the system.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iPhaseAdd       An integer scalar representing the absolute index of a solution phase that is 
    !!                               to be added to the system.
    !> \param[in]   iPhaseRem       An integer scalar representing the index of a solutoin phase that is 
    !!                               to be removed from the system.
    !> \param[out]  lPhasePass          A logical variable indicatin whether the phase assemblage is appropriate 
    !!                                   for consdieration (TRUE) or not (FALSE).
    !
    ! INFO                  An integer scalar indicating a successful exit (zero) or an error from the call
    !                        to CheckPhaseChange.f90. 
    ! nElements             An integer scalar representing the number of elements in the system.
    ! nSolnPhases           An integer scalar representing the number of solution phases in the assemblage.
    ! iAssemblage           Integer vector containing the indices of all phases currently estimated to contribute
    !                       to the equilibrium phase assemblage.  
    ! iPhaseChange          
    ! iSolnPhaseLast        Index of the last solution phase to be added to, withdrawn from, or exchanged, from
    !                       the estimated phase assemblage.
    ! dMolesPhase           A double real vector representing the number of moles of each phase.
    ! dTemp                 A double real scalar temporarily used for work purposes.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine SwapSolnPhaseSpecific(iPhaseAdd,iPhaseRem,lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                       :: i, j, iPhaseAdd, iPhaseRem, INFO, iterBack
    integer, dimension(nElements) :: iAssemblageTest
    real(8)                       :: dTemp
    logical                       :: lPhasePass, lCompEverything, lSwapLater


    ! Initialize variables:
    INFO            = 0
    iterBack        = 200
    lCompEverything = .FALSE.
    lPhasePass      = .FALSE.
    lSwapLater      = .FALSE.
    j               = nElements - iPhaseRem + 1  ! Relative index of phase in iAssemblage/dMolesPhase vectors
    
    ! Check the iteration history:
    iAssemblageTest    = iAssemblage 
    iAssemblageTest(j) = -iPhaseAdd
    
    ! Check whether this particular phase assemblage has been previously considered:
    call CheckIterHistory(iAssemblageTest,iterBack,lSwapLater)

    ! Return if this phase assemblage fails:
    if (lSwapLater .EQV. .TRUE.) return
    
    ! Compute the stoichiometry of the solution phase that is to be added to the system:
    call CompStoichSolnPhase(iPhaseAdd)
                      
    iSolnPhaseLast = -iAssemblage(j)            ! Absolute index of phase to be removed.
    dTemp          = dMolesPhase(j)
    
    ! Swap solution phases:
    iAssemblage(j) = -iPhaseAdd 
    dMolesPhase(j) = 0D0
    
    ! Compute the number of moles of the solution phase that has been added to the system:
    call CompMolSolnPhase    
    
    ! Compute the number of moles of solution species:
    do i = nSpeciesPhase(iPhaseAdd-1) + 1, nSpeciesPhase(iPhaseAdd)
        dMolesSpecies(i) = dMolFraction(i) * dMolesPhase(j)
    end do    
      
    ! Compute the chemical potentials:
    call CompChemicalPotential(lCompEverything)    
            
    ! Check that this phase change is acceptable:                    
    call CheckPhaseChange(lPhasePass,INFO)
        
    if (lPhasePass .EQV. .TRUE.) then
        ! This phase assemblage can be considered.
        iterLastSoln                 = iterGlobal
        iterLast                     = iterGlobal
        iterSwap                     = iterGlobal
        dMolesPhaseLast(j)           = 0D0
        dDrivingForceSoln(iPhaseAdd) = 0D0
        lSolnPhases(iPhaseAdd)       = .TRUE.
        lSolnPhases(iSolnPhaseLast)  = .FALSE.
    else
        ! This phase assemblage cannot be considered.  Revert back to the previous values:
        dMolesPhase(j)               = dTemp
        iAssemblage(j)               = -iSolnPhaseLast
    end if
         
    return
           
end subroutine SwapSolnPhaseSpecific