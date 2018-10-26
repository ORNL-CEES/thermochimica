
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    PostLevelingSolver.f90
    !> \brief   Improve initial estimates from LevelingSolver.f90
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      Thermochimica.f90
    !> \sa      LevelingSolver.f90
    !> \sa      GetNewAssemblage.f90
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   05/02/2011      M.H.A. Piro         Original code
    !   07/31/2011      M.H.A. Piro         Clean up code: remove unnecessary variables, update variable names
    !   10/25/2011      M.H.A. Piro         Clean up code: modules, remove unnecessary variables.
    !   05/25/2012      M.H.A. Piro         Fix bug: Relocate update to the element potentials to after the 
    !                                        call to LAPACK.
    !   
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to improve the initial estimates provided by the Leveling.f90 
    !! subroutine.  The premise of the Leveling algorithm is to temporarily treat all species and phases as pure 
    !! stoichiometric phases, which is mathematically equivalent to assigning a unit activity to all species and 
    !! phases considered in the assemblage.  Similar to the Leveling subroutine, the PostLeveling subroutine only 
    !! considers nElements species in the system; however, the activity of those species is permitted to depart 
    !! from unity for solution phases.  The activity is taken to be equal to the mole fraction and is computed 
    !! using the estimated number of moles of that species.  In some situations, in particular when the number of 
    !! moles of the elements varies by many orders of magnitude, the PostLeveling subroutine can provide 
    !! significantly better initial estimates for optimization. 
    !!
    !! For more details, refer to the following literature:
    !!
    !!      M.H.A. Piro and S. Simunovic, "Performance Enhancing Algorithms for Computing Thermodynamic
    !!      Equilibria," CALPHAD, 39 (2012) 104-110.
    !
    !
    ! Pertinent variables:
    ! ====================
    ! 
    ! iAssemblage               An integer vector containing the indices of phases at equilibrium
    ! dMolesPhase               A double real vector representing the number of moles of a phase at 
    !                            equilibrium
    ! dChemicalPotential        A double real vector represening the chemical potential of each species.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine PostLevelingSolver
  
    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo
  
    implicit none
       
    integer::                                 i, j, INFO,  iterpl
    integer,dimension(nElements)::            IPIV, iAssemblageOld
    real(8),dimension(0:nSolnPhasesSys)::     dTempVec
    real(8),dimension(nElements)::            dTempVecE, dMolesPhaseOld, dElementPotentialOld
    real(8),dimension(nSpecies)::             dChemicalPotentialOld
    real(8),dimension(nElements,nElements)::  A
    
    
    ! Initialize variables:
    iterHistoryLevel = 0
    
    ! Store previous data in case if PostLeveling fails:
    iAssemblageOld        = iAssemblage
    dMolesPhaseOld        = dMolesPhase
    dChemicalPotentialOld = dChemicalPotential 
    dElementPotentialOld  = dElementPotential

    LOOP_POSTLEVEL: do iterpl = 1, nElements
                
        dTempVec  = 0D0
        dTempVecE = 1D0
        dLevel    = 0D0   
                
        do i = 1, nElements
            if (iPhase(iAssemblage(i)) < 0) exit LOOP_POSTLEVEL
        end do

        ! Calculate the total number of moles for each solution phase:
        do i = 1, nElements
            dTempVec(iPhase(iAssemblage(i))) = dTempVec(iPhase(iAssemblage(i))) + dMolesPhase(i)
        end do
        
        ! Compute the estimated mole fractions of species in solution phases 
        ! Note that pure condensed phases are set to 1.
        do i = 1,nElements
            if (iPhase(iAssemblage(i)) /= 0) then
                if (dMolesPhase(i) <= 0D0) exit LOOP_POSTLEVEL
                dTempVecE(i) = dMolesPhase(i) / dTempVec(iPhase(iAssemblage(i)))
            end if
        end do

        ! Establish matrix A and vector B for DGESV:
        do j = 1,nElements
            do i = 1,nElements
                A(i,j) = dAtomFractionSpecies(iAssemblage(i),j)
            end do
            dLevel(j) = dChemicalPotential(iAssemblage(j))  + DLOG(dTempVecE(j)) / & 
                !DFLOAT(iSpeciesTotalAtoms(iAssemblage(j)))
                dSpeciesTotalAtoms(iAssemblage(j))
        end do
                
        ! Call the linear equation solver to solve the adjustments applied to the Gibbs Plane:
        call DGESV( nElements, 1, A, nElements, IPIV, dLevel, nElements, INFO ) 
         
        ! Update the element potentials:
        dElementPotential = dElementPotential + dLevel

        ! Verify that the co-ordinates are real, otherwise return to the main program:       
        if (INFO /= 0) exit LOOP_POSTLEVEL 
        
        ! Readjust the Relative Gibbs Energies (dRelGibbsEnergy) of all species and phases in the system:
        dChemicalPotential = dChemicalPotential - MATMUL(dAtomFractionSpecies,dLevel)
    
        i = MAXVAL(MINLOC(dChemicalPotential))
        
        if ((dChemicalPotential(i) > dTolerance(4)).OR.(iterpl == nElements)) exit LOOP_POSTLEVEL
        
        ! Establish a new phase assemblage:
        call GetNewAssemblage(iterpl)
    
    end do LOOP_POSTLEVEL
    
    
    if (MINVAL(dChemicalPotential) < dTolerance(4)) then
        ! PostLevelling was unable to improve estimates from Leveling.  Return the estimates to their
        ! previous state.
        iAssemblage        = iAssemblageOld
        dMolesPhase        = dMolesPhaseOld
        dChemicalPotential = dChemicalPotentialOld
        dElementPotential  = dElementPotentialOld
    end if
    
    ! Verify that none of the predicted stable phases are dummy phases:
    do i = 1, nElements
        ! Cycle if it is an electron...this is the only exception.
        if (cSpeciesName(iAssemblage(i)) == 'e-') cycle
        if (iPhase(iAssemblage(i)) < 0) then
            INFOThermo = 9
            return
        end if
    end do
    
    return
    
end subroutine PostLevelingSolver