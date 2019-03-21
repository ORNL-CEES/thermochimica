
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ResetThermo.f90
    !> \brief   Deallocate allocatable variables used by the ModuleThermo.f90, ModulePGESolver.f90 modules.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      ModuleThermo.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer         Description of change
    !    ----          ----------         ---------------------
    !    11/04/2011    M.H.A. Piro        Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to attempt to gracefully exit Thermochimica.  Allocatable
    !! arrays are deallocated and memory is stored for output to external packages.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO                  An error is returned if deallocation is unsuccessful.
    ! INFOThermo            An integer scalar identifying whether the program exits successfully or if
    !                       it encounters an error.  A description for each error is given in ThermoDebug.f90.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ResetThermo

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleGEMSolver
    USE ModuleSubMin

    implicit none

    integer::   i, INFO


    ! Initialize variables:
    i = 0

    if (allocated(dStoichSpecies)) then
        ! Deallocate integer arrays from ModuleThermo:
        deallocate (dStoichSpecies,dSpeciesTotalAtoms,iPhase,nSpeciesPhase,iAssemblage,nParamPhase, &
            iRegularParam,iElementSystem,iSpeciesPass,iParticlesPerMole, dMagGibbsEnergy, &
            dCoeffGibbsMagnetic, STAT = INFO)
        i = i + INFO
    end if

    if (allocated(dChemicalPotential)) then
        ! Deallocate real arrays from ModuleThermo:
        deallocate (dChemicalPotential,dExcessGibbsParam,dMolesSpecies,dElementPotential,dMolesPhase,&
            dMolesPhaseLast,dMolesElement,dAtomFractionSpecies,dStdGibbsEnergy,STAT = INFO)
        i = i + INFO
    end if

    if (allocated(cElementName)) then
        ! Deallocate character arrays from ModuleThermo:
        deallocate (cElementName,cSpeciesName,cSolnPhaseType,cSolnPhaseName, STAT = INFO)
        i = i + INFO
    end if

    if (allocated(iterHistoryLevel)) then
        ! Deallocate variables from ModuleThermo specific to Leveling:
        deallocate (iterHistoryLevel,dLevel, STAT = INFO)
        i = i + INFO
    end if

    if (allocated(iterHistory)) then
        ! Deallocate integer arrays from ModulePGESolver:
        deallocate (iterHistory,STAT = INFO)
        i = i + INFO
    end if

    if (allocated(dSumMolFractionSoln)) then
        ! Deallocate real arrays from ModulePGESolver:
        deallocate (dSumMolFractionSoln,dMolFraction,dUpdateVar,dPartialExcessGibbs, &
            dPartialExcessGibbsLast,dEffStoichSolnPhase,lSolnPhases,dGibbsSolnPhase, &
            lMiscibility,dDrivingForceSoln, STAT = INFO)
        i = i + INFO
    end if

    if (allocated(nConstituentSublattice)) then
        ! Deallocate arrays used for sublattices:
        deallocate(iConstituentPass, iPhaseSublattice, dStoichSublattice, dSiteFraction, &
            cConstituentNameSUB, iConstituentSublattice, iPhaseElectronID,nSublatticePhase,&
            nConstituentSublattice,  STAT = INFO)
        i = i + INFO
    end if

    if (allocated(dHessian)) then
        ! Deallocate allocatable arrays used by ModuleSubMin
        deallocate(iHessian, dChemicalPotentialStar, dRHS, dHessian, STAT = INFO)
        i = i + INFO
    end if

    if (allocated(iPairID)) then
        ! Deallocate pair IDs
        i = i + INFO
        deallocate(iPairID,dCoordinationNumber)
    end if

   if (allocated(nPairsSRO)) then
         deallocate(nPairsSRO, STAT = INFO)
      i = i + INFO
   end if 
      

!    if (allocated(lSpeciesStable)) then
!        deallocate(lSpeciesStable, STAT = INFO)
!        i = i + INFO
!    end if

    ! Return an INFOThermo if deallocation of any of the allocatable variables failed:
    if (i > 0) then
        INFOThermo = 15
    end if

    return

end subroutine ResetThermo
