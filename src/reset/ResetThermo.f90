
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
    USE ModuleThermoIO, ONLY: INFOThermo, lRetryAttempted
    USE ModuleGEMSolver
    USE ModuleSubMin

    implicit none

    integer :: i, INFO

    ! Initialize variables:
    i = 0

    if (allocated(dStoichSpecies)) deallocate(dStoichSpecies, STAT = INFO)
    i = i + INFO
    if (allocated(dSpeciesTotalAtoms)) deallocate(dSpeciesTotalAtoms, STAT = INFO)
    i = i + INFO
    if (allocated(iPhase)) deallocate(iPhase, STAT = INFO)
    i = i + INFO
    if (allocated(nSpeciesPhase)) deallocate(nSpeciesPhase, STAT = INFO)
    i = i + INFO
    if (allocated(nParamPhase)) deallocate(nParamPhase, STAT = INFO)
    i = i + INFO
    if (allocated(iElementSystem)) deallocate(iElementSystem, STAT = INFO)
    i = i + INFO
    if (allocated(iSpeciesPass)) deallocate(iSpeciesPass, STAT = INFO)
    i = i + INFO
    if (allocated(iParticlesPerMole)) deallocate(iParticlesPerMole, STAT = INFO)
    i = i + INFO
    if (allocated(dMagGibbsEnergy)) deallocate(dMagGibbsEnergy, STAT = INFO)
    i = i + INFO
    if (allocated(dCoeffGibbsMagnetic)) deallocate(dCoeffGibbsMagnetic, STAT = INFO)
    i = i + INFO
    if (allocated(iAssemblage)) deallocate(iAssemblage, STAT = INFO)
    i = i + INFO
    if (allocated(dChemicalPotential)) deallocate(dChemicalPotential, STAT = INFO)
    i = i + INFO
    if (allocated(dMolesElement)) deallocate(dMolesElement, STAT = INFO)
    i = i + INFO
    if (allocated(dAtomFractionSpecies)) deallocate(dAtomFractionSpecies, STAT = INFO)
    i = i + INFO
    if (allocated(dStdGibbsEnergy)) deallocate(dStdGibbsEnergy, STAT = INFO)
    i = i + INFO
    if (allocated(nMagParamPhase)) deallocate(nMagParamPhase, STAT = INFO)
    i = i + INFO
    if (allocated(iMagneticParam)) deallocate(iMagneticParam, STAT = INFO)
    i = i + INFO
    if (allocated(dMagneticParam)) deallocate(dMagneticParam, STAT = INFO)
    i = i + INFO
    if (allocated(dMolesSpecies)) deallocate(dMolesSpecies, STAT = INFO)
    i = i + INFO
    if (allocated(dElementPotential)) deallocate(dElementPotential, STAT = INFO)
    i = i + INFO
    if (allocated(dMolesPhase)) deallocate(dMolesPhase, STAT = INFO)
    i = i + INFO
    if (allocated(dMolesPhaseLast)) deallocate(dMolesPhaseLast, STAT = INFO)
    i = i + INFO
    if (allocated(cElementName)) deallocate(cElementName, STAT = INFO)
    i = i + INFO
    if (allocated(cSpeciesName)) deallocate(cSpeciesName, STAT = INFO)
    i = i + INFO
    if (allocated(cSolnPhaseType)) deallocate(cSolnPhaseType, STAT = INFO)
    i = i + INFO
    if (allocated(cSolnPhaseName)) deallocate(cSolnPhaseName, STAT = INFO)
    i = i + INFO
    if (allocated(dAtomicMass)) deallocate(dAtomicMass, STAT = INFO)
    i = i + INFO
    if (allocated(iterHistoryLevel)) deallocate(iterHistoryLevel, STAT = INFO)
    i = i + INFO
    if (allocated(dLevel)) deallocate(dLevel, STAT = INFO)
    i = i + INFO
    if (allocated(iterHistory)) deallocate(iterHistory, STAT = INFO)
    i = i + INFO
    if (allocated(dSumMolFractionSoln)) deallocate(dSumMolFractionSoln, STAT = INFO)
    i = i + INFO
    if (allocated(dMolFraction)) deallocate(dMolFraction, STAT = INFO)
    i = i + INFO
    if (allocated(dUpdateVar)) deallocate(dUpdateVar, STAT = INFO)
    i = i + INFO
    if (allocated(dPartialExcessGibbs)) deallocate(dPartialExcessGibbs, STAT = INFO)
    i = i + INFO
    if (allocated(dPartialExcessGibbsLast)) deallocate(dPartialExcessGibbsLast, STAT = INFO)
    i = i + INFO
    if (allocated(dEffStoichSolnPhase)) deallocate(dEffStoichSolnPhase, STAT = INFO)
    i = i + INFO
    if (allocated(lSolnPhases)) deallocate(lSolnPhases, STAT = INFO)
    i = i + INFO
    if (allocated(dGibbsSolnPhase)) deallocate(dGibbsSolnPhase, STAT = INFO)
    i = i + INFO
    if (allocated(lMiscibility)) deallocate(lMiscibility, STAT = INFO)
    i = i + INFO
    if (allocated(dDrivingForceSoln)) deallocate(dDrivingForceSoln, STAT = INFO)
    i = i + INFO
    if (allocated(lSolnPhases)) deallocate(lSolnPhases, STAT = INFO)
    i = i + INFO
    if (allocated(dGibbsSolnPhase)) deallocate(dGibbsSolnPhase, STAT = INFO)
    i = i + INFO
    if (allocated(lMiscibility)) deallocate(lMiscibility, STAT = INFO)
    i = i + INFO
    if (allocated(iConstituentPass)) deallocate(iConstituentPass, STAT = INFO)
    i = i + INFO
    if (allocated(iPhaseSublattice)) deallocate(iPhaseSublattice, STAT = INFO)
    i = i + INFO
    if (allocated(dStoichSublattice)) deallocate(dStoichSublattice, STAT = INFO)
    i = i + INFO
    if (allocated(dSiteFraction)) deallocate(dSiteFraction, STAT = INFO)
    i = i + INFO
    if (allocated(cConstituentNameSUB)) deallocate(cConstituentNameSUB, STAT = INFO)
    i = i + INFO
    if (allocated(iConstituentSublattice)) deallocate(iConstituentSublattice, STAT = INFO)
    i = i + INFO
    if (allocated(iPhaseElectronID)) deallocate(iPhaseElectronID, STAT = INFO)
    i = i + INFO
    if (allocated(nSublatticePhase)) deallocate(nSublatticePhase, STAT = INFO)
    i = i + INFO
    if (allocated(nConstituentSublattice)) deallocate(nConstituentSublattice, STAT = INFO)
    i = i + INFO
    if (allocated(dZetaSpecies)) deallocate(dZetaSpecies, STAT = INFO)
    i = i + INFO
    if (allocated(dConstituentCoefficients)) deallocate(dConstituentCoefficients, STAT = INFO)
    i = i + INFO
    if (allocated(dSublatticeCharge)) deallocate(dSublatticeCharge, STAT = INFO)
    i = i + INFO
    if (allocated(iChemicalGroup)) deallocate(iChemicalGroup, STAT = INFO)
    i = i + INFO
    if (allocated(cPairName)) deallocate(cPairName, STAT = INFO)
    i = i + INFO
    if (allocated(dStoichPairs)) deallocate(dStoichPairs, STAT = INFO)
    i = i + INFO
    if (allocated(iHessian)) deallocate(iHessian, STAT = INFO)
    i = i + INFO
    if (allocated(dChemicalPotentialStar)) deallocate(dChemicalPotentialStar, STAT = INFO)
    i = i + INFO
    if (allocated(dRHS)) deallocate(dRHS, STAT = INFO)
    i = i + INFO
    if (allocated(dHessian)) deallocate(dHessian, STAT = INFO)
    i = i + INFO
    if (allocated(iPairID)) deallocate(iPairID, STAT = INFO)
    i = i + INFO
    if (allocated(dCoordinationNumber)) deallocate(dCoordinationNumber, STAT = INFO)
    i = i + INFO
    if (allocated(nPairsSRO)) deallocate(nPairsSRO, STAT = INFO)
    i = i + INFO
    if (allocated(iRegularParam)) deallocate(iRegularParam, STAT = INFO)
    i = i + INFO
    if (allocated(dExcessGibbsParam)) deallocate(dExcessGibbsParam, STAT = INFO)
    i = i + INFO
    if (allocated(cRegularParam)) deallocate(cRegularParam, STAT = INFO)
    i = i + INFO
    if (allocated(iSUBLParamData)) deallocate(iSUBLParamData, STAT = INFO)
    i = i + INFO
    if (allocated(iSUBIMixType)) deallocate(iSUBIMixType, STAT = INFO)
    i = i + INFO
    if (allocated(dQKTOParams)) deallocate(dQKTOParams, STAT = INFO)
    i = i + INFO

    lRetryAttempted = .FALSE.

    ! Return an INFOThermo if deallocation of any of the allocatable variables failed:
    if (i > 0) then
        INFOThermo = 15
    end if

    return

end subroutine ResetThermo
