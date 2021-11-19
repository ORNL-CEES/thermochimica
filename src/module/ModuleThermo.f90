
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file        ModuleThermo.f90
    !> \brief       Fortran module for internal use of Thermochimica
    !> \details     The purpose of this module is to provide the means to share information amongst the
    !!              various subroutines used by Thermochimica
    !> \author      M.H.A. Piro
    !
    !> \param       nElements           The number of elements in the system. nElements >= nConPhases + nSolnPhases.
    !> \param       nElementsPT         The number of chemical elements on the periodic table.
    !> \param       nSpecies            The number of species in the system.
    !> \param       nSpeciesPhase       A cumulative integer vector representing the number of species in each solution phase.
    !!                                   0-indexed with nSpeciesPhase(0)=0 for determining # of species in phase i via
    !!                                   nSpeciesPhase(i) - nSpeciesPhase(i-1).
    !> \param       nParam              The number of mixing parameters in the system.
    !> \param       nParamPhase         The number of mixing parameters in each solution phase. Indexing follows same scheme as
    !!                                   nSpeciesPhase.
    !> \param       nMaxParam           The maximum number of parameters that are allowed.
    !> \param       nDummySpecies       The number of dummy species that have been added to a ChemSage data-file from FactSage.
    !> \param       nConPhases          The number of pure condensed phases predicted to be stable at equilibrium.
    !> \param       nConPhasesSys       The number of pure condensed phases in the system (not necessarily stable at equilibrium).
    !!                                   nConPhasesSys >= nConPhases.
    !> \param       nSolnPhases         The number of solution phases predicted to be stable at equilibrium.
    !> \param       nSolnPhasesSys      The number of solution phases in the system (not necessarily stable at equilibrium).
    !!                                   nSolnPhasesSys >= nSolnPhases.
    !> \param       iTolNum             The number of numerical tolerances (used only for allocation purposes).
    !> \param       dSpeciesTotalAtoms  The total number of atoms per formula mass of a species.
    !> \param       iPhase              An integer vector of length nSpecies representing the phase type (0: pure condensed phase;
    !!                                   > 0: solution phase index; -1: dummy species).
    !> \param       iParticlesPerMole   The number of particles per mole of constituent.
    !> \param       iAssemblage         Integer vector of length nElements due to Gibbs' phase rule which contains the indices of
    !!                                   phases estimated to be part of the equilibrium phase assemblage.
    !!                                   iAssemblage(1:nConPhases) contains species indices of condensed phases in the assemblage.
    !!                                   iAssemblage(nElements-nSolnPhases+1:nElements) contain negatives of indices of solution
    !!                                   phases in the assemblage. Will contain nElements-nConPhases-nSolnPhases 0 entries.
    !!                                   Example: nElements=4, nConPhases=1, nSolnPhases=2 might result in iAssemblage={8,0,-3,-2}
    !!                                   The name of the stable pure condensed phase is at cSpeciesName(8), and the names of the
    !!                                   stable solution phases are at cSolnPhaseName(3) and cSolnPhaseName(2).
    !> \param       iSpeciesPass        An integer vector that is used soley to determine whether a particular
    !!                                   species will be considered in the system.
    !> \param       dStoichSpecies      A matrix of size nSpecies by nElements representing the stoichiometry coefficients for
    !!                                   each element in each species.
    !> \param       iRegularParam       An integer matrix representing information pertient to regular solution models. The first
    !!                                   coefficient represents the number of components in the sub-system and the other
    !!                                   coefficients represent the indices of components in the sub-system.
    !> \param       iterHistoryLevel    An integer matrix representing all of the indices of phases that contribute to the
    !!                                   equilibrium phase assemblage at each stage in the iteration history during Leveling.
    !> \param       dIdealConstant      The ideal gas constant.
    !> \param       dChemicalPotential  A double real vector of length nSpecies representing the chemical potential of each
    !!                                   species. To be precise, this is defined as the difference between the standard molar
    !!                                   Gibbs energy and the chemical potential defined by the element potentials (represented
    !!                                   in dimensionless units and per formula mass).
    !> \param       dElementPotential   A double real vector of length nElements representing the element potentials.
    !> \param       dExcessGibbsParam   A double real vector of length nParam representing excess Gibbs energy of mixing
    !!                                   parameters.
    !> \param       dMolesElement       A double real vector of length nElements representing the total number of moles of each
    !!                                   element.
    !> \param       dMolesPhase         A double real vector of length nElements representing the moles of each phase in the
    !!                                   assemblage. The ordering of the phases matches iAssemblage.
    !> \param       dMolesSpecies       A double real vector of length nSpecies representing the number of moles each species in
    !!                                   the system.
    !> \param       dMolFraction        A double real vector of length nSpecies representing the mole fraction of each species in
    !!                                   the system.
    !> \param       dLevel              A double real vector of length nElements representing the adjustment applied to the
    !!                                   element potentials.
    !> \param       dAtomFractionSpecies  A double real matrix representing tha atom fraction of each element in each species.
    !> \param       dTolerance          A double real vector representing numerical tolerances (defined in InitThermo.f90).
    !> \param       cElementName        A character vector of length nElements representing the names of elements in the system.
    !> \param       cSpeciesName        A character vector of length nSpecies representing the name of each species in short-form.
    !> \param       cSolnPhaseName      A character vector of length nSolnPhasesSys representing the name of each solution phase.
    !> \param       cSolnPhaseType      A character vector of length nSolnPhasesSys representing the type of each solution phase.
    !> \param       nCountSublattice    An integer scalar representing the number of sublattice phases (e.g., Compound Energy
    !!                                   Formalism, aqueous, MQM, etc.).
    !> \param       iPhaseSublattice    An integer vector representing the sublattice # for each soluton phase.
    !> \param       dSiteFraction       A double real array with 3 dimensions representing the site fractions on each sublattice
    !!                                   for each ionic phase.  The first dimension corresponds to the ionic phase index, the
    !!                                   second dimension corresponds to the sublattice index and the third dimension corresponds
    !!                                   to the constituent index.
    !
    !-------------------------------------------------------------------------------------------------------------


module ModuleThermo

    implicit none

    SAVE

    integer::                                     nElements, nSpecies, nParam, nMaxParam, nDummySpecies, nElemOrComp, nMagParam
    integer::                                     nConPhases, nSolnPhases, nSolnPhasesSys, nChargedConstraints, nConPhasesSys
    integer::                                     nMaxSublatticeSys, nMaxConstituentSys, nCountSublattice
    integer,       parameter::                    iTolNum = 15, nElementsPT = 118, nMaxCompounds = 50
    integer,       dimension(:),   allocatable::  iPhase, nSpeciesPhase, iParticlesPerMole, iSUBIMixType
    integer,       dimension(:),   allocatable::  nParamPhase, iElementSystem, iSpeciesPass, nMagParamPhase
    integer,       dimension(:),   allocatable::  iAssemblage, iAssemblageBest
    integer,       dimension(:),   allocatable::  nSublatticePhase, iPhaseSublattice, iPhaseElectronID
    integer,       dimension(:,:), allocatable::  iRegularParam, iterHistoryLevel, nConstituentSublattice, nPairsSRO, iMagneticParam
    integer,       dimension(:,:), allocatable::  iSUBLParamData
    integer,       dimension(:,:,:),allocatable:: iConstituentPass, iConstituentSublattice, iPairID
    integer,       dimension(:,:,:),allocatable:: iChemicalGroup

    logical::                                     lHeatCapacityCurrent = .FALSE.

    real(8)::                                     dIdealConstant, dNormalizeSum, dNormalizeInput, dMassScale, dTemperatureForLimits
    real(8),       dimension(iTolNum)::           dTolerance
    real(8),       dimension(:),   allocatable::  dStdGibbsEnergy,    dGibbsSolnPhase, dMolesSpecies, dMagGibbsEnergy
    real(8),       dimension(:),   allocatable::  dChemicalPotential, dExcessGibbsParam, dLevel, dSpeciesTotalAtoms
    real(8),       dimension(:),   allocatable::  dElementPotential, dElementPotentialBest, dMolesElement, dAtomicMass
    real(8),       dimension(:),   allocatable::  dMolesPhase, dMolesPhaseBest, dMolFraction, dMolFractionBest
    real(8),       dimension(:,:), allocatable::  dAtomFractionSpecies, dStoichSublattice, dStoichSpecies, dQKTOParams
    real(8),       dimension(:,:), allocatable::  dCoeffGibbsMagnetic, dZetaSpecies, dMagneticParam

    real(8),      dimension(:,:,:),allocatable::  dSiteFraction, dCoordinationNumber, dSublatticeCharge, dStoichPairs
    real(8),      dimension(:,:,:),allocatable::  dConstituentCoefficients

    character(12), dimension(:),   allocatable::  cElementName
    character(30), dimension(:),   allocatable::  cSpeciesName
    character(8),  dimension(:),   allocatable::  cSolnPhaseType
    character(25), dimension(:),   allocatable::  cSolnPhaseName
    character(8),  dimension(:,:,:),allocatable:: cConstituentNameSUB
    character,     dimension(:),   allocatable::  cRegularParam
    character(30),  dimension(:,:), allocatable :: cPairName

end module ModuleThermo
