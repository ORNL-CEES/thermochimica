
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ModuleParseCS.f90
    !> \brief   A Fortran module used to store data for the ChemSage parser.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !   04/24/2012      M.H.A Piro      Original code
    !   03/25/2019      P. Bajpai       Added SUBQ to supported models
    !
    !> \param INFO                      A scalar integer that indicates a successful exit or identifies
    !!                                    an error.
    !> \param nElementsCS               Number of elements in the system.
    !> \param nSolnPhasesSysCS          Number of solution phases in the system.
    !> \param nSolnPhasesSysMax         Maximum number of solution phases in a system that can be considered.
    !> \param nSpeciesPhaseCS           Number of species in a solution phase.
    !> \param nSpeciesCS                Number of species in the system (combined solution species and pure
    !!                                    separate phases).
    !> \param nGibbsEqSpecies           Number of Gibbs energy equations for a particular species.
    !> \param nGibbsCoeff               Number of coefficients for a Gibbs energy equation.
    !> \param nMaxGibbsEqs              Maximum number of Gibbs energy equations per species.
    !> \param nParamMax                 Maximum number of parameters in a sub-system of a non-ideal solution
    !!                                    phase.  The default is set to 4 (i.e., quaternary).
    !> \param dStoichSpeciesCS          A double real matrix representing the number of atoms of a particular
    !!                                    elements in a species (i.e., stoichiometry matrix).
    !> \param iRegularParam             Integer matrix representing the component numbers and exponents in
    !!                                    a regular solution phase.
    !> \param iPhase                    Integer vector containing the index of the phase that a particular
    !!                                    species belongs to.  iPhase = 0 for a pure separate phase, iPhase = -1
    !!                                    for a "dummy species", iPhase > 0 for a solution species, where the
    !!                                    number corresponds to the solution phase index.
    !> \param iParticlesPerMoleCS       An integer vector containing the number of particles per mole of the
    !!                                    constituent species formula mass.  The default value is 1.
    !> \param cSystemTitle              A character string representing the name of the system.
    !> \param cDummy                    A dummy character variable.
    !> \param cElementName              The name of a chemical element.
    !> \param cSpeciesNameCS            The name of a species (short hand).  Note that this can be a solution
    !!                                    species or pure condensed phase.
    !> \param cSolnPhaseName            The name of a solution phase.
    !> \param cSolnPhaseType            The type of a solution phase.
    !> \param dGibbsDummy               An abritrary value for the molar standard Gibbs energy that is applied
    !                                    to dummy species.
    !> \param dAtomicMass               Atomic mass of an element.
    !> \param dGibbsCoeffSpeciesTemp    Temporary double array of coefficients for a Gibbs energy equation.
    !> \param nCountSublatticeCS        An integer scalar representing the number of phases with sublattices
    !!                                  (e.g., SUBL, SUBG, SUBQ).
    !> \param nSolnTypeSupport          An integer scalar representing the number of supported solution phase types.
    !> \param cSolnPhaseTypeSupport     A character array representing the solution phase types that are supported.
    !> \param nSublatticePhaseCS        An integer vector representing the number of sublattices per phase.  This
    !!                                   only applies to charged phases.
    !> \param nConstituentSublattice    An integer matrix representing the number of constituents per sublattice
    !!                                   for each CEF phase.
    !> \param dStoichSublattice         A double real matrix representing the stoichiometry of each sublattice
    !!                                   for each CEF Phase.
    !> \param iConstituentSublatticeCS  An integer matrix representing the coefficients of constituents on each
    !!                                   sublattice. The first dimension corresponds to the phase, the second
    !!                                   corresponds to the sublattice and the third corresponds to the constituent.
    !> \param iPhaseSublatticeCS        An integer vector representing the sublattice ID for each solution phase.
    !!                                   This value is equal to zero for a solution phase that does not contain
    !!                                   any sublattices.
    !
    !-------------------------------------------------------------------------------------------------------------


module ModuleParseCS

    implicit none

    SAVE

    integer                                     :: nElementsCS, nSpeciesCS, nSolnPhasesSysCS, INFO, iMiscSUBI
    integer                                     :: nParamCS, nCountSublatticeCS, nMaxSpeciesPhaseCS, nMagParamCS
    integer,        parameter                   :: nSolnPhasesSysMax = 100, nMaxSublatticeCS = 5
    integer,        parameter                   :: nSolnTypeSupport = 9
    integer,        parameter                   :: nGibbsCoeff = 13, nMaxGibbsEqs = 6, nParamMax = 4
    integer,        dimension(:),   allocatable :: nSpeciesPhaseCS, nGibbsEqSpecies, iPhaseCS, iParticlesPerMoleCS
    integer,        dimension(:),   allocatable :: nParamPhaseCS, iParamPassCS, nSublatticePhaseCS, iPhaseSublatticeCS
    integer,        dimension(:),   allocatable :: iMagParamPassCS, nMagParamPhaseCS, iSUBIMixTypeCS
    integer,        dimension(:,:), allocatable :: iRegularParamCS, nConstituentSublatticeCS, nPairsSROCS, iMagneticParamCS
    integer,        dimension(:,:), allocatable :: iSUBIParamDataCS
    integer,        dimension(:,:,:), allocatable :: iConstituentSublatticeCS, iPairIDCS, iChemicalGroupCS

    real(8),        dimension(:),   allocatable :: dAtomicMassCS
    real(8),        dimension(:,:), allocatable :: dGibbsCoeffSpeciesTemp, dRegularParamCS, dGibbsMagneticCS, dMagneticParamCS
    real(8),        dimension(:,:), allocatable :: dStoichSublatticeCS, dStoichSpeciesCS, dZetaSpeciesCS, dStoichConstituentCS
    real(8),        dimension(:,:,:), allocatable :: dSublatticeChargeCS, dStoichPairsCS, dConstituentCoefficientsCS
    real(8),        dimension(:,:,:), allocatable :: dCoordinationNumberCS

    character(3),   dimension(:),   allocatable :: cElementNameCS
    character(8),   dimension(:),   allocatable :: cSolnPhaseTypeCS
    character(25),  dimension(:),   allocatable :: cSolnPhaseNameCS
    character(30),  dimension(:),   allocatable :: cSpeciesNameCS
    character(30),  dimension(:,:), allocatable :: cPairNameCS
    character(8),   dimension(:,:,:), allocatable :: cConstituentNameSUBCS
    character,      dimension(:),   allocatable :: cRegularParamCS

    character(8),   dimension(nSolnTypeSupport), parameter :: cSolnPhaseTypeSupport = &
                                                    ['IDMX    ','QKTO    ','SUBL    ','RKMP    ','RKMPM   ','SUBLM   ','SUBG    ', &
                                                    'SUBQ    ','SUBI    ']

end module ModuleParseCS
