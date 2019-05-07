
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSHeader.f90
    !> \brief   Parse the header section of a ChemSage data-file.
    !> \author  M.H.A. Piro
    !> \date    Oct. 25, 2018
    !> \sa      ParseCSDataFile.f90
    !> \todo    Do a better job allocating variables iPairIDCS and dCoordinationNumberCS.
    !
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified and is completely
    ! independent of ChemApp and related products, including Solgas, Solgasmix, Fact, FactSage
    ! and ChemSage.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   10/06/2011      M.H.A. Piro         Original code
    !   02/01/2012      M.H.A. Piro         Correction: ChemSage files always include the gas phase in
    !                                        indexing (always treated as the first phase).  If the gas
    !                                        phase is not considered in the system, the indexing scheme
    !                                        assigned zero constituents to this phase, but indexing must
    !                                        be consistent.
    !   02/07/2012      M.H.A. Piro         Correction: read # of species when there aren't any solution
    !                                        phases.
    !   12/20/2012      M.H.A. Piro         Added capability of interpreting electrons as system components.
    !   01/15/2013      M.H.A. Piro         Added variables specific to ionic phases.
    !   10/25/2018      M.H.A. Piro         Fixed bug: nGibbsEqSpecies was not initialized, resulting in
    !                                        memory leaks.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to parse the "header section" of a ChemSage data-file.
    !! The "Header Section" of this datafile contains important information for allocating
    !! memory, such as the number of elements in the system, number of solution phases in the
    !! system, the number of chemical species etc.  Details about the INFO error codes are given
    !! in ParseCSDataFile.f90.  Details about line numbers are explained in the "ChemApp
    !! Programmer's manual".
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO                      A scalar integer that indicates a successful exit or identifies an error.
    ! nElementsPT               Number of elements in the periodic table.
    ! nElements                 Number of elements in the system.
    ! nCountSublatticeCS        Number of sublattice phases in the system.
    ! nSolnPhasesSys            Number of solution phases in the system.
    ! nSolnPhasesSysMax         Maximum number of solution phases in a system that can be considered.
    ! nSpeciesPhase             Number of species in a solution phase.
    ! nSpecies                  Number of species in the system (combined solution species and pure separate
    !                           phases).
    ! nGibbsEqSpecies           Number of Gibbs energy equations for a particular species.
    ! nGibbsCoeff               Number of coefficients for a Gibbs energy equation.
    ! nMaxGibbsEqs              Maximum number of Gibbs energy equations per species.
    ! iSpeciesAtoms             Integer matrix representing the number of atoms of a particular elements in a
    !                           species.
    ! cSystemTitle              A character string representing the name of the system.
    ! cDummy                    A dummy character variable.
    ! cElementName              The name of a chemical element.
    ! cSolnPhaseName            The name of a solution phase.
    ! cSpeciesName              The name of a species (short hand).
    ! dAtomicMass               Atomic mass of an element.
    ! dGibbsCoeffSpeciesTemp    Temporary double array of coefficients for a Gibbs energy equation.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ParseCSHeader

    USE ModuleParseCS

    implicit none

    integer              :: i, j, k, iGasPhase, nMaxSpeciesPhaseCS
    integer,dimension(7) :: iGequation
    integer,dimension(2) :: iDummy
    character(8)         :: cDummy
    character(78)        :: cSystemTitle

    ! Line 1: Read title of system:
    read (1,*,IOSTAT = INFO) cDummy, cSystemTitle
    ! Note: the title of the system is useless and does not need to be stored. Also, ChemSage files crop the
    ! title after a certain number of characters making it useless if you intend to read the string of elements.

    if (INFO /= 0) then
        INFO = 101
        return
    end if

    ! Line 2: Read the number of elements in the system and the number of solution phases:
    read (1,*,IOSTAT = INFO) nElementsCS, nSolnPhasesSysCS, iGasPhase

    if (INFO /= 0) then
        ! Report an error if there is an issue parsing this line.
        INFO = 102
        return
    elseif (nSolnPhasesSysCS > nSolnPhasesSysMax) then
        ! Ensure that the maximum number of solution phases has not exceeded the maximum allowed:
        INFO = 8
        return
    end if

    ! Check to see if there are any gaseous species:
    if (iGasPhase == 0) then
        ! There is no gas phase, the number of solution phases and indexing must be adjusted.
        nSolnPhasesSysCS = nSolnPhasesSysCS - 1
        iGasPhase = 1
    else
        ! There is a gas phase.
        iGasPhase = 0
    end if

    ! Allocate memory for pertinent variables:
    j = MAX(1,nSolnPhasesSysCS)
    allocate(nSpeciesPhaseCS(0:j),nParamPhaseCS(0:j))
    allocate(cElementNameCS(nElementsCS),dAtomicMass(nElementsCS))
    allocate(cSolnPhaseNameCS(j),cSolnPhaseTypeCS(j))
    allocate(iPhaseSublatticeCS(j))
    allocate(nPairsSROCS(j,2))

    ! TEMPORARY:

    ! Initialize variables:
    nCountSublatticeCS    = 0
    nSpeciesPhaseCS       = 0
    iPhaseSublatticeCS    = 0
    nPairsSROCS           = 0

    ! Go back to Line 2:
    backspace(UNIT = 1)

    ! Continue reading the rest of Line 2:
    if (nSolnPhasesSysCS == 0) then
        ! There aren't any solution phases in the system.
        read (1,*,IOSTAT = INFO) nElementsCS, iDummy(1:iGasPhase+1), nSpeciesCS
    else
        read (1,*,IOSTAT = INFO) nElementsCS, iDummy(1:iGasPhase+1), nSpeciesPhaseCS(1:nSolnPhasesSysCS), nSpeciesCS
    end if

    if ((nSolnPhasesSysCS == 1).AND.(nSpeciesPhaseCS(1) == 0)) then
        ! There are no solution phases in the system
        nSolnPhasesSysCS = 0
    end if

    if (INFO /= 0) then
        INFO = 102
        return
    end if

    ! Store the maximum number of species per solution phase in the system:
    nMaxSpeciesPhaseCS = MAXVAL(nSpeciesPhaseCS)

    ! Compute the number of species per solution phase:
    do i = 1, nSolnPhasesSysCS
        nSpeciesPhaseCS(i) = nSpeciesPhaseCS(i) + nSpeciesPhaseCS(i-1)
    end do

    ! Compute the total number of species in the system:
    nSpeciesCS = nSpeciesCS + nSpeciesPhaseCS(nSolnPhasesSysCS)

    allocate(iPairIDCS(nSolnPhasesSysCS,nMaxSpeciesPhaseCS,4))
    allocate(dCoordinationNumberCS(nSolnPhasesSysCS,nMaxSpeciesPhaseCS,4))
    iPairIDCS             = 0
    dCoordinationNumberCS = 0D0

    ! Allocate allocatable arrays:
    allocate(dStoichSpeciesCS(nSpeciesCS,nElementsCS))
    allocate(cSpeciesNameCS(nSpeciesCS),nGibbsEqSpecies(nSpeciesCS))
    allocate(dGibbsCoeffSpeciesTemp(nGibbsCoeff,nSpeciesCS*nMaxGibbsEqs))
    allocate(iRegularParamCS(1000,nParamMax*2+1),dRegularParamCS(1000,6))
    allocate(iPhaseCS(nSpeciesCS),iParticlesPerMoleCS(nSpeciesCS))
    allocate(dGibbsMagneticCS(nSpeciesCS,4))

    ! Initialize variables:
    dGibbsCoeffSpeciesTemp = 0D0
    dStoichSpeciesCS       = 0D0
    nGibbsEqSpecies        = 0

    ! Line 3: List of system components:
    read (1,*,IOSTAT = INFO) cElementNameCS(1:nElementsCS)

    if (INFO /= 0) then
        INFO = 103
        return
    end if

    ! Count the number of charged phases and change the name of the system component if it is an electron:
    do i = 1, nElementsCS
        cDummy = cElementNameCS(i)
        if (cDummy(1:2) == 'e(') then
            cElementNameCS(i) = 'e-'
        end if
    end do

    ! Allocate variables specific to sublattice phases:
    allocate(nSublatticePhaseCS(nSolnPhasesSysCS))
    allocate(dStoichSublatticeCS(nSolnPhasesSysCS,nMaxSublatticeCS))
    allocate(nConstituentSublatticeCS(nSolnPhasesSysCS,nMaxSublatticeCS))
    allocate(cConstituentNameSUBCS(nSolnPhasesSysCS,nMaxSublatticeCS,nMaxSpeciesPhaseCS))
    allocate(iConstituentSublatticeCS(nSolnPhasesSysCS,nMaxSublatticeCS,nMaxSpeciesPhaseCS))
    allocate(iSublatticeElementsCS(nSolnPhasesSysCS,nMaxSublatticeCS,nElementsCS))
    allocate(nSublatticeElementsCS(nSolnPhasesSysCS,nMaxSublatticeCS))

    ! Initialize variables:
    nSublatticePhaseCS       = 0
    nConstituentSublatticeCS = 0
    iConstituentSublatticeCS = 0
    dStoichSublatticeCS      = 0D0
    iSublatticeElementsCS    = 0
    nSublatticeElementsCS    = 0

    ! Line 4: List of atomic masses of the elements:
    read (1,*,IOSTAT = INFO) dAtomicMass(1:nElementsCS)

    if (INFO /= 0) then
        !INFO = 14
        INFO = 104
        return
    end if

    ! Line 5: Definition of the temperature dependence terms:
    read (1,*,IOSTAT = INFO) iGequation(1:7)
    read (1,*,IOSTAT = INFO) iGequation(1:7)

    ! Temporary (the capability of handling different types of temperature dependence terms is probably
    ! not needed):
    iGequation(1:7) = iGequation(1:7) - (/6,1,2,3,4,5,6/)
    if ((SUM(iGequation) /= 0).OR.(INFO /= 0)) then
        !INFO = 15
        INFO = 105
        return
    end if

    return

end subroutine ParseCSHeader
