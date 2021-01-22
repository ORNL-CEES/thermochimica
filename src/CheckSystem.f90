
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckSystem.f90
    !> \brief   Check for consistency between the system and the data-file.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      Thermochimica.f90
    !> \sa      GetElementName.f90
    !> \todo    Do a better job allocating variables iPairID and dCoordinationNumber.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   10/17/2011      M.H.A. Piro         Original code
    !   05/01/2012      M.H.A. Piro         Corrected storage of mixing component indices when solution species
    !                                       are no longer considered.  This affected the variable iSpeciesPass.
    !                                       This also affected the way that pure condensed phases are represented
    !                                       by iSpeciesPass.
    !   02/17/2012      M.H.A. Piro         Check if allocatable arrays have been allocated and if so, have they
    !                                       changed in dimension.
    !   05/24/2012      M.H.A. Piro         Correct the minimum number of moles of an element allowed =
    !                                        (normalizer) * (machine precision) / (mass balance tolerance)
    !   08/21/2012      M.H.A. Piro         Redefine the tolerance for the minimum number of moles of a solution
    !                                        phase that is introduced to the system by the minimum number of
    !                                        moles of an element in the system.
    !   01/14/2013      M.H.A. Piro         Move code relavent to excess terms into a separate subroutine.
    !   02/05/2013      M.H.A. Piro         The check for the minimum number of system components now includes
    !                                        constraints imposed by charge neutrality.
    !   05/14/2013      M.H.A. Piro         Fixed bug in allocating arrays specific to ionic phases.  This was
    !                                        done when the size of cSolnPhaseName changes. This should be independent.
    !   10/23/2018      M.H.A. Piro         Fixed a bug in allocating variables for SUBG phases.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to ensure that the selection of system components in the
    !! parsed ChemSage data-file and the data provided to Thermochimica are consistent.  System components are
    !! always taken to be chemical elements in Thermochimica, or "elements" for sake of brevity.  An element will
    !! only be considered if thermodynamic data is provided by the data-file and if the mass of that particular
    !! element is provided.  If the mass of a particular element is provided but there isn't any data for it
    !! (via the data-file), that element will not be considered.  Similarly, if thermodynamic data (via the
    !! data-file) is provided for a particular element, but the mass of that element is not available, that
    !! element will not be considered.
    !!
    !! This subroutine will select all species and phases that are relevant to this system from the data
    !! parsed from the ChemSage data-file.
    !!
    !! The variable cInputThermo(3) can accept the following values:
    !!
    ! cThermoInputUnits(3)      Description
    ! --------------------      -----------
    !> \details
    !!
    !! <table border="1" width="800">
    !! <tr>
    !!    <td> <b> Units </td> <td> Description </b> </td>
    !! </tr>
    !! <tr>
    !!    <td> "mass fraction" </td>
    !!    <td> Mass fraction in dimensionless units (e.g., gram/gram, kilogram/kilogram, pound/pound, wt%).  </td>
    !! </tr>
    !! <tr>
    !!    <td> "mole fraction" </td>
    !!    <td> Mole fraction in dimensionless units (e.g., mole/mole, mol%).  </td>
    !! </tr>
    !! <tr>
    !!    <td> "atom fraction" </td>
    !!    <td> Atom fraction in dimensionless units (e.g., atom/atom, at%).  </td>
    !! </tr>
    !! <tr>
    !!    <td> "kilograms" </td>
    !!    <td> All quantities are in kilograms.  </td>
    !! </tr>
    !! <tr>
    !!    <td> "grams" </td>
    !!    <td> All quantities are in grams.  </td>
    !! </tr>
    !! <tr>
    !!    <td> "pounds" or "lbs" </td>
    !!    <td> All quantities are in pounds.  </td>
    !! </tr>
    !! <tr>
    !!    <td> "moles" </td>
    !!    <td> All quantities are in moles.  </td>
    !! </tr>
    !! <tr>
    !!    <td> "gram-atoms" </td>
    !!    <td> All quantities are in gram-atoms (same as moles for the pure elements);  </td>
    !! </tr>
    !! <tr>
    !!    <td> "atoms"   </td>
    !!    <td>  All quantities are in atoms.  </td>
    !! </tr>
    !! </table>
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dNormalizeInput                   A double scalar that normalizes the mass units.  This is performed to
    !                                    minimize numerical error when evaluating the Jacobian/Broyden matrix in
    !                                    the GEMSolver.
    !
    ! dTemperature                      Temperature [K]
    ! dPressure                         Absolute hydrostatic pressure [atm]
    ! dElementMass                      Total mass of each element, where the coefficient corresponds to the
    !                                    atomic number (e.g., dMolesElement(92) refers to uranium).
    ! dElementMoleFractionMin           A minimum allowable value for the number of moles of an element.
    ! INFOThermo                        A scalar integer that indicates a successful exit or identifies an error.
    ! cInputThermo                      A character vector containing the units for temperature, pressure and
    !                                    elemental quantity.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckSystem

    USE ModuleParseCS
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver, ONLY: lSolnPhases, lMiscibility

    implicit none

    integer                                 :: i, j, k, l, m, n, nMaxSpeciesPhase, nCountSublatticeTemp, iCon1, iCon2, iCon3, iCon4
    integer,dimension(0:nSolnPhasesSysCS+1) :: iTempVec
    real(8)                                 :: dSum, dElementMoleFractionMin
    character(3),dimension(0:nElementsPT)   :: cElementNamePT
    character(12)                           :: cDummy
    logical                                 :: lPos, lNeg


    ! Check to see if the allocatable arrays have already been allocated
    if (allocated(iElementSystem)) then
        ! Do nothing.
    else
        ! Allocate memory:
        allocate(iElementSystem(1:nElementsCS))
        allocate(iSpeciesPass(nSpeciesCS))
    end if

    ! Check if there are any solution phases with a sublattice:
    if (nCountSublatticeCS > 0) then
        ! Allocate array to check if a constituent passes:
        j = SIZE(iConstituentSublatticeCS, DIM=3)
        if (allocated(iConstituentPass)) deallocate(iConstituentPass)
        allocate(iConstituentPass(nCountSublatticeCS,nMaxSublatticeCS,j))
        iConstituentPass = 0
    end if

    ! Initialize variables:
    iElementSystem      = 0
    iTempVec            = 0
    iSpeciesPass        = 0
    nCountSublattice    = 0
    nMaxSublatticeSys   = 0
    nMaxConstituentSys  = 0
    nMaxSpeciesPhase    = 0
    nCountSublatticeTemp   = 0
    nChargedConstraints = 0
    n                   = 0
    dSum                = 0D0
    dElementMoleFractionMin = dTolerance(6)

    ! Get the name of all the elements on the periodic table:
    call GetElementName(cElementNamePT)

    ! Perform a mass conversion:
    LOOP_Small: do j = 1, nElementsCS
        ! Map the character string of this element to its atomic number:
        LOOP_Big: do i = 0, nElementsPT
            if (cElementNameCS(j) == cElementNamePT(i)) then
                ! Convert the mass of each element to moles:
                select case (cInputUnitMass)
                case ('mass fraction','kilograms','grams','pounds','lbs','g','kg')
                        ! Convert mass unit to moles:
                        dElementMass(i) = dElementMass(i) / dAtomicMassCS(j)
                    case ('mole fraction','atom fraction','atoms','moles','mol','gram-atoms')
                        ! Do nothing
                    case default
                        ! The character string representing input units is not recognized.
                        INFOThermo = 4
                        return
                end select
                iElementSystem(j) = i
                dSum = dSum + dElementMass(i)
                cycle LOOP_Small
            elseif (cElementNameCS(j) == 'e-') then
                ! Electron
                iElementSystem(j) = i
                dElementMass(i)   = 0D0
                cycle LOOP_Small
            end if
        end do LOOP_Big
        ! The chemical element stored by cElementNameCS(j) does not correspond to an
        ! element in cElementNamePT.  Report an error and exit:
        INFOThermo = 31
        return
    end do LOOP_Small

    cInputUnitMass = 'moles'

    if (nCompounds > 0) then
        call CheckCompounds
        dSum = 0
        do i = 1, nElemOrComp
            dSum = dSum + dElementMass(iElementSystem(i))
        end do
    else
        nElemOrComp = nElementsCS
    end if

    if (INFOThermo > 0) then
        return
    end if

    ! Make sure that the sum of element masses is not zero:
    if (dSum == 0D0) then
        INFOThermo = 5
        return
    end if

    ! Make dSum multiplicative and normalize:
    dNormalizeInput = dNormalizeInput / dSum

    k = 0
    ! Normalize dElementMass to mole fraction and establish the elements of the system:
    do j = 1, nElemOrComp
        dElementMass(iElementSystem(j)) = dElementMass(iElementSystem(j)) * dNormalizeInput
        if (dElementMass(iElementSystem(j)) < dElementMoleFractionMin) then
            ! Element j should not be considered.
            dElementMass(iElementSystem(j)) = 0
            iElementSystem(j) = 0
            if (cElementNameCS(j) == 'e-') then
                nElements         = nElements + 1
                iElementSystem(j) = -1
                k = k + 1
            end if
        else
            nElements = nElements + 1
        end if
    end do

    LOOP_checkElements: do i = 1, nElementsPT
        if (dElementMass(i) > 0) then
            do j = 1, nElementsCS
                if (iElementSystem(j) == i) cycle LOOP_checkElements
            end do
            print *, "WARNING: Element ", cElementNamePT(i), " not in database and therefore omitted from calculation"
        end if
    end do LOOP_checkElements

    ! The system requires a minimum of two elements in order to be considered:
    if (nElements < 2 + k) then
        INFOThermo = 5
        return
    end if

    ! Check to see if the system has to be re-adjusted:
    nSpecies = 0
    LOOP_SolnPhases: do i = 1, nSolnPhasesSysCS
        ! Loop through species in solution phases:
        m = 0
        LOOP_SpeciesInSolnPhase: do j = nSpeciesPhaseCS(i-1) + 1, nSpeciesPhaseCS(i)
            ! if (SUM(DABS(dStoichSpeciesCS(j,1:nElemOrComp))) == 0) cycle LOOP_SpeciesInSolnPhase
            do k = 1, nElemOrComp
                if ((dStoichSpeciesCS(j,k) > 0).AND.(iElementSystem(k) == 0)) then
                    ! This species should not be considered
                    cycle LOOP_SpeciesInSolnPhase
                end if
            end do
            if (cSolnPhaseTypeCS(i) == 'SUBI') then
                k = iPhaseSublatticeCS(i)
                iCon1 = iConstituentSublatticeCS(k,2,j-nSpeciesPhaseCS(i-1))
                ! If this is an ionic liquid and there are identical neutrals, skip the later ones
                if (dSublatticeChargeCS(k,2,iCon1) == 0D0) then
                    do l = nSpeciesPhaseCS(i-1) + 1, j - 1
                        iCon2 = iConstituentSublatticeCS(k,2,l-nSpeciesPhaseCS(i-1))
                        if ((dSublatticeChargeCS(k,2,iCon2) == 0D0) .AND. (iCon1 == iCon2)) then
                              cycle LOOP_SpeciesInSolnPhase
                        end if
                    end do
                end if
            end if
            nSpecies = nSpecies + 1
            m = m + 1
            iSpeciesPass(j) = m
            if (cSolnPhaseTypeCS(i) == 'SUBG' .OR. cSolnPhaseTypeCS(i) == 'SUBQ') then
                k = iPhaseSublatticeCS(i)
                iCon1 = iPairIDCS(k,j-nSpeciesPhaseCS(i-1),1)
                iCon2 = iPairIDCS(k,j-nSpeciesPhaseCS(i-1),2)
                iCon3 = iPairIDCS(k,j-nSpeciesPhaseCS(i-1),3) - nSublatticeElementsCS(k,1)
                iCon4 = iPairIDCS(k,j-nSpeciesPhaseCS(i-1),4) - nSublatticeElementsCS(k,1)
                if (iCon1 > 0) iConstituentPass(k,1,iCon1) = 1
                if (iCon2 > 0) iConstituentPass(k,1,iCon2) = 1
                if (iCon3 > 0) iConstituentPass(k,2,iCon3) = 1
                if (iCon4 > 0) iConstituentPass(k,2,iCon4) = 1
            end if
        end do LOOP_SpeciesInSolnPhase

        ! For electrons, there have to be species with positive and negative stoichiometries still in the system.
        ! Otherwise, remove the species that use that electron.
        do k = 1, nElementsCS
            if (cElementNameCS(k) == 'e-') then
                lPos = .FALSE.
                lNeg = .FALSE.
                do j = nSpeciesPhaseCS(i-1) + 1, nSpeciesPhaseCS(i)
                    if ((iSpeciesPass(j) > 0) .AND. dStoichSpeciesCS(j,k) > 0D0) lPos = .TRUE.
                    if ((iSpeciesPass(j) > 0) .AND. dStoichSpeciesCS(j,k) < 0D0) lNeg = .TRUE.
                end do
                if (lPos .NEQV. lNeg) then
                    do j = nSpeciesPhaseCS(i-1) + 1, nSpeciesPhaseCS(i)
                        if ((iSpeciesPass(j) > 0) .AND. DABS(dStoichSpeciesCS(j,k)) > 0D0) then
                            iSpeciesPass(j) = 0
                            nSpecies = nSpecies - 1
                        end if
                    end do
                end if
            end if
        end do

        ! Store temporary counter for the number of charged phases from the CS data-file:
        if ((cSolnPhaseTypeCS(i) == 'SUBL').OR.(cSolnPhaseTypeCS(i) == 'SUBLM').OR. &
             (cSolnPhaseTypeCS(i) == 'SUBG').OR.(cSolnPhaseTypeCS(i) == 'SUBQ').OR. &
             (cSolnPhaseTypeCS(i) == 'SUBI')) then
            nCountSublatticeTemp = nCountSublatticeTemp + 1
        end if

        ! Count the number of solution phases in the system:
        iTempVec(nSolnPhasesSys+1) = nSpecies
        if (iTempVec(nSolnPhasesSys+1) - iTempVec(nSolnPhasesSys) > 1) then
            nSolnPhasesSys = nSolnPhasesSys + 1
            nMaxSpeciesPhase = MAX(nMaxSpeciesPhase, iTempVec(nSolnPhasesSys) - iTempVec(nSolnPhasesSys-1))
            ! Check if this is a charged phase:
            if ((cSolnPhaseTypeCS(i) == 'SUBL').OR.(cSolnPhaseTypeCS(i) == 'SUBLM').OR. &
                 (cSolnPhaseTypeCS(i) == 'SUBG').OR.(cSolnPhaseTypeCS(i) == 'SUBQ').OR. &
                 (cSolnPhaseTypeCS(i) == 'SUBI')) then
                ! Count the number of charged phases:
                nCountSublattice = nCountSublattice + 1
                ! Determine the maximum number of sublattice of any stable phase:
                nMaxSublatticeSys  = MAX(nMaxSublatticeSys,nSublatticePhaseCS(nCountSublatticeTemp))
                m = MAXVAL(nConstituentSublatticeCS(nCountSublatticeTemp,1:nMaxSublatticeCS))
                nMaxConstituentSys = MAX(nMaxConstituentSys,m)
            end if
        elseif (iTempVec(nSolnPhasesSys+1) - iTempVec(nSolnPhasesSys) == 1) then
            ! There is only one species in this solution phase.  This solution phase should not be considered.
            iSpeciesPass(nSpeciesPhaseCS(i-1) + 1 : nSpeciesPhaseCS(i))          = 0
            nSpecies                 = nSpecies - 1
            iTempVec(nSolnPhasesSys) = nSpecies
        else
            ! Do nothing.  The number of species in this solution phase is zero and this phase will not be
            ! considered.
        end if
    end do LOOP_SolnPhases

    ! Loop through pure condensed phases:
    LOOP_PureConPhases: do j = nSpeciesPhaseCS(nSolnPhasesSysCS) + 1, nSpeciesCS
        if (SUM(dStoichSpeciesCS(j,1:nElemOrComp)) == 0) cycle LOOP_PureConPhases
        do k = 1, nElemOrComp
            if ((dStoichSpeciesCS(j,k) > 0).AND.(iElementSystem(k) == 0)) then
                ! This species should not be considered
                cycle LOOP_PureConPhases
            end if
        end do
        nSpecies        = nSpecies + 1
        iSpeciesPass(j) = 1
    end do LOOP_PureConPhases

    ! Re-establish the character vector representing the element names:
    j = 0
    do i = 1, nElemOrComp
        if (iElementSystem(i) /= 0) then
            cDummy           = cElementNameCS(i)
            if (cDummy == 'e-') nChargedConstraints = nChargedConstraints + 1
        end if
    end do

    ! Add dummy species representing electrons to the number of species in the system:
    nSpecies = nSpecies + nChargedConstraints

    ! Check if these variables have already been allocated:
    if (allocated(dChemicalPotential)) then
        ! Check to see if the number of species has changed:

        i = SIZE(dChemicalPotential)
        if (i /= nSpecies) then
            ! The number of species has changed.
            deallocate(dChemicalPotential,iPhase,dSpeciesTotalAtoms,cSpeciesName,&
                iParticlesPerMole,dStdGibbsEnergy,dCoeffGibbsMagnetic,dMagGibbsEnergy, STAT = n)
            if (n /= 0) then
                INFOThermo = 19
                return
            end if
            ! Allocate memory for variables:
            allocate(dChemicalPotential(nSpecies),iPhase(nSpecies),dSpeciesTotalAtoms(nSpecies))
            allocate(cSpeciesName(nSpecies),dStdGibbsEnergy(nSpecies))
            allocate(iParticlesPerMole(nSpecies), dCoeffGibbsMagnetic(nSpecies,4), dMagGibbsEnergy(nSpecies))
        end if

        ! Check to see if the number of elements has changed:
        j = SIZE(cElementName)
        if (j /= nElements) then
            ! The number of elements has changed.
            deallocate(cElementName,dMolesElement,dAtomicMass, STAT = n)
            if (n /= 0) then
                INFOThermo = 19
                return
            end if
            ! Allocate memory:
            allocate(cElementName(nElements),dMolesElement(nElements),dAtomicMass(nElements))
        end if

        ! Check to see if either the number of species or the number of elements has changed:
        if ((i /= nSpecies).OR.(j /= nElements)) then
            deallocate(dAtomFractionSpecies,dStoichSpecies, STAT = n)
            if (n /= 0) then
                INFOThermo = 19
                return
            end if
            ! Allocate memory:
            allocate(dAtomFractionSpecies(nSpecies,nElements),dStoichSpecies(nSpecies,nElements))
        end if

        ! Check to see if the number of solution phases in the system has changed:
        k = SIZE(cSolnPhaseName)
        if (k /= nSolnPhasesSys) then
            ! The number of solution phases in the system has changed.
            deallocate(nSpeciesPhase,nParamPhase,cSolnPhaseType,cSolnPhaseName,lSolnPhases,dGibbsSolnPhase, &
                lMiscibility,nMagParamPhase, STAT = n)
            if (n /= 0) then
                INFOThermo = 19
                return
            end if
            ! Allocate memory:
            allocate(nSpeciesPhase(0:nSolnPhasesSys),nParamPhase(0:nSolnPhasesSys),nMagParamPhase(0:nSolnPhasesSys))
            allocate(cSolnPhaseType(nSolnPhasesSys),cSolnPhaseName(nSolnPhasesSys))
            allocate(lSolnPhases(nSolnPhasesSys),dGibbsSolnPhase(nSolnPhasesSys),lMiscibility(nSolnPhasesSys))

        end if

        ! Only allocate if there are sublattice phases:
        if (nCountSublattice > 0) then

            deallocate(iPhaseSublattice,nSublatticePhase,nConstituentSublattice,dStoichSublattice, &
                    dSiteFraction,cConstituentNameSUB,iConstituentSublattice,nSublatticeElements, &
                     nPairsSRO,iPairID,dCoordinationNumber, dZetaSpecies, &
                     dSublatticeCharge,iChemicalGroup,dStoichPairs,cPairName,dConstituentCoefficients, STAT = n)

            allocate(iPhaseSublattice(nSolnPhasesSys),nSublatticePhase(nCountSublattice))
            allocate(nConstituentSublattice(nCountSublattice,nMaxSublatticeSys))
            allocate(dStoichSublattice(nCountSublattice,nMaxSublatticeSys))
            allocate(dSiteFraction(nCountSublattice,nMaxSublatticeSys,nMaxConstituentSys))
            allocate(cConstituentNameSUB(nCountSublattice,nMaxSublatticeSys,MAX(nMaxConstituentSys,nMaxSpeciesPhase)))
            allocate(iConstituentSublattice(nCountSublattice,nMaxSublatticeSys,nMaxSpeciesPhase))
            allocate(nSublatticeElements(nCountSublattice,nMaxSublatticeSys))
            j = MAXVAL(nSublatticeElementsCS)
            allocate(nPairsSRO(nCountSublattice,2))
            allocate(iPairID(nCountSublattice,nMaxSpeciesPhase,4))
            allocate(dCoordinationNumber(nCountSublattice,nMaxSpeciesPhase,4))
            allocate(dZetaSpecies(nCountSublattice,nMaxSpeciesPhase))
            allocate(dConstituentCoefficients(nCountSublattice,nMaxSpeciesPhase,5))
            allocate(dSublatticeCharge(nCountSublattice,nMaxSublatticeSys,j))
            allocate(iChemicalGroup(nCountSublattice,nMaxSublatticeSys,j))
            allocate(dStoichPairs(nCountSublattice,MAXVAL(nPairsSROCS(:,1)),nElements))
            allocate(cPairName(nCountSublattice,MAXVAL(nPairsSROCS(:,1))))
        end if

    else

        ! Allocate memory for variables:
        allocate(dChemicalPotential(nSpecies),iPhase(nSpecies),dSpeciesTotalAtoms(nSpecies))
        allocate(cSpeciesName(nSpecies),dStdGibbsEnergy(nSpecies))
        allocate(iParticlesPerMole(nSpecies),dCoeffGibbsMagnetic(nSpecies,4),dMagGibbsEnergy(nSpecies))
        allocate(cElementName(nElements),dMolesElement(nElements),dAtomicMass(nElements))
        allocate(dAtomFractionSpecies(nSpecies,nElements),dStoichSpecies(nSpecies,nElements))
        allocate(nSpeciesPhase(0:nSolnPhasesSys),nParamPhase(0:nSolnPhasesSys),nMagParamPhase(0:nSolnPhasesSys))
        allocate(cSolnPhaseType(nSolnPhasesSys),cSolnPhaseName(nSolnPhasesSys))
        allocate(lSolnPhases(nSolnPhasesSys),dGibbsSolnPhase(nSolnPhasesSys),lMiscibility(nSolnPhasesSys))

        ! Only allocate if there are charged phases:
        if (nCountSublattice > 0) then
            allocate(iPhaseSublattice(nSolnPhasesSys),nSublatticePhase(nCountSublattice))
            allocate(nConstituentSublattice(nCountSublattice,nMaxSublatticeSys))
            allocate(dStoichSublattice(nCountSublattice,nMaxSublatticeSys))
            allocate(dSiteFraction(nCountSublattice,nMaxSublatticeSys,nMaxConstituentSys))
            allocate(cConstituentNameSUB(nCountSublattice,nMaxSublatticeSys,MAX(nMaxConstituentSys,nMaxSpeciesPhase)))
            allocate(iConstituentSublattice(nCountSublattice,nMaxSublatticeSys,nMaxSpeciesPhase))
            allocate(nSublatticeElements(nCountSublattice,nMaxSublatticeSys))
            j = MAXVAL(nSublatticeElementsCS)
            allocate(nPairsSRO(nCountSublattice,2))
            allocate(iPairID(nCountSublattice,nMaxSpeciesPhase,4))
            allocate(dCoordinationNumber(nCountSublattice,nMaxSpeciesPhase,4))
            allocate(dZetaSpecies(nCountSublattice,nMaxSpeciesPhase))
            allocate(dConstituentCoefficients(nCountSublattice,nMaxSpeciesPhase,5))
            allocate(dSublatticeCharge(nCountSublattice,nMaxSublatticeSys,j))
            allocate(iChemicalGroup(nCountSublattice,nMaxSublatticeSys,j))
            allocate(dStoichPairs(nCountSublattice,MAXVAL(nPairsSROCS(:,1)),nElements))
            allocate(cPairName(nCountSublattice,MAXVAL(nPairsSROCS(:,1))))
        end if

    end if

    ! Initialize variables:
    iPhase               = 0
    dSpeciesTotalAtoms   = 0
    iParticlesPerMole    = 0
    dStoichSpecies       = 0
    nSpeciesPhase        = 0
    nParamPhase          = 0
    nMagParamPhase       = 0
    dChemicalPotential   = 0D0
    dStdGibbsEnergy      = 0D0
    dMolesElement        = 0D0
    dAtomFractionSpecies = 0D0
    dGibbsSolnPhase      = 0D0
    dCoeffGibbsMagnetic  = 0D0
    dMagGibbsEnergy      = 0D0
    lSolnPhases          = .FALSE.
    lMiscibility         = .FALSE.

    ! Initialize arrays (if necessary) for sublattice phases:
    if (nCountSublattice > 0) then
        dSiteFraction          = 0D0
        dStoichSublattice      = 0D0
        iConstituentSublattice = 0
        nSublatticePhase       = 0
        nConstituentSublattice = 0
        nSublatticeElements  = 0
        nPairsSRO            = 0
        iPairID              = 0
        dCoordinationNumber  = 0D0
        dZetaSpecies         = 0D0
        dConstituentCoefficients = 0D0
        dSublatticeCharge    = 0D0
        iChemicalGroup       = 0
        dStoichPairs         = 0D0
        cConstituentNameSUB  = ' '
        cPairName            = ' '
    end if

    ! Re-establish the character vector representing the element names:
    j = 0
    do i = 1, nElemOrComp
        if (iElementSystem(i) /= 0) then
            j = j + 1
            if (nCompounds == 0) then
                cElementName(j)  = cElementNameCS(i)
                cDummy           = cElementName(j)
                dAtomicMass(j)   = dAtomicMassCS(i)
            else
                cElementName(i)  = cCompoundNames(i)
            end if
            if (iElementSystem(i) > 0) dMolesElement(j) = dElementMass(iElementSystem(i))
        end if
    end do

    ! Redefine the tolerance for the minimum number of moles of a solution phase that is introduced to the system:
    dTolerance(9) = DMIN1(1000D0 * MINVAL(dMolesElement, MASK = dMolesElement > 0D0),dTolerance(9))

    ! Re-establish the nSpeciesPhase vector:
    nSpeciesPhase(0:nSolnPhasesSys) = iTempVec(0:nSolnPhasesSys)

    ! Check the excess terms:
    call CheckSystemExcess

    return

end subroutine CheckSystem
