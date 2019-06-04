
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompThermoData.f90
    !> \brief   Compute thermodynamic data
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      Thermochimica.f90
    !> \sa      CompGibbsMagnetic.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   10/11/2011      M.H.A. Piro         Original code
    !   05/01/2012      M.H.A. Piro         Corrected storage of mixing component indices when solution species
    !                                        are no longer considered.  This affected the variables iSpeciesPass
    !                                        and iRegularParam.
    !   01/31/2012      M.H.A. Piro         The standard molar Gibbs energy coefficient of gaseous species is
    !                                        adjusted by an arbitrary quantity (see comments below).
    !   01/17/2013      M.H.A. Piro         Additional dummy species are included to the system for each electron
    !                                        representing a system component.
    !   02/11/2013      M.H.A. Piro         Mixing parameters for SUBL phases have been incorporated.  These are
    !                                        different than other phases because the constituent indices on each
    !                                        sublattice are stored, rather than the species indices.
    !   07/16/2013      M.H.A. Piro         Fix bug in constructing the iRegularParam vector when dealing with SUBL
    !                                        phases when the mixing constituents are not on the first sublattice.
    !   03/19/2018      M.H.A. Piro         Added capability to handle SUBG phases. This includes a special case
    !                                        to define the stoichiometry of pairs, rather than pure species. Otherwise
    !                                        they would be undefined.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute thermodynamic data for all substances in the
    !! system.  The specific variables that are computed include:
    !!   - the standard molar Gibbs energy for each pure substance using the specified temperature and pressure,
    !!   - the excess molar Gibbs energy of mixing for each sub-system,
    !!   - the total number of atoms per formula mass for each compound,
    !!   - the atomic fraction of each element in a particular compound, and
    !!   - construction of the Hessian matrix, which will be used later by the GEMNewton subroutine.
    !!
    !! The coefficients for the standard molar Gibbs energy equations for pure species originate from a ChemSage
    !! data-file that was parsed from the ParseCSDataFile program.  The format for the coefficients follow:
    !!
    !! \f$ g_i^{\circ} = A + BT + CTln(T) + DT^2 + ET^3 + F/T + (GT^U + HT^V + IT^W + JT^X + KT^Y + LT^Z) \f$
    !!
    !! Note that the terms in paranteses in the above equation are additional terms.  Some of the
    !! exponents used in the additional terms may be 99.  This corresponds to the natural logarithm.
    !!
    !! For some very strange reason, the B coefficient in the above equation is modified by ~0.10945 J/mol for
    !! only gaseous species when FactSage generates a ChemSage data-file.  Another very peculiar observation is
    !! that if a database is constructed in FactSage using the new ChemSage data-file, FactSage will automatically
    !! remove this quantity from the B coefficient.  The B coefficient for gaseous species is corrected in this
    !! subroutine.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nSpecies              The number of species in the system.
    ! nSpeciesCS            The number of species in the original ChemSage data-file.  Note: this may be
    !                        higher than what may be considered in this system.
    ! nSpeciesPhase         An integer vector representing the number of species in each solution phase.
    ! nSpeciesPhaseCS       An integer vector representing the number of species in each solution phase
    !                        from the ChemSage data-file (note: this may be larger than nSpeciesPhase).
    ! dChemicalPotential    A double real vector representing the chemical potential of each species.
    !                        To be precise, this is defined as the difference between the standard molar
    !                        Gibbs energy and the chemical potential defined by the element potentials
    !                        (represented in dimensionless units and per formula mass).
    ! dExcessGibbsParam     A double real vector representing the excess Gibbs energy of mixing
    !                        parameters.
    ! dAtomFractionSpecies  A double real matrix representing tha atom fraction of each element in each
    !                        species.
    ! dJacobianLong         A double real matrix that is used to construct the Jacobian matrix in
    !                        GEMNewton.f90. This matrix is constructed once and never modified throughout
    !                        a calculation.
    ! iParticlesPerMole     The number of particles per mole of constituent.
    ! iSpeciesAtoms         An integer matrix representing the number of atoms of a particular element
    !                        for a particular species.
    ! iSpeciesAtomsCS       An integer matrix representing the number of atoms of a particular element
    !                        for a particular species from the database (note: this may be larger than
    !                        iSpecesAtoms).
    ! iSpeciesTotalAtoms    The total number of atoms per formula mass of a species.
    ! iRegularParam         An integer vector representing important mixing terms for a regular solution
    !                        model.
    ! iRegularParamCS       An integer vector representing important mixing terms for a regular solution
    !                        model from the ChemSage data-file (Note: this may be larger than
    !                        iRegularParam).
    ! dGibbsCoeffSpeciesTemp  A double real array containing the coefficients of Gibbs energy equations
    !                        of species in the database.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CompThermoData

    USE ModuleParseCS
    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer                            :: i, j, k, l, m, n, s, iCounterGibbsEqn, nCounter
    integer                            :: ii, jj, kk, ll, ka, la, iax, iay, ibx, iby
    integer                            :: iSublPhaseIndex, iFirst, nRemove, nA2X2
    integer                            :: iaaxx, ibbxx, iaayy, ibbyy
    integer, dimension(nElementsCS)    :: iRemove
    real(8)                            :: dLogT, dLogP, dTemp, dQx, dQy, dZa, dZb, dZx, dZy
    real(8)                            :: dZaxa, dZbxb, dZaya, dZbyb
    real(8), dimension(6)              :: dGibbsCoeff
    real(8), dimension(nSpeciesCS)     :: dChemicalPotentialTemp


    ! Initialize variables:
    j                = 0
    iCounterGibbsEqn = 0
    nDummySpecies    = 0
    nCounter         = 0
    dTemp            = 1D0 / (dIdealConstant * dTemperature)

    ! Compute Gibbs energy coefficients:
    dGibbsCoeff(1)   = 1D0                             ! A
    dGibbsCoeff(2)   = dTemperature                    ! B
    dGibbsCoeff(3)   = dTemperature*DLOG(dTemperature) ! C
    dGibbsCoeff(4)   = dTemperature**2                 ! D
    dGibbsCoeff(5)   = dTemperature**3                 ! E
    dGibbsCoeff(6)   = 1D0 / dTemperature              ! F
    dLogT            = DLOG(dTemperature)              ! ln(T)
    dLogP            = DLOG(dPressure)                 ! ln(P)

    ! Loop through all species in the system:
    LOOP_nPhasesCS: do n = 1, nSolnPhasesSysCS
        if (SUM(iSpeciesPass(nSpeciesPhaseCS(n-1)+1:nSpeciesPhaseCS(n))) == 0) cycle LOOP_nPhasesCS
        if ((cSolnPhaseTypeCS(n) == 'SUBG') .OR. (cSolnPhaseTypeCS(n) == 'SUBQ')) then
            iSublPhaseIndex = iPhaseSublatticeCS(n)
            iFirst = nSpeciesPhaseCS(n - 1) + 1
            dChemicalPotentialTemp = 0D0
            jj = 0
            LOOP_SROPairs: do i = iFirst, iFirst - 1 + nPairsSROCS(iSublPhaseIndex,1)
                l = 0
                ! Loop through the Gibbs energy equations to figure out which one to use:
                do k = 1, nGibbsEqSpecies(i)
                    iCounterGibbsEqn = iCounterGibbsEqn + 1
                    if ((dTemperature <= dGibbsCoeffSpeciesTemp(1,iCounterGibbsEqn)).AND.(l == 0)) then
                        l = k
                    end if
                end do

                if (l == 0) l = nGibbsEqSpecies(i)

                l = l + iCounterGibbsEqn - nGibbsEqSpecies(i)

                do k = 2, 7
                    dChemicalPotentialTemp(i) = dChemicalPotentialTemp(i) + dGibbsCoeffSpeciesTemp(k,l) * dGibbsCoeff(k-1)
                end do

                ! Compute additional standard molar Gibbs energy terms:
                if (dGibbsCoeffSpeciesTemp(9,l) .EQ. 99) then
                    dChemicalPotentialTemp(i) = dChemicalPotentialTemp(i) + dGibbsCoeffSpeciesTemp(8,l) * dLogT
                else
                    dChemicalPotentialTemp(i) = dChemicalPotentialTemp(i) + dGibbsCoeffSpeciesTemp(8,l) &
                        *dTemperature**dGibbsCoeffSpeciesTemp(9,l)
                end if

                if (dGibbsCoeffSpeciesTemp(11,l) .EQ. 99) then
                    dChemicalPotentialTemp(i) = dChemicalPotentialTemp(i) + dGibbsCoeffSpeciesTemp(10,l) * dLogT
                else
                    dChemicalPotentialTemp(i) = dChemicalPotentialTemp(i) + dGibbsCoeffSpeciesTemp(10,l) &
                        * dTemperature**dGibbsCoeffSpeciesTemp(11,l)
                end if

                if (dGibbsCoeffSpeciesTemp(13,l) .EQ. 99) then
                    dChemicalPotentialTemp(i) = dChemicalPotentialTemp(i) + dGibbsCoeffSpeciesTemp(12,l) * dLogT
                else
                    dChemicalPotentialTemp(i) = dChemicalPotentialTemp(i) + dGibbsCoeffSpeciesTemp(12,l) &
                        * dTemperature**dGibbsCoeffSpeciesTemp(13,l)
                end if

                ! Convert chemical potentials to dimensionless units:
                dChemicalPotentialTemp(i) = dChemicalPotentialTemp(i) * dTemp * DFLOAT(iParticlesPerMoleCS(i))

                do k = 1, nElemOrComp
                    if ((dStoichPairsCS(iSublPhaseIndex,i - iFirst + 1,k) > 0).AND.(iElementSystem(k) == 0)) then
                        cycle LOOP_SROPairs
                    end if
                end do

                ! Check if pair should be saved - only worry about zeta here, LOOP_nSUBGQCS will take care
                ! of only using the necessary reference energy terms.
                ! SUBG/Q internal consituent indices are in iConstituentSublatticeCS
                ! iSublatticeElementsCS converts these to global element indices
                ! The rule to check is are elements A (ii) and X (kk) both in the system:
                ii = iSublatticeElementsCS(iSublPhaseIndex,1,iConstituentSublatticeCS(iSublPhaseIndex,1,i - iFirst + 1))
                kk = iSublatticeElementsCS(iSublPhaseIndex,2,iConstituentSublatticeCS(iSublPhaseIndex,2,i - iFirst + 1))
                if ((ii > 0) .AND. (kk > 0)) then
                    jj = jj + 1
                    dZetaSpecies(iSublPhaseIndex,jj) = dZetaSpeciesCS(iSublPhaseIndex,i - iFirst + 1)
                end if

                ! I'm like pretty sure that these are g_A2/X2 and not g_A/X,
                ! but only because it doesn't work the other way.
                ! Also the seem to be multiplied by the relevant coordination already.
                dChemicalPotentialTemp(i) = dChemicalPotentialTemp(i) !* 4 / dZetaSpecies(iSublPhaseIndex,jj)

            end do LOOP_SROPairs

            LOOP_nSUBGQCS: do i = nSpeciesPhaseCS(n - 1) + 1, nSpeciesPhaseCS(n)
                if (iSpeciesPass(i) == 0) cycle LOOP_nSUBGQCS

                j = j + 1   ! New species index

                cSpeciesName(j)            = cSpeciesNameCS(i)
                iPhase(j)                  = iPhaseCS(i)
                iParticlesPerMole(j)       = iParticlesPerMoleCS(i)
                dCoeffGibbsMagnetic(j,1:4) = dGibbsMagneticCS(i,1:4)

                m = 0
                do k = 1, nElemOrComp
                    if (iElementSystem(k) /= 0) then
                        m = m + 1
                        dStoichSpecies(j,m) = dStoichSpeciesCS(i,k)
                    end if
                end do

                dZa = dCoordinationNumberCS(iSublPhaseIndex,i - iFirst + 1,1)
                dZb = dCoordinationNumberCS(iSublPhaseIndex,i - iFirst + 1,2)
                dZx = dCoordinationNumberCS(iSublPhaseIndex,i - iFirst + 1,3)
                dZy = dCoordinationNumberCS(iSublPhaseIndex,i - iFirst + 1,4)

                ii = iPairIDCS(iSublPhaseIndex,i - iFirst + 1,1)
                jj = iPairIDCS(iSublPhaseIndex,i - iFirst + 1,2)
                kk = iPairIDCS(iSublPhaseIndex,i - iFirst + 1,3)
                ll = iPairIDCS(iSublPhaseIndex,i - iFirst + 1,4)
                ! Anion indices adjusted to start from 1
                ka = kk - nSublatticeElementsCS(iSublPhaseIndex,1)
                la = ll - nSublatticeElementsCS(iSublPhaseIndex,1)

                dQx = dSublatticeChargeCS(iSublPhaseIndex,2,ka)
                dQy = dSublatticeChargeCS(iSublPhaseIndex,2,la)

                ! Indices of reference energy equations for A/X, B/X, A/Y, B/Y
                iax = 0
                ibx = 0
                iay = 0
                iby = 0
                nA2X2 = nSublatticeElementsCS(iSublPhaseIndex,1) * nSublatticeElementsCS(iSublPhaseIndex,2)
                do k = 1, nA2X2
                    if   ((iConstituentSublatticeCS(iSublPhaseIndex,1,k) == ii) &
                    .AND. (iConstituentSublatticeCS(iSublPhaseIndex,2,k) == ka)) then
                        iax = k
                    end if
                    if   ((iConstituentSublatticeCS(iSublPhaseIndex,1,k) == jj) &
                    .AND. (iConstituentSublatticeCS(iSublPhaseIndex,2,k) == ka)) then
                        ibx = k
                    end if
                    if   ((iConstituentSublatticeCS(iSublPhaseIndex,1,k) == ii) &
                    .AND. (iConstituentSublatticeCS(iSublPhaseIndex,2,k) == la)) then
                        iay = k
                    end if
                    if   ((iConstituentSublatticeCS(iSublPhaseIndex,1,k) == jj) &
                    .AND. (iConstituentSublatticeCS(iSublPhaseIndex,2,k) == la)) then
                        iby = k
                    end if
                end do

                iaaxx = ii + (ka - 1) * (nSublatticeElements(iSublPhaseIndex,1) * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2)
                ibbxx = jj + (ka - 1) * (nSublatticeElements(iSublPhaseIndex,1) * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2)
                iaayy = ii + (la - 1) * (nSublatticeElements(iSublPhaseIndex,1) * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2)
                ibbyy = jj + (la - 1) * (nSublatticeElements(iSublPhaseIndex,1) * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2)
                dZaxa = dCoordinationNumberCS(iSublPhaseIndex,iaaxx,1)
                dZbxb = dCoordinationNumberCS(iSublPhaseIndex,ibbxx,2)
                dZaya = dCoordinationNumberCS(iSublPhaseIndex,iaayy,1)
                dZbyb = dCoordinationNumberCS(iSublPhaseIndex,ibbyy,2)

                dChemicalPotential(j) = ((dQx * dChemicalPotentialTemp(iax + iFirst - 1) / (dZa * dZx)) &
                      + (dQx * dChemicalPotentialTemp(ibx + iFirst - 1) / (dZb * dZx)) &
                      + (dQy * dChemicalPotentialTemp(iay + iFirst - 1) / (dZa * dZy)) &
                      + (dQy * dChemicalPotentialTemp(iby + iFirst - 1) / (dZb * dZy))) &
                      / ((dQx/dZx) + (dQy/dZy))
                ! dChemicalPotential(j) = ((dQx * dZaxa * dChemicalPotentialTemp(iax + iFirst - 1) / (dZa * dZx)) &
                !       + (dQx * dZbxb * dChemicalPotentialTemp(ibx + iFirst - 1) / (dZb * dZx)) &
                !       + (dQy * dZaya * dChemicalPotentialTemp(iay + iFirst - 1) / (dZa * dZy)) &
                !       + (dQy * dZbyb * dChemicalPotentialTemp(iby + iFirst - 1) / (dZb * dZy))) &
                !       / (2 * ((dQx/dZx) + (dQy/dZy)))
            end do LOOP_nSUBGQCS
        else
            LOOP_nSpeciesCS: do i = nSpeciesPhaseCS(n - 1) + 1, nSpeciesPhaseCS(n)
                l = 0
                ! Loop through the Gibbs energy equations to figure out which one to use:
                do k = 1, nGibbsEqSpecies(i)
                    iCounterGibbsEqn = iCounterGibbsEqn + 1
                    if ((dTemperature <= dGibbsCoeffSpeciesTemp(1,iCounterGibbsEqn)).AND.(l == 0)) then
                        l = k
                    end if
                end do

                ! This species will not be considered part of the system.
                if (iSpeciesPass(i) == 0) cycle LOOP_nSpeciesCS

                if (l == 0) l = nGibbsEqSpecies(i)

                l = l + iCounterGibbsEqn - nGibbsEqSpecies(i)
                j = j + 1   ! New species index

                cSpeciesName(j)            = cSpeciesNameCS(i)
                iPhase(j)                  = iPhaseCS(i)
                iParticlesPerMole(j)       = iParticlesPerMoleCS(i)
                dCoeffGibbsMagnetic(j,1:4) = dGibbsMagneticCS(i,1:4)

                m = 0
                do k = 1, nElemOrComp
                    if (iElementSystem(k) /= 0) then
                        m = m + 1
                        dStoichSpecies(j,m) = dStoichSpeciesCS(i,k)
                    end if
                end do

                do k = 2, 7
                    dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(k,l) * dGibbsCoeff(k-1)
                end do

                ! Compute additional standard molar Gibbs energy terms:
                if (dGibbsCoeffSpeciesTemp(9,l) .EQ. 99) then
                    dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(8,l) * dLogT
                else
                    dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(8,l) &
                        *dTemperature**dGibbsCoeffSpeciesTemp(9,l)
                end if

                if (dGibbsCoeffSpeciesTemp(11,l) .EQ. 99) then
                    dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(10,l) * dLogT
                else
                    dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(10,l) &
                        * dTemperature**dGibbsCoeffSpeciesTemp(11,l)
                end if

                if (dGibbsCoeffSpeciesTemp(13,l) .EQ. 99) then
                    dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(12,l) * dLogT
                else
                    dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(12,l) &
                        * dTemperature**dGibbsCoeffSpeciesTemp(13,l)
                end if

                ! Compute the magnetic terms (if applicable):
                if ((dGibbsMagneticCS(i,1) /= 0D0).AND.(iPhase(j) == 0)) then
                    call CompGibbsMagnetic(i,j)
                end if

                ! Convert chemical potentials to dimensionless units:
                dChemicalPotential(j) = dChemicalPotential(j) * dTemp * DFLOAT(iParticlesPerMoleCS(i))

                ! Add pressure dependence term to the chemical potential term:
                if (iPhaseCS(i) == 1) then
                    if (cSolnPhaseNameCS(1) == 'gas_ideal') then
                        ! Note: If an ideal gas is included in a ChemSage data-file, then it is always
                        ! the first solution phase in the data-file.
                        dChemicalPotential(j) = dChemicalPotential(j) + dLogP + (0.10945D0 / dIdealConstant)
                    end if
                else if (iPhaseCS(i) == -1) then
                    ! Explicitly set dummy species chemical potentials
                    dChemicalPotential(j) = 100000
                end if

            end do LOOP_nSpeciesCS ! End loop of species (i)
        end if
    end do LOOP_nPhasesCS

    LOOP_nPureConSpeciesCS: do i = nSpeciesPhaseCS(nSolnPhasesSysCS) + 1, nSpeciesCS
        l = 0
        ! Loop through the Gibbs energy equations to figure out which one to use:
        do k = 1, nGibbsEqSpecies(i)
            iCounterGibbsEqn = iCounterGibbsEqn + 1
            if ((dTemperature <= dGibbsCoeffSpeciesTemp(1,iCounterGibbsEqn)).AND.(l == 0)) then
                l = k
            end if
        end do

        ! This species will not be considered part of the system.
        if (iSpeciesPass(i) == 0) cycle LOOP_nPureConSpeciesCS

        if (l == 0) l = nGibbsEqSpecies(i)

        l = l + iCounterGibbsEqn - nGibbsEqSpecies(i)
        j = j + 1   ! New species index

        cSpeciesName(j)            = cSpeciesNameCS(i)
        iPhase(j)                  = iPhaseCS(i)
        iParticlesPerMole(j)       = iParticlesPerMoleCS(i)
        dCoeffGibbsMagnetic(j,1:4) = dGibbsMagneticCS(i,1:4)

        m = 0
        do k = 1, nElemOrComp
            if (iElementSystem(k) /= 0) then
                m = m + 1
                dStoichSpecies(j,m) = dStoichSpeciesCS(i,k)
            end if
        end do

        do k = 2, 7
            dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(k,l) * dGibbsCoeff(k-1)
        end do

        ! Compute additional standard molar Gibbs energy terms:
        if (dGibbsCoeffSpeciesTemp(9,l) .EQ. 99) then
            dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(8,l) * dLogT
        else
            dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(8,l) &
                *dTemperature**dGibbsCoeffSpeciesTemp(9,l)
        end if

        if (dGibbsCoeffSpeciesTemp(11,l) .EQ. 99) then
            dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(10,l) * dLogT
        else
            dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(10,l) &
                * dTemperature**dGibbsCoeffSpeciesTemp(11,l)
        end if

        if (dGibbsCoeffSpeciesTemp(13,l) .EQ. 99) then
            dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(12,l) * dLogT
        else
            dChemicalPotential(j) = dChemicalPotential(j) + dGibbsCoeffSpeciesTemp(12,l) &
                * dTemperature**dGibbsCoeffSpeciesTemp(13,l)
        end if

        ! Compute the magnetic terms (if applicable):
        if ((dGibbsMagneticCS(i,1) /= 0D0).AND.(iPhase(j) == 0)) then
            call CompGibbsMagnetic(i,j)
        end if

        ! Convert chemical potentials to dimensionless units:
        dChemicalPotential(j) = dChemicalPotential(j) * dTemp * DFLOAT(iParticlesPerMoleCS(i))

        ! Add pressure dependence term to the chemical potential term:
        if (iPhaseCS(i) == 1) then
            if (cSolnPhaseNameCS(1) == 'gas_ideal') then
                ! Note: If an ideal gas is included in a ChemSage data-file, then it is always
                ! the first solution phase in the data-file.
                dChemicalPotential(j) = dChemicalPotential(j) + dLogP + (0.10945D0 / dIdealConstant)
            end if
        else if (iPhaseCS(i) == -1) then
            ! Explicitly set dummy species chemical potentials
            dChemicalPotential(j) = 100000
        end if
    end do LOOP_nPureConSpeciesCS ! End loop of species (i)

    ! Update all the excess properties:
    n = 0
    LOOP_SolnPhases: do i = 1, nSolnPhasesSysCS

        if ((cSolnPhaseTypeCS(i) == 'SUBL').OR.(cSolnPhaseTypeCS(i) == 'SUBLM')) nCounter = nCounter + 1

        LOOP_Param: do j = nParamPhaseCS(i-1) + 1, nParamPhaseCS(i)

            ! Proceed if the parameter passed:
            IF_ParamPass: if (iParamPassCS(j) /= 0) then
                n = n + 1

                ! Convert the excess Gibbs energy parameters to dimensionless units:
                iRegularParam(n,1:nParamMax*2+1) = iRegularParamCS(j,1:nParamMax*2+1)

                select case (cSolnPhaseTypeCS(i))
                    case ('QKTO', 'RKMP', 'RKMPM')

                        do k = 1, 6
                            dExcessGibbsParam(n) = dExcessGibbsParam(n) + dRegularParamCS(j,k) * dGibbsCoeff(k)
                        end do
                        dExcessGibbsParam(n) = dExcessGibbsParam(n) * dTemp

                        ! Loop through species involved in mixing parameter:
                        do k = 1, iRegularParamCS(j,1)
                            m                    = iRegularParamCS(j,k+1) + nSpeciesPhaseCS(i-1)
                            iRegularParam(n,k+1) = iSpeciesPass(m)
                        end do

                    case ('SUBG','SUBQ')

                        ! Must remove unused elements from iRegularParam
                        iSublPhaseIndex = iPhaseSublatticeCS(i)
                        nRemove = 0
                        iRemove = 0
                        do k = nSublatticePhaseCS(iSublPhaseIndex), 1, -1
                            do l = nSublatticeElementsCS(iSublPhaseIndex,k), 1, -1
                                if (iSublatticeElementsCS(iSublPhaseIndex,k,l) <= 0) then
                                    nRemove = nRemove + 1
                                    iRemove(nRemove) = l + ((k - 1) * nSublatticeElementsCS(iSublPhaseIndex,1))
                                elseif (iElementSystem(iSublatticeElementsCS(iSublPhaseIndex,k,l)) == 0) then
                                    nRemove = nRemove + 1
                                    iRemove(nRemove) = l + ((k - 1) * nSublatticeElementsCS(iSublPhaseIndex,1))
                                end if
                            end do
                        end do

                        do k = 1, nRemove
                            do l = 2, 5
                                if (iRegularParam(n,l) > iRemove(k)) then
                                    iRegularParam(n,l) = iRegularParam(n,l) - 1
                                end if
                            end do
                        end do

                        ! Note that this is different for SUBG phases than QKTO, RKMP, or SUBL phases:
                        do k = 1, 4
                            dExcessGibbsParam(n) = dExcessGibbsParam(n) + dRegularParamCS(j,k+2) * dGibbsCoeff(k)
                        end do
                        dExcessGibbsParam(n) = dExcessGibbsParam(n) * dTemp

                    case ('SUBL', 'SUBLM')

                        do k = 1, 6
                            dExcessGibbsParam(n) = dExcessGibbsParam(n) + dRegularParamCS(j,k) * dGibbsCoeff(k)
                        end do
                        dExcessGibbsParam(n) = dExcessGibbsParam(n) * dTemp

                        ! Loop through constituents involved in mixing parameter:
                        do k = 1, iRegularParamCS(j,1)

                            ! The constituent numbering scheme from ChemSage does not consider the sublattice #, but just
                            ! a continuing count of the constituents.
                            m = iRegularParamCS(j,k+1)

                            ! Figure out the sublattice (l) and constituent (m) indices:
                            LOOP_SUBL: do s = 1, nSublatticePhaseCS(nCounter)
                                l = s
                                if (m > nConstituentSublatticeCS(nCounter,s)) then
                                    m = m - nConstituentSublatticeCS(nCounter,s)
                                else
                                    exit LOOP_SUBL
                                end if

                            end do LOOP_SUBL

                            ! Apply indexing scheme (l is sublattice index, iCounstituentPass is constituent index on
                            ! sublattice l):
                            iRegularParam(n,k+1) = (10000 * l) + iConstituentPass(nCounter,s,m)

                        end do

                        ! Loop through constituents involved in mixing parameter to see if they need to be shuffled.
                        ! ChemSage files do not order the constituents based on which ones mix.  For instance, there
                        ! may be three constituents mixing where the first constituent is on the first sublattice and
                        ! the second and third are on the second sublattice.
                        LOOP_SUBL_Check: do k = 2, iRegularParamCS(j,1)

                            l = MOD(iRegularParam(n,k), 10000)
                            l = (iRegularParam(n,k) - l) / 10000

                            m = MOD(iRegularParam(n,k+1), 10000)
                            m = (iRegularParam(n,k+1) - m) / 10000

                            ! If these two constituents are on the same sublattice and they correspond to the
                            ! first two mixing constituents, then exit:
                            if (l == m) then

                                ! Check the order of constituents:
                                if (k == 2) then
                                    ! These constituent indices are correctly placed:
                                    exit LOOP_SUBL_Check
                                elseif (k == 3) then
                                    ! Shuffle the vector:
                                    l = iRegularParam(n,2)
                                    iRegularParam(n,2) = iRegularParam(n,3)
                                    iRegularParam(n,3) = iRegularParam(n,4)
                                    iRegularParam(n,4) = l

                                    exit LOOP_SUBL_Check
                                else
                                    ! Report an error and exit:
                                    INFOThermo = 36
                                    exit LOOP_SUBL_Check
                                end if
                            end if
                        end do LOOP_SUBL_Check
                end select
            end if IF_ParamPass
        end do LOOP_Param
    end do LOOP_SolnPhases

    ! Update the phase index vector (iPhase):
    do i = 1, nSolnPhasesSys
        j = nSpeciesPhase(i-1) + 1
        k = nSpeciesPhase(i)
        iPhase(j:k) = i
    end do

    ! If necessary, add dummy species for each ionic phase:
    if (nChargedConstraints > 0) then
        j = 0

        do i = nSpecies, nSpecies - nChargedConstraints + 1, -1
            j = j + 1
            cSpeciesName(i)       = 'e-'
            iPhase(i)             = -1
            dChemicalPotential(i) = 5D6 / (dIdealConstant * dTemperature)
            iParticlesPerMole(i)  = 1
            dStoichSpecies(i,nElements - j + 1) = 1D0
        end do
    end if

    i = nSpeciesPhase(nSolnPhasesSys) + 1
    nDummySpecies = ABS(SUM(iPhase(i:nSpecies)))

    ! Compute the total number of atoms per formula mass of each species:
    dSpeciesTotalAtoms = SUM(ABS(dStoichSpecies),DIM=2)

    ! Compute the atomic fraction of each element in a species:
    do i = 1, nSpecies

        ! SUBG phases will be allocated with more species than appear in the data-file
        ! because it'll work with pair fractions rather than species. Therefore, dummy
        ! entries will be created with no stoichiometric values. To avoid dividing by
        ! zero, the following is included:
        !dSpeciesTotalAtoms(i) = DMAX1(dSpeciesTotalAtoms(i),0.1D0)
        if (dSpeciesTotalAtoms(i) == 0) then
            dSpeciesTotalAtoms(i) = 0.00001D0
        end if

        ! Convert chemical potentials from [J/mol] to [J/g-at]:
        dTemp                 = 1D0 / dSpeciesTotalAtoms(i)
        dChemicalPotential(i) = dChemicalPotential(i) * dTemp

        do j = 1, nElements
            dAtomFractionSpecies(i,j) = dStoichSpecies(i,j) * dTemp
        end do
    end do

    ! Store the standard molar Gibbs energies:
    do i = 1, nSpecies
        dStdGibbsEnergy(i) = dChemicalPotential(i) * dSpeciesTotalAtoms(i) / DFLOAT(iParticlesPerMole(i))
    end do

    ! Store an integer vector representing the component index when the phase is ionic.
    if (allocated(iPhaseElectronID)) deallocate(iPhaseElectronID)
    allocate(iPhaseElectronID(nSolnPhasesSys))
    iPhaseElectronID = 0

    ! Loop through all solution phases:
    LOOP_EID_Phase: do i = 1, nSolnPhasesSys
        ! Check if this phase may be ionic:
        if ((cSolnPhaseType(i) == 'SUBL').OR.(cSolnPhaseType(i) == 'SUBLM')) then
            ! Loop through species in phase:
            LOOP_EID_Species: do j = nSpeciesPhase(i-1) + 1, nSpeciesPhase(i)
                ! Loop through system components represented by electrons:
                do k = nElements - nChargedConstraints + 1, nElements
                    ! If a species in this phase is represented by this electron, store in the index and cycle:

                    if (dStoichSpecies(j,k) /= 0D0) then
                        iPhaseElectronID(i) = k
                        cycle LOOP_EID_Phase
                    end if
                end do
            end do LOOP_EID_Species
        end if
    end do LOOP_EID_Phase

    return

end subroutine CompThermoData
