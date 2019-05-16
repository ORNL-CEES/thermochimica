
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompExcessGibbsEnergySUBG.f90
    !> \brief   Compute the partial molar excess Gibbs energy of mixing of solution phase species in a SUBG
    !!           phase.
    !> \author  M.H.A. Piro
    !> \date    December 10, 2018
    !> \sa      CompExcessGibbsEnergy.f90
    !> \sa      CompExcessGibbsEnergyRKMP.f90
    !> \sa      CompExcessGibbsEnergyQKTO.f90
    !> \sa      CompExcessGibbsEnergySUBL.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/01/2018      M.H.A. Piro         Original code.
    !   12/10/2018      M.H.A. Piro         Fixed a bug in the partial molar excess Gibbs energy
    !                                        expression for BB and AB.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the chemical potentials of pairs of species
    !! (short range order) in a non-ideal solution phase designated as 'SUBG', which is a modified
    !! quasichemical model (MQM). A unique characteristic of the MQM model is that the focus is not
    !! placed on the species (aka 'compound end members'), but rather the pairs of nearest neigbours.
    !! For example, if one had a binary solution phase A-B, this model consideres A-A, B-B, and A-B
    !! as pairs of species, which are distributed about a quasi-lattice. Since the focus is on pairs
    !! of species rather than the species themselves, this considerably changes the calculation in
    !! comparison to other models (e.g., QKTO, RBMK, SUBL).
    !!
    !! For more information on the SUBG model and the derivation of equations, the reader is referred
    !! to the following literature:
    !!
    !!      A.D. Pelton, S.A. Degterov, G. Eriksson, C. Robelin, Y. Dessureault, ``The Modified
    !!       Quasichemical Model I -- Binary Solutions'', Metallurgical and Materials Transactions B,
    !!       31B (2000) 651-659.
    !!
    !!      A.D. Pelton, P. Chartrand, ``The Modified Quasi-Chemical Model: Part II. Multicomponent
    !!       Solutions'', Metallurgical and Materials Transactions B, 32A (2001) 1355-1360.
    !!
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in] iSolnIndex    Absolute index of a solution phase
    !
    ! nSpeciesPhase             An integer vector representing the last index of a species in a particular
    !                            solution phase
    ! iPairID(:,:)              An integer array representing the indices of pairs.
    !
    ! dCoordinationNumber(:,:)  A double array representing the coordination number of pairs.
    ! dX(:)                     A temporary double real vector used to represent the mole fractions of the
    !                            species.
    ! dY(:)                     A temporary double vector used to represent the coordinate equivalent
    !                            fractions of the species.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CompExcessGibbsEnergySUBG(iSolnIndex)

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, l, m, p, q, r, ii, jj, kk ,ll, ka, la
    integer :: iSolnIndex, iSublPhaseIndex, nPhaseElements
    integer :: iFirst, iLast, nA, nX, iWeight
    real(8) :: dTemp, dSum, dEntropy, dPowXij, dPowYi
    real(8) :: dZAAA, dZABA, dZBAB
    real(8) :: x, y, z
    real(8), allocatable, dimension(:) :: dXi, dYi, dNi
    real(8), allocatable, dimension(:,:) :: dXij
    ! X_ij/kl corresponds to dMolFracion


    ! Only proceed if the correct phase type is selected:
    IF_SUBG: if (cSolnPhaseType(iSolnIndex) == 'SUBG' .OR. cSolnPhaseType(iSolnIndex) == 'SUBQ') then

        ! Define temporary variables for sake of convenience:
        iFirst = nSpeciesPhase(iSolnIndex-1) + 1
        iLast  = nSpeciesPhase(iSolnIndex)
        iSublPhaseIndex = iPhaseSublattice(iSolnIndex)

        ! Allocate allocatable arrays:
        if (allocated(dXi)) deallocate(dXi)
        if (allocated(dYi)) deallocate(dYi)
        if (allocated(dNi)) deallocate(dNi)
        if (allocated(dXij)) deallocate(dXij)
        j = iLast - iFirst + 1
        nPhaseElements = nSublatticeElements(iSublPhaseIndex,1) + nSublatticeElements(iSublPhaseIndex,2)
        allocate(dXi(nPhaseElements),dYi(nPhaseElements),dNi(nPhaseElements))
        allocate(dXij(nSublatticeElements(iSublPhaseIndex,1),nSublatticeElements(iSublPhaseIndex,2)))

        ! Initialize variables:
        dXi                               = 0D0
        dYi                               = 0D0
        dXij                              = 0D0
        dChemicalPotential(iFirst:iLast)  = 0D0
        dPartialExcessGibbs(iFirst:iLast) = 0D0

        ! Compute X_i and Y_i
        ! Do cations first:
        dSum = 0D0
        do i = 1, nSublatticeElements(iSublPhaseIndex,1)
            do k = 1, nPairsSRO(iSublPhaseIndex,2)
                l = iFirst + k - 1
                if (i == iPairID(iSublPhaseIndex,k,1))  then
                    dNi(i) = dNi(i) + (dMolFraction(l) / dCoordinationNumber(iSublPhaseIndex,k,1))
                    dYi(i) = dYi(i) + (dMolFraction(l) / 2)
                end if
                if (i == iPairID(iSublPhaseIndex,k,2))  then
                    dNi(i) = dNi(i) + (dMolFraction(l) / dCoordinationNumber(iSublPhaseIndex,k,2))
                    dYi(i) = dYi(i) + (dMolFraction(l) / 2)
                end if
            end do
            dSum = dSum + dNi(i)
        end do
        do i = 1, nSublatticeElements(iSublPhaseIndex,1)
            dXi(i) = dNi(i) / dSum
        end do
        ! Do anions now:
        dSum = 0D0
        do i = 1, nSublatticeElements(iSublPhaseIndex,2)
            j = i + nSublatticeElements(iSublPhaseIndex,1)
            do k = 1, nPairsSRO(iSublPhaseIndex,2)
                l = iFirst + k - 1
                if (j == iPairID(iSublPhaseIndex,k,3))  then
                    dNi(j) = dNi(j) + (dMolFraction(l) / dCoordinationNumber(iSublPhaseIndex,k,3))
                    dYi(j) = dYi(j) + (dMolFraction(l) / 2)
                end if
                if (j == iPairID(iSublPhaseIndex,k,4))  then
                    dNi(j) = dNi(j) + (dMolFraction(l) / dCoordinationNumber(iSublPhaseIndex,k,4))
                    dYi(j) = dYi(j) + (dMolFraction(l) / 2)
                end if
            end do
            dSum = dSum + dNi(j)
        end do
        do i = 1, nSublatticeElements(iSublPhaseIndex,2)
            j = i + nSublatticeElements(iSublPhaseIndex,1)
            dXi(j) = dNi(j) / dSum
        end do

        ! Compute X_i/j
        do i = 1, nSublatticeElements(iSublPhaseIndex,1)
            do j = 1, nSublatticeElements(iSublPhaseIndex,2)
                do k = 1, nPairsSRO(iSublPhaseIndex,2)
                    l = iFirst + k - 1
                    nA = 0
                    if (i == iPairID(iSublPhaseIndex,k,1))  then
                        nA = nA + 1
                    end if
                    if (i == iPairID(iSublPhaseIndex,k,2))  then
                        nA = nA + 1
                    end if
                    nX = 0
                    if ((j + nSublatticeElements(iSublPhaseIndex,1)) == iPairID(iSublPhaseIndex,k,3))  then
                        nX = nX + 1
                    end if
                    if ((j + nSublatticeElements(iSublPhaseIndex,1)) == iPairID(iSublPhaseIndex,k,4))  then
                        nX = nX + 1
                    end if
                    dXij(i,j) = dXij(i,j) + (dMolFraction(l) * nA * nX / 4)
                end do
            end do
        end do

        ! ---------------------------------------------------------------
        ! COMPUTE REFERENCE GIBBS ENERGY AND IDEAL MIXING TERMS
        ! ---------------------------------------------------------------

        do k = 1, nPairsSRO(iSublPhaseIndex,2)
            ! Calculate entropic contributionrs to chemical potentials
            dEntropy = 0D0
            l = iFirst + k - 1

            ! Loop over n_i contributions to entropy
            ! Cations first
            do i = 1, nSublatticeElements(iSublPhaseIndex,1)
                if (i == iPairID(iSublPhaseIndex,k,1))  then
                    dEntropy = dEntropy + (DLOG(dXi(i)) / dCoordinationNumber(iSublPhaseIndex,k,1))
                end if
                if (i == iPairID(iSublPhaseIndex,k,2))  then
                    dEntropy = dEntropy + (DLOG(dXi(i)) / dCoordinationNumber(iSublPhaseIndex,k,2))
                end if
            end do
            ! Now anions
            do i = 1, nSublatticeElements(iSublPhaseIndex,2)
                j = i + nSublatticeElements(iSublPhaseIndex,1)
                if (j == iPairID(iSublPhaseIndex,k,3))  then
                    dEntropy = dEntropy + (DLOG(dXi(j)) / dCoordinationNumber(iSublPhaseIndex,k,3))
                end if
                if (j == iPairID(iSublPhaseIndex,k,4))  then
                    dEntropy = dEntropy + (DLOG(dXi(j)) / dCoordinationNumber(iSublPhaseIndex,k,4))
                end if
            end do

            ! Loop over n_i/j contributions to entropy
            m = 0
            do i = 1, nSublatticeElements(iSublPhaseIndex,1)
                do j = 1, nSublatticeElements(iSublPhaseIndex,2)
                    m = m + 1
                    nA = 0
                    if (i == iPairID(iSublPhaseIndex,k,1))  then
                        nA = nA + 1
                    end if
                    if (i == iPairID(iSublPhaseIndex,k,2))  then
                        nA = nA + 1
                    end if
                    nX = 0
                    if ((j + nSublatticeElements(iSublPhaseIndex,1)) == iPairID(iSublPhaseIndex,k,3))  then
                        nX = nX + 1
                    end if
                    if ((j + nSublatticeElements(iSublPhaseIndex,1)) == iPairID(iSublPhaseIndex,k,4))  then
                        nX = nX + 1
                    end if
                    dEntropy = dEntropy + (DLOG(dXij(i,j) / (dYi(i) * dYi(j + nSublatticeElements(iSublPhaseIndex,1)))) &
                                          * (nA * nX / dZetaSpecies(iSublPhaseIndex,m)))
                end do
            end do

            ! pair indices:
            ii = iPairID(iSublPhaseIndex,k,1)
            jj = iPairID(iSublPhaseIndex,k,2)
            kk = iPairID(iSublPhaseIndex,k,3)
            ll = iPairID(iSublPhaseIndex,k,4)
            ka = kk - nSublatticeElements(iSublPhaseIndex,1)
            la = ll - nSublatticeElements(iSublPhaseIndex,1)

            ! Add n_ij/kl contribution
            iWeight = 1
            if (ii /= jj) iWeight = iWeight * 2
            if (kk /= ll) iWeight = iWeight * 2

            ! SUBG and SUBQ differ in entropy calculation by the powers to which X_i/j and Y_i are raised
            if (cSolnPhaseType(iSolnIndex) == 'SUBG') then
                dPowXij = 1D0
                dPowYi  = 1D0
            else if (cSolnPhaseType(iSolnIndex) == 'SUBQ') then
                dPowXij = 0.75D0
                dPowYi  = 0.5D0
            end if
            dEntropy = dEntropy + DLOG(dMolFraction(l) / (iWeight * (dXij(ii,ka)**dPowXij) * (dXij(ii,la)**dPowXij) &
                                                                  * (dXij(jj,ka)**dPowXij) * (dXij(jj,la)**dPowXij) &
                                            / ((dYi(ii)**dPowYi) * (dYi(jj)**dPowYi) * (dYi(kk)**dPowYi) * (dYi(ll)**dPowYi))))

            dChemicalPotential(l) = dEntropy

        end do

        ! Loop through all pairs:
        LOOP_C: do i = 1, nPairsSRO(iSublPhaseIndex,2)
            j = iFirst + i - 1
            ! AA pairs:
            if (iPairID(iSublPhaseIndex,i,1) == iPairID(iSublPhaseIndex,i,2)) then
                ! Store coordination numbers:
                dZAAA = dCoordinationNumber(iSublPhaseIndex,i,1)

                ! Compute standard reference Gibbs energy (Eq [15]):
                ! dChemicalPotential(j) = dStdGibbsEnergy(j) * 2D0 / dZAAA

                ! Compute ideal mixing component:
                ! dChemicalPotential(j) = dChemicalPotential(j) + (2D0 / dZAAA) * DLOG(dXi(i)) + DLOG(dMolFraction(j) / dYi(i)**2)
                dChemicalPotential(j) = (2D0 / dZAAA) * DLOG(dXi(i)) + DLOG(dMolFraction(j) / dYi(i)**2)
            ! AB pairs:
            else
                k = iPairID(iSublPhaseIndex,i,1)            ! Index of AA
                l = iPairID(iSublPhaseIndex,i,2)            ! Index of BB

                ! Store coordination numbers:
                dZABA = dCoordinationNumber(iSublPhaseIndex,i,1)
                dZBAB = dCoordinationNumber(iSublPhaseIndex,i,2)

                ! Compute standard reference Gibbs energy:
                ! NOTE: I did not include $\Delta g_{AB}^{\circ}$ in this equation because I think it makes more sense
                ! to include it in the excess mixing section.
                ! Eq [16]:
                ! dChemicalPotential(j) = dStdGibbsEnergy(k + iFirst - 1) / dZABA + dStdGibbsEnergy(l + iFirst - 1) / dZBAB

                ! Compute ideal mixing component:
                ! dChemicalPotential(j) = dChemicalPotential(j) + DLOG(dMolFraction(j) / (2D0 * dYi(k) * dYi(l))) &
                !                                               + DLOG(dXi(k)) / dZABA + DLOG(dXi(l)) / dZBAB
                dChemicalPotential(j) = DLOG(dMolFraction(j) / (2D0 * dYi(k) * dYi(l))) &
                                        + DLOG(dXi(k)) / dZABA + DLOG(dXi(l)) / dZBAB
            end if

        end do LOOP_C

        ! ---------------------------------------------------------
        ! COMPUTE PARTIAL MOLAR EXCESS GIBBS ENERGY OF MIXING TERMS
        ! ---------------------------------------------------------

        ! Loop through excess mixing parameters:
        LOOP_Param: do m = nParamPhase(iSolnIndex-1) + 1, nParamPhase(iSolnIndex)
            i = iRegularParam(m,2)              ! Index of AA
            j = iRegularParam(m,3)              ! Index of BB
            k = 0
            ! Find which (if any) position AB is stored at:
            LOOP_FindPair: do l = 1, nPairsSRO(iSublPhaseIndex,2)
                if (((iPairID(iSublPhaseIndex,l,1) == i) .AND. (iPairID(iSublPhaseIndex,l,2) == j)) .OR. &
                    ((iPairID(iSublPhaseIndex,l,1) == j) .AND. (iPairID(iSublPhaseIndex,l,2) == i)))  then
                    k = l
                    EXIT LOOP_FindPair
                end if
            end do LOOP_FindPair
            p = iRegularParam(m,6)              ! Exponent of AA
            q = iRegularParam(m,7)              ! Exponent of BB
            r = iRegularParam(m,8)              ! Exponent of AB
            x = dMolFraction(iFirst + i - 1)    ! x_AA
            y = dMolFraction(iFirst + j - 1)    ! x_BB
            if (k > 0) then
                z = dMolFraction(iFirst + k - 1)! x_AB
            else
                z = 0
            end if
            dSum = x + y + z

            ! Contribution to AA:
            dTemp = z * x**(p-1) * y**(q) * (DBLE(p)*(y+z) - DBLE(q)*x)
            dTemp = dTemp * dExcessGibbsParam(m) / 2D0
            dPartialExcessGibbs(iFirst + i - 1) = dPartialExcessGibbs(iFirst + i - 1) + dTemp

            ! Contribution to BB:
            dTemp = z * x**(p) * y**(q-1) * (DBLE(q)*(x+z) - DBLE(p)*y)
            dTemp = dTemp * dExcessGibbsParam(m) / 2D0
            dPartialExcessGibbs(iFirst + j - 1) = dPartialExcessGibbs(iFirst + j - 1) + dTemp

            ! Contribution to AB (only if pair exists):
            if (k > 0) then
                dTemp = x**(p) * y**(q) * (z*(1D0 - DBLE(p) - DBLE(q)) + x + y)
                dTemp = dTemp * dExcessGibbsParam(m) / 2D0
                dPartialExcessGibbs(iFirst + k - 1) = dPartialExcessGibbs(iFirst + k - 1) + dTemp
            end if

        end do LOOP_Param

    end if IF_SUBG

    ! Deallocate allocatable arrays:
    deallocate(dXi,dYi,dNi,dXij)

    return

end subroutine CompExcessGibbsEnergySUBG
