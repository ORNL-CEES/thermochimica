
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

    integer :: i, j, k, l, m, p, q, r, s, ii, jj, kk ,ll, ka, la
    integer :: a, b, c, x, y, z, e, f, ijkl, abxy, exp, iexp
    integer :: iSolnIndex, iSublPhaseIndex, nPhaseElements
    integer :: iFirst, iLast, nA, nX, iWeight, iBlock, iQuad, iQuad2
    integer :: iA2X2, iB2X2, iA2Y2, iB2Y2
    real(8) :: dSum, dEntropy, dRef, dPowXij, dPowYi
    real(8) :: dZa, dZb, dZx, dZy, dGex, dDgex, dDgexBase
    real(8) :: dXA2X2, dXB2X2, dXA2Y2, dXB2Y2, dXexp
    real(8), allocatable, dimension(:) :: dXi, dYi, dNi
    real(8), allocatable, dimension(:,:) :: dXij, dNij
    ! X_ij/kl corresponds to dMolFracion


    ! Only proceed if the correct phase type is selected:
    if (.NOT. (cSolnPhaseType(iSolnIndex) == 'SUBG' .OR. cSolnPhaseType(iSolnIndex) == 'SUBQ')) return

    ! Define temporary variables for sake of convenience:
    iFirst = nSpeciesPhase(iSolnIndex-1) + 1
    iLast  = nSpeciesPhase(iSolnIndex)
    iSublPhaseIndex = iPhaseSublattice(iSolnIndex)

    ! Allocate allocatable arrays:
    if (allocated(dXi)) deallocate(dXi)
    if (allocated(dYi)) deallocate(dYi)
    if (allocated(dNi)) deallocate(dNi)
    if (allocated(dXij)) deallocate(dXij)
    if (allocated(dNij)) deallocate(dNij)
    j = iLast - iFirst + 1
    nPhaseElements = nSublatticeElements(iSublPhaseIndex,1) + nSublatticeElements(iSublPhaseIndex,2)
    allocate(dXi(nPhaseElements),dYi(nPhaseElements),dNi(nPhaseElements))
    allocate(dXij(nSublatticeElements(iSublPhaseIndex,1),nSublatticeElements(iSublPhaseIndex,2)))
    allocate(dNij(nSublatticeElements(iSublPhaseIndex,1),nSublatticeElements(iSublPhaseIndex,2)))

    ! Initialize variables:
    dXi                               = 0D0
    dYi                               = 0D0
    dNi                               = 0D0
    dXij                              = 0D0
    dNij                              = 0D0
    dChemicalPotential(iFirst:iLast)  = 0D0
    dPartialExcessGibbs(iFirst:iLast) = 0D0

    ! Compute X_i and Y_i
    ! Do cations first:
    dSum = 0D0
    do i = 1, nSublatticeElements(iSublPhaseIndex,1)
        do k = 1, nPairsSRO(iSublPhaseIndex,2)
            l = iFirst + k - 1
            dZa = dCoordinationNumber(iSublPhaseIndex,k,1)
            dZb = dCoordinationNumber(iSublPhaseIndex,k,2)
            if (i == iPairID(iSublPhaseIndex,k,1))  then
                dNi(i) = dNi(i) + (dMolFraction(l) / dZa)
                dYi(i) = dYi(i) + (dMolFraction(l) / 2)
            end if
            if (i == iPairID(iSublPhaseIndex,k,2))  then
                dNi(i) = dNi(i) + (dMolFraction(l) / dZb)
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
            dZx = dCoordinationNumber(iSublPhaseIndex,k,3)
            dZy = dCoordinationNumber(iSublPhaseIndex,k,4)
            if (j == iPairID(iSublPhaseIndex,k,3))  then
                dNi(j) = dNi(j) + (dMolFraction(l) / dZx)
                dYi(j) = dYi(j) + (dMolFraction(l) / 2)
            end if
            if (j == iPairID(iSublPhaseIndex,k,4))  then
                dNi(j) = dNi(j) + (dMolFraction(l) / dZy)
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
    dSum = 0D0
    do i = 1, nSublatticeElements(iSublPhaseIndex,1)
        do j = 1, nSublatticeElements(iSublPhaseIndex,2)
            m = iConstituentSublattice(iSublPhaseIndex,1,i) + &
            ((iConstituentSublattice(iSublPhaseIndex,2,j) - 1) * nSublatticeElements(iSublPhaseIndex,1))
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
                dNij(i,j) = dNij(i,j) + (dMolFraction(l) * nA * nX / dZetaSpecies(iSublPhaseIndex,m))
            end do
            dSum = dSum + dNij(i,j)
        end do
    end do

    do i = 1, nSublatticeElements(iSublPhaseIndex,1)
        do j = 1, nSublatticeElements(iSublPhaseIndex,2)
            dXij(i,j) = dNij(i,j) / dSum
        end do
    end do

    ! ---------------------------------------------------------------
    ! COMPUTE REFERENCE GIBBS ENERGY AND IDEAL MIXING TERMS
    ! ---------------------------------------------------------------

    do k = 1, nPairsSRO(iSublPhaseIndex,2)
        ! Calculate entropic contributions to chemical potentials
        dEntropy = 0D0
        dRef = 0D0
        l = iFirst + k - 1

        ! Coordination numbers for this quadruplet:
        dZa = dCoordinationNumber(iSublPhaseIndex,k,1)
        dZb = dCoordinationNumber(iSublPhaseIndex,k,2)
        dZx = dCoordinationNumber(iSublPhaseIndex,k,3)
        dZy = dCoordinationNumber(iSublPhaseIndex,k,4)

        ! Loop over n_i contributions to entropy
        ! Cations first
        do i = 1, nSublatticeElements(iSublPhaseIndex,1)
            if (i == iPairID(iSublPhaseIndex,k,1))  then
                dEntropy = dEntropy + (DLOG(dXi(i)) / dZa)
            end if
            if (i == iPairID(iSublPhaseIndex,k,2))  then
                dEntropy = dEntropy + (DLOG(dXi(i)) / dZb)
            end if
        end do
        ! Now anions
        do i = 1, nSublatticeElements(iSublPhaseIndex,2)
            j = i + nSublatticeElements(iSublPhaseIndex,1)
            if (j == iPairID(iSublPhaseIndex,k,3))  then
                dEntropy = dEntropy + (DLOG(dXi(j)) / dZx)
            end if
            if (j == iPairID(iSublPhaseIndex,k,4))  then
                dEntropy = dEntropy + (DLOG(dXi(j)) / dZy)
            end if
        end do

        ! Loop over n_i/j contributions to entropy
        m = 0
        do i = 1, nSublatticeElements(iSublPhaseIndex,1)
            do j = 1, nSublatticeElements(iSublPhaseIndex,2)
                m = iConstituentSublattice(iSublPhaseIndex,1,i) + &
                ((iConstituentSublattice(iSublPhaseIndex,2,j) - 1) * nSublatticeElements(iSublPhaseIndex,1))
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

        ! Pair indices:
        ii = iPairID(iSublPhaseIndex,k,1)
        jj = iPairID(iSublPhaseIndex,k,2)
        kk = iPairID(iSublPhaseIndex,k,3)
        ll = iPairID(iSublPhaseIndex,k,4)
        ! Anion indices adjusted to start from 1
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
        dSum = (iWeight * (dXij(ii,ka)**dPowXij) * (dXij(ii,la)**dPowXij) &
                        * (dXij(jj,ka)**dPowXij) * (dXij(jj,la)**dPowXij) &
                        / ((dYi(ii)**dPowYi) * (dYi(jj)**dPowYi) &
                        *  (dYi(kk)**dPowYi) * (dYi(ll)**dPowYi)))
        if (dSum == 0) then
            dEntropy = 100D0
        else
            dEntropy = dEntropy + DLOG(dMolFraction(l) / dSum)
        end if

        dRef = dStdGibbsEnergy(l)

        ! Calculate chemical potential of quadruplet
        dChemicalPotential(l) = dRef + dEntropy
    end do

    ! Loop through excess mixing parameters:
    LOOP_Param: do abxy = nParamPhase(iSolnIndex-1) + 1, nParamPhase(iSolnIndex)
        ! AB/XY parametrization
        a = iRegularParam(abxy,2)              ! Index of A
        b = iRegularParam(abxy,3)              ! Index of B
        x = iRegularParam(abxy,4) - nSublatticeElements(iSublPhaseIndex,1)             ! Index of X
        y = iRegularParam(abxy,5) - nSublatticeElements(iSublPhaseIndex,1)             ! Index of Y
        p = iRegularParam(abxy,6)              ! Exponent of A2/X2
        q = iRegularParam(abxy,7)              ! Exponent of B2/X2
        r = iRegularParam(abxy,8)              ! Exponent of A2/Y2
        s = iRegularParam(abxy,9)              ! Exponent of B2/Y2
        if (x == y) then
            iBlock = (x - 1) * (nSublatticeElements(iSublPhaseIndex,1) &
                             * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2)
        else if (x > y) then
            cycle LOOP_Param
        else
            iBlock = (nSublatticeElements(iSublPhaseIndex,2) + (x - 1) + ((y-2)*(y-1)/2)) &
                   * (nSublatticeElements(iSublPhaseIndex,1) * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2)
        end if
        if (a == b) then
            iBlock = iBlock + a
        else if (a > b) then
            cycle LOOP_Param
        else
            iBlock = iBlock + nSublatticeElements(iSublPhaseIndex,1) + a + ((b-2)*(b-1)/2)
        end if
        iBlock = iBlock + iFirst - 1
        iA2X2 = (x - 1) * (nSublatticeElements(iSublPhaseIndex,1) &
                        * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                        + a + iFirst - 1
        iB2X2 = (x - 1) * (nSublatticeElements(iSublPhaseIndex,1) &
                        * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                        + b + iFirst - 1
        iA2Y2 = (y - 1) * (nSublatticeElements(iSublPhaseIndex,1) &
                        * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                        + a + iFirst - 1
        iB2Y2 = (y - 1) * (nSublatticeElements(iSublPhaseIndex,1) &
                        * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                        + b + iFirst - 1
        dXA2X2 = dMolFraction(iA2X2)
        dXB2X2 = dMolFraction(iB2X2)
        dXA2Y2 = dMolFraction(iA2Y2)
        dXB2Y2 = dMolFraction(iB2Y2)

        dGex = 0D0
        dDgexBase = 0D0
        if ((p == 0) .AND. (q == 0) .AND. (r == 0) .AND. (s == 0)) then
            dGex = dExcessGibbsParam(abxy)
        else
            if (p > 0) then
                dGex = dGex + dExcessGibbsParam(abxy) * dXA2X2**p
                dDgexBase = dDgexBase - dExcessGibbsParam(abxy) * p*dXA2X2**p
            end if
            if (q > 0) then
                dGex = dGex + dExcessGibbsParam(abxy) * dXB2X2**q
                dDgexBase = dDgexBase - dExcessGibbsParam(abxy) * q*dXB2X2**q
            end if
            if (r > 0) then
                dGex = dGex + dExcessGibbsParam(abxy) * dXA2Y2**r
                dDgexBase = dDgexBase - dExcessGibbsParam(abxy) * r*dXA2Y2**r
            end if
            if (s > 0) then
                dGex = dGex + dExcessGibbsParam(abxy) * dXB2Y2**s
                dDgexBase = dDgexBase - dExcessGibbsParam(abxy) * s*dXB2Y2**s
            end if
        end if

        ! First add g^ex contribution to quad corresponding to block AB/XY
        if ((a /= b) .AND. (x /= y)) then
            dPartialExcessGibbs(iBlock) = dPartialExcessGibbs(iBlock) + (dGex / 2)

        ! If A = B add g^ex contribution to quads AC/XY
        else if ((a == b) .AND. (x /= y)) then
            dPartialExcessGibbs(iBlock) = dPartialExcessGibbs(iBlock) + (dGex / 2)
            LOOP_AC1: do c = 1, nSublatticeElements(iSublPhaseIndex,1)
                if (c == a) cycle LOOP_AC1
                e = MIN(a,c)
                f = MAX(a,c)
                iQuad = (nSublatticeElements(iSublPhaseIndex,2) + (x - 1) + ((y-2)*(y-1)/2)) &
                      * (nSublatticeElements(iSublPhaseIndex,1) * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                      +  nSublatticeElements(iSublPhaseIndex,1) + e + ((f-2)*(f-1)/2)
                iQuad = iQuad + iFirst - 1
                dPartialExcessGibbs(iQuad) = dPartialExcessGibbs(iQuad) + ((dGex / 4) &
                                           * (dCoordinationNumber(iSublPhaseIndex,iBlock - iFirst + 1,a) &
                                           /  dCoordinationNumber(iSublPhaseIndex,iQuad - iFirst + 1, a)))
            end do LOOP_AC1

        ! If X = Y add g^ex contribution to quads AB/XZ
        else if ((a /= b) .AND. (x == y)) then
            dPartialExcessGibbs(iBlock) = dPartialExcessGibbs(iBlock) + (dGex / 2)
            LOOP_XZ1: do z = 1, nSublatticeElements(iSublPhaseIndex,2)
                if (z == x) cycle LOOP_XZ1
                e = MIN(x,z)
                f = MAX(x,z)
                iQuad = (nSublatticeElements(iSublPhaseIndex,2) + (e - 1) + ((f-2)*(f-1)/2)) &
                      * (nSublatticeElements(iSublPhaseIndex,1) * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                      +  nSublatticeElements(iSublPhaseIndex,1) + a + ((b-2)*(b-1)/2)
                iQuad = iQuad + iFirst - 1
                dPartialExcessGibbs(iQuad) = dPartialExcessGibbs(iQuad) + ((dGex / 4) &
                                           * (dCoordinationNumber(iSublPhaseIndex,iBlock - iFirst + 1,x) &
                                           /  dCoordinationNumber(iSublPhaseIndex,iQuad - iFirst + 1, x)))
            end do LOOP_XZ1
        end if

        ! Now loop over all quads IJ/KL to add dg^ex contributions
        do ijkl = 1, nPairsSRO(iSublPhaseIndex,2)
            i = iPairID(iSublPhaseIndex,ijkl,1)
            j = iPairID(iSublPhaseIndex,ijkl,2)
            k = iPairID(iSublPhaseIndex,ijkl,3) - nSublatticeElements(iSublPhaseIndex,1)
            l = iPairID(iSublPhaseIndex,ijkl,4) - nSublatticeElements(iSublPhaseIndex,1)
            iQuad2 = ijkl + iFirst - 1
            ! Calculate dg^ex/dn_ij/kl
            dDgex = dDgexBase
            if ((i == j) .AND. (k == l)) then
                ! Should probably fix this so only active exponent gets used, but 0s won't hurt for now
                if ((i == a) .AND. (k == x))                dDgex = dDgex + (dExcessGibbsParam(abxy) * p * (dXA2X2**(p-1)))
                if ((i == b) .AND. (k == x) .AND. (a /= b)) dDgex = dDgex + (dExcessGibbsParam(abxy) * q * (dXB2X2**(q-1)))
                if ((i == a) .AND. (k == y) .AND. (x /= y)) dDgex = dDgex + (dExcessGibbsParam(abxy) * r * (dXA2Y2**(r-1)))
                if ((i == b) .AND. (k == y) .AND. (a /= b) .AND. (x /= y)) &
                                                            dDgex = dDgex + (dExcessGibbsParam(abxy) * s * (dXB2Y2**(s-1)))
            end if

            dPartialExcessGibbs(iQuad2) = dPartialExcessGibbs(iQuad2) + (dMolFraction(iBlock) * dDgex / 2)

            if ((a == b) .AND. (x /= y)) then
                LOOP_AC2: do c = 1, nSublatticeElements(iSublPhaseIndex,1)
                    if (c == a) cycle LOOP_AC2
                    e = MIN(a,c)
                    f = MAX(a,c)
                    iQuad = (nSublatticeElements(iSublPhaseIndex,2) + (x - 1) + ((y-2)*(y-1)/2)) &
                          * (nSublatticeElements(iSublPhaseIndex,1) * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                          +  nSublatticeElements(iSublPhaseIndex,1) + e + ((f-2)*(f-1)/2)
                    iQuad = iQuad + iFirst - 1
                    dPartialExcessGibbs(iQuad2) = dPartialExcessGibbs(iQuad2) + ((dMolFraction(iQuad) * dDgex / 4) &
                                              * (dCoordinationNumber(iSublPhaseIndex,iBlock - iFirst + 1,a) &
                                              /  dCoordinationNumber(iSublPhaseIndex,iQuad - iFirst + 1, a)))
                end do LOOP_AC2
            end if
            ! If X = Y add g^ex contribution to quads AB/XZ
            if ((a /= b) .AND. (x == y)) then
                LOOP_XZ2: do z = 1, nSublatticeElements(iSublPhaseIndex,2)
                    if (z == x) cycle LOOP_XZ2
                    e = MIN(x,z)
                    f = MAX(x,z)
                    iQuad = (nSublatticeElements(iSublPhaseIndex,2) + (e - 1) + ((f-2)*(f-1)/2)) &
                          * (nSublatticeElements(iSublPhaseIndex,1) * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                          +  nSublatticeElements(iSublPhaseIndex,1) + a + ((b-2)*(b-1)/2)
                    iQuad = iQuad + iFirst - 1
                    dPartialExcessGibbs(iQuad2) = dPartialExcessGibbs(iQuad2) + ((dMolFraction(iQuad) * dDgex / 4) &
                                              * (dCoordinationNumber(iSublPhaseIndex,iBlock - iFirst + 1,x) &
                                              /  dCoordinationNumber(iSublPhaseIndex,iQuad - iFirst + 1, x)))
                end do LOOP_XZ2
            end if
        end do

    end do LOOP_Param

    ! Deallocate allocatable arrays:
    deallocate(dXi,dYi,dNi,dXij,dNij)

    return

end subroutine CompExcessGibbsEnergySUBG
