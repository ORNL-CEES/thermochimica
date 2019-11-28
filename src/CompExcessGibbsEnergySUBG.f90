
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
    integer :: a, b, c, d, w, x, y, z, e, f, ijkl, abxy, xx, yy
    integer :: iSolnIndex, iSublPhaseIndex, nPhaseElements
    integer :: iFirst, iLast, nA, nX, iWeight, iBlock, iQuad, iQuad2
    integer :: iA2X2, iB2X2, iA2Y2, iADX2, iD2X2, i2
    integer :: iGroupA, iGroupB, iGroupD
    real(8) :: dSum, dEntropy, dRef, dPowXij, dPowYi
    real(8) :: dZa, dZb, dZx, dZy, dGex, dDgex, dDgexBase, dXtot, dYtot
    real(8) :: dXA2X2, dXB2X2, dXA2Y2, dXADX2, dXD2X2, dX2
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
        xx = iRegularParam(abxy,4)             ! Index of X, unadjusted
        yy = iRegularParam(abxy,5)             ! Index of Y, unadjusted
        x = xx - nSublatticeElements(iSublPhaseIndex,1)   ! Index of X
        y = yy - nSublatticeElements(iSublPhaseIndex,1)   ! Index of Y
        p = iRegularParam(abxy,6)              ! Exponent 1
        q = iRegularParam(abxy,7)              ! Exponent 2
        r = iRegularParam(abxy,8)              ! Exponent 3
        s = iRegularParam(abxy,9)              ! Exponent 4
        d = iRegularParam(abxy,10)             ! Index of ternary constituent on 1st sublattice
        w = iRegularParam(abxy,11)             ! Index of ternary constituent on 2nd sublattice

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
        dXA2X2 = dMolFraction(iA2X2)
        dXB2X2 = dMolFraction(iB2X2)
        dXA2Y2 = dMolFraction(iA2Y2)

        ! Calculate energy for this term
        ! G-type binary terms
        if ((cRegularParam(abxy) == 'G') .AND. (d == 0) .AND. (w == 0)) then
            if ((a /= b) .AND. (x == y)) then
                i2 = iB2X2
                dX2 = dXB2X2
            else if ((a == b) .AND. (x /= y)) then
                i2 = iA2Y2
                dX2 = dXA2Y2
            else
                INFOThermo = 42
            end if
            dXtot = dXA2X2 + dX2 + dMolFraction(iBlock)
            dGex = dExcessGibbsParam(abxy) * dXA2X2**p * dX2**q / (dXtot**(p + q))
            dDgexBase = -dGex * (p + q) / dXtot
        ! Q-type binary terms
        else if (cRegularParam(abxy) == 'Q') then
            dYtot = dYi(a) + dYi(b)
            dGex = dExcessGibbsParam(abxy) * dYi(a)**p * dYi(b)**q * dYi(xx)**r * dYi(yy)**s / (dYtot**(p + q + r + s))
            dDgexBase = -dGex * (p + q + r + s) / dYtot
        ! G-type ternary terms
        else if ((cRegularParam(abxy) == 'G') .AND. (d > 0)) then
            iGroupA = iChemicalGroup(iSublPhaseIndex,1,a)
            iGroupB = iChemicalGroup(iSublPhaseIndex,1,b)
            iGroupD = iChemicalGroup(iSublPhaseIndex,1,d)
            ! Symmetric case
            if ((iGroupA == iGroupB) .OR. ((iGroupA /= iGroupB) .AND. (iGroupA /= iGroupD) .AND. (iGroupB /= iGroupD))) then
                ! Assume this is an AB/XX quadruplet
                dGex = dExcessGibbsParam(abxy) * (dXA2X2 / (dXA2X2 + dXB2X2 + dMolFraction(iBlock)))**p &
                                               * (dXB2X2 / (dXA2X2 + dXB2X2 + dMolFraction(iBlock)))**q &
                                               * dYi(d)**r
                dDgexBase = -dGex * r
            ! Asymmetric case
            else
                if (iGroupA == iGroupD) then
                    iD2X2 = (x - 1) * (nSublatticeElements(iSublPhaseIndex,1) &
                                    * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                                    + d + iFirst - 1
                    if (a < d) then
                        iADX2 = (x - 1) * (nSublatticeElements(iSublPhaseIndex,1) &
                                        * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                                          + nSublatticeElements(iSublPhaseIndex,1) + a + ((d-2)*(d-1)/2)
                    else if (a > d) then
                        iADX2 = (x - 1) * (nSublatticeElements(iSublPhaseIndex,1) &
                                        * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                                          + nSublatticeElements(iSublPhaseIndex,1) + d + ((a-2)*(a-1)/2)
                    end if
                    dXD2X2 = dMolFraction(iD2X2)
                    dXADX2 = dMolFraction(iADX2)
                    dGex = dExcessGibbsParam(abxy) * (dXA2X2 + dXD2X2 + dXADX2)**p * dXB2X2**q &
                           * (dYi(d) / (dYi(a) + dYi(d)))**r
                    dDgexBase = -dGex * (p + q)
                else if (iGroupB == iGroupD) then
                    ! Use same variable names but switch A to B in equations
                    iD2X2 = (x - 1) * (nSublatticeElements(iSublPhaseIndex,1) &
                                    * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                                    + d + iFirst - 1
                    if (a < d) then
                        iADX2 = (x - 1) * (nSublatticeElements(iSublPhaseIndex,1) &
                                        * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                                          + nSublatticeElements(iSublPhaseIndex,1) + b + ((d-2)*(d-1)/2)
                    else if (a > d) then
                        iADX2 = (x - 1) * (nSublatticeElements(iSublPhaseIndex,1) &
                                        * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                                          + nSublatticeElements(iSublPhaseIndex,1) + d + ((b-2)*(b-1)/2)
                    end if
                    dXD2X2 = dMolFraction(iD2X2)
                    dXADX2 = dMolFraction(iADX2)
                    dGex = dExcessGibbsParam(abxy) * dXA2X2**p * (dXB2X2 + dXD2X2 + dXADX2)**q &
                           * (dYi(d) / (dYi(b) + dYi(d)))**r
                    dDgexBase = -dGex * (p + q)
                end if
            end if
        ! Reciprocal terms
        else if (cRegularParam(abxy) == 'R') then
            dGex = dExcessGibbsParam(abxy)
            dDgexBase = 0D0
        else
            INFOThermo = 42
        end if


        ! First add g^ex contribution to quad corresponding to block AB/XY
        dPartialExcessGibbs(iBlock) = dPartialExcessGibbs(iBlock) + (dGex / 2)

        ! If A = B add g^ex contribution to quads AC/XY
        if ((a == b) .AND. (x /= y)) then
            LOOP_AC1: do c = 1, nSublatticeElements(iSublPhaseIndex,1)
                if (c == a) cycle LOOP_AC1
                e = MIN(a,c)
                f = MAX(a,c)
                iQuad = (nSublatticeElements(iSublPhaseIndex,2) + (x - 1) + ((y-2)*(y-1)/2)) &
                      * (nSublatticeElements(iSublPhaseIndex,1) * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                      +  nSublatticeElements(iSublPhaseIndex,1) + e + ((f-2)*(f-1)/2)
                iQuad = iQuad + iFirst - 1
                dPartialExcessGibbs(iQuad) = dPartialExcessGibbs(iQuad) + ((dGex / 4) &
                                           * (dCoordinationNumber(iSublPhaseIndex,iBlock - iFirst + 1,1) &
                                           /  dCoordinationNumber(iSublPhaseIndex,iQuad  - iFirst + 1,1)))
            end do LOOP_AC1
        end if

        ! If X = Y add g^ex contribution to quads AB/XZ
        if ((a /= b) .AND. (x == y)) then
            LOOP_XZ1: do z = 1, nSublatticeElements(iSublPhaseIndex,2)
                if (z == x) cycle LOOP_XZ1
                e = MIN(x,z)
                f = MAX(x,z)
                iQuad = (nSublatticeElements(iSublPhaseIndex,2) + (e - 1) + ((f-2)*(f-1)/2)) &
                      * (nSublatticeElements(iSublPhaseIndex,1) * (nSublatticeElements(iSublPhaseIndex,1) + 1) / 2) &
                      +  nSublatticeElements(iSublPhaseIndex,1) + a + ((b-2)*(b-1)/2)
                iQuad = iQuad + iFirst - 1
                dPartialExcessGibbs(iQuad) = dPartialExcessGibbs(iQuad) + ((dGex / 4) &
                                           * (dCoordinationNumber(iSublPhaseIndex,iBlock - iFirst + 1,3) &
                                           /  dCoordinationNumber(iSublPhaseIndex,iQuad  - iFirst + 1,3)))
            end do LOOP_XZ1
        end if

        ! Now loop over all quads IJ/KL to add dg^ex contributions
        LOOP_ijkl: do ijkl = 1, nPairsSRO(iSublPhaseIndex,2)
            iQuad2 = ijkl + iFirst - 1
            i = iPairID(iSublPhaseIndex,ijkl,1)
            j = iPairID(iSublPhaseIndex,ijkl,2)
            k = iPairID(iSublPhaseIndex,ijkl,3) - nSublatticeElements(iSublPhaseIndex,1)
            l = iPairID(iSublPhaseIndex,ijkl,4) - nSublatticeElements(iSublPhaseIndex,1)
            ! Calculate d(g^ex_ab/xy)/d(n_ij/kl)
            ! G-type binary terms
            if ((cRegularParam(abxy) == 'G') .AND. (d == 0) .AND. (w == 0)) then
                if ((iQuad2 /= iA2X2) .AND. (iQuad2 /= i2) .AND. (iQuad2 /= iBlock)) cycle LOOP_ijkl
                dDgex = dDgexBase
                if (iQuad2 == iA2X2) dDgex = dDgex + dGex * p / dXA2X2
                if (iQuad2 == i2)    dDgex = dDgex + dGex * q / dX2
            ! Q-type binary terms
            else if (cRegularParam(abxy) == 'Q') then
                dDgex = 0D0
                if (i == a) dDgex = dDgex + dDgexBase / 2 + dGex * p / (2 * dYi(a))
                if (j == a) dDgex = dDgex + dDgexBase / 2 + dGex * p / (2 * dYi(a))
                if (i == b) dDgex = dDgex + dDgexBase / 2 + dGex * q / (2 * dYi(b))
                if (j == b) dDgex = dDgex + dDgexBase / 2 + dGex * q / (2 * dYi(b))
            ! G-type ternary terms
            else if ((cRegularParam(abxy) == 'G') .AND. (d > 0)) then
                ! Symmetric case
                dDgex = dDgexBase
                if ((iGroupA == iGroupB) .OR. ((iGroupA /= iGroupB) .AND. (iGroupA /= iGroupD) .AND. (iGroupB /= iGroupD))) then
                    ! Assume this is an AB/XX quadruplet
                    if (i == d) dDgex = dDgex + dGex * r / (2 * dYi(d))
                    if (j == d) dDgex = dDgex + dGex * r / (2 * dYi(d))
                    if ((k == x) .AND. (l == x)) then
                        if (((i == a) .AND. ((j == a) .OR. (j == b))) .OR. ((i == b) .AND. (j == b))) &
                            dDgex = dDgex - dGex * (p + q) / (dXA2X2 + dXB2X2 + dMolFraction(iBlock))
                        if ((i == a) .AND. (j == a)) dDgex = dDgex  + dGex * p / dXA2X2
                        if ((i == b) .AND. (j == b)) dDgex = dDgex  + dGex * q / dXB2X2
                    end if
                ! Asymmetric case
                else
                    if (iGroupA == iGroupD) then
                        if (i == a) dDgex = dDgex - dGex * r / (2 * (dYi(a) + dYi(d)))
                        if (j == a) dDgex = dDgex - dGex * r / (2 * (dYi(a) + dYi(d)))
                        if (i == d) dDgex = dDgex + dGex * r * dYi(a) / (2 * dYi(d) * (dYi(a) + dYi(d)))
                        if (j == d) dDgex = dDgex + dGex * r * dYi(a) / (2 * dYi(d) * (dYi(a) + dYi(d)))
                        if (((i == a) .AND. ((j == a) .OR. (j == d))) .OR. ((i == d) .AND. (j == d))) &
                            dDgex = dDgex + dGex * p / (dXA2X2 + dXD2X2 + dXADX2)
                        if ((i == b) .AND. (j == b) .AND. (k == x) .AND. (l == x)) &
                            dDgex = dDgex + dGex * q / dXB2X2
                    else if (iGroupB == iGroupD) then
                        ! Use same variable names but switch A to B in equations
                        if (i == b) dDgex = dDgex - dGex * r / (2 * (dYi(b) + dYi(d)))
                        if (j == b) dDgex = dDgex - dGex * r / (2 * (dYi(b) + dYi(d)))
                        if (i == d) dDgex = dDgex + dGex * r * dYi(b) / (2 * dYi(d) * (dYi(b) + dYi(d)))
                        if (j == d) dDgex = dDgex + dGex * r * dYi(b) / (2 * dYi(d) * (dYi(b) + dYi(d)))
                        if (((i == b) .AND. ((j == b) .OR. (j == d))) .OR. ((i == d) .AND. (j == d))) &
                            dDgex = dDgex + dGex * q / (dXB2X2 + dXD2X2 + dXADX2)
                        if ((i == a) .AND. (j == a) .AND. (k == x) .AND. (l == x)) &
                            dDgex = dDgex + dGex * p / dXA2X2
                    end if
                end if
            end if

            ! print *, cSpeciesName(iQuad2), dPartialExcessGibbs(iQuad2), dDgex
            dPartialExcessGibbs(iQuad2) = dPartialExcessGibbs(iQuad2) + (dMolFraction(iBlock) * dDgex / 2)

            ! If A = B add dg^ex contribution to quads AC/XY
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
                                              * (dCoordinationNumber(iSublPhaseIndex,iBlock - iFirst + 1,1) &
                                              /  dCoordinationNumber(iSublPhaseIndex,iQuad  - iFirst + 1,1)))
                end do LOOP_AC2
            end if
            ! If X = Y add dg^ex contribution to quads AB/XZ
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
                                              * (dCoordinationNumber(iSublPhaseIndex,iBlock - iFirst + 1,3) &
                                              /  dCoordinationNumber(iSublPhaseIndex,iQuad  - iFirst + 1,3)))
                end do LOOP_XZ2
            end if
        end do LOOP_ijkl

    end do LOOP_Param

    ! Deallocate allocatable arrays:
    deallocate(dXi,dYi,dNi,dXij,dNij)

    return

end subroutine CompExcessGibbsEnergySUBG
