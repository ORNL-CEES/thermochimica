
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

    integer :: i, j, k, l, m, ii, jj, kk ,ll, ka, la
    integer :: a, b, c, d, w, x, y, z, e, f, ijkl, abxy, xx, yy
    integer :: iSolnIndex, iSPI, nPhaseElements
    integer :: iFirst, iLast, nA, nX, iWeight, iBlock, iQuad, iQuad2
    integer :: ia, ix
    ! integer :: nAsymmetric1, nAsymmetric2
    logical, allocatable, dimension(:) :: lAsymmetric1, lAsymmetric2
    real(8) :: dSum, dEntropy, dRef, dPowXij, dPowYi, dSumNij, dSumNsij, p, q, r, s
    real(8) :: dZa, dZb, dZx, dZy, dGex, dDgex, dDgexBase, dXtot
    real(8) :: dXi1, dXi2, dChi1, dChi2, dXiDen, dChiDen, dTernaryFactorG, dTernaryFactorDG, dYik, dYjk, dYdk
    real(8) :: dTernarySum1, dTernarySum2
    real(8), allocatable, dimension(:) :: dXi, dYi, dNi, dFi
    real(8), allocatable, dimension(:,:) :: dXij, dNij, dXsij, dNsij
    ! X_ij/kl corresponds to dMolFraction


    ! Only proceed if the correct phase type is selected:
    if (.NOT. (cSolnPhaseType(iSolnIndex) == 'SUBG' .OR. cSolnPhaseType(iSolnIndex) == 'SUBQ')) return

    ! Define temporary variables for sake of convenience:
    iFirst = nSpeciesPhase(iSolnIndex-1) + 1
    iLast  = nSpeciesPhase(iSolnIndex)
    iSPI = iPhaseSublattice(iSolnIndex)

    ! Allocate allocatable arrays:
    if (allocated(dXi)) deallocate(dXi)
    if (allocated(dYi)) deallocate(dYi)
    if (allocated(dFi)) deallocate(dFi)
    if (allocated(dNi)) deallocate(dNi)
    if (allocated(dXij)) deallocate(dXij)
    if (allocated(dNij)) deallocate(dNij)
    if (allocated(dXsij)) deallocate(dXsij)
    if (allocated(dNsij)) deallocate(dNsij)
    if (allocated(lAsymmetric1)) deallocate(lAsymmetric1)
    if (allocated(lAsymmetric2)) deallocate(lAsymmetric2)
    j = iLast - iFirst + 1
    nPhaseElements = nSublatticeElements(iSPI,1) + nSublatticeElements(iSPI,2)
    allocate(dXi(nPhaseElements),dYi(nPhaseElements),dFi(nPhaseElements),dNi(nPhaseElements))
    allocate(dXij(nSublatticeElements(iSPI,1),nSublatticeElements(iSPI,2)))
    allocate(dNij(nSublatticeElements(iSPI,1),nSublatticeElements(iSPI,2)))
    allocate(dXsij(nSublatticeElements(iSPI,1),nSublatticeElements(iSPI,2)))
    allocate(dNsij(nSublatticeElements(iSPI,1),nSublatticeElements(iSPI,2)))
    allocate(lAsymmetric1(MAX(nSublatticeElements(iSPI,1),nSublatticeElements(iSPI,2))))
    allocate(lAsymmetric2(MAX(nSublatticeElements(iSPI,1),nSublatticeElements(iSPI,2))))

    ! Initialize variables:
    dXi                               = 0D0
    dYi                               = 0D0
    dFi                               = 0D0
    dNi                               = 0D0
    dXij                              = 0D0
    dNij                              = 0D0
    dXsij                             = 0D0
    dNsij                             = 0D0
    dChemicalPotential(iFirst:iLast)  = 0D0
    dPartialExcessGibbs(iFirst:iLast) = 0D0

    ! Compute X_i and Y_i
    ! Do cations first:
    dSum = 0D0
    do i = 1, nSublatticeElements(iSPI,1)
        do k = 1, nPairsSRO(iSPI,2)
            l = iFirst + k - 1
            dZa = dCoordinationNumber(iSPI,k,1)
            dZb = dCoordinationNumber(iSPI,k,2)
            if (i == iPairID(iSPI,k,1))  then
                dNi(i) = dNi(i) + (dMolFraction(l) / dZa)
                dYi(i) = dYi(i) + (dMolFraction(l) / 2D0)
            end if
            if (i == iPairID(iSPI,k,2))  then
                dNi(i) = dNi(i) + (dMolFraction(l) / dZb)
                dYi(i) = dYi(i) + (dMolFraction(l) / 2D0)
            end if
        end do
        dSum = dSum + dNi(i)
    end do
    do i = 1, nSublatticeElements(iSPI,1)
        dXi(i) = dNi(i) / dSum
    end do
    ! Do anions now:
    dSum = 0D0
    do i = 1, nSublatticeElements(iSPI,2)
        j = i + nSublatticeElements(iSPI,1)
        do k = 1, nPairsSRO(iSPI,2)
            l = iFirst + k - 1
            dZx = dCoordinationNumber(iSPI,k,3)
            dZy = dCoordinationNumber(iSPI,k,4)
            if (j == iPairID(iSPI,k,3))  then
                dNi(j) = dNi(j) + (dMolFraction(l) / dZx)
                dYi(j) = dYi(j) + (dMolFraction(l) / 2)
            end if
            if (j == iPairID(iSPI,k,4))  then
                dNi(j) = dNi(j) + (dMolFraction(l) / dZy)
                dYi(j) = dYi(j) + (dMolFraction(l) / 2)
            end if
        end do
        dSum = dSum + dNi(j)
    end do
    do i = 1, nSublatticeElements(iSPI,2)
        j = i + nSublatticeElements(iSPI,1)
        dXi(j) = dNi(j) / dSum
    end do

    ! Compute X_i/j
    dSumNij = 0D0
    dSumNsij = 0D0
    do i = 1, nSublatticeElements(iSPI,1)
        do j = 1, nSublatticeElements(iSPI,2)
            m = iConstituentSublattice(iSPI,1,i) + &
            ((iConstituentSublattice(iSPI,2,j) - 1) * nSublatticeElements(iSPI,1))
            do k = 1, nPairsSRO(iSPI,2)
                l = iFirst + k - 1
                nA = 0
                if (i == iPairID(iSPI,k,1))  then
                    nA = nA + 1
                end if
                if (i == iPairID(iSPI,k,2))  then
                    nA = nA + 1
                end if
                nX = 0
                if ((j + nSublatticeElements(iSPI,1)) == iPairID(iSPI,k,3))  then
                    nX = nX + 1
                end if
                if ((j + nSublatticeElements(iSPI,1)) == iPairID(iSPI,k,4))  then
                    nX = nX + 1
                end if
                dNij(i,j)  = dNij(i,j)  + (dMolFraction(l) * nA * nX)
                dNsij(i,j) = dNsij(i,j) + (dMolFraction(l) * nA * nX / dZetaSpecies(iSPI,m))
            end do
            dSumNij  = dSumNij  + dNij(i,j)
            dSumNsij = dSumNsij + dNsij(i,j)
        end do
    end do

    do i = 1, nSublatticeElements(iSPI,1)
        do j = 1, nSublatticeElements(iSPI,2)
            dXij(i,j)  = dNij(i,j)  / dSumNij
            dXsij(i,j) = dNsij(i,j) / dSumNsij
        end do
    end do

    ! For updated implementation, calculate F_i
    do i = 1, nSublatticeElements(iSPI,1)
        do j = 1, nSublatticeElements(iSPI,2)
            dFi(i) = dFi(i) + dXsij(i,j)
            k = j + nSublatticeElements(iSPI,1)
            dFi(k) = dFi(k) + dXsij(i,j)
        end do
    end do

    ! ---------------------------------------------------------------
    ! COMPUTE REFERENCE GIBBS ENERGY AND IDEAL MIXING TERMS
    ! ---------------------------------------------------------------

    do k = 1, nPairsSRO(iSPI,2)
        ! Calculate entropic contributions to chemical potentials
        dEntropy = 0D0
        dRef = 0D0
        l = iFirst + k - 1

        ! Coordination numbers for this quadruplet:
        dZa = dCoordinationNumber(iSPI,k,1)
        dZb = dCoordinationNumber(iSPI,k,2)
        dZx = dCoordinationNumber(iSPI,k,3)
        dZy = dCoordinationNumber(iSPI,k,4)

        ! Loop over n_i contributions to entropy
        ! Cations first
        do i = 1, nSublatticeElements(iSPI,1)
            if (i == iPairID(iSPI,k,1))  then
                dEntropy = dEntropy + (DLOG(dXi(i)) / dZa)
            end if
            if (i == iPairID(iSPI,k,2))  then
                dEntropy = dEntropy + (DLOG(dXi(i)) / dZb)
            end if
        end do
        ! Now anions
        do i = 1, nSublatticeElements(iSPI,2)
            j = i + nSublatticeElements(iSPI,1)
            if (j == iPairID(iSPI,k,3))  then
                dEntropy = dEntropy + (DLOG(dXi(j)) / dZx)
            end if
            if (j == iPairID(iSPI,k,4))  then
                dEntropy = dEntropy + (DLOG(dXi(j)) / dZy)
            end if
        end do

        ! Loop over n_i/j contributions to entropy
        m = 0
        do i = 1, nSublatticeElements(iSPI,1)
            do j = 1, nSublatticeElements(iSPI,2)
                m = iConstituentSublattice(iSPI,1,i) + &
                ((iConstituentSublattice(iSPI,2,j) - 1) * nSublatticeElements(iSPI,1))
                nA = 0
                if (i == iPairID(iSPI,k,1))  then
                    nA = nA + 1
                end if
                if (i == iPairID(iSPI,k,2))  then
                    nA = nA + 1
                end if
                nX = 0
                if ((j + nSublatticeElements(iSPI,1)) == iPairID(iSPI,k,3))  then
                    nX = nX + 1
                end if
                if ((j + nSublatticeElements(iSPI,1)) == iPairID(iSPI,k,4))  then
                    nX = nX + 1
                end if
                dEntropy = dEntropy + (DLOG(dXsij(i,j) / (dFi(i) * dFi(j + nSublatticeElements(iSPI,1)))) &
                                    * (nA * nX / dZetaSpecies(iSPI,m)))
            end do
        end do

        ! Pair indices:
        ii = iPairID(iSPI,k,1)
        jj = iPairID(iSPI,k,2)
        kk = iPairID(iSPI,k,3)
        ll = iPairID(iSPI,k,4)
        ! Anion indices adjusted to start from 1
        ka = kk - nSublatticeElements(iSPI,1)
        la = ll - nSublatticeElements(iSPI,1)

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

        if (dExcessGibbsParam(abxy) == 0D0) cycle LOOP_Param

        ! AB/XY parametrization
        a = iRegularParam(abxy,2)              ! Index of A
        b = iRegularParam(abxy,3)              ! Index of B
        xx = iRegularParam(abxy,4)             ! Index of X, unadjusted
        yy = iRegularParam(abxy,5)             ! Index of Y, unadjusted
        x = xx - nSublatticeElements(iSPI,1)   ! Index of X
        y = yy - nSublatticeElements(iSPI,1)   ! Index of Y
        p = iRegularParam(abxy,6)              ! Exponent 1
        q = iRegularParam(abxy,7)              ! Exponent 2
        r = iRegularParam(abxy,8)              ! Exponent 3
        s = iRegularParam(abxy,9)              ! Exponent 4
        d = iRegularParam(abxy,10)             ! Index of ternary constituent on 1st sublattice
        w = iRegularParam(abxy,11)             ! Index of ternary constituent on 2nd sublattice

        if (x == y) then
            iBlock = (x - 1) * (nSublatticeElements(iSPI,1) &
                             * (nSublatticeElements(iSPI,1) + 1) / 2)
        else if (x > y) then
            cycle LOOP_Param
        else
            iBlock = (nSublatticeElements(iSPI,2) + (x - 1) + ((y-2)*(y-1)/2)) &
                   * (nSublatticeElements(iSPI,1) * (nSublatticeElements(iSPI,1) + 1) / 2)
        end if
        if (a == b) then
            iBlock = iBlock + a
        else if (a > b) then
            cycle LOOP_Param
        else
            iBlock = iBlock + nSublatticeElements(iSPI,1) + a + ((b-2)*(b-1)/2)
        end if
        iBlock = iBlock + iFirst - 1

        dXi1 = 0D0
        dXi2 = 0D0
        dXiDen = 0D0
        dChi1 = 0D0
        dChi2 = 0D0
        dChiDen = 0D0
        lAsymmetric1 = .FALSE.
        lAsymmetric2 = .FALSE.
        if (x == y) then
            lAsymmetric1(a) = .TRUE.
            lAsymmetric2(b) = .TRUE.
            ! First make a list of which constituents make asymmetric ternaries
            if (iChemicalGroup(iSPI,1,a) /= iChemicalGroup(iSPI,1,b)) then
                do i = 1, nSublatticeElements(iSPI,1)
                    if (iChemicalGroup(iSPI,1,i) == iChemicalGroup(iSPI,1,a)) then
                        lAsymmetric1(i) = .TRUE.
                    else if (iChemicalGroup(iSPI,1,i) == iChemicalGroup(iSPI,1,b)) then
                        lAsymmetric2(i) = .TRUE.
                    end if
                end do
            end if
            ! Now use lists to generate xi and chi
            do i = 1, nSublatticeElements(iSPI,1)
                do j = i, nSublatticeElements(iSPI,1)
                    k = (x - 1) * (nSublatticeElements(iSPI,1) &
                                    * (nSublatticeElements(iSPI,1) + 1) / 2)
                    if (i == j) then
                        k = k + i
                    else
                        k = k + nSublatticeElements(iSPI,1) + i + ((j-2)*(j-1)/2)
                    end if
                    iQuad = k + iFirst - 1
                    if (lAsymmetric1(i) .AND. lAsymmetric1(j)) then
                        dChi1 = dChi1 + dMolFraction(iQuad)
                    end if
                    if (lAsymmetric2(i) .AND. lAsymmetric2(j)) then
                        dChi2 = dChi2 + dMolFraction(iQuad)
                    end if
                    if ((lAsymmetric1(i) .OR. lAsymmetric2(i)) .AND. (lAsymmetric1(j) .OR. lAsymmetric2(j))) then
                        dChiDen = dChiDen + dMolFraction(iQuad)
                    end if
                end do
                ! Below is xi with counting of x /= y quads
                do k = 1, nPairsSRO(iSPI,2)
                    l = k + iFirst - 1
                    if (lAsymmetric1(i)) then
                        if (i == iPairID(iSPI,k,1) .AND. xx == iPairID(iSPI,k,3)) dXi1 = dXi1 + (dMolFraction(l) / 4)
                        if (i == iPairID(iSPI,k,1) .AND. xx == iPairID(iSPI,k,4)) dXi1 = dXi1 + (dMolFraction(l) / 4)
                        if (i == iPairID(iSPI,k,2) .AND. xx == iPairID(iSPI,k,3)) dXi1 = dXi1 + (dMolFraction(l) / 4)
                        if (i == iPairID(iSPI,k,2) .AND. xx == iPairID(iSPI,k,4)) dXi1 = dXi1 + (dMolFraction(l) / 4)
                    end if
                    if (lAsymmetric2(i)) then
                        if (i == iPairID(iSPI,k,1) .AND. xx == iPairID(iSPI,k,3)) dXi2 = dXi2 + (dMolFraction(l) / 4)
                        if (i == iPairID(iSPI,k,1) .AND. xx == iPairID(iSPI,k,4)) dXi2 = dXi2 + (dMolFraction(l) / 4)
                        if (i == iPairID(iSPI,k,2) .AND. xx == iPairID(iSPI,k,3)) dXi2 = dXi2 + (dMolFraction(l) / 4)
                        if (i == iPairID(iSPI,k,2) .AND. xx == iPairID(iSPI,k,4)) dXi2 = dXi2 + (dMolFraction(l) / 4)
                    end if
                end do
            end do
        else if (a == b) then
            lAsymmetric1(x) = .TRUE.
            lAsymmetric2(y) = .TRUE.
            ! First make a list of which constituents make asymmetric ternaries
            if (iChemicalGroup(iSPI,2,x) /= iChemicalGroup(iSPI,2,y)) then
                do i = 1, nSublatticeElements(iSPI,2)
                    if (iChemicalGroup(iSPI,2,i) == iChemicalGroup(iSPI,2,x)) then
                        lAsymmetric1(i) = .TRUE.
                    else if (iChemicalGroup(iSPI,2,i) == iChemicalGroup(iSPI,2,y)) then
                        lAsymmetric2(i) = .TRUE.
                    end if
                end do
            end if
            ! Now use lists to generate xi and chi
            do i = 1, nSublatticeElements(iSPI,2)
                do j = i, nSublatticeElements(iSPI,2)
                    k = a
                    if (i == j) then
                        k = k + (i - 1) * (nSublatticeElements(iSPI,1) &
                                        * (nSublatticeElements(iSPI,1) + 1) / 2)
                    else
                        k = k + (nSublatticeElements(iSPI,2) + (x - 1) + ((y-2)*(y-1)/2)) &
                              * (nSublatticeElements(iSPI,1) * (nSublatticeElements(iSPI,1) + 1) / 2)
                    end if
                    iQuad = k + iFirst - 1
                    if (lAsymmetric1(i) .AND. lAsymmetric1(j)) then
                        dChi1 = dChi1 + dMolFraction(iQuad)
                    end if
                    if (lAsymmetric2(i) .AND. lAsymmetric2(j)) then
                        dChi2 = dChi2 + dMolFraction(iQuad)
                    end if
                    if ((lAsymmetric1(i) .OR. lAsymmetric2(i)) .AND. (lAsymmetric1(j) .OR. lAsymmetric2(j))) then
                        dChiDen = dChiDen + dMolFraction(iQuad)
                    end if
                end do
                ! Below is xi with counting of x /= y quads
                ii = i + nSublatticeElements(iSPI,1)
                do k = 1, nPairsSRO(iSPI,2)
                    l = k + iFirst - 1
                    if (lAsymmetric1(i)) then
                        if (a == iPairID(iSPI,k,1) .AND. ii == iPairID(iSPI,k,3)) dXi1 = dXi1 + (dMolFraction(l) / 4)
                        if (a == iPairID(iSPI,k,1) .AND. ii == iPairID(iSPI,k,4)) dXi1 = dXi1 + (dMolFraction(l) / 4)
                        if (a == iPairID(iSPI,k,2) .AND. ii == iPairID(iSPI,k,3)) dXi1 = dXi1 + (dMolFraction(l) / 4)
                        if (a == iPairID(iSPI,k,2) .AND. ii == iPairID(iSPI,k,4)) dXi1 = dXi1 + (dMolFraction(l) / 4)
                    end if
                    if (lAsymmetric2(i)) then
                        if (a == iPairID(iSPI,k,1) .AND. ii == iPairID(iSPI,k,3)) dXi2 = dXi2 + (dMolFraction(l) / 4)
                        if (a == iPairID(iSPI,k,1) .AND. ii == iPairID(iSPI,k,4)) dXi2 = dXi2 + (dMolFraction(l) / 4)
                        if (a == iPairID(iSPI,k,2) .AND. ii == iPairID(iSPI,k,3)) dXi2 = dXi2 + (dMolFraction(l) / 4)
                        if (a == iPairID(iSPI,k,2) .AND. ii == iPairID(iSPI,k,4)) dXi2 = dXi2 + (dMolFraction(l) / 4)
                    end if
                end do
            end do
        end if

        dXiDen = dXi1 + dXi2
        dChi1 = dChi1 / dChiDen
        dChi2 = dChi2 / dChiDen

        dTernaryFactorG = 1D0
        if (d > 0) then
            dYdk = 0D0
            do k = 1, nPairsSRO(iSPI,2)
                l = k + iFirst - 1
                if (d == iPairID(iSPI,k,1) .AND. xx == iPairID(iSPI,k,3)) dYdk = dYdk + (dMolFraction(l) / 4)
                if (d == iPairID(iSPI,k,1) .AND. xx == iPairID(iSPI,k,4)) dYdk = dYdk + (dMolFraction(l) / 4)
                if (d == iPairID(iSPI,k,2) .AND. xx == iPairID(iSPI,k,3)) dYdk = dYdk + (dMolFraction(l) / 4)
                if (d == iPairID(iSPI,k,2) .AND. xx == iPairID(iSPI,k,4)) dYdk = dYdk + (dMolFraction(l) / 4)
            end do
            if (lAsymmetric2(d)) then
                dYjk = 0D0
                do k = 1, nPairsSRO(iSPI,2)
                    l = k + iFirst - 1
                    if (b == iPairID(iSPI,k,1) .AND. xx == iPairID(iSPI,k,3)) dYjk = dYjk + (dMolFraction(l) / 4)
                    if (b == iPairID(iSPI,k,1) .AND. xx == iPairID(iSPI,k,4)) dYjk = dYjk + (dMolFraction(l) / 4)
                    if (b == iPairID(iSPI,k,2) .AND. xx == iPairID(iSPI,k,3)) dYjk = dYjk + (dMolFraction(l) / 4)
                    if (b == iPairID(iSPI,k,2) .AND. xx == iPairID(iSPI,k,4)) dYjk = dYjk + (dMolFraction(l) / 4)
                end do
                dTernaryFactorG = (dYdk / dXi2) * (1 - (dYjk / dXi2))**(r-1)
            else if (lAsymmetric1(d)) then
                dYik = 0D0
                do k = 1, nPairsSRO(iSPI,2)
                    l = k + iFirst - 1
                    if (a == iPairID(iSPI,k,1) .AND. xx == iPairID(iSPI,k,3)) dYik = dYik + (dMolFraction(l) / 4)
                    if (a == iPairID(iSPI,k,1) .AND. xx == iPairID(iSPI,k,4)) dYik = dYik + (dMolFraction(l) / 4)
                    if (a == iPairID(iSPI,k,2) .AND. xx == iPairID(iSPI,k,3)) dYik = dYik + (dMolFraction(l) / 4)
                    if (a == iPairID(iSPI,k,2) .AND. xx == iPairID(iSPI,k,4)) dYik = dYik + (dMolFraction(l) / 4)
                end do
                dTernaryFactorG = (dYdk / dXi1) * (1 - (dYik / dXi1))**(r-1)
            else
                dTernaryFactorG = dYdk * (1D0 - dXi1 - dXi2)**(r-1D0)
            end if
        end if

        ! Calculate energy for this term
        ! G-type terms
        if (cRegularParam(abxy) == 'G') then
            ! G-type binary terms
            dGex = (dExcessGibbsParam(abxy) * dChi1**p * dChi2**q) * dTernaryFactorG
            dDgexBase = -dGex * (p + q) / dChiDen
        ! Q-type binary terms
        else if (cRegularParam(abxy) == 'Q') then
            dGex = (dExcessGibbsParam(abxy) * dXi1**p * dXi2**q / (dXiDen**(p + q))) * dTernaryFactorG
            dDgexBase = -dGex * (p + q) / dXiDen
        ! B-type binary terms
        else if ((cRegularParam(abxy) == 'B') .OR. (cRegularParam(abxy) == 'H')) then
            ! Start by calculating mixing energy term
            dXtot = dXsij(a,x) + dXsij(b,y)
            dGex = dExcessGibbsParam(abxy) * dXsij(a,x)**(1D0+p) * dXsij(b,y)**(1D0+q) / dXtot**(1D0+p+q)
            dDgexBase = - dGex / dSumNsij
            LOOP_Bder: do i = 1, nSublatticeElements(iSPI,1)
                do j = 1, nSublatticeElements(iSPI,2)
                    dDgex = dDgexBase
                    if ((i == a) .AND. (j == x)) then
                        ! the if below is just to prove that there are numerical issues with this mixing scheme
                        ! if (iterGlobal > 10) then
                        dDgex = dDgex + dGex*(dNsij(b,y) + dNsij(b,y)*p - dNsij(a,x)*q) / (dNsij(a,x) * (dNsij(b,y) + dNsij(a,x)))
                        ! end if
                    else if ((i == b) .AND. (j == y)) then
                        dDgex = dDgex + dGex*(dNsij(a,x) - dNsij(b,y)*p + dNsij(a,x)*q) / (dNsij(b,y) * (dNsij(b,y) + dNsij(a,x)))
                    end if
                    m = iConstituentSublattice(iSPI,1,i) + ((iConstituentSublattice(iSPI,2,j) - 1) * nSublatticeElements(iSPI,1))
                    do k = 1, nPairsSRO(iSPI,2)
                        l = iFirst + k - 1
                        nA = 0
                        if (i == iPairID(iSPI,k,1))  then
                            nA = nA + 1
                        end if
                        if (i == iPairID(iSPI,k,2))  then
                            nA = nA + 1
                        end if
                        nX = 0
                        if ((j + nSublatticeElements(iSPI,1)) == iPairID(iSPI,k,3))  then
                            nX = nX + 1
                        end if
                        if ((j + nSublatticeElements(iSPI,1)) == iPairID(iSPI,k,4))  then
                            nX = nX + 1
                        end if
                        ! Add derivative contribution
                        dPartialExcessGibbs(l) = dPartialExcessGibbs(l) + (dDgex * nA * nX / dZetaSpecies(iSPI,m))
                        ! Add dGex to every quadruplet chemical potential
                        dPartialExcessGibbs(l) = dPartialExcessGibbs(l) + (dGex  * nA * nX / 4D0)
                    end do
                end do
            end do LOOP_Bder
            cycle LOOP_Param
        ! Reciprocal terms
        else if (cRegularParam(abxy) == 'R') then
            dGex = dExcessGibbsParam(abxy)
            dDgexBase = 0D0
        else
            INFOThermo = 42
        end if

        dPartialExcessGibbs(iBlock) = dPartialExcessGibbs(iBlock) + (dGex / 2)

        ! If A = B add g^ex contribution to quads AC/XY
        if ((a == b) .AND. (x /= y)) then
            LOOP_AC1: do c = 1, nSublatticeElements(iSPI,1)
                if (c == a) cycle LOOP_AC1
                e = MIN(a,c)
                f = MAX(a,c)
                ia = MINLOC((/a,c/),1)
                iQuad = (nSublatticeElements(iSPI,2) + (x - 1) + ((y-2)*(y-1)/2)) &
                      * (nSublatticeElements(iSPI,1) * (nSublatticeElements(iSPI,1) + 1) / 2) &
                      +  nSublatticeElements(iSPI,1) + e + ((f-2)*(f-1)/2)
                iQuad = iQuad + iFirst - 1
                ! if (cRegularParam(abxy) /= 'B') then
                    dPartialExcessGibbs(iQuad) = dPartialExcessGibbs(iQuad) + ((dGex / 4) &
                                               * (dCoordinationNumber(iSPI,iBlock - iFirst + 1,1) &
                                               /  dCoordinationNumber(iSPI,iQuad  - iFirst + 1,ia)))
                ! end if
            end do LOOP_AC1
        end if

        ! If X = Y add g^ex contribution to quads AB/XZ
        if ((a /= b) .AND. (x == y)) then
            LOOP_XZ1: do z = 1, nSublatticeElements(iSPI,2)
                if (z == x) cycle LOOP_XZ1
                e = MIN(x,z)
                f = MAX(x,z)
                ix = MINLOC((/x,z/),1) + 2
                iQuad = (nSublatticeElements(iSPI,2) + (e - 1) + ((f-2)*(f-1)/2)) &
                      * (nSublatticeElements(iSPI,1) * (nSublatticeElements(iSPI,1) + 1) / 2) &
                      +  nSublatticeElements(iSPI,1) + a + ((b-2)*(b-1)/2)
                iQuad = iQuad + iFirst - 1
                dPartialExcessGibbs(iQuad) = dPartialExcessGibbs(iQuad) + ((dGex / 4) &
                                           * (dCoordinationNumber(iSPI,iBlock - iFirst + 1,3) &
                                           /  dCoordinationNumber(iSPI,iQuad  - iFirst + 1,ix)))
            end do LOOP_XZ1
        end if

        ! Now loop over all quads IJ/KL to add dg^ex contributions
        LOOP_ijkl: do ijkl = 1, nPairsSRO(iSPI,2)
            iQuad2 = ijkl + iFirst - 1
            i = iPairID(iSPI,ijkl,1)
            j = iPairID(iSPI,ijkl,2)
            k = iPairID(iSPI,ijkl,3) - nSublatticeElements(iSPI,1)
            l = iPairID(iSPI,ijkl,4) - nSublatticeElements(iSPI,1)

            dTernaryFactorDG = 0D0
            if (d > 0) then
                if (lAsymmetric2(d)) then
                    nA = 0
                    if (i == d) nA = nA + 1
                    if (j == d) nA = nA + 1
                    nX = 0
                    if (k == x) nX = nX + 1
                    if (l == x) nX = nX + 1
                    dTernaryFactorDG = dTernaryFactorDG + (nA * nX) / (4D0 * dYdk)

                    dTernarySum2 = 0D0
                    do e = 1, nSublatticeElements(iSPI,1)
                        nA = 0
                        if (i == e) nA = nA + 1
                        if (j == e) nA = nA + 1
                        if (lAsymmetric2(e)) then
                            dTernarySum2 = dTernarySum2 + (nA * nX) / (4D0 * dXi2)
                        end if
                    end do

                    dTernaryFactorDG = dTernaryFactorDG - dTernarySum2

                    nA = 0
                    if (i == b) nA = nA + 1
                    if (j == b) nA = nA + 1

                    dTernaryFactorDG = dTernaryFactorDG - (r - 1D0) * ((nA * nX / 4D0) - dYjk * dTernarySum2) &
                                        / (dXi2 * (1D0 - dYjk / dXi2))
                else if (lAsymmetric1(d)) then
                    nA = 0
                    if (i == d) nA = nA + 1
                    if (j == d) nA = nA + 1
                    nX = 0
                    if (k == x) nX = nX + 1
                    if (l == x) nX = nX + 1
                    dTernaryFactorDG = dTernaryFactorDG + (nA * nX) / (4D0 * dYdk)

                    dTernarySum1 = 0D0
                    do e = 1, nSublatticeElements(iSPI,1)
                        nA = 0
                        if (i == e) nA = nA + 1
                        if (j == e) nA = nA + 1
                        if (lAsymmetric1(e)) then
                            dTernarySum1 = dTernarySum1 + (nA * nX) / (4D0 * dXi1)
                        end if
                    end do

                    dTernaryFactorDG = dTernaryFactorDG - dTernarySum1

                    nA = 0
                    if (i == a) nA = nA + 1
                    if (j == a) nA = nA + 1

                    dTernaryFactorDG = dTernaryFactorDG - (r - 1D0) * ((nA * nX / 4D0) - dYik * dTernarySum1) &
                                        / (dXi1 * (1D0 - dYik / dXi1))
                else
                    dTernaryFactorDG = -r
                    nA = 0
                    if (i == d) nA = nA + 1
                    if (j == d) nA = nA + 1
                    nX = 0
                    if (k == x) nX = nX + 1
                    if (l == x) nX = nX + 1
                    dTernaryFactorDG = dTernaryFactorDG + (nA * nX) / (4D0 * dYdk)

                    dTernarySum1 = 0D0
                    dTernarySum2 = 0D0
                    do e = 1, nSublatticeElements(iSPI,1)
                        nA = 0
                        if (i == e) nA = nA + 1
                        if (j == e) nA = nA + 1
                        if (lAsymmetric1(e)) then
                            dTernarySum1 = dTernarySum1 + (nA * nX) / 4D0
                        end if
                        if (lAsymmetric2(e)) then
                            dTernarySum2 = dTernarySum2 + (nA * nX) / 4D0
                        end if
                    end do
                    dTernaryFactorDG = dTernaryFactorDG + (r - 1D0) * (1D0 - dTernarySum1 - dTernarySum2) / (1D0 - dXi1 - dXi2)
                end if
            end if

            ! Calculate d(g^ex_ab/xy)/d(n_ij/kl)
            ! G-type terms
            if (cRegularParam(abxy) == 'G') then
                dDgex = 0D0
                if (.NOT.((x == y) .AND. ((k /= x) .OR. (l /= x)))) then
                    if (lAsymmetric1(i) .AND. lAsymmetric1(j)) then
                        dDgex = dDgex + dGex * p / dChi1 / dChiDen
                    end if
                    if (lAsymmetric2(i) .AND. lAsymmetric2(j)) then
                        dDgex = dDgex + dGex * q / dChi2 / dChiDen
                    end if
                    if ((lAsymmetric1(i) .OR. lAsymmetric2(i)) .AND. (lAsymmetric1(j) .OR. lAsymmetric2(j))) then
                        dDgex = dDgex + dDgexBase
                    end if
                end if

                dDgex = dDgex + dGex * dTernaryFactorDG
            ! Q-type terms
            else if (cRegularParam(abxy) == 'Q') then
                dDgex = 0D0
                do ii = 1, nSublatticeElements(iSPI,1)
                    ! Below is xi with counting of x /= y quads
                    if (lAsymmetric1(ii)) then
                        if (ii == i .AND. x == k) dDgex = dDgex + dDgexBase / 4 + dGex * p / (4 * dXi1)
                        if (ii == i .AND. x == l) dDgex = dDgex + dDgexBase / 4 + dGex * p / (4 * dXi1)
                        if (ii == j .AND. x == k) dDgex = dDgex + dDgexBase / 4 + dGex * p / (4 * dXi1)
                        if (ii == j .AND. x == l) dDgex = dDgex + dDgexBase / 4 + dGex * p / (4 * dXi1)
                    end if
                    if (lAsymmetric2(ii)) then
                        if (ii == i .AND. x == k) dDgex = dDgex + dDgexBase / 4 + dGex * q / (4 * dXi2)
                        if (ii == i .AND. x == l) dDgex = dDgex + dDgexBase / 4 + dGex * q / (4 * dXi2)
                        if (ii == j .AND. x == k) dDgex = dDgex + dDgexBase / 4 + dGex * q / (4 * dXi2)
                        if (ii == j .AND. x == l) dDgex = dDgex + dDgexBase / 4 + dGex * q / (4 * dXi2)
                    end if
                end do

                dDgex = dDgex + dGex * dTernaryFactorDG
            end if

            dPartialExcessGibbs(iQuad2) = dPartialExcessGibbs(iQuad2) + (dMolFraction(iBlock) * dDgex / 2)

            ! If A = B add dg^ex contribution from quads AC/XY to IJ/KL
            if ((a == b) .AND. (x /= y)) then
                LOOP_AC2: do c = 1, nSublatticeElements(iSPI,1)
                    if (c == a) cycle LOOP_AC2
                    e = MIN(a,c)
                    f = MAX(a,c)
                    ia = MINLOC((/a,c/),1)
                    iQuad = (nSublatticeElements(iSPI,2) + (x - 1) + ((y-2)*(y-1)/2)) &
                          * (nSublatticeElements(iSPI,1) * (nSublatticeElements(iSPI,1) + 1) / 2) &
                          +  nSublatticeElements(iSPI,1) + e + ((f-2)*(f-1)/2)
                    iQuad = iQuad + iFirst - 1
                    dPartialExcessGibbs(iQuad2) = dPartialExcessGibbs(iQuad2) + ((dMolFraction(iQuad) * dDgex / 4) &
                                              * (dCoordinationNumber(iSPI,iBlock - iFirst + 1,1) &
                                              /  dCoordinationNumber(iSPI,iQuad  - iFirst + 1,ia)))
                end do LOOP_AC2
            end if
            ! If X = Y add dg^ex contribution from quads AB/XZ to IJ/KL
            if ((a /= b) .AND. (x == y)) then
                ! if (.NOT.((l == x).AND.(k == x))) cycle LOOP_ijkl
                LOOP_XZ2: do z = 1, nSublatticeElements(iSPI,2)
                    if (z == x) cycle LOOP_XZ2
                    e = MIN(x,z)
                    f = MAX(x,z)
                    ix = MINLOC((/x,z/),1) + 2
                    iQuad = (nSublatticeElements(iSPI,2) + (e - 1) + ((f-2)*(f-1)/2)) &
                          * (nSublatticeElements(iSPI,1) * (nSublatticeElements(iSPI,1) + 1) / 2) &
                          +  nSublatticeElements(iSPI,1) + a + ((b-2)*(b-1)/2)
                    iQuad = iQuad + iFirst - 1
                    dPartialExcessGibbs(iQuad2) = dPartialExcessGibbs(iQuad2) + ((dMolFraction(iQuad) * dDgex / 4) &
                                              * (dCoordinationNumber(iSPI,iBlock - iFirst + 1,3) &
                                              /  dCoordinationNumber(iSPI,iQuad  - iFirst + 1,ix)))
                end do LOOP_XZ2
            end if
        end do LOOP_ijkl

    end do LOOP_Param

    ! Deallocate allocatable arrays:
    deallocate(dXi,dYi,dFi,dNi,dXij,dNij,dXsij,dNsij,lAsymmetric1,lAsymmetric2)

    return

end subroutine CompExcessGibbsEnergySUBG
