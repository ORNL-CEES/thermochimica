
!-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompExcessGibbsEnergySUBM.f90
    !> \brief   Compute the partial molar excess Gibbs energy of mixing of 2-sublattice solution phase.
    !> \author  M. Poschmann
    !> \date    October 25, 2022
    !> \sa      CompExcessGibbsEnergy.f90
    !> \sa      CompExcessGibbsEnergyRKMP.f90
    !> \sa      CompExcessGibbsEnergyQKTO.f90
    !> \sa      CompExcessGibbsEnergySUBL.f90
    !> \sa      CompExcessGibbsEnergySUBG.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   10/25/2022      M. Poschmann        Original code.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the chemical potentials of species in the
    !! 2-sublattice solution phase model. This model is equivalent to IS2L (SUBI), except has no neutral
    !! or vacancy species.
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
    ! dX(:)                     A temporary double real vector used to represent the mole fractions of the
    !                            species.
    ! dY(:)                     A temporary double vector used to represent the coordinate equivalent
    !                            fractions of the species.
    !
!-------------------------------------------------------------------------------------------------------------


subroutine CompExcessGibbsEnergySUBM(iSolnIndex)

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, l, n
    integer :: a, b, c, x, xx, order, iGroup1, iGroup2, iGroupTemp
    integer :: iSolnIndex, iSPI, nPhaseElements, nSub1, nSub2, nA2X2
    integer :: iFirst, iLast, iStartCon, iEndCon, iSub, iOffset
    logical :: lIsException
    logical, allocatable, dimension(:) :: lAsymmetric1, lAsymmetric2
    real(8) :: dSum1, dSum2, q, p, ea, eb, ec, ex, yc, dYcfac
    real(8) :: gref, gideal, gex, natom, dMol, lc1, lc2, dMolAtoms, gexTemp, dgexTemp
    real(8) :: dXi1, dXi2, dXiDen
    real(8), allocatable, dimension(:) :: dXi, dYi, dNi, dgdc, dMolDerivatives, dSumY

    ! Only proceed if the correct phase type is selected:
    if (.NOT. (cSolnPhaseType(iSolnIndex) == 'SUBM')) return

    ! Define temporary variables for sake of convenience:
    iFirst = nSpeciesPhase(iSolnIndex-1) + 1
    iLast  = nSpeciesPhase(iSolnIndex)
    iSPI = iPhaseSublattice(iSolnIndex)
    nSub1 = nConstituentSublattice(iSPI,1)
    nSub2 = nConstituentSublattice(iSPI,2)
    nA2X2 = nSub1 * nSub2

    ! Allocate allocatable arrays:
    if (allocated(dXi)) deallocate(dXi)
    if (allocated(dYi)) deallocate(dYi)
    if (allocated(dNi)) deallocate(dNi)
    if (allocated(lAsymmetric1)) deallocate(lAsymmetric1)
    if (allocated(lAsymmetric2)) deallocate(lAsymmetric2)
    if (allocated(dgdc)) deallocate(dgdc)
    if (allocated(dMolDerivatives)) deallocate(dMolDerivatives)
    if (allocated(dSumY)) deallocate(dSumY)
    nPhaseElements = nSub1 + nSub2
    allocate(dXi(nPhaseElements),dYi(nPhaseElements),dNi(nPhaseElements))
    allocate(lAsymmetric1(nPhaseElements))
    allocate(lAsymmetric2(nPhaseElements))
    allocate(dgdc(nPhaseElements))
    allocate(dMolDerivatives(nA2X2))
    allocate(dSumY(2))

    ! Initialize variables:
    dSiteFraction(iSPI,1:2,1:nMaxConstituentSys) = 0D0
    dXi                               = 0D0
    dYi                               = 0D0
    dNi                               = 0D0
    dChemicalPotential(iFirst:iLast)  = 0D0
    dPartialExcessGibbs(iFirst:iLast) = 0D0
    dgdc                              = 0D0
    dMolDerivatives                   = 0D0

    ! Compute X_i and Y_i
    dSum1 = 0D0
    dSum2 = 0D0
    dSumY = 0D0
    do i = iFirst, iLast
        k = i + 1 - iFirst
        a = iConstituentSublattice(iSPI,1,k)
        x = iConstituentSublattice(iSPI,2,k) + nSub1
        dNi(a)   = dNi(a)   + dMolFraction(i) * dConstituentCoefficients(iSPI,k,1)
        dSum1    = dSum1    + dMolFraction(i) * dConstituentCoefficients(iSPI,k,1)
        dSumY(1) = dSumY(1) + dMolFraction(i) * dConstituentCoefficients(iSPI,k,1) * dSublatticeCharge(iSPI,1,a)
        dNi(x)   = dNi(x)   + dMolFraction(i) * dConstituentCoefficients(iSPI,k,2)
        dSum2    = dSum2    + dMolFraction(i) * dConstituentCoefficients(iSPI,k,2)
        dSumY(2) = dSumY(2) + dMolFraction(i) * dConstituentCoefficients(iSPI,k,2) * dSublatticeCharge(iSPI,2,x-nSub1)
    end do

    ! q and p are charge-weighted sums on each sublattice (used to match SUBI)
    q = 0
    do i = 1, nSub1
        dXi(i) = dNi(i) / dSum1
        dYi(i) = dNi(i) * dSublatticeCharge(iSPI,1,i) / dSumY(1)
        q = q + dXi(i) * dSublatticeCharge(iSPI,1,i)
        dSiteFraction(iSPI,1,i) = dXi(i)
    end do
    p = 0
    do i = 1, nSub2
        k = nSub1 + i
        dXi(k) = dNi(k) / dSum2
        dYi(k) = dNi(k) * dSublatticeCharge(iSPI,2,i) / dSumY(2)
        p = p + dXi(k) * dSublatticeCharge(iSPI,2,i)
        dSiteFraction(iSPI,2,i) = dXi(k)
    end do

    ! Model enforces random mixing: no ordering allowed!
    ! Must set X_a,x = X_a * X_x
    do i = iFirst, iLast
        k = i + 1 - iFirst
        a = iConstituentSublattice(iSPI,1,k)
        x = iConstituentSublattice(iSPI,2,k) + nSub1
        dMolFraction(i) = dXi(a) * dXi(x)
    end do

    ! Compute number of moles and its derivatives
    dMol = q + p
    dMolAtoms = 0D0
    do j = iFirst, iLast
        ! Relative species index:
        n = j - iFirst + 1

        ! Store constituent indices:
        a = iConstituentSublattice(iSPI,1,n)
        x = iConstituentSublattice(iSPI,2,n)

        ! Allocating the correct constituent charges
        lc1 = dSublatticeCharge(iSPI,1,a)
        lc2 = dSublatticeCharge(iSPI,2,x)

        dMolDerivatives(n) = (q-lc1)*lc2/(dSum1*dMol**2)
        dMolDerivatives(n) = dMolDerivatives(n) + (p-lc2)*lc1/(dSum2*dMol**2)

        dMolAtoms = dMolAtoms + dMolFraction(j) * (dSublatticeCharge(iSPI,1,a) + dSublatticeCharge(iSPI,2,x))
    end do

    ! ---------------------------------------------------------------
    ! COMPUTE REFERENCE GIBBS ENERGY AND IDEAL MIXING TERMS
    ! ---------------------------------------------------------------

    gref = 0D0
    do j = iFirst, iLast
        ! Relative species index:
        n = j - iFirst + 1

        ! Store constituent indices:
        a = iConstituentSublattice(iSPI,1,n)
        x = iConstituentSublattice(iSPI,2,n) + nSub1
        
        gref = gref + dXi(a) * dXi(x) * dStdGibbsEnergy(j)
        do i = 1, nSub1
            dgdc(i) = dgdc(i) - dXi(a) * dXi(x) * dStdGibbsEnergy(j)
            if (i == a) dgdc(i) = dgdc(i) + dXi(x) * dStdGibbsEnergy(j)
        end do

        do i = 1, nSub2
            k = nSub1 + i
            dgdc(k) = dgdc(k) - dXi(a) * dXi(x) * dStdGibbsEnergy(j)
            if (k == x) dgdc(k) = dgdc(k) + dXi(a) * dStdGibbsEnergy(j)
        end do
    end do

    ! ---------------------------------------------------------------
    ! COMPUTE IDEAL MIXING ENERGY
    ! ---------------------------------------------------------------

    gideal = 0D0
    ! Calculate entropy derivatives on first sublattice
    do i = 1, nSub1
        gideal = gideal + p * dXi(i) * DLOG(dXi(i))
        dgdc(i) = dgdc(i) + p * (1 - dXi(i)) * DLOG(dXi(i))
        do j = 1, nSub1
            if (.NOT. (i == j)) then
                dgdc(i) = dgdc(i) - p * dXi(j) * DLOG(dXi(j))
            end if
        end do
        do j = 1, nSub2
            l = nSub1 + j
            dgdc(i) = dgdc(i) + (dSublatticeCharge(iSPI,1,i) - q) * dXi(l) * DLOG(dXi(l))
        end do
    end do
    ! Calculate entropy derivatives on second sublattice
    do i = 1, nSub2
        k = nSub1 + i
        gideal = gideal + q * dXi(k) * DLOG(dXi(k))
        dgdc(k) = dgdc(k) + q * (1 - dXi(k)) * DLOG(dXi(k))
        do j = 1, nSub2
            l = nSub1 + j
            if (.NOT. (i == j)) then 
                dgdc(k) = dgdc(k) - q * dXi(l) * DLOG(dXi(l))
            end if
        end do
        do j = 1, nSub1
            dgdc(k) = dgdc(k) + (dSublatticeCharge(iSPI,2,i) - p) * dXi(j) * DLOG(dXi(j))
        end do
    end do

    !---------------------------------------------------------------
    ! COMPUTE EXCESS GIBBS ENERGY
    ! --------------------------------------------------------------

    gex = 0D0
    ! Loop through parameters:
    LOOP_Param: do l = nParamPhase(iSolnIndex-1) + 1, nParamPhase(iSolnIndex)
        if (dExcessGibbsParam(l) == 0D0) cycle LOOP_Param

        ! Get constituents and exponents from parameter
        order = iRegularParam(l,1)          ! Mixing order (3 for binary, 4 for ternary)
        a = iRegularParam(l,2)              ! Index of A
        b = iRegularParam(l,3)              ! Index of B
        xx = iRegularParam(l,order + 1)     ! Index of X, unadjusted
        ! x = xx - nSub1                      ! Index of X

        ea = iRegularParam(l,order + 2)     ! Exponent of a
        eb = iRegularParam(l,order + 3)     ! Exponent of b
        ! ex = iRegularParam(l,order*2 + 1)   ! Exponent of x
        ex = 1D0                            ! Exponent of x (seems to have to be 1)

        c = -1                              ! Index of C (default to -1 so references will break)
        yc = 0D0                            ! Charge-equivalent site fraction of constituent c (default to 0)
        ec = 0D0                            ! Exponent of C (default to 0)
        dYcfac = 1D0                        ! yc**ec, define explicitly for 0**0 = 1
        if (order == 4) then
            ! Ternary mixing, must define parameters for constituent C
            c = iRegularParam(l,4)          ! Index of C
            yc = dYi(c)                     ! Charge-equivalent site fraction of constituent c 
            ec = iRegularParam(l,8)         ! Exponent of C
            dYcfac = yc**ec                 ! yc**ec
        end if

        ! Check which sublattice mixing is on, try to sort out all differences between 1 and 2 here
        if (a <= nSub1) then
            ! Mixing on first sublattice
            iSub      = 1
            iOffset   = 0
            iStartCon = 1
            iEndCon   = nSub1
        else
            ! Mixing on second sublattice
            iSub      = 2
            iOffset   = nSub1
            iStartCon = nSub1 + 1
            iEndCon   = nSub1 + nSub2
        end if

        ! Get chemical groups of first two constituents
        iGroup1 = iChemicalGroup(iSPI,iSub,a-iOffset)
        iGroup2 = iChemicalGroup(iSPI,iSub,b-iOffset)            

        ! Calculate symmetry of all species with these first two
        ! Count them as asymmetric with respect to themselves
        lAsymmetric1 = .FALSE.
        lAsymmetric2 = .FALSE.
        lAsymmetric1(a) = .TRUE.
        lAsymmetric2(b) = .TRUE.
        LOOP_checkSymmetry: do i = iStartCon, iEndCon
            ! First check if this ternary is an exception
            lIsException = .FALSE.
            LOOP_overrides: do k = 1, nInterpolationOverride(iSolnIndex)
                ! Check if this override applies to this ternary
                do j = 1, 3
                    if (.NOT.((iInterpolationOverride(iSolnIndex,k,j) == a) .OR. &
                              (iInterpolationOverride(iSolnIndex,k,j) == b) .OR. &
                              (iInterpolationOverride(iSolnIndex,k,j) == i))) &
                        cycle LOOP_overrides
                end do
                ! If we get here, this one is an exception
                lIsException = .TRUE.
                if (iInterpolationOverride(iSolnIndex,k,5) == b) lAsymmetric1(i) = .TRUE.
                if (iInterpolationOverride(iSolnIndex,k,5) == a) lAsymmetric2(i) = .TRUE.
                exit LOOP_overrides
            end do LOOP_overrides

            if (lIsException) cycle LOOP_checkSymmetry
            ! If they are equal, everything else will be symmetric
            if (iGroup1 /= iGroup2) then
                iGroupTemp = iChemicalGroup(iSPI,iSub,i-iOffset)
                if      (iGroupTemp == iGroup1) then
                    lAsymmetric1(i) = .TRUE.
                else if (iGroupTemp == iGroup2) then
                    lAsymmetric2(i) = .TRUE.
                end if
            end if
        end do LOOP_checkSymmetry

        ! Now calculate xis and sigma
        dXi1 = 0D0
        dXi2 = 0D0
        do i = iStartCon, iEndCon
            if      (lAsymmetric1(i)) then
                dXi1 = dXi1 + dYi(i)
            else if (lAsymmetric2(i)) then
                dXi2 = dXi2 + dYi(i)
            end if
        end do

        dXiDen = dXi1 + dXi2

        ! Same as QKTO
        gexTemp = dExcessGibbsParam(l) * (dXi1**(ea - 1)) * (dXi2**(eb - 1)) / (dXiDen ** (ea + eb - 2))
        gexTemp = gexTemp * dYi(a) * dYi(b) * dSumY(iSub)
        ! Include factor from opposite sublattice
        gexTemp = gexTemp * dYi(xx)
        gex = gex + gexTemp

        ! Contribute to derivatives with respect to constituents
        ! Derivatives on mixing sublattice
        do i = iStartCon, iEndCon
            dgexTemp = 0D0
            ! I can't explain the following line, but hey, it works
            dgexTemp = dgexTemp - (2D0 + dSumY(iSub) - dSumY(iSub)**2D0) / dSumY(iSub)
            if      (lAsymmetric1(i)) then
                dgexTemp = dgexTemp + (ea - 1) / dXi1 - (ea + eb - 2) / dXiDen
            else if (lAsymmetric2(i)) then
                dgexTemp = dgexTemp + (eb - 1) / dXi2 - (ea + eb - 2) / dXiDen
            end if
            if      (i == a) then
                dgexTemp = dgexTemp + 1D0 / dYi(a) / dSumY(iSub)
            else if (i == b) then
                dgexTemp = dgexTemp + 1D0 / dYi(b) / dSumY(iSub)
            else if (i == c) then
                dgexTemp = dgexTemp + 1D0 / dYi(c) / dSumY(iSub)
            end if
            ! Multiply by energy of term and charge of constituent
            dgdc(i) = dgdc(i) + dgexTemp * gexTemp * dSublatticeCharge(iSPI,iSub,i-iOffset)
        end do
        ! Derivatives on opposite sublattice
        if (iSub == 1) then
            do i = 1, nSub2
                k = nSub1 + i
                dgdc(k) = dgdc(k) - ex * dSublatticeCharge(iSPI,2,i) * dSumY(2) * gexTemp
                if (k == xx) dgdc(k) = dgdc(k) + ex * dSublatticeCharge(iSPI,2,i) / dYi(xx) * gexTemp
            end do
        else
            do i = 1, nSub1
                dgdc(i) = dgdc(i) - ex * dSublatticeCharge(iSPI,1,i) * dSumY(1) * gexTemp
                if (k == xx) dgdc(i) = dgdc(i) + ex * dSublatticeCharge(iSPI,1,i) / dYi(xx) * gexTemp
            end do
        end if
    end do LOOP_Param


    ! ---------------------------------------------------------------
    ! SUM ENERGY CONTRIBUTIONS TO CHEMICAL POTENTIALS
    ! ---------------------------------------------------------------

    do i = iFirst, iLast
        ! Get constituents
        k = i + 1 - iFirst
        a = iConstituentSublattice(iSPI,1,k)
        x = iConstituentSublattice(iSPI,2,k) + nSub1

        natom = (dSublatticeCharge(iSPI,1,a) + dSublatticeCharge(iSPI,2,x-nSub1)) / dMol + dMolDerivatives(k) * dMolAtoms

        dChemicalPotential(i) = dChemicalPotential(i) + (gideal + gref + gex) * natom

        ! Calculate entropic contributions from derivatives
        do j = 1, nSub1
            if (j == a) then
                dChemicalPotential(i) = dChemicalPotential(i) + dConstituentCoefficients(iSPI,k,1) * dgdc(j)
            end if
        end do

        do j = 1, nSub2
            l = nSub1 + j
            if (l == x) then
                dChemicalPotential(i) = dChemicalPotential(i) + dConstituentCoefficients(iSPI,k,2) * dgdc(l)
            end if
        end do
    end do

    ! Deallocate allocatable arrays:
    deallocate(dXi,dYi,dNi,lAsymmetric1,lAsymmetric2)
    deallocate(dgdc,dMolDerivatives)

    return

end subroutine CompExcessGibbsEnergySUBM
