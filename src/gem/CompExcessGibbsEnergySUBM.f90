
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
    integer :: a, b, c, x, xx, order
    integer :: iSolnIndex, iSPI, nPhaseElements, nSub1, nSub2, nA2X2
    integer :: iFirst, iLast
    logical, allocatable, dimension(:) :: lAsymmetric1, lAsymmetric2
    real(8) :: dSum1, dSum2, dSumY1, dSumY2, q, p, ea, eb, ec, ex
    real(8) :: gref, gideal, gex, natom, dMol, lc1, lc2, dMolAtoms
    real(8), allocatable, dimension(:) :: dXi, dYi, dNi, dgdc, dMolDerivatives

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
    allocate(dMolDerivatives(nA2X2))
    nPhaseElements = nSub1 + nSub2
    allocate(dXi(nPhaseElements),dYi(nPhaseElements),dNi(nPhaseElements))
    allocate(lAsymmetric1(MAX(nSub1,nSub2)))
    allocate(lAsymmetric2(MAX(nSub1,nSub2)))
    allocate(dgdc(nPhaseElements))

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
    dSumY1 = 0D0
    dSumY2 = 0D0
    do i = iFirst, iLast
        k = i + 1 - iFirst
        a = iConstituentSublattice(iSPI,1,k)
        x = iConstituentSublattice(iSPI,2,k) + nSub1
        dNi(a) = dNi(a) + dMolFraction(i) * dConstituentCoefficients(iSPI,k,1)
        dSum1 = dSum1 + dMolFraction(i) * dConstituentCoefficients(iSPI,k,1)
        dSumY1 = dSumY1 + dMolFraction(i) * dConstituentCoefficients(iSPI,k,1) * dSublatticeCharge(iSPI,1,a)
        dNi(x) = dNi(x) + dMolFraction(i) * dConstituentCoefficients(iSPI,k,2)
        dSum2 = dSum2 + dMolFraction(i) * dConstituentCoefficients(iSPI,k,2)
        dSumY2 = dSumY2 + dMolFraction(i) * dConstituentCoefficients(iSPI,k,2) * dSublatticeCharge(iSPI,2,x-nSub1)
    end do

    ! q and p are charge-weighted sums on each sublattice (used to match SUBI)
    q = 0
    do i = 1, nSub1
        dXi(i) = dNi(i) / dSum1
        dYi(i) = dNi(i) * dSublatticeCharge(iSPI,1,i) / dSumY1
        q = q + dXi(i) * dSublatticeCharge(iSPI,1,i)
        dSiteFraction(iSPI,1,i) = dXi(i)
    end do
    p = 0
    do i = 1, nSub2
        k = nSub1 + i
        dXi(k) = dNi(k) / dSum2
        dYi(k) = dNi(k) * dSublatticeCharge(iSPI,2,i) / dSumY2
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
        x = xx - nSub1                      ! Index of X

        ea = iRegularParam(l,order + 2)     ! Exponent of a
        eb = iRegularParam(l,order + 3)     ! Exponent of b
        ex = iRegularParam(l,order*2 + 1)   ! Exponent of b

        if (order == 4) then
            c = iRegularParam(l,4)          ! Index of C
            ec = iRegularParam(l,8)         ! Index of C
        end if

        gex = gex + dExcessGibbsParam(l) * dYi(a)**ea * dYi(b)**eb * dYi(xx)
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
