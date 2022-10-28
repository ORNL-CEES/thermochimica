
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
    !! 2-sublattice solution phase model.
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
    integer :: a, x, b, y!c, abxy, xx
    integer :: iSolnIndex, iSPI, nPhaseElements, nSub1, nSub2, nA2X2
    integer :: iFirst, iLast
    logical, allocatable, dimension(:) :: lAsymmetric1, lAsymmetric2
    ! logical :: lIsException
    real(8) :: dConfEntropy, dRef, dSum1, dSum2, dSumY1, dSumY2, dCWS1, dCWS2, dTemp!, p, q, r
    real(8), allocatable, dimension(:) :: dXi, dYi, dNi, dgdc

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
    nPhaseElements = nSub1 + nSub2
    allocate(dXi(nPhaseElements),dYi(nPhaseElements),dNi(nPhaseElements))
    allocate(lAsymmetric1(MAX(nSub1,nSub2)))
    allocate(lAsymmetric2(MAX(nSub1,nSub2)))
    allocate(dgdc(nPhaseElements))

    ! Initialize variables:
    dXi                               = 0D0
    dYi                               = 0D0
    dNi                               = 0D0
    dChemicalPotential(iFirst:iLast)  = 0D0
    dPartialExcessGibbs(iFirst:iLast) = 0D0
    dgdc                              = 0D0

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
    ! CWS = charge-weighted sum (needed in entropy calculation)
    dCWS1 = 0
    do i = 1, nSub1
        dXi(i) = dNi(i) / dSum1
        dYi(i) = dNi(i) * dSublatticeCharge(iSPI,1,i) / dSumY1
        dCWS1 = dCWS1 + dXi(i) * dSublatticeCharge(iSPI,1,i)
    end do
    dCWS2 = 0
    do i = 1, nSub2
        k = nSub1 + i
        dXi(k) = dNi(k) / dSum2
        dYi(k) = dNi(k) * dSublatticeCharge(iSPI,2,i) / dSumY2
        dCWS2 = dCWS2 + dXi(k) * dSublatticeCharge(iSPI,2,i)
    end do

    ! Model enforces random mixing: no ordering allowed!
    ! Must set X_a,x = X_a * X_x
    do i = iFirst, iLast
        k = i + 1 - iFirst
        a = iConstituentSublattice(iSPI,1,k)
        x = iConstituentSublattice(iSPI,2,k) + nSub1
        dMolFraction(i) = dXi(a) * dXi(x)
    end do

    ! ---------------------------------------------------------------
    ! COMPUTE REFERENCE GIBBS ENERGY AND IDEAL MIXING TERMS
    ! ---------------------------------------------------------------

    dConfEntropy = 0D0
    ! Calculate entropy derivatives on first sublattice
    do i = 1, nSub1
        dgdc(i) = dgdc(i) + (1 - dXi(i)) + LOG(dXi(i))
        do j = 1, nSub1
            if (.NOT. (i == j)) then
                dgdc(i) = dgdc(i) - dXi(j)
            end if
        end do
    end do
    ! Calculate entropy derivatives on second sublattice
    do i = 1, nSub2
        k = nSub1 + i
        dgdc(k) = dgdc(k) + (1 - dXi(k)) + LOG(dXi(k))
        do j = 1, nSub2
            l = nSub1 + j
            if (.NOT. (i == j)) then 
                dgdc(k) = dgdc(k) - dXi(k)
            end if
        end do
    end do

    do i = iFirst, iLast
        ! Get constituents
        k = i + 1 - iFirst
        a = iConstituentSublattice(iSPI,1,k)
        x = iConstituentSublattice(iSPI,2,k) + nSub1

        ! Loop through species to get reference energy contributions:
        LOOP_species: do j = iFirst, iLast
            ! Relative species index:
            n = j + 1 - iFirst
            b = iConstituentSublattice(iSPI,1,n)
            y = iConstituentSublattice(iSPI,2,n) + nSub1

            ! Compute pre-factor term:
            dTemp = -1D0
            if (a == b) dTemp = dTemp + 1D0 / dXi(a)
            if (x == y) dTemp = dTemp + 1D0 / dXi(x)

            ! Update the reference molar Gibbs energy:
            dChemicalPotential(i) = dChemicalPotential(i) + dTemp * dMolFraction(j) * dStdGibbsEnergy(j)
        end do LOOP_species

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

    ! Loop through excess mixing parameters:
    ! LOOP_Param: do abxy = nParamPhase(iSolnIndex-1) + 1, nParamPhase(iSolnIndex)

    !     if (dExcessGibbsParam(abxy) == 0D0) cycle LOOP_Param

    !     ! AB/XY parametrization
    !     n = iRegularParam(abxy,1)
    !     a = iRegularParam(abxy,2)              ! Index of A
    !     b = iRegularParam(abxy,3)              ! Index of B
    !     xx = iRegularParam(abxy,n+1)           ! Index of X, unadjusted
    !     x = xx - nSub1                         ! Index of X
    !     p = iRegularParam(abxy,n+2)            ! Exponent 1
    !     q = iRegularParam(abxy,n+3)            ! Exponent 2
    !     if (n == 4) then
    !         c = iRegularParam(abxy,4)          ! Index of C
    !         r = iRegularParam(abxy,8)          ! Exponent 3
    !     end if


    !     lAsymmetric1 = .FALSE.
    !     lAsymmetric2 = .FALSE.
    !     lAsymmetric1(a) = .TRUE.
    !     lAsymmetric2(b) = .TRUE.
    !     ! First check if this ternary is an exception
    !     LOOP_checkSymmetry: do i = 1, nSub1
    !         lIsException = .FALSE.
    !         LOOP_overrides: do k = 1, nInterpolationOverride(iSolnIndex)
    !             ! Check if this override applies to this ternary
    !             do l = 1, 3
    !                 if (.NOT.((iInterpolationOverride(iSolnIndex,k,l) == a) .OR. &
    !                         (iInterpolationOverride(iSolnIndex,k,l) == b) .OR. &
    !                         (iInterpolationOverride(iSolnIndex,k,l) == i))) &
    !                     cycle LOOP_overrides
    !             end do
    !             ! If we get here, this one is an exception
    !             lIsException = .TRUE.
    !             if (iInterpolationOverride(iSolnIndex,k,5) == b) lAsymmetric1(i) = .TRUE.
    !             if (iInterpolationOverride(iSolnIndex,k,5) == a) lAsymmetric2(i) = .TRUE.
    !             exit LOOP_overrides
    !         end do LOOP_overrides

    !         if (lIsException) cycle LOOP_checkSymmetry
    !         ! First make a list of which constituents make asymmetric ternaries
    !         if (iChemicalGroup(iSPI,1,a) /= iChemicalGroup(iSPI,1,b)) then
    !                 if (iChemicalGroup(iSPI,1,i) == iChemicalGroup(iSPI,1,a)) then
    !                     lAsymmetric1(i) = .TRUE.
    !                 else if (iChemicalGroup(iSPI,1,i) == iChemicalGroup(iSPI,1,b)) then
    !                     lAsymmetric2(i) = .TRUE.
    !                 end if
    !         end if
    !     end do LOOP_checkSymmetry

    ! end do LOOP_Param

    ! Deallocate allocatable arrays:
    deallocate(dXi,dYi,dNi,lAsymmetric1,lAsymmetric2)
    deallocate(dgdc)

    return

end subroutine CompExcessGibbsEnergySUBM
