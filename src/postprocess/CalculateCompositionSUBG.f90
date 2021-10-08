module ModulePrintUtils
contains
subroutine CalculateCompositionSUBG(iSolnIndex,dMolesPairs,lPrint,cPair,dPair)

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, l, m
    integer :: iSolnIndex, iSPI, nPhaseElements
    integer :: iFirst, iLast, nSub1, nSub2, iMax, nMax
    real(8) :: dSum, dMax, dSumElementQuads, dSumElementPairs,dMolesPairs
    real(8) :: dZa, dZb, dZx, dZy!, dXtot, dYtot
    real(8), allocatable, dimension(:) :: dXi, dYi, dNi
    real(8), allocatable, dimension(:,:) :: dXij, dNij
    logical :: lPrint
    character(30), dimension(*), intent(out),optional :: cPair
    real(8), dimension(*), intent(out),optional       :: dPair
    character(30)                           :: cDummyB
    character(2)                            :: cDummy
    ! X_ij/kl corresponds to dMolFraction


    ! Only proceed if the correct phase type is selected:
    if (.NOT. (cSolnPhaseType(iSolnIndex) == 'SUBG' .OR. cSolnPhaseType(iSolnIndex) == 'SUBQ')) return

    ! Define temporary variables for sake of convenience:
    iFirst = nSpeciesPhase(iSolnIndex-1) + 1
    iLast  = nSpeciesPhase(iSolnIndex)
    iSPI = iPhaseSublattice(iSolnIndex)
    nSub1 = nConstituentSublattice(iSPI,1)
    nSub2 = nConstituentSublattice(iSPI,2)

    ! Allocate allocatable arrays:
    if (allocated(dXi)) deallocate(dXi)
    if (allocated(dYi)) deallocate(dYi)
    if (allocated(dNi)) deallocate(dNi)
    if (allocated(dXij)) deallocate(dXij)
    if (allocated(dNij)) deallocate(dNij)
    j = iLast - iFirst + 1
    nPhaseElements = nSub1 + nSub2
    allocate(dXi(nPhaseElements),dYi(nPhaseElements),dNi(nPhaseElements))
    allocate(dXij(nSub1,nSub2))
    allocate(dNij(nSub1,nSub2))

    ! Initialize variables:
    dXi                               = 0D0
    dYi                               = 0D0
    dNi                               = 0D0
    dXij                              = 0D0
    dNij                              = 0D0

    ! Compute X_i and Y_i
    ! Do cations first:
    dSum = 0D0
    do i = 1, nSub1
        do k = 1, nPairsSRO(iSPI,2)
            l = iFirst + k - 1
            dZa = dCoordinationNumber(iSPI,k,1)
            dZb = dCoordinationNumber(iSPI,k,2)
            if (i == iPairID(iSPI,k,1))  then
                dNi(i) = dNi(i) + (dMolFraction(l) / dZa)
                dYi(i) = dYi(i) + (dMolFraction(l) / 2)
            end if
            if (i == iPairID(iSPI,k,2))  then
                dNi(i) = dNi(i) + (dMolFraction(l) / dZb)
                dYi(i) = dYi(i) + (dMolFraction(l) / 2)
            end if
        end do
        dSum = dSum + dNi(i)
    end do
    if (lPrint) print *, '   Cation fractions:'
    do i = 1, nSub1
        dXi(i) = dNi(i) / dSum
    end do
    nMax = 0
    if (lPrint) then
        do i = 1, nSub1
            cDummyB  = TRIM(cConstituentNameSUB(iSPI,1,i))
            k    = LEN_TRIM(cConstituentNameSUB(iSPI,1,i)) - 1
            nMax = MAX(k, nMax)
            if (i == 1) then
                cDummy = '{ '
            else
                cDummy = '+ '
            end if
            if (i == nSub1) cDummyB(nMax+2:nMax+3) = ' }'
            if (dXi(i) >= 1D-1) then
                print '(A20,F7.5,A3,A35)', cDummy, dXi(i), ' ', cDummyB
            else
                print '(A20,ES10.4,A35)', cDummy, dXi(i), cDummyB
            end if
        end do
    end if

    ! Do anions now:
    dSum = 0D0
    do i = 1, nSub2
        j = i + nSub1
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

    if (lPrint) print *, '   Anion fractions:'
    do i = 1, nSub2
        j = i + nSub1
        dXi(j) = dNi(j) / dSum
    end do
    if (lPrint) then
        do i = 1, nSub2
            j = i + nSub1
            cDummyB  = TRIM(cConstituentNameSUB(iSPI,2,i))
            k    = LEN_TRIM(cConstituentNameSUB(iSPI,2,i)) - 1
            nMax = MAX(k, nMax)
            if (i == 1) then
                cDummy = '{ '
            else
                cDummy = '+ '
            end if
            if (i == nSub2) cDummyB(nMax+2:nMax+3) = ' }'
            if (dXi(i) >= 1D-1) then
                print '(A20,F7.5,A3,A35)', cDummy, dXi(j), ' ', cDummyB
            else
                print '(A20,ES10.4,A35)', cDummy, dXi(j), cDummyB
            end if
        end do
    end if

    dSum = 0D0
    do m = 1, nPairsSRO(iSPI,1)
        i = iConstituentSublattice(iSPI,1,m)
        j = iConstituentSublattice(iSPI,2,m)
        do k = 1, nPairsSRO(iSPI,2)
            l = iFirst + k - 1
            dZa = dCoordinationNumber(iSPI,k,1)
            dZb = dCoordinationNumber(iSPI,k,2)
            if ((i == iPairID(iSPI,k,1)) .AND. ((j + nSub1) == iPairID(iSPI,k,3)))  then
                dNij(i,j) = dNij(i,j) + (dMolFraction(l) / dZa) / dConstituentCoefficients(iSPI,m,1)
            end if
            if ((i == iPairID(iSPI,k,1)) .AND. ((j + nSub1) == iPairID(iSPI,k,4)))  then
                dNij(i,j) = dNij(i,j) + (dMolFraction(l) / dZa) / dConstituentCoefficients(iSPI,m,1)
            end if
            if ((i == iPairID(iSPI,k,2)) .AND. ((j + nSub1) == iPairID(iSPI,k,3)))  then
                dNij(i,j) = dNij(i,j) + (dMolFraction(l) / dZb) / dConstituentCoefficients(iSPI,m,1)
            end if
            if ((i == iPairID(iSPI,k,2)) .AND. ((j + nSub1) == iPairID(iSPI,k,4)))  then
                dNij(i,j) = dNij(i,j) + (dMolFraction(l) / dZb) / dConstituentCoefficients(iSPI,m,1)
            end if
        end do
        dSum = dSum + dNij(i,j)
    end do

    ! Use most abundant element in phase to normalize
    dMax = 0D0
    do i = 1, nElements
        dSumElementQuads = 0D0
        do k = 1, nPairsSRO(iSPI,2)
            l = iFirst + k - 1
            dSumElementQuads = dSumElementQuads + dStoichSpecies(l,i)*dMolesSpecies(l)
        end do
        if (dSumElementQuads > dMax) then
            dMax = dSumElementQuads
            iMax = i
        end if
    end do
    dSumElementQuads = dMax

    dSumElementPairs = 0D0
    do m = 1, nPairsSRO(iSPI,1)
        i = iConstituentSublattice(iSPI,1,m)
        j = iConstituentSublattice(iSPI,2,m)
        dSumElementPairs = dSumElementPairs + dStoichPairs(iSPI,m,iMax)*dNij(i,j)
    end do

    dMolesPairs = dSum*dSumElementQuads/dSumElementPairs
    if (lPrint) print *
    if ((dMolesPairs >= 999.95).OR.(dMolesPairs <= 1D-1)) then
        if (lPrint) print '(A4,ES10.4,A15)', '    ', dMolesPairs, ' Moles of pairs'
    else
        if (lPrint) print '(A7,F7.5,A15)', '       ', dMolesPairs, ' Moles of pairs'
    end if
    if (lPrint) print *, '   Pair fractions:'
    do m = 1, nPairsSRO(iSPI,1)
        i = iConstituentSublattice(iSPI,1,m)
        j = iConstituentSublattice(iSPI,2,m)
        dXij(i,j) = dNij(i,j) / dSum
        if (PRESENT(cPair)) then
            cPair(m) = ADJUSTL(cPairName(iSPI,m))
        end if
        if (PRESENT(dPair)) dPair(m) = dXij(i,j)
    end do


    nMax = 0
    if (lPrint) then
        ! First pair:
        m = 1
        i = iConstituentSublattice(iSPI,1,m)
        j = iConstituentSublattice(iSPI,2,m)
        if (dXij(i,j) >= 1D-1) then
            print '(A20,F7.5,A3,A35)', '{ ', dXij(i,j), ' ', cPairName(iSPI,m)
        else
            print '(A20,ES10.4,A35)', '{ ', dXij(i,j), cPairName(iSPI,m)
        end if
        k    = LEN_TRIM(cPairName(iSPI,m)) - 1
        nMax = MAX(k, nMax)

        ! Print middle pairs:
        do m = 2, nPairsSRO(iSPI,1) - 1
            i = iConstituentSublattice(iSPI,1,m)
            j = iConstituentSublattice(iSPI,2,m)
            if (dXij(i,j) >= 1D-1) then
                print '(A20,F7.5,A3,A35)', '+ ', dXij(i,j), ' ', cPairName(iSPI,m)
            else
                print '(A20,ES10.4,A35)', '+ ', dXij(i,j), cPairName(iSPI,m)
            end if
            k        = LEN_TRIM(cPairName(iSPI,m)) - 1
            nMax = MAX(k, nMax)
        end do

        ! Print last pair:
        m = nPairsSRO(iSPI,1)
        i = iConstituentSublattice(iSPI,1,m)
        j = iConstituentSublattice(iSPI,2,m)
        k        = LEN_TRIM(cPairName(iSPI,m)) - 1
        nMax     = MAX(k, nMax) + 1
        nMax     = MIN(nMax,27)
        cDummyB  = TRIM(cPairName(iSPI,m))
        cDummyB(nMax+2:nMax+3) = '}'

        if (dXij(i,j) >= 1D-1) then
            print '(A20,F7.5,A3,A35)', '+ ', dXij(i,j), ' ', cDummyB
        else
            print '(A20,ES10.4,A35)', '+ ', dXij(i,j), cDummyB
        end if
        print *
    end if

    ! Deallocate allocatable arrays:
    deallocate(dXi,dYi,dNi,dXij,dNij)

end subroutine CalculateCompositionSUBG
end module ModulePrintUtils
