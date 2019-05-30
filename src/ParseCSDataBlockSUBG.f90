

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseCSDataBlockSUBG.f90
    !> \brief   Parse the data block section corresponding to a SUBG phase of a ChemSage data-file.
    !> \author  M.H.A. Piro
    !> \date    Mar. 4, 2018
    !> \sa      ParseCSDataFile.f90
    !> \sa      ParseCSDataBlock.f90
    !> \sa      ParseCSDataBlockGibbs.f90
    !> \todo    There are a number of lines in SUBG phases that I do not yet understand.
    !!           I've asked some experts and they don't know either, which tells me that
    !!           they're not important. Once I
    !!           gain more experience with these models, this will likely become more clear.
    !
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified and is completely
    ! independent of ChemApp and related products, including Solgas, Solgasmix, Fact, FactSage
    ! and ChemSage.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !   03/04/2018      M.H.A. Piro     Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to parse the "data block" section of a ChemSage data-file
    !! containing a "SUBG" phase, which represents the modified quasichemical model. This phase differs
    !! from many other types of thermodynamic models in that it attempts to capture Short Range Order (SRO)
    !! in liquid or solid solutions. This is achieved by focusing on pairs of species, rather than the species
    !! themselves. For more information, see the following paper:
    !!
    !! A.D. Pelton, S.A. Degterov, G. Eriksson, C. Roberlin, Y. Dessureault, "The Modified Quasichemical
    !! Model I -- Binary Solutions", Metallurgical and Materials Transactions B, 31B (2000) 651-659.
    !!
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO                      A scalar integer that indicates a successful exit or identifies an error.
    ! nSpeciesCS                Number of species in the system (combined solution species and pure
    !                            separate phases).
    ! nGibbsEqSpecies           Number of Gibbs energy equations for a particular species.
    ! iSpeciesAtomsCS           Integer matrix representing the number of atoms of a particular
    !                            elements in a species.
    ! iParticlesPerMoleCS       An integer vector containing the number of particles per mole of the
    !                            constituent species formula mass.  The default value is 1.
    ! cSolnPhaseNameCS          The name of a solution phase.
    ! cSolnPhaseTypeCS          The type of a solution phase.
    ! cSolnPhaseTypeSupport     A character array representing solution phase types that are supported.
    ! iRegularParamCS           An integer matrix representing the parameter index for the first dimension
    !                            and the mixing terms on the second dimension.  For the second dimension, the
    !                            first coefficient indicates whether the parameter is a binary or ternary term (n),
    !                            the next n coefficients correspond to the constituent indices, and the last
    !                            coefficient corresponds to the exponent.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ParseCSDataBlockSUBG( i )

    USE ModuleParseCS

    implicit none

    integer                     :: i, j, k, l, n, x, y, p, a, b, nA2X2, nChar
    integer                     :: ia2x2, ib2x2, ia2y2, iabx2, iaby2, ia2xy, ib2xy
    integer,     dimension(10)  :: iTempVec
    integer,     dimension(13)  :: iNumPos
    real(8)                     :: dAnionCoordTemp, dF, qa, qb, qx, qy
    real(8)                     :: dZAa2xy, dZXa2xy, dZYa2xy, dZBb2xy, dZXb2xy, dZYb2xy
    real(8)                     :: dZAabx2, dZBabx2, dZXabx2, dZAaby2, dZBaby2, dZYaby2
    real(8)                     :: dZAa2x2, dZBb2x2, dZXa2x2, dZXb2x2, dZAa2y2, dZYa2y2
    real(8),     dimension(20)  :: dTempVec
    character(8),dimension(20)  :: cDummyVec, cConstituentNames1, cConstituentNames2
    logical, dimension(:), allocatable :: lPairSet

    real(8), dimension(nSpeciesCS,nElementsCS) :: dStoichSpeciesOld

    ! Initialize variables:
    dTempVec = 0D0
    iTempVec = 0

    ! SUBG phases appear to be represented as multi-sublattice phases; however,
    ! they don't appear to make use of any sublattice information. I'm going to
    ! to read these lines for now, but it may need to be revised at a later time.

    ! This line contains N integers (where N is the number of sublattices)
    ! where each integer represents the number of constituents on the respective
    ! sublattice. I think there are always two sublattices for SUBG phases.
    read (1,*,IOSTAT = INFO) nSublatticeElementsCS(nCountSublatticeCS,1:2)
    nSublatticePhaseCS(nCountSublatticeCS) = 2


    ! Read in names of constituents on first sublattice:
    ! NOTE: THIS LINE MAY NEED TO BE REVISED IF THERE ARE A LARGE # OF CONSTITUENTS:
    read (1,*,IOSTAT = INFO) cConstituentNames1(1:nSublatticeElementsCS(nCountSublatticeCS,1))
    ! Match elements on 1st sublattice with elements in dat file order
    do k = 1, nSublatticeElementsCS(nCountSublatticeCS,1)
        ! Find numbers or +/- in name if they are there
        iNumPos = 3
        if (INDEX(cConstituentNames1(k),'1') > 0) iNumPos(1)  = INDEX(cConstituentNames1(k),'1')
        if (INDEX(cConstituentNames1(k),'2') > 0) iNumPos(2)  = INDEX(cConstituentNames1(k),'2')
        if (INDEX(cConstituentNames1(k),'3') > 0) iNumPos(3)  = INDEX(cConstituentNames1(k),'3')
        if (INDEX(cConstituentNames1(k),'4') > 0) iNumPos(4)  = INDEX(cConstituentNames1(k),'4')
        if (INDEX(cConstituentNames1(k),'5') > 0) iNumPos(5)  = INDEX(cConstituentNames1(k),'5')
        if (INDEX(cConstituentNames1(k),'6') > 0) iNumPos(6)  = INDEX(cConstituentNames1(k),'6')
        if (INDEX(cConstituentNames1(k),'7') > 0) iNumPos(7)  = INDEX(cConstituentNames1(k),'7')
        if (INDEX(cConstituentNames1(k),'8') > 0) iNumPos(8)  = INDEX(cConstituentNames1(k),'8')
        if (INDEX(cConstituentNames1(k),'9') > 0) iNumPos(9)  = INDEX(cConstituentNames1(k),'9')
        if (INDEX(cConstituentNames1(k),'0') > 0) iNumPos(10) = INDEX(cConstituentNames1(k),'0')
        if (INDEX(cConstituentNames1(k),'+') > 0) iNumPos(11) = INDEX(cConstituentNames1(k),'+')
        if (INDEX(cConstituentNames1(k),'-') > 0) iNumPos(12) = INDEX(cConstituentNames1(k),'-')
        nChar = MINVAL(iNumPos) - 1
        LOOP_Sublattice1Elements: do j = 1, nElementsCS
            if (cConstituentNames1(k)(1:nChar) == cElementNameCS(j)(1:2)) then
                iSublatticeElementsCS(nCountSublatticeCS, 1, k) = j
                exit LOOP_Sublattice1Elements
            end if
        end do LOOP_Sublattice1Elements
    end do

    ! Read in names of constituents on second sublattice: (ignore for now):
    read (1,*,IOSTAT = INFO) cConstituentNames2(1:nSublatticeElementsCS(nCountSublatticeCS,2))
    ! Match elements on 2nd sublattice with elements in dat file order
    do k = 1, nSublatticeElementsCS(nCountSublatticeCS,2)
        ! Find numbers or +/- in name if they are there
        iNumPos = 3
        if (INDEX(cConstituentNames2(k),'1') > 0) iNumPos(1)  = INDEX(cConstituentNames2(k),'1')
        if (INDEX(cConstituentNames2(k),'2') > 0) iNumPos(2)  = INDEX(cConstituentNames2(k),'2')
        if (INDEX(cConstituentNames2(k),'3') > 0) iNumPos(3)  = INDEX(cConstituentNames2(k),'3')
        if (INDEX(cConstituentNames2(k),'4') > 0) iNumPos(4)  = INDEX(cConstituentNames2(k),'4')
        if (INDEX(cConstituentNames2(k),'5') > 0) iNumPos(5)  = INDEX(cConstituentNames2(k),'5')
        if (INDEX(cConstituentNames2(k),'6') > 0) iNumPos(6)  = INDEX(cConstituentNames2(k),'6')
        if (INDEX(cConstituentNames2(k),'7') > 0) iNumPos(7)  = INDEX(cConstituentNames2(k),'7')
        if (INDEX(cConstituentNames2(k),'8') > 0) iNumPos(8)  = INDEX(cConstituentNames2(k),'8')
        if (INDEX(cConstituentNames2(k),'9') > 0) iNumPos(9)  = INDEX(cConstituentNames2(k),'9')
        if (INDEX(cConstituentNames2(k),'0') > 0) iNumPos(10) = INDEX(cConstituentNames2(k),'0')
        if (INDEX(cConstituentNames2(k),'+') > 0) iNumPos(11) = INDEX(cConstituentNames2(k),'+')
        if (INDEX(cConstituentNames2(k),'-') > 0) iNumPos(12) = INDEX(cConstituentNames2(k),'-')
        nChar = MINVAL(iNumPos) - 1
        LOOP_Sublattice2Elements: do j = 1, nElementsCS
            if (cConstituentNames2(k)(1:nChar) == cElementNameCS(j)(1:2)) then
                iSublatticeElementsCS(nCountSublatticeCS, 2, k) = j
                exit LOOP_Sublattice2Elements
            end if
        end do LOOP_Sublattice2Elements
    end do

    ! Read in the charge of each constituent on the first sublattice.
    ! This seems unnecessary so I'm going to ignore it for now:
    read (1,*,IOSTAT = INFO) dSublatticeChargeCS(nCountSublatticeCS,1,1:nSublatticeElementsCS(nCountSublatticeCS,1))

    ! I think that this entry represents the constituent IDs on the first sublattice (ignore for now):
    read (1,*,IOSTAT = INFO) dTempVec(1:nSublatticeElementsCS(nCountSublatticeCS,1))

    ! Read in the charge of each constituent on the second sublattice.
    ! This seems unnecessary so I'm going to ignore it for now:
    read (1,*,IOSTAT = INFO) dSublatticeChargeCS(nCountSublatticeCS,2,1:nSublatticeElementsCS(nCountSublatticeCS,2))

    ! I think that this entry represents the constituent IDs on the second sublattice (ignore for now):
    read (1,*,IOSTAT = INFO) dTempVec(1:nSublatticeElementsCS(nCountSublatticeCS,2))

    ! This entry appears to represent the IDs matching constituents on the first sublattice to species:
    nA2X2 = nSublatticeElementsCS(nCountSublatticeCS,1) * nSublatticeElementsCS(nCountSublatticeCS,2)
    read (1,*,IOSTAT = INFO) iConstituentSublatticeCS(nCountSublatticeCS, 1, 1:nA2X2)

    ! This entry appears to represent the IDs matching constituents on the second sublattice to species:
    read (1,*,IOSTAT = INFO) iConstituentSublatticeCS(nCountSublatticeCS, 2, 1:nA2X2)

    ! Set up default pair IDs and coordination numbers
    ! dCoordinationNumberCS(nCountSublatticeCS,1:nMaxSpeciesPhaseCS,1:4) = 6D0
    dCoordinationNumberCS(nCountSublatticeCS,1:nMaxSpeciesPhaseCS,1:4) = 0D0
    do y = 1, nSublatticeElementsCS(nCountSublatticeCS,2)
        LOOP_sroPairsOuter: do x = 1, nSublatticeElementsCS(nCountSublatticeCS,2)
            if (x == y) then
                p = (x - 1) * (nSublatticeElementsCS(nCountSublatticeCS,1) * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)
            else if (x > y) then
                cycle LOOP_sroPairsOuter
            else
                p = (nSublatticeElementsCS(nCountSublatticeCS,2) + (x - 1) + ((y-2)*(y-1)/2)) &
                  * (nSublatticeElementsCS(nCountSublatticeCS,1) * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)
            end if
            do k = 1, nSublatticeElementsCS(nCountSublatticeCS,1)
                LOOP_sroPairsInner: do j = 1, nSublatticeElementsCS(nCountSublatticeCS,1)
                    if (j == k) then
                        l = j
                    else if (j > k) then
                        cycle LOOP_sroPairsInner
                    else
                        l = nSublatticeElementsCS(nCountSublatticeCS,1) + j + ((k-2)*(k-1)/2)
                    end if
                    iPairIDCS(nCountSublatticeCS, l + p, 1) = j
                    iPairIDCS(nCountSublatticeCS, l + p, 2) = k
                    iPairIDCS(nCountSublatticeCS, l + p, 3) = x + nSublatticeElementsCS(nCountSublatticeCS,1)
                    iPairIDCS(nCountSublatticeCS, l + p, 4) = y + nSublatticeElementsCS(nCountSublatticeCS,1)
                    end do LOOP_sroPairsInner
            end do
        end do LOOP_sroPairsOuter
    end do

    ! Parse the co-ordination numbers corresponding to all pairs in the phase.
    ! Note that since these lines correspond to pairs, there will always be the same number of
    ! integers and reals on a line, but the number of lines corresponds to the number of pairs.
    ! The SUBG model considers quadruplets, which is why there are four sets.
    ! Note that a quadruplet must satisfy the following constraint:
    ! q(i)/Z(i) + q(j)/Z(j) =  q(x)/Z(x) + q(y)/Z(y)
    allocate(lPairSet(nSpeciesPhaseCS(i) - nSpeciesPhaseCS(i-1)))
    lPairSet = .FALSE.
    LOOP_readPairs: do n = 1, nPairsSROCS(nCountSublatticeCS,2)
        read (1,*,IOSTAT = INFO) j, k, x, y, dTempVec(1:4)
        x = x - nSublatticeElementsCS(nCountSublatticeCS,1)
        y = y - nSublatticeElementsCS(nCountSublatticeCS,1)
        if (x == y) then
            p = (x - 1) * (nSublatticeElementsCS(nCountSublatticeCS,1) * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)
        else if (x > y) then
            cycle LOOP_readPairs
        else
            p = (nSublatticeElementsCS(nCountSublatticeCS,2) + (x - 1) + ((y-2)*(y-1)/2)) &
              * (nSublatticeElementsCS(nCountSublatticeCS,1) * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)
        end if
        if (j == k) then
            l = j
        else if (j > k) then
            cycle LOOP_readPairs
        else
            l = nSublatticeElementsCS(nCountSublatticeCS,1) + j + ((k-2)*(k-1)/2)
        end if
        dCoordinationNumberCS(nCountSublatticeCS, l + p, 1) = dTempVec(1)
        dCoordinationNumberCS(nCountSublatticeCS, l + p, 2) = dTempVec(2)
        dCoordinationNumberCS(nCountSublatticeCS, l + p, 3) = dTempVec(3)
        dCoordinationNumberCS(nCountSublatticeCS, l + p, 4) = dTempVec(4)
        lPairSet(l + p) = .TRUE.
    end do LOOP_readPairs

    ! Increase pairs counter to include default pairs
    nPairsSROCS(nCountSublatticeCS,2) = nSpeciesPhaseCS(i) - nSpeciesPhaseCS(i-1)

    LOOP_allSROPairs: do k = 1, nPairsSROCS(nCountSublatticeCS,2)

        ! If coordinations already set, skip rest
        if (lPairSet(k)) cycle LOOP_allSROPairs

        ! Constituent indices:
        a = iPairIDCS(nCountSublatticeCS,k,1)
        b = iPairIDCS(nCountSublatticeCS,k,2)
        x = iPairIDCS(nCountSublatticeCS,k,3) - nSublatticeElementsCS(nCountSublatticeCS,1)
        y = iPairIDCS(nCountSublatticeCS,k,4) - nSublatticeElementsCS(nCountSublatticeCS,1)

        ! Constituent charges
        qa = dSublatticeChargeCS(nCountSublatticeCS,1,a)
        qb = dSublatticeChargeCS(nCountSublatticeCS,1,b)
        qx = dSublatticeChargeCS(nCountSublatticeCS,1,x)
        qy = dSublatticeChargeCS(nCountSublatticeCS,1,y)

        if (a == b) then
            ia2x2 = a + (x - 1) * (nSublatticeElementsCS(nCountSublatticeCS,1) &
                                    * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)
            ia2y2 = a + (y - 1) * (nSublatticeElementsCS(nCountSublatticeCS,1) &
                                    * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)

            dZAa2x2 = dCoordinationNumberCS(nCountSublatticeCS, ia2x2, 1)
            dZXa2x2 = dCoordinationNumberCS(nCountSublatticeCS, ia2x2, 3)
            dZAa2y2 = dCoordinationNumberCS(nCountSublatticeCS, ia2y2, 1)
            dZYa2y2 = dCoordinationNumberCS(nCountSublatticeCS, ia2y2, 4)

            dCoordinationNumberCS(nCountSublatticeCS, k, 1) = ((dZAa2x2 * qx) + (dZAa2y2 * qy)) / (qx + qy)
            dCoordinationNumberCS(nCountSublatticeCS, k, 2) = ((dZAa2x2 * qx) + (dZAa2y2 * qy)) / (qx + qy)
            dCoordinationNumberCS(nCountSublatticeCS, k, 3) = ((dZXa2x2 * qx) + (dZYa2y2 * qx)) / (qx + qy)
            dCoordinationNumberCS(nCountSublatticeCS, k, 4) = ((dZXa2x2 * qy) + (dZYa2y2 * qy)) / (qx + qy)

            cycle LOOP_allSROPairs
        end if

        if (x == y) then
            ia2x2 = a + (x - 1) * (nSublatticeElementsCS(nCountSublatticeCS,1) &
                                    * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)
            ib2x2 = b + (x - 1) * (nSublatticeElementsCS(nCountSublatticeCS,1) &
                                    * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)

            dZAa2x2 = dCoordinationNumberCS(nCountSublatticeCS, ia2x2, 1)
            dZXa2x2 = dCoordinationNumberCS(nCountSublatticeCS, ia2x2, 3)
            dZBb2x2 = dCoordinationNumberCS(nCountSublatticeCS, ib2x2, 1)
            dZXb2x2 = dCoordinationNumberCS(nCountSublatticeCS, ib2x2, 4)

            dCoordinationNumberCS(nCountSublatticeCS, k, 1) = ((dZAa2x2 * qa) + (dZBb2x2 * qa)) / (qa + qb)
            dCoordinationNumberCS(nCountSublatticeCS, k, 2) = ((dZAa2x2 * qb) + (dZBb2x2 * qb)) / (qa + qb)
            dCoordinationNumberCS(nCountSublatticeCS, k, 3) = ((dZXa2x2 * qa) + (dZXb2x2 * qb)) / (qa + qb)
            dCoordinationNumberCS(nCountSublatticeCS, k, 4) = ((dZXa2x2 * qa) + (dZXb2x2 * qb)) / (qa + qb)

            cycle LOOP_allSROPairs
        end if

        ! Find triplet indices
        l = nSublatticeElementsCS(nCountSublatticeCS,1) + a + ((b-2)*(b-1)/2)
        p = (nSublatticeElementsCS(nCountSublatticeCS,2) + (x - 1) + ((y-2)*(y-1)/2)) &
              * (nSublatticeElementsCS(nCountSublatticeCS,1) * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)

        ia2xy = a + p
        ib2xy = b + p
        iabx2 = l + (x - 1) * (nSublatticeElementsCS(nCountSublatticeCS,1) &
                                * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)
        iaby2 = l + (y - 1) * (nSublatticeElementsCS(nCountSublatticeCS,1) &
                                * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)

        dZAa2xy = dCoordinationNumberCS(nCountSublatticeCS, ia2xy, 1)
        dZXa2xy = dCoordinationNumberCS(nCountSublatticeCS, ia2xy, 3)
        dZYa2xy = dCoordinationNumberCS(nCountSublatticeCS, ia2xy, 4)

        dZBb2xy = dCoordinationNumberCS(nCountSublatticeCS, ib2xy, 2)
        dZXb2xy = dCoordinationNumberCS(nCountSublatticeCS, ib2xy, 3)
        dZYb2xy = dCoordinationNumberCS(nCountSublatticeCS, ib2xy, 4)

        dZAabx2 = dCoordinationNumberCS(nCountSublatticeCS, iabx2, 1)
        dZBabx2 = dCoordinationNumberCS(nCountSublatticeCS, iabx2, 2)
        dZXabx2 = dCoordinationNumberCS(nCountSublatticeCS, iabx2, 3)

        dZAaby2 = dCoordinationNumberCS(nCountSublatticeCS, iaby2, 1)
        dZBaby2 = dCoordinationNumberCS(nCountSublatticeCS, iaby2, 2)
        dZYaby2 = dCoordinationNumberCS(nCountSublatticeCS, iaby2, 4)

        dF = (1D0/8D0) * ((qa / dZAa2xy) + (qb / dZBb2xy) + (qx / dZXabx2) + (qy / dZYaby2))

        ! Set coordinations
        dCoordinationNumberCS(nCountSublatticeCS, k, 1) = 1 / (((dZXabx2 / (qx * dZAabx2)) + (dZYaby2 / (qy * dZAaby2))) * dF)
        dCoordinationNumberCS(nCountSublatticeCS, k, 2) = 1 / (((dZXabx2 / (qx * dZBabx2)) + (dZYaby2 / (qy * dZBaby2))) * dF)
        dCoordinationNumberCS(nCountSublatticeCS, k, 3) = 1 / (((dZAa2xy / (qa * dZXa2xy)) + (dZBb2xy / (qb * dZXb2xy))) * dF)
        dCoordinationNumberCS(nCountSublatticeCS, k, 4) = 1 / (((dZAa2xy / (qa * dZYa2xy)) + (dZBb2xy / (qb * dZYb2xy))) * dF)

    end do LOOP_allSROPairs

    dStoichSpeciesOld = dStoichSpeciesCS(1:nSpeciesCS,1:nElementsCS)
    dStoichPairsCS(nCountSublatticeCS,1:nPairsSROCS(nCountSublatticeCS,2),1:nElementsCS) &
                  = dStoichSpeciesCS((nSpeciesPhaseCS(i-1) + 1):nSpeciesPhaseCS(i),1:nElementsCS)
    dStoichSpeciesCS((nSpeciesPhaseCS(i-1) + 1):nSpeciesPhaseCS(i),1:nElementsCS) = 0D0

    ! Loop through all pairs to calculate stoichiometry entries for quadruplets:
    do j = 1, nPairsSROCS(nCountSublatticeCS,2)
        a = iSublatticeElementsCS(nCountSublatticeCS,1,iPairIDCS(nCountSublatticeCS, j, 1))
        b = iSublatticeElementsCS(nCountSublatticeCS,1,iPairIDCS(nCountSublatticeCS, j, 2))
        x = iSublatticeElementsCS(nCountSublatticeCS,2,iPairIDCS(nCountSublatticeCS, j, 3) &
          - nSublatticeElementsCS(nCountSublatticeCS,1))
        y = iSublatticeElementsCS(nCountSublatticeCS,2,iPairIDCS(nCountSublatticeCS, j, 4) &
          - nSublatticeElementsCS(nCountSublatticeCS,1))

        ! Matching 4 pairs version
        l = j + nSpeciesPhaseCS(i-1)

        ! Just get the quads directly version
        dStoichSpeciesCS(l,a) = dStoichSpeciesCS(l,a) + (1 / dCoordinationNumberCS(nCountSublatticeCS, j, 1))
        dStoichSpeciesCS(l,b) = dStoichSpeciesCS(l,b) + (1 / dCoordinationNumberCS(nCountSublatticeCS, j, 2))
        dStoichSpeciesCS(l,x) = dStoichSpeciesCS(l,x) + (1 / dCoordinationNumberCS(nCountSublatticeCS, j, 3))
        dStoichSpeciesCS(l,y) = dStoichSpeciesCS(l,y) + (1 / dCoordinationNumberCS(nCountSublatticeCS, j, 4))

        ! Create quadruplet names
        cSpeciesNameCS(j + nSpeciesPhaseCS(i-1)) = TRIM(cConstituentNames1(iPairIDCS(nCountSublatticeCS, j, 1))) // '-' &
                                                // TRIM(cConstituentNames1(iPairIDCS(nCountSublatticeCS, j, 2))) // '-' &
                                                // TRIM(cConstituentNames2(iPairIDCS(nCountSublatticeCS, j, 3) &
                                                - nSublatticeElementsCS(nCountSublatticeCS,1)))                  // '-' &
                                                // TRIM(cConstituentNames2(iPairIDCS(nCountSublatticeCS, j, 4) &
                                                - nSublatticeElementsCS(nCountSublatticeCS,1)))

    end do

    ! Loop through excess mixing parameters:
    LOOP_ExcessMixingSUBG: do

        ! Read in number of constituents involved in parameter:
        read (1,*,IOSTAT = INFO) iRegularParamCS(nParamCS+1,1)

        ! The end of the parameter listing is marked by "0":
        if (iRegularParamCS(nParamCS+1,1) == 0) exit LOOP_ExcessMixingSUBG

        ! Check if the parameter is binary or ternary:
        if ((iRegularParamCS(nParamCS+1,1) == 3) .OR. (iRegularParamCS(nParamCS+1,1) == 4)) then

            ! Count the number of parameters:
            nParamCS = nParamCS + 1

            ! Mixing terms:
            read (1,*,IOSTAT = INFO) cDummyVec(1), iRegularParamCS(nParamCS,2:9)

            ! According to Patrice Chartrand, he has no idea what these two lines mean. Ignore.
            read (1,*,IOSTAT = INFO) dTempVec(1:6)
            read (1,*,IOSTAT = INFO) dTempVec(1:6)

            ! Read in the first line of the excess gibbs energy of mixing terms.
            read (1,*,IOSTAT = INFO) dTempVec(1:6)
            dRegularParamCS(nParamCS,1:6) = dTempVec(1:6)

            ! I HAVE NO IDEA IF THIS LINE IS NEEDED OR NOT. DON'T DO ANYTHING WITH IT FOR NOW.
            read (1,*,IOSTAT = INFO) dTempVec(1:2)

        else
            !! This parameter is not recognized; record an error.
            INFO = 1600 + i
            return
        end if

    end do LOOP_ExcessMixingSUBG

    ! Report an error if necessary:
    if (INFO /= 0) INFO = 1600 + i

    deallocate(lPairSet)

    return

end subroutine ParseCSDataBlockSUBG
