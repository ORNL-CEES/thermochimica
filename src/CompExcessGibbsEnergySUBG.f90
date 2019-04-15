
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

    integer :: i, j, k, l, m, p, q, r
    integer :: iSolnIndex
    integer :: iFirst, iLast
    real(8) :: dTemp, dSum
    real(8) :: dZAAA, dZABA, dZBAB
    real(8) :: x, y, z
    real(8), allocatable, dimension(:) :: dX, dY, dN
    ! dX is X_i in Pelton et al., while X_ij in that paper corresponds to dMolesSpecies


    ! Only proceed if the correct phase type is selected:
    IF_SUBG: if (cSolnPhaseType(iSolnIndex) == 'SUBG' .OR. cSolnPhaseType(iSolnIndex) == 'SUBQ') then

        ! Define temporary variables for sake of convenience:
        iFirst = nSpeciesPhase(iSolnIndex-1) + 1
        iLast  = nSpeciesPhase(iSolnIndex)

        ! Allocate allocatable arrays:
        if (allocated(dX)) deallocate(dX)
        if (allocated(dY)) deallocate(dY)
        if (allocated(dN)) deallocate(dN)
        j = iLast - iFirst + 1
        allocate(dX(j),dY(j),dN(j))

        ! Initialize variables:
        dX                                = 0D0
        dY                                = 0D0
        dN                                = 0D0
        dSum                              = 0D0
        dChemicalPotential(iFirst:iLast)  = 0D0
        dPartialExcessGibbs(iFirst:iLast) = 0D0

        ! Compute moles of compound end members
        ! (See eq. [14] from Pelton, Chartrand, Met. Mat. Trans., 32 (2001) 1355):
        LOOP_A: do i = 1, nPairsSRO(iSolnIndex,2)
            j = iFirst + i - 1
            ! Skip AB pairs:
            if (iPairID(j,1) /= iPairID(j,2)) CYCLE LOOP_A
            dZAAA = dCoordinationNumber(j,1)

            ! Loop through i-j pairs to compute
            dTemp = 2D0 * dMolesSpecies(j) / dZAAA
            LOOP_A_interior: do k = 1, nPairsSRO(iSolnIndex,2)
                l = iFirst + k - 1
                ! Skip AA pairs:
                if (iPairID(l,1) == iPairID(l,2)) CYCLE LOOP_A_interior
                ! Verify that AB is indeed paired with AA or BB:
                if (iPairID(j,1) == iPairID(l,1))  then
                    dZABA = dCoordinationNumber(l,1)
                    dTemp = dTemp + dMolesSpecies(l) / dZABA
                elseif (iPairID(j,1) == iPairID(l,2))  then
                    dZBAB = dCoordinationNumber(l,2)
                    dTemp = dTemp + dMolesSpecies(l) / dZBAB
                end if
            end do LOOP_A_interior

            dN(i) = dTemp
            dSum = dSum + dTemp
        end do LOOP_A

        ! Now, compute mole fractions of compound end members and the coordination equivalent fractions:
        LOOP_B: do i = 1, nPairsSRO(iSolnIndex,2)
            j = iFirst + i - 1
            ! Skip AB pairs:
            if (iPairID(j,1) /= iPairID(j,2)) CYCLE LOOP_B
            ! Eq [4]:
            dX(i) = dN(i) / dSum
            ! Eqs [6]:
            dTemp = dMolFraction(j)
            ! Verify that AB is indeed paired with AA or BB:
            LOOP_B_interior: do k = 1, nPairsSRO(iSolnIndex,2)
                l = iFirst + k - 1
                ! Skip AA pairs:
                if (iPairID(l,1) == iPairID(l,2)) CYCLE LOOP_B_interior
                if ((iPairID(j,1) == iPairID(l,1)).OR.(iPairID(j,1) == iPairID(l,2)))  then
                    dTemp = dTemp + dMolFraction(l) / 2D0
                end if
            end do LOOP_B_interior
            dY(i) = dTemp
        end do LOOP_B

        ! ---------------------------------------------------------------
        ! COMPUTE REFERENCE GIBBS ENERGY AND IDEAL MIXING TERMS
        ! ---------------------------------------------------------------

        ! Loop through all pairs:
        LOOP_C: do i = 1, nPairsSRO(iSolnIndex,2)
            j = iFirst + i - 1
            if (iPairID(j,1) == iPairID(j,2)) then
                ! These are AA pairs
                ! Store coordination numbers:
                dZAAA = dCoordinationNumber(j,1)

                ! Compute standard reference Gibbs energy (Eq [15]):
                dChemicalPotential(j) = dStdGibbsEnergy(j) * 2D0 / dZAAA

                ! Compute ideal mixing component (Eq [34]?):
                dChemicalPotential(j) = dChemicalPotential(j) + (2D0 / dZAAA) * DLOG(dX(i)) + DLOG(dMolFraction(j) / dY(i)**2)
            else
                ! These are AB pairs
                k = iPairID(j,1)            ! Index of AA
                l = iPairID(j,2)            ! Index of BB

                ! Store coordination numbers:
                dZABA = dCoordinationNumber(j,1)
                dZBAB = dCoordinationNumber(j,2)

                ! Compute standard reference Gibbs energy:
                ! NOTE: I did not include $\Delta g_{AB}^{\circ}$ in this equation because I think it makes more sense
                ! to include it in the excess mixing section.
                ! Eq [16]:
                dChemicalPotential(j) = dStdGibbsEnergy(k + iFirst - 1) / dZABA + dStdGibbsEnergy(l + iFirst - 1) / dZBAB

                ! Compute ideal mixing component:
                dChemicalPotential(j) = dChemicalPotential(j) + DLOG(dMolFraction(j) / (2D0 * dY(k) * dY(l))) &
                                                              + DLOG(dX(k)) / dZABA + DLOG(dX(l)) / dZBAB
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
            LOOP_FindPair: do l = 1, nPairsSRO(iSolnIndex,2)
                if (((iPairID(l + iFirst - 1,1) == i) .AND. (iPairID(l + iFirst - 1,2) == j)) .OR. &
                    ((iPairID(l + iFirst - 1,1) == j) .AND. (iPairID(l + iFirst - 1,2) == i)))  then
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
    deallocate(dX,dY,dN)

    return

end subroutine CompExcessGibbsEnergySUBG
