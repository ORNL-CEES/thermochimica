
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
        
    integer :: i, j, k, m, p, q, r
    integer :: iSolnIndex
    integer :: iFirst, iLast
    real(8) :: dTemp, dSum
    real(8) :: dZAAA, dZABA, dZBAB
    real(8) :: x, y, z
    real(8), allocatable, dimension(:) :: dX, dY


    ! Only proceed if the correct phase type is selected:
    IF_SUBG: if (cSolnPhaseType(iSolnIndex) == 'SUBG') then

        ! Define temporary variables for sake of convenience:
        iFirst = nSpeciesPhase(iSolnIndex-1) + 1
        iLast  = nSpeciesPhase(iSolnIndex)

        ! Allocate allocatable arrays:
        if (allocated(dX)) deallocate(dX)
        if (allocated(dY)) deallocate(dY)
        j = iLast - iFirst + 1
        allocate(dX(j),dY(j))

        ! Initialize variables:
        dX                                = 0D0
        dY                                = 0D0
        dSum                              = 0D0
        dChemicalPotential(iFirst:iLast)  = 0D0
        dPartialExcessGibbs(iFirst:iLast) = 0D0

        ! Compute moles of compound end members 
        ! (See eq. [14] from Pelton, Chartrand, Met. Mat. Trans., 32 (2001) 1355):
        LOOP_A: do i = iFirst, iFirst + nPairsSRO(1,1) - 1
            j     = i - iFirst + 1
            dZAAA = dCoordinationNumber(j,1)

            ! Loop through i-j pairs to compute
            dTemp = 2D0 * dMolesSpecies(i) / dZAAA
            do k = iFirst + nPairsSRO(1,1), iLast
                ! Index of i-j pair (amoung pairs):
                m = k - iFirst + 1

                ! Verify that AB is indeed paired with AA or BB:
                if (iPairID(j,1) == iPairID(m,1))  then
                    dZABA = dCoordinationNumber(m,1)
                    dTemp = dTemp + dMolesSpecies(k) / dZABA
                elseif (iPairID(j,1) == iPairID(m,2))  then
                    dZBAB = dCoordinationNumber(m,2)
                    dTemp = dTemp + dMolesSpecies(k) / dZBAB
                end if
            end do

            ! Temporarily use the dPartialExcessGibbs for storing the number of moles
            ! this species (to reduce memory requirements):
            dPartialExcessGibbs(i) = dTemp
            dSum = dSum + dTemp
        end do LOOP_A

        ! Now, compute mole fractions of compound end members and the coordination equivalent fractions:
        LOOP_B: do i = iFirst, iFirst + nPairsSRO(1,1) - 1
            j = i - iFirst + 1
            dX(j) = dPartialExcessGibbs(i) / dSum

            dTemp = dMolFraction(i)
            ! Verify that AB is indeed paired with AA or BB:
            do k = iFirst + nPairsSRO(1,1), iLast
                ! Index of i-j pair (amoung pairs):
                m = k - iFirst + 1
                if ((iPairID(j,1) == iPairID(m,1)).OR.(iPairID(j,1) == iPairID(m,2)))  then
                    dTemp = dTemp + dMolFraction(k) / 2D0
                end if
            end do
            dY(j) = dTemp
        end do LOOP_B

        ! -----------------------------------------------------
        ! COMPUTE REFERENCE GIBBS ENERGY AND IDEAL MIXING TERMS
        ! -----------------------------------------------------

        ! Loop through A-A pairs:
        LOOP_C: do i = iFirst, iFirst + nPairsSRO(1,1) - 1
            ! Store indices:
            j = i - iFirst + 1                      ! Index of pairs for coordination number
            k = iFirst + nPairsSRO(1,1) + j         ! Index of AB pair

            ! Store coordination numbers:
            dZAAA = dCoordinationNumber(j,1)

            ! Compute standard reference Gibbs energy:
            dChemicalPotential(i) = dStdGibbsEnergy(i) * 2D0 / dZAAA

            ! Compute ideal mixing component:
            dChemicalPotential(i) = dChemicalPotential(i) + (2D0 / dZAAA) * DLOG(dX(j)) + DLOG(dMolFraction(i) / dY(j)**2)

        end do LOOP_C

        ! Loop through A-B pairs:
        LOOP_D: do i = iFirst + nPairsSRO(1,1), iLast
            m = i - iFirst + 1          ! Index of AB
            j = iPairID(m,1)            ! Index of AA
            k = iPairID(m,2)            ! Index of BB

            ! Store coordination numbers:
            dZABA = dCoordinationNumber(m,1)
            dZBAB = dCoordinationNumber(m,2)

            ! Compute standard reference Gibbs energy:
            ! NOTE: I did not include $\Delta g_{AB}^{\circ}$ in this equation because I think it makes more sense
            ! to include it in the excess mixing section.
            dChemicalPotential(i) = dStdGibbsEnergy(j + iFirst - 1) / dZABA + dStdGibbsEnergy(k + iFirst - 1) / dZBAB

            ! Compute ideal mixing component:
            dChemicalPotential(i) = dChemicalPotential(i) + DLOG(dMolFraction(i) / (2D0 * dY(j) * dY(k))) &
                + DLOG(dX(j)) / dZABA + DLOG(dX(k)) / dZBAB

        end do LOOP_D

        ! ---------------------------------------------------------
        ! COMPUTE PARTIAL MOLAR EXCESS GIBBS ENERGY OF MIXING TERMS
        ! ---------------------------------------------------------

        ! Re-initialize variables:
        dPartialExcessGibbs(iFirst:iLast) = 0D0

        ! Loop through excess mixing parameters: 
        LOOP_Param: do m = nParamPhase(iSolnIndex-1) + 1, nParamPhase(iSolnIndex)
            i = iRegularParam(m,2)              ! Index of AA
            j = iRegularParam(m,3)              ! Index of BB
            k = iRegularParam(m,4)              ! Index of AB
            p = iRegularParam(m,6)              ! Exponent of AA
            q = iRegularParam(m,7)              ! Exponent of BB
            r = iRegularParam(m,8)              ! Exponent of AB
            x = dMolFraction(iFirst + i - 1)    ! x_AA
            y = dMolFraction(iFirst + j - 1)    ! x_BB
            z = dMolFraction(iFirst + k - 1)    ! x_AB

            ! Contribution to AA:
            dTemp = z * x**(p-1) * y**(q) * (DBLE(p)*(y+z) - DBLE(q)*x)
            dTemp = dTemp * dExcessGibbsParam(m) / 2D0
            dPartialExcessGibbs(iFirst + i - 1) = dPartialExcessGibbs(iFirst + i - 1) + dTemp

            ! Contribution to BB:
            dTemp = z * x**(p+1) * y**(q) * (DBLE(q)*(x+z) - DBLE(p)*y)
            dTemp = dTemp * dExcessGibbsParam(m) / (2D0 * x * y)
            dPartialExcessGibbs(iFirst + j - 1) = dPartialExcessGibbs(iFirst + j - 1) + dTemp

            ! Contribution to AB:
            dTemp = x**(p+1) * y**(q) * (z*(1D0 - DBLE(p) - DBLE(q)) + x + y)
            dTemp = dTemp * dExcessGibbsParam(m) / (2D0 * x)
            dPartialExcessGibbs(iFirst + k - 1) = dPartialExcessGibbs(iFirst + k - 1) + dTemp

        end do LOOP_Param

    end if IF_SUBG

    ! Deallocate allocatable arrays:
    deallocate(dX,dY)

    return

end subroutine CompExcessGibbsEnergySUBG
