
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompExcessGibbsEnergyQKTO.f90
    !> \brief   Compute the partial molar excess Gibbs energy of mixing of solution phase constituents in a QKTO
    !!           solution phase.
    !> \author  M.H.A. Piro
    !> \sa      CompExcessGibbsEnergy.f90
    !> \sa      CompExcessGibbsEnergyRKMP.f90
    !> \sa      CompExcessGibbsEnergySUBL.f90
    !> \date    June 13, 2012
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   06/13/2012      M.H.A. Piro         Original code.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the partial molar excess Gibbs energy of mixing
    !! (dPartialExcessGibbs) of all constituents in a non-ideal solution phase designated as 'QKTO'
    !! (Quasi-chemical Kohlter-TOop).  The PolyRegular subroutine computes the excess Gibbs energy of mixing of
    !! a regular solution sub-system (see PolyRegular for a definition) and the KohlerInterpolate subroutine
    !! performs a Kohler interpolation of a sub-system to a phase.
    !!
    !! The molar excess Gibbs energy of mixing of a binary sub-system for a QKTO model is:
    !!
    !! \f$ g_{\lambda,z}^{ex} = L_z x_1^a x_2^b \f$
    !!
    !! where \f$ L_z \f$ is the mixing parameter, \f$ x_1 \f$ and \f$ x_2 \f$ are the mole fractions for
    !! constituents 1 and 2 in the binary term and \f$ a \f$ and \f$ b \f$ are the exponents for constituents
    !! 1 and 2, respectively.
    !!
    !! The molar excess Gibbs energy of mixing for solution phase \f$ \lambda \f$ using Kohler's interpolation
    !! scheme gives
    !!
    !! \f$ g_{\lambda}^{ex} = (x_1 + x_2)^2 g_{\lambda,z}^{ex} \f$
    !!
    !! Similarly, the molar excess Gibbs energy of mixing of a ternary sub-system for a QKTO model is:
    !!
    !! \f$ g_{\lambda,z}^{ex} = L_z x_1^a x_2^b x_3^c \f$
    !!
    !! which is related to the molar excess Gibbs energy of mixing of the phase via Kohler's interpolation:
    !!
    !! \f$ g_{\lambda}^{ex} = (x_1 + x_2 + x_3)^2 g_{\lambda,z}^{ex} \f$
    !!
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in] iSolnIndex    Absolute index of a solution phase
    !
    ! nSpeciesPhase             Highest index number of a species in a particular solution phase
    ! nParam                    Number of parameters
    ! iParam                    Index number of a parameter.
    ! dChemicalPotential        The estimated chemical potential vector.  To be precise, this is defined as the
    !                           molar Gibbs energy of the pure species minus the proper chemical potential
    !                           defined by the element potentials.
    ! dPartialExcessGibbs       Partial molar excess Gibbs energy of mixing of species.
    ! dPartialExcessGibbsLast   Partial molar excess Gibbs energy of mixing of species from the last iteration.
    ! dMolFraction              Current estimated mole fraction.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CompExcessGibbsEnergyQKTO(iSolnIndex)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                       :: iParam, iSolnIndex

    ! Return control to the parent subroutine if there aren't any interaction parameters for this phase:
    if ((nParamPhase(iSolnIndex) - nParamPhase(iSolnIndex-1) == 0).OR. &
        (cSolnPhaseType(iSolnIndex) /= 'QKTO')) return

    ! Loop through all interaction parameters in this phase:
    do iParam = nParamPhase(iSolnIndex-1)+1, nParamPhase(iSolnIndex)
        ! Compute the partial molar excess Gibbs energy of each sub-system in a regular solution phase:
        call PolyRegularQKTO(iSolnIndex,iParam)
    end do

    return

end subroutine CompExcessGibbsEnergyQKTO

    !-------------------------------------------------------------------------------------------------------------
    !
    ! Purpose
    ! =======
    !
    !> \details The purpose of this subroutine is to compute the partial molar excess Gibbs energy of
    !! mixing of species in a regular sub-system.  Note that the effective quantity (i.e., "y") is represented
    !! relative to the particular parameter of interest and differs from that of the phase as a whole.
    !!
    !! For more information on the derivation of the thermodynamic equations used in this subroutine,
    !! refer to the following paper:
    !!
    !!        A.D. Pelton and C.W. Bale, "Computational Techniques for the Treatment
    !!        of Thermodynamic Data in Multicomponent Systems and the Calculation of
    !!        Phase Equilibria," CALPHAD, V. 1, N. 3 (1977) 253-273.
    !!
    !
    ! Pertinent Variables
    ! ===================
    !
    !> \param[in]  iSolnIndex       Integer scalar of the solution index.
    !> \param[in]  iParam           Integer scalar of the parameter index.
    !> \param[out] xT               Sum of mole fractions of actual constituents in solution phase
    !> \param[out] dGParam          Excess Gibbs energy of sub-system
    !> \param[out] dPartialGParam   Partial excess Gibbs energy of a constituent in the sub-system
    !
    ! y                             A double real vector representing the equivalent mole fractions of each
    !                                constituent in the sub-system (parameter).
    ! zT                            A double real scalar representing the sum of exponents in the sub-system.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine PolyRegularQKTO(iSolnIndex,iParam)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit None

    integer                      :: i, j, k, m, zT, iParam, iSolnIndex
    real(8)                      :: xT, dGParam
    real(8),dimension(nMaxParam) :: y, dPartialGParam

    ! Initialize variables:
    xT             = 0D0
    y              = 0D0
    zT             = 0
    dGParam        = dExcessGibbsParam(iParam)
    dPartialGParam = 0D0

    ! Compute the sum of mole fractions of real components and the sum of their exponents in the sub-system:
    do i = 1, iRegularParam(iParam,1)
        j  = nSpeciesPhase(iSolnIndex-1) + iRegularParam(iParam,i+1)
        k  = iRegularParam(iParam,1) + 1 + i
        xT = xT + dMolFraction(j)
        zT = zT + iRegularParam(iParam,k)
    end do

    ! Compute the equivalent mole fractions of components and the integral excess Gibbs energy of the sub-system:
    do i = 1, iRegularParam(iParam,1)
        j       = nSpeciesPhase(iSolnIndex-1) + iRegularParam(iParam,i+1)
        k       = iRegularParam(iParam,1) + 1 + i
        y(i)    = dMolFraction(j) / xT
        dGParam = dGParam * (y(i) ** iRegularParam(iParam,k))
    end do

    ! Compute the partial excess Gibbs energy of mixing per equivalent mole in the sub-system:
    do i = 1, iRegularParam(iParam,1)
        k = iRegularParam(iParam,iRegularParam(iParam,1) + 1 + i)
        dPartialGParam(i) = dExcessGibbsParam(iParam) * (DFLOAT(k) * y(i)**(k - 1) + (1D0 - DFLOAT(zT)) * y(i)**k)

        ! Loop through parameters:
        do j = 1, iRegularParam(iParam,1)
            ! Cycle if it is the same parameter:
            if (j == i) cycle
            m = iRegularParam(iParam,iRegularParam(iParam,1) + 1 + j)
            dPartialGParam(i) = dPartialGParam(i) * y(j)**m
        end do

    end do

    ! Store the number of components in the sub-system (i.e., binary, ternary or quaternary):
    k = iRegularParam(iParam,1)

    ! Compute the partial molar excess Gibbs energy of mixing using the Kohler interpolation scheme:
    do i = 1, iRegularParam(iParam,1)
        j = nSpeciesPhase(iSolnIndex-1) + iRegularParam(iParam,i+1)
        dPartialExcessGibbs(j) = dPartialExcessGibbs(j) + (xT**(k-1)) * (dPartialGParam(i) + &
            dGParam * (1D0 - xT) * (DFLOAT(k)-1D0))
    end do

    ! Compute the equivalent excess Gibbs energy of mixing of species that are not part of the sub-system:
    LOOP_A: do i = nSpeciesPhase(iSolnIndex-1)+1, nSpeciesPhase(iSolnIndex)
        j = i - nSpeciesPhase(iSolnIndex-1)
        do m = 1,k
            if (j == iRegularParam(iParam,m+1)) cycle LOOP_A
        end do
        dPartialExcessGibbs(i) = dPartialExcessGibbs(i) - (xT**(k)) * (dGParam * (DFLOAT(k)-1D0))
    end do LOOP_A

    return

end subroutine PolyRegularQKTO
