
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
        ! Now also performs interpolation
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
    USE ModuleThermoIO

    implicit None

    integer                       :: i, j, zT, iParam, iSolnIndex, iGroup1, iGroup2, iGroupTemp, nSpeciesParam, p, q, r, a, b, c
    real(8)                       :: xT, dGex, dXi1, dXi2, dXiDen, dDgex
    real(8), dimension(nSpeciesPhase(iSolnIndex) - nSpeciesPhase(iSolnIndex-1)) :: y, dPartialGParam
    logical, dimension(nSpeciesPhase(iSolnIndex) - nSpeciesPhase(iSolnIndex-1)) :: lAsymmetric1, lAsymmetric2

    ! Initialize variables:
    xT             = 0D0
    y              = 0D0
    zT             = 0
    dPartialGParam = 0D0
    dXi1           = 0D0
    dXi2           = 0D0
    lAsymmetric1   = .FALSE.
    lAsymmetric2   = .FALSE.
    nSpeciesParam  = iRegularParam(iParam,1)
    a = iRegularParam(iParam,2)
    b = iRegularParam(iParam,3)
    p = iRegularParam(iParam,nSpeciesParam+2)
    q = iRegularParam(iParam,nSpeciesParam+3)
    ! Set ternary parameters initially to 0, so that they may be used for any case
    c = 0
    r = 0
    if (nSpeciesParam == 3) then
        c = iRegularParam(iParam,4)
        r = iRegularParam(iParam,nSpeciesParam+4)
    end if

    ! Get chemical groups of first two species
    iGroup1 = INT(dQKTOParams(nSpeciesPhase(iSolnIndex-1) + a,2))
    iGroup2 = INT(dQKTOParams(nSpeciesPhase(iSolnIndex-1) + b,2))

    ! Calculate symmetry of all species with these first two
    ! Count them as asymmetric with respect to themselves
    lAsymmetric1(a) = .TRUE.
    lAsymmetric2(b) = .TRUE.
    ! If they are equal, everything else will be symmetric
    if (iGroup1 /= iGroup2) then
        do j = nSpeciesPhase(iSolnIndex-1) + 1, nSpeciesPhase(iSolnIndex)
            i = j - nSpeciesPhase(iSolnIndex-1)
            iGroupTemp = INT(dQKTOParams(j,2))
            if      (iGroupTemp == iGroup1) then
                lAsymmetric1(i) = .TRUE.
            else if (iGroupTemp == iGroup2) then
                lAsymmetric2(i) = .TRUE.
            end if
        end do
    end if

    ! Compute sum of equivalent fractions
    do j = nSpeciesPhase(iSolnIndex-1) + 1, nSpeciesPhase(iSolnIndex)
        i = j - nSpeciesPhase(iSolnIndex-1)
        xT = xT + dMolFraction(j) * dQKTOParams(j,1)
    end do

    ! Compute the equivalent mole fractions
    do j = nSpeciesPhase(iSolnIndex-1) + 1, nSpeciesPhase(iSolnIndex)
        i = j - nSpeciesPhase(iSolnIndex-1)
        y(i) = dMolFraction(j) * dQKTOParams(j,1) / xT
    end do

    ! Now calculate xis and sigma
    do j = nSpeciesPhase(iSolnIndex-1) + 1, nSpeciesPhase(iSolnIndex)
        i = j - nSpeciesPhase(iSolnIndex-1)
        if      (lAsymmetric1(i)) then
            dXi1 = dXi1 + y(i)
        else if (lAsymmetric2(i)) then
            dXi2 = dXi2 + y(i)
        end if
    end do

    dXiDen = dXi1 + dXi2

    ! ChemSage guide says only symmetric ternaries are allowed, but I'm not sure if this is still true
    ! Anyway, we haven't seen any others yet
    ! Also the textbook equations for ternaries don't reproduce the symmetry in FactSage, so for now
    ! I'm just implementing symmetric ternary explicitly
    ! I will still try to make this easy to change later
    if (nSpeciesParam == 3) then
        dXi1 = y(a)
        dXi2 = y(b)
        dXiDen = y(a) + y(b) + y(c)
    end if

    ! Calculate g^excess for binary part
    dGex = dExcessGibbsParam(iParam) * (dXi1**(p - 1)) * (dXi2**(q - 1)) / (dXiDen ** (p + q + r - nSpeciesParam))
    dGex = dGex * y(a) * y(b) * xT

    ! Include ternary factor if present
    if (nSpeciesParam == 3) then
        ! Trying to write this parallel to binary
        dGex = dGex * y(c)**(r - 1)
        dGex = dGex * y(c)
    end if

    do j = nSpeciesPhase(iSolnIndex-1) + 1, nSpeciesPhase(iSolnIndex)
        i = j - nSpeciesPhase(iSolnIndex-1)
        dDgex = -1D0 / xT
        if      (lAsymmetric1(i)) then
            dDgex = dDgex + (p - 1) / dXi1 - (p + q + r - nSpeciesParam) / dXiDen
        else if (lAsymmetric2(i)) then
            dDgex = dDgex + (q - 1) / dXi2 - (p + q + r - nSpeciesParam) / dXiDen
        end if
        if      (i == a) then
            dDgex = dDgex + 1D0 / y(a) / xT
        else if (i == b) then
            dDgex = dDgex + 1D0 / y(b) / xT
        else if (i == c) then
            dDgex = dDgex + 1D0 / y(c) / xT
        end if
        ! Ternary part of derivative
        if (nSpeciesParam == 3) then
            dDgex = dDgex - 1D0 / xT
            if (i == c) then
                dDgex = dDgex + (r - 1) / y(c) - (p + q + r - nSpeciesParam) / dXiDen
            end if
        end if
        dPartialExcessGibbs(j) = dPartialExcessGibbs(j) + dDgex * dGex * dQKTOParams(j,1)
    end do

    return

end subroutine PolyRegularQKTO
