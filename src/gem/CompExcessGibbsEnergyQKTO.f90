
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
    real(8)                       :: dGParam, xT
    real(8), dimension(nMaxParam) :: dPartialGParam


    ! Return control to the parent subroutine if there aren't any interaction parameters for this phase:
    if ((nParamPhase(iSolnIndex) - nParamPhase(iSolnIndex-1) == 0).OR. &
        (cSolnPhaseType(iSolnIndex) /= 'QKTO')) return

    ! Loop through all interaction parameters in this phase:
    do iParam = nParamPhase(iSolnIndex-1)+1, nParamPhase(iSolnIndex)

        ! Compute the partial molar excess Gibbs energy of each sub-system in a regular solution phase:
        call PolyRegularQKTO(iSolnIndex,iParam,xT,dGParam,dPartialGParam)

        ! Perform a Kohler interpolation amoungst the sub-systems to compute the overall partial molar
        ! excess Gibbs energy of mixing:
        call KohlerInterpolate(iSolnIndex,iParam,xT,dGParam,dPartialGParam)

    end do

    return

end subroutine CompExcessGibbsEnergyQKTO
