
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    KohlerInterpolate.f90
    !> \brief   Perform a Kohler interpolation for excess mixing terms.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      PolyRegular.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/28/2011      M.H.A. Piro         Original code
    !   07/31/2011      M.H.A. Piro         Clean up code: remove unnecessary variables, update variable names
    !   10/27/2011      M.H.A. Piro         Clean up code: modules, remove unnecessary variables
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to perform a Kohler interpolation of binary/ternary/quaternary
    !! model parameters (provided by PolyRegular.f90) in multi-component phases and return the partial molar
    !! excess Gibbs energy of mixing of a species in a non-ideal solution phase (QKTO).
    !
    !
    ! References:
    ! ===========
    !
    !> \details For more information regarding the Kohler interpolation method and the derivation of the
    !! equations used in this subroutine, refer to the following paper:
    !!
    !!    A.D. Pelton and C.W. Bale, "Computational Techniques for the Treatment
    !!    of Thermodynamic Data in Multicomponent Systems and the Calculation of
    !!    Phase Equilibria," CALPHAD, V. 1, N. 3 (1977) 253-273.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in] iSolnIndex        An integer scalar representing the index of a solution phase.
    !> \param[in] iParam            An integer scalar representing the mixing parameter index.
    !> \param[in] xT                Sum of mole fractions of actual species in solution phase
    !> \param[in] dGParam           Excess Gibbs energy of sub-system
    !> \param[in] dPartialGParam    Partial excess Gibbs energy of species in sub-system
    !
    ! dPartialExcessGibbs           Partial molar excess Gibbs energy of mixing of a species.
    ! nSpeciesPhase                 An integer vector representing the number of species in each solution phase.
    ! iRegularParam                 An integer matrix representing information pertient to regular solution
    !                                models.  The first coefficient represents the number of components in the
    !                                sub-system and the other coefficients represent the indices of components
    !                                in the sub-system.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine KohlerInterpolate(iSolnIndex,iParam,xT,dGParam,dPartialGParam)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                       :: i, j, k, m, iSolnIndex, iParam
    real(8)                       :: xT, dGParam
    real(8), dimension(nMaxParam) :: dPartialGParam


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

end subroutine KohlerInterpolate
