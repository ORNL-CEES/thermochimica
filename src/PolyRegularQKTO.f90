
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    PolyRegularQKTO.f90
    !> \brief   Compute the partial molar excess Gibbs energy of a polynomial regular solution model (QKTO).
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      KohlerInterpolate.f90
    !> \sa      CompExcessGibbsEnergyQKTO.f90
    !> \sa      KohlerInterpolate.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   06/25/2011      M.H.A. Piro         Original code
    !   07/31/2011      M.H.A. Piro         Clean up code: remove unnecessary variables, update variable names
    !   10/27/2011      M.H.A. Piro         Clean up code: modules, remove unnecessary variables
    !
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


subroutine PolyRegularQKTO(iSolnIndex,iParam,xT,dGParam,dPartialGParam)

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

    return

end subroutine PolyRegularQKTO

