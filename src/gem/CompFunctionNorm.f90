
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompFunctionNorm.f90
    !> \brief   Compute the functional norm for the line search.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      GEMLineSearch.f90
    !> \sa      CompStoichSolnPhase.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the functional norm for the line search algorithm to
    !! determine whether the system is converging sufficinetly or diverging. The functional vector is not directly
    !! computed because it is not needed. Thus, the functional norm is computed directly.  This term incorporates
    !! the residuals of the mass balance equations, the average residual between the chemical potential of each
    !! species and the corresponding value computed from the element potentials.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dGEMFunctionNorm      A double real scalar representing the functional norm of the functional  vector.  The
    !                        functional vector represents the relative errors of the mass balance equations,
    !                        the departure of chemical potentials of stable phases from the corresponding values
    !                        computed from the element potentials.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CompFunctionNorm

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModulePhaseConstraints

    implicit none

    integer :: i, j, k, l, c, cIdx, nReal, infoConstraint
    real(8) :: dNormComponent, dSumStoich, dTemp, dAchieved

    ! Initialize variables:
    dGEMFunctionNorm    = 0D0
    dEffStoichSolnPhase = 0D0

    ! Compute the residuals of the mass balance equations for each element:
    do j = 1, nElements
        dNormComponent = dMolesElement(j)
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements - i + 1)       ! Absolute solution phase index
            ! Compute the stoichiometry of this solution phase:
            call CompStoichSolnPhase(k)
            dNormComponent = dNormComponent - dEffStoichSolnPhase(k,j) * dMolesPhase(nElements - i + 1)
        end do
        do i = 1,nConPhases
            dNormComponent = dNormComponent - dMolesPhase(i) * dStoichSpecies(iAssemblage(i),j)
        end do
        dGEMFunctionNorm = dGEMFunctionNorm + (dNormComponent)**(2)
    end do

    ! Compute the residuals of the Gibbs energy difference between each solution phase and the element potentials:
    do l = 1, nSolnPhases
        k      = -iAssemblage(nElements - l + 1)       ! Absolute solution phase index
        cIdx = 0
        if (nPhaseConstraints > 0) then
            if (lPhaseConstrainedSoln(k)) then
                do c = 1, nPhaseConstraints
                    if ((iPhaseConstraintKind(c) == 0) .AND. (iPhaseConstraintID(c) == k)) then
                        cIdx = c
                        exit
                    end if
                end do
            end if
        end if
        do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
            dNormComponent = 0D0
            do j = 1, nElements
                dNormComponent = dNormComponent + dElementPotential(j) * dStoichSpecies(i,j)
            end do
            if (cIdx > 0) then
                nReal = nElements - nChargedConstraints
                if (nReal < 1) nReal = nElements
                dSumStoich = 0D0
                do j = 1, nReal
                    dSumStoich = dSumStoich + dStoichSpecies(i,j)
                end do
                dNormComponent = dNormComponent + dPhaseConstraintLambda(cIdx) * dSumStoich
            end if
            ! Normalize the residual term by the number of particles per formula mass:
            dNormComponent = dNormComponent / DFLOAT(iParticlesPerMole(i))
            ! Compute the residual term weighted by the mole fraction:
            dNormComponent = DABS(dChemicalPotential(i) - dNormComponent) * dMolFraction(i)
            dGEMFunctionNorm = dGEMFunctionNorm + (dNormComponent)**(2)
        end do
    end do

    ! Compute the residuals of the chemical potentials of pure condensed phases and the element potentials:
    do i = 1, nConPhases
        k     = iAssemblage(i)
        dNormComponent = 0D0
        do j = 1, nElements
            dNormComponent = dNormComponent + dElementPotential(j) * dStoichSpecies(k,j)
        end do
        if (nPhaseConstraints > 0) then
            if (lPhaseConstrainedCon(k)) then
                cIdx = 0
                do c = 1, nPhaseConstraints
                    if ((iPhaseConstraintKind(c) == 1) .AND. (iPhaseConstraintID(c) == k)) then
                        cIdx = c
                        exit
                    end if
                end do
                if (cIdx > 0) then
                    call GetPhaseConstraintCoefficient(cIdx, dSumStoich, infoConstraint)
                    if (infoConstraint == 0) then
                        dNormComponent = dNormComponent + dPhaseConstraintLambda(cIdx) * dSumStoich
                    else
                        dGEMFunctionNorm = dGEMFunctionNorm + 1D0
                    end if
                end if
            end if
        end if
        dNormComponent            = dNormComponent - dStdGibbsEnergy(k)
        dGEMFunctionNorm = dGEMFunctionNorm + (dNormComponent)**(2)
    end do

    if (nPhaseConstraints > 0) then
        do c = 1, nPhaseConstraints
            call GetPhaseConstraintResult(c, dAchieved, dTemp, infoConstraint)
            if (infoConstraint == 0) then
                dGEMFunctionNorm = dGEMFunctionNorm + dTemp**2
            else
                dGEMFunctionNorm = dGEMFunctionNorm + 1D0
            end if
        end do
    end if

    ! Finally, the functional norm:
    dGEMFunctionNorm = dGEMFunctionNorm**(0.5)

    if (lDebugMode) print *, 'dGEMFunctionNorm = ', dGEMFunctionNorm

    return

end subroutine CompFunctionNorm
