
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    PostProcess.f90
    !> \brief   Perform post-processing of results.
    !> \author  M.H.A. Piro
    !> \date    January 14, 2013
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/14/2013      M.H.A. Piro         Original code.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to perform post-procesing of results.
    !
    !
    ! Pertinent variables:
    ! ====================

    !
    !-------------------------------------------------------------------------------------------------------------


subroutine PostProcess

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: i, j


    ! Initiate variables:
    dGibbsEnergySys = 0d0

    ! Multiply the number of moles of all phases by the normalizing constant:
    dNormalizeInput = 1D0 / dNormalizeInput
    dMolesPhase     = dMolesPhase   * dNormalizeInput * dMassScale
    dMolesElement   = dMolesElement * dNormalizeInput * dMassScale
    dMolesSpecies   = dMolesSpecies * dNormalizeInput * dMassScale
    dElementMass    = dElementMass  * dMassScale

    ! Write the moles of condensed phases to corresponding species
    conCheck: do i = 1, nConPhases
        do while ((dMolesPhase(i) == 0D0) .AND. (nConPhases > 0))
            do j = i + 1, nConPhases
                iAssemblage(j - 1) = iAssemblage(j)
                dMolesPhase(j - 1) = dMolesPhase(j)
            end do
            iAssemblage(nConPhases) = 0
            dMolesPhase(nConPhases) = 0D0
            nConPhases = nConPhases - 1
            if (i > nConPhases) exit conCheck
        end do
        if (dMolesPhase(i) > 0D0) dMolesSpecies(iAssemblage(i)) = dMolesPhase(i)
    end do conCheck

    ! Compute the integral Gibbs energy of the system:
    do i = 1, nElements
        dGibbsEnergySys = dGibbsEnergySys + dElementPotential(i) * dMolesElement(i)
    end do
    dGibbsEnergySys = dGibbsEnergySys * dIdealConstant * dTemperature

    dNormalizeInput = 1D0 / dNormalizeInput
    return

end subroutine PostProcess
