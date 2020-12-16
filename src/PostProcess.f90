
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

    integer :: i


    ! Initiate variables:
    dGibbsEnergySys = 0d0

    ! Multiply the number of moles of all phases by the normalizing constant:
    dNormalizeInput = 1D0 / dNormalizeInput
    dMolesPhase     = dMolesPhase   * dNormalizeInput * dMassScale
    dMolesElement   = dMolesElement * dNormalizeInput * dMassScale
    dMolesSpecies   = dMolesSpecies * dNormalizeInput * dMassScale
    dElementMass    = dElementMass  * dNormalizeInput * dMassScale

    ! Write the moles of condensed phases to corresponding species
    do i = 1, nConPhases
        dMolesSpecies(iAssemblage(i)) = dMolesPhase(i)
    end do

    ! Compute the integral Gibbs energy of the system:
    do i = 1, nElements
        dGibbsEnergySys = dGibbsEnergySys + dElementPotential(i) * dMolesElement(i)
    end do
    dGibbsEnergySys = dGibbsEnergySys * dIdealConstant * dTemperature

    return

end subroutine PostProcess
