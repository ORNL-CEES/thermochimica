
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GibbsEnergy.f90
    !> \brief   Calculates the Gibbs energy of a given input system.
    !> \author  M. Poschmann
    !> \date    July 24, 2019
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/24/2019      M. Poschmann         Original code.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to print the results to screen in a style similar to FactSage.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nConPhasesIn        An integer scalar representing the number of stable pure condensed phases input.
    ! nSolnPhasesIn       An integer scalar representing the number of stable solution phases input.
    ! iAssemblageIn       An integer vector representing the absolute indices of stable phases input.
    ! dMolesPhaseIn       A double real vector representing the number of moles for each stable phase in the system.
    ! dMolFractionIn      A double real vector representing the mole fraction of each species in the system.
    ! dGibbsEnergySysOut  A double real scalar representing the integral Gibbs energy of the system. This is the output variable.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine GibbsEnergy(nConPhasesIn, nSolnPhasesIn, iAssemblageIn, dMolesPhaseIn, dMolesSpeciesIn, dGibbsEnergySysOut)

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer, intent(in)  :: nConPhasesIn, nSolnPhasesIn
    integer, intent(in)  :: iAssemblageIn(nElements)
    real(8), intent(in)  :: dMolesPhaseIn(nElements), dMolesSpeciesIn(nSpecies)
    real(8), intent(out) :: dGibbsEnergySysOut
    integer              :: i, j, k, iPhaseIndex

    dGibbsEnergySysOut = 0D0

    ! Only proceed for a successful calculation:
    IF_PASS: if (INFOThermo == 0) then

        do i = 1, nSolnPhasesIn
            k = nElements - i + 1
            iPhaseIndex = -iAssemblageIn(k)
            do j = nSpeciesPhase(iPhaseIndex - 1) + 1, nSpeciesPhase(iPhaseIndex)
                dGibbsEnergySysOut = dGibbsEnergySysOut + dChemicalPotential(j) * dMolesSpeciesIn(j)
            end do
        end do
        do i = 1, nConPhasesIn
            iPhaseIndex = iAssemblageIn(i)
            do j = 1, nElements
                dGibbsEnergySysOut = dGibbsEnergySysOut + dElementPotential(j) * dStoichSpecies(iPhaseIndex,j) * dMolesPhaseIn(i)
            end do
        end do
        dGibbsEnergySysOut = dGibbsEnergySysOut * dIdealConstant * dTemperature
        print '(A26,ES12.5,A4)', ' Integral Gibbs energy = ', dGibbsEnergySysOut, ' [J]'

    end if IF_PASS

    return

end subroutine GibbsEnergy


    !---------------------------------------------------------------------------
    !                      END - GibbsEnergy.f90
    !---------------------------------------------------------------------------
