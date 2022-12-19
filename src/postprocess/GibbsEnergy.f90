
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


subroutine GibbsEnergy(nConPhasesIn, nSolnPhasesIn, iAssemblageIn, dMolesPhaseIn, dMolFractionIn, dGibbsEnergySysOut)

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer, intent(in)  :: nConPhasesIn, nSolnPhasesIn
    integer, intent(in)  :: iAssemblageIn(nElements)
    real(8), intent(in)  :: dMolesPhaseIn(nElements), dMolFractionIn(nSpecies)
    real(8), intent(out) :: dGibbsEnergySysOut
    integer              :: i, j, k, iPhaseIndex
    real(8)              :: dGibbsPhase, dMolFractionPhase, dGibbsEnergySysTemp

    ! Temporary variables to store things that must be overwritten
    integer  :: nConPhasesTemp, nSolnPhasesTemp
    integer  :: iAssemblageTemp(nElements)
    real(8)  :: dMolesPhaseTemp(nElements), dMolFractionTemp(nSpecies), dMolesSpeciesTemp(nSpecies)

    nConPhasesTemp = nConPhases
    nSolnPhasesTemp = nSolnPhases
    iAssemblageTemp = iAssemblage
    dMolesPhaseTemp = dMolesPhase
    dMolFractionTemp = dMolFraction
    dMolesSpeciesTemp = dMolesSpecies
    dGibbsEnergySysTemp = dGibbsEnergySys

    nConPhases = nConPhasesIn
    nSolnPhases = nSolnPhasesIn
    iAssemblage = iAssemblageIn
    dMolesPhase = dMolesPhaseIn
    dMolFraction = dMolFractionIn

    dMolesSpecies = 0D0
    do i =1, nSolnPhasesIn
        k = nElements - i + 1
        iPhaseIndex = -iAssemblageIn(k)
        do j = nSpeciesPhase(iPhaseIndex - 1) + 1, nSpeciesPhase(iPhaseIndex)
            dMolesSpecies(j) = dMolesPhaseIn(k) * dMolFractionIn(j)
        end do
    end do

    call CompChemicalPotential(.TRUE.)

    dGibbsEnergySysOut = 0D0
    ! Only proceed for a successful calculation:
    IF_PASS: if (INFOThermo == 0) then

        do i = 1, nSolnPhasesIn
            k = nElements - i + 1
            iPhaseIndex = -iAssemblageIn(k)
            dGibbsPhase = 0D0
            dMolFractionPhase = 0D0
            do j = nSpeciesPhase(iPhaseIndex - 1) + 1, nSpeciesPhase(iPhaseIndex)
                dGibbsPhase = dGibbsPhase + dChemicalPotential(j) * dMolesPhase(k) * dMolFraction(j)
                dMolFractionPhase = dMolFractionPhase + dMolFraction(j)
            end do
            dGibbsEnergySysOut = dGibbsEnergySysOut + dGibbsPhase / dMolFractionPhase
        end do
        do i = 1, nConPhasesIn
            iPhaseIndex = iAssemblageIn(i)
            dGibbsEnergySysOut = dGibbsEnergySysOut + dStdGibbsEnergy(iPhaseIndex) * dMolesPhaseIn(i)
        end do
        dGibbsEnergySysOut = dGibbsEnergySysOut * dIdealConstant * dTemperature
        dGibbsEnergySys = dGibbsEnergySysOut
        ! call PrintResults
        dGibbsEnergySys = dGibbsEnergySysTemp
        print '(A26,ES12.5,A4)', ' Gibbs energy of input  = ', dGibbsEnergySysOut, ' [J]'
        print '(A26,ES12.5,A4)', ' Gibbs energy diff      = ', dGibbsEnergySysOut - dGibbsEnergySys, ' [J]'

    end if IF_PASS

    ! Return original values to variables
    nConPhases = nConPhasesTemp
    nSolnPhases = nSolnPhasesTemp
    iAssemblage = iAssemblageTemp
    dMolesPhase = dMolesPhaseTemp
    dMolFraction = dMolFractionTemp
    dMolesSpecies = dMolesSpeciesTemp

    return

end subroutine GibbsEnergy


    !---------------------------------------------------------------------------
    !                      END - GibbsEnergy.f90
    !---------------------------------------------------------------------------
