
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompChemicalPotential.f90
    !> \brief   Compute the chemical potentials of all solution phase constituents.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      CompMolFraction.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code
    !   06/13/2012      M.H.A. Piro         Include excess parameters
    !   01/18/2013      M.H.A. Piro         The chemical potential terms for each phase component are now
    !                                        computed in the CompExcessGibbsEnergy.f90 subroutine because
    !                                        the equation is model dependent.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the chemical potentials of all solution phase
    !! constituents expected to be stable at equilibrium.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   lCompEverything     A logical scalar indicating whether everything should be computed, or
    !!                                   only what is necessary.  For the most part, it is only necessary to
    !!                                   compute chemical potentials of solution species that are expected to
    !!                                   be stable.
    !
    ! dChemicalPotential                A double real vector representing the chemical potential of a substance.
    ! dGibbsSolnPhase                   A double real vector representing the Gibbs energy of a solution phase.
    ! dPartialExcessGibbs               A double real vector representing the partial molar excess Gibbs energy
    !                                    of mixing of each solution phase constituent.
    ! nSpeciesPhase                     An integer vector representing the index of the last species in a
    !                                   particular solution phase.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CompChemicalPotential(lCompEverything)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k
    real(8) :: dTemp
    logical :: lCompEverything


    ! Initialize variables:
    dChemicalPotential  = 0D0
    dGibbsSolnPhase     = 0D0
    dPartialExcessGibbs = 0D0

    ! Compute the mole fractions of species in solution phases expected to be stable:
    do j = 1, nSolnPhases

        k      = nElements - j + 1          ! Index of phase in iAssemblage
        dTemp  = 1D0 / dMolesPhase(k)
        k      = -iAssemblage(k)            ! Absolute index of solution phase

        ! Compute the mole fractions of all solution phase constituents in the phase:
        do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
            dMolFraction(i) = dMolesSpecies(i) * dTemp
        end do

        ! Compute excess terms:
        call CompExcessGibbsEnergy(k)

    end do

    ! Check if the chemical potentials for everything should be computed:
    if (lCompEverything .EQV. .TRUE.) then

        ! Compute the chemical potentials for every solution phase in the system:
        LOOP_SolnPhasesSys: do j = 1, nSolnPhasesSys

            ! Skip this phase if it is already part of the system:
            if (lSolnPhases(j) .EQV. .TRUE.) cycle LOOP_SolnPhasesSys

            ! Compute the mole fractions of solution phase constituents:
            call CompMolFraction(j)

        end do LOOP_SolnPhasesSys
    end if

    return

end subroutine CompChemicalPotential
