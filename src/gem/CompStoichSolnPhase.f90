
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompStoichSolnPhase.f90
    !> \brief   Compute the stoichiometry of a solution phase.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
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
    !> \details The purpose of this subroutine is to compute the effective stoichiometry of a particular
    !! solution phase. This calculation is independent of a solution phase model type.
    !
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   k           An integer scalar representing the absolute solution phase index.
    !
    ! dEffStoichSolnPhase       A double real matrix representing the effective stoichiometry of a
    !                            solution phase (element,solution phase index).
    ! nSpeciesPhase             An integer vector representing the index of the last species in a
    !                            particular solution phase.
    ! iParticlesPerMole         The number of particles per mole of constituent.
    ! dStoichSpecies            The number of atoms of a particular element for a particular species.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CompStoichSolnPhase(k)

    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::   i, j, k, m, n


    ! Initialize variables:
    dEffStoichSolnPhase(k,1:nElements) = 0D0

    m = nSpeciesPhase(k - 1) + 1        ! First species in phase
    n = nSpeciesPhase(k)                ! Last  species in phase

    ! Compute the effective stoichiometry:
    LOOP_A: do i = m, n

        LOOP_B: do j = 1,nElements
            dEffStoichSolnPhase(k,j) = dEffStoichSolnPhase(k,j) + dMolFraction(i) * dStoichSpecies(i,j) &
                / DFLOAT(iParticlesPerMole(i))

            ! Check for a NAN:
            if (dEffStoichSolnPhase(k,j) /= dEffStoichSolnPhase(k,j)) then
                INFOThermo = 24

                exit LOOP_A
            end if

        end do LOOP_B

    end do LOOP_A

    return

end subroutine CompStoichSolnPhase
