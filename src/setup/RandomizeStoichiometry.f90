
!-------------------------------------------------------------------------------------------------------------
    !
    !> \file        RandomizeStoichiometry.f90
    !> \brief       Fuzzy stoichiometry operations.
    !> \details     Randomizes (and resets) stoichiometry to improve performance in cases of underdetermined systems.
    !> \author      Max Poschmann
    !> \sa          Thermochimica.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !   04/21/2022      Max Poschmann   Original code
    !
    !
    !
    !> \param       lFuzzyStoich    A logical that determines if these routines are enabled.
    !> \param       dFuzzMag        A double real scalar representing magnitude of perturbation to the
    !!                                  stoichiometries in thermodynamic database.
    !
!-------------------------------------------------------------------------------------------------------------


subroutine RandomizeStoichiometry

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver, ONLY: lMiscibility

    implicit none

    integer :: i, j, k, n
    integer, allocatable :: lastSeed(:)
    real(8) :: rValue

    ! Fuzzy stoichiometry
    allocate(dStoichSpeciesUnFuzzed(nSpecies,nElements))
    
    call random_seed(size = n)
    if (allocated(lastSeed)) deallocate(lastSeed)
    allocate(lastSeed(n))
    lastSeed = 0

    ! Species in solution phases
    do k = 1, nSolnPhasesSys
        ! Reset the random seed so all phases start the same
        ! This is for miscibility gap phases: could change to only apply to those
        ! call SRAND(7)
        if (lMiscibility(k)) then
            call RANDOM_SEED(put=lastSeed)
        else
            call RANDOM_SEED(get=lastSeed)
        end if
        
        do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
            do j = 1, nElements
                dStoichSpeciesUnFuzzed(i,j) = dStoichSpecies(i,j)
                if (dMolesElement(j) < 1D-5) cycle
                if (dStoichSpecies(i,j) > 0D0) then
                    call RANDOM_NUMBER(rValue)
                    dStoichSpecies(i,j) = dStoichSpecies(i,j) + 2D0 * (rValue - 0.5D0) * dFuzzMag
                end if
                print *, cSpeciesName(i), dStoichSpecies(i,:)
            end do
        end do
    end do

    ! Stoichiometric phases
    do i = nSpeciesPhase(nSolnPhasesSys) + 1, nSpecies
        do j = 1, nElements
            dStoichSpeciesUnFuzzed(i,j) = dStoichSpecies(i,j)
            if (dMolesElement(j) < 1D-5) cycle
            if (dStoichSpecies(i,j) > 0D0) then
                call RANDOM_NUMBER(rValue)
                dStoichSpecies(i,j) = dStoichSpecies(i,j) + 2D0 * (rValue - 0.5D0) * dFuzzMag
            end if
        end do
    end do

    return

end subroutine RandomizeStoichiometry

subroutine ResetStoichiometry

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: i, j
    logical :: lPhasePass

    ! Un-fuzz stoichiometry
    do i = 1, nSpecies
        do j = 1, nElements
            dStoichSpecies(i,j) = dStoichSpeciesUnFuzzed(i,j)
        end do
    end do
    deallocate(dStoichSpeciesUnFuzzed)

    ! Recompute with reset stoichiometry
    do j = nSolnPhases, 1, -1
        if (dMolesPhase(nElements - j + 1) < dFuzzMag * 1D3) then
            call RemSolnPhase(j,lPhasePass)
        end if
    end do
    do j = nConPhases, 1, -1
        if (dMolesPhase(j) < dFuzzMag * 1D3) then
            call RemPureConPhase(j,lPhasePass,lPhasePass)
        end if
    end do

    return

end subroutine ResetStoichiometry
