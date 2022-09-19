
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

    implicit none

    integer :: i, j

    ! Fuzzy stoichiometry
    if (INFOThermo == 0 .AND. lFuzzyStoich .AND. (.NOT. lRetryAttempted)) then
        allocate(dStoichSpeciesUnFuzzed(nSpecies,nElements))
        do i = 1, nSpecies
            do j = 1, nElements
                if (dMolesElement(j) < 1D-5) cycle
                dStoichSpeciesUnFuzzed(i,j) = dStoichSpecies(i,j)
                if (dStoichSpecies(i,j) > 0D0) then
                    dStoichSpecies(i,j) = dStoichSpecies(i,j) + 2D0 * (RAND(0) - 0.5D0) * dFuzzMag
                end if
            end do
        end do
    end if

    return

end subroutine RandomizeStoichiometry

subroutine ResetStoichiometry

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: i, j
    logical :: lPhasePass

    ! Un-fuzz stoichiometry
    if (INFOThermo == 0 .AND. lFuzzyStoich .AND. (.NOT. lRetryAttempted)) then
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
    end if

    return

end subroutine ResetStoichiometry
