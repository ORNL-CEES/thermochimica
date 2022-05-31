subroutine MakePhaseExclusionList

    USE ModuleThermoIO
    USE ModuleParseCS

    implicit none

    integer :: i, j

    ! If there is an "excluded except" list, then add all other phases to exclusion list
    if (nPhasesExcludedExcept > 0) then
        ! Solution phases
        loop_checkExclusionSolution: do i = 1, nSolnPhasesSysCS
            ! Check if phase is on exception list
            do j = 1, nPhasesExcludedExcept
                if (cSolnPhaseNameCS(i) == cPhasesExcludedExcept(j)) cycle loop_checkExclusionSolution
            end do
            ! Check if phase is on exclusion list
            do j = 1, nPhasesExcluded
                if (cSolnPhaseNameCS(i) == cPhasesExcluded(j)) cycle loop_checkExclusionSolution
            end do
            ! If not, add to exclusion list
            nPhasesExcluded = nPhasesExcluded + 1
            cPhasesExcluded(nPhasesExcluded) = cSolnPhaseNameCS(i)
        end do loop_checkExclusionSolution

        ! Pure condensed phases
        loop_checkExclusionPureCondensed: do i = nSpeciesPhaseCS(nSolnPhasesSysCS) + 1, nSpeciesCS
            ! Check if phase is on exception list
            do j = 1, nPhasesExcludedExcept
                if (cSpeciesNameCS(i) == cPhasesExcludedExcept(j)) cycle loop_checkExclusionPureCondensed
            end do
            ! Check if phase is on exclusion list
            do j = 1, nPhasesExcluded
                if (cSpeciesNameCS(i) == cPhasesExcluded(j)) cycle loop_checkExclusionPureCondensed
            end do
            ! If not, add to exclusion list
            nPhasesExcluded = nPhasesExcluded + 1
            cPhasesExcluded(nPhasesExcluded) = cSpeciesNameCS(i)
        end do loop_checkExclusionPureCondensed
    end if


end subroutine MakePhaseExclusionList
