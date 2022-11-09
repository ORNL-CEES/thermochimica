subroutine GetPhaseNames(cPhaseNamesOut)

    USE ModuleParseCS

    implicit none

    integer :: i, j
    character(30), intent(out), dimension(*) :: cPhaseNamesOut

    do i = 1, nSolnPhasesSysCS
        cPhaseNamesOut(i) = cSolnPhaseNameCS(i)
    end do

    do j = MAXVAL(nSpeciesPhaseCS) + 1, nSpeciesCS
        i = nSolnPhasesSysCS + j - MAXVAL(nSpeciesPhaseCS)
        cPhaseNamesOut(i) = cSpeciesNameCS(j)
    end do

end subroutine GetPhaseNames
