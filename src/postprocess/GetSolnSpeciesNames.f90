subroutine GetSolnSpeciesNames(cSolnSpeciesNamesOut)

    USE ModuleParseCS

    implicit none

    integer :: i, j, k
    character(30), intent(out), dimension(*) :: cSolnSpeciesNamesOut

    k = 0
    do i = 1, nSolnPhasesSysCS
        do j = nSpeciesPhaseCS(i-1) + 1, nSpeciesPhaseCS(i)
            k = k + 1
            cSolnSpeciesNamesOut(k) = cSpeciesNameCS(j)
        end do
    end do

end subroutine GetSolnSpeciesNames
