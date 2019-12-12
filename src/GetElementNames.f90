subroutine GetElementNames(cElementNamesOut)

    USE ModuleParseCS

    implicit none

    integer :: i
    character(3), intent(out), dimension(*) :: cElementNamesOut

    do i = 1, nElementsCS
        cElementNamesOut(i) = cElementNameCS(i)
    end do

end subroutine GetElementNames
