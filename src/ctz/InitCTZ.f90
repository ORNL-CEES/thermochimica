subroutine InitCTZ(inputRange,inputMaxElements,inputMaxAssemblages)

    USE ModuleCTZ

    implicit none

    real(8), intent(in) :: inputRange
    integer, intent(in) :: inputMaxElements,inputMaxAssemblages

    tRange  = inputRange
    maxNorm = 1D-16
    nMaxElements = inputMaxElements
    nMaxAssemblages = inputMaxAssemblages

    if (lCtzInit) call ResetCTZ

    allocate(assemblageHistory(nMaxAssemblages+1,nMaxElements),assemblageTlimits(nMaxAssemblages,2))
    allocate(stoichHistory(nMaxAssemblages,2,nMaxElements,nMaxElements))
    allocate(elementHistory(nMaxAssemblages+1,nMaxElements))
    assemblageHistory      = 0
    assemblageTlimits(:,1) = 1D5
    assemblageTlimits(:,2) = 0D0
    elementHistory         = ''
    stoichHistory          = 0D0
    nAssemblages           = 0
    lCtzInit = .TRUE.

end subroutine InitCTZ
