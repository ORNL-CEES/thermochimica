
!-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ModulePhaseConstraints.f90
    !> \brief   Store and manage phase fraction constraints.
    !> \author  Thermochimica contributors
    !
    !-------------------------------------------------------------------------------------------------------------

module ModulePhaseConstraints

    implicit none

    SAVE

    integer :: nPhaseConstraints = 0
    real(8) :: dPhaseConstraintPenalty = 1D4

    character(25), allocatable :: cPhaseConstraintName(:)
    real(8), allocatable       :: dPhaseConstraintTarget(:)
    real(8), allocatable       :: dPhaseConstraintElemTarget(:)
    real(8), allocatable       :: dPhaseConstraintLambda(:)
    real(8), allocatable       :: dPhaseConstraintLambdaLast(:)
    integer, allocatable       :: iPhaseConstraintKind(:)   ! 0=solution phase, 1=pure condensed
    integer, allocatable       :: iPhaseConstraintID(:)     ! absolute phase index

    logical, allocatable       :: lPhaseConstrainedSoln(:)  ! size nSolnPhasesSys
    logical, allocatable       :: lPhaseConstrainedCon(:)   ! size nSpecies (pure condensed subset)

contains

    subroutine ClearPhaseConstraints()
        implicit none

        if (allocated(cPhaseConstraintName)) deallocate(cPhaseConstraintName)
        if (allocated(dPhaseConstraintTarget)) deallocate(dPhaseConstraintTarget)
        if (allocated(dPhaseConstraintElemTarget)) deallocate(dPhaseConstraintElemTarget)
        if (allocated(dPhaseConstraintLambda)) deallocate(dPhaseConstraintLambda)
        if (allocated(dPhaseConstraintLambdaLast)) deallocate(dPhaseConstraintLambdaLast)
        if (allocated(iPhaseConstraintKind)) deallocate(iPhaseConstraintKind)
        if (allocated(iPhaseConstraintID)) deallocate(iPhaseConstraintID)
        if (allocated(lPhaseConstrainedSoln)) deallocate(lPhaseConstrainedSoln)
        if (allocated(lPhaseConstrainedCon)) deallocate(lPhaseConstrainedCon)

        nPhaseConstraints = 0

    end subroutine ClearPhaseConstraints


    subroutine AddPhaseFractionConstraint(cPhaseIn, dFractionIn)
        implicit none

        character(*), intent(in) :: cPhaseIn
        real(8), intent(in)       :: dFractionIn

        integer :: nNew, i
        character(25), allocatable :: cNamesNew(:)
        real(8), allocatable       :: dTargetsNew(:)

        if (len_trim(cPhaseIn) == 0) return

        nNew = nPhaseConstraints + 1

        allocate(cNamesNew(nNew))
        allocate(dTargetsNew(nNew))

        cNamesNew = ''
        dTargetsNew = 0D0

        if (nPhaseConstraints > 0) then
            do i = 1, nPhaseConstraints
                cNamesNew(i) = cPhaseConstraintName(i)
                dTargetsNew(i) = dPhaseConstraintTarget(i)
            end do
        end if

        cNamesNew(nNew) = trim(adjustl(cPhaseIn(1:min(25,len_trim(cPhaseIn)))))
        dTargetsNew(nNew) = dFractionIn

        if (allocated(cPhaseConstraintName)) deallocate(cPhaseConstraintName)
        if (allocated(dPhaseConstraintTarget)) deallocate(dPhaseConstraintTarget)

        call move_alloc(cNamesNew, cPhaseConstraintName)
        call move_alloc(dTargetsNew, dPhaseConstraintTarget)

        nPhaseConstraints = nNew

    end subroutine AddPhaseFractionConstraint


    subroutine ResolvePhaseConstraints(INFO)
        USE ModuleThermo
        implicit none

        integer, intent(out) :: INFO
        integer :: i, j, k
        character(25) :: cSearch
        logical :: lFound

        INFO = 0

        if (nPhaseConstraints <= 0) return

        if (allocated(iPhaseConstraintKind)) deallocate(iPhaseConstraintKind)
        if (allocated(iPhaseConstraintID)) deallocate(iPhaseConstraintID)

        allocate(iPhaseConstraintKind(nPhaseConstraints))
        allocate(iPhaseConstraintID(nPhaseConstraints))

        iPhaseConstraintKind = -1
        iPhaseConstraintID = -1

        if (allocated(lPhaseConstrainedSoln)) deallocate(lPhaseConstrainedSoln)
        if (allocated(lPhaseConstrainedCon)) deallocate(lPhaseConstrainedCon)

        allocate(lPhaseConstrainedSoln(nSolnPhasesSys))
        allocate(lPhaseConstrainedCon(nSpecies))

        lPhaseConstrainedSoln = .FALSE.
        lPhaseConstrainedCon = .FALSE.

        do i = 1, nPhaseConstraints
            cSearch = trim(adjustl(cPhaseConstraintName(i)))
            lFound = .FALSE.

            ! Check solution phases first
            do j = 1, nSolnPhasesSys
                if (trim(adjustl(cSolnPhaseName(j))) == cSearch) then
                    iPhaseConstraintKind(i) = 0
                    iPhaseConstraintID(i) = j
                    lPhaseConstrainedSoln(j) = .TRUE.
                    lFound = .TRUE.
                    exit
                end if
            end do

            ! Check pure condensed phases
            if (.NOT. lFound) then
                do k = 1, nSpecies
                    if (iPhase(k) == 0) then
                        if (trim(adjustl(cSpeciesName(k))) == cSearch) then
                            iPhaseConstraintKind(i) = 1
                            iPhaseConstraintID(i) = k
                            lPhaseConstrainedCon(k) = .TRUE.
                            lFound = .TRUE.
                            exit
                        end if
                    end if
                end do
            end if

            if (.NOT. lFound) then
                INFO = 1
                return
            end if
        end do

        ! Check for duplicate phase constraints
        do i = 1, nPhaseConstraints
            do j = i + 1, nPhaseConstraints
                if ((iPhaseConstraintKind(i) == iPhaseConstraintKind(j)) .AND. &
                    (iPhaseConstraintID(i) == iPhaseConstraintID(j))) then
                    INFO = 2
                    return
                end if
            end do
        end do

    end subroutine ResolvePhaseConstraints


    subroutine UpdatePhaseConstraintTargets(INFO)
        USE ModuleThermo
        implicit none

        integer, intent(out) :: INFO
        integer :: i, j, nRealElements
        real(8) :: dTotal

        INFO = 0

        if (nPhaseConstraints <= 0) return

        nRealElements = nElements - nChargedConstraints
        if (nRealElements < 1) nRealElements = nElements

        dTotal = 0D0
        do j = 1, nRealElements
            dTotal = dTotal + dMolesElement(j)
        end do

        if (allocated(dPhaseConstraintElemTarget)) deallocate(dPhaseConstraintElemTarget)
        allocate(dPhaseConstraintElemTarget(nPhaseConstraints))

        do i = 1, nPhaseConstraints
            dPhaseConstraintElemTarget(i) = dPhaseConstraintTarget(i) * dTotal
        end do

    end subroutine UpdatePhaseConstraintTargets


    subroutine ValidatePhaseConstraints(INFO)
        USE ModuleThermo
        implicit none

        integer, intent(out) :: INFO
        integer :: i
        real(8) :: dSum

        INFO = 0

        if (nPhaseConstraints <= 0) return

        if (nPhaseConstraints > nElements) then
            INFO = 4
            return
        end if

        dSum = 0D0
        do i = 1, nPhaseConstraints
            dSum = dSum + dPhaseConstraintTarget(i)
        end do

        if (DABS(dSum - 1D0) > 1D-6) then
            INFO = 3
            return
        end if

    end subroutine ValidatePhaseConstraints




end module ModulePhaseConstraints
