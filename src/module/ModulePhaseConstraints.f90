
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
    logical :: lPhaseConstraintResultsAvailable = .FALSE.
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

    pure function NormalizePhaseConstraintName(cPhaseIn) result(cNameOut)
        implicit none

        character(*), intent(in) :: cPhaseIn
        character(25)            :: cNameOut

        cNameOut = ''
        if (len_trim(cPhaseIn) > 0) then
            cNameOut = trim(adjustl(cPhaseIn(1:min(25,len_trim(cPhaseIn)))))
        end if
    end function NormalizePhaseConstraintName


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
        lPhaseConstraintResultsAvailable = .FALSE.

    end subroutine ClearPhaseConstraints


    subroutine AddPhaseFractionConstraint(cPhaseIn, dFractionIn)
        implicit none

        character(*), intent(in) :: cPhaseIn
        real(8), intent(in)       :: dFractionIn

        integer :: nNew, i
        character(25) :: cName
        character(25), allocatable :: cNamesNew(:)
        real(8), allocatable       :: dTargetsNew(:)

        cName = NormalizePhaseConstraintName(cPhaseIn)
        if (len_trim(cName) == 0) return

        lPhaseConstraintResultsAvailable = .FALSE.

        if (nPhaseConstraints > 0) then
            do i = 1, nPhaseConstraints
                if (cPhaseConstraintName(i) == cName) then
                    dPhaseConstraintTarget(i) = dFractionIn
                    return
                end if
            end do
        end if

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

        cNamesNew(nNew) = cName
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
        character(25) :: cCandidate
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
            cSearch = NormalizePhaseConstraintName(cPhaseConstraintName(i))
            lFound = .FALSE.

            ! Check solution phases first
            do j = 1, nSolnPhasesSys
                cCandidate = NormalizePhaseConstraintName(cSolnPhaseName(j))
                if (cCandidate == cSearch) then
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
                        cCandidate = NormalizePhaseConstraintName(cSpeciesName(k))
                        if (cCandidate == cSearch) then
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
        integer :: i
        real(8) :: dTotal

        INFO = 0

        if (nPhaseConstraints <= 0) return

        call GetPhaseConstraintTotalElementMoles(dTotal, INFO)
        if (INFO /= 0) return

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
        integer :: i, nReal
        real(8) :: dSum, dTol

        INFO = 0

        if (nPhaseConstraints <= 0) return

        nReal = nElements - nChargedConstraints
        if (nReal < 1) nReal = nElements
        if (nPhaseConstraints > nReal) then
            INFO = 4
            return
        end if

        dTol = 1D-6
        dSum = 0D0
        do i = 1, nPhaseConstraints
            if (dPhaseConstraintTarget(i) /= dPhaseConstraintTarget(i)) then
                INFO = 5
                return
            end if
            if ((dPhaseConstraintTarget(i) < 0D0) .OR. (dPhaseConstraintTarget(i) > 1D0)) then
                INFO = 5
                return
            end if
            dSum = dSum + dPhaseConstraintTarget(i)
        end do

        if (DABS(dSum - 1D0) > dTol) then
            INFO = 3
            return
        end if

    end subroutine ValidatePhaseConstraints


    subroutine GetPhaseConstraintTotalElementMoles(dTotal, INFO)
        USE ModuleThermo
        implicit none

        real(8), intent(out) :: dTotal
        integer, intent(out) :: INFO
        integer :: j, nReal

        INFO = 0
        dTotal = 0D0
        nReal = nElements - nChargedConstraints
        if (nReal < 1) nReal = nElements

        do j = 1, nReal
            dTotal = dTotal + dMolesElement(j)
        end do

        if ((dTotal <= 0D0) .OR. (dTotal /= dTotal)) INFO = 1

    end subroutine GetPhaseConstraintTotalElementMoles


    subroutine GetPhaseConstraintCoefficient(iConstraint, dCoefficient, INFO)
        USE ModuleThermo
        USE ModuleGEMSolver
        implicit none

        integer, intent(in) :: iConstraint
        real(8), intent(out) :: dCoefficient
        integer, intent(out) :: INFO
        integer :: i, k, nReal

        INFO = 0
        dCoefficient = 0D0

        if (.NOT. allocated(iPhaseConstraintKind) .OR. .NOT. allocated(iPhaseConstraintID)) then
            INFO = 1
            return
        end if
        if ((iConstraint < 1) .OR. (iConstraint > nPhaseConstraints)) then
            INFO = 1
            return
        end if

        nReal = nElements - nChargedConstraints
        if (nReal < 1) nReal = nElements
        k = iPhaseConstraintID(iConstraint)

        select case (iPhaseConstraintKind(iConstraint))
        case (0)
            call CompStoichSolnPhase(k)
            do i = 1, nReal
                dCoefficient = dCoefficient + dEffStoichSolnPhase(k,i)
            end do
        case (1)
            do i = 1, nReal
                dCoefficient = dCoefficient + dStoichSpecies(k,i)
            end do
        case default
            INFO = 1
            return
        end select

        if ((dCoefficient <= 0D0) .OR. (dCoefficient /= dCoefficient)) INFO = 1

    end subroutine GetPhaseConstraintCoefficient


    subroutine GetPhaseConstraintResult(iConstraint, dAchieved, dResidual, INFO)
        USE ModuleThermo
        implicit none

        integer, intent(in) :: iConstraint
        real(8), intent(out) :: dAchieved, dResidual
        integer, intent(out) :: INFO
        integer :: i, k, iSlot
        real(8) :: dCoefficient, dTotal, dPhaseMoles

        INFO = 0
        dAchieved = 0D0
        dResidual = 0D0

        call GetPhaseConstraintTotalElementMoles(dTotal, INFO)
        if (INFO /= 0) return
        call GetPhaseConstraintCoefficient(iConstraint, dCoefficient, INFO)
        if (INFO /= 0) return

        k = iPhaseConstraintID(iConstraint)
        iSlot = 0
        if (iPhaseConstraintKind(iConstraint) == 0) then
            do i = nElements - nSolnPhases + 1, nElements
                if (iAssemblage(i) == -k) then
                    iSlot = i
                    exit
                end if
            end do
        else
            do i = 1, nConPhases
                if (iAssemblage(i) == k) then
                    iSlot = i
                    exit
                end if
            end do
        end if

        dPhaseMoles = 0D0
        if (iSlot > 0) dPhaseMoles = dMolesPhase(iSlot)
        dAchieved = dPhaseMoles * dCoefficient / dTotal
        dResidual = dAchieved - dPhaseConstraintTarget(iConstraint)

    end subroutine GetPhaseConstraintResult


    subroutine GetPhaseConstraintData(iConstraint, dTarget, dAchieved, dResidual, dLambda, INFO)
        USE ModuleThermo, ONLY: dMolesPhase
        USE ModuleThermoIO, ONLY: INFOThermo
        implicit none

        integer, intent(in) :: iConstraint
        real(8), intent(out) :: dTarget, dAchieved, dResidual, dLambda
        integer, intent(out) :: INFO
        integer :: infoResult

        INFO = 0
        dTarget = 0D0
        dAchieved = 0D0
        dResidual = 0D0
        dLambda = 0D0

        if ((iConstraint < 1) .OR. (iConstraint > nPhaseConstraints)) then
            INFO = 1
            return
        end if

        if ((.NOT. lPhaseConstraintResultsAvailable) .OR. (INFOThermo /= 0) .OR. &
            (.NOT. allocated(iPhaseConstraintKind)) .OR. &
            (.NOT. allocated(iPhaseConstraintID)) .OR. (.NOT. allocated(dPhaseConstraintLambda)) .OR. &
            (.NOT. allocated(dMolesPhase))) then
            INFO = 2
            return
        end if

        call GetPhaseConstraintResult(iConstraint, dAchieved, dResidual, infoResult)
        if (infoResult /= 0) then
            dAchieved = 0D0
            dResidual = 0D0
            INFO = 2
            return
        end if

        dTarget = dPhaseConstraintTarget(iConstraint)
        dLambda = dPhaseConstraintLambda(iConstraint)

    end subroutine GetPhaseConstraintData




end module ModulePhaseConstraints
