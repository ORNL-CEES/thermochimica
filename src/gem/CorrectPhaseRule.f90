
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CorrectPhaseRule.f90
    !> \brief   In the event that the Phase Rule has been violated, provide appropriate correction.
    !> \author  M.H.A. Piro
    !> \date    February 4, 2013.
    !> \sa      CheckPhaseAssemblage.f90
    !> \sa      RemSolnPhase.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   02/04/2013      M.H.A. Piro         Original code.
    !
    ! Purpose:
    ! ========
    !
    !> \details In the event that the Phase Rule has been violated, a phase must be removed from the system
    !! in order for the system to achieve true thermodynamic equilibrium.  Generally, the maximum number of
    !! phases considered in Thermochimica is constrained by the number of system components.  This becomes more
    !! convoluted when dealing with phases containing ionic species because the electron is then considered as an
    !! additional system component, which must be accompanied by an additional constraint on the system (i.e.,
    !! charge neutrality).
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[out] lPhaseChange A logical scalar indicating whether a phase change has occured.
    !
    ! nSolnPhases           The number of solution phases in the assemblage.
    ! nElements             The number of elements in the system.
    ! lPhasePass            A logical variable indicating whether the new estimated phase assemblage passed
    !                        (.TRUE.) or failed (.FALSE.).
    ! dMolesPhase           A double real vector representing the number of moles of each phase.
    ! iTempVec              A integer vector that is temporarily used for sorting purposes.
    ! dTempVec              A double real vector that is temporarily used for sorting purposes.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CorrectPhaseRule(lPhaseChange)

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModulePhaseConstraints

    implicit none

    integer                            :: i, j, k
    integer, dimension(:), allocatable :: iTempVec
    real(8), dimension(:), allocatable :: dTempVec
    logical                            :: lPhasePass, lPhaseChange


    ! Only proceed if the Phase Rule has been violated:
    IF_PhaseRule: if (nSolnPhases + nConPhases > nElements - nChargedConstraints) then

        ! Deallocate allocatable arrays if necessary:
        if (allocated(iTempVec)) deallocate(iTempVec)
        if (allocated(dTempVec)) deallocate(dTempVec)

        ! Allocate memory:
        allocate(iTempVec(nSolnPhases), dTempVec(nSolnPhases))

        ! Initialize variables:
        lPhasePass   = .FALSE.
        lPhaseChange = .FALSE.
        dTempVec     = 0D0
        iTempVec     = 0

        ! Assign molar quantities to temporary vector:
        do i = 1, nSolnPhases
            j           = nElements - i + 1
            dTempVec(i) = -dMolesPhase(j)
        end do

        ! Sort the solution phases:
        call SortPick(nSolnPhases,dTempVec,iTempVec)

        ! Loop through stable solution phases:
        LOOP_SolnRem: do i = 1, nSolnPhases

            j = iTempVec(i)

            ! Skip removal attempts for constrained solution phases:
            if (nPhaseConstraints > 0) then
                k = -iAssemblage(nElements - j + 1)
                if (k > 0) then
                    if (lPhaseConstrainedSoln(k)) cycle LOOP_SolnRem
                end if
            end if

            ! Try removing this solution phase from the system:
            call RemSolnPhase(j, lPhasePass)

            ! Exit if the phase has been successfully removed:
            if (lPhasePass) then
                lPhaseChange = .TRUE.
                exit LOOP_SolnRem
            end if

        end do LOOP_SolnRem

        ! Deallocate allocatable arrays:
        deallocate(iTempVec, dTempVec)

    end if IF_PhaseRule

    return

end subroutine CorrectPhaseRule
