
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckPhaseAssemblageID.f90
    !> \brief   Check whether a phase assemblage that is being considered is the same as the current phase
    !!           assemblage.
    !> \author  M.H.A. Piro
    !> \date    March 17, 2013
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/17/2013      M.H.A. Piro         Original code.
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check whether a phase assemblage that is being considered
    !!  is the same as the current phase assemblage.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]  iTempVec    An integer vector representing a phase assemblage that is being considered.
    !> \param[out] lPhasePass  A logical scalar indicating whether the iTempVec and iAssemblage vectors contain
    !!                          the same set of phases (TRUE) or if they differ (FALSE).
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckPhaseAssemblageID(iTempVec, lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                       :: i, j
    integer, dimension(nElements) :: iTempVec
    logical                       :: lPhasePass


    ! Initialize variables:
    lPhasePass = .TRUE.

    ! If the sum of the absolute values of the coefficients are the same, then check each individual coefficient:
    IF_SUM: if (SUM(ABS(iTempVec)) == SUM(ABS(iAssemblage))) then

        ! Loop through coefficients of the iAssemblage vector:
        LOOP_Outer: do j = 1, nElements

            ! Loop through coefficients of the iTempVec vector:
            LOOP_Inner: do i = 1, nElements

                ! Check if the coefficients are the same:
                if (iTempVec(i) == iAssemblage(j)) cycle LOOP_Outer

            end do LOOP_Inner

            ! The iTempVec vector does not contain the phase corresponding to coefficient j in iAssemblage.
            ! Report a failure and exit:
            lPhasePass = .FALSE.
            exit LOOP_Outer

        end do LOOP_Outer

    else
        ! The sum of the absolute values of the coefficients are not the same.
        lPhasePass = .FALSE.

    end if IF_SUM

    return

end subroutine CheckPhaseAssemblageID
