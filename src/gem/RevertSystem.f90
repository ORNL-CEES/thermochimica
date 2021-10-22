
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    RevertSystem.f90
    !> \brief   Revert the system to a previously considered phase assemblage.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   02/28/2012      M.H.A. Piro         Original code
    !   03/08/2012      M.H.A. Piro         Added the capability to revert to any previously considered phase
    !                                       assemblage, instead of just the last phase assemblage.
    !   05/02/2012      M.H.A. Piro         Add the capability to revert to a particular phase assemblage if
    !                                       requested.
    !   06/01/2012      M.H.A. Piro         Fix bug: the number of moles of solution species was not updated.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to revert the system to the last successfully tested phase
    !! assemblage.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]  iterSpecific     An integer scalar representing a specific global iteration to revert to.
    !
    ! nElements             An integer scalar representing the number of elements in the system.
    ! nConPhases            An integer scalar representing the number of pure condensed phases in the system.
    ! nSolnPhases           An integer scalar representing the number of solution phases in the system.
    ! iterRevert            An integer scalar representing the last iteration that the system was reverted to.
    ! iterlast              An integer scalar representing the last iteration when the system was changed.
    ! iAssemblage           An integer vector representing the indices of all phases in the current system.
    ! iterHistory           An integer matrix representing the indices of all phases in the assemblage throughout
    !                        the iteration history.
    ! dMolesPhase           A double real vector representing the number of moles of each phase.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine RevertSystem(iterSpecific)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::                       i, j, k, l, iterStart, iterSpecific
    integer,dimension(nElements)::  iWork, iWorkB
    real(8),dimension(nElements)::  dTempVec
    logical:: lPhasePass


    if (lDebugMode) print *, 'RevertSystem'

    ! Do not revert the system if the system has not yet changed:
    if (iterlast <= 0) return

    ! Initialize variables:
    dTempVec      = dMolesPhase
    dMolesPhase   = 0D0
    lRevertSystem = .FALSE.

   ! dElementPotential = dElementPotentialRevert
    if ((iterSpecific == 2) .AND. iterGlobal > 2000) then
        iAssemblage = iAssemblageBest
        dMolesPhase = dMolesPhaseBest
        dMolFraction = dMolFractionBest
        dElementPotential = dElementPotentialBest

        nConPhases               = 0
        nSolnPhases              = 0
        lSolnPhases              = .FALSE.

        ! Count the number of pure condensed phases and the number of solution phases in the system:
        do j = 1, nElements
            if (iAssemblage(j) > 0) then
                nConPhases = nConPhases + 1
            elseif (iAssemblage(j) < 0) then
                nSolnPhases = nSolnPhases + 1

                k = -iAssemblage(j)

                ! Compute the number of moles of each species:
                do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    dMolesSpecies(i) = dMolesPhase(j) * dMolFraction(i)
                end do

                lSolnPhases(-iAssemblage(j)) = .TRUE.
                dMolesPhaseLast(j) = dMolesPhase(j)
            else
                ! If iAssemblage(j) = 0, then it is an empty placeholder.
            end if

        end do

        iterLast = iterGlobal

    elseif (iterSpecific > 0) then
        ! Revert to a specific iteration:

        iWork(1:nElements) = iterHistory(1:nElements,iterSpecific)

        ! If the row of iterHistory is full of zeros, return:
        if (SUM(ABS(iWork)) == 0) then
            dMolesPhase = dTempVec
            return
        end if

        ! Check if this phase assemblage is the same as the current phase assemblage:
        call CheckPhaseAssemblageID(iWork, lPhasePass)

        ! Return if they are the same:
        if (lPhasePass) then
            dMolesPhase = dTempVec
            return
        end if

        ! Map the number of moles of consistent phases:
        LOOP_Outter: do j = 1, nElements
            LOOP_Inner: do k = 1, nElements
                if (iAssemblage(k) == iWork(j)) then
                    dMolesPhase(j) = dTempVec(k)
                end if
            end do LOOP_Inner
        end do LOOP_Outter

        iAssemblage(1:nElements) = iWork(1:nElements)
        nConPhases               = 0
        nSolnPhases              = 0
        lSolnPhases              = .FALSE.

        ! Count the number of pure condensed phases and the number of solution phases in the system:
        do j = 1, nElements
            if (iAssemblage(j) > 0) then
                nConPhases = nConPhases + 1
            elseif (iAssemblage(j) < 0) then
                nSolnPhases = nSolnPhases + 1
                if (dMolesPhase(j) == 0D0) dMolesPhase(j) = 1D0
                dMolesPhase(j) = MAX(dMolesPhase(j), 1D0)

                k = -iAssemblage(j)

                ! Compute the mole fractions of solution species in this phase:
                call CompMolFraction(k)

                ! Compute the number of moles of each species:
                do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    dMolesSpecies(i) = dMolesPhase(j) * dMolFraction(i)
                end do

                lSolnPhases(-iAssemblage(j)) = .TRUE.
                dMolesPhaseLast(j) = dMolesPhase(j)
            else
                ! If iAssemblage(j) = 0, then it is an empty placeholder.
            end if

        end do

        iterLast = iterGlobal

    else
        ! Revert to the last iteration:

        ! Determine the iteration # to start checking:
        if (iterRevert == 0) then
            iterStart = iterlast
        elseif (iterRevert == 1) then
            iterStart = iterlast
        else
            iterStart = iterRevert - 1
        end if

        ! Search for the iteration to revert to:
        LOOP_A: do i = iterStart, 1, - 1
            iWork(1:nElements) = iterHistory(1:nElements,i)
            iWorkB(1:nElements) = iterHistory(1:nElements,i+1) - iterHistory(1:nElements,i)

            if (SUM(iWorkB) == 0) then
                cycle LOOP_A
            end if

            if ((SUM(iWork) /= 0).AND.(SUM(iWork(1:nElements) - iAssemblage(1:nElements)) /= 0)) then
                ! The system will revert to this iteration.

                ! Map the number of moles of consistent phases:
                LOOP_B: do j = 1, nElements
                    LOOP_C: do k = 1, nElements
                        if (iAssemblage(k) == iWork(j)) then
                            dMolesPhase(j) = dTempVec(k)
                        end if
                    end do LOOP_C
                end do LOOP_B

                iAssemblage(1:nElements) = iWork(1:nElements)
                nConPhases               = 0
                nSolnPhases              = 0
                lSolnPhases              = .FALSE.

                ! Count the number of pure condensed phases and the number of solution phases in the system:
                do j = 1, nElements
                    if (iAssemblage(j) > 0) then
                        nConPhases = nConPhases + 1
                    elseif (iAssemblage(j) < 0) then
                        nSolnPhases = nSolnPhases + 1
                        if (dMolesPhase(j) == 0D0) dMolesPhase(j) = 1D0
                        dMolesPhase(j) = MAX(dMolesPhase(j), 1D0)

                        k = -iAssemblage(j)

                        ! Compute the mole fractions of solution species in this phase:
                        call CompMolFraction(k)

                        ! Compute the number of moles of each species:
                        do l = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                            dMolesSpecies(l) = dMolesPhase(j) * dMolFraction(l)
                        end do

                        lSolnPhases(-iAssemblage(j)) = .TRUE.
                    else
                        ! If iAssemblage(j) = 0, then it is an empty placeholder.
                    end if

                end do

                iterlast = iterGlobal

                ! Store the iteration # that the system is reverted to:
                iterRevert = i

                return
            end if

            iterRevert = 1

        end do LOOP_A

        ! The system was not reverted, return dMolesPhase to its original value:
        dMolesPhase = dTempVec

    end if

    dMolesPhase = DABS(dMolesPhase)

    return

end subroutine RevertSystem
