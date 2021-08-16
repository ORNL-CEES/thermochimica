
!-------------------------------------------------------------------------------
!
!> \file    GetElementMolesInPhase.f90
!> \brief   Get moles of element in phase.
!> \author  M. Poschmann
!> \date    April 21, 2021
!
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   21/04/2021      M. Poschmann        Original code
!
!
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to get the number of moles of
!! an element in phase.
!
!
!
! Pertinent variables:
! ====================
!
!> \param[in]     cElement              A character string representing the
!!                                       element.
!> \param[in]     cPhase                A character string representing the
!!                                       phase.
!> \param[out]    dMolesOut             A double real scalar representing the
!!                                       number of moles of element in phase.
!> \param[out]    INFO                  An integer scalar indicating a successful
!!                                       exit (== 0) or an error (/= 0).
!!
!
!-------------------------------------------------------------------------------


subroutine GetElementMolesInPhase(cElement, lcElement, cPhase, lcPhase, dMolesOut, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer,       intent(out)   :: INFO
    integer                      :: i, j, k, iPhaseInd
    real(8),       intent(out)   :: dMolesOut
    character(*), intent(in)    :: cPhase
    character(*),  intent(in)    :: cElement
    integer, intent(in)          :: lcPhase, lcElement
    character(25)                :: cTempPhase, cTempElement, cSearchPhase, cSearchElement

    cSearchPhase = cPhase !TRIM(cPhase(1:min(25,lcPhase)))
    cSearchPhase = TRIM(cSearchPhase(1:min(25,lcPhase)))
    cSearchElement = cElement! TRIM(cElement(1:min(3,lcElement)))
    cSearchElement = TRIM(cSearchElement(1:min(3,lcElement)))


    ! Initialize variables:
    INFO            = 0
    dMolesOut = 0D0
    k = 0

    LOOP_Elements: do i = 1, nElements
        cTempElement = ADJUSTL(cElementName(i))
        if (cTempElement == cSearchElement) then
            k = i
            exit LOOP_Elements
        end if
    end do LOOP_Elements

    if (k == 0) then
       ! This solution phase was not found.  Report an error:
       INFO = 1
       return
    end if

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then
        ! Remove trailing blanks:
        LOOP_PHASE: do i = 1, nSolnPhases

            iPhaseInd = -iAssemblage(nElements + 1 - i)
            ! iPhaseInd = iPhaseReq
            ! Remove leading blanks:
            cTempPhase = ADJUSTL(cSolnPhaseName(iPhaseInd))

            ! Loop through species in this phase:
            if (cTempPhase == cSearchPhase) then
                ! Solution phase found.  Now look for species.
                LOOP_Species: do j = nSpeciesPhase(iPhaseInd-1) + 1, nSpeciesPhase(iPhaseInd)
                    dMolesOut = dMolesOut + dMolesSpecies(j) * dStoichSpecies(j,k)
                end do LOOP_Species
                exit LOOP_PHASE
            end if

        end do LOOP_PHASE

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetElementMolesInPhase
