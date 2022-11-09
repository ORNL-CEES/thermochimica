
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


subroutine GetElementMolesInPhase(cElement, lcElement, cPhase, lcPhase, dMolesOut, INFO) &
    bind(C, name="TCAPI_getElementMolesInPhase")

    USE ModuleThermo
    USE ModuleThermoIO
    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer,       intent(out)   :: INFO
    real(8),       intent(out)   :: dMolesOut
    character(kind=c_char,len=1), target, intent(in) :: cPhase(*), cElement(*)
    integer(c_size_t), intent(in), value             :: lcPhase, lcElement
    character(kind=c_char,len=lcPhase), pointer      :: fPhase
    character(kind=c_char,len=lcElement), pointer    :: fElement
    integer                      :: i, j, k, iPhaseInd
    character(30)                :: cTempPhase, cTempElement

    call c_f_pointer(cptr=c_loc(cPhase), fptr=fPhase)
    call c_f_pointer(cptr=c_loc(cElement), fptr=fElement)

    ! Initialize variables:
    INFO            = 0
    dMolesOut = 0D0
    k = 0

    LOOP_Elements: do i = 1, nElements
        cTempElement = ADJUSTL(cElementName(i))
        if (cTempElement == fElement) then
            k = i
            exit LOOP_Elements
        end if
    end do LOOP_Elements

    if (k == 0) then
        ! This element was not found.  Report an error:
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
            if (cTempPhase == fPhase) then
                ! Solution phase found.  Now look for species.
                LOOP_Species: do j = nSpeciesPhase(iPhaseInd-1) + 1, nSpeciesPhase(iPhaseInd)
                    dMolesOut = dMolesOut + dMolesSpecies(j) * dStoichSpecies(j,k)
                end do LOOP_Species
                exit LOOP_PHASE
            end if
            INFO = 2
        end do LOOP_PHASE

        if (INFO == 2) then
            INFO = 0
            LOOP_Stoich: do i = 1, nConPhases
                iPhaseInd = iAssemblage(i)
                cTempPhase = ADJUSTL(cSpeciesName(iPhaseInd))
                if (cTempPhase == fPhase) then
                    dMolesOut = dMolesPhase(i) * dStoichSpecies(iPhaseInd,k)
                    exit LOOP_Stoich
                end if
                INFO = 2
            end do LOOP_Stoich
        end if

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetElementMolesInPhase
