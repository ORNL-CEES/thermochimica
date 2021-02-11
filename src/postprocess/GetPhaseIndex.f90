
!-------------------------------------------------------------------------------
!
!> \file    GetSolnPhaseMol.f90
!> \brief   Get specific thermodynamic output.
!> \author  M.H.A. Piro
!> \date    Sept. 16, 2015
!
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   11/10/2016      S.Simunovic         Original code
!
!
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to get number of moles
!! in a specific solution phase.
!
!
!
! Pertinent variables:
! ====================
!
!> \param[in]     cPhaseName              A character string representing the solution
!!                                       phase name.
!> \param[out]    iIndexOut              An integer representing the index of the phase.
!> \param[out]    INFO                  An integer scalar indicating a successful
!!                                       exit (== 0) or an error (/= 0).
!!
!
!-------------------------------------------------------------------------------


subroutine GetPhaseIndex(cPhaseName, lcPhaseName, iIndexOut, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer,       intent(out)   :: INFO
    integer                      :: i, k
    integer,       intent(out)   :: iIndexOut
    character(*),  intent(in)    :: cPhaseName
    integer                      :: lcPhaseName
    character(25)                :: cTemp, cTempPhase
    cTemp = cPhaseName(1:min(25,lcPhaseName))

    ! Initialize variables:
    INFO          = 0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Remove trailing blanks:
        ! cSolnOut    = TRIM(cSolnOut)
        cTemp    = TRIM(cTemp)
        cTemp    = ADJUSTL(cTemp)

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        iIndexOut = 0
        LOOP_ASSEMBLAGE: do i = 1, nElements
            k = iAssemblage(i)
            if (k < 0) then
                cTempPhase = TRIM(cSolnPhaseName(-k))
                cTempPhase = ADJUSTL(cTempPhase)
                if (cTemp == cTempPhase) then
                    ! Phase found.  Record integer index and exit loop.
                    iIndexOut = i
                    exit LOOP_ASSEMBLAGE
                end if
            elseif (k > 0) then
                cTempPhase = TRIM(cSpeciesName(k))
                cTempPhase = ADJUSTL(cTempPhase)
                if (cTemp == cTempPhase) then
                    ! Pure condensed phase found.  Record integer index and exit loop.
                    iIndexOut = i
                    exit LOOP_ASSEMBLAGE
                end if
           end if

        end do LOOP_ASSEMBLAGE
    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetPhaseIndex
