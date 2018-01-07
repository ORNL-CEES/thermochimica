
!-------------------------------------------------------------------------------
!
!> \file    GetPureConPhaseMol.f90
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
!! in a specific pure condensed phase.
!
!
!
! Pertinent variables:
! ====================
!
!> \param[in]     cPureConOut           A character string represnting the pure
!!                                       condensed phase name.
!> \param[out]    dPureConMolOut           A double real scalar representing the
!!                                       number of moles of said phase.
!> \param[out]    INFO                  An integer scalar indicating a successful
!!                                       exit (== 0) or an error (/= 0).
!!
!
!-------------------------------------------------------------------------------


subroutine GetPureConPhaseMol(cPureConOut, dPureConMolOut, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer,       intent(out)   :: INFO
    integer                      :: i, j, k
    real(8),       intent(out)   :: dPureConMolOut
    character(*),  intent(in)    :: cPureConOut
    character(25)                :: cTemp
    character(25)                :: cPureTemp
    cPureTemp = cPureConOut(1:min(25,len(cPureConOut)))

    ! Initialize variables:
    INFO            = 0
    dPureConMolOut   = 0D0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Remove trailing blanks:
        ! cPureConOut    = TRIM(cPureConOut)
        cPureTemp    = TRIM(cPureTemp)

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        j = 0
        LOOP_PureStable: do i = 1, nConPhases
            k = iAssemblage(i)
            ! Remove leading blanks:
            cTemp = ADJUSTL(cSpeciesName(k))
            ! if (cPureConOut == cTemp) then
            if (cPureTemp == cTemp) then
                ! Pure condensed phase found.  Record integer index and exit loop.
               j = i
               exit LOOP_PureStable
            end if

        end do LOOP_PureStable

        ! Check to make sure that the solution phase was found:
        IF_PureStable: if (j /= 0) then
           dPureConMolOut = dMolesPhase(j)
        else
           ! This solution phase was not found.  Report an error:
           INFO = 1
        end if IF_PureStable

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetPureConPhaseMol
