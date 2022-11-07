
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


subroutine GetPureConPhaseMol(cPureConOut, lcPureConOut, dPureConMolOut, INFO) &
    bind(C, name="TCAPI_getPureConPhaseMol")

    USE ModuleThermo
    USE ModuleThermoIO
    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer,       intent(out)   :: INFO
    real(8),       intent(out)   :: dPureConMolOut
    character(kind=c_char,len=1), target, intent(in) :: cPureConOut(*)
    integer(c_size_t), intent(in), value             :: lcPureConOut
    character(kind=c_char,len=lcPureConOut), pointer :: fPureConOut
    integer                      :: i, j, k
    character(30)                :: cTemp

    call c_f_pointer(cptr=c_loc(cPureConOut), fptr=fPureConOut)

    ! Initialize variables:
    INFO            = 0
    dPureConMolOut   = 0D0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        j = 0
        LOOP_PureStable: do i = 1, nConPhases
            k = iAssemblage(i)
            ! Remove leading blanks:
            cTemp = ADJUSTL(cSpeciesName(k))
            ! if (cPureConOut == cTemp) then
            if (fPureConOut == cTemp) then
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
