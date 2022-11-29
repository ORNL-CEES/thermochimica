
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
!> \param[in]     cSolnOut              A character string represnting the solution
!!                                       phase name.
!> \param[out]    dSolnMolOut           A double real scalar representing the
!!                                       number of moles of said solution phase.
!> \param[out]    INFO                  An integer scalar indicating a successful
!!                                       exit (== 0) or an error (/= 0).
!!
!
!-------------------------------------------------------------------------------


subroutine GetSolnPhaseMol(cSolnOut, dSolnMolOut, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer,       intent(out)   :: INFO
    integer                      :: i, j, k, iph
    real(8),       intent(out)   :: dSolnMolOut
    character(*),  intent(in)    :: cSolnOut
    character(25)                :: cTemp
    cTemp = cSolnOut

    ! Initialize variables:
    INFO            = 0
    dSolnMolOut   = 0D0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Remove trailing blanks:
        ! cSolnOut    = TRIM(cSolnOut)
        cTemp    = TRIM(cTemp)

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        j = 0
        LOOP_SOLN: do i = 1, nSolnPhases
           iph = nElements - i + 1
           k = -iAssemblage(iph)
           ! if (cSolnOut == cSolnPhaseName(k)) then
           if (cTemp == cSolnPhaseName(k)) then
              ! Solution phase found.  Record integer index and exit loop.
              j = iph
              exit LOOP_SOLN
           end if

        end do LOOP_SOLN

        ! Check to make sure that the solution phase was found:
        IF_SOLN: if (j /= 0) then
           dSolnMolOut = dMolesPhase(iph)
        else
           ! This solution phase was not found.  Report an error:
           INFO = 1
        end if IF_SOLN

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetSolnPhaseMol
