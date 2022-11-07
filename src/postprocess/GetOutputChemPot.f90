
!-------------------------------------------------------------------------------
!
!> \file    GetOutputChemPot.f90
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
!   09/16/2015      M.H.A. Piro         Original code
!
!
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to get specific thermodynamic
!! output from an equilibrium calculation.
!
!
!
! Pertinent variables:
! ====================
!
!> \param[inout]  cElementNameRequest   A three letter character string
!!                                       representing the name of a chemical
!!                                       element.
!> \param[out]    dElementChemPot       A double real scalar representing the
!!                                       chemical potential of this element.
!> \param[out]    INFO                  An integer scalar indicating a successful
!!                                       exit (== 0) or an error (/= 0).
!
!-------------------------------------------------------------------------------


subroutine GetOutputChemPot(cElementNameRequest, lcElementNameRequest, dElementChemPot, INFO) &
    bind(C, name="TCAPI_getOutputChemPot")

    USE ModuleThermo
    USE ModuleThermoIO
    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer,      intent(out)   :: INFO
    real(8),      intent(out)   :: dElementChemPot
    character(kind=c_char,len=1), target, intent(in)         :: cElementNameRequest(*)
    integer(c_size_t), intent(in), value                     :: lcElementNameRequest
    character(kind=c_char,len=lcElementNameRequest), pointer :: cElementNameUse
    integer                     :: i, j

    call c_f_pointer(cptr=c_loc(cElementNameRequest), fptr=cElementNameUse)

    ! Initialize variables:
    INFO            = 0
    dElementChemPot = 0D0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Loop through elements to find the one corresponding to the element
        ! being requested:
        j = 0
        LOOP_A: do i = 1, nElements
            if (cElementNameUse == cElementName(i)) then
                j = i
                exit LOOP_A
            end if
        end do LOOP_A

        ! Check to make sure that the element was found:
        if (j /= 0) then
            ! The element was found in the list.  Store the chemical potential
            ! of this element and convert back to units of J/g-at:
            dElementChemPot = dElementPotential(j) * dIdealConstant * dTemperature
        else
            ! This element was not found.  Report an error:
            INFO = 1
        end if

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetOutputChemPot
