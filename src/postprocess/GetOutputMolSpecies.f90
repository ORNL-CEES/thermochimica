
!-------------------------------------------------------------------------------
!
!> \file    GetOutputMolSpecies.f90
!> \brief   Get specific thermodynamic output.
!> \author  S. Simunovic
!> \date    Dec. 20, 2017
!
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   12/20/2017      S. Simunovic        Original code
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
!> \param[in]     cSpeciesOut           A character string representing the
!!                                       species.
!> \param[out]    dMolFractionOut       A double real scalar representing the
!!                                       mole fraction of said species.
!> \param[out]    dMolSpecies           A double real scalar representing the
!!                                       moles of said species
!> \param[out]    INFO                  An integer scalar indicating a successful
!!                                       exit (== 0) or an error (/= 0).
!!
!
!-------------------------------------------------------------------------------


subroutine GetOutputMolSpecies(cSpeciesOut, lcSpeciesOut, dMolFractionOut, dMolesOut, INFO) &
    bind(C, name="TCAPI_getOutputMolSpecies")

    USE ModuleThermo
    USE ModuleThermoIO
    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer,       intent(out)   :: INFO
    real(8),       intent(out)   :: dMolFractionOut, dMolesOut
    character(kind=c_char,len=1), target, intent(in) :: cSpeciesOut(*)
    integer(c_size_t), intent(in), value             :: lcSpeciesOut
    character(kind=c_char,len=lcSpeciesOut), pointer :: cSearch
    integer                      :: i, k
    character(30)                :: cTemp

    call c_f_pointer(cptr=c_loc(cSpeciesOut), fptr=cSearch)

    ! Initialize variables:
    INFO            = 0
    dMolFractionOut = 0D0
    dMolesOut       = 0D0
    k=0
    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then
        LOOP_SPECIES: do i = 1, nSpecies

           ! Remove leading blanks:
           cTemp = ADJUSTL(cSpeciesName(i))

           ! Loop through species in this phase:
           if (cTemp == cSearch) then
              ! Solution species found.  Record index and exit loop.
              k=i
              dMolFractionOut = dMolFraction(k)
              dMolesOut       = dMolesSpecies(k)
              exit LOOP_SPECIES
           end if

        end do LOOP_SPECIES

        IF_FOUND: if ( k==0 ) then
           ! This solution phase was not found.  Report an error:
           INFO = 1
        end if IF_FOUND

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetOutputMolSpecies
