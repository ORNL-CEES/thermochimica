
!-------------------------------------------------------------------------------
!
!> \file    GetOutputOinUO2.f90
!> \brief   Get oxygen in UO2+x.
!> \author  M. Poschmann
!> \date    February 6, 2020
!
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   06/02/2020      M. Poschmann        Original code
!
!
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to get specific the amount of
!! oxygen in the UO2+x phase.
!
! Pertinent variables:
! ====================
!
!> \param[out]    dMolO                 A double real scalar representing the
!!                                       number of moles of O in UO2+x.
!> \param[out]    INFO                  An integer scalar indicating a successful
!!                                       exit (== 0) or an error (/= 0).
!
!-------------------------------------------------------------------------------


subroutine GetOutputOinUO2(dMolO, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer,       intent(out)   :: INFO
    integer                      :: i, j, iOxyIndex
    real(8),       intent(out)   :: dMolO

    ! Initialize variables:
    INFO  = 0
    dMolO = 0D0
    iOxyIndex = 0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then
        do i = 1, nElements
            if (cElementName(i) =='O') then
                iOxyIndex = i
            end if
        end do

        if (iOxyIndex == 0) then
           ! This solution phase was not found.  Report an error:
           INFO = 1
        else
            do i = 1, nSolnPhasesSys
                if (cSolnPhaseName(i) == 'O2ZRU_C') then
                    do j = nSpeciesPhase(i - 1) + 1, nSpeciesPhase(i)
                        dMolO = dMolO + dStoichSpecies(j,iOxyIndex) * dMolesSpecies(j)
                    end do
                end if
            end do
        end if

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetOutputOinUO2
