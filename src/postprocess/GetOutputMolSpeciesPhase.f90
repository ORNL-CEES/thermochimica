
!-------------------------------------------------------------------------------
!
!> \file    GetOutputMolSpeciesPhase.f90
!> \brief   Get specific thermodynamic output.
!> \author  M. Poschmann
!> \date    September 9, 2019
!
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   09/09/2019      M. Poschmann        Original code
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


subroutine GetOutputMolSpeciesPhase(cPhase, lcPhase, cSpecies, lcSpecies, dMolFractionOut, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer,       intent(out)   :: INFO
    integer                      :: i, j, k, iPhaseInd
    real(8),       intent(out)   :: dMolFractionOut
    character(*),  intent(in)    :: cPhase, cSpecies
    integer                      :: lcPhase, lcSpecies
    character(30)                :: cTempPhase, cTempSpecies, cSearchPhase, cSearchSpecies

    cSearchPhase = cPhase!TRIM(cPhase(1:min(30,lcPhase)))
    cSearchSpecies = cSpecies!TRIM(cSpecies(1:min(30,lcSpecies)))

    cSearchPhase = TRIM(cSearchPhase(1:min(30,lcPhase)))
    cSearchSpecies = TRIM(cSearchSpecies(1:min(30,lcSpecies)))

    ! Initialize variables:
    INFO            = 0
    dMolFractionOut = 0D0
    k = 0
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
                    cTempSpecies = ADJUSTL(cSpeciesName(j))
                    if (cTempSpecies == cSearchSpecies) then
                        ! Found species in phase, return mole fraction
                        dMolFractionOut = dMolFraction(j)
                        k = j
                        exit LOOP_Species
                    end if
                end do LOOP_Species
                exit LOOP_PHASE
            end if

        end do LOOP_PHASE

        if (k == 0) then
           ! This solution phase was not found.  Report an error:
           INFO = 1
        end if

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetOutputMolSpeciesPhase
