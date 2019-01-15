
!-------------------------------------------------------------------------------
!
!> \file    GetOutputSolnSpecies.f90
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
!> \param[in]     cSolnOut              A character string represnting the solution 
!!                                       phase name.
!> \param[in]     cSpeciesOut           A character string representing the 
!!                                       species.
!> \param[out]    dMolFractionOut       A double real scalar representing the
!!                                       mole fraction of said species.
!> \param[out]    dChemPotSpecies       A double real scalar representing the 
!!                                       chemical potential of said species.
!> \param[out]    INFO                  An integer scalar indicating a successful 
!!                                       exit (== 0) or an error (/= 0).
!!
!
!-------------------------------------------------------------------------------

    
subroutine GetOutputSolnSpecies(cSolnOut, cSpeciesOut, dMolFractionOut, dChemPotSpecies, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer,       intent(out)   :: INFO
    integer                      :: i, j, k
    real(8),       intent(out)   :: dMolFractionOut, dChemPotSpecies
    character(*), intent(inout) :: cSolnOut, cSpeciesOut
    character(25)                :: cTemp


    ! Initialize variables:
    INFO            = 0
    dMolFractionOut = 0D0
    dChemPotSpecies = 0D0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        j = 0
        LOOP_SOLN: do i = 1, nSolnPhases
            k = -iAssemblage(nElements - i + 1)

            if (cSolnOut == cSolnPhaseName(k)) then
                ! Solution phase found.  Record integer index and exit loop.
                j = k
                exit LOOP_SOLN
            end if

        end do LOOP_SOLN

        ! Check to make sure that the solution phase was found:
        IF_SOLN: if (j /= 0) then
            k = 0
            ! Solution phase found.  Now, look for the species in this phase.
            LOOP_SPECIES: do i = nSpeciesPhase(j-1) + 1, nSpeciesPhase(j)

                ! Remove leading blanks:
                cTemp = ADJUSTL(cSpeciesName(i))

                ! Loop through species in this phase:
                if (cTemp == cSpeciesOut) then
                    ! Solution species found.  Record index and exit loop.
                    k = i
                    exit LOOP_SPECIES
                end if

            end do LOOP_SPECIES

            ! Check if the solution species was found:
            if (k /= 0 ) then
                ! Solution species found.
                dMolFractionOut = dMolFraction(k)

                ! Compute the chemical potential of this solution species
                do j = 1, nElements
                    dChemPotSpecies = dChemPotSpecies + dElementPotential(j) * dStoichSpecies(k,j)
                end do

                ! Convert to units of J/mol:
                dChemPotSpecies = dChemPotSpecies * dIdealConstant * dTemperature

            else
                ! Solution species not found.  Record an error:
                INFO = 2
            end if
        else
            ! This solution phase was not found.  Report an error:
            INFO = 1
        end if IF_SOLN

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return
      
end subroutine GetOutputSolnSpecies
