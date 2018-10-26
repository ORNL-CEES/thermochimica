
!-------------------------------------------------------------------------------
!
!> \file    CompOtoMRatio.f90
!> \brief   Compute the oxygent to metal ratio of the UO2 solid solution phase.
!> \author  M.H.A. Piro
!> \date    July 28, 2014
!> \todo    Make the CompOtoMRatio.f90 subroutine more general.
!
!
! Revisions:
! ==========
! 
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   07/28/2014      M.H.A. Piro         Original code
!
!
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to compute the oxygen to metal 
!! ratio of the UO2 solid solution phase.  Note that this value differs from the
!! oxygen to metal ratio of the system, since the system is heterogeneous.  For
!! the latter case, this value necessarily increases as a result of fission; 
!! however, the former is expected to decrease with irradiation since many of the
!! fission products dissolved in this phase have a valence less than 4+.  In 
!! order for the phase to be "stoichiometric", the average valency has to equal
!! 4+ (i.e., for O/M = 2, 2*(O2-) + 1*(M) = 0; thus, M = 4+).
!!
!! The assumption made by this subroutine is that a three sublattice compound
!! energy formalism model is used to represent this phase.  The first sublattice
!! represents the cations on their normal sites, the second sublattice represents
!! oxygen anions on their normal sites (or vacancies) and the third sublattice
!! represents interstitial oxygen anions or vacancies.
!! Of course, this is not entirely general.  For a rainy day...
!
!
!
! Pertinent variables:
! ====================
!
!> \param[in]     cPhaseIn      A character string representing the name of the 
!!                               phase in question.
!> \param[out]    dOtoMRatio    A double real scalar represetning the oxygen to
!!                               metal ratio.
!> \param[out]    INFO          An integer scalar indicating a successful exit (== 0)
!!                               or an error (/= 0).
!!
! Local variables:
! ----------------
! iPhaseID                  An integer scalar representing the absolute index
!                            of the solution phase in quesiton.
! iChargedPhaseID           An integer scalar representing the absolute index
!                            of the sublattice index.
!
!-------------------------------------------------------------------------------

    
subroutine CompOtoMRatio(cPhaseIn, dOtoMRatio, INFO)

    USE ModuleThermo

    implicit none
    
    integer      :: i, j, iPhaseID, iChargedPhaseID, INFO
    character(*) :: cPhaseIn
    real(8)      :: dOtoMRatio


    ! Initialize variables:
    INFO       = 0
    iPhaseID   = 0
    dOtoMRatio = 0D0

    ! Loop through stable solution phases to find the index corresponding to the 
    ! solution phase in question:
    LOOP_CheckPhaseName: do i = 1, nSolnPhases
        j = -iAssemblage(nElements - i + 1)
        ! Check if the solution phase name matches the input variable:
        if (cSolnPhaseName(j) == cPhaseIn) then
            iPhaseID = j
            exit LOOP_CheckPhaseName
        end if
    end do LOOP_CheckPhaseName

    ! Make sure that the
    if ((cSolnPhaseType(j) == 'SUBL').OR.(cSolnPhaseType(j) == 'SUBLM')) then

        ! Store the index # of the charged phase:
        iChargedPhaseID = iPhaseSublattice(j)

        ! Check to make sure that the correct number of sublattices:
        if (nSublatticePhase(iChargedPhaseID) /= 3) then
            INFO = 1
            return
        end if

        dOtoMRatio =  2D0 * dSiteFraction(iChargedPhaseID,2,1) + dSiteFraction(iChargedPhaseID,3,1)

    end if
    
    return
      
end subroutine CompOtoMRatio
