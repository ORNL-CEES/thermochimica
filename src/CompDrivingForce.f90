
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompDrivingForce.f90
    !> \brief   Compute the driving force of all pure condensed phases.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/26/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the driving force of all pure condensed phases in 
    !! the database.  The driving force is defined as the difference between the standard molar Gibbs energy of
    !! a pure condensed phase and the corresponding value computed from the element potentials.  This value is
    !! used to determine whether a pure condensed phase should be added to the system.  For a more thorough 
    !! explanation of the chemical significance of the driving force, refer to the following literature:
    !!
    !! H.L. Lukas, S.G. Fries, B. Sundman, Computational Thermodynamics - The Calphad Method, Cambridge 
    !! University Press, New York, 2007.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[out]  iMaxDrivingForce    An integer scalar representing the index of the pure condensed phase 
    !!                                   with the maximum driving force.  A value of zero is returned by default.
    !> \param[out]  dMaxDrivingForce    A double real scalar representing the maximum driving force of all pure
    !!                                   condensed phases.  A value of zero is returned by default.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine CompDrivingForce(iMaxDrivingForce,dMaxDrivingForce)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none
    
    integer :: i, j, iMaxDrivingForce
    real(8) :: dTemp, dMaxDrivingForce
       
       
    ! Initialize variables:
    iMaxDrivingForce = 0
    dMaxDrivingForce = 0D0
       
    ! Loop through all pure condensed phases:
    do i = nSpeciesPhase(nSolnPhasesSys) + 1, nSpecies - nDummySpecies
    
        ! Compute the chemical potential of this phase as defined by the element potentials:
        dTemp = 0D0
        do j = 1, nElements
            dTemp = dTemp + dElementPotential(j) * dStoichSpecies(i,j)
        end do
    
        ! Compute the driving force:
        dTemp = dStdGibbsEnergy(i) - dTemp
        
        ! Normalize per gram-atom:
        dTemp = dTemp / dSpeciesTotalAtoms(i)
    
        ! Check for most negative value:
        if (dTemp < dMaxDrivingForce) then
            iMaxDrivingForce = i
            dMaxDrivingForce = dTemp
        end if
    
    end do
     
    return 
    
end subroutine CompDrivingForce