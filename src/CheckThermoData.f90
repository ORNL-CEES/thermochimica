
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckThermoData.f90
    !> \brief   Check the thermodynamic database for pure species.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      Thermochimica.f90
    !
    ! References:
    ! ===========
    !
    ! For further information regarding this method, refer to the following material:
    !
    !        M.H.A. Piro, "Computation of Thermodynamic Equilibria Pertinent to Nuclear Materials
    !        in Multi-Physics Codes," PhD Dissertation, Royal Military College of Canada, 2011.
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer       Description of change
    !    ----          ----------       ---------------------
    !    09/02/2011    M.H.A. Piro      Original code
    !    10/20/2011    M.H.A. Piro      Updated modules
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to ensure that the thermodynamic database is appropriate for
    !! use by Thermochimica.  Specifically, at least one pure species of each element must be present
    !! in the database.  A nonzero value of the integer scalar INFOThermo is returned if the thermodynamic
    !! database is inappropriate.  This subroutine is also used as a placeholder for additional checks
    !! in future development.  The following list gives a description of INFOThermo that is relevent
    !! to this subroutine.
    !!
    !   INFOThermo Value        Description
    !   ----------------        -----------
    !> \details
    !!          0          -    Successful exit;
    !!          9          -    A pure chemical species does not exist for at least one element.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nElements                 Number of elements in the system.
    ! nSpecies                  Number of species in the system.
    ! dAtomFractionSpecies      A double real matrix representing the atomic fraction of each element in
    !                            a species.
    ! iAtomFractionSpecies      A temporary integer matrix of dAtomFractionSpecies.
    ! iTempVec                  A temporary integer vector representing the maximum value of each column
    !                            of iAtomFractionSpecies.
    ! INFOThermo                An integer scalar identifying whether the program exits successfully or
    !                            if it encounters an error.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine CheckThermoData

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer                               :: i, j, INFO, nNonDummy
    real(8)                               :: rnorm
    integer,dimension(nSpecies,nElements) :: iAtomFractionSpecies
    integer,dimension(nElements)          :: iTempVec
    integer,dimension(nSpecies - nDummySpecies)          :: indx
    real(8),dimension(nSpecies - nDummySpecies)           :: dMolesElementTemp, x
    real(8), dimension(nElements,nSpecies - nDummySpecies):: dStoichSpeciesTemp
    real(8), dimension(nSpecies - nDummySpecies)  :: work


    ! A pure chemical species will have a value of 1 for iAtomFractionSpecies:
    iAtomFractionSpecies = INT(dAtomFractionSpecies)

    ! The iTempVec vector stores the number of pure species for each element:
    iTempVec = MAXVAL(iAtomFractionSpecies, DIM = 1)

    ! Check to see if at least one pure species exists in the database for each element in the system:
    do i = 1, nElements
        ! Skip if the system component corresponds to an electron:
        if (cElementName(i) == 'e-') cycle
        if (iTempVec(i) == 0) INFOThermo = 9
    end do


    ! Check if the input masses can be represented in terms of the available species
    nNonDummy = nSpecies - nDummySpecies
    dMolesElementTemp(1:nElements) = dMolesElement
    do i = 1, nNonDummy
        do j = 1, nElements
            dStoichSpeciesTemp(j,i) = dStoichSpecies(i,j)
        end do
    end do

    call nnls(dStoichSpeciesTemp, nElements, nNonDummy, dMolesElementTemp, x, rnorm, work, indx, INFO)

    if (rnorm > 1D-12) then
        INFOThermo = 99
    end if


    return

end subroutine CheckThermoData
