
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GetNewAssemblage.f90
    !> \brief   Determine the next phase assemblage to be considered in Leveling.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      ShuffleAssemblage.f90
    !> \sa      LevelingSolver.f90
    !> \sa      PostLevelingSolver.f90
    !
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
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code
    !   07/31/2011      M.H.A. Piro         Clean up code: remove unnecessary variables, update variable names
    !   10/21/2011      M.H.A. Piro         Clean up code: modules, simplify iteration history check.
    !   01/21/2013      M.H.A. Piro         Improved the iteration history check and created the
    !                                        CheckLevelingIterHistory subroutine.  Also, a previous check for
    !                                        the mass balance constraints computed the sum of coefficients of
    !                                        A along the 2nd dimension.  This was changed from an integer
    !                                        vector to a double vector and the absolute quantity of each
    !                                        coefficient is taken to appropriately handle anions.
    !   02/14/2013      M.H.A. Piro         The ShuffleDummySpecies subroutine was created to shuffle dummy
    !                                        species to the front of the iAssemblage vector.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to provide a new estimated phase assemblage to be tested in
    !! the Leveling subroutine.  The phase with the most negative relative Gibbs energy will be introduced
    !! into the previous estimated phase assemblage.  In order to avoid violating the Gibbs Phase Rule,
    !! this new phase will replace a phase from the previous assemblage.  The conditions for considering
    !! a new phase assemblage are:
    !!
    !! <ol>
    !! <li> The number of moles of each phase must be non-negative and real, </li>
    !! <li> The phase assemblage has not been previously tested (this is only
    !!      perfomed after x iterations), </li>
    !! <li> The numerical adjustments to the Gibbs Plane are real.  This last check does
    !!      not add any additional expense because this would have to be computed anyways. </li>
    !! </ol>
    !!
    !! There is an important, yet subtle, verification performed in the third (3) condition.  The Phase Rule
    !! dictates that the maximum number of phases that can coexist in an isobaric-isothermal closed system
    !! cannot exceed the number of system components (elements).  An additional condition to the Phase Rule,
    !! which is normally implied in thermodynamics texts but not explicity stated, is that only one pure
    !! separate phase of one component can exist at equilibrium.  For example, U(BCC) and U(FCC) cannot coexist.
    !! If two pure separate phases of the same component (same X) are included in the phase assemblage, then the
    !! Gibbs Plane is not uniquely defined (A is not a unique matrix) resulting in non-real adjustments to the
    !! Gibbs Plane.  By ensuring that the adjustments applied to the Gibbs Plane are real not only avoids obvious
    !! numerical problems, but also guarentees that the Phase Rule is explicitly satisfied.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO                      An integer scalar used internally by LAPACK that returns 0 for a successful exit
    !                            and a non-zero value when an error has occurred.
    ! iNewPhase                 Integer index representing the new phase that is to be introduced to the system.
    ! iterHistoryLevel          Integer matrix storing the history of the phase assemblages that were tested.
    ! iAssemblage               Integer vector containing the indices of phases in the assemblage
    ! iSpeciesAtoms             Stoichiometry coefficient of a particular species (e.g., UO2: 1 U, 2 O).
    ! nElements                 An integer scalar representing the number of elements in the system.
    ! dMolesPhase               Number of moles of a phase.  These are directly mapped to phases in iAssemblage
    ! dMolesElement             Number of moles of a particular element in the system.
    ! dChemicalPotential        The Relative Gibbs Energy is defined as the difference between the chemical
    !                            potential and standard Gibbs energy of the pure species.
    ! dLevel                    The adjustment applied to the element potentials.
    ! dAtomFractionSpecies      Atomic fraction of a particular element (e.g., UO2: 0.333 U, 0.667 O).
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine GetNewAssemblage(iter)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer                                :: i, j, k, INFO, iNewPhase, iter
    integer,dimension(nElements)           :: iAssemblageLast, IPIV
    real(8),dimension(nElements)           :: dTemp
    real(8),dimension(nElements,nElements) :: A
    logical                                :: lPhasePass


    ! Initialize variables:
    dLevel = 0D0
    IPIV   = 0
    A      = 0D0
    dLevel = 0D0

    ! Determine index of the species with the most negative relative Gibbs energy:
    iNewPhase = MAXVAL(MINLOC(dChemicalPotential))

    ! Shuffle the phase assemblage to make the best candidate phase to be tested first:
    call ShuffleAssemblage(iNewPhase,j)

    ! Shuffle the dummy species:
    call ShuffleDummySpecies

    ! Store the indices of the estimated phase assemblage at each iteration.
    do i = 1, nElements - nChargedConstraints
        iterHistoryLevel(i,iter) = iAssemblage(i)
    end do

    ! Store indices of previous phase assemblage:
    iAssemblageLast = iAssemblage

    ! Loop through all "phases" in the current phase assemblage to determine which one should be
    ! substituted for the new "phase":
    LOOP_NewAssemblage: do k = 1, nElements - nChargedConstraints

        ! The "phase" with the most negative dChemicalPotential replaces one of the phases from the previous
        ! assemblage:
        iAssemblage    = iAssemblageLast
        iAssemblage(k) = iNewPhase

        ! 1) FIRST CHECK: The number of moles of each "phase" are non-negative and real:
        do j = 1,nElements - nChargedConstraints
            do i = 1,nElements - nChargedConstraints
                A(i,j) = dStoichSpecies(iAssemblage(j),i)
            end do
            dMolesPhase(j) = dMolesElement(j)
        end do

        ! Compute the sum of coefficients along each row.  If the sum of coefficients on a row is
        ! equal to zero, this will result in a NaN.  This will quickly reject this phase assemblage
        ! and is intended to reduce computational expense.
        dTemp = SUM(DABS(A), DIM = 2)

        ! Cycle if the A matrix is underdetermined:
        if (MINVAL(dTemp) == 0D0) cycle LOOP_NewAssemblage

        ! Reinitialize variables:
        INFO = 0
        IPIV = 0

        ! Call the linear equation solver to compute molar quantities of the phase assemblage:
        call DGESV( nElements, 1, A, nElements, IPIV, dMolesPhase, nElements, INFO )

        ! Cycle if the number of moles of any "phase" is negative or non-real:
        if ((INFO /= 0).OR.(MINVAL(dMolesPhase) < dTolerance(3))) cycle LOOP_NewAssemblage

        ! 2) SECOND CHECK: Verify that this phase assemblage has not been previously considered.
        call CheckLevelingIterHistory(iter,lPhasePass)

        ! Cycle to the next phase if this assemblage has previsouly been considered:
        if (.NOT.(lPhasePass)) cycle LOOP_NewAssemblage

        ! 3) THIRD CHECK: Verify that the phase assemblage produces real co-ordinates for the Gibbs Plane.
        ! Note: this does not add much computational expense becaue dLevel needs to be computed anyways.
        do j = 1,nElements - nChargedConstraints
            do i = 1,nElements - nChargedConstraints
                A(i,j) = dAtomFractionSpecies(iAssemblage(i),j)
            end do
            dLevel(j) = dChemicalPotential(iAssemblage(j))
        end do

        ! Reinitialize variables:
        INFO = 0
        IPIV = 0

        ! Call linear equation solver to solve the adjustments applied to the Gibbs Plane:
        call DGESV( nElements, 1, A, nElements, IPIV, dLevel, nElements, INFO )

        ! Verify that the co-ordinates are real, otherwise reject the phase assemblage:
        if (INFO == 0) then
            exit LOOP_NewAssemblage
        else
            cycle LOOP_NewAssemblage
        end if

    end do LOOP_NewAssemblage

    ! Write an error code in case if Leveling was unsuccessful:
    if (INFO /= 0) INFOThermo = 10

    return

end subroutine GetNewAssemblage

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to check the iteration history for a
    ! leveling iteration.  The iteration history is scanned backwards and then
    ! every phase in the assemblage in question is compared to the full
    ! assemblage at a particular iteration.  Thus, the order of phases in the
    ! integer array does not affect the calculation.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iter        The iteration index in Leveling.
    !> \param[out]  lPhasePass  A logical scalar indicating whether the phase
    !                            assemblaged has passed (TRUE) or not (FALSE).
    !
    !---------------------------------------------------------------------------


subroutine CheckLevelingIterHistory(iter,lPhasePass)

    USE ModuleThermo

    implicit none

    integer:: i, j, k, iter
    logical:: lPhasePass


    ! Initialize variables:
    lPhasePass = .TRUE.

    ! Loop through iteration history:
    LOOP_Iter: do i = iter, 1, -1

        ! Loop through phases in the assemblage in question:
        LOOP_Outer: do j = 1, nElements - nChargedConstraints

            ! Loop through phases in the assemblage at a particular iteration in the history:
            LOOP_Inner: do k = 1, nElements - nChargedConstraints

                ! If the following statement is true, then the coefficients are consistent and
                ! move on to the next phase in the assemblage in question.
                if (iterHistoryLevel(k,i) == iAssemblage(j)) cycle LOOP_Outer

            end do LOOP_Inner

            ! If this statement is reached, then there is a discrepancy between the two phase assemblages.
            ! Cycle to the next assemblage in the iteration history:
            cycle LOOP_Iter

        end do LOOP_Outer

        ! The phase assemblage in question has previously been considered. The phase assemblage has failed.
        lPhasePass = .FALSE.
        exit LOOP_Iter

    end do LOOP_Iter

    return

end subroutine CheckLevelingIterHistory

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to shuffle dummy species to the front of
    ! the integer vector iAssemblage.  The motivation for doing this is that
    ! dummy species will be considered first to be exchanged for a new phase.
    ! Note that a unique solution is not guaranteed in Leveling in situations
    ! where zero moles is assigned to one of the dummy phases (degrees of freedom
    ! is non-zero).  Therefore, it is possible for G_sys to be at a minimum
    ! but the phase assemblage contains a dummy phase.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! iPhase        An integer vector dimensioned to the total number of
    !                species in the system.  Each coefficient corresponds to
    !                the phase type corresponding to this species.
    !                iPhase < 0: A dummy species.
    !                iPhase = 0: A pure condensed phase.
    !                iPhase > 0: A solution species, where iPhase corresponds
    !                            to the absolute solution phase index.
    !
    !---------------------------------------------------------------------------


subroutine ShuffleDummySpecies

    USE ModuleThermo, ONLY: iAssemblage, nElements, iPhase, nChargedConstraints

    implicit none

    integer                     :: i, j, k, nDummy
    integer,dimension(nElements):: iTempVec


    ! Initialize variables:
    nDummy   = 0
    iTempVec = 0

    ! Count the number of dummy species:
    do i = 1, nElements - nChargedConstraints
        if (iPhase(iAssemblage(i)) < 0) nDummy = nDummy + 1
    end do

    ! Only proceed if a dummy species is currently in the phase assemblage:
    IF_Proceed: if (nDummy > 0) then

        ! Shuffle vector:
        j = 0
        k = nDummy
        do i = 1, nElements - nChargedConstraints
            if (iPhase(iAssemblage(i)) < 0) then
                j = j + 1
                iTempVec(j) = iAssemblage(i)
            else
                k = k + 1
                iTempVec(k) = iAssemblage(i)
            end if
        end do

        ! Update the iAssemblage vector:
        iAssemblage = iTempVec

    end if IF_Proceed

    return

end subroutine ShuffleDummySpecies


    !---------------------------------------------------------------------------
    !                       END - GetNewAssemblage.f90
    !---------------------------------------------------------------------------
