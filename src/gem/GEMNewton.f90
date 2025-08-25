
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GEMNewton.f90
    !> \brief   Compute the direction vector for the GEMSolver using Newton's method.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      GEMSolver.f90
    !> \sa      GEMLineSearch.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code (new GEM solver)
    !   05/25/2012      M.H.A. Piro         Check for a NAN immediately after call to DGESV.
    !   01/31/2013      M.H.A. Piro         Check if a charged phase is contained in the database, but is
    !                                        not represented by the current phase assemblage.
    !   03/04/2013      M.H.A. Piro         Fix bug in correction process when dealing with ionic phases
    !                                        the loop should count back from the number of constraints,
    !                                        not the number of charged phases).
    !   09/06/2021      M. Poschmann        Correct moles of species to be proportional to mole fraction times
    !                                        moles of respective phase before direction vector is computed.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \brief The purpose of this subroutine is to compute the direction vector for the Gibbs energy
    !! minimization (GEM) solver using Newton's method.  The Hessian matrix and its corresponding constraint
    !! vector are first constructed and then the direction vector representing the system parameters is solved
    !! with the DGESV driver routine from LAPACK.  The updated element potentials, adjustments to the number of
    !! moles of solution phases and the number of moles of pure condensed phases are applied in the
    !! GEMLineSearch.f90 subroutine.
    !!
    !! Thermochimica is capable of handling ionic phases, which have an additional charge neutrality
    !! constraint imposed for each ionic phase.  Thus, an electron is added as a system component for every
    !! charged phase in the system.  It may be possible that an ionic phase is not predicted to be stable at
    !! a particular iteration and, thus, there aren't any stable species in the system representing that electron.
    !! To prevent a numerical singularity in the Hessian matrix, a check is performed after the Hessian matrix
    !! has been constructed ensuring that the Hessian does not contain a zero row.  In the event that the
    !! Hessian matrix contains all zeroes in the jth row (and necessarily, the jth column), a unit value is
    !! assigned to A(j,j).  Since the total balance of an electron is necessarily zero (i.e., ensuring charge
    !! neutrality) and there aren't any species for this solution phase, the corresponding value on the b vector
    !! will also be zero.  This procedure effectively ignores the jth row while preventing a numerical
    !! singularity.
    !
    !
    ! References:
    ! ===========
    !
    !> \details For further information regarding this methodology, refer to the following material:
    !! <ul>
    !! <li>  W.B. White, S.M. Johnson, G.B. Dantzig, "Chemical Equilibrium in Complex Mixtures," Journal of
    !!        Chemical Physics, V. 28, N. 5, 1958.
    !!
    !! <li>  G. Eriksson, "Thermodynamic Studies of High Temperature Equilibria," Acta Chemica Scandinavica,
    !!        25, 1971.
    !!
    !! <li>  G. Eriksson, E. Rosen, "General Equations for the Calculation of Equilibria in Multiphase Systems,"
    !!        Chemica Scripta, 4, 1973.
    !! </ul>
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[out]  INFO        An integer scalar used by LAPACK indicating a successful exit or an error.
    !
    ! nVar                      An integer scalar representing the total number of unknowns/linear equations.
    ! nElements                 An integer scalar representing the total number of elements in the system.
    ! nSpeciesPhase             An integer vector representing the number of species in a particular solution
    !                            phase (accumulative indexing)
    ! dStoichSpecies            A double real matrix representing stoichiometry coefficients.
    ! dMolesSpecies             A double real vector representing the number of moles of each species.
    ! dMolesPhase               A double real vector representing the number of moles of each phase.
    ! dMolesElement             A double real vector representing the number of moles of each element.
    ! JacobianLong              A double real matrix representing part of the Jacobian matrix that involves the
    !                            stoichiometry coefficients of solution species.
    ! JacobianShort             A double real vector that incorporates the JacobianLong matrix along with the
    !                            updated number of moles of each solution species.
    ! A                         Hessian matrix
    ! B                         Constraint vector (before call to LAPACK); unknown vector (after call to LAPACK)
    ! dEffStoichSolnPhase       A double real matrix representing the effective stoichiometry of a solution phase.
    ! dUpdateVar                A double real vector represending the updated system variables.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine GEMNewton(INFO)

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleGEMSolver

    implicit none

    integer                              :: i, j, k, l, m, INFO, nVar, iTry, nMaxTry
    integer, dimension(nElements)        :: iErrCol
    real(8)                              :: dTemp

    ! Count phases:
    j = nConPhases
    nConPhases  = 0
    CountCon: do i = 1, j
        if (iAssemblage(i) > 0) then
            nConPhases  = nConPhases  + 1
        else
            exit CountCon
        end if
    end do CountCon

    j = nSolnPhases
    nSolnPhases = 0
    CountSoln: do i = nElements, nElements + 1 - j, -1
        if (iAssemblage(i) < 0) then
            nSolnPhases = nSolnPhases + 1
        else
            exit CountSoln
        end if
    end do CountSoln

    if ((nConPhases + nSolnPhases) <= 0) return

    ! Determine the number of unknowns/linear equations:
    nVar = nElements + nConPhases + nSolnPhases

    iErrCol = 0
    nMaxTry = nElements - (nConPhases + nSolnPhases)
    TryLoop: do iTry = 0, nMaxTry
        ! on retry we are going to use dummy phases
        if (iTry > 0) nVar = nElements * 2

        ! Allocate or resize workspace if necessary:
        if (.NOT. allocated(GEM_A) .OR. SIZE(GEM_A,1) < nVar) then
            if (allocated(GEM_A)) deallocate(GEM_A)
            allocate(GEM_A(nVar, nVar))
        end if
        if (.NOT. allocated(GEM_B) .OR. SIZE(GEM_B) < nVar) then
            if (allocated(GEM_B)) deallocate(GEM_B)
            allocate(GEM_B(nVar))
        end if
        if (.NOT. allocated(GEM_IPIV) .OR. SIZE(GEM_IPIV) < nVar) then
            if (allocated(GEM_IPIV)) deallocate(GEM_IPIV)
            allocate(GEM_IPIV(nVar))
        end if

        ! Initialize variables:
        GEM_IPIV(1:nVar)        = 0
        INFO                    = 0
        GEM_A(1:nVar,1:nVar)    = 0D0
        GEM_B(1:nVar)           = 0D0
        dUpdateVar              = 0D0
        dEffStoichSolnPhase     = 0D0

        do k = 1, nSolnPhases
            ! Absolute solution phase index:
            m = -iAssemblage(nElements - k + 1)
            ! Loop through species in phase:
            do l = nSpeciesPhase(m-1) + 1, nSpeciesPhase(m)
                dMolesSpecies(l) = dMolesPhase(nElements - k + 1) * dMolFraction(l)
                dMolesSpecies(l) = DMAX1(dMolesSpecies(l), dTolerance(8))
            end do
        end do

        ! Construct the Hessian matrix (elements):
        do j = 1, nElements
            do i = j, nElements
                do k = 1, nSolnPhases
                    ! Absolute solution phase index:
                    m = -iAssemblage(nElements - k + 1)
                    ! Loop through species in phase:
                    do l = nSpeciesPhase(m-1) + 1, nSpeciesPhase(m)
                        dTemp  = dStoichSpecies(l,i) * dStoichSpecies(l,j) * dMolesSpecies(l)
                        GEM_A(i,j) = GEM_A(i,j) + dTemp / (DFLOAT(iParticlesPerMole(l))**2)
                    end do
                end do
                ! Apply symmetry:
                GEM_A(j,i) = GEM_A(i,j)
            end do
        end do

        ! Compute the constraint vector (elements):
        do j = 1, nElements
            GEM_B(j) = dMolesElement(j)
            do l = 1, nSolnPhases
                k = -iAssemblage(nElements - l + 1)
                do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    dTemp = dStoichSpecies(i,j) * dMolesSpecies(i) * (dChemicalPotential(i) - 1D0)
                    GEM_B(j)  = GEM_B(j) + dTemp / DFLOAT(iParticlesPerMole(i))
                end do
            end do
        end do

        ! Construct the Hessian matrix and constraint vector (contribution from solution phases):
        do j = nElements + 1, nElements + nSolnPhases
            l = 2 * nElements - j + 1       ! Relative solution phase index (in iAssemblage vector).
            k = -iAssemblage(l)             ! Absolute solution phase index.

            ! Compute the stoichiometry of this phase:
            call CompStoichSolnPhase(k)

            do i = 1,nElements
                GEM_A(i,j) = dEffStoichSolnPhase(k,i) * dMolesPhase(l)
                GEM_A(j,i) = GEM_A(i,j)
            end do
            GEM_B(j) = dGibbsSolnPhase(k)
        end do

        ! Construct the Hessian matrix and constraint vector (contribution from pure condensed phases):
        do j = nElements + nSolnPhases + 1, nElements + nConPhases + nSolnPhases
            k = j - nElements - nSolnPhases
            do i = 1, nElements
                GEM_A(i,j) = dStoichSpecies(iAssemblage(k),i)
                GEM_A(j,i) = GEM_A(i,j)
            end do
            GEM_B(j) = dStdGibbsEnergy(iAssemblage(k))
        end do

        do k = 1, iTry
            i = iErrCol(k)
            j = nElements + nSolnPhases + nConPhases + k
            GEM_A(i,j) = 1D0
            GEM_A(j,i) = GEM_A(i,j)
            GEM_B(j) = 0D0
        end do

        ! Check if the Hessian is properly structured if the system contains any charged phases:
        if (nCountSublattice > 0) then
            ! Loop through elements
            LOOP_SUB: do j = nElements, nElements - nChargedConstraints + 1, -1
                dTemp = 0D0
                ! Loop through coefficients along column:
                do i = 1, nElements
                    dTemp = dTemp + DABS(GEM_A(i,j))
                    if (dTemp > 0D0) cycle LOOP_SUB
                end do
                ! The phase corresponding to this electron is not stable.
                GEM_A(j,j) = 1D0
            end do LOOP_SUB
        end if

        ! Call the linear equation solver:
        if ((nConPhases > 1) .OR. (nSolnPhases > 0)) then
            ! TODO: explore sparse or structured solvers when GEM_A contains many zeros
            call dgesv( nVar, 1, GEM_A, SIZE(GEM_A, 1), GEM_IPIV, GEM_B, SIZE(GEM_B, 1), INFO )
        else
            do i = 1, nElements
                GEM_B(i) = dElementPotential(i)
            end do
            GEM_B(nElements + 1) = dMolesPhase(1)
        end if

        do k = 1, iTry
            j = nElements + nSolnPhases + nConPhases + k
            GEM_B(j) = 0D0
        end do

        ! Check for a NAN:
        LOOP_CheckNan: do i = 1, nVar
            if (GEM_B(i) /= GEM_B(i)) then
                INFO = 1
                exit LOOP_CheckNan
            end if
        end do LOOP_CheckNan

        if (iTry < nMaxTry) then
            if ((INFO <= 0) .OR. (INFO > nElements)) then
                exit TryLoop
            else
                iErrCol(iTry+1) = INFO
                INFO = 0
            end if
        end if
    end do TryLoop

    ! Store the updated variables if LAPACK is successful:
    if (INFO == 0) then
        do j = 1, nVar
            dUpdateVar(j) = GEM_B(j)
        end do

        ! Reset:
        lRevertSystem = .FALSE.
    else
        ! The system failed.  Revert to a previous assemblage.
        lRevertSystem = .TRUE.
        dUpdateVar    = 0D0
    end if

    ! Workspace arrays retained between iterations

    return

end subroutine GEMNewton
