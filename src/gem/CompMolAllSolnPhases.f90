
!-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompMolAllSolnPhases.f90
    !> \brief   Compute the number of moles of all stable pure condensed and solution phases.
    !> \author  M.H.A. Piro
    !> \date    June 20, 2012
    !> \sa      AddSolnPhase.f90
    !> \sa      SwapSolnPhase.f90
    !> \sa      SwapSolnForPureConPhase.f90
    !> \sa      CompMolSolnPhase.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   06/20/2012      M.H.A. Piro         Original code
    !   01/31/2013      M.H.A. Piro         Previously, every row of A was divided by the corresponding value of b
    !                                        (which set b = 1D0) in an attempt to minimize numerical error.  The
    !                                        problem with this approach is that the number of moles of a system
    !                                        component is zero when it is an electron (i.e., a SUBL phase).
    !                                        The A matrix is no longer divided by b and the b vector is not 1.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to estimate the number of moles of each solution and pure
    !! condensed phase in the system.  The stoichiometry of each solution phase and (obviously) pure condensed
    !! phase is fixed and this subroutine determines a particular combination of molar quantities that minimizes
    !! the mass balance residuals using a Linear Least Squares (LLS) technique.
    !!
    !! This subroutine addresses a particular issue that occasionally arises when a pure condensed phase is
    !! driven out of the system (when it should), but it also drives a solution phase to almost be removed when
    !! it shouldn't.  At this point, the system may have sufficiently diverged that it may be impossible for
    !! the system to recover.  The element potentials may be close to the equilibrium values, but the molar
    !! quantities of solution species are way off.  This subroutine therefore fixes the chemical potentials
    !! and determines the best combination of molar quantities of the phases.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! iAssemblage               Integer vector containing the indices of phases at equilibrium
    ! dStoichSpecies            The number of atoms of a particular element for a particular species.
    ! dMolesElement             Total number of moles of an element, where the coefficient corresponds to the
    !                           atomic number (e.g., dMolesElement(92) refers to uranium).
    ! dMolesPhase               Number of moles of a phase at equilibrium
    ! dEffStoichSolnPhase       The effective stoichiometry of a solution phase
    ! dSumMolFractionSoln       Sum of all mole fractions within a solution phase.
    !
!-------------------------------------------------------------------------------------------------------------

subroutine CompMolAllSolnPhases

    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                               :: i, j, k, l, M, N, MN
    integer                               :: NRHS, LDA, LDB, INFO, LWORK
    real(8)                               :: dURF, dTemp
    real(8),  dimension(:,:), allocatable :: A
    real(8),  dimension(:),   allocatable :: B
    real(8),  dimension(:),   allocatable :: WORK
    character                             :: TRANS


    ! Return if there aren't any solution phases:
    if (nSolnPhases == 0) return

    ! Initialize variables:
    TRANS     = 'N'
    INFO      = 0
    NRHS      = 1
    M         = nElements - nChargedConstraints
    N         = nSolnPhases + nConPhases
    LDA       = MAX(1,M)
    LDB       = MAX(1,M,N)
    MN        = MIN(M,N)
    LWORK     = MAX(1,MN + MAX(MN, NRHS))
    dURF      = 0.5D0

    ! Allocate memory:
    allocate(A(LDA,N),B(LDB),WORK(LWORK))

    ! Initialize allocatable variables:
    A        = 0D0
    B        = 0D0
    WORK     = 0D0

    ! Construct matrix A representing the "effective stoichiometry" of each solution phase:
    do j = 1, nSolnPhases
        ! Absolute solution phase index:
        k = -iAssemblage(nElements - j + 1)
        ! Compute the stoichiometry of this solution phase:
        call CompStoichSolnPhase(k)
        do i = 1, M
           A(i,j) = dEffStoichSolnPhase(k,i)
        end do
    end do

    ! Construct matrix A representing the stoichiometry of each pure condensed phase:
    do j = nSolnPhases + 1, nSolnPhases + nConPhases
        k = iAssemblage(j - nSolnPhases)
        do i = 1, M
            A(i,j) = dStoichSpecies(k,i)
        end do
    end do

    ! Apply the right hand side vector:
    do i = 1, M
        B(i) = dMolesElement(i)
    end do

    ! This routine solves the Linear Least Squares (LLS) problem using QR or LQ factorization
    ! in double precision:
    call DGELS(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)

    dMolesPhase = 0D0

    ! Update the estimated number of moles of each solution phase:
    LOOP_SolnPhases: do k = 1, nSolnPhases
        j = nElements - k + 1
        dMolesPhase(j) = B(k)
        if (.NOT. lConverged) dMolesPhase(j) = DMAX1(dMolesPhase(j),dTolerance(9))
        if (.NOT. lConverged) dMolesPhase(j) = DMIN1(dMolesPhase(j),10D0*dTolerance(10))
        l = -iAssemblage(j)
        dTemp = 0D0
        do i = nSpeciesPhase(l-1) + 1, nSpeciesPhase(l)
            dMolesSpecies(i) = dMolFraction(i) * dMolesPhase(j)
            dTemp = dTemp + dMolFraction(i)
        end do
    end do LOOP_SolnPhases

    do k = nSolnPhases + 1, nSolnPhases + nConPhases
        j = k - nSolnPhases
        dMolesPhase(j) = B(k)
        if (.NOT. lConverged) dMolesPhase(j) = DMAX1(dMolesPhase(j),dTolerance(9))
        if (.NOT. lConverged) dMolesPhase(j) = DMIN1(dMolesPhase(j),10D0*dTolerance(10))
    end do

    ! Deallocate memory:
    i = 0
    deallocate(A, B, WORK, STAT = i)

    ! Check for any issues and report an error if appropriate:
    if (i /= 0) INFOThermo = 26

    return

end subroutine CompMolAllSolnPhases
