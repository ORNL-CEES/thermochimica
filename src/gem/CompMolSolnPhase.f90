
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompMolSolnPhase.f90
    !> \brief   Compute the number of moles of all stable solution phases.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      AddSolnPhase.f90
    !> \sa      SwapSolnPhase.f90
    !> \sa      SwapSolnForPureConPhase.f90
    !
    !
    ! References:
    ! ===========
    !
    ! For further information regarding this software, refer to the following material:
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
    !   10/25/2011      M.H.A. Piro         Clean up code: modules, simplify programming
    !   02/10/2012      M.H.A. Piro         Apply a minimum number of moles of a solution phases that is
    !                                       introduced to the system.
    !   04/01/2012      M.H.A. Piro         Constrain the maximum number of moles of a solution phase that
    !                                       is added to the system.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   03/04/2013      M.H.A. Piro         Fix bug when performing LLS when involving ionic phases.  The
    !                                        number of moles of an electron as a system component is zero,
    !                                        and the bug prevented divide by zero.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to estimate the number of moles of each solution phase in the
    !! system when a new solution phase is added.  This can be very important to consider because an initial
    !! estimate of zero moles would mean that the solution phase constituents would not contribute to the
    !! Jacobian matrix.  A more serious issue is if the system is initially estimated to be only comprised of
    !! pure stoichiometric phases and then a solution phase is added, then the upper left quadrant of the
    !! Jacobain matrix would be zero.
    !!
    !! The number of moles of all pure condensed phases are fixed and the number of moles of solution phases
    !! are computed using a Linear Least Squares (LLS) approach.  The linear minimization in this subroutine
    !! involves the relative errors of the mass balances instead of the residuals.  Normalizing each coefficient
    !! of matrix A and vector B by the number of moles of the corresponding element provides greater generality
    !! when the number of moles of different elements varies by many orders of magnitude.
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

subroutine CompMolSolnPhase

    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                               :: i, j, k, M, N, MN
    integer                               :: NRHS, LDA, LDB, INFO, LWORK
    real(8)                               :: dURF
    real(8),  dimension(nElements)        :: dTempVecE
    real(8),  dimension(:,:), allocatable :: A
    real(8),  dimension(:),   allocatable :: B
    real(8),  dimension(:),   allocatable :: WORK
    character                             :: TRANS


    ! Initialize variables:
    TRANS     = 'N'
    NRHS      = 1
    INFO      = 0
    M         = nElements
    N         = nSolnPhases
    LDA       = MAX(1,M)
    LDB       = MAX(1,M,N)
    MN        = MIN(M,N)
    LWORK     = MAX(1,MN + MAX(MN, NRHS))
    dURF      = 0.1D0
    dTempVecE = 0D0

    ! Allocate memory:
    allocate(A(M,N))
    allocate(B(M))
    allocate(WORK(LWORK))

    ! Initialize allocatable variables:
    A    = 0D0
    B    = 0D0
    WORK = 0D0

    ! Construct matrix A representing the "effective stoichiometry" of each solution phase:
    do j = 1, nSolnPhases
        do i = 1,nElements
            k = -iAssemblage(nElements - j + 1)
            A(i,j) = dEffStoichSolnPhase(k,i)
        end do
    end do

    ! Construct vector B, which represents the relative error of each mass balance (with pure condensed
    ! phase, but without solution phases):
    do j = 1,nElements
        B(j) = dMolesElement(j)
        do i = 1, nConPhases
            B(j) = B(j) - dMolesPhase(i) * dStoichSpecies(iAssemblage(i),j)
        end do
    end do

    ! Call the DGELS driver routine from LAPACK to compute the next estimated number of moles for each
    ! solution phase.  This routine solves the Linear Least Squares (LLS) problem using QR or LQ factorization
    ! in double precision.
    call DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )

    ! Store the number of moles of all phases from the previous iteration:
    dTempVecE = dMolesPhase

    ! Update the estimated number of moles of each solution phase:
    do i = 1,nSolnPhases
        j = nElements - i + 1
        k = -iAssemblage(j)
        if (dTempVecE(j) == 0D0) then
            dMolesPhase(j) = dURF * B(i) * dSumMolFractionSoln(k)
            dMolesPhase(j) = DMAX1(dMolesPhase(j),dTolerance(9))
            dMolesPhase(j) = DMIN1(dMolesPhase(j),dTolerance(10))
        end if

        if (dMolesPhase(j) < 0D0) then
            dMolesPhase(j) = DABS(dURF * dMolesPhase(j))
        end if
    end do

    i = 0

    ! Deallocate memory:
    deallocate(A,B,WORK, STAT = i)
    if (i /= 0) then
        INFOThermo = 22
        return
    end if

    return

end subroutine CompMolSolnPhase
