
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file        ModuleGEMSolver.f90
    !> \brief       Fortran module for input/output of the non-linear solver.
    !> \author      M.H.A. Piro
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !
    !> \param iterLast          The last global iteration that the phase assemblage was adjusted.
    !> \param iterHistory       An integer matrix representing all of the indices of phases that contribute
    !!                           to the equilibrium phase assemblage at each stage in the iteration history.
    !> \param dGEMFunctionNorm  A double real scalar representing the norm of the functional vector in the
    !!                          GEMSolver.
    !> \param dGEMFunctionNormLast A double real scalar representing the norm of the functional vector in the
    !!                              GEMSolver from the last iteration.
    !> \param dSumMolFractionSoln  A double real vector representing the sum of mole fractions in each solution
    !!                              phase.
    !> \param dUpdateVar        A double real vector representing the direction vector that updates the function
    !                            vector.
    !> \param dPartialExcessGibbs  A double real vector representing the partial molar excess Gibbs energy of
    !!                              mixing of each species in the system.
    !> \param dEffStoichSolnPhase   A double real matrix representing the effective stoichiometry of each
    !!                               solution phase.
    !> \param lDebugMode        A logical variable used for debugging purposes.  When it is TRUE, a number
    !!                           of print statements are applied.
    !> \param lRevertSystem     A logical variable identifying whether the system should be reverted (TRUE)
    !!                           or not (FALSE).
    !> \param lConverged        A logical variable identifying whether the system has converged (TRUE) or not
    !!                           (FALSE).
    !> \param lSolnPhases       A logical vector indicating whether a particular solution phase is currently
    !!                           assumed to be stable (true) or not (false).
    !> \param lMiscibility      A logical vector indicating whether a particular solution phase has a
    !!                           miscibility gap (true) or not (false).
    !
    !
    ! CONSIDER REMOVING THE FOLLOWING VARIABLES:
    !
    ! iterLastCon           The last iteration that a pure condensed phase was either added to or removed from
    !                        the estimated phase assemblage.
    ! iterLastSoln          The last iteration that a pure condensed phase was either added to or removed from
    !                        the estimated phase assemblage.
    ! iConPhaseLast         The species index of the last pure condensed phase that was either added to, or
    !                        removed from, the estimated phase assemblage.
    ! iSolnPhaseLast        The species index of the last pure condensed phase that was either added to, or
    !                        removed from, the estimated phase assemblage.
    !
    !
    !-------------------------------------------------------------------------------------------------------------

module ModuleGEMSolver

    implicit none

    SAVE

    integer                              ::  iterLast,      iterStep, iterRevert, iterGlobal
    integer                              ::  iterLastCon,   iterLastSoln,         iterSwap,   iterLastMiscGapCheck
    integer                              ::  iConPhaseLast, iSolnPhaseLast,       iSolnSwap,  iPureConSwap
    integer,                 parameter   ::  iterGlobalMax = 3000
    integer, dimension(:,:), allocatable ::  iterHistory

    real(8)                              ::  dGEMFunctionNorm,    dGEMFunctionNormLast, dMaxSpeciesChange, dMinGibbs
    real(8), dimension(:),   allocatable ::  dSumMolFractionSoln, dMolesPhaseLast, dUpdateVar, dDrivingForceSoln
    real(8), dimension(:),   allocatable ::  dPartialExcessGibbs, dPartialExcessGibbsLast
    real(8), dimension(:,:), allocatable ::  dEffStoichSolnPhase

    logical                              ::  lDebugMode, lRevertSystem, lConverged
    logical, dimension(:),   allocatable ::  lSolnPhases, lMiscibility

end module ModuleGEMSolver
