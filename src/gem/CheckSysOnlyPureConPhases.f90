
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckSysOnlyPureConPhases.f90
    !> \brief   Check the system when only pure condensed phases are expected to be stable.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      CompMolFraction.f90
    !> \sa      CheckQKTOSolnPhase.f90
    !> \sa      CheckConvergence.f90
    !> \todo    Recompute the element potentials.
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
    !   02/10/2012      M.H.A. Piro         Original code
    !   02/14/2012      M.H.A. Piro         Perform refinement iteration
    !   02/15/2012      M.H.A. Piro         Call the CompFunction subroutine in case if the functional
    !                                        norm is incorrectly computed.
    !   09/29/2012      M.H.A. Piro         Remove the call to CheckQKTOSolnPhase and rely on the driving
    !                                        force computed by Subminimization.
    !   01/31/2013      M.H.A. Piro         Added correction when charged phases are considered in the system.
    !
    !
    ! Purpose:
    ! ========
    !
    !\details The purpose of this subroutine is to test whether a phase assemblage is appropriate when the
    !! system is comprised of only pure condensed phases.  The number of moles of each phase is computed
    !! twice (by using the residuals from the first call as the input to the second call to the linear equation
    !! solver) to minimize numerical errors since this is a direct calculation (as opposed to an iterative one
    !! when solution phases are included, where the errors are progressively reduced).
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! A                         Stoichiometry matrix
    ! B                         Constraint vector (before call to LAPACK); moles of pure condensed phases
    !                           (after call to LAPACK)
    ! dStoichSpecies            Stoichiometry coefficient (used here for pure condensed phases)
    ! iAssemblage               Integer vector representing the indices of phases included in the current
    !                           estimated phase assemblage.
    ! dMolesElement             Double real vector representing the total number of moles of each element.
    ! lPhasePass                Logical variable indicating whether the phase assemblage has passed or failed.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckSysOnlyPureConPhases

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModulePhaseConstraints

    implicit none

    integer                                 :: i, j, INFO
    integer, dimension(nElements)           :: IPIV
    real(8), dimension(nElements,nElements) :: A, AA
    real(8), dimension(nElements)           :: B, BB
    real(8)                                 :: dTemp


    ! Initialize variables:
    IPIV = 0
    INFO = 0
    A    = 0D0
    AA   = 0D0
    B    = 0D0
    BB   = 0D0

    ! Do not adjust the phase assemblage when constraints are active:
    if (nPhaseConstraints > 0) return

    ! If this comes straight from levelling, phases may have 0 moles
    ! These should be removed
    if (iterGlobal == 0) then
        do i = 1, nConPhases
            if (dMolesPhase(i) < dTolerance(7)) then
                iAssemblage(i) = iAssemblage(nConPhases)
                iAssemblage(nConPhases)   = 0
                dTemp                     = dMolesPhase(i)
                dMolesPhase(i) = dMolesPhase(nConPhases)
                dMolesPhase(nConPhases)   = 0D0
                nConPhases                = nConPhases - 1
            end if
        end do
    end if

    ! Construct the stoichiometry matrix and the constraint vector:
    do j = 1, nElements
        do i = 1,nConPhases
            A(j,i) = dStoichSpecies(iAssemblage(i),j)
        end do
        B(j) = dMolesElement(j)
    end do

    ! Check if there are any charged phases in the database:
    if (nCountSublattice > 0) then
        ! Pure condensed phases are necessarily neutrally charged; thus, this must
        ! be removed from the mass balance calculation.
        do j = nElements, nElements - nChargedConstraints + 1, -1
            A(j,j) = 1D0
        end do
    end if

    if (nConPhases > 1) then
        ! Store the stoichiometry matrix:
        AA = A

        ! Call the linear equation solver to solve the number of moles of each phase:
        call DGESV( nElements, 1, A, nElements, IPIV, B, nElements, INFO )

        ! Compute the residual vector:
        BB = dMolesElement - MATMUL(AA,B)

        ! Perform an additional calculation to refine the number of moles of each phase
        ! (this reduces numerical errors):
        call DGESV( nElements, 1, AA, nElements, IPIV, BB, nElements, INFO )

        ! Update the number of moles of the pure condensed phases:
        if (INFO == 0) dMolesPhase = BB + B

        ! Compute the mole fractions of all solution phase constituents:
        do i = 1, nSolnPhasesSys

            call CompMolFraction(i)

        end do

    end if

    ! Compute functional norm:
    call CompFunctionNorm

    ! Check Convergence:
    call CheckConvergence

    return

end subroutine CheckSysOnlyPureConPhases
