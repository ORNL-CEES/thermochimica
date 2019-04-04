
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    Subminimization.f90
    !> \brief   Determine whether a particular non-ideal solution phase should be added to the system
    !!           by performing a subminimization routine.
    !> \author  M.H.A. Piro
    !> \date    Aug. 21, 2012
    !> \sa      CheckSolnPhaseAdd.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   08/21/2012      M.H.A. Piro         Original code
    !   08/30/2012      M.H.A. Piro         Enforce the Wolfe conditions in the line search algorithm.  Also,
    !                                        a check is performed to ensure that the solution phase is not the
    !                                        same as the other "solution phase."
    !   08/31/2012      M.H.A. Piro         Enforce a minimum mole fraction to prevent the system from
    !                                        continuing to push the mole fraction of a particular constituent
    !                                        to zero.  The tolerance is defined as the quotient of machine
    !                                        precision and the relative error tolerance of the mass balance
    !                                        (i.e., 10**(-15) / 10**(-5)).
    !   09/03/2012      M.H.A. Piro         Revert the last change...do not test for a minimum mole fraction.
    !   09/19/2012      M.H.A. Piro         Replaced the DGESV solver with the ArrowSolver.  This solver
    !                                        exploits the fact that the Hessian is a symmetric arrow matrix.
    !   09/24/2012      M.H.A. Piro         Refine convergence criteria: if the driving force doesn't change
    !                                        by a predefined amount (e.g., 1%).
    !   10/22/2012      M.H.A. Piro         Set a minimum constraint for dMolFraction of 1D-10 when initializing
    !                                        Subminimization.  This prevents extremely small values from being
    !                                        considered.
    !   02/28/2013      M.H.A. Piro         The routine was upgraded to allow for ionic phases (i.e., charge
    !                                        neutrality constraints must be considered).
    !   03/05/2013      M.H.A. Piro         Fix bug in SubminNewton in constructing the Hessian matrix: the
    !                                        stoichiometry coefficients of ionic species should be multiplied
    !                                        by (-1).
    !   03/12/2013      M.H.A. Piro         Compute the functional norm (currently computed by the residual
    !                                        of the mass fractions and charge neutrality constraints) and
    !                                        only proceed if this is below tolerance.
    !   03/14/2013      M.H.A. Piro         Constrain the minimum mole fraction to an arbitrarily small value
    !                                        (e.g., 1D-100) to avoid an floating point underflow error.
    !   06/04/2013      M.H.A. Piro         Lax tolerance for dSubMinFunctionNorm.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to determine whether a particular non-ideal solution phase
    !! should be added to the system by performing a subminimization.  The criteria for adding any type of phase
    !! (regardless of whether it is a pure stoichiometric phase, ideal or non-ideal solution phase) is based on
    !! whether the driving force is positive (it should not be added) or negative (it should be added to the
    !! system).  The driving force is defined as the difference between the molar Gibbs energy of a particular
    !! phase and the corresponding value defined by the element potentials.  A graphical interpretation of the
    !! driving force is the difference between the molar Gibbs energy of a particular phase (represented by a
    !! line for a stoichiometric phase or a point on a curve for a solution phase) and the corresponding value
    !! on the Gibbs hyperplane.
    !!
    !! For further information regarding the term "driving force", refer to:
    !!
    !!      H.L. Lukas, S.G. Fries and B. Sundman, "Computational Thermodynamics: The Calphad Method,"
    !!      Cambridge University Press, New York (2007).
    !!
    !! The premise of the subminimization method is to determine a unique combination of mole fractions for a
    !! particular non-ideal solution phase that minimizes the driving force of that phase.  As an additional
    !! constraint, the sum of the mole fractions of this phase must equal unity.  In the subminimization approach,
    !! a new Lagrangian function is defined as:
    !!
    !! \f$ L_{\lambda} = \sum_{i=1}^{N_{\lambda}} x_{i(\lambda)}^n \left( \mu_{i(\lambda)}^n -
    !! \sum_{j=1}^E a_{i,j} \Gamma_j^m\right) - \pi_{\lambda}^n \left( \sum_{i=1}^{N_{\lambda}}
    !! x_{i(\lambda)}^n - 1 \right)  \f$
    !!
    !! which solves for \f$ N_{\lambda} + 1 \f$ variables corresponding to \f$ x_{i(\lambda)}^n \f$ and
    !! \f$ \pi_{\lambda}^n  \f$.  For more information regarding this approach, refer to the following literature:
    !!
    !!      C.E. Harvie, J.P. Greenberg and J.H. Weare, "A Chemical Equilibrium Algorithm for Highly Non-Ideal
    !!      Multiphase Systems: Free Energy Minimization," Geochimica et Cosmochimica Acta, V. 51 (1987)
    !!      1045-1057.
    !!
    !
    !  In a previous version of Thermochimica, a solution phase was added to the system based on a different
    !  approach.  The previous method would compute the mole fractions of solution phase constituents as a
    !  function of the element potentials and that phase would be added to the system if the sum of mole
    !  fractions of all constituents in that phase was greater than unity.   The biggest issue with doing
    !  this is that the mole fraction of a constituent in a non-ideal solution phase IS NOT a function of
    !  the element potentials due to the inclusion of the partial molar excess Gibbs energy of mixing term.
    !  The previous approach relied heavily on a convoluted technique of dampening excess terms due to the highly
    !  non-linear nature of the equations.  Essentially, the mole fraction is written in that aproach as an
    !  expontential function of mole fractions.
    !
    !  The Subminimization approach achieves the same net effect, but due to the simplicity of the approach, the
    !  numerics are much more robust and efficient.  One of the main issues with the previous approach is that
    !  the mole fractions may not be unique to the current estimated element potentials.  In the subminimization
    !  approach, the mole fractions are unique to the current estimated element potentials.  Furthermore, the
    !  computational expense associated with each global iteration is far less.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iSolnPhaseIndex  An integer scalar representing the absolute index of the solution phase that
    !!                                is being considered.
    !> \param[out]  lPhasePass       A logical scalar indicating whether the phase should be added (i.e., TRUE)
    !                                 to the system or not (i.e., FALSE).
    !
    ! lDuplicate                    A logical scalar indicating whether the mole fractions for this "phase" are
    !                                a duplicate of the mole fractions of the corresponding solution phase.
    !
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine Subminimization(iSolnPhaseIndex,lPhasePass)

    USE ModuleThermo
    USE ModuleSubMin
    uSE ModuleGEMSolver!, ONLY: lMiscibility, dDrivingForceSoln
    USE ModuleThermoIO

    implicit none

    integer :: iterSub,    iSolnPhaseIndex, iterSubMax
    logical :: lPhasePass, lDuplicate


    ! Initialize local variables:
    lPhasePass = .FALSE.
    lDuplicate = .FALSE.

    ! Initialize the subminimization subroutine:
    call SubMinInit(iSolnPhaseIndex,iterSubMax)

    ! Subminimization iteration loop:
    LOOP_IterSub: do iterSub = 1, iterSubMax

        ! Compute the Newton direction vector:
        call SubMinNewton(iSolnPhaseIndex)

        ! Compute an appropriate step length:
        call SubMinLineSearch(iSolnPhaseIndex)

        ! Compute the functional norm:
        call SubMinFunctionNorm(iSolnPhaseIndex)

        ! Check if the minimum mole fraction is below a certain tolerance:
        if (MINVAL(dMolFraction(iFirst:iLast)) < dMinMoleFraction) exit LOOP_IterSub

        ! Only check for convergence if the functional norm is below tolerance:
        if (dSubMinFunctionNorm < dTolerance(1)) then

            ! Check convergence:
            if ((dDrivingForce <= dTolerance(4)).AND.(lMiscibility(iSolnPhaseIndex) .EQV. .FALSE.)) lSubMinConverged = .TRUE.

        end if

        ! Exit if the subminimization has converged:
        if ((lSubMinConverged .EQV. .TRUE.).OR.(INFOThermo /= 0)) exit LOOP_IterSub

        ! Check if the solution phases represening the miscibility gap duplicate one another:
        if (lMiscibility(iSolnPhaseIndex) .EQV. .TRUE.) call SubMinCheckDuplicate(lDuplicate)

        ! Exit if the subminimization has converged:
        if (lDuplicate .EQV. .TRUE.) exit LOOP_IterSub

    end do LOOP_IterSub

    ! If the composition of phases representing a miscibility gap duplicate one another (i.e., they have virtually
    ! the same composition), then set the driving force to zero to prevent this phase from being added to the system.
    if (lDuplicate .EQV. .TRUE.) dDrivingForce = 0D0

    ! If the driving force is less than a specified tolerance and the system has converged,
    ! add this solution phase to the system:
    if ((dDrivingForce < dTolerance(4)).AND.(lSubMinConverged .EQV. .TRUE.)) lPhasePass = .TRUE.

    ! If the mole fraction and charge neutrality constraints were not satisfied, then set the driving force to zero:
    if (dSubMinFunctionNorm > 10D0*dTolerance(1)) dDrivingForce = 0D0

    ! Update the solution driving force:
    dDrivingForceSoln(iSolnPhaseIndex) = dDrivingForce

    ! Deallocate allocatable arrays:
    deallocate(dChemicalPotentialStar, dRHS)

    deallocate(dHessian, iHessian)

    return

end subroutine Subminimization


!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to initialize the subminimiation routine.
    ! Variables are initialized, allocatable arrays are allocated and the
    ! chemical potential terms of the solution phase constituents defined by
    ! the element potentials are computed (i.e., dChemicalPotentialStar).
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! iFirst            Absolute index of first species in phase.
    ! iLast             Absolute index of last species in phase.
    ! nVar              Total number of species in phase.
    ! iterSubMax        An integer scalar representing the maximum number
    !                    of iterations allowed in subminimization.
    ! lSubMinConverged  A logical scalar indicating convergence (true).
    ! dDrivingForceLast A double real scalar representing the driving force
    !                    from the last subminimization iteration.
    !
    !---------------------------------------------------------------------------

subroutine SubMinInit(iSolnPhaseIndex,iterSubMax)

    USE ModuleThermo
    USE ModuleSubMin
    USE ModuleGEMSolver, ONLY: lMiscibility, dDrivingForceSoln

    implicit none

    integer::  i, j, k, iSolnPhaseIndex, iterSubMax


    ! Initialize variables:
    iFirst            = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    iLast             = nSpeciesPhase(iSolnPhaseIndex)
    nVar              = iLast - iFirst + 1
    dDrivingForceLast = 10D0
    lSubMinConverged  = .FALSE.

    dDrivingForceSoln(iSolnPhaseIndex) = 0D0

    ! Check if allocatable arrays are already allocated.  Deallocate if necessary:
    if (allocated(dChemicalPotentialStar)) deallocate(dChemicalPotentialStar)
    if (allocated(dRHS))                   deallocate(dRHS)
    if (allocated(dHessian))               deallocate(dHessian)
    if (allocated(iHessian))               deallocate(iHessian)

    ! Allocate allocatable arrays:
    allocate(dChemicalPotentialStar(nVar))

    ! Determine prefactor for allocating arrays (depends on whether the phase is ionic):
    if (iPhaseElectronID(iSolnPhaseIndex) == 0) then
        i = 1
    else
        i = 2
    end if

    ! Allocate allocatable arrays:
    allocate(dRHS(nVar+i), iHessian(nVar+i), dHessian(nVar+i,nVar+i))

    ! Initialize variables:
    dRHS                   = 0D0
    dChemicalPotentialStar = 0D0

    ! Set a default value for iterSubMax if it is not specified:
    iterSubMax = 100

    ! Loop through all constituents in this solution phase:
    do k = 1, nVar

        ! Absolute species index:
        i = nSpeciesPhase(iSolnPhaseIndex-1) + k

        ! Compute the chemical potentials of all constituents defined by the element potentials:
        dChemicalPotentialStar(k) = 0D0
        do j = 1, nElements
            dChemicalPotentialStar(k) = dChemicalPotentialStar(k) + dElementPotential(j) &
                * dStoichSpecies(i,j)
        end do
        dChemicalPotentialStar(k) = dChemicalPotentialStar(k) / DFLOAT(iParticlesPerMole(i))

        ! Initialize the mole fractions:
        dMolFraction(i) = DMAX1(dMolFraction(i), 1D-15)

    end do

    ! Compute the chemical potentials of solution phase constituents:
    call SubMinChemicalPotential(iSolnPhaseIndex)

    ! Compute the driving force of this solution phase:
    call SubMinDrivingForce

    ! If this phase contains a miscibility gap, determine the absolute index of the corresponding
    ! solution phase with the miscibility gap:
    if (lMiscibility(iSolnPhaseIndex) .EQV. .TRUE.) then

        ! Check the name of the solution phases:
        if (cSolnPhaseName(iSolnPhaseIndex) == cSolnPhaseName(iSolnPhaseIndex-1)) then
            ! The last solution phase in the indexing scheme corresponds to the miscibility gap:
            iSolnPhaseIndexOther = iSolnPhaseIndex - 1
        else
            ! Loop through all solution phases in the data-file system:
            LOOP_SolnSys: do i = 1, nSolnPhasesSys

                ! Cycle if this is same phase:
                if (i == iSolnPhaseIndex) cycle LOOP_SolnSys

                ! Check if these phases are the same:
                if (cSolnPhaseName(iSolnPhaseIndex) == cSolnPhaseName(i)) then
                    iSolnPhaseIndexOther = i
                    exit LOOP_SolnSys
                end if
            end do LOOP_SolnSys
        end if
    end if

end subroutine SubMinInit


!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to compute the chemical potentials of
    ! only the constituents of this solution phase.
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dChemicalPotential    A double real vector representing the chemical
    !                        potential of every species in the system.
    ! dPartialExcessGibbs   A double real vector representing the partial
    !                        molar excess Gibbs energy of every species in
    !                        the system.
    !
    !---------------------------------------------------------------------------


subroutine SubMinChemicalPotential(iSolnPhaseIndex)

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin

    implicit none

    integer::  iSolnPhaseIndex


    ! Initialize variables:
    dChemicalPotential(iFirst:iLast)  = 0D0
    dPartialExcessGibbs(iFirst:iLast) = 0D0

    ! Compute excess terms based on solution phase type:
    call CompExcessGibbsEnergy(iSolnPhaseIndex)

end subroutine SubMinChemicalPotential

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to compute the driving force of a
    ! particular solution phase within the Subminimization routine.
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dDrivingForce     A double real scalar representing the driving force of
    !                    a solution phase.  The driving force represents the
    !                    difference between the molar Gibbs energy of that
    !                    particular phase and the corresponding value calculated
    !                    from the chemical potentials of the system components.
    !                    When negative, this term indicates that a phase should
    !                    be added to the system.  When positive, the phase
    !                    should not be added to hte system.
    !
    !---------------------------------------------------------------------------

subroutine SubMinDrivingForce

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin

    implicit none

    integer::  i, j


    ! Initialize variables:
    dDrivingForce = 0D0

    ! Loop through all species in this phase:
    do j = 1, nVar

        ! Absolute species index:
        i = iFirst + j - 1

        ! Update the driving force:
        dDrivingForce = dDrivingForce + dMolfraction(i) * (dChemicalPotential(i) - dChemicalPotentialStar(j))

    end do

    ! Normalize the driving force:
    dDrivingForce = dDrivingForce / DFLOAT(nVar)

end subroutine SubMinDrivingForce

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to compute the direction vector using
    ! an exact Newton method.
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO          An integer scalar used by LAPACK to identify a successful
    !                exit (INFO = 0) or an error (INFO /= 0).
    ! nEqn          The number of equations in the Hessian matrix.
    ! dRHS          A double real vector representing the right hand side (RHS)
    !                of the Hessian matrix and returns the direction vector.
    !
    !---------------------------------------------------------------------------

subroutine SubMinNewton(iSolnPhaseIndex)

    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin

    implicit none

    integer :: i, j, k, iSolnPhaseIndex, INFO, nEqn, m


    ! Initialize variables:
    nEqn     = nVar + 1
    INFO     = 0
    iHessian = 0
    dRHS     = -1D0
    dHessian = 0D0
    m        = 1

    ! Construct diagonal and part of the arrow head:
    do j = 1, nVar
        i                 = nSpeciesPhase(iSolnPhaseIndex-1) + j
        dHessian(j,j)     = 1D0 / dMolfraction(i)
        dHessian(nEqn,j)  = -1D0
        dHessian(j,nEqn)  = -1D0
        dRHS(j)           = dDrivingForce - (dChemicalPotential(i) + 1 - dChemicalPotentialStar(j))
        dRHS(nEqn)        = dRHS(nEqn) + dMolfraction(i)
    end do

    ! Apply an additional row/column if the phase is ionic:
    if (iPhaseElectronID(iSolnPhaseIndex) /= 0) then
        nEqn = nEqn + 1
        k    = iPhaseElectronID(iSolnPhaseIndex)
        m    = 2

        dRHS(nEqn) = 0D0
        do j = 1, nVar
            i                 = nSpeciesPhase(iSolnPhaseIndex-1) + j
            dHessian(nEqn,j)  = -dStoichSpecies(i,k)
            dHessian(j,nEqn)  = dHessian(nEqn,j)
            dRHS(nEqn)        = dRHS(nEqn) + dStoichSpecies(i,k) * dMolfraction(i)
        end do
    end if

    ! Solve the arrow matrix:
    call ArrowSolver(m, nEqn, dHessian, dRHS, INFO)

    ! If the ArrowSolver failed, then try LAPACK:
    if (INFO /= 0) then
        INFO     = 0
        nEqn     = nVar + 1
        dHessian = 0D0
        dRHS     = -1D0

        ! Construct diagonal and part of the arrow head:
        do j = 1, nVar
            i                 = nSpeciesPhase(iSolnPhaseIndex-1) + j
            dHessian(j,j)     = 1D0 / dMolfraction(i)
            dHessian(nEqn,j)  = -1D0
            dHessian(j,nEqn)  = -1D0
            dRHS(j)           = dDrivingForce - (dChemicalPotential(i) + 1 - dChemicalPotentialStar(j))
            dRHS(nEqn)        = dRHS(nEqn) + dMolfraction(i)
        end do

        ! Apply an additional row/column if the phase is ionic:
        if (iPhaseElectronID(iSolnPhaseIndex) /= 0) then
            nEqn = nEqn + 1
            k    = iPhaseElectronID(iSolnPhaseIndex)
            m    = 2

            dRHS(nEqn) = 0D0
            do j = 1, nVar
                i                 = nSpeciesPhase(iSolnPhaseIndex-1) + j
                dHessian(nEqn,j)  = -dStoichSpecies(i,k)
                dHessian(j,nEqn)  = dHessian(nEqn,j)
                dRHS(nEqn)        = dRHS(nEqn) + dStoichSpecies(i,k) * dMolfraction(i)
            end do
        end if

        ! Call the linear equation solver:
        call DGESV( nEqn, 1, dHessian, nEqn, iHessian, dRHS, nEqn, INFO )

    end if

    ! Check for an error from ArrowSolver:
    if (INFO /= 0) then
        ! Return an error and reset dRHS:
        INFOThermo = 28
        dRHS       = 0D0
    end if

end subroutine SubMinNewton


!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to perform a line search on the
    ! direction vector computed by the SubMinNewton subroutine.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dStepLength   A double real scalar representing the steplength that is
    !                applied to updating the mole fractions.
    ! dTemp         A double real scalar representing a temporary variable.
    ! dMaxChange    A double real scalar representing the maximum change in the
    !                mole fraction of a solution phase constituent.
    !
    !---------------------------------------------------------------------------

subroutine SubMinLineSearch(iSolnPhaseIndex)

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin

    implicit none

    integer :: i, j, k, iSolnPhaseIndex
    real(8) :: dStepLength, dTemp, dMaxChange


    ! Initialize variables:
    dStepLength = 1D0
    dMaxChange  = 0.2D0

    ! Initialize steplength (determine a steplength that prevents negative mole fractions
    ! and constrains maximum change):
    do j = 1, nVar

        ! Absolute species index:
        i = iFirst + j - 1

        ! Check if the mole fraction of this constituent is driven to be negative:
        if (dMolfraction(i) + dRHS(j) <= 0D0) then

            ! Determine a step length that reduces the mole fraction by a factor of 100:
            dTemp = -0.99D0 * dMolfraction(i) / dRHS(j)

            ! Update the step length:
            dStepLength = DMIN1(dStepLength, dTemp)

        end if

        ! Compute a step length that constrains the maximum change by dMaxChange:
        dTemp=dStepLength
        if(dRHS(j) /= 0D0)then
           dTemp       = DABS(dMaxChange / dRHS(j))
        end if
        dStepLength = DMIN1(dStepLength, dTemp)

    end do

    ! Update the mole fractions:
    dTemp = 0D0
    do j = 1, nVar

        ! Absolute species index:
        i = iFirst + j - 1

        ! Apply step length:
        dMolfraction(i) = dMolfraction(i) + dStepLength * dRHS(j)

        ! Store maximum change to the mole fraction:
        dTemp = DMAX1(dTemp, DABS(dRHS(j)))

    end do

    ! Iterate to satisfy Wolfe conditions:
    LOOP_WOLFE: do k = 1, 5

        ! Exit if the minimum mole fraction of this phase is below a specified tolerance:
        if (MINVAL(dMolFraction(iFirst:iLast)) < dMinMoleFraction) exit LOOP_WOLFE

        ! Compute the chemical potentials of solution phase constituents:
        call SubMinChemicalPotential(iSolnPhaseIndex)

        ! Compute the driving force of this solution phase:
        call SubMinDrivingForce

        ! Check for a converging/diverging solution:
        if (dDrivingForce < dDrivingForceLast) then

            ! The driving force is less than the last driving force:
            exit LOOP_WOLFE

        else
            ! Divergence has been detected.  Dampen the sub-system:
            dStepLength = dStepLength * 0.5D0
            dTemp       = 0D0

            ! Adjust the mole fractions of species:
            do j = 1, nVar

                ! Absolute species index:
                i = iFirst + j - 1

                ! Apply step length:
                dMolfraction(i) = dMolfraction(i) - dStepLength * dRHS(j)

                ! Constrain the minimum mole fraction to an arbitrarily small value:
                dMolfraction(i) = DMAX1(dMolfraction(i), 1D-100)

                dTemp = DMAX1(dTemp, dStepLength * dRHS(j))

            end do

        end if

    end do LOOP_WOLFE

    ! Check convergence (maximum change to any mole fraction):
    if (dTemp <= dSubMinTolerance) lSubMinConverged = .TRUE.

    ! Store the driving force from the last iteration:
    dDrivingForceLast = dDrivingForce

end subroutine SubMinLineSearch



!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to compute the functional norm to
    ! test for convergence.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO          An integer scalar used by LAPACK to identify a successful
    !                exit (INFO = 0) or an error (INFO /= 0).
    !
    !---------------------------------------------------------------------------

subroutine SubMinFunctionNorm(iSolnPhaseIndex)

    USE ModuleThermo!, ONLY: dMolfraction, dStoichSpecies, iPhaseElectronID
    USE ModuleSubMin

    implicit none

    integer :: i, j
    integer :: iSolnPhaseIndex


    ! Initialize variables:
    dSubMinFunctionNorm = -1D0

    ! Compute residual of mole fractions:
    do i = iFirst, iLast
        dSubMinFunctionNorm = dSubMinFunctionNorm + dMolfraction(i)
    end do

    ! Compute residual of charge neutrality constraints:
    if (iPhaseElectronID(iSolnPhaseIndex) /= 0) then

        j = iPhaseElectronID(iSolnPhaseIndex)

        do i = iFirst, iLast
            dSubMinFunctionNorm = dSubMinFunctionNorm + dMolfraction(i) * dStoichSpecies(i,j)
        end do
    end if

    ! Compute absolute quantity:
    dSubMinFunctionNorm = DABS(dSubMinFunctionNorm)

end subroutine SubMinFunctionNorm



!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to check whether this yields a duplicate
    ! set of mole fractions as the "other solution phase" that corresponds to
    ! the miscibility gap.
    !
    !---------------------------------------------------------------------------

subroutine SubMinCheckDuplicate(lDuplicate)

    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleSubMin

    implicit none

    integer :: i, j, k, iFirstOther, iLastOther
    real(8) :: dTemp
    logical :: lDuplicate


    ! Initialize variables:
    iFirstOther = nSpeciesPhase(iSolnPhaseIndexOther-1) + 1
    iLastOther  = nSpeciesPhase(iSolnPhaseIndexOther)
    dTemp       = 0D0
    lDuplicate  = .FALSE.

    ! Double check to make sure that they have the same number of constituents:
    if (iLastOther - iFirstOther + 1 /= nVar) then
        INFOThermo = 29
        return
    end if

    ! Compute the Euclidean norm between the two mole fraction vectors:
    do k = 1, nVar

        ! Absolute species index:
        i = iFirst      + k - 1    ! Absolute index for first species in iSolnPhaseIndex
        j = iFirstOther + k - 1    ! Absolute index for first species in iSolnPhaseIndexOther

        dTemp = dTemp + (dMolfraction(i) - dMolfraction(j))**2

    end do

    ! Compute the normalized Euclidean norm between the mole fraction vectors between the two "phases":
    dTemp = ( dTemp**(0.5) ) / DFLOAT(nVar)

    ! Check if the normalized Euclidean norm is less than a specified tolerance:
    if (dTemp < dTolEuclideanNorm) lDuplicate = .TRUE.

end subroutine SubMinCheckDuplicate


    !---------------------------------------------------------------------------
    !                       END - Subminimization.f90
    !---------------------------------------------------------------------------
