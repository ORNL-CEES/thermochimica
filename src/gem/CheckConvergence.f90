
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckConvergence.f90
    !> \brief   Check convergence in the non-linear solver.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      GEMSolver.f90
    !> \todo    Consider removing check for residuals of chemical potential terms...this may be redundant.
    !
    !
    ! References:
    ! ===========
    !
    !        M.H.A. Piro, "Computation of Thermodynamic Equilibria Pertinent to Nuclear Materials
    !        in Multi-Physics Codes," PhD Dissertation, Royal Military College of Canada, 2011.
    !
    !        M.H.A. Piro, T.M. Besmann, S. Simunovic, B.J. Lewis and W.T. Thompson, "Numerical
    !        Verification of Equilibrium Thermodynamic Computations in Nuclear Fuel Performance
    !        Codes," Journal of Nuclear Materials, 414 (2011) 399-407.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code
    !   07/31/2011      M.H.A. Piro         Clean up code: remove unnecessary variables, update variable names
    !   09/12/2011      M.H.A. Piro         Added feature: check if dSumMolFractionSoln or dRelError is a NAN
    !   10/25/2011      M.H.A. Piro         Clean up code: modules, simplify programming
    !   10/27/2011      M.H.A. Piro         Modified check for Gibbs' Criteria: exclude dummy species
    !   04/05/2012      M.H.A. Piro         Clean up code: simplify relative error calculation of mass balance
    !                                        equations by using the dFunction vector.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm.
    !   08/22/2012      M.H.A. Piro         Add check for miscibility gap.
    !   01/31/2013      M.H.A. Piro         Apply appropriate check for residual/relative errors of the mass
    !                                        balance equations when the system component is an electron (i.e.,
    !                                        dMolesElement = 0D0).
    !   02/04/2013      M.H.A. Piro         Add a check to ensure that the Phase Rule is satisfied (this was
    !                                        necessarily satisfied in the previous version because charged
    !                                        phases were not considered.
    !   05/12/2014      M.H.A. Piro         Change the 9th check for a global minimum so that all metastable
    !                                        solution phases are checked for a global minimum, not just phases
    !                                        with a known miscibility gap.
    !   08/21/2015      M.H.A. Piro         Change tolerance for residuals of chemical potentials of individual
    !                                        species to only apply when x > 1.66D-24.  Practically, who cares
    !                                        if the residual of the chemical potential is large for a species
    !                                        whose mole fraction is 1D-50?!?  1.66D-24 is the fraction of a
    !                                        single atom in one mole.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check if convergence has been achieved.  The conditions for
    !! thermodynamic equilibrium are (refer to above references for more details):
    !! <ol>
    !! <li> None of the phases in the estimated phase assemblage are "dummy species", which were introduced
    !!      by the ChemSage data-file, </li>
    !! <li> The number of moles of all species and phases are positive and non-zero. </li>
    !! <li> Gibbs' Phase Rule has been satisfied. </li>
    !! <li> The standard Gibbs energy of a pure condensed phase is not below the Gibbs Plane.
    !!      If so, this phase should be added to the assemblage. </li>
    !! <li> The residuals of the chemical potential terms are below a specified tolerance. </li>
    !! <li> The sum of site fractions on each sublattice must equal unity within tolerance. </li>
    !! <li> The driving force for every solution phase must be non-negative. </li>
    !! <li> The relative errors of the mass balance equations are within tolerance. </li>
    !! <li> The integral Gibbs energy of the system must be at a true global minimum and not
    !!       a local minima. </li>
    !! </ol>
    !!
    !! Each criterion listed above is sequentially tested and the system is considered converged
    !! when all are satisfied.  Control is returned to the PGESolver subroutine when any of the
    !! criterions has not been satisfied.  Note that the order of testing is done in a fashion
    !! that progressively increases computational expense.  For example, testing the mass balance
    !! constraints is the most computationally expensive task, which is why it is performed last.
    !!
    !! Note that the Gibbs Phase Rule is necessarily satisfied because iAssemblage is dimensioned by
    !! nElements.  Therefore, it is impossible for nPhases > nElements and the Gibbs Phase Rule is
    !! implicitly satisfied at all times.
    !!
    !
    ! Pertinent variables:
    ! ====================
    !
    ! lConverged                Logical variable: true if convedRelGibbsEnergyd, false if not converged.
    ! dStoichSpecies            Number of atoms of a particular element per formula mass of a particular
    !                           species (double matrix).
    ! iAssemblage               Integer vector storing the indices of phases contributing to the equilibrium
    !                           phase assemblage.
    ! iElement                  Integer vector storing the atomic numbers of elements in the system.
    ! dRelError                 Relative error of mass balance equations.
    ! dChemicalPotential        A double real vector representing the difference between the chemical potential
    !                           (defined by the element potentials) and the standard molar Gibbs energy of the
    !                           pure species divided by the total number of atoms per formula mass.
    ! dMolesPhase               The number of moles of a phase.
    ! dMolesElement             Total number of moles of an element.
    ! dEffStoichSolnPhase       The "effective" stoichiometry of a solution phase.
    ! dSumMolFractionSoln       Sum of mole fractions in a solution phase.
    ! dTolerance                Acceptable numerical tolerance (defined in Init subroutine).
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckConvergence

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, l, c, iMaxDrivingForce
    real(8) :: dResidual, dMaxDrivingForce
    logical :: lCompEverything, lPhaseChange


    ! Initialize variables:
    dResidual       = 0D0
    lConverged      = .FALSE.
    lCompEverything = .TRUE.
    lPhaseChange    = .FALSE.

    ! Reset lSolnPhases to ensure it is correct
    lSolnPhases     = .FALSE.
    do i = 1, nSolnPhases
        if (iAssemblage(nElements+1-i) < 0) lSolnPhases(-iAssemblage(nElements+1-i)) = .TRUE.
    end do

    ! Test if the largest relative change in species mole fraction is large
    ! This is a self-consistency check
    if (lDebugMode) print *, "Test self-consistency ", dMaxSpeciesChange
    if (dMaxSpeciesChange > LOG(2D0)) then
        return
    end if

    ! TEST #1: Check if any of the phases in the assemblage are "dummy" phases:
    ! -------------------------------------------------------------------------
    if (lDebugMode) print *, "Test 1"
    LOOP_TEST1: do i = 1, nElements
        if (iAssemblage(i) > 0) then
            if (iPhase(iAssemblage(i)) < 0) return
        end if
    end do LOOP_TEST1

    ! TEST #2: Check to make sure that the number of moles of all phases are non-negative.
    ! ------------------------------------------------------------------------------------
    dMolesPhase(nConPhases + 1 : nElements - nSolnPhases) = 0D0
    if (lDebugMode) print *, "Test 2 ", MINVAL(dMolesPhase)
    if (MINVAL(dMolesPhase) < 0D0) return

    ! TEST #3: Check that the Phase Rule has been satisfied:
    ! -----------------------------------------------------------------------------------
    if (lDebugMode) print *, "Test 3 "
    if (nSolnPhases + nConPhases > nElements - nChargedConstraints) then
        if (lDebugMode) print *, "nSolnPhases: ", nSolnPhases, " + nConPhases: ", nConPhases, &
                    " > nElements: ", nElements, " - nChargedConstraints: ", nChargedConstraints
        call CorrectPhaseRule(lPhaseChange)
    end if

    if (lPhaseChange) return

    ! TEST #8: Check the relative errors of the mass balance equations:
    ! -----------------------------------------------------------------
    if (lDebugMode) print *, "Test 8 "
    LOOP_TEST8: do j = 1, nElements
        dResidual = -dMolesElement(j)
        ! Loop through solution phases:
        do l = 1, nSolnPhases
            k = -iAssemblage(nElements - l + 1)                 ! Absolute solution phase index.
            do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)     ! Absolute solution species index.
                dResidual = dResidual + dMolesSpecies(i) * dStoichSpecies(i,j) / DFLOAT(iParticlesPerMole(i))
            end do
        end do

        ! Loop through pure condensed phases:
        do i = 1, nConPhases
            k = iAssemblage(i)                                  ! Absolute pure condensed phase index.
            if (k > 0) dResidual = dResidual + dMolesPhase(i) * dStoichSpecies(k,j)
        end do

        ! Compute residual or relative error term:
        if (dMolesElement(j) > 0D0) then
            dResidual = dResidual / dMolesElement(j)
        end if
        if (DABS(dResidual) > dTolerance(1)) then
            if (lDebugMode) print *, cElementName(j), dMolesElement(j), dResidual
            return
        end if
    end do LOOP_TEST8

    ! CONVERGENCE TEST SHORTCUT
    ! -----------------------------------------------------------------------------------
    if (lDebugMode) print *, "Test Shortcuts start ", dGEMFunctionNorm
    ! Now that crucial tests have been done, can check for convergence shortcut
    ! If the functional norm is less than a specified tolerance and the system hasn't changed,
    ! call it a day:
    if ((dGEMFunctionNorm < dTolerance(12)).AND.(iterGlobal - iterLast > 100).AND.(iterGlobal > 4000))  then
        if (lDebugMode) print *, "Took Shortcut 1 after ", iterGlobal
        lConverged = .TRUE.
        return
    else if ((dGEMFunctionNorm < 1D-2).AND.          (iterGlobal - iterLast > 5000)) then
        if (lDebugMode) print *, "Took Shortcut 2 after ", iterGlobal
        lConverged = .TRUE.
        return
    else if ((dGEMFunctionNorm < 1D-5).AND.(iterGlobal > 4000)) then
        if (lDebugMode) print *, "Took Shortcut 3 after ", iterGlobal
        lConverged = .TRUE.
        return
    end if
    if (lDebugMode) print *, "Test Shortcuts end "
    ! Return if the functional norm is too large.  In other words, it's not worth the flops checking.
    if (dGEMFunctionNorm > dTolerance(1)) return

    ! TEST #4: Check whether a pure condensed phase should be added to the phase assemblage:
    ! --------------------------------------------------------------------------------------

    ! Compute the driving force for all pure condensed phases:
    call CompDrivingForce(iMaxDrivingForce, dMaxDrivingForce)
    if (iMaxDrivingForce > 0) then
        if (lDebugMode) print *, "Test 4 ", dMaxDrivingForce, cSpeciesName(iMaxDrivingForce)
    else
        if (lDebugMode) print *, "Test 4 ", dMaxDrivingForce
    end if
    if (dMaxDrivingForce < dTolerance(4)) return

    ! TEST #5: Check the residuals of chemical potential terms:
    ! ---------------------------------------------------------
    if (lDebugMode) print *, "Test 5 "
    LOOP_TEST5: do j = 1, nSolnPhases
        k = -iAssemblage(nElements - j + 1)

        ! Compute the residual terms:
        LOOP_GRESIDUAL: do i = nSpeciesPhase(k-1)+1, nSpeciesPhase(k)

            ! Compute the chemical potential term from the element potentials:
            dResidual = 0D0
            do l = 1, nElements
                dResidual = dResidual + dElementPotential(l) * dStoichSpecies(i,l)
            end do

            ! Normalize this quantity by the number of particles per mole:
            dResidual = dResidual / DFLOAT(iParticlesPerMole(i))

            ! Compute absolute quantity of the residual term:
            dResidual = DABS(dResidual - dChemicalPotential(i))

            ! Note: the mole fraction is checked if it is greater than 1.66D-24 because of
            ! potential numerical errors.  This particular number is chosen because it is
            ! the fraction of a single atom in 1 mole.
            if ((dResidual > dTolerance(5)).AND.(dMolFraction(i) > 1.66D-24)) return

        end do LOOP_GRESIDUAL

    end do LOOP_TEST5

    ! TEST #6: Check the residuals of site fractions on each sublattice of a CEF phase:
    ! ---------------------------------------------------------------------------------
    if (lDebugMode) print *, "Test 6 "
    ! Loop through solution phases that are stable:
    LOOP_TEST6: do j = 1, nSolnPhases

        ! Absolute solution phase index:
        k = -iAssemblage(nElements - j + 1)

        ! Check if this phase is represented by sublattices:
        if ((cSolnPhaseType(k) == 'SUBL').OR.(cSolnPhaseType(k) == 'SUBLM')) then

            ! Relative charged phased index:
            l = iPhaseSublattice(k)

            ! Loop through sublattices:
            do i = 1, nSublatticePhase(l)
                dResidual = -1D0

                ! Loop through constituents on this sublattice:
                do c = 1, nConstituentSublattice(l,i)
                    dResidual = dResidual + dSiteFraction(l,i,c)
                end do

                dResidual = DABS(dResidual)

                ! Return if the residual is greater than the tolerance:
                if (dResidual > dTolerance(2)) return

            end do
        end if
    end do LOOP_TEST6

    ! TEST #7: Check if a solution phase should be added to the phase assemblage:
    ! ---------------------------------------------------------------------------
    if (lDebugMode) print *, "Test 7 "
    LOOP_TEST7: do i = 1, nSolnPhasesSys

        ! Skip this phase if it is already predicted to be stable:
        if (lSolnPhases(i) .EQV. .TRUE.) cycle LOOP_TEST7

        ! Skip this phase if it is not the first "phase" of a phase containing a miscibility gap
        ! (this will be handled in test #8):
        if (lMiscibility(i) .EQV. .TRUE.) cycle LOOP_TEST7

        ! Check if the driving force of this solution phase is less than the tolerance:
        if (dDrivingForceSoln(i) < dTolerance(4)) then
            if (lDebugMode) print *, 'Negative driving force for ', cSolnPhaseName(i)
            return
        end if

    end do LOOP_TEST7

    ! TEST #9: Check for a miscibility gap:
    ! -------------------------------------
    if (lDebugMode) print *, "Test 9 "
    ! Loop through all solution phases in the system to check for a miscibiltiy gap:
    LOOP_TEST9: do i = 1, nSolnPhasesSys

        ! Skip if this phase is already part of the estimated phase assemblage:
        if (lSolnPhases(i) .EQV. .TRUE.) cycle LOOP_TEST9

        ! Check if this phase may contain a miscibility gap:
        if (lMiscibility(i) .EQV. .TRUE.) then

            ! If this phase contains a miscibility gap but none of the corresponding phases are stable, then cycle:
            LOOP_DoubleCheckMG: do j = 1, nSolnPhases
                k = -iAssemblage(nElements - j + 1)

                if (k == i) cycle LOOP_DoubleCheckMG

                if ((cSolnPhaseName(k) == cSolnPhaseName(i)).AND.(lSolnPhases(k) .EQV. .FALSE.)) cycle LOOP_TEST9

            end do LOOP_DoubleCheckMG

            ! Check for a miscibility gap:
            call CheckMiscibilityGap(i,lPhaseChange)

            ! Return if this phase should be added to the system:
            if (lPhaseChange .EQV. .TRUE.) return

            if (dDrivingForceSoln(i) < dTolerance(4)) return

        else

            ! Check solution phase type:
            if (cSolnPhaseType(i) == 'IDMX') then
                ! This is an ideal mixing phase:
                call CompMolFraction(i)
            else
                ! This is a non-ideal solution phase:
                call CheckMiscibilityGap(i,lPhaseChange)
            end if

            ! If the driving force is negative, return:
            if (dDrivingForceSoln(i) < dTolerance(4)) return

        end if

    end do LOOP_TEST9

    ! If all of the above criterions have been satisfied, then the system has converged.
    if (INFOThermo == 0) lConverged = .TRUE.

    return

end subroutine CheckConvergence
