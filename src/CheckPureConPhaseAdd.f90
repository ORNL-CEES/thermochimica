
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckPureConPhaseAdd.f90
    !> \brief   Check whether a pure condensed phase should be added to the system.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      CheckPhaseAssemblage.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/05/2012      M.H.A. Piro         Original code.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   09/23/2012      M.H.A. Piro         The driving force of pure condensed phases is now computed in
    !                                        the CheckPhaseAssemblage subroutine.
    !   09/24/2012      M.H.A. Piro         ShuffleAssemblage now shuffles the entire assemblage and indicates
    !                                        whether a pure condensed or solution phase should be swapped first.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check whether a pure condensed phase should be added to
    !!  the system.  The condition for a pure condensed phase to be added is when the driving force is less
    !!  than a specified tolerance (e.g., ~ -1D-5).  The driving force is defined as the difference between
    !!  the standard molar Gibbs energy of this phase and the corresponding value computed by the element
    !!  potentials.  Refer to the following literature for a more thorough discussion:
    !!
    !!      H.L. Lukas, S.G. Fries and B. Sundman, "Computational Thermodynamics: The Calphad Method,"
    !!      Cambridge University Press, New York, 2007.
    !!
    !! A pure condensed phase can be added to the system when the current number of phases is less than the
    !! number of elements in the system; however, it must replace an existing phase in the event that the number
    !! of phases is equal to the number of elements to prevent contradicting Gibbs' Phase Rule.  If a particular
    !! phase cannot be added directly to the system, then it will try to replace an existing phase.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iMaxDrivingForce    An integer scalar representing the index of the pure condensed phase
    !!                                   with the most negative driving force.
    !> \param[in]   dMaxDrivingForce    A double real scalar representing the driving force.
    ! nConPhases            The number of pure condensed phases in the assemblage.
    ! nSolnPhases           The number of solution phases in the assemblage.
    ! nElements             The number of elements in the system.
    ! iPhase                An integer vector representing the phase type (0: pure condensed phase; > 0: solution
    !                        phase index; -1: dummy species).
    ! iPhaseChange          An integer scalar representing the index of a phase that will be added to, or
    !                        withdrawn from, the current estimated phase assemblage.
    ! iPhaseTypeOut         An integer scalar representing the phase type to be shuffled by the ShuffleAssemblage
    !                        subroutine.  This subroutine shuffles the assemblage to determine the best candidate
    !                        phase to be withdrawn from the assemblage.
    !                        iPhaseTypeOut = 0 for a pure condensed phase,
    !                        iPhaseTypeOut = 1 for a solution phase.
    ! lPhasePass            A logical variable indicating whether the new estimated phase assemblage passed
    !                        (.TRUE.) or failed (.FALSE.).
    ! lSwapLater            A logical variable indicating whether a phase should be swapped when lPhasPass
    !                        indicates that the phase cannot be added directly.
    ! dGEMFunctionNorm      A double real scalar representing the norm of the functional vector in the PGESolver.
    ! dTolerance            A double real vector representing numerical tolerances defined in InitThermo.f90.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckPureConPhaseAdd(iMaxDrivingForce, dMaxDrivingForce)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: iPhaseChange, iMaxDrivingForce
    real(8) :: dMaxDrivingForce
    logical :: lSwapLater, lSwapCheck, lPhasePass


    ! Initialize variables:
    lSwapLater   = .FALSE.
    lSwapCheck   = .FALSE.
    lPhasePass   = .FALSE.

    ! Do not even try adding a pure condensed phase if any phase has a negative molar quantity:
    if (MINVAL(dMolesPhase) < 0D0) return

    ! Do not try adding a phase unless if at least 5 iterations have passed:
    if (iterGlobal - iterlast < 5) return

    ! Determine whether a pure condensed phase can be swapped:
    if ((dGEMFunctionNorm <= dTolerance(11)).OR.(iterGlobal - iterlast >= 50)) lSwapCheck = .TRUE.

    ! Check if a pure condensed phase should be added:
    IF_CheckDrivingForce: if ((iMaxDrivingForce /= 0).AND.(dMaxDrivingForce < dTolerance(4))) then

        ! The pure condensed phase that is to be added has the maximum driving force:
        iPhaseChange = iMaxDrivingForce

        ! Only proceed if the phase is a pure condensed phase and that the Gibbs Phase Rule will not be violated:
        !if ((iPhase(iPhaseChange) == 0).AND.(nConPhases + nSolnPhases < nElements - nChargedConstraints)) then
        if ((iPhase(iPhaseChange) == 0).AND.(nConPhases + nSolnPhases < nElements )) then

            ! Try adding this pure condensed phase to the phase assemblage:
            call AddPureConPhase(iPhaseChange,lSwapLater,lPhasePass)

            ! The pure condensed phase could not be added the system directly, but it might be possible to
            ! swap it for another phase that is currently in the phase assmeblage:
            if ((lSwapLater).AND.((lSwapCheck).OR.(iterGlobal - iterLast >= 50))) then

                ! Check if a particular combination of phases should be swapped:
                call CheckPureConPhaseSwap(iPhaseChange,lSwapLater,lPhasePass)

            end if

        elseif ((iPhase(iPhaseChange) == 0).AND.(nConPhases + nSolnPhases == nElements ).AND. &
       ! elseif ((iPhase(iPhaseChange) == 0).AND.(nConPhases + nSolnPhases == nElements - nChargedConstraints).AND. &
            ((lSwapCheck))) then

            ! A phase cannot be added to the system when the number of phases already equals the number of
            ! elements as a consequence of the Phase Rule.  This phase must be swapped for another phase
            ! already expected to be stable.

            ! Check if a particular combination of phases should be swapped:
            call CheckPureConPhaseSwap(iPhaseChange,lSwapLater,lPhasePass)

        end if

    end if IF_CheckDrivingForce

    return

end subroutine CheckPureConPhaseAdd


!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to check whether a pure condensed phase
    ! should swap for another pure condensed phase or a solution phase.
    !
    !---------------------------------------------------------------------------


subroutine CheckPureConPhaseSwap(iPhaseChange,lSwapLater,lPhasePass)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::   iPhaseTypeOut, iPhaseChange
    logical::   lSwapLater, lPhasePass


    ! Shuffle the phase assemblage so that the best candidate to be swapped is first:
    call ShuffleAssemblage(iPhaseChange,iPhaseTypeOut)

    ! Check if a pure condensed phase or a solution phase should be swapped first:
    if (iPhaseTypeOut == 0) then
        ! A pure condensed phase should be swapped first.

        ! Try swapping this phase for another pure condensed phase in the assemblage:
        call SwapPureConPhase(iPhaseChange,lSwapLater,lPhasePass)

        ! Try swapping a solution phase for a pure condensed phase:
        if (iterGlobal /= iterLast) call SwapPureConForSolnPhase(iPhaseChange,lPhasePass)

    elseif (iPhaseTypeOut == 1) then
        ! A solution phase should be swapped first.

        ! Try swapping a solution phase for a pure condensed phase:
        call SwapPureConForSolnPhase(iPhaseChange,lPhasePass)

        ! Try swapping this phase for another pure condensed phase in the assemblage:
        if (iterGlobal /= iterLast) call SwapPureConPhase(iPhaseChange,lSwapLater,lPhasePass)

    end if

end subroutine CheckPureConPhaseSwap


    !---------------------------------------------------------------------------
    !                      END - CheckPureConPhaseAdd.f90
    !---------------------------------------------------------------------------
