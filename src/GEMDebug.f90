
subroutine GEMDebug(iDebug)


    !-------------------------------------------------------------------------------------------------------------
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
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
    !   04/05/2012      M.H.A. Piro         Original code
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to debug the PGESolver subroutine and all subroutines associated with it.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nConPhases            The number of pure condensed phases in the assemblage
    ! nSolnPhases           The number of solution phases in the assemblage
    ! nSolnPhasesSys        The number of solution phases in the system
    ! iAssemblage           Integer vector containing the indices of phases in the assemblage
    !                        (1:nConphases represent pure condensed phases and (nElements-nSolnPhases:nSolnPhases)
    !                        represent solution phases.
    ! INFOThermo            An integer scalar identifying whether the program exits successfully or if
    !                        it encounters an error.
    ! lConverged            A logical variable indicating whether the code has convered (.TRUE.) or not (.FALSE.).
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer::   iDebug, i, j, k

    if (iDebug == 1) then
        ! PGESolver.f90

        print *, '=================================='
        print *, 'ITERATION', iterGlobal, iterlast, lRevertSystem
        print *, '=================================='
        print *

    elseif (iDebug == 2) then
        ! PGENewton.f90


    elseif (iDebug == 3) then
        ! CompFunction.f90
        print *, 'Function Norm: '
        print *

    elseif (iDebug == 4) then
        ! CheckPhaseAssemblage.f90
        print *, 'CheckPhaseAssemblage: continue?'
        print *
    elseif (iDebug == 5) then
        ! CheckPhaseAssemblage.f90
        print *, 'test 1'
        print *
    elseif (iDebug == 6) then
        ! CheckPhaseAssemblage.f90
        print *, 'test 2'
        print *
    elseif (iDebug == 7) then
        ! CheckPhaseAssemblage.f90
        print *, 'test 3'
        print *
    elseif (iDebug == 8) then
        ! CheckPhaseAssemblage.f90
        print *, 'test 4'
        print *

    elseif (iDebug == 9) then
        ! GEMSolver.f90


        print *
        print *, '# of pure con phases = ', nConphases
        do i = 1, nConphases
            print *, cSpeciesName(iAssemblage(i)), iAssemblage(i), dMolesPhase(i)
        end do
        print *

        print *, '# of soln phases = ', nSolnPhases
        do i = 1, nSolnPhases
            j = nElements - i + 1
            k = -iAssemblage(j)
            print *, cSolnPhaseName(k), k, lMiscibility(k), dMolesPhase(j)
        end do
        print *


        ! Optional: cycle through metastable phases:
        do i = 1, nSolnPhasesSys

            if (lSolnPhases(i) .EQV. .FALSE.) print *, i, cSolnPhaseName(i), dDrivingForceSoln(i), lSolnPhases(i)

        end do
        print *

    end if

end subroutine GEMDebug
