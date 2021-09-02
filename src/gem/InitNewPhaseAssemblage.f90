
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    InitNewPhaseAssemblage.f90
    !> \brief   Initialize a new phase assemblage.
    !> \author  M.H.A. Piro
    !> \date    Apr. 01, 2013
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/01/2013      M.H.A. Piro         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to initialize the system when the phase assemblage has been
    !! updated.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine InitNewPhaseAssemblage

    USE ModuleThermo
    uSE ModuleGEMSolver

    implicit none

    integer                               :: i, j, k, l, m, n, c, d, e, iter, INFO, LWORK, nVar
    real(8)                               :: dSteplength, dTemp
    real(8),  dimension(:,:), allocatable :: A
    real(8),  dimension(:),   allocatable :: B
    real(8),  dimension(:),   allocatable :: WORK
    logical                               :: lPhasePass
    character                             :: TRANS


    ! Initialize variables:
    M           = 0
    N           = 0
    lPhasePass  = .FALSE.
    TRANS       = 'N'
    dSteplength = 1D0

    ! Loop through solution phases to count the # of solution species:
    do j = 1, nSolnPhases
        k = -iAssemblage(nElements - j + 1)
        M = M + (nSpeciesPhase(k) - nSpeciesPhase(k-1))
    end do
    M = M + nConPhases

    ! Loop through solution phases to count the number of ionic phases that are stable:
    i = 0
    if (nChargedConstraints > 0) then
        do j = 1, nSolnPhases
            k = -iAssemblage(nElements - j + 1)
            if (iPhaseElectronID(k) > 0) i = i + 1
        end do
    end if

    ! Define dimension variables:
    nVar  = nChargedConstraints - i
    N     = nElements - nVar
    LWORK = M + N

    ! Deallocate allocatable arrays (if necessary):
    if (allocated(A))    deallocate(A)
    if (allocated(B))    deallocate(B)
    if (allocated(WORK)) deallocate(WORK)

    ! Allocate allocatable arrays:
    allocate(A(M,N), B(M), WORK(LWORK))

    ! Initialize variables:
    INFO = 0
    WORK = 0D0
    A    = 0D0
    B    = 0D0

    ! Iterate:
    LOOP_Iter: do iter = 1, 1

        ! Reinitialize variables:
        c           = 0
        d           = 0
        dSteplength = 1D0

        ! Loop through stable solution phases to compute the lowest driving force for each phase:
        LOOP_Soln: do j = 1, nSolnPhases

            ! Absolute solution phase index:
            k = -iAssemblage(nElements - j + 1)

            ! Electron ID corresponding to phase:
            e = iPhaseElectronID(k)

            ! Store counter:
            if (e > 0) d = d + 1

            ! Compute the driving force of this solution phase:
            call Subminimization(k,lPhasePass)

            ! Loop through solution species:
            LOOP_SolnSpecies: do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                c = c + 1
                do l = 1, nElements - nChargedConstraints
                    A(c,l) = dStoichSpecies(i,l)
                end do

                ! Add coefficient if this is an ionic phase:
                if (e > 0) A(c,d) = dStoichSpecies(i,l)

                ! Add right hand side:
                B(c) = dChemicalPotential(i)

            end do LOOP_SolnSpecies

        end do LOOP_Soln

        ! Construct matrix A representing the stoichiometry of each pure condensed phase:
        LOOP_PureCon: do j = 1, nConPhases
            k = iAssemblage(j)
            c = c + 1
            do l = 1, nElements - nChargedConstraints
                A(c,l) = dStoichSpecies(k,l)
            end do
            B(c) = dStdGibbsEnergy(k)
        end do LOOP_PureCon

        ! Perform linear regression on the element potentials:
        call DGELS( TRANS, M, N, 1, A, M, B, M, WORK, LWORK, INFO )

        ! Compute the steplength:
        do l = 1, nElements - nChargedConstraints
            dTemp = DABS(dElementPotential(l) - B(l))
            if (dTemp > 1D0) then
                dTemp = 1D0 / dTemp
                dSteplength = DMIN1(dSteplength, dTemp)
            end if
        end do

        ! Update the element potentials:
        do l = 1, nElements - nChargedConstraints
            dElementPotential(l) = (1D0 - dSteplength) * dElementPotential(l) + dSteplength * B(l)
        end do

    end do LOOP_Iter

    ! Recompute the molar quantity of each solution phase.
    call CompMolAllSolnPhases

    ! Loop through stable solution phases to update the molar quantity of each solution species:
    do j = 1, nSolnPhases
        l = nElements - j + 1
        k = -iAssemblage(l)

        ! Loop through species in this phase:
        do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
            dMolesSpecies(i) = dMolFraction(i) * dMolesPhase(l)
        end do

    end do

    ! Deallocate allocatable arrays:
    deallocate(A, B, WORK)

    return

end subroutine InitNewPhaseAssemblage
