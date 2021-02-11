
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    InitGEMSolver.f90
    !> \brief   Initialize the GEMSolver.f90 subroutine.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      GEMSolver.f90
    !> \todo    Improve the estimation of mole fractions and mole species by mapping the results from Leveling to
    !!           here.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code
    !   09/13/2012      M.H.A. Piro         Apply subminimization when nSolnPhases = 0 and better organize
    !                                        this subroutine.
    !   10/21/2012      M.H.A. Piro         Previously, dMolFraction is initialized by 0.1; however, this means
    !                                        that the sum of dMolFraction within a solution phase does not equal
    !                                        unity.  Now, the mole fractions are initializes so that the sum
    !                                        equals unity.
    !   04/23/2013      M.H.A. Piro         Changed the initialization scheme for mole fractions to a constant
    !                                        value.  Computing the mole fractions from the element potentials
    !                                        can give very poor initial estimates.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to initialize the GEMSolver.f90 subroutine.  Specifically, this
    !! subroutine determines which pure condensed phases and solution phases are initially estimated to contribute
    !! to the equilibrium phase assemblage.  Initial estimates of the quantity of each phase was determined by
    !! the LevelingSolver.f90 subroutine.  Also, many allocatable arrays are allocated in this subroutine.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFOThermo            An integer scalar used by Thermochimica to identify a successful exit (0) or an error.
    ! nConPhases            The number of pure condensed phases in the assemblage
    ! nSolnPhases           The number of solution phases in the assemblage
    ! nSolnPhasesSys        The number of solution phases in the system
    ! iAssemblage           Integer vector containing the indices of phases in the assemblage
    !                        (1:nConphases represent pure condensed phases and (nElements-nSolnPhases:nElements)
    !                        represent solution phases.
    ! dMolesPhase           The number of moles of a phase.  These are directly mapped to phases in iAssemblage.
    ! dSumMolFractionSoln   A double real vector representing the sum of mole fractions in each solution phase.
    ! lPhasePass            A logical scalar representing whether a particular phase assemblage is appropriate
    !                        for further consideration.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine InitGEMSolver

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo, lReinitLoaded
    USE ModuleGEMSolver
    USE ModuleReinit

    implicit none

    integer::                               i, j, k, l
    integer,dimension(nElements)::          iAssemblageLast
    real(8)::                               dTemp, dSum
    real(8),dimension(-1:nSolnPhasesSys)::  dTempVec
    logical::                               lPhasePass, lCompEverything


    ! Check to see if allocatable arrays have already been allocated:
    if (allocated(dMolesSpecies)) then
        ! Check to see if the number of species has changed:
        i = SIZE(dMolesSpecies)
        if (i /= nSpecies) then
            deallocate(dMolesSpecies,dMolFraction,dPartialExcessGibbs,dPartialExcessGibbsLast, STAT = l)
            if (l /= 0) then
                INFOThermo = 21
                return
            end if
            allocate(dMolFraction(nSpecies),dMolesSpecies(nSpecies))
            allocate(dPartialExcessGibbs(nSpecies),dPartialExcessGibbsLast(nSpecies))
        end if

        ! Check to see if the number of elements has changed:
        j = SIZE(dUpdateVar)
        if (j /= nElements) then
            deallocate(dUpdateVar,&
                iterHistory, STAT = l)
            if (l /= 0) then
                INFOThermo = 21
                return
            end if
            allocate(dUpdateVar(nElements*2))
            allocate(iterHistory(nElements,iterGlobalMax))
        end if

        ! Check to see if the number of solution phases in the system has changed:
        k = SIZE(dSumMolFractionSoln)
        if (k /= nSolnPhasesSys) then
            deallocate(dSumMolFractionSoln, dDrivingForceSoln, STAT = l)
            if (l /= 0) then
                INFOThermo = 21
                return
            end if
            l = MAX(1,nSolnPhasesSys)
            allocate(dSumMolFractionSoln(l),dDrivingForceSoln(l))
        end if

        ! Check to see if either the number of solution phases or the number of elements has changed:
        if ((j /= nElements).OR.(k /= nSolnPhasesSys)) then
            deallocate(dEffStoichSolnPhase, STAT = l)
            if (l /= 0) then
                INFOThermo = 21
                return
            end if
            l = MAX(1,nSolnPhasesSys)
            allocate(dEffStoichSolnPhase(l,nElements))
        end if
    else

        ! Allocate memory:
        l = MAX(1,nSolnPhasesSys)
        allocate(dMolesSpecies(nSpecies))
        allocate(dPartialExcessGibbs(nSpecies),dPartialExcessGibbsLast(nSpecies))
        allocate(dUpdateVar(nElements*2))
        allocate(iterHistory(nElements,iterGlobalMax))
        allocate(dSumMolFractionSoln(l))
        allocate(dDrivingForceSoln(l))
        allocate(dEffStoichSolnPhase(l,nElements))
        ! This is due to reinit
        if (allocated(dMolFraction)) then
            ! Check to see if the number of species has changed:
            i = SIZE(dMolFraction)
            if (i /= nSpecies) then
                deallocate(dMolFraction, STAT = l)
                if (l /= 0) then
                    INFOThermo = 21
                    return
                end if
                allocate(dMolFraction(nSpecies))
            end if
        else
            allocate(dMolFraction(nSpecies))
        end if
    end if

    if (allocated(dMolesPhaseLast)) deallocate(dMolesPhaseLast)
    allocate(dMolesPhaseLast(nElements))

    ! Initialize variables:
    iterLast                = 0
    iterStep                = 5
    iterLastCon             = 0
    iterLastSoln            = 0
    iterHistory             = 0
    iterRevert              = 0
    iterSwap                = 0
    iterLastMiscGapCheck    = 0
    iterGlobal              = 0
    nConPhases              = 0
    nSolnPhases             = 0
    iConPhaseLast           = 0
    iSolnSwap               = 0
    iPureConSwap            = 0
    iSolnPhaseLast          = 0
    ! From LevelingSolver, iAssemblageLast would always be positive
    ! because only pure condensed phases are considered.
    ! However, when reiniting there may be any type of phase included,
    ! and therefore some indices may be negative.
    iAssemblageLast         = iAssemblage
    iAssemblage             = 0
    dMolesPhaseLast         = dMolesPhase
    dMolesPhase             = 0D0
    dMolFraction            = 0.1D0
    dMolesSpecies           = 0D0
    dTempVec                = 0D0
    dPartialExcessGibbs     = 0D0
    dPartialExcessGibbsLast = 0D0
    dDrivingForceSoln       = 0D0
    lCompEverything         = .FALSE.
    lConverged              = .FALSE.
    lRevertSystem           = .FALSE.

    ! Regrettably, branching here depending on whether reinit data has been loaded
    IF_ReinitLoaded: if (lReinitLoaded) then
      iAssemblage = iAssemblageLast
      dMolFraction = dMolFraction_Old
      do i = 1, nElements
        dMolesPhase = dMolesPhaseLast
        if (iAssemblageLast(i) < 0) then
          nSolnPhases = nSolnPhases + 1
          lSolnPhases(-iAssemblageLast(i)) = .TRUE.
        elseif (iAssemblageLast(i) > 0) then
          nConPhases  = nConPhases + 1
        end if
      end do
    ! If there is no reinit data use old methods:
    else
      ! Calculate the total number of moles for each solution phase:
      do i = 1, nElements
          j           = iPhase(iAssemblageLast(i))
          dTempVec(j) = dTempVec(j) + dMolesPhaseLast(i)
      end do

      ! Count the number of pure condensed phases and solution phases are assumed to be part of the phase
      ! assemblage and establish iAssemblage based on the results of Leveling and PostLeveling:
      LOOP_AddPhase: do i = 1, nElements
          if ((iPhase(iAssemblageLast(i)) == 0).AND.(nConPhases + nSolnPhases < nElements)) then
              nConPhases              = nConPhases + 1
              iAssemblage(nConPhases) = iAssemblageLast(i)
              dMolesPhase(nConPhases) = dMolesPhaseLast(i)
          elseif ((iPhase(iAssemblageLast(i)) > 0).AND.(nConPhases + nSolnPhases < nElements)) then
              do j = 1,nSolnPhases
                  k = nElements - j + 1
                  ! Ensure that this solution phase is not already stored:
                  if (iAssemblage(k) == -iPhase(iAssemblageLast(i))) cycle LOOP_AddPhase
              end do

              nSolnPhases = nSolnPhases + 1

              j = nElements - nSolnPhases + 1
              k = iPhase(iAssemblageLast(i))

              iAssemblage(j) = -k
              dMolesPhase(j) = DMAX1(dTempVec(k),dTolerance(9))
              lSolnPhases(k) = .TRUE.
          end if
      end do LOOP_AddPhase

      ! Initialize mole fractions of all solution phase constituents:
      LOOP_CompX: do k = 1, nSolnPhasesSys

          ! Reinitialize temporary variable:
          dSum = 0D0

          ! The default case assumes an ideal solution phase.
              do i = nSpeciesPhase(k - 1) + 1, nSpeciesPhase(k)
                  dTemp = 0D0
                  do j = 1, nElements
                      dTemp = dTemp + dElementPotential(j) * dStoichSpecies(i,j)
                  end do
                  dTemp           = dTemp / DFLOAT(iParticlesPerMole(i))
                  dMolFraction(i) = DEXP(dTemp - dStdGibbsEnergy(i) - dMagGibbsEnergy(i))
                  dMolFraction(i) = DMIN1(dMolFraction(i),1D0)
                  dSum            = dSum + dMolFraction(i)
              end do

              ! Normalize the mole fractions so that their sum equals unity:
              dSum = 1D0 / dSum
              do i = nSpeciesPhase(k - 1) + 1, nSpeciesPhase(k)
                  dMolFraction(i) = dMolFraction(i) * dSum
              end do

      end do LOOP_CompX
    end if IF_ReinitLoaded

    ! If there aren't any solution phases currently predicted to be stable, check if any should be added:
    if (nSolnPhases == 0) call InitGemCheckSolnPhase

    ! It may be possible that a single solution phase may be expected to form (as predicted by Leveling)
    ! but that there are zero moles of this particular phase.  This can cause problems establishing the
    ! Jacobian.
    if ((nSolnPhases == 1).AND.(dMolesPhase(nElements) == 0D0)) then

        ! Adjust the number of moles of pure condensed phases:
        dMolesPhase = dMolesPhase * 0.95D0

        ! Compute the number of moles of each solution phase and establish the Jacobian constraint vector
        call CompMolSolnPhase

    end if

    ! Now that the initial phase assemblage has been established, compute the number of moles of each
    ! solution phase constituent:
    do i = 1, nSolnPhases
        k = nElements - i + 1  ! Relative solution phase index
        l = -iAssemblage(k)    ! Absolute solution phase index

        do j = nSpeciesPhase(l-1) + 1, nSpeciesPhase(l)
            dMolesSpecies(j) = dMolesPhase(k) * dMolFraction(j)
            dMolesSpecies(j) = DMAX1(dMolesSpecies(j), 1D-300)
        end do

        ! Initialize the driving force for this phase:
        dDrivingForceSoln(l) = 0D0

    end do

    ! Compute the chemical potentials:
    call CompChemicalPotential(lCompEverything)

    ! Check the phase assemblage to make sure that the Jacobian matrix is appropriate.  Only do this if there is
    ! at least one solution phase:
    if (nSolnPhases > 0) then

        i = MAX(1,nConPhases)

        LOOP_CheckPhaseAssemblage: do j = 1, i

                ! Check to make sure that the phase can be added:
                call CheckPhaseChange(lPhasePass,k)

                if (k > nElements + nSolnPhases) then
                    ! A pure condensed phase should be removed.
                    k = k - nElements - nSolnPhases
                    iAssemblage(k)          = iAssemblage(nConPhases)
                    dMolesPhase(k)          = dMolesPhase(nConPhases)
                    iAssemblage(nConPhases) = 0
                    dMolesPhase(nConPhases) = 0D0
                    nConPhases = nConPhases - 1
                elseif (k == 0) then
                    exit LOOP_CheckPhaseAssemblage
                else
                    ! Placeholder...the phase assemblage has failed...
                    exit LOOP_CheckPhaseAssemblage
                end if
        end do LOOP_CheckPhaseAssemblage
    end if

    dGEMFunctionNorm = 1D3

    return

end subroutine InitGEMSolver


!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to check whether a solution phase
    ! should be added to the system.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! i             Absolute index of phase to be added to the system.
    ! lPhasePass    A logical scalar indicating whether the new phase assemblage
    !                has passed (TRUE) or not (FALSE).
    !
    !---------------------------------------------------------------------------


subroutine InitGemCheckSolnPhase

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::   i, j, k, l
    real(8)::   dTemp
    logical::   lPhasePass


    ! Initialize variables:
    lPhasePass = .FALSE.

    ! Loop through all solution phases in the system to check for a miscibiltiy gap:
    LOOP_Subminimization: do i = 1, nSolnPhasesSys

        call CheckMiscibilityGap(i,lPhasePass)

        ! If this phase should be added, then add it:
        IF_AddPhase: if (lPhasePass) then

            ! First, check if the solution phase should be added directly or if it should swap another phase:
            IF_PhaseRule: if (nConPhases + nSolnPhases < nElements) then
                ! This solution phase can be added directly:
                nSolnPhases            = nSolnPhases + 1
                k                      = nElements - nSolnPhases + 1
                iAssemblage(k)         = -i
                lSolnPhases(i)         = .TRUE.
                dSumMolFractionSoln(i) = 1D0
            else
                ! This solution phase should swap a pure condensed phase.

                dMolesPhase            = dMolesPhase * 0.95D0
                lSolnPhases(i)         = .TRUE.
                dSumMolFractionSoln(i) = 1D0

                ! Shuffle the phase assmeblage
                call ShuffleAssemblage(-i,j)

                ! Loop through pure condensed phases to see which one can be swapped:
                LOOP_ConPhases: do j = 1, nConPhases
                    dTemp                   = dMolesPhase(j)
                    k                       = iAssemblage(j)
                    iAssemblage(j)          = iAssemblage(nConPhases)
                    iAssemblage(nConPhases) = -i
                    nConPhases              = nConPhases - 1
                    nSolnPhases             = nSolnPhases + 1

                    ! Compute the number of moles of each solution phase and establish the Jacobian constraint vector
                    call CompMolSolnPhase

                    ! Check to make sure that the phase can be added:
                    call CheckPhaseChange(lPhasePass,l)

                    if (lPhasePass) then
                        ! This phase assemblage is appropriate for testing.
                        exit LOOP_ConPhases
                    else
                        ! This phase assemblage is not appropriate for testing.  Return to the previous assemblage.
                        dMolesPhase(j) = dTemp
                        iAssemblage(j) = k
                        nConPhases     = nConPhases  + 1
                        nSolnPhases    = nSolnPhases - 1
                    end if

                end do LOOP_ConPhases

            end if IF_PhaseRule

            ! Exit the loop:
            exit LOOP_Subminimization

        end if IF_AddPhase

    end do LOOP_Subminimization

end subroutine InitGemCheckSolnPhase


    !---------------------------------------------------------------------------
    !                       END - InitGEMSolver.f90
    !---------------------------------------------------------------------------
