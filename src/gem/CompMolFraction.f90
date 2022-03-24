
!-------------------------------------------------------------------------------
!
!> \file    CompMolFraction.f90
!> \brief   Compute the mole fraction of all solution phase constituents of a
!!           particular solution phase.
!> \author  M.H.A. Piro
!> \date    Apr. 26, 2012
!
! This subroutine is called by:
!> \sa      InitGEMSolver.f90
!> \sa      CompChemicalPotential.f90
!> \sa      CheckSolnPhaseAdd.f90
!> \sa      AddSolnPhase.f90
!> \sa      RemPureConAddSolnPhase.f90
!> \sa      SwapSolnPhase.f90
!> \sa      SwapSolnForPureConPhase.f90
!
!
! Revisions:
! ==========
!
!   Date            Programmer      Description of change
!   ----            ----------      ---------------------
!   31/03/2011      M.H.A. Piro     Original code
!   07/31/2011      M.H.A. Piro     Clean up code: remove unnecessary
!                                    variables, update variable names
!   10/25/2011      M.H.A. Piro     Clean up code: modules, simplify
!                                    programming.
!   02/06/2012      M.H.A. Piro     Provided the capability to handle
!                                    constituents with more than one particle
!                                    per mole.
!   04/26/2012      M.H.A. Piro     Implementing Gibbs energy Minimization
!                                    algorithm.
!   09/22/2012      M.H.A. Piro     The CompMolFractionQKTO subroutine was
!                                    replaced with Subminimization.  Solution
!                                    phases are now added based on the driving
!                                    force of that phase. See
!                                    Subminimization.f90 for more details of the
!                                    advantages of this approach.
!
!
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to compute mole fractions of all
!! solution phase constituents, the sum of mole fractions of all constituents
!! in a solution phase (i.e., dSumMolFractionSoln) and the effective
!! stoichiometry of each solution phase (i.e., dEffStoichSolnPhase).  The
!! relationship between the mole fraction of a solution phase constituent and
!! its chemical potential depends on the type of solution phase.  The solution
!! phase types that are supported follow:
!
! Solution phase types (cSolnPhaseType):
! --------------------------------------
!
!!    QKTO     -      Quasichemical Kohler-Toop;
!!    IDMX     -      Ideal mixing (default).
!
!
! Pertinent variables:
! ====================
!
!> \param[in]       k       Absolute solution phase index.
!
! nSolnPhasesSys            Number of solution phases in the system.
! nSpeciesPhase             Index of last species in a particular solution phase.
! dMolFraction              Double real vector representing the current
!                            estimated mole fraction.
! dSumMolFractionSoln       Sum of all mole fractions within a solution phase.
! dEffStoichSolnPhase       Effective stoichiometry of a particular element in
!                            a solution phase.
! dPartialExcessGibbs       The partial molar excess Gibbs energy of mixing of
!                            a phase constituent.
! dChemicalPotential        A double real vector representing the chemical
!                            potential of each species.  To be precise, this is
!                            defined as the difference between the standard
!                            molar Gibbs energy and the chemical potential
!                            defined by the element potentials (represented in
!                            dimensionless units and per formula mass).
! iParticlesPerMole         The number of particles per mole of constituent.
! dStoichSpecies            The number of atoms of a particular element for a
!                            particular species.
! dSpeciesTotalAtoms        The total number of atoms per formula mass of a
!                            species.
! cSolnPhaseType            Character vector representing the type of solution
!                            phase (see above).
!
!-------------------------------------------------------------------------------


subroutine CompMolFraction(k)

    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, l, m, n, nActiveSpecies, INFO, iFirst, iLast
    real(8) :: dTemp, dDrivingForceTemp, dSum, rnorm
    logical :: lPhasePass
    integer,dimension(:), allocatable           :: indx
    real(8),dimension(nElements)                :: dMolesElementTemp
    real(8),dimension(:), allocatable           :: x, work, dComparisonEnergy
    real(8),dimension(:,:), allocatable         :: dStoichAvail, dStoichSpeciesTemp

    ! Although this variable is not used in this subroutine, initialize the
    ! variable for sake of consistency.
    lPhasePass = .FALSE.

    ! Reset sum of mole fractions and the effective stoichiometry:
    dSumMolFractionSoln(k)             = 0D0
    dEffStoichSolnPhase(k,1:nElements) = 0D0
    dDrivingForceTemp                  = 0D0

    ! Check if the input masses can be represented in terms of the available species
    nActiveSpecies = nConPhases
    do i = 1, nSolnPhases
        j = -iAssemblage(nElements - i + 1)
        nActiveSpecies = nActiveSpecies + nSpeciesPhase(j) - nSpeciesPhase(j-1)
    end do

    allocate(indx(nActiveSpecies), x(nActiveSpecies), work(nActiveSpecies))
    allocate(dStoichAvail(nElements,nActiveSpecies),dStoichSpeciesTemp(nElements,nActiveSpecies))
    n = 0
    do i = 1, nSolnPhases
        j = -iAssemblage(nElements - i + 1)
        do m = nSpeciesPhase(j-1) + 1,  nSpeciesPhase(j)
            n = n + 1
            do l = 1, nElements
                dStoichAvail(l,n) = dStoichSpecies(m,l)
            end do
        end do
    end do

    do i = 1, nConPhases
        j = iAssemblage(i)
        n = n + 1
        do l = 1, nElements
            dStoichAvail(l,n) = dStoichSpecies(j,l)
        end do
    end do

    ! Initialize variables:
    iFirst = nSpeciesPhase(k-1) + 1      ! First constituent in phase.
    iLast = nSpeciesPhase(k)            ! Last  constituent in phase.
    allocate(dComparisonEnergy(iLast + 1 - iFirst))

    ! Compute the mole fraction depending on the type of solution phase:
    select case (cSolnPhaseType(k))

        ! Leave this for now until I improve the estimator:

        case ('QKTO','RKMP','RKMPM','SUBL','SUBLM','SUBG','SUBQ')

            ! Perform subminimization:
            call Subminimization(k,lPhasePass)

        case default

            ! The default case assumes an ideal solution phase.
            dSum = 0D0
            do i = iFirst, iLast
                dMolesElementTemp = dStoichSpecies(i,:)
                dStoichSpeciesTemp = dStoichAvail
                call nnls(dStoichSpeciesTemp, nElements, nActiveSpecies, dMolesElementTemp, x, rnorm, work, indx, INFO)

                if (rnorm > 1D-12) then
                    dTemp = -1D10
                else
                    l = 0
                    dTemp = 0D0
                    do n = 1, nSolnPhases
                        j = -iAssemblage(nElements - n + 1)
                        do m = nSpeciesPhase(j-1) + 1,  nSpeciesPhase(j)
                            l = l + 1
                            dTemp = dTemp + x(l) * dChemicalPotential(m)
                        end do
                    end do

                    do n = 1, nConPhases
                        j = iAssemblage(n)
                        l = l + 1
                        dTemp = dTemp + x(l) * dChemicalPotential(j)
                    end do
                end if
                dTemp = dTemp / DFLOAT(iParticlesPerMole(i))
                dComparisonEnergy(i + 1 - iFirst) = dTemp
                dMolFraction(i) = DMIN1(dTemp - dStdGibbsEnergy(i),0D0)
                dMolFraction(i) = DEXP(dMolFraction(i))
                dSum = dSum + dMolFraction(i)
            end do
            ! Get total moles in gas phase
            if (dSum > 1D0) then
                do i = iFirst, iLast
                     dMolFraction(i) = dMolFraction(i) / dSum
                end do
                dSum = 1D0
            end if
            ! Calculate driving forces
            do i = iFirst, iLast
                dTemp = dComparisonEnergy(i + 1 - iFirst)
                dDrivingForceTemp = dDrivingForceTemp + dMolFraction(i)/dSum * &
                                    (dStdGibbsEnergy(i) + DLOG(DMAX1(dMolFraction(i)/dSum, 1D-75)) - dTemp)
            end do

    end select


    ! Compute the sum of all mole fractions in this phase and the effective
    ! stoichiometry:
    do i = iFirst, iLast
        dSumMolFractionSoln(k) = dSumMolFractionSoln(k) + dMolFraction(i)
        do j = 1,nElements
            dEffStoichSolnPhase(k,j) = dEffStoichSolnPhase(k,j) + dMolFraction(i) * &
                dStoichSpecies(i,j) / DFLOAT(iParticlesPerMole(i))
        end do
    end do

    ! Check for a NAN and record an error:
    if (dSumMolFractionSoln(k) /= dSumMolFractionSoln(k)) INFOThermo = 25

    ! Compute the driving force for ideal mixtures (this is computed by
    ! Subminimization for other solution phases):
    if (cSolnPhaseType(k) == 'IDMX') then
        do i = 1, nElements
            ! phase already in assemblage should not have driving force
            if (-iAssemblage(i) == k) dDrivingForceTemp = 0D0
        end do
        dDrivingForceSoln(k) = dDrivingForceTemp
    end if

    deallocate(indx, x, work, dStoichSpeciesTemp, dStoichAvail, dComparisonEnergy)

    return

end subroutine CompMolFraction
