
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

    integer :: i, j, k, m, n
    real(8) :: dTemp, dDrivingForceTemp
    logical :: lPhasePass


    ! Initialize variables:
    m = nSpeciesPhase(k-1) + 1      ! First constituent in phase.
    n = nSpeciesPhase(k)            ! Last  constituent in phase.

    ! Although this variable is not used in this subroutine, initialize the
    ! variable for sake of consistency.
    lPhasePass = .FALSE.

    ! Reset sum of mole fractions and the effective stoichiometry:
    dSumMolFractionSoln(k)             = 0D0
    dEffStoichSolnPhase(k,1:nElements) = 0D0
    dDrivingForceTemp                  = 0D0

    ! Compute the mole fraction depending on the type of solution phase:
    select case (cSolnPhaseType(k))

        ! Leave this for now until I improve the estimator:

        case ('QKTO','RKMP','RKMPM','SUBL','SUBLM','SUBG','SUBQ')

            ! Perform subminimization:
            call Subminimization(k,lPhasePass)

        case default

            ! The default case assumes an ideal solution phase.
            do i = m, n
                dTemp = 0D0
                do j = 1, nElements
                    dTemp = dTemp + dElementPotential(j) * dStoichSpecies(i,j)
                end do
                dTemp           = dTemp / DFLOAT(iParticlesPerMole(i))
                dMolFraction(i) = DEXP(dTemp - dStdGibbsEnergy(i))
                dMolFraction(i) = DMIN1(dMolFraction(i),1D0)
                dDrivingForceTemp = dDrivingForceTemp + dMolFraction(i) * (dStdGibbsEnergy(i) - dTemp)
            end do

    end select


    ! Compute the sum of all mole fractions in this phase and the effective
    ! stocihiometry:
    do i = m, n
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
        dDrivingForceSoln(k) = dDrivingForceTemp
    end if

    return

end subroutine CompMolFraction
