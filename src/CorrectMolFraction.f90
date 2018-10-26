
!-------------------------------------------------------------------------------
!
!> \file    CorrectMolFraction.f90
!> \brief   Correct mole fractions of minor species by defining them from the
!!           element potentials.
!> \author  M.H.A. Piro
!> \date    August 20, 2015
!
!
! Revisions:
! ==========
! 
!   Date            Programmer      Description of change
!   ----            ----------      ---------------------
!   20/08/2015      M.H.A. Piro     Original code
!
!
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to correct the mole fractions
!! of minor species.  At equilibrium, the chemical potential of every stable
!! species and phase must be related to the chemical potentials of the system
!! components.  However, there might be some residual, which often happens 
!! with minor species.  This subroutine corrects for that error.
!
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


subroutine CorrectMolFraction

    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none
    
    integer :: i, j, k, l
    real(8) :: dTemp, dURF
    logical :: lCompEverything


    ! Initialize variables:
    lCompEverything = .FALSE.
    dURF = 0.5D0


    if (dGEMFunctionNorm > 1D-5) return


    print *, 'CorrectMolFraction'
    print *

    LOOP_PHASE: do j = 1, nSolnPhases
        k = -iAssemblage(nElements - j + 1)

        ! Compute the residual terms:
        LOOP_SPECIES: do i = nSpeciesPhase(k-1)+1, nSpeciesPhase(k)

            ! Compute the chemical potential term from the element potentials:
            dTemp = 0D0
            do l = 1, nElements
                dTemp = dTemp + dElementPotential(l) * dStoichSpecies(i,l)
            end do

            ! Normalize this quantity by the number of particles per mole:
            dTemp = dTemp / DFLOAT(iParticlesPerMole(i))

            ! Compute correction factor to mole fraction:
            dTemp = dTemp - dChemicalPotential(i)


print *, csolnPhaseName(k), cSpeciesName(i), dTemp, dMolFraction(i), dMolFraction(i) * DSQRT(DEXP(dTemp))


            dMolFraction(i) = DSQRT(DEXP(dTemp)) * dMolFraction(i)
            !dMolFraction(i)  = DEXP(dTemp) * dMolFraction(i)
            dMolesSpecies(i) = dMolesPhase(nElements-j+1) * dMolFraction(i)
            !dMolesSpecies(i) = dURF * dTemp + (1D0 - dURF) * dMolesSpecies(i)

        end do LOOP_SPECIES

    end do LOOP_PHASE

    call CompChemicalPotential(lCompEverything)
    
    return

end subroutine CorrectMolFraction