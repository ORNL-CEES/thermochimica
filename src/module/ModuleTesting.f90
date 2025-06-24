    
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ModuleTesting.F90
    !> \brief   Module for testing.
    !> \author  A.E.F. Fitzsimmons
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    09/03/2024    A.E.F. Fitzsimmons   Original Code
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to run tests with cleaner code, checking gibbs energy,
    !! heat capacity and another function to test all mole fractions within a stable phase.
    !
    !-------------------------------------------------------------------------------------------------------------

module ModuleTesting

    USE ModuleThermo
    USE ModuleThermoIO
    
    implicit none

    contains

    subroutine testProperties(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)
        USE ModuleThermo

        !Init variables
        integer, intent(in):: nSpeciesTest
        integer, allocatable, intent(in) :: iSpeciesIndexTest(:)
        real(8), allocatable, intent(in) :: dMolFractionTest(:)
        real(8), intent(in) :: dGibbsCheck, dHeatCapacityCheck
        logical, intent(out) :: lPass
        real(8) :: dToleranceCheck
        integer :: i
        
        dToleranceCheck = 1D-3

        !Test Gibbs and HeatCapacity
        if (DABS((dGibbsEnergySys - dGibbsCheck)/dGibbsCheck) >= dToleranceCheck) then
            print *, "Gibbs Energy Test Failed"
            return
        elseif (DABS(dHeatCapacity - dHeatCapacityCheck)/dHeatCapacityCheck >= dToleranceCheck) then 
            print *, "Heat Capacity Test Failed"
            return
        end if 
        
        !Test Mole Fraction
        do i = 1, nSpeciesTest
            if (DABS(dMolFraction(iSpeciesIndexTest(i)) - dMolFractionTest(i))/dMolFractionTest(i) >= dToleranceCheck) then 
                print *, "Mole fraction Test Failed" 
                return
            end if
        end do

        lPass = .TRUE.

    end subroutine testProperties

    subroutine printMolFractions
        USE ModuleThermo
        integer :: i 

        do i = 1, nSpecies
            print *, i, cSpeciesName(i), dMolFraction(i)
        end do
    
    end subroutine printMolFractions

    subroutine checkTransitionTest(dPhaseTransitionTemp, dTestTransitionTemp, dTolerance, lPass)
        implicit none
        real(8), intent(in) :: dPhaseTransitionTemp(:)
        real(8), intent(in) :: dTestTransitionTemp(:)
        real(8), intent(in) :: dTolerance
        logical, intent(out) :: lPass
        integer :: i, n
        real(8) :: dDiff
    
        n = size(dPhaseTransitionTemp)
        lPass = .TRUE.
    
        do i = 1, n
            ! Only check where the expected value is non-zero
            if (dTestTransitionTemp(i) > 0D0) then
                dDiff = abs(dPhaseTransitionTemp(i) - dTestTransitionTemp(i))
                if (dDiff > dTolerance) then
                    print *, 'Mismatch at index ', i, ': calculated = ', dPhaseTransitionTemp(i), &
                             ', expected = ', dTestTransitionTemp(i), ', dDiff = ', dDiff
                    lPass = .FALSE.
                end if
            end if
        end do
    
    end subroutine checkTransitionTest
end module ModuleTesting