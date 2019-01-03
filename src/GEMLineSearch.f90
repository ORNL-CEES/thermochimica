
    !---------------------------------------------------------------------------
    !
    !> \file    GEMLineSearch.f90
    !> \brief   Perform a line search for the GEMSolver.f90
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      GEMSolver.f90
    !> \sa      GEMNewton.f90
    !> \sa      CompChemicalPotential.f90
    !> \sa      CompFunctionNorm.f90
    !> \todo    Figure out a long term solution for the condition to allow the solution
    !!           to only perform 1 iteration if the functional norm is less than 1E-6.
    !! \todo    Figure out a more permanent solution to how dMaxGamma shoudl be handled.
    !!           I think that I need to get away from that and focus on ensuring that the 
    !!           Wolfe conditions have been satisfied.
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !   04/25/2012      M.H.A. Piro     Original code
    !   05/02/2012      M.H.A. Piro     Improved calculation of initial step 
    !                                    length (constrain changes to the Element 
    !                                    potentials).
    !   05/08/2012      M.H.A. Piro     Improved calculation of initial step length
    !                                    (constrain the maximum change to the total
    !                                    number of moles of a solutoin phase).
    !   04/11/2013      M.H.A. Piro     When computing a steplength that constrains
    !                                    the maximum change to the element potential
    !                                    to less than or equal to unity, exclude elements
    !                                    with zero moles.  In other words, exclude 
    !                                    electrons corresponding to ionic phases
    !                                    that are not currently stable.  This is also
    !                                    done when updating the element potentials.
    !   05/08/2013      M.H.A. Piro     Changed tolerance for calculating the steplength
    !                                    when the number of moles of a species tends to zero
    !                                    to 1D-50 from 1D-100.
    !   06/08/2013      M.H.A. Piro     In determining an initial steplength, do not constrain
    !                                    the maximum decrease in the number of moles of a 
    !                                    solution phase by a certain increment (e.g., 50%) if
    !                                    the number of moles of that phase is below a certain
    !                                    value (e.g., 10**(-9)).  The motivation for doing this
    !                                    is that a solution phase may be driving out of the system
    !                                    but it may be inhibited if this condition is not made.
    !   06/08/2013      M.H.A. Piro     Exit the Wolfe loop if the functional norm is below a 
    !                                    certain tolerance (e.g., 10**(-6)).
    !   04/02/2014      M.H.A. Piro     I changed one of the conditions to satisfy the line search.
    !                                    Specifically, if the relative change of the functional 
    !                                    norm is 0.95 < F_norm < 1.0 to 0.97 < F_norm < 1.0.  
    !                                    Consider the following scenario: the functional norm
    !                                    is 1 at global iteration 5 and currently we are at iter 6.
    !                                    The f-norm of the first line search loop is 6.1 and then
    !                                    its 5.86 (dTemp = 0.961).  Before, the system would then
    !                                    exit, but the f-norm on the global scale is not converging.
    !   08/20/2015      M.H.A. Piro     Previously, the maximum change in the element potentials was
    !                                    constrained to 1.  Now, it is a linear function that varies
    !                                    with respect to the functional norm.  The motivation for this
    !                                    is that an initially poor guess may be far from equilibrium,
    !                                    and laxing this constraint accelerates convergence.  Another
    !                                    change is the way that dStepLength is initialized when 
    !                                    dMolesSpecies tends to zero.  Previously, this would compute 
    !                                    the minimum dStepLength that corresponds to reducing dMolesSpecies
    !                                    by 100.  Now, this is NOT applied when dMolesSpecies is incredibly
    !                                    small (e.g., less than numerical tolerance).
    !   12/18/2018      M.H.A. Piro     I effectively removed the constraint applied to the maximum
    !                                    change to the element potentials. I think that this is inefficient
    !                                    because it's probably best to leave it to ensuring that the
    !                                    Wolfe conditions have been satisfied. Furthermore, it seems to 
    !                                    be an issue for SUBG phases, which are effectively on a knife edge.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to perform a line search using 
    !! the direction vector computed by the Newton/Broyden solver.  The system 
    !! is updated using an appropriate step length that satisfies the Wolfe 
    !! conditions.  Specifically, values of dChemicalPotential and dMolesPhase 
    !! are updated.  It is possible for the system of equations to be ill-
    !! behaved and yield inappropriate results.  An initial step-length is 
    !! computed by normalizing the largest change of the system variables by a 
    !! pre-defined value.  The maximum change to the element potentials is 1 and
    !! the maximum change to the number of moles of a solution phase is twice of
    !! the previous value.  For more information, refer to Chapter 6 of the 
    !! above reference.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nElements             The number of elements in the system.
    ! nConPhases            The number of pure condensed phases in the system.
    ! nSolnPhases           The number of solution phases in the system.
    ! dUpdateVar            A double real vector that contains updates to the 
    !                        system variables (element potentials, moles of 
    !                        solution phases and moles of pure condensed phases).
    ! dStepLength           Step length applied to the direction vector 
    !                        (dUpdateVar)
    ! dMolesPhase           A double real vector representing the number of
    !                        moles of phases predicted to be stable.
    ! dLevel                The adjustment applied to the chemical potentials of
    !                        the elements
    ! iPhaseDampen          Integer vector that counts the number of times that 
    !                        the number of moles
    !                        of a solution phase had to be dampened.
    ! dChemicalPotential    A double real vector representing the chemical 
    !                        potential of each species and pure condensed phase.
    ! dGEMFunctionNorm      A double real scalar representing the norm of the 
    !                        functional vector in the PGESolver.
    ! dGEMFunctionNormLast  A double real scalar representing the norm of the 
    !                        functional vector from the previous iteration.
    !
    !---------------------------------------------------------------------------


subroutine GEMLineSearch

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none
    
    integer                       :: iterWolfe 
    real(8)                       :: dStepLength, dTemp, dWolfeFunctionNormLast
    real(8), dimension(nElements) :: dElementPotentialLast
    real(8), dimension(nSpecies)  :: dMolesSpeciesLast
    logical                       :: lCompEverything

        
    ! Initialize variable:
    dWolfeFunctionNormLast  = 1D-12
    dMolesSpeciesLast       = dMolesSpecies
    dElementPotentialLast   = dElementPotential 
    dGEMFunctionNormLast    = dGEMFunctionNorm
    dPartialExcessGibbsLast = dPartialExcessGibbs
    lCompEverything         = .FALSE.

    ! Initialize the line search method:
    call InitGEMLineSearch(dStepLength,dMolesSpeciesLast,dElementPotentialLast)

    ! Commence line search:
    LOOP_WOLFE: do iterWolfe = 1, 5

        ! Compute the fractional change in the functional norm:
        dTemp = dGEMFunctionNorm / dWolfeFunctionNormLast
    
        ! If the functional norm is already small, call it a day:
        if (dGEMFunctionNorm < 1D-6) exit LOOP_WOLFE

        ! If the number of moles of a solution phase is tending to zero, call it a day:
        if (MINVAL(dMolesPhase(nElements - nSolnPhases + 1: nElements)) < 1D-9) exit LOOP_WOLFE

        ! Check if the system is diverging or has sufficiently progressed, otherwise dampen:
        if (dGEMFunctionNorm < 0.999D0 * dGEMFunctionNormLast) then 

            ! The system has sufficiently progressed, exit:                
            exit LOOP_WOLFE
            
        elseif ((iterWolfe > 1).AND.(dTemp >= 0.999D0).AND.(dTemp <= 1D0)) then
         
            ! Return values to previous line search iteration:
            
            ! Update the system variables:
            dStepLength = 2D0
            
            call UpdateSystemVariables(dStepLength,dMolesSpeciesLast,dElementPotentialLast)
            
            ! Compute the chemical potentials of solution species:
            call CompChemicalPotential(lCompEverything)
       
            ! Compute the functional norm
            call CompFunctionNorm
            
            ! Values are returned to their previous values and exit:
            exit LOOP_WOLFE
            
        !elseif ((iterWolfe > 1).AND.(dTemp >= 0.95D0).AND.(dTemp <= 1D0)) then
        elseif ((iterWolfe > 1).AND.(dTemp >= 0.97D0).AND.(dTemp <= 1D0)) then
        
            ! The system has sufficiently progressed, exit:                
            exit LOOP_WOLFE
            
        elseif ((iterWolfe > 1).AND.((dGEMFunctionNorm >= dWolfeFunctionNormLast))) then
            
            ! Return values to previous line search iteration:
            
            ! Update the system variables:
            dStepLength = 2D0
            
            call UpdateSystemVariables(dStepLength,dMolesSpeciesLast,dElementPotentialLast)
            
            ! Compute the chemical potentials of solution species:
            call CompChemicalPotential(lCompEverything)
        
            ! Compute the functional norm
            call CompFunctionNorm
            
            ! Values are returned to their previous values and exit:
            exit LOOP_WOLFE
            
        else
            ! Dampen the system variables:
            
            ! If the functional norm has only increased by a nominal amount (i.e., 1%), then exit:
            dTemp = dGEMFunctionNorm / dWolfeFunctionNormLast
            if ((dTemp > 1D0).AND.(dTemp < 1.05D0)) exit LOOP_WOLFE
                                                            
            dStepLength = 0.5D0
            
            call UpdateSystemVariables(dStepLength,dMolesSpeciesLast,dElementPotentialLast)
            
            ! Compute the chemical potentials of solution species:
            call CompChemicalPotential(lCompEverything)
            
            dWolfeFunctionNormLast = dGEMFunctionNorm
            
            ! Compute the functional norm:
            call CompFunctionNorm
            
            ! Reiterate:
            cycle LOOP_WOLFE

        end if
        
    end do LOOP_WOLFE

    return
    
end subroutine GEMLineSearch


!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    ! 
    ! Purpose:
    ! ========
    ! 
    ! The purpose of this subroutine is to initialize the line search algorithm.
    ! Specifically, the initial step length needs to be determined before the 
    ! line search loop starts.
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !
    !   05/08/2012      M.H.A. Piro     Original Code
    !   07/04/2012      M.H.A. Piro     If the number of moles of any soln
    !                                    phase is to be significantly reduced
    !                                    and the system is not stagnant, then
    !                                    then further dampen the system.
    !   09/29/2012      M.H.A. Piro     Revert the system if the phase assemblage
    !                                    has not changed in 50 iterations and
    !                                    the maximum change to the system 
    !                                    variables is extremely large.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dStepLength       A double real scalar representing the step length.
    ! dMaxIncrease      A double real scalar representing the maximum 
    !                    increase in the number of moles of a solution phase.
    ! dMaxDecrease      A double real scalar representing the maximum 
    !                    decrease in the number of moles of a solution phase.
    ! lCompEverything   A logical scalar indicating whether everything should
    !                    be computed in a particular subroutine (true) or 
    !                    not (false). 
    !
    !---------------------------------------------------------------------------


subroutine InitGEMLineSearch(dStepLength,dMolesSpeciesLast,dElementPotentialLast)

    USE ModuleThermo
    USE ModuleGEMSolver
    
    implicit none
    
    integer                       :: i, j, k, l, nMisciblePhases
    real(8)                       :: dStepLength, dTemp, dMaxIncrease, dMaxDecrease, dMaxChange, dMaxGamma
    real(8), dimension(nElements) :: dElementPotentialLast
    real(8), dimension(nSpecies)  :: dMolesSpeciesLast
    logical                       :: lCompEverything

    
    ! Initialize variables:
    lCompEverything  = .FALSE.
    dStepLength      = 1D0
    dMolesPhaseLast  = dMolesPhase
    nMisciblePhases  = 0
    
    ! Count the number of stable miscible phases:
    do j = 1, nSolnPhases
        k = -iAssemblage(nElements - j + 1)
        if (lMiscibility(k) .EQV. .TRUE.) nMisciblePhases = nMisciblePhases + 1
    end do

    ! Initialize the maximum increase/decrease of functional variables
    ! (this is dapenned if the number of iterations gets large):
    if (iterGlobal < 1000) then
        dMaxIncrease = 2D0
        dMaxDecrease = 0.5D0
    elseif (iterGlobal < 2000) then
        dMaxIncrease = 1.5D0
        dMaxDecrease = 0.75D0
    else
        dMaxIncrease = 1.25D0
        dMaxDecrease = 0.85D0
    end if

    ! TEMPORARY: I should figure out a better and more permanent solution. My gut tells me that
    ! I need to go away from constraining the maximum change to the element potentials and a more
    ! elegant approach would be leaving it to the Wolfe conditions. For now, do this:
    if (iterGlobal > 1000) then
        dMaxGamma = (5D0 - 1D0) / (100D0 - 25D0) * (dGEMFunctionNorm - 100D0) + 5D0
        dMaxGamma = DMIN1(dMaxGamma,5D0)
        dMaxGamma = DMAX1(dMaxGamma,1D0)
    else
        dMaxGamma = 100D0
    end if

    ! Update the number of moles of pure condensed phases:
    do i = 1, nConPhases
        j = nElements + nSolnPhases + i 
        dMolesPhase(i) = dUpdateVar(j)
    end do

    ! Initialize the step length (prevent moles from being negative):
    do l = 1, nSolnPhases
        k = -iAssemblage(nElements - l + 1)     ! Absolute solution phase index.

        do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
            dTemp = 0D0
            do j = 1, nElements
                dTemp = dTemp + dUpdateVar(j) * dStoichSpecies(i,j)
            end do
            
            dTemp = dTemp / DFLOAT(iParticlesPerMole(i))
            
            ! NOTE: The variable dMolFraction is used temporarily to represent 
            ! the fractional update to dMolesSpecies, but it does not replace
            ! dMolesSpecies in the event that further dampening is required.
            dMolFraction(i) = (1D0 + dUpdateVar(nElements + l) + dTemp - &
                dChemicalPotential(i))

            if (dMolFraction(i) /= 1D0) dTemp = 1D0 / (1D0 - dMolFraction(i))
            
            if ((dMolFraction(i) < 0D0).AND.(dTemp < dStepLength)) then
                if (dMolesSpecies(i) < 1D-50) then
                    ! This species is very small and the system is trying to make it negative.
                    ! Reduce the mass of this species by an arbitrary positive value less than
                    ! unity:
                    dMolFraction(i) = 0.5D0

!!!! TEMPORARY   !!!!!
!!!! TEMPORARY   !!!!!
!!!! TEMPORARY   !!!!!

                elseif ((dMolesSpecies(i) / dMolesPhaseLast(nElements - l + 1)) < 1D-10) then

                    dMolFraction(i) = 0.1D0


                else

                    dStepLength = dTemp * 0.99D0
                end if
            end if
            dMolesSpecies(i) = dMolesSpecies(i) * dMolFraction(i)
            
        end do
    end do

    ! Initialize the steplength (constrain the element potentials to only 
    ! change by dMaxGamma):
    do i = 1, nElements
        dTemp = DABS(dElementPotential(i) - dUpdateVar(i))
        
        if (dUpdateVar(i) == 0D0) cycle
        
        if (dTemp > 1D0) then

            dTemp = dMaxGamma / dTemp
            !dTemp = 1D0 / dTemp
            dStepLength = DMIN1(dTemp, dStepLength)
        end if
    end do

	! Update the element potentials:
    do i = 1, nElements
        dElementPotential(i) = dUpdateVar(i)
    end do

    ! Update the system variables:
    call UpdateSystemVariables(dStepLength,dMolesSpeciesLast,dElementPotentialLast)

    ! Constrain the number of moles of any solution phase to only change by a factor of 2:
    i = 0
    dStepLength = 1D0
    do j = 1, nSolnPhases
        k     = nElements - j + 1     ! Solution phase index in iAssemblage and dMolesPhase

        if (dMolesPhase(k) /= dMolesPhaseLast(k)) then
            dTemp = dMolesPhaseLast(k) * (dMaxIncrease - 1D0) / (dMolesPhase(k) - dMolesPhaseLast(k))
        end if

        if (dTemp < 0D0) dTemp = dStepLength
        dStepLength = DMIN1(dTemp, dStepLength)        

        ! TEMPORARY TO AVOID AN INF:
        dMolesPhaseLast(k) = DMAX1(dMolesPhaseLast(k),1D-10)

        ! Do not allow the number of moles of a solution phase to be reduced by less than half if it is
        ! not being pushed to zero:
        dTemp = DABS(dMolesPhase(k) / dMolesPhaseLast(k))
        
        ! Count the number of times the number of moles of a solution phase is trying to increase
        ! by a large margin:
        if (dTemp > dMaxIncrease) i = i + 1


! THIS FOLLOWING SECTION NEEDS TO BE BETTER FIGURED OUT.  WHY WOULD IT DAMPEN IF A SOLUTOIN PHASE IS 
! BECOMING NEGATIVE WITH DTEMP > 0.011 OR DTEMP < 0.009??


        if (dTemp < 1D0) then
        !if (dTemp < 0.1D0) then
            if ((dTemp > 0.05D0).OR.(dTemp < 0.009D0)) then
            !if ((dTemp > 0.011D0).OR.(dTemp < 0.009D0)) then

                ! TEMPORARY:
                if (dMolesPhase(k) < 1D-9) Cycle 

                ! Count the number of solution phases that have molar quantities that are tending to zero:
                if (dTemp <= dMaxDecrease) i = i + 1    
                   
                ! Determine a step length that would result in halving this solution phase:
                dTemp = DABS(dMaxDecrease * dMolesPhaseLast(k) / (dMolesPhaseLast(k) - dMolesPhase(k)))
                dStepLength = DMIN1(dTemp,dStepLength)
            end if
        end if
    end do

    ! Check if there is at least one solution phase that is tending to zero:
    if (i > 0) then
        ! An issue with the Gibbs Energy Minimization method is that one is effectively attempting to 
        ! minimize two objective functions simultaneously: 1) the integral Gibbs energy of the system, 
        ! and 2) the residual vector of the mass balance equations.  A common problem is that there 
        ! may be a phase change that results in significant changes in both the mass balance residuals
        ! and the Gibbs energy function, making the functional norm a poor indicator of convergence.
        ! Generally, if a single solution phase should be forced out of the system, then the number of
        ! moles of that phase alone will change significantly.  If the number of moles of more than 
        ! one solution phase changes, then the system may need to be further dampened.

        dTemp = 0.05D0
        
        ! Count the number of solution phases (i.e., j) that are changing by at least 5%:
        call CheckStagnation(dTemp,dMaxChange,j)

        ! If there is more than one solution phase that is changing by more than 5%, then dampen some more:
        if (j > 1) dStepLength = dStepLength * 0.5D0

    end if

    ! Further dampen the system if a miscible phase was recently added to the system:
    if ((iterGlobal - iterLast < 10).AND.(nMisciblePhases > 0)) dStepLength = dStepLength * 0.5D0

    ! If the number of moles of a solution phase is changing by too large of an amount, then recompute
    if (dStepLength < 1D0) call UpdateSystemVariables(dStepLength,dMolesSpeciesLast,dElementPotentialLast)

    ! Compute the chemical potentials of solution species:
    call CompChemicalPotential(lCompEverything)
    
    ! Compute the functional norm:
    call CompFunctionNorm

end subroutine InitGEMLineSearch


!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    ! 
    ! Purpose:
    ! ========
    ! 
    ! The purpose of this subroutine is to update the system variables using
    ! a specified steplength.
    !
    !---------------------------------------------------------------------------


subroutine UpdateSystemVariables(dStepLength,dMolesSpeciesLast,dElementPotentialLast)

    USE ModuleThermo
    USE ModuleGEMSolver

    integer                       :: i, j, k
    real(8)                       :: dTemp, dStepLength
    real(8), dimension(nSpecies)  :: dMolesSpeciesLast
    real(8), dimension(nElements) :: dElementPotentialLast


    ! Loop through all solution phases expected to be stable:
    do j = 1, nSolnPhases
        k     = -iAssemblage(nElements - j + 1)
        dTemp = 0D0
        do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
            dMolesSpecies(i) = dStepLength * dMolesSpecies(i) + (1D0 - dStepLength) * dMolesSpeciesLast(i)
            dMolesSpecies(i) = DMAX1(dMolesSpecies(i), dTolerance(8))
            dTemp            = dTemp + dMolesSpecies(i)
        end do
        ! Compute the total number of moles of each solution phase:
        dMolesPhase(nElements-j+1) = dTemp
    end do
    
    ! Dampen the number of moles of pure condensed phases:
    do i = 1, nConPhases
        dMolesPhase(i) = dStepLength * dMolesPhase(i) + (1D0 - dStepLength) * dMolesPhaseLast(i)
    end do
    
    ! Dampen the element potentials:
    LOOP_Gamma: do i = 1, nElements
        if (dElementPotential(i) == 0D0) then
            dElementPotential(i) = dElementPotentialLast(i)
            cycle LOOP_Gamma
        end if
        dElementPotential(i) = dStepLength * dElementPotential(i) + (1D0 - dStepLength) * dElementPotentialLast(i)
    end do LOOP_Gamma

end subroutine UpdateSystemVariables


    !---------------------------------------------------------------------------
    !                            END - GEMLineSearch.f90
    !---------------------------------------------------------------------------

