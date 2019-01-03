
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GEMSolver.f90
    !> \brief   Gibbs Energy Minimization solver.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      Thermochimica.f90
    !> \sa      InitGEMSolver.f90
    !> \sa      GEMNewton.f90
    !> \sa      GEMLineSearch.f90
    !> \sa      CheckPhaseAssemblage.f90
    !> \sa      CheckConvergence.f90
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the quantities of species and phases at thermodynamic
    !! equilibrium using the Gibbs Energy Minimization (GEM) method.  This subroutine uses values of
    !! dMolesPhase, dChemicalPotential and iAssemblage from the Leveling and PostLeveling subroutines as initial 
    !! estimates for computation.  
    !!
    !! The main subroutines used by this solver are summarized below:
    !! <table border="1" width="800">
    !! <tr>
    !!    <td> <b> File name </td> <td> Description </b> </td>
    !! </tr>
    !! <tr>
    !!    <td> InitGEMSolver.f90 </td> 
    !!    <td> Initialize the GEMSolver by establishing the initial phase assemblage and composition.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckSysOnlyPureConPhases.f90 </td> 
    !!    <td> Check the system if there are only pure condensed phases. The system may already be converged.</td>
    !! </tr>
    !! <tr>
    !!    <td> GEMNewton.f90 </td> 
    !!    <td> Compute the direction vector using Newton's method.  </td>
    !! </tr> 
    !! <tr>
    !!    <td> GEMLineSearch.f90 </td> 
    !!    <td> Perform a line search along the direction vector.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckPhaseAssemblage.f90 </td> 
    !!    <td> Check if the phase assemblage needs to be adjusted.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckConvergence.f90 </td> 
    !!    <td> Check if the system has converged.  </td>
    !! </tr>
    !! </table>
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
    ! INFO                  An integer scalar identifying an error from LAPACK.  This is used by the GEMNewton
    !                        subroutine to indicate whether there is a singularity in the Hessian matrix.
    ! lConverged            A logical variable indicating whether the code has convered (.TRUE.) or not (.FALSE.).
    ! lRevertSystem         A logical scalar indicating whether the system should be reverted to a previously
    !                        successful phase assemblage.
    ! dTolerance            A double real vector representing numerical tolerances (defined in InitThermo.f90).
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine GEMSolver

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    
    implicit none
    
    integer::   INFO


    ! Initialize the GEM solver:
    call InitGEMSolver

    !!!
    !!! CONSIDER MOVING THIS INTO THE InitGEMSolver SUBROUTINE:
    !!!
    ! The system may be converged if there aren't any solution phases:
    if ((nSolnPhases == 0).AND.(INFOThermo == 0)) then
    
        ! Check the system if only pure condensed phases are expected to appear:
        call CheckSysOnlyPureConPhases
    
        ! Report an error if this failed:
        if (lConverged .EQV. .FALSE.) INFOThermo = 14
        
    end if

    ! Begin the global iteration cycle:
    LOOP_GEMSolver: do iterGlobal = 1, iterGlobalMax

        ! If in debug mode, call the debugger:
        if (lDebugMode .EQV. .TRUE.) call GEMDebug(1) 
    
        ! Construct the Hessian matrix and compute the direction vector:
        call GEMNewton(INFO)

        ! Perform a line search using the direction vector: 
        call GEMLineSearch            

        ! Check if the estimated phase assemblage needs to be adjusted: 
        call CheckPhaseAssemblage

        ! Check convergence:
        if (iterGlobal /= iterLast) call CheckConvergence

        ! If in debug mode, call the debugger:
        if (lDebugMode .EQV. .TRUE.) call GEMDebug(9)

        ! Return control to the main program if an error has occured or if the solution has converged:
        if ((INFOThermo /= 0).OR.(lConverged .EQV. .TRUE.)) exit LOOP_GEMSolver

    end do LOOP_GEMSolver

    ! Report an error if the GEMSolver did not converge but no other errors were encountered:
    if ((lConverged .EQV. .FALSE.).AND.(INFOThermo == 0)) INFOThermo = 12

    return
 
end subroutine GEMSolver
