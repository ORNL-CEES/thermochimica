
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file        ModuleSubMin.f90
    !> \brief       Fortran module for the Subminimization routine.
    !> \details     The purpose of this module is to
    !> \author      M.H.A. Piro
    !
    !> \param       nVar            The number of species in the specified solution phase.
    !> \param       iFirst          The first species in the specified solution phase.
    !> \param       iLast           The last species in the specified solution phase.
    !> \param       iterSubMax      The maximum number of iterations in the subminimization routine.
    !> \param       dDrivingForce   The driving force of the specified solution phase.
    !> \param       dSubMinTolerance  A numerical tolerance of the maximum change in dMolFraction that is used to
    !!                                identify convergence.
    !> \param       dChemicalPotentialStar  The chemical potential of each solution phase constituent defined by
    !!                                  the element potentials.
    !> \param       dRHS            A working vector that is used to represent the functional vector in the
    !!                               SubMinNewton subroutine and the direction vector that is computed.
    !> \param       dHessian        A double real array representing the Hessian matrix.
    !> \param       lSubMinConverged  A logical scalar indicating a converged (i.e., TRUE) or non-converged
    !!                                solution (i.e., FALSE).
    !> \param       dTolEuclideanNorm A double real scalar representing the tolerance of the Euclidean norm
    !!                                  between the mole fraction vectors of two corresponding solution phases.
    !> \param       dTolDrivingForceChange A double real scalar represending the tolerance for the change in
    !!                                 the driving force.
    !
    !-------------------------------------------------------------------------------------------------------------


module ModuleSubMin

    implicit none

    SAVE

    integer::                               nVar, iFirst, iLast, iSolnPhaseIndexOther
    integer, dimension(:),   allocatable::  iHessian

    real(8)::                               dDrivingForce, dDrivingForceLast, dSubMinFunctionNorm
    real(8), parameter::                    dSubMinTolerance = 1D-3, dMinMoleFraction = 1D-100
    real(8), parameter::                    dTolEuclideanNorm = 1D-2, dTolDrivingForceChange = 1D-3
    real(8), dimension(:),   allocatable::  dChemicalPotentialStar, dRHS
    real(8), dimension(:,:), allocatable::  dHessian

    logical::                               lSubMinConverged

end module ModuleSubMin
