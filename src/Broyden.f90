
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    Broyden.f90
    !> \brief   Perform a unit iteration using Broyden's method.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !    Date          Programmer         Description of change
    !    ----          ----------         ---------------------
    !    11/11/2011    M.H.A. Piro        Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to perform a unit iteration in solving a system of 
    !! non-linear equations using the "good" Broyden method.  This subroutine does not perform an iteration 
    !! process, but is intended to be used within a main iteration cycle.  An estimate of the inverse of the 
    !! Broyden matrix is required as input and an improved estimate is provided as output in addition to a 
    !! new direction vector s.  
    !!
    !! One of the advantages of using Broyden's method in comparison to the more popular BFGS method is 
    !! that Broyden's method permits nonsymmetric updates to the inverse Broyden matrix, whereas the BFGS 
    !! method maintains symmetry.  The true Jacobian for some numerical problems may not always be symmetric.
    !!
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]     nVar          The number of unknown variables (and of course the number of non-linear 
    !!                               equations)
    !> \param[in,out] s             The change of the unknown variable vector (nVar)
    !> \param[in,out] y             The change of the function vector (nVar)
    !> \param[in]     f             The function vector (nVar)
    !> \param[in,out] dInvBroyden   The inverse of the Broyden matrix (nVar x nVar)    
    !
    !-------------------------------------------------------------------------------------------------------------

    
subroutine Broyden(nVar,y,s,f,dInvBroyden)
    
    implicit none
    
    integer::                       i, j, nVar   
    real(8)::                       dTemp 
    real(8),dimension(nVar)::       f, y, s, dTempVec
    real(8),dimension(nVar,nVar)::  dInvBroyden

   
    ! Compute the difference vector multiplied by the current inverse Broyden matrix:
    do i = 1, nVar
        dTempVec(i) = 0D0
        do j = 1, nVar
            dTempVec(i) = dTempVec(i) + dInvBroyden(i,j) * y(j)
        end do
    end do
    
    ! Calculate the dot products for denominators:
    dTemp = 0D0
    do i = 1, nVar
        dTemp = dTemp + s(i) * dTempVec(i)
    end do
    
    ! Make the denominator multiplicative:
    dTemp = 1D0 / dTemp

    ! Compute the work vector:
    do i = 1, nVar
        y(i) = 0D0
        do j = 1, nVar
            y(i) = y(i) + s(j) * dInvBroyden(j,i)
        end do
    end do
    
    ! Update the inverse Broyden matrix:
    do i = 1, nVar
        do j = 1, nVar
            dInvBroyden(i,j) = dInvBroyden(i,j) + dTemp * y(j) * (s(i) - dTempVec(i))
        end do
    end do
    
    ! Compute the direction vector:
    do i = 1, nVar
        s(i) = 0D0
        do j = 1, nVar
            s(i) = s(i) - dInvBroyden(i,j) * f(j)
        end do
    end do
    
    return
    
end subroutine Broyden