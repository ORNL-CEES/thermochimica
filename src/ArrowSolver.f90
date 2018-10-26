
!-------------------------------------------------------------------------------
!
!> \file    ArrowSolver.f90
!> \brief   Solve a system of simultaneous linear equations with a 
!!           symmetric arrow matrix.
!> \author  M.H.A. Piro
!> \date    September 19, 2012
!
!
! Revisions:
! ==========
! 
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   09/19/2012      M.H.A. Piro         Original code
!   03/20/2013      M.H.A. Piro         Allow for multiple rows of the arrow
!                                        head.
!
!
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to solve a real system of linear 
!! equations Ax = b where A is an n-by-n symmetric arrow matrix and x and b are 
!! N vectors.  An arrow matrix is sparse everywhere except for the diagonal
!! (i.e., A(i,i) /= 0 except A(n,n) == 0), and the bottom row/extreme right column.
!! In other words, this matrix resembles the shape of an arrow head.  The arrow head
!! can be represented by a single row/column pair or by multiple rows/columns.  This
!! is defined by m.
!!
!! This subroutine effectively performs Gaussian elimination on only the bottom 
!! m rows.  Computational expense is further reduced for this specific problem when
!! performing back subsituation by exploiting the sparsity of the A matrix.  
!! This subroutine returns an error via the INFO variable if an error has been 
!! detected.  The possible values for INFO are summarized below:
!!
!! INFO        Description
!! ----        -----------
!!   0         Successful exit.
!!  -1         m is equal to or less than zero.
!!  -2         n is equal to or less than zero.
!!  -3         m is greater than or equal to n.
!!  -10        The system of linear equations contains an inconsistent row.
!!   i         The row corresponding to i (where i > 0) in B contains a NAN.
!!
!
!
!
! Pertinent variables:
! ====================
!
!> \param[in]     m         An integer scalar representing the number of rows
!!                            corresponding to the arrow head.
!> \param[in]     n         An integer scalar representing the number of linear
!!                            equations.
!> \param[in]     A         A double real array represending the A matrix.
!> \param[in,out] B         A double real vector (dimension n)   representing 
!!                            the functional vector on input and on output it 
!!                            represents the solution vector.
!> \param[out]    INFO      An integer scalar indicating a successful exit or 
!!                            an error.
!!
! Local variables:
! ----------------
! i, j, k                   Integer scalars used for indexing.
! dMultiplier               A double real scalar used for calculating the 
!                            multiplier in Gaussian elimination.
!
!-------------------------------------------------------------------------------

    
subroutine ArrowSolver(m, n, A, B, INFO)

    implicit none
    
    integer                                :: i, j, k
    integer,                 intent(out)   :: INFO
    integer,                 intent(inout) :: m, n
    real(8)                                :: dMultiplier
    real(8), dimension(n,n), intent(inout) :: A
    real(8), dimension(n),   intent(inout) :: B
    
    
    ! Initialize variables:
    INFO = 0

    ! Ensure that input parameters are within bounds:
    if (m <= 0) then
        INFO = -1
        return
    elseif (n <= 0) then
        INFO = -2
        return
    elseif (m >= n) then
        INFO = -3
        return
    end if
    
    ! Loop through number of rows representing the arrow head:
    LOOP_ArrowHead: do k = n - m + 1, n
    
        ! Loop through columns along the diagonal:
        LOOP_Column: do j = 1, k - 1
        
            ! Cycle if the coefficient is zero:
            if (A(k,j) == 0D0) cycle LOOP_Column
            
            ! Compute the multiplier:
            dMultiplier = A(j,j) / A(k,j)
            
            ! Update all columns along row:
            A(k,j) = 0D0
            do i = j + 1, n
                A(k,i) = dMultiplier * A(k,i) - A(j,i)
            end do
            
            ! Update solution vector component:
            B(k) = dMultiplier * B(k) - B(j)

        end do LOOP_Column
        
    end do LOOP_ArrowHead

    ! Check if the system of equations contains an inconsistent row:
    IF_NULL: if (A(n,n) /= 0D0) then
    
        ! Perform back-step calculation:
        B(n) = B(n) / A(n,n)

        ! Loop backwards:
        LOOP_Back: do i = (n - 1), (1), (-1)
    
            ! Reinitialize temporary variable:
            dMultiplier = B(i)
        
            ! Loop along column:
            LOOP_Inner: do j = n, (n - m), (-1)
                ! Exit if the diagonal has been reached:
                if (j <= i) exit LOOP_Inner
                dMultiplier = dMultiplier - A(i,j) * B(j)
            end do LOOP_Inner
            
            ! Update the solution vector component:
            B(i) = dMultiplier / A(i,i)
        
        end do LOOP_Back
    else
    
        ! The system contains an inconsistent row.  Return an error:
        INFO = -10
        
    end if IF_NULL
    
    ! Check for a NAN:
    LOOP_NAN: do i = 1, n
    
        if (B(i) /= B(i)) then
            ! A NAN has been detected; report an error corresponding
            ! to the row of the first NAN and reset the solution vector.
            INFO = i 
            B    = 0D0
            exit LOOP_NAN
        end if
        
    end do LOOP_NAN
    
    return
      
end subroutine ArrowSolver
