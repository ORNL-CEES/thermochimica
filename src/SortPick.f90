
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    SortPick.f90
    !> \brief   Sort a double real vector (this vector is unchanged, the indices of this vector is sorted).
    !> \author  W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery (modified by M.H.A. Piro)
    !> \sa      CheckSolnPhaseAdd.f90
    !
    !
    ! References:
    ! ===========
    !
    ! This subroutine employs a simple sorting algorithm that was modified from:
    !
    !   W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery, "Numerical Recipes in Fortran 90 
    !   (Second Edition)", Cambridge University Press, New York, 1996
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Modified the sort_pick subroutine from Numerical Recipes to 
    !                                        return an integer vector representing the descending order of 
    !                                        the real vector (unchanged).
    !   02/06/2013      M.H.A. Piro         Only proceed if n is positive.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to take a double real vector as input and return an integer vector 
    !! representing the descending order of the real vector (unchanged).  Although this sorting routine is by
    !! no means efficient for large problems, "it is meant to be invoked only for the most trivial sorting jobs,
    !! say, N < 20." 
    !
    !
    ! Pertinent variables:
    ! ====================
    !   
    !> \param[in]   n           An integer scalar represening the dimension of dVec.
    !> \param[in]   dVec        A double real vector that is to be sorted in descending order.
    !> \param[out]  iVec        An integer vector representing the index of coefficients of dVec in
    !!                           descending order.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine SortPick(n,dVec,iVec)

    implicit none   

    integer                            :: i, j, k
    integer,               intent(in)  :: n 
    integer, dimension(n), intent(out) :: iVec
    real(8), dimension(n), intent(in)  :: dVec
    
    
    ! Only proceed if n is greater than zero:
    if (n > 0) then
        ! Initialize variables:
        iVec = 0
        do i = 1, n
            iVec(i) = i
        end do
    
        ! Sort the indices on dVec in the iVec vector:
        LOOP_A: do j = 2, n
            LOOP_B: do i = j-1, 1, -1
                if (dVec(iVec(i)) >= dVec(j)) exit LOOP_B
                k         = iVec(i)
                iVec(i)   = j
                iVec(i+1) = k
            end do LOOP_B 
        end do LOOP_A
    
    end if
    
    return
    
end subroutine SortPick
