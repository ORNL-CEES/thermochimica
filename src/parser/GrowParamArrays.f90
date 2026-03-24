
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GrowParamArrays.f90
    !> \brief   Dynamically grow the excess parameter arrays used by the ChemSage parser.
    !> \author  P. Bajpai
    !
    ! Purpose:
    ! ========
    !
    !> \details The ChemSage format does not encode the total number of excess mixing parameters in the
    !! file header, so the parser cannot pre-allocate exact sizes.  These subroutines double the first
    !! dimension of the relevant arrays whenever the current capacity is exceeded, using Fortran's
    !! move_alloc. They must be called before any write to
    !! iRegularParamCS(nParamCS+1,...) or iMagneticParamCS(nMagParamCS+1,...).
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine GrowRegularParamArrays

    USE ModuleParseCS

    implicit none

    integer :: oldSize, newSize
    integer,   allocatable, dimension(:,:) :: iTmp
    real(8),   allocatable, dimension(:,:) :: dTmp
    character, allocatable, dimension(:)   :: cTmp

    oldSize = size(iRegularParamCS, 1)
    newSize = oldSize * 2

    allocate(iTmp(newSize, size(iRegularParamCS, 2)))
    iTmp = 0
    iTmp(1:oldSize,:) = iRegularParamCS
    call move_alloc(iTmp, iRegularParamCS)

    allocate(dTmp(newSize, size(dRegularParamCS, 2)))
    dTmp = 0D0
    dTmp(1:oldSize,:) = dRegularParamCS
    call move_alloc(dTmp, dRegularParamCS)

    allocate(cTmp(newSize))
    cTmp = ' '
    cTmp(1:oldSize) = cRegularParamCS
    call move_alloc(cTmp, cRegularParamCS)

    return

end subroutine GrowRegularParamArrays


subroutine GrowMagneticParamArrays

    USE ModuleParseCS

    implicit none

    integer :: oldSize, newSize
    integer, allocatable, dimension(:,:) :: iTmp
    real(8), allocatable, dimension(:,:) :: dTmp

    oldSize = size(iMagneticParamCS, 1)
    newSize = oldSize * 2

    allocate(iTmp(newSize, size(iMagneticParamCS, 2)))
    iTmp = 0
    iTmp(1:oldSize,:) = iMagneticParamCS
    call move_alloc(iTmp, iMagneticParamCS)

    allocate(dTmp(newSize, size(dMagneticParamCS, 2)))
    dTmp = 0D0
    dTmp(1:oldSize,:) = dMagneticParamCS
    call move_alloc(dTmp, dMagneticParamCS)

    return

end subroutine GrowMagneticParamArrays
