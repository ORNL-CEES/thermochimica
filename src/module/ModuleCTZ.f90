
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file        ModuleCTZ.f90
    !> \brief       Fortran module for internal use of Thermochimica
    !> \details     Contains data required for using Common Tangent Zone interpolation
    !! calculations.
    !> \author      M. Poschmann
    !
    !-------------------------------------------------------------------------------------------------------------


module ModuleCTZ

    implicit none

    SAVE

    integer :: nMaxAssemblages, nMaxElements
    integer, dimension(:,:), allocatable :: assemblageHistory
    character(12), dimension(:,:), allocatable :: elementHistory
    real(8), dimension(:,:), allocatable :: assemblageTlimits
    real(8), dimension(:,:,:,:), allocatable :: stoichHistory
    logical :: lCtzInit = .FALSE.

end module ModuleCTZ
