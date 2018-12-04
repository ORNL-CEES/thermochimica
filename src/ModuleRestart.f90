
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file        ModuleRestart.f90
    !> \brief       Fortran module for internal use of Thermochimica
    !> \details     The purpose of this module is to provide the means to share information used for restarting
    !! calculations.
    !> \author      M. Poschmann
    !
    !-------------------------------------------------------------------------------------------------------------


module ModuleRestart

    implicit none

    SAVE

    integer,       dimension(:),   allocatable::  iAssemblage_Old
    real(8),       dimension(:),   allocatable::  dChemicalPotential_Old, dMolesPhase_Old, dElementPotential_Old
    real(8),       dimension(:),   allocatable::  dMolFraction_Old

end module ModuleRestart
