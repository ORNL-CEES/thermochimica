
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ResetThermoAll.f90
    !> \brief   Deallocate all allocatable variables.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      ResetThermo.f90
    !> \sa      ResetThermoParser.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer         Description of change
    !    ----          ----------         ---------------------
    !    02/17/2012    M.H.A. Piro        Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to fully destruct all allocatable variables associated with
    !! Thermochimica.  This includes variables used internally by Thermochimica and by the parser that parses
    !! data-files as input.
    !
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine ResetThermoAll

    USE ModulePhaseConstraints

    ! Reset thermochimica:
    call ResetThermo

    ! Reset the parser used by thermochimica:
    call ResetThermoParser

    ! Reset reinit data:
    call ResetReinit

    ! Reset CTZ data:
    call ResetCTZ

    ! Clear phase constraints:
    call ClearPhaseConstraints

    return

end subroutine ResetThermoAll
