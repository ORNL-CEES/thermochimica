
    !-------------------------------------------------------------------------------------------------------------
    !
    ! DISCLAIMER
    ! ===========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    !
    ! Revisions:
    ! ===========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    28/11/2018    M. Poschmann        Creation
    !
    ! Purpose:
    ! =========
    ! Test the new ParseInput subroutine.
    !
    !-------------------------------------------------------------------------------------------------------------

program ParserTest

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none
    integer :: i
    real    :: start, finish

    ! Call input parser
    call ParseInput('trial-input.ti')

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call cpu_time(start)
    LOOP_time: do i = 1,1
      call Thermochimica
    end do LOOP_time
    call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start

    ! Perform post-processing of results:
    if (iPrintResultsMode > 0)  call PrintResults

    ! Reset Thermochimica:
    call ResetThermo

    ! Call the debugger:
    call ThermoDebug

end program ParserTest
