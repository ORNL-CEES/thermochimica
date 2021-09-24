
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ThermoDebug.f90
    !> \brief   Thermochimica debugger
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \param[in] INFOThermo An integer scalar representing a successful exit or an error.
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer       Description of change
    !    ----          ----------       ---------------------
    !    10/21/2011    M.H.A. Piro      Original code.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to print an error message corresponding to a particular value
    !! of INFOThermo if a problem is encountered.  This subroutine is intended to only be used for debugging
    !! purposes.
    !!
    !-------------------------------------------------------------------------------------------------------------


subroutine ThermoDEBUG

    USE ModuleThermoIO
    USE ModuleParseCS, ONLY: nSolnTypeSupport, cSolnPhaseTypeSupport

    implicit none

    integer::   i, j
    print *

    if (INFOThermo == 0) then
        print *, 'DEBUG: Successful exit.'
    elseif (INFOThermo == 1) then
        print *, 'DEBUG: Temperature is out of range or a NAN.'
        !Check the CheckThermoInput.f90 subroutine.
    elseif (INFOThermo == 2) then
        print *, 'DEBUG: Pressure is out of range or a NAN.'
        ! Check the CheckThermoInput.f90 subroutine.
    elseif (INFOThermo == 3) then
        print *, 'DEBUG: Mass is out of range or a NAN.'
        ! Check the CheckThermoInput.f90 subroutine.
    elseif (INFOThermo == 4) then
        print *, 'DEBUG: The input variable cThermoInputUnits is not recognizable.'
        ! Check the CheckThermoInput.f90 subroutine.
    elseif (INFOThermo == 5) then
        print *, 'DEBUG: The number of acceptable system components is out of range.'
    elseif (INFOThermo == 6) then
        print *, 'DEBUG: The specified data-file was not found.'
    elseif (INFOThermo == 7) then
        print *, 'DEBUG: There is an unknown error in reading the data-file.'
    elseif (INFOThermo == 8) then
        print *, 'DEBUG: The number of solution phases in the system exceeds the maximum allowable number.'
    elseif (INFOThermo == 9) then
        print *, 'DEBUG: A pure chemical species is needed to be present in the database for each element.'
        ! Check the CheckThermoData.f90 subroutine.
    elseif (INFOThermo == 10) then
        print *, 'DEBUG: The Leveling subroutine failed to determine an appropriate phase assemblage for'
        print *, 'further consideration.   '
    elseif (INFOThermo == 11) then
        print *, 'DEBUG: The PostLeveling subroutine failed.'
    elseif (INFOThermo == 12) then
        print *, 'DEBUG: The GEMSolver subroutine failed to converge.'
    elseif (INFOThermo == 13) then
        print *, 'DEBUG: The GEMSolver subroutine detected a NAN.'
    elseif (INFOThermo == 14) then
        print *, 'DEBUG: The GEMSolver determined that there are no solution phases, but the system'
        print *, 'cannot be represented by only pure condensed phases. '
        ! Check the InitGEM.f90 and CheckConvergence.f90 subroutines.
    elseif (INFOThermo == 15) then
        print *, 'DEBUG: Failed to deallocate allocatable variables from the ModuleThermo, ModuleThermoIO'
        print *, 'and/or ModuleGEMSolver modules.'
        ! Check the ResetThermo.f90 subroutine.
    elseif (INFOThermo == 16) then
        print *, 'DEBUG: The LAPACK driver routines were not able to invert the Jacobian matrix in the GEMSolver.'
        ! Check the GEMNewton.f90 subroutine.
    elseif (INFOThermo == 17) then
        print *, 'DEBUG: The data-file contains a solution phase type that is not currently supported by Thermochimica.'
        print *
        print *, 'The following solution phase types are currently supported:'
        do i = 1, nSolnTypeSupport
            print *, ' ', cSolnPhaseTypeSupport(i)
        end do
        print *
        ! Check the ParseCSDataBlock.f90 subroutine.
    elseif (INFOThermo == 18) then
        print *, 'DEBUG: Failed to deallocate allocatable array used in the ModuleParseCS module.'
        ! Check the ResetThermoParser subroutine.
    elseif (INFOThermo == 19) then
        print *, 'DEBUG: Failed to deallocate allocatable array used in the CheckSystem subroutine.'
        ! Check the CheckSystem subroutine.
    elseif (INFOThermo == 20) then
        print *, 'DEBUG: Failed to deallocate allocatable array used in the LevelingSolver subroutine.'
        ! Check the LevelingSolver subroutine.
    elseif (INFOThermo == 21) then
        print *, 'DEBUG: Failed to deallocate allocatable array used in the InitGEM subroutine.'
        ! Check the InitGEM subroutine.
    elseif (INFOThermo == 22) then
        print *, 'DEBUG: Failed to deallocate allocatable array used in the CompMolSolnPhase subroutine.'
        ! Check the CompMolSolnPhase subroutine.
    elseif (INFOThermo == 23) then
        print *, 'DEBUG: Failed to deallocate allocatable array used in the GEMBroyden subroutine.'
        ! Check the GEMBroyden subroutine.
    elseif (INFOThermo == 24) then
        print *, 'DEBUG: A NAN was detected in the CompStoichSolnPhase.f90 subroutine.'
    elseif (INFOThermo == 25) then
        print *, 'DEBUG: A NAN was detected in the CompMolFraction.f90 subroutine.'
    elseif (INFOThermo == 26) then
        print *, 'DEBUG: Failed to deallocate allocatable variables used in the CompMolAllSolnPhases subroutine.'
    elseif (INFOThermo == 27) then
        print *, 'DEBUG: The CheckQKTOSolnPhase subroutine failed to converge.'
    elseif (INFOThermo == 28) then
        print *, 'DEBUG: LAPACK returned an error in the SubMinNewton subroutine.'
    elseif (INFOThermo == 29) then
        print *, 'DEBUG: There is an unknown error in the Subminimization subroutine.'
    elseif (INFOThermo == 30) then
        print *, 'DEBUG: LAPACK reported an error in InitMolFractionMiscGap.'
    elseif (INFOThermo == 31) then
        print *, 'DEBUG: The name of a chemical element in the parsed data-file does not correspond to a known '
        print *, 'element on the periodic table.'
    elseif (INFOThermo == 32) then
        print *, 'DEBUG: An unsupported number of mixing parameters is included in the data-file. Thermochimica '
        print *, 'is currently capable of handling: binary, ternary and quaternary terms.'
    elseif (INFOThermo == 33) then
        print *, 'DEBUG: Parser: The number of sublattices per phase is not supported. '
    elseif (INFOThermo == 34) then
        print *, 'DEBUG: Parser: Magnetic mixing terms are not supported yet. '
    elseif (INFOThermo == 35) then
        print *, 'DEBUG: A coefficient in the iTempVec vector in the CheckPureConPhaseRem subroutine is zero,'
        print *, 'but must be positive.'
    elseif (INFOThermo == 36) then
        print *, 'DEBUG: An error occured in interpreting the parsed data for a SUBL phase.'
    elseif (INFOThermo == 40) then
        print *, 'DEBUG: There is an element in a compound that is not in the dat file.'
    elseif (INFOThermo == 41) then
        print *, 'DEBUG: Error finding stoichiometry in terms of compounds.'
    elseif (INFOThermo == 42) then
        print *, 'DEBUG: Excess mixing term in SUBG/SUBQ phase not supported.'
    elseif (INFOThermo == 43) then
        print *, 'DEBUG: An unsupported type of magnetic mixing parameter is included in the data-file.'
        print *, ' Thermochimica is currently capable of handling: binary terms.'
    elseif (INFOThermo == 44) then
        print *, 'DEBUG: An unsupported type of excess mixing term in the SUBI phase is included in the data-file.'
        print *, 'Please contact the developers and share this data-file.'
    elseif (INFOThermo == 99) then
        print *, 'DEBUG: The input element masses are not representable in terms of the available species.'
        ! Check CompThermoData.f90
    elseif ((INFOThermo >= 100).AND.(INFOThermo < 1000)) then
        i = INFOThermo - 100
        print '(A31,I1,A41)', ' DEBUG: Error in reading "line ', i, '" of the header section of the data-file.'
    elseif ((INFOThermo >= 1000).AND.(INFOThermo < 10000)) then
        i = (INFOThermo - 1000) / 100
        i = INT(i)
        j = (INFOThermo - 1000) - i * 100
        print '(A32,I2,A20,I2,A44)', ' DEBUG: Error in reading "entry ', i, '" of solution phase ', j, &
            ' of the data-block section of the data-file.'
    elseif ((INFOThermo >= 10000).AND.(INFOThermo < 1000000)) then
        i = (INFOThermo - 10000) / 1000
        i = INT(i)
        j = (INFOThermo - 10000) - i * 1000
        print '(A44,I2,A19,I2,A44)', ' DEBUG: Error in reading excess mixing term ', i, ' of solution phase ', j, &
            ' of the data-block section of the data-file.'
    else
        print *, 'DEBUG: An unknown error has occured. Error code ', INFOThermo
    end if

    print *

    return

end subroutine ThermoDEBUG
