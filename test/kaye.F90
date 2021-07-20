

    !-------------------------------------------------------------------------------------------------------------
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified.
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer        Description of change
    !    ----          ----------        ---------------------
    !    10/05/2011    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this program is to provide an example Application Programming Interface (API) for
    ! Thermochimica to simulate the interaction with another code.  The module "ModuleThermoIO" is intended to
    ! be used by another code that calls Thermochimica for input/output operations.
    !
    !   Subroutine                  Brief Description
    !   ----------                  -----------------
    !
    !   ParseCSDataFile             Parse ChemSage data-file.
    !   Thermochimica               Compute thermodynamic equilibrium with provided input.
    !   PrintResults                Print results to screen (optional).
    !   ResetThermoAll              Destructor: deallocate all allocatable arrays (optional).
    !   ThermoDebug                 Debugger: provide an error message if an error has been encountered (optional).
    !
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !   Variable                    Brief Description
    !   --------                    -----------------
    !
    !   cInputUnitTemperature       A character variable representing the input temperature units.
    !   cInputUnitPressure          A character variable representing the input hydrostatic pressure units.
    !   cInputUnitMass              A character variable representing the input mass units.
    !   cThermoFileName             A character variable representing the path and file name of the data-file.
    !   dTemperature                A double real scalar representing temperature.
    !   dPressure                   A double real scalar representing hydrostatic pressure.
    !   dElementMass                A double real vector (0:118) representing the mass of each chemical element.
    !                                Coefficients corespond to the atomic number on the periodic table (e.g.,
    !                                dElementMass(92) corresponds to uranium).
    !   iPrintResultsMode           An integer scalar indicating the print level.
    !                                0 = do not print output to screen.
    !                                1 = print simple details to screen on exit.
    !                                2 = print advanced details to screen on exit.
    !    INFOThermo                 An integer scalar indicating a successful calculation (==0) or an error (/=0).
    !                                See ThermoDebug.f90 for a full explanation of all error codes.
    !
    !-------------------------------------------------------------------------------------------------------------

program thermo

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo

    implicit none

    integer :: i, j, k

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'Kaye_NobleMetals.dat'

    ! Specify values:
    dPressure              = 1D0
    ! dTemperature           = 1723D0
    ! dElementMass(42)       = 0.3D0            ! Mo
    ! dElementMass(44)       = 0.5D0            ! Ru
    ! dElementMass(46)       = 0.2D0            ! Pd

    dTemperature           = 400D0
    dElementMass(42)       = 0.9D0            ! Mo
    dElementMass(44)       = 0.1D0            ! Ru
    ! dElementMass(46)       = 1.1D0            ! Pd

    ! Specify output mode:
    iPrintResultsMode     = 2

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    ! Perform post-processing of results:
    if (iPrintResultsMode > 0)  call PrintResults

    ! do i = 1, nSolnPhases
    !   k = -iAssemblage(nElements + 1 - i)
    !   do j = nSpeciesPhase(k-1)+1, nSpeciesPhase(k)
    !     print *, cSpeciesName(j), ' ', dMolFraction(j), ' ', dPartialExcessGibbs(j) * dIdealConstant * dTemperature
    !   end do
    ! end do

    ! print *, nDummySpecies, nChargedConstraints

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

    call RunKaye

end program thermo
