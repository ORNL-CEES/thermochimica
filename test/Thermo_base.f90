

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
    !    11/10/2016    S.Simunovic       Changed input model, new test
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
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = '../data/CEF4.dat'

    ! Specify values:
    dTemperature            = 1023D0
    dPressure               = 3D0
    dElementMass            = 0D0
    dElementMass(8)     = 8.4030D+03
    dElementMass(92)    = 3.6820D+03
    dElementMass(54)    = 1.3150D+02
    dElementMass(42)    = 1.0080D+02
    dElementMass(40)    = 9.8930D+01
    dElementMass(44)    = 8.9420D+01
    dElementMass(60)    = 8.0156D+01
    dElementMass(46)    = 6.7630D+01
    dElementMass(94)    = 6.5847D+01
    dElementMass(56)    = 5.7420D+01
    dElementMass(55)    = 5.7380D+01
    dElementMass(58)    = 5.4720D+01
    dElementMass(64)    = 3.0632D+01
    dElementMass(57)    = 2.5200D+01
    dElementMass(59)    = 2.1986D+01
    dElementMass(43)    = 1.8530D+01
    dElementMass(52)    = 1.2290D+01
    dElementMass(39)    = 1.1550D+01
    dElementMass(36)    = 1.0110D+01
    dElementMass(37)    = 9.3000D+00
    dElementMass(45)    = 7.5310D+00
    dElementMass(48)    = 7.3460D+00
    dElementMass(53)    = 5.3120D+00
    dElementMass(2)     = 3.4670D+00
    ! dElementMass(8)         = 0.6D0
    ! dElementMass(90)        = 0.4D0

    ! Specify output and debug modes:
    iPrintResultsMode = 2
    lDebugMode        = .FALSE.
    !lDebugMode        = .TRUE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Perform post-processing of results:
    if (iPrintResultsMode > 0) call PrintResults

    ! Destruct everything:
    if (INFOThermo == 0) call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

    

    
end program thermo
