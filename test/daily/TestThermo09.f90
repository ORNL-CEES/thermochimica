
program ThermoTest09

    !-------------------------------------------------------------------------------------------------------------
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer        Description of change
    !    ----          ----------        ---------------------
    !    02/07/2012    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this unit test is to ensure that Thermochimica exits gracefully when the system cannot be
    ! represented by the selection of phases provided by the data-file.  This particular data-file only contains
    ! pure stoichiometric solid phases from the uranium-oxygen system, but without any gas or other phases.  An
    ! oxygen-rich system cannot be represented by this particular collection of phases.
    !
    !-------------------------------------------------------------------------------------------------------------


    USE ModuleThermoIO

    implicit none

    integer :: STATUS = 0

    ! Initialize variables:
    dTemperature            = 300D0
    dPressure               = 1D0
    dElementMass(8)         = 1D0       ! O
    dElementMass(92)        = 0.001D0   ! U
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/Example2.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if (INFOThermo == 14) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo09: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo09: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program ThermoTest09
