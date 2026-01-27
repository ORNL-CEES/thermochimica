!-------------------------------------------------------------------------------------------------------------
!
!> \file    test_error_handling.F90
!> \brief   Test suite for error handling in Thermochimica
!> \details Uses test-drive framework to validate error detection and reporting
!
!-------------------------------------------------------------------------------------------------------------

module test_error_handling
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use ModuleThermoIO
    implicit none
    private

    public :: collect_error_handling

contains

    !> Collect all error handling tests
    subroutine collect_error_handling(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("no_data_file", test_no_data_file), &
            new_unittest("nonexistent_data_file", test_nonexistent_data_file), &
            new_unittest("no_input_units", test_no_input_units), &
            new_unittest("no_temperature", test_no_temperature), &
            new_unittest("no_pressure", test_no_pressure), &
            new_unittest("no_mass", test_no_mass), &
            new_unittest("temperature_out_of_range", test_temperature_out_of_range), &
            new_unittest("pressure_out_of_range", test_pressure_out_of_range), &
            new_unittest("mass_out_of_range", test_mass_out_of_range), &
            new_unittest("pressure_nan", test_pressure_nan), &
            new_unittest("temperature_infinite", test_temperature_infinite), &
            new_unittest("all_elements_zero", test_all_elements_zero), &
            new_unittest("unsupported_mass_unit", test_unsupported_mass_unit) &
            ]
    end subroutine collect_error_handling

    !> Test that Thermochimica detects when no data file is specified
    !> Converted from TestThermo01.F90
    subroutine test_no_data_file(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables:
        dTemperature            = 300D0
        dPressure               = 1D0
        dElementMass            = 1D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        ! Intentionally not setting cThermoFileName

        ! Parse the ChemSage data-file (should fail):
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 6, &
            "Expected error code 6 (no data file), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica (full cleanup including parser):
        call ResetThermoAll

    end subroutine test_no_data_file

    !> Test that Thermochimica detects when mass is out of range
    !> Converted from TestThermo09.F90
    subroutine test_mass_out_of_range(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables:
        dTemperature            = 300D0
        dPressure               = 1D0
        dElementMass            = -1D0  ! Invalid negative mass
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 3, &
            "Expected error code 3 (mass out of range), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica (full cleanup including parser):
        call ResetThermoAll

    end subroutine test_mass_out_of_range

    !> Test that Thermochimica detects when a non-existent data file is specified
    !> Converted from TestThermo02.F90
    subroutine test_nonexistent_data_file(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables:
        dTemperature            = 300D0
        dPressure               = 1D0
        dElementMass            = 1D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'dsts/C-O.dat'  ! Typo in path

        ! Parse the ChemSage data-file (should fail):
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 6, &
            "Expected error code 6 (file not found), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_nonexistent_data_file

    !> Test that Thermochimica detects when input units are not specified
    !> Converted from TestThermo03.F90
    subroutine test_no_input_units(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables:
        dTemperature            = 300D0
        dPressure               = 1D0
        dElementMass            = 1D0
        ! Explicitly set unit strings to empty to simulate "not set"
        cInputUnitTemperature   = ''
        cInputUnitPressure      = ''
        cInputUnitMass          = ''
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 4, &
            "Expected error code 4 (no input units), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_no_input_units

    !> Test that Thermochimica detects when temperature is not specified
    !> Converted from TestThermo04.F90
    subroutine test_no_temperature(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables - set temperature to 0 (triggers "not specified"):
        dTemperature            = 0D0
        dPressure               = 1D0
        dElementMass            = 1D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 1, &
            "Expected error code 1 (temperature not specified), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_no_temperature

    !> Test that Thermochimica detects when pressure is not specified
    !> Converted from TestThermo05.F90
    subroutine test_no_pressure(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables - set pressure to 0 (triggers "not specified"):
        dTemperature            = 300D0
        dPressure               = 0D0
        dElementMass            = 1D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 2, &
            "Expected error code 2 (pressure not specified), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_no_pressure

    !> Test that Thermochimica detects when mass is not specified
    !> Converted from TestThermo06.F90
    subroutine test_no_mass(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables - set all masses to 0 (triggers "not specified"):
        dTemperature            = 2000D0
        dPressure               = 1D0
        dElementMass            = 0D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 5, &
            "Expected error code 5 (mass not specified), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_no_mass

    !> Test that Thermochimica detects when temperature is out of range (negative)
    !> Converted from TestThermo07.F90
    subroutine test_temperature_out_of_range(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables:
        dTemperature            = -300D0  ! Invalid negative temperature
        dPressure               = 1D0
        dElementMass            = 1D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 1, &
            "Expected error code 1 (temperature out of range), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_temperature_out_of_range

    !> Test that Thermochimica detects when pressure is out of range (negative)
    !> Converted from TestThermo08.F90
    subroutine test_pressure_out_of_range(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables:
        dTemperature            = 300D0
        dPressure               = -1D0  ! Invalid negative pressure
        dElementMass            = 1D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 2, &
            "Expected error code 2 (pressure out of range), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_pressure_out_of_range

    !> Test that Thermochimica detects when pressure is NaN
    !> Converted from TestThermo10.F90
    subroutine test_pressure_nan(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables:
        dTemperature            = 300D0
        dPressure               = -1D0
        dPressure               = sqrt(dPressure)  ! Creates NaN
        dElementMass            = 1D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 2, &
            "Expected error code 2 (pressure NaN), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_pressure_nan

    !> Test that Thermochimica detects when temperature is infinite
    !> Converted from TestThermo11.F90
    subroutine test_temperature_infinite(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables:
        dTemperature            = 1E37
        dTemperature            = dTemperature * dTemperature  ! Creates INF
        dPressure               = 0D0
        dElementMass            = 1D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 1, &
            "Expected error code 1 (temperature infinite), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_temperature_infinite

    !> Test that Thermochimica detects when all element masses are zero
    !> Converted from TestThermo12.F90
    subroutine test_all_elements_zero(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables:
        dTemperature            = 300D0
        dPressure               = 1D0
        dElementMass            = 0D0  ! All elements zero
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 5, &
            "Expected error code 5 (no elements specified), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_all_elements_zero

    !> Test that Thermochimica detects when an unsupported mass unit is specified
    !> Converted from TestThermo13.F90
    subroutine test_unsupported_mass_unit(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize variables:
        dTemperature            = 300D0
        dPressure               = 1D0
        dElementMass            = 1D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'elephants'  ! Unsupported unit
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! Check that the correct error code was reported
        call check(error, INFOThermo == 4, &
            "Expected error code 4 (unsupported mass unit), got: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_unsupported_mass_unit

    !-------------------------------------------------------------------------------------------------------------
    ! NOTE: Tests must use ResetThermoAll instead of ResetThermo to fully clean up between tests.
    ! ResetThermoAll deallocates both calculation state AND parser state, preventing allocation errors
    ! when running multiple tests in sequence with different data files.
    !-------------------------------------------------------------------------------------------------------------

    !> Helper function to convert integer to string
    function int_to_str(i) result(str)
        integer, intent(in) :: i
        character(len=20) :: str
        write(str, '(I0)') i
    end function int_to_str

end module test_error_handling
