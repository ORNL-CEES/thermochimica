!-------------------------------------------------------------------------------------------------------------
!
!> \file    test_systems.F90
!> \brief   Test suite for system-specific calculations in Thermochimica
!> \details Uses test-drive framework to validate thermodynamic calculations for various chemical systems
!
!-------------------------------------------------------------------------------------------------------------

module test_systems
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use ModuleThermoIO
    use ModuleThermo
    use ModuleGEMSolver
    implicit none
    private

    public :: collect_systems

contains

    !> Test that Thermochimica can parse a data file with many solution phases
    !> Converted from TestThermo14.F90
    subroutine test_maximum_solution_phases(error)
        type(error_type), allocatable, intent(out) :: error

        ! Initialize element masses
        dElementMass = 0D0

        ! Specify data file with 42 solution phases
        cThermoFileName = DATA_DIRECTORY // 'PdRuTcMo.dat'

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Check that parsing succeeded
        call check(error, INFOThermo == 0, &
            "Failed to parse file with many solution phases, error: " // trim(adjustl(int_to_str(INFOThermo))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_maximum_solution_phases

    !> Collect all system validation tests
    subroutine collect_systems(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("maximum_solution_phases", test_maximum_solution_phases), &
            new_unittest("W_Au_Ar_O_system", test_w_au_ar_o_system), &
            new_unittest("W_Au_Ar_O_system_02", test_w_au_ar_o_system_02), &
            new_unittest("W_Au_Ar_Ne_O_high_temp", test_w_au_ar_ne_o_high_temp), &
            new_unittest("W_Au_Ar_Ne_O_low_temp", test_w_au_ar_ne_o_low_temp) &
            ]
    end subroutine collect_systems

    !> Test W-Au-Ar-O system equilibrium calculation
    !> Converted from TestThermo30.F90
    subroutine test_w_au_ar_o_system(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        ! Initialize element masses to zero (important for test suite sequential execution)
        dElementMass           = 0D0

        ! Specify units:
        cInputUnitTemperature  = 'K'
        cInputUnitPressure     = 'atm'
        cInputUnitMass         = 'moles'
        cThermoFileName        = DATA_DIRECTORY // 'WAuArO-1.dat'

        ! Specify values:
        dPressure              = 1D0
        dTemperature           = 1455D0
        dElementMass(74)       = 1.95D0        ! W
        dElementMass(79)       = 1D0           ! Au
        dElementMass(18)       = 2D0           ! Ar
        dElementMass(8)        = 10D0          ! O

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! First check: ensure calculation succeeded
        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermo
            return
        end if

        ! Second check: validate Gibbs energy
        expected_gibbs = -4.620D5
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))) // &
            ", Relative error: " // trim(adjustl(real_to_str(relative_error))))

        ! Reset Thermochimica (full cleanup including parser):
        call ResetThermoAll

    end subroutine test_w_au_ar_o_system

    !> Test W-Au-Ar-O system 02 equilibrium calculation
    !> Converted from TestThermo31.F90
    subroutine test_w_au_ar_o_system_02(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        ! Initialize element masses
        dElementMass = 0D0

        ! Specify units:
        cInputUnitTemperature  = 'K'
        cInputUnitPressure     = 'atm'
        cInputUnitMass         = 'moles'
        cThermoFileName        = DATA_DIRECTORY // 'WAuArO-2.dat'

        ! Specify values:
        dPressure              = 1D0
        dTemperature           = 1000D0
        dElementMass(74)       = 1D0        ! W
        dElementMass(79)       = 3D0        ! Au
        dElementMass(18)       = 5D0        ! Ar
        dElementMass(8)        = 2D0        ! O

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! First check: ensure calculation succeeded
        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ! Second check: validate Gibbs energy
        expected_gibbs = 6.769D5
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))) // &
            ", Relative error: " // trim(adjustl(real_to_str(relative_error))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_w_au_ar_o_system_02

    !> Test W-Au-Ar-Ne-O system at high temperature (2452K)
    !> Converted from TestThermo32.F90
    subroutine test_w_au_ar_ne_o_high_temp(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        ! Initialize element masses
        dElementMass = 0D0

        ! Specify units:
        cInputUnitTemperature  = 'K'
        cInputUnitPressure     = 'atm'
        cInputUnitMass         = 'moles'
        cThermoFileName        = DATA_DIRECTORY // 'WAuArNeO-1.dat'

        ! Specify values:
        dPressure              = 1D0
        dTemperature           = 2452D0
        dElementMass(74)       = 1.95D0     ! W
        dElementMass(79)       = 1D0        ! Au
        dElementMass(18)       = 2D0        ! Ar
        dElementMass(8)        = 10D0       ! O
        dElementMass(10)       = 10D0       ! Ne

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! First check: ensure calculation succeeded
        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ! Second check: validate Gibbs energy
        expected_gibbs = 1.672D7
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))) // &
            ", Relative error: " // trim(adjustl(real_to_str(relative_error))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_w_au_ar_ne_o_high_temp

    !> Test W-Au-Ar-Ne-O system at low temperature (900K)
    !> Converted from TestThermo33.F90
    subroutine test_w_au_ar_ne_o_low_temp(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_molfrac1, computed_molfrac1, rel_error_mf1
        real(8) :: expected_molfrac9, computed_molfrac9, rel_error_mf9

        ! Initialize element masses
        dElementMass = 0D0

        ! Specify units:
        cInputUnitTemperature  = 'K'
        cInputUnitPressure     = 'atm'
        cInputUnitMass         = 'moles'
        cThermoFileName        = DATA_DIRECTORY // 'WAuArNeO-2.dat'

        ! Specify values:
        dPressure              = 2D0
        dTemperature           = 900D0
        dElementMass(74)       = 20D0       ! W
        dElementMass(79)       = 2D0        ! Au
        dElementMass(18)       = 7D0        ! Ar
        dElementMass(8)        = 5D0        ! O
        dElementMass(10)       = 1D0        ! Ne

        ! Parse the ChemSage data-file:
        call ParseCSDataFile(cThermoFileName)

        ! Call Thermochimica:
        call Thermochimica

        ! First check: ensure calculation succeeded
        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ! Check mole fraction 1
        expected_molfrac1 = 0.75306881663786374D0
        computed_molfrac1 = dMolFraction(1)
        rel_error_mf1 = DABS(computed_molfrac1 - expected_molfrac1) / expected_molfrac1

        call check(error, rel_error_mf1 < 1D-3, &
            "Mole fraction 1 mismatch. Expected: " // trim(adjustl(real_to_str(expected_molfrac1))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_molfrac1))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ! Check mole fraction 9
        expected_molfrac9 = 3.0917444033201544D-002
        computed_molfrac9 = dMolFraction(9)
        rel_error_mf9 = DABS(computed_molfrac9 - expected_molfrac9) / expected_molfrac9

        call check(error, rel_error_mf9 < 1D-3, &
            "Mole fraction 9 mismatch. Expected: " // trim(adjustl(real_to_str(expected_molfrac9))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_molfrac9))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ! Check Gibbs energy
        expected_gibbs = 3.06480D6
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        ! Reset Thermochimica:
        call ResetThermoAll

    end subroutine test_w_au_ar_ne_o_low_temp

    !> Helper function to convert integer to string
    function int_to_str(i) result(str)
        integer, intent(in) :: i
        character(len=20) :: str
        write(str, '(I0)') i
    end function int_to_str

    !> Helper function to convert real to string
    function real_to_str(r) result(str)
        real(8), intent(in) :: r
        character(len=40) :: str
        write(str, '(ES15.6)') r
    end function real_to_str

end module test_systems
