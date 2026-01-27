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
            new_unittest("W_Au_Ar_Ne_O_low_temp", test_w_au_ar_ne_o_low_temp), &
            new_unittest("subm_ternary_one_coeff", test_subm_ternary_one_coeff), &
            new_unittest("subm_first_sublattice_binary", test_subm_first_sublattice_binary), &
            new_unittest("subm_first_sublattice_binary_oxide", test_subm_first_sublattice_binary_oxide), &
            new_unittest("subm_second_sublattice_binary", test_subm_second_sublattice_binary), &
            new_unittest("subm_first_sublattice_ternary", test_subm_first_sublattice_ternary), &
            new_unittest("subm_second_sublattice_ternary", test_subm_second_sublattice_ternary), &
            new_unittest("subm_all_mixing", test_subm_all_mixing), &
            new_unittest("subm_mismatch_coeff_charges_first", test_subm_mismatch_coeff_charges_first), &
            new_unittest("subm_mismatch_coeff_charges_second", test_subm_mismatch_coeff_charges_second), &
            new_unittest("subm_mismatch_coeff_charges_both", test_subm_mismatch_coeff_charges_both), &
            new_unittest("subm_worst_case", test_subm_worst_case), &
            new_unittest("noble_metals_mo_ru_2250k", test_noble_metals_mo_ru_2250k), &
            new_unittest("noble_metals_mo_ru_mix_2250k", test_noble_metals_mo_ru_mix_2250k), &
            new_unittest("noble_metals_mo_tc_liquid_2310k", test_noble_metals_mo_tc_liquid_2310k), &
            new_unittest("noble_metals_pd_ru_low_temp", test_noble_metals_pd_ru_low_temp), &
            new_unittest("noble_metals_pd_tc_1000k", test_noble_metals_pd_tc_1000k), &
            new_unittest("noble_metals_tc_ru_1234k", test_noble_metals_tc_ru_1234k), &
            new_unittest("noble_metals_tc_ru_2250k", test_noble_metals_tc_ru_2250k), &
            new_unittest("noble_metals_ternary_fccn_1973k", test_noble_metals_ternary_fccn_1973k), &
            new_unittest("noble_metals_ternary_liquid_1973k", test_noble_metals_ternary_liquid_1973k), &
            new_unittest("noble_metals_tc_pd_liquid_1900k", test_noble_metals_tc_pd_liquid_1900k), &
            new_unittest("noble_metals_tc_ru_1234k_dup", test_noble_metals_tc_ru_1234k_dup), &
            new_unittest("noble_metals_tc_ru_2250k_dup", test_noble_metals_tc_ru_2250k_dup), &
            new_unittest("noble_metals_ternary_fccn_1973k_dup", test_noble_metals_ternary_fccn_1973k_dup), &
            new_unittest("noble_metals_ternary_hcpn_973k", test_noble_metals_ternary_hcpn_973k), &
            new_unittest("noble_metals_quaternary_bccn_1800k", test_noble_metals_quaternary_bccn_1800k), &
            new_unittest("magnetic_ni_cr_fe_h_300k", test_magnetic_ni_cr_fe_h_300k), &
            new_unittest("magnetic_fe_cu_c_1400k", test_magnetic_fe_cu_c_1400k), &
            new_unittest("subq_fe_ti_v_o_2000k", test_subq_fe_ti_v_o_2000k), &
            new_unittest("subq_fe_ti_o_2000k", test_subq_fe_ti_o_2000k), &
            new_unittest("subq_fe_v_o_2000k", test_subq_fe_v_o_2000k), &
            new_unittest("subq_ti_v_o_2000k", test_subq_ti_v_o_2000k), &
            new_unittest("subl_vacancy_ti_o_1000k", test_subl_vacancy_ti_o_1000k), &
            new_unittest("subl_cl_al_2000k", test_subl_cl_al_2000k), &
            new_unittest("subl_cl_na_al_solid_1000k", test_subl_cl_na_al_solid_1000k), &
            new_unittest("ternary_miscibility_pd_ru_tc_mo_400k", test_ternary_miscibility_pd_ru_tc_mo_400k), &
            new_unittest("csi_low_pressure_673k", test_csi_low_pressure_673k), &
            new_unittest("subi_sn_o_1500k_miscibility", test_subi_sn_o_1500k_miscibility), &
            new_unittest("subi_cr_zr_o_2000k", test_subi_cr_zr_o_2000k), &
            new_unittest("subi_cs_te_700k", test_subi_cs_te_700k), &
            new_unittest("subi_ca_mn_s_2500k_miscibility", test_subi_ca_mn_s_2500k_miscibility), &
            new_unittest("subi_fe_mn_ca_s_1900k_miscibility", test_subi_fe_mn_ca_s_1900k_miscibility), &
            new_unittest("subi_fe_mn_ca_1500k", test_subi_fe_mn_ca_1500k), &
            new_unittest("subi_cs_te_2500k", test_subi_cs_te_2500k), &
            new_unittest("subi_nb_sn_o_2500k_miscibility", test_subi_nb_sn_o_2500k_miscibility), &
            new_unittest("subi_cr_sn_o_2000k", test_subi_cr_sn_o_2000k) &
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

    !> Test Nb-Zr-O-H system (ternary SUBM with one coefficient)
    !> Converted from TestThermo64.F90
    subroutine test_subm_ternary_one_coeff(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        ! Initialize element masses
        dElementMass = 0D0

        ! Specify units and file
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZIRC-test64.dat'

        ! Specify conditions
        dTemperature   = 600D0
        dPressure      = 100D0
        dElementMass(41) = 1D0      ! Nb
        dElementMass(40) = 1D0      ! Zr
        dElementMass(8)  = 1D0      ! O
        dElementMass(1)  = 0.1D0    ! H

        ! Parse and run
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        ! Check success
        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ! Check Gibbs energy
        expected_gibbs = -5.24838D+05
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subm_ternary_one_coeff

    !> Test SUBM first sublattice binary mixing
    !> Converted from TestThermo80.F90
    subroutine test_subm_first_sublattice_binary(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZrFeKClNaFOLi.dat'

        dTemperature   = 1000D0
        dPressure      = 1.0D0
        dElementMass(3)  = 0.1D0    ! Li
        dElementMass(11) = 0.4D0    ! Na
        dElementMass(9)  = 1.6D0    ! F
        dElementMass(26) = 0.3D0    ! Fe
        dElementMass(19) = 0.2D0    ! K

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -4.89601D+04
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subm_first_sublattice_binary

    !> Test SUBM first sublattice binary mixing with oxide
    !> Converted from TestThermo81.F90
    subroutine test_subm_first_sublattice_binary_oxide(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZrFeKClNaFOLi.dat'

        dTemperature   = 1000D0
        dPressure      = 1.0D0
        dElementMass(3)  = 0.1D0    ! Li
        dElementMass(11) = 0.4D0    ! Na
        dElementMass(26) = 0.3D0    ! Fe
        dElementMass(8)  = 0.8D0    ! O
        dElementMass(19) = 0.2D0    ! K

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -5.12663D+04
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subm_first_sublattice_binary_oxide

    !> Test SUBM second sublattice binary mixing
    !> Converted from TestThermo82.F90
    subroutine test_subm_second_sublattice_binary(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZrFeKClNaFOLi.dat'

        dTemperature   = 1000D0
        dPressure      = 1.0D0
        dElementMass(3)  = 1.2D0    ! Li
        dElementMass(17) = 0.5D0    ! Cl
        dElementMass(9)  = 0.3D0    ! F
        dElementMass(8)  = 0.2D0    ! O

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -4.52110D+04
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subm_second_sublattice_binary

    !> Test SUBM first sublattice ternary mixing
    !> Converted from TestThermo83.F90
    subroutine test_subm_first_sublattice_ternary(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZrFeKClNaFOLi.dat'

        dTemperature   = 1000D0
        dPressure      = 1.0D0
        dElementMass(3)  = 0.1D0    ! Li
        dElementMass(11) = 0.4D0    ! Na
        dElementMass(17) = 1.6D0    ! Cl
        dElementMass(26) = 0.3D0    ! Fe
        dElementMass(19) = 0.2D0    ! K

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -5.89226D+04
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subm_first_sublattice_ternary

    !> Test SUBM second sublattice ternary mixing
    !> Converted from TestThermo84.F90
    subroutine test_subm_second_sublattice_ternary(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZrFeKClNaFOLi.dat'

        dTemperature   = 1000D0
        dPressure      = 1.0D0
        dElementMass(17) = 0.5D0    ! Cl
        dElementMass(9)  = 0.3D0    ! F
        dElementMass(8)  = 0.2D0    ! O
        dElementMass(19) = 1.2D0    ! K

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.76944D+04
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subm_second_sublattice_ternary

    !> Test SUBM all mixing terms
    !> Converted from TestThermo85.F90
    subroutine test_subm_all_mixing(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZrFeKClNaFOLi.dat'

        dTemperature   = 1000D0
        dPressure      = 1.0D0
        dElementMass(3)  = 1D0      ! Li
        dElementMass(11) = 2D0      ! Na
        dElementMass(17) = 2D0      ! Cl
        dElementMass(9)  = 12D0     ! F
        dElementMass(26) = 3D0      ! Fe
        dElementMass(8)  = 1D0      ! O
        dElementMass(19) = 4D0      ! K

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -5.64332D+05
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subm_all_mixing

    !> Test SUBM mismatch coefficients/charges, two constituents on first sublattice
    !> Converted from TestThermo86.F90
    subroutine test_subm_mismatch_coeff_charges_first(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZrFeKClNaFOLi.dat'

        dTemperature   = 1000D0
        dPressure      = 1.0D0
        dElementMass(3)  = 2D0      ! Li
        dElementMass(8)  = 3D0      ! O
        dElementMass(40) = 1D0      ! Zr

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.01677D+05
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subm_mismatch_coeff_charges_first

    !> Test SUBM mismatch coefficients/charges, two constituents on second sublattice
    !> Converted from TestThermo87.F90
    subroutine test_subm_mismatch_coeff_charges_second(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZrFeKClNaFOLi.dat'

        dTemperature   = 1000D0
        dPressure      = 1.0D0
        dElementMass(17) = 2D0      ! Cl
        dElementMass(8)  = 1D0      ! O
        dElementMass(40) = 1D0      ! Zr

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.50268D+04
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subm_mismatch_coeff_charges_second

    !> Test SUBM mismatch coefficients/charges, two constituents on both sublattices
    !> Converted from TestThermo88.F90
    subroutine test_subm_mismatch_coeff_charges_both(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZrFeKClNaFOLi.dat'

        dTemperature   = 1000D0
        dPressure      = 1.0D0
        dElementMass(3)  = 3D0      ! Li
        dElementMass(9)  = 5D0      ! F
        dElementMass(8)  = 1D0      ! O
        dElementMass(40) = 1D0      ! Zr

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.90265D+05
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subm_mismatch_coeff_charges_both

    !> Test SUBM worst case scenario
    !> Converted from TestThermo89.F90
    subroutine test_subm_worst_case(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZrFeKClNaFOLi.dat'

        dTemperature   = 1000D0
        dPressure      = 1.0D0
        dElementMass(3)  = 1D0      ! Li
        dElementMass(11) = 2D0      ! Na
        dElementMass(17) = 2D0      ! Cl
        dElementMass(9)  = 12D0     ! F
        dElementMass(26) = 3D0      ! Fe
        dElementMass(8)  = 11D0     ! O
        dElementMass(19) = 4D0      ! K
        dElementMass(40) = 5D0      ! Zr

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -7.61628D+05
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subm_worst_case

    !-------------------------------------------------------------------------------------------------------------
    ! Noble Metals (Pd-Ru-Tc-Mo) System Tests
    ! Converted from TestThermo40-54
    ! These tests validate Gibbs energy, mole fractions, and heat capacity for the noble metals system
    !-------------------------------------------------------------------------------------------------------------

    !> Test Pd-Ru-Tc-Mo system: 2250K with 80% Mo, 20% Ru (BCCN phase)
    !> Converted from TestThermo40.F90
    subroutine test_noble_metals_mo_ru_2250k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: ru_found, mo_found

        ! Initialize
        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 2250D0
        dPressure      = 1D0
        dElementMass(42) = 0.8D0    ! Mo
        dElementMass(44) = 0.2D0    ! Ru

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        ! Check calculation succeeded
        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ! Check Gibbs energy
        expected_gibbs = -1.44373D5
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ! Check mole fractions in BCCN phase
        ru_found = .FALSE.
        mo_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'BCCN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Ru') then
                        relative_error = DABS(dMolFraction(j) - 0.2D0) / 0.2D0
                        if (relative_error < 1D-3) ru_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Mo') then
                        relative_error = DABS(dMolFraction(j) - 0.8D0) / 0.8D0
                        if (relative_error < 1D-3) mo_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, ru_found, "Ru mole fraction not found or incorrect in BCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, mo_found, "Mo mole fraction not found or incorrect in BCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ! Check heat capacity
        expected_cp = 40.2724D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_mo_ru_2250k

    !> Test Pd-Ru-Tc-Mo system: 2250K with Mo:4.3, Ru:4.5 (HCPN phase)
    !> Converted from TestThermo41.F90
    subroutine test_noble_metals_mo_ru_mix_2250k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: ru_found, mo_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 2250D0
        dPressure      = 1D0
        dElementMass(42) = 4.3D0    ! Mo
        dElementMass(44) = 4.5D0    ! Ru

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.30624D6
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ru_found = .FALSE.
        mo_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'HCPN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Ru') then
                        relative_error = DABS(dMolFraction(j) - 0.51136D0) / 0.51136D0
                        if (relative_error < 1D-3) ru_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Mo') then
                        relative_error = DABS(dMolFraction(j) - 0.48864D0) / 0.48864D0
                        if (relative_error < 1D-3) mo_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, ru_found, "Ru mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, mo_found, "Mo mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 352.351D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_mo_ru_mix_2250k

    !> Test Pd-Ru-Tc-Mo system: 2310K with Mo:0.22, Tc:0.78 (LiqN phase)
    !> Converted from TestThermo42.F90
    subroutine test_noble_metals_mo_tc_liquid_2310k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: tc_found, mo_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 2310D0
        dPressure      = 1D0
        dElementMass(42) = 0.22D0   ! Mo
        dElementMass(43) = 0.78D0   ! Tc

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.5588D5
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        tc_found = .FALSE.
        mo_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'LiqN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Tc') then
                        relative_error = DABS(dMolFraction(j) - 0.66535D0) / 0.66535D0
                        if (relative_error < 1D-3) tc_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Mo') then
                        relative_error = DABS(dMolFraction(j) - 0.33465D0) / 0.33465D0
                        if (relative_error < 1D-3) mo_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, tc_found, "Tc mole fraction not found or incorrect in LiqN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, mo_found, "Mo mole fraction not found or incorrect in LiqN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 111.198D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_mo_tc_liquid_2310k

    !> Test Pd-Ru-Tc-Mo system: 400K with Pd:0.4, Ru:0.6 (FCCN phase)
    !> Converted from TestThermo43.F90
    subroutine test_noble_metals_pd_ru_low_temp(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: ru_found, pd_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 400D0
        dPressure      = 1D0
        dElementMass(46) = 0.4D0    ! Pd
        dElementMass(44) = 0.6D0    ! Ru

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.33770D4
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ru_found = .FALSE.
        pd_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'FCCN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Ru') then
                        relative_error = DABS(dMolFraction(j) - 5.4695E-03) / 5.4695E-03
                        if (relative_error < 1D-3) ru_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Pd') then
                        relative_error = DABS(dMolFraction(j) - 0.99453D0) / 0.99453D0
                        if (relative_error < 1D-3) pd_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, ru_found, "Ru mole fraction not found or incorrect in FCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, pd_found, "Pd mole fraction not found or incorrect in FCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 25.7620D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_pd_ru_low_temp

    !> Test Pd-Ru-Tc-Mo system: 1000K with Pd:40, Tc:60 (HCPN phase)
    !> Converted from TestThermo44.F90
    subroutine test_noble_metals_pd_tc_1000k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: pd_found, tc_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 1000D0
        dPressure      = 1D0
        dElementMass(46) = 40D0     ! Pd
        dElementMass(43) = 60D0     ! Tc

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -5.04309D6
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        pd_found = .FALSE.
        tc_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'HCPN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Pd') then
                        relative_error = DABS(dMolFraction(j) - 0.4D0) / 0.4D0
                        if (relative_error < 1D-3) pd_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Tc') then
                        relative_error = DABS(dMolFraction(j) - 0.6D0) / 0.6D0
                        if (relative_error < 1D-3) tc_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, pd_found, "Pd mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, tc_found, "Tc mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 3023.31D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_pd_tc_1000k

    !> Test Pd-Ru-Tc-Mo system: 1234K with Tc:1, Ru:99 (HCPN phase)
    !> Converted from TestThermo45.F90
    subroutine test_noble_metals_tc_ru_1234k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: ru_found, tc_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 1234D0
        dPressure      = 1D0
        dElementMass(43) = 1D0      ! Tc
        dElementMass(44) = 99D0     ! Ru

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -5.64282D6
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ru_found = .FALSE.
        tc_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'HCPN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Ru') then
                        relative_error = DABS(dMolFraction(j) - 0.99D0) / 0.99D0
                        if (relative_error < 1D-3) ru_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Tc') then
                        relative_error = DABS(dMolFraction(j) - 0.01D0) / 0.01D0
                        if (relative_error < 1D-3) tc_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, ru_found, "Ru mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, tc_found, "Tc mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 2983.20D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_tc_ru_1234k

    !> Test Pd-Ru-Tc-Mo system: 2250K with Tc:0.55, Ru:0.45 (HCPN phase)
    !> Converted from TestThermo46.F90
    subroutine test_noble_metals_tc_ru_2250k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: ru_found, tc_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 2250D0
        dPressure      = 1D0
        dElementMass(43) = 0.55D0   ! Tc
        dElementMass(44) = 0.45D0   ! Ru

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.54452D5
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ru_found = .FALSE.
        tc_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'HCPN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Ru') then
                        relative_error = DABS(dMolFraction(j) - 0.45D0) / 0.45D0
                        if (relative_error < 1D-3) ru_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Tc') then
                        relative_error = DABS(dMolFraction(j) - 0.55D0) / 0.55D0
                        if (relative_error < 1D-3) tc_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, ru_found, "Ru mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, tc_found, "Tc mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 38.5394D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_tc_ru_2250k

    !> Test Pd-Ru-Tc-Mo system: 1973K ternary with Mo:0.3, Pd:0.4, Ru:0.3 (FCCN phase)
    !> Converted from TestThermo47.F90
    subroutine test_noble_metals_ternary_fccn_1973k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: ru_found, mo_found, pd_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 1973D0
        dPressure      = 1D0
        dElementMass(42) = 0.3D0    ! Mo
        dElementMass(46) = 0.4D0    ! Pd
        dElementMass(44) = 0.3D0    ! Ru

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.38528D5
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ru_found = .FALSE.
        mo_found = .FALSE.
        pd_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'FCCN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Ru') then
                        relative_error = DABS(dMolFraction(j) - 0.13193D0) / 0.13193D0
                        if (relative_error < 1D-3) ru_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Mo') then
                        relative_error = DABS(dMolFraction(j) - 0.24756D0) / 0.24756D0
                        if (relative_error < 1D-3) mo_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Pd') then
                        relative_error = DABS(dMolFraction(j) - 0.62051D0) / 0.62051D0
                        if (relative_error < 1D-3) pd_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, ru_found, "Ru mole fraction not found or incorrect in FCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, mo_found, "Mo mole fraction not found or incorrect in FCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, pd_found, "Pd mole fraction not found or incorrect in FCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 37.9399D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_ternary_fccn_1973k

    !> Test Pd-Ru-Tc-Mo system: 1973K ternary with Mo:0.1, Pd:0.3, Ru:0.6 (LiqN phase)
    !> Converted from TestThermo48.F90
    subroutine test_noble_metals_ternary_liquid_1973k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: ru_found, mo_found, pd_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 1973D0
        dPressure      = 1D0
        dElementMass(42) = 0.1D0    ! Mo
        dElementMass(46) = 0.3D0    ! Pd
        dElementMass(44) = 0.6D0    ! Ru

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.27255D5
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ru_found = .FALSE.
        mo_found = .FALSE.
        pd_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'LiqN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Ru') then
                        relative_error = DABS(dMolFraction(j) - 0.13768D0) / 0.13768D0
                        if (relative_error < 1D-3) ru_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Mo') then
                        relative_error = DABS(dMolFraction(j) - 0.12624D0) / 0.12624D0
                        if (relative_error < 1D-3) mo_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Pd') then
                        relative_error = DABS(dMolFraction(j) - 0.73608D0) / 0.73608D0
                        if (relative_error < 1D-3) pd_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, ru_found, "Ru mole fraction not found or incorrect in LiqN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, mo_found, "Mo mole fraction not found or incorrect in LiqN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, pd_found, "Pd mole fraction not found or incorrect in LiqN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 78.2758D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_ternary_liquid_1973k

    !> Test Pd-Ru-Tc-Mo system: 1900K with Tc:0.125, Pd:0.874 (LiqN phase)
    !> Converted from TestThermo49.F90
    subroutine test_noble_metals_tc_pd_liquid_1900k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: pd_found, tc_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 1900D0
        dPressure      = 1D0
        dElementMass(43) = 0.125D0  ! Tc
        dElementMass(46) = 0.874D0  ! Pd

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.28092D5
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        pd_found = .FALSE.
        tc_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'LiqN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Pd') then
                        relative_error = DABS(dMolFraction(j) - 0.90358D0) / 0.90358D0
                        if (relative_error < 1D-3) pd_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Tc') then
                        relative_error = DABS(dMolFraction(j) - 9.6420E-02) / 9.6420E-02
                        if (relative_error < 1D-3) tc_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, pd_found, "Pd mole fraction not found or incorrect in LiqN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, tc_found, "Tc mole fraction not found or incorrect in LiqN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 679.082D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_tc_pd_liquid_1900k

    !> Test Pd-Ru-Tc-Mo system: 1234K with Tc:1, Ru:99 (HCPN phase) - Duplicate of Test 45
    !> Converted from TestThermo50.F90
    subroutine test_noble_metals_tc_ru_1234k_dup(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: ru_found, tc_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 1234D0
        dPressure      = 1D0
        dElementMass(43) = 1D0      ! Tc
        dElementMass(44) = 99D0     ! Ru

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -5.64282D6
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ru_found = .FALSE.
        tc_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'HCPN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Ru') then
                        relative_error = DABS(dMolFraction(j) - 0.99D0) / 0.99D0
                        if (relative_error < 1D-3) ru_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Tc') then
                        relative_error = DABS(dMolFraction(j) - 0.01D0) / 0.01D0
                        if (relative_error < 1D-3) tc_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, ru_found, "Ru mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, tc_found, "Tc mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 2983.20D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_tc_ru_1234k_dup

    !> Test Pd-Ru-Tc-Mo system: 2250K with Tc:0.55, Ru:0.45 (HCPN phase) - Duplicate of Test 46
    !> Converted from TestThermo51.F90
    subroutine test_noble_metals_tc_ru_2250k_dup(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: ru_found, tc_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 2250D0
        dPressure      = 1D0
        dElementMass(43) = 0.55D0   ! Tc
        dElementMass(44) = 0.45D0   ! Ru

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.54452D5
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ru_found = .FALSE.
        tc_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'HCPN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Ru') then
                        relative_error = DABS(dMolFraction(j) - 0.45D0) / 0.45D0
                        if (relative_error < 1D-3) ru_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Tc') then
                        relative_error = DABS(dMolFraction(j) - 0.55D0) / 0.55D0
                        if (relative_error < 1D-3) tc_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, ru_found, "Ru mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, tc_found, "Tc mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 38.5394D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_tc_ru_2250k_dup

    !> Test Pd-Ru-Tc-Mo system: 1973K ternary with Mo:0.3, Pd:0.4, Ru:0.3 (FCCN phase) - Duplicate of Test 47
    !> Converted from TestThermo52.F90
    subroutine test_noble_metals_ternary_fccn_1973k_dup(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: ru_found, mo_found, pd_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 1973D0
        dPressure      = 1D0
        dElementMass(42) = 0.3D0    ! Mo
        dElementMass(46) = 0.4D0    ! Pd
        dElementMass(44) = 0.3D0    ! Ru

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.38528D5
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ru_found = .FALSE.
        mo_found = .FALSE.
        pd_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'FCCN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Ru') then
                        relative_error = DABS(dMolFraction(j) - 0.13193D0) / 0.13193D0
                        if (relative_error < 1D-3) ru_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Mo') then
                        relative_error = DABS(dMolFraction(j) - 0.24756D0) / 0.24756D0
                        if (relative_error < 1D-3) mo_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Pd') then
                        relative_error = DABS(dMolFraction(j) - 0.62051D0) / 0.62051D0
                        if (relative_error < 1D-3) pd_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, ru_found, "Ru mole fraction not found or incorrect in FCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, mo_found, "Mo mole fraction not found or incorrect in FCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, pd_found, "Pd mole fraction not found or incorrect in FCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 37.9399D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_ternary_fccn_1973k_dup

    !> Test Pd-Ru-Tc-Mo system: 973K ternary with Mo:0.1, Pd:0.3, Ru:0.6 (HCPN phase)
    !> Converted from TestThermo53.F90
    subroutine test_noble_metals_ternary_hcpn_973k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: ru_found, mo_found, pd_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 973D0
        dPressure      = 1D0
        dElementMass(42) = 0.1D0    ! Mo
        dElementMass(46) = 0.3D0    ! Pd
        dElementMass(44) = 0.6D0    ! Ru

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -4.90922D4
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ru_found = .FALSE.
        mo_found = .FALSE.
        pd_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'HCPN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Ru') then
                        relative_error = DABS(dMolFraction(j) - 0.80273D0) / 0.80273D0
                        if (relative_error < 1D-3) ru_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Mo') then
                        relative_error = DABS(dMolFraction(j) - 0.10757D0) / 0.10757D0
                        if (relative_error < 1D-3) mo_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Pd') then
                        relative_error = DABS(dMolFraction(j) - 8.9693E-02) / 8.9693E-02
                        if (relative_error < 1D-3) pd_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, ru_found, "Ru mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, mo_found, "Mo mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, pd_found, "Pd mole fraction not found or incorrect in HCPN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 30.1581D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_ternary_hcpn_973k

    !> Test Pd-Ru-Tc-Mo system: 1800K quaternary with Tc:0.01, Pd:0.09, Mo:0.9 (BCCN phase)
    !> Converted from TestThermo54.F90
    subroutine test_noble_metals_quaternary_bccn_1800k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp
        integer :: i, j, k
        logical :: tc_found, mo_found, pd_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'NobleMetals-Kaye.dat'

        dTemperature   = 1800D0
        dPressure      = 1D0
        dElementMass(43) = 0.01D0   ! Tc
        dElementMass(46) = 0.09D0   ! Pd
        dElementMass(42) = 0.9D0    ! Mo

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.05595D5
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        tc_found = .FALSE.
        mo_found = .FALSE.
        pd_found = .FALSE.
        do i = 1, nSolnPhases
            k = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(k) == 'BCCN') then
                do j = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
                    if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Tc') then
                        relative_error = DABS(dMolFraction(j) - 9.7175E-03) / 9.7175E-03
                        if (relative_error < 1D-3) tc_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Mo') then
                        relative_error = DABS(dMolFraction(j) - 0.93597D0) / 0.93597D0
                        if (relative_error < 1D-3) mo_found = .TRUE.
                    else if (TRIM(ADJUSTL(cSpeciesName(j))) == 'Pd') then
                        relative_error = DABS(dMolFraction(j) - 5.4310E-02) / 5.4310E-02
                        if (relative_error < 1D-3) pd_found = .TRUE.
                    end if
                end do
            end if
        end do
        call check(error, tc_found, "Tc mole fraction not found or incorrect in BCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, mo_found, "Mo mole fraction not found or incorrect in BCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, pd_found, "Pd mole fraction not found or incorrect in BCCN phase")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 40.0333D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_noble_metals_quaternary_bccn_1800k

    !-------------------------------------------------------------------------------------------------------------
    ! Magnetic Tests
    ! Converted from TestThermo55-56
    ! These tests validate magnetic solution species and excess magnetic terms
    !-------------------------------------------------------------------------------------------------------------

    !> Test Ni-Cr-Fe-H system at 300K with magnetic properties (FCC_A1 and BCC_A2 phases)
    !> Converted from TestThermo55.F90
    subroutine test_magnetic_ni_cr_fe_h_300k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: T, B, StructureFactor
        integer :: i, j, iFirst
        logical :: fcc_found, bcc_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZIRC-noSUBI.dat'

        dPressure      = 1000D0
        dTemperature   = 300D0
        dElementMass(23) = 0.1D0    ! V
        dElementMass(24) = 1D0      ! Cr
        dElementMass(26) = 2D0      ! Fe
        dElementMass(28) = 1D0      ! Ni
        dElementMass(50) = 0.1D0    ! Sn
        dElementMass(1)  = 1D0      ! H

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -4.96674D04
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ! Check magnetic properties for FCC_A1 and BCC_A2 phases
        fcc_found = .FALSE.
        bcc_found = .FALSE.
        do i = 1, nSolnPhases
            j = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(j) == 'FCC_A1') then
                call CompMagneticTemperatureMoment(j, T, B)
                iFirst = nSpeciesPhase(j-1) + 1
                StructureFactor = dCoeffGibbsMagnetic(iFirst, 3)
                T = -T * StructureFactor
                B = -B * StructureFactor
                if ((DABS((T - 97.24D0) / 97.24D0) < 1D-3) .AND. &
                    (DABS((B - 0.13862D0) / 0.13862D0) < 1D-3)) then
                    fcc_found = .TRUE.
                end if
            else if (cSolnPhaseName(j) == 'BCC_A2') then
                call CompMagneticTemperatureMoment(j, T, B)
                if ((DABS((T - 641.38D0) / 641.38D0) < 1D-3) .AND. &
                    (DABS((B - 1.5248D0) / 1.5248D0) < 1D-3)) then
                    bcc_found = .TRUE.
                end if
            end if
        end do

        call check(error, fcc_found, &
            "FCC_A1 magnetic properties not found or incorrect (T=97.24, B=0.13862)")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, bcc_found, &
            "BCC_A2 magnetic properties not found or incorrect (T=641.38, B=1.5248)")

        call ResetThermoAll
    end subroutine test_magnetic_ni_cr_fe_h_300k

    !> Test Fe-Cu-C system at 1400K with magnetic properties and SUBG mixing (FCC_A1 and Liquid phases)
    !> Converted from TestThermo56.F90
    subroutine test_magnetic_fe_cu_c_1400k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: T, B, StructureFactor, expected_val, computed_val
        integer :: i, j, iFirst
        logical :: fcc_found, liquid_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'CuFeC-Kang.dat'

        dTemperature   = 1400D0
        dPressure      = 1.0D0
        dElementMass(6)  = 1.0D0    ! C
        dElementMass(26) = 1.0D0    ! Fe
        dElementMass(29) = 1.0D0    ! Cu

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.73325D05
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        ! Check magnetic properties for FCC_A1 and mole properties for Liquid phase
        fcc_found = .FALSE.
        liquid_found = .FALSE.
        do i = 1, nSolnPhases
            j = -iAssemblage(nElements + 1 - i)
            if (cSolnPhaseName(j) == 'FCC_A1') then
                call CompMagneticTemperatureMoment(j, T, B)
                iFirst = nSpeciesPhase(j-1) + 1
                StructureFactor = dCoeffGibbsMagnetic(iFirst, 3)
                T = -T * StructureFactor
                B = -B * StructureFactor
                if ((DABS((T - 59.4D0) / 59.4D0) < 1D-3) .AND. &
                    (DABS((B - 0.62064D0) / 0.62064D0) < 1D-3)) then
                    fcc_found = .TRUE.
                end if
            else if (cSolnPhaseName(j) == 'Liquid') then
                expected_val = 2.9666D0
                computed_val = dMolesPhase(nElements + 1 - i)
                if (DABS((computed_val - expected_val) / expected_val) < 1D-3) then
                    expected_val = 4.4975D-2
                    computed_val = dMolFraction(nSpeciesPhase(j))
                    if (DABS((computed_val - expected_val) / expected_val) < 1D-3) then
                        liquid_found = .TRUE.
                    end if
                end if
            end if
        end do

        call check(error, fcc_found, &
            "FCC_A1 magnetic properties not found or incorrect (T=59.4, B=0.62064)")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, liquid_found, &
            "Liquid phase properties not found or incorrect (moles=2.9666, molfrac=0.044975)")

        call ResetThermoAll
    end subroutine test_magnetic_fe_cu_c_1400k

    !-------------------------------------------------------------------------------------------------------------
    ! SUBQ (Sub-Quadrilateral) System Tests - FeTiVO Database
    ! Converted from TestThermo57-60
    ! These tests validate SUBQ solution phases with compound constituents
    !-------------------------------------------------------------------------------------------------------------

    !> Test Fe-Ti-V-O system at 2000K - SUBQ with SlagBsoln and Hemasoln phases
    !> Converted from TestThermo57.F90
    subroutine test_subq_fe_ti_v_o_2000k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp, expected_val, computed_val
        integer :: i, j, k
        logical :: slag_found, hema_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'FeTiVO.dat'

        dPressure      = 1D0
        dTemperature   = 2000D0
        dElementMass(8)  = 2D0      ! O
        dElementMass(22) = 0.5D0    ! Ti
        dElementMass(23) = 0.5D0    ! V
        dElementMass(26) = 0.5D0    ! Fe

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.21336D06
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        slag_found = .FALSE.
        hema_found = .FALSE.
        do i = 1, nSolnPhases
            k = nElements + 1 - i
            j = -iAssemblage(k)
            if (cSolnPhaseName(j) == 'SlagBsoln') then
                expected_val = 0.36803D0
                computed_val = dMolesPhase(k)
                if (DABS((computed_val - expected_val) / expected_val) < 1D-3) then
                    expected_val = 1.0885D-2
                    computed_val = dMolFraction(nSpeciesPhase(j-1) + 1)
                    if (DABS((computed_val - expected_val) / expected_val) < 1D-3) slag_found = .TRUE.
                end if
            else if (cSolnPhaseName(j) == 'Hemasoln') then
                expected_val = 0.47842D0
                computed_val = dMolesPhase(k)
                if (DABS((computed_val - expected_val) / expected_val) < 1D-3) then
                    expected_val = 5.6813D-2
                    computed_val = dMolFraction(nSpeciesPhase(j-1) + 1)
                    if (DABS((computed_val - expected_val) / expected_val) < 1D-3) hema_found = .TRUE.
                end if
            end if
        end do

        call check(error, slag_found, "SlagBsoln phase not found or incorrect")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, hema_found, "Hemasoln phase not found or incorrect")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 186.933D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_subq_fe_ti_v_o_2000k

    !> Test Fe-Ti-O system at 2000K - SUBQ with SlagBsoln and gas_ideal phases
    !> Converted from TestThermo58.F90
    subroutine test_subq_fe_ti_o_2000k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp, expected_val, computed_val
        integer :: i, j, k
        logical :: slag_found, gas_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'FeTiVO.dat'

        dPressure      = 1D0
        dTemperature   = 2000D0
        dElementMass(8)  = 2D0      ! O
        dElementMass(22) = 0.5D0    ! Ti
        dElementMass(26) = 0.5D0    ! Fe

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.00057D06
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        slag_found = .FALSE.
        gas_found = .FALSE.
        do i = 1, nSolnPhases
            k = nElements + 1 - i
            j = -iAssemblage(k)
            if (cSolnPhaseName(j) == 'SlagBsoln') then
                expected_val = 1.1768D0
                computed_val = dMolesPhase(k)
                if (DABS((computed_val - expected_val) / expected_val) < 1D-3) then
                    expected_val = 5.1547D-2
                    computed_val = dMolFraction(nSpeciesPhase(j-1) + 1)
                    if (DABS((computed_val - expected_val) / expected_val) < 1D-3) slag_found = .TRUE.
                end if
            else if (cSolnPhaseName(j) == 'gas_ideal') then
                expected_val = 0.14544D0
                computed_val = dMolesPhase(k)
                if (DABS((computed_val - expected_val) / expected_val) < 1D-3) then
                    expected_val = 2.9406D-7
                    computed_val = dMolFraction(nSpeciesPhase(j-1) + 5)  ! Check 5th species
                    if (DABS((computed_val - expected_val) / expected_val) < 1D-3) gas_found = .TRUE.
                end if
            end if
        end do

        call check(error, slag_found, "SlagBsoln phase not found or incorrect")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, gas_found, "gas_ideal phase not found or incorrect")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 105.954D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_subq_fe_ti_o_2000k

    !> Test Fe-V-O system at 2000K - SUBQ with SlagBsoln and gas_ideal phases
    !> Converted from TestThermo59.F90
    subroutine test_subq_fe_v_o_2000k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp, expected_val, computed_val
        integer :: i, j, k
        logical :: slag_found, gas_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'FeTiVO.dat'

        dPressure      = 1D0
        dTemperature   = 2000D0
        dElementMass(8)  = 2D0      ! O
        dElementMass(23) = 0.5D0    ! V
        dElementMass(26) = 0.5D0    ! Fe

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -8.92127D05
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        slag_found = .FALSE.
        gas_found = .FALSE.
        do i = 1, nSolnPhases
            k = nElements + 1 - i
            j = -iAssemblage(k)
            if (cSolnPhaseName(j) == 'SlagBsoln') then
                expected_val = 0.80434D0
                computed_val = dMolesPhase(k)
                if (DABS((computed_val - expected_val) / expected_val) < 1D-3) then
                    expected_val = 0.17025D0
                    computed_val = dMolFraction(nSpeciesPhase(j-1) + 1)
                    if (DABS((computed_val - expected_val) / expected_val) < 1D-3) slag_found = .TRUE.
                end if
            else if (cSolnPhaseName(j) == 'gas_ideal') then
                expected_val = 0.0406D0
                computed_val = dMolesPhase(k)
                if (DABS((computed_val - expected_val) / expected_val) < 1D-3) then
                    expected_val = 0.99973D0
                    computed_val = dMolFraction(nSpeciesPhase(j-1) + 1)
                    if (DABS((computed_val - expected_val) / expected_val) < 1D-3) gas_found = .TRUE.
                end if
            end if
        end do

        call check(error, slag_found, "SlagBsoln phase not found or incorrect")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, gas_found, "gas_ideal phase not found or incorrect")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 133.351D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_subq_fe_v_o_2000k

    !> Test Ti-V-O system at 2000K - SUBQ with SlagBsoln and Hemasoln phases
    !> Converted from TestThermo60.F90
    subroutine test_subq_ti_v_o_2000k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        real(8) :: expected_cp, computed_cp, expected_val, computed_val
        integer :: i, j, k
        logical :: slag_found, hema_found

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'FeTiVO.dat'

        dPressure      = 1D0
        dTemperature   = 2000D0
        dElementMass(8)  = 2D0      ! O
        dElementMass(22) = 0.5D0    ! Ti
        dElementMass(23) = 0.5D0    ! V

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call HeatCapacity

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.09209D06
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        slag_found = .FALSE.
        hema_found = .FALSE.
        do i = 1, nSolnPhases
            k = nElements + 1 - i
            j = -iAssemblage(k)
            if (cSolnPhaseName(j) == 'SlagBsoln') then
                expected_val = 0.60817D0
                computed_val = dMolesPhase(k)
                if (DABS((computed_val - expected_val) / expected_val) < 1D-3) then
                    expected_val = 0.25897D0
                    computed_val = dMolFraction(nSpeciesPhase(j-1) + 1)
                    if (DABS((computed_val - expected_val) / expected_val) < 1D-3) slag_found = .TRUE.
                end if
            else if (cSolnPhaseName(j) == 'Hemasoln') then
                expected_val = 0.18563D0
                computed_val = dMolesPhase(k)
                if (DABS((computed_val - expected_val) / expected_val) < 1D-3) then
                    expected_val = 0.34107D0
                    computed_val = dMolFraction(nSpeciesPhase(j-1) + 1)
                    if (DABS((computed_val - expected_val) / expected_val) < 1D-3) hema_found = .TRUE.
                end if
            end if
        end do

        call check(error, slag_found, "SlagBsoln phase not found or incorrect")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        call check(error, hema_found, "Hemasoln phase not found or incorrect")
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_cp = 100.960D0
        computed_cp = dHeatCapacity
        relative_error = DABS(computed_cp - expected_cp) / DABS(expected_cp)
        call check(error, relative_error < 1D-3, &
            "Heat capacity mismatch. Expected: " // trim(adjustl(real_to_str(expected_cp))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_cp))))

        call ResetThermoAll
    end subroutine test_subq_ti_v_o_2000k

    !-------------------------------------------------------------------------------------------------------------
    ! Remaining Simple Tests (61-63, 65, 90)
    !-------------------------------------------------------------------------------------------------------------

    !> Test Ti-O with vacancy handling at 1000K
    !> Converted from TestThermo61.F90
    subroutine test_subl_vacancy_ti_o_1000k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'FeTiVO.dat'

        dTemperature   = 1000D0
        dPressure      = 1D0
        dElementMass(8)  = 2D0      ! O
        dElementMass(22) = 0.5D0    ! Ti

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -6.24557D05
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subl_vacancy_ti_o_1000k

    !> Test Cl-Al system at 2000K
    !> Converted from TestThermo62.F90
    subroutine test_subl_cl_al_2000k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ClAlNa.dat'

        dTemperature   = 2000D0
        dPressure      = 1D0
        dElementMass(17) = 2D0      ! Cl
        dElementMass(13) = 1D0      ! Al

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -9.64834D+05
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subl_cl_al_2000k

    !> Test Cl-Na-Al system with solid phase at 1000K
    !> Converted from TestThermo63.F90
    subroutine test_subl_cl_na_al_solid_1000k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ClAlNa.dat'

        dTemperature   = 1000D0
        dPressure      = 1D0
        dElementMass(17) = 3D0      ! Cl
        dElementMass(11) = 1D0      ! Na
        dElementMass(13) = 1D0      ! Al

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -1.17685D+06
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_subl_cl_na_al_solid_1000k

    !> Test Pd-Ru-Tc-Mo ternary miscibility gap at 400K
    !> Converted from TestThermo65.F90
    subroutine test_ternary_miscibility_pd_ru_tc_mo_400k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ternaryMiscibility-Kaye.dat'

        dTemperature   = 400D0
        dPressure      = 1D0
        dElementMass(42) = 1D0      ! Mo
        dElementMass(43) = 1D0      ! Tc
        dElementMass(44) = 1D0      ! Ru
        dElementMass(46) = 1D0      ! Pd

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -4.48928E04
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_ternary_miscibility_pd_ru_tc_mo_400k

    !> Test CsI at low pressure (1E-5 atm) and 673K
    !> Converted from TestThermo90.F90
    subroutine test_csi_low_pressure_673k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error

        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'CsI-Pham.dat'

        dTemperature   = 673D0
        dPressure      = 1D-5
        dElementMass(53) = 1D0      ! I
        dElementMass(55) = 1D0      ! Cs

        call ParseCSDataFile(cThermoFileName)
        call Thermochimica

        call check(error, INFOThermo == 0, &
            "Thermochimica failed with error code: " // trim(adjustl(int_to_str(INFOThermo))))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if

        expected_gibbs = -4.41869E+05
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, &
            "Gibbs energy mismatch. Expected: " // trim(adjustl(real_to_str(expected_gibbs))) // &
            ", Got: " // trim(adjustl(real_to_str(computed_gibbs))))

        call ResetThermoAll
    end subroutine test_csi_low_pressure_673k

    !-------------------------------------------------------------------------------------------------------------
    ! SUBI Ionic Liquid Tests (70-78)
    ! These tests validate ionic liquid systems with sublattice ionic model
    !-------------------------------------------------------------------------------------------------------------

    !> Test Sn-O miscibility gap at 1500K - Converted from TestThermo70.F90
    subroutine test_subi_sn_o_1500k_miscibility(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZIRC_no_liq.dat'
        dTemperature   = 1500D0
        dPressure      = 1D0
        dElementMass(8) = 0.3D0
        dElementMass(50) = 0.7D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call check(error, INFOThermo == 0, "Thermochimica failed: " // int_to_str(INFOThermo))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        expected_gibbs = -181212D0
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, "Gibbs mismatch. Expected: " // real_to_str(expected_gibbs) // ", Got: " // real_to_str(computed_gibbs))
        call ResetThermoAll
    end subroutine test_subi_sn_o_1500k_miscibility

    !> Test Cr-Zr-O ternary at 2000K - Converted from TestThermo71.F90
    subroutine test_subi_cr_zr_o_2000k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZIRC_no_liq.dat'
        dTemperature   = 2000D0
        dPressure      = 1D0
        dElementMass(8) = 0.1D0
        dElementMass(24) = 0.2D0
        dElementMass(40) = 0.7D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call check(error, INFOThermo == 0, "Thermochimica failed: " // int_to_str(INFOThermo))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        expected_gibbs = -194030D0
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, "Gibbs mismatch. Expected: " // real_to_str(expected_gibbs) // ", Got: " // real_to_str(computed_gibbs))
        call ResetThermoAll
    end subroutine test_subi_cr_zr_o_2000k

    !> Test Cs-Te binary at 700K - Converted from TestThermo72.F90
    subroutine test_subi_cs_te_700k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'CsTe-1.dat'
        dTemperature   = 700D0
        dPressure      = 1D0
        dElementMass(52) = 0.8D0
        dElementMass(55) = 0.2D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call check(error, INFOThermo == 0, "Thermochimica failed: " // int_to_str(INFOThermo))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        expected_gibbs = -90213.3D0
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, "Gibbs mismatch. Expected: " // real_to_str(expected_gibbs) // ", Got: " // real_to_str(computed_gibbs))
        call ResetThermoAll
    end subroutine test_subi_cs_te_700k

    !> Test Ca-Mn-S miscibility at 2500K - Converted from TestThermo73.F90
    subroutine test_subi_ca_mn_s_2500k_miscibility(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'CaMnS.dat'
        dTemperature   = 2500D0
        dPressure      = 1D0
        dElementMass(20) = 1D0
        dElementMass(25) = 1D0
        dElementMass(16) = 1D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call check(error, INFOThermo == 0, "Thermochimica failed: " // int_to_str(INFOThermo))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        expected_gibbs = -910619D0
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, "Gibbs mismatch. Expected: " // real_to_str(expected_gibbs) // ", Got: " // real_to_str(computed_gibbs))
        call ResetThermoAll
    end subroutine test_subi_ca_mn_s_2500k_miscibility

    !> Test Fe-Mn-Ca-S quaternary at 1900K - Converted from TestThermo74.F90
    subroutine test_subi_fe_mn_ca_s_1900k_miscibility(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'FeMnCaS-1.dat'
        dTemperature   = 1900D0
        dPressure      = 1D0
        dElementMass(26) = 2D0
        dElementMass(25) = 4D0
        dElementMass(16) = 3D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call check(error, INFOThermo == 0, "Thermochimica failed: " // int_to_str(INFOThermo))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        expected_gibbs = -1944880D0
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, "Gibbs mismatch. Expected: " // real_to_str(expected_gibbs) // ", Got: " // real_to_str(computed_gibbs))
        call ResetThermoAll
    end subroutine test_subi_fe_mn_ca_s_1900k_miscibility

    !> Test Fe-Mn-Ca ternary at 1500K - Converted from TestThermo75.F90
    subroutine test_subi_fe_mn_ca_1500k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'FeMnCaS-2.dat'
        dTemperature   = 1500D0
        dPressure      = 1D0
        dElementMass(26) = 0.4D0
        dElementMass(25) = 0.2D0
        dElementMass(20) = 0.5D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call check(error, INFOThermo == 0, "Thermochimica failed: " // int_to_str(INFOThermo))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        expected_gibbs = -116594D0
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, "Gibbs mismatch. Expected: " // real_to_str(expected_gibbs) // ", Got: " // real_to_str(computed_gibbs))
        call ResetThermoAll
    end subroutine test_subi_fe_mn_ca_1500k

    !> Test Cs-Te high temp at 2500K - Converted from TestThermo76.F90
    subroutine test_subi_cs_te_2500k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'CsTe-2.dat'
        dTemperature   = 2500D0
        dPressure      = 1D0
        dElementMass(52) = 0.5D0
        dElementMass(55) = 1.0D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call check(error, INFOThermo == 0, "Thermochimica failed: " // int_to_str(INFOThermo))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        expected_gibbs = -568905D0
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, "Gibbs mismatch. Expected: " // real_to_str(expected_gibbs) // ", Got: " // real_to_str(computed_gibbs))
        call ResetThermoAll
    end subroutine test_subi_cs_te_2500k

    !> Test Nb-Sn-O miscibility at 2500K - Converted from TestThermo77.F90
    subroutine test_subi_nb_sn_o_2500k_miscibility(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZIRC_no_liq_mod1.dat'
        dTemperature   = 2500D0
        dPressure      = 1D0
        dElementMass(41) = 1D0
        dElementMass(8) = 0.3D0
        dElementMass(50) = 0.7D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call check(error, INFOThermo == 0, "Thermochimica failed: " // int_to_str(INFOThermo))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        expected_gibbs = -535982D0
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, "Gibbs mismatch. Expected: " // real_to_str(expected_gibbs) // ", Got: " // real_to_str(computed_gibbs))
        call ResetThermoAll
    end subroutine test_subi_nb_sn_o_2500k_miscibility

    !> Test Cr-Sn-O ternary at 2000K - Converted from TestThermo78.F90
    subroutine test_subi_cr_sn_o_2000k(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: expected_gibbs, computed_gibbs, relative_error
        dElementMass = 0D0
        cInputUnitTemperature = 'K'
        cInputUnitPressure    = 'atm'
        cInputUnitMass        = 'moles'
        cThermoFileName       = DATA_DIRECTORY // 'ZIRC_no_liq_mod2.dat'
        dTemperature   = 2000D0
        dPressure      = 1D0
        dElementMass(24) = 0.5D0
        dElementMass(50) = 0.8D0
        dElementMass(8) = 0.6D0
        call ParseCSDataFile(cThermoFileName)
        call Thermochimica
        call check(error, INFOThermo == 0, "Thermochimica failed: " // int_to_str(INFOThermo))
        if (allocated(error)) then
            call ResetThermoAll
            return
        end if
        expected_gibbs = -1266220D0
        computed_gibbs = dGibbsEnergySys
        relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)
        call check(error, relative_error < 1D-3, "Gibbs mismatch. Expected: " // real_to_str(expected_gibbs) // ", Got: " // real_to_str(computed_gibbs))
        call ResetThermoAll
    end subroutine test_subi_cr_sn_o_2000k

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
