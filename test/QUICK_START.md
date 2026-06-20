# Test-Drive Quick Start Guide

## Run Tests

```bash
# Build and run test-drive tests
make runtests

# Or separately:
make testdrive        # Build
./bin/test_main       # Run
```

## Add a New Test (3 Simple Steps)

### 1. Add test function to appropriate suite module

**For error tests:** Edit `test/test_error_handling.F90`
**For system tests:** Edit `test/test_systems.F90`

```fortran
subroutine test_your_test_name(error)
    type(error_type), allocatable, intent(out) :: error

    ! Initialize (required!)
    dElementMass = 0D0

    ! Setup
    dTemperature = 300D0
    dPressure = 1D0
    cThermoFileName = DATA_DIRECTORY // 'yourdata.dat'

    ! Run
    call ParseCSDataFile(cThermoFileName)
    call Thermochimica

    ! Check
    call check(error, INFOThermo == 0, &
        "Calculation failed with error: " // int_to_str(INFOThermo))

    ! Cleanup (required!)
    call ResetThermoAll
end subroutine
```

### 2. Register test in collect subroutine

```fortran
subroutine collect_error_handling(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
        new_unittest("existing_test", test_existing), &
        new_unittest("your_test_name", test_your_test_name) &  ! Add this
        ]
end subroutine
```

### 3. Rebuild and run

```bash
make runtests
```

Done! Your test is now part of the suite.

## Common Patterns

### Check error code
```fortran
call check(error, INFOThermo == expected_code, &
    "Expected code " // int_to_str(expected_code) // &
    ", got: " // int_to_str(INFOThermo))
```

### Check numerical value
```fortran
relative_error = DABS(computed - expected) / DABS(expected)
call check(error, relative_error < tolerance, &
    "Value mismatch. Expected: " // real_to_str(expected) // &
    ", Got: " // real_to_str(computed))
```

### Exit early on failure
```fortran
call check(error, condition1, "First check failed")
if (allocated(error)) return  ! Stop if this critical check failed

call check(error, condition2, "Second check failed")
```

## CRITICAL Rules

1. ✅ **Always** initialize `dElementMass = 0D0` at test start
2. ✅ **Always** use `ResetThermoAll` for cleanup (not `ResetThermo`)
3. ✅ Provide descriptive error messages
4. ✅ Test after each conversion

## Get Help

- Full documentation: `test/README_TEST_DRIVE.md`
- Migration status: `test/MIGRATION_STATUS.md`
- Working examples: `test/test_error_handling.F90`, `test/test_systems.F90`
