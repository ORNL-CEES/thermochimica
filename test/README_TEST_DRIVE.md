# Thermochimica Test-Drive Framework

## Overview

This directory contains the new test-drive based testing infrastructure for Thermochimica. Test-drive is a lightweight, modern Fortran testing framework that provides better organization, error reporting, and test management compared to the legacy standalone test programs.

## Quick Start

### Building and Running Tests

```bash
# Build the test-drive test suite
make testdrive

# Build and run the test-drive tests
make runtests

# Build everything including both legacy and test-drive tests
make test
```

### Running Tests Directly

```bash
# After building with 'make testdrive', you can run tests directly:
./bin/test_main
```

## Test Organization

The test-drive framework organizes tests into **suites** and individual **test cases**:

```
test_main.F90                    # Main test driver
├── test_error_handling.F90      # Error handling test suite
│   ├── test_no_data_file        # Individual test case
│   └── test_mass_out_of_range   # Individual test case
└── test_systems.F90             # System validation test suite
    └── test_w_au_ar_o_system    # Individual test case
```

## Converting Legacy Tests to Test-Drive

### Pattern 1: Simple Error Validation Tests

**Legacy pattern** (e.g., TestThermo01.F90, TestThermo09.F90):
```fortran
program TestThermo01
    USE ModuleThermoIO

    ! Set up test conditions
    dTemperature = 300D0
    ! ... more setup

    ! Run calculation
    call Thermochimica

    ! Check error code
    if (INFOThermo == 6) then
        print *, 'TestThermo01: PASS'
        call EXIT(0)
    else
        print *, 'TestThermo01: FAIL'
        call EXIT(1)
    end if
end program
```

**Test-drive pattern**:
```fortran
subroutine test_no_data_file(error)
    type(error_type), allocatable, intent(out) :: error

    ! Set up test conditions
    dTemperature = 300D0
    ! ... more setup

    ! Run calculation
    call Thermochimica

    ! Check error code using test-drive's check() function
    call check(error, INFOThermo == 6, &
        "Expected error code 6, got: " // trim(adjustl(int_to_str(INFOThermo))))

    ! Clean up
    call ResetThermoAll
end subroutine
```

**Key differences:**
- No `program` wrapper needed - just a subroutine
- Takes `error` as output parameter (test-drive uses this to report results)
- Use `check()` function instead of manual pass/fail logic
- No `EXIT()` calls - test-drive handles this
- Can provide descriptive error messages

### Pattern 2: Numerical Validation Tests

**Legacy pattern** (e.g., TestThermo30.F90):
```fortran
program TestThermo30
    ! Setup and calculation
    call Thermochimica

    ! Check results with tolerance
    if (INFOThermo == 0) then
        if ((DABS(dGibbsEnergySys - (-4.620D5))/((-4.620D5))) < 1D-3) then
            print *, 'TestThermo30: PASS'
            call EXIT(0)
        else
            print *, 'TestThermo30: FAIL'
            call EXIT(1)
        end if
    end if
end program
```

**Test-drive pattern**:
```fortran
subroutine test_w_au_ar_o_system(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: expected_gibbs, computed_gibbs, relative_error

    ! Setup and calculation
    call Thermochimica

    ! First check: ensure calculation succeeded
    call check(error, INFOThermo == 0, &
        "Calculation failed with error code: " // int_to_str(INFOThermo))
    if (allocated(error)) return  ! Exit early if this check fails

    ! Second check: validate numerical results
    expected_gibbs = -4.620D5
    computed_gibbs = dGibbsEnergySys
    relative_error = DABS(computed_gibbs - expected_gibbs) / DABS(expected_gibbs)

    call check(error, relative_error < 1D-3, &
        "Gibbs energy mismatch. Expected: " // real_to_str(expected_gibbs) // &
        ", Got: " // real_to_str(computed_gibbs))

    call ResetThermoAll
end subroutine
```

**Key differences:**
- Multiple `check()` calls for different validation steps
- Can exit early if a critical check fails
- Rich error messages with actual vs expected values
- Cleaner separation of validation logic

## Step-by-Step Conversion Guide

### 1. Choose or Create a Test Suite Module

Decide which suite your test belongs to:
- **Error handling tests** → `test_error_handling.F90`
- **System calculations** → `test_systems.F90`
- **New category?** → Create a new `test_<category>.F90` file

### 2. Convert the Test Logic

```fortran
! In the appropriate test suite module (e.g., test_error_handling.F90):

subroutine test_your_new_test(error)
    type(error_type), allocatable, intent(out) :: error

    ! 1. Copy setup code from legacy test
    dTemperature = ...
    dPressure = ...

    ! 3. Copy calculation calls
    call ParseCSDataFile(...)
    call Thermochimica

    ! 4. Replace pass/fail logic with check() calls
    call check(error, <condition>, "<descriptive error message>")

    ! 5. Clean up (IMPORTANT: Use ResetThermoAll, not ResetThermo)
    call ResetThermoAll
end subroutine
```

### 3. Register the Test

Add the new test to the suite's `collect` subroutine:

```fortran
subroutine collect_error_handling(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
        new_unittest("no_data_file", test_no_data_file), &
        new_unittest("mass_out_of_range", test_mass_out_of_range), &
        new_unittest("your_new_test", test_your_new_test) &  ! Add this line
        ]
end subroutine
```

### 4. Update the Makefile (if you created a new suite module)

If you created a new test suite file, add it to the Makefile:

```makefile
# In the TEST-DRIVE TESTS section:
TESTSUITE_SRC = test_error_handling.F90 test_systems.F90 test_yournew.F90

# Add compilation rule:
$(OBJ_DIR)/test_yournew.o: $(TST_DIR)/test_yournew.F90 $(TESTDRIVE_OBJ) $(MODS_LNK) | $(OBJ_DIR)
	$(FC) -I$(OBJ_DIR) -J$(OBJ_DIR) $(FCFLAGS) -c $< -o $@
```

And register it in `test_main.F90`:

```fortran
use test_yournew, only : collect_yournew

testsuites = [ &
    new_testsuite("error_handling", collect_error_handling), &
    new_testsuite("systems", collect_systems), &
    new_testsuite("yournew", collect_yournew) &
    ]
```

### 5. Test Your Conversion

```bash
make testdrive
make runtests
```

## Important Notes

### Test Initialization

**CRITICAL:** Always initialize `dElementMass = 0D0` at the start of each test!

In standalone programs, variables start at their default values. But in test-drive, tests run sequentially in the same process. Values set by previous tests persist, so you must explicitly initialize variables.

```fortran
subroutine test_your_test(error)
    type(error_type), allocatable, intent(out) :: error

    ! REQUIRED: Initialize element masses
    dElementMass = 0D0

    ! Then set specific elements as needed
    dElementMass(8) = 1.0D0  ! Oxygen
    ! ...
end subroutine
```

### Test Cleanup

**CRITICAL:** Always use `ResetThermoAll` instead of `ResetThermo` in tests!

- `ResetThermo` - Only resets calculation state
- `ResetThermoAll` - Resets calculation state AND parser state (including deallocating parser arrays)

When running multiple tests in sequence (which test-drive does), using only `ResetThermo` will cause allocation errors when the next test tries to parse a data file, because parser variables like `nSpeciesPhaseCS` remain allocated.

```fortran
! WRONG - will cause errors in test suites
call ResetThermo

! CORRECT - full cleanup for test suites
call ResetThermoAll
```

## Helper Functions

Both `test_error_handling.F90` and `test_systems.F90` include helper functions for formatting values in error messages:

```fortran
! Convert integer to string
int_to_str(i)

! Convert real to string
real_to_str(r)
```

Use these in your error messages to provide actual values:

```fortran
call check(error, computed == expected, &
    "Expected: " // int_to_str(expected) // ", Got: " // int_to_str(computed))
```

## Test-Drive Features

### Basic Assertions

```fortran
! Simple boolean check
call check(error, condition, "Error message if false")

! Check with early exit
call check(error, condition, "Critical check")
if (allocated(error)) return
```

### Tests That Should Fail

Mark tests that are expected to fail (useful for documenting known issues):

```fortran
testsuite = [ &
    new_unittest("known_issue", test_known_issue, should_fail=.true.) &
    ]
```

## Migration Strategy

**Current Status:**
- ✅ test-drive framework integrated
- ✅ 3 example tests converted (2 error tests + 1 system test)
- ⏳ ~87 legacy tests remaining

**Recommended Approach:**
1. Convert tests incrementally by category
2. Run both legacy and test-drive tests in parallel during migration
3. Remove legacy test once its test-drive equivalent is validated
4. Continue until all tests are converted

**Priority Conversion Order:**
1. Error handling tests (simple, establish pattern)
2. Common system tests (frequently used validation cases)
3. Edge cases and regression tests
4. Weekly/specialized tests

## Benefits of Test-Drive

1. **Better organization** - Tests grouped by functionality
2. **Improved error reporting** - Descriptive messages with actual values
3. **Partial failures** - One test failure doesn't stop the entire suite
4. **Test discovery** - Automatic collection and registration
5. **Modern Fortran** - Uses modern language features and patterns
6. **Active maintenance** - Part of the Fortran-lang ecosystem

## Further Reading

- [test-drive GitHub repository](https://github.com/fortran-lang/test-drive)
- `testdrive.F90` - The framework source (includes comprehensive documentation)
- Existing test examples in this directory

## Questions?

See the example test modules (`test_error_handling.F90` and `test_systems.F90`) for working examples of different test patterns.
