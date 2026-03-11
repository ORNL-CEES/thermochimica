# Repository Guidelines

## Project Structure & Module Organization
- `src/`: core Fortran library (solver, setup, parser, postprocess, GEM, API). `src/exec/` contains executable drivers; `src/module/` holds `Module*.f90`.
- `test/`: Fortran driver tests; daily tests live in `test/daily/`.
- `data/`: thermodynamic databases; `inputs/`: input-script examples (`.ti`); `outputs/`: runtime outputs (for example `outputs/thermoout.json`).
- `scripts/`: GUI launchers; `python/`: GUI support and `python/requirements.txt`.
- `doc/`: Doxygen and GUI docs; `cmake/`: CMake configuration; root `Makefile` drives standard builds.

## Build, Test, and Development Commands
- `make`: compile libraries and executables into `lib/` and `bin/`.
- `make test`: build all standard tests; `make dailytest`: build daily tests.
- `./run_tests`: run the compiled `bin/TestThermo*` test suite.
- `./bin/InputScriptMode inputs/demo.ti`: run an input-script workflow.
- `make doc`: build Doxygen HTML/LaTeX docs under `doc/`.

## Coding Style & Naming Conventions
- Primary language is Fortran free-form (`.f90`); preprocessed sources use `.F90`. C/C++ bindings live in `src/*.C`.
- Keep naming consistent with existing files: modules as `Module*.f90`, test drivers as `TestThermoNN.F90` or descriptive names in `test/`.
- Follow local formatting/indentation in the file you are editing; avoid reformatting unrelated lines.

## Testing Guidelines
- Add or update Fortran driver tests in `test/` (or `test/daily/` when appropriate).
- Ensure new tests build via `make test` and run via `./run_tests` (or execute the generated `bin/<test>` directly).
- No explicit coverage tool is wired in; focus on deterministic, scriptable test cases.

## Commit & Pull Request Guidelines
- Commit messages are short, imperative, and specific (for example “Fix sublattice indexing”, “Add output path option”).
- PRs should describe the change, list test commands run, and update `README.md`/docs when behavior or user workflows change.
- Include screenshots when modifying GUI behavior.
