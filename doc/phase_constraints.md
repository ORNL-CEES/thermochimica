# Phase Fraction Constraints

## Goal
Add hard phase-fraction constraints (mole-based, element-moles definition) to Thermochimica. When constraints are active, only constrained phases are allowed to be stable and their fractions must sum to 1.0. When no constraints are specified, behavior must remain unchanged.

## Definition (Mole-Based, Element-Moles)
For constrained phase p with target fraction f_p:

- Total element moles in system:
  - T = sum_{j=1..(nElements - nChargedConstraints)} dMolesElement(j)
- Element moles in phase p:
  - Solution phase:
    - dEffStoichSolnPhase(p,j) = sum_{i in p} x_i * stoich(i,j) / iParticlesPerMole(i)
    - S_p = sum_j dEffStoichSolnPhase(p,j)
    - moles of elements in phase = dMolesPhase(p) * S_p
  - Pure condensed phase:
    - S_p = dSpeciesTotalAtoms(phase)
    - moles of elements in phase = dMolesPhase(p) * S_p

Constraint equation:
  dMolesPhase(p) * S_p = f_p * T

## Hard Constraints (KKT / Lagrange Multipliers)
Expand Newton system with Lagrange multipliers for constraints:

Unknown vector:
  [ element potentials | solution phase moles | condensed phase moles | lambda ]

Block system:
  | A_gg  A_gs  A_gc   0  | |gamma|   | B_g |
  | A_sg   0     0   C_s^T| |n_s | = | B_s |
  | A_cg   0     0   C_c^T| |n_c |   | B_c |
  |  0    C_s   C_c   0  | |lambda| | b_c |

Constraint rows:
- For each constrained solution phase p:
  - C_s(r,p) = S_p
  - b_c(r) = f_p * T
- For each constrained condensed phase p:
  - C_c(r,p) = S_p
  - b_c(r) = f_p * T

## Assemblage Rules
- If constraints exist:
  - Only constrained phases are allowed to be stable.
  - The phase assemblage is fixed to those phases.
  - Add/remove/swap logic is disabled.
  - This is a complete fixed-assemblage partition. Partial constraints that leave an
    unconstrained equilibrium remainder are not currently supported.
- Validation:
  - Sum of fractions must equal 1.0 (within tolerance).
  - Each fraction must be between 0 and 1.
  - Count of constrained phases <= the number of non-charge element constraints.
  - All phase names must resolve.

## Lifetime and Reset Behavior
- `ResetThermo` preserves phase constraints for repeated calculations and sweeps.
- `ResetThermoAllPreservePhaseConstraints` rebuilds solver, parser, reinit, and CTZ
  state while retaining the requested constraints. Sweep drivers use this after a
  failed calculation.
- `ResetThermoAll` clears phase constraints along with all other state.

## Input + API
- Input scripts and run lists support:
  - phase fraction(PhaseName) = 0.25
- TCAPI supports programmatic add/clear:
  - addPhaseFractionConstraint(name, fraction)
  - clearPhaseConstraints()
- TCAPI, C, and C++ callers can query the constraint count and retrieve each solved
  constraint by index. TCAPI and C indices are one-based; C++ indices are zero-based.
  Each result contains the phase name, target and achieved fractions, normalized
  residual, and Lagrange multiplier.
- The C and C++ add functions return status 1 for an empty name and status 2 for a
  non-finite or out-of-range fraction. Phase-name resolution remains a solve-time check.
- Result queries return status 1 for an invalid index and status 2 when a successful
  solved state is unavailable. The returned C name pointer is library-owned and remains
  valid until the next phase-constraint name query.

## Output
Successful JSON output includes a `phase constraints` object containing each requested
phase's target fraction, achieved fraction, normalized residual, and Lagrange multiplier.
The achieved fraction and residual are computed by the same shared routines used by the
solver convergence checks.

## Key Code Touch Points
- New module: src/module/ModulePhaseConstraints.f90
- Parsing: src/parser/ParseInput*.f90, src/exec/RunCalculationList.F90
- Solver:
  - src/gem/InitGEMSolver.f90
  - src/gem/GEMNewton.f90
  - src/gem/CompFunctionNorm.f90
  - src/gem/GEMLineSearch.f90
- Assemblage control: src/gem/CheckPhaseAssemblage.f90 and phase add/remove/swap helpers
- Reset: src/reset/ResetThermoAll.f90
- TCAPI: src/api/CouplingUtilities.f90, src/Thermochimica.h, src/Thermochimica-c.C, src/Thermochimica-cxx.*
