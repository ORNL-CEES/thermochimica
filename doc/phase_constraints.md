# Phase Fraction Constraints Design (Work-In-Progress)

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

## Implementation Phases

### Phase A: Penalty Method (fast path)
- Add constraint residuals to the functional norm:
  - g_p = dMolesPhase(p) * S_p - target_p
  - dGEMFunctionNorm += w_penalty * g_p^2
- Add penalty “force” to GEMNewton RHS rows for phase-moles:
  - B(row_p) += 2 * w_penalty * S_p * g_p
- No change to matrix size or unknown count.

### Phase B: Hard Constraints (KKT / Lagrange multipliers)
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
- Validation:
  - Sum of fractions must equal 1.0 (within tolerance).
  - Count of constrained phases <= nElements.
  - All phase names must resolve.

## Input + API
- Input scripts and run lists support:
  - phase fraction(PhaseName) = 0.25
- TCAPI supports programmatic add/clear:
  - addPhaseFractionConstraint(name, fraction)
  - clearPhaseConstraints()

## Key Code Touch Points
- New module: src/module/ModulePhaseConstraints.f90
- Parsing: src/parser/ParseInput*.f90, src/exec/RunCalculationList.F90
- Solver:
  - src/gem/InitGEMSolver.f90
  - src/gem/GEMNewton.f90
  - src/gem/CompFunctionNorm.f90
  - src/gem/GEMLineSearch.f90 (ignore lambda tail in Phase B)
- Assemblage control: src/gem/CheckPhaseAssemblage.f90 and phase add/remove/swap helpers
- Reset: src/reset/ResetThermo.f90
- TCAPI: src/api/CouplingUtilities.f90, src/Thermochimica.h, src/Thermochimica-c.C, src/Thermochimica-cxx.*

